##
# File: NmrVrptUtility.py
# Date: 19-Apr-2023
#
# Updates:
# 19-Jul-2023  M. Yokochi - add trustPdbxAuthAtomName() for OneDep validation package (DAOTHER-8705)
# 19-Jul-2023  M. Yokochi - fix distance/dihedral angle/RDC averaging when lower/upper bounds are different in a restraint
#                           (DAOTER-8705)
# 18-Dec-2023  M. Yokochi - retrieve non-leaving hydrogens independent of MolProbity (DAOTHER-8945)
# 20-Dec-2023  M. Yokochi - add support for case 'Member_logic_code' value equals 'AND'
# 21-Feb-2024  M. Yokochi - add support for discontinuous model_id (NMR restraint remediation, 2n6j)
# 28-Feb-2024  M. Yokochi - collect atom_ids dictionary for both auth_atom_id and pdbx_auth_atom_name tags
#                           to prevent MISSING ATOM IN MODEL KeyError in restraintsanalysis.py (DAOTHER-9200)
# 22-May-2026  M. Yokochi - transplant BMRB chemical shift analysis except for PANAV support,
#                           add 'nmr-chemical-shift-validation' workflow operation (DAOTHER-9785)
# 16-Jun-2025  M. Yokochi - enable to set working directory and cache file directory (DAOTHER-9785)
# 24-Jul-2026  M. Yokochi - map comma-separated Auth_asym_IDs for calculation of chemical shift completeness (DAOTHER-10898)
# 31-Jul-2026  M. Yokochi - fix unexpected missing of dihedral angle restraint validation
# 03-Aug-2026  M. Yokochi - add RDC restraint analysis (DAOTHER-9785, v1.3.0)
# 18-Aug-2026  M. Yokochi - implement Monte Carlo simulation to estimate uncertainty of calculated RDC values
#                           (DAOTHER-9785, 10893, v1.3.1)
# 09-Sep-2026  M. Yokochi - calculate number of chemical shifts mapped to unmodeled residues as unmapped warning
#                           (DAOTHER-9785, 10987, v1.3.2)
# 16-Sep-2026  M. Yokochi - allow to run chemical shift analysis without coordinates (DAOTHER-9785, v1.3.3)
# 30-Sep-2026  M. Yokochi - import NEXT_STAR_FILE_PATH_KEY in standalone mode too, and remove the temporary CIF file of
#                           a pynmrstar object only when its path was not given (DAOTHER-9785)
# 03-Oct-2026  M. Yokochi - detect duplicated chemical shifts within each list, which the check against the dictionary
#                           of all lists never did, and look up the shifts of each residue for the random coil index
#                           instead of scanning all shifts of the list per residue (DAOTHER-7829, 9785)
# 03-Oct-2026  M. Yokochi - keep the workflow tasks as functions instead of bound methods, so that an instance is not
#                           a reference cycle and is freed as soon as it is dropped (DAOTHER-7829, 9785)
# 13-Oct-2026  M. Yokochi - hand the data blocks parsed from the NMR-STAR file to CifReader instead of writing them to
#                           a CIF file and parsing it back, which cost about 100 s on 6cc9; write the CIF file only when
#                           its path is given (DAOTHER-7829, 9785)
##
""" Wrapper class for NMR chemical shifts and restraints analysis.
    @author: Masashi Yokochi
    @note: This class is alternative implementation of wwpdb.apps.validation.src.ChemicalShiftsValidation.BMRBChemicalShiftsAnalysis
    @note: This class is alternative implementation of wwpdb.apps.validation.src.RestraintsValidation.BMRBRestraintsAnalysis
"""
__docformat__ = "restructuredtext en"
__author__ = "Masashi Yokochi, Kumaran Baskaran"
__email__ = "yokochi@protein.osaka-u.ac.jp, baskaran@uchc.edu"
__license__ = "Apache License 2.0"
__version__ = "v1.3.3"

import collections
import copy
import gc
import gzip
import math
import os
import pickle
import sys
import tempfile
import time
import warnings
from operator import itemgetter
from typing import Any, IO, List, Optional, Tuple

from mmcif.io.IoAdapterPy import IoAdapterPy

import numpy

try:
    from wwpdb.utils.nmr.NmrDpConstant import (MODEL_FILE_PATH_KEY,
                                               NMR_CIF_FILE_PATH_KEY,
                                               NMR_STR_FILE_PATH_KEY,
                                               NEXT_STAR_FILE_PATH_KEY,
                                               RESULT_PKL_FILE_PATH_KEY,
                                               REPORT_FILE_PATH_KEY,
                                               CIF_READER_OBJ_KEY,
                                               NMR_CIF_READER_OBJ_KEY,
                                               PYNMRSTAR_OBJ_KEY,
                                               VRPT_INPUT_PARAM_KEYS,
                                               VRPT_INPUT_FILE_KEYS,
                                               VRPT_OUTPUT_FILE_KEYS,
                                               VRPT_WORKFLOW_OPS,
                                               SUB_DIR_NAME_FOR_CACHE,
                                               LARGE_ASYM_ID,
                                               LEN_MAJOR_ASYM_ID,
                                               EMPTY_VALUE,
                                               STD_MON_DICT,
                                               PROTON_BEGIN_CODE,
                                               AMINO_PROTON_CODE,
                                               REPRESENTATIVE_MODEL_ID,
                                               REPRESENTATIVE_ALT_ID,
                                               REPRESENTATIVE_ASYM_ID,
                                               CARTN_DATA_ITEMS,
                                               THRESHOLD_FOR_CIRCULAR_SHIFT,
                                               DIST_ERROR_MAX,
                                               ANGLE_ERROR_MAX,
                                               RDC_ERROR_MAX,
                                               ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS,
                                               ALLOWED_AMBIGUITY_CODES,
                                               PARAMAGNETIC_ELEMENTS,
                                               FERROMAGNETIC_ELEMENTS,
                                               CUTOFF_AROMATIC,
                                               CUTOFF_PARAMAGNETIC,
                                               VICINITY_AROMATIC,
                                               VICINITY_PARAMAGNETIC,
                                               MAGIC_ANGLE,
                                               GYROMAGNETIC_RATIOS,
                                               PERMEABILITY_0,
                                               PLANCK_CONSTANT,
                                               REDUCED_PLANCK_CONSTANT,
                                               text_md5)
    from wwpdb.utils.nmr.ChemCompUtil import ChemCompUtil
    from wwpdb.utils.nmr.BmrbChemShiftStat import BmrbChemShiftStat
    from wwpdb.utils.nmr.NmrDpReport import NmrDpReport
    from wwpdb.utils.nmr.io.CifReader import (CifReader,
                                              SharingPdbxReader,
                                              calculate_uninstanced_coord,
                                              to_np_array)
    from wwpdb.utils.nmr.mr.ParserListenerUtil import (coordAssemblyChecker,
                                                       getDistConstraintType,
                                                       getRdcCode)
    from wwpdb.utils.nmr.rci.RCI import RCI
except ImportError:
    from nmr.NmrDpConstant import (MODEL_FILE_PATH_KEY,
                                   NMR_CIF_FILE_PATH_KEY,
                                   NMR_STR_FILE_PATH_KEY,
                                   NEXT_STAR_FILE_PATH_KEY,
                                   RESULT_PKL_FILE_PATH_KEY,
                                   REPORT_FILE_PATH_KEY,
                                   CIF_READER_OBJ_KEY,
                                   NMR_CIF_READER_OBJ_KEY,
                                   PYNMRSTAR_OBJ_KEY,
                                   VRPT_INPUT_PARAM_KEYS,
                                   VRPT_INPUT_FILE_KEYS,
                                   VRPT_OUTPUT_FILE_KEYS,
                                   VRPT_WORKFLOW_OPS,
                                   SUB_DIR_NAME_FOR_CACHE,
                                   LARGE_ASYM_ID,
                                   LEN_MAJOR_ASYM_ID,
                                   EMPTY_VALUE,
                                   STD_MON_DICT,
                                   PROTON_BEGIN_CODE,
                                   AMINO_PROTON_CODE,
                                   REPRESENTATIVE_MODEL_ID,
                                   REPRESENTATIVE_ALT_ID,
                                   REPRESENTATIVE_ASYM_ID,
                                   CARTN_DATA_ITEMS,
                                   THRESHOLD_FOR_CIRCULAR_SHIFT,
                                   DIST_ERROR_MAX,
                                   ANGLE_ERROR_MAX,
                                   RDC_ERROR_MAX,
                                   ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS,
                                   ALLOWED_AMBIGUITY_CODES,
                                   PARAMAGNETIC_ELEMENTS,
                                   FERROMAGNETIC_ELEMENTS,
                                   CUTOFF_AROMATIC,
                                   CUTOFF_PARAMAGNETIC,
                                   VICINITY_AROMATIC,
                                   VICINITY_PARAMAGNETIC,
                                   MAGIC_ANGLE,
                                   GYROMAGNETIC_RATIOS,
                                   PERMEABILITY_0,
                                   PLANCK_CONSTANT,
                                   REDUCED_PLANCK_CONSTANT,
                                   text_md5)
    from nmr.ChemCompUtil import ChemCompUtil
    from nmr.BmrbChemShiftStat import BmrbChemShiftStat
    from nmr.NmrDpReport import NmrDpReport
    from nmr.io.CifReader import (CifReader,
                                  SharingPdbxReader,
                                  calculate_uninstanced_coord,
                                  to_np_array)
    from nmr.mr.ParserListenerUtil import (coordAssemblyChecker,
                                           getDistConstraintType,
                                           getRdcCode)
    from nmr.rci.RCI import RCI


# throw RuntimeWarning as error for bug tracking, any runtimewarning should be handled
warnings.filterwarnings('error', category=RuntimeWarning, module='NmrVrptUtility')

NMR_VTF_DIST_VIOL_CUTOFF = 0.1
NMR_VTF_DIHED_VIOL_CUTOFF = 1.0
NMR_VTF_RDC_VIOL_CUTOFF = 1.0  # to be decided

NMR_VTF_DIST_ERR_BINS = (0.1, 0.2, 0.5)
NMR_VTF_DIHED_ERR_BINS = (1.0, 10.0, 20.0)
NMR_VTF_RDC_ERR_BINS = (1.0, 2.0, 5.0)  # to be decided


# effective Monte Carlo simulation cycles for estimating uncertainty of calculated RDC values
# note that the real cycle will be scaled by the effective models
RDC_EFF_MC_CYCLES = 1000
RDC_MIN_MC_CYCLES = RDC_EFF_MC_CYCLES * 0.9


def uncompress_gzip_file(inPath: str, outPath: str) -> None:
    """ Uncompress a given gzip file.
    """

    with gzip.open(inPath, mode='rt') as ifh, \
            open(outPath, 'w', encoding='utf-8') as ofh:
        for line in ifh:
            ofh.write(line)


def compress_as_gzip_file(inPath: str, outPath: str) -> None:
    """ Compress a given file as a gzip file.
    """

    with open(inPath, mode='r', encoding='utf-8') as ifh, \
            gzip.open(outPath, 'wt') as ofh:
        for line in ifh:
            ofh.write(line)


def load_from_pickle(file_name: str, default: Any = None) -> Any:
    """ Load object from pickle file.
    """

    if os.path.exists(file_name):

        try:

            with open(file_name, 'rb') as ifh:
                obj = pickle.load(ifh)
                return obj if obj is not None else default

        except IOError:
            try:
                os.remove(file_name)
            except OSError:
                pass

    return default


def write_as_pickle(obj: Any, file_name: str) -> None:
    """ Write a given object as pickle file.
    """

    if obj is not None:
        with open(file_name, 'wb') as ofh:
            pickle.dump(obj, ofh)


def get_temp_path(src_path: str = None, suffix: str = '~'):
    """ Return tempfile path based on given hints.
    """

    if src_path is None:
        return os.path.join(tempfile.gettempdir(), f"{next(tempfile._get_candidate_names())}{suffix}")  # noqa: E501, pylint: disable=protected-access,line-too-long
    return os.path.join(tempfile.gettempdir(), f"{os.path.basename(src_path)}{suffix}")


def distance(p0: list, p1: list) -> float:
    """ Return distance between two points.
    """

    # math.hypot on the 3 components is ~2x faster than numpy.linalg.norm for a
    # single 3-vector (this is the most-called geometry primitive) and numerically
    # identical. p0/p1 are numpy arrays (their difference is a length-3 array).
    d = p0 - p1

    try:
        return math.hypot(d[0], d[1], d[2])
    except TypeError:
        # Backward compatibility with Python 3.7 and earlier
        return numpy.linalg.norm(d)


def to_unit_vector(a: list) -> list:
    """ Return unit vector of a given vector.
    """

    try:
        return a / math.hypot(a[0], a[1], a[2])
    except TypeError:
        # Backward compatibility with Python 3.7 and earlier
        return a / numpy.linalg.norm(a)


def dist_inv_6_summed(r_list: List[float]) -> float:
    """ Return r^−6-summed distance for a given list of distances for ambiguous restraints as recommended by NMR VTF.
        Reference:
          Calculation of Protein Structures with Ambiguous Distance Restraints. Automated Assignment of
          Ambiguous NOE Crosspeaks and Disulphide Connectivities.
          Michael Nilges.
          J. Mol. Biol. (1995) 245, 645–660.
          DOI: 10.1006/jmbi.1994.0053
        @author: Kumaran Baskaran
        @see: wwpdb.apps.validation.src.RestraintsValidation.BMRBRestraintsAnalysis.r6sum
    """

    if len(r_list) == 1:
        return r_list[0]

    return sum(r ** (-6.0) for r in r_list) ** (-1.0 / 6.0)


def dist_target_values(target_value: Optional[float], target_value_uncertainty: Optional[float],
                       lower_limit: Optional[float], upper_limit: Optional[float],
                       lower_linear_limit: Optional[float], upper_linear_limit: Optional[float]
                       ) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    """ Return estimated distance target value, lower_bound, upper_bound.
        @author: Masashi Yokochi
    """

    lower_bound, upper_bound = lower_limit, upper_limit

    if lower_bound is None and upper_bound is None and lower_linear_limit is None and upper_linear_limit is None:
        if target_value is None:
            return None, lower_bound, upper_bound
        if target_value is not None:
            if target_value_uncertainty is not None:
                lower_bound = target_value - target_value_uncertainty
                upper_bound = target_value + target_value_uncertainty
            else:
                lower_bound = upper_bound = target_value

    if lower_bound is None and lower_linear_limit is not None:
        lower_bound = lower_linear_limit
    if upper_bound is None and upper_linear_limit is not None:
        upper_bound = upper_linear_limit

    if lower_bound is not None:
        lower_bound = max(lower_bound, 0.0)
    if upper_bound is not None:
        upper_bound = max(upper_bound, 0.0)

    if target_value is None:
        if None not in (lower_bound, upper_bound):
            target_value = (lower_bound + upper_bound) / 2.0
        elif lower_bound is not None:
            target_value = lower_bound
        elif upper_bound is not None:
            target_value = upper_bound

    return target_value, lower_bound, upper_bound


def dist_error(lower_bound: Optional[float], upper_bound: Optional[float], dist: Optional[float]) -> float:
    """ Return distance outlier for given lower_bound and upper_bound.
        @author: Masashi Yokochi
    """

    error = 0.0

    try:

        if None not in (lower_bound, upper_bound):
            if lower_bound <= dist <= upper_bound:
                pass
            elif dist > upper_bound:
                error = abs(dist - upper_bound)
            else:
                error = abs(dist - lower_bound)

        elif upper_bound is not None:
            if dist <= upper_bound:
                pass
            elif dist > upper_bound:
                error = abs(dist - upper_bound)

        elif lower_bound is not None:
            if lower_bound <= dist:
                pass
            else:
                error = abs(dist - lower_bound)

    except TypeError:
        pass

    return error


def dihedral_angle(p0: list, p1: list, p2: list, p3: list) -> float:
    """ Return dihedral angle from a series of four points.
    """

    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 = to_unit_vector(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - numpy.dot(b0, b1) * b1
    w = b2 - numpy.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = numpy.dot(v, w)
    y = numpy.dot(numpy.cross(b1, v), w)

    return numpy.degrees(numpy.arctan2(y, x))


def dihedral_angles(p0: numpy.ndarray, p1: numpy.ndarray, p2: numpy.ndarray, p3: numpy.ndarray
                    ) -> numpy.ndarray:
    """ Vectorized dihedral_angle(): p0..p3 are (M, 3) arrays (one row per model);
        returns an (M,) array of angles in degrees. Same Praxeolitic formula as
        dihedral_angle(), applied row-wise. (dihedral_angle() itself is kept for the
        many scalar callers across the package.)
    """

    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2

    b1 = b1 / numpy.linalg.norm(b1, axis=1, keepdims=True)

    v = b0 - numpy.sum(b0 * b1, axis=1, keepdims=True) * b1
    w = b2 - numpy.sum(b2 * b1, axis=1, keepdims=True) * b1

    x = numpy.sum(v * w, axis=1)
    y = numpy.sum(numpy.cross(b1, v) * w, axis=1)

    return numpy.degrees(numpy.arctan2(y, x))


def angle_target_values(target_value: Optional[float], target_value_uncertainty: Optional[float],
                        lower_limit: Optional[float], upper_limit: Optional[float],
                        lower_linear_limit: Optional[float], upper_linear_limit: Optional[float]
                        ) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    """ Return estimated angle target value, lower_bound, upper_bound.
        @author: Masashi Yokochi
        @note: support for the case target_value is not set, but upper/lower_limit and upper/lower_linear_limit are set
                       (i.e. AMBER restraint format, decide target_value by testing anti-clockwise and clockwise mean),
               support for the case target_value is not set, but upper/lower_limit are set
                       (i.e. CYANA restraint format, decide target_value by comparing lower_limit and upper_limit values),
               support for the case upper/lower_linear_limit are set, but missing upper/lower_limit
                       (i.e. XPLOR-NIH/CNS exponent parameter (ed) equals 1)
    """

    lower_bound, upper_bound = lower_limit, upper_limit

    if lower_bound is None and upper_bound is None and lower_linear_limit is None and upper_linear_limit is None:

        if target_value is None:
            return None, lower_bound, upper_bound

        if target_value_uncertainty is not None:
            lower_bound = target_value - target_value_uncertainty
            upper_bound = target_value + target_value_uncertainty
        else:
            lower_bound = upper_bound = target_value

    if lower_bound is None and lower_linear_limit is not None:
        lower_bound = lower_linear_limit
    if upper_bound is None and upper_linear_limit is not None:
        upper_bound = upper_linear_limit

    if target_value is None:  # target values are not always filled (e.g. AMBER/CYANA dihedral angle restraints)
        has_valid_lower_linear_limit =\
            lower_bound is not None and lower_linear_limit is not None and lower_bound != lower_linear_limit
        has_valid_upper_linear_limit =\
            upper_bound is not None and upper_linear_limit is not None and upper_bound != upper_linear_limit

        try:

            target_value_aclock = (lower_bound + upper_bound) / 2.0
            target_value_clock = target_value_aclock + 180.0
            if target_value_clock >= 360.0:
                target_value_clock -= 360.0

            # decide target value from upper/lower_limit and upper/lower_linear_limit (AMBER)
            if has_valid_lower_linear_limit or has_valid_upper_linear_limit:
                target_value_vote_aclock = target_value_vote_clock = 0

                if has_valid_lower_linear_limit:
                    if angle_diff(lower_bound, target_value_aclock) < angle_diff(lower_linear_limit, target_value_aclock):
                        target_value_vote_aclock += 1
                    elif angle_diff(lower_bound, target_value_clock) < angle_diff(lower_linear_limit, target_value_clock):
                        target_value_vote_clock += 1
                if has_valid_upper_linear_limit:
                    if angle_diff(upper_bound, target_value_aclock) < angle_diff(upper_linear_limit, target_value_aclock):
                        target_value_vote_aclock += 1
                    elif angle_diff(upper_bound, target_value_clock) < angle_diff(upper_linear_limit, target_value_clock):
                        target_value_vote_clock += 1

                if target_value_vote_aclock + target_value_vote_clock == 0\
                   or target_value_vote_aclock * target_value_vote_clock != 0:
                    if angle_diff(upper_bound, target_value_aclock) > angle_diff(upper_bound, target_value_clock)\
                       and angle_diff(lower_bound, target_value_aclock) > angle_diff(lower_bound, target_value_clock):
                        return target_value_aclock, lower_bound, upper_bound
                    if angle_diff(upper_bound, target_value_clock) > angle_diff(upper_bound, target_value_aclock)\
                       and angle_diff(lower_bound, target_value_clock) > angle_diff(lower_bound, target_value_aclock):
                        return target_value_clock, lower_bound, upper_bound
                    return None, lower_bound, upper_bound

                target_value = target_value_aclock if target_value_vote_aclock > target_value_vote_clock else target_value_clock

            else:  # estimate target value by comparing lower_limit and upper_limit value (CYANA)
                target_value = target_value_aclock if lower_bound <= upper_bound else target_value_clock

        except TypeError:
            return None, lower_bound, upper_bound

    return target_value, lower_bound, upper_bound


def angle_diff(x: float, y: float) -> float:
    """ Return normalized angular difference.
        @author: Kumaran Baskaran
        @see: wwpdb.apps.validation.src.RestraintsValidation.BMRBRestraintsAnalysis.angle_diff.ac
    """

    if x < 0.0:
        x += 360.0
    if y < 0.0:
        y += 360.0

    a, b = sorted([x, y])

    d = b - a
    if d > 180.0:
        d = 360.0 - d

    return d


def angle_error(lower_bound: Optional[float], upper_bound: Optional[float], target_value: Optional[float], angle: float
                ) -> float:
    """ Return angle outlier for given lower_bound, upper_bound, and target_value.
        @author: Kumaran Baskaran
        @see: wwpdb.apps.validation.src.RestraintsValidation.BMRBRestraintsAnalysis.angle_diff
    """

    def check_angle_range_overlap(x, y, c, g, t: float = 0.5):
        """ Return whether angular range formed by (x, c) and (c, y) matches to a given range (g) with tolerance (t).
            @author: Kumaran Baskaran
            @see: wwpdb.apps.validation.src.RestraintsValidation.BMRBRestraintsAnalysis.angle_diff.check_ac
        """

        l = angle_diff(x, c)  # noqa: E741
        r = angle_diff(y, c)

        return abs(g - (l + r)) < t

    try:

        ld = angle_diff(lower_bound, target_value)
        rd = angle_diff(upper_bound, target_value)

        if check_angle_range_overlap(lower_bound, target_value, angle, ld, 0.5)\
           or check_angle_range_overlap(upper_bound, target_value, angle, rd, 0.5):
            return 0.0

    except TypeError:
        return 0.0

    return min(angle_diff(upper_bound, angle), angle_diff(lower_bound, angle))


def select_realistic_bond_constraint(listener: Any, atom1: dict, atom2: dict,
                                     alt_atom_id1: Optional[str], alt_atom_id2: Optional[str], dst_func: dict,
                                     verbose: bool, log: IO) -> Tuple[dict, dict]:
    """ Return realistic bond constraint taking into account the current coordinates.
        @param listener: MR parser listener providing hasCoord, cR and the atom_site item names
    """

    if not listener.hasCoord:
        return atom1, atom2

    try:

        cR = listener.cR
        authAsymId, authSeqId, authAtomId = listener.authAsymId, listener.authSeqId, listener.authAtomId
        modelNumName = listener.modelNumName
        representativeModelId, representativeAltId = listener.representativeModelId, listener.representativeAltId

        _p1 =\
            cR.getDictListWithFilter('atom_site',
                                     CARTN_DATA_ITEMS,
                                     [{'name': authAsymId, 'type': 'str', 'value': atom1['chain_id']},
                                      {'name': authSeqId, 'type': 'int', 'value': atom1['seq_id']},
                                      {'name': authAtomId, 'type': 'str', 'value': atom1['atom_id']},
                                      {'name': modelNumName, 'type': 'int',
                                       'value': representativeModelId},
                                      {'name': 'label_alt_id', 'type': 'enum',
                                       'enum': (representativeAltId,)}
                                      ])

        if len(_p1) != 1:
            return atom1, atom2

        p1 = to_np_array(_p1[0])

        _p2 =\
            cR.getDictListWithFilter('atom_site',
                                     CARTN_DATA_ITEMS,
                                     [{'name': authAsymId, 'type': 'str', 'value': atom2['chain_id']},
                                      {'name': authSeqId, 'type': 'int', 'value': atom2['seq_id']},
                                      {'name': authAtomId, 'type': 'str', 'value': atom2['atom_id']},
                                      {'name': modelNumName, 'type': 'int',
                                       'value': representativeModelId},
                                      {'name': 'label_alt_id', 'type': 'enum',
                                       'enum': (representativeAltId,)}
                                      ])

        if len(_p2) != 1:
            return atom1, atom2

        p2 = to_np_array(_p2[0])

        d_org = distance(p1, p2)

        lower_limit = dst_func.get('lower_limit')
        if lower_limit is not None:
            lower_limit = float(lower_limit)
        upper_limit = dst_func.get('upper_limit')
        if upper_limit is not None:
            upper_limit = float(upper_limit)

        if alt_atom_id1 is not None:

            _p1 =\
                cR.getDictListWithFilter('atom_site',
                                         CARTN_DATA_ITEMS,
                                         [{'name': authAsymId, 'type': 'str', 'value': atom1['chain_id']},
                                          {'name': authSeqId, 'type': 'int', 'value': atom1['seq_id']},
                                          {'name': authAtomId, 'type': 'str', 'value': alt_atom_id1},
                                          {'name': modelNumName, 'type': 'int',
                                           'value': representativeModelId},
                                          {'name': 'label_alt_id', 'type': 'enum',
                                           'enum': (representativeAltId,)}
                                          ])

            if len(_p1) != 1:
                return atom1, atom2

            p1_alt = to_np_array(_p1[0])

            d_alt = distance(p1_alt, p2)

            if dist_error(lower_limit, upper_limit, d_org) > dist_error(lower_limit, upper_limit, d_alt):
                if 'auth_atom_id' not in atom1:
                    atom1['auth_atom_id'] = atom1['atom_id']
                atom1['atom_id'] = alt_atom_id1

        elif alt_atom_id2 is not None:

            _p2 =\
                cR.getDictListWithFilter('atom_site',
                                         CARTN_DATA_ITEMS,
                                         [{'name': authAsymId, 'type': 'str', 'value': atom2['chain_id']},
                                          {'name': authSeqId, 'type': 'int', 'value': atom2['seq_id']},
                                          {'name': authAtomId, 'type': 'str', 'value': alt_atom_id2},
                                          {'name': modelNumName, 'type': 'int',
                                           'value': representativeModelId},
                                          {'name': 'label_alt_id', 'type': 'enum',
                                           'enum': (representativeAltId,)}
                                          ])

            if len(_p2) != 1:
                return atom1, atom2

            p2_alt = to_np_array(_p2[0])

            d_alt = distance(p1, p2_alt)

            if dist_error(lower_limit, upper_limit, d_org) > dist_error(lower_limit, upper_limit, d_alt):
                if 'auth_atom_id' not in atom2:
                    atom2['auth_atom_id'] = atom2['atom_id']
                atom2['atom_id'] = alt_atom_id2

    except Exception as e:  # pylint: disable=broad-exception-caught
        if verbose:
            log.write(f"+{listener.__class_name__}.selectRealisticBondConstraint() ++ Error  - {str(e)}")

    return atom1, atom2


def select_realistic_chi2_angle_constraint(listener: Any, atom1: dict, atom2: dict, atom3: dict, atom4: dict,
                                           dst_func: dict, verbose: bool, log: IO,
                                           correct_circular_shift: bool = True) -> dict:
    """ Return realistic chi2 angle constraint taking into account the current coordinates.
        @param listener: MR parser listener providing hasCoord, cR and the atom_site item names
    """

    if not listener.hasCoord:
        return dst_func

    try:

        cR = listener.cR
        authAsymId, authSeqId, authAtomId = listener.authAsymId, listener.authSeqId, listener.authAtomId
        modelNumName = listener.modelNumName
        representativeModelId, representativeAltId = listener.representativeModelId, listener.representativeAltId

        _p1 =\
            cR.getDictListWithFilter('atom_site',
                                     CARTN_DATA_ITEMS,
                                     [{'name': authAsymId, 'type': 'str', 'value': atom1['chain_id']},
                                      {'name': authSeqId, 'type': 'int', 'value': atom1['seq_id']},
                                      {'name': authAtomId, 'type': 'str', 'value': atom1['atom_id']},
                                      {'name': modelNumName, 'type': 'int',
                                       'value': representativeModelId},
                                      {'name': 'label_alt_id', 'type': 'enum',
                                       'enum': (representativeAltId,)}
                                      ])

        if len(_p1) != 1:
            return dst_func

        p1 = to_np_array(_p1[0])

        _p2 =\
            cR.getDictListWithFilter('atom_site',
                                     CARTN_DATA_ITEMS,
                                     [{'name': authAsymId, 'type': 'str', 'value': atom2['chain_id']},
                                      {'name': authSeqId, 'type': 'int', 'value': atom2['seq_id']},
                                      {'name': authAtomId, 'type': 'str', 'value': atom2['atom_id']},
                                      {'name': modelNumName, 'type': 'int',
                                       'value': representativeModelId},
                                      {'name': 'label_alt_id', 'type': 'enum',
                                       'enum': (representativeAltId,)}
                                      ])

        if len(_p2) != 1:
            return dst_func

        p2 = to_np_array(_p2[0])

        _p3 =\
            cR.getDictListWithFilter('atom_site',
                                     CARTN_DATA_ITEMS,
                                     [{'name': authAsymId, 'type': 'str', 'value': atom3['chain_id']},
                                      {'name': authSeqId, 'type': 'int', 'value': atom3['seq_id']},
                                      {'name': authAtomId, 'type': 'str', 'value': atom3['atom_id']},
                                      {'name': modelNumName, 'type': 'int',
                                       'value': representativeModelId},
                                      {'name': 'label_alt_id', 'type': 'enum',
                                       'enum': (representativeAltId,)}
                                      ])

        if len(_p3) != 1:
            return dst_func

        p3 = to_np_array(_p3[0])

        _p4 =\
            cR.getDictListWithFilter('atom_site',
                                     CARTN_DATA_ITEMS,
                                     [{'name': authAsymId, 'type': 'str', 'value': atom4['chain_id']},
                                      {'name': authSeqId, 'type': 'int', 'value': atom4['seq_id']},
                                      {'name': authAtomId, 'type': 'str', 'value': 'CD1'},
                                      {'name': modelNumName, 'type': 'int',
                                       'value': representativeModelId},
                                      {'name': 'label_alt_id', 'type': 'enum',
                                       'enum': (representativeAltId,)}
                                      ])

        if len(_p4) != 1:
            return dst_func

        p4 = to_np_array(_p4[0])

        chi2 = dihedral_angle(p1, p2, p3, p4)

        _p4 =\
            cR.getDictListWithFilter('atom_site',
                                     CARTN_DATA_ITEMS,
                                     [{'name': authAsymId, 'type': 'str', 'value': atom4['chain_id']},
                                      {'name': authSeqId, 'type': 'int', 'value': atom4['seq_id']},
                                      {'name': authAtomId, 'type': 'str', 'value': 'CD2'},
                                      {'name': modelNumName, 'type': 'int',
                                       'value': representativeModelId},
                                      {'name': 'label_alt_id', 'type': 'enum',
                                       'enum': (representativeAltId,)}
                                      ])

        if len(_p4) != 1:
            return dst_func

        alt_p4 = to_np_array(_p4[0])

        alt_chi2 = dihedral_angle(p1, p2, p3, alt_p4)

        target_value = dst_func.get('target_value')
        if target_value is not None:
            target_value = float(target_value)
        target_value_uncertainty = dst_func.get('target_value_uncertainty')
        if target_value_uncertainty is not None:
            target_value_uncertainty = float(target_value_uncertainty)

        lower_limit = dst_func.get('lower_limit')
        if lower_limit is not None:
            lower_limit = float(lower_limit)
        upper_limit = dst_func.get('upper_limit')
        if upper_limit is not None:
            upper_limit = float(upper_limit)

        lower_linear_limit = dst_func.get('lower_linear_limit')
        if lower_linear_limit is not None:
            lower_linear_limit = float(lower_linear_limit)
        upper_linear_limit = dst_func.get('upper_linear_limit')
        if upper_linear_limit is not None:
            upper_linear_limit = float(upper_linear_limit)

        target_value, lower_bound, upper_bound =\
            angle_target_values(target_value, target_value_uncertainty,
                                lower_limit, upper_limit,
                                lower_linear_limit, upper_linear_limit)

        if target_value is None:
            return dst_func

        if angle_error(lower_bound, upper_bound, target_value, chi2) >\
           angle_error(lower_bound, upper_bound, target_value, alt_chi2):
            target_value = dst_func.get('target_value')
            if target_value is not None:
                target_value = float(target_value) + 180.0
            lower_limit = dst_func.get('lower_limit')
            if lower_limit is not None:
                lower_limit = float(lower_limit) + 180.0
            upper_limit = dst_func.get('upper_limit')
            if upper_limit is not None:
                upper_limit = float(upper_limit) + 180.0

            if lower_linear_limit is not None:
                lower_linear_limit += 180.0
            if upper_linear_limit is not None:
                upper_linear_limit += 180.0

            _array = numpy.array([target_value, lower_limit, upper_limit, lower_linear_limit, upper_linear_limit],
                                 dtype=float)

            shift = 0.0
            if correct_circular_shift:
                if numpy.nanmin(_array) >= THRESHOLD_FOR_CIRCULAR_SHIFT:
                    shift = -(numpy.nanmax(_array) // 360) * 360
                elif numpy.nanmax(_array) <= -THRESHOLD_FOR_CIRCULAR_SHIFT:
                    shift = -(numpy.nanmin(_array) // 360) * 360
            if target_value is not None:
                dst_func['target_value'] = str(target_value + shift)
            if lower_limit is not None:
                dst_func['lower_limit'] = str(lower_limit + shift)
            if upper_limit is not None:
                dst_func['upper_limit'] = str(upper_limit + shift)
            if lower_linear_limit is not None:
                dst_func['lower_linear_limit'] = str(lower_linear_limit + shift)
            if upper_linear_limit is not None:
                dst_func['upper_linear_limit'] = str(upper_linear_limit + shift)

    except Exception as e:  # pylint: disable=broad-exception-caught
        if verbose:
            log.write(f"+{listener.__class_name__}.selectRealisticChi2AngleConstraint() ++ Error  - {str(e)}")

    return dst_func


def rdc_dmax(atom_type_1: str, atom_type_2: str, bond_distance: float, hz_unit: bool = True):
    """ Return scale factor of residual coupling constant for a given nuclei vector.
    """

    return - PERMEABILITY_0 * (PLANCK_CONSTANT if hz_unit else REDUCED_PLANCK_CONSTANT)\
        * GYROMAGNETIC_RATIOS[atom_type_1] * GYROMAGNETIC_RATIOS[atom_type_2]\
        / math.pow(2.0e-10 * numpy.pi * bond_distance, 3.0)


def rdc_target_values(target_value: Optional[float], target_value_uncertainty: Optional[float],
                      value: Optional[float], value_uncertainty: Optional[float],
                      lower_limit: Optional[float], upper_limit: Optional[float],
                      lower_linear_limit: Optional[float], upper_linear_limit: Optional[float]
                      ) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    """ Return estimated RDC target value, lower_bound, upper_bound.
        @author: Masashi Yokochi
    """

    lower_bound, upper_bound = lower_limit, upper_limit

    if lower_bound is None and upper_bound is None and lower_linear_limit is None and upper_linear_limit is None:
        if target_value is None and value is None:
            return None, lower_bound, upper_bound
        if target_value is not None:
            if target_value_uncertainty is not None:
                lower_bound = target_value - target_value_uncertainty
                upper_bound = target_value + target_value_uncertainty
            else:
                lower_bound = upper_bound = target_value
        else:
            target_value = value
            target_value_uncertainty = value_uncertainty
            if target_value_uncertainty is not None:
                lower_bound = target_value - target_value_uncertainty
                upper_bound = target_value + target_value_uncertainty
            else:
                lower_bound = upper_bound = target_value

    if lower_bound is None and lower_linear_limit is not None:
        lower_bound = lower_linear_limit
    if upper_bound is None and upper_linear_limit is not None:
        upper_bound = upper_linear_limit

    if target_value is None:
        if None not in (lower_bound, upper_bound):
            target_value = (lower_bound + upper_bound) / 2.0
        elif lower_bound is not None:
            target_value = lower_bound
        elif upper_bound is not None:
            target_value = upper_bound

    return target_value, lower_bound, upper_bound


def rdc_error(lower_bound: Optional[float], upper_bound: Optional[float], rdc: float) -> float:
    """ Return RDC outlier for given lower_bound and upper_bound.
        @author: Masashi Yokochi
    """

    error = 0.0

    try:

        if None not in (lower_bound, upper_bound):
            if lower_bound <= rdc <= upper_bound:
                pass
            elif rdc > upper_bound:
                error = abs(rdc - upper_bound)
            else:
                error = abs(rdc - lower_bound)

        elif upper_bound is not None:
            if rdc <= upper_bound:
                pass
            elif rdc > upper_bound:
                error = abs(rdc - upper_bound)

        elif lower_bound is not None:
            if lower_bound <= rdc:
                pass
            else:
                error = abs(rdc - lower_bound)

    except TypeError:
        pass

    return error


def get_violated_model_ids(viol_per_model: List[dict]
                           ) -> List[int]:
    """ Return the list of violated model IDs.
    """

    return [m for m, err in viol_per_model.items() if err is not None and err > 0.0]


def get_violation_statistics_for_each_bin(beg_err_bin: Optional[float], end_err_bin: Optional[float],
                                          total_models: int, eff_model_ids: List[int], viol_dict: dict
                                          ) -> Tuple[Optional[float], Optional[float], Optional[int],
                                                     Optional[float], Optional[float], Optional[int], Optional[float]]:
    """ Calculate violation statistics for each bin.
        The results were summarized against all models and per model, respectively.
    """

    viol_stat_per_model, all_err_list = [], []

    for m in eff_model_ids:
        err_list = []

        for viol_per_model in viol_dict.values():
            err = viol_per_model[m]

            if err is None or err == 0.0:
                continue

            if (end_err_bin is None and beg_err_bin < err)\
               or beg_err_bin < err <= end_err_bin\
               or (beg_err_bin is None and err < end_err_bin):
                err_list.append(err)
                all_err_list.append(err)

        if len(err_list) == 0:
            viol_stat = [None, None, None]
        else:
            viol_stat = [round(min(err_list), 2),
                         round(max(err_list), 2),
                         len(err_list)]

        viol_stat_per_model.append(viol_stat)

    if len(all_err_list) == 0:
        viol_stat_all_model = [None, None, None, None]
    else:
        viol_stat_all_model = [round(min(all_err_list), 2),
                               round(max(all_err_list), 2),
                               len(all_err_list),
                               round(float(len(all_err_list)) / float(total_models), 1)]

    return viol_stat_all_model, viol_stat_per_model


def probability_density(value: float, mean: float, stddev: float) -> float:
    """ Return probability density.
    """

    stddev2 = stddev ** 2.0

    return math.exp(-((value - mean) ** 2.0) / (2.0 * stddev2)) / math.sqrt(2.0 * math.pi * stddev2)


def predict_redox_state_of_cysteine(ca_chem_shift: Optional[float], cb_chem_shift: Optional[float]
                                    ) -> Tuple[float, float]:
    """ Return prediction of redox state of Cysteine using assigned CA, CB chemical shifts.
        @return: probability of oxidized state, probability of reduced state
        Reference:
          13C NMR chemical shifts can predict disulfide bond formation.
          Sharma, D., Rajarathnam, K.
          J Biomol NMR 18, 165–171 (2000).
          DOI: 10.1023/A:1008398416292
    """

    oxi_ca = {'avr': 55.5, 'std': 2.5}
    oxi_cb = {'avr': 40.7, 'std': 3.8}

    red_ca = {'avr': 59.3, 'std': 3.2}
    red_cb = {'avr': 28.3, 'std': 2.2}

    oxi = 1.0
    red = 1.0

    if ca_chem_shift is not None:
        oxi *= probability_density(ca_chem_shift, oxi_ca['avr'], oxi_ca['std'])
        red *= probability_density(ca_chem_shift, red_ca['avr'], red_ca['std'])

    if cb_chem_shift is not None:
        if cb_chem_shift < 32.0:
            oxi = 0.0
        else:
            oxi *= probability_density(cb_chem_shift, oxi_cb['avr'], oxi_cb['std'])
        if cb_chem_shift > 35.0:
            red = 0.0
        else:
            red *= probability_density(cb_chem_shift, red_cb['avr'], red_cb['std'])

    total = oxi + red

    if total in (0.0, 2.0):
        return 0.0, 0.0

    return oxi / total, red / total


def predict_cis_trans_peptide_of_proline(cb_chem_shift, cg_chem_shift
                                         ) -> Tuple[float, float]:
    """ Return prediction of cis-trans peptide bond of Proline using assigned CB, CG chemical shifts.
        @return: probability of cis-peptide bond, probability of trans-peptide bond
        Reference:
          A software tool for the prediction of Xaa-Pro peptide bond conformations in proteins
          based on 13C chemical shift statistics.
          Schubert, M., Labudde, D., Oschkinat, H. et al.
          J Biomol NMR 24, 149–154 (2002)
          DOI: 10.1023/A:1020997118364
    """

    cis_cb = {'avr': 34.16, 'std': 1.15, 'max': 36.23, 'min': 30.74}
    cis_cg = {'avr': 24.52, 'std': 1.09, 'max': 27.01, 'min': 22.10}
    cis_dl = {'avr': 9.64, 'std': 1.27}

    trs_cb = {'avr': 31.75, 'std': 0.98, 'max': 35.83, 'min': 26.30}
    trs_cg = {'avr': 27.26, 'std': 1.05, 'max': 33.39, 'min': 19.31}
    trs_dl = {'avr': 4.51, 'std': 1.17}

    cis = 1.0
    trs = 1.0

    if cb_chem_shift is not None:
        if cb_chem_shift < cis_cb['min'] - cis_cb['std'] or cb_chem_shift > cis_cb['max'] + cis_cb['std']:
            cis = 0.0
        else:
            cis *= probability_density(cb_chem_shift, cis_cb['avr'], cis_cb['std'])
        if cb_chem_shift < trs_cb['min'] - trs_cb['std'] or cb_chem_shift > trs_cb['max'] + trs_cb['std']:
            trs = 0.0
        else:
            trs *= probability_density(cb_chem_shift, trs_cb['avr'], trs_cb['std'])

    if cg_chem_shift is not None:
        if cg_chem_shift < cis_cg['min'] - cis_cg['std'] or cg_chem_shift > cis_cg['max'] + cis_cg['std']:
            cis = 0.0
        else:
            cis *= probability_density(cg_chem_shift, cis_cg['avr'], cis_cg['std'])
        if cg_chem_shift < trs_cg['min'] - trs_cg['std'] or cg_chem_shift > trs_cg['max'] + trs_cg['std']:
            trs = 0.0
        else:
            trs *= probability_density(cg_chem_shift, trs_cg['avr'], trs_cg['std'])

    if (cb_chem_shift is not None) and (cg_chem_shift is not None):
        delta_shift = cb_chem_shift - cg_chem_shift

        cis *= probability_density(delta_shift, cis_dl['avr'], cis_dl['std'])
        trs *= probability_density(delta_shift, trs_dl['avr'], trs_dl['std'])

    total = cis + trs

    if total in (0.0, 2.0):
        return 0.0, 0.0

    return cis / total, trs / total


def predict_tautomer_state_of_histidine(cg_chem_shift: Optional[float], cd2_chem_shift: Optional[float],
                                        nd1_chem_shift: Optional[float], ne2_chem_shift: Optional[float]
                                        ) -> Tuple[float, float]:
    """ Return prediction of tautomeric state of Histidine using assigned CG, CD2, ND1, and NE2 chemical shifts.
        @return: probability of biprotonated, probability of tau tautomer, probability of pi tautomer
        Reference:
          Protonation, Tautomerization, and Rotameric Structure of Histidine:
          A Comprehensive Study by Magic-Angle-Spinning Solid-State NMR.
          Shenhui Li and Mei Hong.
          Journal of the American Chemical Society 2011 133 (5), 1534-1544
          DOI: 10.1021/ja108943n
    """

    bip_cg = {'avr': 131.2, 'std': 0.7}
    bip_cd2 = {'avr': 120.6, 'std': 1.3}
    bip_nd1 = {'avr': 190.0, 'std': 1.9}
    bip_ne2 = {'avr': 176.3, 'std': 1.9}

    tau_cg = {'avr': 135.7, 'std': 2.2}
    tau_cd2 = {'avr': 116.9, 'std': 2.1}
    tau_nd1 = {'avr': 249.4, 'std': 1.9}
    tau_ne2 = {'avr': 171.1, 'std': 1.9}

    pi_cg = {'avr': 125.7, 'std': 2.2}
    pi_cd2 = {'avr': 125.6, 'std': 2.1}
    pi_nd1 = {'avr': 171.8, 'std': 1.9}
    pi_ne2 = {'avr': 248.2, 'std': 1.9}

    bip = 1.0
    tau = 1.0
    pi = 1.0

    if cg_chem_shift is not None:
        bip *= probability_density(cg_chem_shift, bip_cg['avr'], bip_cg['std'])
        tau *= probability_density(cg_chem_shift, tau_cg['avr'], tau_cg['std'])
        pi *= probability_density(cg_chem_shift, pi_cg['avr'], pi_cg['std'])

    if cd2_chem_shift is not None:
        bip *= probability_density(cd2_chem_shift, bip_cd2['avr'], bip_cd2['std'])
        tau *= probability_density(cd2_chem_shift, tau_cd2['avr'], tau_cd2['std'])
        pi *= probability_density(cd2_chem_shift, pi_cd2['avr'], pi_cd2['std'])

    if nd1_chem_shift is not None:
        bip *= probability_density(nd1_chem_shift, bip_nd1['avr'], bip_nd1['std'])
        tau *= probability_density(nd1_chem_shift, tau_nd1['avr'], tau_nd1['std'])
        pi *= probability_density(nd1_chem_shift, pi_nd1['avr'], pi_nd1['std'])

    if ne2_chem_shift is not None:
        bip *= probability_density(ne2_chem_shift, bip_ne2['avr'], bip_ne2['std'])
        tau *= probability_density(ne2_chem_shift, tau_ne2['avr'], tau_ne2['std'])
        pi *= probability_density(ne2_chem_shift, pi_ne2['avr'], pi_ne2['std'])

    total = bip + tau + pi

    if total in (0.0, 3.0):
        return 0.0, 0.0, 0.0

    return bip / total, tau / total, pi / total


def predict_rotamer_state_of_leucine(cd1_chem_shift: Optional[float], cd2_chem_shift: Optional[float]
                                     ) -> Tuple[float, float]:
    """ Return prediction of rotameric state of Leucine using assigned CD1 and CD2 chemical shifts.
        @return: probability of gauche+, trans, gauche-
        Reference:
          Dependence of Amino Acid Side Chain 13C Shifts on Dihedral Angle: Application to Conformational Analysis.
          Robert E. London, Brett D. Wingad, and Geoffrey A. Mueller.
          Journal of the American Chemical Society 2008 130 (33), 11097-11105
          DOI: 10.1021/ja802729t
    """

    if None not in (cd1_chem_shift, cd2_chem_shift):

        delta = cd1_chem_shift - cd2_chem_shift

        pt = (delta + 5.0) / 10.0

        if 0.0 <= pt <= 1.0:
            return 1.0 - pt, pt, 0.0

    gp_cd1 = {'avr': 24.45, 'std': 1.58}
    gp_cd2 = {'avr': 25.79, 'std': 1.68}

    t_cd1 = {'avr': 25.17, 'std': 1.58}
    t_cd2 = {'avr': 23.84, 'std': 1.68}

    gp = 1.0
    t = 1.0

    if cd1_chem_shift is not None:
        gp *= probability_density(cd1_chem_shift, gp_cd1['avr'], gp_cd1['std'])
        t *= probability_density(cd1_chem_shift, t_cd1['avr'], t_cd1['std'])

    if cd2_chem_shift is not None:
        gp *= probability_density(cd2_chem_shift, gp_cd2['avr'], gp_cd2['std'])
        t *= probability_density(cd2_chem_shift, t_cd2['avr'], t_cd2['std'])

    total = gp + t

    if total in (0.0, 2.0):
        return 0.0, 0.0, 0.0

    return gp / total, t / total, 0.0


def predict_rotamer_state_of_valine(cg1_chem_shift: Optional[float], cg2_chem_shift: Optional[float]
                                    ) -> Tuple[float, float]:
    """ Return prediction of rotameric state of Valine using assigned CG1 and CG2 chemical shifts.
        @return: probability of gauche+, trans, gauche-
        Reference:
          Dependence of Amino Acid Side Chain 13C Shifts on Dihedral Angle: Application to Conformational Analysis.
          Robert E. London, Brett D. Wingad, and Geoffrey A. Mueller.
          Journal of the American Chemical Society 2008 130 (33), 11097-11105
          DOI: 10.1021/ja802729t
    """

    gm_cg1 = {'avr': 22.05, 'std': 1.36}
    gm_cg2 = {'avr': 20.1, 'std': 1.55}

    gp_cg1 = {'avr': 20.87, 'std': 1.36}
    gp_cg2 = {'avr': 21.23, 'std': 1.55}

    t_cg1 = {'avr': 21.74, 'std': 1.36}
    t_cg2 = {'avr': 21.97, 'std': 1.55}

    gm = 1.0
    gp = 1.0
    t = 1.0

    if cg1_chem_shift is not None:
        gm *= probability_density(cg1_chem_shift, gm_cg1['avr'], gm_cg1['std'])
        gp *= probability_density(cg1_chem_shift, gp_cg1['avr'], gp_cg1['std'])
        t *= probability_density(cg1_chem_shift, t_cg1['avr'], t_cg1['std'])

    if cg2_chem_shift is not None:
        gm *= probability_density(cg2_chem_shift, gm_cg2['avr'], gm_cg2['std'])
        gp *= probability_density(cg2_chem_shift, gp_cg2['avr'], gp_cg2['std'])
        t *= probability_density(cg2_chem_shift, t_cg2['avr'], t_cg2['std'])

    total = gm + gp + t

    if total in (0.0, 3.0):
        return 0.0, 0.0, 0.0

    return gp / total, t / total, gm / total


def predict_rotamer_state_of_isoleucine(cd1_chem_shift: Optional[float]
                                        ) -> Tuple[float, float, float]:
    """ Return prediction of rotameric state of Isoleucine using assigned CD1 chemical shift.
        @return: probability of gauche+, trans, gauche-
        Reference:
          Determination of Isoleucine Side-Chain Conformations in Ground and Excited States of Proteins from Chemical Shifts.
          D. Flemming Hansen, Philipp Neudecker, and Lewis E. Kay.
          Journal of the American Chemical Society 2010 132 (22), 7589-7591
          DOI: 10.1021/ja102090z
    """

    if cd1_chem_shift is None:
        return 0.0, 0.0, 0.0

    if cd1_chem_shift < 9.3:
        return 0.0, 0.0, 1.0

    if cd1_chem_shift > 14.8:
        return 1.0 * (4.0 / 85.0), 1.0 * (81.0 / 85.0), 0.0

    pgm = (14.8 - cd1_chem_shift) / 5.5

    return (1.0 - pgm) * (4.0 / 85.0), (1.0 - pgm) * (81.0 / 85.0), pgm


class NmrVrptUtility:
    """ Wrapper class for NMR restraint analysis.
    """
    __slots__ = ('__class_name__',
                 '__version__',
                 '__verbose',
                 '__log',
                 '__debug',
                 '__trust_pdbx_auth_atom_name',
                 '__permit_missing_chem_shift',
                 '__use_cache',
                 '__inputParamDict',
                 '__outputParamDict',
                 '__cifPath',
                 '__nmrDataPath',
                 '__dirPath',
                 '__cacheDirPath',
                 '__cifHashCode',
                 '__nmrDataHashCode',
                 '__resultsCacheName',
                 '__inputReport',
                 '__cR',
                 '__rR',
                 '__ccU',
                 '__csStat',
                 '__caC',
                 '__representative_model_id',
                 '__representative_alt_id',
                 '__total_models',
                 '__eff_model_ids',
                 '__entityInstance',
                 '__entityUninstance',
                 '__atomIdList',
                 '__coordinates',
                 '__chemShiftMeta',
                 '__chemShiftTotal',
                 '__chemShiftDict',
                 '__chemShiftUniqDict',
                 '__chemShiftPrimary',
                 '__chemShiftOutlier',
                 '__chemShiftDuplicated',
                 '__chemShiftUnparsed',
                 '__chemShiftUnmapped',
                 '__chemShiftUnmodeled',
                 '__distRestDict',
                 '__distRestDictWithCombKey',
                 '__distRestSeqDict',
                 '__dihedRestDict',
                 '__dihedRestDictWithCombKey',
                 '__dihedRestSeqDict',
                 '__rdcRestDict',
                 '__rdcRestDictWithCombKey',
                 '__rdcRestSeqDict',
                 '__rdcSaupeOrderMatrix',
                 '__rdcCalcDict',
                 '__rdcSyntCalcDict',
                 '__rdcCorrPlotDict',
                 '__distRestViolDict',
                 '__distRestUnmapped',
                 '__distRestViolCombKeyDict',
                 '__dihedRestViolDict',
                 '__dihedRestUnmapped',
                 '__dihedRestViolCombKeyDict',
                 '__rdcRestViolDict',
                 '__rdcRestUnmapped',
                 '__rdcRestViolCombKeyDict',
                 '__results',
                 '__has_prev_results',
                 '__is_diamagnetic',
                 '__procTasksDict')

    def __init__(self, verbose: bool = False, log: IO = sys.stderr,
                 cR: Optional[CifReader] = None, caC: Optional[dict] = None, ccU: Optional[ChemCompUtil] = None,
                 csStat: Optional[BmrbChemShiftStat] = None) -> None:
        self.__class_name__ = self.__class__.__name__
        self.__version__ = __version__

        self.__verbose = verbose
        self.__log = log

        self.__debug = False
        self.__trust_pdbx_auth_atom_name = False
        self.__permit_missing_chem_shift = False
        self.__use_cache = cR is not None and caC is not None

        # auxiliary input resource
        self.__inputParamDict = {}

        # auxiliary output resource
        self.__outputParamDict = {}

        # CIF file path
        self.__cifPath = None

        # NMR data file path
        self.__nmrDataPath = None

        # current working directory
        self.__dirPath = None

        # directory for cache files
        self.__cacheDirPath = None

        # hash code of the coordinate file
        self.__cifHashCode = None

        # hash code of the NMR data file
        self.__nmrDataHashCode = None

        # cache file name for results
        self.__resultsCacheName = None

        # input NMR data processing report
        self.__inputReport = None

        # CIF reader
        self.__cR = cR

        # NMR data reader
        self.__rR = None

        # CCD accessing utility
        self.__ccU = ChemCompUtil(verbose, log) if ccU is None else ccU

        # BMRB chemical shift statistics
        self.__csStat = BmrbChemShiftStat(verbose, log, self.__ccU) if csStat is None else csStat

        # ParserListerUtil.coordAssemblyChecker()
        self.__caC = caC

        # representative model id
        self.__representative_model_id = REPRESENTATIVE_MODEL_ID
        # representative_alt_id
        self.__representative_alt_id = REPRESENTATIVE_ALT_ID
        # total number of models
        self.__total_models = 0
        # list of effective model_id
        self.__eff_model_ids = None

        # entity instances
        self.__entityInstance = None
        # entity unmodeled (DAOTHER-10987)
        self.__entityUninstance = None

        # atom id list for each model_id and atom key (auth_asym_id, auth_seq_id, auth_comp_id, auth_atom_id, PDB_ins_code)
        self.__atomIdList = None
        # coordinates for each model_id and atom key
        self.__coordinates = None

        # metadata of assigned chemical shift lists
        self.__chemShiftMeta = None
        # total cs row size of assigned chemical lists
        self.__chemShiftTotal = None
        # chemical shifts for each list_id
        self.__chemShiftDict = None

        # distance restraints for each restraint key (list_id, restraint_id)
        self.__distRestDict = None
        # distance restraints for each restraint key (list_id, restraint_id, tuple(combination_id, member_id))
        self.__distRestDictWithCombKey = None
        # distance restraint keys for each sequence key (auth_asym_id, auth_seq_id, auth_comp_id)
        self.__distRestSeqDict = None

        # dihedral angle restraints for each restraint key (list_id, restraint_id)
        self.__dihedRestDict = None
        # dihedral angle restraints for each restraint key (list_id, restraint_id, tuple(combination_id))
        self.__dihedRestDictWithCombKey = None
        # dihedral angle restraint keys for each sequence key (auth_asym_id, auth_seq_id, auth_comp_id)
        self.__dihedRestSeqDict = None

        # RDC restraints for each restraint key (list_id, restraint_id)
        self.__rdcRestDict = None
        # RDC restraints for each restraint key (list_id, restraint_id, tuple(combination_id))
        self.__rdcRestDictWithCombKey = None
        # RDC restraint keys for each sequence key (auth_asym_id, auth_seq_id, auth_comp_id)
        self.__rdcRestSeqDict = None

        # Saupe order matrix (Molecular alignment tensor)
        self.__rdcSaupeOrderMatrix = None
        # calculated RDC values based on molecular alignment tensor
        self.__rdcCalcDict = None
        # calculated RDC values for synthetic RDC values derived from Monte Carlo simulations
        self.__rdcSyntCalcDict = None
        # RDC correlation plot of observed and calculated RDCs for each list
        self.__rdcCorrPlotDict = None

        # unique chemical shifts for each list_id and atom_key
        self.__chemShiftUniqDict = None
        # list of chemical shift outliers for each list_id
        self.__chemShiftOutlier = None
        # list of duplicated chemical shifts for each list_id
        self.__chemShiftDuplicated = None
        # list of unparsed chemical shifts for each list_id
        self.__chemShiftUnparsed = None
        # list of unmapped chemical shifts for each list_id
        self.__chemShiftUnmapped = None
        # list of mapped chemical shifts to unmodeled residues for each list_id
        self.__chemShiftUnmodeled = None

        # distance restraint violations for each restraint key
        self.__distRestViolDict = None
        # list of restraint key of unmapped distance restraints
        self.__distRestUnmapped = None
        # combination keys (Combination_ID/Member_ID) of violated distance restraints for each restraint key
        self.__distRestViolCombKeyDict = None

        # dihedral angle restraint violations for each restraint key
        self.__dihedRestViolDict = None
        # list of restraint key of unmapped dihedral angle restraints
        self.__dihedRestUnmapped = None
        # combination key (Combination_ID) of violated dihedral angle restraints for each restraint key
        self.__dihedRestViolCombKeyDict = None

        # RDC restraint violations for each restraint key
        self.__rdcRestViolDict = None
        # list of restraint key of unmapped RDC restraints
        self.__rdcRestUnmapped = None
        # combination key (Combination_ID) of violated RDC restraints for each restraint key
        self.__rdcRestViolCombKeyDict = None

        # summarized restraint analysis results
        self.__results = None

        # whether the previous results have been retrieved
        self.__has_prev_results = False

        # whether molecular assembly is diamagnetic.
        self.__is_diamagnetic = True

        # the tasks are functions, called with the instance, rather than bound methods, which would make each instance a
        # reference cycle (instance -> task list -> bound method -> instance) that holds its coordinates, NMR data and
        # report until a full garbage collection, e.g. the previous NMR unified data in 'nmr-str2str-deposit'
        __csValidTasks = [NmrVrptUtility.__parseCoordinate,
                          NmrVrptUtility.__parseNmrData,
                          NmrVrptUtility.__parseNmrDpReport,
                          NmrVrptUtility.__checkPreviousCsAnalysis,
                          NmrVrptUtility.__retrieveCoordAssemblyChecker,
                          NmrVrptUtility.__extractCoordAtomSites,
                          NmrVrptUtility.__extractEntityInstances,
                          NmrVrptUtility.__extractChemicalShifts,
                          NmrVrptUtility.__validateChemicalShifts,
                          NmrVrptUtility.__summarizeCommonCsAnalysis,
                          NmrVrptUtility.__outputResultsAsPickleFile]

        __mrValidTasks = [NmrVrptUtility.__parseCoordinate,
                          NmrVrptUtility.__parseNmrData,
                          NmrVrptUtility.__checkPreviousMrAnalysis,
                          NmrVrptUtility.__retrieveCoordAssemblyChecker,
                          NmrVrptUtility.__extractCoordAtomSites,
                          NmrVrptUtility.__extractGenDistConstraints,
                          NmrVrptUtility.__extractTorsionAngleConstraints,
                          NmrVrptUtility.__extractRdcConstraints,
                          NmrVrptUtility.__validateDistanceRestraints,
                          NmrVrptUtility.__validateDihedralAngleRestraints,
                          NmrVrptUtility.__validateRdcRestraints,
                          NmrVrptUtility.__summarizeCommonMrAnalysis,
                          NmrVrptUtility.__summarizeDistanceRestraintAnalysis,
                          NmrVrptUtility.__summarizeDihedralAngleRestraintAnalysis,
                          NmrVrptUtility.__summarizeRdcRestraintAnalysis,
                          NmrVrptUtility.__outputResultsAsPickleFile]

        # dictionary of processing tasks of each workflow operation
        self.__procTasksDict = {'nmr-cs-validation': __csValidTasks,
                                'nmr-mr-validation': __mrValidTasks,
                                'nmr-chemical-shift-validation': __csValidTasks,
                                'nmr-restraint-validation': __mrValidTasks}  # for backward compatibility

        # Pre-run GC to avoid OOM because NmrVrptUtility requires a mount of memory
        gc.collect(2)

    @property
    def version(self) -> str:
        """ Retrieve software version.
        """

        return __version__

    @property
    def dirPath(self) -> Optional[str]:
        """ Retrieve the current working directory.
        """

        return self.__dirPath

    @dirPath.setter
    def dirPath(self, dirPath: Optional[str]) -> None:
        """ Set cache directory path.
        """

        self.__dirPath = dirPath

        if dirPath is not None:
            if not os.path.isdir(dirPath):
                os.makedirs(dirPath)

    @property
    def cacheDirPath(self) -> Optional[str]:
        """ Retrieve cache directory path.
        """

        return self.__cacheDirPath

    @cacheDirPath.setter
    def cacheDirPath(self, cacheDirPath: Optional[str]) -> None:
        """ Set cache directory path.
        """

        self.__cacheDirPath = cacheDirPath

        if self.__use_cache and cacheDirPath is not None:
            if not os.path.isdir(cacheDirPath):
                os.makedirs(cacheDirPath)

    def getNextPath(self, src_path: str, suffix: str = '~') -> str:
        """ Return candidate next file path.
        """
        assert len(suffix) > 0

        src_path_next = src_path + suffix

        if self.__dirPath is not None:
            src_path_next = os.path.join(self.__dirPath, os.path.basename(src_path_next))

        return src_path_next

    def setVerbose(self, verbose: bool) -> None:
        """ Set verbose mode.
        """

        self.__verbose = verbose

    def setDebugMode(self, debug: bool) -> None:
        """ Set debug mode.
        """

        self.__debug = debug

    def trustPdbxAuthAtomName(self, trust_pdbx_auth_atom_name: bool) -> None:
        """ Whether to trust _atom_site.pdbx_auth_atom_name rather than _atom_site.auth_atom_id.
            @note: Set True for OneDep validation package.
        """

        self.__trust_pdbx_auth_atom_name = trust_pdbx_auth_atom_name

    def permitMissingChemShift(self, permit_missing_chem_shift: bool) -> None:
        """ Whether permit missing assigned chemical shift.
        """

        self.__permit_missing_chem_shift = permit_missing_chem_shift

    def useCache(self, use_cache: bool) -> None:
        """ Use cache file(s) of the previous run.
            Do not enable this for generation of wwPDB validation report because of no performance improvement.
        """

        self.__use_cache = use_cache

    def getResults(self) -> Optional[dict]:
        """ Return NMR restraint analysis result.
        """

        return self.__results

    def addInput(self, name: Optional[str] = None, value: Any = None, type: str = 'file') -> None:  # noqa: E501, pylint: disable=redefined-builtin,line-too-long
        """ Add a named input and value to the dictionary of input parameters.
        """

        try:

            if type == 'param':
                if name not in VRPT_INPUT_PARAM_KEYS:
                    raise KeyError(f"+{self.__class_name__}.addInput() ++ Error  - Unknown input param {name!r}.")
                self.__inputParamDict[name] = value
            elif type == 'file':
                if name not in VRPT_INPUT_FILE_KEYS:
                    raise KeyError(f"+{self.__class_name__}.addInput() ++ Error  - Unknown input file {name!r}.")
                self.__inputParamDict[name] = os.path.abspath(value)
            else:
                raise KeyError(f"+{self.__class_name__}.addInput() ++ Error  - Unknown input type {type!r}.")

        except Exception as e:  # pylint: disable=broad-exception-caught
            raise ValueError(f"+{self.__class_name__}.addInput() ++ Error  - " + str(e)) from e

    def addOutput(self, name: Optional[str] = None, value: Any = None, type: str = 'file') -> None:  # noqa: E501, pylint: disable=redefined-builtin,line-too-long
        """ Add a named input and value to the dictionary of output parameters.
        """

        try:

            if type == 'file':
                if name not in VRPT_OUTPUT_FILE_KEYS:
                    raise KeyError(f"+{self.__class_name__}.addOutput() ++ Error  - Unknown output file {name!r}.")
                self.__outputParamDict[name] = os.path.abspath(value)
                if name == RESULT_PKL_FILE_PATH_KEY:
                    self.__use_cache = True
            else:
                raise KeyError(f"+{self.__class_name__}.addOutput() ++ Error  - Unknown output type {type!r}.")

        except Exception as e:  # pylint: disable=broad-exception-caught
            raise ValueError(f"+{self.__class_name__}.addOutput() ++ Error  - " + str(e)) from e

    def op(self, op: str
           ) -> Optional[dict]:
        """ Perform a series of tasks for a given workflow operation.
        """

        if op not in VRPT_WORKFLOW_OPS:
            raise KeyError(f"+{self.__class_name__}.op() ++ Error  - Unknown workflow operation {op!r}.")

        if self.__verbose:
            self.__log.write(f"+{self.__class_name__}.op() starting op {op}, use_cache {self.__use_cache}\n")

        if op in self.__procTasksDict:

            for task in self.__procTasksDict[op]:

                if self.__verbose:
                    self.__log.write(f"+{self.__class_name__}.op() starting op {op} - task {task.__name__}\n")

                start_time = time.time()

                if not task(self):
                    break

                if self.__verbose:
                    end_time = time.time()
                    if end_time - start_time > 1.0:
                        self.__log.write(f"op: {op}, task: {task.__name__}, elapsed time: {end_time - start_time:.1f} sec\n")

        return self.__results

    def release(self) -> None:
        """ Release the memory held for the coordinates and NMR data.
        """

        # back to the initial state, which the next op expects, e.g. __extractEntityInstances() tests for None
        self.__coordinates = None
        self.__entityInstance = None

        self.__rR.release()
        self.__nmrDataPath = None

        # check self.__cR is not borrowed from NmrDpUtility before releasing memory for the cooridnates
        if not self.__use_cache:
            self.__cR.release()

        self.__cifPath = None

        gc.collect(2)

    def __parseCoordinate(self) -> bool:
        """ Parse coordinates.
        """

        def extract_ensemble():

            try:

                self.__total_models = 0
                self.__eff_model_ids = []

                ensemble = self.__cR.getDictList('pdbx_nmr_ensemble')

                if len(ensemble) > 0 and 'conformers_submitted_total_number' in ensemble[0]:

                    try:
                        self.__total_models = int(ensemble[0]['conformers_submitted_total_number'])
                    except ValueError:
                        pass

                if len(ensemble) == 0 or self.__total_models == 0:

                    ensemble = self.__cR.getDictList('rcsb_nmr_ensemble')

                    if len(ensemble) > 0 and 'conformers_submitted_total_number' in ensemble[0]:

                        try:
                            self.__total_models = int(ensemble[0]['conformers_submitted_total_number'])
                        except ValueError:
                            pass

                    else:

                        try:

                            model_ids = self.__cR.getDistinctValues('atom_site', 'pdbx_PDB_model_num', 'int')

                            if len(model_ids) > 0:
                                self.__representative_model_id = min(model_ids)
                                self.__total_models = len(model_ids)
                                self.__eff_model_ids = sorted(model_ids)

                        except Exception as e:  # pylint: disable=broad-exception-caught

                            if self.__verbose:
                                self.__log.write(f"+{self.__class_name__}.__parseCoordinate() ++ Error  - {str(e)}\n")

                if len(ensemble) > 0 and 'representative_conformer' in ensemble[0]:

                    try:

                        rep_model_id = int(ensemble[0]['representative_conformer'])

                        if 1 <= rep_model_id <= self.__total_models:
                            self.__representative_model_id = rep_model_id

                    except ValueError:
                        pass

                if len(self.__eff_model_ids) == 0:

                    try:

                        model_ids = self.__cR.getDistinctValues('atom_site', 'pdbx_PDB_model_num', 'int')

                        if len(model_ids) > 0:
                            self.__total_models = len(model_ids)
                            self.__eff_model_ids = sorted(model_ids)

                    except Exception as e:  # pylint: disable=broad-exception-caught

                        if self.__verbose:
                            self.__log.write(f"+{self.__class_name__}.__parseCoordinate() ++ Error  - {str(e)}\n")

                if self.__cR.hasItem('atom_site', 'label_alt_id'):
                    alt_id = self.__cR.getFirstValue('atom_site', 'label_alt_id')

                    if alt_id is not None:
                        self.__representative_alt_id = alt_id

                if self.__cR.hasCategory('chem_comp_atom'):

                    try:

                        type_symbols = self.__cR.getDictListWithFilter('chem_comp_atom',
                                                                       [{'name': 'type_symbol', 'type': 'str'}
                                                                        ])

                        if len(type_symbols) > 0:
                            for elem in type_symbols:
                                if elem in PARAMAGNETIC_ELEMENTS or elem in FERROMAGNETIC_ELEMENTS:
                                    self.__is_diamagnetic = False
                                    break

                    except Exception as e:  # pylint: disable=broad-exception-caught

                        if self.__verbose:
                            self.__log.write(f"+{self.__class_name__}.__parseCoordinate() ++ Error  - {str(e)}\n")

            except Exception:  # pylint: disable=broad-exception-caught
                pass

        if self.__cR is not None:

            self.__cifPath = self.__cR.getFilePath()

            if self.__cifPath is None:
                return True

            if self.__use_cache:

                if self.__dirPath is None:
                    self.__dirPath = os.path.dirname(self.__cifPath)

                if self.__cacheDirPath is None:
                    self.__cacheDirPath = os.path.join(self.__dirPath, SUB_DIR_NAME_FOR_CACHE)

                if not os.path.isdir(self.__cacheDirPath):
                    os.makedirs(self.__cacheDirPath)

                self.__cifHashCode = self.__cR.getHashCode()

            extract_ensemble()

            return True

        if not self.__checkCoordInputSource():

            if MODEL_FILE_PATH_KEY in self.__inputParamDict:

                err = f"No such {self.__inputParamDict[MODEL_FILE_PATH_KEY]!r} file."

                if self.__verbose:
                    self.__log.write(f"+{self.__class_name__}.__parseCoordinate() ++ Error  - {err}\n")

            return False

        file_name = os.path.basename(self.__cifPath)

        if self.__cifPath is None:

            err = f"{file_name!r} is invalid PDBx/mmCIF file."

            if self.__verbose:
                self.__log.write(f"+{self.__class_name__}.__parseCoordinate() ++ Error  - {err}\n")

            return False

        extract_ensemble()

        return True

    def __checkCoordInputSource(self) -> bool:
        """ Check input source of the coordinates.
        """

        if self.__cifPath is not None:
            return True

        if MODEL_FILE_PATH_KEY in self.__inputParamDict:

            fPath = self.__inputParamDict[MODEL_FILE_PATH_KEY]

            if fPath.endswith('.gz'):

                _fPath = os.path.splitext(fPath)[0]

                if not os.path.exists(_fPath):

                    try:

                        uncompress_gzip_file(fPath, _fPath)

                    except Exception as e:  # pylint: disable=broad-exception-caught

                        if self.__verbose:
                            self.__log.write(f"+{self.__class_name__}.__checkCoordInputSource() ++ Error  - {str(e)}\n")

                        return False

                fPath = _fPath

            try:

                if self.__use_cache:

                    if self.__dirPath is None:
                        self.__dirPath = os.path.dirname(fPath)

                    if self.__cacheDirPath is None:
                        self.__cacheDirPath = os.path.join(self.__dirPath, SUB_DIR_NAME_FOR_CACHE)

                    if not os.path.isdir(self.__cacheDirPath):
                        os.makedirs(self.__cacheDirPath)

                self.__cR = CifReader(self.__verbose, self.__log,
                                      use_cache=self.__use_cache)

                if self.__cR.parse(fPath):
                    self.__cifPath = fPath

                    if self.__use_cache:
                        self.__cifHashCode = self.__cR.getHashCode()

                    return True

            except Exception:  # pylint: disable=broad-exception-caught
                pass

        if CIF_READER_OBJ_KEY in self.__inputParamDict:

            self.__cR = self.__inputParamDict[CIF_READER_OBJ_KEY]

            self.__cifPath = self.__cR.getFilePath()

            if self.__use_cache:

                if self.__dirPath is None:
                    self.__dirPath = os.path.dirname(self.__cifPath)

                if self.__cacheDirPath is None:
                    self.__cacheDirPath = os.path.join(self.__dirPath, SUB_DIR_NAME_FOR_CACHE)

                if not os.path.isdir(self.__cacheDirPath):
                    os.makedirs(self.__cacheDirPath)

                self.__cifHashCode = self.__cR.getHashCode()

            return True

        return False

    def __parseNmrData(self) -> bool:
        """ Parse NMR data.
        """

        if not self.__checkNmrDataInputSource():

            if NMR_CIF_FILE_PATH_KEY in self.__inputParamDict:

                err = f"No such {self.__inputParamDict[NMR_CIF_FILE_PATH_KEY]!r} file."

                if self.__verbose:
                    self.__log.write(f"+{self.__class_name__}.__parseNmrData() ++ Error  - {err}\n")

            return False

        return True

    def __parseNmrStarAsCif(self, strPath: str, cifPath: Optional[str] = None) -> bool:
        """ Parse the saveframes of an NMR-STAR file into self.__rR as CIF data blocks, and write them to 'cifPath' too,
            unless it is None. getHashCode() of the reader then returns the digest of the NMR-STAR file.
            The data blocks used to be written to a CIF file in any case, only to be parsed back: that round trip cost
            about 100 s of CPU per call for the 112 MB NMR-STAR file of 6cc9, and gave the data blocks parsed here.
            @return: True for success or False otherwise
        """

        if strPath.endswith(('.gz', '.bz2', '.xz', '.zip')):
            containerList = IoAdapterPy(False, self.__log).readFile(strPath)

        else:
            containerList = []

            # as IoAdapterPy.readFile() did, keep the data blocks parsed before an error, if any
            try:
                with open(strPath, 'r', encoding='utf-8', errors='ignore') as ifh:
                    SharingPdbxReader(ifh).read(containerList)
            except Exception:  # pylint: disable=broad-exception-caught
                pass

        if containerList is None or len(containerList) < 2:
            return False

        if self.__verbose:
            self.__log.write(f"Input container list is {[(c.getName(), c.getType()) for c in containerList]!r}\n")

        for c in containerList:
            c.setType('data')

        dBlockList = containerList[1:]

        # the CIF writer skipped a category without rows, as the last one parsed before a syntax error can be
        for dBlock in dBlockList:
            for catName in [catName for catName in dBlock.getObjNameList() if not dBlock.getObj(catName).getRowList()]:
                dBlock.remove(catName)

        if cifPath is not None:
            IoAdapterPy(False, self.__log).writeFile(cifPath, containerList=dBlockList)

        self.__rR = CifReader(self.__verbose, self.__log,
                              use_cache=False)

        return self.__rR.setDataBlockList(dBlockList, strPath, text_md5(strPath) if self.__use_cache else None)

    def __checkNmrDataInputSource(self) -> bool:
        """ Check input source of NMR data.
        """

        if self.__nmrDataPath is not None:
            return True

        if NMR_CIF_FILE_PATH_KEY in self.__inputParamDict:

            fPath = self.__inputParamDict[NMR_CIF_FILE_PATH_KEY]

            if fPath.endswith('.gz'):

                _fPath = os.path.splitext(fPath)[0]

                if not os.path.exists(_fPath):

                    try:

                        uncompress_gzip_file(fPath, _fPath)

                    except Exception as e:  # pylint: disable=broad-exception-caught

                        if self.__verbose:
                            self.__log.write(f"+{self.__class_name__}.__checkNmrDataInputSource() ++ Error  - {str(e)}\n")

                        return False

                fPath = _fPath

            try:

                self.__rR = CifReader(self.__verbose, self.__log,
                                      use_cache=False)

                if self.__rR.parse(fPath):
                    self.__nmrDataPath = fPath

                    if self.__use_cache:
                        self.__nmrDataHashCode = self.__rR.getHashCode()

                    return True

            except Exception:  # pylint: disable=broad-exception-caught
                pass

        if NMR_CIF_READER_OBJ_KEY in self.__inputParamDict:

            self.__rR = self.__inputParamDict[NMR_CIF_READER_OBJ_KEY]

            self.__nmrDataPath = self.__rR.getFilePath()

            if self.__use_cache:
                self.__nmrDataHashCode = self.__rR.getHashCode()

            return True

        if NMR_STR_FILE_PATH_KEY in self.__inputParamDict:

            fPath = self.__inputParamDict[NMR_STR_FILE_PATH_KEY]

            has_safe_cif_path = False

            if NMR_CIF_FILE_PATH_KEY in self.__outputParamDict and self.__outputParamDict[NMR_CIF_FILE_PATH_KEY] is not None:
                _fPath = self.__outputParamDict[NMR_CIF_FILE_PATH_KEY]
                has_safe_cif_path = True
            else:
                _fPath = self.getNextPath(fPath, '.str2cif')

            try:

                if self.__parseNmrStarAsCif(fPath, _fPath if has_safe_cif_path else None):
                    self.__nmrDataPath = fPath

                    if self.__use_cache:
                        self.__nmrDataHashCode = self.__rR.getHashCode()

                    return True

            except Exception as e:  # pylint: disable=broad-exception-caught
                self.__log.write(f"+{self.__class_name__}.__checkNmrDataInputSource() ++ Error  - {str(e)}\n")

            finally:
                try:
                    if not has_safe_cif_path:
                        if os.path.exists(_fPath):
                            os.remove(_fPath)
                except OSError:
                    pass

        if PYNMRSTAR_OBJ_KEY in self.__inputParamDict:

            master_entry = self.__inputParamDict[PYNMRSTAR_OBJ_KEY]

            has_safe_star_path = has_safe_cif_path = False

            if NEXT_STAR_FILE_PATH_KEY in self.__outputParamDict and self.__outputParamDict[NEXT_STAR_FILE_PATH_KEY] is not None:
                _fPath = self.__outputParamDict[NEXT_STAR_FILE_PATH_KEY]
                has_safe_star_path = True
            else:
                _fPath = get_temp_path(None, '.str')

            if NMR_CIF_FILE_PATH_KEY in self.__outputParamDict and self.__outputParamDict[NMR_CIF_FILE_PATH_KEY] is not None:
                __fPath = self.__outputParamDict[NMR_CIF_FILE_PATH_KEY]
                has_safe_cif_path = True
            else:
                __fPath = f'{_fPath}.str2cif'

            try:

                master_entry.write_to_file(_fPath, show_comments=False, skip_empty_loops=True, skip_empty_tags=False)

                if self.__parseNmrStarAsCif(_fPath, __fPath if has_safe_cif_path else None):
                    self.__nmrDataPath = _fPath

                    if self.__use_cache:
                        self.__nmrDataHashCode = self.__rR.getHashCode()

                    return True

            except Exception as e:  # pylint: disable=broad-exception-caught
                self.__log.write(f"+{self.__class_name__}.__checkNmrDataInputSource() ++ Error  - {str(e)}\n")

            finally:
                try:
                    if not has_safe_star_path:
                        if os.path.exists(_fPath):
                            os.remove(_fPath)
                    if not has_safe_cif_path:
                        if os.path.exists(__fPath):
                            os.remove(__fPath)
                except OSError:
                    pass

        return False

    def __parseNmrDpReport(self) -> bool:
        """ Parse report file, which include information about well-defined region and random coil index.
        """

        self.__inputReport = None

        if REPORT_FILE_PATH_KEY not in self.__inputParamDict:
            return True

        fPath = self.__inputParamDict[REPORT_FILE_PATH_KEY]

        try:

            self.__inputReport = NmrDpReport(self.__verbose, self.__log)

            if self.__inputReport.loadFile(fPath):
                return True

        except Exception:  # pylint: disable=broad-exception-caught
            pass

        self.__inputReport = None

        return False

    def __checkPreviousCsAnalysis(self) -> bool:
        """ Retrieve the previous chemical shift analysis using the identical data sources, if available.
        """

        if None not in (self.__cifHashCode, self.__nmrDataHashCode):
            self.__resultsCacheName = f"{self.__cifHashCode}_{self.__nmrDataHashCode}_vrpt_cs_analysis.pkl"
            cache_path = os.path.join(self.__cacheDirPath, self.__resultsCacheName)

            if self.__debug:
                if os.path.exists(cache_path):
                    os.remove(cache_path)

            self.__results = load_from_pickle(cache_path)
            self.__has_prev_results = self.__results is not None

        return True

    def __checkPreviousMrAnalysis(self) -> bool:
        """ Retrieve the previous restraint analysis using the identical data sources, if available.
        """

        if None not in (self.__cifHashCode, self.__nmrDataHashCode):
            self.__resultsCacheName = f"{self.__cifHashCode}_{self.__nmrDataHashCode}_vrpt_mr_analysis.pkl"
            cache_path = os.path.join(self.__cacheDirPath, self.__resultsCacheName)

            if self.__debug:
                if os.path.exists(cache_path):
                    os.remove(cache_path)

            self.__results = load_from_pickle(cache_path)
            self.__has_prev_results = self.__results is not None

        return True

    def __retrieveCoordAssemblyChecker(self) -> bool:
        """ Wrapper function for ParserListenerUtil.coordAssemblyChecker.
        """

        if self.__has_prev_results:
            return True

        if self.__caC is not None:
            return True

        if self.__use_cache:

            cache_path = None
            if self.__cifHashCode is not None:
                cache_path = os.path.join(self.__cacheDirPath, f"{self.__cifHashCode}_asm_chk.pkl")

                self.__caC = load_from_pickle(cache_path)

                if self.__caC is not None:
                    return True

            self.__caC = coordAssemblyChecker(self.__verbose, self.__log,
                                              self.__representative_model_id,
                                              self.__representative_alt_id,
                                              self.__cR, self.__ccU, None, None)

            if self.__caC is not None and cache_path:
                write_as_pickle(self.__caC, cache_path)

        else:
            self.__caC = coordAssemblyChecker(self.__verbose, self.__log,
                                              self.__representative_model_id,
                                              self.__representative_alt_id,
                                              self.__cR, self.__ccU, None, None, False)

        return True

    def __extractEntityInstances(self) -> bool:
        """ Extract entity instances of coordinate file.
            @author: Masashi Yokochi
        """

        if self.__has_prev_results:
            return True

        if self.__coordinates is None:
            return self.__extractEntityInstancesFromNmr()

        vrpt_entity_instance_cache_path = vrpt_entity_uninstance_cache_path = None

        if self.__cifHashCode is not None:
            vrpt_entity_instance_cache_path =\
                os.path.join(self.__cacheDirPath, f"{self.__cifHashCode}_vrpt_entity_instance.pkl")
            self.__entityInstance = load_from_pickle(vrpt_entity_instance_cache_path)
            vrpt_entity_uninstance_cache_path =\
                os.path.join(self.__cacheDirPath, f"{self.__cifHashCode}_vrpt_entity_uninstance.pkl")
            self.__entityUninstance = load_from_pickle(vrpt_entity_uninstance_cache_path)

            if None not in (self.__entityInstance, self.__entityUninstance):
                return True

        _auth_atom_id = 'pdbx_auth_atom_name'\
            if self.__trust_pdbx_auth_atom_name and self.__cR.hasItem('atom_site', 'pdbx_auth_atom_name')\
            else 'auth_atom_id'

        data_items = [{'name': 'auth_asym_id', 'type': 'str', 'alt_name': 'auth_chain_id'},
                      {'name': 'auth_seq_id', 'type': 'int'},
                      {'name': 'auth_comp_id', 'type': 'str'},
                      {'name': _auth_atom_id, 'type': 'str', 'alt_name': 'auth_atom_id'},
                      # {'name': 'label_asym_id', 'type': 'str'},
                      {'name': 'label_seq_id', 'type': 'int', 'alt_name': 'seq_id'},
                      {'name': 'label_comp_id', 'type': 'str'},
                      # {'name': 'label_atom_id', 'type': 'str'},
                      # {'name': 'label_entity_id', 'type': 'int'},
                      {'name': 'pdbx_PDB_ins_code', 'type': 'str', 'alt_name': 'ins_code', 'default': '?'}
                      ]

        NMR_OBS_NUCS = set(ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS.keys())

        # DAOTHER-8705, 8817
        if _auth_atom_id == 'pdbx_auth_atom_name':
            data_items.append({'name': 'auth_atom_id', 'type': 'str', 'alt_name': 'alt_auth_atom_id'})
            data_items.append({'name': 'pdbx_auth_comp_id', 'type': 'str', 'alt_name': 'alt_auth_comp_id'})

        _filter_items = [{'name': 'type_symbol', 'type': 'enum', 'enum': NMR_OBS_NUCS},
                         # remove label_alt_id filter because of missing of ayth_asym_id: B (2n2k)
                         # {'name': 'label_alt_id', 'type': 'enum', 'enum': (self.__representative_alt_id,)},
                         {'name': 'pdbx_PDB_model_num', 'type': 'int', 'value': self.__eff_model_ids[0]}]

        if len(self.__caC['polymer_sequence']) >= LEN_MAJOR_ASYM_ID:
            _filter_items.append({'name': 'auth_asym_id', 'type': 'enum', 'enum': LARGE_ASYM_ID,
                                  'fetch_first_match': True})  # to process large assembly avoiding forced timeout

        self.__entityInstance, self.__entityUninstance = {}, {}

        atom_name_unchecked = comp_name_unchecked = True
        _auth_atom_id_ = 'auth_atom_id'

        try:

            coord = self.__cR.getDictListWithFilter('atom_site', data_items, _filter_items)

            # DAOTHER-8705
            if atom_name_unchecked:
                atom_name_unchecked = False
                if _auth_atom_id == 'pdbx_auth_atom_name':
                    for c in coord:
                        if None not in (c['auth_atom_id'], c['alt_auth_atom_id'])\
                           and c['auth_atom_id'][0].isdigit() and c['alt_auth_atom_id'][0] == 'H':
                            _auth_atom_id_ = 'alt_auth_atom_id'
                            break

            # DAOTHER-8817
            if comp_name_unchecked:
                comp_name_unchecked = False
                if _auth_atom_id == 'pdbx_auth_atom_name':
                    for c in coord:
                        if c['auth_comp_id'] is not None and c['auth_comp_id'] not in STD_MON_DICT\
                           and c['alt_auth_comp_id'] is not None and c['auth_comp_id'] != c['alt_auth_comp_id']:
                            _auth_atom_id_ = 'alt_auth_atom_id'
                            break

            for c in coord:
                auth_chain_id = c['auth_chain_id']
                if auth_chain_id not in self.__entityInstance:
                    self.__entityInstance[auth_chain_id] = {}

                seq_id = str(c['auth_seq_id'])
                if c['ins_code'] not in EMPTY_VALUE:
                    seq_id += c['ins_code']

                seq_key = (seq_id, c['label_comp_id'])

                if seq_key not in self.__entityInstance[auth_chain_id]:
                    self.__entityInstance[auth_chain_id][seq_key] = {'seq_id': c['seq_id'], 'atoms': []}

                atom_id = c[_auth_atom_id_]

                if atom_id not in self.__entityInstance[auth_chain_id][seq_key]['atoms']:
                    self.__entityInstance[auth_chain_id][seq_key]['atoms'].append(atom_id)

            # include unmodeled entity

            def update_entity_instance(_ps):
                auth_chain_id = _ps['auth_chain_id']
                if auth_chain_id not in self.__entityInstance:
                    return

                if 'ins_code' in _ps:
                    for auth_seq_id, ins_code, seq_id, comp_id\
                            in zip(_ps['auth_seq_id'], _ps['ins_code'], _ps['seq_id'], _ps['comp_id']):
                        auth_seq_id = str(auth_seq_id)
                        if ins_code not in EMPTY_VALUE:
                            auth_seq_id += ins_code

                        seq_key = (auth_seq_id, comp_id)

                        if seq_key in self.__entityInstance[auth_chain_id]:
                            continue

                        self.__entityInstance[auth_chain_id][seq_key] = {'seq_id': seq_id, 'atoms': []}

                        if self.__ccU.updateChemCompDict(comp_id):
                            peptide_like = self.__csStat.peptideLike(comp_id)
                            for cca in self.__ccU.lastAtomDictList:
                                if cca['type_symbol'] in NMR_OBS_NUCS\
                                   and (cca['leaving_atom_flag'] != 'Y'
                                        or (peptide_like
                                            and cca['n_terminal_atom_flag'] == 'N'
                                            and cca['c_terminal_atom_flag'] == 'N')):
                                    self.__entityInstance[auth_chain_id][seq_key]['atoms'].append(cca['atom_id'])

                else:
                    for auth_seq_id, seq_id, comp_id\
                            in zip(_ps['auth_seq_id'], _ps['seq_id'], _ps['comp_id']):
                        seq_key = (str(auth_seq_id), comp_id)

                        if seq_key in self.__entityInstance[auth_chain_id]:
                            continue

                        self.__entityInstance[auth_chain_id][seq_key] = {'seq_id': seq_id, 'atoms': []}

                        if self.__ccU.updateChemCompDict(comp_id):
                            peptide_like = self.__csStat.peptideLike(comp_id)
                            for cca in self.__ccU.lastAtomDictList:
                                if cca['type_symbol'] in NMR_OBS_NUCS\
                                   and (cca['leaving_atom_flag'] != 'Y'
                                        or (peptide_like
                                            and cca['n_terminal_atom_flag'] == 'N'
                                            and cca['c_terminal_atom_flag'] == 'N')):
                                    self.__entityInstance[auth_chain_id][seq_key]['atoms'].append(cca['atom_id'])

            for ps in self.__caC['polymer_sequence']:
                update_entity_instance(ps)

            if 'branched' in self.__caC and self.__caC['branched'] is not None:
                for br in self.__caC['branched']:
                    update_entity_instance(br)

            if 'non_poly' in self.__caC and self.__caC['non_poly'] is not None:
                for np in self.__caC['non_poly']:
                    update_entity_instance(np)

            def get_label_seq_id(auth_chain_id, auth_seq_id):
                ps = next((ps for ps in self.__caC['polymer_sequence']
                           if ps['auth_chain_id'] == auth_chain_id and auth_seq_id in ps['auth_seq_id']), None)
                if ps is not None:
                    return ps['seq_id'][ps['auth_seq_id'].index(auth_seq_id)]
                if 'branched' in self.__caC and self.__caC['branched'] is not None:
                    br = next((br for br in self.__caC['branched']
                               if br['auth_chain_id'] == auth_chain_id and auth_seq_id in br['auth_seq_id']), None)
                    if br is not None:
                        return br['seq_id'][br['auth_seq_id'].index(auth_seq_id)]
                if 'non_poly' in self.__caC and self.__caC['non_poly'] is not None:
                    np = next((np for np in self.__caC['non_poly']
                               if np['auth_chain_id'] == auth_chain_id and auth_seq_id in np['auth_seq_id']), None)
                    if np is not None:
                        return np['seq_id'][np['auth_seq_id'].index(auth_seq_id)]
                return None

            # DAOTHER-10987
            if 'coord_unobs_res' in self.__caC and self.__caC['coord_unobs_res'] is not None:
                for k, v in self.__caC['coord_unobs_res'].items():
                    auth_chain_id, auth_seq_id, comp_id = k[0], k[1], v['comp_id']
                    if comp_id in EMPTY_VALUE:
                        continue
                    if auth_chain_id not in self.__entityUninstance:
                        self.__entityUninstance[auth_chain_id] = {}
                    seq_key = (str(auth_seq_id), comp_id)
                    ps = next((ps for ps in self.__caC['polymer_sequence']
                               if ps['auth_seq_id'] == auth_chain_id and auth_seq_id in ps['auth_seq_id']), None)
                    seq_id = get_label_seq_id(auth_chain_id, auth_seq_id)
                    if seq_id is None or not self.__ccU.updateChemCompDict(comp_id):
                        continue
                    self.__entityUninstance[auth_chain_id][seq_key] =\
                        {'seq_id': seq_id,
                         'atoms': [cca['atom_id'] for cca in self.__ccU.lastAtomDictList
                                   if cca['leaving_atom_flag'] != 'Y'
                                   or (self.__csStat.peptideLike(comp_id)
                                       and cca['n_terminal_atom_flag'] == 'N'
                                       and cca['c_terminal_atom_flag'] == 'N')]}

            # DAOTHER-10987
            if 'coord_unobs_atom' in self.__caC and self.__caC['coord_unobs_atom'] is not None:
                for k, v in self.__caC['coord_unobs_atom'].items():
                    auth_chain_id, auth_seq_id, comp_id, atoms = k[0], k[1], v['comp_id'], v['atom_ids']
                    if comp_id in EMPTY_VALUE or not isinstance(atoms, list) or len(atoms) == 0:
                        continue
                    if auth_chain_id not in self.__entityUninstance:
                        self.__entityUninstance[auth_chain_id] = {}
                    seq_key = (str(auth_seq_id), comp_id)
                    seq_id = get_label_seq_id(auth_chain_id, auth_seq_id)
                    if seq_id is None or not self.__ccU.updateChemCompDict(comp_id):
                        continue
                    if seq_key not in self.__entityUninstance[auth_chain_id]:
                        self.__entityUninstance[auth_chain_id][seq_key] =\
                            {'seq_id': seq_id, 'atoms': atoms}
                    else:
                        self.__entityUninstance[auth_chain_id][seq_key]['atoms'].extend(atoms)

            if self.__cifHashCode is not None:
                write_as_pickle(self.__entityInstance, vrpt_entity_instance_cache_path)
                write_as_pickle(self.__entityUninstance, vrpt_entity_uninstance_cache_path)

            return True

        except Exception as e:  # pylint: disable=broad-exception-caught
            self.__log.write(f"Exception occurred while processing {os.path.basename(self.__cifPath)} "
                             f"and {os.path.basename(self.__nmrDataPath)}\n")
            self.__log.write(f"+{self.__class_name__}.__extractEntityInstances() ++ Error  - {str(e)}\n")

            self.__entityInstance = self.__entityUninstance = None

        return False

    def __extractEntityInstancesFromNmr(self) -> bool:
        """ Extract entity instances of NMR data file.
            @author: Masashi Yokochi
        """

        self.__entityInstance, self.__entityUninstance = {}, {}

        NMR_OBS_NUCS = set(ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS.keys())

        try:

            has_chem_comp_assembly = False

            lp_category = 'Chem_comp_assembly'

            for datablock_name in self.__rR.getDataBlockNameList():

                if not self.__rR.hasCategory(lp_category, datablock_name):
                    continue

                has_chem_comp_assembly = True

                for row in self.__rR.getDictList(lp_category, datablock_name):
                    auth_chain_id = row['Entity_assembly_ID']
                    seq_id = row['Comp_index_ID']
                    # NEF sequence translation provides canonical Comp_index_ID
                    # even when the redundant Seq_ID column is absent.
                    auth_seq_id = row.get('Seq_ID', seq_id)
                    if auth_seq_id in EMPTY_VALUE:
                        auth_seq_id = seq_id
                    comp_id = row['Comp_ID']

                    seq_key = (auth_seq_id, comp_id)

                    if auth_chain_id not in self.__entityInstance:
                        self.__entityInstance[auth_chain_id] = {}

                    if seq_key not in self.__entityInstance[auth_chain_id]:
                        self.__entityInstance[auth_chain_id][seq_key] = {'seq_id': seq_id, 'atoms': []}

                    if self.__ccU.updateChemCompDict(comp_id):
                        peptide_like = self.__csStat.peptideLike(comp_id)
                        for cca in self.__ccU.lastAtomDictList:
                            if cca['type_symbol'] in NMR_OBS_NUCS\
                               and (cca['leaving_atom_flag'] != 'Y'
                                    or (peptide_like
                                        and cca['n_terminal_atom_flag'] == 'N'
                                        and cca['c_terminal_atom_flag'] == 'N')):
                                self.__entityInstance[auth_chain_id][seq_key]['atoms'].append(cca['atom_id'])

                    if 'Auth_variant_ID' in row:
                        variant_id = row['Auth_variant_ID']

                        if variant_id in EMPTY_VALUE:
                            continue

                        variants = variant_id.split(',')

                        for v in variants:
                            v = v.strip()

                            if not v.startswith('-'):
                                continue

                            if auth_chain_id not in self.__entityUninstance:
                                self.__entityUninstance[auth_chain_id] = {}

                            if seq_key not in self.__entityUninstance[auth_chain_id]:
                                self.__entityUninstance[auth_chain_id][seq_key] = {'seq_id': seq_id, 'atoms': []}

                            self.__entityUninstance[auth_chain_id][seq_key]['atoms'].append(v[1:])

            if not has_chem_comp_assembly:

                lp_category = 'Entity_assembly'

                entity_id_mapping = {}

                for datablock_name in self.__rR.getDataBlockNameList():

                    if not self.__rR.hasCategory(lp_category, datablock_name):
                        continue

                    for row in self.__rR.getDictList(lp_category, datablock_name):
                        auth_chain_id = row['ID']
                        entity_id = row['Entity_ID']

                        if auth_chain_id in EMPTY_VALUE or entity_id in EMPTY_VALUE:
                            continue

                        entity_id_mapping[entity_id] = auth_chain_id

                if len(entity_id_mapping) > 0:

                    lp_category = 'Entity_comp_index'

                    for datablock_name in self.__rR.getDataBlockNameList():

                        if not self.__rR.hasCategory(lp_category, datablock_name):
                            continue

                        for row in self.__rR.getDictList(lp_category, datablock_name):
                            seq_id = row['ID']
                            auth_seq_id = row['Auth_seq_ID']
                            comp_id = row['Comp_ID']
                            entity_id = row['Entity_ID']

                            if auth_seq_id in EMPTY_VALUE or not auth_seq_id.isdigit():
                                auth_seq_id = seq_id

                            if entity_id in entity_id_mapping:
                                auth_chain_id = entity_id_mapping[entity_id]

                                if auth_chain_id not in self.__entityInstance:
                                    self.__entityInstance[auth_chain_id] = {}

                                seq_key = (auth_seq_id, comp_id)

                                if seq_key not in self.__entityInstance[auth_chain_id]:
                                    self.__entityInstance[auth_chain_id][seq_key] = {'seq_id': seq_id, 'atoms': []}

                                if self.__ccU.updateChemCompDict(comp_id):
                                    peptide_like = self.__csStat.peptideLike(comp_id)
                                    for cca in self.__ccU.lastAtomDictList:
                                        if cca['type_symbol'] in NMR_OBS_NUCS\
                                           and (cca['leaving_atom_flag'] != 'Y'
                                                or (peptide_like
                                                    and cca['n_terminal_atom_flag'] == 'N'
                                                    and cca['c_terminal_atom_flag'] == 'N')):
                                            self.__entityInstance[auth_chain_id][seq_key]['atoms'].append(cca['atom_id'])

            lp_category = 'Entity_deleted_atom'

            for datablock_name in self.__rR.getDataBlockNameList():

                if not self.__rR.hasCategory(lp_category, datablock_name):
                    continue

                for row in self.__rR.getDictList(lp_category, datablock_name):
                    auth_chain_id = row['Entity_assembly_ID']
                    seq_id = row['Comp_index_ID']
                    auth_seq_id = row['Seq_ID']
                    comp_id = row['Comp_ID']
                    atom_id = row['Atom_ID']

                    if auth_seq_id in EMPTY_VALUE or not auth_seq_id.isdigit():
                        auth_seq_id = seq_id

                    if auth_chain_id not in self.__entityUninstance:
                        self.__entityUninstance[auth_chain_id] = {}

                    seq_key = (auth_seq_id, comp_id)

                    if seq_key not in self.__entityUninstance[auth_chain_id]:
                        self.__entityUninstance[auth_chain_id][seq_key] = {'seq_id': seq_id, 'atoms': []}

                    self.__entityUninstance[auth_chain_id][seq_key]['atoms'].append(atom_id)

            return len(self.__entityInstance) > 0

        except Exception as e:  # pylint: disable=broad-exception-caught
            self.__log.write(f"Exception occurred while processing {os.path.basename(self.__nmrDataPath)}\n")
            self.__log.write(f"+{self.__class_name__}.__extractEntityInstancesFromNmr() ++ Error  - {str(e)}\n")

            self.__entityInstance = self.__entityUninstance = None

        return False

    def __extractCoordAtomSites(self) -> bool:
        """ Extract atom_site of coordinate file.
            @author: Masashi Yokochi
            @note: Derived from wwpdb.apps.validation.src.RestraintsValidation.BMRBRestraintsAnalysis.get_coordinates,
                   written by Kumaran Baskaran
            @change: class method, use of wwpdb.utils.nmr.io.CifReader, use PDB_ins_code for atom identification,
                     performance optimization
        """

        if self.__has_prev_results:
            return True

        if self.__cifPath is None:
            return True

        vrpt_atom_id_cache_path = vrpt_atom_site_cache_path = None

        if self.__cifHashCode is not None:
            vrpt_atom_id_cache_path = os.path.join(self.__cacheDirPath, f"{self.__cifHashCode}_vrpt_atom_id.pkl")
            self.__atomIdList = load_from_pickle(vrpt_atom_id_cache_path)

            vrpt_atom_site_cache_path = os.path.join(self.__cacheDirPath, f"{self.__cifHashCode}_vrpt_atom_site.pkl")
            self.__coordinates = load_from_pickle(vrpt_atom_site_cache_path)

            if None not in (self.__atomIdList, self.__coordinates):
                return True

        _auth_atom_id = 'pdbx_auth_atom_name'\
            if self.__trust_pdbx_auth_atom_name and self.__cR.hasItem('atom_site', 'pdbx_auth_atom_name')\
            else 'auth_atom_id'

        data_items = [{'name': 'auth_asym_id', 'type': 'str'},
                      {'name': 'auth_seq_id', 'type': 'int'},
                      {'name': 'auth_comp_id', 'type': 'str'},
                      {'name': _auth_atom_id, 'type': 'str', 'alt_name': 'auth_atom_id'},
                      {'name': 'label_asym_id', 'type': 'str'},
                      {'name': 'label_seq_id', 'type': 'int'},
                      {'name': 'label_comp_id', 'type': 'str'},
                      # {'name': 'label_atom_id', 'type': 'str'},
                      {'name': 'label_entity_id', 'type': 'int'},
                      {'name': 'label_alt_id', 'type': 'str', 'default': '.'},
                      {'name': 'pdbx_PDB_ins_code', 'type': 'str', 'alt_name': 'ins_code', 'default': '?'},
                      {'name': 'Cartn_x', 'type': 'float', 'alt_name': 'x'},
                      {'name': 'Cartn_y', 'type': 'float', 'alt_name': 'y'},
                      {'name': 'Cartn_z', 'type': 'float', 'alt_name': 'z'}
                      ]

        # DAOTHER-8705, 8817
        if _auth_atom_id == 'pdbx_auth_atom_name':
            data_items.append({'name': 'auth_atom_id', 'type': 'str', 'alt_name': 'alt_auth_atom_id'})
            data_items.append({'name': 'pdbx_auth_comp_id', 'type': 'str', 'alt_name': 'alt_auth_comp_id'})

        _filter_items = []  # {'name': 'label_alt_id', 'type': 'enum', 'enum': (self.__representative_alt_id,)}]

        if len(self.__caC['polymer_sequence']) >= LEN_MAJOR_ASYM_ID:
            _filter_items.append({'name': 'auth_asym_id', 'type': 'enum', 'enum': LARGE_ASYM_ID,
                                  'fetch_first_match': True})  # to process large assembly avoiding forced timeout

        self.__atomIdList, self.__coordinates = {}, {}

        atom_name_unchecked = comp_name_unchecked = True
        _auth_atom_id_ = 'auth_atom_id'

        try:

            for model_id in self.__eff_model_ids:
                filter_items = copy.copy(_filter_items)
                filter_items.append({'name': 'pdbx_PDB_model_num', 'type': 'int',
                                     'value': model_id})

                coord = self.__cR.getDictListWithFilter('atom_site', data_items, filter_items)

                # DAOTHER-8705
                if atom_name_unchecked:
                    atom_name_unchecked = False
                    if _auth_atom_id == 'pdbx_auth_atom_name':
                        for c in coord:
                            if None not in (c['auth_atom_id'], c['alt_auth_atom_id'])\
                               and c['auth_atom_id'][0].isdigit() and c['alt_auth_atom_id'][0] == 'H':
                                _auth_atom_id_ = 'alt_auth_atom_id'
                                break

                # DAOTHER-8817
                if comp_name_unchecked:
                    comp_name_unchecked = False
                    if _auth_atom_id == 'pdbx_auth_atom_name':
                        for c in coord:
                            if c['auth_comp_id'] is not None and c['auth_comp_id'] not in STD_MON_DICT\
                               and c['alt_auth_comp_id'] is not None and c['auth_comp_id'] != c['alt_auth_comp_id']:
                                _auth_atom_id_ = 'alt_auth_atom_id'
                                break

                atom_id_list_per_model, coordinates_per_model = {}, {}

                for c in coord:
                    atom_key = (c['auth_asym_id'], c['auth_seq_id'], c['auth_comp_id'],
                                c[_auth_atom_id_], c['ins_code'])

                    # tokens = ("ent_", "said_", "resname_", "seq_", "resnum_", "altcode_", "icode_", "chain_")

                    atom_id_list_per_model[atom_key] =\
                        (c['label_entity_id'], c['label_asym_id'], c['label_comp_id'], c['label_seq_id'],
                         c['auth_seq_id'], c['label_alt_id'], c['ins_code'], c['auth_asym_id'])

                    coordinates_per_model[atom_key] = to_np_array(c)

                    # DAOTHER-9200 for wwpdb.apps.validation.src.wrapper.restraintsanalysis.generate_formated_output
                    # (MISSING ATOM IN MODEL KeyError)
                    if _auth_atom_id == 'pdbx_auth_atom_name' and c['auth_atom_id'] != c['alt_auth_atom_id']:
                        _atom_key = (c['auth_asym_id'], c['auth_seq_id'], c['auth_comp_id'],
                                     c['alt_auth_atom_id'], c['ins_code'])

                        atom_id_list_per_model[_atom_key] =\
                            (c['label_entity_id'], c['label_asym_id'], c['label_comp_id'], c['label_seq_id'],
                             c['auth_seq_id'], c['label_alt_id'], c['ins_code'], c['auth_asym_id'])

                        coordinates_per_model[_atom_key] = to_np_array(c)

                self.__atomIdList[model_id] = atom_id_list_per_model
                self.__coordinates[model_id] = coordinates_per_model

            if self.__cifHashCode is not None:
                write_as_pickle(self.__atomIdList, vrpt_atom_id_cache_path)
                write_as_pickle(self.__coordinates, vrpt_atom_site_cache_path)

            return True

        except Exception as e:  # pylint: disable=broad-exception-caught
            self.__log.write(f"Exception occurred while processing {os.path.basename(self.__cifPath)} "
                             f"and {os.path.basename(self.__nmrDataPath)}\n")
            self.__log.write(f"+{self.__class_name__}.__extractCoordAtomSites() ++ Error  - {str(e)}\n")

            self.__atomIdList = self.__coordinates = None

        return False

    def __extractChemicalShifts(self) -> bool:
        """ Extract Atom_chem_shift category of NMR data file.
            @author: Masashi Yokochi
            @note: Derived from wwpdb.apps.validation.src.ChemicalShiftsValidation.BMRBChemicalShiftAnalysis.get_chemical_shifts,
                   written by Aleksandras Gutmanas, Kumaran Baskaran
            @change: class method, use of wwpdb.utils.nmr.io.CifReader, improve readability of chemical shifts
        """

        if self.__has_prev_results:
            return True

        if self.__nmrDataPath is None:
            return False

        has_coord = self.__coordinates is not None

        self.__chemShiftMeta = []
        self.__chemShiftTotal, self.__chemShiftDict, self.__chemShiftUnparsed = {}, {}, {}

        lp_category = 'Atom_chem_shift'
        sf_category = 'Assigned_chem_shift_list'

        try:

            for idx, datablock_name in enumerate(self.__rR.getDataBlockNameList()):

                if not self.__rR.hasCategory(lp_category, datablock_name):
                    continue

                sf_tag = self.__rR.getDictList(sf_category, datablock_name)

                list_id = int(sf_tag[0]['ID'])
                sf_framecode = sf_tag[0]['Sf_framecode']

                self.__chemShiftMeta.append((idx, sf_framecode, list_id))
                self.__chemShiftTotal[list_id] = self.__rR.getRowLength(lp_category, datablock_name)

                if has_coord:
                    data_items = [{'name': 'ID', 'type': 'int', 'alt_name': 'id'},
                                  {'name': 'Auth_asym_ID', 'type': 'str', 'alt_name': 'auth_chain_id'},
                                  {'name': 'Auth_seq_ID', 'type': 'int', 'alt_name': 'auth_seq_id'},
                                  {'name': 'Comp_ID', 'type': 'str', 'alt_name': 'comp_id'},
                                  {'name': 'Atom_ID', 'type': 'str', 'alt_name': 'atom_id'},
                                  {'name': 'Val', 'type': 'float', 'alt_name': 'value'},
                                  {'name': 'Val_err', 'type': 'float', 'alt_name': 'error'},
                                  {'name': 'Ambiguity_code', 'type': 'enum-int', 'alt_name': 'ambig_code',
                                   'enum': ALLOWED_AMBIGUITY_CODES},
                                  ]
                else:
                    data_items = [{'name': 'ID', 'type': 'int', 'alt_name': 'id'},
                                  {'name': 'Entity_assembly_ID', 'type': 'str', 'alt_name': 'auth_chain_id'},
                                  {'name': 'Comp_index_ID', 'type': 'int', 'alt_name': 'auth_seq_id'},
                                  {'name': 'Comp_ID', 'type': 'str', 'alt_name': 'comp_id'},
                                  {'name': 'Atom_ID', 'type': 'str', 'alt_name': 'atom_id'},
                                  {'name': 'Val', 'type': 'float', 'alt_name': 'value'},
                                  {'name': 'Val_err', 'type': 'float', 'alt_name': 'error'},
                                  {'name': 'Ambiguity_code', 'type': 'enum-int', 'alt_name': 'ambig_code',
                                   'enum': ALLOWED_AMBIGUITY_CODES},
                                  ]

                tags = self.__rR.getAttributeList(lp_category, datablock_name)

                if 'PDB_ins_code' in tags:
                    data_items.append({'name': 'PDB_ins_code', 'type': 'str', 'alt_name': 'ins_code', 'default': '?'})

                filter_items = [{'name': 'Assigned_chem_shift_list_ID', 'type': 'int', 'value': list_id}]

                data = self.__rR.getDictListWithFilter(lp_category,
                                                       data_items,
                                                       filter_items,
                                                       datablock_name)

                if list_id not in self.__chemShiftDict:

                    if len(data) == 0:
                        self.__log.write(f"Failed in parsing assigned chemical shifts of list ID {list_id}\n")
                        continue

                    self.__chemShiftDict[list_id] = data
                    self.__chemShiftUnparsed[list_id] = []

                    if self.__chemShiftTotal[list_id] > len(data):

                        if has_coord:
                            data_items = [{'name': 'ID', 'type': 'int', 'alt_name': 'id'},
                                          {'name': 'Auth_asym_ID', 'type': 'str', 'alt_name': 'auth_chain_id'},
                                          {'name': 'Auth_seq_ID', 'type': 'int', 'alt_name': 'auth_seq_id'},
                                          {'name': 'Comp_ID', 'type': 'str', 'alt_name': 'comp_id'},
                                          {'name': 'Atom_ID', 'type': 'str', 'alt_name': 'atom_id'},
                                          {'name': 'Val', 'type': 'str', 'alt_name': 'value'},
                                          {'name': 'Val_err', 'type': 'str', 'alt_name': 'error'},
                                          {'name': 'Ambiguity_code', 'type': 'str', 'alt_name': 'ambig_code'}
                                          ]
                        else:
                            data_items = [{'name': 'ID', 'type': 'int', 'alt_name': 'id'},
                                          {'name': 'Entity_assembly_ID', 'type': 'str', 'alt_name': 'auth_chain_id'},
                                          {'name': 'Comp_index_ID', 'type': 'int', 'alt_name': 'auth_seq_id'},
                                          {'name': 'Comp_ID', 'type': 'str', 'alt_name': 'comp_id'},
                                          {'name': 'Atom_ID', 'type': 'str', 'alt_name': 'atom_id'},
                                          {'name': 'Val', 'type': 'str', 'alt_name': 'value'},
                                          {'name': 'Val_err', 'type': 'str', 'alt_name': 'error'},
                                          {'name': 'Ambiguity_code', 'type': 'str', 'alt_name': 'ambig_code'}
                                          ]

                        if 'PDB_ins_code' in tags:
                            data_items.append({'name': 'PDB_ins_code', 'type': 'str', 'alt_name': 'ins_code', 'default': '?'})

                        filter_items = [{'name': 'Assigned_chem_shift_list_ID', 'type': 'int', 'value': list_id}]

                        _data = self.__rR.getDictListWithFilter(lp_category,
                                                                data_items,
                                                                filter_items,
                                                                datablock_name)

                        # Prebuild the parsed (chain, seq, comp, atom) keys once so the
                        # "already parsed?" test below is O(1) instead of a full linear
                        # scan of `data` per unparsed row (was O(len(_data) * len(data))).
                        data_keys = {(cs['auth_chain_id'], cs['auth_seq_id'],
                                      cs['comp_id'], cs['atom_id']) for cs in data}

                        for _cs in _data:
                            cs_auth_chain_id = _cs['auth_chain_id']
                            cs_auth_seq_id = _cs['auth_seq_id']
                            cs_comp_id = _cs['comp_id']
                            cs_atom_id = _cs['atom_id']

                            # DAOTHER-10898
                            if ',' in cs_auth_chain_id:
                                for cs_auth_chain_id in _cs['auth_chain_id'].split(','):
                                    cs_auth_chain_id = cs_auth_chain_id.strip()
                                    if (cs_auth_chain_id, cs_auth_seq_id, cs_comp_id, cs_atom_id) not in data_keys:
                                        cs_value = _cs['value']
                                        cs_error = _cs['error']
                                        ambig_code = _cs['ambig_code']

                                        if 'PDB_ins_code' not in tags or _cs['ins_code'] in EMPTY_VALUE:
                                            cs_key = (cs_auth_chain_id, str(cs_auth_seq_id), cs_comp_id, cs_atom_id)
                                        else:
                                            cs_key = (cs_auth_chain_id, str(cs_auth_seq_id) + _cs['ins_code'],
                                                      cs_comp_id, cs_atom_id)

                                        err_cs_values = list(cs_key)
                                        err_cs_values.extend([cs_value, cs_error, ambig_code])
                                        self.__chemShiftUnparsed[list_id].append(err_cs_values)

                            else:
                                if (cs_auth_chain_id, cs_auth_seq_id, cs_comp_id, cs_atom_id) not in data_keys:
                                    cs_value = _cs['value']
                                    cs_error = _cs['error']
                                    ambig_code = _cs['ambig_code']

                                    if 'PDB_ins_code' not in tags or _cs['ins_code'] in EMPTY_VALUE:
                                        cs_key = (cs_auth_chain_id, str(cs_auth_seq_id), cs_comp_id, cs_atom_id)
                                    else:
                                        cs_key = (cs_auth_chain_id, str(cs_auth_seq_id) + _cs['ins_code'],
                                                  cs_comp_id, cs_atom_id)

                                    err_cs_values = list(cs_key)
                                    err_cs_values.extend([cs_value, cs_error, ambig_code])
                                    self.__chemShiftUnparsed[list_id].append(err_cs_values)

            return True

        except Exception as e:  # pylint: disable=broad-exception-caught
            self.__log.write(f"Exception occurred while processing {os.path.basename(self.__cifPath)} "
                             f"and {os.path.basename(self.__nmrDataPath)}\n")
            self.__log.write(f"+{self.__class_name__}.__extractChemicalShiftgs() ++ Error  - {str(e)}\n")

            self.__chemShiftMeta = self.__chemShiftTotal = self.__chemShiftDict = self.__chemShiftUnparsed = None

        return False

    def __extractGenDistConstraints(self) -> bool:
        """ Extract Gen_dist_constraint category of NMR data file.
            @author: Masashi Yokochi
            @note: Derived from wwpdb.apps.validation.src.RestraintsValidation.BMRBRestraintsAnalysis.get_restraints2,
                   written by Kumaran Baskaran
            @change: class method, use of wwpdb.utils.nmr.io.CifReader, improve readability of restraints,
                     support combinational restraints (_Gen_dist_constraint.Combination_ID, Member_ID)
                     define bond types for diselenide bond and metal coordination (sebond, metal)
                     ignore obvious lower-bound-* potential restraint (e.g. CYAYA .lol file)
        """

        if self.__has_prev_results:
            return True

        if self.__nmrDataPath is None:
            return False

        self.__distRestDict, self.__distRestSeqDict = {}, {}

        # {comp_id: frozenset(backbone atom_ids)} cache — getBackBoneAtoms() rebuilds
        # a list (and may touch the CCD/CSV) on every call, and comp_id repeats across
        # restraints; cache once and use set membership below (was O(restraints * bb_len)).
        bb_atoms_cache = {}

        lp_category = 'Gen_dist_constraint'
        sf_category = 'Gen_dist_constraint_list'

        try:

            has_upper_bound = False

            for datablock_name in self.__rR.getDataBlockNameList():

                if not self.__rR.hasCategory(lp_category, datablock_name):
                    continue

                sf_tag = self.__rR.getDictList(sf_category, datablock_name)

                try:
                    potential_type = sf_tag[0]['Potential_type']
                    if potential_type not in EMPTY_VALUE:
                        if 'upper-bound' in potential_type\
                           or 'square-well' in potential_type\
                           or 'log-harmonic' in potential_type\
                           or potential_type == 'parabolic':
                            has_upper_bound = True
                        break
                except KeyError:
                    pass

            skipped = False

            for datablock_name in self.__rR.getDataBlockNameList():

                if not self.__rR.hasCategory(lp_category, datablock_name):
                    continue

                sf_tag = self.__rR.getDictList(sf_category, datablock_name)

                list_id = int(sf_tag[0]['ID'])

                try:
                    constraint_type = sf_tag[0]['Constraint_type']
                    if constraint_type == 'covalent bond':
                        continue
                except KeyError:
                    pass

                try:
                    potential_type = sf_tag[0]['Potential_type']
                    if potential_type not in EMPTY_VALUE and 'lower-bound' in potential_type and has_upper_bound:
                        continue
                except KeyError:
                    pass

                data_items = [{'name': 'ID', 'type': 'int', 'alt_name': 'id'},
                              {'name': 'Entity_assembly_ID_1', 'type': 'str', 'alt_name': 'entity_asm_id_1'},
                              {'name': 'Auth_asym_ID_1', 'type': 'str', 'alt_name': 'auth_asym_id_1'},
                              {'name': 'Auth_seq_ID_1', 'type': 'int', 'alt_name': 'auth_seq_id_1'},
                              {'name': 'Comp_ID_1', 'type': 'str', 'alt_name': 'comp_id_1'},
                              {'name': 'Atom_ID_1', 'type': 'str', 'alt_name': 'atom_id_1'},
                              {'name': 'Entity_assembly_ID_2', 'type': 'str', 'alt_name': 'entity_asm_id_2'},
                              {'name': 'Auth_asym_ID_2', 'type': 'str', 'alt_name': 'auth_asym_id_2'},
                              {'name': 'Auth_seq_ID_2', 'type': 'int', 'alt_name': 'auth_seq_id_2'},
                              {'name': 'Comp_ID_2', 'type': 'str', 'alt_name': 'comp_id_2'},
                              {'name': 'Atom_ID_2', 'type': 'str', 'alt_name': 'atom_id_2'},
                              {'name': 'Distance_lower_bound_val', 'type': 'float', 'alt_name': 'lower_limit'},
                              {'name': 'Distance_upper_bound_val', 'type': 'float', 'alt_name': 'upper_limit'}
                              ]

                tags = self.__rR.getAttributeList(lp_category, datablock_name)

                has_combination_id = 'Combination_ID' in tags
                has_member_logic_code = 'Member_logic_code' in tags
                has_member_id = 'Member_ID' in tags
                has_pdb_ins_code_1 = 'PDB_ins_code_1' in tags
                has_pdb_ins_code_2 = 'PDB_ins_code_2' in tags
                has_target_val = 'Target_val' in tags
                has_target_val_uncertainty = 'Target_val_uncertainty' in tags
                has_lower_linear_limit = 'Lower_linear_limit' in tags
                has_upper_linear_limit = 'Upper_linear_limit' in tags

                if has_combination_id:
                    data_items.append({'name': 'Combination_ID', 'type': 'int', 'alt_name': 'combination_id'})
                if has_member_logic_code:
                    data_items.append({'name': 'Member_logic_code', 'type': 'enum', 'alt_name': 'member_logic_code',
                                       'enum': ('OR', 'AND')})
                if has_member_id:
                    data_items.append({'name': 'Member_ID', 'type': 'int', 'alt_name': 'member_id'})
                if has_pdb_ins_code_1:
                    data_items.append({'name': 'PDB_ins_code_1', 'type': 'str', 'alt_name': 'ins_code_1', 'default': '?'})
                if has_pdb_ins_code_2:
                    data_items.append({'name': 'PDB_ins_code_2', 'type': 'str', 'alt_name': 'ins_code_2', 'default': '?'})
                if has_target_val:
                    data_items.append({'name': 'Target_val', 'type': 'float', 'alt_name': 'target_value'})
                if has_target_val_uncertainty:
                    data_items.append({'name': 'Target_val_uncertainty', 'type': 'abs-float',
                                       'alt_name': 'target_value_uncertainty'})
                if has_lower_linear_limit:
                    data_items.append({'name': 'Lower_linear_limit', 'type': 'float', 'alt_name': 'lower_linear_limit'})
                if has_upper_linear_limit:
                    data_items.append({'name': 'Upper_linear_limit', 'type': 'float', 'alt_name': 'upper_linear_limit'})

                filter_items = [{'name': 'Gen_dist_constraint_list_ID', 'type': 'int', 'value': list_id}]

                rest = self.__rR.getDictListWithFilter(lp_category,
                                                       data_items,
                                                       filter_items,
                                                       datablock_name)

                for r in rest:
                    rest_key = (list_id, r['id'])

                    if rest_key not in self.__distRestDict:
                        self.__distRestDict[rest_key] = []

                    auth_asym_id_1 = r['auth_asym_id_1']
                    auth_asym_id_2 = r['auth_asym_id_2']
                    auth_seq_id_1 = r['auth_seq_id_1']
                    auth_seq_id_2 = r['auth_seq_id_2']
                    comp_id_1 = r['comp_id_1']
                    comp_id_2 = r['comp_id_2']
                    atom_id_1 = r['atom_id_1']
                    atom_id_2 = r['atom_id_2']
                    ins_code_1 = r.get('ins_code_1', '?')
                    ins_code_2 = r.get('ins_code_2', '?')
                    target_value = r.get('target_value')
                    target_value_uncertainty = r.get('target_value_uncertainty')
                    lower_limit = r['lower_limit']
                    upper_limit = r['upper_limit']
                    lower_linear_limit = r.get('lower_linear_limit')
                    upper_linear_limit = r.get('upper_linear_limit')

                    if None in (atom_id_1, atom_id_2)\
                       or not isinstance(auth_seq_id_1, int) or not isinstance(auth_seq_id_2, int):
                        if 'HOH' not in (comp_id_1, comp_id_2):
                            self.__log.write(f"+{self.__class_name__}.__extractGenDistConstraints() ++ Error  - "
                                             f"distance restraint {rest_key} {r} is not interpretable, "
                                             f"{os.path.basename(self.__nmrDataPath)}.\n")
                        skipped = True
                        continue

                    offset = abs(auth_seq_id_1 - auth_seq_id_2)

                    if r['entity_asm_id_1'] != r['entity_asm_id_2']:
                        distance_type = 'interchain'
                    elif offset == 0:
                        distance_type = 'intraresidue'
                    elif offset == 1:
                        distance_type = 'sequential'
                    elif 1 < offset < 5:
                        distance_type = 'medium'
                    else:
                        distance_type = 'long'

                    if comp_id_1 not in bb_atoms_cache:
                        bb_atoms_cache[comp_id_1] = frozenset(self.__csStat.getBackBoneAtoms(comp_id_1, incl_nstd_bb_atom=True))
                    if comp_id_2 not in bb_atoms_cache:
                        bb_atoms_cache[comp_id_2] = frozenset(self.__csStat.getBackBoneAtoms(comp_id_2, incl_nstd_bb_atom=True))
                    bb_atoms_1 = bb_atoms_cache[comp_id_1]
                    bb_atoms_2 = bb_atoms_cache[comp_id_2]

                    if atom_id_1 in bb_atoms_1 and atom_id_2 in bb_atoms_2:
                        distance_sub_type = 'backbone-backbone'
                    elif atom_id_1 in bb_atoms_1 or atom_id_2 in bb_atoms_2:
                        distance_sub_type = 'backbone-sidechain'
                    else:
                        distance_sub_type = 'sidechain-sidechain'

                    target_value, lower_bound, upper_bound =\
                        dist_target_values(target_value, target_value_uncertainty,
                                           lower_limit, upper_limit, lower_linear_limit, upper_linear_limit)

                    if target_value is None:
                        self.__log.write(f"+{self.__class_name__}.__extractGenDistConstraints() ++ Error  - "
                                         f"distance restraint {rest_key} {r} is not interpretable, "
                                         f"{os.path.basename(self.__nmrDataPath)}.\n")
                        skipped = True
                        continue

                    atom_sels = [[{'chain_id': auth_asym_id_1,
                                   'seq_id': auth_seq_id_1,
                                   'comp_id': comp_id_1,
                                   'atom_id': atom_id_1}],
                                 [{'chain_id': auth_asym_id_2,
                                   'seq_id': auth_seq_id_2,
                                   'comp_id': comp_id_2,
                                   'atom_id': atom_id_2}]]

                    dst_func = {'lower_limit': str(lower_bound) if lower_bound is not None else None,
                                'upper_limit': str(upper_bound) if upper_bound is not None else None}

                    const_type = getDistConstraintType(atom_sels, dst_func, self.__csStat)

                    bond_flag = None
                    if const_type is not None:
                        if 'hydrogen bond' in const_type:
                            bond_flag = 'hbond'
                        elif 'disulfide bond' in const_type:
                            bond_flag = 'sbond'
                        elif 'diselenide bond' in const_type:
                            bond_flag = 'sebond'
                        elif const_type == 'metal coordination':
                            bond_flag = 'metal'

                    or_member = r.get('member_logic_code') != 'AND'

                    self.__distRestDict[rest_key].append({'atom_key_1': (auth_asym_id_1, auth_seq_id_1, comp_id_1,
                                                                         atom_id_1, ins_code_1),
                                                          'atom_key_2': (auth_asym_id_2, auth_seq_id_2, comp_id_2,
                                                                         atom_id_2, ins_code_2),
                                                          'combination_id': r.get('combination_id'),
                                                          'member_id': r.get('member_id'),
                                                          'or_member': or_member,
                                                          'distance_type': distance_type,
                                                          'distance_sub_type': distance_sub_type,
                                                          'bond_flag': bond_flag,
                                                          'lower_bound': lower_bound,
                                                          'upper_bound': upper_bound,
                                                          'target_value': target_value})

                    seq_key_1 = (auth_asym_id_1, auth_seq_id_1, comp_id_1)
                    seq_key_2 = (auth_asym_id_2, auth_seq_id_2, comp_id_2)

                    seq_keys = dict.fromkeys([seq_key_1, seq_key_2])

                    for seq_key in seq_keys:
                        if seq_key not in self.__distRestSeqDict:
                            self.__distRestSeqDict[seq_key] = []
                        self.__distRestSeqDict[seq_key].append(rest_key)

            if skipped:
                __distRestDict__ = copy.copy(self.__distRestDict)
                for k, v in __distRestDict__.items():
                    if len(v) == 0:
                        del self.__distRestDict[k]

            for k, v in self.__distRestSeqDict.items():
                self.__distRestSeqDict[k] = list(set(v))

            if len(self.__distRestDict) == 0:
                self.__distRestDict = self.__distRestSeqDict = None

            return True

        except Exception as e:  # pylint: disable=broad-exception-caught
            self.__log.write(f"Exception occurred while processing {os.path.basename(self.__cifPath)} "
                             f"and {os.path.basename(self.__nmrDataPath)}\n")
            self.__log.write(f"+{self.__class_name__}.__extractGenDistConstraints() ++ Error  - {str(e)}\n")

            self.__distRestDict = self.__distRestSeqDict = None

        return False

    def __extractTorsionAngleConstraints(self) -> bool:
        """ Extract Torsion_angle_constraint category of NMR data file.
            @author: Masashi Yokochi
            @note: Derived from wwpdb.apps.validation.src.RestraintsValidation.BMRBRestraintsAnalysis.get_restraints2,
                   written by Kumaran Baskaran
            @change: class method, use of wwpdb.utils.nmr.io.CifReader, improve readability of restraints,
                     support combinational restraints (_Torsion_angle_constraint.Combination_ID),
                     support for the case target_value is not set, but upper/lower_limit and upper/lower_linear_limit are set
                             (i.e. AMBER restraint format, decide target_value by testing anti-clockwise and clockwise mean),
                     support for the case target_value is not set, but upper/lower_limit are set
                             (i.e. CYANA restraint format, decide target_value by comparing lower_limit and upper_limit values),
                     support for the case upper/lower_linear_limit are set, but missing upper/lower_limit
                             (i.e. XPLOR-NIH/CNS exponent parameter (ed) equals 1)
        """

        if self.__has_prev_results:
            return True

        if self.__nmrDataPath is None:
            return False

        self.__dihedRestDict, self.__dihedRestSeqDict = {}, {}

        lp_category = 'Torsion_angle_constraint'
        sf_category = 'Torsion_angle_constraint_list'

        try:

            skipped = False

            for datablock_name in self.__rR.getDataBlockNameList():

                if not self.__rR.hasCategory(lp_category, datablock_name):
                    continue

                sf_tag = self.__rR.getDictList(sf_category, datablock_name)

                list_id = int(sf_tag[0]['ID'])

                data_items = [{'name': 'ID', 'type': 'int', 'alt_name': 'id'},
                              {'name': 'Auth_asym_ID_1', 'type': 'str', 'alt_name': 'auth_asym_id_1'},
                              {'name': 'Auth_seq_ID_1', 'type': 'int', 'alt_name': 'auth_seq_id_1'},
                              {'name': 'Comp_ID_1', 'type': 'str', 'alt_name': 'comp_id_1'},
                              {'name': 'Atom_ID_1', 'type': 'str', 'alt_name': 'atom_id_1'},
                              {'name': 'Auth_asym_ID_2', 'type': 'str', 'alt_name': 'auth_asym_id_2'},
                              {'name': 'Auth_seq_ID_2', 'type': 'int', 'alt_name': 'auth_seq_id_2'},
                              {'name': 'Comp_ID_2', 'type': 'str', 'alt_name': 'comp_id_2'},
                              {'name': 'Atom_ID_2', 'type': 'str', 'alt_name': 'atom_id_2'},
                              {'name': 'Auth_asym_ID_3', 'type': 'str', 'alt_name': 'auth_asym_id_3'},
                              {'name': 'Auth_seq_ID_3', 'type': 'int', 'alt_name': 'auth_seq_id_3'},
                              {'name': 'Comp_ID_3', 'type': 'str', 'alt_name': 'comp_id_3'},
                              {'name': 'Atom_ID_3', 'type': 'str', 'alt_name': 'atom_id_3'},
                              {'name': 'Auth_asym_ID_4', 'type': 'str', 'alt_name': 'auth_asym_id_4'},
                              {'name': 'Auth_seq_ID_4', 'type': 'int', 'alt_name': 'auth_seq_id_4'},
                              {'name': 'Comp_ID_4', 'type': 'str', 'alt_name': 'comp_id_4'},
                              {'name': 'Atom_ID_4', 'type': 'str', 'alt_name': 'atom_id_4'},
                              {'name': 'Angle_lower_bound_val', 'type': 'float', 'alt_name': 'lower_limit'},
                              {'name': 'Angle_upper_bound_val', 'type': 'float', 'alt_name': 'upper_limit'},
                              {'name': 'Angle_target_val', 'type': 'float', 'alt_name': 'target_value'}
                              ]

                tags = self.__rR.getAttributeList(lp_category, datablock_name)

                has_torsion_angle_name = 'Torsion_angle_name' in tags
                has_combination_id = 'Combination_ID' in tags
                has_pdb_ins_code_1 = 'PDB_ins_code_1' in tags
                has_pdb_ins_code_2 = 'PDB_ins_code_2' in tags
                has_pdb_ins_code_3 = 'PDB_ins_code_3' in tags
                has_pdb_ins_code_4 = 'PDB_ins_code_4' in tags
                has_lower_linear_limit = 'Angle_lower_linear_limit' in tags
                has_upper_linear_limit = 'Angle_upper_linear_limit' in tags
                has_target_val_err = 'Angle_target_val_err' in tags

                if has_torsion_angle_name:
                    data_items.append({'name': 'Torsion_angle_name', 'type': 'str', 'alt_name': 'angle_type', 'default': 'UNNAMED'})
                if has_combination_id:
                    data_items.append({'name': 'Combination_ID', 'type': 'int', 'alt_name': 'combination_id'})
                if has_pdb_ins_code_1:
                    data_items.append({'name': 'PDB_ins_code_1', 'type': 'str', 'alt_name': 'ins_code_1', 'default': '?'})
                if has_pdb_ins_code_2:
                    data_items.append({'name': 'PDB_ins_code_2', 'type': 'str', 'alt_name': 'ins_code_2', 'default': '?'})
                if has_pdb_ins_code_3:
                    data_items.append({'name': 'PDB_ins_code_3', 'type': 'str', 'alt_name': 'ins_code_3', 'default': '?'})
                if has_pdb_ins_code_4:
                    data_items.append({'name': 'PDB_ins_code_4', 'type': 'str', 'alt_name': 'ins_code_4', 'default': '?'})
                if has_lower_linear_limit:
                    data_items.append({'name': 'Angle_lower_linear_limit', 'type': 'float', 'alt_name': 'lower_linear_limit'})
                if has_upper_linear_limit:
                    data_items.append({'name': 'Angle_upper_linear_limit', 'type': 'float', 'alt_name': 'upper_linear_limit'})
                if has_target_val_err:
                    data_items.append({'name': 'Angle_target_val_err', 'type': 'abs-float', 'alt_name': 'target_value_uncertainty'})

                filter_items = [{'name': 'Torsion_angle_constraint_list_ID', 'type': 'int', 'value': list_id}]

                rest = self.__rR.getDictListWithFilter(lp_category,
                                                       data_items,
                                                       filter_items,
                                                       datablock_name)

                for r in rest:
                    rest_key = (list_id, r['id'])

                    if rest_key not in self.__dihedRestDict:
                        self.__dihedRestDict[rest_key] = []

                    angle_type = r.get('angle_type', 'UNNAMED')
                    auth_asym_id_1 = r['auth_asym_id_1']
                    auth_asym_id_2 = r['auth_asym_id_2']
                    auth_asym_id_3 = r['auth_asym_id_3']
                    auth_asym_id_4 = r['auth_asym_id_4']
                    auth_seq_id_1 = r['auth_seq_id_1']
                    auth_seq_id_2 = r['auth_seq_id_2']
                    auth_seq_id_3 = r['auth_seq_id_3']
                    auth_seq_id_4 = r['auth_seq_id_4']
                    comp_id_1 = r['comp_id_1']
                    comp_id_2 = r['comp_id_2']
                    comp_id_3 = r['comp_id_3']
                    comp_id_4 = r['comp_id_4']
                    atom_id_1 = r['atom_id_1']
                    atom_id_2 = r['atom_id_2']
                    atom_id_3 = r['atom_id_3']
                    atom_id_4 = r['atom_id_4']
                    ins_code_1 = r.get('ins_code_1', '?')
                    ins_code_2 = r.get('ins_code_2', '?')
                    ins_code_3 = r.get('ins_code_3', '?')
                    ins_code_4 = r.get('ins_code_4', '?')

                    if None in (atom_id_1, atom_id_2, atom_id_3, atom_id_4)\
                       or not isinstance(auth_seq_id_1, int) or not isinstance(auth_seq_id_2, int)\
                       or not isinstance(auth_seq_id_3, int) or not isinstance(auth_seq_id_4, int):
                        if angle_type not in ('PPA', 'UNNAMED'):
                            self.__log.write(f"+{self.__class_name__}.__extractTorsionAngleConstraints() ++ Error  - "
                                             f"dihedral angle restraint {rest_key} {r} is not interpretable, "
                                             f"{os.path.basename(self.__nmrDataPath)}.\n")
                        skipped = True
                        continue

                    lower_limit = r['lower_limit']
                    upper_limit = r['upper_limit']
                    target_value = r['target_value']
                    lower_linear_limit = r.get('lower_linear_limit')
                    upper_linear_limit = r.get('upper_linear_limit')
                    target_value_uncertainty = r.get('target_value_uncertainty')

                    target_value, lower_bound, upper_bound =\
                        angle_target_values(target_value, target_value_uncertainty,
                                            lower_limit, upper_limit, lower_linear_limit, upper_linear_limit)

                    if target_value is None:
                        self.__log.write(f"+{self.__class_name__}.__extractTorsionAngleConstraints() ++ Error  - "
                                         f"dihedral angle restraint {rest_key} {r} is not interpretable, "
                                         f"{os.path.basename(self.__nmrDataPath)}.\n")
                        skipped = True
                        continue

                    self.__dihedRestDict[rest_key].append({'atom_key_1': (auth_asym_id_1, auth_seq_id_1, comp_id_1,
                                                                          atom_id_1, ins_code_1),
                                                           'atom_key_2': (auth_asym_id_2, auth_seq_id_2, comp_id_2,
                                                                          atom_id_2, ins_code_2),
                                                           'atom_key_3': (auth_asym_id_3, auth_seq_id_3, comp_id_3,
                                                                          atom_id_3, ins_code_3),
                                                           'atom_key_4': (auth_asym_id_4, auth_seq_id_4, comp_id_4,
                                                                          atom_id_4, ins_code_4),
                                                           'combination_id': r.get('combination_id'),
                                                           'angle_type': angle_type,
                                                           'lower_bound': lower_bound,
                                                           'upper_bound': upper_bound,
                                                           'target_value': target_value})

                    seq_key_1 = (auth_asym_id_1, auth_seq_id_1, comp_id_1)
                    seq_key_2 = (auth_asym_id_2, auth_seq_id_2, comp_id_2)
                    seq_key_3 = (auth_asym_id_3, auth_seq_id_3, comp_id_3)
                    seq_key_4 = (auth_asym_id_4, auth_seq_id_4, comp_id_4)

                    seq_keys = dict.fromkeys([seq_key_1, seq_key_2, seq_key_3, seq_key_4])

                    for seq_key in seq_keys:
                        if seq_key not in self.__dihedRestSeqDict:
                            self.__dihedRestSeqDict[seq_key] = []
                        self.__dihedRestSeqDict[seq_key].append(rest_key)

            if skipped:
                __dihedRestDict__ = copy.copy(self.__dihedRestDict)
                for k, v in __dihedRestDict__.items():
                    if len(v) == 0:
                        del self.__dihedRestDict[k]

            for k, v in self.__dihedRestSeqDict.items():
                self.__dihedRestSeqDict[k] = list(set(v))

            if len(self.__dihedRestDict) == 0:
                self.__dihedRestDict = self.__dihedRestSeqDict = None

            return True

        except Exception as e:  # pylint: disable=broad-exception-caught
            self.__log.write(f"Exception occurred while processing {os.path.basename(self.__cifPath)} "
                             f"and {os.path.basename(self.__nmrDataPath)}\n")
            self.__log.write(f"+{self.__class_name__}.__extractTorsionAngleConstraints() ++ Error  - {str(e)}\n")

            self.__dihedRestDict = self.__dihedRestSeqDict = None

        return False

    def __extractRdcConstraints(self) -> bool:
        """ Extract RDC_constraint category of NMR data file.
        """

        if self.__has_prev_results:
            return True

        if self.__nmrDataPath is None:
            return False

        self.__rdcRestDict, self.__rdcRestSeqDict = {}, {}

        lp_category = 'RDC_constraint'
        sf_category = 'RDC_constraint_list'

        try:

            skipped = False

            for datablock_name in self.__rR.getDataBlockNameList():

                if not self.__rR.hasCategory(lp_category, datablock_name):
                    continue

                sf_tag = self.__rR.getDictList(sf_category, datablock_name)

                list_id = int(sf_tag[0]['ID'])

                try:
                    # e.g. RDC_HNC, RDC_NH, RDC_CN_i_1, RDC_CAHA, RDC_HNHA, RDC_HNHA_i_1,
                    # RDC_CAC, RDC_CAN, RDC_HH, RDC_CC, RDC_other
                    rdc_type = sf_tag[0]['Details']
                    if rdc_type in (None, '', '.', '?', 'null'):
                        rdc_type = 'UNNAMED'
                except KeyError:
                    rdc_type = 'UNNAMED'

                data_items = [{'name': 'ID', 'type': 'int', 'alt_name': 'id'},
                              {'name': 'Auth_asym_ID_1', 'type': 'str', 'alt_name': 'auth_asym_id_1'},
                              {'name': 'Auth_seq_ID_1', 'type': 'int', 'alt_name': 'auth_seq_id_1'},
                              {'name': 'Comp_ID_1', 'type': 'str', 'alt_name': 'comp_id_1'},
                              {'name': 'Atom_ID_1', 'type': 'str', 'alt_name': 'atom_id_1'},
                              {'name': 'Auth_asym_ID_2', 'type': 'str', 'alt_name': 'auth_asym_id_2'},
                              {'name': 'Auth_seq_ID_2', 'type': 'int', 'alt_name': 'auth_seq_id_2'},
                              {'name': 'Comp_ID_2', 'type': 'str', 'alt_name': 'comp_id_2'},
                              {'name': 'Atom_ID_2', 'type': 'str', 'alt_name': 'atom_id_2'},
                              {'name': 'Target_value', 'type': 'float', 'alt_name': 'target_value'},
                              ]

                tags = self.__rR.getAttributeList(lp_category, datablock_name)

                has_combination_id = 'Combination_ID' in tags
                has_pdb_ins_code_1 = 'PDB_ins_code_1' in tags
                has_pdb_ins_code_2 = 'PDB_ins_code_2' in tags
                has_val = 'RDC_val' in tags
                has_val_err = 'RDC_val_err' in tags
                has_lower_bound = 'RDC_lower_bound' in tags
                has_upper_bound = 'RDC_upper_bound' in tags
                has_lower_linear_limit = 'RDC_lower_linear_limit' in tags
                has_upper_linear_limit = 'RDC_upper_linear_limit' in tags
                has_target_val_uncertainty = 'Target_value_uncertainty' in tags
                has_weight = 'Weight' in tags
                has_scale_factor = 'RDC_val_scale_factor' in tags

                if has_combination_id:
                    data_items.append({'name': 'Combination_ID', 'type': 'int', 'alt_name': 'combination_id'})
                if has_pdb_ins_code_1:
                    data_items.append({'name': 'PDB_ins_code_1', 'type': 'str', 'alt_name': 'ins_code_1', 'default': '?'})
                if has_pdb_ins_code_2:
                    data_items.append({'name': 'PDB_ins_code_2', 'type': 'str', 'alt_name': 'ins_code_2', 'default': '?'})
                if has_val:
                    data_items.append({'name': 'RDC_val', 'type': 'float', 'alt_name': 'value'})
                if has_val_err:
                    data_items.append({'name': 'RDC_val_err', 'type': 'abs-float', 'alt_name': 'value_uncertainty'})
                if has_lower_bound:
                    data_items.append({'name': 'RDC_lower_bound', 'type': 'float', 'alt_name': 'lower_limit'})
                if has_upper_bound:
                    data_items.append({'name': 'RDC_upper_bound', 'type': 'float', 'alt_name': 'upper_limit'})
                if has_lower_linear_limit:
                    data_items.append({'name': 'RDC_lower_linear_limit', 'type': 'float', 'alt_name': 'lower_linear_limit'})
                if has_upper_linear_limit:
                    data_items.append({'name': 'RDC_upper_linear_limit', 'type': 'float', 'alt_name': 'upper_linear_limit'})
                if has_target_val_uncertainty:
                    data_items.append({'name': 'Target_value_uncertainty', 'type': 'abs-float',
                                       'alt_name': 'target_value_uncertainty'})
                if has_weight:
                    data_items.append({'name': 'Weight', 'type': 'float'})
                if has_scale_factor:
                    data_items.append({'name': 'RDC_val_scale_factor', 'type': 'float', 'alt_name': 'scale_factor'})

                filter_items = [{'name': 'RDC_constraint_list_ID', 'type': 'int', 'value': list_id}]

                rest = self.__rR.getDictListWithFilter(lp_category,
                                                       data_items,
                                                       filter_items,
                                                       datablock_name)

                for r in rest:
                    rest_key = (list_id, r['id'])

                    if rest_key not in self.__rdcRestDict:
                        self.__rdcRestDict[rest_key] = []

                    auth_asym_id_1 = r['auth_asym_id_1']
                    auth_asym_id_2 = r['auth_asym_id_2']
                    auth_seq_id_1 = r['auth_seq_id_1']
                    auth_seq_id_2 = r['auth_seq_id_2']
                    comp_id_1 = r['comp_id_1']
                    comp_id_2 = r['comp_id_2']
                    atom_id_1 = r['atom_id_1']
                    atom_id_2 = r['atom_id_2']
                    ins_code_1 = r.get('ins_code_1', '?')
                    ins_code_2 = r.get('ins_code_2', '?')

                    if None in (atom_id_1, atom_id_2)\
                       or not isinstance(auth_seq_id_1, int) or not isinstance(auth_seq_id_2, int):
                        self.__log.write(f"+{self.__class_name__}.__extractRdcConstraints() "
                                         f"++ Error  - RDC restraint {rest_key} {r} is not interpretable, "
                                         f"{os.path.basename(self.__nmrDataPath)}.\n")
                        skipped = True
                        continue

                    _rdc_type = rdc_type
                    if rdc_type == 'UNNAMED':
                        try:
                            atom1 = {'chain_id': auth_asym_id_1,
                                     'seq_id': int(auth_seq_id_1),
                                     'comp_id': comp_id_1,
                                     'atom_id': atom_id_1}
                            atom2 = {'chain_id': auth_asym_id_2,
                                     'seq_id': int(auth_seq_id_2),
                                     'comp_id': comp_id_2,
                                     'atom_id': atom_id_2}
                            _rdc_type = getRdcCode([atom1, atom2])
                        except (ValueError, TypeError):
                            pass

                    target_value = r['target_value']
                    target_value_uncertainty = r.get('target_value_uncertainty')
                    value = r.get('value')
                    value_uncertainty = r.get('value_uncertainty')
                    lower_limit = r.get('lower_limit')
                    upper_limit = r.get('upper_limit')
                    lower_linear_limit = r.get('lower_linear_limit')
                    upper_linear_limit = r.get('upper_linear_limit')
                    weight = r.get('weight')
                    if weight is None or weight < 0.0:
                        weight = 1.0
                    scale_factor = r.get('scale_factor')
                    if scale_factor is None or scale_factor < 0.0:
                        scale_factor = 1.0

                    target_value, lower_bound, upper_bound =\
                        rdc_target_values(target_value, target_value_uncertainty, value, value_uncertainty,
                                          lower_limit, upper_limit, lower_linear_limit, upper_linear_limit)

                    if target_value is None:
                        self.__log.write(f"+{self.__class_name__}.__extractRdcConstraints() "
                                         f"++ Error  - RDC restraint {rest_key} {r} is not interpretable, "
                                         f"{os.path.basename(self.__nmrDataPath)}.\n")
                        skipped = True
                        continue

                    self.__rdcRestDict[rest_key].append({'atom_key_1': (auth_asym_id_1, auth_seq_id_1, comp_id_1,
                                                                        atom_id_1, ins_code_1),
                                                         'atom_key_2': (auth_asym_id_2, auth_seq_id_2, comp_id_2,
                                                                        atom_id_2, ins_code_2),
                                                         'combination_id': r.get('combination_id'),
                                                         'rdc_type': _rdc_type,
                                                         'lower_bound': lower_bound,
                                                         'upper_bound': upper_bound,
                                                         'target_value': target_value,
                                                         'weight': weight,
                                                         'scale_factor': scale_factor})

                    seq_key_1 = (auth_asym_id_1, auth_seq_id_1, comp_id_1)
                    seq_key_2 = (auth_asym_id_2, auth_seq_id_2, comp_id_2)

                    seq_keys = dict.fromkeys([seq_key_1, seq_key_2])

                    for seq_key in seq_keys:
                        if seq_key not in self.__rdcRestSeqDict:
                            self.__rdcRestSeqDict[seq_key] = []
                        self.__rdcRestSeqDict[seq_key].append(rest_key)

            if skipped:
                __rdcRestDict__ = copy.copy(self.__rdcRestDict)
                for k, v in __rdcRestDict__.items():
                    if len(v) == 0:
                        del self.__rdcRestDict[k]

            for k, v in self.__rdcRestSeqDict.items():
                self.__rdcRestSeqDict[k] = list(set(v))

            if len(self.__rdcRestDict) != 0:
                # settle molecular alignment tensor for each rdc list
                self.__settleMolecularAlignmentTensors()

            else:
                self.__rdcRestDict = self.__rdcRestSeqDict = None

            return True

        except Exception as e:  # pylint: disable=broad-exception-caught
            self.__log.write(f"Exception occurred while processing {os.path.basename(self.__cifPath)} "
                             f"and {os.path.basename(self.__nmrDataPath)}\n")
            self.__log.write(f"+{self.__class_name__}.__extractRdcConstraints() ++ Error  - {str(e)}\n")

            self.__rdcRestDict = self.__rdcRestSeqDict = None

        return False

    def __settleMolecularAlignmentTensors(self) -> bool:
        """ Settle molecular alignment tensors for each RDC list based on coordinates and experimental RDC values.
            @author: Masashi Yokochi
            Reference:
              Order Matrix Analysis of Residual Dipolar Couplings Using Singular Value Decomposition
              J. A. Losonczi, M. Andrec, M. W. F. Fischer, and J. H. Prestegard.
              Journal of Magnetic Resonance, 138, 334-342 (1999)
              DOI: 10.1006/jmre.1999.1754
        """

        self.__rdcSaupeOrderMatrix, self.__rdcCalcDict, self.__rdcSyntCalcDict = {}, {}, {}

        try:

            cycles = max(RDC_EFF_MC_CYCLES // len(self.__eff_model_ids), 1)
            rng = numpy.random.default_rng(999)

            list_ids = set()
            for rest_key in self.__rdcRestDict:
                list_ids.add(rest_key[0])

            for list_id in sorted(list(list_ids)):

                dmax_dict = {}

                self.__rdcSaupeOrderMatrix[list_id] = {}

                for model_id in self.__coordinates:

                    target_rest_keys = []

                    A = numpy.empty((0, 5), dtype=float)
                    b_exp_list, dmax_list = [], []

                    for rest_key, restraints in self.__rdcRestDict.items():

                        if rest_key[0] != list_id:
                            continue

                        for r in restraints:

                            if 0.0 in (r['weight'], r['scale_factor']):
                                continue

                            atom_key_1 = r['atom_key_1']
                            atom_key_2 = r['atom_key_2']

                            atom_present = True

                            try:
                                pos_1 = self.__coordinates[model_id][atom_key_1]
                            except KeyError:
                                atom_present = False

                            try:
                                pos_2 = self.__coordinates[model_id][atom_key_2]
                            except KeyError:
                                atom_present = False

                            if not atom_present:
                                continue

                            atom_type_1 = atom_key_1[3][0]
                            atom_type_2 = atom_key_2[3][0]

                            if atom_type_1 not in GYROMAGNETIC_RATIOS\
                               or atom_type_2 not in GYROMAGNETIC_RATIOS:
                                continue

                            target_rest_keys.append(rest_key)

                            vec_key = (atom_type_1, atom_type_2)

                            if vec_key not in dmax_dict:
                                dmax_dict[vec_key] = rdc_dmax(atom_type_1, atom_type_2, distance(pos_1, pos_2), hz_unit=True)

                            _dmax = dmax_dict[vec_key]

                            vector = to_unit_vector(pos_2 - pos_1)
                            cos_x, cos_y, cos_z = vector[0], vector[1], vector[2]
                            cos_x_2 = cos_x * cos_x
                            cos_y_2 = cos_y * cos_y
                            cos_z_2 = cos_z * cos_z

                            A = numpy.append(A, numpy.array([[cos_y_2 - cos_x_2, cos_z_2 - cos_x_2,
                                                              2.0 * cos_x * cos_y,
                                                              2.0 * cos_x * cos_z,
                                                              2.0 * cos_y * cos_z]]), axis=0)

                            # input observed RDC values calibrated with a given scale factor
                            # so, calculated RDC values will be scaled by default,
                            # this effect should be taken into consideration when comparing (raw) observed RDCs and calculated RDCs
                            b_exp_list.append(r['scale_factor'] * r['target_value'] / _dmax)

                            dmax_list.append(_dmax)

                    b_size = len(b_exp_list)

                    if b_size < 5:
                        continue

                    b_exp = numpy.array(b_exp_list, dtype=float)
                    dmax = numpy.array(dmax_list, dtype=float)

                    U, S, Vh = numpy.linalg.svd(A, full_matrices=False)

                    Si = numpy.diag(numpy.array([1.0 / s for s in list(S)], dtype=float))

                    Ai = Vh.T @ Si @ U.T

                    x = Ai @ b_exp

                    Syy, Szz, Sxy, Sxz, Syz = x[0], x[1], x[2], x[3], x[4]

                    if -0.5 <= Syy <= 1.0 and -0.5 <= Szz <= 1.0 and abs(Sxy) <= 0.75 and abs(Sxz) <= 0.75 and abs(Syz) <= 0.75:
                        Sxx = -(Syy + Szz)
                        Sorder = sorted([abs(Sxx), abs(Syy), abs(Szz)], reverse=True)
                        Szz_, Syy_, Sxx_ = Sorder[0], Sorder[1], Sorder[2]

                        if Szz_ == abs(Szz):
                            Szz_ = Szz
                            if Syy_ == abs(Syy):  # zz > yy > xx (as is)
                                Sxx_, Syy_ = Sxx, Syy
                                Sxy_, Sxz_, Syz_ = Sxy, Sxz, Syz
                            else:  # zz > xx > yy : x <-> y, z -> -z
                                Sxx_, Syy_ = Syy, Sxx
                                Sxy_, Sxz_, Syz_ = Sxy, -Syz, -Sxz
                        elif Szz_ == abs(Syy):
                            Szz_ = Syy
                            if Syy_ == abs(Szz):  # yy > zz > xx : y <-> z, x -> -x
                                Sxx_, Syy_ = Sxx, Szz
                                Sxy_, Sxz_, Syz_ = -Sxz, -Sxy, Syz
                            else:  # yy > xx > zz : y->z, x->y, z->x (rotation)
                                Sxx_, Syy_ = Szz, Sxx
                                Sxy_, Sxz_, Syz_ = Syz, Sxy, Sxz
                        else:
                            Szz_ = Sxx
                            if Syy_ == abs(Szz):  # xx > zz > yy : x->z, z->y, y->x (rotation)
                                Sxx_, Syy_ = Syy, Szz
                                Sxy_, Sxz_, Syz_ = Sxz, Syz, Sxy
                            else:  # xx > yy > zz : x <-> z, y -> -y
                                Sxx_, Syy_ = Szz, Syy
                                Sxy_, Sxz_, Syz_ = -Syz, Sxz, -Sxy

                        eta = (Syy_ - Sxx_) / Szz_

                        self.__rdcSaupeOrderMatrix[list_id][model_id] = {'Sxx': f'{Sxx_:.4e}', 'Syy': f'{Syy_:.4e}',
                                                                         'Szz': f'{Szz_:.4e}',
                                                                         'Sxy': f'{Sxy_:.4e}', 'Sxz': f'{Sxz_:.4e}',
                                                                         'Syz': f'{Syz_:.4e}',
                                                                         'eta': f'{eta:.4e}',
                                                                         'Dmax': f'{numpy.mean(dmax):.4e}'}

                        assert abs(Szz_) >= abs(Syy_) >= abs(Sxx_)
                        assert 0 <= abs(eta) <= 1.0

                        b_calc = A @ x

                        for rest_key, val in zip(target_rest_keys, list(b_calc * dmax)):
                            if rest_key not in self.__rdcCalcDict:
                                self.__rdcCalcDict[rest_key] = {}
                            self.__rdcCalcDict[rest_key][model_id] = val

                        b_std = numpy.std(b_calc - b_exp)

                        for _ in range(cycles):
                            b_noise = rng.normal(loc=0.0, scale=b_std, size=b_size)

                            b_syn = b_exp + b_noise

                            x_syn = Ai @ b_syn

                            _Syy, _Szz, _Sxy, _Sxz, _Syz = x_syn[0], x_syn[1], x_syn[2], x_syn[3], x_syn[4]

                            if -0.5 <= _Syy <= 1.0 and -0.5 <= _Szz <= 1.0\
                               and abs(_Sxy) <= 0.75 and abs(_Sxz) <= 0.75 and abs(_Syz) <= 0.75:

                                b_syn_calc = A @ x_syn

                                for rest_key, val in zip(target_rest_keys, list(b_syn_calc * dmax)):
                                    if rest_key not in self.__rdcSyntCalcDict:
                                        self.__rdcSyntCalcDict[rest_key] = {}
                                    if model_id not in self.__rdcSyntCalcDict[rest_key]:
                                        self.__rdcSyntCalcDict[rest_key][model_id] = []
                                    self.__rdcSyntCalcDict[rest_key][model_id].append(val)

            return True

        except Exception as e:  # pylint: disable=broad-exception-caught
            self.__log.write(f"Exception occurred while processing {os.path.basename(self.__cifPath)} "
                             f"and {os.path.basename(self.__nmrDataPath)}\n")
            self.__log.write(f"+{self.__class_name__}.__settleMolecularAlignmentTensors() ++ Error  - {str(e)}\n")

        return False

    def __validateChemicalShifts(self) -> bool:
        """ Validate assigned chemical shifts.
            @author: Masashi Yokochi
            @note: Derived from wwpdb.apps.validation.src.ChemicalShiftsValidation.BMRBChemicalShiftsAnalysis.get_chemical_shifts
                   written by Aleksandras Gutmanas, Kumaran Baskaran
            @change: class method
        """

        if self.__chemShiftDict is None or len(self.__chemShiftDict) == 0 or self.__has_prev_results:
            return True

        self.__chemShiftUniqDict, self.__chemShiftOutlier, self.__chemShiftDuplicated, \
            self.__chemShiftUnmapped, self.__chemShiftUnmodeled = {}, {}, {}, {}, {}

        has_coord = self.__coordinates is not None

        try:

            def check_entity_instance(list_id, _cs):
                auth_chain_id = _cs['auth_chain_id']

                if auth_chain_id not in self.__entityInstance and auth_chain_id not in self.__entityUninstance:
                    if list_id not in self.__chemShiftUnmapped:
                        self.__chemShiftUnmapped[list_id] = []
                    self.__chemShiftUnmapped[list_id].append(_cs)
                    return

                auth_seq_id = _cs['auth_seq_id']
                comp_id = _cs['comp_id']
                atom_id = _cs['atom_id']

                if 'ins_code' not in _cs or _cs['ins_code'] in EMPTY_VALUE:
                    seq_key = (str(auth_seq_id), comp_id)
                else:
                    seq_key = (str(auth_seq_id) + _cs['ins_code'], comp_id)

                # DAOTHER-10987
                if auth_chain_id in self.__entityUninstance and seq_key in self.__entityUninstance[auth_chain_id]\
                   and atom_id in self.__entityUninstance[auth_chain_id][seq_key]['atoms']:
                    if list_id not in self.__chemShiftUnmodeled:
                        self.__chemShiftUnmodeled[list_id] = []
                    self.__chemShiftUnmodeled[list_id].append(_cs)
                    return

                if seq_key not in self.__entityInstance[auth_chain_id]:
                    # DAOTHER-9875, 10987 in case missing of caC['coord_unobs_res']
                    if self.__inputReport is not None:
                        unmapped_seq_keys = self.__inputReport.getUnmappedSeqKeyOfModelChainId(auth_chain_id, False)
                        if unmapped_seq_keys is not None:
                            if seq_key in unmapped_seq_keys:
                                if list_id not in self.__chemShiftUnmodeled:
                                    self.__chemShiftUnmodeled[list_id] = []
                                self.__chemShiftUnmodeled[list_id].append(_cs)
                                return

                    if list_id not in self.__chemShiftUnmapped:
                        self.__chemShiftUnmapped[list_id] = []
                    self.__chemShiftUnmapped[list_id].append(_cs)
                    return

                atoms = self.__entityInstance[auth_chain_id][seq_key]['atoms']

                if atom_id in atoms:
                    return

                if self.__entityInstance[auth_chain_id][seq_key]['seq_id'] == 1\
                   and self.__ccU.peptideLike(comp_id)\
                   and atom_id in AMINO_PROTON_CODE\
                   and any(True for a_proton in AMINO_PROTON_CODE if a_proton in atoms):
                    return

                self.__chemShiftUnmapped[list_id].append(_cs)

            for list_id, cs_data in self.__chemShiftDict.items():

                if list_id in self.__chemShiftUniqDict:
                    continue

                has_ins_code = any(True for r in cs_data if 'ins_code' in r and r['ins_code'] not in EMPTY_VALUE)

                self.__chemShiftUniqDict[list_id] = {}
                self.__chemShiftOutlier[list_id], self.__chemShiftDuplicated[list_id], \
                    self.__chemShiftUnmapped[list_id], self.__chemShiftUnmodeled[list_id] = [], [], [], []

                for cs in cs_data:
                    cs_auth_chain_id = cs['auth_chain_id']
                    cs_auth_seq_id = cs['auth_seq_id']
                    cs_comp_id = cs['comp_id']
                    cs_atom_id = cs['atom_id']
                    cs_value = cs['value']
                    cs_error = cs['error']
                    ambig_code = cs['ambig_code']

                    if not has_ins_code or cs['ins_code'] in EMPTY_VALUE:
                        cs_key = (cs_auth_chain_id, str(cs_auth_seq_id), cs_comp_id, cs_atom_id)
                    else:
                        cs_key = (cs_auth_chain_id, str(cs_auth_seq_id) + cs['ins_code'], cs_comp_id, cs_atom_id)

                    if cs_key not in self.__chemShiftUniqDict[list_id]:
                        self.__chemShiftUniqDict[list_id][cs_key] = {'value': cs_value,
                                                                     'error': cs_error,
                                                                     'ambig_code': ambig_code}

                        check_entity_instance(list_id, cs)

                        cs_stat = next((cs_stat for cs_stat in self.__csStat.get(cs['comp_id'])
                                        if cs_stat['atom_id'] == cs['atom_id'] and cs_stat['count'] > 0), None)

                        if cs_stat is None:
                            continue

                        avg_value = cs_stat['avg']
                        std_value = cs_stat['std']

                        if None in (avg_value, std_value):
                            continue

                        self.__chemShiftUniqDict[list_id][cs_key]['primary'] = True

                        z_score = round((cs['value'] - avg_value) / std_value, 2)

                        if abs(z_score) > 5.0:
                            details = None

                            if has_coord:
                                na = self.__getNearestAromaticRing(cs_auth_chain_id, cs_auth_seq_id, cs_atom_id)
                                pa = self.__getNearestParaFerroMagneticAtom(cs_auth_chain_id, cs_auth_seq_id, cs_atom_id)

                                if na is None and pa is None:
                                    pass
                                elif pa is None:
                                    if (na['ring_angle'] - MAGIC_ANGLE) * z_score < 0.0 or na['ring_distance'] > VICINITY_AROMATIC:
                                        pass
                                    else:
                                        details = "The nearest aromatic ring "\
                                            f"({na['auth_chain_id']}:{na['auth_seq_id']}:{na['comp_id']}:{na['ring_atoms']}) "\
                                            f"is located at a distance of {na['ring_distance']}Å, "\
                                            f"and has an elevation angle of {na['ring_angle']}° with the ring plane. "\
                                            "It may explain a primary factor behind the outlier."
                                else:
                                    if pa['distance'] > VICINITY_PARAMAGNETIC:
                                        pass
                                    else:
                                        details = "The nearest paramagnetic/ferromagnetic atom "\
                                            f"({pa['auth_chain_id']}:{pa['auth_seq_id']}:{pa['comp_id']}:{pa['atom_id']}) "\
                                            f"is located at a distance of {pa['distance']}Å. "\
                                            "It may explain a primary factor behind the outlier."

                            out_cs_values = list(cs_key)
                            out_cs_values.extend([cs_value, ambig_code, z_score,
                                                  round(avg_value - 5.0 * std_value, 2),
                                                  round(avg_value + 5.0 * std_value, 2),
                                                  details])

                            self.__chemShiftOutlier[list_id].append(out_cs_values)

                    else:
                        dup_cs_values = list(cs_key)
                        dup_cs_values.extend([cs_value, cs_error, ambig_code])
                        self.__chemShiftDuplicated[list_id].append(dup_cs_values)

            return True

        except Exception as e:  # pylint: disable=broad-exception-caught
            self.__log.write(f"Exception occurred while processing {os.path.basename(self.__cifPath)} "
                             f"and {os.path.basename(self.__nmrDataPath)}\n")
            self.__log.write(f"+{self.__class_name__}.__validateChemicalShifts() ++ Error  - {str(e)}\n")

            self.__chemShiftUniqDict = self.__chemShiftOutlier = \
                self.__chemShiftDuplicated = self.__chemShiftUnmapped = self.__chemShiftUnmodeled = None

        return False

    def __getNearestAromaticRing(self, auth_chain_id: str, auth_seq_id: int, atom_id: str
                                 ) -> Optional[dict]:
        """ Return the nearest aromatic ring around a given atom.
            @return: the nearest aromatic ring
        """

        try:

            _origin = self.__cR.getDictListWithFilter('atom_site',
                                                      [{'name': 'Cartn_x', 'type': 'float', 'alt_name': 'x'},
                                                       {'name': 'Cartn_y', 'type': 'float', 'alt_name': 'y'},
                                                       {'name': 'Cartn_z', 'type': 'float', 'alt_name': 'z'}
                                                       ],
                                                      [{'name': 'auth_asym_id', 'type': 'str', 'value': auth_chain_id},
                                                       {'name': 'auth_seq_id', 'type': 'int', 'value': auth_seq_id},
                                                       {'name': 'label_atom_id', 'type': 'str', 'value': atom_id},
                                                       {'name': 'pdbx_PDB_model_num', 'type': 'int',
                                                        'value': self.__representative_model_id},
                                                       {'name': 'label_alt_id', 'type': 'enum',
                                                        'enum': (self.__representative_alt_id,)}
                                                       ])

        except Exception as e:  # pylint: disable=broad-exception-caught

            if self.__verbose:
                self.__log.write(f"+{self.__class_name__}.__getNearestAromaticRing() "
                                 f"++ Error  - {str(e)}\n")

            return None

        if len(_origin) != 1:
            return None

        o = to_np_array(_origin[0])

        try:

            _neighbor = self.__cR.getDictListWithFilter('atom_site',
                                                        [{'name': 'auth_asym_id', 'type': 'str', 'alt_name': 'auth_chain_id'},
                                                         {'name': 'auth_seq_id', 'type': 'int', 'alt_name': 'auth_seq_id'},
                                                         {'name': 'label_comp_id', 'type': 'starts-with-alnum',
                                                          'alt_name': 'comp_id'},
                                                         {'name': 'label_atom_id', 'type': 'starts-with-alnum',
                                                          'alt_name': 'atom_id'},
                                                         {'name': 'Cartn_x', 'type': 'float', 'alt_name': 'x'},
                                                         {'name': 'Cartn_y', 'type': 'float', 'alt_name': 'y'},
                                                         {'name': 'Cartn_z', 'type': 'float', 'alt_name': 'z'},
                                                         {'name': 'type_symbol', 'type': 'str'}
                                                         ],
                                                        [{'name': 'Cartn_x', 'type': 'range-float',
                                                          'range': {'min_exclusive': (o[0] - CUTOFF_AROMATIC),
                                                                    'max_exclusive': (o[0] + CUTOFF_AROMATIC)}},
                                                         {'name': 'Cartn_y', 'type': 'range-float',
                                                          'range': {'min_exclusive': (o[1] - CUTOFF_AROMATIC),
                                                                    'max_exclusive': (o[1] + CUTOFF_AROMATIC)}},
                                                         {'name': 'Cartn_z', 'type': 'range-float',
                                                          'range': {'min_exclusive': (o[2] - CUTOFF_AROMATIC),
                                                                    'max_exclusive': (o[2] + CUTOFF_AROMATIC)}},
                                                         {'name': 'pdbx_PDB_model_num', 'type': 'int',
                                                          'value': self.__representative_model_id},
                                                         {'name': 'label_alt_id', 'type': 'enum',
                                                          'enum': (self.__representative_alt_id,)}
                                                         ])

        except Exception as e:  # pylint: disable=broad-exception-caught

            if self.__verbose:
                self.__log.write(f"+{self.__class_name__}.__getNearestAromaticRing() "
                                 f"++ Error  - {str(e)}\n")

            return None

        if len(_neighbor) == 0:
            return None

        neighbor = [n for n in _neighbor
                    if n['auth_seq_id'] != auth_seq_id
                    and n['type_symbol'] not in PROTON_BEGIN_CODE
                    and distance(to_np_array(n), o) < CUTOFF_AROMATIC
                    and n['atom_id'] in self.__csStat.getAromaticAtoms(n['comp_id'])]

        if len(neighbor) == 0:
            return None

        atom_list = []

        for n in neighbor:
            atom_list.append({'auth_chain_id': n['auth_chain_id'],
                              'auth_seq_id': n['auth_seq_id'],
                              'comp_id': n['comp_id'],
                              'atom_id': n['atom_id'],
                              'distance': distance(to_np_array(n), o)})

        if len(atom_list) == 0:
            return None

        na = sorted(atom_list, key=itemgetter('distance'))[0]

        na_atom_id = na['atom_id']

        if not self.__ccU.updateChemCompDict(na['comp_id']):
            return None

        # matches with comp_id in CCD

        half_ring_traces = []

        for b1 in self.__ccU.lastBondDictList:

            if b1['aromatic_flag'] != 'Y':
                continue

            if b1['atom_id_1'] == na_atom_id and b1['atom_id_2'][0] not in PROTON_BEGIN_CODE:
                na_ = b1['atom_id_2']

            elif b1['atom_id_2'] == na_atom_id and b1['atom_id_1'][0] not in PROTON_BEGIN_CODE:
                na_ = b1['atom_id_1']

            else:
                continue

            for b2 in self.__ccU.lastBondDictList:

                if b2['aromatic_flag'] != 'Y':
                    continue

                if b2['atom_id_1'] == na_ and b2['atom_id_2'][0] not in PROTON_BEGIN_CODE\
                        and b2['atom_id_2'] != na_atom_id:
                    na__ = b2['atom_id_2']

                elif b2['atom_id_2'] == na_ and b2['atom_id_1'][0] not in PROTON_BEGIN_CODE\
                        and b2['atom_id_1'] != na_atom_id:
                    na__ = b2['atom_id_1']

                else:
                    continue

                for b3 in self.__ccU.lastBondDictList:

                    if b3['aromatic_flag'] != 'Y':
                        continue

                    if b3['atom_id_1'] == na__ and b3['atom_id_2'][0] not in PROTON_BEGIN_CODE\
                            and b3['atom_id_2'] != na_:
                        na___ = b3['atom_id_2']

                    elif b3['atom_id_2'] == na__ and b3['atom_id_1'][0] not in PROTON_BEGIN_CODE\
                            and b3['atom_id_1'] != na_:
                        na___ = b3['atom_id_1']

                    else:
                        continue

                    half_ring_traces.append(f"{na_atom_id}:{na_}:{na__}:{na___}")

        len_half_ring_traces = len(half_ring_traces)

        if len_half_ring_traces < 2:
            return None

        ring_traces = []

        for i in range(len_half_ring_traces - 1):

            half_ring_trace_1 = half_ring_traces[i].split(':')

            for j in range(i + 1, len_half_ring_traces):

                half_ring_trace_2 = half_ring_traces[j].split(':')

                # hexagonal ring
                if half_ring_trace_1[3] == half_ring_trace_2[3]:
                    ring_traces.append(f'{half_ring_traces[i]}:{half_ring_trace_2[2]}:{half_ring_trace_2[1]}')

                # pentagonal ring
                elif half_ring_trace_1[3] == half_ring_trace_2[2] and half_ring_trace_1[2] == half_ring_trace_2[3]:
                    ring_traces.append(f'{half_ring_traces[i]}:{half_ring_trace_2[1]}')

        if len(ring_traces) == 0:
            return None

        ring_atoms = None
        ring_trace_score = 0

        for ring_trace in ring_traces:

            _ring_atoms = ring_trace.split(':')

            score = 0

            for a in atom_list:

                if a['auth_chain_id'] != na['auth_chain_id'] or a['auth_seq_id'] != na['auth_seq_id']\
                   or a['comp_id'] != na['comp_id']:
                    continue

                if a['atom_id'] in _ring_atoms:
                    score += 1

            if score > ring_trace_score:
                ring_atoms = _ring_atoms
                ring_trace_score = score

        try:

            _na = self.__cR.getDictListWithFilter('atom_site',
                                                  [{'name': 'label_atom_id', 'type': 'starts-with-alnum',
                                                    'alt_name': 'atom_id'},
                                                   {'name': 'Cartn_x', 'type': 'float', 'alt_name': 'x'},
                                                   {'name': 'Cartn_y', 'type': 'float', 'alt_name': 'y'},
                                                   {'name': 'Cartn_z', 'type': 'float', 'alt_name': 'z'},
                                                   {'name': 'pdbx_PDB_model_num', 'type': 'int', 'alt_name': 'model_id'}
                                                   ],
                                                  [{'name': 'auth_asym_id', 'type': 'str', 'value': na['auth_chain_id']},
                                                   {'name': 'auth_seq_id', 'type': 'int', 'value': na['auth_seq_id']},
                                                   {'name': 'label_comp_id', 'type': 'str', 'value': na['comp_id']},
                                                   {'name': 'label_atom_id', 'type': 'enum', 'enum': ring_atoms},
                                                   {'name': 'label_alt_id', 'type': 'enum',
                                                    'enum': (self.__representative_alt_id,)}
                                                   ])

        except Exception as e:  # pylint: disable=broad-exception-caught

            if self.__verbose:
                self.__log.write(f"+{self.__class_name__}.__getNearestAromaticRing() "
                                 f"++ Error  - {str(e)}\n")

            return None

        if len(_na) == 0:
            return None

        model_ids = set(a['model_id'] for a in _na)

        len_model_ids = 0

        dist = ring_dist = ring_angle = 0.0

        for model_id in model_ids:

            rc = numpy.array([0.0] * 3, dtype=float)

            total = 0

            for a in _na:

                if a['model_id'] == model_id:

                    _a = to_np_array(a)

                    if a['atom_id'] == na_atom_id:
                        dist += distance(_a, o)

                    rc = numpy.add(rc, _a)

                    total += 1

            if total == len(ring_atoms):

                rc = rc / total

                ring_dist += distance(rc, o)

                na_ = next(to_np_array(na_) for na_ in _na if na_['atom_id'] == ring_atoms[0])
                na__ = next(to_np_array(na__) for na__ in _na if na__['atom_id'] == ring_atoms[1])
                na___ = next(to_np_array(na___) for na___ in _na if na___['atom_id'] == ring_atoms[-1])

                ring_vector = numpy.cross(na__ - na_, na___ - na_)

                ring_angle += math.acos(abs(numpy.dot(to_unit_vector(o - rc), to_unit_vector(ring_vector))))

                len_model_ids += 1

        if len_model_ids == 0:  # DAOTHER-8840
            return None

        na['ring_atoms'] = ring_atoms
        na['distance'] = round(dist / len_model_ids, 1)
        na['ring_distance'] = round(ring_dist / len_model_ids, 1)
        na['ring_angle'] = round(numpy.degrees(ring_angle / len_model_ids), 1)

        return na

    def __getNearestParaFerroMagneticAtom(self, auth_chain_id: str, auth_seq_id: int, atom_id: str
                                          ) -> Optional[dict]:
        """ Return the nearest paramagnetic/ferromagnetic atom around a given atom.
            @return: the nearest paramagnetic/ferromagnetic atom
        """

        if self.__is_diamagnetic:
            return None

        try:

            _origin = self.__cR.getDictListWithFilter('atom_site',
                                                      [{'name': 'Cartn_x', 'type': 'float', 'alt_name': 'x'},
                                                       {'name': 'Cartn_y', 'type': 'float', 'alt_name': 'y'},
                                                       {'name': 'Cartn_z', 'type': 'float', 'alt_name': 'z'}
                                                       ],
                                                      [{'name': 'auth_asym_id', 'type': 'str', 'value': auth_chain_id},
                                                       {'name': 'auth_seq_id', 'type': 'int', 'value': auth_seq_id},
                                                       {'name': 'label_atom_id', 'type': 'str', 'value': atom_id},
                                                       {'name': 'pdbx_PDB_model_num', 'type': 'int',
                                                        'value': self.__representative_model_id},
                                                       {'name': 'label_alt_id', 'type': 'enum',
                                                        'enum': (self.__representative_alt_id,)}
                                                       ])

        except Exception as e:  # pylint: disable=broad-exception-caught

            if self.__verbose:
                self.__log.write(f"+{self.__class_name__}.__getNearestParaFerroMagneticAtom() "
                                 f"++ Error  - {str(e)}\n")

            return None

        if len(_origin) != 1:
            return None

        o = to_np_array(_origin[0])

        try:

            _neighbor = self.__cR.getDictListWithFilter('atom_site',
                                                        [{'name': 'auth_asym_id', 'type': 'str',
                                                          'alt_name': 'auth_chain_id', 'default': REPRESENTATIVE_ASYM_ID},
                                                         {'name': 'auth_seq_id', 'type': 'int', 'alt_name': 'auth_seq_id'},
                                                         {'name': 'label_comp_id', 'type': 'starts-with-alnum',
                                                          'alt_name': 'comp_id'},
                                                         {'name': 'label_atom_id', 'type': 'starts-with-alnum',
                                                          'alt_name': 'atom_id'},
                                                         {'name': 'Cartn_x', 'type': 'float', 'alt_name': 'x'},
                                                         {'name': 'Cartn_y', 'type': 'float', 'alt_name': 'y'},
                                                         {'name': 'Cartn_z', 'type': 'float', 'alt_name': 'z'},
                                                         {'name': 'type_symbol', 'type': 'str'}
                                                         ],
                                                        [{'name': 'Cartn_x', 'type': 'range-float',
                                                          'range': {'min_exclusive': (o[0] - CUTOFF_PARAMAGNETIC),
                                                                    'max_exclusive': (o[0] + CUTOFF_PARAMAGNETIC)}},
                                                         {'name': 'Cartn_y', 'type': 'range-float',
                                                          'range': {'min_exclusive': (o[1] - CUTOFF_PARAMAGNETIC),
                                                                    'max_exclusive': (o[1] + CUTOFF_PARAMAGNETIC)}},
                                                         {'name': 'Cartn_z', 'type': 'range-float',
                                                          'range': {'min_exclusive': (o[2] - CUTOFF_PARAMAGNETIC),
                                                                    'max_exclusive': (o[2] + CUTOFF_PARAMAGNETIC)}},
                                                         {'name': 'pdbx_PDB_model_num', 'type': 'int',
                                                          'value': self.__representative_model_id},
                                                         {'name': 'label_alt_id', 'type': 'enum',
                                                          'enum': (self.__representative_alt_id,)}
                                                         ])

        except Exception as e:  # pylint: disable=broad-exception-caught

            if self.__verbose:
                self.__log.write(f"+{self.__class_name__}.__getNearestParaFerroMagneticAtom() "
                                 f"++ Error  - {str(e)}\n")

            return None

        if len(_neighbor) == 0:
            return None

        neighbor = [n for n in _neighbor
                    if n['auth_seq_id'] != auth_seq_id
                    and distance(to_np_array(n), o) < CUTOFF_PARAMAGNETIC
                    and (n['type_symbol'] in PARAMAGNETIC_ELEMENTS
                         or n['type_symbol'] in FERROMAGNETIC_ELEMENTS)]

        if len(neighbor) == 0:
            return None

        atom_list = []

        for n in neighbor:
            atom_list.append({'auth_chain_id': n['auth_chain_id'], 'auth_seq_id': n['auth_seq_id'],
                              'comp_id': n['comp_id'], 'atom_id': n['atom_id'],
                              'distance': distance(to_np_array(n), o)})

        if len(atom_list) == 0:
            return None

        p = sorted(atom_list, key=itemgetter('distance'))[0]

        try:

            _p = self.__cR.getDictListWithFilter('atom_site',
                                                 [{'name': 'Cartn_x', 'type': 'float', 'alt_name': 'x'},
                                                  {'name': 'Cartn_y', 'type': 'float', 'alt_name': 'y'},
                                                  {'name': 'Cartn_z', 'type': 'float', 'alt_name': 'z'}
                                                  ],
                                                 [{'name': 'auth_asym_id', 'type': 'str', 'value': p['auth_chain_id']},
                                                  {'name': 'auth_seq_id', 'type': 'int', 'value': p['auth_seq_id']},
                                                  {'name': 'label_comp_id', 'type': 'str', 'value': p['comp_id']},
                                                  {'name': 'label_atom_id', 'type': 'str', 'value': p['atom_id']},
                                                  {'name': 'label_alt_id', 'type': 'enum',
                                                   'enum': (self.__representative_alt_id,)}
                                                  ])

        except Exception as e:  # pylint: disable=broad-exception-caught

            if self.__verbose:
                self.__log.write(f"+{self.__class_name__}.__getNearestParaFerroMagneticAtom() "
                                 f"++ Error  - {str(e)}\n")

            return None

        if len(_p) == 0:
            return None

        dist = 0.0

        for __p in _p:
            dist += distance(to_np_array(__p), o)

        p['distance'] = round(dist / len(_p), 1)

        return p

    def __validateDistanceRestraints(self) -> bool:
        """ Validate distance restraints.
            @author: Masashi Yokochi
            @note: Derived from wwpdb.apps.validation.src.RestraintsValidation.BMRBRestraintsAnalysis.calculate_distance_violations,
                   written by Kumaran Baskaran
            @change: class method,
                     support combinational restraints (_Gen_dist_constraint.Combination_ID, Member_ID)
        """

        if self.__coordinates is None:
            return False

        self.__distRestDictWithCombKey, self.__distRestViolDict, self.__distRestViolCombKeyDict = {}, {}, {}
        self.__distRestUnmapped = []

        if self.__distRestDict is None or self.__has_prev_results:
            return True

        try:

            def get_uninstanced_hydrogen_coord(model_id, atom_key):

                if atom_key in self.__coordinates[model_id]:
                    return None

                auth_asym_id, auth_seq_id, comp_id, atom_id, ins_code = atom_key

                if atom_id[0] not in PROTON_BEGIN_CODE:
                    return None

                if not self.__ccU.updateChemCompDict(comp_id):
                    return None

                if self.__ccU.lastChemCompDict['release_status'] != 'REL':
                    return None

                cca = next((cca for cca in self.__ccU.lastAtomDictList if cca['atom_id'] == atom_id), None)

                if cca is None:
                    return None

                bonded_to = self.__ccU.getBondedAtoms(comp_id, atom_id)

                if len(bonded_to) == 0:
                    return None

                ref_atom_id = bonded_to[0]

                _atom_key = (auth_asym_id, auth_seq_id, comp_id, ref_atom_id, ins_code)

                if _atom_key not in self.__coordinates[model_id]:
                    return None

                if not (cca['leaving_atom_flag'] != 'Y'
                        or (self.__csStat.peptideLike(comp_id)
                            and cca['n_terminal_atom_flag'] == 'N'
                            and cca['c_terminal_atom_flag'] == 'N')):
                    return None

                neighbor_to = self.__ccU.getBondedAtoms(comp_id, ref_atom_id, exclProton=True)

                if len(neighbor_to) == 0:
                    return None

                if len(neighbor_to) < 2:
                    for _ref_atom_id in neighbor_to:
                        _neighbor_to = self.__ccU.getBondedAtoms(comp_id, _ref_atom_id, exclProton=True)
                        for __ref_atom_id in _neighbor_to:
                            if __ref_atom_id not in neighbor_to:
                                neighbor_to.append(__ref_atom_id)
                                break
                        if len(neighbor_to) >= 2:
                            break

                if len(neighbor_to) < 2:
                    return None

                ref_atom_ids = [ref_atom_id]
                ref_atoms_xyz = [self.__coordinates[model_id][_atom_key]]

                for _ref_atom_id in neighbor_to:
                    __atom_key = (auth_asym_id, auth_seq_id, comp_id, _ref_atom_id, ins_code)
                    if __atom_key in self.__coordinates[model_id]:
                        ref_atom_ids.append(_ref_atom_id)
                        ref_atoms_xyz.append(self.__coordinates[model_id][__atom_key])

                src_ccd_xyz = numpy.asarray([cca['x'], cca['y'], cca['z']], dtype=float)

                ccd_atoms_xyz = []
                for _ref_atom_id in ref_atom_ids:
                    _cca = next((_cca for _cca in self.__ccU.lastAtomDictList if _cca['atom_id'] == _ref_atom_id), None)
                    if _cca is not None:
                        ccd_atoms_xyz.append(numpy.asarray([_cca['x'], _cca['y'], _cca['z']], dtype=float))

                if len(ref_atoms_xyz) != len(ccd_atoms_xyz):
                    return None

                dst_ccd_xyz, rmsd = calculate_uninstanced_coord(numpy.asarray(ccd_atoms_xyz),
                                                                numpy.asarray(ref_atoms_xyz),
                                                                numpy.asarray([src_ccd_xyz]))

                if rmsd > 0.1:
                    return None

                return dst_ccd_xyz[0]

            def calc_dist_rest_viol(rest_key, restraints):

                model_ids = list(self.__coordinates)              # stable order for stacking
                coords_by_model = [self.__coordinates[m] for m in model_ids]

                # per model: {bound_key: [d, ...]}
                dist_list_set_per_model = [{} for _ in model_ids]

                for r in restraints:
                    atom_key_1 = r['atom_key_1']
                    atom_key_2 = r['atom_key_2']
                    bound_key = (r['lower_bound'], r['upper_bound'],
                                 'or' if r['or_member'] else r['member_id'])

                    # per-model coordinates (direct, or reconstructed for uninstanced
                    # protons); collect the models where both atoms resolve
                    present_idx, p1_rows, p2_rows = [], [], []
                    for i, (model_id, coord) in enumerate(zip(model_ids, coords_by_model)):
                        pos_1 = coord.get(atom_key_1)
                        if pos_1 is None:
                            pos_1 = get_uninstanced_hydrogen_coord(model_id, atom_key_1)
                            if pos_1 is None and self.__verbose:
                                self.__log.write(f"Atom (auth_asym_id: {atom_key_1[0]}, auth_seq_id: {atom_key_1[1]}, "
                                                 f"comp_id: {atom_key_1[2]}, atom_id: {atom_key_1[3]}) "
                                                 f"not found in the coordinates for distance restraint {rest_key}.\n")
                        pos_2 = coord.get(atom_key_2)
                        if pos_2 is None:
                            pos_2 = get_uninstanced_hydrogen_coord(model_id, atom_key_2)
                            if pos_2 is None and self.__verbose:
                                self.__log.write(f"Atom (auth_asym_id: {atom_key_2[0]}, auth_seq_id: {atom_key_2[1]}, "
                                                 f"comp_id: {atom_key_2[2]}, atom_id: {atom_key_2[3]}) "
                                                 f"not found in the coordinates for distance restraint {rest_key}.\n")

                        if pos_1 is None or pos_2 is None:
                            self.__distRestUnmapped.append(rest_key)
                            continue

                        present_idx.append(i)
                        p1_rows.append(pos_1)
                        p2_rows.append(pos_2)

                    if len(present_idx) == 0:
                        continue

                    # one vectorized set of distances over present models (was a
                    # per-model scalar distance() call, i.e. models x restraints calls)
                    d_arr = numpy.linalg.norm(numpy.array(p1_rows) - numpy.array(p2_rows), axis=1)

                    for k, i in enumerate(present_idx):
                        d = d_arr[k]
                        if d == 0.0:
                            self.__log.write(f"+{self.__class_name__}.__validateDistanceRestraints() ++ Error  - "
                                             f"distance restraint {rest_key} {r} does not make sense, "
                                             f"{os.path.basename(self.__nmrDataPath)}.\n")
                        dist_list_set_per_model[i].setdefault(bound_key, []).append(d)

                error_per_model = {}

                for i, model_id in enumerate(model_ids):
                    error = None

                    for bound_key, dist_list in dist_list_set_per_model[i].items():
                        lower_bound, upper_bound, _ = bound_key
                        avr_d = dist_inv_6_summed(dist_list)

                        _error = dist_error(lower_bound, upper_bound, avr_d)

                        if error is None or error > _error:
                            error = _error

                    error_per_model[model_id] = error

                return error_per_model

            def fill_smaller_error_for_each_model(error_per_model, min_error_per_model,
                                                  combination_id, member_id, min_comb_key_per_model):
                for model_id in self.__eff_model_ids:
                    error = error_per_model[model_id]

                    if error is not None and error < min_error_per_model[model_id]:
                        min_error_per_model[model_id] = error
                        min_comb_key_per_model[model_id]['combination_id'] = combination_id
                        min_comb_key_per_model[model_id]['member_id'] = member_id

            def get_viol_per_model(min_error_per_model, min_comb_key_per_model):
                viol_per_model, comb_key_per_model = {}, {}

                for model_id in self.__eff_model_ids:
                    error = min_error_per_model[model_id]
                    comb_key = min_comb_key_per_model[model_id]

                    if NMR_VTF_DIST_VIOL_CUTOFF < error < DIST_ERROR_MAX:
                        viol_per_model[model_id] = round(error, 2)
                        comb_key_per_model[model_id] = (comb_key['combination_id'], comb_key['member_id'])
                    else:
                        viol_per_model[model_id] = None
                        comb_key_per_model[model_id] = None

                return viol_per_model, comb_key_per_model

            for rest_key, restraints in self.__distRestDict.items():

                has_combination_id = any(r['combination_id'] is not None for r in restraints)
                has_member_id = any(r['member_id'] is not None for r in restraints)

                self.__distRestDictWithCombKey[rest_key] = {}

                min_error_per_model = {model_id: DIST_ERROR_MAX for model_id in self.__eff_model_ids}
                min_comb_key_per_model = {model_id: {'combination_id': None, 'member_id': None}
                                          for model_id in self.__eff_model_ids}

                if not has_combination_id and not has_member_id:
                    error_per_model = calc_dist_rest_viol(rest_key, restraints)

                    fill_smaller_error_for_each_model(error_per_model, min_error_per_model,
                                                      None, None, min_comb_key_per_model)

                    self.__distRestDictWithCombKey[rest_key][(None, None)] = restraints

                elif not has_combination_id and has_member_id:
                    r_by_mem = {}
                    for r in restraints:
                        member_id = r['member_id']
                        if member_id not in r_by_mem:
                            r_by_mem[member_id] = []
                        r_by_mem[member_id].append(r)

                    for member_id, _restraints in r_by_mem.items():
                        _error_per_model = calc_dist_rest_viol(rest_key, _restraints)

                        fill_smaller_error_for_each_model(_error_per_model, min_error_per_model,
                                                          None, member_id, min_comb_key_per_model)

                        self.__distRestDictWithCombKey[rest_key][(None, member_id)] = _restraints

                elif has_combination_id and not has_member_id:
                    r_by_comb = {}
                    for r in restraints:
                        combination_id = r['combination_id']
                        if combination_id not in r_by_comb:
                            r_by_comb[combination_id] = []
                        r_by_comb[combination_id].append(r)

                    for combination_id, _restraints in r_by_comb.items():
                        _error_per_model = calc_dist_rest_viol(rest_key, _restraints)

                        fill_smaller_error_for_each_model(_error_per_model, min_error_per_model,
                                                          combination_id, None, min_comb_key_per_model)

                        self.__distRestDictWithCombKey[rest_key][(combination_id, None)] = _restraints

                else:
                    r_by_comb_and_mem = {}
                    for r in restraints:
                        combination_id = r['combination_id']
                        member_id = r['member_id']
                        if combination_id not in r_by_comb_and_mem:
                            r_by_comb_and_mem[combination_id] = {}
                        if member_id not in r_by_comb_and_mem[combination_id]:
                            r_by_comb_and_mem[combination_id][member_id] = []
                        r_by_comb_and_mem[combination_id][member_id].append(r)

                    for combination_id, _r_by_mem in r_by_comb_and_mem.items():
                        for member_id, _restraints in _r_by_mem.items():
                            _error_per_model = calc_dist_rest_viol(rest_key, _restraints)

                            fill_smaller_error_for_each_model(_error_per_model, min_error_per_model,
                                                              combination_id, member_id,
                                                              min_comb_key_per_model)

                            self.__distRestDictWithCombKey[rest_key][(combination_id, member_id)] = _restraints

                self.__distRestViolDict[rest_key], self.__distRestViolCombKeyDict[rest_key] =\
                    get_viol_per_model(min_error_per_model, min_comb_key_per_model)

            self.__distRestUnmapped = list(set(self.__distRestUnmapped))

            return True

        except Exception as e:  # pylint: disable=broad-exception-caught
            self.__log.write(f"Exception occurred while processing {os.path.basename(self.__cifPath)} "
                             f"and {os.path.basename(self.__nmrDataPath)}\n")
            self.__log.write(f"+{self.__class_name__}.__validateDistanceRestraints() ++ Error  - {str(e)}\n")

            self.__distRestViolDict = self.__distRestUnmapped = None

        return False

    def __validateDihedralAngleRestraints(self) -> bool:
        """ Validate dihedral angle restraints.
            @author: Masashi Yokochi
            @note: Derived from wwpdb.apps.validation.src.RestraintsValidation.BMRBRestraintsAnalysis.calculate_angle_violations,
                   written by Kumaran Baskaran
            @change: class method,
                     support combinational restraints (_Torsion_angle_constraint.Combination_ID)
        """

        if self.__coordinates is None:
            return False

        self.__dihedRestDictWithCombKey, self.__dihedRestViolDict, self.__dihedRestViolCombKeyDict = {}, {}, {}
        self.__dihedRestUnmapped = []

        if self.__dihedRestDict is None or self.__has_prev_results:
            return True

        try:

            def calc_dihed_rest_viol(rest_key, restraints):
                model_ids = list(self.__coordinates)              # stable order for stacking
                coords_by_model = [self.__coordinates[m] for m in model_ids]

                # per model: {bound_key: [a = dihedral + 180, ...]}
                angle_list_set_per_model = [{} for _ in model_ids]

                for r in restraints:
                    atom_keys = (r['atom_key_1'], r['atom_key_2'], r['atom_key_3'], r['atom_key_4'])
                    bound_key = (r['lower_bound'], r['upper_bound'], r['target_value'])

                    # per-model presence; collect the model indices where all four atoms exist
                    present_idx = []
                    for i, coord in enumerate(coords_by_model):
                        atom_present = True
                        for atom_key in atom_keys:
                            if atom_key not in coord:
                                atom_present = False
                                if self.__verbose:
                                    self.__log.write(f"Atom (auth_asym_id: {atom_key[0]}, auth_seq_id: {atom_key[1]}, "
                                                     f"comp_id: {atom_key[2]}, atom_id: {atom_key[3]}) "
                                                     f"not found in the coordinates for dihedral angle restraint {rest_key}.\n")
                        if atom_present:
                            present_idx.append(i)
                        else:
                            self.__dihedRestUnmapped.append(rest_key)

                    if len(present_idx) == 0:
                        continue

                    # one vectorized dihedral over all present models (was a per-model
                    # scalar dihedral_angle() call, i.e. models x restraints calls)
                    p0, p1, p2, p3 = (numpy.array([coords_by_model[i][ak] for i in present_idx])
                                      for ak in atom_keys)
                    angles = dihedral_angles(p0, p1, p2, p3) + 180.0

                    for k, i in enumerate(present_idx):
                        angle_list_set_per_model[i].setdefault(bound_key, []).append(angles[k])

                error_per_model = {}

                for i, model_id in enumerate(model_ids):
                    error = None

                    for bound_key, angle_list in angle_list_set_per_model[i].items():
                        lower_bound, upper_bound, target_value = bound_key
                        avr_a = numpy.mean(numpy.array(angle_list, dtype=float)) - 180.0

                        _error = angle_error(lower_bound, upper_bound, target_value, avr_a)

                        if error is None or error > _error:
                            error = _error

                    error_per_model[model_id] = error

                return error_per_model

            def fill_smaller_error_for_each_model(error_per_model, min_error_per_model,
                                                  combination_id, min_comb_key_per_model):
                for model_id in self.__eff_model_ids:
                    error = error_per_model[model_id]

                    if error is not None and error < min_error_per_model[model_id]:
                        min_error_per_model[model_id] = error
                        min_comb_key_per_model[model_id]['combination_id'] = combination_id

            def get_viol_per_model(min_error_per_model, min_comb_key_per_model):
                viol_per_model, comb_key_per_model = {}, {}

                for model_id in self.__eff_model_ids:
                    error = min_error_per_model[model_id]
                    comb_key = min_comb_key_per_model[model_id]

                    if NMR_VTF_DIHED_VIOL_CUTOFF < error < ANGLE_ERROR_MAX:
                        viol_per_model[model_id] = round(error, 2)
                        comb_key_per_model[model_id] = (comb_key['combination_id'],)
                    else:
                        viol_per_model[model_id] = None
                        comb_key_per_model[model_id] = None

                return viol_per_model, comb_key_per_model

            for rest_key, restraints in self.__dihedRestDict.items():

                has_combination_id = any(r['combination_id'] is not None for r in restraints)

                self.__dihedRestDictWithCombKey[rest_key] = {}

                min_error_per_model = {model_id: ANGLE_ERROR_MAX for model_id in self.__eff_model_ids}
                min_comb_key_per_model = {model_id: {'combination_id': None}
                                          for model_id in self.__eff_model_ids}

                if not has_combination_id:
                    error_per_model = calc_dihed_rest_viol(rest_key, restraints)

                    fill_smaller_error_for_each_model(error_per_model, min_error_per_model,
                                                      None, min_comb_key_per_model)

                    self.__dihedRestDictWithCombKey[rest_key][(None,)] = restraints

                else:
                    r_by_comb = {}
                    for r in restraints:
                        combination_id = r['combination_id']
                        if combination_id not in r_by_comb:
                            r_by_comb[combination_id] = []
                        r_by_comb[combination_id].append(r)

                    for combination_id, _restraints in r_by_comb.items():
                        _error_per_model = calc_dihed_rest_viol(rest_key, _restraints)

                        fill_smaller_error_for_each_model(_error_per_model, min_error_per_model,
                                                          combination_id, min_comb_key_per_model)

                        self.__dihedRestDictWithCombKey[rest_key][(combination_id,)] = _restraints

                self.__dihedRestViolDict[rest_key], self.__dihedRestViolCombKeyDict[rest_key] =\
                    get_viol_per_model(min_error_per_model, min_comb_key_per_model)

            self.__dihedRestUnmapped = list(set(self.__dihedRestUnmapped))

            return True

        except Exception as e:  # pylint: disable=broad-exception-caught
            self.__log.write(f"Exception occurred while processing {os.path.basename(self.__cifPath)} "
                             f"and {os.path.basename(self.__nmrDataPath)}\n")
            self.__log.write(f"+{self.__class_name__}.__validateDihedralAngleRestraints() ++ Error  - {str(e)}\n")

            self.__dihedRestViolDict = self.__dihedRestUnmapped = None

        return False

    def __validateRdcRestraints(self) -> bool:
        """ Validate RDC restraints.
            @author: Masashi Yokochi
        """

        if self.__coordinates is None:
            return False

        self.__rdcRestDictWithCombKey, self.__rdcRestViolDict, self.__rdcRestViolCombKeyDict, \
            self.__rdcCorrPlotDict = {}, {}, {}, {}  # pylint: disable='unreachable'
        self.__rdcRestUnmapped = []

        if self.__rdcRestDict is None or self.__has_prev_results:
            return True

        try:

            def calc_rdc_rest_viol(rest_key, restraints):
                error_per_model = {}

                # __rdcCalcDict[rest_key][model_id] is one calculated RDC per model,
                # independent of the individual restraint r, so averaging a bound_key
                # list (copies of the same value) == that value; the pos_1/pos_2
                # fetches were only atom-presence checks. Hoist the model-independent
                # rest_key lookup and drop the per-(model, bound_key) numpy.mean.
                calc_per_model = self.__rdcCalcDict.get(rest_key, {})

                for model_id in self.__coordinates:

                    coords = self.__coordinates[model_id]
                    has_calc = model_id in calc_per_model
                    r_calc = calc_per_model[model_id] if has_calc else None

                    error = None
                    seen_bound_keys = set()

                    for r in restraints:
                        atom_key_1 = r['atom_key_1']
                        atom_key_2 = r['atom_key_2']

                        atom_present = True

                        if atom_key_1 not in coords:
                            if self.__verbose:
                                self.__log.write(f"Atom (auth_asym_id: {atom_key_1[0]}, auth_seq_id: {atom_key_1[1]}, "
                                                 f"comp_id: {atom_key_1[2]}, atom_id: {atom_key_1[3]}) "
                                                 f"not found in the coordinates for RDC restraint {rest_key}.\n")
                            atom_present = False

                        if atom_key_2 not in coords:
                            if self.__verbose:
                                self.__log.write(f"Atom (auth_asym_id: {atom_key_2[0]}, auth_seq_id: {atom_key_2[1]}, "
                                                 f"comp_id: {atom_key_2[2]}, atom_id: {atom_key_2[3]}) "
                                                 f"not found in the coordinates for RDC restraint {rest_key}.\n")
                            atom_present = False

                        if not atom_present:
                            self.__rdcRestUnmapped.append(rest_key)
                            continue

                        if not has_calc:
                            continue

                        bound_key = (r['lower_bound'], r['upper_bound'])

                        if bound_key in seen_bound_keys:  # same r_calc -> same error
                            continue
                        seen_bound_keys.add(bound_key)

                        # experimental
                        if bound_key[0] != bound_key[1] or rest_key not in self.__rdcSyntCalcDict:
                            _error = rdc_error(bound_key[0], bound_key[1], r_calc)

                        # theoretical
                        else:
                            rdc_calc_min = rdc_calc_max = None
                            _rdc_synt_calcs = []
                            for v in self.__rdcSyntCalcDict[rest_key].values():
                                _rdc_synt_calcs.extend(v)
                            if len(_rdc_synt_calcs) > RDC_MIN_MC_CYCLES:
                                rdc_synt_calcs = numpy.array(_rdc_synt_calcs, dtype=float) / r['scale_factor']
                                rdc_calc_min = round(numpy.min(rdc_synt_calcs), 2)
                                rdc_calc_max = round(numpy.max(rdc_synt_calcs), 2)
                            if None in (rdc_calc_min, rdc_calc_max):
                                _error = rdc_error(rdc_calc_min, rdc_calc_max, r_calc)
                            else:
                                _error = rdc_error(bound_key[0], bound_key[1], r_calc)

                        if error is None or error > _error:
                            error = _error

                    error_per_model[model_id] = error

                return error_per_model

            def fill_smaller_error_for_each_model(error_per_model, min_error_per_model,
                                                  combination_id, min_comb_key_per_model):
                for model_id in self.__eff_model_ids:
                    error = error_per_model[model_id]

                    if error is not None and error < min_error_per_model[model_id]:
                        min_error_per_model[model_id] = error
                        min_comb_key_per_model[model_id]['combination_id'] = combination_id

            def get_viol_per_model(min_error_per_model, min_comb_key_per_model):
                viol_per_model, comb_key_per_model = {}, {}

                for model_id in self.__eff_model_ids:
                    error = min_error_per_model[model_id]
                    comb_key = min_comb_key_per_model[model_id]

                    if NMR_VTF_RDC_VIOL_CUTOFF < error < RDC_ERROR_MAX:
                        viol_per_model[model_id] = round(error, 2)
                        comb_key_per_model[model_id] = (comb_key['combination_id'],)
                    else:
                        viol_per_model[model_id] = None
                        comb_key_per_model[model_id] = None

                return viol_per_model, comb_key_per_model

            for rest_key, restraints in self.__rdcRestDict.items():

                has_combination_id = any(r['combination_id'] is not None for r in restraints)

                self.__rdcRestDictWithCombKey[rest_key] = {}

                min_error_per_model = {model_id: RDC_ERROR_MAX for model_id in self.__eff_model_ids}
                min_comb_key_per_model = {model_id: {'combination_id': None}
                                          for model_id in self.__eff_model_ids}

                if not has_combination_id:
                    error_per_model = calc_rdc_rest_viol(rest_key, restraints)

                    fill_smaller_error_for_each_model(error_per_model, min_error_per_model,
                                                      None, min_comb_key_per_model)

                    self.__rdcRestDictWithCombKey[rest_key][(None,)] = restraints

                else:
                    r_by_comb = {}
                    for r in restraints:
                        combination_id = r['combination_id']
                        if combination_id not in r_by_comb:
                            r_by_comb[combination_id] = []
                        r_by_comb[combination_id].append(r)

                    for combination_id, _restraints in r_by_comb.items():
                        _error_per_model = calc_rdc_rest_viol(rest_key, _restraints)

                        fill_smaller_error_for_each_model(_error_per_model, min_error_per_model,
                                                          combination_id, min_comb_key_per_model)

                        self.__rdcRestDictWithCombKey[rest_key][(combination_id,)] = _restraints

                self.__rdcRestViolDict[rest_key], self.__rdcRestViolCombKeyDict[rest_key] =\
                    get_viol_per_model(min_error_per_model, min_comb_key_per_model)

            self.__rdcRestUnmapped = list(set(self.__rdcRestUnmapped))

            list_ids = set()
            for rest_key in self.__rdcCalcDict:
                list_ids.add(rest_key[0])

            for list_id in sorted(list(list_ids)):
                rdc_values, rdc_errors, rdc_viols, q_scores = {}, {}, {}, {}

                for rest_key, rdc_calc in self.__rdcCalcDict.items():

                    if rest_key[0] != list_id:
                        continue

                    restraints = self.__rdcRestDict[rest_key]

                    for r in restraints:
                        rdc_type = r['rdc_type']

                        if rdc_type not in rdc_values:
                            rdc_values[rdc_type], rdc_errors[rdc_type], rdc_viols[rdc_type] = [], [], []
                            q_scores[rdc_type] = {'rdc_exp': [], 'rdc_calc': []}

                        ak1 = r['atom_key_1']
                        ak2 = r['atom_key_2']

                        rdc_exp_center = r['target_value']  # raw observed RDC value (without scale factor)

                        # since calculated RDC values have been scaled through scaled RDC input data at SVD operation
                        # cancel scale factor of the calculated RDC values to compare with the raw observed RDC values
                        rdc_calcs = numpy.array(list(rdc_calc.values()), dtype=float) / r['scale_factor']
                        rdc_calc_mean = numpy.mean(rdc_calcs)
                        rdc_calc_center = round(rdc_calc_mean, 2)

                        # uncertainty of calculated RDC values is estimated by Monte Carlo simulation
                        rdc_calc_min = rdc_calc_max = None
                        if rest_key in self.__rdcSyntCalcDict:
                            _rdc_synt_calcs = []
                            for v in self.__rdcSyntCalcDict[rest_key].values():
                                _rdc_synt_calcs.extend(v)
                            if len(_rdc_synt_calcs) > RDC_MIN_MC_CYCLES:
                                rdc_synt_calcs = numpy.array(_rdc_synt_calcs, dtype=float) / r['scale_factor']

                                rdc_calc_min = round(numpy.min(rdc_synt_calcs), 2)
                                rdc_calc_max = round(numpy.max(rdc_synt_calcs), 2)

                        q_scores[rdc_type]['rdc_exp'].append(rdc_exp_center)
                        q_scores[rdc_type]['rdc_calc'].append(rdc_calc_mean)

                        if ak1[0] == ak2[0]:
                            if ak1[1] == ak2[1]:
                                vector_name = f"{ak1[0]}:{ak1[1]}:{ak1[2]}:{ak1[3]}-{ak2[3]}"
                            else:
                                vector_name = f"{ak1[0]}:{ak1[1]}:{ak1[2]}:{ak1[3]}-"\
                                    f"{ak2[1]}:{ak2[2]}:{ak2[3]}"
                        else:
                            vector_name = f"{ak1[0]}:{ak1[1]}:{ak1[2]}:{ak1[3]}-"\
                                f"{ak2[0]}:{ak2[1]}:{ak2[2]}:{ak2[3]}"

                        rdc_values[rdc_type].append([rdc_exp_center, rdc_calc_center, vector_name])
                        rdc_errors[rdc_type].append([rdc_exp_center, rdc_calc_center,
                                                     r['lower_bound'], r['upper_bound'],
                                                     rdc_calc_min,
                                                     rdc_calc_max])

                        if None not in (rdc_calc_min, rdc_calc_max)\
                           and (rdc_exp_center < rdc_calc_min - NMR_VTF_RDC_ERR_BINS[1]
                                or rdc_exp_center > rdc_calc_max + NMR_VTF_RDC_ERR_BINS[1]):
                            rdc_viols[rdc_type].append([rdc_exp_center, rdc_calc_center, vector_name])

                da_array = numpy.array([abs(float(v['Szz'])) * float(v['Dmax'])
                                        for v in self.__rdcSaupeOrderMatrix[list_id].values()], dtype=float)
                eta_array = numpy.array([abs(float(v['eta'])) for v in self.__rdcSaupeOrderMatrix[list_id].values()], dtype=float)
                denominator_unit = 0.4 * da_array.mean() ** 2 * (4.0 + 3.0 * eta_array.mean() ** 2)

                for k, v in copy.copy(q_scores).items():
                    rdc_exp_array = numpy.array(v['rdc_exp'], dtype=float)
                    rdc_calc_array = numpy.array(v['rdc_calc'], dtype=float)
                    rdc_exp_mean = numpy.mean(rdc_exp_array)
                    total_sum_of_square = ((rdc_exp_array - rdc_exp_mean) ** 2).sum()
                    sum_of_squared_errors = ((rdc_exp_array - rdc_calc_array) ** 2).sum()
                    sum_of_squared_values = (rdc_exp_array ** 2).sum()

                    if rdc_exp_array.std() == 0.0 or rdc_calc_array.std() == 0.0:
                        q_scores[k]['r'] = None
                    else:
                        m = numpy.corrcoef(rdc_exp_array, rdc_calc_array)
                        q_scores[k]['r'] = round(m[0, 1], 2)
                    q_scores[k]['r2'] =\
                        round(1.0 - sum_of_squared_errors / total_sum_of_square, 2)\
                        if total_sum_of_square > 0.0 else 1.0
                    q_scores[k]['cornilescu_q'] =\
                        round(math.sqrt(sum_of_squared_errors / sum_of_squared_values), 2)\
                        if sum_of_squared_values > 0.0 else 0.0
                    q_scores[k]['clore_q'] =\
                        round(math.sqrt(sum_of_squared_errors / (rdc_exp_array.shape[0] * denominator_unit)), 2)\
                        if denominator_unit > 0.0 else 0.0

                    del q_scores[k]['rdc_exp']
                    del q_scores[k]['rdc_calc']

                self.__rdcCorrPlotDict[list_id] = {'values': rdc_values, 'errors': rdc_errors,
                                                   'q_scores': q_scores}

                if len(rdc_viols) > 0:
                    self.__rdcCorrPlotDict[list_id]['violations'] = rdc_viols

            return True

        except Exception as e:  # pylint: disable=broad-exception-caught
            self.__log.write(f"Exception occurred while processing {os.path.basename(self.__cifPath)} "
                             f"and {os.path.basename(self.__nmrDataPath)}\n")
            self.__log.write(f"+{self.__class_name__}.__validateRdcRestraints() ++ Error  - {str(e)}\n")

            self.__rdcRestViolDict = self.__rdcRestUnmapped = None

        return False

    def __summarizeCommonCsAnalysis(self) -> bool:
        """ Summarize common chemical shift analysis.
            @author: Masashi Yokochi
            @note: Derived from wwpdb.apps.validation.src.ChemicalShiftsValidation.BMRBChemicalShiftsAnalysis.getCompleteness,
                   written by Aleksandras Gutmanas, Kumaran Baskaran
            @change: Add support for Phosphorus chemical shift, completeness of assigned chemical shifts in favorable ranges
        """

        if self.__has_prev_results:
            return True

        if self.__chemShiftDict is None or len(self.__chemShiftDict) == 0:
            self.__results = {'cs_no_file': True}
            if not self.__permit_missing_chem_shift:
                self.__log.write("There are no assigned chemical shifts that need to be analyzed.\n")
            return False

        has_coord = self.__coordinates is not None

        # metadata, bookkeeping

        self.__results = {'meta_data': (os.path.basename(self.__nmrDataPath), self.__chemShiftMeta),
                          'book_keeping': {'total_shifts': self.__chemShiftTotal,
                                           'cs_error': {'CS_OUTLIER': self.__chemShiftOutlier,
                                                        'CS_DUPLICATE': self.__chemShiftDuplicated,
                                                        'CS_VALUE': self.__chemShiftUnparsed,
                                                        'NO_MAP': self.__chemShiftUnmapped,
                                                        'NO_MODEL': self.__chemShiftUnmodeled
                                                        }
                                           }
                          }

        # solid-state NMR

        if has_coord:
            exptl = self.__cR.getDictList('exptl')

            self.__results['ssnmr'] =\
                exptl[0]['method'] == 'SOLID-STATE NMR' if len(exptl) > 0 and 'method' in exptl[0] else False

        else:
            self.__results['ssnmr'] = False

            if self.__rR.hasItem('Entry', 'Experimental_method_subtype'):
                entry = self.__rR.getDictList('Entry')

                self.__results['ssnmr'] =\
                    entry[0]['Experimental_method_subtype'] == 'solid-state'\
                    if len(entry) > 0 and 'Experimental_method_subtype' in entry[0] else False

        # completeness of assigned chemical shifts

        any_type = 'Total'

        task_keys = ('well_defined', 'full_length')
        atom_types = (any_type, 'H', 'C', 'N', 'P')
        atom_categories = ('overall', 'favorable', 'backbone', 'sidechain', 'aromatic', 'sugar', 'base')

        def do_calc_completeness(list_ids, task_key, auth_chain_id, well_defined_region, seq_key, seq_id, coord_atoms, output):
            if task_key == 'well_defined':
                try:
                    auth_seq_id = int(seq_key[0])
                    if auth_seq_id not in well_defined_region:
                        return
                except ValueError:
                    return

            res_num = seq_key[0]
            comp_id = seq_key[1]

            cs_atoms = []
            for list_id in list_ids:
                if list_id not in self.__chemShiftUniqDict:
                    continue
                cs_atoms.extend([cs_key[3] for cs_key, cs_val in self.__chemShiftUniqDict[list_id].items()
                                 if cs_key[0] == auth_chain_id and cs_key[1] == res_num and cs_key[2] == comp_id
                                 and 'primary' in cs_val])
            if len(list_ids) > 1:
                cs_atoms = list(set(cs_atoms))

            categorized_atom_ids = self.__csStat.getCategorizedAtomIds(comp_id)

            if categorized_atom_ids is None:
                return

            n_terminal_aa = seq_id == 1 and self.__csStat.peptideLike(comp_id)

            for atom_category, atom_ids in categorized_atom_ids.items():
                for atom_id in atom_ids:

                    if n_terminal_aa and atom_id == 'H':
                        continue

                    atom_type = atom_id[0]
                    if atom_type not in ('H', 'C', 'N', 'P'):
                        continue

                    if atom_id in cs_atoms:
                        output[task_key][atom_type][atom_category][0] += 1
                        output[task_key][atom_type]['overall'][0] += 1
                        output[task_key][any_type][atom_category][0] += 1
                        output[task_key][any_type]['overall'][0] += 1

                        if not any(any(True for out_cs in self.__chemShiftOutlier[list_id]
                                   if out_cs[0] == auth_chain_id and out_cs[1] == res_num
                                   and out_cs[2] == comp_id and out_cs[3] == atom_id) for list_id in list_ids):
                            output[task_key][atom_type]['favorable'][0] += 1
                            output[task_key][any_type]['favorable'][0] += 1

                    # count only nitrogens if not quaternary
                    # use coordinates to check whether hydrogen attached
                    elif atom_type == 'N' and atom_id in self.__csStat.getQuaternaryNitrogensIfNoProtonIsBonded(comp_id):
                        bonded_protons = self.__ccU.getBondedAtoms(comp_id, atom_id, onlyProton=True)
                        if len(bonded_protons) == 1 and bonded_protons[0] not in coord_atoms:
                            continue

                    # count only phosphorus having chemical shift
                    elif atom_type == 'P':
                        continue

                    output[task_key][atom_type][atom_category][1] += 1
                    output[task_key][atom_type]['overall'][1] += 1
                    output[task_key][atom_type]['favorable'][1] += 1
                    output[task_key][any_type][atom_category][1] += 1
                    output[task_key][any_type]['overall'][1] += 1
                    output[task_key][any_type]['favorable'][1] += 1

            # Check for stereo-assignment of methyl groups (VAL, LEU)
            for list_id in list_ids:
                methyl_carbons = [a for a in self.__csStat.getMethylAtoms(comp_id) if a.startswith('C')]

                if len(methyl_carbons) < 2:
                    continue

                geminal_methyl_carbons = [a for a in methyl_carbons if self.__csStat.getMaxAmbigCodeWoSetId(comp_id, a) == 2]

                if len(geminal_methyl_carbons) < 2:
                    continue

                has_all_shifts = True
                for a in geminal_methyl_carbons:
                    _cs_key = (auth_chain_id, res_num, comp_id, a)
                    if _cs_key in self.__chemShiftUniqDict[list_id]\
                       and self.__chemShiftUniqDict[list_id][_cs_key]['ambig_code'] == 1:
                        continue
                    has_all_shifts = False
                    break

                if has_all_shifts:
                    output[task_key]['stereomethyl'][0] += 1

                output[task_key]['stereomethyl'][1] += 1

        if self.__inputReport is None and has_coord:
            polySeq = []

            polySeqPdbMonIdName = 'pdb_mon_id' if self.__cR.hasItem('pdbx_poly_seq_scheme', 'pdb_mon_id') else 'mon_id'

            lpCategory = 'pdbx_poly_seq_scheme'
            keyItems = [{'name': 'asym_id', 'type': 'str', 'alt_name': 'chain_id',
                         'default': REPRESENTATIVE_ASYM_ID},
                        {'name': 'seq_id', 'type': 'int', 'alt_name': 'seq_id'},
                        {'name': 'mon_id', 'type': 'starts-with-alnum', 'alt_name': 'comp_id'},
                        {'name': 'pdb_strand_id', 'type': 'str', 'alt_name': 'auth_chain_id',
                         'default': REPRESENTATIVE_ASYM_ID},
                        {'name': 'pdb_seq_num', 'type': 'int', 'alt_name': 'auth_seq_id'},
                        {'name': polySeqPdbMonIdName, 'type': 'str', 'alt_name': 'auth_comp_id',
                         'default-from': 'mon_id'}
                        ]

            try:

                if self.__cR.hasItem(lpCategory, 'pdb_ins_code'):
                    keyItems.append({'name': 'pdb_ins_code', 'type': 'str', 'alt_name': 'ins_code', 'default': '.'})

                if self.__cR.hasItem(lpCategory, 'auth_mon_id'):
                    keyItems.append({'name': 'auth_mon_id', 'type': 'str', 'alt_name': 'alt_comp_id', 'default-from': 'mon_id'})

                try:
                    polySeq, _ = self.__cR.getPolymerSequence(lpCategory, keyItems,
                                                              withRmsd=True,
                                                              totalModels=self.__total_models,
                                                              effModelIds=self.__eff_model_ids,
                                                              repAltId=self.__representative_alt_id)
                except KeyError:  # pdbx_PDB_ins_code throws KeyError
                    pass

                if len(polySeq) == 0:
                    lpCategory = 'atom_site'
                    keyItems = [{'name': 'auth_asym_id', 'type': 'str', 'alt_name': 'auth_chain_id',
                                 'default': REPRESENTATIVE_ASYM_ID},
                                {'name': 'label_asym_id', 'type': 'str', 'alt_name': 'chain_id'},
                                {'name': 'auth_seq_id', 'type': 'int', 'alt_name': 'auth_seq_id'},
                                {'name': 'label_seq_id', 'type': 'str', 'alt_name': 'seq_id',
                                 'default-from': 'auth_seq_id'},
                                {'name': 'auth_comp_id', 'type': 'int', 'alt_name': 'auth_comp_id'},
                                {'name': 'label_comp_id', 'type': 'starts-with-alnum', 'alt_name': 'comp_id'}
                                ]

                    if self.__cR.hasItem(lpCategory, 'pdbx_PDB_ins_code'):
                        keyItems.append({'name': 'pdbx_PDB_ins_code', 'type': 'str', 'alt_name': 'ins_code', 'default': '.'})

                    try:
                        polySeq, _ = self.__cR.getPolymerSequence(lpCategory, keyItems,
                                                                  withRmsd=True,
                                                                  totalModels=self.__total_models,
                                                                  effModelIds=self.__eff_model_ids,
                                                                  repAltId=self.__representative_alt_id)
                    except (KeyError, ValueError):
                        pass

                for ps in polySeq:
                    if 'ins_code' in ps and len(collections.Counter(ps['ins_code']).most_common()) == 1:
                        del ps['ins_code']

            except Exception as e:  # pylint: disable=broad-exception-caught
                if self.__verbose:
                    self.__log.write(f"{self.__class_name__}.__summarizeCommonCsAnalysis() ++ Error  - {str(e)}\n")

        def calc_completeness(list_ids):
            output = {}
            for task_key in task_keys:
                output[task_key] = {'stereomethyl': [0, 0]}
                for atom_type in atom_types:
                    output[task_key][atom_type] = {}
                    for atom_category in atom_categories:
                        output[task_key][atom_type][atom_category] = [0, 0]

            for auth_chain_id, v in self.__entityInstance.items():
                well_defined_region = []

                if self.__inputReport is None and has_coord:
                    cif_ps = next((ps for ps in polySeq if ps['auth_chain_id'] == auth_chain_id), None)

                    if cif_ps is not None:
                        if 'well_defined_region' in cif_ps:
                            for item in cif_ps['well_defined_region']:
                                well_defined_region.extend(item['seq_id'])

                else:
                    cif_ps = self.__inputReport.getModelPolymerSequenceOf(auth_chain_id, label_scheme=False)\
                        if has_coord else self.__inputReport.getNmrPolymerSequenceOf(auth_chain_id)

                    if cif_ps is not None:
                        if 'well_defined_region' in cif_ps:
                            for item in cif_ps['well_defined_region']:
                                well_defined_region.extend(item['seq_id'])

                for seq_key, w in v.items():
                    for task_key in task_keys:
                        do_calc_completeness(list_ids, task_key, auth_chain_id, well_defined_region,
                                             seq_key, w['seq_id'], w['atoms'], output)

            return output

        list_ids = list(self.__chemShiftDict.keys())

        _completeness = calc_completeness(list_ids)

        completeness = {task_key: _completeness[task_key][any_type]['overall'] for task_key in task_keys}

        completeness['chemical_shift_completeness'] =\
            round(100.0 * completeness['well_defined'][0] / completeness['well_defined'][1], 2)\
            if completeness['well_defined'][1] > 0 else 0.0

        completeness['chemical_shift_completeness_full_length'] =\
            round(100.0 * completeness['full_length'][0] / completeness['full_length'][1], 2)\
            if completeness['full_length'][1] > 0 else 0.0

        favor_completeness = {task_key: _completeness[task_key][any_type]['favorable'] for task_key in task_keys}

        completeness['favor_well_defined'] = _completeness['well_defined'][any_type]['favorable']
        completeness['favor_full_length'] = _completeness['full_length'][any_type]['favorable']

        completeness['favor_chem_shift_completeness'] =\
            round(100.0 * favor_completeness['well_defined'][0] / favor_completeness['well_defined'][1], 2)\
            if favor_completeness['well_defined'][1] > 0 else 0.0

        completeness['favor_chem_shift_completeness_full_length'] =\
            round(100.0 * favor_completeness['full_length'][0] / favor_completeness['full_length'][1], 2)\
            if favor_completeness['full_length'][1] > 0 else 0.0

        self.__results['completeness'] = completeness

        self.__results['completeness_items'] = {list_id: calc_completeness([list_id]) for list_id in list_ids}

        # shift summary table

        shift_summary_table = {}

        for idx, list_id in enumerate(list_ids):
            parsed = len(self.__chemShiftDict[list_id])
            unparsed = len(self.__chemShiftUnparsed[list_id])
            unmapped_error = len(self.__chemShiftUnmapped[list_id])
            unmapped_warning = len(self.__chemShiftUnmodeled[list_id])
            mapped = parsed - unmapped_error - unmapped_warning
            shift_summary_table[list_id] = {'block_id': self.__chemShiftMeta[idx][0],
                                            'block_name': self.__chemShiftMeta[idx][1],
                                            'list_id': self.__chemShiftMeta[idx][2],
                                            'file_name': os.path.basename(self.__nmrDataPath),
                                            # 'type': 'full',
                                            'total_number_of_shifts': parsed + unparsed,
                                            'number_of_parsed_shifts': parsed,
                                            'number_of_unparsed_shifts': unparsed,
                                            'number_of_mapped_shifts': mapped,
                                            'number_of_errors_while_mapping': unmapped_error,
                                            'number_of_warnings_while_mapping': unmapped_warning
                                            }

        self.__results['shift_summary_table'] = shift_summary_table

        # random coil index

        rci_result = {}

        rci = RCI()

        rci_atom_ids = ('HA', 'HA1', 'HA2', 'HA3', 'H', 'HN', 'NH', 'C', 'CO', 'N', 'CA', 'CB')

        has_rci_results = False

        for list_id, cs_data in self.__chemShiftUniqDict.items():
            rci_result[list_id] = {}

            auth_chain_ids = list(dict.fromkeys(cs_key[0] for cs_key in cs_data))

            # assigned shifts of each residue, {(auth_chain_id, auth_seq_id, comp_id): {cs_key: cs_vals}}, in the order of
            # cs_data, so that each residue looks its shifts up instead of scanning all shifts of the list
            cs_data_by_residue = {}
            for cs_key, cs_vals in cs_data.items():
                if cs_key[1].isdigit() and cs_vals['value'] not in EMPTY_VALUE:
                    residue_key = (cs_key[0], int(cs_key[1]), cs_key[2])
                    if residue_key not in cs_data_by_residue:
                        cs_data_by_residue[residue_key] = {}
                    cs_data_by_residue[residue_key][cs_key] = cs_vals

            for auth_chain_id in auth_chain_ids:
                if has_coord:
                    ps = next((ps for ps in self.__caC['polymer_sequence'] if ps['auth_chain_id'] == auth_chain_id), None)
                elif self.__inputReport is not None:
                    ps = self.__inputReport.getNmrPolymerSequenceOf(auth_chain_id)

                if ps is None:
                    continue

                rci_residues, rci_assignments, seq_ids_wo_assign, oxidized_cys_seq_ids = [], [], [], []

                for auth_seq_id, comp_id in zip(ps['auth_seq_id'] if 'auth_seq_id' in ps else ps['seq_id'], ps['comp_id']):

                    if comp_id not in EMPTY_VALUE:
                        if comp_id not in STD_MON_DICT:
                            continue
                        if not self.__csStat.peptideLike(comp_id):
                            continue
                        rci_residues.append([comp_id, auth_seq_id])
                    else:
                        continue

                    _cs_data = cs_data_by_residue.get((auth_chain_id, auth_seq_id, comp_id), {})

                    if len(_cs_data) == 0:
                        continue

                    has_bb_atoms = False

                    for cs_key, cs_vals in _cs_data.items():
                        atom_id = cs_key[3]

                        if atom_id not in rci_atom_ids:
                            continue

                        rci_assignments.append([comp_id, auth_seq_id, atom_id, atom_id[0], cs_vals['value']])

                        has_bb_atoms = True

                    if has_bb_atoms:

                        if comp_id in ('CYS', 'DCY'):

                            ca_chem_shift = cb_chem_shift = None

                            for cs_key, cs_vals in _cs_data.items():
                                atom_id = cs_key[3]

                                if atom_id == 'CA':
                                    ca_chem_shift = cs_vals['value']
                                elif atom_id == 'CB':
                                    cb_chem_shift = cs_vals['value']

                            ambig_redox_state = False

                            if cb_chem_shift is not None:
                                if cb_chem_shift < 32.0:
                                    pass
                                elif cb_chem_shift > 35.0:
                                    oxidized_cys_seq_ids.append(auth_seq_id)
                                else:
                                    ambig_redox_state = True
                            elif ca_chem_shift is not None:
                                ambig_redox_state = True

                            if ambig_redox_state:
                                oxi, red = predict_redox_state_of_cysteine(ca_chem_shift, cb_chem_shift)
                                if oxi < 0.001:
                                    pass
                                elif red < 0.001 or oxi > 0.5:
                                    oxidized_cys_seq_ids.append(auth_seq_id)

                    else:
                        seq_ids_wo_assign.append(auth_seq_id)

                if self.__chemShiftUnmapped is not None and list_id in self.__chemShiftUnmapped:
                    unmap_cs_data = self.__chemShiftUnmapped[list_id]

                    unmap_cs_data_ = [cs for cs in unmap_cs_data
                                      if cs['auth_chain_id'] == auth_chain_id
                                      and cs['value'] not in EMPTY_VALUE]

                    if len(unmap_cs_data_) > 0:
                        unmap_rci_residues, unmap_rci_residue_set = [], set()

                        unmap_cs_data_by_residue = {}
                        for cs in unmap_cs_data_:
                            residue_key = (cs['auth_seq_id'], cs['comp_id'])
                            if residue_key not in unmap_cs_data_by_residue:
                                unmap_cs_data_by_residue[residue_key] = []
                            unmap_cs_data_by_residue[residue_key].append(cs)

                        for cs in unmap_cs_data_:
                            auth_seq_id = int(cs['auth_seq_id'])
                            comp_id = cs['comp_id']
                            atom_id = cs['atom_id']

                            if comp_id not in EMPTY_VALUE:
                                if comp_id not in STD_MON_DICT:
                                    continue
                                if not self.__csStat.peptideLike(comp_id):
                                    continue

                                if (comp_id, auth_seq_id) not in unmap_rci_residue_set:
                                    unmap_rci_residue_set.add((comp_id, auth_seq_id))
                                    unmap_rci_residues.append([comp_id, auth_seq_id])

                        if len(unmap_rci_residues) > 0:
                            rci_residues.extend(unmap_rci_residues)
                            rci_residues = sorted(rci_residues, key=itemgetter(1))

                            for comp_id, auth_seq_id in unmap_rci_residues:

                                _unmap_cs_data = unmap_cs_data_by_residue.get((auth_seq_id, comp_id), [])

                                if len(_unmap_cs_data) == 0:
                                    continue

                                has_bb_atoms = False

                                for cs in _unmap_cs_data:
                                    atom_id = cs['atom_id']

                                    if atom_id not in rci_atom_ids:
                                        continue

                                    rci_assignments.append([comp_id, auth_seq_id, atom_id, atom_id[0], cs['value']])

                                    has_bb_atoms = True

                                if has_bb_atoms:

                                    if comp_id in ('CYS', 'DCY'):

                                        ca_chem_shift = cb_chem_shift = None

                                        for cs in _unmap_cs_data:
                                            atom_id = cs['atom_id']

                                            if atom_id == 'CA':
                                                ca_chem_shift = cs['value']
                                            elif atom_id == 'CB':
                                                cb_chem_shift = cs['value']

                                        ambig_redox_state = False

                                        if cb_chem_shift is not None:
                                            if cb_chem_shift < 32.0:
                                                pass
                                            elif cb_chem_shift > 35.0:
                                                oxidized_cys_seq_ids.append(auth_seq_id)
                                            else:
                                                ambig_redox_state = True
                                        elif ca_chem_shift is not None:
                                            ambig_redox_state = True

                                        if ambig_redox_state:
                                            oxi, red = predict_redox_state_of_cysteine(ca_chem_shift, cb_chem_shift)
                                            if oxi < 0.001:
                                                pass
                                            elif red < 0.001 or oxi > 0.5:
                                                oxidized_cys_seq_ids.append(auth_seq_id)

                                else:
                                    seq_ids_wo_assign.append(auth_seq_id)

                if len(rci_assignments) > 0:
                    rci_result[list_id][auth_chain_id] =\
                        rci.calculate(rci_residues, rci_assignments, oxidized_cys_seq_ids, seq_ids_wo_assign)

                    has_rci_results = True

        self.__results['rci'] = rci_result

        if has_rci_results:
            self.__results['rci_version'] = rci.version

        self.__chemShiftDict = None

        return True

    def __summarizeCommonMrAnalysis(self) -> bool:
        """ Summarize common restraint analysis.
        """

        if self.__has_prev_results:
            return True

        self.__results = {'total_models': self.__total_models, 'eff_model_ids': self.__eff_model_ids,
                          'atom_ids': self.__atomIdList, 'key_lists': {}}

        for datablock_name in self.__rR.getDataBlockNameList():

            if not self.__rR.hasCategory('Chem_comp_assembly', datablock_name):
                continue

            self.__results['seq_length'] = self.__rR.getRowLength('Chem_comp_assembly', datablock_name)

            break

        return True

    def __summarizeDistanceRestraintAnalysis(self) -> bool:
        """ Summarize distance restraint analysis results.
            @author: Masashi Yokochi
            @note: Derived from wwpdb.apps.validation.src.RestraintsValidation.BMRBRestraintsAnalysis.generate_output,
                   written by Kumaran Baskaran
            @change: class method, improve readability of restraints, support combinational restraints, performance optimization
        """

        if self.__has_prev_results:
            return True

        if self.__distRestDict is None or len(self.__distRestDict) == 0:
            self.__results['distance'] = False
            return True

        try:

            self.__results['distance'] = self.__distRestViolDict is not None and len(self.__distRestViolDict) > 0

            self.__results['dist_seq_dict'] = self.__distRestSeqDict
            self.__results['unmapped_dist'] = self.__distRestUnmapped

            any_type = 'total'

            distance_type = ('intraresidue', 'sequential', 'medium', 'long', 'interchain', any_type)
            distance_sub_type = ('backbone-backbone', 'backbone-sidechain', 'sidechain-sidechain')
            bond_flag = ('hbond', 'sbond', 'sebond', 'metal', None)

            self.__results['key_lists']['distance_type'] = distance_type
            self.__results['key_lists']['distance_sub_type'] = distance_sub_type
            self.__results['key_lists']['bond_flag'] = bond_flag

            distance_summary, distance_violation = {}, {}

            consistent_distance_violation, distance_violations_vs_models, distance_violations_in_models =\
                {}, {}, {}

            # Non-per-model accumulators: build once (were re-created and overwritten
            # inside the per-model loop below, i.e. total_models times).
            for t in distance_type:
                distance_summary[t], distance_violation[t], \
                    consistent_distance_violation[t], distance_violations_vs_models[t] = {}, {}, {}, {}

                for s in distance_sub_type:
                    distance_summary[t][s], distance_violation[t][s], \
                        consistent_distance_violation[t][s], distance_violations_vs_models[t][s] = {}, {}, {}, {}

                    for b in bond_flag:
                        distance_summary[t][s][b] = distance_violation[t][s][b] =\
                            consistent_distance_violation[t][s][b] = 0
                        distance_violations_vs_models[t][s][b] = [0] * (self.__total_models + 1)

            for m in self.__eff_model_ids:
                distance_violations_in_models[m] = {}

                for t in distance_type:
                    distance_violations_in_models[m][t] = {}

                    for s in distance_sub_type:
                        distance_violations_in_models[m][t][s] = {}

                        for b in bond_flag:
                            distance_violations_in_models[m][t][s][b] = []

            for rest_key, restraints in self.__distRestDict.items():
                r = restraints[0]

                t = r['distance_type']
                s = r['distance_sub_type']
                b = r['bond_flag']

                distance_summary[t][s][b] += 1
                distance_summary[any_type][s][b] += 1

                len_vm = len(get_violated_model_ids(self.__distRestViolDict[rest_key]))

                if len_vm > 0:
                    distance_violation[t][s][b] += 1
                    distance_violation[any_type][s][b] += 1

                if len_vm == self.__total_models:
                    consistent_distance_violation[t][s][b] += 1
                    consistent_distance_violation[any_type][s][b] += 1

                distance_violations_vs_models[t][s][b][len_vm] += 1
                distance_violations_vs_models[any_type][s][b][len_vm] += 1

            self.__results['distance_summary'] = distance_summary
            self.__results['distance_violation'] = distance_violation
            self.__results['consistent_distance_violation'] = consistent_distance_violation
            self.__results['distance_violations_vs_models'] = distance_violations_vs_models

            for rest_key, viol_per_model in self.__distRestViolDict.items():
                for m in self.__eff_model_ids:
                    err = viol_per_model[m]

                    if err is None or err == 0.0:
                        continue

                    comb_key = self.__distRestViolCombKeyDict[rest_key][m]

                    if comb_key is None:
                        continue

                    r = self.__distRestDictWithCombKey[rest_key][comb_key][0]
                    t = r['distance_type']
                    s = r['distance_sub_type']
                    b = r['bond_flag']

                    distance_violations_in_models[m][t][s][b].append(err)
                    distance_violations_in_models[m][any_type][s][b].append(err)

            self.__results['distance_violations_in_models'] = distance_violations_in_models

            dist_range = []

            residual_distance_violation = {}
            for idx, beg_err_bin in enumerate(NMR_VTF_DIST_ERR_BINS):
                if idx == len(NMR_VTF_DIST_ERR_BINS) - 1:
                    end_err_bin = None
                    dist_err_range = f">{beg_err_bin}"
                else:
                    end_err_bin = NMR_VTF_DIST_ERR_BINS[idx + 1]
                    dist_err_range = f"{beg_err_bin}-{end_err_bin}"

                dist_range.append(dist_err_range)
                residual_distance_violation[dist_err_range] =\
                    get_violation_statistics_for_each_bin(beg_err_bin, end_err_bin,
                                                          self.__total_models,
                                                          self.__eff_model_ids,
                                                          self.__distRestViolDict)

            self.__results['key_lists']['dist_range'] = dist_range
            self.__results['residual_distance_violation'] = residual_distance_violation

            most_violated_distance = []
            for rest_key, viol_per_model in self.__distRestViolDict.items():
                vm = get_violated_model_ids(viol_per_model)

                if len(vm) > 1:
                    e = numpy.array([err for err in viol_per_model.values() if err is not None and err > 0.0], dtype=float)
                    # e is fixed per rest_key; compute the five statistics once
                    # (were recomputed for every r below).
                    e_min, e_max, e_mean, e_std, e_median = \
                        numpy.min(e), numpy.max(e), numpy.mean(e), numpy.std(e), numpy.median(e)

                    comb_keys = []
                    for _m in set(vm):
                        comb_key = self.__distRestViolCombKeyDict[rest_key][_m]

                        if comb_key is None:
                            continue

                        if comb_key not in comb_keys:
                            comb_keys.append(comb_key)

                    for comb_key in comb_keys:
                        for r in self.__distRestDictWithCombKey[rest_key][comb_key]:
                            most_violated_distance.append([rest_key,
                                                           r['atom_key_1'],
                                                           r['atom_key_2'],
                                                           r['distance_type'],
                                                           r['distance_sub_type'],
                                                           r['bond_flag'],
                                                           len(vm),
                                                           vm,
                                                           e_min,
                                                           e_max,
                                                           e_mean,
                                                           e_std,
                                                           e_median])

            self.__results['most_violated_distance'] =\
                sorted(most_violated_distance, reverse=True, key=itemgetter(6, 10))

            all_distance_violations = []
            for rest_key, viol_per_model in self.__distRestViolDict.items():
                for m in self.__eff_model_ids:
                    err = viol_per_model[m]

                    if err is None or err == 0.0:
                        continue

                    comb_key = self.__distRestViolCombKeyDict[rest_key][m]

                    if comb_key is None:
                        continue

                    for r in self.__distRestDictWithCombKey[rest_key][comb_key]:
                        all_distance_violations.append([rest_key,
                                                        r['atom_key_1'],
                                                        r['atom_key_2'],
                                                        m,
                                                        r['distance_type'],
                                                        r['distance_sub_type'],
                                                        r['bond_flag'],
                                                        err])

            self.__results['all_distance_violations'] =\
                sorted(all_distance_violations, reverse=True, key=itemgetter(7, 0))

            dist_violation_seq = {}
            for seq_key, rest_keys in self.__distRestSeqDict.items():
                for m in self.__eff_model_ids:
                    _seq_key = (seq_key[0], seq_key[1], seq_key[2], m)

                    if _seq_key not in dist_violation_seq:
                        dist_violation_seq[_seq_key] = []

                    for rest_key in rest_keys:
                        err = self.__distRestViolDict[rest_key][m]

                        if err is None or err == 0.0:
                            continue

                        comb_key = self.__distRestViolCombKeyDict[rest_key][m]

                        if comb_key is None:
                            continue

                        atom_ids = set()
                        distance_type = distance_sub_type = bond_flag = None

                        for r in self.__distRestDictWithCombKey[rest_key][comb_key]:
                            seq_key_1 = (r['atom_key_1'][0], r['atom_key_1'][1], r['atom_key_1'][2])
                            seq_key_2 = (r['atom_key_2'][0], r['atom_key_2'][1], r['atom_key_2'][2])

                            if seq_key not in (seq_key_1, seq_key_2):
                                if self.__verbose and comb_key[0] is None and comb_key[1] is None:
                                    self.__log.write(f"Nothing matches with sequence, {seq_key_1}, {seq_key_2}, {seq_key}\n")
                                continue

                            if seq_key_1 == seq_key_2:
                                atom_ids.add(r['atom_key_1'][3])
                                atom_ids.add(r['atom_key_2'][3])
                            elif seq_key_1 == seq_key:
                                atom_ids.add(r['atom_key_1'][3])
                            else:
                                atom_ids.add(r['atom_key_2'][3])

                            if distance_type is None:
                                distance_type = r['distance_type']
                                distance_sub_type = r['distance_sub_type']
                                bond_flag = r['bond_flag']

                        dist_violation_seq[_seq_key].append([rest_key[0],
                                                             rest_key[1],
                                                             sorted(list(atom_ids)),
                                                             distance_type,
                                                             distance_sub_type,
                                                             bond_flag,
                                                             err])

            self.__results['dist_violation_seq'] = dist_violation_seq

            self.__distRestDict = None

            return True

        except Exception as e:  # pylint: disable=broad-exception-caught
            self.__log.write(f"Exception occurred while processing {os.path.basename(self.__cifPath)} "
                             f"and {os.path.basename(self.__nmrDataPath)}\n")
            self.__log.write(f"+{self.__class_name__}.__summarizeDistanceRestraintAnalysis() ++ Error  - {str(e)}\n")

        return False

    def __summarizeDihedralAngleRestraintAnalysis(self) -> bool:
        """ Summarize dihedral angle restraint analysis results.
            @author: Masashi Yokochi
            @note: Derived from wwpdb.apps.validation.src.RestraintsValidation.BMRBRestraintsAnalysis.generate_output,
                   written by Kumaran Baskaran
            @change: class method, improve readability of restraints, support combinational restraints, performance optimization
        """

        if self.__has_prev_results:
            return True

        if self.__dihedRestDict is None or len(self.__dihedRestDict) == 0:
            self.__results['angle'] = False
            self.__results['error_message_angle'] = None
            return True

        try:

            self.__results['angle'] = self.__dihedRestViolDict is not None and len(self.__dihedRestViolDict) > 0
            self.__results['error_message_angle'] = None

            self.__results['angle_seq_dict'] = self.__dihedRestSeqDict
            self.__results['unmapped_angle'] = self.__dihedRestUnmapped

            any_type = 'Total'

            try:
                angle_type_set = set()
                for r_list in self.__dihedRestDict.values():
                    for r in r_list:
                        angle_type_set.add(r['angle_type'])
                angle_type = sorted(list(angle_type_set)) + [any_type]
            except IndexError:
                self.__log.write("Dihedral angle analysis failed due to data error in the dihedral angle restraints. "
                                 f"{self.__dihedRestDict.values()}\n")
                self.__results['error_message_angle'] =\
                    'Dihedral angle analysis failed due to data error in the dihedral angle restraints, '\
                    'possibly missing target value'
                self.__results['angle'] = False
                return True

            self.__results['key_lists']['angle_type'] = angle_type

            angle_summary, angle_violation = {}, {}

            consistent_angle_violation, angle_violations_vs_models, angle_violations_in_models =\
                {}, {}, {}

            # Non-per-model accumulators: build once (were re-created and overwritten
            # inside the per-model loop below, i.e. total_models times).
            for t in angle_type:
                angle_summary[t] = 0
                angle_violation[t] = 0
                consistent_angle_violation[t] = 0
                angle_violations_vs_models[t] = [0] * (self.__total_models + 1)

            for m in self.__eff_model_ids:
                angle_violations_in_models[m] = {}

                for t in angle_type:
                    angle_violations_in_models[m][t] = []

            for rest_key, restraints in self.__dihedRestDict.items():
                r = restraints[0]

                t = r['angle_type']

                angle_summary[t] += 1
                angle_summary[any_type] += 1

                len_vm = len(get_violated_model_ids(self.__dihedRestViolDict[rest_key]))

                if len_vm > 0:
                    angle_violation[t] += 1
                    angle_violation[any_type] += 1

                if len_vm == self.__total_models:
                    consistent_angle_violation[t] += 1
                    consistent_angle_violation[any_type] += 1

                angle_violations_vs_models[t][len_vm] += 1
                angle_violations_vs_models[any_type][len_vm] += 1

            self.__results['angle_summary'] = angle_summary
            self.__results['angle_violation'] = angle_violation
            self.__results['consistent_angle_violation'] = consistent_angle_violation
            self.__results['angle_violations_vs_models'] = angle_violations_vs_models

            for rest_key, viol_per_model in self.__dihedRestViolDict.items():
                for m in self.__eff_model_ids:
                    err = viol_per_model[m]

                    if err is None or err == 0.0:
                        continue

                    comb_key = self.__dihedRestViolCombKeyDict[rest_key][m]

                    if comb_key is None:
                        continue

                    r = self.__dihedRestDictWithCombKey[rest_key][comb_key][0]

                    t = r['angle_type']

                    angle_violations_in_models[m][t].append(err)
                    angle_violations_in_models[m][any_type].append(err)

            self.__results['angle_violations_in_models'] = angle_violations_in_models

            angle_range = []

            residual_angle_violation = {}
            for idx, beg_err_bin in enumerate(NMR_VTF_DIHED_ERR_BINS):
                if idx == len(NMR_VTF_DIHED_ERR_BINS) - 1:
                    end_err_bin = None
                    dihed_err_range = f">{beg_err_bin}"
                else:
                    end_err_bin = NMR_VTF_DIHED_ERR_BINS[idx + 1]
                    dihed_err_range = f"{beg_err_bin}-{end_err_bin}"

                angle_range.append(dihed_err_range)
                residual_angle_violation[dihed_err_range] =\
                    get_violation_statistics_for_each_bin(beg_err_bin, end_err_bin,
                                                          self.__total_models,
                                                          self.__eff_model_ids,
                                                          self.__dihedRestViolDict)

            self.__results['key_lists']['angle_range'] = angle_range
            self.__results['residual_angle_violation'] = residual_angle_violation

            most_violated_angle = []
            for rest_key, viol_per_model in self.__dihedRestViolDict.items():
                vm = get_violated_model_ids(viol_per_model)

                if len(vm) > 1:
                    e = numpy.array([err for err in viol_per_model.values() if err is not None and err > 0.0], dtype=float)
                    e_min, e_max, e_mean, e_std, e_median = \
                        numpy.min(e), numpy.max(e), numpy.mean(e), numpy.std(e), numpy.median(e)

                    comb_keys = []
                    for _m in set(vm):
                        comb_key = self.__dihedRestViolCombKeyDict[rest_key][_m]

                        if comb_key is None:
                            continue

                        if comb_key not in comb_keys:
                            comb_keys.append(comb_key)

                    for comb_key in comb_keys:
                        for r in self.__dihedRestDictWithCombKey[rest_key][comb_key]:
                            most_violated_angle.append([rest_key,
                                                        r['atom_key_1'],
                                                        r['atom_key_2'],
                                                        r['atom_key_3'],
                                                        r['atom_key_4'],
                                                        r['angle_type'],
                                                        len(vm),
                                                        vm,
                                                        e_min,
                                                        e_max,
                                                        e_mean,
                                                        e_std,
                                                        e_median])

            self.__results['most_violated_angle'] =\
                sorted(most_violated_angle, reverse=True, key=itemgetter(6, 10))

            all_angle_violations = []
            for rest_key, viol_per_model in self.__dihedRestViolDict.items():
                for m in self.__eff_model_ids:
                    err = viol_per_model[m]

                    if err is None or err == 0.0:
                        continue

                    comb_key = self.__dihedRestViolCombKeyDict[rest_key][m]

                    if comb_key is None:
                        continue

                    for r in self.__dihedRestDictWithCombKey[rest_key][comb_key]:
                        all_angle_violations.append([rest_key,
                                                    r['atom_key_1'],
                                                    r['atom_key_2'],
                                                    r['atom_key_3'],
                                                    r['atom_key_4'],
                                                    m,
                                                    r['angle_type'],
                                                    err])

            self.__results['all_angle_violations'] =\
                sorted(all_angle_violations, reverse=True, key=itemgetter(7, 0))

            angle_violation_seq = {}
            for seq_key, rest_keys in self.__dihedRestSeqDict.items():
                for m in self.__eff_model_ids:
                    _seq_key = (seq_key[0], seq_key[1], seq_key[2], m)

                    if _seq_key not in angle_violation_seq:
                        angle_violation_seq[_seq_key] = []

                    for rest_key in rest_keys:
                        err = self.__dihedRestViolDict[rest_key][m]

                        if err is None or err == 0.0:
                            continue

                        comb_key = self.__dihedRestViolCombKeyDict[rest_key][m]

                        if comb_key is None:
                            continue

                        atom_ids_1, atom_ids_2, atom_ids_3, atom_ids_4 = [], [], [], []
                        angle_type = None

                        for r in self.__dihedRestDictWithCombKey[rest_key][comb_key]:
                            seq_key_1 = (r['atom_key_1'][0], r['atom_key_1'][1], r['atom_key_1'][2])
                            if seq_key_1 == seq_key:
                                atom_ids_1.append(r['atom_key_1'][3])
                            seq_key_2 = (r['atom_key_2'][0], r['atom_key_2'][1], r['atom_key_2'][2])
                            if seq_key_2 == seq_key:
                                atom_ids_2.append(r['atom_key_2'][3])
                            seq_key_3 = (r['atom_key_3'][0], r['atom_key_3'][1], r['atom_key_3'][2])
                            if seq_key_3 == seq_key:
                                atom_ids_3.append(r['atom_key_3'][3])
                            seq_key_4 = (r['atom_key_4'][0], r['atom_key_4'][1], r['atom_key_4'][2])
                            if seq_key_4 == seq_key:
                                atom_ids_4.append(r['atom_key_4'][3])

                            if angle_type is None:
                                angle_type = r['angle_type']

                        atom_ids = list(dict.fromkeys(atom_ids_1))
                        atom_ids.extend(list(dict.fromkeys(atom_ids_2)))
                        atom_ids.extend(list(dict.fromkeys(atom_ids_3)))
                        atom_ids.extend(list(dict.fromkeys(atom_ids_4)))

                        angle_violation_seq[_seq_key].append([rest_key[0],
                                                              rest_key[1],
                                                              atom_ids,
                                                              angle_type,
                                                              err])

            self.__results['angle_violation_seq'] = angle_violation_seq

            self.__dihedRestDict = None

            return True

        except Exception as e:  # pylint: disable=broad-exception-caught
            self.__log.write(f"Exception occurred while processing {os.path.basename(self.__cifPath)} "
                             f"and {os.path.basename(self.__nmrDataPath)}\n")
            self.__log.write(f"+{self.__class_name__}.__summarizeDihedralAngleRestraintAnalysis() ++ Error  - {str(e)}\n")

        return False

    def __summarizeRdcRestraintAnalysis(self) -> bool:
        """ Summarize RDC restraint analysis results.
            @author: Masashi Yokochi
        """

        if self.__has_prev_results:
            return True

        if self.__rdcRestDict is None or len(self.__rdcRestDict) == 0:
            self.__results['rdc'] = False
            self.__results['error_message_rdc'] = None
            return True

        try:

            self.__results['rdc'] = self.__rdcRestViolDict is not None and len(self.__rdcRestViolDict) > 0
            self.__results['error_message_rdc'] = None

            self.__results['rdc_seq_dict'] = self.__rdcRestSeqDict
            self.__results['unmapped_rdc'] = self.__rdcRestUnmapped

            any_type = 'Total'

            try:
                rdc_type_set = set()
                for r_list in self.__rdcRestDict.values():
                    for r in r_list:
                        rdc_type_set.add(r['rdc_type'])
                rdc_type = sorted(list(rdc_type_set)) + [any_type]
            except IndexError:
                self.__log.write(f"RDC analysis failed due to data error in the RDC restraints. {self.__rdcRestDict.values()}\n")
                self.__results['error_message_rdc'] = 'RDC analysis failed due to data error in the RDC angle restraints'
                self.__results['rdc'] = False
                return True

            self.__results['key_lists']['rdc_type'] = rdc_type

            rdc_summary, rdc_violation = {}, {}

            consistent_rdc_violation, rdc_violations_vs_models, rdc_violations_in_models =\
                {}, {}, {}

            # Non-per-model accumulators: build once (were re-created and overwritten
            # inside the per-model loop below, i.e. total_models times).
            for t in rdc_type:
                rdc_summary[t] = 0
                rdc_violation[t] = 0
                consistent_rdc_violation[t] = 0
                rdc_violations_vs_models[t] = [0] * (self.__total_models + 1)

            for m in self.__eff_model_ids:
                rdc_violations_in_models[m] = {}

                for t in rdc_type:
                    rdc_violations_in_models[m][t] = []

            for rest_key, restraints in self.__rdcRestDict.items():
                r = restraints[0]

                t = r['rdc_type']

                rdc_summary[t] += 1
                rdc_summary[any_type] += 1

                len_vm = len(get_violated_model_ids(self.__rdcRestViolDict[rest_key]))

                if len_vm > 0:
                    rdc_violation[t] += 1
                    rdc_violation[any_type] += 1

                if len_vm == self.__total_models:
                    consistent_rdc_violation[t] += 1
                    consistent_rdc_violation[any_type] += 1

                rdc_violations_vs_models[t][len_vm] += 1
                rdc_violations_vs_models[any_type][len_vm] += 1

            self.__results['rdc_summary'] = rdc_summary
            self.__results['rdc_violation'] = rdc_violation
            self.__results['consistent_rdc_violation'] = consistent_rdc_violation
            self.__results['rdc_violations_vs_models'] = rdc_violations_vs_models

            for rest_key, viol_per_model in self.__rdcRestViolDict.items():
                for m in self.__eff_model_ids:
                    err = viol_per_model[m]

                    if err is None or err == 0.0:
                        continue

                    comb_key = self.__rdcRestViolCombKeyDict[rest_key][m]

                    if comb_key is None:
                        continue

                    r = self.__rdcRestDictWithCombKey[rest_key][comb_key][0]

                    t = r['rdc_type']

                    rdc_violations_in_models[m][t].append(err)
                    rdc_violations_in_models[m][any_type].append(err)

            self.__results['rdc_violations_in_models'] = rdc_violations_in_models

            rdc_range = []

            residual_rdc_violation = {}
            for idx, beg_err_bin in enumerate(NMR_VTF_RDC_ERR_BINS):
                if idx == len(NMR_VTF_RDC_ERR_BINS) - 1:
                    end_err_bin = None
                    rdc_err_range = f">{beg_err_bin}"
                else:
                    end_err_bin = NMR_VTF_RDC_ERR_BINS[idx + 1]
                    rdc_err_range = f"{beg_err_bin}-{end_err_bin}"

                rdc_range.append(rdc_err_range)
                residual_rdc_violation[rdc_err_range] =\
                    get_violation_statistics_for_each_bin(beg_err_bin, end_err_bin,
                                                          self.__total_models,
                                                          self.__eff_model_ids,
                                                          self.__rdcRestViolDict)

            self.__results['key_lists']['rdc_range'] = rdc_range
            self.__results['residual_rdc_violation'] = residual_rdc_violation

            most_violated_rdc = []
            for rest_key, viol_per_model in self.__rdcRestViolDict.items():
                vm = get_violated_model_ids(viol_per_model)

                if len(vm) > 1:
                    e = numpy.array([err for err in viol_per_model.values() if err is not None and err > 0.0], dtype=float)
                    e_min, e_max, e_mean, e_std, e_median = \
                        numpy.min(e), numpy.max(e), numpy.mean(e), numpy.std(e), numpy.median(e)

                    comb_keys = []
                    for _m in set(vm):
                        comb_key = self.__rdcRestViolCombKeyDict[rest_key][_m]

                        if comb_key is None:
                            continue

                        if comb_key not in comb_keys:
                            comb_keys.append(comb_key)

                    for comb_key in comb_keys:
                        for r in self.__rdcRestDictWithCombKey[rest_key][comb_key]:
                            most_violated_rdc.append([rest_key,
                                                      r['atom_key_1'],
                                                      r['atom_key_2'],
                                                      r['rdc_type'],
                                                      len(vm),
                                                      vm,
                                                      e_min,
                                                      e_max,
                                                      e_mean,
                                                      e_std,
                                                      e_median])

            self.__results['most_violated_rdc'] =\
                sorted(most_violated_rdc, reverse=True, key=itemgetter(4, 8))

            all_rdc_violations = []
            for rest_key, viol_per_model in self.__rdcRestViolDict.items():
                for m in self.__eff_model_ids:
                    err = viol_per_model[m]

                    if err is None or err == 0.0:
                        continue

                    comb_key = self.__rdcRestViolCombKeyDict[rest_key][m]

                    if comb_key is None:
                        continue

                    for r in self.__rdcRestDictWithCombKey[rest_key][comb_key]:
                        all_rdc_violations.append([rest_key,
                                                  r['atom_key_1'],
                                                  r['atom_key_2'],
                                                  m,
                                                  r['rdc_type'],
                                                  err])

            self.__results['all_rdc_violations'] =\
                sorted(all_rdc_violations, reverse=True, key=itemgetter(5, 0))

            rdc_violation_seq = {}
            for seq_key, rest_keys in self.__rdcRestSeqDict.items():
                for m in self.__eff_model_ids:
                    _seq_key = (seq_key[0], seq_key[1], seq_key[2], m)

                    if _seq_key not in rdc_violation_seq:
                        rdc_violation_seq[_seq_key] = []

                    for rest_key in rest_keys:
                        err = self.__rdcRestViolDict[rest_key][m]

                        if err is None or err == 0.0:
                            continue

                        comb_key = self.__rdcRestViolCombKeyDict[rest_key][m]

                        if comb_key is None:
                            continue

                        atom_ids_1, atom_ids_2 = [], []
                        rdc_type = None

                        for r in self.__rdcRestDictWithCombKey[rest_key][comb_key]:
                            seq_key_1 = (r['atom_key_1'][0], r['atom_key_1'][1], r['atom_key_1'][2])
                            if seq_key_1 == seq_key:
                                atom_ids_1.append(r['atom_key_1'][3])
                            seq_key_2 = (r['atom_key_2'][0], r['atom_key_2'][1], r['atom_key_2'][2])
                            if seq_key_2 == seq_key:
                                atom_ids_2.append(r['atom_key_2'][3])

                            if rdc_type is None:
                                rdc_type = r['rdc_type']

                        atom_ids = list(dict.fromkeys(atom_ids_1))
                        atom_ids.extend(list(dict.fromkeys(atom_ids_2)))

                        rdc_violation_seq[_seq_key].append([rest_key[0],
                                                            rest_key[1],
                                                            atom_ids,
                                                            rdc_type,
                                                            err])

            self.__results['rdc_violation_seq'] = rdc_violation_seq

            if len(self.__rdcCorrPlotDict) > 0:
                self.__results['rdc_correlation_plot'] = self.__rdcCorrPlotDict

            self.__rdcRestDict = None

            return True

        except Exception as e:  # pylint: disable=broad-exception-caught
            self.__log.write(f"Exception occurred while processing {os.path.basename(self.__cifPath)} "
                             f"and {os.path.basename(self.__nmrDataPath)}\n")
            self.__log.write(f"+{self.__class_name__}.__summarizeRdcRestraintAnalysis() ++ Error  - {str(e)}\n")

        return False

    def __outputResultsAsPickleFile(self) -> bool:
        """ Output results if 'result_pickle_file_path' was set in output parameter.
        """

        if self.__use_cache:

            cache_path = None
            if self.__resultsCacheName is not None:
                cache_path = os.path.join(self.__cacheDirPath, self.__resultsCacheName)

                if not os.path.exists(cache_path):
                    write_as_pickle(self.__results, cache_path)

            if RESULT_PKL_FILE_PATH_KEY in self.__outputParamDict:
                pickle_file = self.__outputParamDict[RESULT_PKL_FILE_PATH_KEY]

                if cache_path is None:
                    write_as_pickle(self.__results, pickle_file)
                    return True

                try:

                    if os.path.exists(pickle_file):
                        os.remove(pickle_file)

                    os.symlink(cache_path, pickle_file)

                except OSError:
                    pass

        elif RESULT_PKL_FILE_PATH_KEY in self.__outputParamDict:
            pickle_file = self.__outputParamDict[RESULT_PKL_FILE_PATH_KEY]

            write_as_pickle(self.__results, pickle_file)

        return True
