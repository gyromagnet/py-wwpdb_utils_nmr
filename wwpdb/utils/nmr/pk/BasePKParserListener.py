##
# File: BasePKParserListener.py
# Date: 03-Dec-2024
#
# Updates:
# 09-Oct-2026  M. Yokochi - assignCoordPolymerSequenceWithChainIdWithoutCompId() without a chain returns no assignment again,
#                           so that extractPeakAssignment() falls back to ...WithoutCompId() without pinning the default segment id
#                           (AttributeError on the bare base instance of XEASY PROT) (DAOTHER-7829)
# 09-Oct-2026  M. Yokochi - keep the chain id in assignCoordPolymerSequence{WithChainId,}WithoutCompId() instead of
#                           overwriting it with in-loop sequence remap lookups, which recorded the assignment under
#                           chain None or spread it to every chain with that residue number; the chain-assignment index
#                           of extractPeakAssignment() no longer clobbers the term index of its loop; J-coupling transfer
#                           remediation converts an integer chain id as relayed transfer remediation does (DAOTHER-7829)
""" ParserLister base class for any peak list file.
    @author: Masashi Yokochi
"""
__docformat__ = "restructuredtext en"
__author__ = "Masashi Yokochi"
__email__ = "yokochi@protein.osaka-u.ac.jp"
__license__ = "Apache License 2.0"
__version__ = "1.2.0"

import collections
import copy
import functools
import itertools
import re
import sys
from typing import IO, List, Optional, Tuple, Union

import numpy

import pynmrstar

try:
    from wwpdb.utils.nmr.NmrDpConstant import (LARGE_ASYM_ID,
                                               EMPTY_VALUE,
                                               ELEMENT_SYMBOLS,
                                               STD_MON_DICT,
                                               PROTON_BEGIN_CODE,
                                               PSE_PRO_BEGIN_CODE,
                                               AMINO_PROTON_CODE,
                                               CARBOXYL_CODE,
                                               ZINC_ION_CODE,
                                               CALCIUM_ION_CODE,
                                               ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS,
                                               REPRESENTATIVE_MODEL_ID,
                                               REPRESENTATIVE_ALT_ID,
                                               MAX_PREF_LABEL_SCHEME_COUNT,
                                               MAX_ALLOWED_EXT_SEQ,
                                               UNREAL_AUTH_SEQ_NUM,
                                               MAX_CONFLICT_ATTEMPT,
                                               CS_RESTRAINT_RANGE,
                                               CS_RESTRAINT_ERROR,
                                               WEIGHT_RANGE,
                                               CS_RANGE_MIN,
                                               CS_RANGE_MAX,
                                               CS_ERROR_MIN,
                                               CS_ERROR_MAX,
                                               WEIGHT_RANGE_MIN,
                                               WEIGHT_RANGE_MAX,
                                               HEME_LIKE_RES_NAMES,
                                               SPECTRAL_DIM_TEMPLATE,
                                               ASSIGNMENT_SEPARATOR_PAT,
                                               ASSIGNMENT_RESID_PAT,
                                               ASSIGNMENT_HALF_SPIN_NUCLEUS,
                                               C_CARBONYL_CENTER_MAX_TOR,
                                               C_CARBONYL_CENTER_MAX,
                                               C_CARBONYL_CENTER_MIN,
                                               C_AROMATIC_CENTER_MAX,
                                               C_AROMATIC_CENTER_MIN,
                                               C_AROMATIC_CENTER_MIN_TOR,
                                               N_AMIDE_CENTER_MAX,
                                               N_AMIDE_CENTER_MIN,
                                               C_ALL_CENTER_MAX,
                                               C_ALL_CENTER_MIN,
                                               C_ALIPHATIC_CENTER_MAX,
                                               C_ALIPHATIC_CENTER_MIN,
                                               C_METHYL_CENTER_MAX,
                                               C_METHYL_CENTER_MIN,
                                               HN_AROMATIC_CENTER_MAX,
                                               HN_AROMATIC_CENTER_MIN,
                                               H_IMIDE_CENTER_MAX,
                                               H_IMIDE_CENTER_MIN,
                                               H_ALL_CENTER_MAX,
                                               H_ALL_CENTER_MIN,
                                               H_ALIPHATIC_CENTER_MAX,
                                               H_ALIPHATIC_CENTER_MIN,
                                               H_METHYL_CENTER_MAX,
                                               H_METHYL_CENTER_MIN,
                                               MIN_CORRCOEF_FOR_ONE_BOND_TRANSFER,
                                               INSTRUCTION_FOR_FULL_SEQUENCE)
    from wwpdb.utils.nmr.AlignUtil import (getOneLetterCode,
                                           updatePolySeqRst,
                                           revertPolySeqRst,
                                           sortPolySeqRst,
                                           syncCompIdOfPolySeqRst,
                                           alignPolymerSequence,
                                           alignPolymerSequenceWithConflicts,
                                           assignPolymerSequence,
                                           trimSequenceAlignment,
                                           updatePolySeqRstAmbig,
                                           retrieveAtomIdentFromMRMap,
                                           retrieveAtomIdFromMRMap,
                                           retrieveRemappedSeqId,
                                           retrieveRemappedSeqIdAndCompId,
                                           splitPolySeqRstForMultimers,
                                           splitPolySeqRstForExactNoes,
                                           retrieveRemappedChainId,
                                           splitPolySeqRstForNonPoly,
                                           retrieveRemappedNonPoly,
                                           splitPolySeqRstForBranched,
                                           retrieveOriginalSeqIdFromMRMap)
    from wwpdb.utils.nmr.CifToNmrStar import (get_first_sf_tag,
                                              set_sf_tag)
    from wwpdb.utils.nmr.nef.NefTranslator import NefTranslator
    from wwpdb.utils.nmr.io.CifReader import CifReader
    from wwpdb.utils.nmr.mr.ParserListenerUtil import (coordAssemblyChecker,
                                                       translateToStdResName,
                                                       translateToStdAtomName,
                                                       translateToStdAtomNameNoRef,
                                                       translateToStdAtomNameWithRef,
                                                       translateToLigandName,
                                                       backTranslateFromStdResName,
                                                       hasInterChainRestraint,
                                                       isLongRangeRestraint,
                                                       isCyclicPolymer,
                                                       isStructConn,
                                                       guessCompIdFromAtomId,
                                                       getMetalCoordOf,
                                                       getRestraintName,
                                                       contentSubtypeOf,
                                                       incListIdCounter,
                                                       decListIdCounter,
                                                       getSfDictOf,
                                                       getSaveframe,
                                                       getPkLoop,
                                                       getAltLoops,
                                                       getAuxLoops,
                                                       getPkRow,
                                                       getAltPkRow,
                                                       getPkGenCharRow,
                                                       getPkCharRow,
                                                       getPkChemShiftRow,
                                                       getSpectralDimRow,
                                                       getSpectralDimTransferRow,
                                                       getMaxEffDigits,
                                                       getStarAtom,
                                                       roundString)
except ImportError:
    from nmr.NmrDpConstant import (LARGE_ASYM_ID,
                                   EMPTY_VALUE,
                                   ELEMENT_SYMBOLS,
                                   STD_MON_DICT,
                                   PROTON_BEGIN_CODE,
                                   PSE_PRO_BEGIN_CODE,
                                   AMINO_PROTON_CODE,
                                   CARBOXYL_CODE,
                                   ZINC_ION_CODE,
                                   CALCIUM_ION_CODE,
                                   ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS,
                                   REPRESENTATIVE_MODEL_ID,
                                   REPRESENTATIVE_ALT_ID,
                                   MAX_PREF_LABEL_SCHEME_COUNT,
                                   MAX_ALLOWED_EXT_SEQ,
                                   UNREAL_AUTH_SEQ_NUM,
                                   MAX_CONFLICT_ATTEMPT,
                                   CS_RESTRAINT_RANGE,
                                   CS_RESTRAINT_ERROR,
                                   WEIGHT_RANGE,
                                   CS_RANGE_MIN,
                                   CS_RANGE_MAX,
                                   CS_ERROR_MIN,
                                   CS_ERROR_MAX,
                                   WEIGHT_RANGE_MIN,
                                   WEIGHT_RANGE_MAX,
                                   HEME_LIKE_RES_NAMES,
                                   SPECTRAL_DIM_TEMPLATE,
                                   ASSIGNMENT_SEPARATOR_PAT,
                                   ASSIGNMENT_RESID_PAT,
                                   ASSIGNMENT_HALF_SPIN_NUCLEUS,
                                   C_CARBONYL_CENTER_MAX_TOR,
                                   C_CARBONYL_CENTER_MAX,
                                   C_CARBONYL_CENTER_MIN,
                                   C_AROMATIC_CENTER_MAX,
                                   C_AROMATIC_CENTER_MIN,
                                   C_AROMATIC_CENTER_MIN_TOR,
                                   N_AMIDE_CENTER_MAX,
                                   N_AMIDE_CENTER_MIN,
                                   C_ALL_CENTER_MAX,
                                   C_ALL_CENTER_MIN,
                                   C_ALIPHATIC_CENTER_MAX,
                                   C_ALIPHATIC_CENTER_MIN,
                                   C_METHYL_CENTER_MAX,
                                   C_METHYL_CENTER_MIN,
                                   HN_AROMATIC_CENTER_MAX,
                                   HN_AROMATIC_CENTER_MIN,
                                   H_IMIDE_CENTER_MAX,
                                   H_IMIDE_CENTER_MIN,
                                   H_ALL_CENTER_MAX,
                                   H_ALL_CENTER_MIN,
                                   H_ALIPHATIC_CENTER_MAX,
                                   H_ALIPHATIC_CENTER_MIN,
                                   H_METHYL_CENTER_MAX,
                                   H_METHYL_CENTER_MIN,
                                   MIN_CORRCOEF_FOR_ONE_BOND_TRANSFER,
                                   INSTRUCTION_FOR_FULL_SEQUENCE)
    from nmr.AlignUtil import (getOneLetterCode,
                               updatePolySeqRst,
                               revertPolySeqRst,
                               sortPolySeqRst,
                               syncCompIdOfPolySeqRst,
                               alignPolymerSequence,
                               alignPolymerSequenceWithConflicts,
                               assignPolymerSequence,
                               trimSequenceAlignment,
                               updatePolySeqRstAmbig,
                               retrieveAtomIdentFromMRMap,
                               retrieveAtomIdFromMRMap,
                               retrieveRemappedSeqId,
                               retrieveRemappedSeqIdAndCompId,
                               splitPolySeqRstForMultimers,
                               splitPolySeqRstForExactNoes,
                               retrieveRemappedChainId,
                               splitPolySeqRstForNonPoly,
                               retrieveRemappedNonPoly,
                               splitPolySeqRstForBranched,
                               retrieveOriginalSeqIdFromMRMap)
    from nmr.CifToNmrStar import (get_first_sf_tag,
                                  set_sf_tag)
    from nmr.nef.NefTranslator import NefTranslator
    from nmr.io.CifReader import CifReader
    from nmr.mr.ParserListenerUtil import (coordAssemblyChecker,
                                           translateToStdResName,
                                           translateToStdAtomName,
                                           translateToStdAtomNameNoRef,
                                           translateToStdAtomNameWithRef,
                                           translateToLigandName,
                                           backTranslateFromStdResName,
                                           hasInterChainRestraint,
                                           isLongRangeRestraint,
                                           isCyclicPolymer,
                                           isStructConn,
                                           guessCompIdFromAtomId,
                                           getMetalCoordOf,
                                           getRestraintName,
                                           contentSubtypeOf,
                                           incListIdCounter,
                                           decListIdCounter,
                                           getSfDictOf,
                                           getSaveframe,
                                           getPkLoop,
                                           getAltLoops,
                                           getAuxLoops,
                                           getPkRow,
                                           getAltPkRow,
                                           getPkGenCharRow,
                                           getPkCharRow,
                                           getPkChemShiftRow,
                                           getSpectralDimRow,
                                           getSpectralDimTransferRow,
                                           getMaxEffDigits,
                                           getStarAtom,
                                           roundString)


DIM_TRANSFER_PAT_2D = ((0, 1), )
DIM_TRANSFER_PAT_3D = ((0, 1), (1, 2), (2, 0))
DIM_TRANSFER_PAT_4D = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))


def guess_primary_dim_transfer_type(solid_state_nmr: bool, data_file_name: str, d: int, cur_spectral_dim: dict) -> str:
    """ Return expected primary dimensional transfer type from given frequencies.
    """

    file_name = data_file_name.lower()

    is_noesy = 'noe' in file_name or 'roe' in file_name
    no_aromatic = all(_dict['spectral_region'] != 'C-aromatic' for _dict in cur_spectral_dim.values())

    acq_dim_id = 1

    for __d, __v in cur_spectral_dim.items():
        if 'freq_hint' in __v:
            __v['freq_hint'] = numpy.array(__v['freq_hint'], dtype=float)  # list -> numpy array
            if __v['freq_hint'].size > 0:
                center = numpy.mean(__v['freq_hint'])
                max_ppm = __v['freq_hint'].max()
                min_ppm = __v['freq_hint'].min()

                if __v['atom_isotope_number'] is None:
                    if (C_AROMATIC_CENTER_MIN_TOR if 'aro' in file_name or 'anoe' in file_name else C_AROMATIC_CENTER_MIN)\
                       < center <= C_AROMATIC_CENTER_MAX and min_ppm > C_ALL_CENTER_MAX:
                        __v['atom_type'] = 'C'
                        __v['atom_isotope_number'] = 13
                        __v['axis_code'] = 'C-aromatic'
                    elif N_AMIDE_CENTER_MIN < center <= N_AMIDE_CENTER_MAX and min_ppm > C_ALL_CENTER_MIN\
                            and max_ppm < C_CARBONYL_CENTER_MIN:
                        __v['atom_type'] = 'N'
                        __v['atom_isotope_number'] = 15
                        __v['axis_code'] = 'N'
                    elif C_CARBONYL_CENTER_MIN <= center <= C_CARBONYL_CENTER_MAX:
                        __v['atom_type'] = 'C'
                        __v['atom_isotope_number'] = 13
                        __v['axis_code'] = 'CO'
                    elif HN_AROMATIC_CENTER_MIN < center <= HN_AROMATIC_CENTER_MAX:
                        __v['atom_type'] = 'H'
                        __v['atom_isotope_number'] = 1
                        __v['axis_code'] = 'HN/H-aromatic'
                    elif H_ALL_CENTER_MIN < center <= H_ALL_CENTER_MAX:
                        __v['atom_type'] = 'H'
                        __v['atom_isotope_number'] = 1
                        __v['axis_code'] = 'H'  # all
                    elif H_ALIPHATIC_CENTER_MIN < center <= H_ALIPHATIC_CENTER_MAX:
                        __v['atom_type'] = 'H'
                        __v['atom_isotope_number'] = 1
                        __v['axis_code'] = 'H-aliphatic' if max_ppm < 7 else 'H'
                    elif H_METHYL_CENTER_MIN < center <= H_METHYL_CENTER_MAX:
                        __v['atom_type'] = 'H'
                        __v['atom_isotope_number'] = 1
                        __v['axis_code'] = 'H-methyl' if max_ppm < 3 else 'H-aliphatic'
                    elif C_ALL_CENTER_MIN < center <= C_ALL_CENTER_MAX:
                        __v['atom_type'] = 'C'
                        __v['atom_isotope_number'] = 13
                        __v['axis_code'] = 'C'  # all
                    elif C_ALIPHATIC_CENTER_MIN < center <= C_ALIPHATIC_CENTER_MAX:
                        __v['atom_type'] = 'C'
                        __v['atom_isotope_number'] = 13
                        __v['axis_code'] = 'C-aliphatic'
                    elif C_METHYL_CENTER_MIN < center <= C_METHYL_CENTER_MAX:
                        __v['atom_type'] = 'C'
                        __v['atom_isotope_number'] = 13
                        __v['axis_code'] = 'C-methyl'
                    elif C_METHYL_CENTER_MIN < min_ppm and max_ppm < C_CARBONYL_CENTER_MAX_TOR:
                        __v['atom_type'] = 'C'
                        __v['atom_isotope_number'] = 13
                        __v['axis_code'] = 'C'  # all

                isotope_number = __v['atom_isotope_number']

                if isotope_number is not None:
                    __v['acquisition'] = 'yes' if __d == acq_dim_id\
                        and (isotope_number == 1 or (isotope_number == 13 and solid_state_nmr)) else 'no'

                    if __d == 1 and __v['acquisition'] == 'no':
                        acq_dim_id = d

                    __v['under_sampling_type'] = 'not observed' if __v['acquisition'] == 'yes' else 'aliased'

            if __v['spectral_region'] is None and __v['freq_hint'].size > 0:
                atom_type = __v['atom_type']
                if C_AROMATIC_CENTER_MIN_TOR < center <= C_AROMATIC_CENTER_MAX and min_ppm > C_ALL_CENTER_MAX and atom_type == 'C':
                    __v['spectral_region'] = 'C-aromatic'
                elif N_AMIDE_CENTER_MIN < center <= N_AMIDE_CENTER_MAX and min_ppm > C_ALL_CENTER_MAX:
                    if atom_type == 'N':
                        __v['spectral_region'] = 'N'
                    if atom_type == 'C':
                        __v['spectral_region'] = 'C-aromatic'
                elif C_CARBONYL_CENTER_MIN <= center <= C_CARBONYL_CENTER_MAX and atom_type == 'C':
                    __v['spectral_region'] = 'CO'
                elif HN_AROMATIC_CENTER_MIN < center <= HN_AROMATIC_CENTER_MAX and atom_type == 'H':
                    __v['spectral_region'] = 'HN/H-aromatic'
                elif H_ALL_CENTER_MIN < center <= H_ALL_CENTER_MAX and atom_type == 'H':
                    __v['spectral_region'] = 'H'  # all
                elif H_ALIPHATIC_CENTER_MIN < center <= H_ALIPHATIC_CENTER_MAX and atom_type == 'H':
                    __v['spectral_region'] = 'H-aliphatic' if max_ppm < 7 else 'H'
                elif H_METHYL_CENTER_MIN < center <= H_METHYL_CENTER_MAX and atom_type == 'H':
                    __v['spectral_region'] = 'H-methyl' if max_ppm < 3 else 'H-aliphatic'
                elif C_ALL_CENTER_MIN < center <= C_ALL_CENTER_MAX and atom_type == 'C':
                    __v['spectral_region'] = 'C'  # all
                elif C_ALIPHATIC_CENTER_MIN < center <= C_ALIPHATIC_CENTER_MAX and atom_type == 'C':
                    __v['spectral_region'] = 'C-aliphatic'
                elif C_METHYL_CENTER_MIN < center <= C_METHYL_CENTER_MAX and atom_type == 'C':
                    __v['spectral_region'] = 'C-methyl'
                elif __v['atom_isotope_number'] == 1:
                    __v['spectral_region'] = 'H'
                elif __v['atom_isotope_number'] == 13:
                    __v['spectral_region'] = 'C'
                elif __v['atom_isotope_number'] == 15:
                    __v['spectral_region'] = 'N'
                elif __v['atom_isotope_number'] == 19:
                    __v['spectral_region'] = 'F'
                elif __v['atom_isotope_number'] == 31:
                    __v['spectral_region'] = 'P'

            if __v['freq_hint'].size > 0 and d > 2 and __d >= 2\
               and not solid_state_nmr and __v['atom_isotope_number'] == 13:
                min_ppm = __v['freq_hint'].min()
                width = max_ppm - min_ppm
                if center < 100.0 and width < 50.0:
                    __v['under_sampling_type'] = 'folded'

            if __v['spectrometer_frequency'] is None and 'obs_freq_hint' in __v and len(__v['obs_freq_hint']) > 0:
                __v['spectrometer_frequency'] = collections.Counter(__v['obs_freq_hint']).most_common()[0][0]

            if 'obs_freq_hint' in __v:
                del __v['obs_freq_hint']

            if __v['spectrometer_frequency'] is not None and __v['sweep_width_units'] == 'ppm':
                row = [str(__v['sweep_width']), str(__v['spectrometer_frequency'])]
                max_eff_digits = getMaxEffDigits(row)

                __v['sweep_width'] = float(roundString(str(__v['sweep_width'] * __v['spectrometer_frequency']),
                                                       max_eff_digits))
                __v['sweep_width_units'] = 'Hz'

    for __v in cur_spectral_dim.values():
        if __v['axis_code'] == 'HN/H-aromatic':
            has_a = any(___v['spectral_region'] == 'C-aromatic' for ___v in cur_spectral_dim.values())
            __v['axis_code'] = 'H-aromatic' if has_a else 'H'
        if __v['spectral_region'] == 'HN/H-aromatic':
            has_a = any(___v['spectral_region'] == 'C-aromatic' for ___v in cur_spectral_dim.values())
            __v['spectral_region'] = 'H-aromatic' if has_a else 'H'
            __v['_spectral_region'] = 'H-aromatic' if has_a else 'HN'
        else:
            __v['_spectral_region'] = __v['spectral_region']

    cur_spectral_dim_transfer = []

    dim_to_code = {1: 'x', 2: 'y', 3: 'z', 4: 'a'}

    # onebond: 'Any transfer that connects only directly bonded atoms in this experiment'
    for _dim_id1, _dict1 in cur_spectral_dim.items():
        _region1 = _dict1['_spectral_region']
        if _region1 in ('HN', 'H-aliphatic', 'H-aromatic', 'H-methyl'):
            cases = 0
            max_corrcoef = 0.0
            for _dim_id2, _dict2 in cur_spectral_dim.items():
                _region2 = _dict2['_spectral_region']
                if (_region1 == 'HN' and _region2 == 'N')\
                   or ((_region1 == 'H-aliphatic' or (is_noesy and no_aromatic and _region1 == 'H'))
                       and (_region2 == 'C-aliphatic' or (is_noesy and _region2 == 'C')))\
                   or (_region1 == 'H-aromatic' and _region2 == 'C-aromatic')\
                   or (_region1 == 'H-methyl' and _region2 == 'C-methyl'):
                    if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                        if not any(True for _transfer in cur_spectral_dim_transfer
                                   if _transfer['type'] == 'onebond'
                                   and (_dim_id1 in [_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']]
                                        or _dim_id2 in [_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']])):
                            if _dict1['freq_hint'].size < 2\
                               or numpy.max(_dict1['freq_hint']) == numpy.min(_dict1['freq_hint'])\
                               or numpy.max(_dict2['freq_hint']) == numpy.min(_dict2['freq_hint']):
                                continue
                            _corrcoef = numpy.corrcoef(_dict1['freq_hint'], _dict2['freq_hint'])[0][1]
                            if _corrcoef < 0.0 or (_corrcoef < MIN_CORRCOEF_FOR_ONE_BOND_TRANSFER and no_aromatic):
                                continue
                            cases += 1
                            max_corrcoef = max(max_corrcoef, _corrcoef)

                            if d == 3:
                                for _dim_id3, _dict3 in cur_spectral_dim.items():
                                    if _dim_id3 in (_dim_id1, _dim_id2):
                                        continue
                                    if _dict3['atom_type'] != _dict1['atom_type']:
                                        continue
                                    if not any(True for _transfer in cur_spectral_dim_transfer
                                       if _transfer['type'] == 'onebond'
                                       and (_dim_id3 in [_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']]
                                            or _dim_id2 in [_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']])):
                                        if _dict3['freq_hint'].size < 2\
                                           or numpy.max(_dict3['freq_hint']) == numpy.min(_dict3['freq_hint'])\
                                           or numpy.max(_dict2['freq_hint']) == numpy.min(_dict2['freq_hint']):
                                            continue
                                        _corrcoef = numpy.corrcoef(_dict3['freq_hint'], _dict2['freq_hint'])[0][1]
                                        if _corrcoef < 0.0 or (_corrcoef < MIN_CORRCOEF_FOR_ONE_BOND_TRANSFER and no_aromatic):
                                            continue
                                        cases += 1
                                        max_corrcoef = max(max_corrcoef, _corrcoef)

            if cases == 1:
                for _dim_id2, _dict2 in cur_spectral_dim.items():
                    _region2 = _dict2['_spectral_region']
                    if (_region1 == 'HN' and _region2 == 'N')\
                       or ((_region1 == 'H-aliphatic' or (is_noesy and no_aromatic and _region1 == 'H'))
                           and (_region2 == 'C-aliphatic' or (is_noesy and _region2 == 'C')))\
                       or (_region1 == 'H-aromatic' and _region2 == 'C-aromatic')\
                       or (_region1 == 'H-methyl' and _region2 == 'C-methyl'):
                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                            if not any(True for _transfer in cur_spectral_dim_transfer
                                       if _transfer['type'] == 'onebond'
                                       and (_dim_id1 in [_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']]
                                            or _dim_id2 in [_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']])):
                                transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                            'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                            'type': 'onebond',
                                            'indirect': 'no'}
                                if transfer in cur_spectral_dim_transfer:
                                    continue
                                cur_spectral_dim_transfer.append(transfer)

            elif cases > 1:
                for _dim_id2, _dict2 in cur_spectral_dim.items():
                    _region2 = _dict2['_spectral_region']
                    if (_region1 == 'HN' and _region2 == 'N')\
                       or ((_region1 == 'H-aliphatic' or (is_noesy and no_aromatic and _region1 == 'H'))
                           and (_region2 == 'C-aliphatic' or (is_noesy and _region2 == 'C')))\
                       or (_region1 == 'H-aromatic' and _region2 == 'C-aromatic')\
                       or (_region1 == 'H-methyl' and _region2 == 'C-methyl'):
                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                            if not any(True for _transfer in cur_spectral_dim_transfer
                                       if _transfer['type'] == 'onebond'
                                       and (_dim_id1 in [_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']]
                                            or _dim_id2 in [_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']])):
                                if _dict1['freq_hint'].size < 2\
                                   or numpy.max(_dict1['freq_hint']) == numpy.min(_dict1['freq_hint'])\
                                   or numpy.max(_dict2['freq_hint']) == numpy.min(_dict2['freq_hint']):
                                    continue
                                _corrcoef = numpy.corrcoef(_dict1['freq_hint'], _dict2['freq_hint'])[0][1]
                                if _corrcoef < 0.0 or (_corrcoef < max_corrcoef and no_aromatic):
                                    continue
                                transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                            'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                            'type': 'onebond',
                                            'indirect': 'no'}
                                if transfer in cur_spectral_dim_transfer:
                                    continue
                                cur_spectral_dim_transfer.append(transfer)

    for _dim_id1, _dict1 in cur_spectral_dim.items():
        _region1 = _dict1['_spectral_region']
        if _region1 in ('HN', 'H-aliphatic', 'H-aromatic', 'H-methyl'):
            for _dim_id2, _dict2 in cur_spectral_dim.items():
                _region2 = _dict2['_spectral_region']
                if (_region1 == 'HN' and _region2 == 'N')\
                   or ((_region1 == 'H-aliphatic' or (is_noesy and no_aromatic and _region1 == 'H'))
                       and (_region2 == 'C-aliphatic' or (is_noesy and _region2 == 'C')))\
                   or (_region1 == 'H-aromatic' and _region2 == 'C-aromatic')\
                   or (_region1 == 'H-methyl' and _region2 == 'C-methyl'):
                    if _dict1['acquisition'] == 'no' and _dict2['acquisition'] == 'no':
                        if not any(True for _transfer in cur_spectral_dim_transfer
                                   if _transfer['type'] == 'onebond'
                                   and (_dim_id1 in [_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']]
                                        or _dim_id2 in [_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']])):
                            if _dict1['freq_hint'].size < 2\
                               or numpy.max(_dict1['freq_hint']) == numpy.min(_dict1['freq_hint'])\
                               or numpy.max(_dict2['freq_hint']) == numpy.min(_dict2['freq_hint']):
                                continue
                            _corrcoef = numpy.corrcoef(_dict1['freq_hint'], _dict2['freq_hint'])[0][1]
                            if _corrcoef < 0.0 or (_corrcoef < max_corrcoef and no_aromatic):
                                continue
                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                        'type': 'onebond',
                                        'indirect': 'no'}
                            if transfer in cur_spectral_dim_transfer:
                                continue
                            cur_spectral_dim_transfer.append(transfer)

    # jcoupling: 'Transfer via direct J coupling over one or more bonds'
    if 'cosy' in file_name:
        if d == 2:
            for _dim_id1, _dict1 in cur_spectral_dim.items():
                _iso_num1 = _dict1['atom_isotope_number']
                if _iso_num1 == 1:
                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                        if _dim_id1 == _dim_id2 or _iso_num1 != _dict2['atom_isotope_number']:
                            continue
                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                        'type': 'jcoupling',
                                        'indirect': 'no'}
                        if transfer in cur_spectral_dim_transfer:
                            continue
                        cur_spectral_dim_transfer.append(transfer)

        elif d == 3:
            for _dim_id1, _dict1 in cur_spectral_dim.items():
                _region1 = _dict1['_spectral_region']
                if _region1 == 'HN':
                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                        if _dim_id1 == _dim_id2 or _dict2['atom_isotope_number'] not in (1, 13):
                            continue
                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                        'type': 'jcoupling',
                                        'indirect': 'yes'}
                            if transfer in cur_spectral_dim_transfer:
                                continue
                            cur_spectral_dim_transfer.append(transfer)
            for _dim_id1, _dict1 in cur_spectral_dim.items():
                _region1 = _dict1['_spectral_region']
                if _region1 == 'H-aliphatic':
                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                        _isotope2 = _dict2['atom_isotope']
                        if _dim_id1 == _dim_id2 or _isotope2 not in (1, 13):
                            continue
                        if _isotope2 == 13\
                           and not any(True for _transfer in cur_spectral_dim_transfer
                                       if _transfer['type'] == 'onebond'
                                       and {_dim_id1, _dim_id2} == {_transfer['spectral_dim_id_1'],
                                                                    _transfer['spectral_dim_id_2']}):
                            if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                            'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                            'type': 'jcoupling',
                                            'indirect': 'yes' if _isotope2 == 1 else 'no'}
                                if transfer in cur_spectral_dim_transfer:
                                    continue
                                cur_spectral_dim_transfer.append(transfer)

        elif d == 4:
            for _dim_id1, _dict1 in cur_spectral_dim.items():
                _region1 = _dict1['_spectral_region']
                if _region1 == 'HN':
                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                        if _dim_id1 == _dim_id2 or _dict2['atom_isotope_number'] != 1:
                            continue
                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                        'type': 'jcoupling',
                                        'indirect': 'yes'}
                            if transfer in cur_spectral_dim_transfer:
                                continue
                            cur_spectral_dim_transfer.append(transfer)
            for _dim_id1, _dict1 in cur_spectral_dim.items():
                _region1 = _dict1['_spectral_region']
                if _region1 == 'H-aliphatic':
                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                        if _dim_id1 == _dim_id2 or _dict2['atom_isotope_number'] != 1:
                            continue
                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                        'type': 'jcoupling',
                                        'indirect': 'no'}
                            if transfer in cur_spectral_dim_transfer:
                                continue
                            cur_spectral_dim_transfer.append(transfer)

    # jmultibond: 'Transfer via direct J coupling over multiple bonds'

    # relayed: 'Transfer via multiple successive J coupling steps (TOCSY relay)'
    if 'tocsy' in file_name:
        if d == 2:
            for _dim_id1, _dict1 in cur_spectral_dim.items():
                _iso_num1 = _dict1['atom_isotope_number']
                if _iso_num1 == 1:
                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                        if _dim_id1 == _dim_id2 or _iso_num1 != _dict2['atom_isotope_number']:
                            continue
                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                        'type': 'relayed',
                                        'indirect': 'no'}
                        if transfer in cur_spectral_dim_transfer:
                            continue
                        cur_spectral_dim_transfer.append(transfer)

        elif d == 3:
            for _dim_id1, _dict1 in cur_spectral_dim.items():
                _region1 = _dict1['_spectral_region']
                if _region1 == 'HN':
                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                        if _dim_id1 == _dim_id2 or _dict2['atom_isotope_number'] not in (1, 13):
                            continue
                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                        'type': 'relayed',
                                        'indirect': 'yes'}
                            if transfer in cur_spectral_dim_transfer:
                                continue
                            cur_spectral_dim_transfer.append(transfer)
            for _dim_id1, _dict1 in cur_spectral_dim.items():
                _region1 = _dict1['_spectral_region']
                if _region1 == 'H-aliphatic':
                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                        _isotope2 = _dict2['atom_isotope']
                        if _dim_id1 == _dim_id2 or _isotope2 not in (1, 13):
                            continue
                        if _isotope2 == 13\
                           and not any(True for _transfer in cur_spectral_dim_transfer
                                       if _transfer['type'] == 'onebond'
                                       and {_dim_id1, _dim_id2} == {_transfer['spectral_dim_id_1'],
                                                                    _transfer['spectral_dim_id_2']}):
                            if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                            'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                            'type': 'relayed',
                                            'indirect': 'yes' if _isotope2 == 1 else 'no'}
                                if transfer in cur_spectral_dim_transfer:
                                    continue
                                cur_spectral_dim_transfer.append(transfer)

        elif d == 4:
            for _dim_id1, _dict1 in cur_spectral_dim.items():
                _region1 = _dict1['_spectral_region']
                if _region1 == 'HN':
                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                        if _dim_id1 == _dim_id2 or _dict2['atom_isotope_number'] != 1:
                            continue
                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                        'type': 'relayed',
                                        'indirect': 'yes'}
                            if transfer in cur_spectral_dim_transfer:
                                continue
                            cur_spectral_dim_transfer.append(transfer)
            for _dim_id1, _dict1 in cur_spectral_dim.items():
                _region1 = _dict1['_spectral_region']
                if _region1 == 'H-aliphatic':
                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                        if _dim_id1 == _dim_id2 or _dict2['atom_isotope_number'] != 1:
                            continue
                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                        'type': 'relayed',
                                        'indirect': 'no'}
                            if transfer in cur_spectral_dim_transfer:
                                continue
                            cur_spectral_dim_transfer.append(transfer)

    # relayed-alternate: 'Relayed transfer where peaks from an odd resp. even number of transfer steps have opposite sign'

    # through-space: 'Any transfer that does not go through the covalent bonded skeleton'
    if 'noe' in file_name or 'roe' in file_name:
        for _dim_id1, _dict1 in cur_spectral_dim.items():
            _region1 = _dict1['_spectral_region']
            if _region1 in ('H', 'HN', 'H-aliphatic', 'H-aromatic', 'H-methyl'):
                for _dim_id2, _dict2 in cur_spectral_dim.items():
                    if _dim_id1 == _dim_id2 or _dict1['atom_isotope_number'] != _dict2['atom_isotope_number']:
                        continue
                    if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                        transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                    'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                    'type': 'through-space',
                                    'indirect': 'yes'}
                        if transfer in cur_spectral_dim_transfer:
                            continue
                        cur_spectral_dim_transfer.append(transfer)
                        if d == 2 and _region1 == 'H-aliphatic':
                            _dict1['spectral_region'] = _dict2['spectral_region'] = 'H'  # all
                            if 'H-aliphatic' in (_dict1['axis_code'], _dict2['axis_code']):
                                _dict1['axis_code'] = f'H{dim_to_code[_dim_id1]}'
                                _dict2['axis_code'] = f'H{dim_to_code[_dim_id2]}'
                        if d == 3:
                            _transfer = next((_transfer for _transfer in cur_spectral_dim_transfer
                                              if _transfer['type'] == 'onebond'), None)
                            if _transfer is not None:
                                if _dim_id1 not in (_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']):
                                    if _region1 == 'H-aliphatic':
                                        _dict1['spectral_region'] = 'H'  # all
                                        if _dict1['axis_code'] == 'H-aliphatic':
                                            _dict1['axis_code'] = 'H'
                                elif _region1 == 'HN':
                                    _dict1['spectral_region'] = _dict1['axis_code'] = 'HN'
                                if _dim_id2 not in (_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']):
                                    if _region1 == 'H-aliphatic':
                                        _dict2['spectral_region'] = 'H'  # all
                                        if _dict2['axis_code'] == 'H-aliphatic':
                                            _dict2['axis_code'] = 'H'
                                elif _dict2['_spectral_region'] == 'HN':
                                    _dict2['spectral_region'] = _dict2['axis_code'] = 'HN'

    if solid_state_nmr and d == 2:
        if 'rfdr' in file_name or 'darr' in file_name:
            for _dim_id1, _dict1 in cur_spectral_dim.items():
                _iso_num1 = _dict1['atom_isotope_number']
                if _iso_num1 in (1, 13):
                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                        if _dim_id1 == _dim_id2 or _iso_num1 != _dict2['atom_isotope_number']:
                            continue
                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                        'type': 'through-space',
                                        'indirect': 'yes'}
                            if transfer in cur_spectral_dim_transfer:
                                continue
                            cur_spectral_dim_transfer.append(transfer)
                            if d == 2:  # and _dict1['spectral_region'] == _dict2['spectral_region']:
                                nuc = _dict1[1]['spectral_regison'][0]
                                _dict1['spectral_region'] = _dict2['spectral_region'] = nuc
                                if _dict1['axis_code'] == _dict2['axis_code']:
                                    _dict1['axis_code'] = f'{nuc}{dim_to_code[_dim_id1]}'
                                    _dict2['axis_code'] = f'{nuc}{dim_to_code[_dim_id2]}'

        elif 'redor' in file_name or 'tedor' in file_name:
            for _dim_id1, _dict1 in cur_spectral_dim.items():
                _iso_num1 = _dict1['atom_isotope_number']
                if _iso_num1 in (13, 15, 19, 31):
                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                        _iso_num2 = _dict2['atom_isotope_number']
                        if _dim_id1 == _dim_id2 or _iso_num2 not in (13, 15, 19, 31) or _iso_num1 == _iso_num2:
                            continue
                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                        'type': 'through-space',
                                        'indirect': 'yes'}
                            if transfer in cur_spectral_dim_transfer:
                                continue
                            cur_spectral_dim_transfer.append(transfer)

    for _dim_id1, _dict1 in cur_spectral_dim.items():
        _region1 = _dict1['_spectral_region']
        if _region1 in ('H', 'HN', 'H-aliphatic', 'H-aromatic', 'H-methyl') and d > 2:
            for _dim_id2, _dict2 in cur_spectral_dim.items():
                if _dim_id1 == _dim_id2 or _dict1['atom_isotope_number'] != _dict2['atom_isotope_number']:
                    continue
                if not any(True for _transfer in cur_spectral_dim_transfer
                           if {_dim_id1, _dim_id2} == {_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']}):
                    if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                        transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                    'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                    'type':
                                    'through-space' if 'long_range' in cur_spectral_dim[1] else 'through-space?',
                                    'indirect': 'yes'}
                        if transfer in cur_spectral_dim_transfer:
                            continue
                        cur_spectral_dim_transfer.append(transfer)

    for _dim_id1, _dict1 in cur_spectral_dim.items():
        _region1 = _dict1['_spectral_region']
        if _region1 == 'H' and d == 2:  # all
            for _dim_id2, _dict2 in cur_spectral_dim.items():
                if _dim_id1 == _dim_id2:
                    continue
                if 'long_range' not in cur_spectral_dim[1]:
                    if _dict2['_spectral_region'] != _region1:
                        continue
                if not any(True for _transfer in cur_spectral_dim_transfer
                           if {_dim_id1, _dim_id2} == {_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']}):
                    if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                        transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                    'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                    'type':
                                    'through-space' if 'long_range' in cur_spectral_dim[1] else 'through-space?',
                                    'indirect': 'yes'}
                        if transfer in cur_spectral_dim_transfer:
                            continue
                        cur_spectral_dim_transfer.append(transfer)

    if solid_state_nmr and d == 2:
        for _dim_id1, _dict1 in cur_spectral_dim.items():
            _region1 = _dict1['_spectral_region']
            if _region1 == 'C':  # all
                for _dim_id2, _dict2 in cur_spectral_dim.items():
                    if _dim_id1 == _dim_id2:
                        continue
                    if 'long_range' not in cur_spectral_dim[1]:
                        if _dict2['_spectral_region'] != _region1:
                            continue
                    if not any(True for _transfer in cur_spectral_dim_transfer
                               if {_dim_id1, _dim_id2} == {_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']}):
                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                        'type':
                                        'through-space' if 'long_range' in cur_spectral_dim[1] else 'through-space?',
                                        'indirect': 'yes'}
                            if transfer in cur_spectral_dim_transfer:
                                continue
                            cur_spectral_dim_transfer.append(transfer)

        for _dim_id1, _dict1 in cur_spectral_dim.items():
            _iso_num1 = _dict1['atom_isotope_number']
            if _iso_num1 in (13, 15, 19, 31):
                for _dim_id2, _dict2 in cur_spectral_dim.items():
                    _iso_num2 = _dict2['atom_isotope_number']
                    if _dim_id1 == _dim_id2 or _iso_num2 not in (13, 15, 19, 31) or _iso_num1 == _iso_num2:
                        continue
                    if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                        transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                    'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                    'type':
                                    'through-space' if 'long_range' in cur_spectral_dim[1] else 'through-space?',
                                    'indirect': 'yes'}
                        if transfer in cur_spectral_dim_transfer:
                            continue
                        cur_spectral_dim_transfer.append(transfer)

    for __v in cur_spectral_dim.values():
        if 'freq_hint' in __v:
            del __v['freq_hint']

    primary_dim_transfer = ''
    for transfer in cur_spectral_dim_transfer:
        if transfer['type'] == 'through-space':
            primary_dim_transfer = transfer['type']
            break
        if transfer['type'] == 'through-space?'\
                and primary_dim_transfer != 'through-space':
            primary_dim_transfer = transfer['type']
        elif transfer['type'].startswith('relayed')\
                and primary_dim_transfer not in ('through-space', 'through-space?'):
            primary_dim_transfer = transfer['type']
        elif transfer['type'] == 'jmultibond'\
                and primary_dim_transfer not in ('through-space', 'through-space?', 'relayed', 'relayed-alternate'):
            primary_dim_transfer = transfer['type']
        elif transfer['type'] == 'jcoupling'\
                and primary_dim_transfer not in ('through-space', 'through-space?', 'relayed', 'relayed-alternate', 'jmultibond'):
            primary_dim_transfer = transfer['type']
        elif transfer['type'] == 'onebond'\
                and primary_dim_transfer == '':
            primary_dim_transfer = transfer['type']

    return primary_dim_transfer


class BasePKParserListener():
    """ ParserLister base class for any peak list file.
    """
    __slots__ = ('representativeModelId',
                 'representativeAltId',
                 '__mrAtomNameMapping',
                 'cR',
                 '__hasCoord',
                 'ccU',
                 'modelNumName',
                 '__authSeqId',
                 'polySeq',
                 '__altPolySeq',
                 'nonPoly',
                 'branched',
                 '__coordAtomSite',
                 '__coordUnobsRes',
                 '__coordUnobsAtom',
                 '__labelToAuthSeq',
                 '__authToLabelSeq',
                 '__authToStarSeq',
                 '__authToOrigSeq',
                 '__modResidue',
                 '__splitLigand',
                 '__entityAssembly',
                 'exptlMethod',
                 '__offsetHolder',
                 'hasPolySeq',
                 'hasNonPoly',
                 'hasBranched',
                 'hasNonPolySeq',
                 'nonPolySeq',
                 'authAsymIdSet',
                 'compIdSet',
                 'altCompIdSet',
                 'polyPeptide',
                 'polyDeoxyribonucleotide',
                 'polyRibonucleotide',
                 'isFirstResidueAla',
                 'cyanaCompIdSet',
                 '__uniqAtomIdToSeqKey',
                 'csStat',
                 'nefT',
                 'pA',
                 'reasons',
                 '__preferAuthSeqCount',
                 '__preferLabelSeqCount',
                 'reasonsForReParsing',
                 'peaks2D',
                 'peaks3D',
                 'peaks4D',
                 'sfDict',
                 '__cachedDictForStarAtom',
                 'num_of_dim',
                 'acq_dim_id',
                 'cur_spectral_dim',
                 'spectral_dim',
                 'spectral_dim_transfer',
                 'atom_type_history',
                 'onebond_idx_history',
                 'jcoupling_idx_history',
                 'relayed_idx_history',
                 'listIdInternal',
                 'chainNumberDict',
                 'extResKey',
                 'polySeqRst',
                 'polySeqRstFailed',
                 'polySeqRstFailedAmbig',
                 'compIdMap',
                 'f',
                 '__allow_ext_seq')

    file_type = None
    software_name = None

    __debug = False
    __verbose_debug = False
    __internal = False
    __ignore_diagonal = True

    __createSfDict = False
    __enforcePeakRowFormat = False

    __shiftNonPosSeq = None
    __defaultSegId = None
    __defaultSegId__ = None

    __preferAuthSeq = True
    __extendAuthSeq = False

    seqAlign = None
    chainAssign = None

    cur_subtype = ''
    cur_list_id = -1
    use_peak_row_format = True
    null_value = None
    null_string = None

    spectrum_name = None

    # whether current assignment derived not from unreliable extra comment
    no_extra_comment = False

    # collection of atom selection set for multiple assignments to a peak
    atomSelectionSets = []
    asIsSets = []

    # collection of atom selection
    atomSelectionSet = []

    # collection of auxiliary atom selection
    auxAtomSelectionSet = []

    # collection of position selection
    positionSelection = []

    # collection of number selection
    numberSelection = []
    originalNumberSelection = []

    # collection of assignment (XEASY)
    assignmentSelection = []

    warningMessage = None

    # original source MR file name
    __originalFileName = '.'

    # list id counter
    __listIdCounter = {}

    # reserved list ids for NMR data remediation Phase 2
    __reservedListIds = {}

    # entry ID
    __entryId = '.'

    # list of assigned chemical shift loops
    __csLoops = None

    # tentative chemical shift values
    __tempCsValues = []

    # default saveframe name for error handling
    __def_err_sf_framecode = None

    def __init__(self, verbose: bool = True, log: IO = sys.stdout,
                 representativeModelId: int = REPRESENTATIVE_MODEL_ID,
                 representativeAltId: str = REPRESENTATIVE_ALT_ID,
                 mrAtomNameMapping: Optional[List[dict]] = None,
                 cR: Optional[CifReader] = None, caC: Optional[dict] = None,
                 nefT: NefTranslator = None,
                 reasons: Optional[dict] = None) -> None:

        self.exptlMethod = '.'
        self.representativeModelId = representativeModelId
        self.representativeAltId = representativeAltId
        # atom name mapping of public MR file between the archive coordinates and submitted ones
        self.__mrAtomNameMapping = None if mrAtomNameMapping is None or len(mrAtomNameMapping) == 0 else mrAtomNameMapping

        # CIF reader
        self.cR = cR
        self.__hasCoord = cR is not None

        self.nefT = nefT
        self.ccU = nefT.ccU
        self.csStat = nefT.csStat
        self.pA = nefT.pA

        if self.__hasCoord:
            ret = coordAssemblyChecker(verbose, log, representativeModelId, representativeAltId,
                                       cR, self.ccU, caC)
            self.modelNumName = ret['model_num_name']
            self.__authSeqId = ret['auth_seq_id']
            self.polySeq = ret['polymer_sequence']
            self.__altPolySeq = ret['alt_polymer_sequence']
            self.nonPoly = ret['non_polymer']
            self.branched = ret['branched']
            self.__coordAtomSite = ret['coord_atom_site']
            self.__coordUnobsRes = ret['coord_unobs_res']
            self.__coordUnobsAtom = ret['coord_unobs_atom'] if 'coord_unobs_atom' in ret else {}
            self.__labelToAuthSeq = ret['label_to_auth_seq']
            self.__authToLabelSeq = ret['auth_to_label_seq']
            self.__authToStarSeq = ret['auth_to_star_seq']
            self.__authToOrigSeq = ret['auth_to_orig_seq']
            self.__modResidue = ret['mod_residue']
            self.__splitLigand = ret['split_ligand']
            self.__entityAssembly = ret['entity_assembly']

            exptl = cR.getDictList('exptl')
            if len(exptl) > 0:
                for item in exptl:
                    if 'method' in item:
                        if 'NMR' in item['method']:
                            self.exptlMethod = item['method']
        else:
            self.exptlMethod = ''
            self.modelNumName = None
            self.__authSeqId = None
            self.polySeq = None
            self.__altPolySeq = None
            self.nonPoly = None
            self.branched = None
            self.nonPolySeq = None
            self.__coordAtomSite = None
            self.__coordUnobsRes = None
            self.__coordUnobsAtom = None
            self.__labelToAuthSeq = None
            self.__authToLabelSeq = None
            self.__authToStarSeq = None
            self.__authToOrigSeq = None
            self.__modResidue = None
            self.__splitLigand = None
            self.__entityAssembly = None

            self.exptlMethod = '.'

        self.polyPeptide = False
        self.polyDeoxyribonucleotide = False
        self.polyRibonucleotide = False

        self. __uniqAtomIdToSeqKey = None

        self.__offsetHolder = {}

        self.hasPolySeq = self.polySeq is not None and len(self.polySeq) > 0
        self.hasNonPoly = self.nonPoly is not None and len(self.nonPoly) > 0
        self.hasBranched = self.branched is not None and len(self.branched) > 0
        if self.hasNonPoly or self.hasBranched:
            self.hasNonPolySeq = True
            if self.hasNonPoly and self.hasBranched:
                self.nonPolySeq = self.nonPoly
                self.nonPolySeq.extend(self.branched)
            elif self.hasNonPoly:
                self.nonPolySeq = self.nonPoly
            else:
                self.nonPolySeq = self.branched

        else:
            self.hasNonPolySeq = False
            self.nonPolySeq = None

        if self.hasPolySeq:
            self.authAsymIdSet = set(ps['auth_chain_id'] for ps in self.polySeq)

            self.compIdSet = set()
            self.altCompIdSet = set()

            def is_data(array: list) -> bool:
                return not any(True for d in array if d in EMPTY_VALUE)

            for ps in self.polySeq:
                self.compIdSet.update(set(filter(is_data, ps['comp_id'])))
                if 'auth_comp_id' in ps and ps['comp_id'] != ps['auth_comp_id']:
                    self.altCompIdSet.update(set(filter(is_data, ps['auth_comp_id'])))
                if 'alt_comp_id' in ps and ps['comp_id'] != ps['alt_comp_id']:
                    self.altCompIdSet.update(set(filter(is_data, ps['alt_comp_id'])))

            if self.hasNonPolySeq:
                for np in self.nonPolySeq:
                    self.compIdSet.update(set(filter(is_data, np['comp_id'])))
                    if 'auth_comp_id' in np and np['comp_id'] != np['auth_comp_id']:
                        self.altCompIdSet.update(set(filter(is_data, np['auth_comp_id'])))
                    if 'alt_comp_id' in np and np['comp_id'] != np['alt_comp_id']:
                        self.altCompIdSet.update(set(filter(is_data, np['alt_comp_id'])))
                self.authAsymIdSet.update(set(np['auth_chain_id'] for np in self.nonPolySeq))

            for entity in self.__entityAssembly:
                if 'entity_poly_type' in entity:
                    poly_type = entity['entity_poly_type']
                    if poly_type.startswith('polypeptide'):
                        self.polyPeptide = True
                    elif poly_type == 'polydeoxyribonucleotide':
                        self.polyDeoxyribonucleotide = True
                    elif poly_type == 'polyribonucleotide':
                        self.polyRibonucleotide = True
            if 'ALA' in self.compIdSet:
                self.isFirstResidueAla = any(ps['comp_id'][0] == 'ALA' for ps in self.polySeq)

        else:
            self.authAsymIdSet = set()
            self.compIdSet = self.altCompIdSet = set(STD_MON_DICT.keys())
            self.isFirstResidueAla = False

        self.cyanaCompIdSet = set()
        for compId in self.compIdSet:
            self.cyanaCompIdSet |= backTranslateFromStdResName(compId)

        self.__uniqAtomIdToSeqKey = {}
        if self.hasNonPoly:
            atom_list = []
            for v in self.__coordAtomSite.values():
                atom_list.extend(v['atom_id'])
            common_atom_list = collections.Counter(atom_list).most_common()
            uniq_atom_ids = [atom_id for atom_id, count in common_atom_list if count == 1]
            if len(uniq_atom_ids) > 0:
                for k, v in self.__coordAtomSite.items():
                    if any(True for np in self.nonPoly if np['comp_id'][0] == v['comp_id']):
                        for atom_id in v['atom_id']:
                            if atom_id in uniq_atom_ids:
                                self.__uniqAtomIdToSeqKey[atom_id] = k

        # reasons for re-parsing request from the previous trial
        self.reasons = reasons
        self.__preferAuthSeqCount = 0
        self.__preferLabelSeqCount = 0

        self.reasonsForReParsing = {}  # reset to prevent interference from the previous run

        self.peaks2D = 0
        self.peaks3D = 0
        self.peaks4D = 0

        self.sfDict = {}  # dictionary of pynmrstar saveframes

        self.__cachedDictForStarAtom = {}

        self.num_of_dim = -1
        self.acq_dim_id = 1
        self.cur_spectral_dim = {}
        self.spectral_dim = {}
        self.spectral_dim_transfer = {}

        self.atom_type_history = {}
        self.onebond_idx_history = {}
        self.jcoupling_idx_history = {}
        self.relayed_idx_history = {}

        self.listIdInternal = {}
        self.chainNumberDict = {}
        self.extResKey = []

        self.polySeqRst = []
        self.polySeqRstFailed = []
        self.polySeqRstFailedAmbig = []
        self.compIdMap = {}

        self.f = []

        self.__allow_ext_seq = False

    @property
    def debug(self) -> bool:
        """ Retrieve debug mode.
        """

        return self.__debug

    @debug.setter
    def debug(self, debug: bool) -> None:
        self.__debug = debug

    @property
    def verbose_debug(self) -> bool:
        """ Retrieve verbose debug mode.
        """

        return self.__verbose_debug

    @verbose_debug.setter
    def verbose_debug(self, verbose_debug: bool) -> None:
        self.__verbose_debug = verbose_debug

    @property
    def internal(self) -> bool:
        """ Retrieve internal mode.
        """

        return self.__internal

    @internal.setter
    def internal(self, internal: bool) -> None:
        self.__internal = internal

    @property
    def createSfDict(self) -> bool:
        """ Whether to create saveframe dictionary.
        """

        return self.__createSfDict

    @createSfDict.setter
    def createSfDict(self, createSfDict: bool) -> None:
        self.__createSfDict = createSfDict

    @property
    def enforcePeakRowFormat(self) -> bool:
        """ Whether to enforce peak row format or not.
        """

        return self.__enforcePeakRowFormat

    @enforcePeakRowFormat.setter
    def enforcePeakRowFormat(self, enforcePeakRowFormat: bool) -> None:
        self.__enforcePeakRowFormat = enforcePeakRowFormat

    @property
    def originalFileName(self) -> str:
        """ Retrieve the original file name.
        """

        return self.__originalFileName

    @originalFileName.setter
    def originalFileName(self, originalFileName: str) -> None:
        self.__originalFileName = originalFileName

    @property
    def listIdCounter(self) -> dict:
        """ Retrieve list ID counter dictionary.
        """

        return self.__listIdCounter

    @listIdCounter.setter
    def listIdCounter(self, listIdCounter: dict) -> None:
        self.__listIdCounter = listIdCounter

    @property
    def reservedListIds(self) -> dict:
        """ Retrieve reserved list IDs.
        """

        return self.__reservedListIds

    @reservedListIds.setter
    def reservedListIds(self, reservedListIds: dict) -> None:
        self.__reservedListIds = reservedListIds

    @property
    def entryId(self) -> str:
        """ Retrieve entry ID.
        """

        return self.__entryId

    @entryId.setter
    def entryId(self, entryId: str) -> None:
        self.__entryId = entryId

    @property
    def csLoops(self) -> List[dict]:
        """ Retrieve CS loops as reference for peak assignment.
        """

        return self.__csLoops

    @csLoops.setter
    def csLoops(self, csLoops: List[dict]) -> None:
        self.__csLoops = csLoops

        if self.__csLoops is None or len(self.__csLoops) == 0:
            return

        segIds = set()
        for lp in self.__csLoops:
            for row in lp['data']:
                if 'Auth_asym_ID' in row:
                    seg_id = row['Auth_asym_ID']
                    if seg_id in EMPTY_VALUE:
                        continue
                    segIds.add(seg_id)
        if len(segIds) > 0:
            if len(self.polySeq) > 1:
                for ps1, ps2 in itertools.combinations(self.polySeq, 2):
                    if 'identical_auth_chain_id' in ps1 and ps2['auth_chain_id'] in ps1['identical_auth_chain_id']:
                        continue
                    if len(set(ps1['auth_seq_id']) & set(ps2['auth_seq_id'])) != len(ps1['auth_seq_id']) + len(ps2['auth_seq_id']):
                        return
            self.__defaultSegId__ = collections.Counter(segIds).most_common()[0][0]

    def exit(self, spectrum_names: Optional[dict] = None) -> None:
        """ Common function to exit a parse tree.
        """

        self.fillPkAuxLoops(spectrum_names)

        try:

            if self.hasPolySeq and self.polySeqRst is not None:
                sortPolySeqRst(self.polySeqRst,
                               None if self.reasons is None else self.reasons.get('non_poly_remap'))

                self.seqAlign, _ = alignPolymerSequence(self.pA, self.polySeq, self.polySeqRst,
                                                        resolvedMultimer=self.reasons is not None)
                self.chainAssign, message = assignPolymerSequence(self.pA, self.ccU, self.file_type,
                                                                  self.polySeq, self.polySeqRst, self.seqAlign)

                # allow conflict to detect sequence mismatch with sequence number shift (8dhz)
                if len(self.seqAlign) == 0 and not self.hasNonPolySeq:
                    for c in range(1, MAX_CONFLICT_ATTEMPT):
                        self.seqAlign, _ = alignPolymerSequenceWithConflicts(self.pA, self.polySeq, self.polySeqRst, c)
                        if len(self.seqAlign) > 0:
                            self.chainAssign, message = assignPolymerSequence(self.pA, self.ccU, self.file_type,
                                                                              self.polySeq, self.polySeqRst, self.seqAlign)
                            break

                if len(message) > 0:
                    self.f.extend(message)

                if self.chainAssign is not None:

                    if len(self.polySeq) == len(self.polySeqRst):

                        chain_mapping = {}

                        for ca in self.chainAssign:
                            ref_chain_id = ca['ref_chain_id']
                            test_chain_id = ca['test_chain_id']

                            if ref_chain_id != test_chain_id:
                                chain_mapping[test_chain_id] = ref_chain_id

                        if len(chain_mapping) == len(self.polySeq):

                            for ps in self.polySeqRst:
                                if ps['chain_id'] in chain_mapping:
                                    ps['chain_id'] = chain_mapping[ps['chain_id']]

                            self.seqAlign, _ = alignPolymerSequence(self.pA, self.polySeq, self.polySeqRst,
                                                                    resolvedMultimer=self.reasons is not None)
                            self.chainAssign, _ = assignPolymerSequence(self.pA, self.ccU, self.file_type,
                                                                        self.polySeq, self.polySeqRst, self.seqAlign)

                    trimSequenceAlignment(self.seqAlign, self.chainAssign)

                    if self.reasons is None and any(True for f in self.f
                                                    if '[Atom not found]' in f or '[Sequence mismatch]' in f
                                                    or '[Invalid atom nomenclature]' in f):

                        seqIdRemap = []

                        cyclicPolymer = {}

                        for ca in self.chainAssign:
                            ref_chain_id = ca['ref_chain_id']
                            test_chain_id = ca['test_chain_id']

                            sa = next(sa for sa in self.seqAlign
                                      if sa['ref_chain_id'] == ref_chain_id
                                      and sa['test_chain_id'] == test_chain_id)

                            poly_seq_model = next(ps for ps in self.polySeq
                                                  if ps['auth_chain_id'] == ref_chain_id)
                            poly_seq_rst = next(ps for ps in self.polySeqRst
                                                if ps['chain_id'] == test_chain_id)

                            seq_id_mapping = {}
                            offset = None
                            for ref_seq_id, mid_code, test_seq_id in zip(sa['ref_seq_id'], sa['mid_code'], sa['test_seq_id']):
                                if test_seq_id is None:
                                    continue
                                if mid_code == '|':
                                    try:
                                        seq_id_mapping[test_seq_id] =\
                                            next(auth_seq_id for auth_seq_id, seq_id
                                                 in zip(poly_seq_model['auth_seq_id'], poly_seq_model['seq_id'])
                                                 if seq_id == ref_seq_id and isinstance(auth_seq_id, int))
                                        if offset is None:
                                            offset = seq_id_mapping[test_seq_id] - test_seq_id
                                    except StopIteration:
                                        pass
                                elif mid_code == ' ' and test_seq_id in poly_seq_rst['seq_id']:
                                    idx = poly_seq_rst['seq_id'].index(test_seq_id)
                                    if poly_seq_rst['comp_id'][idx] == '.' and poly_seq_rst['auth_comp_id'][idx] not in EMPTY_VALUE:
                                        try:
                                            seq_id_mapping[test_seq_id] =\
                                                next(auth_seq_id for auth_seq_id, seq_id
                                                     in zip(poly_seq_model['auth_seq_id'], poly_seq_model['seq_id'])
                                                     if seq_id == ref_seq_id and isinstance(auth_seq_id, int))
                                        except StopIteration:
                                            pass

                            if offset is not None and all(v - k == offset for k, v in seq_id_mapping.items()):
                                test_seq_id_list = list(seq_id_mapping.keys())
                                min_test_seq_id = min(test_seq_id_list)
                                max_test_seq_id = max(test_seq_id_list)
                                for test_seq_id in range(min_test_seq_id + 1, max_test_seq_id):
                                    if test_seq_id not in seq_id_mapping:
                                        seq_id_mapping[test_seq_id] = test_seq_id + offset

                            if ref_chain_id not in cyclicPolymer:
                                cyclicPolymer[ref_chain_id] =\
                                    isCyclicPolymer(self.cR, self.polySeq, ref_chain_id,
                                                    self.representativeModelId, self.representativeAltId, self.modelNumName)

                            if cyclicPolymer[ref_chain_id]:

                                poly_seq_model = next(ps for ps in self.polySeq
                                                      if ps['auth_chain_id'] == ref_chain_id)

                                offset = None
                                for seq_id, comp_id in zip(poly_seq_rst['seq_id'], poly_seq_rst['comp_id']):
                                    if seq_id is not None and seq_id not in seq_id_mapping:
                                        _seq_id = next((_seq_id for _seq_id, _comp_id
                                                        in zip(poly_seq_model['seq_id'], poly_seq_model['comp_id'])
                                                        if _seq_id not in seq_id_mapping.values() and _comp_id == comp_id), None)
                                        if _seq_id is not None:
                                            offset = seq_id - _seq_id
                                            break

                                if offset is not None:
                                    for seq_id in poly_seq_rst['seq_id']:
                                        if seq_id is not None and seq_id not in seq_id_mapping:
                                            seq_id_mapping[seq_id] = seq_id - offset

                            if any(True for k, v in seq_id_mapping.items() if k != v)\
                               and not any(True for k, v in seq_id_mapping.items()
                                           if v in poly_seq_model['seq_id']
                                           and k == poly_seq_model['auth_seq_id'][poly_seq_model['seq_id'].index(v)]):
                                offsets = [v - k for k, v in seq_id_mapping.items()]
                                offsets = collections.Counter(offsets).most_common()
                                if len(offsets) == 1:
                                    offset = offsets[0][0]
                                    seq_id_mapping = {ref_seq_id - offset: ref_seq_id
                                                      for ref_seq_id in poly_seq_model['auth_seq_id']}
                                item = {'chain_id': test_chain_id, 'seq_id_dict': seq_id_mapping}
                                if item not in seqIdRemap:
                                    seqIdRemap.append(item)

                        if len(seqIdRemap) > 0:
                            if 'seq_id_remap' not in self.reasonsForReParsing:
                                self.reasonsForReParsing['seq_id_remap'] = seqIdRemap

                        if any(True for ps in self.polySeq if 'identical_chain_id' in ps):
                            polySeqRst, chainIdMapping =\
                                splitPolySeqRstForMultimers(self.pA, self.polySeq, self.polySeqRst,
                                                            self.chainAssign)

                            if polySeqRst is not None and (not self.hasNonPoly or len(self.polySeq) // len(self.nonPoly) in (1, 2)):
                                self.polySeqRst = polySeqRst
                                if 'chain_id_remap' not in self.reasonsForReParsing:
                                    self.reasonsForReParsing['chain_id_remap'] = chainIdMapping

                        if len(self.polySeq) == 1 and len(self.polySeqRst) == 1:
                            polySeqRst, chainIdMapping, _ =\
                                splitPolySeqRstForExactNoes(self.pA, self.polySeq, self.polySeqRst,
                                                            self.chainAssign)

                            if polySeqRst is not None:
                                self.polySeqRst = polySeqRst
                                if 'chain_id_clone' not in self.reasonsForReParsing:
                                    self.reasonsForReParsing['chain_id_clone'] = chainIdMapping

                        if self.hasNonPoly:
                            polySeqRst, nonPolyMapping =\
                                splitPolySeqRstForNonPoly(self.ccU, self.nonPoly, self.polySeqRst,
                                                          self.seqAlign, self.chainAssign)

                            if polySeqRst is not None:
                                self.polySeqRst = polySeqRst
                                if 'non_poly_remap' not in self.reasonsForReParsing:
                                    self.reasonsForReParsing['non_poly_remap'] = nonPolyMapping
                                else:
                                    for k, v in nonPolyMapping.items():
                                        if k not in self.reasonsForReParsing['non_poly_remap']:
                                            self.reasonsForReParsing['non_poly_remap'][k] = v
                                        else:
                                            for k2, v2 in v.items():
                                                if k2 not in self.reasonsForReParsing['non_poly_remap'][k]:
                                                    self.reasonsForReParsing['non_poly_remap'][k][k2] = v2

                        if self.hasBranched:
                            polySeqRst, branchedMapping =\
                                splitPolySeqRstForBranched(self.pA, self.polySeq, self.branched, self.polySeqRst,
                                                           self.chainAssign)

                            if polySeqRst is not None:
                                self.polySeqRst = polySeqRst
                                if 'branched_remap' not in self.reasonsForReParsing:
                                    self.reasonsForReParsing['branched_remap'] = branchedMapping

                        if len(self.polySeqRstFailed) > 0:
                            sortPolySeqRst(self.polySeqRstFailed)
                            if not any(True for f in self.f if '[Sequence mismatch]' in f):  # 2n6y
                                syncCompIdOfPolySeqRst(self.polySeqRstFailed, self.compIdMap)  # 2mx9

                            seqAlignFailed, _ = alignPolymerSequence(self.pA, self.polySeq, self.polySeqRstFailed)
                            chainAssignFailed, _ = assignPolymerSequence(self.pA, self.ccU, self.file_type,
                                                                         self.polySeq, self.polySeqRstFailed, seqAlignFailed)

                            if chainAssignFailed is not None:
                                seqIdRemapFailed = []

                                uniq_ps = not any(True for ps in self.polySeq if 'identical_chain_id' in ps)

                                for ca in chainAssignFailed:
                                    if ca['conflict'] > 0:
                                        continue
                                    ref_chain_id = ca['ref_chain_id']
                                    test_chain_id = ca['test_chain_id']

                                    sa = next((sa for sa in seqAlignFailed
                                               if sa['ref_chain_id'] == ref_chain_id
                                               and sa['test_chain_id'] == test_chain_id), None)

                                    if sa is None:
                                        continue

                                    poly_seq_model = next(ps for ps in self.polySeq
                                                          if ps['auth_chain_id'] == ref_chain_id)

                                    seq_id_mapping = {}
                                    for ref_seq_id, mid_code, test_seq_id\
                                            in zip(sa['ref_seq_id'], sa['mid_code'], sa['test_seq_id']):
                                        if test_seq_id is None:
                                            continue
                                        if mid_code == '|':
                                            try:
                                                seq_id_mapping[test_seq_id] =\
                                                    next(auth_seq_id for auth_seq_id, seq_id
                                                         in zip(poly_seq_model['auth_seq_id'], poly_seq_model['seq_id'])
                                                         if seq_id == ref_seq_id and isinstance(auth_seq_id, int))
                                            except StopIteration:
                                                if uniq_ps:
                                                    seq_id_mapping[test_seq_id] = ref_seq_id

                                    offset = None
                                    offsets = [v - k for k, v in seq_id_mapping.items()]
                                    if len(offsets) > 0\
                                       and ('gap_in_auth_seq' not in poly_seq_model or not poly_seq_model['gap_in_auth_seq']):
                                        offsets = collections.Counter(offsets).most_common()
                                        if len(offsets) > 1:
                                            offset = offsets[0][0]
                                            for k, v in seq_id_mapping.items():
                                                if v - k != offset:
                                                    seq_id_mapping[k] = k + offset

                                    if uniq_ps and offset is not None and len(seq_id_mapping) > 0\
                                       and ('gap_in_auth_seq' not in poly_seq_model or not poly_seq_model['gap_in_auth_seq']):
                                        for ref_seq_id, mid_code, test_seq_id, ref_code, test_code\
                                                in zip(sa['ref_seq_id'], sa['mid_code'], sa['test_seq_id'],
                                                       sa['ref_code'], sa['test_code']):
                                            if test_seq_id is None:
                                                continue
                                            if mid_code == '|' and test_seq_id not in seq_id_mapping:
                                                seq_id_mapping[test_seq_id] = test_seq_id + offset
                                            elif ref_code != '.' and test_code == '.':
                                                seq_id_mapping[test_seq_id] = test_seq_id + offset

                                    if any(True for k, v in seq_id_mapping.items() if k != v)\
                                       and not any(True for k, v in seq_id_mapping.items()
                                                   if v in poly_seq_model['seq_id']
                                                   and k == poly_seq_model['auth_seq_id'][poly_seq_model['seq_id'].index(v)]):
                                        offsets = [v - k for k, v in seq_id_mapping.items()]
                                        offsets = collections.Counter(offsets).most_common()
                                        if len(offsets) == 1:
                                            offset = offsets[0][0]
                                            seq_id_mapping = {ref_seq_id - offset: ref_seq_id
                                                              for ref_seq_id in poly_seq_model['auth_seq_id']}
                                        item = {'chain_id': ref_chain_id,
                                                'seq_id_dict': seq_id_mapping,
                                                'comp_id_set': list(set(poly_seq_model['comp_id']))}
                                        if item not in seqIdRemap:
                                            seqIdRemapFailed.append(item)

                                if len(seqIdRemapFailed) > 0:
                                    if 'chain_seq_id_remap' not in self.reasonsForReParsing:
                                        seqIdRemap =\
                                            self.reasonsForReParsing['seq_id_remap'] if 'seq_id_remap' in self.reasonsForReParsing\
                                            else []
                                        if len(seqIdRemap) != len(seqIdRemapFailed)\
                                           or seqIdRemap[0]['chain_id'] != seqIdRemapFailed[0]['chain_id']\
                                           or not all(src_seq_id in seqIdRemap[0] for src_seq_id in seqIdRemapFailed[0]):
                                            self.reasonsForReParsing['chain_seq_id_remap'] = seqIdRemapFailed

                                else:
                                    for ps in self.polySeqRstFailed:
                                        for ca in self.chainAssign:
                                            ref_chain_id = ca['ref_chain_id']
                                            test_chain_id = ca['test_chain_id']

                                            if test_chain_id != ps['chain_id']:
                                                continue

                                            sa = next(sa for sa in self.seqAlign
                                                      if sa['ref_chain_id'] == ref_chain_id
                                                      and sa['test_chain_id'] == test_chain_id)

                                            if len(sa['test_seq_id']) != len(sa['ref_seq_id']):
                                                continue

                                            poly_seq_model = next(ps for ps in self.polySeq
                                                                  if ps['auth_chain_id'] == ref_chain_id)

                                            seq_id_mapping, comp_id_mapping = {}, {}

                                            for seq_id, comp_id in zip(ps['seq_id'], ps['comp_id']):
                                                if seq_id in sa['test_seq_id']:
                                                    idx = sa['test_seq_id'].index(seq_id)
                                                    auth_seq_id = sa['ref_seq_id'][idx]
                                                    seq_id_mapping[seq_id] = auth_seq_id
                                                    comp_id_mapping[seq_id] = comp_id
                                            if any(True for k, v in seq_id_mapping.items() if k != v)\
                                               or ('label_seq_scheme' not in self.reasonsForReParsing
                                                   and all(v not in poly_seq_model['auth_seq_id']
                                                           for v in seq_id_mapping.values())):
                                                seqIdRemapFailed.append({'chain_id': ref_chain_id, 'seq_id_dict': seq_id_mapping,
                                                                         'comp_id_dict': comp_id_mapping})

                                    if len(seqIdRemapFailed) > 0:
                                        if 'ext_chain_seq_id_remap' not in self.reasonsForReParsing:
                                            seqIdRemap = []
                                            if 'seq_id_remap' in self.reasonsForReParsing:
                                                seqIdRemap = self.reasonsForReParsing['seq_id_remap']
                                            if len(seqIdRemap) != len(seqIdRemapFailed)\
                                               or seqIdRemap[0]['chain_id'] != seqIdRemapFailed[0]['chain_id']\
                                               or not all(src_seq_id in seqIdRemap[0] for src_seq_id in seqIdRemapFailed[0]):
                                                self.reasonsForReParsing['ext_chain_seq_id_remap'] = seqIdRemapFailed

            if 'local_seq_scheme' in self.reasonsForReParsing:
                if 'non_poly_remap' in self.reasonsForReParsing or 'branched_remap' in self.reasonsForReParsing:
                    del self.reasonsForReParsing['local_seq_scheme']
                elif 'seq_id_remap' in self.reasonsForReParsing:
                    del self.reasonsForReParsing['local_seq_scheme']
                elif 'chain_seq_id_remap' in self.reasonsForReParsing:
                    del self.reasonsForReParsing['local_seq_scheme']
                elif 'ext_chain_seq_id_remap' in self.reasonsForReParsing:
                    del self.reasonsForReParsing['local_seq_scheme']

            if 'local_seq_scheme' in self.reasonsForReParsing and len(self.reasonsForReParsing) == 1:
                sortPolySeqRst(self.polySeqRstFailed)
                if len(self.polySeqRstFailed) > 0:
                    self.reasonsForReParsing['extend_seq_scheme'] = self.polySeqRstFailed
                del self.reasonsForReParsing['local_seq_scheme']

        finally:
            self.warningMessage = sorted(list(set(self.f)), key=self.f.index)
            if all('(list_id=1, ' in f for f in self.warningMessage):
                self.warningMessage = [f.replace('(list_id=1, ', '(') for f in self.warningMessage]

            self.translateToStdResNameWrapper.cache_clear()
            self.__getCoordAtomSiteOf.cache_clear()

            translateToStdAtomNameNoRef.cache_clear()
            translateToStdAtomNameWithRef.cache_clear()

    def initSpectralDim(self) -> None:
        """ Initialize spectral_dim dictionary.
        """

        if self.num_of_dim not in (2, 3, 4):
            return
        self.cur_subtype = f'peak{self.num_of_dim}d'
        if self.num_of_dim not in self.listIdInternal:
            self.listIdInternal[self.num_of_dim] = 0
        self.listIdInternal[self.num_of_dim] += 1
        self.cur_list_id = self.listIdInternal[self.num_of_dim]
        if self.num_of_dim not in self.spectral_dim:
            self.spectral_dim[self.num_of_dim] = {}
        if self.cur_list_id not in self.spectral_dim[self.num_of_dim]:
            self.spectral_dim[self.num_of_dim][self.cur_list_id] = {}
        for _dim_id in range(1, self.num_of_dim + 1):
            self.spectral_dim[self.num_of_dim][self.cur_list_id][_dim_id] =\
                copy.copy(SPECTRAL_DIM_TEMPLATE
                          if len(self.cur_spectral_dim) == 0
                          or _dim_id not in self.cur_spectral_dim
                          else self.cur_spectral_dim[_dim_id])
            self.spectral_dim[self.num_of_dim][self.cur_list_id][_dim_id]['freq_hint'] = []  # list -> numpy array before exit()
            if self.file_type == 'nm-pea-pip':
                self.spectral_dim[self.num_of_dim][self.cur_list_id][_dim_id]['obs_freq_hint'] = []
        if self.num_of_dim not in self.spectral_dim_transfer:
            self.spectral_dim_transfer[self.num_of_dim] = {}
        if self.cur_list_id not in self.spectral_dim_transfer[self.num_of_dim]:
            self.spectral_dim_transfer[self.num_of_dim][self.cur_list_id] = []
        if self.num_of_dim not in self.atom_type_history:
            self.atom_type_history[self.num_of_dim] = {}
        if self.cur_list_id not in self.atom_type_history[self.num_of_dim]:
            self.atom_type_history[self.num_of_dim][self.cur_list_id] = []
        if self.num_of_dim not in self.onebond_idx_history:
            self.onebond_idx_history[self.num_of_dim] = {}
        if self.cur_list_id not in self.onebond_idx_history[self.num_of_dim]:
            self.onebond_idx_history[self.num_of_dim][self.cur_list_id] = []
        if self.num_of_dim not in self.jcoupling_idx_history:
            self.jcoupling_idx_history[self.num_of_dim] = {}
        if self.cur_list_id not in self.jcoupling_idx_history[self.num_of_dim]:
            self.jcoupling_idx_history[self.num_of_dim][self.cur_list_id] = -1
        if self.num_of_dim not in self.relayed_idx_history:
            self.relayed_idx_history[self.num_of_dim] = {}
        if self.cur_list_id not in self.relayed_idx_history[self.num_of_dim]:
            self.relayed_idx_history[self.num_of_dim][self.cur_list_id] = -1
        if self.num_of_dim == 2:
            self.peaks2D = 0
        if self.num_of_dim == 3:
            self.peaks3D = 0
        if self.num_of_dim == 4:
            self.peaks4D = 0
        self.use_peak_row_format = True

    def testAssignment(self, _dim_id: int, _assign: List[dict], _label: str
                       ) -> Tuple[bool, Optional[str]]:
        """ Perform consistency test between peak assignment and spectral_dim.
        """

        cur_spectral_dim = self.spectral_dim[self.num_of_dim][self.cur_list_id][_dim_id]
        if cur_spectral_dim['atom_type'] is not None:
            if len(_assign) > 0 and 'atom_id' in _assign[0]:
                atom_type = cur_spectral_dim['atom_type']
                _atom_type = _assign[0]['atom_id'][0]
                if atom_type == _atom_type:
                    return True, None
                concat_type = _atom_type + atom_type
                if concat_type in _label:
                    _label = _label.replace(concat_type, atom_type)
                    return False, _label
                return False, None
        return True, None

    def validateAtomType(self, _dim_id: int, atom_type: str, position: str) -> bool:
        """ Validate peak position based on atom type.
        """

        if self.reasons is None and self.software_name == 'PIPP':
            return True

        cur_spectral_dim = self.spectral_dim[self.num_of_dim][self.cur_list_id][_dim_id]

        def predict_spectral_region_from_position():
            _position = float(position)

            if C_AROMATIC_CENTER_MIN_TOR < _position <= C_AROMATIC_CENTER_MAX and atom_type == 'C':
                return 'C-aromatic'
            if N_AMIDE_CENTER_MIN < _position <= N_AMIDE_CENTER_MAX:
                if atom_type == 'N':
                    return 'N'
                if atom_type == 'C':
                    return 'C-aromatic'
            if C_CARBONYL_CENTER_MIN <= _position <= C_CARBONYL_CENTER_MAX and atom_type == 'C':
                return 'CO'
            if HN_AROMATIC_CENTER_MIN < _position <= HN_AROMATIC_CENTER_MAX and atom_type == 'H':
                return 'HN/H-aromatic'
            if H_ALL_CENTER_MIN < _position <= H_ALL_CENTER_MAX and atom_type == 'H':
                return 'H'  # all
            if H_IMIDE_CENTER_MIN < _position <= H_IMIDE_CENTER_MAX and atom_type == 'H':
                return 'H-imide'
            max_ppm = max(cur_spectral_dim['freq_hint']) if len(cur_spectral_dim['freq_hint']) > 0 else _position
            if H_ALIPHATIC_CENTER_MIN < _position <= H_ALIPHATIC_CENTER_MAX and atom_type == 'H':
                return 'H-aliphatic' if max_ppm < 7 else 'H'
            if H_METHYL_CENTER_MIN < _position <= H_METHYL_CENTER_MAX and atom_type == 'H':
                return 'H-methyl' if max_ppm < 3 else 'H-aliphatic'
            if C_ALL_CENTER_MIN < _position <= C_ALL_CENTER_MAX and atom_type == 'C':
                return 'C'  # all
            if C_ALIPHATIC_CENTER_MIN < _position <= C_ALIPHATIC_CENTER_MAX and atom_type == 'C':
                return 'C-aliphatic'
            if C_METHYL_CENTER_MIN < _position <= C_METHYL_CENTER_MAX and atom_type == 'C':
                return 'C-methyl'
            return ''

        if self.reasons is not None and 'atom_type_history' in self.reasons:
            atom_types = self.reasons['atom_type_history'][self.num_of_dim][self.cur_list_id]
            if len(atom_types) == self.num_of_dim:
                _atom_type = atom_types[_dim_id - 1]
                if atom_type in PROTON_BEGIN_CODE:
                    return _atom_type == 'H'
                return _atom_type == atom_type
        if cur_spectral_dim['atom_type'] is not None:
            if cur_spectral_dim['atom_type'] == atom_type:
                cur_spectral_dim['fixed'] = True  # be robust against interference of unreliable assignments (bmr36675)
                return True
            if 'fixed' in cur_spectral_dim:  # XEASY INNAME label is not reliable (2kj5)
                return False
            if atom_type in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS:
                _atom_type = predict_spectral_region_from_position()
                if len(_atom_type) == 0:
                    if self.num_of_dim == 2 and self.peaks2D < 20:
                        return True
                    if self.num_of_dim == 3 and self.peaks3D < 20:
                        return True
                    if self.num_of_dim == 4 and self.peaks4D < 20:
                        return True
                    return False
                if not _atom_type.startswith(atom_type):
                    return False
                cur_spectral_dim['fixed'] = True
                cur_spectral_dim['atom_type'] = None
                cur_spectral_dim['axis_code'] = None
                cur_spectral_dim['atom_isotope_number'] = None
                return True
        if atom_type in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS:
            _atom_type = predict_spectral_region_from_position()
            if len(_atom_type) == 0:
                if self.num_of_dim == 2 and self.peaks2D < 20:
                    return True
                if self.num_of_dim == 3 and self.peaks3D < 20:
                    return True
                if self.num_of_dim == 4 and self.peaks4D < 20:
                    return True
                return False
            if not _atom_type.startswith(atom_type):
                return False
            cur_spectral_dim['atom_type'] = atom_type
            cur_spectral_dim['axis_code'] = atom_type
            cur_spectral_dim['atom_isotope_number'] = ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS[atom_type][0]
            return True
        return False

    def fillPkAuxLoops(self, spectrum_names: Optional[dict]) -> None:
        """ Create auxiliary loops of spectral peak list.
        """

        if len(self.spectral_dim) > 0:
            dim_to_code = {1: 'x', 2: 'y', 3: 'z', 4: 'a'}

            for d, v in self.spectral_dim.items():
                for _id, _v in v.items():

                    file_name = self.__originalFileName.lower()
                    alt_file_name = ''
                    try:
                        if spectrum_names is not None:
                            alt_file_name = spectrum_names[d][_id].lower()
                    except (KeyError, AttributeError):
                        pass

                    _file_names = (file_name, alt_file_name)

                    self.acq_dim_id = 1

                    for __d, __v in _v.items():
                        if 'freq_hint' in __v:
                            __v['freq_hint'] = numpy.array(__v['freq_hint'], dtype=float)  # list -> numpy array
                            if __v['freq_hint'].size > 0:
                                center = numpy.mean(__v['freq_hint'])
                                max_ppm = __v['freq_hint'].max()
                                min_ppm = __v['freq_hint'].min()

                                if __v['atom_isotope_number'] is None:
                                    if (C_AROMATIC_CENTER_MIN_TOR if any('aro' in n for n in _file_names)
                                        or any('anoe' in n for n in _file_names) else C_AROMATIC_CENTER_MIN)\
                                       < center <= C_AROMATIC_CENTER_MAX and min_ppm > C_ALL_CENTER_MAX:
                                        __v['atom_type'] = 'C'
                                        __v['atom_isotope_number'] = 13
                                        __v['axis_code'] = 'C-aromatic'
                                    elif N_AMIDE_CENTER_MIN < center <= N_AMIDE_CENTER_MAX\
                                            and min_ppm > C_ALL_CENTER_MIN and max_ppm < C_CARBONYL_CENTER_MIN:
                                        __v['atom_type'] = 'N'
                                        __v['atom_isotope_number'] = 15
                                        __v['axis_code'] = 'N'
                                    elif C_CARBONYL_CENTER_MIN <= center <= C_CARBONYL_CENTER_MAX:
                                        __v['atom_type'] = 'C'
                                        __v['atom_isotope_number'] = 13
                                        __v['axis_code'] = 'CO'
                                    elif HN_AROMATIC_CENTER_MIN < center <= HN_AROMATIC_CENTER_MAX:
                                        __v['atom_type'] = 'H'
                                        __v['atom_isotope_number'] = 1
                                        __v['axis_code'] = 'HN/H-aromatic'
                                    elif H_ALL_CENTER_MIN < center <= H_ALL_CENTER_MAX:
                                        __v['atom_type'] = 'H'
                                        __v['atom_isotope_number'] = 1
                                        __v['axis_code'] = 'H'  # all
                                    elif H_ALIPHATIC_CENTER_MIN < center <= H_ALIPHATIC_CENTER_MAX:
                                        __v['atom_type'] = 'H'
                                        __v['atom_isotope_number'] = 1
                                        __v['axis_code'] = 'H-aliphatic' if max_ppm < 7 else 'H'
                                    elif H_METHYL_CENTER_MIN < center <= H_METHYL_CENTER_MAX:
                                        __v['atom_type'] = 'H'
                                        __v['atom_isotope_number'] = 1
                                        __v['axis_code'] = 'H-methyl' if max_ppm < 3 else 'H-aliphatic'
                                    elif C_ALL_CENTER_MIN < center <= C_ALL_CENTER_MAX:
                                        __v['atom_type'] = 'C'
                                        __v['atom_isotope_number'] = 13
                                        __v['axis_code'] = 'C'  # all
                                    elif C_ALIPHATIC_CENTER_MIN < center <= C_ALIPHATIC_CENTER_MAX:
                                        __v['atom_type'] = 'C'
                                        __v['atom_isotope_number'] = 13
                                        __v['axis_code'] = 'C-aliphatic'
                                    elif C_METHYL_CENTER_MIN < center <= C_METHYL_CENTER_MAX:
                                        __v['atom_type'] = 'C'
                                        __v['atom_isotope_number'] = 13
                                        __v['axis_code'] = 'C-methyl'
                                    elif C_METHYL_CENTER_MIN < min_ppm and max_ppm < C_CARBONYL_CENTER_MAX_TOR:
                                        __v['atom_type'] = 'C'
                                        __v['atom_isotope_number'] = 13
                                        __v['axis_code'] = 'C'  # all

                                isotope_number = __v['atom_isotope_number']

                                if isotope_number is not None:
                                    __v['acquisition'] =\
                                        'yes' if __d == self.acq_dim_id\
                                        and (isotope_number == 1
                                             or (isotope_number == 13 and self.exptlMethod == 'SOLID-STATE NMR')) else 'no'

                                    if __d == 1 and __v['acquisition'] == 'no':
                                        self.acq_dim_id = d

                                    __v['under_sampling_type'] = 'not observed' if __v['acquisition'] == 'yes' else 'aliased'

                            if __v['spectral_region'] is None and __v['freq_hint'].size > 0:
                                atom_type = __v['atom_type']
                                if C_AROMATIC_CENTER_MIN_TOR < center <= C_AROMATIC_CENTER_MAX and min_ppm > C_ALL_CENTER_MAX\
                                   and atom_type == 'C':
                                    __v['spectral_region'] = 'C-aromatic'
                                elif N_AMIDE_CENTER_MIN < center <= N_AMIDE_CENTER_MAX and min_ppm > C_ALL_CENTER_MAX:
                                    if atom_type == 'N':
                                        __v['spectral_region'] = 'N'
                                    if atom_type == 'C':
                                        __v['spectral_region'] = 'C-aromatic'
                                elif C_CARBONYL_CENTER_MIN <= center <= C_CARBONYL_CENTER_MAX and atom_type == 'C':
                                    __v['spectral_region'] = 'CO'
                                elif HN_AROMATIC_CENTER_MIN < center <= HN_AROMATIC_CENTER_MAX and atom_type == 'H':
                                    __v['spectral_region'] = 'HN/H-aromatic'
                                elif H_ALL_CENTER_MIN < center <= H_ALL_CENTER_MAX and atom_type == 'H':
                                    __v['spectral_region'] = 'H'  # all
                                elif H_ALIPHATIC_CENTER_MIN < center <= H_ALIPHATIC_CENTER_MAX and atom_type == 'H':
                                    __v['spectral_region'] = 'H-aliphatic' if max_ppm < 7 else 'H'
                                elif H_METHYL_CENTER_MIN < center <= H_METHYL_CENTER_MAX and atom_type == 'H':
                                    __v['spectral_region'] = 'H-methyl' if max_ppm < 3 else 'H-aliphatic'
                                elif C_ALL_CENTER_MIN < center <= C_ALL_CENTER_MAX and atom_type == 'C':
                                    __v['spectral_region'] = 'C'  # all
                                elif C_ALIPHATIC_CENTER_MIN < center <= C_ALIPHATIC_CENTER_MAX and atom_type == 'C':
                                    __v['spectral_region'] = 'C-aliphatic'
                                elif C_METHYL_CENTER_MIN < center <= C_METHYL_CENTER_MAX and atom_type == 'C':
                                    __v['spectral_region'] = 'C-methyl'
                                elif __v['atom_isotope_number'] == 1:
                                    __v['spectral_region'] = 'H'
                                elif __v['atom_isotope_number'] == 13:
                                    __v['spectral_region'] = 'C'
                                elif __v['atom_isotope_number'] == 15:
                                    __v['spectral_region'] = 'N'
                                elif __v['atom_isotope_number'] == 19:
                                    __v['spectral_region'] = 'F'
                                elif __v['atom_isotope_number'] == 31:
                                    __v['spectral_region'] = 'P'
                                elif self.reasons is not None and 'atom_type_history' in self.reasons:
                                    atom_type = self.reasons['atom_type_history'][d][_id][__d - 1]
                                    if atom_type == 'H':
                                        __v['atom_type'] = __v['axis_code'] = __v['spectral_region'] = 'H'
                                        __v['atom_isotope_number'] = 1
                                    elif atom_type == 'C':
                                        __v['atom_type'] = __v['axis_code'] = __v['spectral_region'] = 'C'
                                        __v['atom_isotope_number'] = 13
                                    elif atom_type == 'N':
                                        __v['atom_type'] = __v['axis_code'] = __v['spectral_region'] = 'N'
                                        __v['atom_isotope_number'] = 15
                                    elif atom_type == 'F':
                                        __v['atom_type'] = __v['axis_code'] = __v['spectral_region'] = 'F'
                                        __v['atom_isotope_number'] = 19
                                    elif atom_type == 'P':
                                        __v['atom_type'] = __v['axis_code'] = __v['spectral_region'] = 'P'
                                        __v['atom_isotope_number'] = 31

                            if __v['freq_hint'].size > 0 and d > 2 and __d >= 2\
                               and self.exptlMethod != 'SOLID-STATE NMR' and __v['atom_isotope_number'] == 13:
                                min_ppm = __v['freq_hint'].min()
                                width = max_ppm - min_ppm
                                if center < 100.0 and width < 50.0:
                                    __v['under_sampling_type'] = 'folded'

                            if __v['spectrometer_frequency'] is None and 'obs_freq_hint' in __v and len(__v['obs_freq_hint']) > 0:
                                __v['spectrometer_frequency'] = collections.Counter(__v['obs_freq_hint']).most_common()[0][0]

                            if 'obs_freq_hint' in __v:
                                del __v['obs_freq_hint']

                            if __v['spectrometer_frequency'] is not None and __v['sweep_width_units'] == 'ppm':
                                row = [str(__v['sweep_width']), str(__v['spectrometer_frequency'])]
                                max_eff_digits = getMaxEffDigits(row)

                                __v['sweep_width'] = float(roundString(str(__v['sweep_width'] * __v['spectrometer_frequency']),
                                                                       max_eff_digits))
                                __v['sweep_width_units'] = 'Hz'

                    for __v in _v.values():
                        if __v['axis_code'] == 'HN/H-aromatic':
                            has_a = any(___v['spectral_region'] == 'C-aromatic' for ___v in _v.values())
                            __v['axis_code'] = 'H-aromatic' if has_a else 'H'
                        if __v['spectral_region'] == 'HN/H-aromatic':
                            has_a = any(___v['spectral_region'] == 'C-aromatic' for ___v in _v.values())
                            __v['spectral_region'] = 'H-aromatic' if has_a else 'H'
                            __v['_spectral_region'] = 'H-aromatic' if has_a else 'HN'
                        else:
                            __v['_spectral_region'] = __v['spectral_region']
                        self.atom_type_history[d][_id].append(__v['atom_type'])

            for d, v in self.spectral_dim.items():
                for _id, cur_spectral_dim in v.items():

                    if self.__debug:
                        print("original file name: "
                              f'{self.__originalFileName}{", spectrum name: " + str(spectrum_names[d][_id]) if spectrum_names is not None else ""}')  # noqa: E501, pylint: disable=line-too-long

                    file_name = self.__originalFileName.lower()
                    alt_file_name = ''
                    try:
                        if spectrum_names is not None:
                            alt_file_name = spectrum_names[d][_id].lower()
                    except (KeyError, AttributeError):
                        pass

                    _file_names = (file_name, alt_file_name)

                    cur_spectral_dim_transfer = self.spectral_dim_transfer[d][_id]

                    is_noesy = any('noe' in n for n in _file_names) or any('roe' in n for n in _file_names)
                    no_aromatic = all(_dict['_spectral_region'] != 'C-aromatic' for _dict in cur_spectral_dim.values())

                    # determine 'onebond' coherence transfer based on capitalized axis names
                    if 'axis_order' in cur_spectral_dim[1]:
                        upper_count = lower_count = 0
                        for _dict in cur_spectral_dim.values():
                            if _dict['axis_order'].isupper():
                                upper_count += 1
                            if _dict['axis_order'].islower():
                                lower_count += 1
                        if 0 not in (upper_count, lower_count) and 2 in (upper_count, lower_count):
                            for _dim_id1, _dict1 in cur_spectral_dim.items():
                                for _dim_id2, _dict2 in cur_spectral_dim.items():
                                    if _dim_id1 >= _dim_id2:
                                        continue
                                    if (upper_count == 2 and _dict1['axis_order'].isupper() and _dict2['axis_order'].isupper())\
                                       or (lower_count == 2 and _dict1['axis_order'].islower() and _dict2['axis_order'].islower()):
                                        transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                    'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                    'type': 'onebond',
                                                    'indirect': 'no'}
                                        if transfer in cur_spectral_dim_transfer:
                                            continue
                                        cur_spectral_dim_transfer.append(transfer)

                    # determine 'onebond' coherence transfer based on assigned chemical shifts
                    history = self.onebond_idx_history[d][_id]
                    if len(history) > 0:
                        onebond_idx = self.onebond_idx_history[d][_id] = collections.Counter(history).most_common()[0][0]
                        if d == 2:
                            _dim_id1, _dim_id2 = DIM_TRANSFER_PAT_2D[onebond_idx]
                            _dim_id1, _dim_id2 = _dim_id1 + 1, _dim_id2 + 1
                            if not any(True for _transfer in cur_spectral_dim_transfer
                                       if _transfer['type'] == 'onebond'
                                       and (_dim_id1 in [_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']]
                                            or _dim_id2 in [_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']])):
                                transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                            'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                            'type': 'onebond',
                                            'indirect': 'no'}
                                cur_spectral_dim_transfer.append(transfer)
                        elif d == 3:
                            _dim_id1, _dim_id2 = DIM_TRANSFER_PAT_3D[onebond_idx]
                            _dim_id1, _dim_id2 = _dim_id1 + 1, _dim_id2 + 1
                            if not any(True for _transfer in cur_spectral_dim_transfer
                                       if _transfer['type'] == 'onebond'
                                       and (_dim_id1 in [_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']]
                                            or _dim_id2 in [_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']])):
                                transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                            'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                            'type': 'onebond',
                                            'indirect': 'no'}
                                cur_spectral_dim_transfer.append(transfer)
                        elif d == 4:
                            (_dim_id1, _dim_id2), (_dim_id3, _dim_id4) = DIM_TRANSFER_PAT_4D[onebond_idx]
                            _dim_id1, _dim_id2, _dim_id3, _dim_id4 = _dim_id1 + 1, _dim_id2 + 1, _dim_id3 + 1, _dim_id4 + 1
                            if not any(True for _transfer in cur_spectral_dim_transfer
                                       if _transfer['type'] == 'onebond'
                                       and (_dim_id1 in [_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']]
                                            or _dim_id2 in [_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']])):
                                transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                            'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                            'type': 'onebond',
                                            'indirect': 'no'}
                            if not any(True for _transfer in cur_spectral_dim_transfer
                                       if _transfer['type'] == 'onebond'
                                       and (_dim_id3 in [_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']]
                                            or _dim_id4 in [_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']])):
                                transfer = {'spectral_dim_id_1': min([_dim_id3, _dim_id4]),
                                            'spectral_dim_id_2': max([_dim_id3, _dim_id4]),
                                            'type': 'onebond',
                                            'indirect': 'no'}
                                cur_spectral_dim_transfer.append(transfer)

                    # onebond: 'Any transfer that connects only directly bonded atoms in this experiment'
                    for _dim_id1, _dict1 in cur_spectral_dim.items():
                        _region1 = _dict1['_spectral_region']
                        if _region1 in ('HN', 'H-aliphatic', 'H-aromatic', 'H-methyl'):
                            cases = 0
                            max_corrcoef = 0.0
                            for _dim_id2, _dict2 in cur_spectral_dim.items():
                                _region2 = _dict2['_spectral_region']
                                if (_region1 == 'HN' and _region2 == 'N')\
                                   or ((_region1 == 'H-aliphatic' or (is_noesy and no_aromatic and _region1 == 'H'))
                                       and (_region2 == 'C-aliphatic' or (is_noesy and _region2 == 'C')))\
                                   or (_region1 == 'H-aromatic' and _region2 == 'C-aromatic')\
                                   or (_region1 == 'H-methyl' and _region2 == 'C-methyl'):
                                    if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                        if not any(True for _transfer in cur_spectral_dim_transfer
                                                   if _transfer['type'] == 'onebond'
                                                   and (_dim_id1 in [_transfer['spectral_dim_id_1'],
                                                                     _transfer['spectral_dim_id_2']]
                                                        or _dim_id2 in [_transfer['spectral_dim_id_1'],
                                                                        _transfer['spectral_dim_id_2']])):
                                            if _dict1['freq_hint'].size < 2\
                                               or numpy.max(_dict1['freq_hint']) == numpy.min(_dict1['freq_hint'])\
                                               or numpy.max(_dict2['freq_hint']) == numpy.min(_dict2['freq_hint']):
                                                continue
                                            _corrcoef = numpy.corrcoef(_dict1['freq_hint'], _dict2['freq_hint'])[0][1]
                                            if _corrcoef < 0.0 or (_corrcoef < MIN_CORRCOEF_FOR_ONE_BOND_TRANSFER and no_aromatic):
                                                continue
                                            cases += 1
                                            max_corrcoef = max(max_corrcoef, _corrcoef)

                                            if d == 3:
                                                for _dim_id3, _dict3 in cur_spectral_dim.items():
                                                    if _dim_id3 in (_dim_id1, _dim_id2):
                                                        continue
                                                    if _dict3['atom_type'] != _dict1['atom_type']:
                                                        continue
                                                    if not any(True for _transfer in cur_spectral_dim_transfer
                                                       if _transfer['type'] == 'onebond'
                                                       and (_dim_id3 in [_transfer['spectral_dim_id_1'],
                                                                         _transfer['spectral_dim_id_2']]
                                                            or _dim_id2 in [_transfer['spectral_dim_id_1'],
                                                                            _transfer['spectral_dim_id_2']])):
                                                        if _dict3['freq_hint'].size < 2\
                                                           or numpy.max(_dict3['freq_hint']) == numpy.min(_dict3['freq_hint'])\
                                                           or numpy.max(_dict2['freq_hint']) == numpy.min(_dict2['freq_hint']):
                                                            continue
                                                        _corrcoef = numpy.corrcoef(_dict3['freq_hint'], _dict2['freq_hint'])[0][1]
                                                        if _corrcoef < 0.0\
                                                           or (_corrcoef < MIN_CORRCOEF_FOR_ONE_BOND_TRANSFER and no_aromatic):
                                                            continue
                                                        cases += 1
                                                        max_corrcoef = max(max_corrcoef, _corrcoef)

                            if cases == 1:
                                for _dim_id2, _dict2 in cur_spectral_dim.items():
                                    _region2 = _dict2['_spectral_region']
                                    if (_region1 == 'HN' and _region2 == 'N')\
                                       or ((_region1 == 'H-aliphatic' or (is_noesy and no_aromatic and _region1 == 'H'))
                                       and (_region2 == 'C-aliphatic' or (is_noesy and _region2 == 'C')))\
                                       or (_region1 == 'H-aromatic' and _region2 == 'C-aromatic')\
                                       or (_region1 == 'H-methyl' and _region2 == 'C-methyl'):
                                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                            if not any(True for _transfer in cur_spectral_dim_transfer
                                                       if _transfer['type'] == 'onebond'
                                                       and (_dim_id1 in [_transfer['spectral_dim_id_1'],
                                                                         _transfer['spectral_dim_id_2']]
                                                            or _dim_id2 in [_transfer['spectral_dim_id_1'],
                                                                            _transfer['spectral_dim_id_2']])):
                                                transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                            'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                            'type': 'onebond',
                                                            'indirect': 'no'}
                                                if transfer in cur_spectral_dim_transfer:
                                                    continue
                                                cur_spectral_dim_transfer.append(transfer)

                            elif cases > 1:
                                for _dim_id2, _dict2 in cur_spectral_dim.items():
                                    _region2 = _dict2['_spectral_region']
                                    if (_region1 == 'HN' and _region2 == 'N')\
                                       or ((_region1 == 'H-aliphatic' or (is_noesy and no_aromatic and _region1 == 'H'))
                                       and (_region2 == 'C-aliphatic' or (is_noesy and _region2 == 'C')))\
                                       or (_region1 == 'H-aromatic' and _region2 == 'C-aromatic')\
                                       or (_region1 == 'H-methyl' and _region2 == 'C-methyl'):
                                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                            if not any(True for _transfer in cur_spectral_dim_transfer
                                                       if _transfer['type'] == 'onebond'
                                                       and (_dim_id1 in [_transfer['spectral_dim_id_1'],
                                                                         _transfer['spectral_dim_id_2']]
                                                            or _dim_id2 in [_transfer['spectral_dim_id_1'],
                                                                            _transfer['spectral_dim_id_2']])):
                                                if _dict1['freq_hint'].size < 2\
                                                   or numpy.max(_dict1['freq_hint']) == numpy.min(_dict1['freq_hint'])\
                                                   or numpy.max(_dict2['freq_hint']) == numpy.min(_dict2['freq_hint']):
                                                    continue
                                                _corrcoef = numpy.corrcoef(_dict1['freq_hint'], _dict2['freq_hint'])[0][1]
                                                if _corrcoef < 0.0 or (_corrcoef < max_corrcoef and no_aromatic):
                                                    continue
                                                transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                            'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                            'type': 'onebond',
                                                            'indirect': 'no'}
                                                if transfer in cur_spectral_dim_transfer:
                                                    continue
                                                cur_spectral_dim_transfer.append(transfer)

                    for _dim_id1, _dict1 in cur_spectral_dim.items():
                        _region1 = _dict1['_spectral_region']
                        if _region1 in ('HN', 'H-aliphatic', 'H-aromatic', 'H-methyl'):
                            for _dim_id2, _dict2 in cur_spectral_dim.items():
                                _region2 = _dict2['_spectral_region']
                                if (_region1 == 'HN' and _region2 == 'N')\
                                   or ((_region1 == 'H-aliphatic' or (is_noesy and no_aromatic and _region1 == 'H'))
                                       and (_region2 == 'C-aliphatic' or (is_noesy and _region2 == 'C')))\
                                   or (_region1 == 'H-aromatic' and _region2 == 'C-aromatic')\
                                   or (_region1 == 'H-methyl' and _region2 == 'C-methyl'):
                                    if _dict1['acquisition'] == 'no' and _dict2['acquisition'] == 'no':
                                        if not any(True for _transfer in cur_spectral_dim_transfer
                                                   if _transfer['type'] == 'onebond'
                                                   and (_dim_id1 in [_transfer['spectral_dim_id_1'],
                                                                     _transfer['spectral_dim_id_2']]
                                                        or _dim_id2 in [_transfer['spectral_dim_id_1'],
                                                                        _transfer['spectral_dim_id_2']])):
                                            if _dict1['freq_hint'].size < 2\
                                               or numpy.max(_dict1['freq_hint']) == numpy.min(_dict1['freq_hint'])\
                                               or numpy.max(_dict2['freq_hint']) == numpy.min(_dict2['freq_hint']):
                                                continue
                                            _corrcoef = numpy.corrcoef(_dict1['freq_hint'], _dict2['freq_hint'])[0][1]
                                            if _corrcoef < 0.0 or (_corrcoef < max_corrcoef and no_aromatic):
                                                continue
                                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                        'type': 'onebond',
                                                        'indirect': 'no'}
                                            if transfer in cur_spectral_dim_transfer:
                                                continue
                                            cur_spectral_dim_transfer.append(transfer)

                    # jcoupling: 'Transfer via direct J coupling over one or more bonds'
                    if any('cosy' in n for n in _file_names):
                        if d == 2:
                            for _dim_id1, _dict1 in cur_spectral_dim.items():
                                _iso_num1 = _dict1['atom_isotope_number']
                                if _iso_num1 == 1:
                                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                                        if _dim_id1 == _dim_id2 or _iso_num1 != _dict2['atom_isotope_number']:
                                            continue
                                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                        'type': 'jcoupling',
                                                        'indirect': 'no'}
                                        if transfer in cur_spectral_dim_transfer:
                                            continue
                                        cur_spectral_dim_transfer.append(transfer)

                        elif d == 3:
                            for _dim_id1, _dict1 in cur_spectral_dim.items():
                                _region1 = _dict1['_spectral_region']
                                if _region1 == 'HN':
                                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                                        if _dim_id1 == _dim_id2 or _dict2['atom_isotope_number'] not in (1, 13):
                                            continue
                                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                        'type': 'jcoupling',
                                                        'indirect': 'yes'}
                                            if transfer in cur_spectral_dim_transfer:
                                                continue
                                            cur_spectral_dim_transfer.append(transfer)
                            for _dim_id1, _dict1 in cur_spectral_dim.items():
                                _region1 = _dict1['_spectral_region']
                                if _region1 == 'H-aliphatic':
                                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                                        _isotope2 = _dict2['atom_isotope_number']
                                        if _dim_id1 == _dim_id2 or _isotope2 not in (1, 13):
                                            continue
                                        if _isotope2 == 13\
                                           and not any(True for _transfer in cur_spectral_dim_transfer
                                                       if _transfer['type'] == 'onebond'
                                                       and {_dim_id1, _dim_id2} == {_transfer['spectral_dim_id_1'],
                                                                                    _transfer['spectral_dim_id_2']}):
                                            if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                                transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                            'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                            'type': 'jcoupling',
                                                            'indirect': 'yes' if _isotope2 == 1 else 'no'}
                                                if transfer in cur_spectral_dim_transfer:
                                                    continue
                                                cur_spectral_dim_transfer.append(transfer)

                        elif d == 4:
                            for _dim_id1, _dict1 in cur_spectral_dim.items():
                                _region1 = _dict1['_spectral_region']
                                if _region1 == 'HN':
                                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                                        if _dim_id1 == _dim_id2 or _dict2['atom_isotope_number'] != 1:
                                            continue
                                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                        'type': 'jcoupling',
                                                        'indirect': 'yes'}
                                            if transfer in cur_spectral_dim_transfer:
                                                continue
                                            cur_spectral_dim_transfer.append(transfer)
                            for _dim_id1, _dict1 in cur_spectral_dim.items():
                                _region1 = _dict1['_spectral_region']
                                if _region1 == 'H-aliphatic':
                                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                                        if _dim_id1 == _dim_id2 or _dict2['atom_isotope_number'] != 1:
                                            continue
                                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                        'type': 'jcoupling',
                                                        'indirect': 'no'}
                                            if transfer in cur_spectral_dim_transfer:
                                                continue
                                            cur_spectral_dim_transfer.append(transfer)

                    # jmultibond: 'Transfer via direct J coupling over multiple bonds'

                    # relayed: 'Transfer via multiple successive J coupling steps (TOCSY relay)'
                    if any('tocsy' in n for n in _file_names):
                        if d == 2:
                            for _dim_id1, _dict1 in cur_spectral_dim.items():
                                _iso_num1 = _dict1['atom_isotope_number']
                                if _iso_num1 == 1:
                                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                                        if _dim_id1 == _dim_id2 or _iso_num1 != _dict2['atom_isotope_number']:
                                            continue
                                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                        'type': 'relayed',
                                                        'indirect': 'no'}
                                        if transfer in cur_spectral_dim_transfer:
                                            continue
                                        cur_spectral_dim_transfer.append(transfer)

                        elif d == 3:
                            for _dim_id1, _dict1 in cur_spectral_dim.items():
                                _region1 = _dict1['_spectral_region']
                                if _region1 == 'HN':
                                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                                        if _dim_id1 == _dim_id2 or _dict2['atom_isotope_number'] not in (1, 13):
                                            continue
                                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                        'type': 'relayed',
                                                        'indirect': 'yes'}
                                            if transfer in cur_spectral_dim_transfer:
                                                continue
                                            cur_spectral_dim_transfer.append(transfer)
                            for _dim_id1, _dict1 in cur_spectral_dim.items():
                                _region1 = _dict1['_spectral_region']
                                if _region1 == 'H-aliphatic':
                                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                                        _isotope2 = _dict2['atom_isotope_number']
                                        if _dim_id1 == _dim_id2 or _isotope2 not in (1, 13):
                                            continue
                                        if _isotope2 == 13\
                                           and not any(True for _transfer in cur_spectral_dim_transfer
                                                       if _transfer['type'] == 'onebond'
                                                       and {_dim_id1, _dim_id2} == {_transfer['spectral_dim_id_1'],
                                                                                    _transfer['spectral_dim_id_2']}):
                                            if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                                transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                            'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                            'type': 'relayed',
                                                            'indirect': 'yes' if _isotope2 == 1 else 'no'}
                                                if transfer in cur_spectral_dim_transfer:
                                                    continue
                                                cur_spectral_dim_transfer.append(transfer)

                        elif d == 4:
                            for _dim_id1, _dict1 in cur_spectral_dim.items():
                                _region1 = _dict1['_spectral_region']
                                if _region1 == 'HN':
                                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                                        if _dim_id1 == _dim_id2 or _dict2['atom_isotope_number'] != 1:
                                            continue
                                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                        'type': 'relayed',
                                                        'indirect': 'yes'}
                                            if transfer in cur_spectral_dim_transfer:
                                                continue
                                            cur_spectral_dim_transfer.append(transfer)
                            for _dim_id1, _dict1 in cur_spectral_dim.items():
                                _region1 = _dict1['_spectral_region']
                                if _region1 == 'H-aliphatic':
                                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                                        if _dim_id1 == _dim_id2 or _dict2['atom_isotope_number'] != 1:
                                            continue
                                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                        'type': 'relayed',
                                                        'indirect': 'no'}
                                            if transfer in cur_spectral_dim_transfer:
                                                continue
                                            cur_spectral_dim_transfer.append(transfer)

                    # relayed-alternate: 'Relayed transfer where peaks from an odd resp.
                    # even number of transfer steps have opposite sign'

                    # through-space: 'Any transfer that does not go through the covalent bonded skeleton'
                    if any('noe' in n for n in _file_names) or any('roe' in n for n in _file_names):
                        for _dim_id1, _dict1 in cur_spectral_dim.items():
                            _region1 = _dict1['_spectral_region']
                            if _region1 in ('H', 'HN', 'H-aliphatic', 'H-aromatic', 'H-methyl'):
                                for _dim_id2, _dict2 in cur_spectral_dim.items():
                                    if _dim_id1 == _dim_id2 or _dict1['atom_isotope_number'] != _dict2['atom_isotope_number']:
                                        continue
                                    if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                        transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                    'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                    'type': 'through-space',
                                                    'indirect': 'yes'}
                                        if transfer in cur_spectral_dim_transfer:
                                            continue
                                        cur_spectral_dim_transfer.append(transfer)
                                        if d == 2 and _region1 == 'H-aliphatic':
                                            _dict1['spectral_region'] = _dict2['spectral_region'] = 'H'  # all
                                            if 'H-aliphatic' in (_dict1['axis_code'], _dict2['axis_code']):
                                                _dict1['axis_code'] = f'H{dim_to_code[_dim_id1]}'
                                                _dict2['axis_code'] = f'H{dim_to_code[_dim_id2]}'
                                        if d == 3:
                                            _transfer = next((_transfer for _transfer in cur_spectral_dim_transfer
                                                              if _transfer['type'] == 'onebond'), None)
                                            if _transfer is not None:
                                                if _dim_id1 not in (_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']):
                                                    if _region1 == 'H-aliphatic':
                                                        _dict1['spectral_region'] = 'H'  # all
                                                        if _dict1['axis_code'] == 'H-aliphatic':
                                                            _dict1['axis_code'] = 'H'
                                                elif _region1 == 'HN':
                                                    _dict1['spectral_region'] = _dict1['axis_code'] = 'HN'
                                                if _dim_id2 not in (_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']):
                                                    if _region1 == 'H-aliphatic':
                                                        _dict2['spectral_region'] = 'H'  # all
                                                        if _dict2['axis_code'] == 'H-aliphatic':
                                                            _dict2['axis_code'] = 'H'
                                                elif _dict2['_spectral_region'] == 'HN':
                                                    _dict2['spectral_region'] = _dict2['axis_code'] = 'HN'

                    if self.exptlMethod == 'SOLID-STATE NMR':
                        if any('rfdr' in n for n in _file_names) or any('darr' in n for n in _file_names):
                            for _dim_id1, _dict1 in cur_spectral_dim.items():
                                _iso_num1 = _dict1['atom_isotope_number']
                                if _iso_num1 in (1, 13):
                                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                                        if _dim_id1 == _dim_id2 or _iso_num1 != _dict2['atom_isotope_number']:
                                            continue
                                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                        'type': 'through-space',
                                                        'indirect': 'yes'}
                                            if transfer in cur_spectral_dim_transfer:
                                                continue
                                            cur_spectral_dim_transfer.append(transfer)
                                            if d == 2:  # and _dict1['spectral_region'] == _dict2['spectral_region']:
                                                nuc = _dict1[1]['spectral_regison'][0]
                                                _dict1['spectral_region'] = _dict2['spectral_region'] = nuc
                                                if _dict1['axis_code'] == _dict2['axis_code']:
                                                    _dict1['axis_code'] = f'{nuc}{dim_to_code[_dim_id1]}'
                                                    _dict2['axis_code'] = f'{nuc}{dim_to_code[_dim_id2]}'

                        elif any('redor' in n for n in _file_names) or any('tedor' in n for n in _file_names):
                            for _dim_id1, _dict1 in cur_spectral_dim.items():
                                _iso_num1 = _dict1['atom_isotope_number']
                                if _iso_num1 in (13, 15, 19, 31):
                                    for _dim_id2, _dict2 in cur_spectral_dim.items():
                                        _iso_num2 = _dict2['atom_isotope_number']
                                        if _dim_id1 == _dim_id2 or _iso_num2 not in (13, 15, 19, 31) or _iso_num1 == _iso_num2:
                                            continue
                                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                        'type': 'through-space',
                                                        'indirect': 'yes'}
                                            if transfer in cur_spectral_dim_transfer:
                                                continue
                                            cur_spectral_dim_transfer.append(transfer)

                    for _dim_id1, _dict1 in cur_spectral_dim.items():
                        _region1 = _dict1['_spectral_region']
                        if _region1 in ('H', 'HN', 'H-aliphatic', 'H-aromatic', 'H-methyl') and d > 2:
                            for _dim_id2, _dict2 in cur_spectral_dim.items():
                                if _dim_id1 == _dim_id2 or _dict1['atom_isotope_number'] != _dict2['atom_isotope_number']:
                                    continue
                                if not any(True for _transfer in cur_spectral_dim_transfer
                                           if {_dim_id1, _dim_id2} == {_transfer['spectral_dim_id_1'],
                                                                       _transfer['spectral_dim_id_2']}):
                                    if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                        transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                    'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                    'type':
                                                    'through-space' if 'long_range' in cur_spectral_dim[1] else 'through-space?',
                                                    'indirect': 'yes'}
                                        if transfer in cur_spectral_dim_transfer:
                                            continue
                                        cur_spectral_dim_transfer.append(transfer)

                    for _dim_id1, _dict1 in cur_spectral_dim.items():
                        _region1 = _dict1['_spectral_region']
                        if _region1 == 'H' and d == 2:  # all
                            for _dim_id2, _dict2 in cur_spectral_dim.items():
                                if _dim_id1 == _dim_id2:
                                    continue
                                if 'long_range' not in cur_spectral_dim[1]:
                                    if _dict2['_spectral_region'] != _region1:
                                        continue
                                if not any(True for _transfer in cur_spectral_dim_transfer
                                           if {_dim_id1, _dim_id2} == {_transfer['spectral_dim_id_1'],
                                                                       _transfer['spectral_dim_id_2']}):
                                    if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                        transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                    'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                    'type':
                                                    'through-space' if 'long_range' in cur_spectral_dim[1] else 'through-space?',
                                                    'indirect': 'yes'}
                                        if transfer in cur_spectral_dim_transfer:
                                            continue
                                        cur_spectral_dim_transfer.append(transfer)

                    if self.exptlMethod == 'SOLID-STATE NMR' and d == 2:
                        for _dim_id1, _dict1 in cur_spectral_dim.items():
                            _region1 = _dict1['_spectral_region']
                            if _region1 == 'C':  # all
                                for _dim_id2, _dict2 in cur_spectral_dim.items():
                                    if _dim_id1 == _dim_id2:
                                        continue
                                    if 'long_range' not in cur_spectral_dim[1]:
                                        _region2 = _dict2['_spectral_region']
                                        if _region2 != _region1:
                                            if {_region1, _region2} != {'C', 'C-aliphatic'}:  # darr
                                                continue
                                    if not any(True for _transfer in cur_spectral_dim_transfer
                                               if {_dim_id1, _dim_id2} == {_transfer['spectral_dim_id_1'],
                                                                           _transfer['spectral_dim_id_2']}):
                                        if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                            transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                        'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                        'type':
                                                        'through-space' if 'long_range' in cur_spectral_dim[1]
                                                        else 'through-space?',
                                                        'indirect': 'yes'}
                                            if transfer in cur_spectral_dim_transfer:
                                                continue
                                            cur_spectral_dim_transfer.append(transfer)

                        for _dim_id1, _dict1 in cur_spectral_dim.items():
                            _iso_num1 = _dict1['atom_isotope_number']
                            if _iso_num1 in (13, 15, 19, 31):
                                for _dim_id2, _dict2 in cur_spectral_dim.items():
                                    _iso_num2 = _dict2['atom_isotope_number']
                                    if _dim_id1 == _dim_id2 or _iso_num2 not in (13, 15, 19, 31) or _iso_num1 == _iso_num2:
                                        continue
                                    if 'yes' in (_dict1['acquisition'], _dict2['acquisition']):
                                        transfer = {'spectral_dim_id_1': min([_dim_id1, _dim_id2]),
                                                    'spectral_dim_id_2': max([_dim_id1, _dim_id2]),
                                                    'type':
                                                    'through-space' if 'long_range' in cur_spectral_dim[1] else 'through-space?',
                                                    'indirect': 'yes'}
                                        if transfer in cur_spectral_dim_transfer:
                                            continue
                                        cur_spectral_dim_transfer.append(transfer)

                    for __v in cur_spectral_dim.values():
                        if 'freq_hint' in __v:
                            del __v['freq_hint']

                    if self.__debug:
                        print(f'num_of_dim: {d}, list_id: {_id}')
                        print('spectral_dim')
                        for __d, __v in cur_spectral_dim.items():
                            print(f'{__d} {__v}')
                        print('spectral_dim_transfer')
                        for transfer in cur_spectral_dim_transfer:
                            print(transfer)

                    primary_dim_transfer = ''
                    for transfer in cur_spectral_dim_transfer:
                        if transfer['type'] == 'through-space':
                            primary_dim_transfer = transfer['type']
                            break
                        if transfer['type'] == 'through-space?'\
                                and primary_dim_transfer != 'through-space':
                            primary_dim_transfer = transfer['type']
                        elif transfer['type'].startswith('relayed')\
                                and primary_dim_transfer not in ('through-space', 'through-space?'):
                            primary_dim_transfer = transfer['type']
                        elif transfer['type'] == 'jmultibond'\
                                and primary_dim_transfer not in ('through-space', 'through-space?',
                                                                 'relayed', 'relayed-alternate'):
                            primary_dim_transfer = transfer['type']
                        elif transfer['type'] == 'jcoupling'\
                                and primary_dim_transfer not in ('through-space', 'through-space?',
                                                                 'relayed', 'relayed-alternate', 'jmultibond'):
                            primary_dim_transfer = transfer['type']
                        elif transfer['type'] == 'onebond'\
                                and primary_dim_transfer == '':
                            primary_dim_transfer = transfer['type']

                    exp_class = '.'

                    if primary_dim_transfer != 'onebond':
                        onebonds, onebond_codes = [], []
                        for transfer in cur_spectral_dim_transfer:
                            if transfer['type'] == 'onebond':
                                dim_id_1 = transfer['spectral_dim_id_1']
                                dim_id_2 = transfer['spectral_dim_id_2']
                                atom_type_1 = next(v['atom_type'] for k, v in cur_spectral_dim.items() if k == dim_id_1)
                                atom_type_2 = next(v['atom_type'] for k, v in cur_spectral_dim.items() if k == dim_id_2)
                                for _transfer in cur_spectral_dim_transfer:
                                    if _transfer['type'] != 'onebond':
                                        if dim_id_1 in (_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']):
                                            onebonds.append((dim_id_1, dim_id_2))
                                            onebond_codes.append(f'{atom_type_1}[{atom_type_2}]')
                                        if dim_id_2 in (_transfer['spectral_dim_id_1'], _transfer['spectral_dim_id_2']):
                                            onebonds.append((dim_id_2, dim_id_1))
                                            onebond_codes.append(f'{atom_type_2}[{atom_type_1}]')

                        for transfer in cur_spectral_dim_transfer:
                            if transfer['type'] == primary_dim_transfer:
                                dim_id_1 = transfer['spectral_dim_id_1']
                                dim_id_2 = transfer['spectral_dim_id_2']
                                onebond_1 = next((onebond for onebond in onebonds if dim_id_1 in onebond), None)
                                v_1 = next(v for k, v in cur_spectral_dim.items() if k == dim_id_1)
                                onebond_2 = next((onebond for onebond in onebonds if dim_id_2 in onebond), None)
                                v_2 = next(v for k, v in cur_spectral_dim.items() if k == dim_id_2)
                                if v_1['acquisition'] != v_2['acquisition'] and v_1['acquisition'] == 'no':
                                    onebond_1, onebond_2 = onebond_2, onebond_1
                                    v_1, v_2 = v_2, v_1
                                # pylint: disable=line-too-long
                                exp_class = f'{onebond_codes[onebonds.index(onebond_1)] if onebond_1 is not None else v_1["atom_type"]}_'\
                                    f'{onebond_codes[onebonds.index(onebond_2)] if onebond_2 is not None else v_2["atom_type"]}.'\
                                    f'{primary_dim_transfer}'  # noqa: E501
                                break

                    else:
                        for transfer in cur_spectral_dim_transfer:
                            if transfer['type'] == primary_dim_transfer:
                                dim_id_1 = transfer['spectral_dim_id_1']
                                dim_id_2 = transfer['spectral_dim_id_2']
                                v_1 = next(v for k, v in cur_spectral_dim.items() if k == dim_id_1)
                                v_2 = next(v for k, v in cur_spectral_dim.items() if k == dim_id_2)
                                if v_1['acquisition'] != v_2['acquisition'] and v_1['acquisition'] == 'no':
                                    v_1, v_2 = v_2, v_1
                                exp_class = f'{v_1["atom_type"]}_{v_2["atom_type"]}.{primary_dim_transfer}'
                                break

                    if self.__debug:
                        print(f'experiment class: {exp_class}')

                    if self.software_name == 'PIPP'\
                       and any(transfer['type'] == 'onebond' for transfer in cur_spectral_dim_transfer):
                        if d == 2:
                            transfer = next(transfer for transfer in cur_spectral_dim_transfer if transfer['type'] == 'onebond')
                            pro_axis = hvy_axis = -1
                            for _dim_id, _dict in cur_spectral_dim.items():
                                if _dim_id in (transfer['spectral_dim_id_1'], transfer['spectral_dim_id_2']):
                                    if _dict['atom_isotope_number'] == 1:
                                        pro_axis = _dim_id
                                    else:
                                        hvy_axis = _dim_id
                            if pro_axis != -1 and hvy_axis != -1:
                                self.reasonsForReParsing['onebond_resolved'] = {0: hvy_axis - 1, 1: pro_axis - 1}
                        elif d == 3:
                            transfer = next(transfer for transfer in cur_spectral_dim_transfer if transfer['type'] == 'onebond')
                            pro_axis = hvy_axis = -1
                            for _dim_id, _dict in cur_spectral_dim.items():
                                if _dim_id in (transfer['spectral_dim_id_1'], transfer['spectral_dim_id_2']):
                                    if _dict['atom_isotope_number'] == 1:
                                        pro_axis = _dim_id
                                    else:
                                        hvy_axis = _dim_id
                            if pro_axis != -1 and hvy_axis != -1:
                                self.reasonsForReParsing['onebond_resolved'] = {0: {0: hvy_axis - 1,
                                                                                    1: pro_axis - 1,
                                                                                    2: 5 - hvy_axis - pro_axis},
                                                                                1: {1: hvy_axis - 1,
                                                                                    2: pro_axis - 1,
                                                                                    0: 5 - hvy_axis - pro_axis}}
                        elif d == 4 and len([transfer for transfer in cur_spectral_dim_transfer
                                             if transfer['type'] == 'onebond']) == 2:
                            r = self.reasonsForReParsing['onebond_resolved'] = {}
                            for offset, transfer in enumerate([transfer for transfer in cur_spectral_dim_transfer
                                                               if transfer['type'] == 'onebond']):
                                pro_axis = hvy_axis = -1
                                for _dim_id, _dict in cur_spectral_dim.items():
                                    if _dim_id in (transfer['spectral_dim_id_1'], transfer['spectral_dim_id_2']):
                                        if _dict['atom_isotope_number'] == 1:
                                            pro_axis = _dim_id
                                        else:
                                            hvy_axis = _dim_id
                                if pro_axis != -1 and hvy_axis != -1:
                                    r[offset * 2] = hvy_axis - 1
                                    r[offset * 2 + 1] = pro_axis - 1

                    if self.__createSfDict:
                        self.cur_subtype = f'peak{d}d'
                        self.cur_list_id = _id

                        sf = self.getSf()

                        if sf['peak_row_format']:
                            sf['saveframe'].add_loop(sf['loop'])
                            del sf['alt_loops']
                        else:
                            for alt_loop in sf['alt_loops']:
                                sf['saveframe'].add_loop(alt_loop)
                            del sf['loop']

                        list_id = sf['list_id']
                        sf['aux_loops'] = getAuxLoops('spectral_peak')

                        aux_lp = next((aux_lp for aux_lp in sf['aux_loops'] if aux_lp.category == '_Spectral_dim'), None)

                        if aux_lp is None:
                            continue

                        for _dim_id, _dict in cur_spectral_dim.items():
                            aux_lp.add_data(getSpectralDimRow(_dim_id, list_id, self.__entryId, _dict))

                        sf['saveframe'].add_loop(aux_lp)

                        aux_lp = next((aux_lp for aux_lp in sf['aux_loops'] if aux_lp.category == '_Spectral_dim_transfer'), None)

                        if aux_lp is None:
                            continue

                        for _dict in cur_spectral_dim_transfer:
                            aux_lp.add_data(getSpectralDimTransferRow(list_id, self.__entryId, _dict))

                        sf['saveframe'].add_loop(aux_lp)

                        if exp_class not in EMPTY_VALUE:
                            tags = [t[0] for t in sf['saveframe'].tags]
                            if 'Experiment_class' in tags:
                                sf['saveframe'].tags[tags.index('Experiment_class')][1] = exp_class

                        if spectrum_names is not None:
                            try:
                                alt_file_name = spectrum_names[d][_id]
                                if alt_file_name not in EMPTY_VALUE:
                                    exp_type = get_first_sf_tag(sf['saveframe'], 'Experiment_type')
                                    if exp_type in EMPTY_VALUE:
                                        set_sf_tag(sf['saveframe'], 'Experiment_type', alt_file_name)
                            except (KeyError, AttributeError):
                                pass

                        has_assign = False

                        if sf['peak_row_format']:
                            lp = sf['loop']
                            for row in lp.get_tag('Auth_asym_ID_1'):
                                if row not in EMPTY_VALUE:
                                    has_assign = True
                                    break

                        else:
                            lp = next((alt_lp for alt_lp in sf['alt_loops']
                                       if alt_lp.category == '_Assigned_peak_chem_shift'), None)
                            if lp is not None and len(lp) > 0:
                                for row in lp.get_tag('Auth_entity_ID'):
                                    if row not in EMPTY_VALUE:
                                        has_assign = True
                                        break

                        if not has_assign:
                            continue

                        self.__remediatePeakAssignmentForAtomType(d, self.atom_type_history[d][_id], sf['peak_row_format'], lp)

                        self.__remediateIncompletePeakAssignment(d, sf['peak_row_format'], lp)

                        if any(transfer['type'] == 'onebond' for transfer in cur_spectral_dim_transfer):
                            onebond_dim_transfers = [[transfer['spectral_dim_id_1'], transfer['spectral_dim_id_2']]
                                                     for transfer in cur_spectral_dim_transfer
                                                     if transfer['type'] == 'onebond']

                            self.__remediatePeakAssignmentForOneBondTransfer(d, onebond_dim_transfers, sf['peak_row_format'], lp)

                        if any(transfer['type'] == 'jcoupling' for transfer in cur_spectral_dim_transfer):
                            jcoupling_dim_transfers = [[transfer['spectral_dim_id_1'], transfer['spectral_dim_id_2']]
                                                       for transfer in cur_spectral_dim_transfer
                                                       if transfer['type'] == 'jcoupling']

                            for (_dim_id_1, _dim_id_2) in jcoupling_dim_transfers:
                                _pair = set([_dim_id_1 - 1, _dim_id_2 - 1])
                                if d == 2:
                                    for jcoupling_idx, (__dim_id_1, __dim_id_2) in enumerate(DIM_TRANSFER_PAT_2D):
                                        __pair = set([__dim_id_1, __dim_id_2])
                                        if _pair == __pair:
                                            self.jcoupling_idx_history[d][_id] = jcoupling_idx
                                            break
                                elif d == 3:
                                    for jcoupling_idx, (__dim_id_1, __dim_id_2) in enumerate(DIM_TRANSFER_PAT_3D):
                                        __pair = set([__dim_id_1, __dim_id_2])
                                        if _pair == __pair:
                                            self.jcoupling_idx_history[d][_id] = jcoupling_idx
                                            break
                                elif d == 4:
                                    for jcoupling_idx, ((__dim_id_1, __dim_id_2), (__dim_id_3, __dim_id_4))\
                                            in enumerate(DIM_TRANSFER_PAT_3D):
                                        __pair_1 = set([__dim_id_1, __dim_id_2])
                                        __pair_2 = set([__dim_id_3, __dim_id_4])
                                        if _pair in (__pair_1, __pair_2):
                                            self.jcoupling_idx_history[d][_id] = jcoupling_idx
                                            break

                            self.__remediatePeakAssignmentForJcouplingTransfer(d, jcoupling_dim_transfers, sf['peak_row_format'],
                                                                               lp)

                        if any(transfer['type'] == 'relayed' for transfer in cur_spectral_dim_transfer)\
                           and self.__csLoops is not None and len(self.__csLoops) > 0:
                            relayed_dim_transfers = [[transfer['spectral_dim_id_1'], transfer['spectral_dim_id_2']]
                                                     for transfer in cur_spectral_dim_transfer
                                                     if transfer['type'] == 'relayed']

                            for (_dim_id_1, _dim_id_2) in relayed_dim_transfers:
                                _pair = set([_dim_id_1 - 1, _dim_id_2 - 1])
                                if d == 2:
                                    for relayed_idx, (__dim_id_1, __dim_id_2) in enumerate(DIM_TRANSFER_PAT_2D):
                                        __pair = set([__dim_id_1, __dim_id_2])
                                        if _pair == __pair:
                                            self.relayed_idx_history[d][_id] = relayed_idx
                                            break
                                elif d == 3:
                                    for relayed_idx, (__dim_id_1, __dim_id_2) in enumerate(DIM_TRANSFER_PAT_3D):
                                        __pair = set([__dim_id_1, __dim_id_2])
                                        if _pair == __pair:
                                            self.relayed_idx_history[d][_id] = relayed_idx
                                            break
                                elif d == 4:
                                    for relayed_idx, ((__dim_id_1, __dim_id_2), (__dim_id_3, __dim_id_4))\
                                            in enumerate(DIM_TRANSFER_PAT_3D):
                                        __pair_1 = set([__dim_id_1, __dim_id_2])
                                        __pair_2 = set([__dim_id_3, __dim_id_4])
                                        if _pair in (__pair_1, __pair_2):
                                            self.relayed_idx_history[d][_id] = relayed_idx
                                            break

                            self.__remediatePeakAssignmentForRelayedTransfer(d, relayed_dim_transfers, sf['peak_row_format'], lp)

    def __canRemediatePeakAssignmentForOneBondTransfer(self, atom1: dict, atom2: dict, position: float, position2: float) -> bool:
        """ Return whether peak assignment remediation based on onebond transfer is possible.
        """

        chain_id, seq_id, comp_id, atom_id, chain_id2, seq_id2, comp_id2, atom_id2 =\
            atom1['chain_id'], atom1['seq_id'], atom1['comp_id'], atom1['atom_id'], \
            atom2['chain_id'], atom2['seq_id'], atom2['comp_id'], atom2['atom_id']

        if chain_id == chain_id2 and seq_id == seq_id2 and atom_id != atom_id2\
           and self.ccU.updateChemCompDict(comp_id):
            _atom_ids = self.nefT.get_valid_star_atom(comp_id, atom_id, leave_unmatched=False)[0]
            _atom_ids2 = self.nefT.get_valid_star_atom(comp_id, atom_id2, leave_unmatched=False)[0]

            if any(True for b in self.ccU.lastBondDictList
                   if ((b['atom_id_1'] in _atom_ids and b['atom_id_2'] in _atom_ids2)
                       or (b['atom_id_1'] in _atom_ids2 and b['atom_id_2'] in _atom_ids))):
                return True

            _atom_id, _atom_id2 = _atom_ids[0], _atom_ids2[0]

            _atom_id2_ = self.ccU.getBondedAtoms(comp_id, _atom_id,
                                                 exclProton=_atom_id[0] in PROTON_BEGIN_CODE,
                                                 onlyProton=_atom_id[0] not in PROTON_BEGIN_CODE)
            _atom_id_ = self.ccU.getBondedAtoms(comp_id, _atom_id2,
                                                exclProton=_atom_id2[0] in PROTON_BEGIN_CODE,
                                                onlyProton=_atom_id2[0] not in PROTON_BEGIN_CODE)

            len_atom_id_ = len(_atom_id_)
            len_atom_id2_ = len(_atom_id2_)

            shift, weight = self.__getCsValue(chain_id, seq_id, comp_id, _atom_id)
            shift2, weight2 = self.__getCsValue(chain_id, seq_id, comp_id, _atom_id2)

            shift_ = shift2_ = None
            if len_atom_id_ > 0 and _atom_id_[0][0] == _atom_id[0]:
                shift_, _ = self.__getCsValue(chain_id, seq_id, comp_id, _atom_id_[0])
            if len_atom_id2_ > 0 and _atom_id2_[0][0] == _atom_id2[0]:
                shift2_, _ = self.__getCsValue(chain_id, seq_id, comp_id, _atom_id2_[0])

            if shift_ is not None and not self.ccU.hasBond(comp_id, _atom_id_[0], _atom_id2):
                shift_ = None

            if shift2_ is not None and not self.ccU.hasBond(comp_id, _atom_id, _atom_id2_[0]):
                shift2_ = None

            if None not in (shift, shift2):
                diff = ((position - shift) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                diff *= 2.0

                diff_ = diff2_ = None
                if shift_ is not None:
                    diff_ = ((position - shift_) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                if shift2_ is not None:
                    diff2_ = ((position - shift) * weight) ** 2 + ((position2 - shift2_) * weight2) ** 2

                if diff_ is not None and diff2_ is not None:
                    if diff_ < diff and diff2_ < diff:
                        if diff_ < diff2_ and len_atom_id_ > 0:
                            return True
                        if diff_ > diff2_ and len_atom_id2_ > 0:
                            return True
                    elif diff_ < diff and len_atom_id_ > 0:
                        return True
                    elif diff2_ < diff and len_atom_id2_ > 0:
                        return True
                    else:
                        if diff_ < diff2_ and diff_ < 1.0 and len_atom_id_ > 0:
                            return True
                        if diff_ > diff2_ and diff2_ < 1.0 and len_atom_id2_ > 0:
                            return True

                elif diff_ is not None and (diff_ < diff or diff_ < 1.0) and len_atom_id_ > 0:
                    return True

                elif diff2_ is not None and (diff2_ < diff or diff2_ < 1.0) and len_atom_id2_ > 0:
                    return True

            return False

        if chain_id == chain_id2 and seq_id != seq_id2:
            shift, weight = self.__getCsValue(chain_id, seq_id, comp_id, atom_id)
            shift2, weight2 = self.__getCsValue(chain_id, seq_id2, comp_id2, atom_id2)

            shift_, _ = self.__getCsValue(chain_id, seq_id2, comp_id2, atom_id)
            shift2_, _ = self.__getCsValue(chain_id, seq_id, comp_id, atom_id2)

            if shift_ is not None and not self.ccU.hasBond(comp_id2, atom_id, atom_id2):
                shift_ = None

            if shift2_ is not None and not self.ccU.hasBond(comp_id, atom_id, atom_id2):
                shift2_ = None

            if None in (shift, shift2):

                if shift is None and shift2 is None:
                    return False

                if shift is None and shift_ is not None:
                    return True

                if shift2 is None and shift2_ is not None:
                    return True

                if shift is None and shift_ is None and None not in (shift2, shift2_):
                    return True

                if shift2 is None and shift2_ is None and None not in (shift, shift_):
                    return True

            else:

                diff = ((position - shift) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                diff *= 2.0

                diff_ = diff2_ = None
                if shift_ is not None:
                    diff_ = ((position - shift_) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                if shift2_ is not None:
                    diff2_ = ((position - shift) * weight) ** 2 + ((position2 - shift2_) * weight2) ** 2

                if diff_ is not None and diff2_ is not None:
                    if diff_ < diff and diff2_ < diff:
                        if diff_ < diff2_:
                            return True
                        if diff_ > diff2_:
                            return True
                    elif diff_ < diff:
                        return True
                    elif diff2_ < diff:
                        return True
                    else:
                        if diff_ < diff2_ and diff_ < 1.0:
                            return True
                        if diff_ > diff2_ and diff2_ < 1.0:
                            return True

                    return False

                if diff_ is not None and (diff_ < diff or diff_ < 1.0):
                    return True

                if diff2_ is not None and (diff2_ < diff or diff2_ < 1.0):
                    return True

                _atom_ids = self.nefT.get_valid_star_atom(comp_id, atom_id, leave_unmatched=False)[0]
                _atom_ids2 = self.nefT.get_valid_star_atom(comp_id2, atom_id2, leave_unmatched=False)[0]

                _atom_id, _atom_id2 = _atom_ids[0], _atom_ids2[0]

                _atom_id2_ = self.ccU.getBondedAtoms(comp_id, _atom_id,
                                                     exclProton=_atom_id[0] in PROTON_BEGIN_CODE,
                                                     onlyProton=_atom_id[0] not in PROTON_BEGIN_CODE)
                _atom_id_ = self.ccU.getBondedAtoms(comp_id2, _atom_id2,
                                                    exclProton=_atom_id2[0] in PROTON_BEGIN_CODE,
                                                    onlyProton=_atom_id2[0] not in PROTON_BEGIN_CODE)

                len_atom_id_ = len(_atom_id_)
                len_atom_id2_ = len(_atom_id2_)

                shift_ = shift2_ = None
                if len_atom_id_ > 0 and _atom_id_[0][0] == _atom_id[0]:
                    shift_, _ = self.__getCsValue(chain_id, seq_id2, comp_id2, _atom_id_[0])
                if len_atom_id2_ > 0 and _atom_id2_[0][0] == _atom_id2[0]:
                    shift2_, _ = self.__getCsValue(chain_id, seq_id, comp_id, _atom_id2_[0])

                diff_ = diff2_ = None
                if shift_ is not None:
                    diff_ = ((position - shift_) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                if shift2_ is not None:
                    diff2_ = ((position - shift) * weight) ** 2 + ((position2 - shift2_) * weight2) ** 2

                if diff_ is not None and diff2_ is not None:
                    if diff_ < diff and diff2_ < diff:
                        if diff_ < diff2_ and len_atom_id_ > 0:
                            return True
                        if diff_ > diff2_ and len_atom_id2_ > 0:
                            return True
                    elif diff_ < diff and len_atom_id_ > 0:
                        return True
                    elif diff2_ < diff and len_atom_id2_ > 0:
                        return True
                    else:
                        if diff_ < diff2_ and diff_ < 1.0 and len_atom_id_ > 0:
                            return True
                        if diff_ > diff2_ and diff2_ < 1.0 and len_atom_id2_ > 0:
                            return True

                elif diff_ is not None and (diff_ < diff or diff_ < 1.0) and len_atom_id_ > 0:
                    return True

                elif diff2_ is not None and (diff2_ < diff or diff2_ < 1.0) and len_atom_id2_ > 0:
                    return True

            return False

        if chain_id != chain_id2:
            shift, weight = self.__getCsValue(chain_id, seq_id, comp_id, atom_id)
            shift2, weight2 = self.__getCsValue(chain_id2, seq_id2, comp_id2, atom_id2)

            shift_, _ = self.__getCsValue(chain_id2, seq_id2, comp_id2, atom_id)
            shift2_, _ = self.__getCsValue(chain_id, seq_id, comp_id, atom_id2)

            if shift_ is not None and not self.ccU.hasBond(comp_id2, atom_id, atom_id2):
                shift_ = None

            if shift2_ is not None and not self.ccU.hasBond(comp_id, atom_id, atom_id2):
                shift2_ = None

            if None in (shift, shift2):

                if shift is None and shift2 is None:
                    return False

                if shift is None and shift_ is not None:
                    return True

                if shift2 is None and shift2_ is not None:
                    return True

                if shift is None and shift_ is None and None not in (shift2, shift2_):
                    return True

                if shift2 is None and shift2_ is None and None not in (shift, shift_):
                    return True

                return False

            diff = ((position - shift) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
            diff *= 2.0

            diff_ = diff2_ = None
            if shift_ is not None:
                diff_ = ((position - shift_) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
            if shift2_ is not None:
                diff2_ = ((position - shift) * weight) ** 2 + ((position2 - shift2_) * weight2) ** 2

            if diff_ is not None and diff2_ is not None:
                if diff_ < diff and diff2_ < diff:
                    if diff_ < diff2_:
                        return True
                    if diff_ > diff2_:
                        return True
                elif diff_ < diff:
                    return True
                elif diff2_ < diff:
                    return True
                else:
                    if diff_ < diff2_ and diff_ < 1.0:
                        return True
                    if diff_ > diff2_ and diff2_ < 1.0:
                        return True

            elif diff_ is not None and (diff_ < diff or diff_ < 1.0):
                return True

            elif diff2_ is not None and (diff2_ < diff or diff2_ < 1.0):
                return True

            return False

        return False

    def __canRemediatePeakAssignmentForJcouplingTransfer(self, atom1: dict, atom2: dict, position: float, position2: float) -> bool:
        """ Return whether peak assignment remediation based on J-coupling transfer is possible.
        """

        return self.__canRemediatePeakAssignmentForTransfer(atom1, atom2, position, position2, False)

    def __canRemediatePeakAssignmentForRelayedTransfer(self, atom1: dict, atom2: dict, position: float, position2: float) -> bool:
        """ Return whether peak assignment remediation based on relayed transfer is possible.
        """

        return self.__canRemediatePeakAssignmentForTransfer(atom1, atom2, position, position2, True)

    def __canRemediatePeakAssignmentForTransfer(self, atom1: dict, atom2: dict, position: float, position2: float, relayed: bool
                                                ) -> bool:
        """ Return whether peak assignment remediation based on relayed transfer (relayed), or J-coupling transfer otherwise,
            is possible.
        """

        chain_id, seq_id, comp_id, atom_id, chain_id2, seq_id2, comp_id2, atom_id2 =\
            atom1['chain_id'], atom1['seq_id'], atom1['comp_id'], atom1['atom_id'], \
            atom2['chain_id'], atom2['seq_id'], atom2['comp_id'], atom2['atom_id']

        if chain_id == chain_id2 and (abs(seq_id - seq_id2) < 2 if relayed else seq_id == seq_id2):
            return True

        if chain_id == chain_id2 and (relayed or seq_id != seq_id2):
            shift, weight = self.__getCsValue(chain_id, seq_id, comp_id, atom_id)
            shift2, weight2 = self.__getCsValue(chain_id, seq_id2, comp_id2, atom_id2)

            if None in (shift, shift2):
                return False

            shift_, _ = self.__getCsValue(chain_id, seq_id2, comp_id2, atom_id)
            shift2_, _ = self.__getCsValue(chain_id, seq_id, comp_id, atom_id2)

            diff = ((position - shift) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
            diff *= 2.0

            diff_ = diff2_ = None
            if shift_ is not None:
                diff_ = ((position - shift_) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
            if shift2_ is not None:
                diff2_ = ((position - shift) * weight) ** 2 + ((position2 - shift2_) * weight2) ** 2

            if diff_ is not None and diff2_ is not None:
                if diff_ < diff and diff2_ < diff:
                    if diff_ < diff2_:
                        return True
                    if diff_ > diff2_:
                        return True
                elif diff_ < diff:
                    return True
                elif diff2_ < diff:
                    return True
                else:
                    if diff_ < diff2_ and diff_ < 1.0:
                        return True
                    if diff_ > diff2_ and diff2_ < 1.0:
                        return True

            elif diff_ is not None and (diff_ < diff or diff_ < 1.0):
                return True

            elif diff2_ is not None and (diff2_ < diff or diff2_ < 1.0):
                return True

        return False

    def __remediatePeakAssignmentForAtomType(self, num_of_dim: int, atom_type: List[str],
                                             use_peak_row_format: bool, loop: pynmrstar.Loop) -> None:
        """ Remediate peak assignment for decided atom type (PIPP specific)
        """

        is_reparsable = self.reasons is None and self.software_name != 'PIPP'

        if not is_reparsable:
            return

        if use_peak_row_format:

            tags = [f'Atom_ID_{dim_id}' for dim_id in range(1, num_of_dim + 1)]

            dat = loop.get_tag(tags)

            for row in dat:

                for col in range(num_of_dim):

                    atom_id = row[col]

                    if atom_id in EMPTY_VALUE or atom_type[col] not in ASSIGNMENT_HALF_SPIN_NUCLEUS:
                        continue

                    if atom_id[0] in PROTON_BEGIN_CODE:

                        if atom_type[col] == 'H':
                            continue

                        self.reasonsForReParsing['atom_type_history'] = self.atom_type_history

                        return

                    if atom_id[0] != atom_type[col]:

                        self.reasonsForReParsing['atom_type_history'] = self.atom_type_history

                        return

        else:

            tags = ['Spectral_dim_ID', 'Atom_ID']

            dat = loop.get_tag(tags)

            for row in dat:

                col, atom_id = row[0] - 1, row[1]

                if atom_id in EMPTY_VALUE or atom_type[col] not in ASSIGNMENT_HALF_SPIN_NUCLEUS:
                    continue

                if atom_id[0] in PROTON_BEGIN_CODE:

                    if atom_type[col] == 'H':
                        continue

                    self.reasonsForReParsing['atom_type_history'] = self.atom_type_history

                    return

                if atom_id[0] != atom_type[col]:

                    self.reasonsForReParsing['atom_type_history'] = self.atom_type_history

                    return

    def __remediateIncompletePeakAssignment(self, num_of_dim: int, use_peak_row_format: bool, loop: pynmrstar.Loop  # noqa: E501, pylint: disable=no-self-use,line-too-long
                                            ) -> None:
        """ Remediate incomplete peak assignment.
        """

        del_idx_list = []

        if use_peak_row_format:
            tags = [f'Atom_ID_{_dim_id}' for _dim_id in range(1, num_of_dim + 1)]

            dat = loop.get_tag(tags)

            for idx, row in enumerate(dat):
                if all(row[_dim_id] in EMPTY_VALUE for _dim_id in range(num_of_dim))\
                   or all(row[_dim_id] not in EMPTY_VALUE for _dim_id in range(num_of_dim)):
                    continue
                del_idx_list.append(idx)

        else:
            tags = ['Peak_ID', 'Set_ID']

            dat = loop.get_tag(tags)

            peak_set = {}

            for row in dat:
                key = (row[0], row[1])
                if key not in peak_set:
                    peak_set[key] = 1
                else:
                    peak_set[key] += 1

            incomplete_peak_idx = []

            for k, v in peak_set.items():
                if v == num_of_dim:
                    continue
                incomplete_peak_idx.append(k)

            if len(incomplete_peak_idx) > 0:

                for idx, row in enumerate(dat):
                    key = (row[0], row[1])
                    if key in incomplete_peak_idx:
                        del_idx_list.append(idx)

        if len(del_idx_list) > 0:
            for _idx in reversed(del_idx_list):
                del loop.data[_idx]

    def __remediatePeakAssignmentForOneBondTransfer(self, num_of_dim: int, onebond_transfers: List[List[int]],
                                                    use_peak_row_format: bool, loop: pynmrstar.Loop) -> None:
        """ Remediate peak assignment based on onebond transfer.
        """

        is_reparsable = self.reasons is None and self.software_name != 'PIPP'

        details_col = loop.tags.index('Details')

        del_idx_list = []

        for dim_id_1, dim_id_2 in onebond_transfers:

            if use_peak_row_format:

                tags = [f'Entity_assembly_ID_{dim_id_1}', f'Comp_index_ID_{dim_id_1}',
                        f'Comp_ID_{dim_id_1}', f'Atom_ID_{dim_id_1}', f'Position_{dim_id_1}',
                        f'Entity_assembly_ID_{dim_id_2}', f'Comp_index_ID_{dim_id_2}',
                        f'Comp_ID_{dim_id_2}', f'Atom_ID_{dim_id_2}', f'Position_{dim_id_2}']

                dat = loop.get_tag(tags)

                for idx, row in enumerate(dat):

                    if any(True for col in range(10) if row[col] in EMPTY_VALUE):
                        continue

                    chain_id, seq_id, comp_id, atom_id, chain_id2, seq_id2, comp_id2, atom_id2 =\
                        row[0], row[1], row[2], row[3], row[5], row[6], row[7], row[8]

                    if isinstance(chain_id, int):
                        chain_id = str(chain_id)

                    if isinstance(chain_id2, int):
                        chain_id2 = str(chain_id2)

                    if isinstance(seq_id, str):
                        seq_id = int(seq_id)

                    if isinstance(seq_id2, str):
                        seq_id2 = int(seq_id2)

                    if chain_id == chain_id2 and seq_id == seq_id2 and atom_id != atom_id2\
                       and self.ccU.updateChemCompDict(comp_id):
                        _atom_ids = self.nefT.get_valid_star_atom(comp_id, atom_id, leave_unmatched=False)[0]
                        _atom_ids2 = self.nefT.get_valid_star_atom(comp_id, atom_id2, leave_unmatched=False)[0]

                        if any(True for b in self.ccU.lastBondDictList
                               if ((b['atom_id_1'] in _atom_ids and b['atom_id_2'] in _atom_ids2)
                                   or (b['atom_id_1'] in _atom_ids2 and b['atom_id_2'] in _atom_ids))):
                            continue

                        if self.software_name == 'PIPP':
                            del_idx_list.append(idx)
                            continue

                        _atom_id, _atom_id2 = _atom_ids[0], _atom_ids2[0]

                        _atom_id2_ = self.ccU.getBondedAtoms(comp_id, _atom_id,
                                                             exclProton=_atom_id[0] in PROTON_BEGIN_CODE,
                                                             onlyProton=_atom_id[0] not in PROTON_BEGIN_CODE)
                        _atom_id_ = self.ccU.getBondedAtoms(comp_id, _atom_id2,
                                                            exclProton=_atom_id2[0] in PROTON_BEGIN_CODE,
                                                            onlyProton=_atom_id2[0] not in PROTON_BEGIN_CODE)

                        len_atom_id_ = len(_atom_id_)
                        len_atom_id2_ = len(_atom_id2_)

                        position, position2 = row[4], row[9]

                        if isinstance(position, str):
                            position = float(position)

                        if isinstance(position2, str):
                            position2 = float(position2)

                        shift, weight = self.__getCsValue(chain_id, seq_id, comp_id, _atom_id)
                        shift2, weight2 = self.__getCsValue(chain_id, seq_id, comp_id, _atom_id2)

                        _common_atom_id_ = _common_atom_id2_ = shift_ = shift2_ = None
                        if len_atom_id_ > 0:
                            _common_atom_id_ = _atom_id_[0]
                            if _common_atom_id_[0] == _atom_id[0]:
                                shift_, _ = self.__getCsValue(chain_id, seq_id, comp_id, _common_atom_id_)
                        if len_atom_id2_ > 0:
                            _common_atom_id2_ = _atom_id2_[0]
                            if _common_atom_id2_[0] == _atom_id2[0]:
                                shift2_, _ = self.__getCsValue(chain_id, seq_id, comp_id, _common_atom_id2_)

                        if shift_ is not None and not self.ccU.hasBond(comp_id, _common_atom_id_, _atom_id2):
                            shift_ = None

                        if shift2_ is not None and not self.ccU.hasBond(comp_id, _atom_id, _common_atom_id2_):
                            shift2_ = None

                        # pylint: disable=cell-var-from-loop
                        def swap_atom_id_1(atom_id_list):
                            common_atom_id = atom_id_list[0]
                            if len(atom_id_list) > 1 or len(_atom_ids2) > 1:
                                atom_sel = [{'atom_id': _a_} for _a_ in atom_id_list]
                                if len(_atom_ids2) > 1:
                                    for _atom_id2 in _atom_ids2[1:]:
                                        __atom_id__ = self.ccU.getBondedAtoms(comp_id, _atom_id2,
                                                                              exclProton=_atom_id2[0] in PROTON_BEGIN_CODE,
                                                                              onlyProton=_atom_id2[0] not in PROTON_BEGIN_CODE)
                                        for __a__ in __atom_id__:
                                            atom_sel.append({'atom_id': __a__})
                                try:
                                    common_atom_id = self.__extractCommonAtom(atom_sel)['auth_atom_id']
                                except KeyError:
                                    pass
                            if loop.data[idx][details_col] in EMPTY_VALUE:
                                loop.data[idx][details_col] = f'{atom_id} -> {common_atom_id}'
                            loop.data[idx][loop.tags.index(f'Atom_ID_{dim_id_1}')] =\
                                loop.data[idx][loop.tags.index(f'Auth_atom_ID_{dim_id_1}')] = common_atom_id

                        # pylint: disable=cell-var-from-loop
                        def swap_atom_id_2(atom_id2_list):
                            common_atom_id2 = atom_id2_list[0]
                            if len(atom_id2_list) > 1 or len(_atom_ids) > 1:
                                atom_sel = [{'atom_id': _a_} for _a_ in atom_id2_list]
                                if len(_atom_ids) > 1:
                                    for _atom_id in _atom_ids[1:]:
                                        __atom_id2__ = self.ccU.getBondedAtoms(comp_id, _atom_id,
                                                                               exclProton=_atom_id[0] in PROTON_BEGIN_CODE,
                                                                               onlyProton=_atom_id[0] not in PROTON_BEGIN_CODE)
                                        for __a__ in __atom_id2__:
                                            atom_sel.append({'atom_id': __a__})
                                try:
                                    common_atom_id2 = self.__extractCommonAtom(atom_sel)['auth_atom_id']
                                except KeyError:
                                    pass
                            if loop.data[idx][details_col] in EMPTY_VALUE:
                                loop.data[idx][details_col] = f'{atom_id2} -> {common_atom_id2}'
                            loop.data[idx][loop.tags.index(f'Atom_ID_{dim_id_2}')] =\
                                loop.data[idx][loop.tags.index(f'Auth_atom_ID_{dim_id_2}')] = common_atom_id2

                        if None not in (shift, shift2):
                            diff = ((position - shift) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                            diff *= 2.0

                            diff_ = diff2_ = None
                            if shift_ is not None:
                                diff_ = ((position - shift_) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                            if shift2_ is not None:
                                diff2_ = ((position - shift) * weight) ** 2 + ((position2 - shift2_) * weight2) ** 2

                            if diff_ is not None and diff2_ is not None:
                                if diff_ < diff and diff2_ < diff:
                                    if diff_ < diff2_ and len_atom_id_ > 0:
                                        swap_atom_id_1(_atom_id_)
                                        continue
                                    if diff_ > diff2_ and len_atom_id2_ > 0:
                                        swap_atom_id_2(_atom_id2_)
                                        continue
                                elif diff_ < diff and len_atom_id_ > 0:
                                    swap_atom_id_1(_atom_id_)
                                    continue
                                elif diff2_ < diff and len_atom_id2_ > 0:
                                    swap_atom_id_2(_atom_id2_)
                                    continue
                                else:
                                    if diff_ < diff2_ and diff_ < 1.0 and len_atom_id_ > 0:
                                        swap_atom_id_1(_atom_id_)
                                        continue
                                    if diff_ > diff2_ and diff2_ < 1.0 and len_atom_id2_ > 0:
                                        swap_atom_id_2(_atom_id2_)
                                        continue

                            elif diff_ is not None and (diff_ < diff or diff_ < 1.0) and len_atom_id_ > 0:
                                swap_atom_id_1(_atom_id_)
                                continue

                            elif diff2_ is not None and (diff2_ < diff or diff2_ < 1.0) and len_atom_id2_ > 0:
                                swap_atom_id_2(_atom_id2_)
                                continue

                            if is_reparsable:
                                self.reasonsForReParsing['onebond_idx_history'] = self.onebond_idx_history
                                is_reparsable = False

                        if len_atom_id_ > 0 and len_atom_id2_ > 0\
                           and atom_id[0] != _common_atom_id_[0] and atom_id2[0] == _common_atom_id2_[0]:
                            swap_atom_id_2(_atom_id2_)
                        elif len_atom_id_ > 0 and len_atom_id2_ > 0\
                                and atom_id[0] == _common_atom_id_[0] and atom_id2[0] != _common_atom_id2_[0]:
                            swap_atom_id_1(_atom_id_)
                        elif 0 < len_atom_id2_ < len_atom_id_:
                            swap_atom_id_2(_atom_id2_)
                        elif 0 < len_atom_id_ < len_atom_id2_:
                            swap_atom_id_1(_atom_id_)
                        elif len(atom_id2) < len(atom_id) and len_atom_id2_ > 0:
                            swap_atom_id_2(_atom_id2_)
                        elif len(atom_id) < len(atom_id2) and len_atom_id_ > 0:
                            swap_atom_id_1(_atom_id_)
                        else:

                            if is_reparsable:
                                self.reasonsForReParsing['onebond_idx_history'] = self.onebond_idx_history
                                is_reparsable = False

                            if _atom_id2[0] in PROTON_BEGIN_CODE and len_atom_id2_ > 0:
                                swap_atom_id_2(_atom_id2_)
                            elif _atom_id[0] in PROTON_BEGIN_CODE and len_atom_id_ > 0:
                                swap_atom_id_1(_atom_id_)
                            else:
                                self.f.append("[Inconsistent peak assignment] "
                                              f"[Check row of Index_ID {loop.data[idx][loop.tags.index('Index_ID')]}] "
                                              "Inconsistent assignments of spectral peak with onebond coherence transfer type, "
                                              f"({chain_id}:{seq_id}:{comp_id}:{atom_id}) vs "
                                              f"({chain_id}:{seq_id}:{comp_id}:{atom_id2}) have been cleared.")

                                if loop.data[idx][details_col] in EMPTY_VALUE:
                                    loop.data[idx][details_col] = f'{atom_id}, {atom_id2} -> cleared'

                                loop.data[idx][loop.tags.index(f'Atom_ID_{dim_id_1}')] =\
                                    loop.data[idx][loop.tags.index(f'Auth_atom_ID_{dim_id_1}')] =\
                                    loop.data[idx][loop.tags.index(f'Atom_ID_{dim_id_2}')] =\
                                    loop.data[idx][loop.tags.index(f'Auth_atom_ID_{dim_id_2}')] = None

                    elif chain_id == chain_id2 and seq_id != seq_id2:
                        position, position2 = row[4], row[9]

                        if isinstance(position, str):
                            position = float(position)

                        if isinstance(position2, str):
                            position2 = float(position2)

                        shift, weight = self.__getCsValue(chain_id, seq_id, comp_id, atom_id)
                        shift2, weight2 = self.__getCsValue(chain_id, seq_id2, comp_id2, atom_id2)

                        shift_, _ = self.__getCsValue(chain_id, seq_id2, comp_id2, atom_id)
                        shift2_, _ = self.__getCsValue(chain_id, seq_id, comp_id, atom_id2)

                        if shift_ is not None and not self.ccU.hasBond(comp_id2, atom_id, atom_id2):
                            shift_ = None

                        if shift2_ is not None and not self.ccU.hasBond(comp_id, atom_id, atom_id2):
                            shift2_ = None

                        # pylint: disable=cell-var-from-loop
                        def swap_seq_id_1():
                            if loop.data[idx][details_col] in EMPTY_VALUE:
                                loop.data[idx][details_col] = f'{seq_id}:{comp_id} -> {seq_id2}:{comp_id2}'
                            loop.data[idx][loop.tags.index(f'Comp_index_ID_{dim_id_1}')] =\
                                loop.data[idx][loop.tags.index(f'Comp_index_ID_{dim_id_2}')]
                            loop.data[idx][loop.tags.index(f'Comp_ID_{dim_id_1}')] =\
                                loop.data[idx][loop.tags.index(f'Comp_ID_{dim_id_2}')]
                            loop.data[idx][loop.tags.index(f'Auth_seq_ID_{dim_id_1}')] =\
                                loop.data[idx][loop.tags.index(f'Auth_seq_ID_{dim_id_2}')]
                            loop.data[idx][loop.tags.index(f'Auth_comp_ID_{dim_id_1}')] =\
                                loop.data[idx][loop.tags.index(f'Auth_comp_ID_{dim_id_2}')]

                        # pylint: disable=cell-var-from-loop
                        def swap_seq_id_2():
                            if loop.data[idx][details_col] in EMPTY_VALUE:
                                loop.data[idx][details_col] = f'{seq_id2}:{comp_id2} -> {seq_id}:{comp_id}'
                            loop.data[idx][loop.tags.index(f'Comp_index_ID_{dim_id_2}')] =\
                                loop.data[idx][loop.tags.index(f'Comp_index_ID_{dim_id_1}')]
                            loop.data[idx][loop.tags.index(f'Comp_ID_{dim_id_2}')] =\
                                loop.data[idx][loop.tags.index(f'Comp_ID_{dim_id_1}')]
                            loop.data[idx][loop.tags.index(f'Auth_seq_ID_{dim_id_2}')] =\
                                loop.data[idx][loop.tags.index(f'Auth_seq_ID_{dim_id_1}')]
                            loop.data[idx][loop.tags.index(f'Auth_comp_ID_{dim_id_2}')] =\
                                loop.data[idx][loop.tags.index(f'Auth_comp_ID_{dim_id_1}')]

                        if None in (shift, shift2):

                            if shift is None and shift2 is None:

                                if is_reparsable:
                                    self.reasonsForReParsing['onebond_idx_history'] = self.onebond_idx_history
                                    is_reparsable = False

                                continue

                            if shift is None and shift_ is not None:
                                swap_seq_id_1()
                                continue

                            if shift2 is None and shift2_ is not None:
                                swap_seq_id_2()
                                continue

                            if shift is None and shift_ is None and None not in (shift2, shift2_):
                                swap_seq_id_1()
                                continue

                            if shift2 is None and shift2_ is None and None not in (shift, shift_):
                                swap_seq_id_2()
                                continue

                        else:

                            diff = ((position - shift) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                            diff *= 2.0

                            diff_ = diff2_ = None
                            if shift_ is not None:
                                diff_ = ((position - shift_) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                            if shift2_ is not None:
                                diff2_ = ((position - shift) * weight) ** 2 + ((position2 - shift2_) * weight2) ** 2

                            if diff_ is not None and diff2_ is not None:
                                if diff_ < diff and diff2_ < diff:
                                    if diff_ < diff2_:
                                        swap_seq_id_1()
                                        continue
                                    if diff_ > diff2_:
                                        swap_seq_id_2()
                                        continue
                                elif diff_ < diff:
                                    swap_seq_id_1()
                                    continue
                                elif diff2_ < diff:
                                    swap_seq_id_2()
                                    continue
                                else:
                                    if diff_ < diff2_ and diff_ < 1.0:
                                        swap_seq_id_1()
                                        continue
                                    if diff_ > diff2_ and diff2_ < 1.0:
                                        swap_seq_id_2()
                                        continue

                                if is_reparsable:
                                    self.reasonsForReParsing['onebond_idx_history'] = self.onebond_idx_history
                                    is_reparsable = False

                                continue

                            if diff_ is not None and (diff_ < diff or diff_ < 1.0):
                                swap_seq_id_1()
                                continue

                            if diff2_ is not None and (diff2_ < diff or diff2_ < 1.0):
                                swap_seq_id_2()
                                continue

                            _atom_ids = self.nefT.get_valid_star_atom(comp_id, atom_id, leave_unmatched=False)[0]
                            _atom_ids2 = self.nefT.get_valid_star_atom(comp_id2, atom_id2, leave_unmatched=False)[0]

                            _atom_id, _atom_id2 = _atom_ids[0], _atom_ids2[0]

                            _atom_id2_ = self.ccU.getBondedAtoms(comp_id, _atom_id,
                                                                 exclProton=_atom_id[0] in PROTON_BEGIN_CODE,
                                                                 onlyProton=_atom_id[0] not in PROTON_BEGIN_CODE)
                            _atom_id_ = self.ccU.getBondedAtoms(comp_id2, _atom_id2,
                                                                exclProton=_atom_id2[0] in PROTON_BEGIN_CODE,
                                                                onlyProton=_atom_id2[0] not in PROTON_BEGIN_CODE)

                            len_atom_id_ = len(_atom_id_)
                            len_atom_id2_ = len(_atom_id2_)

                            _common_atom_id_ = _common_atom_id2_ = shift_ = shift2_ = None
                            if len_atom_id_ > 0:
                                _common_atom_id_ = _atom_id_[0]
                                if _common_atom_id_[0] == _atom_id[0]:
                                    shift_, _ = self.__getCsValue(chain_id, seq_id2, comp_id2, _common_atom_id_)
                            if len_atom_id2_ > 0:
                                _common_atom_id2_ = _atom_id2_[0]
                                if _common_atom_id2_[0] == _atom_id2[0]:
                                    shift2_, _ = self.__getCsValue(chain_id, seq_id, comp_id, _common_atom_id2_)

                            # pylint: disable=cell-var-from-loop
                            def swap_seq_atom_id_1(atom_id_list):
                                common_atom_id = atom_id_list[0]
                                if len(atom_id_list) > 1 or len(_atom_ids2) > 1:
                                    atom_sel = [{'atom_id': _a_} for _a_ in atom_id_list]
                                    if len(_atom_ids2) > 1:
                                        for _atom_id2 in _atom_ids2[1:]:
                                            __atom_id__ = self.ccU.getBondedAtoms(comp_id, _atom_id2,
                                                                                  exclProton=_atom_id2[0] in PROTON_BEGIN_CODE,
                                                                                  onlyProton=_atom_id2[0] not in PROTON_BEGIN_CODE)
                                            for __a__ in __atom_id__:
                                                atom_sel.append({'atom_id': __a__})
                                    try:
                                        common_atom_id = self.__extractCommonAtom(atom_sel)['auth_atom_id']
                                    except KeyError:
                                        pass
                                if loop.data[idx][details_col] in EMPTY_VALUE:
                                    loop.data[idx][details_col] = f'{seq_id}:{comp_id}:{atom_id} -> '\
                                        f'{seq_id2}:{comp_id2}:{common_atom_id}'
                                loop.data[idx][loop.tags.index(f'Comp_index_ID_{dim_id_1}')] =\
                                    loop.data[idx][loop.tags.index(f'Comp_index_ID_{dim_id_2}')]
                                loop.data[idx][loop.tags.index(f'Comp_ID_{dim_id_1}')] =\
                                    loop.data[idx][loop.tags.index(f'Comp_ID_{dim_id_2}')]
                                loop.data[idx][loop.tags.index(f'Auth_seq_ID_{dim_id_1}')] =\
                                    loop.data[idx][loop.tags.index(f'Auth_seq_ID_{dim_id_2}')]
                                loop.data[idx][loop.tags.index(f'Auth_comp_ID_{dim_id_1}')] =\
                                    loop.data[idx][loop.tags.index(f'Auth_comp_ID_{dim_id_2}')]
                                loop.data[idx][loop.tags.index(f'Atom_ID_{dim_id_1}')] =\
                                    loop.data[idx][loop.tags.index(f'Auth_atom_ID_{dim_id_1}')] = common_atom_id

                            # pylint: disable=cell-var-from-loop
                            def swap_seq_atom_id_2(atom_id2_list):
                                common_atom_id2 = atom_id2_list[0]
                                if len(atom_id2_list) > 1 or len(_atom_ids) > 1:
                                    atom_sel = [{'atom_id': _a_} for _a_ in atom_id2_list]
                                    if len(_atom_ids) > 1:
                                        for _atom_id in _atom_ids[1:]:
                                            __atom_id2__ = self.ccU.getBondedAtoms(comp_id, _atom_id,
                                                                                   exclProton=_atom_id[0] in PROTON_BEGIN_CODE,
                                                                                   onlyProton=_atom_id[0] not in PROTON_BEGIN_CODE)
                                            for __a__ in __atom_id2__:
                                                atom_sel.append({'atom_id': __a__})
                                    try:
                                        common_atom_id2 = self.__extractCommonAtom(atom_sel)['auth_atom_id']
                                    except KeyError:
                                        pass
                                if loop.data[idx][details_col] in EMPTY_VALUE:
                                    loop.data[idx][details_col] = f'{seq_id2}:{comp_id2}:{atom_id2} -> '\
                                        f'{seq_id}:{comp_id}:{common_atom_id2}'
                                loop.data[idx][loop.tags.index(f'Comp_index_ID_{dim_id_2}')] =\
                                    loop.data[idx][loop.tags.index(f'Comp_index_ID_{dim_id_1}')]
                                loop.data[idx][loop.tags.index(f'Comp_ID_{dim_id_2}')] =\
                                    loop.data[idx][loop.tags.index(f'Comp_ID_{dim_id_1}')]
                                loop.data[idx][loop.tags.index(f'Auth_seq_ID_{dim_id_2}')] =\
                                    loop.data[idx][loop.tags.index(f'Auth_seq_ID_{dim_id_1}')]
                                loop.data[idx][loop.tags.index(f'Auth_comp_ID_{dim_id_2}')] =\
                                    loop.data[idx][loop.tags.index(f'Auth_comp_ID_{dim_id_1}')]
                                loop.data[idx][loop.tags.index(f'Atom_ID_{dim_id_2}')] =\
                                    loop.data[idx][loop.tags.index(f'Auth_atom_ID_{dim_id_2}')] = common_atom_id2

                            diff_ = diff2_ = None
                            if shift_ is not None:
                                diff_ = ((position - shift_) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                            if shift2_ is not None:
                                diff2_ = ((position - shift) * weight) ** 2 + ((position2 - shift2_) * weight2) ** 2

                            if diff_ is not None and diff2_ is not None:
                                if diff_ < diff and diff2_ < diff:
                                    if diff_ < diff2_ and len_atom_id_ > 0:
                                        swap_seq_atom_id_1(_atom_id_)
                                        continue
                                    if diff_ > diff2_ and len_atom_id2_ > 0:
                                        swap_seq_atom_id_2(_atom_id2_)
                                        continue
                                elif diff_ < diff and len_atom_id_ > 0:
                                    swap_seq_atom_id_1(_atom_id_)
                                    continue
                                elif diff2_ < diff and len_atom_id2_ > 0:
                                    swap_seq_atom_id_2(_atom_id2_)
                                    continue
                                else:
                                    if diff_ < diff2_ and diff_ < 1.0 and len_atom_id_ > 0:
                                        swap_seq_atom_id_1(_atom_id_)
                                        continue
                                    if diff_ > diff2_ and diff2_ < 1.0 and len_atom_id2_ > 0:
                                        swap_seq_atom_id_2(_atom_id2_)
                                        continue

                            elif diff_ is not None and (diff_ < diff or diff_ < 1.0) and len_atom_id_ > 0:
                                swap_seq_atom_id_1(_atom_id_)
                                continue

                            elif diff2_ is not None and (diff2_ < diff or diff2_ < 1.0) and len_atom_id2_ > 0:
                                swap_seq_atom_id_2(_atom_id2_)
                                continue

                            if is_reparsable:
                                self.reasonsForReParsing['onebond_idx_history'] = self.onebond_idx_history
                                is_reparsable = False

                    elif chain_id != chain_id2:
                        position, position2 = row[4], row[9]

                        if isinstance(position, str):
                            position = float(position)

                        if isinstance(position2, str):
                            position2 = float(position2)

                        shift, weight = self.__getCsValue(chain_id, seq_id, comp_id, atom_id)
                        shift2, weight2 = self.__getCsValue(chain_id2, seq_id2, comp_id2, atom_id2)

                        shift_, _ = self.__getCsValue(chain_id2, seq_id2, comp_id2, atom_id)
                        shift2_, _ = self.__getCsValue(chain_id, seq_id, comp_id, atom_id2)

                        if shift_ is not None and not self.ccU.hasBond(comp_id2, atom_id, atom_id2):
                            shift_ = None

                        if shift2_ is not None and not self.ccU.hasBond(comp_id, atom_id, atom_id2):
                            shift2_ = None

                        # pylint: disable=cell-var-from-loop
                        def swap_chain_seq_id_1():
                            if loop.data[idx][details_col] in EMPTY_VALUE:
                                loop.data[idx][details_col] = f'{chain_id}:{seq_id}:{comp_id} -> {chain_id2}:{seq_id2}:{comp_id2}'
                            loop.data[idx][loop.tags.index(f'Entity_assembly_ID_{dim_id_1}')] =\
                                loop.data[idx][loop.tags.index(f'Entity_assembly_ID_{dim_id_2}')]
                            loop.data[idx][loop.tags.index(f'Entity_ID_{dim_id_1}')] =\
                                loop.data[idx][loop.tags.index(f'Entity_ID_{dim_id_2}')]
                            loop.data[idx][loop.tags.index(f'Comp_index_ID_{dim_id_1}')] =\
                                loop.data[idx][loop.tags.index(f'Comp_index_ID_{dim_id_2}')]
                            loop.data[idx][loop.tags.index(f'Comp_ID_{dim_id_1}')] =\
                                loop.data[idx][loop.tags.index(f'Comp_ID_{dim_id_2}')]
                            loop.data[idx][loop.tags.index(f'Auth_asym_ID_{dim_id_1}')] =\
                                loop.data[idx][loop.tags.index(f'Auth_asym_ID_{dim_id_2}')]
                            loop.data[idx][loop.tags.index(f'Auth_seq_ID_{dim_id_1}')] =\
                                loop.data[idx][loop.tags.index(f'Auth_seq_ID_{dim_id_2}')]
                            loop.data[idx][loop.tags.index(f'Auth_comp_ID_{dim_id_1}')] =\
                                loop.data[idx][loop.tags.index(f'Auth_comp_ID_{dim_id_2}')]

                        # pylint: disable=cell-var-from-loop
                        def swap_chain_seq_id_2():
                            if loop.data[idx][details_col] in EMPTY_VALUE:
                                loop.data[idx][details_col] = f'{chain_id2}:{seq_id2}:{comp_id2} -> {chain_id}:{seq_id}:{comp_id}'
                            loop.data[idx][loop.tags.index(f'Entity_assembly_ID_{dim_id_2}')] =\
                                loop.data[idx][loop.tags.index(f'Entity_assembly_ID_{dim_id_1}')]
                            loop.data[idx][loop.tags.index(f'Entity_ID_{dim_id_2}')] =\
                                loop.data[idx][loop.tags.index(f'Entity_ID_{dim_id_1}')]
                            loop.data[idx][loop.tags.index(f'Comp_index_ID_{dim_id_2}')] =\
                                loop.data[idx][loop.tags.index(f'Comp_index_ID_{dim_id_1}')]
                            loop.data[idx][loop.tags.index(f'Comp_ID_{dim_id_2}')] =\
                                loop.data[idx][loop.tags.index(f'Comp_ID_{dim_id_1}')]
                            loop.data[idx][loop.tags.index(f'Auth_asym_ID_{dim_id_2}')] =\
                                loop.data[idx][loop.tags.index(f'Auth_asym_ID_{dim_id_1}')]
                            loop.data[idx][loop.tags.index(f'Auth_seq_ID_{dim_id_2}')] =\
                                loop.data[idx][loop.tags.index(f'Auth_seq_ID_{dim_id_1}')]
                            loop.data[idx][loop.tags.index(f'Auth_comp_ID_{dim_id_2}')] =\
                                loop.data[idx][loop.tags.index(f'Auth_comp_ID_{dim_id_1}')]

                        if None in (shift, shift2):

                            if shift is None and shift2 is None:

                                if is_reparsable:
                                    self.reasonsForReParsing['onebond_idx_history'] = self.onebond_idx_history
                                    is_reparsable = False

                                continue

                            if shift is None and shift_ is not None:
                                swap_chain_seq_id_1()
                                continue

                            if shift2 is None and shift2_ is not None:
                                swap_chain_seq_id_2()
                                continue

                            if shift is None and shift_ is None and None not in (shift2, shift2_):
                                swap_chain_seq_id_1()
                                continue

                            if shift2 is None and shift2_ is None and None not in (shift, shift_):
                                swap_chain_seq_id_2()
                                continue

                            if is_reparsable:
                                self.reasonsForReParsing['onebond_idx_history'] = self.onebond_idx_history
                                is_reparsable = False

                            continue

                        diff = ((position - shift) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                        diff *= 2.0

                        diff_ = diff2_ = None
                        if shift_ is not None:
                            diff_ = ((position - shift_) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                        if shift2_ is not None:
                            diff2_ = ((position - shift) * weight) ** 2 + ((position2 - shift2_) * weight2) ** 2

                        if diff_ is not None and diff2_ is not None:
                            if diff_ < diff and diff2_ < diff:
                                if diff_ < diff2_:
                                    swap_chain_seq_id_1()
                                    continue
                                if diff_ > diff2_:
                                    swap_chain_seq_id_2()
                                    continue
                            elif diff_ < diff:
                                swap_chain_seq_id_1()
                                continue
                            elif diff2_ < diff:
                                swap_chain_seq_id_2()
                                continue
                            else:
                                if diff_ < diff2_ and diff_ < 1.0:
                                    swap_chain_seq_id_1()
                                    continue
                                if diff_ > diff2_ and diff2_ < 1.0:
                                    swap_chain_seq_id_2()
                                    continue

                        elif diff_ is not None and (diff_ < diff or diff_ < 1.0):
                            swap_chain_seq_id_1()
                            continue

                        elif diff2_ is not None and (diff2_ < diff or diff2_ < 1.0):
                            swap_chain_seq_id_2()
                            continue

                        if is_reparsable:
                            self.reasonsForReParsing['onebond_idx_history'] = self.onebond_idx_history
                            is_reparsable = False

            else:

                tags = ['Peak_ID', 'Spectral_dim_ID', 'Entity_assembly_ID', 'Comp_index_ID', 'Comp_ID', 'Atom_ID', 'Val']

                chain_id_col = loop.tags.index('Entity_assembly_ID')
                entity_id_col = loop.tags.index('Entity_ID')
                seq_id_col = loop.tags.index('Comp_index_ID')
                comp_id_col = loop.tags.index('Comp_ID')
                atom_id_col = loop.tags.index('Atom_ID')
                auth_chain_id_col = loop.tags.index('Auth_entity_ID')
                auth_seq_id_col = loop.tags.index('Auth_seq_ID')
                auth_comp_id_col = loop.tags.index('Auth_comp_ID')
                auth_atom_id_col = loop.tags.index('Auth_atom_ID')

                dat = loop.get_tag(tags)

                peak_id = None

                for idx, row in enumerate(dat):
                    dim_id = row[1]

                    if any(True for col in range(7) if row[col] in EMPTY_VALUE):
                        continue

                    if dim_id == 1:
                        peak_id = row[0]
                        chain_ids, seq_ids, comp_ids, atom_ids, positions = [], [], [], [], []
                    else:
                        if peak_id != row[0]:
                            continue

                    chain_ids.append(row[2] if isinstance(row[2], str) else str(row[2]))  # noqa: E501, pylint: disable=possibly-used-before-assignment,line-too-long
                    seq_ids.append(int(row[3]) if isinstance(row[3], str) else row[3])  # noqa: E501, pylint: disable=possibly-used-before-assignment,line-too-long
                    comp_ids.append(row[4])  # pylint: disable=possibly-used-before-assignment
                    atom_ids.append(row[5])  # pylint: disable=possibly-used-before-assignment
                    positions.append(float(row[6]) if isinstance(row[6], str) else row[6])  # noqa: E501, pylint: disable=possibly-used-before-assignment,line-too-long

                    if dim_id < num_of_dim:
                        continue

                    if len(atom_ids) < num_of_dim:
                        continue

                    _dim_id_1 = dim_id_1 - 1
                    _dim_id_2 = dim_id_2 - 1

                    chain_id, seq_id, comp_id, atom_id, chain_id2, seq_id2, comp_id2, atom_id2 =\
                        chain_ids[_dim_id_1], seq_ids[_dim_id_1], comp_ids[_dim_id_1], atom_ids[_dim_id_1], \
                        chain_ids[_dim_id_2], seq_ids[_dim_id_2], comp_ids[_dim_id_2], atom_ids[_dim_id_2]

                    if chain_id == chain_id2 and seq_id == seq_id2 and atom_id != atom_id2\
                       and self.ccU.updateChemCompDict(comp_id):
                        _atom_ids = self.nefT.get_valid_star_atom(comp_id, atom_id, leave_unmatched=False)[0]
                        _atom_ids2 = self.nefT.get_valid_star_atom(comp_id, atom_id2, leave_unmatched=False)[0]
                        if any(True for b in self.ccU.lastBondDictList
                           if ((b['atom_id_1'] in _atom_ids and b['atom_id_2'] in _atom_ids2)
                               or (b['atom_id_1'] in _atom_ids2 and b['atom_id_2'] in _atom_ids))):
                            continue

                        if self.software_name == 'PIPP':
                            del_idx_list.extend(list(range(idx - num_of_dim + 1, idx + 1)))
                            continue

                        _atom_id, _atom_id2 = _atom_ids[0], _atom_ids2[0]

                        _atom_id2_ = self.ccU.getBondedAtoms(comp_id, _atom_id,
                                                             exclProton=_atom_id[0] in PROTON_BEGIN_CODE,
                                                             onlyProton=_atom_id[0] not in PROTON_BEGIN_CODE)
                        _atom_id_ = self.ccU.getBondedAtoms(comp_id, _atom_id2,
                                                            exclProton=_atom_id2[0] in PROTON_BEGIN_CODE,
                                                            onlyProton=_atom_id2[0] not in PROTON_BEGIN_CODE)

                        len_atom_id_ = len(_atom_id_)
                        len_atom_id2_ = len(_atom_id2_)

                        position, position2 = positions[_dim_id_1], positions[_dim_id_2]

                        shift, weight = self.__getCsValue(chain_id, seq_id, comp_id, _atom_id)
                        shift2, weight2 = self.__getCsValue(chain_id, seq_id, comp_id, _atom_id2)

                        _common_atom_id_ = _common_atom_id2_ = shift_ = shift2_ = None
                        if len_atom_id_ > 0:
                            _common_atom_id_ = _atom_id_[0]
                            if _common_atom_id_[0] == _atom_id[0]:
                                shift_, _ = self.__getCsValue(chain_id, seq_id, comp_id, _common_atom_id_)
                        if len_atom_id2_ > 0:
                            _common_atom_id2_ = _atom_id2_[0]
                            if _common_atom_id2_[0] == _atom_id2[0]:
                                shift2_, _ = self.__getCsValue(chain_id, seq_id, comp_id, _common_atom_id2_)

                        if shift_ is not None and not self.ccU.hasBond(comp_id, _common_atom_id_, _atom_id2):
                            shift_ = None

                        if shift2_ is not None and not self.ccU.hasBond(comp_id, _atom_id, _common_atom_id2_):
                            shift2_ = None

                        # pylint: disable=cell-var-from-loop
                        def alt_swap_atom_id_1(atom_id_list):
                            common_atom_id = atom_id_list[0]
                            if len(atom_id_list) > 1 or len(_atom_ids2) > 1:
                                atom_sel = [{'atom_id': _a_} for _a_ in atom_id_list]
                                if len(_atom_ids2) > 1:
                                    for _atom_id2 in _atom_ids2[1:]:
                                        __atom_id__ = self.ccU.getBondedAtoms(comp_id, _atom_id2,
                                                                              exclProton=_atom_id2[0] in PROTON_BEGIN_CODE,
                                                                              onlyProton=_atom_id2[0] not in PROTON_BEGIN_CODE)
                                        for __a__ in __atom_id__:
                                            atom_sel.append({'atom_id': __a__})
                                try:
                                    common_atom_id = self.__extractCommonAtom(atom_sel)['auth_atom_id']
                                except KeyError:
                                    pass
                            if loop.data[idx - num_of_dim + dim_id_1][details_col] in EMPTY_VALUE:
                                loop.data[idx - num_of_dim + dim_id_1][details_col] = f'{atom_id} -> {_common_atom_id_}'
                            loop.data[idx - num_of_dim + dim_id_1][atom_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_1][auth_atom_id_col] = common_atom_id

                        # pylint: disable=cell-var-from-loop
                        def alt_swap_atom_id_2(atom_id2_list):
                            common_atom_id2 = atom_id2_list[0]
                            if len(atom_id2_list) > 1 or len(_atom_ids) > 1:
                                atom_sel = [{'atom_id': _a_} for _a_ in atom_id2_list]
                                if len(_atom_ids) > 1:
                                    for _atom_id in _atom_ids[1:]:
                                        __atom_id2__ = self.ccU.getBondedAtoms(comp_id, _atom_id,
                                                                               exclProton=_atom_id[0] in PROTON_BEGIN_CODE,
                                                                               onlyProton=_atom_id[0] not in PROTON_BEGIN_CODE)
                                        for __a__ in __atom_id2__:
                                            atom_sel.append({'atom_id': __a__})
                                try:
                                    common_atom_id2 = self.__extractCommonAtom(atom_sel)['auth_atom_id']
                                except KeyError:
                                    pass
                            if loop.data[idx - num_of_dim + dim_id_2][details_col] in EMPTY_VALUE:
                                loop.data[idx - num_of_dim + dim_id_2][details_col] = f'{atom_id2} -> {common_atom_id2}'
                            loop.data[idx - num_of_dim + dim_id_2][atom_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_2][auth_atom_id_col] = common_atom_id2

                        if None not in (shift, shift2):
                            diff = ((position - shift) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                            diff *= 2.0

                            diff_ = diff2_ = None
                            if shift_ is not None:
                                diff_ = ((position - shift_) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                            if shift2_ is not None:
                                diff2_ = ((position - shift) * weight) ** 2 + ((position2 - shift2_) * weight2) ** 2

                            if diff_ is not None and diff2_ is not None:
                                if diff_ < diff and diff2_ < diff:
                                    if diff_ < diff2_ and len_atom_id_ > 0:
                                        alt_swap_atom_id_1(_atom_id_)
                                        continue
                                    if diff_ > diff2_ and len_atom_id2_ > 0:
                                        alt_swap_atom_id_2(_atom_id2_)
                                        continue
                                elif diff_ < diff and len_atom_id_ > 0:
                                    alt_swap_atom_id_1(_atom_id_)
                                    continue
                                elif diff2_ < diff and len_atom_id2_ > 0:
                                    alt_swap_atom_id_2(_atom_id2_)
                                    continue
                                else:
                                    if diff_ < diff2_ and diff_ < 1.0 and len_atom_id_ > 0:
                                        alt_swap_atom_id_1(_atom_id_)
                                        continue
                                    if diff_ > diff2_ and diff2_ < 1.0 and len_atom_id2_ > 0:
                                        alt_swap_atom_id_2(_atom_id2_)
                                        continue

                                if is_reparsable:
                                    self.reasonsForReParsing['onebond_idx_history'] = self.onebond_idx_history
                                    is_reparsable = False

                                continue

                            if diff_ is not None and (diff_ < diff or diff_ < 1.0) and len_atom_id_ > 0:
                                alt_swap_atom_id_1(_atom_id_)
                                continue

                            if diff2_ is not None and (diff2_ < diff or diff2_ < 1.0) and len_atom_id2_ > 0:
                                alt_swap_atom_id_2(_atom_id2_)
                                continue

                        if len_atom_id_ > 0 and len_atom_id2_ > 0\
                           and atom_id[0] != _common_atom_id_[0] and atom_id2[0] == _common_atom_id2_[0]:
                            alt_swap_atom_id_2(_atom_id2_)
                        elif len_atom_id_ > 0 and len_atom_id2_ > 0\
                                and atom_id[0] == _common_atom_id_[0] and atom_id2[0] != _common_atom_id2_[0]:
                            alt_swap_atom_id_1(_atom_id_)
                        elif 0 < len_atom_id2_ < len_atom_id_:
                            alt_swap_atom_id_2(_atom_id2_)
                        elif 0 < len_atom_id_ < len_atom_id2_:
                            alt_swap_atom_id_1(_atom_id_)
                        elif len(atom_id2) < len(atom_id) and len_atom_id2_ > 0:
                            alt_swap_atom_id_2(_atom_id2_)
                        elif len(atom_id) < len(atom_id2) and len_atom_id_ > 0:
                            alt_swap_atom_id_1(_atom_id_)
                        else:

                            if is_reparsable:
                                self.reasonsForReParsing['onebond_idx_history'] = self.onebond_idx_history
                                is_reparsable = False

                            if _atom_id2[0] in PROTON_BEGIN_CODE and len_atom_id2_ > 0:
                                alt_swap_atom_id_2(_atom_id2_)
                            elif _atom_id[0] in PROTON_BEGIN_CODE and len_atom_id_ > 0:
                                alt_swap_atom_id_1(_atom_id_)
                            else:
                                self.f.append("[Inconsistent peak assignment] "
                                              f"[Check row of Peak_ID {loop.data[idx][loop.tags.index('Peak_ID')]}] "
                                              "Inconsistent assignments of spectral peak with onebond coherence transfer type, "
                                              f"({chain_id}:{seq_id}:{comp_id}:{atom_id}) vs "
                                              f"({chain_id}:{seq_id}:{comp_id}:{atom_id2}) have been cleared.")

                                if loop.data[idx - num_of_dim + dim_id_1][details_col] in EMPTY_VALUE:
                                    loop.data[idx - num_of_dim + dim_id_1][details_col] = f'{atom_id} -> cleared'
                                if loop.data[idx - num_of_dim + dim_id_2][details_col] in EMPTY_VALUE:
                                    loop.data[idx - num_of_dim + dim_id_2][details_col] = f'{atom_id2} -> cleared'

                                loop.data[idx - num_of_dim + dim_id_1][atom_id_col] =\
                                    loop.data[idx - num_of_dim + dim_id_1][auth_atom_id_col] =\
                                    loop.data[idx - num_of_dim + dim_id_2][atom_id_col] =\
                                    loop.data[idx - num_of_dim + dim_id_2][auth_atom_id_col] = None

                    elif chain_id == chain_id2 and seq_id != seq_id2:
                        position, position2 = positions[_dim_id_1], positions[_dim_id_2]

                        shift, weight = self.__getCsValue(chain_id, seq_id, comp_id, atom_id)
                        shift2, weight2 = self.__getCsValue(chain_id, seq_id2, comp_id2, atom_id2)

                        shift_, _ = self.__getCsValue(chain_id, seq_id2, comp_id2, atom_id)
                        shift2_, _ = self.__getCsValue(chain_id, seq_id, comp_id, atom_id2)

                        if shift_ is not None and not self.ccU.hasBond(comp_id2, atom_id, atom_id2):
                            shift_ = None

                        if shift2_ is not None and not self.ccU.hasBond(comp_id, atom_id, atom_id2):
                            shift2_ = None

                        # pylint: disable=cell-var-from-loop
                        def alt_swap_seq_id_1():
                            if loop.data[idx][details_col] in EMPTY_VALUE:
                                loop.data[idx - num_of_dim + dim_id_1][details_col] =\
                                    f'{seq_id}:{comp_id} -> {seq_id2}:{comp_id2}'
                            loop.data[idx - num_of_dim + dim_id_1][seq_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_2][seq_id_col]
                            loop.data[idx - num_of_dim + dim_id_1][comp_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_2][comp_id_col]
                            loop.data[idx - num_of_dim + dim_id_1][auth_seq_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_2][auth_seq_id_col]
                            loop.data[idx - num_of_dim + dim_id_1][auth_comp_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_2][auth_comp_id_col]

                        # pylint: disable=cell-var-from-loop
                        def alt_swap_seq_id_2():
                            if loop.data[idx - num_of_dim + dim_id_2][details_col] in EMPTY_VALUE:
                                loop.data[idx - num_of_dim + dim_id_2][details_col] =\
                                    f'{seq_id2}:{comp_id2} -> {seq_id}:{comp_id}'
                            loop.data[idx - num_of_dim + dim_id_2][seq_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_1][seq_id_col]
                            loop.data[idx - num_of_dim + dim_id_2][comp_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_1][comp_id_col]
                            loop.data[idx - num_of_dim + dim_id_2][auth_seq_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_1][auth_seq_id_col]
                            loop.data[idx - num_of_dim + dim_id_2][auth_comp_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_1][auth_comp_id_col]

                        if None in (shift, shift2):

                            if shift is None and shift2 is None:

                                if is_reparsable:
                                    self.reasonsForReParsing['onebond_idx_history'] = self.onebond_idx_history
                                    is_reparsable = False

                                else:
                                    del_idx_list.extend(list(range(idx - num_of_dim + 1, idx + 1)))

                                continue

                            if shift is None and shift_ is not None:
                                alt_swap_seq_id_1()
                                continue

                            if shift2 is None and shift2_ is not None:
                                alt_swap_seq_id_2()
                                continue

                            if shift is None and shift_ is None and None not in (shift2, shift2_):
                                alt_swap_seq_id_1()
                                continue

                            if shift2 is None and shift2_ is None and None not in (shift, shift_):
                                alt_swap_seq_id_2()
                                continue

                        else:

                            diff = ((position - shift) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                            diff *= 2.0

                            diff_ = diff2_ = None
                            if shift_ is not None:
                                diff_ = ((position - shift_) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                            if shift2_ is not None:
                                diff2_ = ((position - shift) * weight) ** 2 + ((position2 - shift2_) * weight2) ** 2

                            if diff_ is not None and diff2_ is not None:
                                if diff_ < diff and diff2_ < diff:
                                    if diff_ < diff2_:
                                        alt_swap_seq_id_1()
                                        continue
                                    if diff_ > diff2_:
                                        alt_swap_seq_id_2()
                                        continue
                                elif diff_ < diff:
                                    alt_swap_seq_id_1()
                                    continue
                                elif diff2_ < diff:
                                    alt_swap_seq_id_2()
                                    continue
                                else:
                                    if diff_ < diff2_ and diff_ < 1.0:
                                        alt_swap_seq_id_1()
                                        continue
                                    if diff_ > diff2_ and diff2_ < 1.0:
                                        alt_swap_seq_id_2()
                                        continue

                                if is_reparsable:
                                    self.reasonsForReParsing['onebond_idx_history'] = self.onebond_idx_history
                                    is_reparsable = False

                                else:
                                    del_idx_list.extend(list(range(idx - num_of_dim + 1, idx + 1)))

                                continue

                            if diff_ is not None and (diff_ < diff or diff_ < 1.0):
                                alt_swap_seq_id_1()
                                continue

                            if diff2_ is not None and (diff2_ < diff or diff2_ < 1.0):
                                alt_swap_seq_id_2()
                                continue

                            _atom_ids = self.nefT.get_valid_star_atom(comp_id, atom_id, leave_unmatched=False)[0]
                            _atom_ids2 = self.nefT.get_valid_star_atom(comp_id2, atom_id2, leave_unmatched=False)[0]

                            _atom_id, _atom_id2 = _atom_ids[0], _atom_ids2[0]

                            _atom_id2_ = self.ccU.getBondedAtoms(comp_id, _atom_id,
                                                                 exclProton=_atom_id[0] in PROTON_BEGIN_CODE,
                                                                 onlyProton=_atom_id[0] not in PROTON_BEGIN_CODE)
                            _atom_id_ = self.ccU.getBondedAtoms(comp_id2, _atom_id2,
                                                                exclProton=_atom_id2[0] in PROTON_BEGIN_CODE,
                                                                onlyProton=_atom_id2[0] not in PROTON_BEGIN_CODE)

                            len_atom_id_ = len(_atom_id_)
                            len_atom_id2_ = len(_atom_id2_)

                            _common_atom_id_ = _common_atom_id2_ = shift_ = shift2_ = None
                            if len_atom_id_ > 0:
                                _common_atom_id_ = _atom_id_[0]
                                if _common_atom_id_[0] == _atom_id[0]:
                                    shift_, _ = self.__getCsValue(chain_id, seq_id2, comp_id2, _common_atom_id_)
                            if len_atom_id2_ > 0:
                                _common_atom_id2_ = _atom_id2_[0]
                                if _common_atom_id2_[0] == _atom_id2[0]:
                                    shift2_, _ = self.__getCsValue(chain_id, seq_id, comp_id, _common_atom_id2_)

                            # pylint: disable=cell-var-from-loop
                            def alt_swap_seq_atom_id_1(atom_id_list):
                                common_atom_id = atom_id_list[0]
                                if len(atom_id_list) > 1 or len(_atom_ids2) > 1:
                                    atom_sel = [{'atom_id': _a_} for _a_ in atom_id_list]
                                    if len(_atom_ids2) > 1:
                                        for _atom_id2 in _atom_ids2[1:]:
                                            __atom_id__ = self.ccU.getBondedAtoms(comp_id, _atom_id2,
                                                                                  exclProton=_atom_id2[0] in PROTON_BEGIN_CODE,
                                                                                  onlyProton=_atom_id2[0] not in PROTON_BEGIN_CODE)
                                            for __a__ in __atom_id__:
                                                atom_sel.append({'atom_id': __a__})
                                    try:
                                        common_atom_id = self.__extractCommonAtom(atom_sel)['auth_atom_id']
                                    except KeyError:
                                        pass
                                if loop.data[idx][details_col] in EMPTY_VALUE:
                                    loop.data[idx - num_of_dim + dim_id_1][details_col] =\
                                        f'{seq_id}:{comp_id}:{atom_id} -> {seq_id2}:{comp_id2}:{common_atom_id}'
                                loop.data[idx - num_of_dim + dim_id_1][seq_id_col] =\
                                    loop.data[idx - num_of_dim + dim_id_2][seq_id_col]
                                loop.data[idx - num_of_dim + dim_id_1][comp_id_col] =\
                                    loop.data[idx - num_of_dim + dim_id_2][comp_id_col]
                                loop.data[idx - num_of_dim + dim_id_1][auth_seq_id_col] =\
                                    loop.data[idx - num_of_dim + dim_id_2][auth_seq_id_col]
                                loop.data[idx - num_of_dim + dim_id_1][auth_comp_id_col] =\
                                    loop.data[idx - num_of_dim + dim_id_2][auth_comp_id_col]
                                loop.data[idx - num_of_dim + dim_id_1][atom_id_col] = \
                                    loop.data[idx - num_of_dim + dim_id_1][auth_atom_id_col] = common_atom_id

                            # pylint: disable=cell-var-from-loop
                            def alt_swap_seq_atom_id_2(atom_id2_list):
                                common_atom_id2 = atom_id2_list[0]
                                if len(atom_id2_list) > 1 or len(_atom_ids) > 1:
                                    atom_sel = [{'atom_id': _a_} for _a_ in atom_id2_list]
                                    if len(_atom_ids) > 1:
                                        for _atom_id in _atom_ids[1:]:
                                            __atom_id2__ = self.ccU.getBondedAtoms(comp_id, _atom_id,
                                                                                   exclProton=_atom_id[0] in PROTON_BEGIN_CODE,
                                                                                   onlyProton=_atom_id[0] not in PROTON_BEGIN_CODE)
                                            for __a__ in __atom_id2__:
                                                atom_sel.append({'atom_id': __a__})
                                    try:
                                        common_atom_id2 = self.__extractCommonAtom(atom_sel)['auth_atom_id']
                                    except KeyError:
                                        pass
                                if loop.data[idx - num_of_dim + dim_id_2][details_col] in EMPTY_VALUE:
                                    loop.data[idx - num_of_dim + dim_id_2][details_col] =\
                                        f'{seq_id2}:{comp_id2}:{atom_id2} -> {seq_id}:{comp_id}:{common_atom_id2}'
                                loop.data[idx - num_of_dim + dim_id_2][seq_id_col] =\
                                    loop.data[idx - num_of_dim + dim_id_1][seq_id_col]
                                loop.data[idx - num_of_dim + dim_id_2][comp_id_col] =\
                                    loop.data[idx - num_of_dim + dim_id_1][comp_id_col]
                                loop.data[idx - num_of_dim + dim_id_2][auth_seq_id_col] =\
                                    loop.data[idx - num_of_dim + dim_id_1][auth_seq_id_col]
                                loop.data[idx - num_of_dim + dim_id_2][auth_comp_id_col] =\
                                    loop.data[idx - num_of_dim + dim_id_1][auth_comp_id_col]
                                loop.data[idx - num_of_dim + dim_id_2][atom_id_col] =\
                                    loop.data[idx - num_of_dim + dim_id_2][auth_atom_id_col] = common_atom_id2

                            diff_ = diff2_ = None
                            if shift_ is not None:
                                diff_ = ((position - shift_) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                            if shift2_ is not None:
                                diff2_ = ((position - shift) * weight) ** 2 + ((position2 - shift2_) * weight2) ** 2

                            if diff_ is not None and diff2_ is not None:
                                if diff_ < diff and diff2_ < diff:
                                    if diff_ < diff2_ and len_atom_id_ > 0:
                                        alt_swap_seq_atom_id_1(_atom_id_)
                                        continue
                                    if diff_ > diff2_ and len_atom_id2_ > 0:
                                        alt_swap_seq_atom_id_2(_atom_id2_)
                                        continue
                                elif diff_ < diff and len_atom_id_ > 0:
                                    alt_swap_seq_atom_id_1(_atom_id_)
                                    continue
                                elif diff2_ < diff and len_atom_id2_ > 0:
                                    alt_swap_seq_atom_id_2(_atom_id2_)
                                    continue
                                else:
                                    if diff_ < diff2_ and diff_ < 1.0 and len_atom_id_ > 0:
                                        alt_swap_seq_atom_id_1(_atom_id_)
                                        continue
                                    if diff_ > diff2_ and diff2_ < 1.0 and len_atom_id2_ > 0:
                                        alt_swap_seq_atom_id_2(_atom_id2_)
                                        continue

                            elif diff_ is not None and (diff_ < diff or diff_ < 1.0) and len_atom_id_ > 0:
                                alt_swap_seq_atom_id_1(_atom_id_)
                                continue

                            elif diff2_ is not None and (diff2_ < diff or diff2_ < 1.0) and len_atom_id2_ > 0:
                                alt_swap_seq_atom_id_2(_atom_id2_)
                                continue

                        if is_reparsable:
                            self.reasonsForReParsing['onebond_idx_history'] = self.onebond_idx_history
                            is_reparsable = False

                        else:
                            del_idx_list.extend(list(range(idx - num_of_dim + 1, idx + 1)))

                    elif chain_id != chain_id2:
                        position, position2 = positions[_dim_id_1], positions[_dim_id_2]

                        shift, weight = self.__getCsValue(chain_id, seq_id, comp_id, atom_id)
                        shift2, weight2 = self.__getCsValue(chain_id2, seq_id2, comp_id2, atom_id2)

                        shift_, _ = self.__getCsValue(chain_id2, seq_id2, comp_id2, atom_id)
                        shift2_, _ = self.__getCsValue(chain_id, seq_id, comp_id, atom_id2)

                        if shift_ is not None and not self.ccU.hasBond(comp_id2, atom_id, atom_id2):
                            shift_ = None

                        if shift2_ is not None and not self.ccU.hasBond(comp_id, atom_id, atom_id2):
                            shift2_ = None

                        # pylint: disable=cell-var-from-loop
                        def alt_swap_chain_seq_id_1():
                            if loop.data[idx][details_col] in EMPTY_VALUE:
                                loop.data[idx - num_of_dim + dim_id_1][details_col] =\
                                    f'{chain_id}:{seq_id}:{comp_id} -> {chain_id2}:{seq_id2}:{comp_id2}'
                            loop.data[idx - num_of_dim + dim_id_1][chain_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_2][chain_id_col]
                            loop.data[idx - num_of_dim + dim_id_1][entity_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_2][entity_id_col]
                            loop.data[idx - num_of_dim + dim_id_1][seq_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_2][seq_id_col]
                            loop.data[idx - num_of_dim + dim_id_1][comp_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_2][comp_id_col]
                            loop.data[idx - num_of_dim + dim_id_1][auth_chain_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_2][auth_chain_id_col]
                            loop.data[idx - num_of_dim + dim_id_1][auth_seq_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_2][auth_seq_id_col]
                            loop.data[idx - num_of_dim + dim_id_1][auth_comp_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_2][auth_comp_id_col]

                        # pylint: disable=cell-var-from-loop
                        def alt_swap_chain_seq_id_2():
                            if loop.data[idx - num_of_dim + dim_id_2][details_col] in EMPTY_VALUE:
                                loop.data[idx - num_of_dim + dim_id_2][details_col] =\
                                    f'{chain_id2}:{seq_id2}:{comp_id2} -> {chain_id}:{seq_id}:{comp_id}'
                            loop.data[idx - num_of_dim + dim_id_2][chain_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_1][chain_id_col]
                            loop.data[idx - num_of_dim + dim_id_2][entity_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_1][entity_id_col]
                            loop.data[idx - num_of_dim + dim_id_2][seq_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_1][seq_id_col]
                            loop.data[idx - num_of_dim + dim_id_2][comp_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_1][comp_id_col]
                            loop.data[idx - num_of_dim + dim_id_2][auth_chain_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_1][auth_chain_id_col]
                            loop.data[idx - num_of_dim + dim_id_2][auth_seq_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_1][auth_seq_id_col]
                            loop.data[idx - num_of_dim + dim_id_2][auth_comp_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_1][auth_comp_id_col]

                        if None in (shift, shift2):

                            if shift is None and shift2 is None:

                                if is_reparsable:
                                    self.reasonsForReParsing['onebond_idx_history'] = self.onebond_idx_history
                                    is_reparsable = False

                                else:
                                    del_idx_list.extend(list(range(idx - num_of_dim + 1, idx + 1)))

                                continue

                            if shift is None and shift_ is not None:
                                alt_swap_chain_seq_id_1()
                                continue

                            if shift2 is None and shift2_ is not None:
                                alt_swap_chain_seq_id_2()
                                continue

                            if shift is None and shift_ is None and None not in (shift2, shift2_):
                                alt_swap_chain_seq_id_1()
                                continue

                            if shift2 is None and shift2_ is None and None not in (shift, shift_):
                                alt_swap_chain_seq_id_2()
                                continue

                            if is_reparsable:
                                self.reasonsForReParsing['onebond_idx_history'] = self.onebond_idx_history
                                is_reparsable = False

                            else:
                                del_idx_list.extend(list(range(idx - num_of_dim + 1, idx + 1)))

                            continue

                        diff = ((position - shift) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                        diff *= 2.0

                        diff_ = diff2_ = None
                        if shift_ is not None:
                            diff_ = ((position - shift_) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                        if shift2_ is not None:
                            diff2_ = ((position - shift) * weight) ** 2 + ((position2 - shift2_) * weight2) ** 2

                        if diff_ is not None and diff2_ is not None:
                            if diff_ < diff and diff2_ < diff:
                                if diff_ < diff2_:
                                    alt_swap_chain_seq_id_1()
                                    continue
                                if diff_ > diff2_:
                                    alt_swap_chain_seq_id_2()
                                    continue
                            elif diff_ < diff:
                                alt_swap_chain_seq_id_1()
                                continue
                            elif diff2_ < diff:
                                alt_swap_chain_seq_id_2()
                                continue
                            else:
                                if diff_ < diff2_ and diff_ < 1.0:
                                    alt_swap_chain_seq_id_1()
                                    continue
                                if diff_ > diff2_ and diff2_ < 1.0:
                                    alt_swap_chain_seq_id_2()
                                    continue

                        elif diff_ is not None and (diff_ < diff or diff_ < 1.0):
                            alt_swap_chain_seq_id_1()
                            continue

                        elif diff2_ is not None and (diff2_ < diff or diff2_ < 1.0):
                            alt_swap_chain_seq_id_2()
                            continue

                        if is_reparsable:
                            self.reasonsForReParsing['onebond_idx_history'] = self.onebond_idx_history
                            is_reparsable = False

                        else:
                            del_idx_list.extend(list(range(idx - num_of_dim + 1, idx + 1)))

        if len(del_idx_list) > 0:
            for _idx in reversed(del_idx_list):
                del loop.data[_idx]

    def __remediatePeakAssignmentForJcouplingTransfer(self, num_of_dim: int, jcoupling_transfers: List[List[int]],
                                                      use_peak_row_format: bool, loop: pynmrstar.Loop) -> None:
        """ Remediate peak assignment based on J-coupling transfer.
        """

        self.__remediatePeakAssignmentForTransfer(num_of_dim, jcoupling_transfers, use_peak_row_format, loop, False)

    def __remediatePeakAssignmentForRelayedTransfer(self, num_of_dim: int, relayed_transfers: List[List[int]],
                                                    use_peak_row_format: bool, loop: pynmrstar.Loop) -> None:
        """ Remediate peak assignment based on relayed transfer.
        """

        self.__remediatePeakAssignmentForTransfer(num_of_dim, relayed_transfers, use_peak_row_format, loop, True)

    def __remediatePeakAssignmentForTransfer(self, num_of_dim: int, transfers: List[List[int]],
                                             use_peak_row_format: bool, loop: pynmrstar.Loop, relayed: bool) -> None:
        """ Remediate peak assignment based on relayed transfer (relayed), or J-coupling transfer otherwise.
        """

        idx_history_key = 'relayed_idx_history' if relayed else 'jcoupling_idx_history'

        is_reparsable = self.reasons is None and self.software_name != 'PIPP'

        details_col = loop.tags.index('Details')

        for dim_id_1, dim_id_2 in transfers:

            if use_peak_row_format:

                tags = [f'Entity_assembly_ID_{dim_id_1}', f'Comp_index_ID_{dim_id_1}',
                        f'Comp_ID_{dim_id_1}', f'Atom_ID_{dim_id_1}', f'Position_{dim_id_1}',
                        f'Entity_assembly_ID_{dim_id_2}', f'Comp_index_ID_{dim_id_2}',
                        f'Comp_ID_{dim_id_2}', f'Atom_ID_{dim_id_2}', f'Position_{dim_id_2}']

                dat = loop.get_tag(tags)

                for idx, row in enumerate(dat):

                    if any(True for col in range(10) if row[col] in EMPTY_VALUE):
                        continue

                    chain_id, seq_id, comp_id, atom_id, chain_id2, seq_id2, comp_id2, atom_id2 =\
                        row[0], row[1], row[2], row[3], row[5], row[6], row[7], row[8]

                    if isinstance(chain_id, int):
                        chain_id = str(chain_id)

                    if isinstance(chain_id2, int):
                        chain_id2 = str(chain_id2)

                    if isinstance(seq_id, str):
                        seq_id = int(seq_id)

                    if isinstance(seq_id2, str):
                        seq_id2 = int(seq_id2)

                    if chain_id == chain_id2 and (abs(seq_id - seq_id2) < 2 if relayed else seq_id == seq_id2):
                        continue

                    if chain_id == chain_id2 and (relayed or seq_id != seq_id2):
                        position, position2 = row[4], row[9]

                        if isinstance(chain_id, int):
                            chain_id = str(chain_id)

                        if isinstance(position, str):
                            position = float(position)

                        if isinstance(position2, str):
                            position2 = float(position2)

                        shift, weight = self.__getCsValue(chain_id, seq_id, comp_id, atom_id)
                        shift2, weight2 = self.__getCsValue(chain_id, seq_id2, comp_id2, atom_id2)

                        # pylint: disable=cell-var-from-loop
                        def swap_seq_id_1():
                            if loop.data[idx][details_col] in EMPTY_VALUE:
                                loop.data[idx][details_col] = f'{seq_id}:{comp_id} -> {seq_id2}:{comp_id2}'
                            loop.data[idx][loop.tags.index(f'Comp_index_ID_{dim_id_1}')] =\
                                loop.data[idx][loop.tags.index(f'Comp_index_ID_{dim_id_2}')]
                            loop.data[idx][loop.tags.index(f'Comp_ID_{dim_id_1}')] =\
                                loop.data[idx][loop.tags.index(f'Comp_ID_{dim_id_2}')]
                            loop.data[idx][loop.tags.index(f'Auth_seq_ID_{dim_id_1}')] =\
                                loop.data[idx][loop.tags.index(f'Auth_seq_ID_{dim_id_2}')]
                            loop.data[idx][loop.tags.index(f'Auth_comp_ID_{dim_id_1}')] =\
                                loop.data[idx][loop.tags.index(f'Auth_comp_ID_{dim_id_2}')]

                        # pylint: disable=cell-var-from-loop
                        def swap_seq_id_2():
                            if loop.data[idx][details_col] in EMPTY_VALUE:
                                loop.data[idx][details_col] = f'{seq_id2}:{comp_id2} -> {seq_id}:{comp_id}'
                            loop.data[idx][loop.tags.index(f'Comp_index_ID_{dim_id_2}')] =\
                                loop.data[idx][loop.tags.index(f'Comp_index_ID_{dim_id_1}')]
                            loop.data[idx][loop.tags.index(f'Comp_ID_{dim_id_2}')] =\
                                loop.data[idx][loop.tags.index(f'Comp_ID_{dim_id_1}')]
                            loop.data[idx][loop.tags.index(f'Auth_seq_ID_{dim_id_2}')] =\
                                loop.data[idx][loop.tags.index(f'Auth_seq_ID_{dim_id_1}')]
                            loop.data[idx][loop.tags.index(f'Auth_comp_ID_{dim_id_2}')] =\
                                loop.data[idx][loop.tags.index(f'Auth_comp_ID_{dim_id_1}')]

                        if None in (shift, shift2):

                            if is_reparsable:
                                self.reasonsForReParsing[idx_history_key] = getattr(self, idx_history_key)
                                is_reparsable = False

                            continue

                        shift_, _ = self.__getCsValue(chain_id, seq_id2, comp_id2, atom_id)
                        shift2_, _ = self.__getCsValue(chain_id, seq_id, comp_id, atom_id2)

                        diff = ((position - shift) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                        diff *= 2.0

                        diff_ = diff2_ = None
                        if shift_ is not None:
                            diff_ = ((position - shift_) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                        if shift2_ is not None:
                            diff2_ = ((position - shift) * weight) ** 2 + ((position2 - shift2_) * weight2) ** 2

                        if diff_ is not None and diff2_ is not None:
                            if diff_ < diff and diff2_ < diff:
                                if diff_ < diff2_:
                                    swap_seq_id_1()
                                    continue
                                if diff_ > diff2_:
                                    swap_seq_id_2()
                                    continue
                            elif diff_ < diff:
                                swap_seq_id_1()
                                continue
                            elif diff2_ < diff:
                                swap_seq_id_2()
                                continue
                            else:
                                if diff_ < diff2_ and diff_ < 1.0:
                                    swap_seq_id_1()
                                    continue
                                if diff_ > diff2_ and diff2_ < 1.0:
                                    swap_seq_id_2()
                                    continue

                        elif diff_ is not None and (diff_ < diff or diff_ < 1.0):
                            swap_seq_id_1()
                            continue

                        elif diff2_ is not None and (diff2_ < diff or diff2_ < 1.0):
                            swap_seq_id_2()
                            continue

                        if is_reparsable:
                            self.reasonsForReParsing[idx_history_key] = getattr(self, idx_history_key)
                            is_reparsable = False

                    else:

                        if is_reparsable:
                            self.reasonsForReParsing[idx_history_key] = getattr(self, idx_history_key)
                            is_reparsable = False

            else:

                tags = ['Peak_ID', 'Spectral_dim_ID', 'Entity_assembly_ID', 'Comp_index_ID', 'Comp_ID', 'Atom_ID', 'Val']

                seq_id_col = loop.tags.index('Comp_index_ID')
                comp_id_col = loop.tags.index('Comp_ID')
                auth_seq_id_col = loop.tags.index('Auth_seq_ID')
                auth_comp_id_col = loop.tags.index('Auth_comp_ID')

                dat = loop.get_tag(tags)

                peak_id = None

                for idx, row in enumerate(dat):
                    dim_id = row[1]

                    if any(True for col in range(7) if row[col] in EMPTY_VALUE):
                        continue

                    if dim_id == 1:
                        peak_id = row[0]
                        chain_ids, seq_ids, comp_ids, atom_ids, positions = [], [], [], [], []
                    else:
                        if peak_id != row[0]:
                            continue

                    chain_ids.append(row[2] if isinstance(row[2], str) else str(row[2]))  # noqa: E501, pylint: disable=possibly-used-before-assignment,line-too-long
                    seq_ids.append(int(row[3]) if isinstance(row[3], str) else row[3])  # noqa: E501, pylint: disable=possibly-used-before-assignment,line-too-long
                    comp_ids.append(row[4])  # pylint: disable=possibly-used-before-assignment
                    atom_ids.append(row[5])  # pylint: disable=possibly-used-before-assignment
                    positions.append(float(row[6]) if isinstance(row[6], str) else row[6])  # noqa: E501, pylint: disable=possibly-used-before-assignment,line-too-long

                    if dim_id < num_of_dim:
                        continue

                    if len(atom_ids) < num_of_dim:
                        continue

                    _dim_id_1 = dim_id_1 - 1
                    _dim_id_2 = dim_id_2 - 1

                    chain_id, seq_id, comp_id, atom_id, position, chain_id2, seq_id2, comp_id2, atom_id2, position2 =\
                        chain_ids[_dim_id_1], seq_ids[_dim_id_1], comp_ids[_dim_id_1], atom_ids[_dim_id_1], positions[_dim_id_1], \
                        chain_ids[_dim_id_2], seq_ids[_dim_id_2], comp_ids[_dim_id_2], atom_ids[_dim_id_2], positions[_dim_id_2]

                    if chain_id == chain_id2 and (abs(seq_id - seq_id2) < 2 if relayed else seq_id == seq_id2):
                        continue

                    if chain_id == chain_id2 and (relayed or seq_id != seq_id2):
                        chain_id = chain_ids[_dim_id_1]

                        if isinstance(chain_id, int):
                            chain_id = str(chain_id)

                        shift, weight = self.__getCsValue(chain_id, seq_id, comp_id, atom_id)
                        shift2, weight2 = self.__getCsValue(chain_id, seq_id2, comp_id2, atom_id2)

                        # pylint: disable=cell-var-from-loop
                        def alt_swap_seq_id_1():
                            if loop.data[idx][details_col] in EMPTY_VALUE:
                                loop.data[idx - num_of_dim + dim_id_1][details_col] = f'{seq_id}:{comp_id} -> {seq_id2}:{comp_id2}'
                            loop.data[idx - num_of_dim + dim_id_1][seq_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_2][seq_id_col]
                            loop.data[idx - num_of_dim + dim_id_1][comp_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_2][comp_id_col]
                            loop.data[idx - num_of_dim + dim_id_1][auth_seq_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_2][auth_seq_id_col]
                            loop.data[idx - num_of_dim + dim_id_1][auth_comp_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_2][auth_comp_id_col]

                        # pylint: disable=cell-var-from-loop
                        def alt_swap_seq_id_2():
                            if loop.data[idx - num_of_dim + dim_id_2][details_col] in EMPTY_VALUE:
                                loop.data[idx - num_of_dim + dim_id_2][details_col] = f'{seq_id2}:{comp_id2} -> {seq_id}:{comp_id}'
                            loop.data[idx - num_of_dim + dim_id_2][seq_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_1][seq_id_col]
                            loop.data[idx - num_of_dim + dim_id_2][comp_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_1][comp_id_col]
                            loop.data[idx - num_of_dim + dim_id_2][auth_seq_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_1][auth_seq_id_col]
                            loop.data[idx - num_of_dim + dim_id_2][auth_comp_id_col] =\
                                loop.data[idx - num_of_dim + dim_id_1][auth_comp_id_col]

                        if None in (shift, shift2):

                            if is_reparsable:
                                self.reasonsForReParsing[idx_history_key] = getattr(self, idx_history_key)
                                is_reparsable = False

                            continue

                        shift_, _ = self.__getCsValue(chain_id, seq_id2, comp_id2, atom_id)
                        shift2_, _ = self.__getCsValue(chain_id, seq_id, comp_id, atom_id2)

                        diff = ((position - shift) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                        diff *= 2.0

                        diff_ = diff2_ = None
                        if shift_ is not None:
                            diff_ = ((position - shift_) * weight) ** 2 + ((position2 - shift2) * weight2) ** 2
                        if shift2_ is not None:
                            diff2_ = ((position - shift) * weight) ** 2 + ((position2 - shift2_) * weight2) ** 2

                        if diff_ is not None and diff2_ is not None:
                            if diff_ < diff and diff2_ < diff:
                                if diff_ < diff2_:
                                    alt_swap_seq_id_1()
                                    continue
                                if diff_ > diff2_:
                                    alt_swap_seq_id_2()
                                    continue
                            elif diff_ < diff:
                                alt_swap_seq_id_1()
                                continue
                            elif diff2_ < diff:
                                alt_swap_seq_id_2()
                                continue
                            else:
                                if diff_ < diff2_ and diff_ < 1.0:
                                    alt_swap_seq_id_1()
                                    continue
                                if diff_ > diff2_ and diff2_ < 1.0:
                                    alt_swap_seq_id_2()
                                    continue

                        elif diff_ is not None and (diff_ < diff or diff_ < 1.0):
                            alt_swap_seq_id_1()
                            continue

                        elif diff2_ is not None and (diff2_ < diff or diff2_ < 1.0):
                            alt_swap_seq_id_2()
                            continue

                        if is_reparsable:
                            self.reasonsForReParsing[idx_history_key] = getattr(self, idx_history_key)
                            is_reparsable = False

                    else:

                        if is_reparsable:
                            self.reasonsForReParsing[idx_history_key] = getattr(self, idx_history_key)
                            is_reparsable = False

    def __getCsValue(self, chain_id: str, seq_id: int, comp_id: str, atom_id: str
                     ) -> Tuple[Optional[float], Optional[float]]:
        """ Retrieve CS value for a given atom as reference for peak assignment.
        """

        if self.__csLoops is None or len(self.__csLoops) == 0:
            return None, None

        _atom_ids = self.nefT.get_valid_star_atom(comp_id, atom_id, leave_unmatched=False)[0]

        for lp in self.__csLoops:
            row = next((row for row in lp['data']
                        if row['Entity_assembly_ID'] == chain_id and row['Comp_index_ID'] == seq_id
                        and row['Comp_ID'] == comp_id and row['Atom_ID'] in _atom_ids), None)

            if row is not None:
                isotope_num = row['Atom_isotope_number']
                weight = 1.0
                if isotope_num == 13:
                    weight = 0.251449530
                elif isotope_num == 15:
                    weight = 0.101329118
                elif isotope_num == 2:
                    weight = 0.153506088
                return row['Val'], weight

            try:

                row = next(row for row in self.__tempCsValues
                           if row['chain_id'] == chain_id
                           and row['seq_id'] == seq_id
                           and row['comp_id'] == comp_id
                           and row['atom_id'] == atom_id)

                return row['value'], row['weight']

            except StopIteration:
                pass

        return None, None

    def __setTempCsValue(self, star_atom: dict, values: List[float]) -> None:
        """ Set reference CS value for a given atom as reference for peak assignment.
        """

        try:

            next(row for row in self.__tempCsValues
                 if row['chain_id'] == star_atom['chain_id']
                 and row['seq_id'] == star_atom['seq_id']
                 and row['comp_id'] == star_atom['comp_id']
                 and row['atom_id'] == star_atom['atom_id'])

        except StopIteration:
            atom_type = star_atom['atom_id'][0]

            weight = 1.0
            if atom_type == 'C':
                weight = 0.251449530
            elif atom_type == 'N':
                weight = 0.101329118

            own_values = [row['value'] for row in self.__tempCsValues
                          if row['chain_id'] == star_atom['chain_id']
                          and row['seq_id'] == star_atom['seq_id']
                          and row['comp_id'] == star_atom['comp_id']]

            for value in values:

                unique = True
                for _value in own_values:
                    if ((value - _value) * weight) ** 2 < 0.02:
                        unique = False
                        break

                if unique:
                    self.__tempCsValues.append({'chain_id': star_atom['chain_id'],
                                                'seq_id': star_atom['seq_id'],
                                                'comp_id': star_atom['comp_id'],
                                                'atom_id': star_atom['atom_id'],
                                                'weight': weight,
                                                'value': value})
                    return

    def __validatePeak(self, index: int, positions: tuple, pos_uncs: tuple, lws: tuple, lw_hzs: tuple,
                       height: Optional[str], height_uncertainty: Optional[str],
                       volume: Optional[str], volume_uncertainty: Optional[str],
                       figure_of_merit: Optional[Union[float, int]], num_of_dim: int
                       ) -> Optional[dict]:
        """ Validate value range of peak of a given dimension (positions, pos_uncs, lws, lw_hzs are given per dimension).
        """

        validRange = True
        dstFunc = {}

        for dim, pos in enumerate(positions, start=1):
            if CS_ERROR_MIN < pos < CS_ERROR_MAX:
                dstFunc[f'position_{dim}'] = str(pos)
            else:
                validRange = False
                self.f.append(f"[Range value error] {self.getCurrentSpectralPeak(n=index)}"
                              f"The position_{dim}='{pos}' must be within range {CS_RESTRAINT_ERROR}.")

        if not validRange:
            return None

        for dim, pos in enumerate(positions, start=1):
            if CS_RANGE_MIN <= pos <= CS_RANGE_MAX:
                pass
            else:
                self.f.append(f"[Range value warning] {self.getCurrentSpectralPeak(n=index)}"
                              f"The position_{dim}='{pos}' should be within range {CS_RESTRAINT_RANGE}.")

        try:
            if height is not None and float(height) != 0.0:
                dstFunc['height'] = height
        except ValueError:
            pass
        try:
            if volume is not None and float(volume) != 0.0:
                dstFunc['volume'] = volume
        except ValueError:
            pass
        try:
            if height_uncertainty is not None and float(height_uncertainty) != 0.0:
                dstFunc['height_uncertainty'] = height_uncertainty
        except ValueError:
            pass
        try:
            if volume_uncertainty is not None and float(volume_uncertainty) != 0.0:
                dstFunc['volume_uncertainty'] = volume_uncertainty
        except ValueError:
            pass

        if 'height' not in dstFunc and 'volume' not in dstFunc and not self.__internal:
            self.f.append(f"[Missing data] {self.getCurrentSpectralPeak(n=index)}"
                          "Neither peak height nor peak volume value is set. Please re-upload the NMR spectral peak list file.")
            return None

        for dim, pos_unc in enumerate(pos_uncs, start=1):
            if pos_unc is not None and pos_unc != 0.0:
                dstFunc[f'position_uncertainty_{dim}'] = str(pos_unc) if pos_unc > 0.0 else str(abs(pos_unc))

        for dim, (lw, lw_hz) in enumerate(zip(lws, lw_hzs), start=1):
            if lw_hz is not None and lw_hz != 0.0:
                dstFunc[f'line_width_{dim}'] = str(lw_hz) if lw_hz > 0.0 else str(abs(lw_hz))
            elif lw is not None and lw != 0.0:
                dstFunc[f'line_width_{dim}'] = str(lw) if lw > 0.0 else str(abs(lw))

        if figure_of_merit is not None:
            if WEIGHT_RANGE_MIN <= figure_of_merit <= WEIGHT_RANGE_MAX:
                dstFunc['figure_of_merit'] = str(figure_of_merit)
            else:
                self.f.append(f"[Range value warning] {self.getCurrentSpectralPeak(n=index)}"
                              f"The figure_of_merit='{figure_of_merit}' should be within range {WEIGHT_RANGE}.")

        if getattr(self, f'peaks{num_of_dim}D') == 1 and self.__defaultSegId__ is not None:
            self.__defaultSegId = self.__defaultSegId__
            if self.reasons is not None and 'default_seg_id' in self.reasons:
                try:
                    self.__defaultSegId = self.reasons['default_seg_id'][num_of_dim][self.cur_list_id]
                except KeyError:
                    pass

        return dstFunc

    def validatePeak2D(self, index: int, pos_1: float, pos_2: float,
                       pos_unc_1: Optional[float], pos_unc_2: Optional[float],
                       lw_1: Optional[float], lw_2: Optional[float],
                       pos_hz_1: Optional[float], pos_hz_2: Optional[float],  # pylint: disable=unused-argument
                       lw_hz_1: Optional[float], lw_hz_2: Optional[float],
                       height: Optional[str], height_uncertainty: Optional[str],
                       volume: Optional[str], volume_uncertainty: Optional[str],
                       figure_of_merit: Optional[Union[float, int]] = None
                       ) -> Optional[dict]:
        """ Validate value range of 2D peak.
        """

        return self.__validatePeak(index, (pos_1, pos_2), (pos_unc_1, pos_unc_2),
                                   (lw_1, lw_2), (lw_hz_1, lw_hz_2),
                                   height, height_uncertainty, volume, volume_uncertainty, figure_of_merit, 2)

    def validatePeak3D(self, index: int, pos_1: float, pos_2: float, pos_3: float,
                       pos_unc_1: Optional[float], pos_unc_2: Optional[float], pos_unc_3: Optional[float],
                       lw_1: Optional[float], lw_2: Optional[float], lw_3: Optional[float],
                       pos_hz_1: Optional[float], pos_hz_2: Optional[float],  # pylint: disable=unused-argument
                       pos_hz_3: Optional[float],  # pylint: disable=unused-argument
                       lw_hz_1: Optional[float], lw_hz_2: Optional[float], lw_hz_3: Optional[float],
                       height: Optional[str], height_uncertainty: Optional[str],
                       volume: Optional[str], volume_uncertainty: Optional[str],
                       figure_of_merit: Optional[Union[float, int]] = None
                       ) -> Optional[dict]:
        """ Validate value range of 3D peak.
        """

        return self.__validatePeak(index, (pos_1, pos_2, pos_3), (pos_unc_1, pos_unc_2, pos_unc_3),
                                   (lw_1, lw_2, lw_3), (lw_hz_1, lw_hz_2, lw_hz_3),
                                   height, height_uncertainty, volume, volume_uncertainty, figure_of_merit, 3)

    def validatePeak4D(self, index: int, pos_1: float, pos_2: float, pos_3: float, pos_4: float,
                       pos_unc_1: Optional[float], pos_unc_2: Optional[float],
                       pos_unc_3: Optional[float], pos_unc_4: Optional[float],
                       lw_1: Optional[float], lw_2: Optional[float], lw_3: Optional[float], lw_4: Optional[float],
                       pos_hz_1: Optional[float], pos_hz_2: Optional[float],  # pylint: disable=unused-argument
                       pos_hz_3: Optional[float], pos_hz_4: Optional[float],  # pylint: disable=unused-argument
                       lw_hz_1: Optional[float], lw_hz_2: Optional[float], lw_hz_3: Optional[float], lw_hz_4: Optional[float],
                       height: Optional[str], height_uncertainty: Optional[str],
                       volume: Optional[str], volume_uncertainty: Optional[str],
                       figure_of_merit: Optional[Union[float, int]] = None
                       ) -> Optional[dict]:
        """ Validate value range of 4D peak.
        """

        return self.__validatePeak(index, (pos_1, pos_2, pos_3, pos_4), (pos_unc_1, pos_unc_2, pos_unc_3, pos_unc_4),
                                   (lw_1, lw_2, lw_3, lw_4), (lw_hz_1, lw_hz_2, lw_hz_3, lw_hz_4),
                                   height, height_uncertainty, volume, volume_uncertainty, figure_of_merit, 4)

    def selectProbablePosition(self, index: int, label: str, positions: List[float],
                               with_segid: Optional[str] = None, with_compid: Optional[str] = None) -> float:
        """ Select the most realistic peak position with given conditions.
        """

        position = positions[0]

        if label is None:
            return position

        ext = self.extractPeakAssignment(1, label, index, with_segid, with_compid)

        if ext is None:
            return position

        assignment = None

        status, _label = self.testAssignment(1, ext, label)

        if status:
            assignment = ext
        elif _label is not None:
            ext = self.extractPeakAssignment(1, _label, index, with_segid, with_compid)
            if ext is not None:
                assignment = ext

        if assignment is None:
            return position

        self.checkAssignment(index, assignment)

        if len(self.atomSelectionSet) == 0:
            return position

        _diff = None

        for atom in self.atomSelectionSet[0]:
            star_atom = getStarAtom(self.__authToStarSeq, self.__authToOrigSeq, self.__offsetHolder, atom)

            if star_atom is None:
                continue

            star_atom['chain_id'] = str(star_atom['chain_id'])

            shift, weight =\
                self.__getCsValue(star_atom['chain_id'], star_atom['seq_id'], star_atom['comp_id'], star_atom['atom_id'])

            if shift is None:
                self.__setTempCsValue(star_atom, positions)

                shift, weight =\
                    self.__getCsValue(star_atom['chain_id'], star_atom['seq_id'], star_atom['comp_id'], star_atom['atom_id'])

                if shift is None:
                    continue

            _position = None

            for _position_ in positions:
                diff = ((_position_ - shift) * weight) ** 2
                if _diff is None or diff < _diff:
                    _position = _position_
                    _diff = diff

            if _position is not None and _diff < 1.0:
                position = _position

        self.atomSelectionSet.clear()

        return position

    def checkAssignment(self, index: int, assignment: List[dict]) -> None:
        """ Evaluate peak assignment for selecting the most realistic peak position.
        """

        if assignment is not None:

            self.retrieveLocalSeqScheme()

            try:

                hasChainId = all(a['chain_id'] is not None for a in assignment)
                hasCompId = all(a['comp_id'] is not None for a in assignment)

                for a1 in assignment:

                    self.atomSelectionSet.clear()

                    if hasChainId and hasCompId:
                        chainAssign1, _ =\
                            self.assignCoordPolymerSequenceWithChainId(a1['chain_id'], a1['seq_id'], a1['comp_id'], a1['atom_id'],
                                                                       index)

                    elif hasChainId:
                        chainAssign1 =\
                            self.assignCoordPolymerSequenceWithChainIdWithoutCompId(a1['chain_id'], a1['seq_id'], a1['atom_id'],
                                                                                    index)

                    elif hasCompId:
                        chainAssign1, _ =\
                            self.assignCoordPolymerSequence(a1['chain_id'], a1['seq_id'], a1['comp_id'], a1['atom_id'], index)

                    else:
                        chainAssign1 =\
                            self.assignCoordPolymerSequenceWithoutCompId(a1['seq_id'], a1['atom_id'], index)

                    if len(chainAssign1) > 0:
                        self.selectCoordAtoms(chainAssign1, a1.get('auth_chain_id'),
                                              a1['seq_id'], a1['comp_id'], a1['atom_id'], index)

            except (KeyError, TypeError):
                pass

    def checkAssignments2D(self, index: int, assignments: List[List[dict]], dstFunc: dict
                           ) -> Tuple[bool, bool, Optional[bool], Optional[bool]]:
        """ Evaluate 2D peak assignment.
        """

        has_assignments = has_multiple_assignments = False
        asis1 = asis2 = None

        primary_specral_dim = self.spectral_dim[self.num_of_dim][self.cur_list_id][1]
        has_long_range = 'long_range' in primary_specral_dim

        if 'is_noesy' not in primary_specral_dim:
            file_name = self.__originalFileName.lower()
            alt_file_name = '' if self.spectrum_name is None else self.spectrum_name.lower()
            _file_names = (file_name, alt_file_name)
            primary_specral_dim['is_noesy'] = any('noe' in n for n in _file_names) or any('roe' in n for n in _file_names)
            if primary_specral_dim['is_noesy']:
                primary_specral_dim['long_range'] = True

        if all(assignment is not None for assignment in assignments):

            if self.reasons is not None and 'onebond_resolved' in self.reasons:
                _assignments = [None] * 2
                for k, v in self.reasons['onebond_resolved'].items():
                    _assignments[v] = assignments[k]
                assignments = _assignments

            self.retrieveLocalSeqScheme()

            try:

                hasChainId = all(all(_a['chain_id'] is not None for _a in a) for a in assignments)
                hasCompId = all(all(_a['comp_id'] is not None for _a in a) for a in assignments)

                has_multiple_assignments = any(True for assignment in assignments if len(assignment) > 1)

                pairs = []
                if len(assignments[0]) == len(assignments[1]):
                    for a1, a2 in zip(assignments[0], assignments[1]):
                        pairs.append((a1, a2))
                else:
                    for a1, a2 in itertools.product(assignments[0], assignments[1]):
                        pairs.append((a1, a2))

                for a1, a2 in pairs:

                    self.atomSelectionSet.clear()
                    asis1 = asis2 = None

                    if hasChainId and hasCompId:
                        chainAssign1, asis1 =\
                            self.assignCoordPolymerSequenceWithChainId(a1['chain_id'], a1['seq_id'], a1['comp_id'], a1['atom_id'],
                                                                       index)
                        chainAssign2, asis2 =\
                            self.assignCoordPolymerSequenceWithChainId(a2['chain_id'], a2['seq_id'], a2['comp_id'], a2['atom_id'],
                                                                       index)

                    elif hasChainId:
                        chainAssign1 =\
                            self.assignCoordPolymerSequenceWithChainIdWithoutCompId(a1['chain_id'], a1['seq_id'], a1['atom_id'],
                                                                                    index)
                        chainAssign2 =\
                            self.assignCoordPolymerSequenceWithChainIdWithoutCompId(a2['chain_id'], a2['seq_id'], a2['atom_id'],
                                                                                    index)

                    elif hasCompId:
                        chainAssign1, asis1 =\
                            self.assignCoordPolymerSequence(a1['chain_id'], a1['seq_id'], a1['comp_id'], a1['atom_id'], index)
                        chainAssign2, asis2 =\
                            self.assignCoordPolymerSequence(a2['chain_id'], a2['seq_id'], a2['comp_id'], a2['atom_id'], index)

                    else:
                        chainAssign1 = self.assignCoordPolymerSequenceWithoutCompId(a1['seq_id'], a1['atom_id'], index)
                        chainAssign2 = self.assignCoordPolymerSequenceWithoutCompId(a2['seq_id'], a2['atom_id'], index)

                    if len(chainAssign1) > 0 and len(chainAssign2) > 0:
                        self.selectCoordAtoms(chainAssign1, a1.get('auth_chain_id'),
                                              a1['seq_id'], a1['comp_id'], a1['atom_id'], index)
                        self.selectCoordAtoms(chainAssign2, a2.get('auth_chain_id'),
                                              a2['seq_id'], a2['comp_id'], a2['atom_id'], index)

                        if len(self.atomSelectionSet) == self.num_of_dim:
                            has_assignments = True
                            has_assignments &=\
                                self.validateAtomType(1, self.atomSelectionSet[0][0]['atom_id'][0], dstFunc['position_1'])
                            has_assignments &=\
                                self.validateAtomType(2, self.atomSelectionSet[1][0]['atom_id'][0], dstFunc['position_2'])
                            if has_assignments:
                                self.atomSelectionSets.append([[dict(a) if a.__class__ is dict else a for a in sel]
                                                               for sel in self.atomSelectionSet])
                                self.asIsSets.append([asis1, asis2])
                                if self.reasons is not None:
                                    if 'onebond_idx_history' in self.reasons:
                                        onebond_idx = self.reasons['onebond_idx_history'][self.num_of_dim][self.cur_list_id]
                                        _atom1, _atom2 =\
                                            self.atomSelectionSet[DIM_TRANSFER_PAT_2D[onebond_idx][0]][0], \
                                            self.atomSelectionSet[DIM_TRANSFER_PAT_2D[onebond_idx][1]][0]
                                        if _atom1['chain_id'] != _atom2['chain_id']\
                                           or _atom1['seq_id'] != _atom2['seq_id']\
                                           or _atom1['atom_id'][0] == _atom2['atom_id'][0]\
                                           or not self.ccU.hasBond(_atom1['comp_id'], _atom1['atom_id'], _atom2['atom_id']):
                                            if not self.__canRemediatePeakAssignmentForOneBondTransfer(
                                                    _atom1, _atom2,
                                                    float(dstFunc[f'position_{DIM_TRANSFER_PAT_2D[onebond_idx][0] + 1}']),
                                                    float(dstFunc[f'position_{DIM_TRANSFER_PAT_2D[onebond_idx][1] + 1}'])):
                                                has_assignments = False
                                    if 'jcoupling_idx_history' in self.reasons:
                                        jcoupling_idx = self.reasons['jcoupling_idx_history'][self.num_of_dim][self.cur_list_id]
                                        if jcoupling_idx != -1:
                                            _atom1, _atom2 =\
                                                self.atomSelectionSet[DIM_TRANSFER_PAT_2D[jcoupling_idx][0]][0], \
                                                self.atomSelectionSet[DIM_TRANSFER_PAT_2D[jcoupling_idx][1]][0]
                                            if _atom1['chain_id'] != _atom2['chain_id']\
                                               or _atom1['seq_id'] != _atom2['seq_id']:
                                                if not self.__canRemediatePeakAssignmentForJcouplingTransfer(
                                                        _atom1, _atom2,
                                                        float(dstFunc[f'position_{DIM_TRANSFER_PAT_2D[jcoupling_idx][0] + 1}']),
                                                        float(dstFunc[f'position_{DIM_TRANSFER_PAT_2D[jcoupling_idx][1] + 1}'])):
                                                    has_assignments = False
                                    if 'relayed_idx_history' in self.reasons:
                                        relayed_idx = self.reasons['relayed_idx_history'][self.num_of_dim][self.cur_list_id]
                                        if relayed_idx != -1:
                                            _atom1, _atom2 =\
                                                self.atomSelectionSet[DIM_TRANSFER_PAT_2D[relayed_idx][0]][0], \
                                                self.atomSelectionSet[DIM_TRANSFER_PAT_2D[relayed_idx][1]][0]
                                            if _atom1['chain_id'] != _atom2['chain_id']\
                                               or _atom1['seq_id'] != _atom2['seq_id']:
                                                if not self.__canRemediatePeakAssignmentForRelayedTransfer(
                                                        _atom1, _atom2,
                                                        float(dstFunc[f'position_{DIM_TRANSFER_PAT_2D[relayed_idx][0] + 1}']),
                                                        float(dstFunc[f'position_{DIM_TRANSFER_PAT_2D[relayed_idx][1] + 1}'])):
                                                    has_assignments = False
                                if not has_long_range and has_assignments:
                                    if hasInterChainRestraint(self.atomSelectionSet):
                                        has_long_range = True
                                        if 'long_range' not in primary_specral_dim:
                                            primary_specral_dim['long_range'] = True
                                    else:
                                        for atom1, atom2 in itertools.product(self.atomSelectionSet[0],
                                                                              self.atomSelectionSet[1]):
                                            has_long_range |= isLongRangeRestraint([atom1, atom2], self.polySeq)
                                            if has_long_range:
                                                if 'long_range' not in primary_specral_dim:
                                                    primary_specral_dim['long_range'] = True
                                                break
                            else:
                                break
                        else:
                            has_assignments = False
                            break

            except (KeyError, TypeError):
                pass

        return has_assignments, has_multiple_assignments, asis1, asis2

    def checkAssignments3D(self, index: int, assignments: List[List[dict]], dstFunc: dict, onebondOrder: int = 0
                           ) -> Tuple[bool, bool, Optional[bool], Optional[bool], Optional[bool]]:
        """ Evaluate 3D peak assignment.
        """

        has_assignments = has_multiple_assignments = False
        asis1 = asis2 = asis3 = None

        primary_specral_dim = self.spectral_dim[self.num_of_dim][self.cur_list_id][1]
        has_long_range = 'long_range' in primary_specral_dim

        if 'is_noesy' not in primary_specral_dim:
            file_name = self.__originalFileName.lower()
            alt_file_name = '' if self.spectrum_name is None else self.spectrum_name.lower()
            _file_names = (file_name, alt_file_name)
            primary_specral_dim['is_noesy'] = any('noe' in n for n in _file_names) or any('roe' in n for n in _file_names)
            if primary_specral_dim['is_noesy']:
                primary_specral_dim['long_range'] = True

        if all(assignment is not None for assignment in assignments):

            if self.reasons is not None and 'onebond_resolved' in self.reasons:
                _assignments = [None] * 3
                for k, v in self.reasons['onebond_resolved'][onebondOrder].items():
                    _assignments[v] = assignments[k]
                assignments = _assignments

            self.retrieveLocalSeqScheme()

            try:

                hasChainId = all(all(_a['chain_id'] is not None for _a in a) for a in assignments)
                hasCompId = all(all(_a['comp_id'] is not None for _a in a) for a in assignments)

                has_multiple_assignments = any(True for assignment in assignments if len(assignment) > 1)

                pairs = []
                if len(assignments[0]) == len(assignments[1]) == len(assignments[2]):
                    for a1, a2, a3 in zip(assignments[0], assignments[1], assignments[2]):
                        pairs.append((a1, a2, a3))
                else:
                    for a1, a2, a3 in itertools.product(assignments[0], assignments[1], assignments[2]):
                        pairs.append((a1, a2, a3))

                for a1, a2, a3 in pairs:

                    self.atomSelectionSet.clear()
                    asis1 = asis2 = asis3 = None

                    if hasChainId and hasCompId:
                        chainAssign1, asis1 =\
                            self.assignCoordPolymerSequenceWithChainId(a1['chain_id'], a1['seq_id'], a1['comp_id'], a1['atom_id'],
                                                                       index)
                        chainAssign2, asis2 =\
                            self.assignCoordPolymerSequenceWithChainId(a2['chain_id'], a2['seq_id'], a2['comp_id'], a2['atom_id'],
                                                                       index)
                        chainAssign3, asis3 =\
                            self.assignCoordPolymerSequenceWithChainId(a3['chain_id'], a3['seq_id'], a3['comp_id'], a3['atom_id'],
                                                                       index)

                    elif hasChainId:
                        chainAssign1 =\
                            self.assignCoordPolymerSequenceWithChainIdWithoutCompId(a1['chain_id'], a1['seq_id'], a1['atom_id'],
                                                                                    index)
                        chainAssign2 =\
                            self.assignCoordPolymerSequenceWithChainIdWithoutCompId(a2['chain_id'], a2['seq_id'], a2['atom_id'],
                                                                                    index)
                        chainAssign3 =\
                            self.assignCoordPolymerSequenceWithChainIdWithoutCompId(a3['chain_id'], a3['seq_id'], a3['atom_id'],
                                                                                    index)

                    elif hasCompId:
                        chainAssign1, asis1 =\
                            self.assignCoordPolymerSequence(a1['chain_id'], a1['seq_id'], a1['comp_id'], a1['atom_id'], index)
                        chainAssign2, asis2 =\
                            self.assignCoordPolymerSequence(a2['chain_id'], a2['seq_id'], a2['comp_id'], a2['atom_id'], index)
                        chainAssign3, asis3 =\
                            self.assignCoordPolymerSequence(a3['chain_id'], a3['seq_id'], a3['comp_id'], a3['atom_id'], index)

                    else:
                        chainAssign1 = self.assignCoordPolymerSequenceWithoutCompId(a1['seq_id'], a1['atom_id'], index)
                        chainAssign2 = self.assignCoordPolymerSequenceWithoutCompId(a2['seq_id'], a2['atom_id'], index)
                        chainAssign3 = self.assignCoordPolymerSequenceWithoutCompId(a3['seq_id'], a3['atom_id'], index)

                    if len(chainAssign1) > 0 and len(chainAssign2) > 0 and len(chainAssign3) > 0:
                        self.selectCoordAtoms(chainAssign1, a1.get('auth_chain_id'),
                                              a1['seq_id'], a1['comp_id'], a1['atom_id'], index)
                        self.selectCoordAtoms(chainAssign2, a2.get('auth_chain_id'),
                                              a2['seq_id'], a2['comp_id'], a2['atom_id'], index)
                        self.selectCoordAtoms(chainAssign3, a3.get('auth_chain_id'),
                                              a3['seq_id'], a3['comp_id'], a3['atom_id'], index)

                        if len(self.atomSelectionSet) == self.num_of_dim:
                            has_assignments = True
                            has_assignments &=\
                                self.validateAtomType(1, self.atomSelectionSet[0][0]['atom_id'][0], dstFunc['position_1'])
                            has_assignments &=\
                                self.validateAtomType(2, self.atomSelectionSet[1][0]['atom_id'][0], dstFunc['position_2'])
                            has_assignments &=\
                                self.validateAtomType(3, self.atomSelectionSet[2][0]['atom_id'][0], dstFunc['position_3'])
                            if has_assignments:
                                self.atomSelectionSets.append([[dict(a) if a.__class__ is dict else a for a in sel]
                                                               for sel in self.atomSelectionSet])
                                self.asIsSets.append([asis1, asis2, asis3])
                                if self.reasons is not None:
                                    if 'onebond_idx_history' in self.reasons:
                                        onebond_idx = self.reasons['onebond_idx_history'][self.num_of_dim][self.cur_list_id]
                                        _atom1, _atom2 =\
                                            self.atomSelectionSet[DIM_TRANSFER_PAT_3D[onebond_idx][0]][0], \
                                            self.atomSelectionSet[DIM_TRANSFER_PAT_3D[onebond_idx][1]][0]
                                        if _atom1['chain_id'] != _atom2['chain_id']\
                                           or _atom1['seq_id'] != _atom2['seq_id']\
                                           or _atom1['atom_id'][0] == _atom2['atom_id'][0]\
                                           or not self.ccU.hasBond(_atom1['comp_id'], _atom1['atom_id'], _atom2['atom_id']):
                                            if not self.__canRemediatePeakAssignmentForOneBondTransfer(
                                                    _atom1, _atom2,
                                                    float(dstFunc[f'position_{DIM_TRANSFER_PAT_3D[onebond_idx][0] + 1}']),
                                                    float(dstFunc[f'position_{DIM_TRANSFER_PAT_3D[onebond_idx][1] + 1}'])):
                                                has_assignments = False
                                    if 'jcoupling_idx_history' in self.reasons:
                                        jcoupling_idx = self.reasons['jcoupling_idx_history'][self.num_of_dim][self.cur_list_id]
                                        if jcoupling_idx != -1:
                                            _atom1, _atom2 =\
                                                self.atomSelectionSet[DIM_TRANSFER_PAT_3D[jcoupling_idx][0]][0], \
                                                self.atomSelectionSet[DIM_TRANSFER_PAT_3D[jcoupling_idx][1]][0]
                                            if _atom1['chain_id'] != _atom2['chain_id']\
                                               or _atom1['seq_id'] != _atom2['seq_id']:
                                                if not self.__canRemediatePeakAssignmentForJcouplingTransfer(
                                                        _atom1, _atom2,
                                                        float(dstFunc[f'position_{DIM_TRANSFER_PAT_3D[jcoupling_idx][0] + 1}']),
                                                        float(dstFunc[f'position_{DIM_TRANSFER_PAT_3D[jcoupling_idx][1] + 1}'])):
                                                    has_assignments = False
                                    if 'relayed_idx_history' in self.reasons:
                                        relayed_idx = self.reasons['relayed_idx_history'][self.num_of_dim][self.cur_list_id]
                                        if jcoupling_idx != -1:
                                            _atom1, _atom2 =\
                                                self.atomSelectionSet[DIM_TRANSFER_PAT_3D[relayed_idx][0]][0], \
                                                self.atomSelectionSet[DIM_TRANSFER_PAT_3D[relayed_idx][1]][0]
                                            if _atom1['chain_id'] != _atom2['chain_id']\
                                               or _atom1['seq_id'] != _atom2['seq_id']:
                                                if not self.__canRemediatePeakAssignmentForRelayedTransfer(
                                                        _atom1, _atom2,
                                                        float(dstFunc[f'position_{DIM_TRANSFER_PAT_3D[relayed_idx][0] + 1}']),
                                                        float(dstFunc[f'position_{DIM_TRANSFER_PAT_3D[relayed_idx][1] + 1}'])):
                                                    has_assignments = False
                                if not has_long_range and has_assignments:
                                    if hasInterChainRestraint(self.atomSelectionSet):
                                        has_long_range = True
                                        if 'long_range' not in primary_specral_dim:
                                            primary_specral_dim['long_range'] = True
                                    else:
                                        for atom1, atom2, atom3 in itertools.product(self.atomSelectionSet[0],
                                                                                     self.atomSelectionSet[1],
                                                                                     self.atomSelectionSet[2]):
                                            has_long_range |= isLongRangeRestraint([atom1, atom2, atom3], self.polySeq)
                                            if has_long_range:
                                                if 'long_range' not in primary_specral_dim:
                                                    primary_specral_dim['long_range'] = True
                                                break
                            else:
                                break
                        else:
                            has_assignments = False
                            break

            except (KeyError, TypeError):
                pass

        return has_assignments, has_multiple_assignments, asis1, asis2, asis3

    def checkAssignments4D(self, index: int, assignments: List[List[dict]], dstFunc: dict
                           ) -> Tuple[bool, bool, Optional[bool], Optional[bool], Optional[bool], Optional[bool]]:
        """ Evaluate 4D peak assignment.
        """

        has_assignments = has_multiple_assignments = False
        asis1 = asis2 = asis3 = asis4 = None

        primary_specral_dim = self.spectral_dim[self.num_of_dim][self.cur_list_id][1]
        has_long_range = 'long_range' in primary_specral_dim

        if 'is_noesy' not in primary_specral_dim:
            file_name = self.__originalFileName.lower()
            alt_file_name = '' if self.spectrum_name is None else self.spectrum_name.lower()
            _file_names = (file_name, alt_file_name)
            primary_specral_dim['is_noesy'] = any('noe' in n for n in _file_names) or any('roe' in n for n in _file_names)
            if primary_specral_dim['is_noesy']:
                primary_specral_dim['long_range'] = True

        if all(assignment is not None for assignment in assignments):

            if self.reasons is not None and 'onebond_resolved' in self.reasons:
                _assignments = [None] * 4
                for k, v in self.reasons['onebond_resolved'].items():
                    _assignments[v] = assignments[k]
                assignments = _assignments

            self.retrieveLocalSeqScheme()

            try:

                hasChainId = all(all(_a['chain_id'] is not None for _a in a) for a in assignments)
                hasCompId = all(all(_a['comp_id'] is not None for _a in a) for a in assignments)

                has_multiple_assignments = any(True for assignment in assignments if len(assignment) > 1)

                pairs = []
                if len(assignments[0]) == len(assignments[1]) == len(assignments[2]) == len(assignments[3]):
                    for a1, a2, a3, a4 in zip(assignments[0], assignments[1], assignments[2], assignments[3]):
                        pairs.append((a1, a2, a3, a4))
                else:
                    for a1, a2, a3, a4 in itertools.product(assignments[0], assignments[1], assignments[2], assignments[3]):
                        pairs.append((a1, a2, a3, a4))

                for a1, a2, a3, a4 in pairs:

                    self.atomSelectionSet.clear()
                    asis1 = asis2 = asis3 = asis4 = None

                    if hasChainId and hasCompId:
                        chainAssign1, asis1 =\
                            self.assignCoordPolymerSequenceWithChainId(a1['chain_id'], a1['seq_id'], a1['comp_id'], a1['atom_id'],
                                                                       index)
                        chainAssign2, asis2 =\
                            self.assignCoordPolymerSequenceWithChainId(a2['chain_id'], a2['seq_id'], a2['comp_id'], a2['atom_id'],
                                                                       index)
                        chainAssign3, asis3 =\
                            self.assignCoordPolymerSequenceWithChainId(a3['chain_id'], a3['seq_id'], a3['comp_id'], a3['atom_id'],
                                                                       index)
                        chainAssign4, asis4 =\
                            self.assignCoordPolymerSequenceWithChainId(a4['chain_id'], a4['seq_id'], a4['comp_id'], a4['atom_id'],
                                                                       index)

                    elif hasChainId:
                        chainAssign1 =\
                            self.assignCoordPolymerSequenceWithChainIdWithoutCompId(a1['chain_id'], a1['seq_id'], a1['atom_id'],
                                                                                    index)
                        chainAssign2 =\
                            self.assignCoordPolymerSequenceWithChainIdWithoutCompId(a2['chain_id'], a2['seq_id'], a2['atom_id'],
                                                                                    index)
                        chainAssign3 =\
                            self.assignCoordPolymerSequenceWithChainIdWithoutCompId(a3['chain_id'], a3['seq_id'], a3['atom_id'],
                                                                                    index)
                        chainAssign4 =\
                            self.assignCoordPolymerSequenceWithChainIdWithoutCompId(a4['chain_id'], a4['seq_id'], a4['atom_id'],
                                                                                    index)

                    elif hasCompId:
                        chainAssign1, asis1 =\
                            self.assignCoordPolymerSequence(a1['chain_id'], a1['seq_id'], a1['comp_id'], a1['atom_id'], index)
                        chainAssign2, asis2 =\
                            self.assignCoordPolymerSequence(a2['chain_id'], a2['seq_id'], a2['comp_id'], a2['atom_id'], index)
                        chainAssign3, asis3 =\
                            self.assignCoordPolymerSequence(a3['chain_id'], a3['seq_id'], a3['comp_id'], a3['atom_id'], index)
                        chainAssign4, asis4 =\
                            self.assignCoordPolymerSequence(a4['chain_id'], a4['seq_id'], a4['comp_id'], a4['atom_id'], index)

                    else:
                        chainAssign1 = self.assignCoordPolymerSequenceWithoutCompId(a1['seq_id'], a1['atom_id'], index)
                        chainAssign2 = self.assignCoordPolymerSequenceWithoutCompId(a2['seq_id'], a2['atom_id'], index)
                        chainAssign3 = self.assignCoordPolymerSequenceWithoutCompId(a3['seq_id'], a3['atom_id'], index)
                        chainAssign4 = self.assignCoordPolymerSequenceWithoutCompId(a4['seq_id'], a4['atom_id'], index)

                    if len(chainAssign1) > 0 and len(chainAssign2) > 0 and len(chainAssign3) > 0 and len(chainAssign4) > 0:
                        self.selectCoordAtoms(chainAssign1, a1.get('auth_chain_id'),
                                              a1['seq_id'], a1['comp_id'], a1['atom_id'], index)
                        self.selectCoordAtoms(chainAssign2, a2.get('auth_chain_id'),
                                              a2['seq_id'], a2['comp_id'], a2['atom_id'], index)
                        self.selectCoordAtoms(chainAssign3, a3.get('auth_chain_id'),
                                              a3['seq_id'], a3['comp_id'], a3['atom_id'], index)
                        self.selectCoordAtoms(chainAssign4, a4.get('auth_chain_id'),
                                              a4['seq_id'], a4['comp_id'], a4['atom_id'], index)

                        if len(self.atomSelectionSet) == self.num_of_dim:
                            has_assignments = True
                            has_assignments &=\
                                self.validateAtomType(1, self.atomSelectionSet[0][0]['atom_id'][0], dstFunc['position_1'])
                            has_assignments &=\
                                self.validateAtomType(2, self.atomSelectionSet[1][0]['atom_id'][0], dstFunc['position_2'])
                            has_assignments &=\
                                self.validateAtomType(3, self.atomSelectionSet[2][0]['atom_id'][0], dstFunc['position_3'])
                            has_assignments &=\
                                self.validateAtomType(4, self.atomSelectionSet[3][0]['atom_id'][0], dstFunc['position_4'])
                            if has_assignments:
                                self.atomSelectionSets.append([[dict(a) if a.__class__ is dict else a for a in sel]
                                                               for sel in self.atomSelectionSet])
                                self.asIsSets.append([asis1, asis2, asis3, asis4])
                                if self.reasons is not None:
                                    if 'onebond_idx_history' in self.reasons:
                                        onebond_idx = self.reasons['onebond_idx_history'][self.num_of_dim][self.cur_list_id]
                                        _atom1, _atom2, _atom3, _atom4 =\
                                            self.atomSelectionSet[DIM_TRANSFER_PAT_4D[onebond_idx][0][0]][0], \
                                            self.atomSelectionSet[DIM_TRANSFER_PAT_4D[onebond_idx][0][1]][0], \
                                            self.atomSelectionSet[DIM_TRANSFER_PAT_4D[onebond_idx][1][0]][0], \
                                            self.atomSelectionSet[DIM_TRANSFER_PAT_4D[onebond_idx][1][1]][0]
                                        if _atom1['chain_id'] != _atom2['chain_id']\
                                           or _atom1['seq_id'] != _atom2['seq_id']\
                                           or _atom1['atom_id'][0] == _atom2['atom_id'][0]\
                                           or _atom1['atom_id'][0] == _atom2['atom_id'][0]\
                                           or not self.ccU.hasBond(_atom1['comp_id'], _atom1['atom_id'], _atom2['atom_id']):
                                            if not self.__canRemediatePeakAssignmentForOneBondTransfer(
                                                    _atom1, _atom2,
                                                    float(dstFunc[f'position_{DIM_TRANSFER_PAT_4D[onebond_idx][0][0] + 1}']),
                                                    float(dstFunc[f'position_{DIM_TRANSFER_PAT_4D[onebond_idx][0][1] + 1}'])):
                                                has_assignments = False
                                        if _atom3['chain_id'] != _atom4['chain_id']\
                                           or _atom3['seq_id'] != _atom4['seq_id']\
                                           or _atom3['atom_id'][0] == _atom4['atom_id'][0]\
                                           or _atom3['atom_id'][0] == _atom4['atom_id'][0]\
                                           or not self.ccU.hasBond(_atom3['comp_id'], _atom3['atom_id'], _atom4['atom_id']):
                                            if not self.__canRemediatePeakAssignmentForOneBondTransfer(
                                                    _atom3, _atom4,
                                                    float(dstFunc[f'position_{DIM_TRANSFER_PAT_4D[onebond_idx][1][0] + 1}']),
                                                    float(dstFunc[f'position_{DIM_TRANSFER_PAT_4D[onebond_idx][1][1] + 1}'])):
                                                has_assignments = False
                                    if 'jcoupling_idx_history' in self.reasons:
                                        jcoupling_idx = self.reasons['jcoupling_idx_history'][self.num_of_dim][self.cur_list_id]
                                        if jcoupling_idx != -1:
                                            _atom1, _atom2, _atom3, _atom4 =\
                                                self.atomSelectionSet[DIM_TRANSFER_PAT_4D[jcoupling_idx][0][0]][0], \
                                                self.atomSelectionSet[DIM_TRANSFER_PAT_4D[jcoupling_idx][0][1]][0], \
                                                self.atomSelectionSet[DIM_TRANSFER_PAT_4D[jcoupling_idx][1][0]][0], \
                                                self.atomSelectionSet[DIM_TRANSFER_PAT_4D[jcoupling_idx][1][1]][0]
                                            if _atom1['chain_id'] != _atom2['chain_id']\
                                               or _atom1['seq_id'] != _atom2['seq_id']:
                                                if not self.__canRemediatePeakAssignmentForJcouplingTransfer(
                                                        _atom1, _atom2,
                                                        float(dstFunc[f'position_{DIM_TRANSFER_PAT_4D[jcoupling_idx][0][0] + 1}']),
                                                        float(dstFunc[f'position_{DIM_TRANSFER_PAT_4D[jcoupling_idx][0][1] + 1}'])):
                                                    has_assignments = False
                                            if _atom3['chain_id'] != _atom4['chain_id']\
                                               or _atom3['seq_id'] != _atom4['seq_id']:
                                                if not self.__canRemediatePeakAssignmentForJcouplingTransfer(
                                                        _atom3, _atom4,
                                                        float(dstFunc[f'position_{DIM_TRANSFER_PAT_4D[jcoupling_idx][1][0] + 1}']),
                                                        float(dstFunc[f'position_{DIM_TRANSFER_PAT_4D[jcoupling_idx][1][1] + 1}'])):
                                                    has_assignments = False
                                    if 'relayed_idx_history' in self.reasons:
                                        relayed_idx = self.reasons['relayed_idx_history'][self.num_of_dim][self.cur_list_id]
                                        if jcoupling_idx != -1:
                                            _atom1, _atom2, _atom3, _atom4 =\
                                                self.atomSelectionSet[DIM_TRANSFER_PAT_4D[relayed_idx][0][0]][0], \
                                                self.atomSelectionSet[DIM_TRANSFER_PAT_4D[relayed_idx][0][1]][0], \
                                                self.atomSelectionSet[DIM_TRANSFER_PAT_4D[relayed_idx][1][0]][0], \
                                                self.atomSelectionSet[DIM_TRANSFER_PAT_4D[relayed_idx][1][1]][0]
                                            if _atom1['chain_id'] != _atom2['chain_id']\
                                               or _atom1['seq_id'] != _atom2['seq_id']:
                                                if not self.__canRemediatePeakAssignmentForRelayedTransfer(
                                                        _atom1, _atom2,
                                                        float(dstFunc[f'position_{DIM_TRANSFER_PAT_4D[relayed_idx][0][0] + 1}']),
                                                        float(dstFunc[f'position_{DIM_TRANSFER_PAT_4D[relayed_idx][0][1] + 1}'])):
                                                    has_assignments = False
                                            if _atom3['chain_id'] != _atom4['chain_id']\
                                               or _atom3['seq_id'] != _atom4['seq_id']:
                                                if not self.__canRemediatePeakAssignmentForRelayedTransfer(
                                                        _atom3, _atom4,
                                                        float(dstFunc[f'position_{DIM_TRANSFER_PAT_4D[relayed_idx][1][0] + 1}']),
                                                        float(dstFunc[f'position_{DIM_TRANSFER_PAT_4D[relayed_idx][1][1] + 1}'])):
                                                    has_assignments = False
                                if not has_long_range and has_assignments:
                                    if hasInterChainRestraint(self.atomSelectionSet):
                                        has_long_range = True
                                        if 'long_range' not in primary_specral_dim:
                                            primary_specral_dim['long_range'] = True
                                    else:
                                        for atom1, atom2, atom3, atom4 in itertools.product(self.atomSelectionSet[0],
                                                                                            self.atomSelectionSet[1],
                                                                                            self.atomSelectionSet[2],
                                                                                            self.atomSelectionSet[3]):
                                            has_long_range |= isLongRangeRestraint([atom1, atom2, atom3, atom4], self.polySeq)
                                            if has_long_range:
                                                if 'long_range' not in primary_specral_dim:
                                                    primary_specral_dim['long_range'] = True
                                                break
                            else:
                                break
                        else:
                            has_assignments = False
                            break

            except (KeyError, TypeError):
                pass

        return has_assignments, has_multiple_assignments, asis1, asis2, asis3, asis4

    def __extractCommonAtom(self, atom_sel: List[dict]) -> dict:
        """ Return standardized pseudoatom name from given atom selections.
        """

        if len(atom_sel) == 0:
            return {}

        if len(atom_sel) == 1:
            return atom_sel[0]

        strings = [a['atom_id'] for a in atom_sel]

        min_str = min(strings, key=len)
        max_str = max(strings, key=len)

        len_min_str = len(min_str)
        len_max_str = len(max_str)
        longest_substr = ''

        for i in range(len_min_str):
            for j in range(i + 1, len_min_str + 1):
                substr = min_str[i:j]
                if all(substr in s for s in strings):
                    if len(substr) > len(longest_substr):
                        longest_substr = substr

        if len(longest_substr) == 0:
            return atom_sel[0]

        self.pA.clear()
        self.pA.setReferenceSequence(list(longest_substr), 'REFNAME')
        self.pA.addTestSequence(list(max_str), 'NAME')
        self.pA.doAlign()

        myAlign = self.pA.getAlignment('NAME')

        length = len(myAlign)

        if length == 0:
            return atom_sel[0]

        common_name = []

        for i in range(length):
            myPr = myAlign[i]
            myPr0 = str(myPr[0])
            myPr1 = str(myPr[1])
            if myPr0 == myPr1:
                if myPr0 not in EMPTY_VALUE:
                    common_name.append(myPr0)
            elif myPr0 in EMPTY_VALUE:
                if myPr1 not in EMPTY_VALUE:
                    common_name.append(('#' if myPr1.isdigit() else '%') if len_min_str == len_max_str else '*')

        if len(common_name) == 0:
            return atom_sel[0]

        common_name = ''.join(common_name)

        while '##' in common_name:
            common_name = common_name.replace('##', '*')

        while '%%' in common_name:
            common_name = common_name.replace('%%', '*')

        while '*%' in common_name:
            common_name = common_name.replace('*%', '*')

        while '%*' in common_name:
            common_name = common_name.replace('%*', '*')

        while '**' in common_name:
            common_name = common_name.replace('**', '*')

        ambig_code = 1
        if any(atom1['chain_id'] != atom2['chain_id'] for atom1, atom2 in itertools.combinations(atom_sel, 2)):
            ambig_code = 6
        elif any(atom1['seq_id'] != atom2['seq_id'] for atom1, atom2 in itertools.combinations(atom_sel, 2)):
            ambig_code = 5
        else:
            comp_id = atom_sel[0]['comp_id']
            for atom1, atom2 in itertools.combinations(atom_sel, 2):
                atom_id1 = atom1['atom_id']
                atom_id2 = atom2['atom_id']
                ambig_code1 = self.csStat.getMaxAmbigCodeWoSetId(comp_id, atom_id1)
                ambig_code2 = self.csStat.getMaxAmbigCodeWoSetId(comp_id, atom_id2)
                if ambig_code1 != ambig_code2:
                    ambig_code = 4
                    break
                if ambig_code1 == 1:
                    if atom_id2 in self.csStat.getProtonsInSameGroup(comp_id, atom_id1, excl_self=True):
                        continue
                    ambig_code = 4
                    break
                if ambig_code1 == 2\
                   and atom_id2 in self.csStat.getProtonsInSameGroup(comp_id, atom_id1, excl_self=True):
                    continue
                if ambig_code1 in (2, 3):
                    _atom_id2 = self.csStat.getGeminalAtom(comp_id, atom_id1)
                    if _atom_id2 is None:
                        ambig_code = 4
                        break
                    if _atom_id2 == atom_id2\
                       or (ambig_code1 == 2 and atom_id2 in self.csStat.getProtonsInSameGroup(comp_id, _atom_id2, excl_self=True)):
                        continue
                    ambig_code = 4
                    break

        if len(common_name) == 2 and common_name.endswith('*'):  # avoid 'H*' for amide proton
            if self.csStat.getTypeOfCompId(atom_sel[0]['comp_id'])[0]:
                return atom_sel[0]

        _atom_sel = copy.copy(atom_sel[0])
        if 'auth_atom_id' in _atom_sel:
            _atom_sel['orig_atom_id'] = _atom_sel['auth_atom_id']
        _atom_sel['auth_atom_id'] = common_name
        if ambig_code != 1:
            _atom_sel['ambig_code'] = ambig_code

        return _atom_sel

    def addAssignedPkRow2D(self, index: int, dstFunc: dict, has_assignments: bool, has_multiple_assignments: bool,
                           asis1: Optional[bool], asis2: Optional[bool],
                           debug_label: Optional[str], details: Optional[str]) -> None:
        """ Add assigned peak list row (2D).
        """

        if self.__debug:
            if not has_assignments:
                print(f"subtype={self.cur_subtype} id={self.peaks2D} (index={index}) "
                      f"{debug_label}None None {dstFunc}")
            for idx, atomSelectionSet in enumerate(self.atomSelectionSets, start=1):
                if has_multiple_assignments:
                    print(f"subtype={self.cur_subtype} id={self.peaks2D} (index={index}) combination_id={idx} "
                          f"{debug_label}{atomSelectionSet[0]} "
                          f"{atomSelectionSet[1]} {dstFunc}")
                else:
                    print(f"subtype={self.cur_subtype} id={self.peaks2D} (index={index}) "
                          f"{debug_label}{atomSelectionSet[0]} "
                          f"{atomSelectionSet[1]} {dstFunc}")

        if self.__createSfDict:
            sf = self.getSf()

            if sf is not None:
                sf['id'] = index

                if has_assignments and has_multiple_assignments:
                    history = self.onebond_idx_history[self.num_of_dim][self.cur_list_id]
                    if self.software_name != 'PIPP':
                        for atomSelectionSet in self.atomSelectionSets:
                            for onebond_idx, (dim1, dim2) in enumerate(DIM_TRANSFER_PAT_2D):
                                _atom1, _atom2 = atomSelectionSet[dim1][0], atomSelectionSet[dim2][0]
                                if _atom1['chain_id'] != _atom2['chain_id']\
                                   or _atom1['seq_id'] != _atom2['seq_id']\
                                   or _atom1['atom_id'][0] == _atom2['atom_id'][0]:
                                    continue
                                if self.ccU.hasBond(_atom1['comp_id'], _atom1['atom_id'], _atom2['atom_id']):
                                    history.append(onebond_idx)
                    atom_set1, atom_set2 = [], []
                    for atomSelectionSet in self.atomSelectionSets:
                        valid = True
                        if len(history) > 0:
                            onebond_idx = collections.Counter(history).most_common()[0][0]
                            _atom1, _atom2 =\
                                atomSelectionSet[DIM_TRANSFER_PAT_2D[onebond_idx][0]][0], \
                                atomSelectionSet[DIM_TRANSFER_PAT_2D[onebond_idx][1]][0]
                            if _atom1['chain_id'] != _atom2['chain_id']\
                               or _atom1['seq_id'] != _atom2['seq_id']\
                               or _atom1['atom_id'][0] == _atom2['atom_id'][0]\
                               or not self.ccU.hasBond(_atom1['comp_id'], _atom1['atom_id'], _atom2['atom_id']):
                                valid = False
                        if not valid:
                            continue
                        atom_set1.extend(atomSelectionSet[0])
                        atom_set2.extend(atomSelectionSet[1])
                    common_atom1 = self.__extractCommonAtom([dict(s) for s in dict.fromkeys(frozenset(a.items())
                                                                                            for a in atom_set1
                                                                                            if isinstance(a, dict))]
                                                            if len(atom_set1) > 1 else atom_set1)
                    common_atom2 = self.__extractCommonAtom([dict(s) for s in dict.fromkeys(frozenset(a.items())
                                                                                            for a in atom_set2
                                                                                            if isinstance(a, dict))]
                                                            if len(atom_set2) > 1 else atom_set2)
                    for atomSelectionSet, asIsSet in zip(self.atomSelectionSets, self.asIsSets):
                        ambig_code1 = ambig_code2 = None
                        atom1 = self.__extractCommonAtom(atomSelectionSet[0])
                        atom2 = self.__extractCommonAtom(atomSelectionSet[1])
                        asis1, asis2 = asIsSet
                        if len(atomSelectionSet[0]) > 1:
                            ambig_code1 = self.csStat.getMaxAmbigCodeWoSetId(atom1['comp_id'], atom1['atom_id'])
                            if ambig_code1 == 0:
                                ambig_code1 = None
                        if 'ambig_code' in common_atom1:
                            ambig_code1 = common_atom1['ambig_code']
                        if len(atomSelectionSet[1]) > 1:
                            ambig_code2 = self.csStat.getMaxAmbigCodeWoSetId(atom2['comp_id'], atom2['atom_id'])
                            if ambig_code2 == 0:
                                ambig_code2 = None
                        if 'ambig_code' in common_atom2:
                            ambig_code2 = common_atom2['ambig_code']

                        sf['row_index_id'] += 1

                        row = getPkRow(self.cur_subtype, sf['id'], sf['row_index_id'],
                                       sf['list_id'], self.__entryId, dstFunc,
                                       self.__authToStarSeq, self.__authToOrigSeq, self.__offsetHolder,
                                       atom1, atom2, asis1=asis1, asis2=asis2,
                                       ambig_code1=ambig_code1, ambig_code2=ambig_code2,
                                       details=details)
                        sf['loop'].add_data(row)

                else:

                    sf['row_index_id'] += 1

                    ambig_code1 = ambig_code2 = None
                    if has_assignments:
                        atom1 = self.__extractCommonAtom(self.atomSelectionSet[0])
                        atom2 = self.__extractCommonAtom(self.atomSelectionSet[1])
                        if len(self.atomSelectionSet[0]) > 1:
                            ambig_code1 = self.csStat.getMaxAmbigCodeWoSetId(atom1['comp_id'], atom1['atom_id'])
                            if ambig_code1 == 0:
                                ambig_code1 = None
                        if 'ambig_code' in atom1:
                            ambig_code1 = atom1['ambig_code']
                        if len(self.atomSelectionSet[1]) > 1:
                            ambig_code2 = self.csStat.getMaxAmbigCodeWoSetId(atom2['comp_id'], atom2['atom_id'])
                            if ambig_code2 == 0:
                                ambig_code2 = None
                        if 'ambig_code' in atom2:
                            ambig_code2 = atom2['ambig_code']
                    else:
                        atom1 = atom2 = None

                    row = getPkRow(self.cur_subtype, sf['id'], sf['row_index_id'],
                                   sf['list_id'], self.__entryId, dstFunc,
                                   self.__authToStarSeq, self.__authToOrigSeq, self.__offsetHolder,
                                   atom1, atom2, asis1=asis1, asis2=asis2,
                                   ambig_code1=ambig_code1, ambig_code2=ambig_code2,
                                   details=details)
                    sf['loop'].add_data(row)

                if not has_assignments and details is not None:
                    self.f.append(f"[Conflicted peak assignment] {self.getCurrentSpectralPeak(n=index)}"
                                  f"Peak assignments {details!r} could not map to the coordinates.")

                sf['index_id'] += 1

                row = getAltPkRow(self.cur_subtype, sf['index_id'], sf['id'], sf['list_id'], self.__entryId, dstFunc)
                if row is not None:
                    sf['alt_loops'][0].add_data(row)

                row = getPkGenCharRow(self.cur_subtype, sf['id'], sf['list_id'], self.__entryId, dstFunc, 'volume')
                if row is not None:
                    sf['alt_loops'][1].add_data(row)
                row = getPkGenCharRow(self.cur_subtype, sf['id'], sf['list_id'], self.__entryId, dstFunc, 'height')
                if row is not None:
                    sf['alt_loops'][1].add_data(row)

                for idx in range(self.num_of_dim):
                    row = getPkCharRow(self.cur_subtype, sf['id'], sf['list_id'], self.__entryId, dstFunc, idx + 1)
                    sf['alt_loops'][2].add_data(row)
                if has_assignments:
                    history = self.onebond_idx_history[self.num_of_dim][self.cur_list_id]
                    if self.software_name != 'PIPP':
                        for atomSelectionSet in self.atomSelectionSets:
                            for onebond_idx, (dim1, dim2) in enumerate(DIM_TRANSFER_PAT_2D):
                                _atom1, _atom2 = atomSelectionSet[dim1][0], atomSelectionSet[dim2][0]
                                if _atom1['chain_id'] != _atom2['chain_id']\
                                   or _atom1['seq_id'] != _atom2['seq_id']\
                                   or _atom1['atom_id'][0] == _atom2['atom_id'][0]:
                                    continue
                                if self.ccU.hasBond(_atom1['comp_id'], _atom1['atom_id'], _atom2['atom_id']):
                                    history.append(onebond_idx)
                    set_id = None if len(self.atomSelectionSets) < 2 else 0
                    if set_id is None:
                        for atomSelectionSet, asIsSet in zip(self.atomSelectionSets, self.asIsSets):
                            uniqAtoms = []
                            for idx in range(self.num_of_dim):
                                atom_sel = atomSelectionSet[idx]
                                common_atom = self.__extractCommonAtom(atom_sel)
                                if common_atom not in uniqAtoms:
                                    asis = asIsSet[idx]
                                    ambig_code = None
                                    if len(atom_sel) > 1:
                                        comp_id = common_atom['comp_id']
                                        rep_atom_id = atom_sel[0]['atom_id']
                                        ambig_code = self.csStat.getMaxAmbigCodeWoSetId(comp_id, rep_atom_id)
                                        if ambig_code == 0:
                                            ambig_code = None
                                        elif ambig_code == 2:
                                            hvy_grp_atoms, _ = self.nefT.get_group(comp_id, rep_atom_id)
                                            hvy_gem_atoms, pro_gem_atoms = self.nefT.get_geminal_group(comp_id, rep_atom_id)
                                            if None not in (hvy_grp_atoms, hvy_gem_atoms):
                                                gem_atom_ids = pro_gem_atoms if rep_atom_id[0] in PROTON_BEGIN_CODE\
                                                    else hvy_gem_atoms
                                                if not any(a['atom_id'] in gem_atom_ids for a in atom_sel):
                                                    ambig_code = 1
                                    if 'ambig_code' in common_atom:
                                        ambig_code = common_atom['ambig_code']
                                    row = getPkChemShiftRow(self.cur_subtype, sf['id'], sf['list_id'], self.__entryId,
                                                            dstFunc, set_id, idx + 1,
                                                            self.__authToStarSeq, self.__authToOrigSeq, self.__offsetHolder,
                                                            common_atom, asis, ambig_code)
                                    if row is None:
                                        continue
                                    sf['alt_loops'][3].add_data(row)
                                    uniqAtoms.append(common_atom)
                    else:
                        atom_set1, atom_set2 = [], []
                        for atomSelectionSet in self.atomSelectionSets:
                            valid = True
                            if len(history) > 0:
                                onebond_idx = collections.Counter(history).most_common()[0][0]
                                _atom1, _atom2 =\
                                    atomSelectionSet[DIM_TRANSFER_PAT_2D[onebond_idx][0]][0], \
                                    atomSelectionSet[DIM_TRANSFER_PAT_2D[onebond_idx][1]][0]
                                if _atom1['chain_id'] != _atom2['chain_id']\
                                   or _atom1['seq_id'] != _atom2['seq_id']\
                                   or _atom1['atom_id'][0] == _atom2['atom_id'][0]\
                                   or not self.ccU.hasBond(_atom1['comp_id'], _atom1['atom_id'], _atom2['atom_id']):
                                    valid = False
                            if not valid:
                                continue
                            atom_set1.extend(atomSelectionSet[0])
                            atom_set2.extend(atomSelectionSet[1])
                        common_atoms = []
                        common_atoms.append(self.__extractCommonAtom([dict(s) for s in dict.fromkeys(frozenset(a.items())
                                                                                                     for a in atom_set1
                                                                                                     if isinstance(a, dict))]
                                                                     if len(atom_set1) > 1 else atom_set1))
                        common_atoms.append(self.__extractCommonAtom([dict(s) for s in dict.fromkeys(frozenset(a.items())
                                                                                                     for a in atom_set2
                                                                                                     if isinstance(a, dict))]
                                                                     if len(atom_set2) > 1 else atom_set2))
                        set_id = 1
                        for atomSelectionSet, asIsSet in itertools.zip_longest(self.atomSelectionSets, self.asIsSets):
                            valid = True
                            if len(history) > 0:
                                onebond_idx = collections.Counter(history).most_common()[0][0]
                                _atom1, _atom2 =\
                                    atomSelectionSet[DIM_TRANSFER_PAT_2D[onebond_idx][0]][0], \
                                    atomSelectionSet[DIM_TRANSFER_PAT_2D[onebond_idx][1]][0]
                                if _atom1['chain_id'] != _atom2['chain_id']\
                                   or _atom1['seq_id'] != _atom2['seq_id']\
                                   or _atom1['atom_id'][0] == _atom2['atom_id'][0]\
                                   or not self.ccU.hasBond(_atom1['comp_id'], _atom1['atom_id'], _atom2['atom_id']):
                                    valid = False
                            if not valid:
                                continue
                            uniqAtoms = []
                            for idx in range(self.num_of_dim):
                                atom_sel = atomSelectionSet[idx]
                                common_atom = self.__extractCommonAtom(atom_sel)
                                if common_atom not in uniqAtoms:
                                    asis = asIsSet[idx]
                                    ambig_code = None
                                    if len(atom_sel) > 1:
                                        comp_id = common_atom['comp_id']
                                        rep_atom_id = atom_sel[0]['atom_id']
                                        ambig_code = self.csStat.getMaxAmbigCodeWoSetId(comp_id, rep_atom_id)
                                        if ambig_code == 0:
                                            ambig_code = None
                                        elif ambig_code == 2:
                                            hvy_grp_atoms, _ = self.nefT.get_group(comp_id, rep_atom_id)
                                            hvy_gem_atoms, pro_gem_atoms = self.nefT.get_geminal_group(comp_id, rep_atom_id)
                                            if None not in (hvy_grp_atoms, hvy_gem_atoms):
                                                gem_atom_ids = pro_gem_atoms if rep_atom_id[0] in PROTON_BEGIN_CODE\
                                                    else hvy_gem_atoms
                                                if not any(a['atom_id'] in gem_atom_ids for a in atom_sel):
                                                    ambig_code = 1
                                    if 'ambig_code' in common_atoms[idx]:
                                        ambig_code = common_atoms[idx]['ambig_code']
                                    row = getPkChemShiftRow(self.cur_subtype, sf['id'], sf['list_id'], self.__entryId,
                                                            dstFunc, set_id, idx + 1,
                                                            self.__authToStarSeq, self.__authToOrigSeq, self.__offsetHolder,
                                                            common_atom, asis, ambig_code)
                                    if row is None:
                                        continue
                                    sf['alt_loops'][3].add_data(row)
                                    uniqAtoms.append(common_atom)
                            set_id += 1

    def addAssignedPkRow3D(self, index: int, dstFunc: dict, has_assignments: bool, has_multiple_assignments: bool,
                           asis1: Optional[bool], asis2: Optional[bool], asis3: Optional[bool],
                           debug_label: Optional[str], details: Optional[str]) -> None:
        """ Add assigned peak list row (3D).
        """

        if self.__debug:
            if not has_assignments:
                print(f"subtype={self.cur_subtype} id={self.peaks3D} (index={index}) "
                      f"{debug_label}None None None {dstFunc}")
            for idx, atomSelectionSet in enumerate(self.atomSelectionSets, start=1):
                if has_multiple_assignments:
                    print(f"subtype={self.cur_subtype} id={self.peaks3D} (index={index}) combination_id={idx} "
                          f"{debug_label}{atomSelectionSet[0]} "
                          f"{atomSelectionSet[1]} "
                          f"{atomSelectionSet[2]} {dstFunc}")
                else:
                    print(f"subtype={self.cur_subtype} id={self.peaks3D} (index={index}) "
                          f"{debug_label}{atomSelectionSet[0]} "
                          f"{atomSelectionSet[1]} "
                          f"{atomSelectionSet[2]} {dstFunc}")

        if self.__createSfDict:
            sf = self.getSf()

            if sf is not None:
                sf['id'] = index

                if has_assignments and has_multiple_assignments:
                    history = self.onebond_idx_history[self.num_of_dim][self.cur_list_id]
                    if self.software_name != 'PIPP':
                        for atomSelectionSet in self.atomSelectionSets:
                            for onebond_idx, (dim1, dim2) in enumerate(DIM_TRANSFER_PAT_3D):
                                _atom1, _atom2 = atomSelectionSet[dim1][0], atomSelectionSet[dim2][0]
                                if _atom1['chain_id'] != _atom2['chain_id']\
                                   or _atom1['seq_id'] != _atom2['seq_id']\
                                   or _atom1['atom_id'][0] == _atom2['atom_id'][0]:
                                    continue
                                if self.ccU.hasBond(_atom1['comp_id'], _atom1['atom_id'], _atom2['atom_id']):
                                    history.append(onebond_idx)
                    atom_set1, atom_set2, atom_set3 = [], [], []
                    for atomSelectionSet in self.atomSelectionSets:
                        valid = True
                        if len(history) > 0:
                            onebond_idx = collections.Counter(history).most_common()[0][0]
                            _atom1, _atom2 =\
                                atomSelectionSet[DIM_TRANSFER_PAT_3D[onebond_idx][0]][0], \
                                atomSelectionSet[DIM_TRANSFER_PAT_3D[onebond_idx][1]][0]
                            if _atom1['chain_id'] != _atom2['chain_id']\
                               or _atom1['seq_id'] != _atom2['seq_id']\
                               or _atom1['atom_id'][0] == _atom2['atom_id'][0]\
                               or not self.ccU.hasBond(_atom1['comp_id'], _atom1['atom_id'], _atom2['atom_id']):
                                valid = False
                        if not valid:
                            continue
                        atom_set1.extend(atomSelectionSet[0])
                        atom_set2.extend(atomSelectionSet[1])
                        atom_set3.extend(atomSelectionSet[2])
                    common_atom1 = self.__extractCommonAtom([dict(s) for s in dict.fromkeys(frozenset(a.items())
                                                                                            for a in atom_set1
                                                                                            if isinstance(a, dict))]
                                                            if len(atom_set1) > 1 else atom_set1)
                    common_atom2 = self.__extractCommonAtom([dict(s) for s in dict.fromkeys(frozenset(a.items())
                                                                                            for a in atom_set2
                                                                                            if isinstance(a, dict))]
                                                            if len(atom_set2) > 1 else atom_set2)
                    common_atom3 = self.__extractCommonAtom([dict(s) for s in dict.fromkeys(frozenset(a.items())
                                                                                            for a in atom_set3
                                                                                            if isinstance(a, dict))]
                                                            if len(atom_set3) > 1 else atom_set3)
                    for atomSelectionSet, asIsSet in zip(self.atomSelectionSets, self.asIsSets):
                        ambig_code1 = ambig_code2 = ambig_code3 = None
                        atom1 = self.__extractCommonAtom(atomSelectionSet[0])
                        atom2 = self.__extractCommonAtom(atomSelectionSet[1])
                        atom3 = self.__extractCommonAtom(atomSelectionSet[2])
                        asis1, asis2, asis3 = asIsSet
                        if len(atomSelectionSet[0]) > 1:
                            ambig_code1 = self.csStat.getMaxAmbigCodeWoSetId(atom1['comp_id'], atom1['atom_id'])
                            if ambig_code1 == 0:
                                ambig_code1 = None
                        if 'ambig_code' in common_atom1:
                            ambig_code1 = common_atom1['ambig_code']
                        if len(atomSelectionSet[1]) > 1:
                            ambig_code2 = self.csStat.getMaxAmbigCodeWoSetId(atom2['comp_id'], atom2['atom_id'])
                            if ambig_code2 == 0:
                                ambig_code2 = None
                        if 'ambig_code' in common_atom2:
                            ambig_code2 = common_atom2['ambig_code']
                        if len(atomSelectionSet[2]) > 1:
                            ambig_code3 = self.csStat.getMaxAmbigCodeWoSetId(atom3['comp_id'], atom3['atom_id'])
                            if ambig_code3 == 0:
                                ambig_code3 = None
                        if 'ambig_code' in common_atom3:
                            ambig_code3 = common_atom3['ambig_code']

                        sf['row_index_id'] += 1

                        row = getPkRow(self.cur_subtype, sf['id'], sf['row_index_id'],
                                       sf['list_id'], self.__entryId, dstFunc,
                                       self.__authToStarSeq, self.__authToOrigSeq, self.__offsetHolder,
                                       atom1, atom2, atom3, asis1=asis1, asis2=asis2, asis3=asis3,
                                       ambig_code1=ambig_code1, ambig_code2=ambig_code2,
                                       ambig_code3=ambig_code3,
                                       details=details)
                        sf['loop'].add_data(row)

                else:

                    sf['row_index_id'] += 1

                    ambig_code1 = ambig_code2 = ambig_code3 = None
                    if has_assignments:
                        atom1 = self.__extractCommonAtom(self.atomSelectionSet[0])
                        atom2 = self.__extractCommonAtom(self.atomSelectionSet[1])
                        atom3 = self.__extractCommonAtom(self.atomSelectionSet[2])
                        if len(self.atomSelectionSet[0]) > 1:
                            ambig_code1 = self.csStat.getMaxAmbigCodeWoSetId(atom1['comp_id'], atom1['atom_id'])
                            if ambig_code1 == 0:
                                ambig_code1 = None
                        if 'ambig_code' in atom1:
                            ambig_code1 = atom1['ambig_code']
                        if len(self.atomSelectionSet[1]) > 1:
                            ambig_code2 = self.csStat.getMaxAmbigCodeWoSetId(atom2['comp_id'], atom2['atom_id'])
                            if ambig_code2 == 0:
                                ambig_code2 = None
                        if 'ambig_code' in atom2:
                            ambig_code2 = atom2['ambig_code']
                        if len(self.atomSelectionSet[2]) > 1:
                            ambig_code3 = self.csStat.getMaxAmbigCodeWoSetId(atom3['comp_id'], atom3['atom_id'])
                            if ambig_code3 == 0:
                                ambig_code3 = None
                        if 'ambig_code' in atom3:
                            ambig_code3 = atom3['ambig_code']
                    else:
                        atom1 = atom2 = atom3 = None

                    row = getPkRow(self.cur_subtype, sf['id'], sf['row_index_id'],
                                   sf['list_id'], self.__entryId, dstFunc,
                                   self.__authToStarSeq, self.__authToOrigSeq, self.__offsetHolder,
                                   atom1, atom2, atom3, asis1=asis1, asis2=asis2, asis3=asis3,
                                   ambig_code1=ambig_code1, ambig_code2=ambig_code2,
                                   ambig_code3=ambig_code3,
                                   details=details)
                    sf['loop'].add_data(row)

                if not has_assignments and details is not None:
                    self.f.append(f"[Conflicted peak assignment] {self.getCurrentSpectralPeak(n=index)}"
                                  f"Peak assignments {details!r} could not map to the coordinates.")

                sf['index_id'] += 1

                row = getAltPkRow(self.cur_subtype, sf['index_id'], sf['id'], sf['list_id'], self.__entryId, dstFunc)
                if row is not None:
                    sf['alt_loops'][0].add_data(row)

                row = getPkGenCharRow(self.cur_subtype, sf['id'], sf['list_id'], self.__entryId, dstFunc, 'volume')
                if row is not None:
                    sf['alt_loops'][1].add_data(row)
                row = getPkGenCharRow(self.cur_subtype, sf['id'], sf['list_id'], self.__entryId, dstFunc, 'height')
                if row is not None:
                    sf['alt_loops'][1].add_data(row)

                for idx in range(self.num_of_dim):
                    row = getPkCharRow(self.cur_subtype, sf['id'], sf['list_id'], self.__entryId, dstFunc, idx + 1)
                    sf['alt_loops'][2].add_data(row)
                if has_assignments:
                    history = self.onebond_idx_history[self.num_of_dim][self.cur_list_id]
                    if self.software_name != 'PIPP':
                        for atomSelectionSet in self.atomSelectionSets:
                            for onebond_idx, (dim1, dim2) in enumerate(DIM_TRANSFER_PAT_3D):
                                _atom1, _atom2 = atomSelectionSet[dim1][0], atomSelectionSet[dim2][0]
                                if _atom1['chain_id'] != _atom2['chain_id']\
                                   or _atom1['seq_id'] != _atom2['seq_id']\
                                   or _atom1['atom_id'][0] == _atom2['atom_id'][0]:
                                    continue
                                if self.ccU.hasBond(_atom1['comp_id'], _atom1['atom_id'], _atom2['atom_id']):
                                    history.append(onebond_idx)
                    set_id = None if len(self.atomSelectionSets) < 2 else 0
                    if set_id is None:
                        for atomSelectionSet, asIsSet in zip(self.atomSelectionSets, self.asIsSets):
                            uniqAtoms = []
                            for idx in range(self.num_of_dim):
                                atom_sel = atomSelectionSet[idx]
                                common_atom = self.__extractCommonAtom(atom_sel)
                                if common_atom not in uniqAtoms:
                                    asis = asIsSet[idx]
                                    ambig_code = None
                                    if len(atom_sel) > 1:
                                        comp_id = common_atom['comp_id']
                                        rep_atom_id = atom_sel[0]['atom_id']
                                        ambig_code = self.csStat.getMaxAmbigCodeWoSetId(comp_id, rep_atom_id)
                                        if ambig_code == 0:
                                            ambig_code = None
                                        elif ambig_code == 2:
                                            hvy_grp_atoms, _ = self.nefT.get_group(comp_id, rep_atom_id)
                                            hvy_gem_atoms, pro_gem_atoms = self.nefT.get_geminal_group(comp_id, rep_atom_id)
                                            if None not in (hvy_grp_atoms, hvy_gem_atoms):
                                                gem_atom_ids = pro_gem_atoms if rep_atom_id[0] in PROTON_BEGIN_CODE\
                                                    else hvy_gem_atoms
                                                if not any(a['atom_id'] in gem_atom_ids for a in atom_sel):
                                                    ambig_code = 1
                                    if 'ambig_code' in common_atom:
                                        ambig_code = common_atom['ambig_code']
                                    row = getPkChemShiftRow(self.cur_subtype, sf['id'], sf['list_id'], self.__entryId,
                                                            dstFunc, set_id, idx + 1,
                                                            self.__authToStarSeq, self.__authToOrigSeq, self.__offsetHolder,
                                                            common_atom, asis, ambig_code)
                                    if row is None:
                                        continue
                                    sf['alt_loops'][3].add_data(row)
                                    uniqAtoms.append(common_atom)
                    else:
                        atom_set1, atom_set2, atom_set3 = [], [], []
                        for atomSelectionSet in self.atomSelectionSets:
                            valid = True
                            if len(history) > 0:
                                onebond_idx = collections.Counter(history).most_common()[0][0]
                                _atom1, _atom2 =\
                                    atomSelectionSet[DIM_TRANSFER_PAT_3D[onebond_idx][0]][0], \
                                    atomSelectionSet[DIM_TRANSFER_PAT_3D[onebond_idx][1]][0]
                                if _atom1['chain_id'] != _atom2['chain_id']\
                                   or _atom1['seq_id'] != _atom2['seq_id']\
                                   or _atom1['atom_id'][0] == _atom2['atom_id'][0]\
                                   or not self.ccU.hasBond(_atom1['comp_id'], _atom1['atom_id'], _atom2['atom_id']):
                                    valid = False
                            if not valid:
                                continue
                            atom_set1.extend(atomSelectionSet[0])
                            atom_set2.extend(atomSelectionSet[1])
                            atom_set3.extend(atomSelectionSet[2])
                        common_atoms = []
                        common_atoms.append(self.__extractCommonAtom([dict(s) for s in dict.fromkeys(frozenset(a.items())
                                                                                                     for a in atom_set1
                                                                                                     if isinstance(a, dict))]
                                                                     if len(atom_set1) > 1 else atom_set1))
                        common_atoms.append(self.__extractCommonAtom([dict(s) for s in dict.fromkeys(frozenset(a.items())
                                                                                                     for a in atom_set2
                                                                                                     if isinstance(a, dict))]
                                                                     if len(atom_set2) > 1 else atom_set2))
                        common_atoms.append(self.__extractCommonAtom([dict(s) for s in dict.fromkeys(frozenset(a.items())
                                                                                                     for a in atom_set3
                                                                                                     if isinstance(a, dict))]
                                                                     if len(atom_set3) > 1 else atom_set3))
                        set_id = 1
                        for atomSelectionSet, asIsSet in itertools.zip_longest(self.atomSelectionSets, self.asIsSets):
                            valid = True
                            if len(history) > 0:
                                onebond_idx = collections.Counter(history).most_common()[0][0]
                                _atom1, _atom2 =\
                                    atomSelectionSet[DIM_TRANSFER_PAT_3D[onebond_idx][0]][0], \
                                    atomSelectionSet[DIM_TRANSFER_PAT_3D[onebond_idx][1]][0]
                                if _atom1['chain_id'] != _atom2['chain_id']\
                                   or _atom1['seq_id'] != _atom2['seq_id']\
                                   or _atom1['atom_id'][0] == _atom2['atom_id'][0]\
                                   or not self.ccU.hasBond(_atom1['comp_id'], _atom1['atom_id'], _atom2['atom_id']):
                                    valid = False
                            if not valid:
                                continue
                            uniqAtoms = []
                            for idx in range(self.num_of_dim):
                                atom_sel = atomSelectionSet[idx]
                                common_atom = self.__extractCommonAtom(atom_sel)
                                if common_atom not in uniqAtoms:
                                    asis = asIsSet[idx]
                                    ambig_code = None
                                    if len(atom_sel) > 1:
                                        comp_id = common_atom['comp_id']
                                        rep_atom_id = atom_sel[0]['atom_id']
                                        ambig_code = self.csStat.getMaxAmbigCodeWoSetId(comp_id, rep_atom_id)
                                        if ambig_code == 0:
                                            ambig_code = None
                                        elif ambig_code == 2:
                                            hvy_grp_atoms, _ = self.nefT.get_group(comp_id, rep_atom_id)
                                            hvy_gem_atoms, pro_gem_atoms = self.nefT.get_geminal_group(comp_id, rep_atom_id)
                                            if None not in (hvy_grp_atoms, hvy_gem_atoms):
                                                gem_atom_ids = pro_gem_atoms if rep_atom_id[0] in PROTON_BEGIN_CODE\
                                                    else hvy_gem_atoms
                                                if not any(a['atom_id'] in gem_atom_ids for a in atom_sel):
                                                    ambig_code = 1
                                    if 'ambig_code' in common_atoms[idx]:
                                        ambig_code = common_atoms[idx]['ambig_code']
                                    row = getPkChemShiftRow(self.cur_subtype, sf['id'], sf['list_id'], self.__entryId,
                                                            dstFunc, set_id, idx + 1,
                                                            self.__authToStarSeq, self.__authToOrigSeq, self.__offsetHolder,
                                                            common_atom, asis, ambig_code)
                                    if row is None:
                                        continue
                                    sf['alt_loops'][3].add_data(row)
                                    uniqAtoms.append(common_atom)
                            set_id += 1

    def addAssignedPkRow4D(self, index: int, dstFunc: dict, has_assignments: bool, has_multiple_assignments: bool,
                           asis1: Optional[bool], asis2: Optional[bool], asis3: Optional[bool], asis4: Optional[bool],
                           debug_label: Optional[str], details: Optional[str]) -> None:
        """ Add assigned peak list row (4D).
        """

        if self.__debug:
            if not has_assignments:
                print(f"subtype={self.cur_subtype} id={self.peaks4D} (index={index}) "
                      f"{debug_label}None None None None {dstFunc}")
            for idx, atomSelectionSet in enumerate(self.atomSelectionSets, start=1):
                if has_multiple_assignments:
                    print(f"subtype={self.cur_subtype} id={self.peaks4D} (index={index}) combination_id={idx} "
                          f"{debug_label}{atomSelectionSet[0]} "
                          f"{atomSelectionSet[1]} "
                          f"{atomSelectionSet[2]} "
                          f"{atomSelectionSet[3]} {dstFunc}")
                else:
                    print(f"subtype={self.cur_subtype} id={self.peaks4D} (index={index}) "
                          f"{debug_label}{atomSelectionSet[0]} "
                          f"{atomSelectionSet[1]} "
                          f"{atomSelectionSet[2]} "
                          f"{atomSelectionSet[3]} {dstFunc}")

        if self.__createSfDict:
            sf = self.getSf()

            if sf is not None:
                sf['id'] = index

                if has_assignments and has_multiple_assignments:
                    history = self.onebond_idx_history[self.num_of_dim][self.cur_list_id]
                    if self.software_name != 'PIPP':
                        for atomSelectionSet in self.atomSelectionSets:
                            for onebond_idx, ((dim1, dim2), (dim3, dim4)) in enumerate(DIM_TRANSFER_PAT_4D):
                                _atom1, _atom2, _atom3, _atom4 =\
                                    atomSelectionSet[dim1][0], atomSelectionSet[dim2][0], \
                                    atomSelectionSet[dim3][0], atomSelectionSet[dim4][0]
                                if _atom1['chain_id'] != _atom2['chain_id']\
                                   or _atom1['seq_id'] != _atom2['seq_id']\
                                   or _atom1['atom_id'][0] == _atom2['atom_id'][0]\
                                   or _atom3['chain_id'] != _atom4['chain_id']\
                                   or _atom3['seq_id'] != _atom4['seq_id']\
                                   or _atom3['atom_id'][0] == _atom4['atom_id'][0]:
                                    continue
                                if self.ccU.hasBond(_atom1['comp_id'], _atom1['atom_id'], _atom2['atom_id'])\
                                   and self.ccU.hasBond(_atom3['comp_id'], _atom3['atom_id'], _atom4['atom_id']):
                                    history.append(onebond_idx)
                    atom_set1, atom_set2, atom_set3, atom_set4 = [], [], [], []
                    for atomSelectionSet in self.atomSelectionSets:
                        valid = True
                        if len(history) > 0:
                            onebond_idx = collections.Counter(history).most_common()[0][0]
                            _atom1, _atom2, _atom3, _atom4 =\
                                atomSelectionSet[DIM_TRANSFER_PAT_4D[onebond_idx][0][0]][0], \
                                atomSelectionSet[DIM_TRANSFER_PAT_4D[onebond_idx][0][1]][0], \
                                atomSelectionSet[DIM_TRANSFER_PAT_4D[onebond_idx][1][0]][0], \
                                atomSelectionSet[DIM_TRANSFER_PAT_4D[onebond_idx][1][1]][0]
                            if _atom1['chain_id'] != _atom2['chain_id']\
                               or _atom1['seq_id'] != _atom2['seq_id']\
                               or _atom1['atom_id'][0] == _atom2['atom_id'][0]\
                               or not self.ccU.hasBond(_atom1['comp_id'], _atom1['atom_id'], _atom2['atom_id'])\
                               or _atom3['chain_id'] != _atom4['chain_id']\
                               or _atom3['seq_id'] != _atom4['seq_id']\
                               or _atom3['atom_id'][0] == _atom4['atom_id'][0]\
                               or not self.ccU.hasBond(_atom3['comp_id'], _atom3['atom_id'], _atom4['atom_id']):
                                valid = False
                        if not valid:
                            continue
                        atom_set1.extend(atomSelectionSet[0])
                        atom_set2.extend(atomSelectionSet[1])
                        atom_set3.extend(atomSelectionSet[2])
                        atom_set4.extend(atomSelectionSet[3])
                    common_atom1 = self.__extractCommonAtom([dict(s) for s in dict.fromkeys(frozenset(a.items())
                                                                                            for a in atom_set1
                                                                                            if isinstance(a, dict))]
                                                            if len(atom_set1) > 1 else atom_set1)
                    common_atom2 = self.__extractCommonAtom([dict(s) for s in dict.fromkeys(frozenset(a.items())
                                                                                            for a in atom_set2
                                                                                            if isinstance(a, dict))]
                                                            if len(atom_set2) > 1 else atom_set2)
                    common_atom3 = self.__extractCommonAtom([dict(s) for s in dict.fromkeys(frozenset(a.items())
                                                                                            for a in atom_set3
                                                                                            if isinstance(a, dict))]
                                                            if len(atom_set3) > 1 else atom_set3)
                    common_atom4 = self.__extractCommonAtom([dict(s) for s in dict.fromkeys(frozenset(a.items())
                                                                                            for a in atom_set4
                                                                                            if isinstance(a, dict))]
                                                            if len(atom_set4) > 1 else atom_set4)
                    for atomSelectionSet, asIsSet in zip(self.atomSelectionSets, self.asIsSets):
                        ambig_code1 = ambig_code2 = ambig_code3 = ambig_code4 = None
                        atom1 = self.__extractCommonAtom(atomSelectionSet[0])
                        atom2 = self.__extractCommonAtom(atomSelectionSet[1])
                        atom3 = self.__extractCommonAtom(atomSelectionSet[2])
                        atom4 = self.__extractCommonAtom(atomSelectionSet[3])
                        asis1, asis2, asis3, asis4 = asIsSet
                        if len(atomSelectionSet[0]) > 1:
                            ambig_code1 = self.csStat.getMaxAmbigCodeWoSetId(atom1['comp_id'], atom1['atom_id'])
                            if ambig_code1 == 0:
                                ambig_code1 = None
                        if 'ambig_code' in common_atom1:
                            ambig_code1 = common_atom1['ambig_code']
                        if len(atomSelectionSet[1]) > 1:
                            ambig_code2 = self.csStat.getMaxAmbigCodeWoSetId(atom2['comp_id'], atom2['atom_id'])
                            if ambig_code2 == 0:
                                ambig_code2 = None
                        if 'ambig_code' in common_atom2:
                            ambig_code2 = common_atom2['ambig_code']
                        if len(atomSelectionSet[2]) > 1:
                            ambig_code3 = self.csStat.getMaxAmbigCodeWoSetId(atom3['comp_id'], atom3['atom_id'])
                            if ambig_code3 == 0:
                                ambig_code3 = None
                        if 'ambig_code' in common_atom3:
                            ambig_code3 = common_atom3['ambig_code']
                        if len(atomSelectionSet[3]) > 1:
                            ambig_code4 = self.csStat.getMaxAmbigCodeWoSetId(atom4['comp_id'], atom4['atom_id'])
                            if ambig_code4 == 0:
                                ambig_code4 = None
                        if 'ambig_code' in common_atom4:
                            ambig_code4 = common_atom4['ambig_code']

                        sf['row_index_id'] += 1

                        row = getPkRow(self.cur_subtype, sf['id'], sf['row_index_id'],
                                       sf['list_id'], self.__entryId, dstFunc,
                                       self.__authToStarSeq, self.__authToOrigSeq, self.__offsetHolder,
                                       atom1, atom2, atom3, atom4,
                                       asis1=asis1, asis2=asis2, asis3=asis3, asis4=asis4,
                                       ambig_code1=ambig_code1, ambig_code2=ambig_code2,
                                       ambig_code3=ambig_code3, ambig_code4=ambig_code4,
                                       details=details)
                        sf['loop'].add_data(row)

                else:

                    sf['row_index_id'] += 1

                    ambig_code1 = ambig_code2 = ambig_code3 = ambig_code4 = None
                    if has_assignments:
                        atom1 = self.__extractCommonAtom(self.atomSelectionSet[0])
                        atom2 = self.__extractCommonAtom(self.atomSelectionSet[1])
                        atom3 = self.__extractCommonAtom(self.atomSelectionSet[2])
                        atom4 = self.__extractCommonAtom(self.atomSelectionSet[3])
                        if len(self.atomSelectionSet[0]) > 1:
                            ambig_code1 = self.csStat.getMaxAmbigCodeWoSetId(atom1['comp_id'], atom1['atom_id'])
                            if ambig_code1 == 0:
                                ambig_code1 = None
                        if 'ambig_code' in atom1:
                            ambig_code1 = atom1['ambig_code']
                        if len(self.atomSelectionSet[1]) > 1:
                            ambig_code2 = self.csStat.getMaxAmbigCodeWoSetId(atom2['comp_id'], atom2['atom_id'])
                            if ambig_code2 == 0:
                                ambig_code2 = None
                        if 'ambig_code' in atom2:
                            ambig_code2 = atom2['ambig_code']
                        if len(self.atomSelectionSet[2]) > 1:
                            ambig_code3 = self.csStat.getMaxAmbigCodeWoSetId(atom3['comp_id'], atom3['atom_id'])
                            if ambig_code3 == 0:
                                ambig_code3 = None
                        if 'ambig_code' in atom3:
                            ambig_code3 = atom3['ambig_code']
                        if len(self.atomSelectionSet[3]) > 1:
                            ambig_code4 = self.csStat.getMaxAmbigCodeWoSetId(atom4['comp_id'], atom4['atom_id'])
                            if ambig_code4 == 0:
                                ambig_code4 = None
                        if 'ambig_code' in atom4:
                            ambig_code4 = atom4['ambig_code']
                    else:
                        atom1 = atom2 = atom3 = atom4 = None

                    row = getPkRow(self.cur_subtype, sf['id'], sf['index_id'],
                                   sf['list_id'], self.__entryId, dstFunc,
                                   self.__authToStarSeq, self.__authToOrigSeq, self.__offsetHolder,
                                   atom1, atom2, atom3, atom4,
                                   asis1=asis1, asis2=asis2, asis3=asis3, asis4=asis4,
                                   ambig_code1=ambig_code1, ambig_code2=ambig_code2,
                                   ambig_code3=ambig_code3, ambig_code4=ambig_code4,
                                   details=details)
                    sf['loop'].add_data(row)

                if not has_assignments and details is not None:
                    self.f.append(f"[Conflicted peak assignment] {self.getCurrentSpectralPeak(n=index)}"
                                  f"Peak assignments {details!r} could not map to the coordinates.")

                sf['index_id'] += 1

                row = getAltPkRow(self.cur_subtype, sf['index_id'], sf['id'], sf['list_id'], self.__entryId, dstFunc)
                if row is not None:
                    sf['alt_loops'][0].add_data(row)

                row = getPkGenCharRow(self.cur_subtype, sf['id'], sf['list_id'], self.__entryId, dstFunc, 'volume')
                if row is not None:
                    sf['alt_loops'][1].add_data(row)
                row = getPkGenCharRow(self.cur_subtype, sf['id'], sf['list_id'], self.__entryId, dstFunc, 'height')
                if row is not None:
                    sf['alt_loops'][1].add_data(row)

                for idx in range(self.num_of_dim):
                    row = getPkCharRow(self.cur_subtype, sf['id'], sf['list_id'], self.__entryId, dstFunc, idx + 1)
                    sf['alt_loops'][2].add_data(row)
                if has_assignments:
                    history = self.onebond_idx_history[self.num_of_dim][self.cur_list_id]
                    if self.software_name != 'PIPP':
                        for atomSelectionSet in self.atomSelectionSets:
                            for onebond_idx, ((dim1, dim2), (dim3, dim4)) in enumerate(DIM_TRANSFER_PAT_4D):
                                _atom1, _atom2, _atom3, _atom4 =\
                                    atomSelectionSet[dim1][0], atomSelectionSet[dim2][0], \
                                    atomSelectionSet[dim3][0], atomSelectionSet[dim4][0]
                                if _atom1['chain_id'] != _atom2['chain_id']\
                                   or _atom1['seq_id'] != _atom2['seq_id']\
                                   or _atom1['atom_id'][0] == _atom2['atom_id'][0]\
                                   or _atom3['chain_id'] != _atom4['chain_id']\
                                   or _atom3['seq_id'] != _atom4['seq_id']\
                                   or _atom3['atom_id'][0] == _atom4['atom_id'][0]:
                                    continue
                                if self.ccU.hasBond(_atom1['comp_id'], _atom1['atom_id'], _atom2['atom_id'])\
                                   and self.ccU.hasBond(_atom3['comp_id'], _atom3['atom_id'], _atom4['atom_id']):
                                    history.append(onebond_idx)
                    set_id = None if len(self.atomSelectionSets) < 2 else 0
                    if set_id is None:
                        for atomSelectionSet, asIsSet in zip(self.atomSelectionSets, self.asIsSets):
                            uniqAtoms = []
                            for idx in range(self.num_of_dim):
                                atom_sel = atomSelectionSet[idx]
                                common_atom = self.__extractCommonAtom(atom_sel)
                                if common_atom not in uniqAtoms:
                                    asis = asIsSet[idx]
                                    ambig_code = None
                                    if len(atom_sel) > 1:
                                        comp_id = common_atom['comp_id']
                                        rep_atom_id = atom_sel[0]['atom_id']
                                        ambig_code = self.csStat.getMaxAmbigCodeWoSetId(comp_id, rep_atom_id)
                                        if ambig_code == 0:
                                            ambig_code = None
                                        elif ambig_code == 2:
                                            hvy_grp_atoms, _ = self.nefT.get_group(comp_id, rep_atom_id)
                                            hvy_gem_atoms, pro_gem_atoms = self.nefT.get_geminal_group(comp_id, rep_atom_id)
                                            if None not in (hvy_grp_atoms, hvy_gem_atoms):
                                                gem_atom_ids = pro_gem_atoms if rep_atom_id[0] in PROTON_BEGIN_CODE\
                                                    else hvy_gem_atoms
                                                if not any(a['atom_id'] in gem_atom_ids for a in atom_sel):
                                                    ambig_code = 1
                                    if 'ambig_code' in common_atom:
                                        ambig_code = common_atom['ambig_code']
                                    row = getPkChemShiftRow(self.cur_subtype, sf['id'], sf['list_id'], self.__entryId,
                                                            dstFunc, set_id, idx + 1,
                                                            self.__authToStarSeq, self.__authToOrigSeq, self.__offsetHolder,
                                                            common_atom, asis, ambig_code)
                                    if row is None:
                                        continue
                                    sf['alt_loops'][3].add_data(row)
                                    uniqAtoms.append(common_atom)
                    else:
                        atom_set1, atom_set2, atom_set3, atom_set4 = [], [], [], []
                        for atomSelectionSet in self.atomSelectionSets:
                            valid = True
                            if len(history) > 0:
                                onebond_idx = collections.Counter(history).most_common()[0][0]
                                _atom1, _atom2, _atom3, _atom4 =\
                                    atomSelectionSet[DIM_TRANSFER_PAT_4D[onebond_idx][0][0]][0], \
                                    atomSelectionSet[DIM_TRANSFER_PAT_4D[onebond_idx][0][1]][0], \
                                    atomSelectionSet[DIM_TRANSFER_PAT_4D[onebond_idx][1][0]][0], \
                                    atomSelectionSet[DIM_TRANSFER_PAT_4D[onebond_idx][1][1]][0]
                                if _atom1['chain_id'] != _atom2['chain_id']\
                                   or _atom1['seq_id'] != _atom2['seq_id']\
                                   or _atom1['atom_id'][0] == _atom2['atom_id'][0]\
                                   or not self.ccU.hasBond(_atom1['comp_id'], _atom1['atom_id'], _atom2['atom_id'])\
                                   or _atom3['chain_id'] != _atom4['chain_id']\
                                   or _atom3['seq_id'] != _atom4['seq_id']\
                                   or _atom3['atom_id'][0] == _atom4['atom_id'][0]\
                                   or not self.ccU.hasBond(_atom3['comp_id'], _atom3['atom_id'], _atom4['atom_id']):
                                    valid = False
                            if not valid:
                                continue
                            atom_set1.extend(atomSelectionSet[0])
                            atom_set2.extend(atomSelectionSet[1])
                            atom_set3.extend(atomSelectionSet[2])
                            atom_set4.extend(atomSelectionSet[3])
                        common_atoms = []
                        common_atoms.append(self.__extractCommonAtom([dict(s) for s in dict.fromkeys(frozenset(a.items())
                                                                                                     for a in atom_set1
                                                                                                     if isinstance(a, dict))]
                                                                     if len(atom_set1) > 1 else atom_set1))
                        common_atoms.append(self.__extractCommonAtom([dict(s) for s in dict.fromkeys(frozenset(a.items())
                                                                                                     for a in atom_set2
                                                                                                     if isinstance(a, dict))]
                                                                     if len(atom_set2) > 1 else atom_set2))
                        common_atoms.append(self.__extractCommonAtom([dict(s) for s in dict.fromkeys(frozenset(a.items())
                                                                                                     for a in atom_set3
                                                                                                     if isinstance(a, dict))]
                                                                     if len(atom_set3) > 1 else atom_set3))
                        common_atoms.append(self.__extractCommonAtom([dict(s) for s in dict.fromkeys(frozenset(a.items())
                                                                                                     for a in atom_set4
                                                                                                     if isinstance(a, dict))]
                                                                     if len(atom_set4) > 1 else atom_set4))
                        set_id = 1
                        for atomSelectionSet, asIsSet in itertools.zip_longest(self.atomSelectionSets, self.asIsSets):
                            valid = True
                            if len(history) > 0:
                                onebond_idx = collections.Counter(history).most_common()[0][0]
                                _atom1, _atom2, _atom3, _atom4 =\
                                    atomSelectionSet[DIM_TRANSFER_PAT_4D[onebond_idx][0][0]][0], \
                                    atomSelectionSet[DIM_TRANSFER_PAT_4D[onebond_idx][0][1]][0], \
                                    atomSelectionSet[DIM_TRANSFER_PAT_4D[onebond_idx][1][0]][0], \
                                    atomSelectionSet[DIM_TRANSFER_PAT_4D[onebond_idx][1][1]][0]
                                if _atom1['chain_id'] != _atom2['chain_id']\
                                   or _atom1['seq_id'] != _atom2['seq_id']\
                                   or _atom1['atom_id'][0] == _atom2['atom_id'][0]\
                                   or not self.ccU.hasBond(_atom1['comp_id'], _atom1['atom_id'], _atom2['atom_id'])\
                                   or _atom3['chain_id'] != _atom4['chain_id']\
                                   or _atom3['seq_id'] != _atom4['seq_id']\
                                   or _atom3['atom_id'][0] == _atom4['atom_id'][0]\
                                   or not self.ccU.hasBond(_atom3['comp_id'], _atom3['atom_id'], _atom4['atom_id']):
                                    valid = False
                            if not valid:
                                continue
                            uniqAtoms = []
                            for idx in range(self.num_of_dim):
                                atom_sel = atomSelectionSet[idx]
                                common_atom = self.__extractCommonAtom(atom_sel)
                                if common_atom not in uniqAtoms:
                                    asis = asIsSet[idx]
                                    ambig_code = None
                                    if len(atom_sel) > 1:
                                        comp_id = common_atom['comp_id']
                                        rep_atom_id = atom_sel[0]['atom_id']
                                        ambig_code = self.csStat.getMaxAmbigCodeWoSetId(comp_id, rep_atom_id)
                                        if ambig_code == 0:
                                            ambig_code = None
                                        elif ambig_code == 2:
                                            hvy_grp_atoms, _ = self.nefT.get_group(comp_id, rep_atom_id)
                                            hvy_gem_atoms, pro_gem_atoms = self.nefT.get_geminal_group(comp_id, rep_atom_id)
                                            if None not in (hvy_grp_atoms, hvy_gem_atoms):
                                                gem_atom_ids = pro_gem_atoms if rep_atom_id[0] in PROTON_BEGIN_CODE\
                                                    else hvy_gem_atoms
                                                if not any(a['atom_id'] in gem_atom_ids for a in atom_sel):
                                                    ambig_code = 1
                                    if 'ambig_code' in common_atoms[idx]:
                                        ambig_code = common_atoms[idx]['ambig_code']
                                    row = getPkChemShiftRow(self.cur_subtype, sf['id'], sf['list_id'], self.__entryId,
                                                            dstFunc, set_id, idx + 1,
                                                            self.__authToStarSeq, self.__authToOrigSeq, self.__offsetHolder,
                                                            common_atom, asis, ambig_code)
                                    if row is None:
                                        continue
                                    sf['alt_loops'][3].add_data(row)
                                    uniqAtoms.append(common_atom)
                            set_id += 1

    def extractPeakAssignment(self, numOfDim: int, string: str, src_index: int,
                              with_segid: Optional[str] = None, with_compid: Optional[str] = None,
                              hint: Optional[List[dict]] = None, dim_id_hint: Optional[int] = None
                              ) -> Optional[List[dict]]:
        """ Extract peak assignment from a given string.
        """

        if numOfDim not in (1, 2, 3, 4) or string is None:
            return None

        _str_ = ASSIGNMENT_SEPARATOR_PAT.sub(' ', string).split()
        _str = ASSIGNMENT_SEPARATOR_PAT.sub(' ', string.upper()).split()
        lenStr = len(_str)

        segIdLike, resIdLike, resNameLike, atomNameLike, _atomNameLike, __atomNameLike, ___atomNameLike, atomNameLike_ =\
            [False] * lenStr, [False] * lenStr, [False] * lenStr, [False] * lenStr, \
            [False] * lenStr, [False] * lenStr, [False] * lenStr, [False] * lenStr

        segIdSpan, resIdSpan, resNameSpan, atomNameSpan, _atomNameSpan, __atomNameSpan, ___atomNameSpan, siblingAtomName =\
            [None] * lenStr, [None] * lenStr, [None] * lenStr, [None] * lenStr, \
            [None] * lenStr, [None] * lenStr, [None] * lenStr, [None] * lenStr

        if not self.__hasCoord:
            if self.compIdSet is None:
                self.compIdSet = self.altCompIdSet = set(STD_MON_DICT.keys())

        oneLetterCodeSet = []
        extMonDict3 = {}
        if self.polyPeptide and not self.polyDeoxyribonucleotide and not self.polyRibonucleotide:
            oneLetterCodeSet = [getOneLetterCode(compId) for compId in self.compIdSet if len(compId) == 3]
            extMonDict3 = {compId: getOneLetterCode(compId) for compId in self.compIdSet if len(compId) == 3}
        elif not self.polyPeptide and self.polyDeoxyribonucleotide and not self.polyRibonucleotide:
            oneLetterCodeSet = [getOneLetterCode(compId) for compId in self.compIdSet if len(compId) == 2]
            extMonDict3 = {compId: getOneLetterCode(compId) for compId in self.compIdSet if len(compId) == 2}
        elif not self.polyPeptide and not self.polyDeoxyribonucleotide and self.polyRibonucleotide:
            oneLetterCodeSet = [getOneLetterCode(compId) for compId in self.compIdSet if len(compId) == 1]
            extMonDict3 = {compId: getOneLetterCode(compId) for compId in self.compIdSet if len(compId) == 1}

        if self.hasNonPolySeq:
            for np in self.nonPoly:
                if np['comp_id'][0][0].isalpha() and np['comp_id'][0][0] not in oneLetterCodeSet:
                    oneLetterCodeSet.append(np['comp_id'][0][0])
                    extMonDict3[np['comp_id'][0]] = np['comp_id'][0][0]
                if np['auth_comp_id'][0][0].isalpha() and np['auth_comp_id'][0][0] not in oneLetterCodeSet:
                    oneLetterCodeSet.append(np['auth_comp_id'][0][0])
                    extMonDict3[np['comp_id'][0]] = np['auth_comp_id'][0][0]
                if np['alt_comp_id'][0][0].isalpha() and np['alt_comp_id'][0][0] not in oneLetterCodeSet:
                    oneLetterCodeSet.append(np['alt_comp_id'][0][0])
                    extMonDict3[np['comp_id'][0]] = np['alt_comp_id'][0][0]

        hasOneLetterCodeSet = len(oneLetterCodeSet) > 0
        useOneLetterCodeSet = forceOneLetterCodeSet = \
            self.polyRibonucleotide and not self.polyPeptide and not self.polyDeoxyribonucleotide
        ligCompId = ligAtomId = None
        _ligSeqId = _ligCompId = _ligAtomId = None

        for idx, term in enumerate(_str):
            for segId in self.authAsymIdSet:
                if term.startswith(segId) and (with_segid is None or term.startswith(with_segid)):
                    segIdLike[idx] = True
                    segIdSpan[idx] = (0, len(segId))
                    break

            resIdTest = ASSIGNMENT_RESID_PAT.search(term)
            if resIdTest:
                if term[0] == 'D' and len(term) == 3 and term[-1] in ('5', '3')\
                   and translateToStdResName(term, ccU=self.ccU) in self.compIdSet:
                    pass
                else:
                    resIdLike[idx] = True
                    resIdSpan[idx] = resIdTest.span()

            minIndex = len(term)

            for compId in self.compIdSet:
                if compId in term:
                    resNameLike[idx] = True
                    index = term.index(compId)
                    if index < minIndex:
                        resNameSpan[idx] = (index, index + len(compId))
                        minIndex = index

            if not resNameLike[idx]:
                for compId in self.altCompIdSet:
                    if compId in term:
                        resNameLike[idx] = True
                        index = term.index(compId)
                        if index < minIndex:
                            resNameSpan[idx] = (index, index + len(compId))
                            minIndex = index

            if not resNameLike[idx] and self.cyanaCompIdSet is not None:
                for compId in self.cyanaCompIdSet:
                    if compId in term:
                        resNameLike[idx] = True
                        index = term.index(compId)
                        if index < minIndex:
                            resNameSpan[idx] = (index, index + len(compId))
                            minIndex = index

            if hasOneLetterCodeSet and not useOneLetterCodeSet and resNameLike[idx]\
               and len(term[resNameSpan[idx][0]:resNameSpan[idx][1]]) > 1:
                hasOneLetterCodeSet = False

            if with_compid is not None and len(with_compid) > 1:
                hasOneLetterCodeSet = False

            if not resNameLike[idx] and hasOneLetterCodeSet:
                if not any(compId in term for compId in STD_MON_DICT if len(compId) == 3):
                    for compId in oneLetterCodeSet:
                        if compId in term:
                            resNameLike[idx] = True
                            index = term.index(compId)
                            if index < minIndex:
                                resNameSpan[idx] = (index, index + len(compId))
                                minIndex = index
                elif resIdLike[idx] and self.hasPolySeq:
                    _compId = next(compId for compId in STD_MON_DICT if len(compId) == 3 and compId in term)
                    resId = int(term[resIdSpan[idx][0]:resIdSpan[idx][1]])
                    for ps in self.polySeq:
                        _, _, compId = self.getRealChainSeqId(ps, resId, None)
                        if len(compId) == 3 and compId in STD_MON_DICT and _compId[0:2] == compId[0:2]:
                            resNameLike[idx] = True
                            index = term.index(_compId)
                            resNameSpan[idx] = (index, index + len(compId))
                            term = _str[idx] = term.replace(_compId, compId)
                            break

            if resNameLike[idx]:
                compId = term[resNameSpan[idx][0]:resNameSpan[idx][1]]
                if compId in self.compIdSet and compId not in STD_MON_DICT:
                    ligCompId = compId

                if len(compId) == 1:
                    _ligCompId = next((k for k, v in extMonDict3.items() if v == compId and k not in STD_MON_DICT), None)

                    if resIdLike[idx] and self.hasPolySeq:
                        resId = int(term[resIdSpan[idx][0]:resIdSpan[idx][1]])
                        for ps in self.polySeq:
                            _, _, _compId = self.getRealChainSeqId(ps, resId, None)
                            if _compId is not None and _compId not in STD_MON_DICT:
                                resNameLike[idx] = False
                                break
                    if resNameLike[idx]:
                        index = resNameSpan[idx][1]
                        if index < len(term):
                            if term[index].isdigit() or term[index] in ASSIGNMENT_HALF_SPIN_NUCLEUS\
                               or term[index] in PSE_PRO_BEGIN_CODE:
                                pass
                            else:
                                resNameLike[idx] = False

            if ligCompId is not None and ligCompId != term:
                _, _, details = self.nefT.get_valid_star_atom_in_xplor(ligCompId, term, leave_unmatched=True)
                if details is None or term[0] in ASSIGNMENT_HALF_SPIN_NUCLEUS:
                    atomNameLike[idx] = True
                    atomNameSpan[idx] = (0, len(term) + 1)
                    ligAtomId = term

            if _ligCompId is not None and resNameLike[idx] and len(term[resNameSpan[idx][1]:]) > 0:
                _, _, details = self.nefT.get_valid_star_atom_in_xplor(_ligCompId, term[resNameSpan[idx][1]:], leave_unmatched=True)
                if details is None or term[resNameSpan[idx][1]] in ASSIGNMENT_HALF_SPIN_NUCLEUS:
                    atomNameLike[idx] = True
                    atomNameSpan[idx] = (resNameSpan[idx][1], len(term) + 1)
                    _ligAtomId = term[resNameSpan[idx][1]:len(term) + 1]
                    for np in self.nonPoly:
                        if np['comp_id'][0] == _ligCompId:
                            _ligSeqId = np['auth_seq_id'][0]

            if resIdLike[idx] and resIdSpan[idx][1] + 1 <= len(term)\
               and _str_[idx][resIdSpan[idx][1]].islower() and _str[idx][resIdSpan[idx][1]].isupper():
                if resIdSpan[idx][1] + 1 < len(term)\
                   and any(_str_[idx][resIdSpan[idx][1] + 1].startswith(elem) for elem in ASSIGNMENT_HALF_SPIN_NUCLEUS):
                    term = _str[idx] = term[0:resIdSpan[idx][1]] + term[resIdSpan[idx][1] + 1:]
                elif resIdSpan[idx][1] + 1 == len(term):
                    term = _str[idx] = term[0:resIdSpan[idx][1]]

            for elem in reversed(ASSIGNMENT_HALF_SPIN_NUCLEUS) if 'NH' in term else ASSIGNMENT_HALF_SPIN_NUCLEUS:
                if len(elem) == 1 and ligAtomId is None and _ligAtomId is None:
                    if elem in term:

                        # handle ambiguous assigned peak '(14Trp/11Trp)Hh2'
                        # seen in 6r28/bmr34380/work/data/D_1292101294_nmr-peaks-upload_P1.dat.V1
                        if idx > 0 and not resNameLike[idx] and any(resIdLike[_idx] and resNameLike[_idx] and not atomNameLike[_idx]
                                                                    for _idx in range(idx)):
                            _index = term.index(elem)
                            _atomId = term[_index:len(term)]
                            for _idx in range(idx):
                                if resIdLike[_idx] and resNameLike[_idx] and not atomNameLike[_idx]:
                                    _compId = _str[_idx][resNameSpan[_idx][0]:resNameSpan[_idx][1]]
                                    if len(_compId) == 1 and hasOneLetterCodeSet:
                                        _compId = next(k for k, v in extMonDict3.items() if k in self.compIdSet and v == _compId)
                                        _, _, details =\
                                            self.nefT.get_valid_star_atom_in_xplor(_compId, _atomId, leave_unmatched=True)
                                        if details is None:
                                            atomNameLike[idx] = atomNameLike_[_idx] = useOneLetterCodeSet = True
                                            atomNameSpan[idx] = (_index, len(term))
                                            if siblingAtomName[_idx] is None:
                                                siblingAtomName[_idx] = []
                                            if _atomId not in siblingAtomName[_idx]:
                                                siblingAtomName[_idx].append(_atomId)
                                            break
                                    if (with_compid is not None and _atomId.startswith(with_compid)) or _atomId.startswith('MET'):
                                        continue
                                    _atomId = translateToStdAtomName(_atomId, _compId, ccU=self.ccU)
                                    _, _, details = self.nefT.get_valid_star_atom_in_xplor(_compId, _atomId, leave_unmatched=True)
                                    if details is None:
                                        atomNameLike[idx] = atomNameLike_[_idx] = True
                                        atomNameSpan[idx] = (_index, len(term))
                                        if siblingAtomName[_idx] is None:
                                            siblingAtomName[_idx] = []
                                        if _atomId not in siblingAtomName[_idx]:
                                            siblingAtomName[_idx].append(_atomId)
                                        break
                            if atomNameLike[idx]:
                                break

                        # prevent to split HH2 -> res_name:'HIS', atom_name:'H2'
                        if resNameLike[idx] and resNameSpan[idx][1] - resNameSpan[idx][0] == 1\
                           and resNameSpan[idx][1] == term.rindex(elem)\
                           and term[resNameSpan[idx][0]] in ASSIGNMENT_HALF_SPIN_NUCLEUS:
                            _index = term.index(elem)
                            _atomId = term[_index:len(term)]
                            for compId in self.compIdSet:
                                _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                                if details is None:
                                    atomNameLike[idx] = True
                                    atomNameSpan[idx] = (_index, len(term))
                                    resNameLike[idx] = False
                                    break
                            if atomNameLike[idx]:
                                break

                        # resolve concatenation of residue number and XPLOR-NIH atom nomenclature of proton,
                        # D1391HB -> res_id:139, res_name:'ASP', atom_name:'1HB'
                        # seen in 8e1d/bmr31038/work/data/D_1000267621_nmr-peaks-upload_P6.dat.V1
                        if self.hasPolySeq and resIdLike[idx] and resNameLike[idx]\
                           and resIdSpan[idx][1] == term.rindex(elem) and elem == 'H'\
                           and term[resIdSpan[idx][1] - 1] in ('1', '2', '3') and resIdSpan[idx][1] - resIdSpan[idx][0] > 3:
                            _resId = int(term[resIdSpan[idx][0]:resIdSpan[idx][1]])
                            _compId = term[resNameSpan[idx][0]:resNameSpan[idx][1]]
                            if len(_compId) == 1 and hasOneLetterCodeSet:
                                _compId = next((k for k, v in extMonDict3.items() if k in self.compIdSet and v == _compId), _compId)
                            valid = False
                            for ps in self.polySeq:
                                _, _, _compId_ = self.getRealChainSeqId(ps, _resId, None)
                                if _compId is not None and _compId == _compId_:
                                    valid = True
                                    break
                            if not valid:
                                _resId = int(term[resIdSpan[idx][0]:resIdSpan[idx][1] - 1])
                                for ps in self.polySeq:
                                    _, _, _compId_ = self.getRealChainSeqId(ps, _resId, None)
                                    if _compId is not None and _compId == _compId_:
                                        valid = True
                                        break
                                if valid:
                                    _index = term.rindex(elem)
                                    _atomId = term[_index:len(term)] + term[resIdSpan[idx][1] - 1]
                                    _, _, details = self.nefT.get_valid_star_atom_in_xplor(_compId, _atomId, leave_unmatched=True)
                                    if details is None:
                                        atomNameLike[idx] = True
                                        atomNameSpan[idx] = (_index - 1, len(term))
                                        resIdSpan[idx] = (resIdSpan[idx][0], resIdSpan[idx][1] - 1)
                                        break

                        index = term.rindex(elem)
                        atomId = term[index:len(term)]
                        if index - 1 >= 0 and term[index - 1] in ASSIGNMENT_HALF_SPIN_NUCLEUS:
                            if not resNameLike[idx]:
                                if hint is not None and 'comp_id' in hint[0]:
                                    compId = hint[0]['comp_id']
                                    _index = term.index(elem)
                                    _atomId = term[_index:len(term)]
                                    _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                                    if details is None:
                                        atomNameLike[idx] = True
                                        atomNameSpan[idx] = (_index, len(term))
                                    else:
                                        _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                                        if details is None:
                                            atomNameLike[idx] = True
                                            atomNameSpan[idx] = (index, len(term))
                                continue
                            if compId[-1] in ASSIGNMENT_HALF_SPIN_NUCLEUS and index == resNameSpan[idx][1]:
                                pass
                            elif len(compId) == 1 and hasOneLetterCodeSet:
                                compId = next(k for k, v in extMonDict3.items() if k in self.compIdSet and v == compId)
                                _atomId = term[index - 1:len(term)]
                                _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                                if details is None:
                                    index -= 1
                                    atomId = _atomId
                                else:
                                    continue
                            elif compId in self.compIdSet:
                                _atomId = term[index - 1:len(term)]
                                _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                                if details is None:
                                    index -= 1
                                    atomId = _atomId
                                else:
                                    continue
                            else:
                                continue

                        if atomId[0] in ('Q', 'M') and index + 1 < len(term) and term[index + 1].isdigit():
                            if resNameLike[idx] and resNameSpan[idx][0] == index:
                                continue
                            if self.csStat.peptideLike(compId):
                                ligand = False
                                if resIdLike[idx] and self.reasons is not None:
                                    resId = int(term[resIdSpan[idx][0]:resIdSpan[idx][1]])
                                    if 'non_poly_remap' in self.reasons\
                                       and compId in self.reasons['non_poly_remap']\
                                       and resId in self.reasons['non_poly_remap'][compId]:
                                        ligand = True
                                if not ligand:
                                    continue
                        if ((with_compid is not None and atomId.startswith(with_compid)) or atomId.startswith('MET'))\
                           and ((index + 3 < len(term) and term[index + 3].isdigit()
                                 or (index + 4 < len(term) and term[index + 4].isdigit()))):
                            continue
                        if resNameLike[idx] and len(compId) > 1 and compId[-1] == elem and index + 1 == resNameSpan[idx][1]:
                            continue
                        if resNameLike[idx]:
                            compId = term[resNameSpan[idx][0]:resNameSpan[idx][1]]
                            if len(compId) == 1 and hasOneLetterCodeSet and not forceOneLetterCodeSet:
                                compId = next(k for k, v in extMonDict3.items() if k in self.compIdSet and v == compId)
                                _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                                if details is None:
                                    if compId in ('DT', 'T') and atomId == 'C7' and self.cur_list_id != -1:
                                        cur_sp_dim = self.spectral_dim[self.num_of_dim][self.cur_list_id]
                                        if all('fixed' in cur_sp_dim[_dim_id] and cur_sp_dim[_dim_id]['atom_type'] == 'H'
                                               for _dim_id in range(1, self.num_of_dim + 1)):
                                            _str[idx] = _str[idx].replace('C7', 'H7')
                                            _string = ' '.join(_str)
                                            return self.extractPeakAssignment(numOfDim, _string, src_index,
                                                                              with_segid, with_compid, hint, dim_id_hint)
                                        if dim_id_hint is not None and 'freq_hint' in cur_sp_dim[dim_id_hint]\
                                           and cur_sp_dim[dim_id_hint]['freq_hint'][-1] < H_METHYL_CENTER_MAX:
                                            _str[idx] = _str[idx].replace('C7', 'H7')
                                            _string = ' '.join(_str)
                                            return self.extractPeakAssignment(numOfDim, _string, src_index,
                                                                              with_segid, with_compid, hint, dim_id_hint)
                                    atomNameLike[idx] = useOneLetterCodeSet = True
                                    atomNameSpan[idx] = (index, len(term))
                                    if resNameSpan[idx][0] == atomNameSpan[idx][0]:
                                        resNameLike[idx] = False
                                    break
                                if self.cur_list_id != -1:
                                    cur_sp_dim = self.spectral_dim[self.num_of_dim][self.cur_list_id]
                                    if all('fixed' in cur_sp_dim[_dim_id] and cur_sp_dim[_dim_id]['atom_type'] == 'H'
                                           for _dim_id in range(1, self.num_of_dim + 1)):
                                        if compId in ('DT', 'T') and atomId == 'CM':
                                            _str[idx] = _str[idx].replace('CM', 'H7')
                                            _string = ' '.join(_str)
                                            return self.extractPeakAssignment(numOfDim, _string, src_index,
                                                                              with_segid, with_compid, hint, dim_id_hint)
                                        if compId in ('DC', 'C') and atomId.startswith("NH"):
                                            if atomId == "NH''":
                                                _str[idx] = _str[idx].replace("NH''", 'H42')
                                            elif atomId == "NH'":
                                                _str[idx] = _str[idx].replace("NH'", 'H41')
                                            else:
                                                _str[idx] = _str[idx].replace("NH", 'H4')
                                            _string = ' '.join(_str)
                                            return self.extractPeakAssignment(numOfDim, _string, src_index,
                                                                              with_segid, with_compid, hint, dim_id_hint)
                                        if compId in ('DA', 'A') and atomId.startswith("NH"):
                                            if atomId == "NH''":
                                                _str[idx] = _str[idx].replace("NH''", 'H62')
                                            elif atomId == "NH'":
                                                _str[idx] = _str[idx].replace("NH'", 'H61')
                                            else:
                                                _str[idx] = _str[idx].replace("NH", 'H6')
                                            _string = ' '.join(_str)
                                            return self.extractPeakAssignment(numOfDim, _string, src_index,
                                                                              with_segid, with_compid, hint, dim_id_hint)
                                    if dim_id_hint is not None and 'freq_hint' in cur_sp_dim[dim_id_hint]\
                                       and cur_sp_dim[dim_id_hint]['freq_hint'][-1] < H_METHYL_CENTER_MAX:
                                        if compId in ('DT', 'T') and atomId == 'CM':
                                            _str[idx] = _str[idx].replace('CM', 'H7')
                                            _string = ' '.join(_str)
                                            return self.extractPeakAssignment(numOfDim, _string, src_index,
                                                                              with_segid, with_compid, hint, dim_id_hint)
                                if with_compid is not None and (atomId.startswith(with_compid)
                                                                or (atomId in with_compid and index < resNameSpan[idx][1])):
                                    continue
                                _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                                _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                                if details is None:
                                    atomNameLike[idx] = True
                                    atomNameSpan[idx] = (index, len(term))
                                    if resNameSpan[idx][0] == atomNameSpan[idx][0]:
                                        resNameLike[idx] = False
                                    break
                            # 'THR' should not be split 'T' and 'HR' (bmr26379)
                            if compId.endswith(atomId) and len(compId) > 2 and compId != atomId:
                                atomNameLike[idx] = False
                                break
                            # need to check the next token in case of 'MET' (bmr26379)
                            if compId == atomId == 'MET' and idx < len(_str) - 1:
                                _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _str[idx + 1], leave_unmatched=True)
                                if details is None:
                                    atomNameLike[idx] = False
                                    break

                        if not atomNameLike[idx] and hint is not None and 'comp_id' in hint[0] and self.cur_list_id != -1:
                            _compId = hint[0]['comp_id']
                            cur_sp_dim = self.spectral_dim[self.num_of_dim][self.cur_list_id]
                            if all('fixed' in cur_sp_dim[_dim_id] and cur_sp_dim[_dim_id]['atom_type'] == 'H'
                                   for _dim_id in range(1, self.num_of_dim + 1)):
                                if _compId in ('DT', 'T') and atomId == 'C7':
                                    _str[idx] = _str[idx].replace('C7', 'H7')
                                    _string = ' '.join(_str)
                                    return self.extractPeakAssignment(numOfDim, _string, src_index,
                                                                      with_segid, with_compid, hint, dim_id_hint)
                                if _compId in ('DT', 'T') and atomId == 'CM':
                                    _str[idx] = _str[idx].replace('CM', 'H7')
                                    _string = ' '.join(_str)
                                    return self.extractPeakAssignment(numOfDim, _string, src_index,
                                                                      with_segid, with_compid, hint, dim_id_hint)
                                if _compId in ('DC', 'C') and atomId.startswith("NH"):
                                    if atomId == "NH''":
                                        _str[idx] = _str[idx].replace("NH''", 'H42')
                                    elif atomId == "NH'":
                                        _str[idx] = _str[idx].replace("NH'", 'H41')
                                    else:
                                        _str[idx] = _str[idx].replace("NH", 'H4')
                                    _string = ' '.join(_str)
                                    return self.extractPeakAssignment(numOfDim, _string, src_index,
                                                                      with_segid, with_compid, hint, dim_id_hint)
                                if _compId in ('DA', 'A') and atomId.startswith("NH"):
                                    if atomId == "NH''":
                                        _str[idx] = _str[idx].replace("NH''", 'H62')
                                    elif atomId == "NH'":
                                        _str[idx] = _str[idx].replace("NH'", 'H61')
                                    else:
                                        _str[idx] = _str[idx].replace("NH", 'H6')
                                    _string = ' '.join(_str)
                                    return self.extractPeakAssignment(numOfDim, _string, src_index,
                                                                      with_segid, with_compid, hint, dim_id_hint)
                            if dim_id_hint is not None and 'freq_hint' in cur_sp_dim[dim_id_hint]\
                               and cur_sp_dim[dim_id_hint]['freq_hint'][-1] < H_METHYL_CENTER_MAX:
                                if _compId in ('DT', 'T') and atomId == 'C7':
                                    _str[idx] = _str[idx].replace('C7', 'H7')
                                    _string = ' '.join(_str)
                                    return self.extractPeakAssignment(numOfDim, _string, src_index,
                                                                      with_segid, with_compid, hint, dim_id_hint)
                                if _compId in ('DT', 'T') and atomId == 'CM':
                                    _str[idx] = _str[idx].replace('CM', 'H7')
                                    _string = ' '.join(_str)
                                    return self.extractPeakAssignment(numOfDim, _string, src_index,
                                                                      with_segid, with_compid, hint, dim_id_hint)

                        for compId in self.compIdSet:
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                            if details is None:
                                atomNameLike[idx] = True
                                atomNameSpan[idx] = (index, len(term))
                                break
                            if with_compid is not None and atomId.startswith(with_compid):
                                continue
                            _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                            if details is None:
                                atomNameLike[idx] = True
                                atomNameSpan[idx] = (index, len(term))
                                break
                        if atomNameLike[idx]:
                            break

            if atomNameLike[idx]:
                _term = term[0:atomNameSpan[idx][0]]
                for elem in reversed(ASSIGNMENT_HALF_SPIN_NUCLEUS) if 'NH' in _term else ASSIGNMENT_HALF_SPIN_NUCLEUS:
                    if len(elem) == 1 and ligAtomId is None and _ligAtomId is None:
                        if elem in _term:
                            index = _term.rindex(elem)
                            atomId = _term[index:len(_term)]
                            if index - 1 >= 0 and _term[index - 1] in ASSIGNMENT_HALF_SPIN_NUCLEUS:
                                if resNameLike[idx] and compId[-1] in ASSIGNMENT_HALF_SPIN_NUCLEUS and index == resNameSpan[idx][1]:
                                    pass
                                else:
                                    continue
                            if atomId[0] in ('Q', 'M') and index + 1 < len(_term) and _term[index + 1].isdigit():
                                continue
                            if ((with_compid is not None and atomId.startswith(with_compid)) or atomId.startswith('MET'))\
                               and ((index + 3 < len(_term) and _term[index + 3].isdigit()
                                     or (index + 4 < len(_term) and _term[index + 4].isdigit()))):
                                continue
                            if resNameLike[idx] and len(compId) > 1 and compId[-1] == elem and index + 1 == resNameSpan[idx][1]:
                                continue
                            if len(_term) == atomNameSpan[idx][0]:
                                continue
                            if resNameLike[idx]:
                                compId = term[resNameSpan[idx][0]:resNameSpan[idx][1]]
                                if len(compId) == 1 and hasOneLetterCodeSet:
                                    compId = next(k for k, v in extMonDict3.items() if k in self.compIdSet and v == compId)
                                    _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                                    if details is None:
                                        _atomNameLike[idx] = useOneLetterCodeSet = True
                                        _atomNameSpan[idx] = (index, len(_term))
                                        if resNameSpan[idx][0] == _atomNameSpan[idx][0]:
                                            resNameLike[idx] = False
                                        break
                                    if with_compid is not None and (atomId.startswith(with_compid)
                                                                    or (atomId in with_compid and index < resNameSpan[idx][1])):
                                        continue
                                    _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                                    _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                                    if details is None:
                                        _atomNameLike[idx] = True
                                        _atomNameSpan[idx] = (index, len(_term))
                                        if resNameSpan[idx][0] == _atomNameSpan[idx][0]:
                                            resNameLike[idx] = False
                                        break
                            for compId in self.compIdSet:
                                _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                                if details is None:
                                    _atomNameLike[idx] = True
                                    _atomNameSpan[idx] = (index, len(_term))
                                    break
                                if with_compid is not None and atomId.startswith(with_compid):
                                    continue
                                _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                                _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                                if details is None:
                                    _atomNameLike[idx] = True
                                    _atomNameSpan[idx] = (index, len(_term))
                                    break
                            if _atomNameLike[idx]:
                                break

            if numOfDim >= 3 and _atomNameLike[idx]:
                __term = term[0:_atomNameSpan[idx][0]]
                for elem in reversed(ASSIGNMENT_HALF_SPIN_NUCLEUS) if 'NH' in __term else ASSIGNMENT_HALF_SPIN_NUCLEUS:
                    if len(elem) == 1 and ligAtomId is None and _ligAtomId is None:
                        if elem in __term:
                            index = __term.rindex(elem)
                            atomId = __term[index:len(__term)]
                            if index - 1 >= 0 and __term[index - 1] in ASSIGNMENT_HALF_SPIN_NUCLEUS:
                                if resNameLike[idx] and compId[-1] in ASSIGNMENT_HALF_SPIN_NUCLEUS and index == resNameSpan[idx][1]:
                                    pass
                                else:
                                    continue
                            if atomId[0] in ('Q', 'M') and index + 1 < len(__term) and __term[index + 1].isdigit():
                                continue
                            if ((with_compid is not None and atomId.startswith(with_compid)) or atomId.startswith('MET'))\
                               and ((index + 3 < len(__term) and __term[index + 3].isdigit()
                                     or (index + 4 < len(__term) and __term[index + 4].isdigit()))):
                                continue
                            if resNameLike[idx] and len(compId) > 1 and compId[-1] == elem and index + 1 == resNameSpan[idx][1]:
                                continue
                            if len(__term) == _atomNameSpan[idx][0]:
                                continue
                            if resNameLike[idx]:
                                compId = term[resNameSpan[idx][0]:resNameSpan[idx][1]]
                                if len(compId) == 1 and hasOneLetterCodeSet:
                                    compId = next(k for k, v in extMonDict3.items() if k in self.compIdSet and v == compId)
                                    _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                                    if details is None:
                                        __atomNameLike[idx] = useOneLetterCodeSet = True
                                        __atomNameSpan[idx] = (index, len(__term))
                                        if resNameSpan[idx][0] == __atomNameSpan[idx][0]:
                                            resNameLike[idx] = False
                                        break
                                    if with_compid is not None and atomId.startswith(with_compid):
                                        continue
                                    _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                                    _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                                    if details is None:
                                        __atomNameLike[idx] = True
                                        __atomNameSpan[idx] = (index, len(__term))
                                        if resNameSpan[idx][0] == __atomNameSpan[idx][0]:
                                            resNameLike[idx] = False
                                        break
                            for compId in self.compIdSet:
                                _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                                if details is None:
                                    __atomNameLike[idx] = True
                                    __atomNameSpan[idx] = (index, len(__term))
                                    break
                                if with_compid is not None and atomId.startswith(with_compid):
                                    continue
                                _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                                _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                                if details is None:
                                    __atomNameLike[idx] = True
                                    __atomNameSpan[idx] = (index, len(__term))
                                    break
                            if __atomNameLike[idx]:
                                break

            if numOfDim >= 4 and __atomNameLike[idx]:
                ___term = term[0:__atomNameSpan[idx][0]]
                for elem in reversed(ASSIGNMENT_HALF_SPIN_NUCLEUS) if 'NH' in ___term else ASSIGNMENT_HALF_SPIN_NUCLEUS:
                    if len(elem) == 1 and ligAtomId is None and _ligAtomId is None:
                        if elem in ___term:
                            index = ___term.rindex(elem)
                            atomId = ___term[index:len(___term)]
                            if index - 1 >= 0 and ___term[index - 1] in ASSIGNMENT_HALF_SPIN_NUCLEUS:
                                if resNameLike[idx] and compId[-1] in ASSIGNMENT_HALF_SPIN_NUCLEUS and index == resNameSpan[idx][1]:
                                    pass
                                else:
                                    continue
                            if atomId[0] in ('Q', 'M') and index + 1 < len(___term) and ___term[index + 1].isdigit():
                                continue
                            if ((with_compid is not None and atomId.startswith(with_compid)) or atomId.startswith('MET'))\
                               and ((index + 3 < len(___term) and ___term[index + 3].isdigit()
                                     or (index + 4 < len(___term) and ___term[index + 4].isdigit()))):
                                continue
                            if resNameLike[idx] and len(compId) > 1 and compId[-1] == elem and index + 1 == resNameSpan[idx][1]:
                                continue
                            if len(___term) == __atomNameSpan[idx][0]:
                                continue
                            if resNameLike[idx]:
                                compId = term[resNameSpan[idx][0]:resNameSpan[idx][1]]
                                if len(compId) == 1 and hasOneLetterCodeSet:
                                    compId = next(k for k, v in extMonDict3.items() if k in self.compIdSet and v == compId)
                                    _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                                    if details is None:
                                        ___atomNameLike[idx] = useOneLetterCodeSet = True
                                        ___atomNameSpan[idx] = (index, len(___term))
                                        if resNameSpan[idx][0] == ___atomNameSpan[idx][0]:
                                            resNameLike[idx] = False
                                        break
                                    if with_compid is not None and atomId.startswith(with_compid):
                                        continue
                                    _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                                    _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                                    if details is None:
                                        ___atomNameLike[idx] = True
                                        ___atomNameSpan[idx] = (index, len(___term))
                                        if resNameSpan[idx][0] == ___atomNameSpan[idx][0]:
                                            resNameLike[idx] = False
                                        break
                            for compId in self.compIdSet:
                                _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                                if details is None:
                                    ___atomNameLike[idx] = True
                                    ___atomNameSpan[idx] = (index, len(___term))
                                    break
                                if with_compid is not None and atomId.startswith(with_compid):
                                    continue
                                _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                                _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                                if details is None:
                                    ___atomNameLike[idx] = True
                                    ___atomNameSpan[idx] = (index, len(___term))
                                    break
                            if ___atomNameLike[idx]:
                                break

            if _atomNameLike[idx] and atomNameLike[idx]:
                concat = False
                if _atomNameSpan[idx][1] == atomNameSpan[idx][0]:
                    atomId = term[_atomNameSpan[idx][0]:atomNameSpan[idx][1]]
                    if resNameLike[idx]:
                        compId = term[resNameSpan[idx][0]:resNameSpan[idx][1]]
                        resId = term[resIdSpan[idx][0]:resIdSpan[idx][1]] if resIdLike[idx] else '.'
                        if len(compId) == 1 and hasOneLetterCodeSet:
                            compId = next(k for k, v in extMonDict3.items() if k in self.compIdSet and v == compId)
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                            if details is None:
                                _atomNameLike[idx] = False
                                atomNameSpan[idx] = (_atomNameSpan[idx][0], atomNameSpan[idx][1])
                                if resNameSpan[idx][0] == atomNameSpan[idx][0]:
                                    resNameLike[idx] = False
                                concat = True
                            else:
                                _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                                _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                                if details is None:
                                    _atomNameLike[idx] = False
                                    atomNameSpan[idx] = (_atomNameSpan[idx][0], atomNameSpan[idx][1])
                                    if resNameSpan[idx][0] == atomNameSpan[idx][0]:
                                        resNameLike[idx] = False
                                    concat = True
                            if details is not None and numOfDim == 1:
                                _atomNameLike[idx] = False
                                for shift in range(_atomNameSpan[idx][0], atomNameSpan[idx][1]):
                                    _atomId = term[_atomNameSpan[idx][0] + shift:atomNameSpan[idx][1]]
                                    if len(_atomId) == 0:
                                        break
                                    if _atomId[0] == resId[0]:
                                        continue
                                    _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                                    if details is None:
                                        _atomNameSpan[idx] = (_atomNameSpan[idx][0] + shift, len(_atomId))
                                        atomNameSpan[idx] = (_atomNameSpan[idx][0], atomNameSpan[idx][1])
                                        concat = True
                                        break
                                    __atomId = translateToStdAtomName(_atomId, compId, ccU=self.ccU)
                                    _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, __atomId, leave_unmatched=True)
                                    if details is None:
                                        _atomNameSpan[idx] = (_atomNameSpan[idx][0] + shift, len(_atomId))
                                        atomNameSpan[idx] = (_atomNameSpan[idx][0], atomNameSpan[idx][1])
                                        concat = True
                                        break
                        else:
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                            if details is None:
                                _atomNameLike[idx] = False
                                atomNameSpan[idx] = (_atomNameSpan[idx][0], atomNameSpan[idx][1])
                                concat = True
                            else:
                                _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                                _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                                if details is None:
                                    _atomNameLike[idx] = False
                                    atomNameSpan[idx] = (_atomNameSpan[idx][0], atomNameSpan[idx][1])
                                    concat = True
                            if details is not None and numOfDim == 1:
                                _atomNameLike[idx] = False
                                for shift in range(_atomNameSpan[idx][0], atomNameSpan[idx][1]):
                                    _atomId = term[_atomNameSpan[idx][0] + shift:atomNameSpan[idx][1]]
                                    if len(_atomId) == 0:
                                        break
                                    if _atomId[0] == resId[0]:
                                        continue
                                    _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                                    if details is None:
                                        _atomNameSpan[idx] = (_atomNameSpan[idx][0] + shift, len(_atomId))
                                        atomNameSpan[idx] = (_atomNameSpan[idx][0], atomNameSpan[idx][1])
                                        concat = True
                                        break
                                    __atomId = translateToStdAtomName(_atomId, compId, ccU=self.ccU)
                                    _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, __atomId, leave_unmatched=True)
                                    if details is None:
                                        _atomNameSpan[idx] = (_atomNameSpan[idx][0] + shift, len(_atomId))
                                        atomNameSpan[idx] = (_atomNameSpan[idx][0], atomNameSpan[idx][1])
                                        concat = True
                                        break
                    if not concat:
                        for compId in self.compIdSet:
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                            if details is None:
                                _atomNameLike[idx] = False
                                atomNameSpan[idx] = (_atomNameSpan[idx][0], atomNameSpan[idx][1])
                                break
                            _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                            if details is None:
                                _atomNameLike[idx] = False
                                atomNameSpan[idx] = (_atomNameSpan[idx][0], atomNameSpan[idx][1])
                                break

            if __atomNameLike[idx] and _atomNameLike[idx] and atomNameLike[idx]:
                concat = False
                if __atomNameSpan[idx][1] == _atomNameSpan[idx][0] and _atomNameSpan[idx][1] == atomNameSpan[idx][0]:
                    atomId = term[__atomNameSpan[idx][0]:atomNameSpan[idx][1]]
                    if resNameLike[idx]:
                        compId = term[resNameSpan[idx][0]:resNameSpan[idx][1]]
                        if len(compId) == 1 and hasOneLetterCodeSet:
                            compId = next(k for k, v in extMonDict3.items() if k in self.compIdSet and v == compId)
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                            if details is None:
                                __atomNameLike[idx] = _atomNameLike[idx] = False
                                atomNameSpan[idx] = (__atomNameSpan[idx][0], atomNameSpan[idx][1])
                                if resNameSpan[idx][0] == atomNameSpan[idx][0]:
                                    resNameLike[idx] = False
                                concat = True
                            else:
                                _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                                _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                                if details is None:
                                    __atomNameLike[idx] = _atomNameLike[idx] = False
                                    atomNameSpan[idx] = (__atomNameSpan[idx][0], atomNameSpan[idx][1])
                                    if resNameSpan[idx][0] == atomNameSpan[idx][0]:
                                        resNameLike[idx] = False
                                    concat = True

                    if not concat:
                        for compId in self.compIdSet:
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                            if details is None:
                                __atomNameLike[idx] = _atomNameLike[idx] = False
                                atomNameSpan[idx] = (__atomNameSpan[idx][0], atomNameSpan[idx][1])
                                break
                            _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                            if details is None:
                                __atomNameLike[idx] = _atomNameLike[idx] = False
                                atomNameSpan[idx] = (__atomNameSpan[idx][0], atomNameSpan[idx][1])
                                break

            if ___atomNameLike[idx] and __atomNameLike[idx] and _atomNameLike[idx] and atomNameLike[idx]:
                concat = False
                if ___atomNameSpan[idx][1] == __atomNameSpan[idx][0] and __atomNameSpan[idx][1] == _atomNameSpan[idx][0]\
                   and _atomNameSpan[idx][1] == atomNameSpan[idx][0]:
                    atomId = term[___atomNameSpan[idx][0]:atomNameSpan[idx][1]]
                    if resNameLike[idx]:
                        compId = term[resNameSpan[idx][0]:resNameSpan[idx][1]]
                        if len(compId) == 1 and hasOneLetterCodeSet:
                            compId = next(k for k, v in extMonDict3.items() if k in self.compIdSet and v == compId)
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                            if details is None:
                                ___atomNameLike[idx] = __atomNameLike[idx] = _atomNameLike[idx] = False
                                atomNameSpan[idx] = (___atomNameSpan[idx][0], atomNameSpan[idx][1])
                                if resNameSpan[idx][0] == atomNameSpan[idx][0]:
                                    resNameLike[idx] = False
                                concat = True
                            else:
                                _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                                _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                                if details is None:
                                    ___atomNameLike[idx] = __atomNameLike[idx] = _atomNameLike[idx] = False
                                    atomNameSpan[idx] = (___atomNameSpan[idx][0], atomNameSpan[idx][1])
                                    if resNameSpan[idx][0] == atomNameSpan[idx][0]:
                                        resNameLike[idx] = False
                                    concat = True
                    if not concat:
                        for compId in self.compIdSet:
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                            if details is None:
                                ___atomNameLike[idx] = __atomNameLike[idx] = _atomNameLike[idx] = False
                                atomNameSpan[idx] = (___atomNameSpan[idx][0], atomNameSpan[idx][1])
                                break
                            _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                            if details is None:
                                ___atomNameLike[idx] = __atomNameLike[idx] = _atomNameLike[idx] = False
                                atomNameSpan[idx] = (___atomNameSpan[idx][0], atomNameSpan[idx][1])
                                break

            if __atomNameLike[idx] and _atomNameLike[idx]:
                concat = False
                if __atomNameSpan[idx][1] == _atomNameSpan[idx][0]:
                    atomId = term[__atomNameSpan[idx][0]:_atomNameSpan[idx][1]]
                    if resNameLike[idx]:
                        compId = term[resNameSpan[idx][0]:resNameSpan[idx][1]]
                        if len(compId) == 1 and hasOneLetterCodeSet:
                            compId = next(k for k, v in extMonDict3.items() if k in self.compIdSet and v == compId)
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                            if details is None:
                                __atomNameLike[idx] = False
                                _atomNameSpan[idx] = (__atomNameSpan[idx][0], _atomNameSpan[idx][1])
                                if resNameSpan[idx][0] == _atomNameSpan[idx][0]:
                                    resNameLike[idx] = False
                                concat = True
                            else:
                                _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                                _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                                if details is None:
                                    __atomNameLike[idx] = False
                                    _atomNameSpan[idx] = (__atomNameSpan[idx][0], _atomNameSpan[idx][1])
                                    if resNameSpan[idx][0] == _atomNameSpan[idx][0]:
                                        resNameLike[idx] = False
                                    concat = True
                    if not concat:
                        for compId in self.compIdSet:
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                            if details is None:
                                __atomNameLike[idx] = False
                                _atomNameSpan[idx] = (__atomNameSpan[idx][0], _atomNameSpan[idx][1])
                                break
                            _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                            if details is None:
                                __atomNameLike[idx] = False
                                _atomNameSpan[idx] = (__atomNameSpan[idx][0], _atomNameSpan[idx][1])
                                break

            if ___atomNameLike[idx] and __atomNameLike[idx] and _atomNameLike[idx]:
                concat = False
                if ___atomNameSpan[idx][1] == __atomNameSpan[idx][0] and __atomNameSpan[idx][1] == _atomNameSpan[idx][0]:
                    atomId = term[___atomNameSpan[idx][0]:_atomNameSpan[idx][1]]
                    if resNameLike[idx]:
                        compId = term[resNameSpan[idx][0]:resNameSpan[idx][1]]
                        if len(compId) == 1 and hasOneLetterCodeSet:
                            compId = next(k for k, v in extMonDict3.items() if k in self.compIdSet and v == compId)
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                            if details is None:
                                ___atomNameLike[idx] = __atomNameLike[idx] = False
                                _atomNameSpan[idx] = (___atomNameSpan[idx][0], _atomNameSpan[idx][1])
                                if resNameSpan[idx][0] == _atomNameSpan[idx][0]:
                                    resNameLike[idx] = False
                                concat = True
                            else:
                                _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                                _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                                if details is None:
                                    ___atomNameLike[idx] = __atomNameLike[idx] = False
                                    _atomNameSpan[idx] = (___atomNameSpan[idx][0], _atomNameSpan[idx][1])
                                    if resNameSpan[idx][0] == _atomNameSpan[idx][0]:
                                        resNameLike[idx] = False
                                    concat = True
                    if not concat:
                        for compId in self.compIdSet:
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                            if details is None:
                                ___atomNameLike[idx] = __atomNameLike[idx] = False
                                _atomNameSpan[idx] = (___atomNameSpan[idx][0], _atomNameSpan[idx][1])
                                break
                            _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                            if details is None:
                                ___atomNameLike[idx] = __atomNameLike[idx] = False
                                _atomNameSpan[idx] = (___atomNameSpan[idx][0], _atomNameSpan[idx][1])
                                break

            if ___atomNameLike[idx] and __atomNameLike[idx]:
                concat = False
                if ___atomNameSpan[idx][1] == __atomNameSpan[idx][0]:
                    atomId = term[___atomNameSpan[idx][0]:__atomNameSpan[idx][1]]
                    if resNameLike[idx]:
                        compId = term[resNameSpan[idx][0]:resNameSpan[idx][1]]
                        if len(compId) == 1 and hasOneLetterCodeSet:
                            compId = next(k for k, v in extMonDict3.items() if k in self.compIdSet and v == compId)
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                            if details is None:
                                ___atomNameLike[idx] = False
                                __atomNameSpan[idx] = (___atomNameSpan[idx][0], __atomNameSpan[idx][1])
                                if resNameSpan[idx][0] == __atomNameSpan[idx][0]:
                                    resNameLike[idx] = False
                                concat = True
                            else:
                                _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                                _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                                if details is None:
                                    ___atomNameLike[idx] = False
                                    __atomNameSpan[idx] = (___atomNameSpan[idx][0], __atomNameSpan[idx][1])
                                    if resNameSpan[idx][0] == __atomNameSpan[idx][0]:
                                        resNameLike[idx] = False
                                    concat = True
                    if not concat:
                        for compId in self.compIdSet:
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, atomId, leave_unmatched=True)
                            if details is None:
                                ___atomNameLike[idx] = False
                                __atomNameSpan[idx] = (___atomNameSpan[idx][0], __atomNameSpan[idx][1])
                                break
                            _atomId = translateToStdAtomName(atomId, compId, ccU=self.ccU)
                            _, _, details = self.nefT.get_valid_star_atom_in_xplor(compId, _atomId, leave_unmatched=True)
                            if details is None:
                                ___atomNameLike[idx] = False
                                __atomNameSpan[idx] = (___atomNameSpan[idx][0], __atomNameSpan[idx][1])
                                break

            if idx == 0 and atomNameLike[0] and not _atomNameLike[0] and self.hasNonPoly:
                if not segIdLike[0] and not resIdLike[0] and not resNameLike[0]:
                    ligands = 0
                    for np in self.nonPoly:
                        if 'alt_comp_id' in np and _str[0][0] == np['alt_comp_id'][0][0]:
                            ligands += 1
                    if ligands == 1:
                        for np in self.nonPoly:
                            if 'alt_comp_id' in np and _str[0][0] == np['alt_comp_id'][0][0]:
                                _string = f"{np['auth_chain_id']} {np['auth_seq_id'][0]} "\
                                    f"{np['comp_id'][0]}{string[atomNameSpan[0][1]:]}"
                                return self.extractPeakAssignment(numOfDim, _string, src_index,
                                                                  with_segid, with_compid, hint, dim_id_hint)

                if _str[0][0] == 'X' and term[atomNameSpan[idx][0]] != 'X' and len(self.nonPoly) == 1:
                    np = self.nonPoly[0]
                    _string = f"{np['auth_chain_id']} {np['auth_seq_id'][0]} {np['comp_id'][0]}{string[atomNameSpan[0][0]:]}"
                    return self.extractPeakAssignment(numOfDim, _string, src_index, with_segid, with_compid, hint, dim_id_hint)

            if self.__verbose_debug:
                print(f'{idx} {term!r} segid:{segIdLike[idx]} {term[segIdSpan[idx][0]:segIdSpan[idx][1]] if segIdLike[idx] else ""}, '  # noqa: E501, pylint: disable=line-too-long
                      f'resid:{resIdLike[idx]} {term[resIdSpan[idx][0]:resIdSpan[idx][1]] if resIdLike[idx] else ""}, '
                      f'resname:{resNameLike[idx]} {term[resNameSpan[idx][0]:resNameSpan[idx][1]] if resNameLike[idx] else ""}, '
                      f'atomname:{atomNameLike[idx]} {term[atomNameSpan[idx][0]:atomNameSpan[idx][1]] if atomNameLike[idx] else ""}, '  # noqa: E501, pylint: disable=line-too-long
                      f'_atomname:{_atomNameLike[idx]} {term[_atomNameSpan[idx][0]:_atomNameSpan[idx][1]] if _atomNameLike[idx] else ""},  '  # noqa: E501, pylint: disable=line-too-long
                      f'__atomname:{__atomNameLike[idx]} {term[__atomNameSpan[idx][0]:__atomNameSpan[idx][1]] if __atomNameLike[idx] else ""}, '  # noqa: E501, pylint: disable=line-too-long
                      f'___atomname:{___atomNameLike[idx]} {term[___atomNameSpan[idx][0]:___atomNameSpan[idx][1]] if ___atomNameLike[idx] else ""}')  # noqa: E501, pylint: disable=line-too-long

        atomNameCount = 0
        for idx in range(lenStr):
            if atomNameLike[idx]:
                atomNameCount += 1
            if _atomNameLike[idx]:
                atomNameCount += 1
            if __atomNameLike[idx]:
                atomNameCount += 1
            if ___atomNameLike[idx]:
                atomNameCount += 1

        if atomNameCount < numOfDim:
            return None

        if atomNameCount > numOfDim > 1:
            atomNameCount = 0
            ignoreBefore = False
            for idx in range(lenStr - 1, 0, -1):
                if ignoreBefore:
                    atomNameLike[idx] = _atomNameLike[idx] = __atomNameLike[idx] = ___atomNameLike[idx] = False
                else:
                    if atomNameLike[idx]:
                        atomNameCount += 1
                    if _atomNameLike[idx]:
                        atomNameCount += 1
                    if __atomNameLike[idx]:
                        atomNameCount += 1
                    if ___atomNameLike[idx]:
                        atomNameCount += 1
                    if atomNameCount >= numOfDim:
                        ignoreBefore = True

        hasResName = hasResId = False
        for idx in range(lenStr):
            if ___atomNameLike[idx]:
                if resNameLike[idx]:
                    if not hasResName and resNameSpan[idx][0] < ___atomNameSpan[idx][0]\
                       and resNameSpan[idx][1] >= ___atomNameSpan[idx][1]:
                        ___atomNameLike[idx] = False
                    elif resNameSpan[idx][1] > ___atomNameSpan[idx][0]:
                        term = _str[idx]
                        if not hasResId and resNameSpan[idx][1] - resNameSpan[idx][0] == 3\
                           and term[resNameSpan[idx][0]:resNameSpan[idx][1]] ==\
                           term[___atomNameSpan[idx][0]:___atomNameSpan[idx][1]]\
                           and any(True for _idx in range(idx + 1, lenStr)
                                   if __atomNameLike[_idx] or _atomNameLike[_idx] or atomNameLike[_idx]):
                            ___atomNameLike[idx] = False
                        else:
                            resNameLike[idx] = False
                if resIdLike[idx]:
                    if resIdSpan[idx][1] > ___atomNameSpan[idx][0]:
                        resIdLike[idx] = False
                if segIdLike[idx]:
                    if segIdSpan[idx][1] > ___atomNameSpan[idx][0]:
                        segIdLike[idx] = False

            elif __atomNameLike[idx]:
                if resNameLike[idx]:
                    if not hasResName and resNameSpan[idx][0] < __atomNameSpan[idx][0]\
                       and resNameSpan[idx][1] >= __atomNameSpan[idx][1]:
                        __atomNameLike[idx] = False
                    elif resNameSpan[idx][1] > __atomNameSpan[idx][0]:
                        term = _str[idx]
                        if not hasResId and resNameSpan[idx][1] - resNameSpan[idx][0] == 3\
                           and term[resNameSpan[idx][0]:resNameSpan[idx][1]] == term[__atomNameSpan[idx][0]:__atomNameSpan[idx][1]]\
                           and any(True for _idx in range(idx + 1, lenStr) if _atomNameLike[_idx] or atomNameLike[_idx]):
                            __atomNameLike[idx] = False
                        else:
                            resNameLike[idx] = False
                if resIdLike[idx]:
                    if resIdSpan[idx][1] > __atomNameSpan[idx][0]:
                        resIdLike[idx] = False
                if segIdLike[idx]:
                    if segIdSpan[idx][1] > __atomNameSpan[idx][0]:
                        segIdLike[idx] = False

            elif _atomNameLike[idx]:
                if resNameLike[idx]:
                    if not hasResName and resNameSpan[idx][0] < _atomNameSpan[idx][0]\
                       and resNameSpan[idx][1] >= _atomNameSpan[idx][1]:
                        _atomNameLike[idx] = False
                    elif resNameSpan[idx][1] > _atomNameSpan[idx][0]:
                        term = _str[idx]
                        if not hasResId and resNameSpan[idx][1] - resNameSpan[idx][0] == 3\
                           and term[resNameSpan[idx][0]:resNameSpan[idx][1]] == term[_atomNameSpan[idx][0]:_atomNameSpan[idx][1]]\
                           and any(True for _idx in range(idx + 1, lenStr) if atomNameLike[_idx]):
                            _atomNameLike[idx] = False
                        else:
                            resNameLike[idx] = False
                if resIdLike[idx]:
                    if resIdSpan[idx][1] > _atomNameSpan[idx][0]:
                        resIdLike[idx] = False
                if segIdLike[idx]:
                    if segIdSpan[idx][1] > _atomNameSpan[idx][0]:
                        segIdLike[idx] = False

            elif atomNameLike[idx]:
                if resNameLike[idx]:
                    if not hasResName and resNameSpan[idx][0] < atomNameSpan[idx][0]\
                       and resNameSpan[idx][1] >= atomNameSpan[idx][1]:
                        atomNameLike[idx] = False
                    elif resNameSpan[idx][1] > atomNameSpan[idx][0]:
                        term = _str[idx]
                        if not hasResId and resNameSpan[idx][1] - resNameSpan[idx][0] == 3\
                           and term[resNameSpan[idx][0]:resNameSpan[idx][1]] == term[atomNameSpan[idx][0]:atomNameSpan[idx][1]]\
                           and any(True for _idx in range(idx + 1, lenStr) if atomNameLike[_idx]):
                            atomNameLike[idx] = False
                        else:
                            if forceOneLetterCodeSet and resIdLike[idx] and len(atomNameLike) > idx + 1 and atomNameLike[idx + 1]:
                                atomNameLike[idx] = False
                            else:
                                resNameLike[idx] = False
                if resIdLike[idx]:
                    if resIdSpan[idx][1] > atomNameSpan[idx][0]:
                        if forceOneLetterCodeSet and resNameLike[idx]:
                            pass
                        else:
                            resIdLike[idx] = False
                if segIdLike[idx]:
                    if segIdSpan[idx][1] > atomNameSpan[idx][0]:
                        segIdLike[idx] = False

            if resNameLike[idx]:
                hasResName = True
                if segIdLike[idx]:
                    if segIdSpan[idx][1] > resNameSpan[idx][0]:
                        if numOfDim > 1 or not any(resNameLike[_idx] for _idx in range(idx + 1, lenStr))\
                           or (resIdLike[idx] and segIdSpan[idx][1] == resIdSpan[idx][0]):
                            _resName = _str[idx][resNameSpan[idx][0]:resNameSpan[idx][1]]
                            _resId = next((int(_str[_idx][resIdSpan[_idx][0]:resIdSpan[_idx][1]]) for _idx in range(idx + 1, lenStr)
                                           if resIdLike[_idx]), None)
                            _atomName = next((_str[_idx][atomNameSpan[_idx][0]:atomNameSpan[_idx][1]]
                                              for _idx in range(idx + 1, lenStr)
                                              if atomNameLike[_idx]), None)
                            checked = True
                            if _resId is not None and _atomName is not None and len(_resName) == 1 and len(self.polySeq) > 1:
                                for ps in self.polySeq:
                                    _, _, _compId_ = self.getRealChainSeqId(ps, _resId, None)
                                    if _resName == ps['auth_chain_id'] and _compId_ in STD_MON_DICT\
                                            and _resName != STD_MON_DICT[_compId_]:
                                        checked = False
                                        break
                            if checked:
                                segIdLike[idx] = False
                            else:
                                resNameLike[idx] = False
                        else:
                            if not hasResId and resNameSpan[idx][1] - resNameSpan[idx][0] == 3\
                               and term[resNameSpan[idx][0]:resNameSpan[idx][1]] == term[atomNameSpan[idx][0]:atomNameSpan[idx][1]]\
                               and any(True for _idx in range(idx + 1, lenStr) if atomNameLike[_idx]):
                                atomNameLike[idx] = False
                            else:
                                resNameLike[idx] = False

            if resIdLike[idx]:
                hasResId = True
                if atomNameLike[idx]:
                    if resIdSpan[idx][1] > atomNameSpan[idx][0]:
                        resIdLike[idx] = False

            if self.__verbose_debug:
                print(f' -> {idx} segid:{segIdLike[idx]}, resid:{resIdLike[idx]}, resname:{resNameLike[idx]}, '
                      f'atomname:{atomNameLike[idx]}, _atomname:{_atomNameLike[idx]}, '
                      f'__atomname:{__atomNameLike[idx]}, ___atomname:{___atomNameLike[idx]}')

        resIdCount = 0
        for idx in range(lenStr):
            if resIdLike[idx]:
                resIdCount += 1

        _resId = [h['seq_id'] for h in hint] if hint is not None else None
        if resIdCount == 0:
            if _resId is None and _ligAtomId is None:
                return None

        _resNameDict = [{h['auth_seq_id']: h['comp_id']} for h in hint] if hint is not None else None

        resIdLater = resIdCount == numOfDim
        if resIdLater:
            atomNameCount = 0
            for idx in range(lenStr):
                if atomNameLike[idx]:
                    atomNameCount += 1
            resIdLater = atomNameCount == numOfDim
            if resIdLater:
                anyResId = False
                for idx in range(lenStr):
                    if resIdLike[idx]:
                        anyResId = True
                    if atomNameLike[idx]:
                        if anyResId:
                            resIdLater = False
                        break

        if self.__verbose_debug:
            print(f'num_of_dim: {numOfDim}, resid_count: {resIdCount}, resid_later:{resIdLater}')

        def is_valid_chain_assign(chain_assign, res_name):
            return len(chain_assign) > 0 and ((res_name in extMonDict3 and any(True for a in chain_assign if a[2] == res_name))
                                              or res_name not in extMonDict3)

        ret = []

        segId = _segId = resId = resName = atomName = _segId_ = _resId_ = authResId = None
        dimId = 1

        def append_atom_assignment() -> bool:
            """ Assign the current atom (atomName) to the coordinates and append it to ret.
                @return: whether the atom is skipped as an ignorable diagonal assignment
            """
            nonlocal _, atomName, details, dimId, resName, segId

            if self.__hasCoord:
                if segId is None and resName is None:
                    chainAssign =\
                        self.assignCoordPolymerSequenceWithChainIdWithoutCompId(self.__defaultSegId, resId, atomName, src_index)
                    if len(chainAssign) > 0:
                        if self.__defaultSegId is None:
                            ca_idx = next((chainAssign.index(a) for a in chainAssign if a[1] == resId), -1)
                        else:
                            ca_idx = next((chainAssign.index(a) for a in chainAssign
                                           if a[0] == self.__defaultSegId and a[1] == resId), -1)
                            if ca_idx == -1:
                                ca_idx = next((chainAssign.index(a) for a in chainAssign if a[1] == resId), -1)
                        if ca_idx != -1:
                            # if self.__defaultSegId is None:
                            self.__defaultSegId = chainAssign[ca_idx][0]
                        else:
                            ca_idx = 0
                        segId, _, resName, _ = chainAssign[ca_idx]
                    else:
                        chainAssign = self.assignCoordPolymerSequenceWithoutCompId(resId, atomName, src_index)
                        if len(chainAssign) > 0:
                            if self.__defaultSegId is None:
                                ca_idx = next((chainAssign.index(a) for a in chainAssign if a[1] == resId), -1)
                            else:
                                ca_idx = next((chainAssign.index(a) for a in chainAssign
                                               if a[0] == self.__defaultSegId and a[1] == resId), -1)
                                if ca_idx == -1:
                                    ca_idx = next((chainAssign.index(a) for a in chainAssign if a[1] == resId), -1)
                            if ca_idx != -1:
                                # if self.__defaultSegId is None:
                                # self.__defaultSegId = chainAssign[ca_idx][0]
                                pass
                            else:
                                ca_idx = 0
                            segId, _, resName, _ = chainAssign[ca_idx]
                elif segId is None:
                    chainAssign, _ = self.assignCoordPolymerSequence(self.__defaultSegId,
                                                                     resId, resName, atomName, src_index)
                    is_valid = is_valid_chain_assign(chainAssign, resName)
                    if not is_valid:
                        if self.__defaultSegId != self.__defaultSegId__:
                            __preferAuthSeq = self.__preferAuthSeq
                            self.__preferAuthSeq = not __preferAuthSeq
                            chainAssign, _ = self.assignCoordPolymerSequenceWithChainId(self.__defaultSegId__,
                                                                                        resId, resName, atomName, src_index)
                            is_valid = is_valid_chain_assign(chainAssign, resName)
                            self.__preferAuthSeq = __preferAuthSeq
                    if is_valid:
                        ca_idx = next((chainAssign.index(a) for a in chainAssign if a[2] == resName), -1)
                        if ca_idx != -1:
                            # if self.__defaultSegId is None:
                            self.__defaultSegId = chainAssign[ca_idx][0]
                        else:
                            ca_idx = 0
                        segId = chainAssign[ca_idx][0]
                    else:
                        chainAssign, _ = self.assignCoordPolymerSequence(None,
                                                                         resId, resName, atomName, src_index)
                        if is_valid_chain_assign(chainAssign, resName):
                            ca_idx = next((chainAssign.index(a) for a in chainAssign if a[2] == resName), -1)
                            if ca_idx != -1:
                                # if self.__defaultSegId is None:
                                self.__defaultSegId = chainAssign[ca_idx][0]
                                if self.reasons is None:
                                    r = self.__getNamedReasonsForReparsing('default_seg_id')
                                    if self.num_of_dim not in r:
                                        r[self.num_of_dim] = {}
                                    if self.cur_list_id != -1 and self.cur_list_id not in r[self.num_of_dim]:
                                        r[self.num_of_dim][self.cur_list_id] = self.__defaultSegId
                            else:
                                ca_idx = 0
                            segId = chainAssign[ca_idx][0]
                elif resName is None:
                    chainAssign = self.assignCoordPolymerSequenceWithChainIdWithoutCompId(segId, resId, atomName, src_index)
                    if len(chainAssign) > 0:
                        ca_idx = next((chainAssign.index(a) for a in chainAssign if a[1] == resId), 0)
                        resName = chainAssign[ca_idx][2]
                _, _, details = self.nefT.get_valid_star_atom_in_xplor(resName, atomName, leave_unmatched=True)
                if details is not None:
                    atomName = translateToStdAtomName(atomName, resName, ccU=self.ccU)
                if any(True for item in ret
                       if item['chain_id'] == segId and item['seq_id'] == resId and item['atom_id'] == atomName):
                    if self.__ignore_diagonal:
                        return True
                ret.append({'dim_id': dimId, 'chain_id': segId, 'auth_chain_id': _segId, 'seq_id': resId,
                            'auth_seq_id': authResId, 'comp_id': resName, 'atom_id': atomName})
            else:
                if any(True for item in ret
                       if (segId is None or item['chain_id'] == segId)
                       and item['seq_id'] == resId and item['atom_id'] == atomName):
                    if self.__ignore_diagonal:
                        return True
                ass = {'dim': dimId, 'atom_id': atomName}
                if _segId is not None:
                    ass['auth_chain_id'] = _segId
                if segId is not None:
                    ass['chain_id'] = segId
                if resId is not None:
                    ass['seq_id'] = resId
                if resName is not None:
                    ass['comp_id'] = resName
                ret.append(ass)
            dimId += 1
            return False

        for idx, term in enumerate(_str):
            if segIdLike[idx]:
                segId = _segId = term[segIdSpan[idx][0]:segIdSpan[idx][1]]
                if _segId.isdigit():
                    _segId = None
                if _segId_ is not None and segId != _segId_:
                    resId = resName = None
                _segId_ = segId
            if resIdLike[idx]:
                resId = authResId = int(term[resIdSpan[idx][0]:resIdSpan[idx][1]])
                if _resId_ is not None and resId != _resId_:
                    resName = None
                _resId_ = resId
            if resNameLike[idx]:
                resName = term[resNameSpan[idx][0]:resNameSpan[idx][1]]
                if len(resName) == 1 and hasOneLetterCodeSet:
                    resName = next(k for k, v in extMonDict3.items() if k in self.compIdSet and v == resName)
            elif _resNameDict is not None and resName is None and len(ret) < len(_resNameDict)\
                    and _resId is not None and len(ret) < len(_resId) and _resId[len(ret)] in _resNameDict[len(ret)]\
                    and (authResId is None or authResId == _resId[len(ret)]):
                resName = _resNameDict[len(ret)][_resId[len(ret)]]
            if ___atomNameLike[idx]:
                if resIdLater:
                    for _idx, _term in enumerate(_str):
                        if _idx > idx and resIdLike[_idx]:
                            resId = int(_term[resIdSpan[_idx][0]:resIdSpan[_idx][1]])
                            segId = resName = None
                            break
                if resId is None:
                    if _resId is None or len(ret) >= len(_resId):
                        return None
                    resId = _resId[len(ret)]
                atomName = term[___atomNameSpan[idx][0]:___atomNameSpan[idx][1]]
                if append_atom_assignment():
                    continue
            if __atomNameLike[idx]:
                if resIdLater:
                    for _idx, _term in enumerate(_str):
                        if _idx > idx and resIdLike[_idx]:
                            resId = int(_term[resIdSpan[_idx][0]:resIdSpan[_idx][1]])
                            segId = resName = None
                            break
                if resId is None:
                    if _resId is None or len(ret) >= len(_resId):
                        return None
                    resId = _resId[len(ret)]
                atomName = term[__atomNameSpan[idx][0]:__atomNameSpan[idx][1]]
                if append_atom_assignment():
                    continue
            if _atomNameLike[idx]:
                if resIdLater:
                    for _idx, _term in enumerate(_str):
                        if _idx > idx and resIdLike[_idx]:
                            resId = int(_term[resIdSpan[_idx][0]:resIdSpan[_idx][1]])
                            segId = resName = None
                            break
                if resId is None:
                    if _resId is None or len(ret) >= len(_resId):
                        return None
                    resId = _resId[len(ret)]
                atomName = term[_atomNameSpan[idx][0]:_atomNameSpan[idx][1]]
                if append_atom_assignment():
                    continue
            if atomNameLike[idx]:
                if resIdLater:
                    for _idx, _term in enumerate(_str):
                        if _idx > idx and resIdLike[_idx]:
                            resId = int(_term[resIdSpan[_idx][0]:resIdSpan[_idx][1]])
                            segId = resName = None
                            break
                if resId is None and _ligAtomId is None:
                    if _resId is None or len(ret) >= len(_resId):
                        return None
                    resId = _resId[len(ret)]
                if _ligAtomId is not None:
                    resId = _ligSeqId
                atomName = term[atomNameSpan[idx][0]:atomNameSpan[idx][1]]
                if append_atom_assignment():
                    continue
            elif atomNameLike_[idx] and siblingAtomName[idx] is not None:
                if resIdLater:
                    for _idx, _term in enumerate(_str):
                        if _idx > idx and resIdLike[_idx]:
                            resId = int(_term[resIdSpan[_idx][0]:resIdSpan[_idx][1]])
                            segId = resName = None
                            break
                if resId is None:
                    if _resId is None or len(ret) >= len(_resId):
                        return None
                    resId = _resId[len(ret)]
                for atomName in siblingAtomName[idx]:
                    if append_atom_assignment():
                        continue

        multiple = len(ret) > numOfDim

        if multiple:
            if self.__createSfDict and self.use_peak_row_format and not self.__enforcePeakRowFormat:
                sf = self.getSf()
                sf['peak_row_format'] = self.use_peak_row_format = False

        return ret if len(ret) > 0 else None

    def getRealChainSeqId(self, ps: dict, seqId: int, compId: Optional[str], isPolySeq: bool = True,
                          isFirstTrial: bool = True
                          ) -> Tuple[str, int, Optional[str]]:
        """ Return a realistic sequence for a given polymer sequence, sequence code, and residue name.
        """

        if compId is not None:
            compId = _compId = translateToStdResName(compId, ccU=self.ccU)
            if len(_compId) == 2 and _compId.startswith('D'):
                _compId = compId[1]
        if not self.__preferAuthSeq:
            seqKey = (ps['auth_chain_id'], seqId)
            if seqKey in self.__labelToAuthSeq:
                _chainId, _seqId = self.__labelToAuthSeq[seqKey]
                if _seqId in ps['auth_seq_id']:
                    return _chainId, _seqId, ps['comp_id'][ps['seq_id'].index(seqId)]
                if seqKey[1] in ps['seq_id']:  # resolve conflict between label/auth sequence schemes of polymer/non-polymer (2l90)
                    idx = ps['seq_id'].index(seqKey[1])
                    return _chainId, ps['auth_seq_id'][idx], ps['comp_id'][idx]
        if seqId in ps['auth_seq_id']:
            if compId is None:
                return ps['auth_chain_id'], seqId, ps['comp_id'][ps['auth_seq_id'].index(seqId)]
            for idx in [_idx for _idx, _seqId in enumerate(ps['auth_seq_id']) if _seqId == seqId]:
                if 'alt_comp_id' in ps and idx < len(ps['alt_comp_id']):
                    if compId in (ps['comp_id'][idx], ps['auth_comp_id'][idx], ps['alt_comp_id'][idx]):
                        return ps['auth_chain_id'], seqId, ps['comp_id'][idx]
                    if compId != _compId and _compId in (ps['comp_id'][idx], ps['auth_comp_id'][idx], ps['alt_comp_id'][idx]):
                        return ps['auth_chain_id'], seqId, ps['comp_id'][idx]
                if compId in (ps['comp_id'][idx], ps['auth_comp_id'][idx])\
                   or (isPolySeq and seqId == 1
                       and ((compId.endswith('-N') and all(c in ps['comp_id'][idx] for c in compId.split('-')[0]))
                            or (ps['comp_id'][idx] == 'PCA' and 'P' == compId[0] and ('GL' in compId or 'N' in compId)))):
                    return ps['auth_chain_id'], seqId, ps['comp_id'][idx]
                if compId != _compId and _compId in (ps['comp_id'][idx], ps['auth_comp_id'][idx]):
                    return ps['auth_chain_id'], seqId, ps['comp_id'][idx]
        if self.reasons is not None and 'extend_seq_scheme' in self.reasons:
            _ps = next((_ps for _ps in self.reasons['extend_seq_scheme'] if _ps['chain_id'] == ps['auth_chain_id']), None)
            if _ps is not None:
                if seqId in _ps['seq_id']:
                    return ps['auth_chain_id'], seqId, _ps['comp_id'][_ps['seq_id'].index(seqId)]
        if 'Check the 1th row of' in self.getCurrentSpectralPeak(-1) and isFirstTrial and isPolySeq\
           and (self.reasons is None
                or not ('seq_id_remap' in self.reasons or 'chain_seq_id_remap' in self.reasons
                        or 'ext_chain_seq_id_remap' in self.reasons)):
            try:
                if not any(_ps['auth_seq_id'][0] - len(_ps['seq_id']) <= seqId <= _ps['auth_seq_id'][-1] + len(_ps['seq_id'])
                           for _ps in self.polySeq):
                    self.__preferAuthSeq = not self.__preferAuthSeq
                    trial = self.getRealChainSeqId(ps, seqId, compId, isPolySeq, False)
                    if trial[2] is not None and compId == trial[2]:
                        return trial
                    self.__preferAuthSeq = not self.__preferAuthSeq
            except ValueError:
                pass
        return ps['auth_chain_id'], seqId, None

    @functools.lru_cache(maxsize=256)
    def translateToStdResNameWrapper(self, seqId: int, compId: str, preferNonPoly: bool = False) -> str:
        """ A wrapper function for ParserListenerUtil.translateToStdResName().
        """

        _compId = compId
        refCompId = None
        for ps in self.polySeq:
            if preferNonPoly:
                continue
            _, _, refCompId = self.getRealChainSeqId(ps, seqId, _compId)
            if refCompId is not None:
                compId = translateToStdResName(_compId, refCompId=refCompId, ccU=self.ccU)
                if compId != _compId and compId in STD_MON_DICT and _compId in STD_MON_DICT:
                    continue
                break
        if refCompId is None and self.hasNonPolySeq:
            for np in self.nonPolySeq:
                _, _, refCompId = self.getRealChainSeqId(np, seqId, _compId, False)
                if refCompId is not None:
                    compId = translateToStdResName(_compId, refCompId=refCompId, ccU=self.ccU)
                    break
        if refCompId is None:
            compId = translateToStdResName(_compId, ccU=self.ccU)
        return compId

    def assignCoordPolymerSequence(self, refChainId: str, seqId: int, compId: str, atomId: str, index: int
                                   ) -> Tuple[List[Tuple[str, int, str, bool]], bool]:
        """ Assign polymer sequences of the coordinates.
        """

        _refChainId = refChainId

        chainAssign = set()
        asis = False

        _seqId = seqId
        _compId = compId

        fixedChainId = fixedSeqId = fixedCompId = None

        preferNonPoly = False

        self.__allow_ext_seq = False

        if self.hasNonPoly:

            resolved = False

            for np in self.nonPoly:
                if 'alt_comp_id' in np and 'alt_auth_seq_id' in np\
                   and compId in np['alt_comp_id'] and seqId in np['alt_auth_seq_id']:
                    npCompId = np['comp_id'][0]
                    npSeqId = np['auth_seq_id'][0]
                    for ps in self.polySeq:
                        if 'ambig_auth_seq_id' in ps and seqId in ps['ambig_auth_seq_id']:
                            psCompId = ps['comp_id'][ps['auth_seq_id'].index(seqId)]
                            _atomId, _, details = self.nefT.get_valid_star_atom_in_xplor(psCompId, atomId, leave_unmatched=True)
                            if details is None:
                                _, _coordAtomSite =\
                                    self.getCoordAtomSiteOf(ps['auth_chain_id'], seqId, psCompId, cifCheck=self.__hasCoord)
                                if _coordAtomSite is not None\
                                   and all(_atomId_ in _coordAtomSite['atom_id'] for _atomId_ in _atomId):
                                    compId = _compId = psCompId
                                    resolved = True
                                    break
                            _, _coordAtomSite =\
                                self.getCoordAtomSiteOf(np['auth_chain_id'], npSeqId, npCompId, cifCheck=self.__hasCoord)
                            if self.__mrAtomNameMapping is not None:
                                atomId = retrieveAtomIdFromMRMap(self.ccU, self.__mrAtomNameMapping,
                                                                 npSeqId, npCompId, atomId, _coordAtomSite)
                            _atomId, _, details = self.nefT.get_valid_star_atom_in_xplor(npCompId, atomId, leave_unmatched=True)
                            if details is None:
                                if _coordAtomSite is not None\
                                   and all(_atomId_ in _coordAtomSite['atom_id'] for _atomId_ in _atomId):
                                    compId = _compId = npCompId
                                    seqId = _seqId = npSeqId
                                    preferNonPoly = resolved = True
                                    break

            if not resolved and compId in ('CYSZ', 'CYZ', 'CYS', 'ION', 'ZN1', 'ZN2')\
               and atomId in ZINC_ION_CODE:
                znCount = 0
                znSeqId = None
                for np in self.nonPoly:
                    if np['comp_id'][0] == 'ZN':
                        znSeqId = np['auth_seq_id'][0]
                        znCount += 1
                if znCount > 0:
                    compId = _compId = 'ZN'
                    if znCount == 1:
                        seqId = _seqId = znSeqId
                        atomId = 'ZN'
                        resolved = True
                    preferNonPoly = True

            if not resolved and compId in ('CYS', 'CYSC', 'CYC', 'CCA', 'CYO', 'ION', 'CA1', 'CA2')\
               and atomId in CALCIUM_ION_CODE:
                caCount = 0
                caSeqId = None
                for np in self.nonPoly:
                    if np['comp_id'][0] == 'CA':
                        caSeqId = np['auth_seq_id'][0]
                        caCount += 1
                if caCount > 0:
                    compId = _compId = 'CA'
                    if caCount == 1:
                        seqId = _seqId = caSeqId
                        atomId = 'CA'
                        resolved = True
                    preferNonPoly = True

            if not resolved and len(atomId) > 1 and atomId in ELEMENT_SYMBOLS:
                elemCount = 0
                for np in self.nonPoly:
                    if np['comp_id'][0] == atomId:
                        elemCount += 1
                if elemCount > 0:
                    _, elemSeqId = getMetalCoordOf(self.cR, seqId, compId, atomId)
                    if elemSeqId is not None:
                        seqId = _seqId = elemSeqId
                        compId = _compId = atomId
                        preferNonPoly = resolved = True
                    elif elemCount == 1:
                        for np in self.nonPoly:
                            if np['comp_id'][0] == atomId:
                                seqId = _seqId = np['auth_seq_id'][0]
                                compId = _compId = atomId
                                preferNonPoly = resolved = True

            if not resolved and len(compId) > 1 and compId in ELEMENT_SYMBOLS:
                elemCount = 0
                for np in self.nonPoly:
                    if np['comp_id'][0] == compId:
                        elemCount += 1
                if elemCount > 0:
                    _, elemSeqId = getMetalCoordOf(self.cR, seqId, compId, compId)
                    if elemSeqId is not None:
                        seqId = _seqId = elemSeqId
                        atomId = _compId = compId
                        preferNonPoly = True
                    elif elemCount == 1:
                        for np in self.nonPoly:
                            if np['comp_id'][0] == compId:
                                seqId = _seqId = np['auth_seq_id'][0]
                                atomId = _compId = compId
                                preferNonPoly = True

        if self.__splitLigand is not None and len(self.__splitLigand):
            found = False
            for (_, _seqId_, _compId_), ligList in self.__splitLigand.items():
                if _seqId_ != seqId or _compId_ != compId:
                    continue
                for idx, lig in enumerate(ligList):
                    _atomId = atomId
                    if self.__mrAtomNameMapping is not None and compId not in STD_MON_DICT:
                        _, _, _atomId = retrieveAtomIdentFromMRMap(self.ccU, self.__mrAtomNameMapping,
                                                                   seqId, compId, atomId)

                    if _atomId in lig['atom_ids']:
                        seqId = _seqId = lig['auth_seq_id']
                        compId = _compId = lig['comp_id']
                        atomId = _atomId
                        preferNonPoly = idx > 0
                        found = True
                        break
                if found:
                    break

        if self.__mrAtomNameMapping is not None and compId not in STD_MON_DICT:
            seqId, compId, _ = retrieveAtomIdentFromMRMap(self.ccU, self.__mrAtomNameMapping,
                                                          seqId, compId, atomId)

        compId = self.translateToStdResNameWrapper(_seqId, _compId, preferNonPoly)

        if len(self.__modResidue) > 0:
            modRes = next((modRes for modRes in self.__modResidue
                           if modRes['auth_comp_id'] == compId
                           and (compId != _compId or seqId in (modRes['auth_seq_id'], modRes['seq_id']))), None)
            if modRes is not None:
                compId = modRes['comp_id']

        if self.reasons is not None:
            if 'non_poly_remap' in self.reasons and _compId in self.reasons['non_poly_remap']\
               and seqId in self.reasons['non_poly_remap'][_compId]:
                fixedChainId, fixedSeqId = retrieveRemappedNonPoly(self.reasons['non_poly_remap'], None, refChainId, seqId, _compId)
                refChainId = fixedChainId
                preferNonPoly = True
            if 'branched_remap' in self.reasons and seqId in self.reasons['branched_remap']:
                fixedChainId, fixedSeqId = retrieveRemappedChainId(self.reasons['branched_remap'], seqId)
                refChainId = fixedChainId
                preferNonPoly = True
            if not preferNonPoly:
                if 'chain_id_remap' in self.reasons and seqId in self.reasons['chain_id_remap']:
                    fixedChainId, fixedSeqId = retrieveRemappedChainId(self.reasons['chain_id_remap'], seqId)
                    refChainId = fixedChainId
                elif 'chain_id_clone' in self.reasons and seqId in self.reasons['chain_id_clone']:
                    fixedChainId, fixedSeqId = retrieveRemappedChainId(self.reasons['chain_id_clone'], seqId)
                    refChainId = fixedChainId
                elif 'seq_id_remap' in self.reasons\
                        or 'chain_seq_id_remap' in self.reasons\
                        or 'ext_chain_seq_id_remap' in self.reasons:
                    if 'ext_chain_seq_id_remap' in self.reasons:
                        fixedChainId, fixedSeqId, fixedCompId =\
                            retrieveRemappedSeqIdAndCompId(self.reasons['ext_chain_seq_id_remap'], refChainId, seqId,
                                                           compId if compId in STD_MON_DICT else None)
                        self.__allow_ext_seq = fixedCompId is not None
                        if fixedSeqId is not None:
                            refChainId = fixedChainId
                    if fixedSeqId is None and 'chain_seq_id_remap' in self.reasons:
                        fixedChainId, fixedSeqId = retrieveRemappedSeqId(self.reasons['chain_seq_id_remap'], refChainId, seqId,
                                                                         compId if compId in STD_MON_DICT else None)
                        if fixedSeqId is not None:
                            refChainId = fixedChainId
                    if fixedSeqId is None and 'seq_id_remap' in self.reasons:
                        _, fixedSeqId = retrieveRemappedSeqId(self.reasons['seq_id_remap'], refChainId, seqId)
            if fixedSeqId is not None:
                _seqId = fixedSeqId

        updatePolySeqRst(self.polySeqRst,
                         self.polySeq[0]['chain_id'] if refChainId is None else refChainId, _seqId, compId, _compId)

        types = self.csStat.getTypeOfCompId(compId)
        if all(not t for t in types) or compId in ('MTS', 'ORI'):
            types = None
        elif compId != _compId:
            if types != self.csStat.getTypeOfCompId(_compId):
                types = None

        def comp_id_unmatched_with(ps, cif_comp_id):
            if 'alt_comp_id' in ps and self.csStat.peptideLike(cif_comp_id) and compId.startswith('D') and len(compId) >= 3\
               and self.ccU.lastChemCompDict['type'].upper() == 'D-PEPTIDE LINKING':
                revertPolySeqRst(self.polySeqRst, ps['chain_id'] if fixedChainId is None else fixedChainId, _seqId, compId)

            if types is None or ('alt_comp_id' in ps and _compId in ps['alt_comp_id']):
                return False
            if compId not in STD_MON_DICT and cif_comp_id not in STD_MON_DICT:
                return False
            return types != self.csStat.getTypeOfCompId(cif_comp_id)

        def comp_id_in_polymer(np):
            return (_seqId == 1
                    and ((compId.endswith('-N') and all(c in np['comp_id'][0] for c in compId.split('-')[0]))
                         or (np['comp_id'][0] == 'PCA' and 'P' == compId[0] and ('GL' in compId or 'N' in compId))))\
                or (compId in STD_MON_DICT
                    and any(compId in ps['comp_id'] for ps in self.polySeq)
                    and compId not in np['comp_id'])

        if refChainId is not None or refChainId != _refChainId:
            if any(True for ps in self.polySeq if ps['auth_chain_id'] == _refChainId):
                fixedChainId = _refChainId
            elif self.hasNonPolySeq:
                if any(True for np in self.nonPolySeq if np['auth_chain_id'] == _refChainId):
                    fixedChainId = _refChainId

        for ps in self.polySeq:
            if preferNonPoly:
                continue
            chainId, seqId, cifCompId = self.getRealChainSeqId(ps, _seqId, compId)
            if fixedChainId is None and refChainId is not None and refChainId != chainId and refChainId in self.chainNumberDict:
                if chainId != self.chainNumberDict[refChainId]:
                    continue
            if self.reasons is not None:
                if fixedChainId is not None:
                    if fixedChainId != chainId:
                        continue
                    if fixedSeqId is not None:
                        seqId = fixedSeqId
                elif fixedSeqId is not None:
                    seqId = fixedSeqId
            if seqId <= 0 and self.__shiftNonPosSeq is not None and chainId in self.__shiftNonPosSeq:
                seqId -= 1
            if seqId in ps['auth_seq_id'] or fixedCompId is not None:
                if fixedCompId is not None:
                    cifCompId = origCompId = fixedCompId
                else:
                    if cifCompId is not None:
                        idx = next((_idx for _idx, (_seqId_, _cifCompId_) in enumerate(zip(ps['auth_seq_id'], ps['comp_id']))
                                    if _seqId_ == seqId and _cifCompId_ == cifCompId), ps['auth_seq_id'].index(seqId))
                    else:
                        idx = ps['auth_seq_id'].index(seqId) if seqId in ps['auth_seq_id'] else ps['seq_id'].index(seqId)
                    cifCompId = ps['comp_id'][idx]
                    origCompId = ps['auth_comp_id'][idx]
                    if comp_id_unmatched_with(ps, cifCompId):
                        continue
                if cifCompId != compId:
                    if (self.__shiftNonPosSeq is None or chainId not in self.__shiftNonPosSeq)\
                       and seqId <= 0 and seqId - 1 in ps['auth_seq_id']\
                       and compId == ps['comp_id'][ps['auth_seq_id'].index(seqId - 1)]:
                        seqId -= 1
                        if self.__shiftNonPosSeq is None:
                            self.__shiftNonPosSeq = {}
                        self.__shiftNonPosSeq[chainId] = True
                    compIds = [_compId for _seqId, _compId in zip(ps['auth_seq_id'], ps['comp_id']) if _seqId == seqId]
                    if compId in compIds:
                        cifCompId = compId
                        origCompId = next(origCompId for _seqId, _compId, origCompId
                                          in zip(ps['auth_seq_id'], ps['comp_id'], ps['auth_comp_id'])
                                          if _seqId == seqId and _compId == compId)
                if self.__mrAtomNameMapping is not None and origCompId not in STD_MON_DICT:
                    _, coordAtomSite = self.getCoordAtomSiteOf(chainId, seqId, cifCompId, cifCheck=self.__hasCoord)
                    atomId = retrieveAtomIdFromMRMap(self.ccU, self.__mrAtomNameMapping,
                                                     seqId, origCompId, atomId, coordAtomSite)
                if compId in (cifCompId, origCompId, 'MTS', 'ORI'):
                    if len(self.nefT.get_valid_star_atom(cifCompId, atomId)[0]) > 0:
                        chainAssign.add((chainId, seqId, cifCompId, True))
                        if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                            self.chainNumberDict[refChainId] = chainId
                else:
                    _atomId, _, details = self.nefT.get_valid_star_atom(cifCompId, atomId)
                    if len(_atomId) > 0 and (details is None or _compId not in STD_MON_DICT):
                        chainAssign.add((chainId, seqId, cifCompId, True))
                        if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                            self.chainNumberDict[refChainId] = chainId

            elif 'gap_in_auth_seq' in ps and ps['gap_in_auth_seq']:
                auth_seq_id_list = list(filter(None, ps['auth_seq_id']))
                if len(auth_seq_id_list) > 0:
                    min_auth_seq_id = min(auth_seq_id_list)
                    max_auth_seq_id = max(auth_seq_id_list)
                    if min_auth_seq_id <= seqId <= max_auth_seq_id:
                        _seqId_ = seqId + 1
                        while _seqId_ <= max_auth_seq_id:
                            if _seqId_ in ps['auth_seq_id']:
                                break
                            _seqId_ += 1
                        if _seqId_ not in ps['auth_seq_id']:
                            _seqId_ = seqId - 1
                            while _seqId_ >= min_auth_seq_id:
                                if _seqId_ in ps['auth_seq_id']:
                                    break
                                _seqId_ -= 1
                        if _seqId_ in ps['auth_seq_id']:
                            idx = ps['auth_seq_id'].index(_seqId_) - (_seqId_ - seqId)
                            try:
                                seqId_ = ps['auth_seq_id'][idx]
                                cifCompId = ps['comp_id'][idx]
                                origCompId = ps['auth_comp_id'][idx]
                                if comp_id_unmatched_with(ps, cifCompId):
                                    continue
                                if cifCompId != compId:
                                    compIds = [_compId for _seqId, _compId in zip(ps['auth_seq_id'], ps['comp_id'])
                                               if _seqId == seqId]
                                    if compId in compIds:
                                        cifCompId = compId
                                        origCompId = next(origCompId for _seqId, _compId, origCompId
                                                          in zip(ps['auth_seq_id'], ps['comp_id'], ps['auth_comp_id'])
                                                          if _seqId == seqId and _compId == compId)
                                if self.__mrAtomNameMapping is not None and origCompId not in STD_MON_DICT:
                                    _, coordAtomSite = self.getCoordAtomSiteOf(chainId, seqId_, cifCompId, cifCheck=self.__hasCoord)
                                    atomId = retrieveAtomIdFromMRMap(self.ccU, self.__mrAtomNameMapping,
                                                                     seqId, origCompId, atomId, coordAtomSite)
                                if compId in (cifCompId, origCompId, 'MTS', 'ORI'):
                                    if len(self.nefT.get_valid_star_atom(cifCompId, atomId)[0]) > 0:
                                        chainAssign.add((chainId, seqId_, cifCompId, True))
                                        if refChainId is not None and refChainId != chainId\
                                           and refChainId not in self.chainNumberDict:
                                            self.chainNumberDict[refChainId] = chainId
                                else:
                                    _atomId, _, details = self.nefT.get_valid_star_atom(cifCompId, atomId)
                                    if len(_atomId) > 0 and (details is None or _compId not in STD_MON_DICT):
                                        chainAssign.add((chainId, seqId_, cifCompId, True))
                                        if refChainId is not None and refChainId != chainId\
                                           and refChainId not in self.chainNumberDict:
                                            self.chainNumberDict[refChainId] = chainId
                            except IndexError:
                                pass

        if self.hasNonPolySeq:
            ligands = 0
            if self.hasNonPoly:
                for np in self.nonPoly:
                    ligands += np['comp_id'].count(_compId)
                if ligands == 0:
                    for np in self.nonPoly:
                        ligands += np['comp_id'].count(compId)
                    if ligands == 1:
                        _compId = compId
                if ligands == 0:
                    for np in self.nonPoly:
                        if 'alt_comp_id' in np:
                            ligands += np['alt_comp_id'].count(_compId)
                if ligands == 0:
                    for np in self.nonPoly:
                        if 'alt_comp_id' in np:
                            ligands += np['alt_comp_id'].count(compId)
                    if ligands == 1:
                        _compId = compId
                if ligands == 0 and len(chainAssign) == 0:
                    __compId = None
                    for np in self.nonPoly:
                        for ligand in np['comp_id']:
                            __compId = translateToLigandName(_compId, ligand, self.ccU)
                            if __compId == ligand:
                                ligands += 1
                    if ligands == 1:
                        compId = _compId = __compId
                    elif len(self.nonPoly) == 1 and self.ccU.updateChemCompDict(_compId, False):
                        if self.ccU.lastChemCompDict['release_status'] == 'OBS':
                            compId = _compId = self.nonPoly[0]['comp_id'][0]
                            ligands = 1
                if self.reasons is None and atomId in self.__uniqAtomIdToSeqKey:
                    seqKey = self.__uniqAtomIdToSeqKey[atomId]
                    if _seqId != seqKey[1]:
                        r = self.__getNamedReasonsForReparsing('non_poly_remap')
                        if _compId not in r:
                            r[_compId] = {}
                        if _seqId not in r[_compId]:
                            r[_compId][_seqId] = {'chain_id': seqKey[0], 'seq_id': seqKey[1], 'original_chain_id': None}
            for np in self.nonPolySeq:
                chainId, seqId, cifCompId = self.getRealChainSeqId(np, _seqId, compId, False)
                if fixedChainId is None and refChainId is not None and refChainId != chainId and refChainId in self.chainNumberDict:
                    if chainId != self.chainNumberDict[refChainId]:
                        continue
                if self.reasons is not None:
                    if fixedChainId is not None:
                        if fixedChainId != chainId:
                            continue
                        if fixedSeqId is not None:
                            seqId = fixedSeqId
                    elif fixedSeqId is not None:
                        seqId = fixedSeqId
                if comp_id_in_polymer(np):
                    continue
                if 'alt_auth_seq_id' in np and seqId not in np['auth_seq_id'] and seqId in np['alt_auth_seq_id']:
                    try:
                        seqId = next(_seqId_ for _seqId_, _altSeqId_ in zip(np['auth_seq_id'], np['alt_auth_seq_id'])
                                     if _altSeqId_ == seqId)
                    except StopIteration:
                        pass
                if seqId in np['auth_seq_id']\
                   or (ligands == 1 and (_compId in np['comp_id'] or ('alt_comp_id' in np and _compId in np['alt_comp_id']))):
                    if ligands == 1 and cifCompId is None:
                        cifCompId = _compId
                    idx = -1
                    try:
                        if cifCompId is not None:
                            idx = next(_idx for _idx, (_seqId_, _cifCompId_) in enumerate(zip(np['auth_seq_id'], np['comp_id']))
                                       if (_seqId_ == seqId or ligands == 1) and _cifCompId_ == cifCompId)
                            if ligands == 1:
                                seqId = np['auth_seq_id'][idx]
                    except StopIteration:
                        pass
                    if idx == -1:
                        idx = np['auth_seq_id'].index(seqId) if seqId in np['auth_seq_id']\
                            else np['seq_id'].index(seqId) if seqId in np['seq_id'] else 0
                    cifCompId = np['comp_id'][idx]
                    origCompId = np['auth_comp_id'][idx]
                    seqId = np['auth_seq_id'][idx]
                    if cifCompId in ('ZN', 'CA') and atomId[0] in PROTON_BEGIN_CODE:  # 2loa
                        continue
                    if self.__mrAtomNameMapping is not None and origCompId not in STD_MON_DICT:
                        _, coordAtomSite = self.getCoordAtomSiteOf(chainId, seqId, cifCompId, cifCheck=self.__hasCoord)
                        atomId = retrieveAtomIdFromMRMap(self.ccU, self.__mrAtomNameMapping,
                                                         _seqId, origCompId, atomId, coordAtomSite)
                    if compId in (cifCompId, origCompId, 'MTS', 'ORI'):
                        if len(self.nefT.get_valid_star_atom(cifCompId, atomId)[0]) > 0:
                            if (ligands == 1 or cifCompId in HEME_LIKE_RES_NAMES) and any(a[3] for a in chainAssign):
                                chainAssign.clear()
                            chainAssign.add((chainId, seqId, cifCompId, False))
                            if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                                self.chainNumberDict[refChainId] = chainId
                    else:
                        _atomId, _, details = self.nefT.get_valid_star_atom(cifCompId, atomId)
                        if len(_atomId) > 0 and (details is None or _compId not in STD_MON_DICT):
                            if (ligands == 1 or cifCompId in HEME_LIKE_RES_NAMES) and any(a[3] for a in chainAssign):
                                chainAssign.clear()
                            chainAssign.add((chainId, seqId, cifCompId, False))
                            if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                                self.chainNumberDict[refChainId] = chainId

        if len(chainAssign) == 0:
            for ps in self.polySeq:
                if preferNonPoly:
                    continue
                chainId = ps['chain_id']
                if refChainId is not None and refChainId != chainId and refChainId in self.chainNumberDict:
                    if chainId != self.chainNumberDict[refChainId]:
                        continue
                if fixedChainId is not None and fixedChainId != chainId:
                    continue
                seqKey = (chainId, _seqId)
                if seqKey in self.__authToLabelSeq:
                    _, seqId = self.__authToLabelSeq[seqKey]
                    if seqId in ps['seq_id']:
                        idx = ps['seq_id'].index(seqId)
                        cifCompId = ps['comp_id'][idx]
                        origCompId = ps['auth_comp_id'][idx]
                        if comp_id_unmatched_with(ps, cifCompId):
                            continue
                        if cifCompId != compId:
                            compIds = [_compId for _seqId, _compId in zip(ps['auth_seq_id'], ps['comp_id']) if _seqId == seqId]
                            if compId in compIds:
                                cifCompId = compId
                                origCompId =\
                                    next(origCompId
                                         for _seqId, _compId, origCompId
                                         in zip(ps['auth_seq_id'], ps['comp_id'], ps['auth_comp_id'])
                                         if _seqId == seqId and _compId == compId)
                        if self.__mrAtomNameMapping is not None and origCompId not in STD_MON_DICT:
                            _, coordAtomSite = self.getCoordAtomSiteOf(chainId, _seqId, cifCompId, cifCheck=self.__hasCoord)
                            atomId = retrieveAtomIdFromMRMap(self.ccU, self.__mrAtomNameMapping,
                                                             seqId, origCompId, atomId, coordAtomSite)
                        if compId in (cifCompId, origCompId, 'MTS', 'ORI'):
                            if len(self.nefT.get_valid_star_atom(cifCompId, atomId)[0]) > 0:
                                chainAssign.add((ps['auth_chain_id'], _seqId, cifCompId, True))
                                if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                                    self.chainNumberDict[refChainId] = chainId
                        else:
                            _atomId, _, details = self.nefT.get_valid_star_atom(cifCompId, atomId)
                            if len(_atomId) > 0 and (details is None or _compId not in STD_MON_DICT):
                                chainAssign.add((ps['auth_chain_id'], _seqId, cifCompId, True))
                                if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                                    self.chainNumberDict[refChainId] = chainId

            if self.hasNonPolySeq:
                for np in self.nonPolySeq:
                    chainId = np['auth_chain_id']
                    if refChainId is not None and refChainId != chainId and refChainId in self.chainNumberDict:
                        if chainId != self.chainNumberDict[refChainId]:
                            continue
                    if fixedChainId is not None and fixedChainId != chainId:
                        continue
                    if comp_id_in_polymer(np):
                        continue
                    seqKey = (chainId, _seqId)
                    if seqKey in self.__authToLabelSeq:
                        _, seqId = self.__authToLabelSeq[seqKey]
                        if seqId in np['seq_id']:
                            idx = np['seq_id'].index(seqId)
                            cifCompId = np['comp_id'][idx]
                            origCompId = np['auth_comp_id'][idx]
                            if self.__mrAtomNameMapping is not None and origCompId not in STD_MON_DICT:
                                _, coordAtomSite = self.getCoordAtomSiteOf(chainId, _seqId, cifCompId, cifCheck=self.__hasCoord)
                                atomId = retrieveAtomIdFromMRMap(self.ccU, self.__mrAtomNameMapping,
                                                                 seqId, origCompId, atomId, coordAtomSite)
                            if compId in (cifCompId, origCompId, 'MTS', 'ORI'):
                                if len(self.nefT.get_valid_star_atom(cifCompId, atomId)[0]) > 0:
                                    chainAssign.add((np['auth_chain_id'], _seqId, cifCompId, False))
                                    if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                                        self.chainNumberDict[refChainId] = chainId
                            else:
                                _atomId, _, details = self.nefT.get_valid_star_atom(cifCompId, atomId)
                                if len(_atomId) > 0 and (details is None or _compId not in STD_MON_DICT):
                                    chainAssign.add((np['auth_chain_id'], _seqId, cifCompId, False))
                                    if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                                        self.chainNumberDict[refChainId] = chainId

        if len(chainAssign) == 0 and self.__altPolySeq is not None:
            for ps in self.__altPolySeq:
                if preferNonPoly:
                    continue
                chainId = ps['auth_chain_id']
                if refChainId is not None and refChainId != chainId and refChainId in self.chainNumberDict:
                    if chainId != self.chainNumberDict[refChainId]:
                        continue
                if fixedChainId is not None and fixedChainId != chainId:
                    continue
                if _seqId in ps['auth_seq_id']:
                    cifCompId = ps['comp_id'][ps['auth_seq_id'].index(_seqId)]
                    if comp_id_unmatched_with(ps, cifCompId):
                        continue
                    if cifCompId != compId:
                        if cifCompId in STD_MON_DICT and compId in STD_MON_DICT:
                            continue
                        compIds = [_compId for _seqId, _compId in zip(ps['auth_seq_id'], ps['comp_id']) if _seqId == seqId]
                        if compId in compIds:
                            cifCompId = compId
                    chainAssign.add((chainId, _seqId, cifCompId, True))
                    if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                        self.chainNumberDict[refChainId] = chainId

        if len(chainAssign) == 0 and (self.__preferAuthSeqCount - self.__preferLabelSeqCount < MAX_PREF_LABEL_SCHEME_COUNT
                                      or len(self.polySeq) > 1):
            for ps in self.polySeq:
                if preferNonPoly:
                    continue
                chainId = ps['chain_id']
                if refChainId is not None and refChainId != chainId and refChainId in self.chainNumberDict:
                    if chainId != self.chainNumberDict[refChainId]:
                        continue
                if fixedChainId is not None and fixedChainId != chainId:
                    continue
                seqKey = (chainId, _seqId)
                if seqKey in self.__labelToAuthSeq:
                    _, seqId = self.__labelToAuthSeq[seqKey]
                    if seqId in ps['auth_seq_id']:
                        idx = ps['seq_id'].index(_seqId)
                        cifCompId = ps['comp_id'][idx]
                        origCompId = ps['auth_comp_id'][idx]
                        if comp_id_unmatched_with(ps, cifCompId):
                            continue
                        if cifCompId != compId:
                            compIds = [_compId for _seqId, _compId in zip(ps['auth_seq_id'], ps['comp_id']) if _seqId == seqId]
                            if compId in compIds:
                                cifCompId = compId
                                origCompId =\
                                    next(origCompId
                                         for _seqId, _compId, origCompId
                                         in zip(ps['auth_seq_id'], ps['comp_id'], ps['auth_comp_id'])
                                         if _seqId == seqId and _compId == compId)
                        if self.__mrAtomNameMapping is not None and origCompId not in STD_MON_DICT:
                            _, coordAtomSite = self.getCoordAtomSiteOf(chainId, seqId, cifCompId, cifCheck=self.__hasCoord)
                            atomId = retrieveAtomIdFromMRMap(self.ccU, self.__mrAtomNameMapping,
                                                             seqId, origCompId, atomId, coordAtomSite)
                        if compId in (cifCompId, origCompId, 'MTS', 'ORI'):
                            if len(self.nefT.get_valid_star_atom(cifCompId, atomId)[0]) > 0:
                                chainAssign.add((ps['auth_chain_id'], seqId, cifCompId, True))
                                if compId in (cifCompId, origCompId):
                                    self.__authSeqId = 'label_seq_id'
                                    self.__setLocalSeqScheme()
                                    if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                                        self.chainNumberDict[refChainId] = chainId
                        else:
                            _atomId, _, details = self.nefT.get_valid_star_atom(cifCompId, atomId)
                            if len(_atomId) > 0 and (details is None or _compId not in STD_MON_DICT):
                                chainAssign.add((ps['auth_chain_id'], seqId, cifCompId, True))
                                self.__authSeqId = 'label_seq_id'
                                self.__setLocalSeqScheme()
                                if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                                    self.chainNumberDict[refChainId] = chainId

        if len(chainAssign) == 0:
            if seqId == 1 or (refChainId, seqId - 1) in self.__coordUnobsRes:
                if atomId in AMINO_PROTON_CODE and atomId != 'H1':
                    return self.assignCoordPolymerSequence(refChainId, seqId, compId, 'H1', index)
            auth_seq_id_list = list(filter(None, self.polySeq[0]['auth_seq_id']))
            min_auth_seq_id = max_auth_seq_id = UNREAL_AUTH_SEQ_NUM
            if len(auth_seq_id_list) > 0:
                min_auth_seq_id = min(auth_seq_id_list)
                max_auth_seq_id = max(auth_seq_id_list)
            if len(self.polySeq) == 1\
               and (seqId < 1
                    or (compId == 'ACE' and seqId == min_auth_seq_id - 1)
                    or (compId == 'NH2' and seqId == max_auth_seq_id + 1)
                    or (compId in STD_MON_DICT
                        and self.__preferAuthSeqCount - self.__preferLabelSeqCount >= MAX_PREF_LABEL_SCHEME_COUNT)):
                refChainId = self.polySeq[0]['auth_chain_id']
                if (compId == 'ACE' and seqId == min_auth_seq_id - 1)\
                   or (compId == 'NH2' and seqId == max_auth_seq_id + 1)\
                   or (compId in STD_MON_DICT
                       and self.__preferAuthSeqCount - self.__preferLabelSeqCount >= MAX_PREF_LABEL_SCHEME_COUNT
                       and (min_auth_seq_id - MAX_ALLOWED_EXT_SEQ <= seqId < min_auth_seq_id
                            or max_auth_seq_id < seqId <= max_auth_seq_id + MAX_ALLOWED_EXT_SEQ)):
                    self.f.append(f"[Sequence mismatch warning] {self.getCurrentSpectralPeak(n=index)}"
                                  f"The residue '{_seqId}:{_compId}' is not present in polymer sequence "
                                  f"of chain {refChainId} of the coordinates. {INSTRUCTION_FOR_FULL_SEQUENCE}")
                    resKey = (_seqId, _compId)
                    if resKey not in self.extResKey:
                        self.extResKey.append(resKey)
                    chainAssign.add((refChainId, _seqId, compId, True))
                    asis = True
                elif compId in STD_MON_DICT\
                        and self.__preferAuthSeqCount - self.__preferLabelSeqCount >= MAX_PREF_LABEL_SCHEME_COUNT:
                    self.f.append(f"[Sequence mismatch warning] {self.getCurrentSpectralPeak(n=index)}"
                                  f"The residue '{_seqId}:{_compId}' is not present in polymer sequence "
                                  f"of chain {refChainId} of the coordinates. {INSTRUCTION_FOR_FULL_SEQUENCE}")
                    resKey = (_seqId, _compId)
                    if resKey not in self.extResKey:
                        self.extResKey.append(resKey)
                elif self.no_extra_comment:
                    self.f.append(f"[Atom not found] {self.getCurrentSpectralPeak(n=index)}"
                                  f"{_seqId}:{_compId}:{atomId} is not present in the coordinates. "
                                  f"The residue number '{_seqId}' is not present in polymer sequence "
                                  f"of chain {refChainId} of the coordinates. {INSTRUCTION_FOR_FULL_SEQUENCE}")
            else:
                ext_seq = False
                if (compId in STD_MON_DICT or compId in ('ACE', 'NH2'))\
                   and self.__preferAuthSeqCount - self.__preferLabelSeqCount >= MAX_PREF_LABEL_SCHEME_COUNT:
                    refChainIds = []
                    _auth_seq_id_list = auth_seq_id_list
                    for idx, ps in enumerate(self.polySeq):
                        if idx > 0:
                            auth_seq_id_list = list(filter(None, ps['auth_seq_id']))
                            _auth_seq_id_list.extend(auth_seq_id_list)
                        if len(auth_seq_id_list) > 0:
                            if idx > 0:
                                min_auth_seq_id = min(auth_seq_id_list)
                                max_auth_seq_id = max(auth_seq_id_list)
                            if min_auth_seq_id - MAX_ALLOWED_EXT_SEQ <= seqId < min_auth_seq_id\
                               and (compId in STD_MON_DICT or (compId == 'ACE' and seqId == min_auth_seq_id - 1)):
                                refChainIds.append(ps['auth_chain_id'])
                                ext_seq = True
                            if max_auth_seq_id < seqId <= max_auth_seq_id + MAX_ALLOWED_EXT_SEQ\
                               and (compId in STD_MON_DICT or (compId == 'NH2' and seqId == max_auth_seq_id + 1)):
                                refChainIds.append(ps['auth_chain_id'])
                                ext_seq = True
                    if ext_seq and seqId in _auth_seq_id_list:
                        ext_seq = False
                if ext_seq:
                    refChainId = refChainIds[0] if len(refChainIds) == 1 else refChainIds
                    self.f.append(f"[Sequence mismatch warning] {self.getCurrentSpectralPeak(n=index)}"
                                  f"The residue '{_seqId}:{_compId}' is not present in polymer sequence "
                                  f"of chain {refChainId} of the coordinates. {INSTRUCTION_FOR_FULL_SEQUENCE}")
                    resKey = (_seqId, _compId)
                    if resKey not in self.extResKey:
                        self.extResKey.append(resKey)
                    if isinstance(refChainId, str):
                        chainAssign.add((refChainId, _seqId, compId, True))
                    else:
                        for _refChainId in refChainIds:
                            chainAssign.add((_refChainId, _seqId, compId, True))
                    asis = True
                elif self.no_extra_comment:
                    self.f.append(f"[Atom not found] {self.getCurrentSpectralPeak(n=index)}"
                                  f"{_seqId}:{_compId}:{atomId} is not present in the coordinates.")
                updatePolySeqRst(self.polySeqRstFailed,
                                 self.polySeq[0]['chain_id'] if refChainId is None else refChainId, _seqId, compId, _compId)

        elif any(True for ca in chainAssign if ca[0] == refChainId) and any(True for ca in chainAssign if ca[0] != refChainId):
            _chainAssign = copy.copy(chainAssign)
            for _ca in _chainAssign:
                if _ca[0] != refChainId:
                    chainAssign.remove(_ca)

        return list(chainAssign), asis

    def assignCoordPolymerSequenceWithChainId(self, refChainId: str, seqId: int, compId: str, atomId: str, index: int
                                              ) -> Tuple[List[Tuple[str, int, str, bool]], bool]:
        """ Assign polymer sequences of the coordinates.
        """

        _refChainId = refChainId

        chainAssign = set()
        asis = False
        _seqId = seqId
        _compId = compId

        fixedChainId = fixedSeqId = fixedCompId = None

        preferNonPoly = False

        if self.hasNonPoly:

            resolved = False

            for np in self.nonPoly:
                if 'alt_comp_id' in np and 'alt_auth_seq_id' in np\
                   and compId in np['alt_comp_id'] and seqId in np['alt_auth_seq_id']:
                    npCompId = np['comp_id'][0]
                    npSeqId = np['auth_seq_id'][0]
                    for ps in self.polySeq:
                        if 'ambig_auth_seq_id' in ps and seqId in ps['ambig_auth_seq_id']:
                            psCompId = ps['comp_id'][ps['auth_seq_id'].index(seqId)]
                            _atomId, _, details = self.nefT.get_valid_star_atom_in_xplor(psCompId, atomId, leave_unmatched=True)
                            if details is None:
                                _, _coordAtomSite =\
                                    self.getCoordAtomSiteOf(ps['auth_chain_id'], seqId, psCompId, cifCheck=self.__hasCoord)
                                if _coordAtomSite is not None\
                                   and all(_atomId_ in _coordAtomSite['atom_id'] for _atomId_ in _atomId):
                                    compId = _compId = psCompId
                                    resolved = True
                                    break
                            _, _coordAtomSite =\
                                self.getCoordAtomSiteOf(np['auth_chain_id'], npSeqId, npCompId, cifCheck=self.__hasCoord)
                            if self.__mrAtomNameMapping is not None:
                                atomId = retrieveAtomIdFromMRMap(self.ccU, self.__mrAtomNameMapping,
                                                                 npSeqId, npCompId, atomId, _coordAtomSite)
                            _atomId, _, details = self.nefT.get_valid_star_atom_in_xplor(npCompId, atomId, leave_unmatched=True)
                            if details is None:
                                if _coordAtomSite is not None\
                                   and all(_atomId_ in _coordAtomSite['atom_id'] for _atomId_ in _atomId):
                                    compId = _compId = npCompId
                                    seqId = _seqId = npSeqId
                                    preferNonPoly = resolved = True
                                    break

            if not resolved and compId in ('CYS', 'CYSZ', 'CYZ', 'CZN', 'CYO', 'ION', 'ZN1', 'ZN2')\
               and atomId in ZINC_ION_CODE:
                znCount = 0
                znSeqId = None
                for np in self.nonPoly:
                    if np['comp_id'][0] == 'ZN':
                        znSeqId = np['auth_seq_id'][0]
                        znCount += 1
                if znCount > 0:
                    compId = _compId = 'ZN'
                    if znCount == 1:
                        seqId = _seqId = znSeqId
                        atomId = 'ZN'
                        resolved = True
                    preferNonPoly = True

            if not resolved and compId in ('CYS', 'CYSC', 'CYC', 'CCA', 'CYO', 'ION', 'CA1', 'CA2')\
               and atomId in CALCIUM_ION_CODE:
                caCount = 0
                caSeqId = None
                for np in self.nonPoly:
                    if np['comp_id'][0] == 'CA':
                        caSeqId = np['auth_seq_id'][0]
                        caCount += 1
                if caCount > 0:
                    compId = _compId = 'CA'
                    if caCount == 1:
                        seqId = _seqId = caSeqId
                        atomId = 'CA'
                        resolved = True
                    preferNonPoly = True

            if not resolved and len(atomId) > 1 and atomId in ELEMENT_SYMBOLS:
                elemCount = 0
                for np in self.nonPoly:
                    if np['comp_id'][0] == atomId:
                        elemCount += 1
                if elemCount > 0:
                    _, elemSeqId = getMetalCoordOf(self.cR, seqId, compId, atomId)
                    if elemSeqId is not None:
                        seqId = _seqId = elemSeqId
                        compId = _compId = atomId
                        preferNonPoly = resolved = True
                    elif elemCount == 1:
                        for np in self.nonPoly:
                            if np['comp_id'][0] == atomId:
                                seqId = _seqId = np['auth_seq_id'][0]
                                compId = _compId = atomId
                                preferNonPoly = resolved = True

            if not resolved and len(compId) > 1 and compId in ELEMENT_SYMBOLS:
                elemCount = 0
                for np in self.nonPoly:
                    if np['comp_id'][0] == compId:
                        elemCount += 1
                if elemCount > 0:
                    _, elemSeqId = getMetalCoordOf(self.cR, seqId, compId, compId)
                    if elemSeqId is not None:
                        seqId = _seqId = elemSeqId
                        atomId = _compId = compId
                        preferNonPoly = True
                    elif elemCount == 1:
                        for np in self.nonPoly:
                            if np['comp_id'][0] == compId:
                                seqId = _seqId = np['auth_seq_id'][0]
                                atomId = _compId = compId
                                preferNonPoly = True

        if self.__splitLigand is not None and len(self.__splitLigand):
            found = False
            for (_, _seqId_, _compId_), ligList in self.__splitLigand.items():
                if _seqId_ != seqId or _compId_ != compId:
                    continue
                for idx, lig in enumerate(ligList):
                    _atomId = atomId
                    if self.__mrAtomNameMapping is not None and compId not in STD_MON_DICT:
                        _, _, _atomId = retrieveAtomIdentFromMRMap(self.ccU, self.__mrAtomNameMapping,
                                                                   seqId, compId, atomId)

                    if _atomId in lig['atom_ids']:
                        seqId = _seqId = lig['auth_seq_id']
                        compId = _compId = lig['comp_id']
                        atomId = _atomId
                        preferNonPoly = idx > 0
                        found = True
                        break
                if found:
                    break

        if self.__mrAtomNameMapping is not None and compId not in STD_MON_DICT:
            seqId, compId, _ = retrieveAtomIdentFromMRMap(self.ccU, self.__mrAtomNameMapping,
                                                          seqId, compId, atomId)

        compId = self.translateToStdResNameWrapper(_seqId, _compId, preferNonPoly)

        if len(self.__modResidue) > 0:
            modRes = next((modRes for modRes in self.__modResidue
                           if modRes['auth_comp_id'] == compId
                           and (compId != _compId or seqId in (modRes['auth_seq_id'], modRes['seq_id']))), None)
            if modRes is not None:
                compId = modRes['comp_id']

        self.__allow_ext_seq = False

        if self.reasons is not None:
            if 'unambig_atom_id_remap' in self.reasons and _compId in self.reasons['unambig_atom_id_remap']\
               and atomId in self.reasons['unambig_atom_id_remap'][_compId]:
                atomId = self.reasons['unambig_atom_id_remap'][_compId][atomId][0]  # select representative one
            if 'non_poly_remap' in self.reasons and _compId in self.reasons['non_poly_remap']\
               and seqId in self.reasons['non_poly_remap'][_compId]:
                fixedChainId, fixedSeqId = retrieveRemappedNonPoly(self.reasons['non_poly_remap'], None,
                                                                   str(refChainId), seqId, _compId)
                refChainId = fixedChainId
                preferNonPoly = True
            if 'branched_remap' in self.reasons and seqId in self.reasons['branched_remap']:
                fixedChainId, fixedSeqId = retrieveRemappedChainId(self.reasons['branched_remap'], seqId)
                refChainId = fixedChainId
                preferNonPoly = True
            if not preferNonPoly:
                if 'chain_id_remap' in self.reasons and seqId in self.reasons['chain_id_remap']:
                    fixedChainId, fixedSeqId = retrieveRemappedChainId(self.reasons['chain_id_remap'], seqId)
                    refChainId = fixedChainId
                elif 'chain_id_clone' in self.reasons and seqId in self.reasons['chain_id_clone']:
                    fixedChainId, fixedSeqId = retrieveRemappedChainId(self.reasons['chain_id_clone'], seqId)
                    refChainId = fixedChainId
                elif 'seq_id_remap' in self.reasons\
                        or 'chain_seq_id_remap' in self.reasons\
                        or 'ext_chain_seq_id_remap' in self.reasons:
                    if 'ext_chain_seq_id_remap' in self.reasons:
                        fixedChainId, fixedSeqId, fixedCompId =\
                            retrieveRemappedSeqIdAndCompId(self.reasons['ext_chain_seq_id_remap'], str(refChainId), seqId,
                                                           compId if compId in STD_MON_DICT else None)
                        self.__allow_ext_seq = fixedCompId is not None
                        if fixedSeqId is not None:
                            refChainId = fixedChainId
                    if fixedSeqId is None and 'chain_seq_id_remap' in self.reasons:
                        fixedChainId, fixedSeqId = retrieveRemappedSeqId(self.reasons['chain_seq_id_remap'], str(refChainId), seqId,
                                                                         compId if compId in STD_MON_DICT else None)
                        if fixedSeqId is not None:
                            refChainId = fixedChainId
                    if fixedSeqId is None and 'seq_id_remap' in self.reasons:
                        _, fixedSeqId = retrieveRemappedSeqId(self.reasons['seq_id_remap'], str(refChainId), seqId)
            if fixedSeqId is not None:
                _seqId = fixedSeqId

        updatePolySeqRst(self.polySeqRst, str(refChainId), _seqId, compId, _compId)

        types = self.csStat.getTypeOfCompId(compId)
        if all(not t for t in types) or compId in ('MTS', 'ORI'):
            types = None
        elif compId != _compId:
            if types != self.csStat.getTypeOfCompId(_compId):
                types = None

        def comp_id_unmatched_with(ps, cif_comp_id):
            if 'alt_comp_id' in ps and self.csStat.peptideLike(cif_comp_id) and compId.startswith('D') and len(compId) >= 3\
               and self.ccU.lastChemCompDict['type'].upper() == 'D-PEPTIDE LINKING':
                revertPolySeqRst(self.polySeqRst, str(refChainId), _seqId, compId)

            if types is None or ('alt_comp_id' in ps and _compId in ps['alt_comp_id']):
                return False
            if compId not in STD_MON_DICT and cif_comp_id not in STD_MON_DICT:
                return False
            return types != self.csStat.getTypeOfCompId(cif_comp_id)

        def comp_id_in_polymer(np):
            return (_seqId == 1
                    and ((compId.endswith('-N') and all(c in np['comp_id'][0] for c in compId.split('-')[0]))
                         or (np['comp_id'][0] == 'PCA' and 'P' == compId[0] and ('GL' in compId or 'N' in compId))))\
                or (compId in STD_MON_DICT
                    and any(compId in ps['comp_id'] for ps in self.polySeq)
                    and compId not in np['comp_id'])

        if refChainId is not None or refChainId != _refChainId:
            if any(True for ps in self.polySeq if ps['auth_chain_id'] == _refChainId):
                fixedChainId = _refChainId
            elif self.hasNonPolySeq:
                if any(True for np in self.nonPolySeq if np['auth_chain_id'] == _refChainId):
                    fixedChainId = _refChainId

        for ps in self.polySeq:
            if preferNonPoly:
                continue
            chainId, seqId, cifCompId = self.getRealChainSeqId(ps, _seqId, compId)
            if fixedChainId is None and refChainId is not None and refChainId != chainId and refChainId in self.chainNumberDict:
                if chainId != self.chainNumberDict[refChainId]:
                    continue
            if self.reasons is not None:
                if fixedChainId is not None:
                    if fixedChainId != chainId:
                        continue
                    if fixedSeqId is not None:
                        seqId = fixedSeqId
                elif fixedSeqId is not None:
                    seqId = fixedSeqId
            if seqId <= 0 and self.__shiftNonPosSeq is not None and chainId in self.__shiftNonPosSeq:
                seqId -= 1
            if seqId in ps['auth_seq_id'] or fixedCompId is not None:
                if fixedCompId is not None:
                    cifCompId = origCompId = fixedCompId
                else:
                    if cifCompId is not None:
                        idx = next((_idx for _idx, (_seqId_, _cifCompId_) in enumerate(zip(ps['auth_seq_id'], ps['comp_id']))
                                    if _seqId_ == seqId and _cifCompId_ == cifCompId), ps['auth_seq_id'].index(seqId))
                    else:
                        idx = ps['auth_seq_id'].index(seqId) if seqId in ps['auth_seq_id'] else ps['seq_id'].index(seqId)
                    cifCompId = ps['comp_id'][idx]
                    origCompId = ps['auth_comp_id'][idx]
                    if comp_id_unmatched_with(ps, cifCompId):
                        continue
                if cifCompId != compId:
                    if (self.__shiftNonPosSeq is None or chainId not in self.__shiftNonPosSeq)\
                       and seqId <= 0 and seqId - 1 in ps['auth_seq_id']\
                       and compId == ps['comp_id'][ps['auth_seq_id'].index(seqId - 1)]:
                        seqId -= 1
                        if self.__shiftNonPosSeq is None:
                            self.__shiftNonPosSeq = {}
                        self.__shiftNonPosSeq[chainId] = True
                    compIds = [_compId for _seqId, _compId in zip(ps['auth_seq_id'], ps['comp_id']) if _seqId == seqId]
                    if compId in compIds:
                        cifCompId = compId
                        origCompId = next(origCompId for _seqId, _compId, origCompId
                                          in zip(ps['auth_seq_id'], ps['comp_id'], ps['auth_comp_id'])
                                          if _seqId == seqId and _compId == compId)
                if self.__mrAtomNameMapping is not None and origCompId not in STD_MON_DICT:
                    _, coordAtomSite = self.getCoordAtomSiteOf(chainId, seqId, cifCompId, cifCheck=self.__hasCoord)
                    atomId = retrieveAtomIdFromMRMap(self.ccU, self.__mrAtomNameMapping,
                                                     seqId, origCompId, atomId, coordAtomSite)
                if compId in (cifCompId, origCompId, 'MTS', 'ORI'):
                    if len(self.nefT.get_valid_star_atom(cifCompId, atomId)[0]) > 0:
                        chainAssign.add((chainId, seqId, cifCompId, True))
                        if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                            self.chainNumberDict[refChainId] = chainId
                else:
                    _atomId, _, details = self.nefT.get_valid_star_atom(cifCompId, atomId)
                    if len(_atomId) > 0 and (details is None or _compId not in STD_MON_DICT):
                        chainAssign.add((chainId, seqId, cifCompId, True))
                        if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                            self.chainNumberDict[refChainId] = chainId

            elif 'gap_in_auth_seq' in ps and ps['gap_in_auth_seq']:
                auth_seq_id_list = list(filter(None, ps['auth_seq_id']))
                if len(auth_seq_id_list) > 0:
                    min_auth_seq_id = min(auth_seq_id_list)
                    max_auth_seq_id = max(auth_seq_id_list)
                    if min_auth_seq_id <= seqId <= max_auth_seq_id:
                        _seqId_ = seqId + 1
                        while _seqId_ <= max_auth_seq_id:
                            if _seqId_ in ps['auth_seq_id']:
                                break
                            _seqId_ += 1
                        if _seqId_ not in ps['auth_seq_id']:
                            _seqId_ = seqId - 1
                            while _seqId_ >= min_auth_seq_id:
                                if _seqId_ in ps['auth_seq_id']:
                                    break
                                _seqId_ -= 1
                        if _seqId_ in ps['auth_seq_id']:
                            idx = ps['auth_seq_id'].index(_seqId_) - (_seqId_ - seqId)
                            try:
                                seqId_ = ps['auth_seq_id'][idx]
                                cifCompId = ps['comp_id'][idx]
                                origCompId = ps['auth_comp_id'][idx]
                                if comp_id_unmatched_with(ps, cifCompId):
                                    continue
                                if cifCompId != compId:
                                    compIds = [_compId for _seqId, _compId in zip(ps['auth_seq_id'], ps['comp_id'])
                                               if _seqId == seqId]
                                    if compId in compIds:
                                        cifCompId = compId
                                        origCompId = next(origCompId for _seqId, _compId, origCompId
                                                          in zip(ps['auth_seq_id'], ps['comp_id'], ps['auth_comp_id'])
                                                          if _seqId == seqId and _compId == compId)
                                if self.__mrAtomNameMapping is not None and origCompId not in STD_MON_DICT:
                                    _, coordAtomSite = self.getCoordAtomSiteOf(chainId, seqId_, cifCompId, cifCheck=self.__hasCoord)
                                    atomId = retrieveAtomIdFromMRMap(self.ccU, self.__mrAtomNameMapping,
                                                                     seqId, origCompId, atomId, coordAtomSite)
                                if compId in (cifCompId, origCompId, 'MTS', 'ORI'):
                                    if len(self.nefT.get_valid_star_atom(cifCompId, atomId)[0]) > 0:
                                        chainAssign.add((chainId, seqId_, cifCompId, True))
                                    else:
                                        _atomId, _, details = self.nefT.get_valid_star_atom(cifCompId, atomId)
                                        if len(_atomId) > 0 and (details is None or _compId not in STD_MON_DICT):
                                            chainAssign.add((chainId, seqId_, cifCompId, True))
                            except IndexError:
                                pass

        if self.hasNonPolySeq:
            ligands = 0
            if self.hasNonPoly:
                for np in self.nonPoly:
                    ligands += np['comp_id'].count(_compId)
                if ligands == 0:
                    for np in self.nonPoly:
                        ligands += np['comp_id'].count(compId)
                    if ligands == 1:
                        _compId = compId
                if ligands == 0:
                    for np in self.nonPoly:
                        if 'alt_comp_id' in np:
                            ligands += np['alt_comp_id'].count(_compId)
                if ligands == 0:
                    for np in self.nonPoly:
                        if 'alt_comp_id' in np:
                            ligands += np['alt_comp_id'].count(compId)
                    if ligands == 1:
                        _compId = compId
                if ligands == 0 and len(chainAssign) == 0:
                    __compId = None
                    for np in self.nonPoly:
                        for ligand in np['comp_id']:
                            __compId = translateToLigandName(_compId, ligand, self.ccU)
                            if __compId == ligand:
                                ligands += 1
                    if ligands == 1:
                        compId = _compId = __compId
                    elif len(self.nonPoly) == 1 and self.ccU.updateChemCompDict(_compId, False):
                        if self.ccU.lastChemCompDict['release_status'] == 'OBS':
                            compId = _compId = self.nonPoly[0]['comp_id'][0]
                            ligands = 1
                if self.reasons is None and atomId in self.__uniqAtomIdToSeqKey:
                    seqKey = self.__uniqAtomIdToSeqKey[atomId]
                    if _seqId != seqKey[1]:
                        r = self.__getNamedReasonsForReparsing('non_poly_remap')
                        if _compId not in r:
                            r[_compId] = {}
                        if _seqId not in r[_compId]:
                            r[_compId][_seqId] = {'chain_id': seqKey[0], 'seq_id': seqKey[1], 'original_chain_id': refChainId}
            for np in self.nonPolySeq:
                chainId, seqId, cifCompId = self.getRealChainSeqId(np, _seqId, compId, False)
                if fixedChainId is None and refChainId is not None and refChainId != chainId and refChainId in self.chainNumberDict:
                    if chainId != self.chainNumberDict[refChainId]:
                        continue
                if self.reasons is not None:
                    if fixedChainId is not None:
                        if fixedChainId != chainId:
                            continue
                        if fixedSeqId is not None:
                            seqId = fixedSeqId
                    elif fixedSeqId is not None:
                        seqId = fixedSeqId
                if comp_id_in_polymer(np):
                    continue
                if 'alt_auth_seq_id' in np and seqId not in np['auth_seq_id'] and seqId in np['alt_auth_seq_id']:
                    try:
                        seqId = next(_seqId_ for _seqId_, _altSeqId_ in zip(np['auth_seq_id'], np['alt_auth_seq_id'])
                                     if _altSeqId_ == seqId)
                    except StopIteration:
                        pass
                if seqId in np['auth_seq_id']\
                   or (ligands == 1 and (_compId in np['comp_id'] or ('alt_comp_id' in np and _compId in np['alt_comp_id']))):
                    if ligands == 1 and cifCompId is None:
                        cifCompId = _compId
                    idx = -1
                    try:
                        if cifCompId is not None:
                            idx = next(_idx for _idx, (_seqId_, _cifCompId_) in enumerate(zip(np['auth_seq_id'], np['comp_id']))
                                       if (_seqId_ == seqId or ligands == 1) and _cifCompId_ == cifCompId)
                            if ligands == 1:
                                seqId = np['auth_seq_id'][idx]
                    except StopIteration:
                        pass
                    if idx == -1:
                        idx = np['auth_seq_id'].index(seqId) if seqId in np['auth_seq_id']\
                            else np['seq_id'].index(seqId) if seqId in np['seq_id'] else 0
                    cifCompId = np['comp_id'][idx]
                    origCompId = np['auth_comp_id'][idx]
                    seqId = np['auth_seq_id'][idx]
                    if cifCompId in ('ZN', 'CA') and atomId[0] in PROTON_BEGIN_CODE:  # 2loa
                        continue
                    if self.__mrAtomNameMapping is not None and origCompId not in STD_MON_DICT:
                        _, coordAtomSite = self.getCoordAtomSiteOf(chainId, seqId, cifCompId, cifCheck=self.__hasCoord)
                        atomId = retrieveAtomIdFromMRMap(self.ccU, self.__mrAtomNameMapping,
                                                         _seqId, origCompId, atomId, coordAtomSite)
                    if compId in (cifCompId, origCompId, 'MTS', 'ORI'):
                        if len(self.nefT.get_valid_star_atom(cifCompId, atomId)[0]) > 0:
                            if (ligands == 1 or cifCompId in HEME_LIKE_RES_NAMES) and any(a[3] for a in chainAssign):
                                chainAssign.clear()
                            chainAssign.add((chainId, seqId, cifCompId, False))
                            if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                                self.chainNumberDict[refChainId] = chainId
                    else:
                        _atomId, _, details = self.nefT.get_valid_star_atom(cifCompId, atomId)
                        if len(_atomId) > 0 and (details is None or _compId not in STD_MON_DICT):
                            if (ligands == 1 or cifCompId in HEME_LIKE_RES_NAMES) and any(a[3] for a in chainAssign):
                                chainAssign.clear()
                            chainAssign.add((chainId, seqId, cifCompId, False))
                            if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                                self.chainNumberDict[refChainId] = chainId

        if len(chainAssign) == 0:
            for ps in self.polySeq:
                if preferNonPoly:
                    continue
                chainId = ps['chain_id']
                if refChainId is not None and refChainId != chainId and refChainId in self.chainNumberDict:
                    if chainId != self.chainNumberDict[refChainId]:
                        continue
                if fixedChainId is not None and fixedChainId != chainId:
                    continue
                seqKey = (chainId, _seqId)
                if seqKey in self.__authToLabelSeq:
                    _, seqId = self.__authToLabelSeq[seqKey]
                    if seqId in ps['seq_id']:
                        idx = ps['seq_id'].index(seqId)
                        cifCompId = ps['comp_id'][idx]
                        origCompId = ps['auth_comp_id'][idx]
                        if comp_id_unmatched_with(ps, cifCompId):
                            continue
                        if cifCompId != compId:
                            compIds = [_compId for _seqId, _compId in zip(ps['auth_seq_id'], ps['comp_id']) if _seqId == seqId]
                            if compId in compIds:
                                cifCompId = compId
                                origCompId =\
                                    next(origCompId
                                         for _seqId, _compId, origCompId
                                         in zip(ps['auth_seq_id'], ps['comp_id'], ps['auth_comp_id'])
                                         if _seqId == seqId and _compId == compId)
                        if self.__mrAtomNameMapping is not None and origCompId not in STD_MON_DICT:
                            _, coordAtomSite = self.getCoordAtomSiteOf(chainId, _seqId, cifCompId, cifCheck=self.__hasCoord)
                            atomId = retrieveAtomIdFromMRMap(self.ccU, self.__mrAtomNameMapping,
                                                             seqId, origCompId, atomId, coordAtomSite)
                        if compId in (cifCompId, origCompId, 'MTS', 'ORI'):
                            if len(self.nefT.get_valid_star_atom(cifCompId, atomId)[0]) > 0:
                                chainAssign.add((ps['auth_chain_id'], _seqId, cifCompId, True))
                                if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                                    self.chainNumberDict[refChainId] = chainId
                        else:
                            _atomId, _, details = self.nefT.get_valid_star_atom(cifCompId, atomId)
                            if len(_atomId) > 0 and (details is None or _compId not in STD_MON_DICT):
                                chainAssign.add((ps['auth_chain_id'], _seqId, cifCompId, True))
                                if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                                    self.chainNumberDict[refChainId] = chainId

            if self.hasNonPolySeq:
                for np in self.nonPolySeq:
                    chainId = np['auth_chain_id']
                    if refChainId is not None and refChainId != chainId and refChainId in self.chainNumberDict:
                        if chainId != self.chainNumberDict[refChainId]:
                            continue
                    if fixedChainId is not None and fixedChainId != chainId:
                        continue
                    if comp_id_in_polymer(np):
                        continue
                    seqKey = (chainId, _seqId)
                    if seqKey in self.__authToLabelSeq:
                        _, seqId = self.__authToLabelSeq[seqKey]
                        if seqId in np['seq_id']:
                            idx = np['seq_id'].index(seqId)
                            cifCompId = np['comp_id'][idx]
                            origCompId = np['auth_comp_id'][idx]
                            if self.__mrAtomNameMapping is not None and origCompId not in STD_MON_DICT:
                                _, coordAtomSite = self.getCoordAtomSiteOf(chainId, _seqId, cifCompId, cifCheck=self.__hasCoord)
                                atomId = retrieveAtomIdFromMRMap(self.ccU, self.__mrAtomNameMapping,
                                                                 seqId, origCompId, atomId, coordAtomSite)
                            if compId in (cifCompId, origCompId, 'MTS', 'ORI'):
                                if len(self.nefT.get_valid_star_atom(cifCompId, atomId)[0]) > 0:
                                    chainAssign.add((np['auth_chain_id'], _seqId, cifCompId, False))
                                    if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                                        self.chainNumberDict[refChainId] = chainId
                            else:
                                _atomId, _, details = self.nefT.get_valid_star_atom(cifCompId, atomId)
                                if len(_atomId) > 0 and (details is None or _compId not in STD_MON_DICT):
                                    chainAssign.add((np['auth_chain_id'], _seqId, cifCompId, False))
                                    if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                                        self.chainNumberDict[refChainId] = chainId

        if len(chainAssign) == 0 and self.__altPolySeq is not None:
            for ps in self.__altPolySeq:
                if preferNonPoly:
                    continue
                chainId = ps['auth_chain_id']
                if fixedChainId is None and refChainId is not None and refChainId != chainId and refChainId in self.chainNumberDict:
                    if chainId != self.chainNumberDict[refChainId]:
                        continue
                if fixedChainId is not None:
                    if fixedChainId != chainId:
                        continue
                    if fixedSeqId is not None:
                        _seqId = fixedSeqId
                elif fixedSeqId is not None:
                    _seqId = fixedSeqId
                if _seqId in ps['auth_seq_id']:
                    cifCompId = ps['comp_id'][ps['auth_seq_id'].index(_seqId)]
                    if comp_id_unmatched_with(ps, cifCompId):
                        continue
                    if cifCompId != compId:
                        if cifCompId in STD_MON_DICT and compId in STD_MON_DICT:
                            continue
                        compIds = [_compId for _seqId, _compId in zip(ps['auth_seq_id'], ps['comp_id']) if _seqId == seqId]
                        if compId in compIds:
                            cifCompId = compId
                    chainAssign.add((chainId, _seqId, cifCompId, True))
                    if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                        self.chainNumberDict[refChainId] = chainId

        if len(chainAssign) == 0 and (self.__preferAuthSeqCount - self.__preferLabelSeqCount < MAX_PREF_LABEL_SCHEME_COUNT
                                      or len(self.polySeq) > 1):
            for ps in self.polySeq:
                if preferNonPoly:
                    continue
                chainId = ps['chain_id']
                if refChainId is not None and refChainId != chainId and refChainId in self.chainNumberDict:
                    if chainId != self.chainNumberDict[refChainId]:
                        continue
                if fixedChainId is not None and fixedChainId != chainId:
                    continue
                seqKey = (chainId, _seqId)
                if seqKey in self.__labelToAuthSeq:
                    _, seqId = self.__labelToAuthSeq[seqKey]
                    if seqId in ps['auth_seq_id']:
                        idx = ps['seq_id'].index(_seqId)
                        cifCompId = ps['comp_id'][idx]
                        origCompId = ps['auth_comp_id'][idx]
                        if comp_id_unmatched_with(ps, cifCompId):
                            continue
                        if cifCompId != compId:
                            compIds = [_compId for _seqId, _compId in zip(ps['auth_seq_id'], ps['comp_id']) if _seqId == seqId]
                            if compId in compIds:
                                cifCompId = compId
                                origCompId =\
                                    next(origCompId
                                         for _seqId, _compId, origCompId
                                         in zip(ps['auth_seq_id'], ps['comp_id'], ps['auth_comp_id'])
                                         if _seqId == seqId and _compId == compId)
                        if self.__mrAtomNameMapping is not None and origCompId not in STD_MON_DICT:
                            _, coordAtomSite = self.getCoordAtomSiteOf(chainId, seqId, cifCompId, cifCheck=self.__hasCoord)
                            atomId = retrieveAtomIdFromMRMap(self.ccU, self.__mrAtomNameMapping,
                                                             seqId, origCompId, atomId, coordAtomSite)
                        if compId in (cifCompId, origCompId, 'MTS', 'ORI'):
                            if len(self.nefT.get_valid_star_atom(cifCompId, atomId)[0]) > 0:
                                chainAssign.add((ps['auth_chain_id'], seqId, cifCompId, True))
                                if compId in (cifCompId, origCompId):
                                    self.__authSeqId = 'label_seq_id'
                                    self.__setLocalSeqScheme()
                                    if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                                        self.chainNumberDict[refChainId] = chainId
                        else:
                            _atomId, _, details = self.nefT.get_valid_star_atom(cifCompId, atomId)
                            if len(_atomId) > 0 and (details is None or _compId not in STD_MON_DICT):
                                chainAssign.add((ps['auth_chain_id'], seqId, cifCompId, True))
                                self.__authSeqId = 'label_seq_id'
                                self.__setLocalSeqScheme()
                                if refChainId is not None and refChainId != chainId and refChainId not in self.chainNumberDict:
                                    self.chainNumberDict[refChainId] = chainId

        if len(chainAssign) == 0:
            if seqId == 1 or (refChainId, seqId - 1) in self.__coordUnobsRes:
                if atomId in AMINO_PROTON_CODE and atomId != 'H1':
                    return self.assignCoordPolymerSequenceWithChainId(refChainId, seqId, compId, 'H1', index)
            else:
                auth_seq_id_list = list(filter(None, self.polySeq[0]['auth_seq_id']))
                min_auth_seq_id = max_auth_seq_id = UNREAL_AUTH_SEQ_NUM
                if len(auth_seq_id_list) > 0:
                    min_auth_seq_id = min(auth_seq_id_list)
                    max_auth_seq_id = max(auth_seq_id_list)
                if len(self.polySeq) == 1\
                   and (seqId < 1
                        or (compId == 'ACE' and seqId == min_auth_seq_id - 1)
                        or (compId == 'NH2' and seqId == max_auth_seq_id + 1)
                        or (compId in STD_MON_DICT
                            and self.__preferAuthSeqCount - self.__preferLabelSeqCount >= MAX_PREF_LABEL_SCHEME_COUNT)):
                    refChainId = self.polySeq[0]['auth_chain_id']
                    if (compId == 'ACE' and seqId == min_auth_seq_id - 1)\
                       or (compId == 'NH2' and seqId == max_auth_seq_id + 1)\
                       or (compId in STD_MON_DICT
                           and self.__preferAuthSeqCount - self.__preferLabelSeqCount >= MAX_PREF_LABEL_SCHEME_COUNT
                           and (min_auth_seq_id - MAX_ALLOWED_EXT_SEQ <= seqId < min_auth_seq_id
                                or max_auth_seq_id < seqId <= max_auth_seq_id + MAX_ALLOWED_EXT_SEQ)):
                        self.f.append(f"[Sequence mismatch warning] {self.getCurrentSpectralPeak(n=index)}"
                                      f"The residue '{_seqId}:{_compId}' is not present in polymer sequence "
                                      f"of chain {refChainId} of the coordinates. {INSTRUCTION_FOR_FULL_SEQUENCE}")
                        resKey = (_seqId, _compId)
                        if resKey not in self.extResKey:
                            self.extResKey.append(resKey)
                        chainAssign.add((refChainId, _seqId, compId, True))
                        asis = True
                    elif compId in STD_MON_DICT\
                            and self.__preferAuthSeqCount - self.__preferLabelSeqCount >= MAX_PREF_LABEL_SCHEME_COUNT:
                        self.f.append(f"[Sequence mismatch warning] {self.getCurrentSpectralPeak(n=index)}"
                                      f"The residue '{_seqId}:{_compId}' is not present in polymer sequence "
                                      f"of chain {refChainId} of the coordinates. {INSTRUCTION_FOR_FULL_SEQUENCE}")
                        resKey = (_seqId, _compId)
                        if resKey not in self.extResKey:
                            self.extResKey.append(resKey)
                    elif self.no_extra_comment:
                        self.f.append(f"[Atom not found] {self.getCurrentSpectralPeak(n=index)}"
                                      f"{_seqId}:{_compId}:{atomId} is not present in the coordinates. "
                                      f"The residue number '{_seqId}' is not present in polymer sequence "
                                      f"of chain {refChainId} of the coordinates. {INSTRUCTION_FOR_FULL_SEQUENCE}")
                else:
                    ext_seq = False
                    if (compId in STD_MON_DICT or compId in ('ACE', 'NH2'))\
                       and self.__preferAuthSeqCount - self.__preferLabelSeqCount >= MAX_PREF_LABEL_SCHEME_COUNT:
                        refChainIds = []
                        _auth_seq_id_list = auth_seq_id_list
                        for idx, ps in enumerate(self.polySeq):
                            if idx > 0:
                                auth_seq_id_list = list(filter(None, ps['auth_seq_id']))
                                _auth_seq_id_list.extend(auth_seq_id_list)
                            if len(auth_seq_id_list) > 0:
                                if idx > 0:
                                    min_auth_seq_id = min(auth_seq_id_list)
                                    max_auth_seq_id = max(auth_seq_id_list)
                                if min_auth_seq_id - MAX_ALLOWED_EXT_SEQ <= seqId < min_auth_seq_id\
                                   and (compId in STD_MON_DICT or (compId == 'ACE' and seqId == min_auth_seq_id - 1)):
                                    refChainIds.append(ps['auth_chain_id'])
                                    ext_seq = True
                                if max_auth_seq_id < seqId <= max_auth_seq_id + MAX_ALLOWED_EXT_SEQ\
                                   and (compId in STD_MON_DICT or (compId == 'NH2' and seqId == max_auth_seq_id + 1)):
                                    refChainIds.append(ps['auth_chain_id'])
                                    ext_seq = True
                        if ext_seq and seqId in _auth_seq_id_list:
                            ext_seq = False
                    if ext_seq:
                        refChainId = refChainIds[0] if len(refChainIds) == 1 else refChainIds
                        self.f.append(f"[Sequence mismatch warning] {self.getCurrentSpectralPeak(n=index)}"
                                      f"The residue '{_seqId}:{_compId}' is not present in polymer sequence "
                                      f"of chain {refChainId} of the coordinates. {INSTRUCTION_FOR_FULL_SEQUENCE}")
                        resKey = (_seqId, _compId)
                        if resKey not in self.extResKey:
                            self.extResKey.append(resKey)
                        if isinstance(refChainId, str):
                            chainAssign.add((refChainId, _seqId, compId, True))
                        else:
                            for _refChainId in refChainIds:
                                chainAssign.add((_refChainId, _seqId, compId, True))
                        asis = True
                    elif self.no_extra_comment:
                        self.f.append(f"[Atom not found] {self.getCurrentSpectralPeak(n=index)}"
                                      f"{_seqId}:{_compId}:{atomId} is not present in the coordinates.")
                    updatePolySeqRst(self.polySeqRstFailed, str(refChainId), _seqId, compId, _compId)

        elif any(True for ca in chainAssign if ca[0] == refChainId) and any(True for ca in chainAssign if ca[0] != refChainId):
            _chainAssign = copy.copy(chainAssign)
            for _ca in _chainAssign:
                if _ca[0] != refChainId:
                    chainAssign.remove(_ca)

        return list(chainAssign), asis

    def assignCoordPolymerSequenceWithoutCompId(self, seqId: int, atomId: str, index: int, fixedChainId: Optional[str] = None
                                                ) -> List[Tuple[str, int, str, bool]]:
        """ Assign polymer sequences of the coordinates, restricted to a given chain (fixedChainId) if any.
        """

        chainAssign = set()
        _seqId = seqId

        _refChainId = fixedChainId
        fixedSeqId = fixedCompId = None

        self.__allow_ext_seq = False

        if self.reasons is not None:
            if 'branched_remap' in self.reasons and seqId in self.reasons['branched_remap']:
                fixedChainId, fixedSeqId = retrieveRemappedChainId(self.reasons['branched_remap'], seqId)
            if 'chain_id_remap' in self.reasons and seqId in self.reasons['chain_id_remap']:
                fixedChainId, fixedSeqId = retrieveRemappedChainId(self.reasons['chain_id_remap'], seqId)
            elif 'chain_id_clone' in self.reasons and seqId in self.reasons['chain_id_clone']:
                fixedChainId, fixedSeqId = retrieveRemappedChainId(self.reasons['chain_id_clone'], seqId)
            if fixedChainId is None:
                fixedChainId = _refChainId
            if fixedSeqId is not None:
                seqId = _seqId = fixedSeqId

        for ps in self.polySeq:
            chainId, seqId, cifCompId = self.getRealChainSeqId(ps, _seqId, None)
            if fixedChainId is not None and chainId != fixedChainId:
                continue
            if self.reasons is not None:
                if 'seq_id_remap' not in self.reasons\
                   and 'chain_seq_id_remap' not in self.reasons\
                   and 'ext_chain_seq_id_remap' not in self.reasons:
                    if fixedChainId is not None and fixedChainId != chainId:
                        continue
                else:
                    if 'ext_chain_seq_id_remap' in self.reasons:
                        remapChainId, fixedSeqId, fixedCompId =\
                            retrieveRemappedSeqIdAndCompId(self.reasons['ext_chain_seq_id_remap'], chainId, seqId)
                        if remapChainId is not None and remapChainId != chainId:
                            continue
                        if fixedSeqId is not None:
                            self.__allow_ext_seq = fixedCompId is not None
                            seqId = _seqId = fixedSeqId
                    if fixedSeqId is None and 'chain_seq_id_remap' in self.reasons:
                        remapChainId, fixedSeqId = retrieveRemappedSeqId(self.reasons['chain_seq_id_remap'], chainId, seqId)
                        if remapChainId is not None and remapChainId != chainId:
                            continue
                        if fixedSeqId is not None:
                            seqId = _seqId = fixedSeqId
                    if fixedSeqId is None and 'seq_id_remap' in self.reasons:
                        _, fixedSeqId = retrieveRemappedSeqId(self.reasons['seq_id_remap'], chainId, seqId)
                        if fixedSeqId is not None:
                            seqId = _seqId = fixedSeqId
            if seqId in ps['auth_seq_id'] or fixedCompId is not None:
                if fixedCompId is not None:
                    cifCompId = fixedCompId
                else:
                    if cifCompId is not None:
                        idx = next((_idx for _idx, (_seqId_, _cifCompId_) in enumerate(zip(ps['auth_seq_id'], ps['comp_id']))
                                    if _seqId_ == seqId and _cifCompId_ == cifCompId), ps['auth_seq_id'].index(seqId))
                    else:
                        idx = ps['auth_seq_id'].index(seqId) if seqId in ps['auth_seq_id'] else ps['seq_id'].index(seqId)
                    cifCompId = ps['comp_id'][idx]
                if self.reasons is not None:
                    if 'non_poly_remap' in self.reasons and cifCompId in self.reasons['non_poly_remap']\
                       and seqId in self.reasons['non_poly_remap'][cifCompId]:
                        remapChainId, fixedSeqId = retrieveRemappedNonPoly(self.reasons['non_poly_remap'], None,
                                                                           chainId, seqId, cifCompId)
                        if fixedSeqId is not None:
                            seqId = _seqId = fixedSeqId
                        if (remapChainId is not None and remapChainId != chainId) or seqId not in ps['auth_seq_id']:
                            continue
                updatePolySeqRst(self.polySeqRst, chainId, _seqId, cifCompId)
                if atomId is None or len(self.nefT.get_valid_star_atom(cifCompId, atomId)[0]) > 0:
                    chainAssign.add((chainId, seqId, cifCompId, True))
            elif 'gap_in_auth_seq' in ps and ps['gap_in_auth_seq']:
                auth_seq_id_list = list(filter(None, ps['auth_seq_id']))
                if len(auth_seq_id_list) > 0:
                    min_auth_seq_id = min(auth_seq_id_list)
                    max_auth_seq_id = max(auth_seq_id_list)
                    if min_auth_seq_id <= seqId <= max_auth_seq_id:
                        _seqId_ = seqId + 1
                        while _seqId_ <= max_auth_seq_id:
                            if _seqId_ in ps['auth_seq_id']:
                                break
                            _seqId_ += 1
                        if _seqId_ not in ps['auth_seq_id']:
                            _seqId_ = seqId - 1
                            while _seqId_ >= min_auth_seq_id:
                                if _seqId_ in ps['auth_seq_id']:
                                    break
                                _seqId_ -= 1
                        if _seqId_ in ps['auth_seq_id']:
                            idx = ps['auth_seq_id'].index(_seqId_) - (_seqId_ - seqId)
                            try:
                                seqId_ = ps['auth_seq_id'][idx]
                                cifCompId = ps['comp_id'][idx]
                                updatePolySeqRst(self.polySeqRst, chainId, _seqId, cifCompId)
                                if atomId is None or len(self.nefT.get_valid_star_atom(cifCompId, atomId)[0]) > 0:
                                    chainAssign.add((chainId, seqId_, cifCompId, True))
                            except IndexError:
                                pass

        if self.hasNonPolySeq:
            for np in self.nonPolySeq:
                chainId, seqId, cifCompId = self.getRealChainSeqId(np, _seqId, None, False)
                if fixedChainId is not None and chainId != fixedChainId:
                    continue
                if self.reasons is not None:
                    if 'seq_id_remap' not in self.reasons and 'chain_seq_id_remap' not in self.reasons:
                        if fixedChainId is not None and fixedChainId != chainId:
                            continue
                    else:
                        if 'chain_seq_id_remap' in self.reasons:
                            remapChainId, fixedSeqId = retrieveRemappedSeqId(self.reasons['chain_seq_id_remap'], chainId, seqId)
                            if remapChainId is not None and remapChainId != chainId:
                                continue
                            if fixedSeqId is not None:
                                seqId = _seqId = fixedSeqId
                        if fixedSeqId is None and 'seq_id_remap' in self.reasons:
                            _, fixedSeqId = retrieveRemappedSeqId(self.reasons['seq_id_remap'], chainId, seqId)
                            if fixedSeqId is not None:
                                seqId = _seqId = fixedSeqId
                if seqId in np['auth_seq_id']:
                    if cifCompId is not None:
                        idx = next((_idx for _idx, (_seqId_, _cifCompId_) in enumerate(zip(np['auth_seq_id'], np['comp_id']))
                                    if _seqId_ == seqId and _cifCompId_ == cifCompId), np['auth_seq_id'].index(seqId))
                    else:
                        idx = np['auth_seq_id'].index(seqId) if seqId in np['auth_seq_id'] else np['seq_id'].index(seqId)
                    cifCompId = np['comp_id'][idx]
                    updatePolySeqRst(self.polySeqRst, chainId, _seqId, cifCompId)
                    if atomId is None or len(self.nefT.get_valid_star_atom(cifCompId, atomId)[0]) > 0:
                        chainAssign.add((chainId, seqId, cifCompId, False))

        if len(chainAssign) == 0:
            for ps in self.polySeq:
                chainId = ps['chain_id']
                if fixedChainId is not None and fixedChainId != chainId:
                    continue
                seqKey = (chainId, _seqId)
                if seqKey in self.__authToLabelSeq:
                    _, seqId = self.__authToLabelSeq[seqKey]
                    if seqId in ps['seq_id']:
                        cifCompId = ps['comp_id'][ps['seq_id'].index(seqId)]
                        updatePolySeqRst(self.polySeqRst, chainId, _seqId, cifCompId)
                        if atomId is None or len(self.nefT.get_valid_star_atom(cifCompId, atomId)[0]) > 0:
                            chainAssign.add((ps['auth_chain_id'], _seqId, cifCompId, True))

            if self.hasNonPolySeq:
                for np in self.nonPolySeq:
                    chainId = np['auth_chain_id']
                    if fixedChainId is not None and fixedChainId != chainId:
                        continue
                    seqKey = (chainId, _seqId)
                    if seqKey in self.__authToLabelSeq:
                        _, seqId = self.__authToLabelSeq[seqKey]
                        if seqId in np['seq_id']:
                            cifCompId = np['comp_id'][np['seq_id'].index(seqId)]
                            updatePolySeqRst(self.polySeqRst, chainId, _seqId, cifCompId)
                            if atomId is None or len(self.nefT.get_valid_star_atom(cifCompId, atomId)[0]) > 0:
                                chainAssign.add((np['auth_chain_id'], _seqId, cifCompId, False))

        if len(chainAssign) == 0 and self.__altPolySeq is not None:
            for ps in self.__altPolySeq:
                chainId = ps['auth_chain_id']
                if fixedChainId is not None and fixedChainId != chainId:
                    continue
                if _seqId in ps['auth_seq_id']:
                    cifCompId = ps['comp_id'][ps['auth_seq_id'].index(_seqId)]
                    updatePolySeqRst(self.polySeqRst, chainId, _seqId, cifCompId)
                    chainAssign.add((chainId, _seqId, cifCompId, True))

        if len(chainAssign) == 0 and (self.__preferAuthSeqCount - self.__preferLabelSeqCount < MAX_PREF_LABEL_SCHEME_COUNT
                                      or len(self.polySeq) > 1):
            for ps in self.polySeq:
                chainId = ps['chain_id']
                if fixedChainId is not None and fixedChainId != chainId:
                    continue
                seqKey = (chainId, _seqId)
                if seqKey in self.__labelToAuthSeq:
                    _, seqId = self.__labelToAuthSeq[seqKey]
                    if seqId in ps['auth_seq_id']:
                        cifCompId = ps['comp_id'][ps['auth_seq_id'].index(seqId)]
                        updatePolySeqRst(self.polySeqRst, chainId, seqId, cifCompId)
                        if atomId is None or len(self.nefT.get_valid_star_atom(cifCompId, atomId)[0]) > 0:
                            chainAssign.add((ps['auth_chain_id'], seqId, cifCompId, True))
                            self.__authSeqId = 'label_seq_id'
                            self.__setLocalSeqScheme()

        if len(chainAssign) == 0:
            _chainId_ = '' if _refChainId is None else f'{fixedChainId}:'  # prefix of the restraint's chain if given
            if seqId == 1 or (chainId if fixedChainId is None else fixedChainId, seqId - 1) in self.__coordUnobsRes:
                if atomId is not None and atomId in AMINO_PROTON_CODE and atomId != 'H1':
                    return self.assignCoordPolymerSequenceWithoutCompId(seqId, 'H1', index,
                                                                        None if _refChainId is None else fixedChainId)
            if atomId is not None and (('-' in atomId and ':' in atomId) or '.' in atomId):
                if self.no_extra_comment:
                    self.f.append(f"[Atom not found] {self.getCurrentSpectralPeak(n=index)}"
                                  f"{_chainId_}{_seqId}:?:{atomId} is not present in the coordinates.")
            elif atomId is not None:
                if len(self.polySeq) == 1 and seqId < 1:
                    refChainId = self.polySeq[0]['auth_chain_id']
                    if self.no_extra_comment:
                        self.f.append(f"[Atom not found] {self.getCurrentSpectralPeak(n=index)}"
                                      f"{_seqId}:?:{atomId} is not present in the coordinates. "
                                      f"The residue number '{_seqId}' is not present in polymer sequence "
                                      f"of chain {refChainId} of the coordinates. {INSTRUCTION_FOR_FULL_SEQUENCE}")
                else:
                    if self.no_extra_comment:
                        self.f.append(f"[Atom not found] {self.getCurrentSpectralPeak(n=index)}"
                                      f"{_chainId_}{_seqId}:{atomId} is not present in the coordinates.")
                    compIds = guessCompIdFromAtomId([atomId], self.polySeq, self.nefT)
                    if compIds is not None:
                        chainId = fixedChainId
                        if chainId is None and len(self.polySeq) == 1:
                            chainId = self.polySeq[0]['chain_id']
                        if chainId is not None:
                            if len(compIds) == 1:
                                updatePolySeqRst(self.polySeqRstFailed, chainId, seqId, compIds[0])
                            else:
                                updatePolySeqRstAmbig(self.polySeqRstFailedAmbig, chainId, seqId, compIds)

        return list(chainAssign)

    def assignCoordPolymerSequenceWithChainIdWithoutCompId(self, fixedChainId: str, seqId: int, atomId: str, index: int
                                                           ) -> List[Tuple[str, int, str, bool]]:
        """ Assign polymer sequences of the coordinates of a given chain.
            @return: no assignment without a chain, so that callers fall back to assignCoordPolymerSequenceWithoutCompId()
        """

        if fixedChainId is None:
            return []

        return self.assignCoordPolymerSequenceWithoutCompId(seqId, atomId, index, fixedChainId)

    def selectCoordAtoms(self, chainAssign: List[Tuple[str, int, str, bool]], authChainId: Optional[str],
                         seqId: int, compId: str, atomId: str,
                         index: int, allowAmbig: bool = True, offset: int = 0) -> None:
        """ Select atoms of the coordinates.
        """

        atomSelection = []

        authAtomId = atomId

        __compId = compId
        __atomId = atomId

        if self.__mrAtomNameMapping is not None and compId not in STD_MON_DICT:
            __atomId = retrieveAtomIdFromMRMap(self.ccU, self.__mrAtomNameMapping, seqId, compId, atomId)

        compId = self.translateToStdResNameWrapper(seqId, __compId)

        for chainId, cifSeqId, cifCompId, isPolySeq in chainAssign:

            if offset != 0:
                cifSeqId += offset
                cifCompId = compId

            seqKey, coordAtomSite = self.getCoordAtomSiteOf(chainId, cifSeqId, cifCompId, asis=self.__preferAuthSeq)

            if self.__mrAtomNameMapping is not None and cifCompId not in STD_MON_DICT:
                __atomId = retrieveAtomIdFromMRMap(self.ccU, self.__mrAtomNameMapping,
                                                   cifSeqId, cifCompId, atomId, coordAtomSite)
                if atomId != __atomId and coordAtomSite is not None\
                   and (__atomId in coordAtomSite['atom_id']
                        or (__atomId.endswith('%') and f'{__atomId[:-1]}2' in coordAtomSite['atom_id'])):
                    atomId = __atomId
                elif self.reasons is not None and 'branched_remap' in self.reasons:
                    _seqId = retrieveOriginalSeqIdFromMRMap(self.reasons['branched_remap'], chainId, cifSeqId)
                    if _seqId != cifSeqId:
                        _, _, atomId = retrieveAtomIdentFromMRMap(self.ccU, self.__mrAtomNameMapping,
                                                                  _seqId, cifCompId, atomId, None, coordAtomSite)

            _atomId = []
            if not isPolySeq and atomId[0] in ('Q', 'M') and coordAtomSite is not None:
                key = (chainId, cifSeqId, compId, atomId)
                if key in self.__cachedDictForStarAtom:
                    _atomId = list(self.__cachedDictForStarAtom[key])
                else:
                    pattern = re.compile(fr'H{atomId[1:]}\d+') if compId in STD_MON_DICT else re.compile(fr'H{atomId[1:]}\S?$')
                    atomIdList = [a for a in coordAtomSite['atom_id'] if re.search(pattern, a) and a[-1] in ('1', '2', '3')]
                    if len(atomIdList) > 1:
                        hvyAtomIdList = [a for a in coordAtomSite['atom_id'] if a[0] in ('C', 'N')]
                        hvyAtomId = None
                        for canHvyAtomId in hvyAtomIdList:
                            if isStructConn(self.cR, chainId, cifSeqId, canHvyAtomId, chainId, cifSeqId, atomIdList[0],
                                            representativeModelId=self.representativeModelId,
                                            representativeAltId=self.representativeAltId,
                                            modelNumName=self.modelNumName):
                                hvyAtomId = canHvyAtomId
                                break
                        if hvyAtomId is not None:
                            for _atomId_ in atomIdList:
                                if isStructConn(self.cR, chainId, cifSeqId, hvyAtomId, chainId, cifSeqId, _atomId_,
                                                representativeModelId=self.representativeModelId,
                                                representativeAltId=self.representativeAltId,
                                                modelNumName=self.modelNumName):
                                    _atomId.append(_atomId_)
                    if len(_atomId) > 1:
                        self.__cachedDictForStarAtom[key] = list(_atomId)
            if len(_atomId) > 1:
                details = None
            else:
                _atomId, _, details = self.nefT.get_valid_star_atom_in_xplor(cifCompId, atomId, leave_unmatched=True)
                if details is not None:
                    if atomId != __atomId:
                        _atomId, _, details = self.nefT.get_valid_star_atom_in_xplor(cifCompId, __atomId, leave_unmatched=True)
                    elif len(atomId) > 1 and not atomId[-1].isalpha()\
                            and (atomId[0] in PSE_PRO_BEGIN_CODE or atomId[0] in ('C', 'N', 'P', 'F')):
                        _atomId, _, details = self.nefT.get_valid_star_atom_in_xplor(cifCompId, atomId[:-1], leave_unmatched=True)
                        if atomId[-1].isdigit() and int(atomId[-1]) <= len(_atomId):
                            _atomId = [_atomId[int(atomId[-1]) - 1]]

            if details is not None or atomId.endswith('"'):
                _atomId_ = translateToStdAtomName(atomId, cifCompId, ccU=self.ccU)
                if _atomId_ != atomId:
                    if atomId.startswith('HT') and len(_atomId_) == 2:
                        _atomId_ = 'H'
                    __atomId__ = self.nefT.get_valid_star_atom_in_xplor(cifCompId, _atomId_)[0]
                    if coordAtomSite is not None:
                        if any(True for _atomId_ in __atomId__ if _atomId_ in coordAtomSite['atom_id']):
                            _atomId = __atomId__
                        elif __atomId__[0][0] in PROTON_BEGIN_CODE:
                            __bondedTo = self.ccU.getBondedAtoms(cifCompId, __atomId__[0])
                            if len(__bondedTo) > 0 and __bondedTo[0] in coordAtomSite['atom_id']:
                                _atomId = __atomId__
                elif coordAtomSite is not None:
                    _atomId = []
            # _atomId = self.nefT.get_valid_star_atom(cifCompId, atomId)[0]

            if coordAtomSite is not None\
               and not any(True for _atomId_ in _atomId if _atomId_ in coordAtomSite['atom_id']):
                if atomId in coordAtomSite['atom_id']:
                    _atomId = [atomId]
                elif seqId == 1 and atomId == 'H1' and self.csStat.peptideLike(compId) and 'H' in coordAtomSite['atom_id']:
                    _atomId = ['H']

            if coordAtomSite is None and not isPolySeq and self.hasNonPolySeq:
                try:
                    for np in self.nonPolySeq:
                        if np['auth_chain_id'] == chainId and cifSeqId in np['auth_seq_id']:
                            cifSeqId = np['seq_id'][np['auth_seq_id'].index(cifSeqId)]
                            seqKey, coordAtomSite = self.getCoordAtomSiteOf(chainId, cifSeqId, cifCompId)
                            if coordAtomSite is not None:
                                break
                except ValueError:
                    pass

            if coordAtomSite is not None:
                atomSiteAtomId = coordAtomSite['atom_id']
                if len(_atomId) == 0 and __atomId in ZINC_ION_CODE and 'ZN' in atomSiteAtomId:
                    compId = atomId = 'ZN'
                    _atomId = [atomId]
                elif len(_atomId) == 0 and __atomId in CALCIUM_ION_CODE and 'CA' in atomSiteAtomId:
                    compId = atomId = 'CA'
                    _atomId = [atomId]
                elif not any(_atomId_ in atomSiteAtomId for _atomId_ in _atomId):
                    pass
                elif atomId[0] not in PSE_PRO_BEGIN_CODE and not all(_atomId in atomSiteAtomId for _atomId in _atomId):
                    _atomId = [_atomId_ for _atomId_ in _atomId if _atomId_ in atomSiteAtomId]

            lenAtomId = len(_atomId)
            if self.reasons is not None and compId != cifCompId and __compId == cifCompId:
                compId = cifCompId
            if compId != cifCompId and compId in STD_MON_DICT and cifCompId in STD_MON_DICT:
                multiChain = insCode = False
                if len(chainAssign) > 0:
                    chainIds = [ca[0] for ca in chainAssign]
                    multiChain = len(collections.Counter(chainIds).most_common()) > 1
                ps = next((ps for ps in self.polySeq if ps['auth_chain_id'] == chainId), None)
                if ps is not None:
                    compIds = [_compId for _seqId, _compId in zip(ps['auth_seq_id'], ps['comp_id']) if _seqId == cifSeqId]
                    if compId in compIds:
                        insCode = True
                        cifCompId = compId
                if not multiChain and not insCode:
                    if self.__preferAuthSeq:
                        _seqKey, _coordAtomSite = self.getCoordAtomSiteOf(chainId, cifSeqId, cifCompId, asis=False)
                        if _coordAtomSite is not None and _coordAtomSite['comp_id'] == compId:
                            if lenAtomId > 0 and _atomId[0] in _coordAtomSite['atom_id']:
                                self.__authSeqId = 'label_seq_id'
                                self.__setLocalSeqScheme()
                                continue
                    self.f.append(f"[Sequence mismatch] {self.getCurrentSpectralPeak(n=index)}"
                                  f"Residue name {__compId!r} of the spectral peak list does not match with "
                                  f"{chainId}:{cifSeqId}:{cifCompId} of the coordinates.")
                    continue

            if compId != cifCompId and cifCompId in STD_MON_DICT and not isPolySeq:
                continue

            if lenAtomId == 0 and not isPolySeq and cifCompId in ELEMENT_SYMBOLS:
                _atomId = [cifCompId]
                lenAtomId = 1

            if lenAtomId == 0:
                if compId != cifCompId and any(True for item in chainAssign if item[2] == compId):
                    continue
                if seqId == 1 and isPolySeq and cifCompId == 'ACE' and cifCompId != compId and offset == 0:
                    self.selectCoordAtoms(chainAssign, authChainId, seqId, compId, atomId, index, allowAmbig, offset=1)
                    return
                self.f.append(f"[Invalid atom nomenclature] {self.getCurrentSpectralPeak(n=index)}"
                              f"{seqId}:{__compId}:{__atomId} is invalid atom nomenclature.")
                continue
            if lenAtomId > 1 and not allowAmbig:
                self.f.append(f"[Invalid atom selection] {self.getCurrentSpectralPeak(n=index)}"
                              f"Ambiguous atom selection '{seqId}:{__compId}:{__atomId}' is not allowed.")
                continue

            if __compId != cifCompId and __compId not in self.compIdMap:
                self.compIdMap[__compId] = cifCompId

            for cifAtomId in _atomId:

                if seqKey in self.__coordUnobsRes and cifCompId in STD_MON_DICT\
                   and self.reasons is not None and 'non_poly_remap' in self.reasons:
                    if self.ccU.updateChemCompDict(cifCompId):
                        try:
                            next(cca for cca in self.ccU.lastAtomDictList
                                 if cca['atom_id'] == cifAtomId and cca['leaving_atom_flag'] != 'Y')
                        except StopIteration:
                            continue
                        try:
                            if len(authAtomId) > len(cifAtomId):
                                next(cca for cca in self.ccU.lastAtomDictList
                                     if cca['atom_id'] == authAtomId and cca['leaving_atom_flag'] != 'Y')
                        except StopIteration:
                            break

                if authAtomId in ('H', 'HN') and cifAtomId in ('HN1', 'HN2', 'HNA') and self.csStat.peptideLike(cifCompId)\
                   and coordAtomSite is not None and cifAtomId not in coordAtomSite['atom_id']:
                    if cifAtomId in ('HN2', 'HNA'):
                        if 'H2' not in coordAtomSite['atom_id']:
                            continue
                        cifAtomId = 'H2'
                    if cifAtomId == 'HN1' and 'H' in coordAtomSite['atom_id']:
                        cifAtomId = 'H'

                atomSelection.append({'chain_id': chainId, 'seq_id': cifSeqId, 'comp_id': cifCompId, 'atom_id': cifAtomId,
                                      'auth_chain_id': authChainId, 'auth_seq_id': seqId,
                                      'auth_comp_id': compId, 'auth_atom_id': authAtomId})

                self.testCoordAtomIdConsistency(chainId, cifSeqId, cifCompId, cifAtomId, seqKey, coordAtomSite, index)

        if len(atomSelection) > 0:
            self.atomSelectionSet.append(atomSelection)

    def testCoordAtomIdConsistency(self, chainId: str, seqId: int, compId: str, atomId: str,
                                   seqKey: Tuple[str, int], coordAtomSite: Optional[dict], index: int
                                   ) -> Tuple[str, bool]:
        """ Perform consistency test for each atom in reference to the coordinates.
        """

        asis = False
        if not self.__hasCoord:
            return atomId, asis

        found = False

        if coordAtomSite is not None:
            if atomId in coordAtomSite['atom_id']:
                found = True
            elif atomId in ('HN1', 'HN2', 'HN3') and (f'{atomId[-1]}HN' in coordAtomSite['atom_id']
                                                      or f'H{atomId[-1]}' in coordAtomSite['atom_id']):
                atomId = f'{atomId[-1]}HN' if f'{atomId[-1]}HN' in coordAtomSite['atom_id'] else f'H{atomId[-1]}'
                found = True
            elif 'alt_atom_id' in coordAtomSite and atomId in coordAtomSite['alt_atom_id']:
                found = True
                # self.__authAtomId = 'auth_atom_id'

            elif self.__preferAuthSeq:
                _seqKey, _coordAtomSite = self.getCoordAtomSiteOf(chainId, seqId, compId, asis=False)
                if _coordAtomSite is not None and _coordAtomSite['comp_id'] == compId and not self.__extendAuthSeq:
                    if atomId in _coordAtomSite['atom_id']:
                        found = True
                        self.__preferAuthSeq = False
                        self.__authSeqId = 'label_seq_id'
                        seqKey = _seqKey
                        self.__setLocalSeqScheme()
                    elif atomId in ('HN1', 'HN2', 'HN3') and (f'{atomId[-1]}HN' in _coordAtomSite['atom_id']
                                                              or f'H{atomId[-1]}' in _coordAtomSite['atom_id']):
                        atomId = f'{atomId[-1]}HN' if f'{atomId[-1]}HN' in _coordAtomSite['atom_id'] else f'H{atomId[-1]}'
                        found = True
                        self.__preferAuthSeq = False
                        self.__authSeqId = 'label_seq_id'
                        seqKey = _seqKey
                        self.__setLocalSeqScheme()
                    elif 'alt_atom_id' in _coordAtomSite and atomId in _coordAtomSite['alt_atom_id']:
                        found = True
                        self.__preferAuthSeq = False
                        self.__authSeqId = 'label_seq_id'
                        # self.__authAtomId = 'auth_atom_id'
                        seqKey = _seqKey
                        self.__setLocalSeqScheme()

            else:
                self.__preferAuthSeq = True
                _seqKey, _coordAtomSite = self.getCoordAtomSiteOf(chainId, seqId, compId)
                if _coordAtomSite is not None and _coordAtomSite['comp_id'] == compId:
                    if atomId in _coordAtomSite['atom_id']:
                        found = True
                        self.__authSeqId = 'auth_seq_id'
                        seqKey = _seqKey
                        self.__setLocalSeqScheme()
                    elif atomId in ('HN1', 'HN2', 'HN3') and (f'{atomId[-1]}HN' in _coordAtomSite['atom_id']
                                                              or f'H{atomId[-1]}' in _coordAtomSite['atom_id']):
                        atomId = f'{atomId[-1]}HN' if f'{atomId[-1]}HN' in _coordAtomSite['atom_id'] else f'H{atomId[-1]}'
                        found = True
                        self.__authSeqId = 'auth_seq_id'
                        seqKey = _seqKey
                        self.__setLocalSeqScheme()
                    elif 'alt_atom_id' in _coordAtomSite and atomId in _coordAtomSite['alt_atom_id']:
                        found = True
                        self.__authSeqId = 'auth_seq_id'
                        # self.__authAtomId = 'auth_atom_id'
                        seqKey = _seqKey
                        self.__setLocalSeqScheme()
                    elif not self.__extendAuthSeq:
                        self.__preferAuthSeq = False
                elif not self.__extendAuthSeq:
                    self.__preferAuthSeq = False

        elif self.__preferAuthSeq and seqKey not in self.__coordUnobsRes:
            _seqKey, _coordAtomSite = self.getCoordAtomSiteOf(chainId, seqId, compId, asis=False)
            if _coordAtomSite is not None and _coordAtomSite['comp_id'] == compId and not self.__extendAuthSeq:
                if atomId in _coordAtomSite['atom_id']:
                    found = True
                    self.__preferAuthSeq = False
                    self.__authSeqId = 'label_seq_id'
                    seqKey = _seqKey
                    self.__setLocalSeqScheme()
                elif atomId in ('HN1', 'HN2', 'HN3') and (f'{atomId[-1]}HN' in _coordAtomSite['atom_id']
                                                          or f'H{atomId[-1]}' in _coordAtomSite['atom_id']):
                    atomId = f'{atomId[-1]}HN' if f'{atomId[-1]}HN' in _coordAtomSite['atom_id'] else f'H{atomId[-1]}'
                    found = True
                    self.__preferAuthSeq = False
                    self.__authSeqId = 'label_seq_id'
                    seqKey = _seqKey
                    self.__setLocalSeqScheme()
                elif 'alt_atom_id' in _coordAtomSite and atomId in _coordAtomSite['alt_atom_id']:
                    found = True
                    self.__preferAuthSeq = False
                    self.__authSeqId = 'label_seq_id'
                    # self.__authAtomId = 'auth_atom_id'
                    seqKey = _seqKey
                    self.__setLocalSeqScheme()

        elif not self.__preferAuthSeq:
            self.__preferAuthSeq = True
            _seqKey, _coordAtomSite = self.getCoordAtomSiteOf(chainId, seqId, compId)
            if _coordAtomSite is not None and _coordAtomSite['comp_id'] == compId:
                if atomId in _coordAtomSite['atom_id']:
                    found = True
                    self.__authSeqId = 'auth_seq_id'
                    seqKey = _seqKey
                    self.__setLocalSeqScheme()
                elif atomId in ('HN1', 'HN2', 'HN3') and (f'{atomId[-1]}HN' in _coordAtomSite['atom_id']
                                                          or f'H{atomId[-1]}' in _coordAtomSite['atom_id']):
                    atomId = f'{atomId[-1]}HN' if f'{atomId[-1]}HN' in _coordAtomSite['atom_id'] else f'H{atomId[-1]}'
                    found = True
                    self.__authSeqId = 'auth_seq_id'
                    seqKey = _seqKey
                    self.__setLocalSeqScheme()
                elif 'alt_atom_id' in _coordAtomSite and atomId in _coordAtomSite['alt_atom_id']:
                    found = True
                    self.__authSeqId = 'auth_seq_id'
                    # self.__authAtomId = 'auth_atom_id'
                    seqKey = _seqKey
                    self.__setLocalSeqScheme()
                elif not self.__extendAuthSeq:
                    self.__preferAuthSeq = False
            elif not self.__extendAuthSeq:
                self.__preferAuthSeq = False

        if found:
            if self.__preferAuthSeq:
                self.__preferAuthSeqCount += 1
            return atomId, asis

        if chainId in self.chainNumberDict.values():

            if self.__preferAuthSeq and seqKey not in self.__coordUnobsRes:
                _seqKey, _coordAtomSite = self.getCoordAtomSiteOf(chainId, seqId, compId, asis=False)
                if _coordAtomSite is not None and _coordAtomSite['comp_id'] == compId and not self.__extendAuthSeq:
                    if atomId in _coordAtomSite['atom_id']:
                        found = True
                        self.__preferAuthSeq = False
                        self.__authSeqId = 'label_seq_id'
                        seqKey = _seqKey
                        self.__setLocalSeqScheme()
                    elif atomId in ('HN1', 'HN2', 'HN3') and (f'{atomId[-1]}HN' in _coordAtomSite['atom_id']
                                                              or f'H{atomId[-1]}' in _coordAtomSite['atom_id']):
                        atomId = f'{atomId[-1]}HN' if f'{atomId[-1]}HN' in _coordAtomSite['atom_id'] else f'H{atomId[-1]}'
                        found = True
                        self.__preferAuthSeq = False
                        self.__authSeqId = 'label_seq_id'
                        seqKey = _seqKey
                        self.__setLocalSeqScheme()
                    elif 'alt_atom_id' in _coordAtomSite and atomId in _coordAtomSite['alt_atom_id']:
                        found = True
                        self.__preferAuthSeq = False
                        self.__authSeqId = 'label_seq_id'
                        # self.__authAtomId = 'auth_atom_id'
                        seqKey = _seqKey
                        self.__setLocalSeqScheme()
            elif not self.__preferAuthSeq:
                self.__preferAuthSeq = True
                _seqKey, _coordAtomSite = self.getCoordAtomSiteOf(chainId, seqId, compId)
                if _coordAtomSite is not None and _coordAtomSite['comp_id'] == compId:
                    if atomId in _coordAtomSite['atom_id']:
                        found = True
                        self.__authSeqId = 'auth_seq_id'
                        seqKey = _seqKey
                        self.__setLocalSeqScheme()
                    elif atomId in ('HN1', 'HN2', 'HN3') and (f'{atomId[-1]}HN' in _coordAtomSite['atom_id']
                                                              or f'H{atomId[-1]}' in _coordAtomSite['atom_id']):
                        atomId = f'{atomId[-1]}HN' if f'{atomId[-1]}HN' in _coordAtomSite['atom_id'] else f'H{atomId[-1]}'
                        found = True
                        self.__authSeqId = 'auth_seq_id'
                        seqKey = _seqKey
                        self.__setLocalSeqScheme()
                    elif 'alt_atom_id' in _coordAtomSite and atomId in _coordAtomSite['alt_atom_id']:
                        found = True
                        self.__authSeqId = 'auth_seq_id'
                        # self.__authAtomId = 'auth_atom_id'
                        seqKey = _seqKey
                        self.__setLocalSeqScheme()
                    elif not self.__extendAuthSeq:
                        self.__preferAuthSeq = False
                elif not self.__extendAuthSeq:
                    self.__preferAuthSeq = False

            if found:
                if self.__preferAuthSeq:
                    self.__preferAuthSeqCount += 1
                return atomId, asis

        if self.ccU.updateChemCompDict(compId):
            cca = next((cca for cca in self.ccU.lastAtomDictList if cca['atom_id'] == atomId), None)
            if cca is not None and seqKey not in self.__coordUnobsRes\
               and self.ccU.lastChemCompDict['release_status'] == 'REL':
                checked = False
                ps = next((ps for ps in self.polySeq if ps['auth_chain_id'] == chainId), None)
                auth_seq_id_list = list(filter(None, ps['auth_seq_id'])) if ps is not None else None
                min_auth_seq_id = max_auth_seq_id = UNREAL_AUTH_SEQ_NUM
                if auth_seq_id_list is not None and len(auth_seq_id_list) > 0:
                    min_auth_seq_id = min(auth_seq_id_list)
                    max_auth_seq_id = max(auth_seq_id_list)
                if seqId == 1 or (chainId, seqId - 1) in self.__coordUnobsRes\
                   or seqId == min_auth_seq_id:
                    if atomId in AMINO_PROTON_CODE and atomId != 'H1':
                        return self.testCoordAtomIdConsistency(chainId, seqId, compId, 'H1', seqKey, coordAtomSite, index)
                    if atomId in AMINO_PROTON_CODE or atomId == 'P' or atomId.startswith('HOP'):
                        checked = True
                if not checked:
                    if atomId[0] in PROTON_BEGIN_CODE:
                        bondedTo = self.ccU.getBondedAtoms(compId, atomId)
                        if len(bondedTo) > 0 and bondedTo[0][0] != 'P':
                            if coordAtomSite is not None and bondedTo[0] in coordAtomSite['atom_id']:
                                if cca['leaving_atom_flag'] != 'Y'\
                                   or (self.csStat.peptideLike(compId)
                                       and cca['n_terminal_atom_flag'] == 'N'
                                       and cca['c_terminal_atom_flag'] == 'N'):
                                    self.f.append(f"[Hydrogen not instantiated] {self.getCurrentSpectralPeak(n=index)}"
                                                  f"{chainId}:{seqId}:{compId}:{atomId} is not properly instantiated "
                                                  "in the coordinates. Please re-upload the model file.")
                                    return atomId, asis
                            if bondedTo[0][0] == 'O':
                                return 'Ignorable hydroxyl group', asis
                    if seqId == max_auth_seq_id\
                       or (chainId, seqId + 1) in self.__coordUnobsRes and self.csStat.peptideLike(compId):
                        if coordAtomSite is not None and atomId in CARBOXYL_CODE\
                           and not isCyclicPolymer(self.cR, self.polySeq, chainId,
                                                   self.representativeModelId, self.representativeAltId, self.modelNumName):
                            self.f.append(f"[Coordinate issue] {self.getCurrentSpectralPeak(n=index)}"
                                          f"{chainId}:{seqId}:{compId}:{atomId} is not properly instantiated in the coordinates. "
                                          "Please re-upload the model file.")
                            return atomId, asis

                    ext_seq = False
                    if auth_seq_id_list is not None and len(auth_seq_id_list) > 0:
                        if (compId == 'ACE' and seqId == min_auth_seq_id - 1)\
                           or (compId == 'NH2' and seqId == max_auth_seq_id + 1)\
                           or (compId in STD_MON_DICT
                               and self.__preferAuthSeqCount - self.__preferLabelSeqCount >= MAX_PREF_LABEL_SCHEME_COUNT
                               and (min_auth_seq_id - MAX_ALLOWED_EXT_SEQ <= seqId < min_auth_seq_id
                                    or max_auth_seq_id < seqId <= max_auth_seq_id + MAX_ALLOWED_EXT_SEQ)):
                            ext_seq = True
                    if chainId in LARGE_ASYM_ID:
                        if ext_seq:
                            return atomId, asis
                        if self.__allow_ext_seq:
                            self.f.append(f"[Sequence mismatch warning] {self.getCurrentSpectralPeak(n=index)}"
                                          f"The residue '{chainId}:{seqId}:{compId}' is not present in polymer sequence "
                                          f"of chain {chainId} of the coordinates. {INSTRUCTION_FOR_FULL_SEQUENCE}")
                            asis = True
                        else:
                            if seqKey in self.__coordUnobsAtom\
                               and (atomId in self.__coordUnobsAtom[seqKey]['atom_ids']
                                    or (atomId[0] in PROTON_BEGIN_CODE
                                        and any(True for bondedTo in self.ccU.getBondedAtoms(compId, atomId, exclProton=True)
                                                if bondedTo in self.__coordUnobsAtom[seqKey]['atom_ids']))):
                                if self.no_extra_comment:
                                    self.f.append(f"[Coordinate issue] {self.getCurrentSpectralPeak(n=index)}"
                                                  f"{chainId}:{seqId}:{compId}:{atomId} is not present in the coordinates.")
                                return atomId, asis
                            if self.no_extra_comment:
                                self.f.append(f"[Atom not found] {self.getCurrentSpectralPeak(n=index)}"
                                              f"{chainId}:{seqId}:{compId}:{atomId} is not present in the coordinates.")
                            updatePolySeqRst(self.polySeqRstFailed, chainId, seqId, compId)
        return atomId, asis

    def getCoordAtomSiteOf(self, chainId: str, seqId: int, compId: Optional[str] = None, cifCheck: bool = True, asis: bool = True
                           ) -> Tuple[Tuple[str, int], Optional[dict]]:
        """ A wrapper function to retrieve sequence key and coordinates for given conditions.
        """

        return self.__getCoordAtomSiteOf(chainId, seqId, compId, cifCheck, asis, self.__preferAuthSeq)

    @functools.lru_cache(maxsize=2048)
    def __getCoordAtomSiteOf(self, chainId: str, seqId: int, compId: Optional[str] = None, cifCheck: bool = True, asis: bool = True,
                             __preferAuthSeq: bool = True
                             ) -> Tuple[Tuple[str, int], Optional[dict]]:
        seqKey = (chainId, seqId)
        if cifCheck:
            preferAuthSeq = __preferAuthSeq if asis else not __preferAuthSeq
            if preferAuthSeq:
                if compId is not None:
                    _seqKey = (chainId, seqId, compId)
                    if _seqKey in self.__coordAtomSite:
                        return seqKey, self.__coordAtomSite[_seqKey]
                if seqKey in self.__coordAtomSite:
                    if compId is None:
                        return seqKey, self.__coordAtomSite[seqKey]
                    _compId = self.__coordAtomSite[seqKey]['comp_id']
                    if compId == _compId:
                        return seqKey, self.__coordAtomSite[seqKey]
                    if self.hasNonPoly:
                        npList = [np for np in self.nonPoly if np['auth_chain_id'] == chainId]
                        for np in npList:
                            if np['comp_id'][0] == compId and np['auth_seq_id'][0] == seqId:
                                _seqKey = (chainId, np['seq_id'][0])
                                if _seqKey in self.__coordAtomSite and self.__coordAtomSite[_seqKey]['comp_id'] == compId:
                                    return _seqKey, self.__coordAtomSite[_seqKey]
                    return seqKey, self.__coordAtomSite[seqKey]
            else:
                if seqKey in self.__labelToAuthSeq:
                    seqKey = self.__labelToAuthSeq[seqKey]
                    if cifCheck and compId is not None:
                        _seqKey = (seqKey[0], seqKey[1], compId)
                        if _seqKey in self.__coordAtomSite:
                            return seqKey, self.__coordAtomSite[_seqKey]
                    if seqKey in self.__coordAtomSite:
                        return seqKey, self.__coordAtomSite[seqKey]
        return seqKey, None

    def getCurrentSpectralPeak(self, n: int) -> str:
        """ Retrieve indicator of the current spectral peak.
        """

        if self.cur_subtype == 'peak2d':
            return f"[Check the {self.peaks2D}th row of 2D spectral peaks "\
                f"(list_id={self.cur_list_id}, index={n}), {self.__def_err_sf_framecode}] "
        if self.cur_subtype == 'peak3d':
            return f"[Check the {self.peaks3D}th row of 3D spectral peaks "\
                f"(list_id={self.cur_list_id}, index={n}), {self.__def_err_sf_framecode}] "
        if self.cur_subtype == 'peak4d':
            return f"[Check the {self.peaks4D}th row of 4D spectral peaks "\
                f"(list_id={self.cur_list_id}, index~{n}), {self.__def_err_sf_framecode}] "
        return ''

    def __getNamedReasonsForReparsing(self, name: str) -> dict:
        """ Retrieve reasons for re-parsing by name. create a dictionary if name does not exist.
        """

        if name in self.reasonsForReParsing:
            return self.reasonsForReParsing[name]
        self.reasonsForReParsing[name] = {}
        return self.reasonsForReParsing[name]

    def __setLocalSeqScheme(self) -> None:
        """ Set sequence scheme for each spectral peak.
        """

        r = self.__getNamedReasonsForReparsing('local_seq_scheme')
        preferAuthSeq = self.__authSeqId == 'auth_seq_id'
        if self.cur_subtype == 'peak2d':
            r[(self.cur_subtype, self.cur_list_id, self.peaks2D)] = preferAuthSeq
        if self.cur_subtype == 'peak3d':
            r[(self.cur_subtype, self.cur_list_id, self.peaks3D)] = preferAuthSeq
        if self.cur_subtype == 'peak4d':
            r[(self.cur_subtype, self.cur_list_id, self.peaks4D)] = preferAuthSeq
        if not preferAuthSeq:
            self.__preferLabelSeqCount += 1
            if self.__preferLabelSeqCount > MAX_PREF_LABEL_SCHEME_COUNT:
                self.reasonsForReParsing['label_seq_scheme'] = True

    def retrieveLocalSeqScheme(self) -> None:
        """ Retrieve sequence scheme for each spectral peak.
        """

        if self.reasons is None\
           or ('label_seq_scheme' not in self.reasons
               and 'local_seq_scheme' not in self.reasons
               and 'extend_seq_scheme' not in self.reasons):
            return
        if 'extend_seq_scheme' in self.reasons:
            self.__preferAuthSeq = self.__extendAuthSeq = True
            return
        if 'label_seq_scheme' in self.reasons and self.reasons['label_seq_scheme']:
            self.__preferAuthSeq = False
            # self.__authSeqId = 'label_seq_id'
            return
        if self.cur_subtype == 'peak2d':
            key = (self.cur_subtype, self.cur_list_id, self.peaks2D)
        if self.cur_subtype == 'peak3d':
            key = (self.cur_subtype, self.cur_list_id, self.peaks3D)
        if self.cur_subtype == 'peak4d':
            key = (self.cur_subtype, self.cur_list_id, self.peaks4D)
        else:
            return

        if key in self.reasons['local_seq_scheme']:
            self.__preferAuthSeq = self.reasons['local_seq_scheme'][key]

    def __addSf(self) -> None:
        """ Add saveframe for given conditions if not exists.
        """

        content_subtype = contentSubtypeOf(self.cur_subtype)

        if content_subtype is None:
            return

        self.__listIdCounter = incListIdCounter(self.cur_subtype, self.__listIdCounter, reservedListIds=self.__reservedListIds)

        key = (self.cur_subtype, self.cur_list_id)

        if key in self.sfDict:
            if len(self.sfDict[key]) > 0:
                decListIdCounter(self.cur_subtype, self.__listIdCounter, reservedListIds=self.__reservedListIds)
                return
        else:
            self.sfDict[key] = []

        list_id = self.__listIdCounter[content_subtype]

        restraint_name = getRestraintName(self.cur_subtype)

        sf_framecode = (f'{self.software_name}_' if self.software_name is not None else '')\
            + restraint_name.replace(' ', '_') + f'_{list_id}'

        sf = getSaveframe(self.cur_subtype, sf_framecode, list_id, self.__entryId, self.__originalFileName,
                          numOfDim=self.num_of_dim, spectrumName=self.spectrum_name)

        lp = getPkLoop(self.cur_subtype)

        alt_loops = getAltLoops(content_subtype)

        item = {'file_type': self.file_type, 'saveframe': sf, 'loop': lp, 'alt_loops': alt_loops, 'list_id': list_id,
                'id': 0, 'index_id': 0, 'row_index_id': 0, 'num_of_dim': self.num_of_dim, 'peak_row_format': True,
                'sf_framecode': sf_framecode}

        self.sfDict[key].append(item)

    def getSf(self) -> dict:
        """ Retrieve saveframe for given conditions.
        """

        key = (self.cur_subtype, self.cur_list_id)

        if key not in self.sfDict:
            self.__addSf()

        cur_sf = self.sfDict[key][-1]

        self.__def_err_sf_framecode = cur_sf['sf_framecode']

        return cur_sf

    def getContentSubtype(self) -> dict:
        """ Return content subtype of PK file.
        """

        n = sum(v for v in self.listIdInternal.values())

        return {'spectral_peak': n} if n > 0 else {}

    def getPolymerSequence(self) -> Optional[List[dict]]:
        """ Return polymer sequence of PK file.
        """

        return None if self.polySeqRst is None or len(self.polySeqRst) == 0 else self.polySeqRst

    def getSequenceAlignment(self) -> Optional[List[dict]]:
        """ Return sequence alignment between coordinates and PK.
        """

        return None if self.seqAlign is None or len(self.seqAlign) == 0 else self.seqAlign

    def getChainAssignment(self) -> Optional[List[dict]]:
        """ Return chain assignment between coordinates and PK.
        """

        return None if self.chainAssign is None or len(self.chainAssign) == 0 else self.chainAssign

    def getReasonsForReparsing(self) -> Optional[dict]:
        """ Return reasons for re-parsing PK file.
        """

        return None if len(self.reasonsForReParsing) == 0 else self.reasonsForReParsing

    def getSfDict(self) -> Tuple[dict, Optional[dict]]:
        """ Return a dictionary of pynmrstar saveframes.
        """

        self.__listIdCounter, sfDict = getSfDictOf(self.sfDict, self.__listIdCounter, self.__reservedListIds)
        return self.__listIdCounter, sfDict
