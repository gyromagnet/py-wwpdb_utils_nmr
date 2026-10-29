##
# File: NmrDpRemediationLegacyMr.py
# Date: 17-Aug-2026
#
# Updates:
##
""" Validation of legacy restraint files during NMR data remediation.
    @author: Masashi Yokochi
"""
__docformat__ = "restructuredtext en"
__author__ = "Masashi Yokochi"
__email__ = "yokochi@protein.osaka-u.ac.jp"
__license__ = "Apache License 2.0"
__version__ = "5.4.4"

import copy
import os
import re
from operator import itemgetter
from typing import Optional, Tuple

try:
    from wwpdb.utils.nmr.NmrDpConstant import (AR_FILE_PATH_LIST_KEY,
                                               EMPTY_VALUE,
                                               WORK_MODEL_FILE_NAME_PAT,
                                               SEQ_MISMATCH_WARNING_PAT,
                                               INCONSISTENT_RESTRAINT_WARNING_PAT)
    from wwpdb.utils.nmr.AlignUtil import (getScoreOfSeqAlign,
                                           updatePolySeqRst,
                                           sortPolySeqRst,
                                           alignPolymerSequence,
                                           assignPolymerSequence,
                                           trimSequenceAlignment)
    from wwpdb.utils.nmr.CifToNmrStar import has_key_value
    from wwpdb.utils.nmr.mr.ParserListenerUtil import (getRestraintName,
                                                       contentSubtypeOf,
                                                       incListIdCounter,
                                                       getSaveframe,
                                                       getLoop,
                                                       getRow)
    from wwpdb.utils.nmr.mr.AmberMRReader import AmberMRReader
    from wwpdb.utils.nmr.mr.AmberPTReader import AmberPTReader
    from wwpdb.utils.nmr.mr.AriaMRReader import AriaMRReader
    from wwpdb.utils.nmr.mr.AriaMRXReader import AriaMRXReader
    from wwpdb.utils.nmr.mr.BareMRReader import BareMRReader
    from wwpdb.utils.nmr.mr.BarePDBReader import BarePDBReader
    from wwpdb.utils.nmr.mr.BiosymMRReader import BiosymMRReader
    from wwpdb.utils.nmr.mr.CharmmCRDReader import CharmmCRDReader
    from wwpdb.utils.nmr.mr.CharmmMRReader import CharmmMRReader
    from wwpdb.utils.nmr.mr.CnsMRReader import CnsMRReader
    from wwpdb.utils.nmr.mr.CyanaMRReader import CyanaMRReader
    from wwpdb.utils.nmr.mr.CyanaNOAReader import CyanaNOAReader
    from wwpdb.utils.nmr.mr.DynamoMRReader import DynamoMRReader
    from wwpdb.utils.nmr.mr.GromacsMRReader import GromacsMRReader
    from wwpdb.utils.nmr.mr.GromacsPTReader import GromacsPTReader
    from wwpdb.utils.nmr.mr.IsdMRReader import IsdMRReader
    from wwpdb.utils.nmr.mr.RosettaMRReader import RosettaMRReader
    from wwpdb.utils.nmr.mr.SchrodingerMRReader import SchrodingerMRReader
    from wwpdb.utils.nmr.mr.SybylMRReader import SybylMRReader
    from wwpdb.utils.nmr.mr.XplorMRReader import XplorMRReader
    from wwpdb.utils.nmr.NmrDpRemediationBase import NmrDpRemediationBase
except ImportError:
    from nmr.NmrDpConstant import (AR_FILE_PATH_LIST_KEY,
                                   EMPTY_VALUE,
                                   WORK_MODEL_FILE_NAME_PAT,
                                   SEQ_MISMATCH_WARNING_PAT,
                                   INCONSISTENT_RESTRAINT_WARNING_PAT)
    from nmr.AlignUtil import (getScoreOfSeqAlign,
                               updatePolySeqRst,
                               sortPolySeqRst,
                               alignPolymerSequence,
                               assignPolymerSequence,
                               trimSequenceAlignment)
    from nmr.CifToNmrStar import has_key_value
    from nmr.mr.ParserListenerUtil import (getRestraintName,
                                           contentSubtypeOf,
                                           incListIdCounter,
                                           getSaveframe,
                                           getLoop,
                                           getRow)
    from nmr.mr.AmberMRReader import AmberMRReader
    from nmr.mr.AmberPTReader import AmberPTReader
    from nmr.mr.AriaMRReader import AriaMRReader
    from nmr.mr.AriaMRXReader import AriaMRXReader
    from nmr.mr.BareMRReader import BareMRReader
    from nmr.mr.BarePDBReader import BarePDBReader
    from nmr.mr.BiosymMRReader import BiosymMRReader
    from nmr.mr.CharmmCRDReader import CharmmCRDReader
    from nmr.mr.CharmmMRReader import CharmmMRReader
    from nmr.mr.CnsMRReader import CnsMRReader
    from nmr.mr.CyanaMRReader import CyanaMRReader
    from nmr.mr.CyanaNOAReader import CyanaNOAReader
    from nmr.mr.DynamoMRReader import DynamoMRReader
    from nmr.mr.GromacsMRReader import GromacsMRReader
    from nmr.mr.GromacsPTReader import GromacsPTReader
    from nmr.mr.IsdMRReader import IsdMRReader
    from nmr.mr.RosettaMRReader import RosettaMRReader
    from nmr.mr.SchrodingerMRReader import SchrodingerMRReader
    from nmr.mr.SybylMRReader import SybylMRReader
    from nmr.mr.XplorMRReader import XplorMRReader
    from nmr.NmrDpRemediationBase import NmrDpRemediationBase

# Reader dispatch for the seven uniform branches of validateLegacyMr(). Every
# other format needs extra reader configuration, a different source of re-parse
# reasons, or its own bookkeeping, and keeps an explicit branch below.
LEGACY_MR_READERS = {
    'nm-res-ari': (AriaMRReader, 'ARIA'),
    'nm-res-arx': (AriaMRXReader, 'ARIA'),
    'nm-res-bar': (BareMRReader, 'Bare WSV/TSV/CSV'),
    'nm-res-bio': (BiosymMRReader, 'BIOSYM'),
    'nm-res-dyn': (DynamoMRReader, 'DYNAMO/PALES/TALOS'),
    'nm-res-isd': (IsdMRReader, 'ISD'),
    'nm-res-syb': (SybylMRReader, 'SYBYL'),
}


class NmrDpRemediationLegacyMr(NmrDpRemediationBase):
    """ Validation of legacy restraint files during NMR data remediation.
    """
    __slots__ = ()

    def _parseLegacyMr(self, spec: tuple, file_path: str, file_name: str,
                       original_file_name: str, create_sf_dict: bool,
                       deal_warn_for_lazy_eval) -> Tuple[bool, Optional[object]]:
        """ Parse a legacy restraint file with the reader for a given file type,
            re-parsing once if the parser listener asks for it.
            @param spec: the LEGACY_MR_READERS entry for the file type
            @return: (whether the initial parse yielded a listener, the final listener),
                     reported separately because the caller's processing was gated on
                     the *initial* parse, exactly as the per-format branches were.
        """
        reader_cls = spec[0]

        def new_reader(*trailing):
            return reader_cls(self._reg.verbose, self._reg.log,
                              self._reg.representative_model_id,
                              self._reg.representative_alt_id,
                              self._reg.mr_atom_name_mapping,
                              self._reg.cR, self._reg.caC,
                              self._reg.ccU, self._reg.csStat, self._reg.nefT,
                              *trailing)

        def parse_with(reader, list_id_counter):
            return reader.parse(file_path, self._reg.cifPath,
                                createSfDict=create_sf_dict, originalFileName=original_file_name,
                                listIdCounter=list_id_counter,
                                entryId=self._reg.entry_id)[0]

        _list_id_counter = copy.copy(self._reg.list_id_counter)

        listener = parse_with(new_reader(), self._reg.list_id_counter)

        if listener is None:
            return False, None

        reasons = listener.getReasonsForReparsing()

        if reasons is not None:
            deal_warn_for_lazy_eval(file_name, listener)

            if 'model_chain_id_ext' in reasons:
                self._reg.auth_asym_ids_with_chem_exch.update(reasons['model_chain_id_ext'])
            if 'chain_id_clone' in reasons:
                self._reg.auth_seq_ids_with_chem_exch.update(reasons['chain_id_clone'])

            listener = parse_with(new_reader(reasons), _list_id_counter)

        return True, listener

    def validateLegacyMr(self) -> bool:
        """ Validate data content of legacy restraint files.
        """

        if self._reg.combined_mode and not self._reg.bmrb_only:
            return True

        if AR_FILE_PATH_LIST_KEY not in self._reg.inputParamDict:
            return True

        src_id = self._reg.report.getInputSourceIdOfCoord()

        if src_id < 0:
            return False

        cif_input_source = self._reg.report.input_sources[src_id]
        cif_input_source_dic = cif_input_source.get()

        has_poly_seq = has_key_value(cif_input_source_dic, 'polymer_sequence')

        if not has_poly_seq:
            return False

        nmr_vs_model = None
        chain_assign_dic = self._reg.report.chain_assignment.get()
        if 'nmr_poly_seq_vs_model_poly_seq' in chain_assign_dic:
            nmr_vs_model = chain_assign_dic['nmr_poly_seq_vs_model_poly_seq']

        if self._reg.versioned_atom_name_mapping is not None:
            for atom_map in self._reg.versioned_atom_name_mapping:
                if atom_map not in self._reg.mr_atom_name_mapping:
                    self._reg.mr_atom_name_mapping.append(atom_map)

        if len(self._reg.internal_atom_name_mapping) > 0:
            for ver in sorted(list(self._reg.internal_atom_name_mapping), reverse=True):
                if self._reg.internal_atom_name_mapping[ver] is None:
                    continue
                for atom_map in self._reg.internal_atom_name_mapping[ver]:
                    if atom_map not in self._reg.mr_atom_name_mapping:
                        self._reg.mr_atom_name_mapping.append(atom_map)

        if self._reg.mr_atom_name_mapping is not None and len(self._reg.mr_atom_name_mapping) > 1:
            self._reg.mr_atom_name_mapping = list(reversed(self._reg.mr_atom_name_mapping))

        amberAtomNumberDict = charmmAtomNumberDict = gromacsAtomNumberDict = pdbAtomNumberDict = None
        _amberAtomNumberDict = {}

        has_aux_amb = has_aux_cha = has_aux_gro = False
        has_res_amb = has_res_cha = has_res_gro = has_res_sch = False

        cyanaUplDistRest = cyanaLolDistRest = 0

        def deal_aux_warn_message(file_name, listener, dry_run=False):

            valid = True

            if listener.warningMessage is not None:

                for warn in listener.warningMessage:

                    msg_dict = {'file_name': file_name, 'description': warn, 'inheritable': True}

                    if warn.startswith('[Concatenated sequence]'):
                        if not dry_run:
                            self._reg.report.warning.appendDescription('concatenated_sequence', msg_dict)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Sequence mismatch]'):
                        if not dry_run and not has_res_sch:
                            self._reg.report.error.appendDescription('sequence_mismatch', msg_dict)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Error  - {warn}\n")

                            valid = False

                    elif warn.startswith('[Unknown atom name]'):
                        if not dry_run:
                            self._reg.report.warning.appendDescription('inconsistent_mr_data', msg_dict)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Unknown residue name]'):
                        if not dry_run:
                            self._reg.report.warning.appendDescription('inconsistent_mr_data', msg_dict)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {warn}\n")

                    else:
                        if not dry_run:
                            self._reg.report.error.appendDescription('internal_error',
                                                                     f"+{self.__class_name__}.validateLegacyMr() "
                                                                     "++ KeyError  - " + warn)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ KeyError  - {warn}\n")

                        valid = False

            else:
                valid = False

            return valid

        fileListId = self._reg.file_path_list_len

        for ar in self._reg.inputParamDict[AR_FILE_PATH_LIST_KEY]:
            file_path = ar['file_name']

            input_source = self._reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_type = input_source_dic['file_type']

            fileListId += 1

            if file_type == 'nm-aux-amb':
                has_aux_amb = True

            if file_type == 'nm-aux-cha':
                has_aux_cha = True

            if file_type == 'nm-aux-gro':
                has_aux_gro = True

            if file_type == 'nm-res-amb':
                has_res_amb = True

            if file_type == 'nm-res-cha':
                has_res_cha = True

            if file_type == 'nm-res-gro':
                has_res_gro = True

            if file_type == 'nm-res-sch':
                has_res_sch = True

        fileListId = self._reg.file_path_list_len

        for ar in self._reg.inputParamDict[AR_FILE_PATH_LIST_KEY]:
            file_path = ar['file_name']

            input_source = self._reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_type = input_source_dic['file_type']
            content_subtype = input_source_dic['content_subtype']

            fileListId += 1

            if file_type == 'nm-aux-amb' and content_subtype is not None and 'topology' in content_subtype:

                if 'is_valid' in ar and ar['is_valid']:

                    file_name = input_source_dic['file_name']

                    original_file_name = None
                    if 'original_file_name' in input_source_dic:
                        if input_source_dic['original_file_name'] is not None:
                            original_file_name = os.path.basename(input_source_dic['original_file_name'])
                        if file_name != original_file_name and original_file_name is not None:
                            file_name = f"{original_file_name} ({file_name})"

                    reader = AmberPTReader(self._reg.verbose, self._reg.log,
                                           self._reg.representative_model_id,
                                           self._reg.representative_alt_id,
                                           self._reg.mr_atom_name_mapping,
                                           self._reg.cR, self._reg.caC,
                                           self._reg.ccU, self._reg.csStat, self._reg.nefT)

                    listener, _, _ = reader.parse(file_path, self._reg.cifPath)

                    if listener is not None:

                        deal_aux_warn_message(file_name, listener)

                        amberAtomNumberDict = listener.getAtomNumberDict()

                        poly_seq = listener.getPolymerSequence()
                        if poly_seq is not None:
                            input_source.setItemValue('polymer_sequence', poly_seq)

                        seq_align = listener.getSequenceAlignment()
                        if seq_align is not None:
                            self._reg.report.sequence_alignment.setItemValue('model_poly_seq_vs_mr_topology', seq_align)

            elif file_type == 'nm-aux-cha' and content_subtype is not None and 'topology' in content_subtype:

                if 'is_valid' in ar and ar['is_valid']:

                    file_name = input_source_dic['file_name']

                    original_file_name = None
                    if 'original_file_name' in input_source_dic:
                        if input_source_dic['original_file_name'] is not None:
                            original_file_name = os.path.basename(input_source_dic['original_file_name'])
                        if file_name != original_file_name and original_file_name is not None:
                            file_name = f"{original_file_name} ({file_name})"

                    reader = CharmmCRDReader(self._reg.verbose, self._reg.log,
                                             self._reg.representative_model_id,
                                             self._reg.representative_alt_id,
                                             self._reg.mr_atom_name_mapping,
                                             self._reg.cR, self._reg.caC,
                                             self._reg.ccU, self._reg.csStat, self._reg.nefT)

                    listener, _, _ = reader.parse(file_path, self._reg.cifPath)

                    if listener is not None:

                        deal_aux_warn_message(file_name, listener)

                        charmmAtomNumberDict = listener.getAtomNumberDict()

                        poly_seq = listener.getPolymerSequence()
                        if poly_seq is not None:
                            input_source.setItemValue('polymer_sequence', poly_seq)

                        seq_align = listener.getSequenceAlignment()
                        if seq_align is not None:
                            self._reg.report.sequence_alignment.setItemValue('model_poly_seq_vs_mr_topology', seq_align)

            elif file_type == 'nm-aux-gro' and content_subtype is not None and 'topology' in content_subtype:

                if 'is_valid' in ar and ar['is_valid']:

                    file_name = input_source_dic['file_name']

                    original_file_name = None
                    if 'original_file_name' in input_source_dic:
                        if input_source_dic['original_file_name'] is not None:
                            original_file_name = os.path.basename(input_source_dic['original_file_name'])
                        if file_name != original_file_name and original_file_name is not None:
                            file_name = f"{original_file_name} ({file_name})"

                    reader = GromacsPTReader(self._reg.verbose, self._reg.log,
                                             self._reg.representative_model_id,
                                             self._reg.representative_alt_id,
                                             self._reg.mr_atom_name_mapping,
                                             self._reg.cR, self._reg.caC,
                                             self._reg.ccU, self._reg.csStat, self._reg.nefT)

                    listener, _, _ = reader.parse(file_path, self._reg.cifPath)

                    if listener is not None:

                        deal_aux_warn_message(file_name, listener)

                        gromacsAtomNumberDict = listener.getAtomNumberDict()

                        poly_seq = listener.getPolymerSequence()
                        if poly_seq is not None:
                            input_source.setItemValue('polymer_sequence', poly_seq)

                        seq_align = listener.getSequenceAlignment()
                        if seq_align is not None:
                            self._reg.report.sequence_alignment.setItemValue('model_poly_seq_vs_mr_topology', seq_align)

            elif file_type == 'nm-aux-pdb' and content_subtype is not None and 'topology' in content_subtype:

                if 'is_valid' in ar and ar['is_valid']:

                    if (has_res_amb and not has_aux_amb)\
                       or (has_res_cha and not has_aux_cha)\
                       or (has_res_gro and not has_aux_gro)\
                       or has_res_sch:

                        file_name = input_source_dic['file_name']

                        original_file_name = None
                        if 'original_file_name' in input_source_dic:
                            if input_source_dic['original_file_name'] is not None:
                                original_file_name = os.path.basename(input_source_dic['original_file_name'])
                            if file_name != original_file_name and original_file_name is not None:
                                file_name = f"{original_file_name} ({file_name})"

                        reader = BarePDBReader(self._reg.verbose, self._reg.log,
                                               self._reg.representative_model_id,
                                               self._reg.representative_alt_id,
                                               self._reg.mr_atom_name_mapping,
                                               self._reg.cR, self._reg.caC,
                                               self._reg.ccU, self._reg.csStat, self._reg.nefT)

                        file_path = self.testPathWithSuffix(file_path, '-corrected')

                        listener, _, _ = reader.parse(file_path, self._reg.cifPath)

                        if listener is not None:

                            valid = deal_aux_warn_message(file_name, listener, True)

                            if not valid:
                                continue

                            poly_seq = listener.getPolymerSequence()

                            if not has_res_sch and any(True for pdb_ps in poly_seq if len(pdb_ps['seq_id']) == 0):
                                continue

                            if poly_seq is not None:
                                input_source.setItemValue('polymer_sequence', poly_seq)

                            seq_align = listener.getSequenceAlignment()
                            if seq_align is not None:
                                self._reg.report.sequence_alignment.setItemValue('model_poly_seq_vs_mr_topology', seq_align)

                                if valid or not has_res_amb:

                                    for sa in seq_align:
                                        ref_chain_id = sa['ref_chain_id']
                                        test_chain_id = sa['test_chain_id']
                                        cif_ps = next((cif_ps for cif_ps in self._reg.caC['polymer_sequence']
                                                       if cif_ps['auth_chain_id'] == ref_chain_id), None)
                                        pdb_ps = next(pdb_ps for pdb_ps in poly_seq if pdb_ps['chain_id'] == test_chain_id)

                                        if cif_ps is None:
                                            valid = False
                                            break

                                        unobs_res_count = 0
                                        for seq_id in cif_ps['seq_id']:
                                            seq_key = (ref_chain_id, seq_id)
                                            if seq_key in self._reg.caC['coord_unobs_res']:
                                                unobs_res_count += 1

                                        if len(cif_ps['seq_id']) - unobs_res_count > len(pdb_ps['seq_id']):
                                            valid = False
                                            break

                                    if valid:
                                        _pdbAtomNumberDict = listener.getAtomNumberDict()

                                        if _pdbAtomNumberDict is not None:
                                            if pdbAtomNumberDict is None:
                                                pdbAtomNumberDict = _pdbAtomNumberDict
                                            else:
                                                for k, v in _pdbAtomNumberDict.items():
                                                    if k not in pdbAtomNumberDict:
                                                        pdbAtomNumberDict[k] = v

                                        deal_aux_warn_message(file_name, listener, True)

                            elif has_res_sch and pdbAtomNumberDict is None:
                                pdbAtomNumberDict = listener.getAtomNumberDict()

            elif file_type == 'nm-res-cya' and content_subtype is not None and 'dist_restraint' in content_subtype:
                if 'is_valid' in ar and ar['is_valid']:
                    if ar['dist_type'] in ('upl', 'both'):
                        cyanaUplDistRest += 1
                    if ar['dist_type'] in ('lol', 'both'):
                        cyanaLolDistRest += 1

        if ((has_res_amb and not has_aux_amb)
            or (has_res_cha and not has_aux_cha)
            or (has_res_gro and not has_aux_gro)
            or has_res_sch)\
           and pdbAtomNumberDict is None and self._reg.internal_mode:
            cif_file_name = os.path.basename(self._reg.cifPath)
            if re.match(r'^([Pp][Dd][Bb]_)?(\d{4})?\d\w{3}.cif$', cif_file_name):
                dep_id = cif_file_name[:-4]

                file_path = os.path.join(self._reg.cR.getDirPath(), f'{dep_id}_model-upload_P1.pdb.V1')

                if os.path.exists(file_path):
                    reader = BarePDBReader(self._reg.verbose, self._reg.log,
                                           self._reg.representative_model_id,
                                           self._reg.representative_alt_id,
                                           self._reg.mr_atom_name_mapping,
                                           self._reg.cR, self._reg.caC,
                                           self._reg.ccU, self._reg.csStat, self._reg.nefT)

                    listener, _, _ = reader.parse(file_path, self._reg.cifPath)

                    if listener is not None:
                        pdbAtomNumberDict = listener.getAtomNumberDict()

        if has_res_sch and pdbAtomNumberDict is None and not self._reg.internal_mode:
            cif_file_name = os.path.basename(self._reg.cifPath)

            if WORK_MODEL_FILE_NAME_PAT.match(cif_file_name):
                dep_id = WORK_MODEL_FILE_NAME_PAT.search(cif_file_name).groups()[0]

                file_path = os.path.join(self._reg.cR.getDirPath(), f'{dep_id}_model-upload_P1.pdb.V1')

                if os.path.exists(file_path):
                    reader = BarePDBReader(self._reg.verbose, self._reg.log,
                                           self._reg.representative_model_id,
                                           self._reg.representative_alt_id,
                                           self._reg.mr_atom_name_mapping,
                                           self._reg.cR, self._reg.caC,
                                           self._reg.ccU, self._reg.csStat, self._reg.nefT)

                    listener, _, _ = reader.parse(file_path, self._reg.cifPath)

                    if listener is not None:
                        pdbAtomNumberDict = listener.getAtomNumberDict()

        if pdbAtomNumberDict is not None:
            if has_res_amb and amberAtomNumberDict is None:
                amberAtomNumberDict = pdbAtomNumberDict
            if has_res_cha and charmmAtomNumberDict is None:
                charmmAtomNumberDict = pdbAtomNumberDict
            if has_res_gro and gromacsAtomNumberDict is None:
                gromacsAtomNumberDict = pdbAtomNumberDict

        fileListId = self._reg.file_path_list_len

        # 6gbm, NOE restraint files must take precedence over other distance constraints such as hydrogen bonds
        ar_file_order, ar_file_any_dist, ar_file_wo_dist = [], [], []

        derived_from_public_mr = False

        hint_for_noe_dist = ['noe', 'roe', 'dist']
        hint_for_any_dist = ['bond', 'disul', 'not', 'seen', 'pre', 'paramag', 'cidnp',
                             'csp', 'perturb', 'mutat', 'protect', 'symm']

        for ar in self._reg.inputParamDict[AR_FILE_PATH_LIST_KEY]:
            file_path = ar['file_name']
            file_size = os.path.getsize(file_path)

            input_source = self._reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_type = input_source_dic['file_type']

            if fileListId == self._reg.file_path_list_len and file_type == 'nm-res-mr':
                derived_from_public_mr = True

            fileListId += 1

            if file_type is None\
               or file_type in ('nm-aux-amb', 'nm-aux-cha', 'nm-aux-gro', 'nm-aux-pdb', 'nm-res-oth', 'nm-res-mr', 'nm-res-sax')\
               or file_type.startswith('nm-pea'):
                continue

            if self._reg.remediation_mode and os.path.exists(self.testPathWithSuffix(file_path, '-ignored', True)):
                continue

            content_subtype = input_source_dic['content_subtype']

            if content_subtype is None or len(content_subtype) == 0:
                continue

            if 'is_valid' not in ar or not ar['is_valid']:
                continue

            file_name = input_source_dic['file_name']

            original_file_name = None
            if 'original_file_name' in input_source_dic:
                if input_source_dic['original_file_name'] is not None:
                    original_file_name = os.path.basename(input_source_dic['original_file_name'])
                if file_name != original_file_name and original_file_name is not None:
                    file_name = original_file_name

            file_name = file_name.lower()

            if content_subtype is not None and 'dist_restraint' in content_subtype.keys():
                if any(k in file_name for k in hint_for_noe_dist) and not any(k in file_name for k in hint_for_any_dist):
                    ar_file_order.append((input_source, ar, file_size))
                else:
                    ar_file_any_dist.append((input_source, ar, file_size))
            else:
                ar_file_wo_dist.append((input_source, ar, file_size))

        ar_file_order = sorted(ar_file_order, key=itemgetter(2), reverse=True)
        ar_file_order.extend(sorted(ar_file_any_dist, key=itemgetter(2), reverse=True))
        ar_file_order.extend(sorted(ar_file_wo_dist, key=itemgetter(2), reverse=True))

        poly_seq_set = []

        create_sf_dict = self._reg.remediation_mode

        if self._reg.list_id_counter is None:
            self._reg.list_id_counter = {}
        if self._reg.mr_sf_dict_holder is None:
            self._reg.mr_sf_dict_holder = {}

        if self._reg.nmr_ext_poly_seq is None:
            self._reg.nmr_ext_poly_seq = []

        if not self._reg.bmrb_only or not self._reg.internal_mode:  # nmrPolySeq is None in __retrieveCoordAssemblyChecker()

            input_source = self._reg.report.input_sources[0]
            input_source_dic = input_source.get()

            has_poly_seq = has_key_value(input_source_dic, 'polymer_sequence')

            nmr_poly_seq = input_source_dic['polymer_sequence']
            cif_poly_seq = self._reg.caC['polymer_sequence']

            seq_align, _ = alignPolymerSequence(self._reg.pA, cif_poly_seq, nmr_poly_seq)
            chain_assign, _ = assignPolymerSequence(self._reg.pA, self._reg.ccU, 'nmr-star',
                                                    cif_poly_seq, nmr_poly_seq, seq_align)

            if chain_assign is not None:

                for ca in chain_assign:
                    ref_chain_id = ca['ref_chain_id']
                    test_chain_id = ca['test_chain_id']

                    sa = next(sa for sa in seq_align
                              if sa['ref_chain_id'] == ref_chain_id
                              and sa['test_chain_id'] == test_chain_id)

                    if sa['conflict'] > 0 or sa['unmapped'] == 0:
                        continue

                    ps1 = next(ps for ps in nmr_poly_seq if ps['chain_id'] == test_chain_id)
                    ps2 = next(ps for ps in cif_poly_seq if ps['auth_chain_id'] == ref_chain_id)

                    self._reg.pA.clear()
                    self._reg.pA.setReferenceSequence(ps1['comp_id'], f'REF{test_chain_id}')
                    self._reg.pA.addTestSequence(ps2['comp_id'], test_chain_id)
                    self._reg.pA.doAlign()

                    myAlign = self._reg.pA.getAlignment(test_chain_id)

                    length = len(myAlign)

                    _matched, unmapped, conflict, offset_1, offset_2 = getScoreOfSeqAlign(myAlign)

                    if conflict == 0 and unmapped > 0:

                        nmr_seq_ids, cif_auth_seq_ids = [], []

                        for i in range(length):
                            if str(myAlign[i][0]) != '.' and i < len(ps1['seq_id']):
                                nmr_seq_ids.append(ps1['seq_id'][i])
                            else:
                                nmr_seq_ids.append(None)

                        for i in range(length):
                            if str(myAlign[i][1]) != '.' and i < len(ps2['seq_id']):
                                cif_auth_seq_ids.append(ps2['auth_seq_id'][i])
                            else:
                                cif_auth_seq_ids.append(None)

                        for i in range(length):
                            nmr_comp_id, cif_comp_id = str(myAlign[i][0]), str(myAlign[i][1])

                            if nmr_comp_id == cif_comp_id:
                                continue

                            if cif_comp_id == '.' and nmr_comp_id != '.':
                                nmr_seq_id = nmr_seq_ids[i] - offset_1 if nmr_seq_ids[i] is not None else None
                                if nmr_seq_id is not None:
                                    offset = None
                                    for _offset in range(1, 20):
                                        if i + _offset < length:
                                            _myPr = myAlign[i + _offset]
                                            if _myPr[0] == _myPr[1]:
                                                offset = _offset
                                                break
                                        if i - _offset >= 0:
                                            _myPr = myAlign[i - _offset]
                                            if _myPr[0] == _myPr[1]:
                                                offset = -_offset
                                                break

                                    if offset is not None and cif_auth_seq_ids[i + offset] is not None:
                                        cif_auth_seq_id = cif_auth_seq_ids[i + offset] - offset - offset_2
                                        self._reg.nmr_ext_poly_seq.append({'auth_chain_id': ps2['auth_chain_id'],
                                                                           'auth_seq_id': cif_auth_seq_id,
                                                                           'auth_comp_id': nmr_comp_id})

        reasons_dict = {}

        suspended_errors_for_lazy_eval = []

        def consume_suspended_message():

            if len(suspended_errors_for_lazy_eval) > 0:
                for msg in suspended_errors_for_lazy_eval:
                    for k, v in msg.items():
                        self._reg.report.error.appendDescription(k, v)
                suspended_errors_for_lazy_eval.clear()

        def deal_res_warn_message(file_name, listener, ignore_error):

            if listener.warningMessage is not None:

                for warn in listener.warningMessage:

                    msg_dict = {'file_name': file_name, 'description': warn, 'inheritable': True}
                    if INCONSISTENT_RESTRAINT_WARNING_PAT.match(warn):
                        g = INCONSISTENT_RESTRAINT_WARNING_PAT.search(warn).groups()
                        if g not in EMPTY_VALUE:
                            msg_dict['sf_framecode'] = g[1]
                            msg_dict['description'] = warn.replace(f', {g[1]}', '')

                    if warn.startswith('[Concatenated sequence]'):
                        self._reg.report.warning.appendDescription('concatenated_sequence', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Sequence mismatch]'):
                        # consume_suspended_message()

                        self._reg.report.error.appendDescription('sequence_mismatch', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Error  - {warn}\n")

                    elif warn.startswith('[Atom not found]'):

                        if not self._reg.remediation_mode or 'Macromolecules page' not in warn or self._reg.conversion_server:
                            consume_suspended_message()

                            self._reg.report.error.appendDescription('atom_not_found', msg_dict)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Error  - {warn}\n")
                        else:
                            self._reg.report.warning.appendDescription('sequence_mismatch', msg_dict)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Hydrogen not instantiated]'):

                        if self._reg.remediation_mode and not self._reg.conversion_server:

                            self._reg.report.warning.appendDescription('hydrogen_not_instantiated', msg_dict)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {warn}\n")

                        else:
                            consume_suspended_message()

                            self._reg.report.error.appendDescription('hydrogen_not_instantiated', msg_dict)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Error  - {warn}\n")

                    elif warn.startswith('[Coordinate issue]'):
                        # consume_suspended_message()

                        if self._reg.internal_mode:  # and not self._reg.conversion_server:

                            self._reg.report.warning.appendDescription('coordinate_issue', msg_dict)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {warn}\n")

                        else:

                            self._reg.report.error.appendDescription('coordinate_issue', msg_dict)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Error  - {warn}\n")

                    elif warn.startswith('[Invalid atom nomenclature]'):
                        consume_suspended_message()

                        self._reg.report.error.appendDescription('invalid_atom_nomenclature', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Error  - {warn}\n")

                    elif warn.startswith('[Invalid atom selection]') or warn.startswith('[Invalid data]'):
                        consume_suspended_message()

                        self._reg.report.error.appendDescription('invalid_data', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ ValueError  - {warn}\n")

                    elif warn.startswith('[Sequence mismatch warning]'):
                        self._reg.report.warning.appendDescription('sequence_mismatch', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {warn}\n")

                        if SEQ_MISMATCH_WARNING_PAT.match(warn):
                            g = SEQ_MISMATCH_WARNING_PAT.search(warn).groups()
                            d = {'auth_chain_id': g[2],
                                 'auth_seq_id': int(g[0]),
                                 'auth_comp_id': g[1]}
                            if d not in self._reg.nmr_ext_poly_seq:
                                self._reg.nmr_ext_poly_seq.append(d)

                    elif warn.startswith('[Missing data]'):
                        # consume_suspended_message()

                        self._reg.report.error.appendDescription('missing_data', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ ValueError  - {warn}\n")

                    elif warn.startswith('[Enum mismatch]'):
                        self._reg.report.warning.appendDescription('enum_mismatch', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Enum mismatch ignorable]'):
                        self._reg.report.warning.appendDescription('enum_mismatch_ignorable', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Unmatched atom type]'):
                        self._reg.report.warning.appendDescription('inconsistent_mr_data', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Inconsistent dihedral angle atoms]'):
                        self._reg.report.warning.appendDescription('inconsistent_mr_data', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Range value error]') and not self._reg.remediation_mode:
                        # consume_suspended_message()

                        self._reg.report.error.appendDescription('anomalous_data', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ ValueError  - {warn}\n")

                    elif warn.startswith('[Range value warning]')\
                            or (warn.startswith('[Range value error]') and self._reg.remediation_mode):
                        self._reg.report.warning.appendDescription('inconsistent_mr_data', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Insufficient atom selection]') or warn.startswith('[Insufficient angle selection]'):

                        if self._reg.conversion_server:
                            self._reg.report.error.appendDescription('unparsed_data', msg_dict)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Error  - {warn}\n")

                        else:
                            self._reg.report.warning.appendDescription('insufficient_mr_data', msg_dict)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Redundant data]'):
                        self._reg.report.warning.appendDescription('redundant_mr_data', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Ambiguous dihedral angle]'):
                        self._reg.report.warning.appendDescription('ambiguous_dihedral_angle', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Anomalous RDC vector]'):
                        self._reg.report.warning.appendDescription('anomalous_rdc_vector', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Anomalous data]'):
                        self._reg.report.warning.appendDescription('anomalous_data', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Unsupported data]'):
                        self._reg.report.warning.appendDescription('unsupported_mr_data', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {warn}\n")

                    elif not ignore_error:
                        self._reg.report.error.appendDescription('internal_error',
                                                                 f"+{self.__class_name__}.validateLegacyMr() "
                                                                 "++ KeyError  - " + warn)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ KeyError  - {warn}\n")

        def deal_res_warn_message_for_lazy_eval(file_name, listener):

            if listener.warningMessage is not None:

                def_sf_framecode = ''
                for warn in listener.warningMessage:

                    if INCONSISTENT_RESTRAINT_WARNING_PAT.match(warn):
                        g = INCONSISTENT_RESTRAINT_WARNING_PAT.search(warn).groups()
                        if g[1] not in EMPTY_VALUE:
                            def_sf_framecode = g[1]
                            break

                for warn in listener.warningMessage:

                    msg_dict = {'file_name': file_name, 'description': warn, 'inheritable': True}
                    if INCONSISTENT_RESTRAINT_WARNING_PAT.match(warn):
                        g = INCONSISTENT_RESTRAINT_WARNING_PAT.search(warn).groups()
                        msg_dict['sf_framecode'] = g[1] if g[1] not in EMPTY_VALUE else def_sf_framecode
                        msg_dict['description'] = warn.replace(f', {g[1]}', '')

                    if warn.startswith('[Sequence mismatch]'):
                        suspended_errors_for_lazy_eval.append({'sequence_mismatch': msg_dict})

                    elif warn.startswith('[Atom not found]'):
                        if not self._reg.remediation_mode or 'Macromolecules page' not in warn:
                            suspended_errors_for_lazy_eval.append({'atom_not_found': msg_dict})

                    elif warn.startswith('[Hydrogen not instantiated]'):
                        if self._reg.remediation_mode:
                            pass
                        else:
                            suspended_errors_for_lazy_eval.append({'hydrogen_not_instantiated': msg_dict})

                    # elif warn.startswith('[Coordinate issue]'):
                    #     suspended_errors_for_lazy_eval.append({'coordinate_issue': msg_dict})

                    # elif warn.startswith('[Invalid atom nomenclature]'):
                    #     suspended_errors_for_lazy_eval.append({'invalid_atom_nomenclature': msg_dict})

                    elif warn.startswith('[Invalid atom selection]') or warn.startswith('[Invalid data]'):
                        suspended_errors_for_lazy_eval.append({'invalid_data': msg_dict})

                    # elif warn.startswith('[Missing data]'):
                    #     suspended_errors_for_lazy_eval.append({'missing_data': msg_dict})

                    # elif warn.startswith('[Range value error]') and not self._reg.remediation_mode:
                    #     suspended_errors_for_lazy_eval.append({'anomalous_data': msg_dict})

        for input_source, ar, _ in ar_file_order:
            file_path = ar['file_name']

            input_source_dic = input_source.get()

            file_type = input_source_dic['file_type']
            content_subtype = input_source_dic['content_subtype']

            ignore_error = False if 'ignore_error' not in input_source_dic else input_source_dic['ignore_error']

            if file_type is None\
               or file_type in ('nm-aux-amb', 'nm-aux-cha', 'nm-aux-gro', 'nm-aux-pdb', 'nm-res-oth', 'nm-res-mr', 'nm-res-sax')\
               or file_type.startswith('nm-pea'):
                continue

            if self._reg.remediation_mode and os.path.exists(self.testPathWithSuffix(file_path, '-ignored', True)):
                continue

            file_name = input_source_dic['file_name']

            original_file_name = None
            if 'original_file_name' in input_source_dic:
                if input_source_dic['original_file_name'] is not None:
                    original_file_name = os.path.basename(input_source_dic['original_file_name'])
                if file_name != original_file_name and original_file_name is not None:
                    file_name = f"{original_file_name} ({file_name})"
            if original_file_name in EMPTY_VALUE and self._reg.internal_mode:
                original_file_name = file_name

            if file_type == 'nm-res-amb' and amberAtomNumberDict is None and 'has_comments' in ar and not ar['has_comments']:

                err = f"To verify AMBER restraint file {file_name!r}, AMBER topology file must be uploaded "\
                    "or Sander comments should be included in the AMBER restraint file."

                if self._reg.internal_mode:

                    self._reg.report.warning.appendDescription('missing_content',
                                                               {'file_name': file_name, 'description': err})

                    if self._reg.verbose:
                        self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {err}\n")

                else:

                    self._reg.report.error.appendDescription('missing_mandatory_content',
                                                             {'file_name': file_name, 'description': err})

                    if self._reg.verbose:
                        self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Error  - {err}\n")

                    continue

            if file_type == 'nm-res-cha' and charmmAtomNumberDict is None and not has_aux_cha:

                err = f"CHARMM topology file (aka. CRD or CHARM CARD) must be uploaded "\
                    f"to verify CHARMM restraint file {file_name!r}."

                if self._reg.internal_mode:

                    self._reg.report.warning.appendDescription('missing_content',
                                                               {'file_name': file_name, 'description': err})

                    if self._reg.verbose:
                        self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {err}\n")

                else:

                    self._reg.report.error.appendDescription('missing_mandatory_content',
                                                             {'file_name': file_name, 'description': err})

                    if self._reg.verbose:
                        self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Error  - {err}\n")

                    continue

            if file_type == 'nm-res-gro' and gromacsAtomNumberDict is None and not has_aux_gro:

                err = f"GROMACS topology file must be uploaded to verify GROMACS restraint file {file_name!r}."

                if self._reg.internal_mode:

                    self._reg.report.warning.appendDescription('missing_content',
                                                               {'file_name': file_name, 'description': err})

                    if self._reg.verbose:
                        self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Warning  - {err}\n")

                else:

                    self._reg.report.error.appendDescription('missing_mandatory_content',
                                                             {'file_name': file_name, 'description': err})

                    if self._reg.verbose:
                        self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Error  - {err}\n")

                    continue

            if content_subtype is None or len(content_subtype) == 0:
                continue

            if 'is_valid' not in ar or not ar['is_valid']:
                continue

            _reasons = reasons_dict.get(file_type)
            if _reasons is not None\
               and 'label_seq_scheme' not in _reasons\
               and 'global_auth_sequence_offset' not in _reasons\
               and 'chain_id_remap' not in _reasons\
               and file_type != 'nm-res-amb'\
               and 'dist_restraint' not in content_subtype:
                _reasons = None

            reasons = _reasons

            suspended_errors_for_lazy_eval.clear()

            if file_type == 'nm-res-amb':
                reader = AmberMRReader(self._reg.verbose, self._reg.log,
                                       self._reg.representative_model_id,
                                       self._reg.representative_alt_id,
                                       self._reg.mr_atom_name_mapping,
                                       self._reg.cR, self._reg.caC,
                                       self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                       amberAtomNumberDict, _amberAtomNumberDict,
                                       reasons)
                reader.setInternalMode(self._reg.internal_mode and derived_from_public_mr)

                _list_id_counter = copy.copy(self._reg.list_id_counter)

                listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                              createSfDict=create_sf_dict, originalFileName=original_file_name,
                                              listIdCounter=self._reg.list_id_counter,
                                              entryId=self._reg.entry_id)

                if listener is not None:
                    reasons = reader.getReasons()

                    if reasons is not None and _reasons is not None and listener.warningMessage is not None\
                       and len(listener.warningMessage) > 0:
                        deal_res_warn_message_for_lazy_eval(file_name, listener)

                        reader = AmberMRReader(self._reg.verbose, self._reg.log,
                                               self._reg.representative_model_id,
                                               self._reg.representative_alt_id,
                                               self._reg.mr_atom_name_mapping,
                                               self._reg.cR, self._reg.caC,
                                               self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                               amberAtomNumberDict, _amberAtomNumberDict,
                                               None)
                        reader.setInternalMode(self._reg.internal_mode and derived_from_public_mr)

                        listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                                      createSfDict=create_sf_dict, originalFileName=original_file_name,
                                                      listIdCounter=_list_id_counter,
                                                      entryId=self._reg.entry_id)

                        if listener is not None:
                            reasons = reader.getReasons()

                    if reasons is not None:

                        if 'dist_restraint' in content_subtype.keys():
                            reasons_dict[file_type] = reasons

                    deal_res_warn_message(file_name, listener, ignore_error)

                    cur_dict = listener.getAtomNumberDict()
                    if cur_dict is not None:
                        if len(_amberAtomNumberDict) == 0:
                            _amberAtomNumberDict = cur_dict
                        else:
                            for k, v in cur_dict.items():
                                if k not in _amberAtomNumberDict:
                                    _amberAtomNumberDict[k] = v

                    poly_seq = listener.getPolymerSequence()
                    if poly_seq is not None:
                        input_source.setItemValue('polymer_sequence', poly_seq)
                        poly_seq_set.append(poly_seq)

                    seq_align = listener.getSequenceAlignment()
                    if seq_align is not None:
                        self._reg.report.sequence_alignment.setItemValue('model_poly_seq_vs_mr_restraint', seq_align)

                    if create_sf_dict:
                        if len(listener.getContentSubtype()) == 0 and not ignore_error:
                            err = f"Failed to validate the restraint file (AMBER) {file_name!r}."

                            self._reg.report.error.appendDescription('internal_error',
                                                                     f"+{self.__class_name__}.validateLegacyMr() "
                                                                     "++ Error  - " + err)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Error  - {err}\n")

                        self._reg.list_id_counter, sf_dict = listener.getSfDict()
                        if sf_dict is not None:
                            for k, v in sf_dict.items():
                                content_subtype = contentSubtypeOf(k[0])
                                if content_subtype not in self._reg.mr_sf_dict_holder:
                                    self._reg.mr_sf_dict_holder[content_subtype] = []
                                for sf in v:
                                    if sf not in self._reg.mr_sf_dict_holder[content_subtype]:
                                        self._reg.mr_sf_dict_holder[content_subtype].append(sf)

            elif file_type in LEGACY_MR_READERS:
                spec = LEGACY_MR_READERS[file_type]

                parsed, listener = self._parseLegacyMr(
                    spec, file_path, file_name, original_file_name, create_sf_dict,
                    deal_res_warn_message_for_lazy_eval)

                if parsed:
                    deal_res_warn_message(file_name, listener, ignore_error)

                    poly_seq = listener.getPolymerSequence()
                    if poly_seq is not None:
                        input_source.setItemValue('polymer_sequence', poly_seq)
                        poly_seq_set.append(poly_seq)

                    seq_align = listener.getSequenceAlignment()
                    if seq_align is not None:
                        self._reg.report.sequence_alignment.setItemValue('model_poly_seq_vs_mr_restraint', seq_align)

                    if create_sf_dict:
                        if len(listener.getContentSubtype()) == 0 and not ignore_error:
                            label = spec[1]
                            err = f"Failed to validate the restraint file ({label}) {file_name!r}."

                            self._reg.report.error.appendDescription('internal_error',
                                                                     f"+{self.__class_name__}.validateLegacyMr() "
                                                                     "++ Error  - " + err)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Error  - {err}\n")

                        self._reg.list_id_counter, sf_dict = listener.getSfDict()
                        if sf_dict is not None:
                            for k, v in sf_dict.items():
                                content_subtype = contentSubtypeOf(k[0])
                                if content_subtype not in self._reg.mr_sf_dict_holder:
                                    self._reg.mr_sf_dict_holder[content_subtype] = []
                                for sf in v:
                                    if sf not in self._reg.mr_sf_dict_holder[content_subtype]:
                                        self._reg.mr_sf_dict_holder[content_subtype].append(sf)

            elif file_type == 'nm-res-cha':
                reader = CharmmMRReader(self._reg.verbose, self._reg.log,
                                        self._reg.representative_model_id,
                                        self._reg.representative_alt_id,
                                        self._reg.mr_atom_name_mapping,
                                        self._reg.cR, self._reg.caC,
                                        self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                        charmmAtomNumberDict,
                                        reasons)
                reader.setInternalMode(self._reg.internal_mode and derived_from_public_mr)
                reader.setNmrVsModel(nmr_vs_model)
                if file_path in self._reg.sll_pred_forced:
                    reader.setSllPredMode(True)

                _list_id_counter = copy.copy(self._reg.list_id_counter)
                __list_id_counter = copy.copy(self._reg.list_id_counter)

                listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                              createSfDict=create_sf_dict, originalFileName=original_file_name,
                                              listIdCounter=self._reg.list_id_counter,
                                              entryId=self._reg.entry_id)

                if listener is not None:
                    reasons = listener.getReasonsForReparsing()

                    if None not in (reasons, _reasons):

                        reader = CharmmMRReader(self._reg.verbose, self._reg.log,
                                                self._reg.representative_model_id,
                                                self._reg.representative_alt_id,
                                                self._reg.mr_atom_name_mapping,
                                                self._reg.cR, self._reg.caC,
                                                self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                                charmmAtomNumberDict,
                                                None)
                        reader.setInternalMode(self._reg.internal_mode and derived_from_public_mr)
                        reader.setNmrVsModel(nmr_vs_model)
                        if file_path in self._reg.sll_pred_forced:
                            reader.setSllPredMode(True)

                        listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                                      createSfDict=create_sf_dict, originalFileName=original_file_name,
                                                      listIdCounter=_list_id_counter,
                                                      entryId=self._reg.entry_id)

                        if listener is not None:
                            reasons = listener.getReasonsForReparsing()

                    if reasons is not None:
                        deal_res_warn_message_for_lazy_eval(file_name, listener)

                        if 'dist_restraint' in content_subtype.keys():
                            reasons_dict[file_type] = reasons

                        if 'model_chain_id_ext' in reasons:
                            self._reg.auth_asym_ids_with_chem_exch.update(reasons['model_chain_id_ext'])
                        if 'chain_id_clone' in reasons:
                            self._reg.auth_seq_ids_with_chem_exch.update(reasons['chain_id_clone'])

                        reader = CharmmMRReader(self._reg.verbose, self._reg.log,
                                                self._reg.representative_model_id,
                                                self._reg.representative_alt_id,
                                                self._reg.mr_atom_name_mapping,
                                                self._reg.cR, self._reg.caC,
                                                self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                                charmmAtomNumberDict,
                                                reasons)
                        reader.setInternalMode(self._reg.internal_mode and derived_from_public_mr)
                        reader.setNmrVsModel(nmr_vs_model)
                        if file_path in self._reg.sll_pred_forced:
                            reader.setSllPredMode(True)

                        listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                                      createSfDict=create_sf_dict, originalFileName=original_file_name,
                                                      listIdCounter=__list_id_counter,
                                                      entryId=self._reg.entry_id)

                    deal_res_warn_message(file_name, listener, ignore_error)

                    poly_seq = listener.getPolymerSequence()
                    if poly_seq is not None:
                        input_source.setItemValue('polymer_sequence', poly_seq)
                        poly_seq_set.append(poly_seq)

                    seq_align = listener.getSequenceAlignment()
                    if seq_align is not None:
                        self._reg.report.sequence_alignment.setItemValue('model_poly_seq_vs_mr_restraint', seq_align)

                    if create_sf_dict:
                        if len(listener.getContentSubtype()) == 0 and not ignore_error:
                            err = f"Failed to validate the restraint file (CHARMM) {file_name!r}."

                            self._reg.report.error.appendDescription('internal_error',
                                                                     f"+{self.__class_name__}.validateLegacyMr() "
                                                                     "++ Error  - " + err)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Error  - {err}\n")

                        self._reg.list_id_counter, sf_dict = listener.getSfDict()
                        if sf_dict is not None:
                            for k, v in sf_dict.items():
                                content_subtype = contentSubtypeOf(k[0])
                                if content_subtype not in self._reg.mr_sf_dict_holder:
                                    self._reg.mr_sf_dict_holder[content_subtype] = []
                                for sf in v:
                                    if sf not in self._reg.mr_sf_dict_holder[content_subtype]:
                                        self._reg.mr_sf_dict_holder[content_subtype].append(sf)

            elif file_type == 'nm-res-cns':
                reader = CnsMRReader(self._reg.verbose, self._reg.log,
                                     self._reg.representative_model_id,
                                     self._reg.representative_alt_id,
                                     self._reg.mr_atom_name_mapping,
                                     self._reg.cR, self._reg.caC,
                                     self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                     reasons)
                reader.setInternalMode(self._reg.internal_mode and derived_from_public_mr)
                reader.setNmrVsModel(nmr_vs_model)
                if file_path in self._reg.sll_pred_forced:
                    reader.setSllPredMode(True)

                _list_id_counter = copy.copy(self._reg.list_id_counter)
                __list_id_counter = copy.copy(self._reg.list_id_counter)

                listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                              createSfDict=create_sf_dict, originalFileName=original_file_name,
                                              listIdCounter=self._reg.list_id_counter,
                                              entryId=self._reg.entry_id)

                if listener is not None:
                    reasons = listener.getReasonsForReparsing()

                    if None not in (reasons, _reasons):

                        reader = CnsMRReader(self._reg.verbose, self._reg.log,
                                             self._reg.representative_model_id,
                                             self._reg.representative_alt_id,
                                             self._reg.mr_atom_name_mapping,
                                             self._reg.cR, self._reg.caC,
                                             self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                             None)
                        reader.setInternalMode(self._reg.internal_mode and derived_from_public_mr)
                        reader.setNmrVsModel(nmr_vs_model)
                        if file_path in self._reg.sll_pred_forced:
                            reader.setSllPredMode(True)

                        listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                                      createSfDict=create_sf_dict, originalFileName=original_file_name,
                                                      listIdCounter=_list_id_counter,
                                                      entryId=self._reg.entry_id)

                        if listener is not None:
                            reasons = listener.getReasonsForReparsing()

                    if reasons is not None:
                        deal_res_warn_message_for_lazy_eval(file_name, listener)

                        if 'dist_restraint' in content_subtype.keys():
                            reasons_dict[file_type] = reasons

                        if 'model_chain_id_ext' in reasons:
                            self._reg.auth_asym_ids_with_chem_exch.update(reasons['model_chain_id_ext'])
                        if 'chain_id_clone' in reasons:
                            self._reg.auth_seq_ids_with_chem_exch.update(reasons['chain_id_clone'])

                        reader = CnsMRReader(self._reg.verbose, self._reg.log,
                                             self._reg.representative_model_id,
                                             self._reg.representative_alt_id,
                                             self._reg.mr_atom_name_mapping,
                                             self._reg.cR, self._reg.caC,
                                             self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                             reasons)
                        reader.setInternalMode(self._reg.internal_mode and derived_from_public_mr)
                        reader.setNmrVsModel(nmr_vs_model)
                        if file_path in self._reg.sll_pred_forced:
                            reader.setSllPredMode(True)

                        listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                                      createSfDict=create_sf_dict, originalFileName=original_file_name,
                                                      listIdCounter=__list_id_counter,
                                                      entryId=self._reg.entry_id)

                    deal_res_warn_message(file_name, listener, ignore_error)

                    poly_seq = listener.getPolymerSequence()
                    if poly_seq is not None:
                        input_source.setItemValue('polymer_sequence', poly_seq)
                        poly_seq_set.append(poly_seq)

                    seq_align = listener.getSequenceAlignment()
                    if seq_align is not None:
                        self._reg.report.sequence_alignment.setItemValue('model_poly_seq_vs_mr_restraint', seq_align)

                    if create_sf_dict:
                        if len(listener.getContentSubtype()) == 0 and not ignore_error:
                            err = f"Failed to validate the restraint file (CNS) {file_name!r}."

                            self._reg.report.error.appendDescription('internal_error',
                                                                     f"+{self.__class_name__}.validateLegacyMr() "
                                                                     "++ Error  - " + err)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Error  - {err}\n")

                        self._reg.list_id_counter, sf_dict = listener.getSfDict()
                        if sf_dict is not None:
                            for k, v in sf_dict.items():
                                content_subtype = contentSubtypeOf(k[0])
                                if content_subtype not in self._reg.mr_sf_dict_holder:
                                    self._reg.mr_sf_dict_holder[content_subtype] = []
                                for sf in v:
                                    if sf not in self._reg.mr_sf_dict_holder[content_subtype]:
                                        self._reg.mr_sf_dict_holder[content_subtype].append(sf)

            elif file_type == 'nm-res-cya':
                has_dist_restraint = 'dist_restraint' in content_subtype

                upl_or_lol = None
                if has_dist_restraint:
                    dist_type = ar['dist_type']
                    if cyanaLolDistRest == 0 and dist_type == 'upl':
                        upl_or_lol = 'upl_only'
                    elif cyanaUplDistRest == 0 and dist_type == 'lol':
                        upl_or_lol = 'lol_only'
                    elif dist_type == 'upl':
                        upl_or_lol = 'upl_w_lol'
                    elif dist_type == 'lol':
                        upl_or_lol = 'lol_w_upl'
                    else:
                        upl_or_lol = None

                cya_file_ext = self._reg.dpS.retrieveOriginalFileExtensionOfCyanaMrFile() if self._reg.dpS is not None else None

                reader = CyanaMRReader(self._reg.verbose, self._reg.log,
                                       self._reg.representative_model_id,
                                       self._reg.representative_alt_id,
                                       self._reg.mr_atom_name_mapping,
                                       self._reg.cR, self._reg.caC,
                                       self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                       reasons, upl_or_lol, cya_file_ext)
                reader.setRemediateMode(self._reg.remediation_mode)

                _list_id_counter = copy.copy(self._reg.list_id_counter)
                __list_id_counter = copy.copy(self._reg.list_id_counter)

                listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                              createSfDict=create_sf_dict, originalFileName=original_file_name,
                                              listIdCounter=self._reg.list_id_counter,
                                              entryId=self._reg.entry_id)

                if listener is not None:
                    reasons = listener.getReasonsForReparsing()

                    if None not in (reasons, _reasons):

                        reader = CyanaMRReader(self._reg.verbose, self._reg.log,
                                               self._reg.representative_model_id,
                                               self._reg.representative_alt_id,
                                               self._reg.mr_atom_name_mapping,
                                               self._reg.cR, self._reg.caC,
                                               self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                               None, upl_or_lol, cya_file_ext)
                        reader.setRemediateMode(self._reg.remediation_mode)

                        listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                                      createSfDict=create_sf_dict, originalFileName=original_file_name,
                                                      listIdCounter=_list_id_counter,
                                                      entryId=self._reg.entry_id)

                        if listener is not None:
                            reasons = listener.getReasonsForReparsing()

                    if reasons is not None:
                        deal_res_warn_message_for_lazy_eval(file_name, listener)

                        if 'dist_restraint' in content_subtype.keys():
                            reasons_dict[file_type] = reasons

                        if 'model_chain_id_ext' in reasons:
                            self._reg.auth_asym_ids_with_chem_exch.update(reasons['model_chain_id_ext'])
                        if 'chain_id_clone' in reasons:
                            self._reg.auth_seq_ids_with_chem_exch.update(reasons['chain_id_clone'])

                        reader = CyanaMRReader(self._reg.verbose, self._reg.log,
                                               self._reg.representative_model_id,
                                               self._reg.representative_alt_id,
                                               self._reg.mr_atom_name_mapping,
                                               self._reg.cR, self._reg.caC,
                                               self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                               reasons, upl_or_lol, cya_file_ext)
                        reader.setRemediateMode(self._reg.remediation_mode)

                        listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                                      createSfDict=create_sf_dict, originalFileName=original_file_name,
                                                      listIdCounter=__list_id_counter,
                                                      entryId=self._reg.entry_id)

                    deal_res_warn_message(file_name, listener, ignore_error)

                    poly_seq = listener.getPolymerSequence()
                    if poly_seq is not None:
                        input_source.setItemValue('polymer_sequence', poly_seq)
                        poly_seq_set.append(poly_seq)

                    seq_align = listener.getSequenceAlignment()
                    if seq_align is not None:
                        self._reg.report.sequence_alignment.setItemValue('model_poly_seq_vs_mr_restraint', seq_align)

                    # support content subtype change during MR validation with the coordinates
                    input_source.setItemValue('content_subtype', listener.getContentSubtype())

                    if create_sf_dict:
                        if len(listener.getContentSubtype()) == 0 and not ignore_error:
                            err = f"Failed to validate the restraint file (CYANA) {file_name!r}."

                            self._reg.report.error.appendDescription('internal_error',
                                                                     f"+{self.__class_name__}.validateLegacyMr() "
                                                                     "++ Error  - " + err)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Error  - {err}\n")

                        self._reg.list_id_counter, sf_dict = listener.getSfDict()
                        if sf_dict is not None:
                            for k, v in sf_dict.items():
                                content_subtype = contentSubtypeOf(k[0])
                                if content_subtype not in self._reg.mr_sf_dict_holder:
                                    self._reg.mr_sf_dict_holder[content_subtype] = []
                                for sf in v:
                                    if sf not in self._reg.mr_sf_dict_holder[content_subtype]:
                                        self._reg.mr_sf_dict_holder[content_subtype].append(sf)

            elif file_type == 'nm-res-gro':
                reader = GromacsMRReader(self._reg.verbose, self._reg.log,
                                         self._reg.representative_model_id,
                                         self._reg.representative_alt_id,
                                         self._reg.mr_atom_name_mapping,
                                         self._reg.cR, self._reg.caC,
                                         self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                         gromacsAtomNumberDict)

                listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                              createSfDict=create_sf_dict, originalFileName=original_file_name,
                                              listIdCounter=self._reg.list_id_counter,
                                              entryId=self._reg.entry_id)

                if listener is not None:
                    deal_res_warn_message(file_name, listener, ignore_error)

                    poly_seq = listener.getPolymerSequence()
                    if poly_seq is not None:
                        input_source.setItemValue('polymer_sequence', poly_seq)
                        poly_seq_set.append(poly_seq)

                    seq_align = listener.getSequenceAlignment()
                    if seq_align is not None:
                        self._reg.report.sequence_alignment.setItemValue('model_poly_seq_vs_mr_restraint', seq_align)

                    if create_sf_dict:
                        if len(listener.getContentSubtype()) == 0 and not ignore_error:
                            err = f"Failed to validate the restraint file (GROMACS) {file_name!r}."

                            self._reg.report.error.appendDescription('internal_error',
                                                                     f"+{self.__class_name__}.validateLegacyMr() "
                                                                     "++ Error  - " + err)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Error  - {err}\n")

                        self._reg.list_id_counter, sf_dict = listener.getSfDict()
                        if sf_dict is not None:
                            for k, v in sf_dict.items():
                                content_subtype = contentSubtypeOf(k[0])
                                if content_subtype not in self._reg.mr_sf_dict_holder:
                                    self._reg.mr_sf_dict_holder[content_subtype] = []
                                for sf in v:
                                    if sf not in self._reg.mr_sf_dict_holder[content_subtype]:
                                        self._reg.mr_sf_dict_holder[content_subtype].append(sf)

            elif file_type == 'nm-res-noa':
                reader = CyanaNOAReader(self._reg.verbose, self._reg.log,
                                        self._reg.representative_model_id,
                                        self._reg.representative_alt_id,
                                        self._reg.mr_atom_name_mapping,
                                        self._reg.cR, self._reg.caC,
                                        self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                        reasons)

                _list_id_counter = copy.copy(self._reg.list_id_counter)
                __list_id_counter = copy.copy(self._reg.list_id_counter)

                listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                              createSfDict=create_sf_dict, originalFileName=original_file_name,
                                              listIdCounter=self._reg.list_id_counter,
                                              entryId=self._reg.entry_id)

                if listener is not None:
                    reasons = listener.getReasonsForReparsing()

                    if None not in (reasons, _reasons):

                        reader = CyanaNOAReader(self._reg.verbose, self._reg.log,
                                                self._reg.representative_model_id,
                                                self._reg.representative_alt_id,
                                                self._reg.mr_atom_name_mapping,
                                                self._reg.cR, self._reg.caC,
                                                self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                                None)

                        listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                                      createSfDict=create_sf_dict, originalFileName=original_file_name,
                                                      listIdCounter=_list_id_counter,
                                                      entryId=self._reg.entry_id)

                        if listener is not None:
                            reasons = listener.getReasonsForReparsing()

                    if reasons is not None:
                        deal_res_warn_message_for_lazy_eval(file_name, listener)

                        if 'dist_restraint' in content_subtype.keys():
                            reasons_dict[file_type] = reasons

                        if 'model_chain_id_ext' in reasons:
                            self._reg.auth_asym_ids_with_chem_exch.update(reasons['model_chain_id_ext'])
                        if 'chain_id_clone' in reasons:
                            self._reg.auth_seq_ids_with_chem_exch.update(reasons['chain_id_clone'])

                        reader = CyanaNOAReader(self._reg.verbose, self._reg.log,
                                                self._reg.representative_model_id,
                                                self._reg.representative_alt_id,
                                                self._reg.mr_atom_name_mapping,
                                                self._reg.cR, self._reg.caC,
                                                self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                                reasons)

                        listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                                      createSfDict=create_sf_dict, originalFileName=original_file_name,
                                                      listIdCounter=__list_id_counter,
                                                      entryId=self._reg.entry_id)

                    deal_res_warn_message(file_name, listener, ignore_error)

                    poly_seq = listener.getPolymerSequence()
                    if poly_seq is not None:
                        input_source.setItemValue('polymer_sequence', poly_seq)
                        poly_seq_set.append(poly_seq)

                    seq_align = listener.getSequenceAlignment()
                    if seq_align is not None:
                        self._reg.report.sequence_alignment.setItemValue('model_poly_seq_vs_mr_restraint', seq_align)

                    # support content subtype change during MR validation with the coordinates
                    input_source.setItemValue('content_subtype', listener.getContentSubtype())

                    if create_sf_dict:
                        if len(listener.getContentSubtype()) == 0 and not ignore_error:
                            err = f"Failed to validate the restraint file (CYANA NOA) {file_name!r}."

                            self._reg.report.error.appendDescription('internal_error',
                                                                     f"+{self.__class_name__}.validateLegacyMr() "
                                                                     "++ Error  - " + err)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Error  - {err}\n")

                        self._reg.list_id_counter, sf_dict = listener.getSfDict()
                        if sf_dict is not None:
                            for k, v in sf_dict.items():
                                content_subtype = contentSubtypeOf(k[0])
                                if content_subtype not in self._reg.mr_sf_dict_holder:
                                    self._reg.mr_sf_dict_holder[content_subtype] = []
                                for sf in v:
                                    if sf not in self._reg.mr_sf_dict_holder[content_subtype]:
                                        self._reg.mr_sf_dict_holder[content_subtype].append(sf)

            elif file_type == 'nm-res-ros':
                reader = RosettaMRReader(self._reg.verbose, self._reg.log,
                                         self._reg.representative_model_id,
                                         self._reg.representative_alt_id,
                                         self._reg.mr_atom_name_mapping,
                                         self._reg.cR, self._reg.caC,
                                         self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                         reasons)
                reader.setRemediateMode(self._reg.remediation_mode)

                _list_id_counter = copy.copy(self._reg.list_id_counter)
                __list_id_counter = copy.copy(self._reg.list_id_counter)

                listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                              createSfDict=create_sf_dict, originalFileName=original_file_name,
                                              listIdCounter=self._reg.list_id_counter,
                                              entryId=self._reg.entry_id)

                if listener is not None:
                    reasons = listener.getReasonsForReparsing()

                    if None not in (reasons, _reasons):

                        reader = RosettaMRReader(self._reg.verbose, self._reg.log,
                                                 self._reg.representative_model_id,
                                                 self._reg.representative_alt_id,
                                                 self._reg.mr_atom_name_mapping,
                                                 self._reg.cR, self._reg.caC,
                                                 self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                                 None)
                        reader.setRemediateMode(self._reg.remediation_mode)

                        listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                                      createSfDict=create_sf_dict, originalFileName=original_file_name,
                                                      listIdCounter=_list_id_counter,
                                                      entryId=self._reg.entry_id)

                        if listener is not None:
                            reasons = listener.getReasonsForReparsing()

                    if reasons is not None:
                        deal_res_warn_message_for_lazy_eval(file_name, listener)

                        if 'dist_restraint' in content_subtype.keys():
                            reasons_dict[file_type] = reasons

                        if 'model_chain_id_ext' in reasons:
                            self._reg.auth_asym_ids_with_chem_exch.update(reasons['model_chain_id_ext'])
                        if 'chain_id_clone' in reasons:
                            self._reg.auth_seq_ids_with_chem_exch.update(reasons['chain_id_clone'])

                        reader = RosettaMRReader(self._reg.verbose, self._reg.log,
                                                 self._reg.representative_model_id,
                                                 self._reg.representative_alt_id,
                                                 self._reg.mr_atom_name_mapping,
                                                 self._reg.cR, self._reg.caC,
                                                 self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                                 reasons)
                        reader.setRemediateMode(self._reg.remediation_mode)

                        listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                                      createSfDict=create_sf_dict, originalFileName=original_file_name,
                                                      listIdCounter=__list_id_counter,
                                                      entryId=self._reg.entry_id)

                    deal_res_warn_message(file_name, listener, ignore_error)

                    poly_seq = listener.getPolymerSequence()
                    if poly_seq is not None:
                        input_source.setItemValue('polymer_sequence', poly_seq)
                        poly_seq_set.append(poly_seq)

                    seq_align = listener.getSequenceAlignment()
                    if seq_align is not None:
                        self._reg.report.sequence_alignment.setItemValue('model_poly_seq_vs_mr_restraint', seq_align)

                    if create_sf_dict:
                        if len(listener.getContentSubtype()) == 0 and not ignore_error:
                            err = f"Failed to validate the restraint file (ROSETTA) {file_name!r}."

                            self._reg.report.error.appendDescription('internal_error',
                                                                     f"+{self.__class_name__}.validateLegacyMr() "
                                                                     "++ Error  - " + err)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Error  - {err}\n")

                        self._reg.list_id_counter, sf_dict = listener.getSfDict()
                        if sf_dict is not None:
                            for k, v in sf_dict.items():
                                content_subtype = contentSubtypeOf(k[0])
                                if content_subtype not in self._reg.mr_sf_dict_holder:
                                    self._reg.mr_sf_dict_holder[content_subtype] = []
                                for sf in v:
                                    if sf not in self._reg.mr_sf_dict_holder[content_subtype]:
                                        self._reg.mr_sf_dict_holder[content_subtype].append(sf)

            elif file_type == 'nm-res-sch':
                reader = SchrodingerMRReader(self._reg.verbose, self._reg.log,
                                             self._reg.representative_model_id,
                                             self._reg.representative_alt_id,
                                             self._reg.mr_atom_name_mapping,
                                             self._reg.cR, self._reg.caC,
                                             self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                             pdbAtomNumberDict, reasons)
                reader.setInternalMode(self._reg.internal_mode and derived_from_public_mr)
                reader.setNmrVsModel(nmr_vs_model)

                _list_id_counter = copy.copy(self._reg.list_id_counter)
                __list_id_counter = copy.copy(self._reg.list_id_counter)

                listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                              createSfDict=create_sf_dict, originalFileName=original_file_name,
                                              listIdCounter=self._reg.list_id_counter,
                                              entryId=self._reg.entry_id)

                if listener is not None:
                    reasons = listener.getReasonsForReparsing()

                    if None not in (reasons, _reasons):

                        reader = SchrodingerMRReader(self._reg.verbose, self._reg.log,
                                                     self._reg.representative_model_id,
                                                     self._reg.representative_alt_id,
                                                     self._reg.mr_atom_name_mapping,
                                                     self._reg.cR, self._reg.caC,
                                                     self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                                     pdbAtomNumberDict, None)
                        reader.setInternalMode(self._reg.internal_mode and derived_from_public_mr)
                        reader.setNmrVsModel(nmr_vs_model)

                        listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                                      createSfDict=create_sf_dict, originalFileName=original_file_name,
                                                      listIdCounter=_list_id_counter,
                                                      entryId=self._reg.entry_id)

                        if listener is not None:
                            reasons = listener.getReasonsForReparsing()

                    if reasons is not None:
                        deal_res_warn_message_for_lazy_eval(file_name, listener)

                        if 'dist_restraint' in content_subtype.keys():
                            reasons_dict[file_type] = reasons

                        if 'model_chain_id_ext' in reasons:
                            self._reg.auth_asym_ids_with_chem_exch.update(reasons['model_chain_id_ext'])
                        if 'chain_id_clone' in reasons:
                            self._reg.auth_seq_ids_with_chem_exch.update(reasons['chain_id_clone'])

                        reader = SchrodingerMRReader(self._reg.verbose, self._reg.log,
                                                     self._reg.representative_model_id,
                                                     self._reg.representative_alt_id,
                                                     self._reg.mr_atom_name_mapping,
                                                     self._reg.cR, self._reg.caC,
                                                     self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                                     pdbAtomNumberDict, reasons)
                        reader.setInternalMode(self._reg.internal_mode and derived_from_public_mr)
                        reader.setNmrVsModel(nmr_vs_model)

                        listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                                      createSfDict=create_sf_dict, originalFileName=original_file_name,
                                                      listIdCounter=__list_id_counter,
                                                      entryId=self._reg.entry_id)

                    deal_res_warn_message(file_name, listener, ignore_error)

                    poly_seq = listener.getPolymerSequence()
                    if poly_seq is not None:
                        input_source.setItemValue('polymer_sequence', poly_seq)
                        poly_seq_set.append(poly_seq)

                    seq_align = listener.getSequenceAlignment()
                    if seq_align is not None:
                        self._reg.report.sequence_alignment.setItemValue('model_poly_seq_vs_mr_restraint', seq_align)

                    if create_sf_dict:
                        if len(listener.getContentSubtype()) == 0 and not ignore_error:
                            err = f"Failed to validate the restraint file (SCHRODINGER/ASL) {file_name!r}."

                            self._reg.report.error.appendDescription('internal_error',
                                                                     f"+{self.__class_name__}.validateLegacyMr() "
                                                                     "++ Error  - " + err)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Error  - {err}\n")

                        self._reg.list_id_counter, sf_dict = listener.getSfDict()
                        if sf_dict is not None:
                            for k, v in sf_dict.items():
                                content_subtype = contentSubtypeOf(k[0])
                                if content_subtype not in self._reg.mr_sf_dict_holder:
                                    self._reg.mr_sf_dict_holder[content_subtype] = []
                                for sf in v:
                                    if sf not in self._reg.mr_sf_dict_holder[content_subtype]:
                                        self._reg.mr_sf_dict_holder[content_subtype].append(sf)

            elif file_type == 'nm-res-xpl':
                reader = XplorMRReader(self._reg.verbose, self._reg.log,
                                       self._reg.representative_model_id,
                                       self._reg.representative_alt_id,
                                       self._reg.mr_atom_name_mapping,
                                       self._reg.cR, self._reg.caC,
                                       self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                       reasons)
                reader.setRemediateMode(self._reg.remediation_mode and derived_from_public_mr)
                reader.setInternalMode(self._reg.internal_mode and derived_from_public_mr)
                reader.setNmrVsModel(nmr_vs_model)
                if file_path in self._reg.sll_pred_forced:
                    reader.setSllPredMode(True)

                _list_id_counter = copy.copy(self._reg.list_id_counter)
                __list_id_counter = copy.copy(self._reg.list_id_counter)

                listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                              createSfDict=create_sf_dict, originalFileName=original_file_name,
                                              listIdCounter=self._reg.list_id_counter,
                                              entryId=self._reg.entry_id)

                if listener is not None:
                    reasons = listener.getReasonsForReparsing()

                    if None not in (reasons, _reasons):

                        reader = XplorMRReader(self._reg.verbose, self._reg.log,
                                               self._reg.representative_model_id,
                                               self._reg.representative_alt_id,
                                               self._reg.mr_atom_name_mapping,
                                               self._reg.cR, self._reg.caC,
                                               self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                               None)
                        reader.setRemediateMode(self._reg.remediation_mode and derived_from_public_mr)
                        reader.setInternalMode(self._reg.internal_mode and derived_from_public_mr)
                        reader.setNmrVsModel(nmr_vs_model)
                        if file_path in self._reg.sll_pred_forced:
                            reader.setSllPredMode(True)

                        listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                                      createSfDict=create_sf_dict, originalFileName=original_file_name,
                                                      listIdCounter=_list_id_counter,
                                                      entryId=self._reg.entry_id)

                        if listener is not None:
                            reasons = listener.getReasonsForReparsing()

                    if reasons is not None:
                        deal_res_warn_message_for_lazy_eval(file_name, listener)

                        if 'dist_restraint' in content_subtype.keys():
                            reasons_dict[file_type] = reasons

                        if 'model_chain_id_ext' in reasons:
                            self._reg.auth_asym_ids_with_chem_exch.update(reasons['model_chain_id_ext'])
                        if 'chain_id_clone' in reasons:
                            self._reg.auth_seq_ids_with_chem_exch.update(reasons['chain_id_clone'])

                        reader = XplorMRReader(self._reg.verbose, self._reg.log,
                                               self._reg.representative_model_id,
                                               self._reg.representative_alt_id,
                                               self._reg.mr_atom_name_mapping,
                                               self._reg.cR, self._reg.caC,
                                               self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                               reasons)
                        reader.setRemediateMode(self._reg.remediation_mode and derived_from_public_mr)
                        reader.setInternalMode(self._reg.internal_mode and derived_from_public_mr)
                        reader.setNmrVsModel(nmr_vs_model)
                        if file_path in self._reg.sll_pred_forced:
                            reader.setSllPredMode(True)

                        listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                                      createSfDict=create_sf_dict, originalFileName=original_file_name,
                                                      listIdCounter=__list_id_counter,
                                                      entryId=self._reg.entry_id)

                    deal_res_warn_message(file_name, listener, ignore_error)

                    poly_seq = listener.getPolymerSequence()
                    if poly_seq is not None:
                        input_source.setItemValue('polymer_sequence', poly_seq)
                        poly_seq_set.append(poly_seq)

                    seq_align = listener.getSequenceAlignment()
                    if seq_align is not None:
                        self._reg.report.sequence_alignment.setItemValue('model_poly_seq_vs_mr_restraint', seq_align)

                    # support content subtype change during MR validation with the coordinates
                    input_source.setItemValue('content_subtype', listener.getContentSubtype())

                    if create_sf_dict:
                        if len(listener.getContentSubtype()) == 0 and not ignore_error:
                            err = f"Failed to validate the restraint file (XPLOR-NIH) {file_name!r}."

                            self._reg.report.error.appendDescription('internal_error',
                                                                     f"+{self.__class_name__}.validateLegacyMr() "
                                                                     "++ Error  - " + err)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyMr() ++ Error  - {err}\n")

                        self._reg.list_id_counter, sf_dict = listener.getSfDict()
                        if sf_dict is not None:
                            for k, v in sf_dict.items():
                                content_subtype = contentSubtypeOf(k[0])
                                if content_subtype not in self._reg.mr_sf_dict_holder:
                                    self._reg.mr_sf_dict_holder[content_subtype] = []
                                for sf in v:
                                    if sf not in self._reg.mr_sf_dict_holder[content_subtype]:
                                        self._reg.mr_sf_dict_holder[content_subtype].append(sf)

        if len(poly_seq_set) > 1:

            poly_seq_rst = None
            for idx, poly_seq in enumerate(poly_seq_set):
                if idx == 0:
                    poly_seq_rst = poly_seq
                    continue
                for ps in poly_seq:
                    chain_id = ps['chain_id']
                    for seq_id, comp_id in zip(ps['seq_id'], ps['comp_id']):
                        updatePolySeqRst(poly_seq_rst, chain_id, seq_id, comp_id)

            poly_seq_model = self._reg.caC['polymer_sequence']

            sortPolySeqRst(poly_seq_rst)

            file_type = 'nm-res-mr'

            seq_align, _ = alignPolymerSequence(self._reg.pA, poly_seq_model, poly_seq_rst, conservative=False)
            chain_assign, _ = assignPolymerSequence(self._reg.pA, self._reg.ccU, file_type,
                                                    poly_seq_model, poly_seq_rst, seq_align)

            if chain_assign is not None:

                if len(poly_seq_model) == len(poly_seq_rst):

                    chain_mapping = {}

                    for ca in chain_assign:
                        ref_chain_id = ca['ref_chain_id']
                        test_chain_id = ca['test_chain_id']

                        if ref_chain_id != test_chain_id:
                            chain_mapping[test_chain_id] = ref_chain_id

                    if len(chain_mapping) == len(poly_seq_model):

                        for ps in poly_seq_rst:
                            if ps['chain_id'] in chain_mapping:
                                ps['chain_id'] = chain_mapping[ps['chain_id']]

                        seq_align, _ = alignPolymerSequence(self._reg.pA, poly_seq_model, poly_seq_rst, conservative=False)
                        chain_assign, _ = assignPolymerSequence(self._reg.pA, self._reg.ccU, file_type,
                                                                poly_seq_model, poly_seq_rst, seq_align)

                    trimSequenceAlignment(seq_align, chain_assign)

            input_source.setItemValue('polymer_sequence', poly_seq_rst)

            self._reg.report.sequence_alignment.setItemValue('model_poly_seq_vs_mr_restraint', seq_align)

        if len(self._reg.nmr_ext_poly_seq) > 0:
            entity_assembly = self._reg.caC['entity_assembly']
            auth_chain_ids = list(set(d['auth_chain_id'] for d in self._reg.nmr_ext_poly_seq))
            for auth_chain_id in auth_chain_ids:
                try:
                    item = next(item for item in entity_assembly if auth_chain_id in item['auth_asym_id'].split(','))
                except StopIteration:
                    continue
                if item['entity_type'] == 'polymer':
                    poly_type = item['entity_poly_type']
                    if poly_type.startswith('polypeptide'):
                        unknown_residue = 'UNK'
                    elif any(True for comp_id in item['comp_id_set'] if comp_id in ('DA', 'DC', 'DG', 'DT'))\
                            and any(True for comp_id in item['comp_id_set'] if comp_id in ('A', 'C', 'G', 'U')):
                        unknown_residue = 'DN'
                    elif poly_type == 'polydeoxyribonucleotide':
                        unknown_residue = 'DN'
                    elif poly_type == 'polyribonucleotide':
                        unknown_residue = 'N'
                    else:
                        continue
                    ps = next(ps for ps in self._reg.caC['polymer_sequence'] if ps['auth_chain_id'] == auth_chain_id)
                    auth_seq_ids = [d['auth_seq_id'] for d in self._reg.nmr_ext_poly_seq if d['auth_chain_id'] == auth_chain_id]
                    auth_seq_ids.extend(list(filter(None, ps['auth_seq_id'])))
                    min_auth_seq_id = min(auth_seq_ids)
                    max_auth_seq_id = max(auth_seq_ids)
                    for auth_seq_id in range(min_auth_seq_id, max_auth_seq_id + 1):
                        if auth_seq_id not in ps['auth_seq_id']\
                           and not any(True for d in self._reg.nmr_ext_poly_seq
                                       if d['auth_chain_id'] == auth_chain_id and d['auth_seq_id'] == auth_seq_id):
                            self._reg.nmr_ext_poly_seq.append({'auth_chain_id': auth_chain_id,
                                                               'auth_seq_id': auth_seq_id,
                                                               'auth_comp_id': unknown_residue})

            self._reg.nmr_ext_poly_seq = sorted(self._reg.nmr_ext_poly_seq, key=itemgetter('auth_chain_id', 'auth_seq_id'))

        return not self._reg.report.isError()

    def validateSaxsMr(self) -> bool:
        """ Validate SAXS restraint files.
        """

        if self._reg.combined_mode:
            return True

        if AR_FILE_PATH_LIST_KEY not in self._reg.inputParamDict:
            return True

        content_subtype = 'saxs_restraint'

        if self._reg.list_id_counter is None:
            self._reg.list_id_counter = {}
        if self._reg.mr_sf_dict_holder is None:
            self._reg.mr_sf_dict_holder = {}

        if content_subtype not in self._reg.mr_sf_dict_holder:
            self._reg.mr_sf_dict_holder[content_subtype] = []

        fileListId = self._reg.file_path_list_len

        for ar in self._reg.inputParamDict[AR_FILE_PATH_LIST_KEY]:
            file_path = ar['file_name']

            input_source = self._reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_type = input_source_dic['file_type']

            fileListId += 1

            if file_type != 'nm-res-sax':
                continue

            file_name = input_source_dic['file_name']

            original_file_name = os.path.basename(file_name).replace('-corrected', '').replace('-selected-as-res-sax', '')
            if '-div_' in original_file_name:
                original_file_name = None
                # if 'original_file_name' in input_source_dic:
                #     if input_source_dic['original_file_name'] is not None:
                #         original_file_name = os.path.basename(input_source_dic['original_file_name'])

            sf_item = {}

            title = _row = None

            _q_value = 0.0

            lp_count = 0

            with open(file_path, 'r', encoding='utf-8') as ifh:
                for line in ifh:

                    line = ' '.join(line.split())

                    _line = line.split()

                    len_line = len(_line)

                    if len_line == 0:
                        continue

                    if line.startswith('#') or line.startswith('!'):

                        if len(line) == 1:
                            continue

                        __line = line[1:].split(' ')

                        if len(__line) == 0 or line[1].startswith('#') or line[1].startswith('!'):
                            continue

                        title_can = __line[0]
                        if len(title_can) == 0 and len(__line) > 1:
                            title_can = __line[1]

                        try:
                            float(title_can)
                            continue
                        except ValueError:
                            if len(title_can) > 1\
                               and '(' not in title_can and ')' not in title_can\
                               and '[' not in title_can and ']' not in title_can:
                                if len(title_can) > 0:
                                    title = title_can
                            _row = None

                        continue

                    if len_line != 3:
                        continue

                    try:

                        q_value = float(_line[0])
                        float(_line[1])
                        float(_line[2])

                        dstFunc = {'weight': '1.0',
                                   'target_value': _line[1].replace('E', 'e'),
                                   'target_value_uncertainty': _line[2].replace('E', 'e')}

                        if _q_value == 0.0:

                            if len(sf_item) > 0 and sf_item['id'] > 0:
                                self._reg.mr_sf_dict_holder[content_subtype].append(sf_item)
                                lp_count += 1

                            self._reg.list_id_counter =\
                                incListIdCounter(content_subtype, self._reg.list_id_counter, reduced=False)

                            list_id = self._reg.list_id_counter[content_subtype]

                            restraint_name = getRestraintName(content_subtype)

                            sf_framecode = restraint_name.replace(' ', '_').lower() + f'_{list_id}'

                            if original_file_name is not None:
                                title = original_file_name

                            sf = getSaveframe(content_subtype, sf_framecode, list_id, self._reg.entry_id, title, reduced=False)

                            _restraint_name = restraint_name.split()

                            sf_item = {'file_type': file_type, 'saveframe': sf, 'list_id': list_id,
                                       'id': 0, 'index_id': 0,
                                       'constraint_type': ' '.join(_restraint_name[:-1])}

                            lp = getLoop(content_subtype, reduced=False)

                            sf.add_loop(lp)
                            sf_item['loop'] = lp

                            if _row is not None:

                                sf_item['loop'].add_data(_row)

                                sf_item['id'] = sf_item['index_id'] = 1

                                _row = None

                        if q_value > _q_value:

                            sf_item['id'] += 1
                            sf_item['index_id'] += 1

                            row = getRow('saxs', sf_item['id'], sf_item['index_id'], None, None, _line[0].replace('E', 'e'),
                                         sf_item['list_id'], self._reg.entry_id, dstFunc, None, None, None, None, None)
                            sf_item['loop'].add_data(row)

                            _q_value = q_value

                        else:

                            _row = getRow('saxs', 1, 1, None, None, _line[0].replace('E', 'e'),
                                          sf_item['list_id'] + 1, self._reg.entry_id, dstFunc, None, None, None, None, None)

                            _q_value = 0.0

                    except ValueError:
                        continue

            if len(sf_item) > 0 and sf_item['id'] > 0:
                self._reg.mr_sf_dict_holder[content_subtype].append(sf_item)

                lp_count += 1
                input_source.setItemValue('content_subtype', {'saxs_restraint': lp_count})

        if len(self._reg.mr_sf_dict_holder[content_subtype]) == 0:
            del self._reg.mr_sf_dict_holder[content_subtype]

        return True
