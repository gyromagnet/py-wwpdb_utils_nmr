##
# File: NmrDpRemediationCsLoop.py
# Date: 17-Aug-2026
#
# Updates:
##
""" Remediation of the assigned chemical shift loop of NMR data.
    @author: Masashi Yokochi
"""
__docformat__ = "restructuredtext en"
__author__ = "Masashi Yokochi"
__email__ = "yokochi@protein.osaka-u.ac.jp"
__license__ = "Apache License 2.0"
__version__ = "5.4.4"

import collections
import copy
import functools
import itertools
import re
from operator import itemgetter
from typing import Union

import pynmrstar

try:
    from wwpdb.utils.nmr.NmrDpConstant import (DATA_ITEMS,
                                               AUX_LP_CATEGORIES,
                                               LOW_SEQ_COVERAGE,
                                               EMPTY_VALUE,
                                               STD_MON_DICT,
                                               PROTON_BEGIN_CODE,
                                               PSE_PRO_BEGIN_CODE,
                                               AMINO_PROTON_CODE,
                                               ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS,
                                               ALLOWED_AMBIGUITY_CODES,
                                               GLOBAL_OFFSET_ATTEMPT,
                                               PERIPH_OFFSET_ATTEMPT,
                                               CONCAT_SEQ_ID_INS_CODE_PAT,
                                               REPRESENTATIVE_ASYM_ID)
    from wwpdb.utils.nmr.AlignUtil import (letterToDigit,
                                           alignPolymerSequence,
                                           assignPolymerSequence)
    from wwpdb.utils.nmr.CifToNmrStar import (has_key_value,
                                              get_first_sf_tag,
                                              set_sf_tag)
    from wwpdb.utils.nmr.NmrVrptUtility import write_as_pickle
    from wwpdb.utils.nmr.mr.ParserListenerUtil import (translateToStdResName,
                                                       translateToStdAtomName)
    from wwpdb.utils.nmr.NmrDpRemediationBase import (NmrDpRemediationBase,
                                                      CCA_ENT_ASM_ID,
                                                      CCA_ENTITY_ID,
                                                      CCA_COMP_IDX,
                                                      CCA_SEQ_ID,
                                                      CCA_COMP_ID,
                                                      CCA_AUTH_ASYM)
except ImportError:
    from nmr.NmrDpConstant import (DATA_ITEMS,
                                   AUX_LP_CATEGORIES,
                                   LOW_SEQ_COVERAGE,
                                   EMPTY_VALUE,
                                   STD_MON_DICT,
                                   PROTON_BEGIN_CODE,
                                   PSE_PRO_BEGIN_CODE,
                                   AMINO_PROTON_CODE,
                                   ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS,
                                   ALLOWED_AMBIGUITY_CODES,
                                   GLOBAL_OFFSET_ATTEMPT,
                                   PERIPH_OFFSET_ATTEMPT,
                                   CONCAT_SEQ_ID_INS_CODE_PAT,
                                   REPRESENTATIVE_ASYM_ID)
    from nmr.AlignUtil import (letterToDigit,
                               alignPolymerSequence,
                               assignPolymerSequence)
    from nmr.CifToNmrStar import (has_key_value,
                                  get_first_sf_tag,
                                  set_sf_tag)
    from nmr.NmrVrptUtility import write_as_pickle
    from nmr.mr.ParserListenerUtil import (translateToStdResName,
                                           translateToStdAtomName)
    from nmr.NmrDpRemediationBase import (NmrDpRemediationBase,
                                          CCA_ENT_ASM_ID,
                                          CCA_ENTITY_ID,
                                          CCA_COMP_IDX,
                                          CCA_SEQ_ID,
                                          CCA_COMP_ID,
                                          CCA_AUTH_ASYM)


# Column positions of the rows remediateCsLoop() builds. They are named so that a
# transposition in a multi-column assignment is visible: _row[CS_ORIG_SEQ] cannot
# be mistaken for _row[CS_ORIG_COMP] the way _row[21] and _row[22] can. The tag
# tuples below are what the loops are actually built from, so the indices and the
# column order have a single definition and cannot drift apart.

CS_TAGS = ('ID',
           'Entity_assembly_ID',
           'Entity_ID',
           'Comp_index_ID',
           'Seq_ID',
           'Comp_ID',
           'Atom_ID',
           'Atom_type',
           'Atom_isotope_number',
           'Val',
           'Val_err',
           'Assign_fig_of_merit',
           'Ambiguity_code',
           'Ambiguity_set_ID',
           'Occupancy',
           'Resonance_ID',
           'Auth_asym_ID',
           'Auth_seq_ID',
           'Auth_comp_ID',
           'Auth_atom_ID',
           'Original_PDB_strand_ID',
           'Original_PDB_residue_no',
           'Original_PDB_residue_name',
           'Original_PDB_atom_name',
           'Details',
           'Entry_ID',
           'Assigned_chem_shift_list_ID')

# positions within CS_TAGS
CS_ID = 0  # ID
CS_ENT_ASM_ID = 1  # Entity_assembly_ID
CS_ENTITY_ID = 2  # Entity_ID
CS_COMP_IDX = 3  # Comp_index_ID
CS_SEQ_ID = 4  # Seq_ID
CS_COMP_ID = 5  # Comp_ID
CS_ATOM_ID = 6  # Atom_ID
CS_ATOM_TYPE = 7  # Atom_type
CS_ISOTOPE = 8  # Atom_isotope_number
CS_VAL = 9  # Val
CS_VAL_ERR = 10  # Val_err
CS_FIG_MERIT = 11  # Assign_fig_of_merit
CS_AMBIG_CODE = 12  # Ambiguity_code
CS_AMBIG_SET = 13  # Ambiguity_set_ID
CS_OCCUPANCY = 14  # Occupancy
CS_RESONANCE = 15  # Resonance_ID
CS_AUTH_ASYM = 16  # Auth_asym_ID
CS_AUTH_SEQ = 17  # Auth_seq_ID
CS_AUTH_COMP = 18  # Auth_comp_ID
CS_AUTH_ATOM = 19  # Auth_atom_ID
CS_ORIG_ASYM = 20  # Original_PDB_strand_ID
CS_ORIG_SEQ = 21  # Original_PDB_residue_no
CS_ORIG_COMP = 22  # Original_PDB_residue_name
CS_ORIG_ATOM = 23  # Original_PDB_atom_name
CS_DETAILS = 24  # Details
CS_ENTRY_ID = 25  # Entry_ID
CS_LIST_ID = 26  # Assigned_chem_shift_list_ID
CS_INS_CODE = 27  # PDB_ins_code, appended only when has_ins_code

NEF_CS_TAGS = ('chain_code',
               'sequence_code',
               'residue_name',
               'atom_name',
               'value',
               'value_uncertainty',
               'element',
               'isotope_number')

# positions within NEF_CS_TAGS
NEF_CHAIN = 0  # chain_code
NEF_SEQ = 1  # sequence_code
NEF_COMP = 2  # residue_name
NEF_ATOM = 3  # atom_name
NEF_VAL = 4  # value
NEF_VAL_ERR = 5  # value_uncertainty
NEF_ELEM = 6  # element
NEF_ISOTOPE = 7  # isotope_number


class NmrDpRemediationCsLoop(NmrDpRemediationBase):
    """ Remediation of the assigned chemical shift loop of NMR data.
    """
    __slots__ = ()

    def remediateCsLoop(self, file_list_id: int, file_type: str, content_subtype: str,
                        sf: Union[pynmrstar.Saveframe, pynmrstar.Loop],
                        list_id: int, sf_framecode: str, lp_category: str) -> bool:
        """ Remediate assigned chemical shift loop based on coordinates.
        """

        def peripheral_loop_rows():
            """ Yield rows of loop at increasing distance from _idx, forward then
                backward, so a scan can inspect the neighbourhood of the current row.
            """
            for offset in range(1, PERIPH_OFFSET_ATTEMPT):
                if _idx + offset < len(loop):
                    yield loop.data[_idx + offset]
                if _idx - offset >= 0:
                    yield loop.data[_idx - offset]

        input_source = self._reg.report.input_sources[file_list_id]
        input_source_dic = input_source.get()

        file_name = input_source_dic['file_name']

        has_poly_seq_in_lp = has_key_value(input_source_dic, 'polymer_sequence_in_loop')

        if not has_poly_seq_in_lp:
            return False

        try:

            if file_type == 'nmr-star':

                _lp_category = '_Systematic_chem_shift_offset'

                _loop = sf.get_loop(_lp_category)

                if 'Type' in _loop.tags:
                    type_col = _loop.tags.index('Type')
                    for _row in _loop:
                        if _row[type_col] in EMPTY_VALUE:
                            continue
                        text = _row[type_col].lower()
                        if 'sail' in text or 'stereo-array isotope labeling' in text:
                            self._reg.sail_flag = True
                            break

                if 'sample' in self._reg.sf_category_list\
                   and '_Sample_component' in self._reg.lp_category_list:

                    _lp_category = '_Sample_component'

                    for _sf in self._reg.star_data[file_list_id].get_saveframes_by_category('sample'):

                        _loop = _sf.get_loop(_lp_category)

                        if 'Isotopic_labeling' in _loop.tags:
                            isotopic_labeling_col = _loop.tags.index('Isotopic_labeling')
                            for _row in _loop:
                                if _row[isotopic_labeling_col] in EMPTY_VALUE:
                                    continue
                                text = _row[isotopic_labeling_col].lower()
                                if 'sail' in text or 'stereo-array isotope labeling' in text:
                                    self._reg.sail_flag = True
                                    break

        except KeyError:
            pass

        poly_seq_in_lp = input_source_dic['polymer_sequence_in_loop']

        poly_seq_common = input_source_dic['polymer_sequence']

        poly_seq = seq_align = chain_assign = br_seq_align = br_chain_assign = np_seq_align = np_chain_assign = None

        if content_subtype in poly_seq_in_lp and self._reg.caC is not None:
            _poly_seq_in_lp = next((_poly_seq_in_lp for _poly_seq_in_lp in poly_seq_in_lp[content_subtype]
                                    if _poly_seq_in_lp['sf_framecode'] == sf_framecode), None)

            if _poly_seq_in_lp is not None:
                list_id = _poly_seq_in_lp['list_id']
                poly_seq = _poly_seq_in_lp['polymer_sequence']

                seq_align, _ =\
                    alignPolymerSequence(self._reg.pA, self._reg.caC['polymer_sequence'], poly_seq, conservative=False)
                chain_assign, _ =\
                    assignPolymerSequence(self._reg.pA, self._reg.ccU, file_type, self._reg.caC['polymer_sequence'],
                                          poly_seq, seq_align)

                if self._reg.caC['branched'] is not None:
                    br_seq_align, _ =\
                        alignPolymerSequence(self._reg.pA, self._reg.caC['branched'], poly_seq, conservative=False)
                    br_chain_assign, _ =\
                        assignPolymerSequence(self._reg.pA, self._reg.ccU, file_type, self._reg.caC['branched'],
                                              poly_seq, br_seq_align)

                if self._reg.caC['non_polymer'] is not None:
                    np_seq_align, _ =\
                        alignPolymerSequence(self._reg.pA, self._reg.caC['non_polymer'], poly_seq, conservative=False)
                    np_chain_assign, _ =\
                        assignPolymerSequence(self._reg.pA, self._reg.ccU, file_type, self._reg.caC['non_polymer'],
                                              poly_seq, np_seq_align)

        @functools.lru_cache()
        def get_auth_seq_scheme(chain_id, seq_id):
            auth_asym_id = auth_seq_id = None

            if seq_id is not None:

                if chain_assign is not None:
                    auth_asym_id = next((ca['ref_chain_id'] for ca in chain_assign if ca['test_chain_id'] == chain_id), None)
                    if auth_asym_id is not None:
                        sa = next((sa for sa in seq_align
                                   if sa['ref_chain_id'] == auth_asym_id and sa['test_chain_id'] == chain_id
                                   and seq_id in sa['test_seq_id'] and sa['sequence_coverage'] >= LOW_SEQ_COVERAGE), None)
                        if sa is not None:
                            _ref_seq_id_name = 'ref_auth_seq_id' if 'ref_auth_seq_id' in sa else 'ref_seq_id'
                            auth_seq_id = next((ref_seq_id for ref_seq_id, test_seq_id
                                                in zip(sa[_ref_seq_id_name], sa['test_seq_id'])
                                                if test_seq_id == seq_id), None)

                if None in (auth_asym_id, auth_seq_id) and br_seq_align is not None:
                    auth_asym_id = next((ca['ref_chain_id'] for ca in br_chain_assign if ca['test_chain_id'] == chain_id), None)
                    if auth_asym_id is not None:
                        sa = next((sa for sa in br_seq_align
                                   if sa['ref_chain_id'] == auth_asym_id and sa['test_chain_id'] == chain_id
                                   and seq_id in sa['test_seq_id'] and sa['sequence_coverage'] >= LOW_SEQ_COVERAGE), None)
                        if sa is not None:
                            _ref_seq_id_name = 'ref_auth_seq_id' if 'ref_auth_seq_id' in sa else 'ref_seq_id'
                            auth_seq_id = next((ref_seq_id for ref_seq_id, test_seq_id
                                                in zip(sa[_ref_seq_id_name], sa['test_seq_id'])
                                                if test_seq_id == seq_id), None)

                if None in (auth_asym_id, auth_seq_id) and np_seq_align is not None:
                    auth_asym_id = next((ca['ref_chain_id'] for ca in np_chain_assign if ca['test_chain_id'] == chain_id), None)
                    if auth_asym_id is not None:
                        sa = next((sa for sa in np_seq_align
                                   if sa['ref_chain_id'] == auth_asym_id and sa['test_chain_id'] == chain_id
                                   and seq_id in sa['test_seq_id'] and sa['sequence_coverage'] >= LOW_SEQ_COVERAGE), None)
                        if sa is not None:
                            _ref_seq_id_name = 'ref_auth_seq_id' if 'ref_auth_seq_id' in sa else 'ref_seq_id'
                            auth_seq_id = next((ref_seq_id for ref_seq_id, test_seq_id
                                                in zip(sa[_ref_seq_id_name], sa['test_seq_id'])
                                                if test_seq_id == seq_id), None)

            return auth_asym_id, auth_seq_id

        @functools.lru_cache()
        def get_label_seq_scheme(chain_id, seq_id):
            auth_asym_id = auth_seq_id = label_seq_id = None

            if seq_id is not None:

                if chain_assign is not None:
                    auth_asym_id = next((ca['ref_chain_id'] for ca in chain_assign if ca['test_chain_id'] == chain_id), None)
                    if auth_asym_id is not None:
                        sa = next((sa for sa in seq_align
                                   if sa['ref_chain_id'] == auth_asym_id and sa['test_chain_id'] == chain_id
                                   and seq_id in sa['ref_auth_seq_id' if 'ref_auth_seq_id' in sa else 'ref_seq_id']
                                   and sa['sequence_coverage'] >= LOW_SEQ_COVERAGE), None)
                        if sa is not None:
                            _ref_seq_id_name = 'ref_auth_seq_id' if 'ref_auth_seq_id' in sa else 'ref_seq_id'
                            auth_seq_id, label_seq_id = next(((ref_seq_id, test_seq_id)
                                                              for ref_seq_id, test_seq_id
                                                              in zip(sa[_ref_seq_id_name], sa['test_seq_id'])
                                                              if ref_seq_id == seq_id), (None, None))

                if None in (auth_asym_id, auth_seq_id) and br_seq_align is not None:
                    auth_asym_id = next((ca['ref_chain_id'] for ca in br_chain_assign if ca['test_chain_id'] == chain_id), None)
                    if auth_asym_id is not None:
                        sa = next((sa for sa in br_seq_align
                                   if sa['ref_chain_id'] == auth_asym_id and sa['test_chain_id'] == chain_id
                                   and seq_id in sa['ref_auth_seq_id' if 'ref_auth_seq_id' in sa else 'ref_seq_id']
                                   and sa['sequence_coverage'] >= LOW_SEQ_COVERAGE), None)
                        if sa is not None:
                            _ref_seq_id_name = 'ref_auth_seq_id' if 'ref_auth_seq_id' in sa else 'ref_seq_id'
                            auth_seq_id, label_seq_id = next(((ref_seq_id, test_seq_id)
                                                              for ref_seq_id, test_seq_id
                                                              in zip(sa[_ref_seq_id_name], sa['test_seq_id'])
                                                              if ref_seq_id == seq_id), (None, None))

                if None in (auth_asym_id, auth_seq_id) and np_seq_align is not None:
                    auth_asym_id = next((ca['ref_chain_id'] for ca in np_chain_assign if ca['test_chain_id'] == chain_id), None)
                    if auth_asym_id is not None:
                        sa = next((sa for sa in np_seq_align
                                   if sa['ref_chain_id'] == auth_asym_id and sa['test_chain_id'] == chain_id
                                   and seq_id in sa['ref_auth_seq_id' if 'ref_auth_seq_id' in sa else 'ref_seq_id']
                                   and sa['sequence_coverage'] >= LOW_SEQ_COVERAGE), None)
                        if sa is not None:
                            _ref_seq_id_name = 'ref_auth_seq_id' if 'ref_auth_seq_id' in sa else 'ref_seq_id'
                            auth_seq_id, label_seq_id = next(((ref_seq_id, test_seq_id)
                                                              for ref_seq_id, test_seq_id
                                                              in zip(sa[_ref_seq_id_name], sa['test_seq_id'])
                                                              if ref_seq_id == seq_id), (None, None))

            return auth_asym_id, auth_seq_id, label_seq_id

        has_ins_code = False

        if poly_seq is not None:

            for ps in poly_seq:

                if has_ins_code:
                    break

                auth_asym_id, _ = get_auth_seq_scheme(ps['chain_id'], ps['seq_id'][0])

                if self._reg.caC['polymer_sequence'] is not None\
                   and any(True for cif_ps in self._reg.caC['polymer_sequence']
                           if cif_ps['auth_chain_id'] == auth_asym_id and 'ins_code' in cif_ps):
                    has_ins_code = True

                if self._reg.caC['branched'] is not None\
                   and any(True for cif_ps in self._reg.caC['branched']
                           if cif_ps['auth_chain_id'] == auth_asym_id and 'ins_code' in cif_ps):
                    has_ins_code = True

                if self._reg.caC['non_polymer'] is not None\
                   and any(True for cif_ps in self._reg.caC['non_polymer']
                           if cif_ps['auth_chain_id'] == auth_asym_id and 'ins_code' in cif_ps):
                    has_ins_code = True

        loop = sf if self._reg.star_data_type[file_list_id] == 'Loop' else sf.get_loop(lp_category)

        # cleanup unnecessary '?'
        item_names = [item['name'] for item in self._reg.key_items[file_type][content_subtype]]
        item_names.extend([item['name'] for item in DATA_ITEMS[file_type][content_subtype]])
        first_row = loop.data[0]
        for item_name in set(loop.tags) - set(item_names):
            item_col = loop.tags.index(item_name)
            if first_row[item_col] == '?':
                for row in loop:
                    if row[item_col] == '?':
                        row[item_col] = None

        aux_lp = None

        aux_lp_category = AUX_LP_CATEGORIES[file_type][content_subtype][0] if file_type == 'nmr-star' else ''

        def delete_aux_loop():

            if file_type == 'nmr-star' and not isinstance(sf, pynmrstar.Loop):

                try:

                    aux_loop = sf.get_loop(aux_lp_category)

                    del sf[aux_loop]

                except KeyError:
                    pass

        if file_type == 'nef':

            items = list(NEF_CS_TAGS)

            mandatory_items = [item['name'] for item in self._reg.key_items[file_type][content_subtype]
                               if 'remove-bad-pattern' in item]

            if not all(tag in loop.tags for tag in mandatory_items):

                err = f"Assigned chemical shifts of {sf_framecode!r} saveframe was not parsed properly. "\
                    "Please fix problems reported."

                self._reg.report.error.appendDescription('missing_mandatory_content',
                                                         {'file_name': file_name, 'description': err})

                if self._reg.verbose:
                    self._reg.log.write(f"+{self.__class_name__}.remediateCsLoop() ++ Error  - {err}\n")

                return False

            mandatory_items = [item['name'] for item in self._reg.key_items[file_type][content_subtype]]
            for item in DATA_ITEMS[file_type][content_subtype]:
                if item['mandatory']:
                    mandatory_items.append(item['name'])

            if not all(tag for tag in mandatory_items if tag in loop.tags):
                return False

            coord_atom_site = self._reg.caC['coord_atom_site'] if self._reg.caC is not None else {}

            chain_id_col = loop.tags.index('chain_code')
            seq_id_col = loop.tags.index('sequence_code')
            comp_id_col = loop.tags.index('residue_name')
            atom_id_col = loop.tags.index('atom_name')
            val_col = loop.tags.index('value')
            val_err_col = loop.tags.index('value_uncertainty') if 'value_uncertainty' in loop.tags else -1

            lp = pynmrstar.Loop.from_scratch(lp_category)

            tags = [f'{lp_category}.{item}' for item in items]

            lp.add_tag(tags)

            for idx, row in enumerate(loop):

                _row = [None] * len(tags)

                try:
                    seq_key = (row[chain_id_col], int(row[seq_id_col]))
                except (ValueError, TypeError):
                    continue

                if seq_key in self._reg.seq_id_map_for_remediation:
                    seq_key = self._reg.seq_id_map_for_remediation[seq_key]

                _row[NEF_CHAIN], _row[NEF_SEQ] = seq_key

                if seq_key in coord_atom_site:
                    _row[NEF_COMP] = coord_atom_site[seq_key]['comp_id']
                else:
                    _row[NEF_COMP] = row[comp_id_col].upper()

                _row[NEF_ATOM] = row[atom_id_col]
                atom_id = row[atom_id_col].upper()

                _row[NEF_VAL] = row[val_col]

                try:
                    float(_row[NEF_VAL])
                except ValueError:
                    continue

                if val_err_col != -1:
                    val_err = row[val_err_col]
                    _row[NEF_VAL_ERR] = val_err

                    if val_err not in EMPTY_VALUE:
                        try:
                            _val_err = float(val_err)
                            if _val_err < 0.0:
                                _row[NEF_VAL_ERR] = abs(_val_err)
                        except ValueError:
                            pass

                _row[NEF_ELEM] = 'H' if _row[NEF_ATOM][0] in PROTON_BEGIN_CODE else atom_id[0]
                if _row[NEF_ELEM] in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS:
                    _row[NEF_ISOTOPE] = ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS[_row[NEF_ELEM]][0]

                lp.add_data(_row)

            key_items = self._reg.key_items[file_type][content_subtype]

            conflict_id = self._reg.nefT.get_conflict_id(lp, lp_category, key_items)[0]

            if len(conflict_id) > 0:
                conflict_id_set = self._reg.nefT.get_conflict_id_set(lp, lp_category, key_items)[0]

                for _id in conflict_id:
                    _id_set = next(id_set for id_set in conflict_id_set if _id in id_set)

                    if len(set(str(lp.data[_id_]) for _id_ in _id_set)) == 1:
                        continue

                    msg = ' vs '.join([str(lp.data[_id_]).replace('None', '.').replace(',', '').replace("'", '')
                                       for _id_ in _id_set])

                    warn = f"Resolved redundancy of assigned chemical shifts ({msg}) by deletion of the latter one."

                    self._reg.report.warning.appendDescription('redundant_data',
                                                               {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                'category': lp_category, 'description': warn})

                    if self._reg.verbose:
                        self._reg.log.write(f"+{self.__class_name__}.remediateCsLoop() ++ Warning  - {warn}\n")

                for _id in conflict_id:
                    del lp.data[_id]

        else:

            if 'original_file_name' in input_source_dic:
                tagNames = [t[0] for t in sf.tags]
                if 'Data_file_name' not in tagNames:
                    sf.add_tag('Data_file_name', input_source_dic['original_file_name'])

            items = list(CS_TAGS)

            if has_ins_code:
                items.append('PDB_ins_code')

            mandatory_items = [item['name'] for item in self._reg.key_items[file_type][content_subtype]
                               if 'remove-bad-pattern' in item]

            if not all(tag in loop.tags for tag in mandatory_items):

                err = f"Assigned chemical shifts of {sf_framecode!r} saveframe was not parsed properly. "\
                    "Please fix problems reported."

                self._reg.report.error.appendDescription('missing_mandatory_content',
                                                         {'file_name': file_name, 'description': err})

                if self._reg.verbose:
                    self._reg.log.write(f"+{self.__class_name__}.remediateCsLoop() ++ Error  - {err}\n")

                return False

            mandatory_items = [item['name'] for item in self._reg.key_items[file_type][content_subtype]]
            for item in DATA_ITEMS[file_type][content_subtype]:
                if item['mandatory']:
                    mandatory_items.append(item['name'])

            if not all(tag for tag in mandatory_items if tag in loop.tags):
                return False

            auth_pdb_tags = ['Auth_asym_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID']
            orig_pdb_tags = ['Original_PDB_strand_ID', 'Original_PDB_residue_no',
                             'Original_PDB_residue_name', 'Original_PDB_atom_name']

            entity_assembly = self._reg.caC['entity_assembly'] if self._reg.caC is not None else []
            auth_to_entity_type = self._reg.caC['auth_to_entity_type'] if self._reg.caC is not None else {}
            auth_to_star_seq = self._reg.caC['auth_to_star_seq'] if self._reg.caC is not None else {}
            auth_to_orig_seq = self._reg.caC['auth_to_orig_seq'] if self._reg.caC is not None else {}
            auth_to_ins_code = self._reg.caC['auth_to_ins_code'] if self._reg.caC is not None else {}
            auth_to_star_seq_ann = self._reg.caC['auth_to_star_seq_ann'] if self._reg.caC is not None else {}
            coord_atom_site = self._reg.caC['coord_atom_site'] if self._reg.caC is not None else {}
            coord_unobs_res = self._reg.caC['coord_unobs_res'] if self._reg.caC is not None else {}
            auth_atom_name_to_id = self._reg.caC['auth_atom_name_to_id'] if self._reg.caC is not None else {}
            auth_atom_name_to_id_ext = self._reg.caC['auth_atom_name_to_id_ext'] if self._reg.caC is not None else {}
            mis_poly_link = self._reg.caC['missing_polymer_linkage'] if self._reg.caC is not None else []

            _auth_to_orig_seq = {}

            has_auth_seq = valid_auth_seq = has_auth_chain = False
            aux_auth_seq_id_col = aux_auth_comp_id_col = aux_auth_atom_id_col = -1
            auth_asym_id_col = loop.tags.index('Auth_asym_ID') if 'Auth_asym_ID' in loop.tags else -1
            auth_seq_id_col = loop.tags.index('Auth_seq_ID') if 'Auth_seq_ID' in loop.tags else -1
            auth_comp_id_col = loop.tags.index('Auth_comp_ID') if 'Auth_comp_ID' in loop.tags else -1
            auth_atom_id_col = loop.tags.index('Auth_atom_ID') if 'Auth_atom_ID' in loop.tags else -1

            # split concatenation of auth_seq_id and ins_code (DAOTHER-10418)
            if auth_to_ins_code is not None and len(auth_to_ins_code) > 0 and auth_seq_id_col != -1:
                auth_dat = loop.get_tag(['Auth_seq_ID'])

                if any(True for row in auth_dat if isinstance(row, str)):

                    ins_code_col = loop.tags.index('PDB_ins_code') if 'PDB_ins_code' in loop.tags else -1
                    if ins_code_col == -1:
                        loop.add_tag('PDB_ins_code', update_data=True)

                    for idx, row in enumerate(auth_dat):
                        if isinstance(row, str) and CONCAT_SEQ_ID_INS_CODE_PAT.match(row):
                            g = CONCAT_SEQ_ID_INS_CODE_PAT.search(row).groups()
                            loop.data[idx][auth_seq_id_col] = g[0]
                            if g[1] not in EMPTY_VALUE:
                                loop.data[idx][ins_code_col] = g[1]

                    if not has_ins_code:
                        items.append('PDB_ins_code')
                        has_ins_code = True

            valid_auth_seq_per_chain = []

            if set(auth_pdb_tags) & set(loop.tags) == set(auth_pdb_tags):
                auth_dat = loop.get_tag(auth_pdb_tags)
                if len(auth_dat) > 0:
                    has_auth_seq = valid_auth_seq = True
                    if not self._reg.annotation_mode or len(coord_unobs_res) > 0:
                        for row in auth_dat:
                            try:
                                seq_key = (row[0], int(row[1]), row[2])
                                if seq_key not in auth_to_star_seq_ann:
                                    valid_auth_seq = False
                                    break
                            except (ValueError, TypeError):
                                has_auth_seq = valid_auth_seq = False
                                break

            # DAOTHER-9281
            elif auth_asym_id_col != -1:
                has_auth_chain = True
                _auth_pdb_tags = ['Auth_asym_ID']
                if 'Auth_seq_ID' in loop.tags:
                    _auth_pdb_tags.append('Auth_seq_ID')
                elif 'Comp_index_ID' in loop.tags:
                    _auth_pdb_tags.append('Comp_index_ID')
                elif 'Seq_ID' in loop.tags:
                    _auth_pdb_tags.append('Seq_ID')

                if 'Auth_comp_ID' in loop.tags:
                    _auth_pdb_tags.append('Auth_comp_ID')
                elif 'Comp_ID' in loop.tags:
                    _auth_pdb_tags.append('Comp_ID')

                if 'Auth_atom_ID' in loop.tags:
                    _auth_pdb_tags.append('Auth_atom_ID')
                elif 'Atom_ID' in loop.tags:
                    _auth_pdb_tags.append('Atom_ID')

                if len(_auth_pdb_tags) == 4:
                    auth_dat = loop.get_tag(_auth_pdb_tags)
                    if len(auth_dat) > 0:
                        aux_auth_seq_id_col = loop.tags.index(_auth_pdb_tags[1])
                        aux_auth_comp_id_col = loop.tags.index(_auth_pdb_tags[2])
                        aux_auth_atom_id_col = loop.tags.index(_auth_pdb_tags[3])

                        valid_auth_seq = True
                        if not self._reg.annotation_mode or len(coord_unobs_res) > 0:
                            for row in auth_dat:
                                try:
                                    seq_key = (row[0], int(row[1]), row[2])
                                    if seq_key not in auth_to_star_seq_ann:
                                        valid_auth_seq = False
                                        break
                                except (ValueError, TypeError):
                                    valid_auth_seq = False
                                    break
                        if not valid_auth_seq:
                            for row in auth_dat:
                                if row[0] not in valid_auth_seq_per_chain:
                                    valid_auth_seq_per_chain.append(row[0])
                            if not self._reg.annotation_mode or len(coord_unobs_res) > 0:
                                for row in auth_dat:
                                    try:
                                        seq_key = (row[0], int(row[1]), row[2])
                                        if seq_key not in auth_to_star_seq_ann:
                                            if row[0] in valid_auth_seq_per_chain:
                                                valid_auth_seq_per_chain.remove(row[0])
                                    except (ValueError, TypeError):
                                        if row[0] in valid_auth_seq_per_chain:
                                            valid_auth_seq_per_chain.remove(row[0])

            has_orig_seq = ch2_name_in_xplor = ch3_name_in_xplor = False

            if self._reg.remediation_mode:
                if set(orig_pdb_tags) & set(loop.tags) == set(orig_pdb_tags):
                    orig_dat = loop.get_tag(orig_pdb_tags)
                    if len(orig_dat) > 0:
                        for row in orig_dat:
                            if all(d not in EMPTY_VALUE for d in row):
                                has_orig_seq = True
                                break
                        if has_orig_seq:
                            orig_pdb_tags.append('Comp_ID')
                            orig_pdb_tags.append('Atom_ID')
                            dat = loop.get_tag(orig_pdb_tags)
                            for row in dat:
                                if row[3] in EMPTY_VALUE:
                                    continue
                                orig_atom_id = row[3].upper()
                                comp_id = row[4]
                                atom_id = row[5]
                                if orig_atom_id == atom_id:
                                    continue
                                ambig_code = self._reg.csStat.getMaxAmbigCodeWoSetId(comp_id, atom_id)
                                if ambig_code == 0 or atom_id[0] not in PROTON_BEGIN_CODE:
                                    continue
                                len_in_grp = len(self._reg.csStat.getProtonsInSameGroup(comp_id, atom_id))
                                if len_in_grp == 2 and ambig_code == 2:
                                    ch2_name_in_xplor = any(True for r, o in zip(atom_id, orig_atom_id) if r == '3' and o == '1')
                                elif len_in_grp == 3 and atom_id[-1] == orig_atom_id[0]:
                                    ch3_name_in_xplor = True
            else:
                if set(orig_pdb_tags) & set(loop.tags) == set(orig_pdb_tags):
                    orig_dat = loop.get_tag(orig_pdb_tags)
                    if len(orig_dat) > 0:
                        for row in orig_dat:
                            if all(d not in EMPTY_VALUE for d in row):
                                has_orig_seq = True
                                break

            entity_assembly_mappping = {}
            if self._reg.bmrb_only and self._reg.internal_mode:
                if isinstance(self._reg.star_data[file_list_id], pynmrstar.Entry):
                    for asm_sf in self._reg.star_data[file_list_id].get_saveframes_by_category('assembly'):
                        try:
                            ea_loop = asm_sf.get_loop('_Entity_assembly')
                            dat = ea_loop.get_tag(['ID', 'Entity_ID'])
                            for row in dat:
                                entity_assembly_id = row[0] if isinstance(row[0], str) else str(row[0])
                                entity_assembly_mappping[entity_assembly_id] = row[1]
                        except KeyError:
                            continue

            chain_id_col = loop.tags.index('Entity_assembly_ID')
            entity_id_col = loop.tags.index('Entity_ID') if 'Entity_ID' in loop.tags else -1
            seq_id_col = loop.tags.index('Comp_index_ID')
            comp_id_col = loop.tags.index('Comp_ID')
            atom_id_col = loop.tags.index('Atom_ID')
            val_col = loop.tags.index('Val') if 'Val' in loop.tags else loop.tags.index('Chem_shift_val')
            val_err_col = loop.tags.index('Val_err') if 'Val_err' in loop.tags\
                else loop.tags.index('Chem_shift_val_err') if 'Chem_shift_val_err' in loop.tags else -1
            fig_of_merit_col = loop.tags.index('Assign_fig_of_merit') if 'Assign_fig_of_merit' in loop.tags\
                else loop.tags.index('Chem_shift_assign_fig_of_merit') if 'Chem_shift_assign_fig_of_merit' in loop.tags else -1
            ambig_code_col = loop.tags.index('Ambiguity_code') if 'Ambiguity_code' in loop.tags\
                else loop.tags.index('Chem_shift_ambiguity_code') if 'Chem_shift_ambiguity_code' in loop.tags else -1
            ambig_set_id_col = loop.tags.index('Ambiguity_set_ID') if 'Ambiguity_set_ID' in loop.tags else -1
            occupancy_col = loop.tags.index('Occupancy') if 'Occupancy' in loop.tags else -1
            reson_id_col = loop.tags.index('Resonance_ID') if 'Resonance_ID' in loop.tags else -1
            details_col = loop.tags.index('Details') if 'Details' in loop.tags else -1

            if self._reg.annotation_mode and details_col != -1:
                for row in loop:
                    if row[details_col] == 'UNMAPPED':
                        row[details_col] = None

            trial = 0

            incomplete_comp_id_annotation = []  # DAOTHER-9286
            truncated_loop_sequence = []  # DAOTHER-9644
            for mis in mis_poly_link:
                auth_chain_id = mis['auth_chain_id']
                auth_seq_id_1 = mis['auth_seq_id_1']
                auth_seq_id_2 = mis['auth_seq_id_2']

                cif_ps = next((cif_ps for cif_ps in self._reg.caC['polymer_sequence']
                               if cif_ps['auth_chain_id'] == auth_chain_id), None)

                if cif_ps is not None and auth_seq_id_1 in cif_ps['auth_seq_id'] and auth_seq_id_2 in cif_ps['auth_seq_id']\
                   and auth_seq_id_1 < auth_seq_id_2:
                    for auth_seq_id in range(auth_seq_id_1 + 1, auth_seq_id_2):
                        _seq_key = (auth_chain_id, auth_seq_id)
                        truncated_loop_sequence.append(_seq_key)

            # DAOTHER-10898
            copied_auth_asym_id_mapping = {}

            def apply_orig_seq_naming(*, std_res_branch=False, aux_cols=False,
                                      guard_has_auth_seq=False, record_can_mapping=False):
                """ Carry the original (author) sequence and atom naming of the current row
                    into _row, then fill the row in.

                    Collapses five sites that ran this same sequence; the keyword flags are
                    the only things that ever differed between them.
                    @param std_res_branch: also accept a comp_id that translates to auth_comp_id
                    @param aux_cols: read the author naming from the auxiliary loop columns
                    @param guard_has_auth_seq: copy the author naming only if the loop has it
                    @param record_can_mapping: remember the canonical author chain mapping
                """
                nonlocal _row, _seq_key, ambig_code, comp_id, index, len_in_grp, reparse, reparse_request, seq_key

                if has_ins_code and seq_key in auth_to_ins_code:
                    _row[CS_INS_CODE] = auth_to_ins_code[seq_key]

                if seq_key in auth_to_orig_seq:
                    if _row[CS_ORIG_ASYM] not in EMPTY_VALUE and seq_key not in _auth_to_orig_seq:
                        orig_seq_id, orig_comp_id = auth_to_orig_seq[seq_key]
                        __seq_key = (_seq_key[0], orig_seq_id, comp_id)
                        if self._reg.csStat.getTypeOfCompId(comp_id)[2]\
                           and seq_key not in coord_atom_site and __seq_key in auth_to_star_seq:
                            _seq_key = __seq_key
                            if _row[CS_ORIG_SEQ] in EMPTY_VALUE or _row[CS_ORIG_COMP] in EMPTY_VALUE:
                                _row[CS_ORIG_SEQ], _row[CS_ORIG_COMP] = orig_seq_id, orig_comp_id
                        else:
                            _auth_to_orig_seq[seq_key] = (_row[CS_ORIG_ASYM], orig_seq_id, orig_comp_id)
                    if not has_orig_seq:
                        orig_seq_id, orig_comp_id = auth_to_orig_seq[seq_key]
                        if orig_seq_id in EMPTY_VALUE:
                            orig_seq_id = auth_seq_id
                        if orig_comp_id in EMPTY_VALUE:
                            orig_comp_id = comp_id
                        _row[CS_ORIG_ASYM], _row[CS_ORIG_SEQ], _row[CS_ORIG_COMP], _row[CS_ORIG_ATOM] =\
                            auth_asym_id, orig_seq_id, orig_comp_id, _orig_atom_id
                    elif any(True for d in orig_dat[idx] if d in EMPTY_VALUE):
                        if seq_key in _auth_to_orig_seq:
                            _row[CS_ORIG_ASYM], _row[CS_ORIG_SEQ], _row[CS_ORIG_COMP] = _auth_to_orig_seq[seq_key]
                        elif std_res_branch and comp_id != auth_comp_id\
                                and translateToStdResName(comp_id, ccU=self._reg.ccU) == auth_comp_id:
                            _row[CS_ORIG_ASYM], _row[CS_ORIG_SEQ], _row[CS_ORIG_COMP] = auth_asym_id, auth_seq_id, comp_id
                            _row[CS_COMP_ID] = comp_id = auth_comp_id
                        if _row[CS_ORIG_ATOM] in EMPTY_VALUE:
                            _row[CS_ORIG_ATOM] = atom_id
                        ambig_code = self._reg.csStat.getMaxAmbigCodeWoSetId(comp_id, atom_id)
                        if ambig_code > 0:
                            orig_seq_id, orig_comp_id = auth_to_orig_seq[seq_key]
                            if orig_seq_id in EMPTY_VALUE:
                                orig_seq_id = auth_seq_id
                            if orig_comp_id in EMPTY_VALUE:
                                orig_comp_id = comp_id
                            _row[CS_ORIG_ASYM], _row[CS_ORIG_SEQ], _row[CS_ORIG_COMP] =\
                                auth_asym_id, orig_seq_id, orig_comp_id
                            if atom_id[0] not in PROTON_BEGIN_CODE:
                                _row[CS_ORIG_ATOM] = atom_id
                            else:
                                len_in_grp = len(self._reg.csStat.getProtonsInSameGroup(comp_id, atom_id))
                                if len_in_grp == 2:
                                    _row[CS_ORIG_ATOM] = f'{atom_id[0:-1]}1'\
                                        if ambig_code == 2 and ch2_name_in_xplor and atom_id[-1] == '3' else atom_id
                                elif len_in_grp == 3:
                                    _row[CS_ORIG_ATOM] = (atom_id[-1] + atom_id[0:-1])\
                                        if ch3_name_in_xplor and atom_id[0] == 'H'\
                                        and atom_id[-1] in ('1', '2', '3')\
                                        else atom_id
                                elif _row[CS_ORIG_ATOM] in EMPTY_VALUE:
                                    _row[CS_ORIG_ATOM] = atom_id

                else:
                    seq_key = next((k for k, v in auth_to_star_seq.items()
                                    if v[0] == entity_assembly_id and v[1] == seq_id and v[2] == entity_id), None)
                    if seq_key is not None:
                        _seq_key = (seq_key[0], seq_key[1])
                        _row[CS_AUTH_ASYM], _row[CS_AUTH_SEQ], _row[CS_AUTH_COMP], _row[CS_AUTH_ATOM] =\
                            seq_key[0], seq_key[1], seq_key[2], atom_id

                        if has_ins_code and seq_key in auth_to_ins_code:
                            _row[CS_INS_CODE] = auth_to_ins_code[seq_key]

                    if not guard_has_auth_seq or has_auth_seq:
                        _row[CS_ORIG_ASYM], _row[CS_ORIG_SEQ], _row[CS_ORIG_COMP], _row[CS_ORIG_ATOM] =\
                            (row[auth_asym_id_col], row[aux_auth_seq_id_col],
                             row[aux_auth_comp_id_col], row[aux_auth_atom_id_col]) if aux_cols\
                            else (row[auth_asym_id_col], row[auth_seq_id_col],
                                  row[auth_comp_id_col], row[auth_atom_id_col])

                index, _row, reparse = fill_cs_row(lp, index, _row, prefer_auth_atom_name, coord_atom_site, _seq_key,
                                                   comp_id, atom_id, loop, idx)
                reparse_request |= reparse

                if record_can_mapping and chain_id not in can_auth_asym_id_mapping:
                    can_auth_asym_id_mapping[chain_id] = {'auth_asym_id': auth_asym_id,
                                                          'ref_auth_seq_id': auth_seq_id
                                                          }

            def fill_cs_row(lp, index, _row, prefer_auth_atom_name, coord_atom_site, _seq_key, comp_id, atom_id, src_lp, src_idx):
                original_index = index

                def peripheral_src_rows():
                    """ Yield rows of src_lp at increasing distance from src_idx, forward then
                        backward, so a scan can inspect the neighbourhood of the current row.
                    """
                    for offset in range(1, PERIPH_OFFSET_ATTEMPT):
                        if src_idx + offset < len(src_lp):
                            yield src_lp.data[src_idx + offset]
                        if src_idx - offset >= 0:
                            yield src_lp.data[src_idx - offset]

                reparse = False
                _src_idx = src_idx
                if src_idx > 0:
                    src_idx -= 1
                fill_auth_atom_id = self._reg.annotation_mode\
                    or (_row[CS_AUTH_ATOM] in EMPTY_VALUE and _row[CS_AUTH_COMP] not in EMPTY_VALUE)
                fill_orig_atom_id = _row[CS_ORIG_ATOM] not in EMPTY_VALUE

                if _seq_key is not None:
                    if _seq_key in truncated_loop_sequence and _row[CS_DETAILS] in EMPTY_VALUE:
                        _row[CS_DETAILS] = 'UNMAPPED'
                    seq_key = (_seq_key[0], _seq_key[1], comp_id)
                    _seq_key = seq_key if seq_key in coord_atom_site else _seq_key
                if _seq_key in coord_atom_site\
                   and (coord_atom_site[_seq_key]['comp_id'] == comp_id
                        or (_seq_key not in coord_unobs_res and coord_atom_site[_seq_key]['comp_id'] not in STD_MON_DICT)
                        or (_seq_key in coord_unobs_res and coord_unobs_res[_seq_key]['comp_id'] != comp_id)):
                    # 8b9r: A:24:VAL (unobserved), A:24:CU
                    _coord_atom_site = coord_atom_site[_seq_key]
                    _atom_site_atom_id = _coord_atom_site['atom_id']
                    # DAOTHER-8817
                    if 'chain_id' in _coord_atom_site:
                        _row[CS_AUTH_ASYM] = _coord_atom_site['chain_id']
                        # DAOTHER-10898
                        if _row[CS_AUTH_ASYM] in copied_auth_asym_id_mapping:
                            _row[CS_AUTH_ASYM] = _row[CS_ORIG_ASYM] = copied_auth_asym_id_mapping[_row[CS_AUTH_ASYM]]
                    _row[CS_COMP_ID] = _row[CS_AUTH_COMP] = comp_id = _coord_atom_site['comp_id']
                    valid = True
                    missing_ch3 = []
                    if not self._reg.annotation_mode and atom_id in self._reg.csStat.getMethylProtons(comp_id):
                        missing_ch3 = self._reg.csStat.getProtonsInSameGroup(comp_id, atom_id, True)
                        valid = self._reg.sail_flag
                        row_src = src_lp.data[src_idx]
                        seq_id_src = row_src[seq_id_col]
                        if 0 <= src_idx < len(src_lp):
                            for row in peripheral_src_rows():
                                if (row[seq_id_col] == str(_row[CS_COMP_IDX])
                                    or (_row[CS_COMP_IDX] != seq_id_src and row[seq_id_col] == seq_id_src))\
                                   and row[comp_id_col].upper() == comp_id\
                                   and row[atom_id_col] in missing_ch3:
                                    valid = True
                                    missing_ch3.remove(row[atom_id_col])
                                    if len(missing_ch3) == 0:
                                        break
                    if atom_id in _atom_site_atom_id and valid and len(missing_ch3) == 0\
                       and (not self._reg.annotation_mode or comp_id not in incomplete_comp_id_annotation):
                        _row[CS_ATOM_ID] = atom_id
                        if fill_auth_atom_id or _row[CS_ATOM_ID] != _row[CS_AUTH_ATOM]:
                            _row[CS_AUTH_ATOM] = _row[CS_ATOM_ID]
                        _row[CS_ATOM_TYPE] = _coord_atom_site['type_symbol'][_atom_site_atom_id.index(atom_id)]
                        if _row[CS_ATOM_TYPE] in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS:
                            _row[CS_ISOTOPE] = ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS[_row[CS_ATOM_TYPE]][0]
                        # """ need to preserve Original_PDB_atom_name for atom name mapping history
                        # if fill_orig_atom_id and _row[CS_ATOM_ID] != _row[CS_ORIG_ATOM] and _row[CS_ORIG_ATOM] in _atom_site_atom_id:  # noqa: E501, pylint: disable=line-too-long
                        #     if _row[CS_ORIG_ATOM] in self._reg.csStat.getProtonsInSameGroup(comp_id, atom_id, True):
                        #         _row[CS_ORIG_ATOM] = copy.copy(atom_id)
                        # """
                    else:
                        if atom_id in ('H1', 'HT1') and 'H' in _atom_site_atom_id\
                           and atom_id not in _atom_site_atom_id:
                            if self._reg.ccU.updateChemCompDict(comp_id):
                                cca = next((cca for cca in self._reg.ccU.lastAtomDictList
                                            if cca['atom_id'] == atom_id
                                            and cca['leaving_atom_flag'] == 'N'), None)
                                if cca is None:
                                    atom_id = 'H'
                                    if fill_auth_atom_id:
                                        _row[CS_AUTH_ATOM] = atom_id
                            else:
                                atom_id = 'H'
                                if fill_auth_atom_id:
                                    _row[CS_AUTH_ATOM] = atom_id
                        elif atom_id in ('H', 'HT1') and 'H1' in _atom_site_atom_id\
                                and atom_id not in _atom_site_atom_id:
                            if self._reg.ccU.updateChemCompDict(comp_id):
                                cca = next((cca for cca in self._reg.ccU.lastAtomDictList
                                            if cca['atom_id'] == atom_id
                                            and cca['leaving_atom_flag'] == 'N'), None)
                                if cca is None:
                                    atom_id = 'H1'
                                    if fill_auth_atom_id:
                                        _row[CS_AUTH_ATOM] = atom_id
                            else:
                                atom_id = 'H1'
                                if fill_auth_atom_id:
                                    _row[CS_AUTH_ATOM] = atom_id
                        elif atom_id in AMINO_PROTON_CODE and f'C{atom_id[1:]}' in _atom_site_atom_id:
                            bonded = self._reg.ccU.getBondedAtoms(comp_id, f'C{atom_id[1:]}', onlyProton=True)
                            if len(bonded) == 1 and bonded[0] in _atom_site_atom_id:
                                atom_id = bonded[0]
                                if fill_auth_atom_id:
                                    _row[CS_AUTH_ATOM] = atom_id
                        if len(missing_ch3) > 0 and (_row[CS_VAL] in EMPTY_VALUE or float(_row[CS_VAL]) >= 4.0):
                            heme = False
                            if _row[CS_VAL] not in EMPTY_VALUE:
                                if self._reg.ccU.updateChemCompDict(comp_id):
                                    heme = comp_id == 'HEM' or 'HEME' in self._reg.ccU.lastChemCompDict['name']
                            if not heme:
                                missing_ch3 = []
                        _atom_id = atom_id
                        if not valid and len(missing_ch3) > 0 and atom_id not in _atom_site_atom_id:
                            atom_id = atom_id[:-1]
                            if _atom_id in self._reg.csStat.getRepMethylProtons(comp_id):
                                atom_id = _atom_id
                        if (valid and atom_id in _atom_site_atom_id)\
                           or ((prefer_auth_atom_name or _row[CS_DETAILS] == 'UNMAPPED') and atom_id[0] not in ('Q', 'M')):
                            atom_ids = [atom_id]
                            # DAOTHER-9286
                            if self._reg.annotation_mode and comp_id in incomplete_comp_id_annotation and trial > 0:
                                atom_ids =\
                                    self._reg.dpV.getAtomIdListInXplorForLigandRemap(comp_id,
                                                                                     _row[CS_ORIG_ATOM] if fill_orig_atom_id
                                                                                     else atom_id,
                                                                                     _coord_atom_site)
                        else:
                            atom_ids = self._reg.dpV.getAtomIdListInXplor(comp_id, atom_id)
                            if len(atom_ids) == 0 or atom_ids[0] not in _atom_site_atom_id:
                                atom_ids =\
                                    self._reg.dpV.getAtomIdListInXplor(comp_id,
                                                                       translateToStdAtomName(atom_id, comp_id, _atom_site_atom_id,
                                                                                              ccU=self._reg.ccU))
                                if len(atom_ids) == 1 and atom_ids[0] in _atom_site_atom_id and atom_id not in _atom_site_atom_id:
                                    atom_id = atom_ids[0]
                            # DAOTHER-9286
                            if self._reg.annotation_mode and (len(atom_ids) == 0 or atom_ids[0] not in _atom_site_atom_id):
                                atom_ids = self._reg.dpV.getAtomIdListInXplorForLigandRemap(comp_id, atom_id, _coord_atom_site)
                                if comp_id not in incomplete_comp_id_annotation:
                                    incomplete_comp_id_annotation.append(comp_id)
                        if valid and len(missing_ch3) > 0:
                            if not fill_orig_atom_id or not any(c in ('x', 'y', 'X', 'Y') for c in _row[CS_ORIG_ATOM])\
                               and len(self._reg.dpV.getAtomIdListInXplor(comp_id, _row[CS_ORIG_ATOM])) > 1\
                               and _row[CS_DETAILS] != 'UNMAPPED':
                                atom_ids = self._reg.dpV.getAtomIdListInXplor(comp_id, _row[CS_ORIG_ATOM])
                            else:
                                missing_ch3.clear()
                        if not valid and len(missing_ch3) > 0 and atom_id in _atom_site_atom_id:
                            atom_ids.extend(missing_ch3)
                        len_atom_ids = len(atom_ids)
                        if len_atom_ids == 0:
                            _row[CS_ATOM_ID] = atom_id
                            if fill_auth_atom_id:
                                _row[CS_AUTH_ATOM] = _row[CS_ATOM_ID]
                            _row[CS_ATOM_TYPE] = 'H' if atom_id[0] in PSE_PRO_BEGIN_CODE else atom_id[0]
                            if _row[CS_ATOM_TYPE] in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS:
                                _row[CS_ISOTOPE] = ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS[_row[CS_ATOM_TYPE]][0]
                        else:
                            methyl_atoms = self._reg.csStat.getMethylAtoms(comp_id)
                            atom_ids = sorted(atom_ids)
                            _row[CS_ATOM_ID] = atom_ids[0]
                            _row[CS_AUTH_ATOM] = None
                            fill_auth_atom_id = _row[CS_AUTH_COMP] not in EMPTY_VALUE
                            if self._reg.ccU.updateChemCompDict(comp_id):
                                cca = next((cca for cca in self._reg.ccU.lastAtomDictList
                                            if cca['atom_id'] == _row[CS_ATOM_ID]), None)
                                if cca is not None:
                                    _row[CS_ATOM_TYPE] = cca['type_symbol']
                                    if _row[CS_ATOM_TYPE] in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS:
                                        _row[CS_ISOTOPE] = ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS[_row[CS_ATOM_TYPE]][0]
                                else:
                                    _row[CS_ATOM_TYPE] = 'H' if _row[CS_ATOM_ID][0] in PROTON_BEGIN_CODE else atom_id[0]
                                    if _row[CS_ATOM_TYPE] in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS:
                                        _row[CS_ISOTOPE] = ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS[_row[CS_ATOM_TYPE]][0]
                            else:
                                _row[CS_ATOM_TYPE] = 'H' if atom_id[0] in PSE_PRO_BEGIN_CODE else atom_id[0]
                                if _row[CS_ATOM_TYPE] in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS:
                                    _row[CS_ISOTOPE] = ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS[_row[CS_ATOM_TYPE]][0]

                            ambig_code = _row[CS_AMBIG_CODE]
                            if ambig_code == 0:
                                _row[CS_AMBIG_CODE] = None

                            elif ambig_code in (2, 3):
                                _ambig_code = self._reg.csStat.getMaxAmbigCodeWoSetId(comp_id, _row[CS_ATOM_ID])
                                if _ambig_code not in (0, ambig_code):
                                    if _ambig_code != 1:
                                        _row[CS_AMBIG_CODE] = _ambig_code
                                    else:
                                        _row[CS_AMBIG_CODE] = ambig_code = 4
                                        if 0 <= _src_idx < len(src_lp):
                                            row_src = src_lp.data[_src_idx]
                                            chain_id_src = row_src[chain_id_col]
                                            seq_id_src = row_src[seq_id_col]
                                            atom_type = row_src[atom_id_col][0]
                                            val = float(row_src[val_col])
                                            sig = self._reg.ccU.getBondSignature(comp_id, atom_id)
                                            for row in peripheral_src_rows():
                                                if row[chain_id_col] == chain_id_src\
                                                   and row[seq_id_col] == seq_id_src\
                                                   and row[comp_id_col] == comp_id\
                                                   and row[atom_id_col][0] == atom_type\
                                                   and abs(float(row[val_col]) - val) < 1.0\
                                                   and self._reg.ccU.getBondSignature(comp_id, row[atom_id_col]) == sig:
                                                    row[ambig_code_col] = '4'
                                                    reparse = True

                            elif ambig_code == 4:
                                if not self._reg.annotation_mode and _row[CS_DETAILS] != 'UNMAPPED':
                                    row_src = src_lp.data[_src_idx]
                                    chain_id_src = row_src[chain_id_col]
                                    atom_id_src = row_src[atom_id_col]
                                    atom_type = atom_id_src[0]
                                    ambig_code_src = row_src[ambig_code_col]
                                    atom_ids_in_group_src = self._reg.ccU.getProtonsInSameGroup(comp_id, atom_id_src)\
                                        if atom_type in PROTON_BEGIN_CODE else []
                                    ambig_code_4_test = hetero_group_test = False
                                    for row in peripheral_src_rows():
                                        if row[comp_id_col] == comp_id\
                                           and row[atom_id_col][0] == atom_type\
                                           and row[ambig_code_col] == str(_row[CS_AMBIG_CODE])\
                                           or (_row[CS_AMBIG_CODE] != ambig_code_src
                                               and row[ambig_code_col] == ambig_code_src):
                                            if not (row[chain_id_col] == str(_row[CS_ENT_ASM_ID])
                                                    or (_row[CS_ENT_ASM_ID] != chain_id_src
                                                        and row[chain_id_col] == chain_id_src)):
                                                break
                                            _seq_id =\
                                                row[seq_id_col] if isinstance(row[seq_id_col], int) else int(row[seq_id_col])
                                            if _seq_id in (_row[CS_COMP_IDX], _row[CS_AUTH_SEQ]):
                                                ambig_code_4_test = True
                                                if row[atom_id_col] not in atom_ids_in_group_src:
                                                    hetero_group_test = True
                                                    break
                                    if not ambig_code_4_test:
                                        ambig_code_5_test = False
                                        for row in peripheral_src_rows():
                                            if row[comp_id_col] == comp_id\
                                               and row[atom_id_col][0] == atom_type\
                                               and row[ambig_code_col] == str(_row[CS_AMBIG_CODE])\
                                               or (_row[CS_AMBIG_CODE] != ambig_code_src
                                                   and row[ambig_code_col] == ambig_code_src):
                                                if not (row[chain_id_col] == str(_row[CS_ENT_ASM_ID])
                                                        or (_row[CS_ENT_ASM_ID] != chain_id_src
                                                            and row[chain_id_col] == chain_id_src)):
                                                    break
                                                _seq_id = row[seq_id_col] if isinstance(row[seq_id_col], int)\
                                                    else int(row[seq_id_col])
                                                if _seq_id in (_row[CS_COMP_IDX], _row[CS_AUTH_SEQ]):
                                                    break
                                                _row[CS_AMBIG_CODE] = ambig_code = 5
                                                ambig_code_5_test = True
                                                break
                                        if not ambig_code_5_test:
                                            for row in peripheral_src_rows():
                                                if row[comp_id_col] == comp_id\
                                                   and row[atom_id_col][0] == atom_type\
                                                   and row[ambig_code_col] == str(_row[CS_AMBIG_CODE])\
                                                   or (_row[CS_AMBIG_CODE] != ambig_code_src
                                                       and row[ambig_code_col] == ambig_code_src):
                                                    if not (row[chain_id_col] == str(_row[CS_ENT_ASM_ID])
                                                            or (_row[CS_ENT_ASM_ID] != chain_id_src
                                                                and row[chain_id_col] == chain_id_src)):
                                                        _row[CS_AMBIG_CODE] = ambig_code = 6
                                                        break
                                            if ambig_code == 4:
                                                _row[CS_AMBIG_CODE] = ambig_code = 1
                                    elif not hetero_group_test:
                                        _row[CS_AMBIG_CODE] = ambig_code = 1

                            elif ambig_code == 5:
                                if not self._reg.annotation_mode and _row[CS_DETAILS] != 'UNMAPPED':
                                    row_src = src_lp.data[_src_idx]
                                    chain_id_src = row_src[chain_id_col]
                                    atom_type = row_src[atom_id_col][0]
                                    ambig_code_src = row_src[ambig_code_col]
                                    ambig_code_5_test = False
                                    for row in peripheral_src_rows():
                                        if row[comp_id_col] == comp_id\
                                           and row[atom_id_col][0] == atom_type\
                                           and row[ambig_code_col] == str(_row[CS_AMBIG_CODE])\
                                           or (_row[CS_AMBIG_CODE] != ambig_code_src
                                               and row[ambig_code_col] == ambig_code_src):
                                            if not (row[chain_id_col] == str(_row[CS_ENT_ASM_ID])
                                                    or (_row[CS_ENT_ASM_ID] != chain_id_src
                                                        and row[chain_id_col] == chain_id_src)):
                                                break
                                            _seq_id =\
                                                row[seq_id_col] if isinstance(row[seq_id_col], int) else int(row[seq_id_col])
                                            if _seq_id in (_row[CS_COMP_IDX], _row[CS_AUTH_SEQ]):
                                                break
                                            _row[CS_AMBIG_CODE] = ambig_code = 5
                                            ambig_code_5_test = True
                                            break
                                    if not ambig_code_5_test:
                                        for row in peripheral_src_rows():
                                            if row[comp_id_col] == comp_id\
                                               and row[atom_id_col][0] == atom_type\
                                               and row[ambig_code_col] == str(_row[CS_AMBIG_CODE])\
                                               or (_row[CS_AMBIG_CODE] != ambig_code_src
                                                   and row[ambig_code_col] == ambig_code_src):
                                                if not (row[chain_id_col] == str(_row[CS_ENT_ASM_ID])
                                                        or (_row[CS_ENT_ASM_ID] != chain_id_src
                                                            and row[chain_id_col] == chain_id_src)):
                                                    _row[CS_AMBIG_CODE] = ambig_code = 6
                                                    break

                            elif ambig_code == 6:
                                if len([item for item in entity_assembly
                                        if item['entity_type'] not in ('non-polymer', 'water')]) == 1\
                                   and len(entity_assembly[0]['label_asym_id'].split(',')) == 1:
                                    _row[CS_AMBIG_CODE] = ambig_code = 5

                            if ambig_code in (1, 2, 3):
                                if _row[CS_AMBIG_SET] is not None:
                                    _row[CS_AMBIG_SET] = None

                            if len_atom_ids > 1:
                                if _row[CS_AMBIG_CODE] == 1 or _row[CS_AMBIG_CODE] in EMPTY_VALUE:
                                    if _row[CS_ATOM_ID] not in methyl_atoms\
                                       or (_row[CS_ATOM_ID] in methyl_atoms
                                           and ((_row[CS_ATOM_TYPE][0] == 'H' and len_atom_ids == 6)
                                                or (_row[CS_ATOM_TYPE][0] == 'C' and len_atom_ids == 2))):
                                        _row[CS_AMBIG_CODE] =\
                                            self._reg.csStat.getMaxAmbigCodeWoSetId(comp_id, _row[CS_ATOM_ID], None)
                                __row = copy.copy(_row)
                                if fill_auth_atom_id:
                                    __row[CS_AUTH_ATOM] = __row[CS_ATOM_ID]

                                lp.add_data(__row)

                                for _atom_id in atom_ids[1:]:
                                    __row = copy.copy(_row)

                                    index += 1

                                    __row[CS_ID] = index
                                    __row[CS_ATOM_ID] = _atom_id
                                    if fill_auth_atom_id:
                                        __row[CS_AUTH_ATOM] = __row[CS_ATOM_ID]
                                    if fill_orig_atom_id and len(missing_ch3) > 0 and __row[CS_ORIG_ATOM] in EMPTY_VALUE:
                                        if _atom_id in methyl_atoms:
                                            if ch3_name_in_xplor and _atom_id[0] in PROTON_BEGIN_CODE:
                                                __row[CS_ORIG_ATOM] = __row[CS_ATOM_ID][-1] + __row[CS_ATOM_ID][:-1]
                                            else:
                                                __row[CS_ORIG_ATOM] = copy.copy(__row[CS_ATOM_ID])

                                    lp.add_data(__row)

                                index += 1

                                _row[CS_ATOM_ID] = atom_ids[-1]

                            if fill_auth_atom_id:
                                _row[CS_AUTH_ATOM] = _row[CS_ATOM_ID]
                            if fill_orig_atom_id and len(missing_ch3) > 0 and _row[CS_ORIG_ATOM] in EMPTY_VALUE:
                                if _row[CS_ATOM_ID] in methyl_atoms:
                                    if ch3_name_in_xplor and _row[CS_ATOM_ID][0] in PROTON_BEGIN_CODE:
                                        _row[CS_ORIG_ATOM] = _row[CS_ATOM_ID][-1] + _row[CS_ATOM_ID][:-1]
                                    else:
                                        _row[CS_ORIG_ATOM] = copy.copy(_row[CS_ATOM_ID])

                else:
                    _row[CS_COMP_ID] = comp_id
                    valid = True
                    missing_ch3 = []
                    if atom_id in self._reg.csStat.getMethylProtons(comp_id):
                        missing_ch3 = self._reg.csStat.getProtonsInSameGroup(comp_id, atom_id, True)
                        valid = self._reg.sail_flag
                        if 0 <= src_idx < len(src_lp):
                            row_src = src_lp.data[src_idx]
                            seq_id_src = row_src[seq_id_col]
                            for row in peripheral_src_rows():
                                if (row[seq_id_col] == str(_row[CS_COMP_IDX])
                                    or (_row[CS_COMP_IDX] != seq_id_src and row[seq_id_col] == seq_id_src)
                                    or (_row[CS_DETAILS] == 'UNMAPPED' and row[seq_id_col] == str(_row[CS_AUTH_SEQ])))\
                                   and row[comp_id_col].upper() == comp_id\
                                   and row[atom_id_col] in missing_ch3:
                                    valid = True
                                    missing_ch3.remove(row[atom_id_col])
                                    if len(missing_ch3) == 0:
                                        break
                    if len(missing_ch3) > 0 and (_row[CS_VAL] in EMPTY_VALUE or float(_row[CS_VAL]) >= 4.0):
                        heme = False
                        if _row[CS_VAL] not in EMPTY_VALUE:
                            if self._reg.ccU.updateChemCompDict(comp_id):
                                heme = comp_id == 'HEM' or 'HEME' in self._reg.ccU.lastChemCompDict['name']
                        if not heme:
                            missing_ch3 = []
                    _atom_id = atom_id
                    if not valid and len(missing_ch3) > 0:
                        atom_id = atom_id[:-1]
                        if _atom_id in self._reg.csStat.getRepMethylProtons(comp_id):
                            atom_id = _atom_id
                    if (valid or prefer_auth_atom_name or _row[CS_DETAILS] == 'UNMAPPED') and atom_id[0] not in ('Q', 'M'):
                        atom_ids = [atom_id]
                    else:
                        atom_ids = self._reg.dpV.getAtomIdListInXplor(comp_id, atom_id)
                        if len(atom_ids) == 0:
                            atom_ids = self._reg.dpV.getAtomIdListInXplor(comp_id,
                                                                          translateToStdAtomName(atom_id, comp_id,
                                                                                                 ccU=self._reg.ccU))
                    if valid and len(missing_ch3) > 0:
                        if not fill_orig_atom_id or not any(c in ('x', 'y', 'X', 'Y') for c in _row[CS_ORIG_ATOM])\
                           and len(self._reg.dpV.getAtomIdListInXplor(comp_id, _row[CS_ORIG_ATOM])) > 1\
                           and _row[CS_DETAILS] != 'UNMAPPED':
                            atom_ids = self._reg.dpV.getAtomIdListInXplor(comp_id, _row[CS_ORIG_ATOM])
                        else:
                            missing_ch3.clear()
                    if not valid and len(missing_ch3) > 0:
                        atom_ids.extend(missing_ch3)
                    len_atom_ids = len(atom_ids)
                    if len_atom_ids == 0:
                        _row[CS_ATOM_ID] = atom_id
                        if fill_auth_atom_id or _row[CS_ATOM_ID] != _row[CS_AUTH_ATOM]:
                            _row[CS_AUTH_ATOM] = _row[CS_ATOM_ID]
                        _row[CS_ATOM_TYPE] = 'H' if atom_id[0] in PSE_PRO_BEGIN_CODE else atom_id[0]
                        if _row[CS_ATOM_TYPE] in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS:
                            _row[CS_ISOTOPE] = ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS[_row[CS_ATOM_TYPE]][0]
                    else:
                        methyl_atoms = self._reg.csStat.getMethylAtoms(comp_id)
                        atom_ids = sorted(atom_ids)
                        _row[CS_ATOM_ID] = atom_ids[0]
                        _row[CS_AUTH_ATOM] = None
                        fill_auth_atom_id = _row[CS_AUTH_COMP] not in EMPTY_VALUE
                        if self._reg.ccU.updateChemCompDict(comp_id):
                            cca = next((cca for cca in self._reg.ccU.lastAtomDictList
                                        if cca['atom_id'] == _row[CS_ATOM_ID]), None)
                            if cca is not None:
                                _row[CS_ATOM_TYPE] = cca['type_symbol']
                                if _row[CS_ATOM_TYPE] in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS:
                                    _row[CS_ISOTOPE] = ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS[_row[CS_ATOM_TYPE]][0]
                            else:
                                _row[CS_ATOM_TYPE] = 'H' if _row[CS_ATOM_ID][0] in PROTON_BEGIN_CODE else atom_id[0]
                                if _row[CS_ATOM_TYPE] in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS:
                                    _row[CS_ISOTOPE] = ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS[_row[CS_ATOM_TYPE]][0]
                        else:
                            _row[CS_ATOM_TYPE] = 'H' if atom_id[0] in PSE_PRO_BEGIN_CODE else atom_id[0]
                            if _row[CS_ATOM_TYPE] in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS:
                                _row[CS_ISOTOPE] = ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS[_row[CS_ATOM_TYPE]][0]

                        ambig_code = _row[CS_AMBIG_CODE]
                        if ambig_code == 0:
                            _row[CS_AMBIG_CODE] = None

                        elif ambig_code in (2, 3):
                            _ambig_code = self._reg.csStat.getMaxAmbigCodeWoSetId(comp_id, _row[CS_ATOM_ID])
                            if _ambig_code not in (0, ambig_code):
                                if _ambig_code != 1:
                                    _row[CS_AMBIG_CODE] = _ambig_code
                                else:
                                    _row[CS_AMBIG_CODE] = ambig_code = 4
                                    if 0 <= _src_idx < len(src_lp):
                                        row_src = src_lp.data[_src_idx]
                                        chain_id_src = row_src[chain_id_col]
                                        seq_id_src = row_src[seq_id_col]
                                        atom_type = row_src[atom_id_col][0]
                                        val = float(row_src[val_col])
                                        sig = self._reg.ccU.getBondSignature(comp_id, atom_id)
                                        for row in peripheral_src_rows():
                                            if row[chain_id_col] == chain_id_src\
                                               and row[seq_id_col] == seq_id_src\
                                               and row[comp_id_col] == comp_id\
                                               and row[atom_id_col][0] == atom_type\
                                               and abs(float(row[val_col]) - val) < 1.0\
                                               and self._reg.ccU.getBondSignature(comp_id, row[atom_id_col]) == sig:
                                                row[ambig_code_col] = '4'
                                                reparse = True

                        elif ambig_code == 4:
                            if not self._reg.annotation_mode and _row[CS_DETAILS] != 'UNMAPPED':
                                row_src = src_lp.data[_src_idx]
                                chain_id_src = row_src[chain_id_col]
                                atom_id_src = row_src[atom_id_col]
                                atom_type = atom_id_src[0]
                                ambig_code_src = row_src[ambig_code_col]
                                atom_ids_in_group = self._reg.ccU.getProtonsInSameGroup(comp_id, atom_id_src)\
                                    if atom_type in PROTON_BEGIN_CODE else []
                                ambig_code_4_test = hetero_group_test = False
                                for row in peripheral_src_rows():
                                    if row[comp_id_col] == comp_id\
                                       and row[atom_id_col][0] == atom_type\
                                       and row[ambig_code_col] == str(_row[CS_AMBIG_CODE])\
                                       or (_row[CS_AMBIG_CODE] != ambig_code_src
                                           and row[ambig_code_col] == ambig_code_src):
                                        if not (row[chain_id_col] == str(_row[CS_ENT_ASM_ID])
                                                or (_row[CS_ENT_ASM_ID] != chain_id_src
                                                    and row[chain_id_col] == chain_id_src)):
                                            break
                                        _seq_id = row[seq_id_col] if isinstance(row[seq_id_col], int) else int(row[seq_id_col])
                                        if _seq_id in (_row[CS_COMP_IDX], _row[CS_AUTH_SEQ]):
                                            ambig_code_4_test = True
                                            if row[atom_id_col] not in atom_ids_in_group:
                                                hetero_group_test = True
                                                break
                                if not ambig_code_4_test:
                                    ambig_code_5_test = False
                                    for row in peripheral_src_rows():
                                        if row[comp_id_col] == comp_id\
                                           and row[atom_id_col][0] == atom_type\
                                           and row[ambig_code_col] == str(_row[CS_AMBIG_CODE])\
                                           or (_row[CS_AMBIG_CODE] != ambig_code_src
                                               and row[ambig_code_col] == ambig_code_src):
                                            if not (row[chain_id_col] == str(_row[CS_ENT_ASM_ID])
                                                    or (_row[CS_ENT_ASM_ID] != chain_id_src
                                                        and row[chain_id_col] == chain_id_src)):
                                                break
                                            _seq_id =\
                                                row[seq_id_col] if isinstance(row[seq_id_col], int) else int(row[seq_id_col])
                                            if _seq_id in (_row[CS_COMP_IDX], _row[CS_AUTH_SEQ]):
                                                break
                                            _row[CS_AMBIG_CODE] = ambig_code = 5
                                            ambig_code_5_test = True
                                            break
                                    if not ambig_code_5_test:
                                        for row in peripheral_src_rows():
                                            if row[comp_id_col] == comp_id\
                                               and row[atom_id_col][0] == atom_type\
                                               and row[ambig_code_col] == str(_row[CS_AMBIG_CODE])\
                                               or (_row[CS_AMBIG_CODE] != ambig_code_src
                                                   and row[ambig_code_col] == ambig_code_src):
                                                if not (row[chain_id_col] == str(_row[CS_ENT_ASM_ID])
                                                        or (_row[CS_ENT_ASM_ID] != chain_id_src
                                                            and row[chain_id_col] == chain_id_src)):
                                                    _row[CS_AMBIG_CODE] = ambig_code = 6
                                                    break
                                        if ambig_code == 4:
                                            _row[CS_AMBIG_CODE] = ambig_code = 1
                                elif not hetero_group_test:
                                    _row[CS_AMBIG_CODE] = ambig_code = 1

                        elif ambig_code == 5:
                            if not self._reg.annotation_mode and _row[CS_DETAILS] != 'UNMAPPED':
                                row_src = src_lp.data[_src_idx]
                                chain_id_src = row_src[chain_id_col]
                                atom_type = row_src[atom_id_col][0]
                                ambig_code_src = row_src[ambig_code_col]
                                ambig_code_5_test = False
                                for row in peripheral_src_rows():
                                    if row[comp_id_col] == comp_id\
                                       and row[atom_id_col][0] == atom_type\
                                       and row[ambig_code_col] == str(_row[CS_AMBIG_CODE])\
                                       or (_row[CS_AMBIG_CODE] != ambig_code_src
                                           and row[ambig_code_col] == ambig_code_src):
                                        if not (row[chain_id_col] == str(_row[CS_ENT_ASM_ID])
                                                or (_row[CS_ENT_ASM_ID] != chain_id_src
                                                    and row[chain_id_col] == chain_id_src)):
                                            break
                                        _seq_id = row[seq_id_col] if isinstance(row[seq_id_col], int) else int(row[seq_id_col])
                                        if _seq_id in (_row[CS_COMP_IDX], _row[CS_AUTH_SEQ]):
                                            break
                                        _row[CS_AMBIG_CODE] = ambig_code = 5
                                        ambig_code_5_test = True
                                        break
                                if not ambig_code_5_test:
                                    for row in peripheral_src_rows():
                                        if row[comp_id_col] == comp_id\
                                           and row[atom_id_col][0] == atom_type\
                                           and row[ambig_code_col] == str(_row[CS_AMBIG_CODE])\
                                           or (_row[CS_AMBIG_CODE] != ambig_code_src
                                               and row[ambig_code_col] == ambig_code_src):
                                            if not (row[chain_id_col] == str(_row[CS_ENT_ASM_ID])
                                                    or (_row[CS_ENT_ASM_ID] != chain_id_src
                                                        and row[chain_id_col] == chain_id_src)):
                                                _row[CS_AMBIG_CODE] = ambig_code = 6
                                                break

                        elif ambig_code == 6:
                            if len([item for item in entity_assembly
                                    if item['entity_type'] not in ('non-polymer', 'water')]) == 1\
                               and len(entity_assembly[0]['label_asym_id'].split(',')) == 1:
                                _row[CS_AMBIG_CODE] = ambig_code = 5

                        if ambig_code in (1, 2, 3):
                            if _row[CS_AMBIG_SET] is not None:
                                _row[CS_AMBIG_SET] = None

                        if len_atom_ids > 1:
                            if _row[CS_AMBIG_CODE] == 1 or _row[CS_AMBIG_CODE] in EMPTY_VALUE:
                                if _row[CS_ATOM_ID] not in methyl_atoms\
                                   or (_row[CS_ATOM_ID] in methyl_atoms
                                       and ((_row[CS_ATOM_TYPE][0] == 'H' and len_atom_ids == 6)
                                            or (_row[CS_ATOM_TYPE][0] == 'C' and len_atom_ids == 2))):
                                    _row[CS_AMBIG_CODE] = self._reg.csStat.getMaxAmbigCodeWoSetId(comp_id, _row[CS_ATOM_ID], None)
                            __row = copy.copy(_row)
                            if fill_auth_atom_id:
                                __row[CS_AUTH_ATOM] = __row[CS_ATOM_ID]
                            lp.add_data(__row)

                            for _atom_id in atom_ids[1:]:
                                __row = copy.copy(_row)

                                index += 1

                                __row[CS_ID] = index
                                __row[CS_ATOM_ID] = _atom_id
                                if fill_auth_atom_id:
                                    __row[CS_AUTH_ATOM] = __row[CS_ATOM_ID]
                                if fill_orig_atom_id and len(missing_ch3) > 0\
                                   and __row[CS_ORIG_ATOM] in EMPTY_VALUE:
                                    if _atom_id in methyl_atoms:
                                        if ch3_name_in_xplor and _atom_id[0] in PROTON_BEGIN_CODE:
                                            __row[CS_ORIG_ATOM] = __row[CS_ATOM_ID][-1] + __row[CS_ATOM_ID][:-1]
                                        else:
                                            __row[CS_ORIG_ATOM] = copy.copy(__row[CS_ATOM_ID])

                                lp.add_data(__row)

                            index += 1

                            _row[CS_ATOM_ID] = atom_ids[-1]

                        if fill_auth_atom_id:
                            _row[CS_AUTH_ATOM] = _row[CS_ATOM_ID]
                        if fill_orig_atom_id and len(missing_ch3) > 0\
                           and _row[CS_ORIG_ATOM] in EMPTY_VALUE:
                            if _row[CS_ATOM_ID] in methyl_atoms:
                                if ch3_name_in_xplor and _row[CS_ATOM_ID][0] in PROTON_BEGIN_CODE:
                                    _row[CS_ORIG_ATOM] = _row[CS_ATOM_ID][-1] + _row[CS_ATOM_ID][:-1]
                                else:
                                    _row[CS_ORIG_ATOM] = copy.copy(_row[CS_ATOM_ID])

                # Canonical atom resolution must not erase a supplied author label.
                # Chain reconstruction and methyl validation may have populated it
                # with a coordinate name. Restore only unchanged, unexpanded source
                # observations; annotation and real atom remapping retain their rules.
                source_row = src_lp.data[_src_idx]
                if (
                    not self._reg.annotation_mode and index == original_index and auth_atom_id_col != -1
                    and source_row[auth_atom_id_col] not in EMPTY_VALUE
                    and source_row[atom_id_col] == _row[CS_ATOM_ID]
                    and source_row[seq_id_col] == str(_row[CS_COMP_IDX])
                    and source_row[comp_id_col] == _row[CS_COMP_ID]
                ):
                    _row[CS_AUTH_ATOM] = source_row[auth_atom_id_col]
                    if not has_orig_seq:
                        _row[CS_ORIG_ATOM] = source_row[auth_atom_id_col]

                return index, _row, reparse

            copied_auth_chain_ids = set()
            copied_chain_ids = set()
            copied_auth_asym_id_mapping = {}

            if has_auth_seq:
                auth_asym_ids = [row[0] for row in auth_dat]

                common_auth_asym_ids = collections.Counter(auth_asym_ids).most_common()

                if len(common_auth_asym_ids) > 1:
                    auth_cs_tags = ['Auth_asym_ID', 'Auth_seq_ID', 'Auth_comp_ID', 'Auth_atom_ID', 'Val']

                    _common_auth_asym_ids = dict(common_auth_asym_ids)

                    for _auth_chain_id_1, _auth_chain_id_2 in itertools.combinations(_common_auth_asym_ids.keys(), 2):

                        if _common_auth_asym_ids[_auth_chain_id_1] != _common_auth_asym_ids[_auth_chain_id_2]:
                            continue

                        try:
                            _auth_seq_id_1, _auth_comp_id_1 =\
                                next((int(row[1]), row[2]) for row in auth_dat if row[0] == _auth_chain_id_1)
                            _auth_seq_id_2, _auth_comp_id_2 =\
                                next((int(row[1]), row[2]) for row in auth_dat if row[0] == _auth_chain_id_2)
                        except (ValueError, TypeError):
                            continue

                        _seq_key_1 = (_auth_chain_id_1, _auth_seq_id_1, _auth_comp_id_1)
                        _seq_key_2 = (_auth_chain_id_2, _auth_seq_id_2, _auth_comp_id_2)

                        if _seq_key_1 not in auth_to_entity_type or _seq_key_2 not in auth_to_entity_type:
                            continue

                        if auth_to_entity_type[_seq_key_1] != auth_to_entity_type[_seq_key_2]\
                           or auth_to_entity_type[_seq_key_1] in ('non-polymer', 'water'):
                            continue

                        _auth_cs_1 = [row[1:] for row in loop.get_tag(auth_cs_tags) if row[0] == _auth_chain_id_1]
                        _auth_cs_2 = [row[1:] for row in loop.get_tag(auth_cs_tags) if row[0] == _auth_chain_id_2]

                        _auth_cs_1 = sorted(_auth_cs_1, key=itemgetter(0, 2))
                        _auth_cs_2 = sorted(_auth_cs_2, key=itemgetter(0, 2))

                        if _auth_cs_1 == _auth_cs_2:
                            copied_auth_chain_ids.add(_auth_chain_id_2)

                            # DAOTHER-10898
                            if _auth_chain_id_1 not in copied_auth_asym_id_mapping:
                                copied_auth_asym_id_mapping[_auth_chain_id_1] = [_auth_chain_id_1]
                            if _auth_chain_id_2 not in copied_auth_asym_id_mapping[_auth_chain_id_1]:
                                copied_auth_asym_id_mapping[_auth_chain_id_1].append(_auth_chain_id_2)

                # DAOTHER-10898
                if len(copied_auth_chain_ids) > 0:
                    for k, v in copied_auth_asym_id_mapping.items():
                        copied_auth_asym_id_mapping[k] = ','.join(sorted(v))

            else:

                tags = ['Entity_assembly_ID', 'Comp_index_ID', 'Comp_ID', 'Atom_ID']
                dat = loop.get_tag(tags)

                chain_ids = [row[0] for row in dat]

                common_chain_ids = collections.Counter(chain_ids).most_common()

                if len(common_chain_ids) > 1:
                    cs_tags = ['Entity_assembly_ID', 'Comp_index_ID', 'Comp_ID', 'Atom_ID', 'Val']

                    _common_chain_ids = dict(common_chain_ids)

                    for _chain_id_1, _chain_id_2 in itertools.combinations(_common_chain_ids, 2):

                        if _common_chain_ids[_chain_id_1] != _common_chain_ids[_chain_id_2]:
                            continue

                        _cs_1 = [row[1:] for row in loop.get_tag(cs_tags) if row[0] == _chain_id_1]
                        _cs_2 = [row[1:] for row in loop.get_tag(cs_tags) if row[0] == _chain_id_2]

                        _cs_1 = sorted(_cs_1, key=itemgetter(0, 2))
                        _cs_2 = sorted(_cs_2, key=itemgetter(0, 2))

                        if _cs_1 == _cs_2:
                            copied_chain_ids.add(_chain_id_2)

            if has_orig_seq:
                orig_asym_id_col = loop.tags.index('Original_PDB_strand_ID')
                orig_seq_id_col = loop.tags.index('Original_PDB_residue_no')
                orig_comp_id_col = loop.tags.index('Original_PDB_residue_name')
                orig_atom_id_col = loop.tags.index('Original_PDB_atom_name')

            lp = pynmrstar.Loop.from_scratch(lp_category)

            tags = [f'{lp_category}.{item}' for item in items]

            lp.add_tag(tags)

            prefer_auth_atom_name = False

            if (self._reg.annotation_mode or self._reg.native_combined) and len(auth_atom_name_to_id) > 0:

                def get_auth_seq_id(val):
                    if isinstance(val, int):
                        return val
                    if val in EMPTY_VALUE:
                        return None
                    if val.isdigit():
                        return int(val)
                    return int(re.findall(r'\d+', val)[0])

                count_auth_name = count_auth_id = 0

                for row in loop:

                    auth_asym_id = row[auth_asym_id_col]
                    auth_seq_id = get_auth_seq_id(row[auth_seq_id_col])
                    auth_comp_id = row[auth_comp_id_col]
                    auth_atom_id = row[auth_atom_id_col]

                    if auth_seq_id is None:
                        continue

                    seq_key = (auth_asym_id, auth_seq_id, auth_comp_id)
                    try:
                        auth_to_star_seq[seq_key]  # pylint: disable=pointless-statement
                    except KeyError:
                        auth_asym_id = next((_auth_asym_id for _auth_asym_id, _auth_seq_id, _auth_comp_id in auth_to_star_seq
                                             if _auth_seq_id == auth_seq_id and _auth_comp_id == auth_comp_id), auth_asym_id)
                        if (auth_asym_id, auth_seq_id, auth_comp_id) not in auth_to_star_seq:
                            auth_comp_id = next((_auth_comp_id for _auth_asym_id, _auth_seq_id, _auth_comp_id in auth_to_star_seq
                                                 if _auth_asym_id == auth_asym_id and _auth_seq_id == auth_seq_id), auth_comp_id)

                    if auth_comp_id in auth_atom_name_to_id:
                        if auth_atom_id in auth_atom_name_to_id[auth_comp_id]:
                            count_auth_name += 1
                        if auth_atom_id in auth_atom_name_to_id[auth_comp_id].values():
                            count_auth_id += 1

                if count_auth_name + count_auth_id == 0:

                    for row in loop:

                        auth_asym_id = row[auth_asym_id_col]
                        auth_seq_id = get_auth_seq_id(row[auth_seq_id_col])
                        auth_comp_id = row[auth_comp_id_col]
                        auth_atom_id = row[auth_atom_id_col]

                        if auth_seq_id is None:
                            continue

                        seq_key = (auth_asym_id, auth_seq_id, auth_comp_id)
                        try:
                            auth_to_star_seq_ann[seq_key]  # pylint: disable=pointless-statement
                            _seq_key = seq_key if seq_key in coord_atom_site else (seq_key[0], seq_key[1])
                            if _seq_key in coord_atom_site:  # DAOTHER-8817
                                auth_comp_id = coord_atom_site[_seq_key]['comp_id']
                        except KeyError:
                            continue

                        if auth_comp_id in auth_atom_name_to_id:
                            if auth_atom_id in auth_atom_name_to_id[auth_comp_id]:
                                count_auth_name += 1
                            if auth_atom_id in auth_atom_name_to_id[auth_comp_id].values():
                                count_auth_id += 1

                prefer_auth_atom_name = count_auth_name > count_auth_id

            has_genuine_ambig_code = False

            can_auth_asym_id_mapping = {}  # DAOTHER-8751
            seq_id_offset_for_unmapped = {}  # DAOTHER-9065
            label_seq_id_offset_for_extended = {}  # D_1300044764

            while True:

                reparse_request = can_auth_asym_id_mapping_failed = False  # DAOTHER-9065, DAOTHER-9158

                lp.clear_data()

                index = 1

                for idx, row in enumerate(loop):

                    _row = [None] * len(tags)

                    comp_id = _orig_comp_id = row[comp_id_col].upper()
                    _orig_atom_id = row[atom_id_col]
                    atom_id = _orig_atom_id.upper()

                    _row[CS_VAL] = row[val_col]

                    try:
                        float(_row[CS_VAL])
                    except ValueError:
                        continue

                    if val_err_col != -1:
                        val_err = row[val_err_col]
                        _row[CS_VAL_ERR] = val_err

                        if val_err not in EMPTY_VALUE:
                            try:
                                _val_err = float(val_err)
                                if _val_err < 0.0:
                                    _row[CS_VAL_ERR] = abs(_val_err)
                            except ValueError:
                                pass

                    if fig_of_merit_col != -1:
                        _row[CS_FIG_MERIT] = row[fig_of_merit_col]

                    if ambig_code_col != -1:
                        ambig_code = row[ambig_code_col]
                        if ambig_code not in EMPTY_VALUE:
                            try:
                                ambig_code = int(ambig_code) if isinstance(ambig_code, str) else ambig_code
                                if ambig_code in ALLOWED_AMBIGUITY_CODES:
                                    _row[CS_AMBIG_CODE] = ambig_code
                                else:
                                    _row[CS_AMBIG_CODE] = None
                            except ValueError:
                                _row[CS_AMBIG_CODE] = None

                    if ambig_set_id_col != -1:
                        ambig_set_id = row[ambig_set_id_col]
                        if ambig_set_id not in EMPTY_VALUE:
                            try:
                                ambig_set_id = int(ambig_set_id)
                                if ambig_set_id > 0:
                                    _row[CS_AMBIG_SET] = ambig_set_id
                            except ValueError:
                                _row[CS_AMBIG_SET] = None

                    if occupancy_col != -1:
                        try:
                            occupancy = row[occupancy_col]
                        except IndexError:
                            occupancy = '.'
                        if occupancy not in EMPTY_VALUE:
                            try:
                                occupancy = float(occupancy)
                                if occupancy >= 0.0:
                                    _row[CS_OCCUPANCY] = occupancy
                            except ValueError:
                                pass

                    if reson_id_col != -1:
                        reson_id = row[reson_id_col]
                        if reson_id not in EMPTY_VALUE:
                            try:
                                reson_id = int(reson_id)
                                if reson_id > 0:
                                    _row[CS_RESONANCE] = reson_id
                            except ValueError:
                                pass

                    if has_auth_seq:

                        if row[auth_asym_id_col] in copied_auth_chain_ids:
                            continue

                        _row[CS_AUTH_ASYM], _row[CS_AUTH_SEQ], _row[CS_AUTH_COMP], _row[CS_AUTH_ATOM] =\
                            row[auth_asym_id_col], row[auth_seq_id_col], \
                            row[auth_comp_id_col], row[auth_atom_id_col]

                    elif self._reg.bmrb_only and self._reg.internal_mode:
                        if auth_seq_id_col != -1 and auth_comp_id_col != -1 and auth_atom_id_col != -1:
                            _row[CS_AUTH_SEQ], _row[CS_AUTH_COMP], _row[CS_AUTH_ATOM] =\
                                row[auth_seq_id_col], row[auth_comp_id_col], row[auth_atom_id_col]

                    if has_orig_seq:
                        _row[CS_ORIG_ASYM], _row[CS_ORIG_SEQ], _row[CS_ORIG_COMP], _row[CS_ORIG_ATOM] =\
                            row[orig_asym_id_col], row[orig_seq_id_col], \
                            row[orig_comp_id_col], row[orig_atom_id_col]

                    if details_col != -1:
                        _row[CS_DETAILS] = row[details_col]

                    _row[CS_ENTRY_ID], _row[CS_LIST_ID] = self._reg.entry_id, list_id

                    resolved = True

                    if has_auth_seq and len(auth_to_star_seq) > 0:
                        auth_asym_id = row[auth_asym_id_col]
                        auth_seq_id = row[auth_seq_id_col]
                        auth_comp_id = row[auth_comp_id_col]

                        if valid_auth_seq and auth_seq_id not in EMPTY_VALUE:
                            auth_seq_id_ = int(auth_seq_id)
                            seq_key = (auth_asym_id, auth_seq_id_, auth_comp_id)
                            _seq_key = (seq_key[0], seq_key[1])
                            try:
                                entity_assembly_id, seq_id, entity_id, _ = auth_to_star_seq[seq_key]
                                if atom_id != _row[CS_AUTH_ATOM]:
                                    if _seq_key in coord_atom_site:
                                        _coord_atom_site = coord_atom_site[_seq_key]
                                        if atom_id in _coord_atom_site['atom_id']:
                                            _row[CS_AUTH_ATOM] = atom_id
                            except KeyError:
                                entity_assembly_id = None
                                if self._reg.annotation_mode or self._reg.native_combined:
                                    auth_asym_id =\
                                        next((_auth_asym_id for _auth_asym_id, _auth_seq_id, _auth_comp_id in auth_to_star_seq
                                              if _auth_seq_id == auth_seq_id_ and _auth_comp_id == auth_comp_id), auth_asym_id)
                                    seq_key = (auth_asym_id, auth_seq_id_, auth_comp_id)
                                    if seq_key in auth_to_star_seq:
                                        _row[CS_AUTH_ASYM] = row[auth_asym_id_col] = auth_asym_id
                                        if has_orig_seq:
                                            _row[CS_ORIG_ASYM] = row[orig_asym_id_col] = auth_asym_id
                                        else:
                                            _row[CS_ORIG_ASYM] = auth_asym_id
                                        _seq_key = (seq_key[0], seq_key[1])
                                        entity_assembly_id, seq_id, entity_id, _ = auth_to_star_seq[seq_key]
                                    else:
                                        auth_asym_id, auth_comp_id =\
                                            next(((_auth_asym_id, _auth_comp_id)
                                                  for _auth_asym_id, _auth_seq_id, _auth_comp_id in auth_to_star_seq
                                                  if _auth_seq_id == auth_seq_id_), (auth_asym_id, auth_comp_id))
                                        seq_key = (auth_asym_id, auth_seq_id_, auth_comp_id)
                                        if seq_key in auth_to_star_seq:
                                            _row[CS_AUTH_ASYM] = row[auth_asym_id_col] = auth_asym_id
                                            if has_orig_seq:
                                                _row[CS_ORIG_ASYM] = row[orig_asym_id_col] = auth_asym_id
                                            else:
                                                _row[CS_ORIG_ASYM] = auth_asym_id
                                            _row[CS_COMP_ID] = row[comp_id_col] = auth_comp_id
                                            _row[CS_AUTH_COMP] = row[auth_comp_id_col] = auth_comp_id
                                            comp_id = auth_comp_id
                                            _seq_key = (seq_key[0], seq_key[1])
                                            entity_assembly_id, seq_id, entity_id, _ = auth_to_star_seq[seq_key]
                                if seq_key not in auth_to_star_seq:
                                    auth_comp_id =\
                                        next((_auth_comp_id for _auth_asym_id, _auth_seq_id, _auth_comp_id in auth_to_star_seq
                                              if _auth_asym_id == auth_asym_id and _auth_seq_id == auth_seq_id_), auth_comp_id)
                                    comp_id = _row[CS_AUTH_COMP] = auth_comp_id
                                    seq_key = (auth_asym_id, auth_seq_id_, auth_comp_id)
                                    if seq_key in auth_to_star_seq:
                                        entity_assembly_id, seq_id, entity_id, _ = auth_to_star_seq[seq_key]
                                    elif seq_key in auth_to_star_seq_ann:
                                        entity_assembly_id, seq_id, entity_id, _ = auth_to_star_seq_ann[seq_key]

                            if entity_assembly_id is not None:
                                self._reg.ent_asym_id_with_exptl_data.add(entity_assembly_id)
                                _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID] = entity_assembly_id, entity_id
                                _row[CS_COMP_IDX] = _row[CS_SEQ_ID] = seq_id

                            if prefer_auth_atom_name:
                                if has_orig_seq:
                                    orig_atom_id = _row[CS_ORIG_ATOM]
                                _atom_id = atom_id
                                _seq_key = seq_key if seq_key in coord_atom_site else (seq_key[0], seq_key[1])
                                if _seq_key in coord_atom_site:
                                    _coord_atom_site = coord_atom_site[_seq_key]
                                    if comp_id in auth_atom_name_to_id and comp_id == _coord_atom_site['comp_id']:
                                        if _atom_id in auth_atom_name_to_id[comp_id]:
                                            if auth_atom_name_to_id[comp_id][_atom_id] in _coord_atom_site['atom_id']:
                                                _row[CS_AUTH_ATOM] = atom_id = auth_atom_name_to_id[comp_id][_atom_id]
                                            elif 'split_comp_id' not in _coord_atom_site and has_orig_seq\
                                                    and orig_atom_id in auth_atom_name_to_id[comp_id]:
                                                _row[CS_AUTH_ATOM] = atom_id = auth_atom_name_to_id[comp_id][orig_atom_id]
                                    if 'alt_atom_id' in _coord_atom_site and _atom_id in _coord_atom_site['alt_atom_id']\
                                       and comp_id == _coord_atom_site['comp_id']:
                                        _row[CS_AUTH_ATOM] = atom_id =\
                                            _coord_atom_site['atom_id'][_coord_atom_site['alt_atom_id'].index(_atom_id)]
                                    # DAOTHER-8751, 8817 (D_1300043061)
                                    elif 'alt_comp_id' in _coord_atom_site and 'alt_atom_id' in _coord_atom_site\
                                         and _atom_id in _coord_atom_site['alt_atom_id']\
                                         and comp_id == _coord_atom_site['alt_comp_id'][_coord_atom_site['alt_atom_id'].index(_atom_id)]:  # noqa: E501, pylint: disable=line-too-long
                                        _row[CS_AUTH_COMP] = comp_id
                                        # Entity_assembly_ID, Entity_ID, Comp_index_ID, Seq_ID, Comp_ID, Auth_asym_ID, Auth_seq_ID
                                        cca_row =\
                                            next((cca_row for cca_row in self._reg.chem_comp_asm_dat
                                                  if cca_row[CCA_COMP_ID] == comp_id and cca_row[CCA_SEQ_ID] == _seq_key[0]
                                                  and cca_row[CCA_AUTH_ASYM] == _seq_key[1]), None)
                                        if cca_row is not None:
                                            _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID], _row[CS_COMP_IDX], _row[CS_SEQ_ID] =\
                                                cca_row[CCA_ENT_ASM_ID], cca_row[CCA_ENTITY_ID], \
                                                cca_row[CCA_COMP_IDX], cca_row[CCA_SEQ_ID]
                                        if comp_id in auth_atom_name_to_id_ext and _atom_id in auth_atom_name_to_id_ext[comp_id]\
                                           and len(set(_coord_atom_site['alt_comp_id'])) > 1:
                                            _row[CS_AUTH_ATOM] = atom_id = auth_atom_name_to_id_ext[comp_id][_atom_id]
                                        else:
                                            _row[CS_AUTH_ATOM] = atom_id =\
                                                _coord_atom_site['atom_id'][_coord_atom_site['alt_atom_id'].index(_atom_id)]
                                    elif 'split_comp_id' in _coord_atom_site:
                                        for _comp_id in _coord_atom_site['split_comp_id']:
                                            if _comp_id == comp_id:
                                                continue
                                            __seq_key = (_seq_key[0], _seq_key[1], _comp_id)
                                            __coord_atom_site = coord_atom_site[__seq_key]
                                            if __coord_atom_site is None:
                                                continue
                                            if 'alt_comp_id' in __coord_atom_site and 'alt_atom_id' in __coord_atom_site\
                                               and _atom_id in __coord_atom_site['alt_atom_id']:
                                                comp_id = _comp_id
                                                _row[CS_AUTH_COMP] = comp_id
                                                # Entity_assembly_ID, Entity_ID, Comp_index_ID,
                                                # Seq_ID, Comp_ID, Auth_asym_ID, Auth_seq_ID
                                                cca_row =\
                                                    next((cca_row for cca_row in self._reg.chem_comp_asm_dat
                                                          if cca_row[CCA_COMP_ID] == comp_id and cca_row[CCA_SEQ_ID] == _seq_key[0]
                                                          and cca_row[CCA_AUTH_ASYM] == _seq_key[1]), None)
                                                if cca_row is not None:
                                                    _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID], _row[CS_COMP_IDX], _row[CS_SEQ_ID] =\
                                                        cca_row[CCA_ENT_ASM_ID], cca_row[CCA_ENTITY_ID], \
                                                        cca_row[CCA_COMP_IDX], cca_row[CCA_SEQ_ID]
                                                _row[CS_AUTH_ATOM] = atom_id =\
                                                    __coord_atom_site['atom_id'][__coord_atom_site['alt_atom_id'].index(_atom_id)]
                                                _seq_key = __seq_key
                                                break
                                            if _atom_id in __coord_atom_site['atom_id']:
                                                comp_id = _comp_id
                                                _row[CS_AUTH_COMP] = comp_id
                                                # Entity_assembly_ID, Entity_ID, Comp_index_ID,
                                                # Seq_ID, Comp_ID, Auth_asym_ID, Auth_seq_ID
                                                cca_row =\
                                                    next((cca_row for cca_row in self._reg.chem_comp_asm_dat
                                                          if cca_row[CCA_COMP_ID] == comp_id and cca_row[CCA_SEQ_ID] == _seq_key[0]
                                                          and cca_row[CCA_AUTH_ASYM] == _seq_key[1]), None)
                                                if cca_row is not None:
                                                    _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID], _row[CS_COMP_IDX], _row[CS_SEQ_ID] =\
                                                        cca_row[CCA_ENT_ASM_ID], cca_row[CCA_ENTITY_ID], \
                                                        cca_row[CCA_COMP_IDX], cca_row[CCA_SEQ_ID]
                                                _seq_key = __seq_key
                                                break

                            apply_orig_seq_naming(std_res_branch=True)

                        elif auth_asym_id not in EMPTY_VALUE and auth_seq_id not in EMPTY_VALUE and auth_comp_id not in EMPTY_VALUE:

                            try:
                                _auth_seq_id = int(auth_seq_id)
                                seq_key = (auth_asym_id, _auth_seq_id, auth_comp_id)
                                _seq_key = (seq_key[0], seq_key[1])
                                if seq_key in auth_to_star_seq:
                                    entity_assembly_id, seq_id, entity_id, _ = auth_to_star_seq[seq_key]
                                    self._reg.ent_asym_id_with_exptl_data.add(entity_assembly_id)
                                    _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID] = entity_assembly_id, entity_id
                                    _row[CS_COMP_IDX] = _row[CS_SEQ_ID] = seq_id

                                    apply_orig_seq_naming()

                                else:
                                    resolved = False

                            except ValueError:
                                resolved = False

                        else:
                            resolved = False

                    # DAOTHER-9281
                    elif has_auth_chain and (valid_auth_seq or row[auth_asym_id_col] in valid_auth_seq_per_chain):
                        auth_asym_id = row[auth_asym_id_col]
                        auth_seq_id = row[aux_auth_seq_id_col]
                        auth_comp_id = row[aux_auth_comp_id_col]

                        _row[CS_AUTH_SEQ] = auth_seq_id

                        auth_seq_id_ = int(auth_seq_id)
                        seq_key = (auth_asym_id, auth_seq_id_, auth_comp_id)
                        _seq_key = (seq_key[0], seq_key[1])
                        try:
                            entity_assembly_id, seq_id, entity_id, _ = auth_to_star_seq[seq_key]
                            if atom_id != _row[CS_AUTH_ATOM]:
                                if _seq_key in coord_atom_site:
                                    _coord_atom_site = coord_atom_site[_seq_key]
                                    if atom_id in _coord_atom_site['atom_id']:
                                        _row[CS_AUTH_ATOM] = atom_id
                        except KeyError:
                            entity_assembly_id = None
                            if self._reg.annotation_mode or self._reg.native_combined:
                                auth_asym_id =\
                                    next((_auth_asym_id for _auth_asym_id, _auth_seq_id, _auth_comp_id in auth_to_star_seq
                                          if _auth_seq_id == auth_seq_id_ and _auth_comp_id == auth_comp_id), auth_asym_id)
                                seq_key = (auth_asym_id, auth_seq_id_, auth_comp_id)
                                if seq_key in auth_to_star_seq:
                                    _row[CS_AUTH_ASYM] = row[auth_asym_id_col] = auth_asym_id
                                    if has_orig_seq:
                                        _row[CS_ORIG_ASYM] = row[orig_asym_id_col] = auth_asym_id
                                    else:
                                        _row[CS_ORIG_ASYM] = auth_asym_id
                                    _seq_key = (seq_key[0], seq_key[1])
                                    entity_assembly_id, seq_id, entity_id, _ = auth_to_star_seq[seq_key]
                                else:
                                    auth_asym_id, auth_comp_id =\
                                        next(((_auth_asym_id, _auth_comp_id)
                                              for _auth_asym_id, _auth_seq_id, _auth_comp_id in auth_to_star_seq
                                              if _auth_seq_id == auth_seq_id_), (auth_asym_id, auth_comp_id))
                                    seq_key = (auth_asym_id, auth_seq_id_, auth_comp_id)
                                    if seq_key in auth_to_star_seq:
                                        _row[CS_AUTH_ASYM] = row[auth_asym_id_col] = auth_asym_id
                                        if has_orig_seq:
                                            _row[CS_ORIG_ASYM] = row[orig_asym_id_col] = auth_asym_id
                                        else:
                                            _row[CS_ORIG_ASYM] = auth_asym_id
                                        _row[CS_COMP_ID] = row[comp_id_col] = auth_comp_id
                                        _row[CS_AUTH_COMP] = row[auth_comp_id_col] = auth_comp_id
                                        comp_id = auth_comp_id
                                        _seq_key = (seq_key[0], seq_key[1])
                                        entity_assembly_id, seq_id, entity_id, _ = auth_to_star_seq[seq_key]
                            if seq_key not in auth_to_star_seq:
                                auth_comp_id =\
                                    next((_auth_comp_id for _auth_asym_id, _auth_seq_id, _auth_comp_id in auth_to_star_seq
                                          if _auth_asym_id == auth_asym_id and _auth_seq_id == auth_seq_id_), auth_comp_id)
                                comp_id = _row[CS_AUTH_COMP] = auth_comp_id
                                seq_key = (auth_asym_id, auth_seq_id_, auth_comp_id)
                                if seq_key in auth_to_star_seq:
                                    entity_assembly_id, seq_id, entity_id, _ = auth_to_star_seq[seq_key]
                                elif seq_key in auth_to_star_seq_ann:
                                    entity_assembly_id, seq_id, entity_id, _ = auth_to_star_seq_ann[seq_key]

                        if entity_assembly_id is not None:
                            self._reg.ent_asym_id_with_exptl_data.add(entity_assembly_id)
                            _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID] = entity_assembly_id, entity_id
                            _row[CS_COMP_IDX] = _row[CS_SEQ_ID] = seq_id

                        if prefer_auth_atom_name:
                            if has_orig_seq:
                                orig_atom_id = _row[CS_ORIG_ATOM]
                            _atom_id = atom_id
                            _seq_key = seq_key if seq_key in coord_atom_site else (seq_key[0], seq_key[1])
                            if _seq_key in coord_atom_site:
                                _coord_atom_site = coord_atom_site[_seq_key]
                                if comp_id in auth_atom_name_to_id and comp_id == _coord_atom_site['comp_id']:
                                    if _atom_id in auth_atom_name_to_id[comp_id]:
                                        if auth_atom_name_to_id[comp_id][_atom_id] in _coord_atom_site['atom_id']:
                                            _row[CS_AUTH_ATOM] = atom_id = auth_atom_name_to_id[comp_id][_atom_id]
                                        elif 'split_comp_id' not in _coord_atom_site and has_orig_seq\
                                                and orig_atom_id in auth_atom_name_to_id[comp_id]:
                                            _row[CS_AUTH_ATOM] = atom_id = auth_atom_name_to_id[comp_id][orig_atom_id]
                                if 'alt_atom_id' in _coord_atom_site and _atom_id in _coord_atom_site['alt_atom_id']\
                                   and comp_id == _coord_atom_site['comp_id']:
                                    _row[CS_AUTH_ATOM] = atom_id =\
                                        _coord_atom_site['atom_id'][_coord_atom_site['alt_atom_id'].index(_atom_id)]
                                # DAOTHER-8751, 8817 (D_1300043061)
                                elif 'alt_comp_id' in _coord_atom_site and 'alt_atom_id' in _coord_atom_site\
                                     and _atom_id in _coord_atom_site['alt_atom_id']\
                                     and comp_id == _coord_atom_site['alt_comp_id'][_coord_atom_site['alt_atom_id'].index(_atom_id)]:  # noqa: E501, pylint: disable=line-too-long
                                    _row[CS_AUTH_COMP] = comp_id
                                    # Entity_assembly_ID, Entity_ID, Comp_index_ID, Seq_ID, Comp_ID, Auth_asym_ID, Auth_seq_ID
                                    cca_row =\
                                        next((cca_row for cca_row in self._reg.chem_comp_asm_dat
                                              if cca_row[CCA_COMP_ID] == comp_id and cca_row[CCA_SEQ_ID] == _seq_key[0]
                                              and cca_row[CCA_AUTH_ASYM] == _seq_key[1]), None)
                                    if cca_row is not None:
                                        _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID], _row[CS_COMP_IDX], _row[CS_SEQ_ID] =\
                                            cca_row[CCA_ENT_ASM_ID], cca_row[CCA_ENTITY_ID], \
                                            cca_row[CCA_COMP_IDX], cca_row[CCA_SEQ_ID]
                                    if comp_id in auth_atom_name_to_id_ext and _atom_id in auth_atom_name_to_id_ext[comp_id]\
                                       and len(set(_coord_atom_site['alt_comp_id'])) > 1:
                                        _row[CS_AUTH_ATOM] = atom_id = auth_atom_name_to_id_ext[comp_id][_atom_id]
                                    else:
                                        _row[CS_AUTH_ATOM] = atom_id =\
                                            _coord_atom_site['atom_id'][_coord_atom_site['alt_atom_id'].index(_atom_id)]
                                elif 'split_comp_id' in _coord_atom_site:
                                    for _comp_id in _coord_atom_site['split_comp_id']:
                                        if _comp_id == comp_id:
                                            continue
                                        __seq_key = (_seq_key[0], _seq_key[1], _comp_id)
                                        __coord_atom_site = coord_atom_site[__seq_key]
                                        if __coord_atom_site is None:
                                            continue
                                        if 'alt_comp_id' in __coord_atom_site and 'alt_atom_id' in __coord_atom_site\
                                           and _atom_id in __coord_atom_site['alt_atom_id']:
                                            comp_id = _comp_id
                                            _row[CS_AUTH_COMP] = comp_id
                                            # Entity_assembly_ID, Entity_ID, Comp_index_ID,
                                            # Seq_ID, Comp_ID, Auth_asym_ID, Auth_seq_ID
                                            cca_row =\
                                                next((cca_row for cca_row in self._reg.chem_comp_asm_dat
                                                      if cca_row[CCA_COMP_ID] == comp_id and cca_row[CCA_SEQ_ID] == _seq_key[0]
                                                      and cca_row[CCA_AUTH_ASYM] == _seq_key[1]), None)
                                            if cca_row is not None:
                                                _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID], _row[CS_COMP_IDX], _row[CS_SEQ_ID] =\
                                                    cca_row[CCA_ENT_ASM_ID], cca_row[CCA_ENTITY_ID], \
                                                    cca_row[CCA_COMP_IDX], cca_row[CCA_SEQ_ID]
                                            _row[CS_AUTH_ATOM] = atom_id =\
                                                __coord_atom_site['atom_id'][__coord_atom_site['alt_atom_id'].index(_atom_id)]
                                            _seq_key = __seq_key
                                            break
                                        if _atom_id in __coord_atom_site['atom_id']:
                                            comp_id = _comp_id
                                            _row[CS_AUTH_COMP] = comp_id
                                            # Entity_assembly_ID, Entity_ID, Comp_index_ID,
                                            # Seq_ID, Comp_ID, Auth_asym_ID, Auth_seq_ID
                                            cca_row =\
                                                next((cca_row for cca_row in self._reg.chem_comp_asm_dat
                                                      if cca_row[CCA_COMP_ID] == comp_id and cca_row[CCA_SEQ_ID] == _seq_key[0]
                                                      and cca_row[CCA_AUTH_ASYM] == _seq_key[1]), None)
                                            if cca_row is not None:
                                                _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID], _row[CS_COMP_IDX], _row[CS_SEQ_ID] =\
                                                    cca_row[CCA_ENT_ASM_ID], cca_row[CCA_ENTITY_ID], \
                                                    cca_row[CCA_COMP_IDX], cca_row[CCA_SEQ_ID]
                                            _seq_key = __seq_key
                                            break

                        apply_orig_seq_naming(std_res_branch=True, aux_cols=True)

                    else:
                        resolved = False

                    if not resolved:

                        chain_id = row[chain_id_col]
                        if chain_id in EMPTY_VALUE:
                            chain_id = REPRESENTATIVE_ASYM_ID

                        if chain_id in copied_chain_ids:
                            continue

                        try:
                            seq_id = int(row[seq_id_col])
                        except (ValueError, TypeError):
                            seq_id = None

                        if auth_asym_id_col != -1 and row[auth_asym_id_col] == 'UNMAPPED':
                            _row[CS_DETAILS] = 'UNMAPPED'

                        auth_asym_id, auth_seq_id = get_auth_seq_scheme(chain_id, seq_id)

                        resolved = True

                        if None not in (auth_asym_id, auth_seq_id):
                            seq_key = (auth_asym_id, auth_seq_id, _orig_comp_id)
                            _seq_key = (seq_key[0], seq_key[1])
                            if seq_key in auth_to_star_seq:
                                entity_assembly_id, seq_id, entity_id, _ = auth_to_star_seq[seq_key]
                                comp_id = next((_v[1] for _k, _v in auth_to_orig_seq.items() if _k == seq_key), _orig_comp_id)
                                self._reg.ent_asym_id_with_exptl_data.add(entity_assembly_id)
                                _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID] = entity_assembly_id, entity_id
                                _row[CS_COMP_IDX] = _row[CS_SEQ_ID] = seq_id

                                _row[CS_AUTH_ASYM], _row[CS_AUTH_SEQ], _row[CS_AUTH_COMP], _row[CS_AUTH_ATOM] =\
                                    auth_asym_id, auth_seq_id, comp_id, atom_id
                                apply_orig_seq_naming(guard_has_auth_seq=True, record_can_mapping=True)

                            else:

                                item = next((item for item in entity_assembly if item['auth_asym_id'] == auth_asym_id), None)

                                if item is not None and poly_seq is not None\
                                   and any(True for _ps in poly_seq
                                           if _ps['chain_id'] == auth_asym_id and auth_seq_id in _ps['seq_id']):
                                    entity_assembly_id = item['entity_assembly_id']
                                    entity_id = item['entity_id']

                                    _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID] = entity_assembly_id, entity_id
                                    _row[CS_COMP_IDX] = _row[CS_SEQ_ID] = seq_id

                                    seq_key = next((k for k, v in auth_to_star_seq.items()
                                                    if v[0] == entity_assembly_id and v[1] == seq_id and v[2] == entity_id), None)
                                    if seq_key is not None and (seq_id == seq_key[1] or comp_id == seq_key[2]):
                                        _seq_key = (seq_key[0], seq_key[1])
                                        _row[CS_AUTH_ASYM], _row[CS_AUTH_SEQ], _row[CS_AUTH_COMP], _row[CS_AUTH_ATOM] =\
                                            seq_key[0], seq_key[1], seq_key[2], atom_id
                                        if has_ins_code and seq_key in auth_to_ins_code:
                                            _row[CS_INS_CODE] = auth_to_ins_code[seq_key]
                                    elif seq_key is not None:  # 5ydx
                                        offset = seq_id - seq_key[1]
                                        seq_id += offset
                                        _row[CS_COMP_IDX] = _row[CS_SEQ_ID] = seq_id
                                        _seq_key = None
                                        _row[CS_DETAILS] = 'UNMAPPED'
                                    else:
                                        resolved = False

                                    if has_auth_seq:
                                        _row[CS_ORIG_ASYM], _row[CS_ORIG_SEQ], _row[CS_ORIG_COMP], _row[CS_ORIG_ATOM] =\
                                            row[auth_asym_id_col], row[auth_seq_id_col], \
                                            row[auth_comp_id_col], row[auth_atom_id_col]

                                    if resolved:
                                        index, _row, reparse = fill_cs_row(lp, index, _row, prefer_auth_atom_name,
                                                                           coord_atom_site, _seq_key,
                                                                           comp_id, atom_id, loop, idx)
                                        reparse_request |= reparse

                                    if chain_id not in can_auth_asym_id_mapping:
                                        can_auth_asym_id_mapping[chain_id] = {'auth_asym_id': auth_asym_id,
                                                                              'ref_auth_seq_id': auth_seq_id
                                                                              }

                                else:
                                    resolved = False

                        else:

                            if has_auth_seq:

                                try:

                                    auth_asym_id = row[auth_asym_id_col]
                                    auth_seq_id = int(row[auth_seq_id_col])

                                    item = next((item for item in entity_assembly if item['auth_asym_id'] == auth_asym_id), None)

                                    if item is not None and poly_seq is not None\
                                       and any(True for _ps in poly_seq
                                               if _ps['chain_id'] == auth_asym_id and auth_seq_id in _ps['seq_id']):
                                        entity_assembly_id = item['entity_assembly_id']
                                        entity_id = item['entity_id']

                                        _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID] = entity_assembly_id, entity_id
                                        _row[CS_COMP_IDX] = _row[CS_SEQ_ID] = seq_id

                                        seq_key = next((k for k, v in auth_to_star_seq.items()
                                                        if v[0] == entity_assembly_id and v[1] == seq_id
                                                        and v[2] == entity_id), None)
                                        if seq_key is not None:

                                            if comp_id != seq_key[2] and comp_id in STD_MON_DICT and seq_key[2] in STD_MON_DICT:
                                                resolved = False

                                            else:
                                                _seq_key = (seq_key[0], seq_key[1])
                                                _row[CS_AUTH_ASYM], _row[CS_AUTH_SEQ], _row[CS_AUTH_COMP], _row[CS_AUTH_ATOM] =\
                                                    seq_key[0], seq_key[1], seq_key[2], atom_id
                                                if has_ins_code and seq_key in auth_to_ins_code:
                                                    _row[CS_INS_CODE] = auth_to_ins_code[seq_key]

                                        if resolved:
                                            _row[CS_ORIG_ASYM], _row[CS_ORIG_SEQ], _row[CS_ORIG_COMP], _row[CS_ORIG_ATOM] =\
                                                row[auth_asym_id_col], row[auth_seq_id_col], \
                                                row[auth_comp_id_col], row[auth_atom_id_col]

                                            index, _row, reparse = fill_cs_row(lp, index, _row, prefer_auth_atom_name,
                                                                               coord_atom_site, _seq_key,
                                                                               comp_id, atom_id, loop, idx)
                                            reparse_request |= reparse

                                    else:
                                        resolved = False

                                except (ValueError, TypeError):
                                    resolved = False

                            else:

                                def retrieve_label_comp_id(can_seq_id, can_comp_id):
                                    for _seq_key in auth_to_star_seq.keys():
                                        if _seq_key[1] == can_seq_id and _seq_key in auth_to_orig_seq:
                                            if can_comp_id == auth_to_orig_seq[_seq_key][1]:
                                                return _seq_key[2]
                                    return can_comp_id

                                can_auth_asym_id = [_auth_asym_id for _auth_asym_id, _auth_seq_id, _comp_id in auth_to_star_seq
                                                    if _auth_seq_id == seq_id and _comp_id in (_orig_comp_id, '.')]

                                if len(can_auth_asym_id) == 0:
                                    can_auth_asym_id = [_auth_asym_id for _auth_asym_id, _auth_seq_id, _comp_id in auth_to_star_seq
                                                        if _auth_seq_id == seq_id
                                                        and _comp_id == retrieve_label_comp_id(seq_id, _orig_comp_id)]
                                    if len(can_auth_asym_id) > 0:
                                        _orig_comp_id = retrieve_label_comp_id(seq_id, _orig_comp_id)

                                if len(can_auth_asym_id) != 1:
                                    resolved = False

                                else:
                                    auth_asym_id, auth_seq_id = can_auth_asym_id[0], seq_id

                                    seq_key = (auth_asym_id, auth_seq_id, _orig_comp_id)
                                    dummy_key = (auth_asym_id, auth_seq_id, '.')
                                    _seq_key = (seq_key[0], seq_key[1])
                                    if seq_key in auth_to_star_seq\
                                       or (dummy_key in auth_to_star_seq and _orig_comp_id in STD_MON_DICT):
                                        try:
                                            entity_assembly_id, seq_id, entity_id, _ = auth_to_star_seq[seq_key]
                                        except KeyError:  # DAOTHER-9644: map residue on truncated loop
                                            if _seq_key not in truncated_loop_sequence:
                                                truncated_loop_sequence.append(_seq_key)
                                            entity_assembly_id, seq_id, entity_id, _ = auth_to_star_seq[dummy_key]
                                            auth_to_star_seq[seq_key] = auth_to_star_seq[dummy_key]
                                            del auth_to_star_seq[dummy_key]
                                            auth_to_orig_seq[seq_key] = (auth_to_orig_seq[dummy_key][0], _orig_comp_id)
                                            del auth_to_orig_seq[dummy_key]
                                            if has_ins_code:
                                                auth_to_ins_code[seq_key] = auth_to_ins_code[dummy_key]
                                                del auth_to_ins_code[dummy_key]
                                            cif_ps = next((cif_ps for cif_ps in self._reg.caC['polymer_sequence']
                                                           if cif_ps['auth_chain_id'] == auth_asym_id
                                                           and auth_seq_id in cif_ps['auth_seq_id']), None)
                                            if cif_ps is not None:
                                                _idx_ = cif_ps['auth_seq_id'].index(auth_seq_id)
                                                if cif_ps['comp_id'][_idx_] in EMPTY_VALUE:
                                                    cif_ps['comp_id'][_idx_] = cif_ps['auth_comp_id'][_idx_] = _orig_comp_id
                                                    if self._reg.asmChkCachePath is not None:
                                                        write_as_pickle(self._reg.caC, self._reg.asmChkCachePath)
                                        comp_id = next((_v[1] for _k, _v in auth_to_orig_seq.items()
                                                        if _k == seq_key), _orig_comp_id)
                                        self._reg.ent_asym_id_with_exptl_data.add(entity_assembly_id)
                                        _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID] = entity_assembly_id, entity_id
                                        _row[CS_COMP_IDX] = _row[CS_SEQ_ID] = seq_id

                                        _row[CS_AUTH_ASYM], _row[CS_AUTH_SEQ], _row[CS_AUTH_COMP], _row[CS_AUTH_ATOM] =\
                                            auth_asym_id, auth_seq_id, comp_id, atom_id
                                        apply_orig_seq_naming(guard_has_auth_seq=True, record_can_mapping=True)

                                    else:

                                        item = next((item for item in entity_assembly
                                                     if item['auth_asym_id'] == auth_asym_id), None)

                                        if item is not None and poly_seq is not None\
                                           and any(True for _ps in poly_seq
                                                   if _ps['chain_id'] == auth_asym_id and auth_seq_id in _ps['seq_id']):
                                            entity_assembly_id = item['entity_assembly_id']
                                            entity_id = item['entity_id']

                                            _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID] = entity_assembly_id, entity_id
                                            _row[CS_COMP_IDX] = _row[CS_SEQ_ID] = seq_id

                                            seq_key = next((k for k, v in auth_to_star_seq.items()
                                                            if v[0] == entity_assembly_id and v[1] == seq_id
                                                            and v[2] == entity_id), None)
                                            if seq_key is not None:
                                                _seq_key = (seq_key[0], seq_key[1])
                                                _row[CS_AUTH_ASYM], _row[CS_AUTH_SEQ], _row[CS_AUTH_COMP], _row[CS_AUTH_ATOM] =\
                                                    seq_key[0], seq_key[1], seq_key[2], atom_id
                                                if has_ins_code and seq_key in auth_to_ins_code:
                                                    _row[CS_INS_CODE] = auth_to_ins_code[seq_key]

                                            if has_auth_seq:
                                                _row[CS_ORIG_ASYM], _row[CS_ORIG_SEQ], _row[CS_ORIG_COMP], _row[CS_ORIG_ATOM] =\
                                                    row[auth_asym_id_col], row[auth_seq_id_col], \
                                                    row[auth_comp_id_col], row[auth_atom_id_col]

                                            index, _row, reparse = fill_cs_row(lp, index, _row, prefer_auth_atom_name,
                                                                               coord_atom_site, _seq_key,
                                                                               comp_id, atom_id, loop, idx)
                                            reparse_request |= reparse

                                            if chain_id not in can_auth_asym_id_mapping:
                                                can_auth_asym_id_mapping[chain_id] = {'auth_asym_id': auth_asym_id,
                                                                                      'ref_auth_seq_id': seq_key[1]
                                                                                      }

                                        else:
                                            resolved = False

                        is_valid, cc_name, _ = self._reg.dpV.getChemCompNameAndStatusOf(comp_id)
                        comp_id_bmrb_only = not is_valid and cc_name is not None and 'processing site' in cc_name

                        if not resolved and has_auth_seq and not comp_id_bmrb_only:
                            try:
                                seq_id = int(row[auth_seq_id_col])
                            except (ValueError, TypeError):
                                seq_id = None

                        if not resolved and seq_id is not None and self._reg.cifChecked:

                            def test_seq_id_offset(lp, index, row, _row, _idx, chain_id, seq_id, comp_id, offset):
                                _found = _resolved = _reparse = False
                                _index = index

                                auth_asym_id, auth_seq_id = get_auth_seq_scheme(chain_id, seq_id + offset)
                                if None not in (auth_asym_id, auth_seq_id):
                                    _found = _resolved = True

                                    item = next((item for item in entity_assembly if item['auth_asym_id'] == auth_asym_id), None)

                                    if item is not None and poly_seq_common is not None\
                                       and any(True for _ps in poly_seq_common
                                               if _ps['chain_id'] == auth_asym_id and auth_seq_id in _ps['seq_id']):
                                        entity_assembly_id = item['entity_assembly_id']
                                        entity_id = item['entity_id']

                                        _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID] = entity_assembly_id, entity_id
                                        _row[CS_COMP_IDX] = _row[CS_SEQ_ID] = seq_id

                                        seq_key = next((k for k, v in auth_to_star_seq.items()
                                                        if v[0] == entity_assembly_id and v[1] == seq_id + offset
                                                        and v[2] == entity_id), None)
                                        _seq_key = None
                                        if seq_key is not None and comp_id == seq_key[2]:
                                            _seq_key = (seq_key[0], seq_key[1])
                                            _row[CS_AUTH_ASYM], _row[CS_AUTH_SEQ], _row[CS_AUTH_COMP], _row[CS_AUTH_ATOM] =\
                                                seq_key[0], seq_key[1] - offset, comp_id, atom_id
                                            if has_ins_code and seq_key in auth_to_ins_code:
                                                _row[CS_INS_CODE] = auth_to_ins_code[seq_key]
                                        else:
                                            if has_orig_seq:  # DAOTHER-8758
                                                try:
                                                    orig_asym_id = row[orig_asym_id_col]
                                                    orig_seq_id = int(row[orig_seq_id_col])
                                                    _item = next((item for item in entity_assembly
                                                                  if item['auth_asym_id'] == orig_asym_id), None)
                                                    if _item is not None:
                                                        _entity_assembly_id = _item['entity_assembly_id']
                                                        _entity_id = _item['entity_id']
                                                        __seq_key = next((k for k, v in auth_to_star_seq.items()
                                                                          if v[0] == _entity_assembly_id
                                                                          and v[1] in (seq_id, orig_seq_id)
                                                                          and v[2] == _entity_id), None)
                                                        if __seq_key is not None:
                                                            comp_id = __seq_key[2]
                                                            _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID], \
                                                                _row[CS_COMP_IDX], _row[CS_SEQ_ID] =\
                                                                _entity_assembly_id, _entity_id, __seq_key[1], __seq_key[1]
                                                            _seq_key = (__seq_key[0], __seq_key[1])
                                                            _row[CS_AUTH_ASYM], _row[CS_AUTH_SEQ], \
                                                                _row[CS_AUTH_COMP], _row[CS_AUTH_ATOM] =\
                                                                __seq_key[0], __seq_key[1], comp_id, atom_id
                                                            if has_ins_code and __seq_key in auth_to_ins_code:
                                                                _row[CS_INS_CODE] = auth_to_ins_code[__seq_key]
                                                        else:
                                                            _seq_key = (auth_asym_id, auth_seq_id + offset)
                                                    else:
                                                        _seq_key = (auth_asym_id, auth_seq_id + offset)
                                                except ValueError:
                                                    _seq_key = (auth_asym_id, auth_seq_id + offset)
                                            else:
                                                _item = next((item for item in entity_assembly
                                                              if item['auth_asym_id'] == chain_id), None)
                                                if _item is not None:
                                                    _entity_assembly_id = _item['entity_assembly_id']
                                                    _entity_id = _item['entity_id']
                                                    __seq_key = next((k for k, v in auth_to_star_seq.items()
                                                                      if v[0] == _entity_assembly_id and v[1] == seq_id + offset
                                                                      and v[2] == _entity_id), None)
                                                    if __seq_key is not None:
                                                        _offset = __seq_key[1] - (seq_id + offset)
                                                        _seq_id = seq_id - _offset
                                                        __seq_key = next((k for k, v in auth_to_star_seq.items()
                                                                          if v[0] == _entity_assembly_id and v[1] == _seq_id
                                                                          and v[2] == _entity_id), None)
                                                        if __seq_key is not None and comp_id == __seq_key[2]:
                                                            comp_id = __seq_key[2]
                                                            _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID], \
                                                                _row[CS_COMP_IDX], _row[CS_SEQ_ID] =\
                                                                _entity_assembly_id, _entity_id, _seq_id, _seq_id
                                                            _seq_key = (__seq_key[0], __seq_key[1])
                                                            _row[CS_AUTH_ASYM], _row[CS_AUTH_SEQ], \
                                                                _row[CS_AUTH_COMP], _row[CS_AUTH_ATOM] =\
                                                                __seq_key[0], __seq_key[1], comp_id, atom_id
                                                            if has_ins_code and __seq_key in auth_to_ins_code:
                                                                _row[CS_INS_CODE] = auth_to_ins_code[__seq_key]
                                                        else:
                                                            _resolved = False
                                                    else:
                                                        _resolved = False
                                                else:
                                                    _resolved = False

                                        if has_auth_seq:
                                            _row[CS_ORIG_ASYM], _row[CS_ORIG_SEQ], _row[CS_ORIG_COMP], _row[CS_ORIG_ATOM] =\
                                                row[auth_asym_id_col], row[auth_seq_id_col], \
                                                row[auth_comp_id_col], row[auth_atom_id_col]
                                        else:
                                            _row[CS_ORIG_ASYM], _row[CS_ORIG_SEQ], _row[CS_ORIG_COMP], _row[CS_ORIG_ATOM] =\
                                                _row[CS_AUTH_ASYM], _row[CS_AUTH_SEQ], _row[CS_AUTH_COMP], _row[CS_AUTH_ATOM]

                                        if _resolved:
                                            _index, _row, _reparse = fill_cs_row(lp, index, _row, prefer_auth_atom_name,
                                                                                 coord_atom_site, _seq_key,
                                                                                 comp_id, atom_id, loop, _idx)

                                    else:
                                        _resolved = False

                                return _found, _resolved, _reparse, _index, _row

                            found = False
                            for offset in range(1, GLOBAL_OFFSET_ATTEMPT):
                                found, resolved, reparse, _index, __row =\
                                    test_seq_id_offset(lp, index, row, _row, idx, chain_id, seq_id, comp_id, offset)

                                if found:
                                    if resolved:
                                        index, _row = _index, __row
                                        reparse_request |= reparse
                                    break

                                found, resolved, reparse, _index, __row =\
                                    test_seq_id_offset(lp, index, row, _row, idx, chain_id, seq_id, comp_id, -offset)

                                if found:
                                    if resolved:
                                        index, _row = _index, __row
                                        reparse_request |= reparse
                                    break

                            if not resolved and chain_id in can_auth_asym_id_mapping:  # DAOTHER-8751, 8755

                                if can_auth_asym_id_mapping_failed and trial == 0:  # DAOTHER-9158
                                    reparse_request = True

                                mapping = can_auth_asym_id_mapping[chain_id]

                                auth_asym_id = mapping['auth_asym_id']
                                ref_auth_seq_id = mapping['ref_auth_seq_id']

                                item = next((item for item in entity_assembly if item['auth_asym_id'] == auth_asym_id), None)

                                if item is not None and poly_seq_common is not None\
                                    and any(True for _ps in poly_seq_common
                                            if _ps['chain_id'] in (auth_asym_id, str(letterToDigit(auth_asym_id)))
                                            and ref_auth_seq_id in _ps['seq_id']):
                                    resolved = True
                                    found = False

                                    entity_assembly_id = item['entity_assembly_id']
                                    entity_id = item['entity_id']

                                    _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID] = entity_assembly_id, entity_id
                                    _row[CS_COMP_IDX] = _row[CS_SEQ_ID] = seq_id

                                    _row[CS_AUTH_ASYM], _row[CS_AUTH_SEQ], _row[CS_AUTH_COMP], _row[CS_AUTH_ATOM] =\
                                        auth_asym_id, seq_id, comp_id, atom_id

                                    if has_auth_seq:
                                        _row[CS_ORIG_ASYM], _row[CS_ORIG_SEQ], _row[CS_ORIG_COMP], _row[CS_ORIG_ATOM] =\
                                            row[auth_asym_id_col], row[auth_seq_id_col], \
                                            row[auth_comp_id_col], row[auth_atom_id_col]
                                    else:
                                        _row[CS_ORIG_ASYM], _row[CS_ORIG_SEQ], _row[CS_ORIG_COMP], _row[CS_ORIG_ATOM] =\
                                            _row[CS_AUTH_ASYM], _row[CS_AUTH_SEQ], _row[CS_AUTH_COMP], _row[CS_AUTH_ATOM]

                                    # DAOTHER-9281
                                    if isinstance(_row[CS_ENT_ASM_ID], int)\
                                       and str(_row[CS_ENT_ASM_ID]) in seq_id_offset_for_unmapped:
                                        __offset = seq_id_offset_for_unmapped[str(_row[CS_ENT_ASM_ID])]
                                    elif isinstance(_row[CS_ENT_ASM_ID], str)\
                                            and _row[CS_ENT_ASM_ID] in seq_id_offset_for_unmapped:
                                        __offset = seq_id_offset_for_unmapped[_row[CS_ENT_ASM_ID]]
                                    else:
                                        __offset = 0

                                    if comp_id not in STD_MON_DICT:
                                        for item in entity_assembly:
                                            if 'comp_id' in item and comp_id == item['comp_id']:
                                                _entity_assembly_id = item['entity_assembly_id']
                                                _entity_id = item['entity_id']

                                                __seq_key = next((k for k, v in auth_to_star_seq.items()
                                                                  if v[0] == _entity_assembly_id and v[1] == seq_id
                                                                  and v[2] == _entity_id), None)
                                                if __seq_key is not None:
                                                    found = True
                                                    comp_id = __seq_key[2]
                                                    _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID], _row[CS_COMP_IDX], _row[CS_SEQ_ID] =\
                                                        _entity_assembly_id, _entity_id, __seq_key[1], __seq_key[1]
                                                    _seq_key = (__seq_key[0], __seq_key[1])
                                                    _row[CS_AUTH_ASYM], _row[CS_AUTH_SEQ], _row[CS_AUTH_COMP], _row[CS_AUTH_ATOM] =\
                                                        __seq_key[0], __seq_key[1], comp_id, atom_id
                                                    if has_ins_code and __seq_key in auth_to_ins_code:
                                                        _row[CS_INS_CODE] = auth_to_ins_code[__seq_key]
                                                    break

                                                if self._reg.caC['non_polymer'] is not None:
                                                    ligands = 0
                                                    for np in self._reg.caC['non_polymer']:
                                                        if comp_id == np['comp_id'][0]:
                                                            ligands += len(np['seq_id'])
                                                    if ligands == 1:  # DAOTHER-9063, 2nd case
                                                        __seq_key = next((k for k, v in auth_to_star_seq.items()
                                                                          if v[0] == _entity_assembly_id
                                                                          and v[2] == _entity_id), None)
                                                        if __seq_key is not None:
                                                            seq_id = auth_to_star_seq[__seq_key][1]
                                                            found = True
                                                            comp_id = __seq_key[2]
                                                            _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID] =\
                                                                _entity_assembly_id, _entity_id
                                                            _row[CS_COMP_IDX] = _row[CS_SEQ_ID] = seq_id
                                                            _seq_key = (__seq_key[0], __seq_key[1])
                                                            _row[CS_AUTH_ASYM], _row[CS_AUTH_SEQ], \
                                                                _row[CS_AUTH_COMP], _row[CS_AUTH_ATOM] =\
                                                                __seq_key[0], __seq_key[1], comp_id, atom_id
                                                            if has_ins_code and __seq_key in auth_to_ins_code:
                                                                _row[CS_INS_CODE] = auth_to_ins_code[__seq_key]
                                                            break

                                    else:

                                        __seq_key = next((k for k, v in auth_to_star_seq.items()
                                                          if v[0] == entity_assembly_id
                                                          and v[1] == seq_id + __offset
                                                          and v[2] == entity_id), None)
                                        if __seq_key is not None:
                                            __comp_id = __seq_key[2]
                                            if self._reg.ccU.updateChemCompDict(comp_id):
                                                cc_type = self._reg.ccU.lastChemCompDict['type']
                                                if self._reg.ccU.updateChemCompDict(__comp_id):
                                                    __cc_type = self._reg.ccU.lastChemCompDict['type']
                                                    if cc_type == __cc_type:  # DAOTHER-9198
                                                        found = True
                                                        comp_id = __seq_key[2]
                                                        _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID], \
                                                            _row[CS_COMP_IDX], _row[CS_SEQ_ID] =\
                                                            entity_assembly_id, entity_id, __seq_key[1], __seq_key[1]
                                                        _seq_key = (__seq_key[0], __seq_key[1])
                                                        _row[CS_AUTH_ASYM], _row[CS_AUTH_SEQ], \
                                                            _row[CS_AUTH_COMP], _row[CS_AUTH_ATOM] =\
                                                            __seq_key[0], __seq_key[1], comp_id, atom_id
                                                        if has_ins_code and __seq_key in auth_to_ins_code:
                                                            _row[CS_INS_CODE] = auth_to_ins_code[__seq_key]

                                    if not found:
                                        _row[CS_DETAILS] = 'UNMAPPED'
                                        # DAOTHER-9065
                                        if __offset != 0:
                                            _row[CS_COMP_IDX] += __offset
                                            _row[CS_SEQ_ID] = _row[CS_COMP_IDX]
                                        elif trial == 0:
                                            reparse_request = True

                                        _seq_key = (auth_asym_id, seq_id)

                                    _index, _row, reparse = fill_cs_row(lp, index, _row, prefer_auth_atom_name,
                                                                        coord_atom_site, _seq_key,
                                                                        comp_id, atom_id, loop, idx)
                                    reparse_request |= reparse

                            if not resolved and seq_id is not None and self._reg.cifChecked:

                                can_auth_asym_id_mapping_failed = True  # DAOTHER-9158

                                def test_seq_id_offset_as_is(lp, index, _row, _idx, chain_id, seq_id, comp_id, offset):
                                    _resolved = _reparse = False
                                    _index, _seq_id = index, seq_id

                                    auth_asym_id, auth_seq_id, label_seq_id = get_label_seq_scheme(chain_id, seq_id + offset)
                                    if None not in (auth_asym_id, auth_seq_id):
                                        _resolved = True

                                        item = next((item for item in entity_assembly
                                                     if item['auth_asym_id'] == auth_asym_id), None)

                                        if item is not None and poly_seq_common is not None\
                                           and any(True for _ps in poly_seq_common
                                                   if _ps['chain_id'] == chain_id and label_seq_id in _ps['seq_id']):
                                            entity_assembly_id = item['entity_assembly_id']
                                            entity_id = item['entity_id']

                                            seq_key = next((k for k, v in auth_to_star_seq.items()
                                                            if k[0] == auth_asym_id and k[1] == auth_seq_id
                                                            and v[0] == entity_assembly_id and v[2] == entity_id), None)

                                            if seq_key is not None:
                                                _, _label_seq_id, _, _ = auth_to_star_seq[seq_key]

                                                if entity_id not in label_seq_id_offset_for_extended\
                                                   or _label_seq_id - label_seq_id == label_seq_id_offset_for_extended[entity_id]:
                                                    seq_id += (label_seq_id - auth_seq_id)
                                                    seq_id += (_label_seq_id - label_seq_id)

                                                    if entity_id not in label_seq_id_offset_for_extended:
                                                        label_seq_id_offset_for_extended[entity_id] = _label_seq_id - label_seq_id

                                                    _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID] = entity_assembly_id, entity_id
                                                    _row[CS_COMP_IDX] = _row[CS_SEQ_ID] = seq_id

                                                    if _row[CS_AUTH_SEQ] in EMPTY_VALUE:
                                                        _row[CS_AUTH_SEQ] = _seq_id
                                                    if _row[CS_AUTH_COMP] in EMPTY_VALUE:
                                                        _row[CS_AUTH_COMP] = comp_id
                                                    if _row[CS_AUTH_ATOM] in EMPTY_VALUE:
                                                        _row[CS_AUTH_ATOM] = atom_id

                                                    _row[CS_AUTH_ASYM] = _row[CS_ORIG_ASYM] = auth_asym_id
                                                    if _row[CS_ORIG_SEQ] in EMPTY_VALUE:
                                                        _row[CS_ORIG_SEQ] = _row[CS_AUTH_SEQ]
                                                    if _row[CS_ORIG_COMP] in EMPTY_VALUE:
                                                        _row[CS_ORIG_COMP] = _row[CS_AUTH_COMP]
                                                    if _row[CS_ORIG_ATOM] in EMPTY_VALUE:
                                                        _row[CS_ORIG_ATOM] = _row[CS_AUTH_ATOM]
                                                    if _row[CS_DETAILS] in EMPTY_VALUE:
                                                        _row[CS_DETAILS] = 'UNMAPPED'

                                                    _index, _row, reparse = fill_cs_row(lp, index, _row, prefer_auth_atom_name,
                                                                                        coord_atom_site, None,
                                                                                        comp_id, atom_id, loop, _idx)
                                                    _reparse |= reparse

                                                else:
                                                    _resolved = False

                                            else:
                                                _resolved = False

                                        else:
                                            _resolved = False

                                    return _resolved, _reparse, _index, _row

                                found = False
                                for offset in range(1, GLOBAL_OFFSET_ATTEMPT):
                                    resolved, reparse, _index, __row =\
                                        test_seq_id_offset_as_is(lp, index, _row, idx, chain_id, seq_id, comp_id, offset)

                                    if resolved:
                                        index, _row = _index, __row
                                        reparse_request |= reparse
                                        break

                                    resolved, reparse, _index, __row =\
                                        test_seq_id_offset_as_is(lp, index, _row, idx, chain_id, seq_id, comp_id, -offset)

                                    if resolved:
                                        index, _row = _index, __row
                                        reparse_request |= reparse
                                        break

                        if not resolved:

                            entity_id = None
                            if (self._reg.combined_mode or (self._reg.bmrb_only and self._reg.internal_mode))\
                               and entity_id_col != -1:
                                try:
                                    entity_id = int(row[entity_id_col])
                                except (ValueError, TypeError):
                                    entity_id = None

                            if not self._reg.cifChecked:
                                seq_id = int(row[seq_id_col])

                            _row[CS_ENT_ASM_ID], _row[CS_ENTITY_ID], _row[CS_COMP_ID] = chain_id, entity_id, comp_id
                            _row[CS_COMP_IDX] = _row[CS_SEQ_ID] = seq_id

                            # DAOTHER-9065
                            if details_col != -1 and row[details_col] == 'UNMAPPED':
                                if isinstance(_row[CS_ENT_ASM_ID], int)\
                                   and str(_row[CS_ENT_ASM_ID]) in seq_id_offset_for_unmapped:
                                    __offset = seq_id_offset_for_unmapped[str(_row[CS_ENT_ASM_ID])]
                                elif isinstance(_row[CS_ENT_ASM_ID], str)\
                                        and _row[CS_ENT_ASM_ID] in seq_id_offset_for_unmapped:
                                    __offset = seq_id_offset_for_unmapped[_row[CS_ENT_ASM_ID]]
                                else:
                                    __offset = None

                                if __offset is not None:
                                    offset = None
                                    if isinstance(_row[CS_AUTH_SEQ], int):
                                        offset = _row[CS_COMP_IDX] - _row[CS_AUTH_SEQ]
                                    elif isinstance(_row[CS_AUTH_SEQ], str) and _row[CS_AUTH_SEQ].isdigit():
                                        offset = _row[CS_COMP_IDX] - int(_row[CS_AUTH_SEQ])
                                    if offset is not None and offset != __offset:
                                        if isinstance(_row[CS_AUTH_SEQ], int):
                                            _row[CS_COMP_IDX] = _row[CS_AUTH_SEQ] + __offset
                                            _row[CS_SEQ_ID] = _row[CS_COMP_IDX]
                                        else:
                                            _row[CS_COMP_IDX] = int(_row[CS_AUTH_SEQ]) + __offset
                                            _row[CS_SEQ_ID] = _row[CS_COMP_IDX]
                                elif trial == 0:
                                    reparse_request = True

                            atom_ids = self._reg.dpV.getAtomIdListInXplor(comp_id, atom_id)
                            if len(atom_ids) == 0 or atom_ids[0] not in self._reg.csStat.getAllAtoms(comp_id):
                                atom_ids = self._reg.dpV.getAtomIdListInXplor(comp_id,
                                                                              translateToStdAtomName(atom_id, comp_id,
                                                                                                     ccU=self._reg.ccU))
                            len_atom_ids = len(atom_ids)
                            if len_atom_ids == 0 or comp_id_bmrb_only or _row[CS_DETAILS] == 'UNMAPPED':
                                _row[CS_ATOM_ID] = atom_id
                                _row[CS_ATOM_TYPE] = 'H' if atom_id[0] in PSE_PRO_BEGIN_CODE else atom_id[0]
                                if _row[CS_ATOM_TYPE] in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS:
                                    _row[CS_ISOTOPE] = ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS[_row[CS_ATOM_TYPE]][0]
                            else:
                                _row[CS_ATOM_ID] = atom_ids[0]
                                _row[CS_AUTH_ATOM] = None
                                fill_auth_atom_id = _row[CS_AUTH_COMP] not in EMPTY_VALUE
                                if self._reg.ccU.updateChemCompDict(comp_id):
                                    cca = next((cca for cca in self._reg.ccU.lastAtomDictList
                                                if cca['atom_id'] == _row[CS_ATOM_ID]), None)
                                    if cca is not None:
                                        _row[CS_ATOM_TYPE] = cca['type_symbol']
                                        if _row[CS_ATOM_TYPE] in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS:
                                            _row[CS_ISOTOPE] = ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS[_row[CS_ATOM_TYPE]][0]
                                    else:
                                        _row[CS_ATOM_TYPE] = 'H' if _row[CS_ATOM_ID][0] in PROTON_BEGIN_CODE else atom_id[0]
                                        if _row[CS_ATOM_TYPE] in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS:
                                            _row[CS_ISOTOPE] = ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS[_row[CS_ATOM_TYPE]][0]
                                else:
                                    _row[CS_ATOM_TYPE] = 'H' if atom_id[0] in PSE_PRO_BEGIN_CODE else atom_id[0]
                                    if _row[CS_ATOM_TYPE] in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS:
                                        _row[CS_ISOTOPE] = ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS[_row[CS_ATOM_TYPE]][0]

                                if len_atom_ids > 1:
                                    __row = copy.copy(_row)
                                    lp.add_data(__row)

                                    for _atom_id in atom_ids[1:]:
                                        __row = copy.copy(_row)

                                        index += 1

                                        __row[CS_ID] = index
                                        __row[CS_ATOM_ID] = _atom_id

                                        lp.add_data(__row)

                                    index += 1

                                    _row[CS_ATOM_ID] = atom_ids[-1]

                                if fill_auth_atom_id:
                                    _row[CS_AUTH_ATOM] =\
                                        _row[CS_ATOM_ID] if self._reg.caC is not None or row[auth_atom_id_col] in EMPTY_VALUE\
                                        else row[auth_atom_id_col]

                    # DAOTHER-9065
                    if isinstance(_row[CS_ENT_ASM_ID], int) and str(_row[CS_ENT_ASM_ID]) not in seq_id_offset_for_unmapped\
                       and _row[CS_DETAILS] in EMPTY_VALUE:
                        if isinstance(_row[CS_COMP_IDX], int):
                            if isinstance(_row[CS_AUTH_SEQ], int):
                                seq_id_offset_for_unmapped[str(_row[CS_ENT_ASM_ID])] = _row[CS_COMP_IDX] - _row[CS_AUTH_SEQ]
                            elif isinstance(_row[CS_AUTH_SEQ], str) and _row[CS_AUTH_SEQ].isdigit():
                                seq_id_offset_for_unmapped[str(_row[CS_ENT_ASM_ID])] = _row[CS_COMP_IDX] - int(_row[CS_AUTH_SEQ])
                    elif isinstance(_row[CS_ENT_ASM_ID], str) and _row[CS_ENT_ASM_ID] not in seq_id_offset_for_unmapped\
                            and _row[CS_DETAILS] in EMPTY_VALUE:
                        if isinstance(_row[CS_COMP_IDX], int):
                            if isinstance(_row[CS_AUTH_SEQ], int):
                                seq_id_offset_for_unmapped[_row[CS_ENT_ASM_ID]] = _row[CS_COMP_IDX] - _row[CS_AUTH_SEQ]
                            elif isinstance(_row[CS_AUTH_SEQ], str) and _row[CS_AUTH_SEQ].isdigit():
                                seq_id_offset_for_unmapped[_row[CS_ENT_ASM_ID]] = _row[CS_COMP_IDX] - int(_row[CS_AUTH_SEQ])

                    if self._reg.bmrb_only and self._reg.internal_mode and _row[CS_ENT_ASM_ID] not in EMPTY_VALUE\
                       and _row[CS_ENTITY_ID] in EMPTY_VALUE:
                        entity_assembly_id =\
                            _row[CS_ENT_ASM_ID] if isinstance(_row[CS_ENT_ASM_ID], str) else str(_row[CS_ENT_ASM_ID])
                        if entity_assembly_id in entity_assembly_mappping:
                            _row[CS_ENTITY_ID] = entity_assembly_mappping[entity_assembly_id]

                    if isinstance(_row[CS_AMBIG_CODE], int):
                        comp_id = _row[CS_COMP_ID]
                        atom_id = _row[CS_ATOM_ID]
                        ambig_code = _row[CS_AMBIG_CODE]

                        if ambig_code == 0:
                            _row[CS_AMBIG_CODE] = None

                        elif ambig_code in (2, 3):
                            _ambig_code = self._reg.csStat.getMaxAmbigCodeWoSetId(comp_id, atom_id)
                            if _ambig_code not in (0, ambig_code):
                                if _ambig_code != 1:
                                    _row[CS_AMBIG_CODE] = _ambig_code
                                else:
                                    ambig_code_4_test = False
                                    _chain_id = row[chain_id_col]
                                    _seq_id = row[seq_id_col]
                                    _atom_type = row[atom_id_col][0]
                                    _ambig_code = str(ambig_code)
                                    _idx = idx
                                    for row_ in peripheral_loop_rows():
                                        if row_[comp_id_col] == comp_id\
                                           and row_[atom_id_col][0] == _atom_type\
                                           and str(row_[ambig_code_col]) == _ambig_code:
                                            if row_[chain_id_col] != _chain_id:
                                                continue
                                            if row_[seq_id_col] == _seq_id:
                                                ambig_code_4_test = True
                                                break
                                    if ambig_code_4_test:
                                        _row[CS_AMBIG_CODE] = ambig_code = 4
                                        val = float(row[val_col])
                                        sig = self._reg.ccU.getBondSignature(comp_id, atom_id)
                                        for row_ in peripheral_loop_rows():
                                            if row_[chain_id_col] == _chain_id and row_[seq_id_col] == _seq_id\
                                               and row_[comp_id_col] == comp_id and row_[atom_id_col][0] == _atom_type\
                                               and abs(float(row_[val_col]) - val) < 1.0\
                                               and self._reg.ccU.getBondSignature(comp_id, row[atom_id_col]) == sig:
                                                row[ambig_code_col] = 4
                                                reparse_request = True
                                    else:
                                        _row[CS_AMBIG_CODE] = ambig_code = 1

                        elif ambig_code == 4:
                            if not self._reg.annotation_mode and _row[CS_DETAILS] != 'UNMAPPED':
                                _chain_id = row[chain_id_col]
                                _seq_id = row[seq_id_col]
                                _atom_id = row[atom_id_col]
                                _atom_type = _atom_id[0]
                                _ambig_code = str(ambig_code)
                                _atom_ids_in_group = self._reg.ccU.getProtonsInSameGroup(comp_id, _atom_id)\
                                    if _atom_type in PROTON_BEGIN_CODE else []
                                ambig_code_4_test = hetero_group_test = False
                                _idx = idx
                                for offset in range(1, PERIPH_OFFSET_ATTEMPT):
                                    if _idx + offset < len(loop):
                                        row_ = loop.data[_idx + offset]
                                        if row_[comp_id_col] == comp_id\
                                           and row_[atom_id_col][0] == _atom_type\
                                           and str(row_[ambig_code_col]) == _ambig_code:
                                            if row_[chain_id_col] != _chain_id:
                                                break
                                            if row_[seq_id_col] == _seq_id:
                                                ambig_code_4_test = True
                                                if row_[atom_id_col] not in _atom_ids_in_group:
                                                    hetero_group_test = True
                                                    break
                                    if _idx - offset >= 0:
                                        row_ = loop.data[_idx - offset]
                                        if row_[comp_id_col] == comp_id\
                                           and row_[atom_id_col][0] == _atom_type\
                                           and str(row_[ambig_code_col]) == _ambig_code:
                                            if row_[chain_id_col] != _chain_id:
                                                break
                                            if row_[seq_id_col] == _seq_id:
                                                ambig_code_4_test = True
                                                ambig_code_4_test = True
                                                if row_[atom_id_col] not in _atom_ids_in_group:
                                                    hetero_group_test = True
                                                    break
                                if not ambig_code_4_test:
                                    ambig_code_5_test = False
                                    for row_ in peripheral_loop_rows():
                                        if row_[comp_id_col] == comp_id\
                                           and row_[atom_id_col][0] == _atom_type\
                                           and str(row_[ambig_code_col]) == _ambig_code:
                                            if row_[chain_id_col] != _chain_id:
                                                break
                                            if row_[seq_id_col] == _seq_id:
                                                break
                                            _row[CS_AMBIG_CODE] = ambig_code = 5
                                            ambig_code_5_test = True
                                            break
                                    if not ambig_code_5_test:
                                        for row_ in peripheral_loop_rows():
                                            if row_[comp_id_col] == comp_id\
                                               and row_[atom_id_col][0] == _atom_type\
                                               and str(row_[ambig_code_col]) == _ambig_code:
                                                if row_[chain_id_col] != _chain_id:
                                                    _row[CS_AMBIG_CODE] = ambig_code = 6
                                                    break
                                        if ambig_code == 4:
                                            _row[CS_AMBIG_CODE] = ambig_code = 1
                                elif not hetero_group_test:
                                    _row[CS_AMBIG_CODE] = ambig_code = 1

                        elif ambig_code == 5:
                            if not self._reg.annotation_mode and _row[CS_DETAILS] != 'UNMAPPED':
                                _chain_id = row[chain_id_col]
                                _seq_id = row[seq_id_col]
                                _atom_type = row[atom_id_col][0]
                                _ambig_code = str(ambig_code)
                                _idx = idx
                                ambig_code_5_test = False
                                for row_ in peripheral_loop_rows():
                                    if row_[comp_id_col] == comp_id\
                                       and row_[atom_id_col][0] == _atom_type\
                                       and str(row_[ambig_code_col]) == _ambig_code:
                                        if row_[chain_id_col] != _chain_id:
                                            continue
                                        if row_[seq_id_col] == _seq_id:
                                            break
                                        ambig_code_5_test = True
                                        break
                                if not ambig_code_5_test:
                                    for row_ in peripheral_loop_rows():
                                        if row_[comp_id_col] == comp_id\
                                           and row_[atom_id_col][0] == _atom_type\
                                           and str(row_[ambig_code_col]) == _ambig_code:
                                            if row_[chain_id_col] != _chain_id:
                                                _row[CS_AMBIG_CODE] = ambig_code = 6
                                                break

                        elif ambig_code == 6:
                            if len([item for item in entity_assembly
                                    if item['entity_type'] not in ('non-polymer', 'water')]) == 1\
                               and len(entity_assembly[0]['label_asym_id'].split(',')) == 1:
                                _row[CS_AMBIG_CODE] = ambig_code = 5

                        if ambig_code in (1, 2, 3):
                            if _row[CS_AMBIG_SET] is not None:
                                _row[CS_AMBIG_SET] = None

                        elif ambig_code in (4, 5, 6, 9):
                            has_genuine_ambig_code = True

                    if _row[CS_ISOTOPE] not in EMPTY_VALUE:  # DAOTHER-9520: Atom_isotppe_number is mandatory
                        lp.add_data(_row)

                        index += 1

                if trial == 0 and len(incomplete_comp_id_annotation) > 0:  # DAOTHER-9286
                    reparse_request = True

                if not reparse_request or trial > 0:
                    break

                trial += 1

            key_items = self._reg.key_items[file_type][content_subtype]

            conflict_id = self._reg.nefT.get_conflict_id(lp, lp_category, key_items)[0]

            if len(conflict_id) > 0:
                conflict_id_set = self._reg.nefT.get_conflict_id_set(lp, lp_category, key_items)[0]
                orig_atom_id_col = lp.tags.index('Original_PDB_atom_name') if 'Original_PDB_atom_name' in lp.tags else -1

                for _id in conflict_id:
                    _id_set = next(id_set for id_set in conflict_id_set if _id in id_set)

                    if len(set(str(lp.data[_id_]) for _id_ in _id_set)) == 1:
                        continue

                    # DAOTHER-9520: suppress known warnings
                    if orig_atom_id_col != -1\
                       and any(lp.data[_id1_][orig_atom_id_col] != lp.data[_id2_][orig_atom_id_col]
                               for (_id1_, _id2_) in itertools.combinations(_id_set, 2)):

                        msg = ' vs '.join([str(lp.data[_id_]).replace('None', '.').replace(',', '').replace("'", '')
                                           for _id_ in _id_set])

                        warn = f"Resolved redundancy of assigned chemical shifts ({msg}) by deletion of the latter one."

                        self._reg.report.warning.appendDescription('redundant_data',
                                                                   {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                    'category': lp_category, 'description': warn})

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.remediateCsLoop() ++ Warning  - {warn}\n")

                for _id in conflict_id:
                    del lp.data[_id]

            if not any(True for _row in lp if _row[CS_ENT_ASM_ID] in EMPTY_VALUE
                       or (isinstance(_row[CS_ENT_ASM_ID], str) and not _row[CS_ENT_ASM_ID].isdigit())):
                try:
                    lp.sort_rows(['Atom_ID', 'Atom_isotope_number', 'Comp_index_ID', 'Entity_assembly_ID'])
                except (TypeError, ValueError):
                    pass

            lp.renumber_rows('ID')

            if has_genuine_ambig_code:

                for _row in lp:

                    if _row[CS_AMBIG_CODE] not in (4, 5):
                        continue

                    ambig_code = _row[CS_AMBIG_CODE]

                    if _row[CS_AMBIG_SET] not in EMPTY_VALUE:
                        ambig_code = copy.copy(_row[CS_AMBIG_CODE])
                        ambig_set_id = copy.copy(_row[CS_AMBIG_SET])
                        chain_id = _row[CS_ENT_ASM_ID]
                        seq_id = _row[CS_COMP_IDX]
                        comp_id = _row[CS_COMP_ID]
                        atom_id = _row[CS_ATOM_ID]
                        atom_type = _row[CS_ATOM_TYPE]

                        if atom_type == 'H':
                            atom_in_same_group = self._reg.csStat.getProtonsInSameGroup(comp_id, atom_id)

                            if not any(((ambig_code == 5 and (__row[CS_ENT_ASM_ID] != chain_id or __row[CS_COMP_IDX] != seq_id))
                                        or (ambig_code == 4 and __row[CS_ATOM_ID] not in atom_in_same_group))
                                       for __row in lp
                                       if __row[CS_AMBIG_CODE] == ambig_code and __row[CS_AMBIG_SET] == ambig_set_id):

                                for __row in lp:
                                    if __row[CS_AMBIG_CODE] == ambig_code and __row[CS_AMBIG_SET] == ambig_set_id:
                                        __row[CS_AMBIG_CODE] = self._reg.csStat.getMaxAmbigCodeWoSetId(comp_id, atom_id, None)
                                        __row[CS_AMBIG_SET] = None

                                if not isinstance(sf, pynmrstar.Loop)\
                                   and any(True for aux_loop in sf if aux_loop.category == aux_lp_category):

                                    aux_loop = sf.get_loop(aux_lp_category)

                                    if 'Ambiguous_shift_set_ID' in aux_loop.tags:
                                        ambig_set_id_col = aux_loop.tags.index('Ambiguous_shift_set_ID')

                                        del_row_idx = []

                                        for idx, __row in enumerate(aux_loop):
                                            if __row[ambig_set_id_col] == ambig_set_id:
                                                del_row_idx.append(idx)

                                        if len(del_row_idx) > 0:
                                            for idx in reversed(del_row_idx):
                                                del aux_loop.data[idx]

                                            if len(aux_loop) == 0:
                                                delete_aux_loop()

                has_genuine_ambig_code = False

                for _row in lp:

                    if _row[CS_AMBIG_CODE] not in (4, 5, 6, 9):
                        continue

                    ambig_code = _row[CS_AMBIG_CODE]
                    _ambig_code = self._reg.csStat.getMaxAmbigCodeWoSetId(comp_id, atom_id)

                    chain_id = _row[CS_ENT_ASM_ID]
                    seq_id = _row[CS_COMP_IDX]
                    comp_id = _row[CS_COMP_ID]
                    atom_id = _row[CS_ATOM_ID]
                    atom_type = _row[CS_ATOM_TYPE]

                    if ambig_code == 4:
                        _atom_id_set_w_same_ambig_code =\
                            set(_row_[CS_ATOM_ID] for _row_ in lp
                                if _row != _row_ and _row_[CS_ENT_ASM_ID] == chain_id and _row_[CS_COMP_IDX] == seq_id
                                and _row_[CS_ATOM_TYPE] == atom_type and _row_[CS_AMBIG_CODE] == ambig_code)

                        if atom_type == 'H':
                            atom_in_same_group = self._reg.csStat.getProtonsInSameGroup(comp_id, atom_id, True)

                            if len(_atom_id_set_w_same_ambig_code - set(atom_in_same_group)) == 0:
                                if _ambig_code > 1:
                                    _row[CS_AMBIG_CODE] = _ambig_code
                                    _row[CS_AMBIG_SET] = None

                        else:
                            geminal_atom = self._reg.csStat.getGeminalAtom(comp_id, atom_id)

                            if geminal_atom is not None and len(_atom_id_set_w_same_ambig_code - set([geminal_atom])) == 0:
                                if _ambig_code > 1:
                                    _row[CS_AMBIG_CODE] = _ambig_code
                                    _row[CS_AMBIG_SET] = None

                    elif ambig_code == 5:
                        _atom_id_set_w_same_ambig_code =\
                            set(_row_[CS_ATOM_ID] for _row_ in lp
                                if _row != _row_ and _row_[CS_ENT_ASM_ID] == chain_id and _row_[CS_COMP_IDX] != seq_id
                                and _row_[CS_ATOM_TYPE] == atom_type and _row_[CS_AMBIG_CODE] == ambig_code)

                        if len(_atom_id_set_w_same_ambig_code) == 0 and _ambig_code != 0:
                            _row[CS_AMBIG_CODE] = _ambig_code
                            _row[CS_AMBIG_SET] = None

                    else:
                        _row[CS_AMBIG_SET] = None

                    if _row[CS_AMBIG_CODE] in (4, 5):
                        has_genuine_ambig_code = True

                if has_genuine_ambig_code:

                    aux_lp = pynmrstar.Loop.from_scratch(aux_lp_category)

                    aux_items = ['Ambiguous_shift_set_ID', 'Atom_chem_shift_ID', 'Entry_ID', 'Assigned_chem_shift_list_ID']

                    aux_tags = [f'{aux_lp_category}.{item}' for item in aux_items]

                    aux_lp.add_tag(aux_tags)

                    inter_residue_seq_id = {}

                    for _row in lp:

                        if _row[CS_AMBIG_CODE] != 5:
                            continue

                        chain_id = _row[CS_ENT_ASM_ID]
                        seq_id = _row[CS_COMP_IDX]

                        if chain_id not in inter_residue_seq_id:
                            inter_residue_seq_id[chain_id] = set()

                        inter_residue_seq_id[chain_id].add(seq_id)

                    if len(inter_residue_seq_id) > 0:

                        for k, v in inter_residue_seq_id.items():
                            if len(v) == 1:
                                chain_id = k
                                seq_id = list(v)[0]

                                for _row in lp:

                                    if _row[CS_AMBIG_CODE] != 5:
                                        continue

                                    if _row[CS_ENT_ASM_ID] == chain_id and _row[CS_COMP_IDX] == seq_id:
                                        _row[CS_AMBIG_CODE] = 4

                    aux_index_id = 0
                    ambig_shift_set_id = {}

                    for _idx, _row in enumerate(lp):

                        if _row[CS_AMBIG_CODE] not in (4, 5):
                            continue

                        ambig_code = _row[CS_AMBIG_CODE]

                        chain_id = _row[CS_ENT_ASM_ID]
                        seq_id = _row[CS_COMP_IDX]
                        atom_type = _row[CS_ATOM_TYPE]

                        if ambig_code == 4:
                            key = (chain_id, str(seq_id), atom_type, ambig_code)
                        else:
                            key = (chain_id, str(inter_residue_seq_id[chain_id]), atom_type, ambig_code)

                        if key not in ambig_shift_set_id:
                            aux_index_id += 1
                            ambig_shift_set_id[key] = aux_index_id

                        lp.data[_idx][CS_AMBIG_SET] = ambig_shift_set_id[key]

                        _aux_row = [None] * 4
                        _aux_row[0], _aux_row[1], _aux_row[2], _aux_row[3] =\
                            ambig_shift_set_id[key], _row[CS_ID], self._reg.entry_id, list_id

                        aux_lp.add_data(_aux_row)

                else:
                    delete_aux_loop()

            else:
                delete_aux_loop()

        del sf[loop]

        sf.add_loop(lp)

        if aux_lp is not None and len(aux_lp) > 0:
            delete_aux_loop()

            sf.add_loop(aux_lp)

        if file_type == 'nmr-star':
            val = get_first_sf_tag(sf, 'ID')
            if (isinstance(val, int) and val == list_id) or (isinstance(val, str) and val.isdigit() and int(val) == list_id):
                pass
            else:
                set_sf_tag(sf, 'ID', list_id)

        if not self._reg.native_combined:
            self._reg.dpV.testDataConsistencyInLoop(file_list_id, file_name, file_type, content_subtype,
                                                    sf, sf_framecode, lp_category, list_id)

        get_auth_seq_scheme.cache_clear()
        get_label_seq_scheme.cache_clear()

        return True
