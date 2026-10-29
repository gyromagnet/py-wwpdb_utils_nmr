##
# File: NmrDpRemediationLegacyPk.py
# Date: 17-Aug-2026
#
# Updates:
##
""" Validation of legacy spectral peak list files during NMR data remediation.
    @author: Masashi Yokochi
"""
__docformat__ = "restructuredtext en"
__author__ = "Masashi Yokochi"
__email__ = "yokochi@protein.osaka-u.ac.jp"
__license__ = "Apache License 2.0"
__version__ = "5.4.4"

import copy
import os
from operator import itemgetter
from typing import Optional, Tuple

import pynmrstar

try:
    from wwpdb.utils.nmr.NmrDpConstant import (AR_FILE_PATH_LIST_KEY,
                                               SF_CATEGORIES,
                                               EMPTY_VALUE,
                                               SEQ_MISMATCH_WARNING_PAT,
                                               INCONSISTENT_RESTRAINT_WARNING_PAT,
                                               MISMATCHED_INPUT_ERR_MSG,
                                               EXTRANEOUS_INPUT_ERR_MSG)
    from wwpdb.utils.nmr.NmrDpMrSplitter import detect_encoding
    from wwpdb.utils.nmr.AlignUtil import (getScoreOfSeqAlign,
                                           getRestraintFormatName,
                                           updatePolySeqRst,
                                           sortPolySeqRst,
                                           alignPolymerSequence,
                                           assignPolymerSequence,
                                           trimSequenceAlignment)
    from wwpdb.utils.nmr.CifToNmrStar import (has_key_value,
                                              get_first_sf_tag,
                                              set_sf_tag)
    from wwpdb.utils.nmr.mr.ParserListenerUtil import contentSubtypeOf
    from wwpdb.utils.nmr.pk.AriaPKReader import AriaPKReader
    from wwpdb.utils.nmr.pk.BarePKReader import BarePKReader
    from wwpdb.utils.nmr.pk.CcpnPKReader import CcpnPKReader
    from wwpdb.utils.nmr.pk.OliviaPKReader import OliviaPKReader
    from wwpdb.utils.nmr.pk.NmrPipePKReader import NmrPipePKReader
    from wwpdb.utils.nmr.pk.PonderosaPKReader import PonderosaPKReader
    from wwpdb.utils.nmr.pk.NmrViewPKReader import NmrViewPKReader
    from wwpdb.utils.nmr.pk.NmrViewNPKReader import NmrViewNPKReader
    from wwpdb.utils.nmr.pk.SparkyPKReader import SparkyPKReader
    from wwpdb.utils.nmr.pk.SparkyNPKReader import SparkyNPKReader
    from wwpdb.utils.nmr.pk.SparkyRPKReader import SparkyRPKReader
    from wwpdb.utils.nmr.pk.SparkySPKReader import SparkySPKReader
    from wwpdb.utils.nmr.pk.TopSpinPKReader import TopSpinPKReader
    from wwpdb.utils.nmr.pk.VnmrPKReader import VnmrPKReader
    from wwpdb.utils.nmr.pk.XeasyPKReader import XeasyPKReader
    from wwpdb.utils.nmr.pk.XeasyPROTReader import XeasyPROTReader
    from wwpdb.utils.nmr.pk.XwinNmrPKReader import XwinNmrPKReader
    from wwpdb.utils.nmr.NmrDpRemediationBase import NmrDpRemediationBase
except ImportError:
    from nmr.NmrDpConstant import (AR_FILE_PATH_LIST_KEY,
                                   SF_CATEGORIES,
                                   EMPTY_VALUE,
                                   SEQ_MISMATCH_WARNING_PAT,
                                   INCONSISTENT_RESTRAINT_WARNING_PAT,
                                   MISMATCHED_INPUT_ERR_MSG,
                                   EXTRANEOUS_INPUT_ERR_MSG)
    from nmr.NmrDpMrSplitter import detect_encoding
    from nmr.AlignUtil import (getScoreOfSeqAlign,
                               getRestraintFormatName,
                               updatePolySeqRst,
                               sortPolySeqRst,
                               alignPolymerSequence,
                               assignPolymerSequence,
                               trimSequenceAlignment)
    from nmr.CifToNmrStar import (has_key_value,
                                  get_first_sf_tag,
                                  set_sf_tag)
    from nmr.mr.ParserListenerUtil import contentSubtypeOf
    from nmr.pk.AriaPKReader import AriaPKReader
    from nmr.pk.BarePKReader import BarePKReader
    from nmr.pk.CcpnPKReader import CcpnPKReader
    from nmr.pk.OliviaPKReader import OliviaPKReader
    from nmr.pk.NmrPipePKReader import NmrPipePKReader
    from nmr.pk.PonderosaPKReader import PonderosaPKReader
    from nmr.pk.NmrViewPKReader import NmrViewPKReader
    from nmr.pk.NmrViewNPKReader import NmrViewNPKReader
    from nmr.pk.SparkyPKReader import SparkyPKReader
    from nmr.pk.SparkyNPKReader import SparkyNPKReader
    from nmr.pk.SparkyRPKReader import SparkyRPKReader
    from nmr.pk.SparkySPKReader import SparkySPKReader
    from nmr.pk.TopSpinPKReader import TopSpinPKReader
    from nmr.pk.VnmrPKReader import VnmrPKReader
    from nmr.pk.XeasyPKReader import XeasyPKReader
    from nmr.pk.XeasyPROTReader import XeasyPROTReader
    from nmr.pk.XwinNmrPKReader import XwinNmrPKReader
    from nmr.NmrDpRemediationBase import NmrDpRemediationBase

# Reader dispatch for the collapsible branches of validateLegacyPk(). The keys
# beyond 'reader' and 'label' record this format's deviations from the common
# parse / re-parse sequence:
#
#   label=None             the failure message carries no format name
#   cs_loops=False         parse() is not given csLoops
#   xml_like               lexer errors are ignored (incomplete XML formats):
#                          the lexer listener is discarded, the parse result is
#                          not used to refine _content_subtype, and
#                          deal_lexer_or_parser_error() is passed None for it
#   lexer_var              name to unpack the lexer listener into when xml_like
#   needs_xeasy_dict       the reader takes the XEASY atom number dictionary
#   reparse_needs_warning  re-parse only if the listener also emitted warnings
#
# nm-pea-bar, nm-pea-spa and nm-pea-vie are deliberately absent: they keep their
# own branches below because their control flow genuinely differs.
LEGACY_PK_READERS = {
    'nm-pea-ari': {'reader': AriaPKReader, 'label': 'ARIA', 'xml_like': True},
    'nm-pea-ccp': {'reader': CcpnPKReader, 'label': 'CCPN'},
    'nm-pea-oli': {'reader': OliviaPKReader, 'label': 'OLIVIA'},
    'nm-pea-pip': {'reader': NmrPipePKReader, 'label': 'NMRPIPE'},
    'nm-pea-pon': {'reader': PonderosaPKReader, 'label': 'PONDEROSA'},
    'nm-pea-sps': {'reader': SparkySPKReader, 'label': 'SPARKY'},
    'nm-pea-top': {'reader': TopSpinPKReader, 'label': 'TOPSPIN', 'xml_like': True,
                   'cs_loops': False},
    'nm-pea-vnm': {'reader': VnmrPKReader, 'label': 'VNMR'},
    'nm-pea-xea': {'reader': XeasyPKReader, 'label': 'XEASY',
                   'needs_xeasy_dict': True, 'reparse_needs_warning': True},
    'nm-pea-xwi': {'reader': XwinNmrPKReader, 'label': 'XWINNMR', 'cs_loops': False},
}


class NmrDpRemediationLegacyPk(NmrDpRemediationBase):
    """ Validation of legacy spectral peak list files during NMR data remediation.
    """
    __slots__ = ()

    def _parseLegacyPk(self, spec: dict, extra_args: tuple, file_path: str, file_name: str,
                       original_file_name: str, create_sf_dict: bool,
                       reserved_list_ids: dict, content_subtype: Optional[dict],
                       a_pk_format_name: str, deal_lexer_or_parser_error,
                       deal_warn_for_lazy_eval) -> Tuple[bool, bool, Optional[object]]:
        """ Parse a legacy spectral peak list file with the reader for a given file type,
            re-parsing once if the parser listener asks for it.
            @param spec: the LEGACY_PK_READERS entry for the file type
            @param extra_args: extra positional arguments for the reader constructor
            @param content_subtype: the content subtype declared by the input source,
                                    used as-is for xml_like formats and otherwise
                                    refined from the initial parse
            @return: (whether to skip this file, whether the initial parse yielded a
                     listener, the final listener). The second and third are reported
                     separately because the caller's processing was gated on the
                     *initial* parse, exactly as the per-format branches were.
        """
        reader_cls = spec['reader']
        xml_like = spec.get('xml_like', False)
        cs_loops = spec.get('cs_loops', True)

        def new_reader(*trailing):
            reader = reader_cls(self._reg.verbose, self._reg.log,
                                self._reg.representative_model_id,
                                self._reg.representative_alt_id,
                                self._reg.mr_atom_name_mapping,
                                self._reg.cR, self._reg.caC,
                                self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                *extra_args, *trailing)
            reader.enforcePeakRowFormat(self._reg.enforce_peak_row_format)
            reader.setInternalMode(self._reg.internal_mode)
            return reader

        def parse_with(reader, list_id_counter):
            kwargs = {'createSfDict': create_sf_dict, 'originalFileName': original_file_name,
                      'listIdCounter': list_id_counter, 'reservedListIds': reserved_list_ids,
                      'entryId': self._reg.entry_id}
            if cs_loops:
                kwargs['csLoops'] = self._reg.lp_data['chem_shift']
            return reader.parse(file_path, self._reg.cifPath, **kwargs)

        _list_id_counter = copy.copy(self._reg.list_id_counter)

        listener, parser_err_listener, lexer_err_listener = parse_with(new_reader(),
                                                                       self._reg.list_id_counter)

        if not xml_like:
            content_subtype = listener.getContentSubtype() if listener is not None else None
            if content_subtype is not None and len(content_subtype) == 0:
                content_subtype = None

        if xml_like:
            checked = None not in (parser_err_listener, listener)\
                and (parser_err_listener.getMessageList() is None or content_subtype is not None)
            lexer_arg = None
        else:
            checked = None not in (lexer_err_listener, parser_err_listener, listener)\
                and ((lexer_err_listener.getMessageList() is None and parser_err_listener.getMessageList() is None)
                     or content_subtype is not None)
            lexer_arg = lexer_err_listener

        if checked and deal_lexer_or_parser_error(a_pk_format_name, file_name,
                                                  lexer_arg, parser_err_listener)[0]:
            return True, False, None

        if listener is None:
            return False, False, None

        reasons = listener.getReasonsForReparsing()

        reparse = reasons is not None
        if reparse and spec.get('reparse_needs_warning', False):
            reparse = listener.warningMessage is not None and len(listener.warningMessage) > 0

        if reparse:
            deal_warn_for_lazy_eval(file_name, listener)

            listener = parse_with(new_reader(reasons), _list_id_counter)[0]

        return False, True, listener

    def validateLegacyPk(self) -> bool:
        """ Validate data content of legacy spectral peak files and merge them if possible.
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

        file_type = 'nmr-star'
        content_subtype = 'spectral_peak'

        sf_category = SF_CATEGORIES[file_type][content_subtype]

        rlist_ids = []
        if len(self._reg.star_data) > 0 and isinstance(self._reg.star_data[0], pynmrstar.Entry):
            for idx, sf in enumerate(self._reg.star_data[0].get_saveframes_by_category(sf_category), start=1):
                list_id = get_first_sf_tag(sf, 'ID')
                rlist_ids.append(int(list_id) if list_id not in EMPTY_VALUE else idx)

        reserved_list_ids = {content_subtype: rlist_ids} if len(rlist_ids) > 0 else None

        xeasyAtomNumberDict = None

        has_aux_xea = has_pea_xea = False

        fileListId = self._reg.file_path_list_len

        for ar in self._reg.inputParamDict[AR_FILE_PATH_LIST_KEY]:
            file_path = ar['file_name']

            input_source = self._reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_type = input_source_dic['file_type']

            if file_type == 'nm-pea-xea':
                has_pea_xea = True
                break

            fileListId += 1

        def deal_lexer_or_parser_error(a_pk_format_name, file_name, lexer_err_listener, parser_err_listener):
            _type = 'default'
            _err = ''
            if lexer_err_listener is not None:
                messageList = lexer_err_listener.getMessageList()

                if messageList is not None:
                    for description in messageList:
                        _err = f"[Syntax error as {a_pk_format_name} file] "\
                               f"line {description['line_number']}:{description['column_position']} {description['message']}\n"
                        if 'input' in description:
                            enc = detect_encoding(description['input'])
                            is_not_ascii = False
                            if enc is not None and enc != 'ascii':
                                _err += f"{description['input']}\n".encode().decode('ascii', 'backslashreplace')
                                is_not_ascii = True
                            else:
                                _err += f"{description['input']}\n"
                            _err += f"{description['marker']}\n"
                            if is_not_ascii:
                                _err += f"[Unexpected text encoding] Encoding used in the above line is {enc!r} "\
                                    "and must be 'ascii'.\n"

            if parser_err_listener is not None and len(_err) == 0:
                messageList = parser_err_listener.getMessageList()

                if messageList is not None:
                    for description in messageList:
                        if 'SPARKY' in a_pk_format_name\
                           and (MISMATCHED_INPUT_ERR_MSG in description['message']
                                or EXTRANEOUS_INPUT_ERR_MSG in description['message'])\
                           and "expecting {Lw1_Hz_LA, Lw2_Hz_LA, Lw3_Hz_LA, Lw4_Hz_LA" in description['message']\
                           and ('Height' in description['message']
                                or 'Data' in description['message']):
                            _type = 'reverse'
                        elif 'SPARKY' in a_pk_format_name\
                                and MISMATCHED_INPUT_ERR_MSG in description['message']\
                                and "'\\n' expecting {Integer, Float, Real, Real_vol}" in description['message']:
                            _type = 'no'
                        elif 'NMRVIEW' in a_pk_format_name\
                                and MISMATCHED_INPUT_ERR_MSG in description['message']\
                                and "' expecting L_brace" in description['message']:
                            _type = 'no_brace'
                        else:
                            _err += f"[Syntax error as {a_pk_format_name} file] "\
                                    f"line {description['line_number']}:{description['column_position']} "\
                                    f"{description['message']}\n"
                            if 'input' in description:
                                _err += f"{description['input']}\n"
                                _err += f"{description['marker']}\n"

            if len(_err) == 0:
                return False, _type

            err = f"The spectral peak list file {file_name!r} looks like {a_pk_format_name} file. "\
                "Please re-upload the spectral peak list file.\n"\
                "The following issues need to be fixed before re-upload.\n" + _err[:-1]

            self._reg.report.error.appendDescription('format_issue',
                                                     {'file_name': file_name, 'description': err})

            self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Error  - {file_name} {err}\n")

            return True, _type

        def deal_aux_warn_message(file_name, listener):

            if listener.warningMessage is not None:

                for warn in listener.warningMessage:

                    msg_dict = {'file_name': file_name, 'description': warn, 'inheritable': True}

                    if warn.startswith('[Concatenated sequence]'):
                        self._reg.report.warning.appendDescription('concatenated_sequence', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Sequence mismatch]'):
                        self._reg.report.warning.appendDescription('sequence_mismatch', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Unknown atom name]'):
                        self._reg.report.warning.appendDescription('inconsistent_peak_list', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Unknown residue name]'):
                        self._reg.report.warning.appendDescription('inconsistent_peak_list', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Warning  - {warn}\n")

                    else:
                        self._reg.report.error.appendDescription('internal_error',
                                                                 f"+{self.__class_name__}.validateLegacyPk() "
                                                                 "++ KeyError  - " + warn)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ KeyError  - {warn}\n")

        if has_pea_xea:

            fileListId = self._reg.file_path_list_len

            for ar in self._reg.inputParamDict[AR_FILE_PATH_LIST_KEY]:
                file_path = self.testPathWithSuffix(ar['file_name'], '-corrected')

                input_source = self._reg.report.input_sources[fileListId]
                input_source_dic = input_source.get()

                file_type = input_source_dic['file_type']

                fileListId += 1

                if file_type == 'nm-aux-xea':
                    has_aux_xea = True

                    file_name = input_source_dic['file_name']

                    original_file_name = None
                    if 'original_file_name' in input_source_dic:
                        if input_source_dic['original_file_name'] is not None:
                            original_file_name = os.path.basename(input_source_dic['original_file_name'])
                        if file_name != original_file_name and original_file_name is not None:
                            file_name = f"{original_file_name} ({file_name})"

                    reader = XeasyPROTReader(self._reg.verbose, self._reg.log,
                                             self._reg.representative_model_id,
                                             self._reg.representative_alt_id,
                                             self._reg.mr_atom_name_mapping,
                                             self._reg.cR, self._reg.caC,
                                             self._reg.ccU, self._reg.csStat, self._reg.nefT)

                    listener, parser_err_listener, lexer_err_listener = reader.parse(file_path, self._reg.cifPath)

                    _content_subtype = listener.getContentSubtype() if listener is not None else None
                    if _content_subtype is not None and len(_content_subtype) == 0:
                        _content_subtype = None

                    if None not in (lexer_err_listener, parser_err_listener, listener)\
                       and ((lexer_err_listener.getMessageList() is None and parser_err_listener.getMessageList() is None)
                            or _content_subtype is not None):
                        _pk_format_name = getRestraintFormatName(file_type)
                        pk_format_name = _pk_format_name.split()[0]
                        a_pk_format_name = ('an ' if pk_format_name[0] in ('AINMX') else 'a ') + _pk_format_name
                        if deal_lexer_or_parser_error(a_pk_format_name, file_name, lexer_err_listener, parser_err_listener)[0]:
                            continue

                    if listener is not None:

                        deal_aux_warn_message(file_name, listener)

                        xeasyAtomNumberDict = listener.getAtomNumberDict()

                    break

        poly_seq_set = []

        create_sf_dict = self._reg.remediation_mode

        if self._reg.list_id_counter is None:
            self._reg.list_id_counter = {}

        pk_sf_dict_holder = {}

        proc_nmr_ext_poly_seq = False

        if self._reg.nmr_ext_poly_seq is None and not self._reg.bmrb_only or not self._reg.internal_mode:
            proc_nmr_ext_poly_seq = True

            self._reg.nmr_ext_poly_seq = []

            input_source = self._reg.report.input_sources[0]
            input_source_dic = input_source.get()

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

        suspended_errors_for_lazy_eval = []

        def consume_suspended_message():

            if len(suspended_errors_for_lazy_eval) > 0:
                for msg in suspended_errors_for_lazy_eval:
                    for k, v in msg.items():
                        self._reg.report.error.appendDescription(k, v)
                suspended_errors_for_lazy_eval.clear()

        def deal_pea_warn_message(file_name, listener, ignore_error):

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
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Sequence mismatch]'):
                        # consume_suspended_message()

                        self._reg.report.error.appendDescription('sequence_mismatch', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Error  - {warn}\n")

                    elif warn.startswith('[Atom not found]'):
                        if not self._reg.remediation_mode or 'Macromolecules page' not in warn:
                            consume_suspended_message()

                            self._reg.report.warning.appendDescription('assigned_peak_atom_not_found', msg_dict)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Warning  - {warn}\n")
                        else:
                            self._reg.report.warning.appendDescription('sequence_mismatch', msg_dict)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Hydrogen not instantiated]'):
                        if (self._reg.remediation_mode or self._reg.internal_mode) and not self._reg.conversion_server:
                            pass
                        else:
                            consume_suspended_message()

                        self._reg.report.warning.appendDescription('hydrogen_not_instantiated', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Coordinate issue]'):
                        # consume_suspended_message()

                        self._reg.report.warning.appendDescription('coordinate_issue', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Invalid atom nomenclature]'):
                        consume_suspended_message()

                        # DAOTHER-8905: change warning level from 'invalid_atom_nomenclature' error
                        # to 'atom_nomenclature_mismatch' warning
                        # because we accept atom nomenclature provided by depositor for peak list
                        self._reg.report.warning.appendDescription('atom_nomenclature_mismatch', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Invalid atom selection]') or warn.startswith('[Invalid data]'):
                        # consume_suspended_message()

                        self._reg.report.warning.appendDescription('inconsistent_peak_list', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Sequence mismatch warning]'):
                        self._reg.report.warning.appendDescription('sequence_mismatch', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Warning  - {warn}\n")

                        if SEQ_MISMATCH_WARNING_PAT.match(warn):
                            g = SEQ_MISMATCH_WARNING_PAT.search(warn).groups()
                            d = {'auth_chain_id': g[2],
                                 'auth_seq_id': int(g[0]),
                                 'auth_comp_id': g[1]}
                            if d not in self._reg.nmr_ext_poly_seq:
                                self._reg.nmr_ext_poly_seq.append(d)

                    elif warn.startswith('[Inconsistent peak assignment]'):
                        self._reg.report.warning.appendDescription('inconsistent_peak_list', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Conflicted peak assignment]'):
                        self._reg.report.warning.appendDescription('conflicted_peak_list', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Warning  - {warn}\n")

                    elif warn.startswith('[Missing data]'):
                        if (self._reg.remediation_mode or self._reg.internal_mode) and not self._reg.conversion_server:
                            pass
                        else:
                            self._reg.report.error.appendDescription('missing_mandatory_item', msg_dict)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Error  - {warn}\n")

                    elif warn.startswith('[Range value error]') and not self._reg.remediation_mode:
                        # consume_suspended_message()

                        self._reg.report.error.appendDescription('anomalous_data', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ ValueError  - {warn}\n")

                    elif warn.startswith('[Range value warning]')\
                            or (warn.startswith('[Range value error]') and self._reg.remediation_mode):
                        self._reg.report.warning.appendDescription('inconsistent_peak_list', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Warning  - {warn}\n")

                    elif not ignore_error:
                        self._reg.report.error.appendDescription('internal_error',
                                                                 f"+{self.__class_name__}.validateLegacyPk() "
                                                                 "++ KeyError  - " + warn)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ KeyError  - {warn}\n")

        def deal_pea_warn_message_for_lazy_eval(file_name, listener):

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

                    # elif warn.startswith('[Atom not found]'):
                    #     if not self._reg.remediation_mode or 'Macromolecules page' not in warn:
                    #         suspended_errors_for_lazy_eval.append({'atom_not_found': msg_dict})

                    # elif warn.startswith('[Hydrogen not instantiated]'):
                    #     if self._reg.remediation_mode:
                    #         pass
                    #     else:
                    #         suspended_errors_for_lazy_eval.append({'hydrogen_not_instantiated': msg_dict})

                    # elif warn.startswith('[Coordinate issue]'):
                    #     suspended_errors_for_lazy_eval.append({'coordinate_issue': msg_dict})

                    # elif warn.startswith('[Invalid atom nomenclature]'):
                    #     suspended_errors_for_lazy_eval.append({'invalid_atom_nomenclature': msg_dict})

                    elif warn.startswith('[Invalid atom selection]') or warn.startswith('[Invalid data]'):
                        suspended_errors_for_lazy_eval.append({'invalid_data': msg_dict})

                    # elif warn.startswith('[Range value error]') and not self._reg.remediation_mode:
                    #     suspended_errors_for_lazy_eval.append({'anomalous_data': msg_dict})

        fileListId = self._reg.file_path_list_len

        for ar in self._reg.inputParamDict[AR_FILE_PATH_LIST_KEY]:
            file_path = ar['file_name']

            input_source = self._reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_type = input_source_dic['file_type']

            ignore_error = False if 'ignore_error' not in input_source_dic else input_source_dic['ignore_error']

            fileListId += 1

            if file_type is None or file_type.startswith('nm-res') or file_type.startswith('nm-aux'):
                continue

            if self._reg.remediation_mode and os.path.exists(self.testPathWithSuffix(file_path, '-ignored', True)):
                continue

            file_path = self.testPathWithSuffix(file_path, '-corrected')

            file_name = input_source_dic['file_name']

            original_file_name = os.path.basename(file_path)

            if file_type == 'nm-pea-any':

                warn = f"We could not identify peak list file format of {file_name!r}. "\
                    "In order to add file format support in the future, "\
                    "the contents is temporarily stored as-is in the _Spectral_peak_list.Text_data tag "\
                    "and will be converted during future data remediation if the data matches a known peak list format."

                self._reg.report.warning.appendDescription('unsupported_peak_list',
                                                           {'file_name': file_name, 'description': warn, 'inheritable': True})

                if self._reg.verbose:
                    self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Warning  - {warn}\n")

                continue

            if file_type == 'nm-pea-xea' and not has_aux_xea and not self._reg.internal_mode:

                err = f"XEASY PROT file should be uploaded to verify XEASY spectral peak list file {file_name!r}."

                suspended_errors_for_lazy_eval.append({'missing_mandatory_content':
                                                       {'file_name': file_name, 'description': err}})

            _content_subtype = input_source_dic['content_subtype']

            if _content_subtype is None or len(_content_subtype) == 0:
                continue

            _pk_format_name = getRestraintFormatName(file_type)
            pk_format_name = _pk_format_name.split()[0]
            a_pk_format_name = ('an ' if pk_format_name[0] in ('AINMX') else 'a ') + _pk_format_name

            suspended_errors_for_lazy_eval.clear()

            spec = LEGACY_PK_READERS.get(file_type)

            if spec is not None:
                extra_args = (xeasyAtomNumberDict,) if spec.get('needs_xeasy_dict') else ()

                skip, parsed, listener = self._parseLegacyPk(
                    spec, extra_args, file_path, file_name, original_file_name,
                    create_sf_dict, reserved_list_ids, _content_subtype,
                    a_pk_format_name, deal_lexer_or_parser_error,
                    deal_pea_warn_message_for_lazy_eval)

                if skip:
                    continue

                if parsed:
                    deal_pea_warn_message(file_name, listener, ignore_error)

                    poly_seq = listener.getPolymerSequence()
                    if poly_seq is not None:
                        input_source.setItemValue('polymer_sequence', poly_seq)
                        poly_seq_set.append(poly_seq)

                    seq_align = listener.getSequenceAlignment()
                    if seq_align is not None:
                        self._reg.report.sequence_alignment.setItemValue(f'model_poly_seq_vs_{content_subtype}', seq_align)

                    if create_sf_dict:
                        if len(listener.getContentSubtype()) == 0 and not ignore_error:
                            label = spec['label']
                            err = f"Failed to validate spectral peak list file ({label}) {file_name!r}."

                            self._reg.report.error.appendDescription('internal_error',
                                                                     f"+{self.__class_name__}.validateLegacyPk() "
                                                                     "++ Error  - " + err)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Error  - {err}\n")

                        self._reg.list_id_counter, sf_dict = listener.getSfDict()
                        if sf_dict is not None:
                            for k, v in sf_dict.items():
                                content_subtype = contentSubtypeOf(k[0])
                                if content_subtype not in pk_sf_dict_holder:
                                    pk_sf_dict_holder[content_subtype] = []
                                for sf in v:
                                    if sf not in pk_sf_dict_holder[content_subtype]:
                                        pk_sf_dict_holder[content_subtype].append(sf)

            elif file_type == 'nm-pea-bar':
                reader = BarePKReader(self._reg.verbose, self._reg.log,
                                      self._reg.representative_model_id,
                                      self._reg.representative_alt_id,
                                      self._reg.mr_atom_name_mapping,
                                      self._reg.cR, self._reg.caC,
                                      self._reg.ccU, self._reg.csStat, self._reg.nefT)
                reader.enforcePeakRowFormat(self._reg.enforce_peak_row_format)
                reader.setInternalMode(self._reg.internal_mode)

                _list_id_counter = copy.copy(self._reg.list_id_counter)

                # ignore lexer error because of incomplete XML file format
                listener, parser_err_listener, lexer_err_listener =\
                    reader.parse(file_path, self._reg.cifPath,
                                 createSfDict=create_sf_dict, originalFileName=original_file_name,
                                 listIdCounter=self._reg.list_id_counter, reservedListIds=reserved_list_ids,
                                 entryId=self._reg.entry_id,
                                 csLoops=self._reg.lp_data['chem_shift'])

                if None not in (lexer_err_listener, parser_err_listener, listener)\
                   and ((lexer_err_listener.getMessageList() is None and parser_err_listener.getMessageList() is None)
                        or _content_subtype is not None):
                    if deal_lexer_or_parser_error(a_pk_format_name, file_name, None, parser_err_listener)[0]:
                        continue

                if listener is not None:
                    reasons = listener.getReasonsForReparsing()

                    if reasons is not None:
                        deal_pea_warn_message_for_lazy_eval(file_name, listener)

                        reader = BarePKReader(self._reg.verbose, self._reg.log,
                                              self._reg.representative_model_id,
                                              self._reg.representative_alt_id,
                                              self._reg.mr_atom_name_mapping,
                                              self._reg.cR, self._reg.caC,
                                              self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                              reasons)
                        reader.enforcePeakRowFormat(self._reg.enforce_peak_row_format)
                        reader.setInternalMode(self._reg.internal_mode)

                        listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                                      createSfDict=create_sf_dict, originalFileName=original_file_name,
                                                      listIdCounter=_list_id_counter, reservedListIds=reserved_list_ids,
                                                      entryId=self._reg.entry_id,
                                                      csLoops=self._reg.lp_data['chem_shift'])

                    deal_pea_warn_message(file_name, listener, ignore_error)

                    poly_seq = listener.getPolymerSequence()
                    if poly_seq is not None:
                        input_source.setItemValue('polymer_sequence', poly_seq)
                        poly_seq_set.append(poly_seq)

                    seq_align = listener.getSequenceAlignment()
                    if seq_align is not None:
                        self._reg.report.sequence_alignment.setItemValue(f'model_poly_seq_vs_{content_subtype}', seq_align)

                    if create_sf_dict:
                        if len(listener.getContentSubtype()) == 0 and not ignore_error:
                            err = f"Failed to validate spectral peak list file {file_name!r}."

                            self._reg.report.error.appendDescription('internal_error',
                                                                     f"+{self.__class_name__}.validateLegacyPk() "
                                                                     "++ Error  - " + err)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Error  - {err}\n")

                        self._reg.list_id_counter, sf_dict = listener.getSfDict()
                        if sf_dict is not None:
                            for k, v in sf_dict.items():
                                content_subtype = contentSubtypeOf(k[0])
                                if content_subtype not in pk_sf_dict_holder:
                                    pk_sf_dict_holder[content_subtype] = []
                                for sf in v:
                                    if sf not in pk_sf_dict_holder[content_subtype]:
                                        pk_sf_dict_holder[content_subtype].append(sf)

            elif file_type == 'nm-pea-spa':
                __list_id_counter = copy.copy(self._reg.list_id_counter)

                reader = SparkyPKReader(self._reg.verbose, self._reg.log,
                                        self._reg.representative_model_id,
                                        self._reg.representative_alt_id,
                                        self._reg.mr_atom_name_mapping,
                                        self._reg.cR, self._reg.caC,
                                        self._reg.ccU, self._reg.csStat, self._reg.nefT)
                reader.enforcePeakRowFormat(self._reg.enforce_peak_row_format)
                reader.setInternalMode(self._reg.internal_mode)

                _list_id_counter = copy.copy(self._reg.list_id_counter)

                listener, parser_err_listener, lexer_err_listener =\
                    reader.parse(file_path, self._reg.cifPath,
                                 createSfDict=create_sf_dict, originalFileName=original_file_name,
                                 listIdCounter=self._reg.list_id_counter, reservedListIds=reserved_list_ids,
                                 entryId=self._reg.entry_id,
                                 csLoops=self._reg.lp_data['chem_shift'])

                _content_subtype = listener.getContentSubtype() if listener is not None else None
                if _content_subtype is not None and len(_content_subtype) == 0:
                    _content_subtype = None

                spa_type = 'default'
                if None not in (lexer_err_listener, parser_err_listener, listener)\
                   and ((lexer_err_listener.getMessageList() is None and parser_err_listener.getMessageList() is None)
                        or _content_subtype is not None):
                    skip, spa_type = deal_lexer_or_parser_error(a_pk_format_name, file_name,
                                                                lexer_err_listener, parser_err_listener)
                    if skip and spa_type == 'default':
                        continue

                if spa_type == 'reverse':
                    self._reg.list_id_counter = copy.copy(__list_id_counter)

                    reader = SparkyRPKReader(self._reg.verbose, self._reg.log,
                                             self._reg.representative_model_id,
                                             self._reg.representative_alt_id,
                                             self._reg.mr_atom_name_mapping,
                                             self._reg.cR, self._reg.caC,
                                             self._reg.ccU, self._reg.csStat, self._reg.nefT)
                    reader.enforcePeakRowFormat(self._reg.enforce_peak_row_format)
                    reader.setInternalMode(self._reg.internal_mode)

                    _list_id_counter = copy.copy(self._reg.list_id_counter)

                    listener, parser_err_listener, lexer_err_listener =\
                        reader.parse(file_path, self._reg.cifPath,
                                     createSfDict=create_sf_dict, originalFileName=original_file_name,
                                     listIdCounter=self._reg.list_id_counter, reservedListIds=reserved_list_ids,
                                     entryId=self._reg.entry_id,
                                     csLoops=self._reg.lp_data['chem_shift'])

                    _content_subtype = listener.getContentSubtype() if listener is not None else None
                    if _content_subtype is not None and len(_content_subtype) == 0:
                        _content_subtype = None

                    if None not in (lexer_err_listener, parser_err_listener, listener)\
                       and ((lexer_err_listener.getMessageList() is None and parser_err_listener.getMessageList() is None)
                            or _content_subtype is not None):
                        if deal_lexer_or_parser_error(a_pk_format_name, file_name, lexer_err_listener, parser_err_listener)[0]:
                            continue

                if spa_type == 'no':
                    if self._reg.internal_mode:
                        self._reg.list_id_counter = copy.copy(__list_id_counter)

                        reader = SparkyNPKReader(self._reg.verbose, self._reg.log,
                                                 self._reg.representative_model_id,
                                                 self._reg.representative_alt_id,
                                                 self._reg.mr_atom_name_mapping,
                                                 self._reg.cR, self._reg.caC,
                                                 self._reg.ccU, self._reg.csStat, self._reg.nefT)
                        reader.enforcePeakRowFormat(self._reg.enforce_peak_row_format)
                        reader.setInternalMode(self._reg.internal_mode)

                        _list_id_counter = copy.copy(self._reg.list_id_counter)

                        listener, parser_err_listener, lexer_err_listener =\
                            reader.parse(file_path, self._reg.cifPath,
                                         createSfDict=create_sf_dict, originalFileName=original_file_name,
                                         listIdCounter=self._reg.list_id_counter, reservedListIds=reserved_list_ids,
                                         entryId=self._reg.entry_id,
                                         csLoops=self._reg.lp_data['chem_shift'])

                        _content_subtype = listener.getContentSubtype() if listener is not None else None
                        if _content_subtype is not None and len(_content_subtype) == 0:
                            _content_subtype = None

                        if None not in (lexer_err_listener, parser_err_listener, listener)\
                           and ((lexer_err_listener.getMessageList() is None and parser_err_listener.getMessageList() is None)
                                or _content_subtype is not None):
                            if deal_lexer_or_parser_error(a_pk_format_name, file_name, lexer_err_listener, parser_err_listener)[0]:
                                continue

                    else:
                        warn = "Neither peak height nor peak volume are included in the file. "\
                            "Please re-upload the spectral peak list file."
                        msg_dict = {'file_name': file_name, 'description': warn, 'inheritable': True}

                        self._reg.report.error.appendDescription('format_issue', msg_dict)

                        if self._reg.verbose:
                            self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Error  - {warn}\n")
                        continue

                if listener is not None:
                    reasons = listener.getReasonsForReparsing()

                    if reasons is not None:
                        deal_pea_warn_message_for_lazy_eval(file_name, listener)

                        if spa_type == 'reverse':
                            reader = SparkyRPKReader(self._reg.verbose, self._reg.log,
                                                     self._reg.representative_model_id,
                                                     self._reg.representative_alt_id,
                                                     self._reg.mr_atom_name_mapping,
                                                     self._reg.cR, self._reg.caC,
                                                     self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                                     reasons)
                        elif spa_type == 'default' or not self._reg.internal_mode:
                            reader = SparkyPKReader(self._reg.verbose, self._reg.log,
                                                    self._reg.representative_model_id,
                                                    self._reg.representative_alt_id,
                                                    self._reg.mr_atom_name_mapping,
                                                    self._reg.cR, self._reg.caC,
                                                    self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                                    reasons)
                        else:
                            reader = SparkyNPKReader(self._reg.verbose, self._reg.log,
                                                     self._reg.representative_model_id,
                                                     self._reg.representative_alt_id,
                                                     self._reg.mr_atom_name_mapping,
                                                     self._reg.cR, self._reg.caC,
                                                     self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                                     reasons)
                            reader.setInternalMode(self._reg.internal_mode)

                        reader.enforcePeakRowFormat(self._reg.enforce_peak_row_format)
                        reader.setInternalMode(self._reg.internal_mode)

                        listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                                      createSfDict=create_sf_dict, originalFileName=original_file_name,
                                                      listIdCounter=_list_id_counter, reservedListIds=reserved_list_ids,
                                                      entryId=self._reg.entry_id,
                                                      csLoops=self._reg.lp_data['chem_shift'])

                    deal_pea_warn_message(file_name, listener, ignore_error)

                    poly_seq = listener.getPolymerSequence()
                    if poly_seq is not None:
                        input_source.setItemValue('polymer_sequence', poly_seq)
                        poly_seq_set.append(poly_seq)

                    seq_align = listener.getSequenceAlignment()
                    if seq_align is not None:
                        self._reg.report.sequence_alignment.setItemValue(f'model_poly_seq_vs_{content_subtype}', seq_align)

                    if create_sf_dict:
                        if len(listener.getContentSubtype()) == 0 and not ignore_error:
                            err = f"Failed to validate spectral peak list file (SPARKY) {file_name!r}."

                            self._reg.report.error.appendDescription('internal_error',
                                                                     f"+{self.__class_name__}.validateLegacyPk() "
                                                                     "++ Error  - " + err)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Error  - {err}\n")

                        self._reg.list_id_counter, sf_dict = listener.getSfDict()
                        if sf_dict is not None:
                            for k, v in sf_dict.items():
                                content_subtype = contentSubtypeOf(k[0])
                                if content_subtype not in pk_sf_dict_holder:
                                    pk_sf_dict_holder[content_subtype] = []
                                for sf in v:
                                    if sf not in pk_sf_dict_holder[content_subtype]:
                                        pk_sf_dict_holder[content_subtype].append(sf)

            elif file_type == 'nm-pea-vie':
                __list_id_counter = copy.copy(self._reg.list_id_counter)

                reader = NmrViewPKReader(self._reg.verbose, self._reg.log,
                                         self._reg.representative_model_id,
                                         self._reg.representative_alt_id,
                                         self._reg.mr_atom_name_mapping,
                                         self._reg.cR, self._reg.caC,
                                         self._reg.ccU, self._reg.csStat, self._reg.nefT)
                reader.enforcePeakRowFormat(self._reg.enforce_peak_row_format)
                reader.setInternalMode(self._reg.internal_mode)

                _list_id_counter = copy.copy(self._reg.list_id_counter)

                listener, parser_err_listener, lexer_err_listener =\
                    reader.parse(file_path, self._reg.cifPath,
                                 createSfDict=create_sf_dict, originalFileName=original_file_name,
                                 listIdCounter=self._reg.list_id_counter, reservedListIds=reserved_list_ids,
                                 entryId=self._reg.entry_id,
                                 csLoops=self._reg.lp_data['chem_shift'])

                _content_subtype = listener.getContentSubtype() if listener is not None else None
                if _content_subtype is not None and len(_content_subtype) == 0:
                    _content_subtype = None

                vie_type = 'default'
                if None not in (lexer_err_listener, parser_err_listener, listener)\
                   and ((lexer_err_listener.getMessageList() is None and parser_err_listener.getMessageList() is None)
                        or _content_subtype is not None):
                    skip, vie_type = deal_lexer_or_parser_error(a_pk_format_name, file_name,
                                                                lexer_err_listener, parser_err_listener)
                    if skip and vie_type == 'default':
                        continue

                if vie_type != 'default':
                    self._reg.list_id_counter = copy.copy(__list_id_counter)

                    reader = NmrViewNPKReader(self._reg.verbose, self._reg.log,
                                              self._reg.representative_model_id,
                                              self._reg.representative_alt_id,
                                              self._reg.mr_atom_name_mapping,
                                              self._reg.cR, self._reg.caC,
                                              self._reg.ccU, self._reg.csStat, self._reg.nefT)
                    reader.enforcePeakRowFormat(self._reg.enforce_peak_row_format)
                    reader.setInternalMode(self._reg.internal_mode)

                    _list_id_counter = copy.copy(self._reg.list_id_counter)

                    listener, parser_err_listener, lexer_err_listener =\
                        reader.parse(file_path, self._reg.cifPath,
                                     createSfDict=create_sf_dict, originalFileName=original_file_name,
                                     listIdCounter=self._reg.list_id_counter, reservedListIds=reserved_list_ids,
                                     entryId=self._reg.entry_id,
                                     csLoops=self._reg.lp_data['chem_shift'])

                    _content_subtype = listener.getContentSubtype() if listener is not None else None
                    if _content_subtype is not None and len(_content_subtype) == 0:
                        _content_subtype = None

                    if None not in (lexer_err_listener, parser_err_listener, listener)\
                       and ((lexer_err_listener.getMessageList() is None and parser_err_listener.getMessageList() is None)
                            or _content_subtype is not None):
                        if deal_lexer_or_parser_error(a_pk_format_name, file_name, lexer_err_listener, parser_err_listener)[0]:
                            continue

                if listener is not None:
                    reasons = listener.getReasonsForReparsing()

                    if reasons is not None:
                        deal_pea_warn_message_for_lazy_eval(file_name, listener)

                        if vie_type == 'default':
                            reader = NmrViewPKReader(self._reg.verbose, self._reg.log,
                                                     self._reg.representative_model_id,
                                                     self._reg.representative_alt_id,
                                                     self._reg.mr_atom_name_mapping,
                                                     self._reg.cR, self._reg.caC,
                                                     self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                                     reasons)
                        else:
                            reader = NmrViewNPKReader(self._reg.verbose, self._reg.log,
                                                      self._reg.representative_model_id,
                                                      self._reg.representative_alt_id,
                                                      self._reg.mr_atom_name_mapping,
                                                      self._reg.cR, self._reg.caC,
                                                      self._reg.ccU, self._reg.csStat, self._reg.nefT,
                                                      reasons)

                        reader.enforcePeakRowFormat(self._reg.enforce_peak_row_format)
                        reader.setInternalMode(self._reg.internal_mode)

                        listener, _, _ = reader.parse(file_path, self._reg.cifPath,
                                                      createSfDict=create_sf_dict, originalFileName=original_file_name,
                                                      listIdCounter=_list_id_counter, reservedListIds=reserved_list_ids,
                                                      entryId=self._reg.entry_id,
                                                      csLoops=self._reg.lp_data['chem_shift'])

                    deal_pea_warn_message(file_name, listener, ignore_error)

                    poly_seq = listener.getPolymerSequence()
                    if poly_seq is not None:
                        input_source.setItemValue('polymer_sequence', poly_seq)
                        poly_seq_set.append(poly_seq)

                    seq_align = listener.getSequenceAlignment()
                    if seq_align is not None:
                        self._reg.report.sequence_alignment.setItemValue(f'model_poly_seq_vs_{content_subtype}', seq_align)

                    if create_sf_dict:
                        if len(listener.getContentSubtype()) == 0 and not ignore_error:
                            err = f"Failed to validate spectral peak list file (NMRVIEW) {file_name!r}."

                            self._reg.report.error.appendDescription('internal_error',
                                                                     f"+{self.__class_name__}.validateLegacyPk() "
                                                                     "++ Error  - " + err)

                            if self._reg.verbose:
                                self._reg.log.write(f"+{self.__class_name__}.validateLegacyPk() ++ Error  - {err}\n")

                        self._reg.list_id_counter, sf_dict = listener.getSfDict()
                        if sf_dict is not None:
                            for k, v in sf_dict.items():
                                content_subtype = contentSubtypeOf(k[0])
                                if content_subtype not in pk_sf_dict_holder:
                                    pk_sf_dict_holder[content_subtype] = []
                                for sf in v:
                                    if sf not in pk_sf_dict_holder[content_subtype]:
                                        pk_sf_dict_holder[content_subtype].append(sf)

        if len(self._reg.star_data) > 0 and isinstance(self._reg.star_data[0], pynmrstar.Entry):
            master_entry = self._reg.star_data[0]

            if content_subtype in pk_sf_dict_holder:

                for sf in pk_sf_dict_holder[content_subtype]:

                    cs_list = get_first_sf_tag(sf['saveframe'], 'Chemical_shift_list')

                    if cs_list in EMPTY_VALUE:
                        sf_category = SF_CATEGORIES['nmr-star']['chem_shift']
                        cs_sf_list = master_entry.get_saveframes_by_category(sf_category)
                        if len(cs_sf_list) == 1:
                            set_sf_tag(sf['saveframe'], 'Chemical_shift_list', get_first_sf_tag(cs_sf_list[0], 'Sf_framecode'))

                    # prevent duplication of spectral peak list
                    data_file_name = get_first_sf_tag(sf['saveframe'], 'Data_file_name')
                    if data_file_name not in EMPTY_VALUE\
                       and len(master_entry.get_saveframes_by_tag_and_value('Data_file_name', data_file_name)) > 0:
                        continue

                    try:
                        master_entry.add_saveframe(sf['saveframe'])
                    except ValueError:
                        pass

                self._reg.pk_sf_holder = pk_sf_dict_holder['spectral_peak']

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

            file_type = 'nm-pea-any'

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

            self._reg.report.sequence_alignment.setItemValue(f'model_poly_seq_vs_{content_subtype}', seq_align)

        if proc_nmr_ext_poly_seq and len(self._reg.nmr_ext_poly_seq) > 0:
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
