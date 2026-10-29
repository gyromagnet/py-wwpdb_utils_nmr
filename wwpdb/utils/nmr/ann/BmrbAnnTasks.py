# File: BmrbAnnTasks.py
# Date: 24-Dec-2023
#
# Updates:
# 09-Jul-2026  M. Yokochi - implement BMRB's data provenance check in standalone NMR data conversion service (DAOTHER-9785)
# 09-Jul-2026  M. Yokochi - update/generate deposit_data_file saveframe in standalone NMR data conversion service (DAOTHER-9785)
##
""" Wrapper class for BMRB annotation tasks.
    @author: Masashi Yokochi
"""
__docformat__ = "restructuredtext en"
__author__ = "Masashi Yokochi"
__email__ = "yokochi@protein.osaka-u.ac.jp"
__license__ = "Apache License 2.0"
__version__ = "1.2.0"

import copy
import itertools
import re
from operator import itemgetter
from typing import Optional

import pynmrstar

try:
    from wwpdb.utils.nmr.NmrDpConstant import (EMPTY_VALUE,
                                               CNV_ID_PAT,
                                               STD_MON_DICT,
                                               ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS,
                                               ISOTOPE_NAMES_OF_NMR_OBS_NUCS,
                                               ALLOWED_AMBIGUITY_CODES,
                                               ALLOWED_ISOTOPE_NUMBERS,
                                               WELL_KNOWN_ISOTOPE_NUMBERS,
                                               CS_UNCERT_MAX,
                                               MAX_ROWS_TO_CHECK_SPECTRAL_PEAK_IDENTITY,
                                               PROTEIN_RELATED_WORDS,
                                               DNA_RELATED_WORDS,
                                               RNA_RELATED_WORDS,
                                               ALLOWED_THIOL_STATES,
                                               SF_CATEGORIES)
    from wwpdb.utils.nmr.NmrDpRegistry import NmrDpRegistry
    from wwpdb.utils.nmr.AlignUtil import (getOneLetterCodeCan,
                                           getScoreOfSeqAlign,
                                           getRestraintFormatName,
                                           getChemShiftFormatName)
    from wwpdb.utils.nmr.CifToNmrStar import (get_first_sf_tag,
                                              set_sf_tag,
                                              has_key_value)
    from wwpdb.utils.nmr.mr.ParserListenerUtil import (getMaxEffDigits,
                                                       roundString,
                                                       retrieveOriginalFileName,
                                                       getPdbxNmrSoftwareName)
except ImportError:
    from nmr.NmrDpConstant import (EMPTY_VALUE,
                                   CNV_ID_PAT,
                                   STD_MON_DICT,
                                   ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS,
                                   ISOTOPE_NAMES_OF_NMR_OBS_NUCS,
                                   ALLOWED_AMBIGUITY_CODES,
                                   ALLOWED_ISOTOPE_NUMBERS,
                                   WELL_KNOWN_ISOTOPE_NUMBERS,
                                   CS_UNCERT_MAX,
                                   MAX_ROWS_TO_CHECK_SPECTRAL_PEAK_IDENTITY,
                                   PROTEIN_RELATED_WORDS,
                                   DNA_RELATED_WORDS,
                                   RNA_RELATED_WORDS,
                                   ALLOWED_THIOL_STATES,
                                   SF_CATEGORIES)
    from nmr.NmrDpRegistry import NmrDpRegistry
    from nmr.AlignUtil import (getOneLetterCodeCan,
                               getScoreOfSeqAlign,
                               getRestraintFormatName,
                               getChemShiftFormatName)
    from nmr.CifToNmrStar import (get_first_sf_tag,
                                  set_sf_tag,
                                  has_key_value)
    from nmr.mr.ParserListenerUtil import (getMaxEffDigits,
                                           roundString,
                                           retrieveOriginalFileName,
                                           getPdbxNmrSoftwareName)


DEFAULT_SF_LABEL_TAGS = ('_Assigned_chem_shift_list.Sample_condition_list_label',
                         '_Assigned_chem_shift_list.Chem_shift_reference_label',
                         '_Spectral_peak_list.Sample_label',
                         '_Spectral_peak_list.Sample_condition_list_label',
                         '_Conformer_stat_list.Conf_family_coord_set_label',
                         '_Conformer_stat_list.Representative_conformer_label',
                         '_Conformer_family_coord_set.Sample_condition_list_label')


def annotate_assembly_experimental_data(loop: pynmrstar.Loop, reported: set[int]) -> None:
    """Mark assemblies with known experimental data; leave other states unknown.

    Extracted input tables may omit the annotation column. A missing column must
    not crash conversion or be interpreted as a set of negative observations.
    """
    if 'ID' not in loop.tags:
        return
    if 'Experimental_data_reported' not in loop.tags:
        loop.add_tag('Experimental_data_reported')
        for row in loop.data:
            row.append(None)
    id_column = loop.tags.index('ID')
    data_column = loop.tags.index('Experimental_data_reported')
    for row in loop.data:
        identifier = row[id_column]
        if isinstance(identifier, str) and identifier.isdigit() and int(identifier) in reported:
            row[data_column] = 'yes'


class BmrbAnnTasks:
    """ Wrapper class for BMRB annotation tasks.
    """
    __slots__ = ('__class_name__',
                 '__version__',
                 '__reg',
                 '__derivedEntryId',
                 '__derivedEntryTitle',
                 '__defSfLabelTag',
                 '__update_related_entries',
                 '__secret_key',
                 '__service_host',
                 '__dep_sys_name')

    def __init__(self, registry: NmrDpRegistry) -> None:
        self.__class_name__ = self.__class__.__name__
        self.__version__ = __version__

        self.__reg = registry

        # provenance information to be included in _Related_entries loop if it does not exist
        self.__derivedEntryId = None
        self.__derivedEntryTitle = None

        self.__update_related_entries = True
        if has_key_value(self.__reg.inputParamDict, 'update_related_entries'):
            if isinstance(self.__reg.inputParamDict['update_related_entries'], bool):
                self.__update_related_entries = self.__reg.inputParamDict['update_related_entries']

        # parameters for BMRB's data provenance check in standalone NMR data conversion service (DAOTHER-9785)
        self.__secret_key = None
        self.__service_host = None
        self.__dep_sys_name = 'unknown'
        if self.__reg.conversion_server and self.__reg.entry_id is not None and CNV_ID_PAT.match(self.__reg.entry_id):
            if has_key_value(self.__reg.inputParamDict, 'secret_key')\
               and has_key_value(self.__reg.inputParamDict, 'service_host'):
                secret_key = self.__reg.inputParamDict['secret_key']
                service_host = self.__reg.inputParamDict['service_host']
                if isinstance(secret_key, str) and secret_key not in EMPTY_VALUE\
                   and isinstance(service_host, str) and service_host not in EMPTY_VALUE:
                    self.__secret_key = secret_key
                    self.__service_host = service_host
                    if has_key_value(self.__reg.inputParamDict, 'dep_sys_name'):
                        dep_sys_name = self.__reg.inputParamDict['dep_sys_name']
                        if isinstance(dep_sys_name, str) and dep_sys_name not in EMPTY_VALUE:
                            self.__dep_sys_name = dep_sys_name

    def setProvenanceInfo(self, derivedEntryId: Optional[str], derivedEntryTitle: Optional[str]) -> None:
        """ Set provenance information.
        """

        self.__derivedEntryId = derivedEntryId
        self.__derivedEntryTitle = derivedEntryTitle

    def perform(self, master_entry: pynmrstar.Entry) -> bool:
        """ Perform a series of BMRB annotation tasks.
        """

        if not isinstance(master_entry, pynmrstar.Entry):
            return False

        def eff_digits(value: str):
            """ Return effective digits.
            """

            return len(value) - 1 - value.index('.') if '.' in value else 0

        def not_title(title: str):
            if title in EMPTY_VALUE:
                return True
            title = title.lower()
            return title.startswith('unidentified') or title.startswith('unclassified')\
                or title.startswith('unknown') or title.startswith('unspecified')\
                or title.startswith('none') or title.startswith('not ')\
                or title.startswith('n/a')\
                or title.startswith('na ') or title in ('na', 'no')\
                or title.startswith('subunit')

        sf_category = 'entry_interview'

        int_sf = None

        if sf_category in self.__reg.sf_category_list:
            int_sf = master_entry.get_saveframes_by_category(sf_category)[0]

        sf_category = 'entry_information'

        if sf_category in self.__reg.sf_category_list:
            ent_sf = master_entry.get_saveframes_by_category(sf_category)[0]

            lp_category = '_Contact_person'

            try:

                lp = ent_sf.get_loop(lp_category)

                tags = ['Email_address', 'Given_name', 'Family_name',
                        'Department_and_institution', 'Country', 'State_province', 'City',
                        'Postal_code', 'Role', 'Organization_type']

                if set(tags) & set(lp.tags) == set(tags):
                    contact = lp.get_tag(tags)
                    len_contact = len(contact)

                    dup_idx = set()
                    for idx, row in enumerate(contact):
                        for _idx in range(idx + 1, len_contact):
                            if contact[_idx] == row:
                                dup_idx.add(_idx)

                    if len(dup_idx) > 0:
                        dup_idx = sorted(dup_idx, reverse=True)

                        for _idx in dup_idx:
                            del lp.data[_idx]

                if 'ID' in lp.tags:
                    lp.renumber_rows('ID')

            except KeyError:
                pass

            lp_category = '_Entry_author'

            try:

                lp = ent_sf.get_loop(lp_category)

                tags = ['Given_name', 'Family_name', 'First_initial', 'Middle_initials', 'Family_title', 'ORCID']

                if set(tags) & set(lp.tags) == set(tags):
                    author = lp.get_tag(tags)
                    len_author = len(author)

                    dup_idx = set()
                    for idx, row in enumerate(author):
                        for _idx in range(idx + 1, len_author):
                            if author[_idx] == row:
                                dup_idx.add(_idx)

                    if len(dup_idx) > 0:
                        dup_idx = sorted(dup_idx, reverse=True)

                        for _idx in dup_idx:
                            del lp.data[_idx]

                if 'Ordinal' in lp.tags:
                    lp.renumber_rows('Ordinal')

            except KeyError:
                pass

            if not self.__reg.internal_mode and self.__derivedEntryId not in EMPTY_VALUE and self.__update_related_entries:

                lp_category = '_Related_entries'

                try:

                    lp = ent_sf.get_loop(lp_category)

                    tags = ['Database_name', 'Database_accession_code']

                    if set(tags) & set(lp.tags) == set(tags):
                        related_entries = lp.get_tag(tags)

                        has_provenance = False
                        for idx, row in enumerate(related_entries):
                            if row == ['BMRB', self.__derivedEntryId]:
                                has_provenance = True
                                if self.__derivedEntryTitle not in EMPTY_VALUE:
                                    lp.data[idx][lp.tags.index('Relationship')] = self.__derivedEntryTitle
                                break

                        if not has_provenance:
                            row = [None] * len(lp.tags)
                            row[lp.tags.index('Database_name')] = 'BMRB'
                            row[lp.tags.index('Database_accession_code')] = self.__derivedEntryId
                            if self.__derivedEntryTitle not in EMPTY_VALUE:
                                row[lp.tags.index('Relationship')] = self.__derivedEntryTitle
                            row[lp.tags.index('Entry_ID')] = self.__reg.entry_id

                            lp.add_data(row)

                            lp.sort_rows('Database_accession_code')

                except KeyError:
                    # items = ['Database_name', 'Database_accession_code', 'Relationship', 'Entry_ID']
                    #
                    # lp = pynmrstar.Loop.from_scratch(lp_category)
                    #
                    # tags = [f'{lp_category}.{item}' for item in items]
                    #
                    # lp.add_tag(tags)
                    #
                    # lp.add_data(['BMRB', self.__derivedEntryId, self.__derivedEntryTitle, self.__reg.entry_id])
                    #
                    # ent_sf.add_loop(lp)
                    pass

            # write signature for BMRB's data provenance check in standalone NMR data conversion service (DAOTHER-9785)
            if self.__reg.conversion_server and None not in (self.__secret_key, self.__service_host, self.__reg.entry_id):

                signed_by = get_first_sf_tag(ent_sf, 'Signed_by')
                signature = get_first_sf_tag(ent_sf, 'Signature')

                if signed_by in EMPTY_VALUE or signature in EMPTY_VALUE:

                    try:

                        from itsdangerous import URLSafeSerializer  # pylint: disable=import-outside-toplevel

                        serializer = URLSafeSerializer(self.__secret_key, self.__service_host)

                        signature = serializer.dumps({'converion_id': self.__reg.entry_id, 'dep_sys_name': self.__dep_sys_name})

                        set_sf_tag(ent_sf, 'Signed_by', self.__service_host)
                        set_sf_tag(ent_sf, 'Signature', signature)

                    except ImportError:
                        pass

        # generate/update deposited_data_files saveframe in standalone NMR data conversion service (DAOTHER-9785)
        if self.__reg.conversion_server and self.__reg.op in ('nmr-cs-mr-merge', 'nmr-str-replace-cs',
                                                              'nmr-nef2str-deposit', 'nmr-str2str-deposit'):

            sf_category = 'deposited_data_files'
            lp_category = '_Upload_data'

            def gen_upload_data_loop(list_id):

                loop = pynmrstar.Loop.from_scratch(lp_category)

                items = ['Data_file_ID', 'Data_file_name', 'Data_file_content_type', 'Data_file_Sf_category',
                         'Data_file_syntax', 'Entry_ID', 'Deposoted_data_files_ID']

                tags = [f'{lp_category}.{item}' for item in items]

                loop.add_tag(tags)

                for fileListId, input_source in enumerate(self.__reg.report.input_sources):

                    input_source_dic = input_source.get()

                    file_name = input_source_dic['file_name']
                    file_type = input_source_dic['file_type']
                    original_file_name = input_source_dic['original_file_name']
                    content_type = input_source_dic['content_type']
                    content_subtype = input_source_dic['content_subtype']

                    row = [None] * len(loop.tags)

                    row[0] = str(fileListId + 1)
                    row[1] = file_name if original_file_name in EMPTY_VALUE else original_file_name
                    row[2] = content_type
                    if content_subtype is not None and len(content_subtype) > 0:
                        sf_cat_list = []
                        if content_type == 'model':
                            sf_cat_list.append('conformer_family_coord_set')
                        # elif content_type == 'nmr-chemical-shifts':
                        #     sf_cat_list.append('assigned_chemical_shifts')
                        elif not file_type.startswith('nm-csp'):
                            for k in content_subtype:
                                if k in SF_CATEGORIES['nmr-star']:
                                    sf_cat_list.append(SF_CATEGORIES['nmr-star'][k])
                        else:
                            sf_cat_list.append('chem_shift_perturbation')
                        if len(sf_cat_list) > 0:
                            row[3] = ','.join(sf_cat_list)
                    if file_type in ('nmr-star', 'nm-shi'):
                        row[4] = 'NMR-STAR V3'
                    elif file_type == 'nef':
                        row[4] = 'NEF'
                    elif file_type == 'pdbx':
                        row[4] = 'PDBx/mmCIF'
                    elif file_type.startswith('nm-res')\
                            or file_type.startswith('nm-aux')\
                            or file_type.startswith('nm-pea'):
                        syntax = getRestraintFormatName(file_type)
                        if not syntax.startswith('other'):
                            syntax = getPdbxNmrSoftwareName(syntax.split()[0])
                        row[4] = syntax
                    elif file_type.startswith('nm-shi') or file_type.startswith('nm-csp'):
                        if file_type == 'nm-shi-st2':
                            syntax = 'NMR-STAR V2'
                        else:
                            syntax = getChemShiftFormatName(file_type)
                            if not syntax.startswith('other'):
                                syntax = getPdbxNmrSoftwareName(syntax.split()[0])
                        row[4] = syntax
                    row[5] = self.__reg.entry_id
                    row[6] = list_id

                    loop.add_data(row)

                return loop

            if sf_category in self.__reg.sf_category_list and self.__reg.op == 'nmr-str-replace-cs':

                for sf in master_entry.get_saveframes_by_category(sf_category):

                    list_id = get_first_sf_tag(sf, 'ID')

                    try:

                        lp = sf.get_loop(lp_category)

                        id_col = lp.tags.index('Data_file_ID') if 'Data_file_ID' in lp.tags else -1
                        name_col = lp.tags.index('Data_file_name') if 'Data_file_name' in lp.tags else -1
                        type_col = lp.tags.index('Data_file_content_type') if 'Data_file_content_type' in lp.tags else -1
                        sf_cat_col = lp.tags.index('Data_file_Sf_category') if 'Data_file_Sf_category' in lp.tags else -1
                        syntax_col = lp.tags.index('Data_file_syntax') if 'Data_file_syntax' in lp.tags else -1
                        entry_id_col = lp.tags.index('Entry_ID') if 'Entry_ID' in lp.tags else -1
                        list_id_col = lp.tags.index('Deposoted_data_files_ID') if 'Deposited_data_files_ID' in lp.tags else -1

                        data_file_names = lp.get_tag(['Data_file_name'])

                        has_model = False
                        for row in lp:
                            if row[syntax_col] in ('PDB', 'mmCIF', 'PDBx/mmCIF'):
                                has_model = True
                                break

                        is_emply = len(lp) > (1 if has_model else 0)

                        file_id = len(lp) + 1

                        for fileListId, input_source in enumerate(self.__reg.report.input_sources):

                            input_source_dic = input_source.get()

                            file_name = input_source_dic['file_name']
                            file_type = input_source_dic['file_type']
                            original_file_name = input_source_dic['original_file_name']
                            content_type = input_source_dic['content_type']
                            content_subtype = input_source_dic['content_subtype']

                            if not is_emply and fileListId == 0:
                                continue

                            if has_model and file_type == 'pdbx':
                                continue

                            row = [None] * len(lp.tags)

                            if id_col != -1:
                                row[id_col] = str(file_id)
                            if name_col != -1:
                                row[name_col] = file_name if original_file_name in EMPTY_VALUE else original_file_name
                                if row[name_col] in data_file_names:
                                    continue
                            if type_col != -1:
                                row[type_col] = content_type
                            if sf_cat_col != -1 and content_subtype is not None and len(content_subtype) > 0:
                                sf_cat_list = []
                                if content_type == 'model':
                                    sf_cat_list.append('conformer_family_coord_set')
                                # elif content_type == 'nmr-chemical-shifts':
                                #     sf_cat_list.append('assigned_chemical_shifts')
                                elif not file_type.startswith('nm-csp'):
                                    for k in content_subtype:
                                        if k in SF_CATEGORIES['nmr-star']:
                                            sf_cat_list.append(SF_CATEGORIES['nmr-star'][k])
                                else:
                                    sf_cat_list.append('chem_shift_perturbation')
                                if len(sf_cat_list) > 0:
                                    row[sf_cat_col] = ','.join(sf_cat_list)
                            if syntax_col != -1:
                                if file_type in ('nmr-star', 'nm-shi'):
                                    row[syntax_col] = 'NMR-STAR V3'
                                elif file_type == 'nef':
                                    row[syntax_col] = 'NEF'
                                elif file_type == 'pdbx':
                                    row[syntax_col] = 'PDBx/mmCIF'
                                elif file_type.startswith('nm-res')\
                                        or file_type.startswith('nm-aux')\
                                        or file_type.startswith('nm-pea'):
                                    syntax = getRestraintFormatName(file_type)
                                    if not syntax.startswith('other'):
                                        syntax = getPdbxNmrSoftwareName(syntax.split()[0])
                                    row[syntax_col] = syntax
                                elif file_type.startswith('nm-shi') or file_type.startswith('nm-csp'):
                                    if file_type == 'nm-shi-st2':
                                        syntax = 'NMR-STAR V2'
                                    else:
                                        syntax = getChemShiftFormatName(file_type)
                                        if not syntax.startswith('other'):
                                            syntax = getPdbxNmrSoftwareName(syntax.split()[0])
                                    row[syntax_col] = syntax
                            if entry_id_col != -1:
                                row[entry_id_col] = self.__reg.entry_id
                            if list_id_col != -1:
                                row[list_id_col] = list_id

                            lp.add_data(row)

                            file_id += 1

                    except KeyError:
                        sf.add_loop(gen_upload_data_loop(list_id))

                    break

            else:

                if sf_category in self.__reg.sf_category_list:
                    for sf in master_entry.get_saveframes_by_category(sf_category):
                        master_entry.remove_saveframe(sf.name)

                list_id = '1'

                sf = pynmrstar.Saveframe.from_scratch(sf_category, '_Deposited_data_files')

                sf.add_tag('Sf_category', sf_category)
                sf.add_tag('Sf_framecode', sf_category)
                sf.add_tag('Entry_ID', self.__reg.entry_id)
                sf.add_tag('ID', list_id)
                sf.add_tag('Atomic_coordinate_file_name', '.')
                sf.add_tag('Atomic_coordinate_file_syntax', '.')
                sf.add_tag('Constraint_file_name', '.')
                sf.add_tag('Constraint_file_syntax', '.')
                sf.add_tag('Precheck_flag', '.')
                sf.add_tag('Validate_flag', '.')

                sf.add_loop(gen_upload_data_loop(list_id))

                master_entry.add_saveframe(sf)

        # strip citation author names

        sf_category = 'citations'

        if sf_category in self.__reg.sf_category_list:

            for sf in master_entry.get_saveframes_by_category(sf_category):

                lp_category = '_Citation_author'

                try:

                    lp = sf.get_loop(lp_category)

                    tags = ['Given_name', 'Family_name', 'First_initial', 'Middle_initials', 'Family_title']

                    dat = lp.get_tag(tags)

                    for idx, row in enumerate(dat):
                        for col in range(5):
                            if row[col].startswith(' ') or row[col].endswith(' '):
                                lp.data[idx][lp.tags.index(tags[col])] = row[col].strip()

                except KeyError:
                    pass

        sf_category = 'experiment_list'

        if sf_category in self.__reg.sf_category_list:
            for sf in master_entry.get_saveframes_by_category(sf_category):

                lp_category = '_Experiment'

                try:

                    lp = sf.get_loop(lp_category)

                    tags = ['Sample_ID', 'Sample_label',
                            'Sample_condition_list_ID', 'Sample_condition_list_label',
                            'NMR_spectrometer_ID', 'NMR_spectrometer_label']

                    if set(tags) & set(lp.tags) == set(tags):
                        exp_list = lp.get_tag(tags)

                        def sync_exp_lp_and_sf(row, idx, id_col, label_col, parent_sf_tag_prefix):
                            if (row[id_col] in EMPTY_VALUE and row[label_col] not in EMPTY_VALUE)\
                               or (row[id_col] not in EMPTY_VALUE and row[label_col] in EMPTY_VALUE):
                                list_id = sf_framecode = None
                                if row[id_col] in EMPTY_VALUE:
                                    _sf_framecode = row[label_col].lstrip('$')
                                    sf_framecode = _sf_framecode.replace(' ', '_')
                                    try:
                                        parent_sf = master_entry.get_saveframe_by_name(sf_framecode)
                                        list_id = get_first_sf_tag(parent_sf, 'ID')
                                    except KeyError:
                                        pass
                                    if list_id in EMPTY_VALUE and sf_framecode.split('_')[-1].isdigit():
                                        list_id = sf_framecode.split('_')[-1]
                                        if _sf_framecode != sf_framecode:
                                            sf_framecode = f'{parent_sf_tag_prefix[1:]}_{list_id}'
                                            for _list_id, parent_sf\
                                                    in enumerate(master_entry.get_saveframes_by_category(parent_sf_tag_prefix[1:]),
                                                                 start=1):
                                                if str(_list_id) == list_id:
                                                    set_sf_tag(parent_sf, 'Sf_framecode', sf_framecode)
                                                    set_sf_tag(parent_sf, 'ID', list_id)
                                                    name = get_first_sf_tag(parent_sf, 'Name')
                                                    if len(name) == 0:
                                                        set_sf_tag(parent_sf, 'Name', _sf_framecode)
                                else:
                                    list_id = row[id_col]
                                    if isinstance(list_id, int):
                                        list_id = str(list_id)
                                    for parent_sf in master_entry.get_saveframes_by_tag_and_value(f'{parent_sf_tag_prefix}.ID', list_id):  # noqa: E501, pylint: disable=line-too-long
                                        sf_framecode = get_first_sf_tag(parent_sf, 'Sf_framecode').replace(' ', '_')
                                        break
                                if None not in (list_id, sf_framecode):
                                    lp.data[idx][lp.tags.index(tags[id_col])] = list_id
                                    lp.data[idx][lp.tags.index(tags[label_col])] = f'${sf_framecode}'

                        for idx, row in enumerate(exp_list):
                            sync_exp_lp_and_sf(row, idx, 0, 1, '_Sample')
                            sync_exp_lp_and_sf(row, idx, 2, 3, '_Sample_condition_list')
                            sync_exp_lp_and_sf(row, idx, 4, 5, '_NMR_spectrometer')

                except KeyError:
                    pass

        # section 11: assigned chemical shifts

        isotope_nums, cs_ref_sf_framecode, smpl_cond_sf_framecode, isotope_nums_per_entity = {}, {}, {}, {}

        sf_category = 'sample'

        if sf_category in self.__reg.sf_category_list:
            for sf in master_entry.get_saveframes_by_category(sf_category):
                smpl_id = get_first_sf_tag(sf, 'ID')
                if not isinstance(smpl_id, int):
                    if len(smpl_id) == 0 or smpl_id in EMPTY_VALUE or not smpl_id.isdigit():
                        smpl_id = 1
                        set_sf_tag(sf, 'ID', smpl_id)
                    else:
                        smpl_id = int(smpl_id)
                    self.__reg.c2S.set_local_sf_id(sf, smpl_id)

        sf_category = 'chem_shift_reference'

        if sf_category in self.__reg.sf_category_list:
            for sf in master_entry.get_saveframes_by_category(sf_category):
                cs_ref_id = get_first_sf_tag(sf, 'ID')
                if not isinstance(cs_ref_id, int):
                    if len(cs_ref_id) == 0 or cs_ref_id in EMPTY_VALUE or not cs_ref_id.isdigit():
                        cs_ref_id = 1
                        set_sf_tag(sf, 'ID', cs_ref_id)
                    else:
                        cs_ref_id = int(cs_ref_id)
                    self.__reg.c2S.set_local_sf_id(sf, cs_ref_id)

                isotope_nums[cs_ref_id] = set()
                cs_ref_sf_framecode[cs_ref_id] = get_first_sf_tag(sf, 'Sf_framecode')

        sf_category = 'sample_conditions'

        if sf_category in self.__reg.sf_category_list:
            for sf in master_entry.get_saveframes_by_category(sf_category):
                smpl_cond_list_id = get_first_sf_tag(sf, 'ID')
                if not isinstance(smpl_cond_list_id, int):
                    if len(smpl_cond_list_id) == 0 or smpl_cond_list_id in EMPTY_VALUE or not smpl_cond_list_id.isdigit():
                        smpl_cond_list_id = 1
                        set_sf_tag(sf, 'ID', smpl_cond_list_id)
                    else:
                        smpl_cond_list_id = int(smpl_cond_list_id)
                    self.__reg.c2S.set_local_sf_id(sf, smpl_cond_list_id)

                smpl_cond_sf_framecode[smpl_cond_list_id] = get_first_sf_tag(sf, 'Sf_framecode')

        empty_val_err_tags, zero_val_err_tags = [], []
        zero_shift_val_err = ''
        label_to_auth_seq = {}
        ent_asym_id_with_exptl_data = set()
        has_cs_lp = True

        sf_category = 'assigned_chemical_shifts'

        if sf_category in self.__reg.sf_category_list:
            for list_id, sf in enumerate(master_entry.get_saveframes_by_category(sf_category), start=1):
                self.__reg.c2S.set_local_sf_id(sf, list_id)
                cs_ref_id = get_first_sf_tag(sf, 'Chem_shift_reference_ID')
                if not isinstance(cs_ref_id, int):
                    if len(cs_ref_id) == 0 or cs_ref_id in EMPTY_VALUE or not cs_ref_id.isdigit():
                        cs_ref_id = 1
                        set_sf_tag(sf, 'Chem_shift_reference_ID', cs_ref_id)
                    else:
                        cs_ref_id = int(cs_ref_id)
                if cs_ref_id in cs_ref_sf_framecode:
                    set_sf_tag(sf, 'Chem_shift_reference_label', f'${cs_ref_sf_framecode[cs_ref_id]}')
                else:
                    isotope_nums[cs_ref_id] = set()
                    if len(cs_ref_sf_framecode) > 0:
                        cs_ref_id = min(list(cs_ref_sf_framecode.keys()))
                        set_sf_tag(sf, 'Chem_shift_reference_ID', cs_ref_id)
                        set_sf_tag(sf, 'Chem_shift_reference_label', f'${cs_ref_sf_framecode[cs_ref_id]}')
                    else:
                        cs_ref_id = 1
                        set_sf_tag(sf, 'Chem_shift_reference_ID', cs_ref_id)

                smpl_cond_list_id = get_first_sf_tag(sf, 'Sample_condition_list_ID')
                if not isinstance(smpl_cond_list_id, int):
                    if len(smpl_cond_list_id) == 0 or smpl_cond_list_id in EMPTY_VALUE or not smpl_cond_list_id.isdigit():
                        smpl_cond_list_id = 1
                        set_sf_tag(sf, 'Sample_condition_list_ID', smpl_cond_list_id)
                    else:
                        smpl_cond_list_id = int(smpl_cond_list_id)
                if smpl_cond_list_id in smpl_cond_sf_framecode:
                    set_sf_tag(sf, 'Sample_condition_list_label', f'${smpl_cond_sf_framecode[smpl_cond_list_id]}')
                else:
                    if len(smpl_cond_sf_framecode) > 0:
                        smpl_cond_list_id = min(list(smpl_cond_sf_framecode.keys()))
                        set_sf_tag(sf, 'Sample_condition_list_ID', smpl_cond_list_id)
                        set_sf_tag(sf, 'Sample_condition_list_label', f'${smpl_cond_sf_framecode[smpl_cond_list_id]}')
                    else:
                        smpl_cond_list_id = 1
                        set_sf_tag(sf, 'Sample_condition_list_ID', smpl_cond_list_id)

                exp_list_sf_category = 'experiment_list'

                if exp_list_sf_category in self.__reg.sf_category_list:
                    exp_list_sf = master_entry.get_saveframes_by_category(exp_list_sf_category)[0]

                    exp_lp_category = '_Experiment'

                    try:

                        exp_lp = exp_list_sf.get_loop(exp_lp_category)
                        exp_lp.sort_rows('ID')

                    except KeyError:
                        exp_lp = None

                    if exp_lp is not None:

                        exp_tags = ['ID', 'Name', 'Sample_ID', 'Sample_label', 'Sample_state']

                        if set(exp_tags) & set(exp_lp.tags) == set(exp_tags):
                            exp_list = exp_lp.get_tag(exp_tags)

                            if len(exp_list) > 0:

                                lp_category = '_Chem_shift_experiment'
                                lp = None

                                try:

                                    lp = sf.get_loop(lp_category)
                                    lp.sort_rows('Experiment_ID')

                                    tags = ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state']

                                    if set(tags) & set(lp.tags) == set(tags):
                                        exp_id_col = lp.tags.index('Experiment_ID')
                                        exp_name_col = lp.tags.index('Experiment_name')
                                        sample_id_col = lp.tags.index('Sample_ID')
                                        sample_label_col = lp.tags.index('Sample_label')
                                        sample_state_col = lp.tags.index('Sample_state')

                                        cs_exp_list = lp.get_tag(tags)
                                        reserved_ids, duplicated_idxs, reserved_names = [], [], []

                                        for idx, cs_exp in enumerate(cs_exp_list):
                                            if cs_exp[0] not in EMPTY_VALUE:
                                                if cs_exp[0] not in reserved_ids:
                                                    reserved_ids.append(cs_exp[0])
                                                else:
                                                    duplicated_idxs.append(idx)
                                                exp = next((exp for exp in exp_list if exp[0] == cs_exp[0]), None)
                                                if exp is not None:
                                                    if cs_exp[2:5] != exp[2:5]:
                                                        lp.data[idx][exp_name_col] = exp[1]
                                                        lp.data[idx][sample_id_col] = exp[2]
                                                        lp.data[idx][sample_label_col] = exp[3]
                                                        lp.data[idx][sample_state_col] = exp[4]
                                                    reserved_names.append(exp[1])
                                            else:
                                                exp = next((exp for exp in exp_list if exp[1] not in reserved_names), None)
                                                if exp is not None:
                                                    lp.data[idx][exp_id_col] = exp_list.index(exp) + 1
                                                    lp.data[idx][exp_name_col] = exp[1]
                                                    lp.data[idx][sample_id_col] = exp[2]
                                                    lp.data[idx][sample_label_col] = exp[3]
                                                    lp.data[idx][sample_state_col] = exp[4]
                                                    reserved_names.append(exp[1])

                                        if len(duplicated_idxs) > 0:
                                            for idx in reversed(duplicated_idxs):
                                                del lp.data[idx]

                                except (KeyError, TypeError, ValueError):

                                    if lp is not None:
                                        del sf[lp]

                                    items = ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state',
                                             'Entry_ID', 'Assigned_chem_shift_list_ID']

                                    cs_list_id = get_first_sf_tag(sf, 'ID')

                                    lp = pynmrstar.Loop.from_scratch(lp_category)

                                    tags = [f'{lp_category}.{item}' for item in items]

                                    lp.add_tag(tags)

                                    for exp in exp_list:
                                        row = exp
                                        row.extend([self.__reg.entry_id, cs_list_id])
                                        lp.add_data(row)

                                    lp.sort_rows('Experiment_ID')

                                    sf.add_loop(lp)

                try:

                    lp_category = '_Atom_chem_shift'

                    lp = sf.get_loop(lp_category)

                    if self.__reg.report.getInputSourceIdOfCoord() < 0:

                        try:

                            sf.get_loop('_Ambiguous_atom_chem_shift')

                        except KeyError:

                            if 'Entity_assembly_ID' in lp.tags:
                                chain_id_col = lp.tags.index('Entity_assembly_ID')

                                if not any(True for _row in lp
                                           if _row[chain_id_col] in EMPTY_VALUE
                                           or (isinstance(_row[chain_id_col], str) and not _row[chain_id_col].isdigit())):
                                    try:
                                        lp.sort_rows(['Atom_ID', 'Atom_isotope_number', 'Comp_index_ID', 'Entity_assembly_ID'])
                                    except (TypeError, ValueError):
                                        pass

                            if 'Index_ID' in lp.tags:
                                lp.renumber_rows('Index_ID')

                    tags = ['Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Auth_seq_ID']

                    if set(tags) & set(lp.tags) == set(tags):

                        dat = lp.get_tag(tags)

                        for row in dat:
                            if isinstance(row[0], str) and row[0] not in EMPTY_VALUE and row[0].isdigit():
                                ent_asym_id_with_exptl_data.add(int(row[0]))
                            if row[1] not in EMPTY_VALUE and row[2] not in EMPTY_VALUE and row[3] not in EMPTY_VALUE:
                                try:
                                    seq_key = (int(row[1]), int(row[2]))
                                    if seq_key not in label_to_auth_seq:
                                        label_to_auth_seq[seq_key] = int(row[3])
                                except ValueError:
                                    continue

                    tags = ['Atom_isotope_number', 'Val_err', 'Ambiguity_code', 'Entity_ID']

                    if set(tags) & set(lp.tags) == set(tags):

                        dat = lp.get_tag(tags)

                        val_err_col = lp.tags.index('Val_err')
                        ambig_code_col = lp.tags.index('Ambiguity_code')

                        isotope_nums_in_loop = set()
                        isotope_nums_with_empty_val_err = set()
                        isotope_nums_with_zero_val_err = set()
                        isotope_nums_with_non_zero_val_err = set()
                        empty_ambig_code = False

                        for idx, row in enumerate(dat):
                            try:
                                isotope_number = int(row[0])
                                isotope_nums[cs_ref_id].add(isotope_number)
                                isotope_nums_in_loop.add(isotope_number)
                                if row[1] in EMPTY_VALUE and isotope_number in ALLOWED_ISOTOPE_NUMBERS:
                                    isotope_nums_with_empty_val_err.add(isotope_number)
                                if row[1] not in EMPTY_VALUE:
                                    try:
                                        val_err = float(row[1])
                                        if val_err > 0.0:
                                            if isotope_number in ALLOWED_ISOTOPE_NUMBERS:
                                                isotope_nums_with_non_zero_val_err.add(isotope_number)
                                        elif val_err == 0.0:
                                            if len(zero_shift_val_err) == 0:
                                                zero_shift_val_err = row[1]
                                            if isotope_number in ALLOWED_ISOTOPE_NUMBERS:
                                                isotope_nums_with_zero_val_err.add(isotope_number)
                                            if row[1][0] == '-':
                                                lp.data[idx][val_err_col] = row[1][1:]
                                        else:
                                            lp.data[idx][val_err_col] = abs(val_err)
                                            if isotope_number in ALLOWED_ISOTOPE_NUMBERS:
                                                isotope_nums_with_non_zero_val_err.add(isotope_number)
                                    except ValueError:
                                        pass
                                if row[2] in EMPTY_VALUE:
                                    empty_ambig_code = True
                                if row[3] not in EMPTY_VALUE:
                                    if isinstance(row[3], str):
                                        entity_id = int(row[3])
                                    else:
                                        entity_id = row[3]
                                    if entity_id not in isotope_nums_per_entity:
                                        isotope_nums_per_entity[entity_id] = set()
                                    isotope_nums_per_entity[entity_id].add(isotope_number)
                            except (ValueError, TypeError):
                                continue

                        if len(isotope_nums_with_empty_val_err) > 0:

                            for isotope_number in isotope_nums_with_empty_val_err:

                                if isotope_number in (1, 2, 13, 15, 19, 31):
                                    type_symbol = next(k for k, v in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS.items() if isotope_number in v)
                                    cs_val_err = get_first_sf_tag(sf, f'Chem_shift_{isotope_number}{type_symbol}_err')
                                else:
                                    continue

                                if len(cs_val_err) == 0 or cs_val_err in EMPTY_VALUE:
                                    empty_val_err_tag = f'_Assigned_chem_shift_list.Chem_shift_{isotope_number}{type_symbol}_err'
                                    if empty_val_err_tag not in empty_val_err_tags:
                                        empty_val_err_tags.append(empty_val_err_tag)
                                    continue

                                try:
                                    cs_val_err = abs(float(cs_val_err))
                                    if cs_val_err > CS_UNCERT_MAX:
                                        continue

                                    for idx, row in enumerate(dat):
                                        try:
                                            if int(row[0]) != isotope_number:
                                                continue
                                            if row[1] in EMPTY_VALUE:
                                                lp.data[idx][val_err_col] = cs_val_err
                                        except (ValueError, TypeError):
                                            continue

                                except ValueError:
                                    continue

                        if len(isotope_nums_with_zero_val_err) > 0:

                            for isotope_number in isotope_nums_with_zero_val_err:

                                if isotope_number in (1, 2, 13, 15, 19, 31):
                                    type_symbol = next(k for k, v in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS.items() if isotope_number in v)
                                    cs_val_err = get_first_sf_tag(sf, f'Chem_shift_{isotope_number}{type_symbol}_err')
                                else:
                                    continue

                                if len(cs_val_err) == 0 or cs_val_err in EMPTY_VALUE:
                                    zero_val_err_tag = f'_Assigned_chem_shift_list.Chem_shift_{isotope_number}{type_symbol}_err'
                                    if zero_val_err_tag not in zero_val_err_tags:
                                        zero_val_err_tags.append(zero_val_err_tag)
                                    continue

                                try:
                                    cs_val_err = abs(float(cs_val_err))
                                    if cs_val_err > CS_UNCERT_MAX:
                                        continue

                                    fill_non_zero_val_err = isotope_number not in isotope_nums_with_non_zero_val_err
                                    for idx, row in enumerate(dat):
                                        try:
                                            if int(row[0]) != isotope_number:
                                                continue
                                            if fill_non_zero_val_err or row[1] in EMPTY_VALUE:
                                                lp.data[idx][val_err_col] = cs_val_err
                                        except (ValueError, TypeError):
                                            continue

                                except ValueError:
                                    continue

                        if empty_ambig_code:

                            tags = ['Comp_ID', 'Atom_ID', 'Ambiguity_code', 'Ambiguity_set_ID']

                            dat = lp.get_tag(tags)

                            # wo coordinates (bmrbdep)
                            if self.__reg.report.getInputSourceIdOfCoord() < 0:

                                for idx, row in enumerate(dat):
                                    comp_id = row[0]
                                    atom_id = row[1]
                                    ambig_code = row[2]
                                    ambig_set_id = row[3]

                                    _ambig_code = 1 if self.__reg.sail_flag\
                                        else self.__reg.csStat.getMaxAmbigCodeWoSetId(comp_id, atom_id)

                                    checked = _ambig_code != 1

                                    if ambig_code in EMPTY_VALUE or (isinstance(ambig_code, str) and not ambig_code.isdigit()):
                                        checked = False
                                    else:
                                        ambig_code = int(ambig_code)
                                        if ambig_code not in ALLOWED_AMBIGUITY_CODES:
                                            checked = False
                                        else:
                                            if _ambig_code == 0:
                                                pass
                                            elif {ambig_code, _ambig_code} == {2, 3}:
                                                checked = False
                                            elif ambig_code in (4, 5) and ambig_set_id not in EMPTY_VALUE:
                                                checked = True

                                    if not checked and _ambig_code > 0:
                                        lp.data[idx][ambig_code_col] = _ambig_code

                            # with coordinates
                            else:

                                for idx, row in enumerate(dat):
                                    comp_id = row[0]
                                    atom_id = row[1]
                                    ambig_code = row[2]

                                    _ambig_code = 1 if self.__reg.sail_flag\
                                        else self.__reg.csStat.getMaxAmbigCodeWoSetId(comp_id, atom_id)

                                    checked = _ambig_code != 1

                                    if ambig_code in EMPTY_VALUE or (isinstance(ambig_code, str) and not ambig_code.isdigit()):
                                        pass
                                    else:
                                        ambig_code = int(ambig_code)
                                        if ambig_code not in ALLOWED_AMBIGUITY_CODES:
                                            checked = False
                                            if _ambig_code != 1:
                                                _ambig_code = '.'
                                        else:
                                            if _ambig_code == 0:
                                                pass
                                            elif {ambig_code, _ambig_code} == {2, 3}:
                                                checked = False

                                    if not checked:
                                        lp.data[idx][ambig_code_col] = _ambig_code

                except KeyError:
                    has_cs_lp = False

        sf_category = 'chem_shift_perturbation'

        if sf_category in self.__reg.sf_category_list:
            for list_id, sf in enumerate(master_entry.get_saveframes_by_category(sf_category), start=1):
                self.__reg.c2S.set_local_sf_id(sf, list_id)
                cs_ref_id = get_first_sf_tag(sf, 'Chem_shift_ref_set_ID')
                if not isinstance(cs_ref_id, int):
                    if len(cs_ref_id) == 0 or cs_ref_id in EMPTY_VALUE or not cs_ref_id.isdigit():
                        cs_ref_id = 1
                        set_sf_tag(sf, 'Chem_shift_ref_set_ID', cs_ref_id)
                    else:
                        cs_ref_id = int(cs_ref_id)
                if cs_ref_id in cs_ref_sf_framecode:
                    set_sf_tag(sf, 'Chem_shift_ref_set_label', f'${cs_ref_sf_framecode[cs_ref_id]}')
                else:
                    isotope_nums[cs_ref_id] = set()
                    if len(cs_ref_sf_framecode) > 0:
                        cs_ref_id = min(list(cs_ref_sf_framecode.keys()))
                        set_sf_tag(sf, 'Chem_shift_ref_set_ID', cs_ref_id)
                        set_sf_tag(sf, 'Chem_shift_ref_set_label', f'${cs_ref_sf_framecode[cs_ref_id]}')
                    else:
                        cs_ref_id = 1
                        set_sf_tag(sf, 'Chem_shift_ref_set_ID', cs_ref_id)

                smpl_cond_list_id = get_first_sf_tag(sf, 'Sample_condition_list_ID')
                if not isinstance(smpl_cond_list_id, int):
                    if len(smpl_cond_list_id) == 0 or smpl_cond_list_id in EMPTY_VALUE or not smpl_cond_list_id.isdigit():
                        smpl_cond_list_id = 1
                        set_sf_tag(sf, 'Sample_condition_list_ID', smpl_cond_list_id)
                    else:
                        smpl_cond_list_id = int(smpl_cond_list_id)
                if smpl_cond_list_id in smpl_cond_sf_framecode:
                    set_sf_tag(sf, 'Sample_condition_list_label', f'${smpl_cond_sf_framecode[smpl_cond_list_id]}')
                else:
                    if len(smpl_cond_sf_framecode) > 0:
                        smpl_cond_list_id = min(list(smpl_cond_sf_framecode.keys()))
                        set_sf_tag(sf, 'Sample_condition_list_ID', smpl_cond_list_id)
                        set_sf_tag(sf, 'Sample_condition_list_label', f'${smpl_cond_sf_framecode[smpl_cond_list_id]}')
                    else:
                        smpl_cond_list_id = 1
                        set_sf_tag(sf, 'Sample_condition_list_ID', smpl_cond_list_id)

                exp_list_sf_category = 'experiment_list'

                if exp_list_sf_category in self.__reg.sf_category_list:
                    exp_list_sf = master_entry.get_saveframes_by_category(exp_list_sf_category)[0]

                    exp_lp_category = '_Experiment'

                    try:

                        exp_lp = exp_list_sf.get_loop(exp_lp_category)
                        exp_lp.sort_rows('ID')

                    except KeyError:
                        exp_lp = None

                    if exp_lp is not None:

                        exp_tags = ['ID', 'Name', 'Sample_ID', 'Sample_label', 'Sample_state']

                        if set(exp_tags) & set(exp_lp.tags) == set(exp_tags):
                            exp_list = exp_lp.get_tag(exp_tags)

                            if len(exp_list) > 0:

                                lp_category = '_Chem_shift_perturbation_experiment'
                                lp = None

                                try:

                                    lp = sf.get_loop(lp_category)
                                    lp.sort_rows('Experiment_ID')

                                    tags = ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state']

                                    if set(tags) & set(lp.tags) == set(tags):
                                        exp_id_col = lp.tags.index('Experiment_ID')
                                        exp_name_col = lp.tags.index('Experiment_name')
                                        sample_id_col = lp.tags.index('Sample_ID')
                                        sample_label_col = lp.tags.index('Sample_label')
                                        sample_state_col = lp.tags.index('Sample_state')

                                        cs_exp_list = lp.get_tag(tags)
                                        reserved_ids, duplicated_idxs, reserved_names = [], [], []

                                        for idx, cs_exp in enumerate(cs_exp_list):
                                            if cs_exp[0] not in EMPTY_VALUE:
                                                if cs_exp[0] not in reserved_ids:
                                                    reserved_ids.append(cs_exp[0])
                                                else:
                                                    duplicated_idxs.append(idx)
                                                exp = next((exp for exp in exp_list if exp[0] == cs_exp[0]), None)
                                                if exp is not None:
                                                    if cs_exp[2:5] != exp[2:5]:
                                                        lp.data[idx][exp_name_col] = exp[1]
                                                        lp.data[idx][sample_id_col] = exp[2]
                                                        lp.data[idx][sample_label_col] = exp[3]
                                                        lp.data[idx][sample_state_col] = exp[4]
                                                    reserved_names.append(exp[1])
                                            else:
                                                exp = next((exp for exp in exp_list if exp[1] not in reserved_names), None)
                                                if exp is not None:
                                                    lp.data[idx][exp_id_col] = exp_list.index(exp) + 1
                                                    lp.data[idx][exp_name_col] = exp[1]
                                                    lp.data[idx][sample_id_col] = exp[2]
                                                    lp.data[idx][sample_label_col] = exp[3]
                                                    lp.data[idx][sample_state_col] = exp[4]
                                                    reserved_names.append(exp[1])

                                        if len(duplicated_idxs) > 0:
                                            for idx in reversed(duplicated_idxs):
                                                del lp.data[idx]

                                except (KeyError, TypeError, ValueError):

                                    if lp is not None:
                                        del sf[lp]

                                    items = ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state',
                                             'Entry_ID', 'Chem_shift_perturbation_list_ID']

                                    cs_list_id = get_first_sf_tag(sf, 'ID')

                                    lp = pynmrstar.Loop.from_scratch(lp_category)

                                    tags = [f'{lp_category}.{item}' for item in items]

                                    lp.add_tag(tags)

                                    for exp in exp_list:
                                        row = exp
                                        row.extend([self.__reg.entry_id, cs_list_id])
                                        lp.add_data(row)

                                    lp.sort_rows('Experiment_ID')

                                    sf.add_loop(lp)

                try:

                    lp_category = '_Chem_shift_perturbation_experiment'

                    lp = sf.get_loop(lp_category)

                    tags = ['Chem_shift_val_err', 'Difference_chem_shift_val_err']

                    if set(tags) & set(lp.tags) == set(tags):

                        dat = lp.get_tag(tags)

                        val_err_col = lp.tags.index('Chem_shift_val_err')
                        diff_val_err_col = lp.tags.index('Difference_chem_shift_val_err')

                        for idx, row in enumerate(dat):
                            try:
                                if row[0] not in EMPTY_VALUE:
                                    try:
                                        val_err = float(row[0])
                                        if val_err > 0.0:
                                            pass
                                        elif val_err == 0.0:
                                            if row[0][0] == '-':
                                                lp.data[idx][val_err_col] = row[0][1:]
                                        else:
                                            lp.data[idx][val_err_col] = abs(val_err)
                                    except ValueError:
                                        pass
                                if row[1] not in EMPTY_VALUE:
                                    try:
                                        diff_val_err = float(row[1])
                                        if diff_val_err > 0.0:
                                            pass
                                        elif diff_val_err == 0.0:
                                            if row[1][0] == '-':
                                                lp.data[idx][diff_val_err_col] = row[1][1:]
                                        else:
                                            lp.data[idx][diff_val_err_col] = abs(diff_val_err)
                                    except ValueError:
                                        pass
                            except (ValueError, TypeError):
                                continue

                except KeyError:
                    pass

        # section 12: anomalous chemical shift assignments

        sf_category = 'assembly'

        assembly_id = assembly_label = '.'

        if sf_category in self.__reg.sf_category_list:
            asm_sf = master_entry.get_saveframes_by_category(sf_category)[0]
            try:
                assembly_id = int(get_first_sf_tag(asm_sf, 'ID'))
                assembly_label = f"${get_first_sf_tag(asm_sf, 'Sf_framecode')}"
            except ValueError:
                pass

        sf_category = 'entity'

        has_non_polymer = has_nstd_monomer = False
        entity_dict = {}
        polymer_common_types, non_polymer_types, nstd_monomers = [], [], []

        if sf_category in self.__reg.sf_category_list:
            for sf in master_entry.get_saveframes_by_category(sf_category):
                try:
                    entity_id = int(get_first_sf_tag(sf, 'ID'))
                    _type = get_first_sf_tag(sf, 'Type')
                    polymer_type = get_first_sf_tag(sf, 'Polymer_type')
                    polymer_common_type = get_first_sf_tag(sf, 'Polymer_common_type')
                    num_of_monomers = get_first_sf_tag(sf, 'Number_of_monomers')

                    if len(_type) == 0 or _type in ('solvent', 'water') or len(polymer_type) == 0:

                        lp_category = '_Entity_comp_index'

                        protein = dna = rna = False

                        try:

                            lp = sf.get_loop(lp_category)

                            dat = lp.get_tag(['Comp_ID'])

                            for row in dat:
                                if row not in EMPTY_VALUE:
                                    if row in ('ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                                               'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                                               'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                                               'SER', 'THR', 'TRP', 'TYR', 'VAL',
                                               'ASX', 'GLX'):
                                        protein = True
                                    elif row in ('DA', 'DC', 'DG', 'DT', 'DU', 'DI'):
                                        dna = True
                                    elif row in ('A', 'C', 'G', 'I', 'T', 'U'):
                                        rna = True
                                    else:
                                        protein = dna = rna = False
                                        break

                            if protein or dna or rna:
                                _type = 'polymer'
                                set_sf_tag(sf, 'Type', _type)

                                if protein:
                                    polymer_type = 'polypeptide(L)'
                                    polymer_common_type = 'protein'
                                elif dna and not rna:
                                    polymer_type = 'polydeoxyribonucleotide'
                                    polymer_common_type = 'DNA'
                                elif rna and not dna:
                                    polymer_type = 'polyribonucleotide'
                                    polymer_common_type = 'RNA'
                                else:
                                    polymer_type = 'polydeoxyribonucleotide/polyribonucleotide hybrid'
                                    polymer_common_type = 'DNA/RNA hybrid'
                                set_sf_tag(sf, 'Polymer_type', polymer_type)
                                set_sf_tag(sf, 'Polymer_common_type', polymer_common_type)

                        except KeyError:
                            pass

                    if polymer_type not in EMPTY_VALUE:
                        if polymer_type == 'polypeptide(L)':
                            polymer_common_type = 'protein'
                        elif polymer_type == 'polydeoxyribonucleotide':
                            polymer_common_type = 'DNA'
                        elif polymer_type == 'polyribonucleotide':
                            polymer_common_type = 'RNA'
                        elif polymer_type == 'polydeoxyribonucleotide/polyribonucleotide hybrid':
                            polymer_common_type = 'DNA/RNA hybrid'
                        if len(polymer_common_type) > 0:
                            set_sf_tag(sf, 'Polymer_common_type', polymer_common_type)

                    if polymer_common_type not in EMPTY_VALUE:
                        polymer_common_types.append({polymer_common_type: entity_id})

                    try:

                        lp_category = '_Entity_comp_index'

                        lp = sf.get_loop(lp_category)

                        lp_category = '_Entity_poly_seq'

                        _lp = sf.get_loop(lp_category)

                        del sf[_lp]

                        _lp = pynmrstar.Loop.from_scratch(lp_category)

                        tags = ['Hetero', 'Mon_ID', 'Num', 'Comp_index_ID', 'Entry_ID', 'Entity_ID']

                        _lp.add_tag(tags)

                        dat = lp.get_tag(['ID', 'Comp_ID'])

                        for row in dat:
                            _row = ['.', row[1], row[0], row[0], self.__reg.entry_id, entity_id]
                            _lp.add_data(_row)

                        sf.add_loop(_lp)

                        if num_of_monomers in EMPTY_VALUE:
                            num_of_monomers = len(dat)

                            lp_category = '_Entity_comp_index'

                            protein = dna = rna = False

                            lp = sf.get_loop(lp_category)

                            dat = lp.get_tag(['Comp_ID'])

                            for row in dat:
                                if row not in EMPTY_VALUE:
                                    if row in ('ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                                               'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                                               'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                                               'SER', 'THR', 'TRP', 'TYR', 'VAL',
                                               'ASX', 'GLX'):
                                        protein = True
                                    elif row in ('DA', 'DC', 'DG', 'DT', 'DU', 'DI'):
                                        dna = True
                                    elif row in ('A', 'C', 'G', 'I', 'T', 'U'):
                                        rna = True
                                    else:
                                        protein = dna = rna = False
                                        break

                            if protein or dna or rna:
                                set_sf_tag(sf, 'Number_of_monomers', num_of_monomers)

                    except KeyError:
                        pass

                    has_non_polymer = _type == 'non-polymer'

                    nstd_monomer = get_first_sf_tag(sf, 'Nstd_monomer')
                    has_nstd_monomer = nstd_monomer == 'yes'

                    lp_category = '_Entity_comp_index'

                    total_cys = 0

                    try:

                        lp = sf.get_loop(lp_category)

                        dat = lp.get_tag(['Comp_ID'])

                        for row in dat:
                            if row not in EMPTY_VALUE:
                                if row not in STD_MON_DICT:
                                    if row not in nstd_monomers:
                                        nstd_monomers.append(row)
                                    if _type != 'non-polymer':
                                        has_nstd_monomer = True
                                if row in ('CYS', 'DCY'):
                                    total_cys += 1

                        id_col = lp.tags.index('ID')
                        auth_seq_id_col = lp.tags.index('Auth_seq_ID')

                        trial = 0

                        while True:

                            offset = None
                            reserved_auth_seq_ids, conflict_auth_seq_ids = [], []
                            for row in reversed(lp):
                                seq_key = (entity_id, int(row[id_col]))
                                if seq_key in label_to_auth_seq:
                                    row[auth_seq_id_col] = label_to_auth_seq[seq_key]
                                    offset = label_to_auth_seq[seq_key] - int(row[id_col])
                                elif offset is not None:
                                    if row[auth_seq_id_col] in EMPTY_VALUE:
                                        row[auth_seq_id_col] = int(row[id_col]) + offset
                                if row[auth_seq_id_col] in EMPTY_VALUE:
                                    pass
                                elif isinstance(row[auth_seq_id_col], int):
                                    if row[auth_seq_id_col] in reserved_auth_seq_ids:
                                        conflict_auth_seq_ids.append(row[auth_seq_id_col])
                                    else:
                                        reserved_auth_seq_ids.append(row[auth_seq_id_col])
                                else:
                                    try:
                                        auth_seq_id = int(row[auth_seq_id_col])
                                        if auth_seq_id in reserved_auth_seq_ids:
                                            conflict_auth_seq_ids.append(auth_seq_id)
                                        else:
                                            reserved_auth_seq_ids.append(auth_seq_id)
                                    except ValueError:
                                        row[auth_seq_id_col] = None

                            if len(conflict_auth_seq_ids) == 0\
                               and (all(row[auth_seq_id_col] not in EMPTY_VALUE for row in lp) or len(lp) == 1):
                                break

                            for row in reversed(lp):
                                if row[auth_seq_id_col] in EMPTY_VALUE:
                                    pass
                                elif isinstance(row[auth_seq_id_col], int):
                                    if row[auth_seq_id_col] in conflict_auth_seq_ids:
                                        row[auth_seq_id_col] = None
                                else:
                                    auth_seq_id = int(row[auth_seq_id_col])
                                    if auth_seq_id in conflict_auth_seq_ids:
                                        row[auth_seq_id_col] = None
                                seq_key = (entity_id, int(row[id_col]))
                                if seq_key in label_to_auth_seq:
                                    row[auth_seq_id_col] = label_to_auth_seq[seq_key]
                                    offset = label_to_auth_seq[seq_key] - int(row[id_col])
                                elif offset is not None:
                                    if row[auth_seq_id_col] in EMPTY_VALUE:
                                        row[auth_seq_id_col] = int(row[id_col]) + offset

                            trial += 1

                            if trial > 4:
                                break

                    except KeyError:
                        pass

                    _lp_category = '_Entity_common_name'
                    common_names = []

                    try:

                        _lp = sf.get_loop(_lp_category)

                        _dat = _lp.get_tag(['Name'])

                        for _row in _dat:
                            if _row not in EMPTY_VALUE:
                                for _name in _row.split(','):
                                    common_names.append(_name.strip())

                    except KeyError:
                        pass

                    _sf_category = 'natural_source'
                    gene_mnemonic = []
                    if _sf_category in self.__reg.sf_category_list:
                        for _sf in master_entry.get_saveframes_by_category(_sf_category):

                            _lp_category = '_Entity_natural_src'

                            try:

                                _lp = _sf.get_loop(_lp_category)

                                if 'Gene_mnemonic' in _lp.tags:
                                    _dat = _lp.get_tag(['Entity_ID', 'Gene_mnemonic'])

                                    for _row in _dat:
                                        if int(_row[0]) == entity_id and _row[1] not in EMPTY_VALUE:
                                            for _name in _row[1].split(','):
                                                gene_mnemonic.append(_name.strip())
                                            break

                            except (KeyError, ValueError):
                                pass

                    if 'peptide' in polymer_type:
                        sample_type = 'protein' if len(lp) >= 24 else 'peptide'
                        entity_dict[entity_id] = {'name': get_first_sf_tag(sf, 'Name'),
                                                  'sf_framecode': get_first_sf_tag(sf, 'Sf_framecode'),
                                                  'sample_type': sample_type,
                                                  'assembly_id': assembly_id,
                                                  'assembly_label': assembly_label,
                                                  'fragment': get_first_sf_tag(sf, 'Fragment'),
                                                  'mutation': get_first_sf_tag(sf, 'Mutation'),
                                                  'details': get_first_sf_tag(sf, 'Details'),
                                                  'common_names': common_names,
                                                  'gene_mnemonic': gene_mnemonic,
                                                  'total_cys': total_cys
                                                  }

                    elif 'ribonucleotide' in polymer_type:
                        if polymer_type == 'polydeoxyribonucleotide':
                            sample_type = 'DNA'
                        elif polymer_type == 'polyribonucleotide':
                            sample_type = 'RNA'
                        elif polymer_type == 'polydeoxyribonucleotide/polyribonucleotide hybrid':
                            sample_type = 'DNA/RNA hybrid'
                        else:
                            sample_type = '?'
                        entity_dict[entity_id] = {'name': get_first_sf_tag(sf, 'Name'),
                                                  'sf_framecode': get_first_sf_tag(sf, 'Sf_framecode'),
                                                  'sample_type': sample_type,
                                                  'assembly_id': assembly_id,
                                                  'assembly_label': assembly_label,
                                                  'fragment': get_first_sf_tag(sf, 'Fragment'),
                                                  'mutation': get_first_sf_tag(sf, 'Mutation'),
                                                  'details': get_first_sf_tag(sf, 'Details'),
                                                  'common_names': common_names,
                                                  'gene_mnemonic': gene_mnemonic
                                                  }

                    elif 'saccharide' in polymer_type or 'carbohydrate' in polymer_type:
                        entity_dict[entity_id] = {'name': get_first_sf_tag(sf, 'Name'),
                                                  'sf_framecode': get_first_sf_tag(sf, 'Sf_framecode'),
                                                  'sample_type':
                                                  'carbohydrate' if polymer_type.startswith('carbo') else 'polysaccharide',
                                                  'assembly_id': assembly_id,
                                                  'assembly_label': assembly_label,
                                                  'fragment': get_first_sf_tag(sf, 'Fragment'),
                                                  'mutation': get_first_sf_tag(sf, 'Mutation'),
                                                  'details': get_first_sf_tag(sf, 'Details'),
                                                  'common_names': common_names,
                                                  'gene_mnemonic': gene_mnemonic
                                                  }

                    if _type == 'non-polymer':
                        non_polymer_type = 'ligand' if 'ION' not in get_first_sf_tag(sf, 'Name') else 'metal ion'
                        non_polymer_types.append({non_polymer_type: entity_id})
                        entity_dict[entity_id] = {'name': get_first_sf_tag(sf, 'Name'),
                                                  'sf_framecode': get_first_sf_tag(sf, 'Sf_framecode'),
                                                  'sample_type': non_polymer_type,
                                                  'assembly_id': assembly_id,
                                                  'assembly_label': assembly_label,
                                                  'fragment': get_first_sf_tag(sf, 'Fragment'),
                                                  'mutation': get_first_sf_tag(sf, 'Mutation'),
                                                  'details': get_first_sf_tag(sf, 'Details'),
                                                  'common_names': common_names,
                                                  'gene_mnemonic': gene_mnemonic
                                                  }

                except ValueError:
                    continue

        # section 9: NMR applied experiments

        spectrometer_dict = {}

        sf_category = 'NMR_spectrometer_list'

        if sf_category in self.__reg.sf_category_list:

            for sf in master_entry.get_saveframes_by_category(sf_category):

                lp_category = '_NMR_spectrometer_view'

                try:

                    lp = sf.get_loop(lp_category)

                    tags = ['ID', 'Name', 'Manufacturer', 'Model', 'Field_strength']

                    has_serial_number = 'Serial_number' in lp.tags
                    has_details = 'Details' in lp.tags

                    if has_serial_number:
                        tags.append('Serial_number')
                    if has_details:
                        tags.append('Details')

                    dat = lp.get_tag(tags)

                    for row in dat:
                        try:
                            spectrometer_id = int(row[0])
                            if spectrometer_id in spectrometer_dict:
                                continue
                            spectrometer_dict[spectrometer_id] = {'label': row[1],
                                                                  'manufacturer': row[2],
                                                                  'model': row[3],
                                                                  'field_strength': row[4]}
                            if has_serial_number:
                                spectrometer_dict[spectrometer_id]['serial_number'] = row[5]
                            if has_details:
                                spectrometer_dict[spectrometer_id]['details'] = row[6 if has_serial_number else 5]
                        except ValueError:
                            continue

                    lp.sort_rows('ID')

                except KeyError:
                    continue

        field_strength_max_didits = 0
        if len(spectrometer_dict) > 0:
            field_strength_max_didits = getMaxEffDigits([str(spectrometer['field_strength'])
                                                         for spectrometer in spectrometer_dict.values()])

        sf_category = 'experiment_list'

        if sf_category in self.__reg.sf_category_list:

            for sf in master_entry.get_saveframes_by_category(sf_category):

                lp_category = '_Experiment'

                try:

                    lp = sf.get_loop(lp_category)

                    tags = ['NMR_spectrometer_ID', 'NMR_spectrometer_label']

                    dat = lp.get_tag(tags)

                    for idx, row in enumerate(dat):
                        if row[1] in EMPTY_VALUE:
                            continue
                        spectrometer_label = row[1][1:] if row[1][0] == '$' else row[1]
                        if row[0] in EMPTY_VALUE:
                            spectrometer_id = next((k for k, v in spectrometer_dict.items()
                                                    if v['label'] == spectrometer_label), None)
                            if spectrometer_id is not None:
                                spectrometer_id_col = lp.tags.index('NMR_spectrometer_ID')
                                lp.data[idx][spectrometer_id_col] = str(spectrometer_id)

                except KeyError:
                    continue

        if len(spectrometer_dict) > 0:

            for spectrometer_id, spectrometer in spectrometer_dict.items():

                sf_framecode = spectrometer['label']
                if sf_framecode not in EMPTY_VALUE and ' ' in sf_framecode:
                    sf_framecode = f'NMR_spectrometer_{spectrometer_id}'

                sf = next((sf for sf in master_entry.frame_list
                           if sf.category == 'NMR_spectrometer' and sf.name == sf_framecode), None)

                if sf is None:
                    sf = pynmrstar.Saveframe.from_scratch(sf_framecode, '_NMR_spectrometer')
                    sf.add_tag('Sf_category', 'NMR_spectrometer')
                    sf.add_tag('Sf_framecode', sf_framecode)
                    sf.add_tag('Entry_ID', self.__reg.entry_id)
                    sf.add_tag('ID', spectrometer_id)
                    sf.add_tag('Name', spectrometer['label'])
                    sf.add_tag('Details', spectrometer['details'] if 'details' in spectrometer else '.')
                    sf.add_tag('Manufacturer', spectrometer['manufacturer'])
                    sf.add_tag('Model', spectrometer['model'])
                    sf.add_tag('Serial_number', spectrometer['serial_number'] if 'serial_nubmer' in spectrometer else '.')
                    sf.add_tag('Field_strength', spectrometer['field_strength'])

                    master_entry.add_saveframe(sf)

        # section 8: sample conditions

        sf_category = 'sample_conditions'

        pH = ionic_strength = None

        if sf_category in self.__reg.sf_category_list:
            for sf in master_entry.get_saveframes_by_category(sf_category):

                lp_category = '_Sample_condition_variable'

                try:

                    lp = sf.get_loop(lp_category)

                    tags = ['Type', 'Val', 'Val_units']

                    dat = lp.get_tag(tags)

                    val_col = lp.tags.index('Val')
                    val_err_col = lp.tags.index('Val_err')
                    val_unit_col = lp.tags.index('Val_units')

                    for idx, row in enumerate(dat):
                        if not isinstance(row[1], str):
                            continue
                        minus_code_count = row[1].count('-')
                        if minus_code_count == 1:
                            if row[1].startswith('+-') or row[1].startswith('-+'):
                                lp.data[idx][val_col] = round(0, eff_digits(row[1][2:]))
                                lp.data[idx][val_err_col] = row[1][2:]
                            elif row[1].startswith('-'):
                                pass
                            else:
                                try:
                                    _val = row[1].split('-')
                                    _val_0 = float(_val[0])
                                    _val_1 = float(_val[1])
                                    _eff_digits_0 = eff_digits(_val[0])
                                    _eff_digits_1 = eff_digits(_val[1])
                                    max_eff_digits = max([_eff_digits_0, _eff_digits_1])
                                    lp.data[idx][val_col] = round((_val_0 + _val_1) / 2.0, max_eff_digits)
                                    lp.data[idx][val_err_col] = round(abs(_val_0 - _val_1) / 2.0, max_eff_digits)
                                except ValueError:
                                    pass
                        elif minus_code_count == 2 and row[1].startswith('-'):
                            try:
                                _val = row[1][1:].split('-')
                                _val_0 = float(f'-{_val[0]}')
                                _val_1 = float(_val[1])
                                _eff_digits_0 = eff_digits(_val[0])
                                _eff_digits_1 = eff_digits(_val[1])
                                max_eff_digits = max([_eff_digits_0, _eff_digits_1])
                                lp.data[idx][val_col] = round((_val_0 + _val_1) / 2.0, max_eff_digits)
                                lp.data[idx][val_err_col] = round(abs(_val_0 - _val_1) / 2.0, max_eff_digits)
                            except ValueError:
                                pass

                        for delimiter in ('/', ':', ',', '~'):
                            if delimiter not in row[1]:
                                continue
                            code_count = row[1].count(delimiter)
                            if code_count == 1:
                                try:
                                    _val = row[1].split(delimiter)
                                    _val_0 = float(_val[0])
                                    _val_1 = float(_val[1])
                                    _eff_digits_0 = eff_digits(_val[0] if isinstance(_val[0], str) else str(_val[0]))
                                    _eff_digits_1 = eff_digits(_val[1] if isinstance(_val[1], str) else str(_val[1]))
                                    max_eff_digits = max([_eff_digits_0, _eff_digits_1])
                                    lp.data[idx][val_col] = round((_val_0 + _val_1) / 2.0, max_eff_digits)
                                    lp.data[idx][val_err_col] = round(abs(_val_0 - _val_1) / 2.0, max_eff_digits)
                                    break
                                except ValueError:
                                    pass

                        if row[0] == 'pressure' and isinstance(row[1], str) and row[1].lower() == 'ambient':  # and row[2] == 'atm':
                            lp.data[idx][val_col] = '1'
                            lp.data[idx][val_unit_col] = 'atm'

                        elif row[0] == 'temperature' and row[2] == 'C':
                            try:
                                celcius = float(row[1])
                                kelvin = round(celcius + 273.15, eff_digits(row[1]))
                                lp.data[idx][val_col] = kelvin
                                lp.data[idx][val_unit_col] = 'K'
                            except (ValueError, TypeError):
                                pass

                        elif row[0] == 'ionic strength':
                            if isinstance(row[1], str) and row[1].lower()[0] in ('n', 'u'):
                                lp.data[idx][val_col] = '.'

                            ionic_strength_pat = re.compile(r'\s*(\d+)\s*(mM|M).*')

                            if isinstance(lp.data[idx][val_col], float) and lp.data[idx][val_col] > 0.0:
                                ionic_strength = f'{lp.data[idx][val_col]} {row[2]}'
                            else:
                                try:
                                    val = float(lp.data[idx][val_col])
                                    if val > 0.0:
                                        ionic_strength = f'{lp.data[idx][val_col]} {row[2]}'
                                except (ValueError, TypeError):
                                    if row[1] in EMPTY_VALUE:
                                        lp.data[idx][val_col] = '?'
                                    else:
                                        if ionic_strength_pat.match(row[1]):
                                            g = ionic_strength_pat.search(row[1]).groups()
                                            lp.data[idx][val_col] = g[0]
                                            lp.data[idx][val_unit_col] = g[1]

                                            ionic_strength = f'{g[0]} {g[1]}'

                                        else:
                                            ionic_strength = f'{lp.data[idx][val_col]} {row[2]}'

                                            _lp = pynmrstar.Loop.from_scratch(lp_category)

                                            tags = ['Type', 'Val', 'Val_err', 'Val_units']

                                            _lp.add_tag(tags)

                                            data = lp.get_tag(tags)

                                            for row in data:
                                                _lp.add_data(row)

                        elif row[0].startswith('pH') or row[0].startswith('pD'):
                            if isinstance(row[1], str) and row[1].lower()[0] in ('n', 'u'):
                                lp.data[idx][val_col] = '.'

                            if isinstance(lp.data[idx][val_col], float) and lp.data[idx][val_col] > 0.0:
                                pH = f'{row[0]}: {lp.data[idx][val_col]!r}'
                            else:
                                try:
                                    val = float(lp.data[idx][val_col])
                                    if val > 0.0:
                                        pH = f'{row[0]}: {lp.data[idx][val_col]!r}'
                                except (ValueError, TypeError):
                                    pass

                except KeyError:
                    continue

        # section 7: sample description

        sf_category = 'sample'

        solvent_with_percent_pat = re.compile(r'(\d+|\d+\.\d+)%\s*([\S ]+)')
        deuterated_pat = re.compile(r'(?:[Dd]2[Oo]|\S+(\s\S+)?-[Dd]\d+|.*deuterate.*)')
        perdeuterated_pat = re.compile(r'.*perdeuterate.*')
        percent_value_pat = re.compile(r'(\d+|\d+\.\d+)\s*%')
        redundant_solvent_pat = re.compile(r'\d?%?\s*([\S]+)\s+and\s+100%\s*([\S]+)')
        and_pat = r'\s[Aa][Nn][Dd]\s'
        has_salt = has_buffer = has_reducing_agent = has_chelating_agent = False
        default_internal_reference = 'DSS'
        is_single_protein_ligand_complex =\
            len(polymer_common_types) == 1 and 'protein' in polymer_common_types[0]\
            and len(non_polymer_types) == 1 and 'ligand' in non_polymer_types[0]

        def is_natural_abundance(isotopic_labeling: str):
            isotopic_labeling = isotopic_labeling.lower()
            return 'abundance' in isotopic_labeling\
                or isotopic_labeling.startswith('non') or isotopic_labeling.startswith('not ')\
                or isotopic_labeling.startswith('n/a') or isotopic_labeling.startswith('natural')\
                or isotopic_labeling.startswith('na ') or isotopic_labeling.startswith('na-')\
                or isotopic_labeling in ('na', 'no', '/', '1h')\
                or isotopic_labeling.startswith('non-lab')\
                or isotopic_labeling.startswith('nonlab')\
                or isotopic_labeling.startswith('unlab')\
                or 'protonated' in isotopic_labeling

        def cross_check_entity(key, entity):
            for f in entity['fragment'].split():
                _f = f.lower()
                if len(f) > 3 and (_f in key or key in _f):
                    return True
            for f in entity['mutation'].split():
                _f = f.lower()
                if len(f) > 2 and (_f in key or key in _f):
                    return True
            for f in entity['details'].split():
                _f = f.lower()
                if len(f) > 3 and (_f in key or key in _f):
                    return True
            for n in entity['common_names']:
                for f in n.split():
                    _f = f.lower()
                    if len(f) > 3 and (_f in key or key in _f):
                        return True
            for n in entity['gene_mnemonic']:
                for f in n.split():
                    _f = f.lower()
                    if len(f) > 2 and (_f in key or key in _f):
                        return True
            for k in key.split():
                if any((k in _f or _f in k) for _f in entity['sf_framecode'].lower().split() if len(_f) > 3)\
                   or any((k in _f or _f in k) for _f in entity['name'].lower().split() if len(_f) > 3):
                    return True

            if key in ('k-pi', 'na-pi', 'na-p', 'kpi', 'napi', 'nap')\
               or 'buffer' in key\
               or 'acetate' in key\
               or 'acetic' in key\
               or 'ch3coo' in key\
               or 'citric' in key\
               or 'hepes' in key\
               or 'pipes' in key\
               or 'mes' in key\
               or 'mops' in key\
               or ('morpholino' in key and 'propane' in key and 'sulfonic' in key)\
               or 'pbs' in key\
               or 'tris' in key\
               or 'phosphate' in key\
               or 'po4' in key\
               or 'kphos' in key\
               or 'naphos' in key:
                return False

            self.__reg.pA.clear()
            self.__reg.pA.setReferenceSequence(list(entity['name'].replace(' ', '').lower()), 'REFNAME')
            self.__reg.pA.addTestSequence(list(key.replace(' ', '').lower()), 'NAME')
            self.__reg.pA.doAlign()

            myAlign = self.__reg.pA.getAlignment('NAME')
            matched, unmapped, conflict, _, _ = getScoreOfSeqAlign(myAlign)

            return (matched > 3 and matched > conflict)\
                or (matched >= 3 and unmapped == conflict == 0)

        if sf_category in self.__reg.sf_category_list:
            sf_list = master_entry.get_saveframes_by_category(sf_category)
            is_single_sample_loop = len(sf_list) == 1
            for sf in sf_list:
                try:
                    sample_id = int(get_first_sf_tag(sf, 'ID'))

                    solvent_system, solvent_isotope = {}, {}
                    _solvent_system = get_first_sf_tag(sf, 'Solvent_system')

                    if '/' in _solvent_system:
                        for _solvent in _solvent_system.split('/'):
                            _solvent = _solvent.strip(',').strip()
                            if solvent_with_percent_pat.match(_solvent):
                                g = solvent_with_percent_pat.search(_solvent).groups()
                                solvent_name = g[1].strip()
                                if solvent_name.lower() == 'water':
                                    solvent_name = 'H2O'
                                if redundant_solvent_pat.match(solvent_name):
                                    h = redundant_solvent_pat.search(solvent_name)
                                    if h[1] == h[2]:
                                        solvent_name = h[1]
                                solvent_system[solvent_name] = float(g[0]) if '.' in g[0] else int(g[0])
                                solvent_isotope[solvent_name] = '[U-2H]' if deuterated_pat.match(solvent_name)\
                                    and not perdeuterated_pat.match(solvent_name)\
                                    else '[U-?% 2H]' if perdeuterated_pat.match(solvent_name) else 'natural abundance'
                            else:
                                if _solvent.lower() == 'water':
                                    _solvent = 'H2O'
                                solvent_system[_solvent] = 0
                                solvent_isotope[_solvent] = '[U-2H]' if deuterated_pat.match(_solvent)\
                                    and not perdeuterated_pat.match(_solvent)\
                                    else '[U-?% 2H]' if perdeuterated_pat.match(_solvent) else 'natural abundance'

                        total_v = sum(v for v in solvent_system.values())
                        if total_v < 100:
                            zero_comps = len([v for v in solvent_system.values() if v == 0])
                            if zero_comps > 0:
                                resid = (100 - total_v) // zero_comps
                                solvent_system = {k: v if v > 0 else resid for k, v in solvent_system.items()}

                    elif re.search(and_pat, _solvent_system):
                        _solvent_systems = re.split(and_pat, _solvent_system)
                        for _solvent in _solvent_systems:
                            _solvent = _solvent.strip(',').strip()
                            if solvent_with_percent_pat.match(_solvent):
                                g = solvent_with_percent_pat.search(_solvent).groups()
                                solvent_name = g[1].strip()
                                if solvent_name.lower() == 'water':
                                    solvent_name = 'H2O'
                                if redundant_solvent_pat.match(solvent_name):
                                    h = redundant_solvent_pat.search(solvent_name)
                                    if h[1] == h[2]:
                                        solvent_name = h[1]
                                solvent_system[solvent_name] = float(g[0]) if '.' in g[0] else int(g[0])
                                solvent_isotope[solvent_name] = '[U-2H]' if deuterated_pat.match(solvent_name)\
                                    and not perdeuterated_pat.match(solvent_name)\
                                    else '[U-?% 2H]' if perdeuterated_pat.match(solvent_name) else 'natural abundance'
                            else:
                                if _solvent.lower() == 'water':
                                    _solvent = 'H2O'
                                solvent_system[_solvent] = 100
                                solvent_isotope[_solvent] = '[U-2H]' if deuterated_pat.match(_solvent)\
                                    and not perdeuterated_pat.match(_solvent)\
                                    else '[U-?% 2H]' if perdeuterated_pat.match(_solvent) else 'natural abundance'

                        total_v = sum(v for v in solvent_system.values())
                        if total_v < 100:
                            zero_comps = len([v for v in solvent_system.values() if v == 0])
                            if zero_comps > 0:
                                resid = (100 - total_v) // zero_comps
                                solvent_system = {k: v if v > 0 else resid for k, v in solvent_system.items()}

                    elif _solvent_system not in EMPTY_VALUE:
                        if solvent_with_percent_pat.match(_solvent_system):
                            g = solvent_with_percent_pat.search(_solvent_system).groups()
                            solvent_name = g[1].strip()
                            if solvent_name.lower() == 'water':
                                solvent_name = 'H2O'
                            if redundant_solvent_pat.match(solvent_name):
                                h = redundant_solvent_pat.search(solvent_name)
                                if h[1] == h[2]:
                                    solvent_name = h[1]
                            solvent_system[solvent_name] = float(g[0]) if '.' in g[0] else int(g[0])
                            solvent_isotope[solvent_name] = '[U-2H]' if deuterated_pat.match(solvent_name)\
                                and not perdeuterated_pat.match(solvent_name)\
                                else '[U-?% 2H]' if perdeuterated_pat.match(solvent_name) else 'natural abundance'
                        else:
                            if _solvent_system.lower() == 'water':
                                _solvent_system = 'H2O'
                            solvent_system[_solvent_system] = 100
                            solvent_isotope[_solvent_system] = '[U-2H]' if deuterated_pat.match(_solvent_system)\
                                and not perdeuterated_pat.match(_solvent_system)\
                                else '[U-?% 2H]' if perdeuterated_pat.match(_solvent_system) else 'natural abundance'

                    lp_category = '_Sample_component'

                    lp = sf.get_loop(lp_category)

                    tags = ['ID',
                            'Mol_common_name',
                            'Isotopic_labeling',
                            'Assembly_ID',
                            'Assembly_label',
                            'Entity_ID',
                            'Entity_label',
                            'Type',
                            'Concentration_val',
                            'Concentration_val_min',
                            'Concentration_val_max',
                            'Concentration_val_units',
                            'Concentration_val_err']

                    dat = lp.get_tag(tags)

                    id_col = lp.tags.index('ID')
                    mol_common_name_col = lp.tags.index('Mol_common_name')
                    isotopic_labeling_col = lp.tags.index('Isotopic_labeling')
                    assembly_id_col = lp.tags.index('Assembly_ID')
                    assembly_label_col = lp.tags.index('Assembly_label')
                    entity_id_col = lp.tags.index('Entity_ID')
                    entity_label_col = lp.tags.index('Entity_label')
                    type_col = lp.tags.index('Type')
                    concentration_val_col = lp.tags.index('Concentration_val')
                    concentration_val_min_col = lp.tags.index('Concentration_val_min')
                    concentration_val_max_col = lp.tags.index('Concentration_val_max')
                    concentration_val_units_col = lp.tags.index('Concentration_val_units')
                    sample_id_col = lp.tags.index('Sample_ID')
                    entry_id_col = lp.tags.index('Entry_ID')

                    has_mand_concentration_val = True
                    has_poly_entity = has_hvy_shifts_with_natural_abundance = has_perdeuteration_sol = False
                    solvent_in_sample_loop, appended_solvent_in_sample_loop, touched_entity_id = [], [], []

                    for idx, row in enumerate(dat):
                        if row[7] == 'solvent' and row[1] not in EMPTY_VALUE and re.search(and_pat, row[1]):
                            del lp.data[idx]
                            break

                    entity_id_remap = {}
                    if len(entity_dict) > 1:
                        for idx, row in enumerate(dat):
                            if row[1] not in EMPTY_VALUE and row[5] not in EMPTY_VALUE:
                                key = row[1].lower()
                                entity_id = int(row[5]) if isinstance(row[5], str) else row[5]
                                if entity_id in entity_dict:
                                    entity = entity_dict[entity_id]
                                    if cross_check_entity(key, entity):
                                        continue
                                    for k, entity in entity_dict.items():
                                        if k != entity_id and cross_check_entity(key, entity):
                                            entity_id_remap[entity_id] = k
                                            if len(entity_dict) == 2:
                                                entity_id_remap[k] = entity_id

                    for idx, row in enumerate(dat):
                        if row[0] in EMPTY_VALUE:
                            lp.data[idx][id_col] = idx + 1
                        entity_id = None
                        if row[5] not in EMPTY_VALUE:
                            entity_id = _entity_id = int(row[5]) if isinstance(row[5], str) else row[5]
                            entity_id = entity_id_remap.get(entity_id, entity_id)
                            if entity_id in entity_dict:
                                entity = entity_dict[entity_id]
                                if row[1] in EMPTY_VALUE:
                                    lp.data[idx][mol_common_name_col] = entity['name']
                                lp.data[idx][entity_id_col] = entity_id
                                lp.data[idx][assembly_id_col] = entity['assembly_id']
                                lp.data[idx][assembly_label_col] = entity['assembly_label']
                                lp.data[idx][entity_label_col] = f"${entity['sf_framecode']}"
                                lp.data[idx][type_col] = entity['sample_type']
                                if row[8] in EMPTY_VALUE and (row[9] in EMPTY_VALUE or row[10] in EMPTY_VALUE):
                                    has_mand_concentration_val = False
                                has_poly_entity = True
                                touched_entity_id.append(entity_id)
                        elif row[1] not in EMPTY_VALUE:
                            if redundant_solvent_pat.match(row[1]):
                                h = redundant_solvent_pat.search(row[1])
                                if h[1] == h[2]:
                                    row[1] = lp.data[idx][mol_common_name_col] = h[1]
                            key = row[1].lower()
                            if key not in ('nacl', 'kcl', 'cacl2', 'zncl2', 'mgcl2',
                                           'ca+2', 'zn+2', 'mg+2', 'ca2+', 'zn2+', 'mg2+')\
                               and 'chloride' not in key:
                                entity_id = next((k for k, v in entity_dict.items()
                                                  if v['sf_framecode'] == row[1] or v['name'] == row[1]
                                                  and k not in touched_entity_id), None)
                                if entity_id is None:
                                    for k, entity in entity_dict.items():
                                        if cross_check_entity(key, entity) and k not in touched_entity_id:
                                            entity_id = k
                                            break
                            if entity_id is None and any(word in key for word in PROTEIN_RELATED_WORDS):
                                can_entity_id = [k for k, v in entity_dict.items()
                                                 if k not in touched_entity_id and v['sample_type'] in ('protein', 'peptide')]
                                if len(can_entity_id) == 1:
                                    entity_id = can_entity_id[0]
                            if entity_id is None and any(word in key for word in DNA_RELATED_WORDS):
                                can_entity_id = [k for k, v in entity_dict.items()
                                                 if k not in touched_entity_id and v['sample_type'] == 'DNA']
                                if len(can_entity_id) == 1:
                                    entity_id = can_entity_id[0]
                            if entity_id is None and any(word in key for word in RNA_RELATED_WORDS):
                                can_entity_id = [k for k, v in entity_dict.items()
                                                 if k not in touched_entity_id and v['sample_type'] == 'RNA']
                                if len(can_entity_id) == 1:
                                    entity_id = can_entity_id[0]
                            if entity_id is not None:
                                if entity_id not in touched_entity_id:
                                    touched_entity_id.append(entity_id)
                                entity = entity_dict[entity_id]
                                lp.data[idx][assembly_id_col] = entity['assembly_id']
                                lp.data[idx][assembly_label_col] = entity['assembly_label']
                                lp.data[idx][entity_id_col] = entity_id
                                lp.data[idx][entity_label_col] = f"${entity['sf_framecode']}"
                                lp.data[idx][type_col] = entity['sample_type']
                                if row[8] in EMPTY_VALUE and (row[9] in EMPTY_VALUE or row[10] in EMPTY_VALUE):
                                    has_mand_concentration_val = False
                                if entity['sample_type'] in ('protein', 'peptide', 'DNA', 'RNA', 'DNA/RNA hybrid'):
                                    has_poly_entity = True
                            elif row[1] in solvent_system:
                                if row[1] in solvent_isotope:
                                    lp.data[idx][isotopic_labeling_col] = solvent_isotope[row[1]]
                                lp.data[idx][type_col] = 'solvent'
                                if row[8] in EMPTY_VALUE or row[11] in EMPTY_VALUE or row[11] == '%':
                                    lp.data[idx][concentration_val_col] = solvent_system[row[1]]
                                if row[11] in EMPTY_VALUE:
                                    lp.data[idx][concentration_val_units_col] = '%'
                                if row[1] not in solvent_in_sample_loop:
                                    solvent_in_sample_loop.append(row[1])
                            else:
                                mol_common_name = row[1].lower()
                                if mol_common_name in ('nan3', 'nano3') or 'azide' in mol_common_name:
                                    lp.data[idx][type_col] = 'cytocide'
                                elif mol_common_name in ('k-pi', 'na-pi', 'na-p', 'kpi', 'napi', 'nap')\
                                        or 'buffer' in mol_common_name\
                                        or 'acetate' in mol_common_name\
                                        or 'acetic' in mol_common_name\
                                        or 'ch3coo' in mol_common_name\
                                        or 'citric' in mol_common_name\
                                        or 'hepes' in mol_common_name\
                                        or 'kphos' in mol_common_name\
                                        or 'mes' in mol_common_name\
                                        or 'mops' in mol_common_name\
                                        or ('morpholino' in mol_common_name and 'propane' in mol_common_name
                                            and 'sulfonic' in mol_common_name)\
                                        or 'naphos' in mol_common_name\
                                        or 'pbs' in mol_common_name\
                                        or 'phosphate' in mol_common_name\
                                        or 'pipes' in mol_common_name\
                                        or 'po4' in mol_common_name\
                                        or 'tricine' in mol_common_name\
                                        or 'tris' in mol_common_name:
                                    lp.data[idx][type_col] = 'buffer'
                                    lp.data[idx][isotopic_labeling_col] = '[U-2H]'\
                                        if mol_common_name.startswith('d-')\
                                        or 'deuterate' in mol_common_name\
                                        or re.match(r'd\d+-.*', mol_common_name) else 'natural abundance'
                                    has_buffer = True
                                elif mol_common_name in ('nacl', 'kcl', 'na2so4',
                                                         'cacl2', 'zncl2', 'mgcl2',
                                                         'ca+2', 'zn+2', 'mg+2',
                                                         'ca2+', 'zn2+', 'mg2+',
                                                         'sodium sulfate')\
                                        or 'chloride' in mol_common_name:
                                    lp.data[idx][type_col] = 'salt'
                                    lp.data[idx][isotopic_labeling_col] = 'natural abundance'
                                    has_salt = True
                                elif 'dtt' in mol_common_name or 'dithiothreitol' in mol_common_name\
                                        or 'dtba' in mol_common_name or 'dithiobutylamine' in mol_common_name\
                                        or 'cysteamine' in mol_common_name or 'betame' in mol_common_name\
                                        or 'mercaptoethanol' in mol_common_name or 'beta-me' in mol_common_name\
                                        or 'betame' in mol_common_name or 'tcep' in mol_common_name:
                                    lp.data[idx][type_col] = 'reducing agent'
                                    lp.data[idx][isotopic_labeling_col] = '[U-2H]'\
                                        if mol_common_name.startswith('d-')\
                                        or 'deuterate' in mol_common_name\
                                        or re.match(r'd\d+-.*', mol_common_name) else 'natural abundance'
                                    has_reducing_agent = True
                                elif mol_common_name in ('edta', 'egta'):
                                    lp.data[idx][type_col] = 'chelating agent'
                                    lp.data[idx][isotopic_labeling_col] = 'natural abundance'
                                    has_chelating_agent = True
                                elif mol_common_name in ('dss', 'tsp'):
                                    lp.data[idx][type_col] = 'internal reference'
                                    lp.data[idx][isotopic_labeling_col] = 'natural abundance'
                                    default_internal_reference = mol_common_name.upper()
                                elif 'pyridostatin' in mol_common_name or mol_common_name in ('pds', 'pypds'):
                                    lp.data[idx][type_col] = 'G-quadruplex stabilizing agent'
                                elif 'lps' in mol_common_name or 'lipopolysaccharide' in mol_common_name:
                                    lp.data[idx][type_col] = 'bacterial outer membrane'
                                elif 'bicelle' in mol_common_name or 'phage' in mol_common_name or mol_common_name == 'pf1':
                                    lp.data[idx][type_col] = 'molecular alignment inductor'
                                elif mol_common_name in ('dpc', 'dpc-d38', 'lpg', 'lppg',
                                                         'sds', 'sodium dodecyl sulfate')\
                                        or 'dodecylphosphocholine' in mol_common_name\
                                        or 'lysophosph' in mol_common_name\
                                        or 'micelle' in mol_common_name:
                                    lp.data[idx][type_col] = 'micelles'
                                elif mol_common_name in ('chaps', 'chapso', 'dhpc', 'popc', 'np-40')\
                                        or ('deuterate' in mol_common_name
                                            and ('dhpc' in mol_common_name) or 'popc' in mol_common_name)\
                                        or mol_common_name.startswith('zwittergent'):
                                    lp.data[idx][type_col] = 'detergent'
                                    lp.data[idx][isotopic_labeling_col] = '[U-2H]'\
                                        if mol_common_name.startswith('d-')\
                                        or 'deuterate' in mol_common_name\
                                        or re.match(r'd\d+-.*', mol_common_name) else 'natural abundance'
                                elif 'phosph' in mol_common_name\
                                        and ('ylethanolamine' in mol_common_name or 'ylcholine' in mol_common_name
                                             or 'ylserine' in mol_common_name or 'ylinositol' in mol_common_name
                                             or 'lipid' in mol_common_name):
                                    lp.data[idx][type_col] = 'phospholipid'
                                elif mol_common_name in ('dmpc', 'dopc', 'dppc', 'popg', 'popg-na') or 'lipids' in mol_common_name:
                                    lp.data[idx][type_col] = 'phospholipid'
                                elif mol_common_name == 'dmso-d6':
                                    lp.data[idx][isotopic_labeling_col] = '[U-2H]'
                                    lp.data[idx][type_col] = 'solvent'
                                    solvent_in_sample_loop.append('DMSO')
                                elif mol_common_name == 'dmso':
                                    lp.data[idx][type_col] = 'solvent'
                                    solvent_in_sample_loop.append('DMSO')
                                elif mol_common_name == 'glycerol':
                                    lp.data[idx][type_col] = 'solvent'
                                    solvent_in_sample_loop.append('Glycerol')
                                elif mol_common_name == 'glycerol-d8':
                                    lp.data[idx][type_col] = 'solvent'
                                    solvent_in_sample_loop.append('Glycerol')
                                    lp.data[idx][isotopic_labeling_col] = '[U-2H]'
                                elif mol_common_name in ('tfe', 'trifluoroethanol'):
                                    lp.data[idx][type_col] = 'solvent'
                                    solvent_in_sample_loop.append('TFE')
                                elif 'deuterate' in mol_common_name and 'd2o' in mol_common_name:
                                    lp.data[idx][mol_common_name_col] = 'D2O'
                                    lp.data[idx][type_col] = 'solvent'
                                    solvent_in_sample_loop.append('D2O')

                        if row[2] not in EMPTY_VALUE:
                            if is_natural_abundance(row[2]):
                                lp.data[idx][isotopic_labeling_col] = 'natural abundance'
                            if ',' in row[2] and ']-' not in row[2]:
                                isotopic_labelings = []
                                effective_labeling = True
                                for txt in row[2].split(','):
                                    txt = txt.strip().lstrip('[').rstrip(']')
                                    if ';' in txt:
                                        for _txt in txt.split(';'):
                                            _txt = _txt.strip()
                                            iso_name = next((name for name in ISOTOPE_NAMES_OF_NMR_OBS_NUCS if name in _txt), None)
                                            if iso_name is None:
                                                effective_labeling = False
                                                break
                                            if '%' in _txt:
                                                try:
                                                    g = percent_value_pat.search(_txt).groups()
                                                    p = float(g[0]) if '.' in g[0] else int(g[0])
                                                    if 0 < p <= 100:
                                                        isotopic_labelings.append(f'U-{p}% {iso_name}')
                                                    else:
                                                        effective_labeling = False
                                                        break
                                                except AttributeError:
                                                    effective_labeling = False
                                                    break
                                            else:
                                                isotopic_labelings.append(f'U-{iso_name}')
                                        if not effective_labeling:
                                            break
                                    else:
                                        iso_name = next((name for name in ISOTOPE_NAMES_OF_NMR_OBS_NUCS if name in txt), None)
                                        if iso_name is None:
                                            effective_labeling = False
                                            break
                                        if '%' in txt:
                                            try:
                                                g = percent_value_pat.search(txt).groups()
                                                p = float(g[0]) if '.' in g[0] else int(g[0])
                                                if 0 < p <= 100:
                                                    isotopic_labelings.append(f'U-{p}% {iso_name}')
                                                else:
                                                    effective_labeling = False
                                                    break
                                            except AttributeError:
                                                effective_labeling = False
                                                break
                                        else:
                                            isotopic_labelings.append(f'U-{iso_name}')
                                if effective_labeling:
                                    lp.data[idx][isotopic_labeling_col] = f'[{"; ".join(isotopic_labelings)}]'

                        else:
                            lp.data[idx][isotopic_labeling_col] = 'natural abundance'

                        if row[8] not in EMPTY_VALUE and isinstance(row[8], str):
                            for delimiter in ('-', '/', ':', ',', '~'):
                                if delimiter not in row[8]:
                                    continue
                                code_count = row[8].count(delimiter)
                                if code_count == 1:
                                    try:
                                        _val = row[8].split(delimiter)
                                        _val_0 = float(_val[0])
                                        _val_1 = float(_val[1])
                                        _val = [float(_val[0]), float(_val[1])]
                                        lp.data[idx][concentration_val_col] = '.'
                                        lp.data[idx][concentration_val_min_col] = min(_val)
                                        lp.data[idx][concentration_val_max_col] = max(_val)
                                        break
                                    except ValueError:
                                        pass

                    if len(lp) == 1 and len(entity_dict) == 1 and not has_poly_entity and dat[0][7] in EMPTY_VALUE:
                        for entity_id, entity in entity_dict.items():
                            if entity['sample_type'] in ('protein', 'peptide', 'DNA', 'RNA', 'DNA/RNA hybrid'):
                                if dat[0][1] in EMPTY_VALUE:
                                    lp.data[0][mol_common_name_col] = entity['name']
                                lp.data[0][entity_id_col] = entity_id
                                lp.data[0][assembly_id_col] = entity['assembly_id']
                                lp.data[0][assembly_label_col] = entity['assembly_label']
                                lp.data[0][entity_label_col] = f"${entity['sf_framecode']}"
                                lp.data[0][type_col] = entity['sample_type']
                                if row[8] in EMPTY_VALUE and (row[9] in EMPTY_VALUE or row[10] in EMPTY_VALUE):
                                    has_mand_concentration_val = False
                                has_poly_entity = True

                    cur_id = len(lp) + 1

                    if not has_poly_entity:
                        for entity_id, entity in entity_dict.items():
                            if entity['sample_type'] in ('protein', 'peptide', 'DNA', 'RNA', 'DNA/RNA hybrid'):
                                row = [None] * len(lp.tags)
                                row[id_col] = cur_id
                                row[mol_common_name_col] = entity['name']
                                row[isotopic_labeling_col] = 'natural abundance'
                                if entity_id in isotope_nums_per_entity:
                                    isotopic_labeling = []
                                    for isotope_number in sorted(list(isotope_nums_per_entity[entity_id])):
                                        if isotope_number == 1:
                                            continue
                                        if isotope_number not in ALLOWED_ISOTOPE_NUMBERS:
                                            continue
                                        type_symbol = next(k for k, v in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS.items()
                                                           if isotope_number in v)
                                        isotopic_labeling.append(f'U-{isotope_number}{type_symbol}')
                                    if len(isotopic_labeling) > 0:
                                        row[isotopic_labeling_col] = f"[{'; '.join(isotopic_labeling)}]"
                                row[assembly_id_col] = entity['assembly_id']
                                row[assembly_label_col] = entity['assembly_label']
                                row[entity_id_col] = entity_id
                                row[entity_label_col] = f"${entity['sf_framecode']}"
                                row[type_col] = entity['sample_type']
                                row[concentration_val_col] = '?'
                                row[concentration_val_units_col] = 'mM'
                                row[sample_id_col] = sample_id
                                row[entry_id_col] = self.__reg.entry_id

                                has_mand_concentration_val = False
                                cur_id += 1

                    for entity_id, entity in entity_dict.items():
                        if entity_id in touched_entity_id:
                            continue
                        if entity['sample_type'] in ('DNA', 'RNA', 'DNA/RNA hybrid'):
                            _entity_id =\
                                next((_entity_id for _entity_id in touched_entity_id
                                      if entity_dict[_entity_id]['sample_type'] in ('DNA', 'RNA', 'DNA/RNA hybrid')), None)
                            if _entity_id is not None:
                                _row = next((_row for _row in lp.data
                                             if _row[entity_id_col] not in EMPTY_VALUE
                                             and ((isinstance(_row[entity_id_col], str)
                                                   and int(_row[entity_id_col]) == _entity_id)
                                                  or (isinstance(_row[entity_id_col], int)
                                                      and _row[entity_id_col] == _entity_id))), None)
                                if _row is not None:
                                    row = copy.copy(_row)
                                    row[id_col] = cur_id
                                    row[mol_common_name_col] = entity['name']
                                    row[entity_id_col] = entity_id
                                    row[entity_label_col] = f"${entity['sf_framecode']}"
                                    if has_mand_concentration_val:
                                        lp.add_data(row)
                                    cur_id += 1

                    dat = lp.get_tag(tags)

                    for idx, row in enumerate(dat):
                        if row[1] == 'H2O':
                            if row[2] != 'natural abundance':
                                lp.data[idx][isotopic_labeling_col] = 'natural abundance'
                        elif row[1] == 'D2O':
                            if not row[2].endswith('2H]'):
                                lp.data[idx][isotopic_labeling_col] = '[U-2H]'

                    if is_single_sample_loop:

                        for row in dat:
                            if row[5] in EMPTY_VALUE:
                                continue
                            if not is_natural_abundance(row[2]):
                                continue
                            entity_id = row[5]
                            if isinstance(entity_id, str):
                                entity_id = int(entity_id)
                            if entity_id not in isotope_nums_per_entity:
                                continue
                            _isotope_nums = isotope_nums_per_entity[entity_id]
                            if 13 in _isotope_nums or 15 in _isotope_nums:
                                has_hvy_shifts_with_natural_abundance = True

                    if is_single_protein_ligand_complex and has_mand_concentration_val:
                        ref_concentration_val = ref_common_name = None
                        for row in lp:
                            if (isinstance(row[entity_id_col], int) and row[entity_id_col] == polymer_common_types[0]['protein'])\
                               or (isinstance(row[entity_id_col], str) and row[entity_id_col].isdigit()
                                   and int(row[entity_id_col]) == polymer_common_types[0]['protein']):
                                if ref_common_name is not None:
                                    ref_concentration_val = ref_common_name = None
                                    break
                                ref_concentration_val = f'{row[concentration_val_col]} '\
                                    f'{row[concentration_val_min_col]} '\
                                    f'{row[concentration_val_max_col]} '\
                                    f'{row[concentration_val_units_col]}'
                                ref_common_name = row[mol_common_name_col]

                        if ref_concentration_val is not None:
                            can_ligand_idx = None
                            for idx, row in enumerate(lp):
                                if row[entity_id_col] in EMPTY_VALUE:
                                    test_concentration_val = f'{row[concentration_val_col]} '\
                                        f'{row[concentration_val_min_col]} '\
                                        f'{row[concentration_val_max_col]} '\
                                        f'{row[concentration_val_units_col]}'
                                    if test_concentration_val == ref_concentration_val\
                                       and not any(word in row[mol_common_name_col] for word in PROTEIN_RELATED_WORDS)\
                                       and not row[mol_common_name_col].startswith(ref_common_name):
                                        if isinstance(can_ligand_idx, int):
                                            can_ligand_idx = None
                                            break
                                        can_ligand_idx = idx
                            if can_ligand_idx is not None:
                                entity = next(entity for entity_id, entity in entity_dict.items()
                                              if entity_id == non_polymer_types[0]['ligand'])
                                row = lp.data[can_ligand_idx]
                                row[assembly_id_col] = entity['assembly_id']
                                row[assembly_label_col] = entity['assembly_label']
                                row[entity_id_col] = non_polymer_types[0]['ligand']
                                row[entity_label_col] = f"${entity['sf_framecode']}"
                                row[type_col] = 'ligand'

                    for solvent in solvent_system:
                        if solvent in solvent_in_sample_loop:
                            continue
                        if solvent not in appended_solvent_in_sample_loop and not not_title(solvent):
                            appended_solvent_in_sample_loop.append(solvent)
                        row = [None] * len(lp.tags)
                        row[id_col] = cur_id
                        row[mol_common_name_col] = solvent
                        if solvent in solvent_isotope:
                            row[isotopic_labeling_col] = solvent_isotope[solvent]
                        row[type_col] = 'solvent'
                        row[concentration_val_col] = solvent_system[solvent]
                        row[concentration_val_units_col] = '%'
                        row[sample_id_col] = sample_id
                        row[entry_id_col] = self.__reg.entry_id
                        if not has_mand_concentration_val:
                            row[id_col] = len(lp) + 1
                        lp.add_data(row)
                        cur_id += 1

                        if solvent in solvent_isotope and '?' in solvent_isotope[solvent]:
                            has_perdeuteration_sol = True

                    if has_hvy_shifts_with_natural_abundance:

                        dat = lp.get_tag(tags)

                        for row in dat:
                            if row[5] in EMPTY_VALUE:
                                continue
                            if not is_natural_abundance(row[2]):
                                continue
                            entity_id = row[5]
                            if isinstance(entity_id, str):
                                entity_id = int(entity_id)
                            if entity_id not in isotope_nums_per_entity:
                                continue

                    if has_perdeuteration_sol:

                        dat = lp.get_tag(tags)

                        for idx, row in enumerate(dat):
                            if row[2] not in EMPTY_VALUE and '?' in row[2] and row[7] == 'solvent':
                                lp.data[idx][isotopic_labeling_col] = row[2].replace('?% ', '')

                    # sort sample loop

                    lp = sf.get_loop(lp_category)

                    _lp = pynmrstar.Loop.from_scratch(lp_category)

                    _lp.add_tag(lp.tags)

                    cur_id = 1

                    tags = ['Entity_ID', 'Type']

                    data = lp.get_tag(tags)

                    for sample_type in ['protein', 'peptide', 'DNA', 'RNA', 'DNA/RNA hybrid', 'carbohydrate', 'polysaccharide']:
                        idx_dict = []
                        for idx, row in enumerate(data):
                            if row[1] not in EMPTY_VALUE and row[1] == sample_type:
                                entity_id = 9999
                                if row[0] not in EMPTY_VALUE:
                                    if isinstance(row[0], str):
                                        if row[0].isdigit():
                                            entity_id = int(row[0])
                                    else:
                                        entity_id = row[0]
                                idx_dict.append({'entity_id': entity_id, 'row_idx': idx})
                        if len(idx_dict) > 0:
                            for d in sorted(idx_dict, key=itemgetter('entity_id')):
                                row = copy.copy(lp.data[d['row_idx']])
                                row[id_col] = cur_id
                                _lp.add_data(row)
                                cur_id += 1

                    tags = ['Entity_ID', 'Type', 'Mol_common_name', 'Isotopic_labeling']

                    data = lp.get_tag(tags)

                    for idx, row in enumerate(data):
                        if row[0] in EMPTY_VALUE and row[1] in EMPTY_VALUE:
                            entity_id =\
                                next((_row[0] for _row in data
                                      if len([term in _row[2] for term in row[2].split()]) > 2 and row[3] != _row[3]), None)
                            if entity_id is None:
                                continue
                            if isinstance(entity_id, str):
                                entity_id = int(entity_id)
                            entity_id = entity_id_remap.get(entity_id, entity_id)
                            if entity_id is not None and entity_id in entity_dict:
                                entity = entity_dict[entity_id]
                                lp.data[idx][assembly_id_col] = entity['assembly_id']
                                lp.data[idx][assembly_label_col] = entity['assembly_label']
                                lp.data[idx][entity_id_col] = entity_id
                                lp.data[idx][entity_label_col] = f"${entity['sf_framecode']}"
                                lp.data[idx][type_col] = entity['sample_type']
                                _lp.add_data(lp.data[idx])
                                cur_id += 1

                    tags = ['Entity_ID', 'Type']

                    data = lp.get_tag(tags)

                    unsorted_sample_types = []
                    for idx, row in enumerate(data):
                        if row[1] in EMPTY_VALUE:
                            if '' not in unsorted_sample_types:
                                unsorted_sample_types.append('')
                            continue
                        if row[1] in ('protein', 'peptide', 'DNA', 'RNA', 'DNA/RNA hybrid',
                                      'carbohydrate', 'polysaccharide',
                                      'reducing agent', 'chelating agent', 'salt', 'buffer',
                                      'phospholipid', 'internal reference', 'solvent'):
                            continue
                        if row[1] not in unsorted_sample_types:
                            unsorted_sample_types.append(row[1])

                    if len(unsorted_sample_types) > 0:
                        for sample_type in sorted(unsorted_sample_types):
                            idx_dict = []
                            for idx, row in enumerate(data):
                                if row[1] not in EMPTY_VALUE and row[1] == sample_type:
                                    entity_id = 9999
                                    if row[0] not in EMPTY_VALUE:
                                        if isinstance(row[0], str):
                                            if row[0].isdigit():
                                                entity_id = int(row[0])
                                        else:
                                            entity_id = row[0]
                                    idx_dict.append({'entity_id': entity_id, 'row_idx': idx})
                            if len(idx_dict) > 0:
                                for d in sorted(idx_dict, key=itemgetter('entity_id')):
                                    row = copy.copy(lp.data[d['row_idx']])
                                    row[id_col] = cur_id
                                    _lp.add_data(row)
                                    cur_id += 1

                    for sample_type in ['reducing agent', 'chelating agent', 'salt', 'buffer',
                                        'phospholipid', 'internal reference', 'solvent']:
                        idx_dict = []
                        for idx, row in enumerate(data):
                            if row[1] not in EMPTY_VALUE and row[1] == sample_type:
                                entity_id = 9999
                                if row[0] not in EMPTY_VALUE:
                                    if isinstance(row[0], str):
                                        if row[0].isdigit():
                                            entity_id = int(row[0])
                                    else:
                                        entity_id = row[0]
                                idx_dict.append({'entity_id': entity_id, 'row_idx': idx})
                        if len(idx_dict) > 0:
                            for d in sorted(idx_dict, key=itemgetter('entity_id')):
                                row = copy.copy(lp.data[d['row_idx']])
                                row[id_col] = cur_id
                                _lp.add_data(row)
                                cur_id += 1
                    # """
                    # for idx, row in enumerate(data):
                    #     if row[1] in EMPTY_VALUE:
                    #         _row = lp.data[idx]
                    #         _row[type_col] = 'na (yet not decided)'
                    #         _row[id_col] = cur_id
                    #         _lp.add_data(_row)
                    #         cur_id += 1
                    # """
                    if has_mand_concentration_val:
                        del sf[lp]

                        sf.add_loop(_lp)

                except (ValueError, KeyError):
                    continue

            if ionic_strength is not None and not has_salt and not has_buffer:
                if not has_chelating_agent:
                    sf = master_entry.get_saveframes_by_category(sf_category)[0]

                    lp_category = '_Sample_component'

                    try:

                        lp = sf.get_loop(lp_category)

                        _lp = pynmrstar.Loop.from_scratch(lp_category)

                        tags = ['ID',
                                'Mol_common_name',
                                'Isotopic_labeling',
                                'Type',
                                'Concentration_val',
                                'Concentration_val_min',
                                'Concentration_val_max',
                                'Concentration_val_units',
                                'Concentration_val_err']

                        _lp.add_tag(tags)

                        data = lp.get_tag(tags)

                        for row in data:
                            _lp.add_data(row)

                        cur_id = len(data) + 1

                        row = ['?'] * len(tags)
                        row[0] = cur_id
                        row[3] = 'salt'
                        row[7] = 'mM'
                        cur_id += 1

                        _lp.add_data(row)

                        if not has_buffer:
                            row = ['?'] * len(tags)
                            row[0] = cur_id
                            row[3] = 'buffer'
                            row[7] = 'mM'

                            _lp.add_data(row)

                    except KeyError:
                        pass

            elif pH is not None and not has_buffer:
                sf = master_entry.get_saveframes_by_category(sf_category)[0]

                lp_category = '_Sample_component'

                try:

                    lp = sf.get_loop(lp_category)

                    _lp = pynmrstar.Loop.from_scratch(lp_category)

                    tags = ['ID',
                            'Mol_common_name',
                            'Isotopic_labeling',
                            'Type',
                            'Concentration_val',
                            'Concentration_val_min',
                            'Concentration_val_max',
                            'Concentration_val_units',
                            'Concentration_val_err']

                    _lp.add_tag(tags)

                    data = lp.get_tag(tags)

                    for row in data:
                        _lp.add_data(row)

                    cur_id = len(data) + 1

                    row = ['?'] * len(tags)
                    row[0] = cur_id
                    row[3] = 'buffer'
                    row[7] = 'mM'

                    _lp.add_data(row)

                except KeyError:
                    pass

        sf_category = 'spectral_peak_list'

        # resolve duplication of spectral peak lists

        if sf_category in self.__reg.sf_category_list and self.__reg.internal_mode:

            update_sf_name = {}

            tags_2d = ['Volume', 'Height', 'Position_1', 'Position_2', 'Auth_asym_ID_1']
            tags_3d = ['Volume', 'Height', 'Position_1', 'Position_2', 'Position_3', 'Auth_asym_ID_1']
            tags_4d = ['Volume', 'Height', 'Position_1', 'Position_2', 'Position_3', 'Position_4', 'Auth_asym_ID_1']

            tags_pk_char = ['Peak_ID', 'Spectral_dim_ID', 'Chem_shift_val']
            tags_pk_gen_char = ['Peak_ID', 'Intensity_val', 'Measurement_method']

            sp_info = {}
            empty_sf_framecodes = []

            truncated = False

            for idx, sf in enumerate(master_entry.get_saveframes_by_category(sf_category)):
                sf_framecode = get_first_sf_tag(sf, 'Sf_framecode')

                sp_dim = {}

                try:

                    _lp_category = '_Spectral_dim'
                    _lp = sf.get_loop(_lp_category)
                    _dat = _lp.get_tag(['ID', 'Atom_type'])
                    for _row in _dat:
                        sp_dim[int(_row[0]) if isinstance(_row[0], str) else _row[0]] = _row[1]

                except KeyError:
                    pass

                lp_category = '_Peak_row_format'

                try:

                    lp = sf.get_loop(lp_category)

                    if len(lp) > 0:
                        sp_info[idx] = {'sf_framecode': sf_framecode,
                                        'num_of_dim': get_first_sf_tag(sf, 'Number_of_spectral_dimensions'),
                                        'class': get_first_sf_tag(sf, 'Experiment_class').replace('?', ''),
                                        'type': get_first_sf_tag(sf, 'Experiment_type'),
                                        'spectral_dim': sp_dim,
                                        'size': len(lp)}

                        num_of_dim = sp_info[idx]['num_of_dim']
                        if isinstance(num_of_dim, str):
                            num_of_dim = sp_info[idx]['num_of_dim'] = int(num_of_dim)

                        has_volume = has_height = has_assign = False
                        assigns = details = 0
                        signature = []

                        if num_of_dim == 2:
                            dat = lp.get_tag(tags_2d)

                            for _idx, row in enumerate(dat, start=1):

                                if _idx <= MAX_ROWS_TO_CHECK_SPECTRAL_PEAK_IDENTITY:
                                    volume = row[0]
                                    if volume in EMPTY_VALUE:
                                        volume = None
                                    else:
                                        volume = float(volume)
                                        has_volume = True

                                    height = row[1]
                                    if height in EMPTY_VALUE:
                                        height = None
                                    else:
                                        height = float(height)
                                        has_height = True

                                    position = {float(row[2]), float(row[3])}

                                    signature.append({'position': position,
                                                      'volume': volume,
                                                      'height': height})

                                if row[4] not in EMPTY_VALUE:
                                    has_assign = True

                            if has_assign:
                                for row in dat:
                                    if row[4] not in EMPTY_VALUE:
                                        assigns += 1

                                if 'Details' in lp.tags:
                                    dat = lp.get_tag(['Details'])
                                    for row in dat:
                                        if row not in EMPTY_VALUE:
                                            if ' -> ' in row:
                                                details -= 1
                                            else:
                                                details += 1

                        elif num_of_dim == 3:
                            dat = lp.get_tag(tags_3d)

                            for _idx, row in enumerate(dat, start=1):

                                if _idx <= MAX_ROWS_TO_CHECK_SPECTRAL_PEAK_IDENTITY:
                                    volume = row[0]
                                    if volume in EMPTY_VALUE:
                                        volume = None
                                    else:
                                        volume = float(volume)

                                    height = row[1]
                                    if height in EMPTY_VALUE:
                                        height = None
                                    else:
                                        height = float(height)

                                    position = {float(row[2]), float(row[3]), float(row[4])}

                                    signature.append({'position': position,
                                                      'volume': volume,
                                                      'height': height})

                                if row[5] not in EMPTY_VALUE:
                                    has_assign = True
                                    if _idx > MAX_ROWS_TO_CHECK_SPECTRAL_PEAK_IDENTITY:
                                        break

                            if has_assign:
                                for row in dat:
                                    if row[5] not in EMPTY_VALUE:
                                        assigns += 1

                            if has_assign and 'Details' in lp.tags:
                                dat = lp.get_tag(['Details'])
                                for row in dat:
                                    if row not in EMPTY_VALUE:
                                        if ' -> ' in row:
                                            details -= 1
                                        else:
                                            details += 1

                        elif num_of_dim == 4:
                            dat = lp.get_tag(tags_4d)

                            for _idx, row in enumerate(dat, start=1):

                                if _idx <= MAX_ROWS_TO_CHECK_SPECTRAL_PEAK_IDENTITY:
                                    volume = row[0]
                                    if volume in EMPTY_VALUE:
                                        volume = None
                                    else:
                                        volume = float(volume)

                                    height = row[1]
                                    if height in EMPTY_VALUE:
                                        height = None
                                    else:
                                        height = float(height)

                                    position = {float(row[2]), float(row[3]), float(row[4]), float(row[5])}

                                    signature.append({'position': position,
                                                      'volume': volume,
                                                      'height': height})

                                if row[6] not in EMPTY_VALUE:
                                    has_assign = True
                                    if _idx > MAX_ROWS_TO_CHECK_SPECTRAL_PEAK_IDENTITY:
                                        break

                            if has_assign:
                                for row in dat:
                                    if row[6] not in EMPTY_VALUE:
                                        assigns += 1

                            if has_assign and 'Details' in lp.tags:
                                dat = lp.get_tag(['Details'])
                                for row in dat:
                                    if row not in EMPTY_VALUE:
                                        if ' -> ' in row:
                                            details -= 1
                                        else:
                                            details += 1

                        sp_info[idx]['has_volume'] = has_volume
                        sp_info[idx]['has_height'] = has_height
                        sp_info[idx]['has_assign'] = has_assign
                        sp_info[idx]['assigns'] = assigns
                        sp_info[idx]['details'] = details
                        sp_info[idx]['signature'] = signature

                except KeyError:

                    try:

                        pk_char_category = '_Peak_char'
                        pk_char = sf.get_loop(pk_char_category)

                        pk_gen_char_category = '_Peak_general_char'
                        pk_gen_char = sf.get_loop(pk_gen_char_category)

                        if len(pk_char) > 0 and len(pk_gen_char) > 0:
                            sp_info[idx] = {'sf_framecode': sf_framecode,
                                            'num_of_dim': get_first_sf_tag(sf, 'Number_of_spectral_dimensions'),
                                            'class': get_first_sf_tag(sf, 'Experiment_class').replace('?', ''),
                                            'type': get_first_sf_tag(sf, 'Experiment_type'),
                                            'spectral_dim': sp_dim,
                                            'size': len(lp)}

                            num_of_dim = sp_info[idx]['num_of_dim']
                            if isinstance(num_of_dim, str):
                                num_of_dim = sp_info[idx]['num_of_dim'] = int(num_of_dim)

                            has_volume = has_height = has_assign = False
                            signature = []

                            dat_pk_char = pk_char.get_tag(tags_pk_char)
                            dat_pk_gen_char = pk_gen_char.get_tag(tags_pk_gen_char)

                            peak_ids = []
                            for row_pk_char in dat_pk_char:
                                if row_pk_char[0] not in peak_ids:
                                    peak_ids.append(row_pk_char[0])
                                    if len(peak_ids) == MAX_ROWS_TO_CHECK_SPECTRAL_PEAK_IDENTITY:
                                        break

                            position = [None] * num_of_dim

                            for row_pk_char in dat_pk_char:
                                if row_pk_char[0] in peak_ids:
                                    sp_dim_id = int(row_pk_char[1])
                                    position[sp_dim_id - 1] = float(row_pk_char[2])
                                    if all(pos is not None for pos in position):
                                        volume = next((float(row_pk_gen_char[1]) for row_pk_gen_char in dat_pk_gen_char
                                                       if row_pk_gen_char[0] == row_pk_char[0]
                                                       and row_pk_gen_char[2] == 'volume'), None)
                                        if volume is not None:
                                            has_volume = True
                                        height = next((float(row_pk_gen_char[1]) for row_pk_gen_char in dat_pk_gen_char
                                                       if row_pk_gen_char[0] == row_pk_char[0]
                                                       and row_pk_gen_char[2] == 'height'), None)
                                        if height is not None:
                                            has_height = True
                                        signature.append({'position': set(position),
                                                          'volume': volume,
                                                          'height': height})

                                        for _idx in range(num_of_dim):
                                            position[_idx] = None
                                else:
                                    break

                            sp_info[idx]['has_volume'] = has_volume
                            sp_info[idx]['has_height'] = has_height
                            sp_info[idx]['has_assign'] = False
                            assigns = details = 0
                            try:
                                pk_cs = sf.get_loop('_Assigned_peak_chem_shift')
                                dat_pk_cs = pk_cs.get_tag(['Auth_entity_ID'])
                                for row_pk_cs in dat_pk_cs:
                                    if row_pk_cs not in EMPTY_VALUE:
                                        sp_info[idx]['has_assign'] = True
                                        break
                                if sp_info[idx]['has_assign']:
                                    for row_pk_cs in dat_pk_cs:
                                        if row_pk_cs not in EMPTY_VALUE:
                                            assigns += 1
                                    if 'Details' in pk_cs.tags:
                                        dat_pk_cs = pk_cs.get_tag(['Details'])
                                        for row_pk_cs in dat_pk_cs:
                                            if row_pk_cs not in EMPTY_VALUE:
                                                if ' -> ' in row_pk_cs:
                                                    details -= 1
                                                else:
                                                    details += 1
                            except KeyError:
                                pass
                            sp_info[idx]['assigns'] = assigns
                            sp_info[idx]['details'] = details
                            sp_info[idx]['signature'] = signature

                    except KeyError:
                        if get_first_sf_tag(sf, 'Text_data') in EMPTY_VALUE:
                            empty_sf_framecodes.append(sf_framecode)

            if len(sp_info) > 1:
                dup_idx = set()
                res_idx = {}

                for idx1, idx2 in itertools.combinations(sp_info, 2):

                    if sp_info[idx1]['num_of_dim'] != sp_info[idx2]['num_of_dim']\
                       or sp_info[idx1]['spectral_dim'] != sp_info[idx2]['spectral_dim']\
                       or (not self.__reg.enforce_peak_row_format and sp_info[idx1]['size'] != sp_info[idx2]['size']):
                        continue

                    if not all(sig1['position'] == sig2['position']
                               for sig1, sig2 in zip(sp_info[idx1]['signature'],
                                                     sp_info[idx2]['signature'])):
                        continue

                    if not all(sig1['volume'] == sig2['volume']
                               for sig1, sig2 in zip(sp_info[idx1]['signature'],
                                                     sp_info[idx2]['signature']))\
                       and not all(sig1['height'] == sig2['height']
                                   for sig1, sig2 in zip(sp_info[idx1]['signature'],
                                                         sp_info[idx2]['signature'])):
                        continue

                    if sp_info[idx1]['has_assign'] != sp_info[idx2]['has_assign']:
                        if sp_info[idx1]['has_assign']:
                            dup_idx.add(idx2)
                            res_idx[idx2] = idx1
                        else:
                            dup_idx.add(idx1)
                            res_idx[idx1] = idx2
                        continue

                    if sp_info[idx1]['has_assign'] and sp_info[idx2]['has_assign']:
                        if sp_info[idx1]['assigns'] > sp_info[idx2]['assigns']:
                            dup_idx.add(idx2)
                            res_idx[idx2] = idx1
                            continue
                        if sp_info[idx1]['assigns'] < sp_info[idx2]['assigns']:
                            dup_idx.add(idx1)
                            res_idx[idx1] = idx2
                            continue
                        if sp_info[idx1]['details'] < sp_info[idx2]['details']:
                            dup_idx.add(idx2)
                            res_idx[idx2] = idx1
                            continue
                        if sp_info[idx1]['details'] > sp_info[idx2]['details']:
                            dup_idx.add(idx1)
                            res_idx[idx1] = idx2
                            continue

                    if sp_info[idx1]['has_volume'] != sp_info[idx2]['has_volume']:
                        if sp_info[idx1]['has_volume']:
                            dup_idx.add(idx2)
                            res_idx[idx2] = idx1
                        else:
                            dup_idx.add(idx1)
                            res_idx[idx1] = idx2
                        continue

                    if sp_info[idx1]['has_height'] != sp_info[idx2]['has_height']:
                        if sp_info[idx1]['has_height']:
                            dup_idx.add(idx2)
                            res_idx[idx2] = idx1
                        else:
                            dup_idx.add(idx1)
                            res_idx[idx1] = idx2
                        continue

                    if len(sp_info[idx1]['class']) != len(sp_info[idx2]['class']):
                        if len(sp_info[idx1]['class']) > len(sp_info[idx2]['class']):
                            dup_idx.add(idx2)
                            res_idx[idx2] = idx1
                        else:
                            dup_idx.add(idx1)
                            res_idx[idx1] = idx2
                        continue

                    if len(sp_info[idx1]['type']) != len(sp_info[idx2]['type']):
                        if len(sp_info[idx1]['type']) > len(sp_info[idx2]['type']):
                            dup_idx.add(idx2)
                            res_idx[idx2] = idx1
                        else:
                            dup_idx.add(idx1)
                            res_idx[idx1] = idx2
                        continue

                    if idx1 in dup_idx or idx2 in dup_idx:
                        continue

                    dup_idx.add(idx2)
                    res_idx[idx2] = idx1

                pk_name_pat = re.compile(r'D_\d{6,10}_nmr-peaks-upload_P(\d+)\.dat\.V(\d+)$')
                mr_name_pat = re.compile(r'D_\d{6,10}_mr-(\S+)_P(\d+)\.(\S+)\.V(\d+)$')

                for idx in sorted(list(dup_idx)):
                    sf_framecode = sp_info[idx]['sf_framecode']
                    sf = master_entry.get_saveframe_by_name(sf_framecode)
                    file_name = get_first_sf_tag(sf, 'Data_file_name')
                    file_name_ = retrieveOriginalFileName(file_name) if file_name not in EMPTY_VALUE else file_name

                    res_sf_framecode = sp_info[res_idx[idx]]['sf_framecode']

                    _sf = master_entry.get_saveframe_by_name(res_sf_framecode)
                    _file_name = get_first_sf_tag(_sf, 'Data_file_name')
                    _file_name_ = retrieveOriginalFileName(_file_name) if _file_name not in EMPTY_VALUE else _file_name

                    if _file_name in EMPTY_VALUE and _file_name_ not in EMPTY_VALUE:
                        set_sf_tag(_sf, 'Data_file_name', _file_name_)
                    elif _file_name_ in EMPTY_VALUE and file_name_ not in EMPTY_VALUE:
                        set_sf_tag(_sf, 'Data_file_name', file_name_)

                    if _file_name_ not in EMPTY_VALUE and file_name_ not in EMPTY_VALUE:
                        if _file_name_ != file_name_ and not mr_name_pat.match(file_name_)\
                           and pk_name_pat.match(_file_name_):
                            set_sf_tag(_sf, 'Data_file_name', file_name_)

                    update_sf_name[res_sf_framecode] = sf_framecode

                    warn = f'There are equivalent NMR spectral peak lists, {res_sf_framecode!r} vs {sf_framecode!r}. '\
                        f'For the sake of simplicity, saveframe {sf_framecode!r} derived from {file_name!r} was ignored.'

                    self.__reg.report.warning.appendDescription('redundant_mr_data',
                                                                {'file_name': file_name,
                                                                 'sf_framecode': sf_framecode,
                                                                 'description': warn})

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.perform() ++ Warning  - {warn}\n")

                    master_entry.remove_saveframe(sf_framecode)

                    truncated = True

            if len(empty_sf_framecodes) > 0:
                for sf_framecode in empty_sf_framecodes:
                    master_entry.remove_saveframe(sf_framecode)

                truncated = True

            if truncated:

                count = 0
                for idx, sf in enumerate(master_entry.get_saveframes_by_category(sf_category), start=1):
                    set_sf_tag(sf, 'ID', idx)
                    self.__reg.c2S.set_local_sf_id(sf, idx)
                    count += 1

                _sf_category = 'entry_information'

                if _sf_category in self.__reg.sf_category_list:
                    ent_sf = master_entry.get_saveframes_by_category(_sf_category)[0]

                    _lp_category = '_Data_set'

                    try:

                        _lp = ent_sf.get_loop(_lp_category)

                        type_col = _lp.tags.index('Type')
                        count_col = _lp.tags.index('Count')

                        for _idx, _row in enumerate(_lp):
                            if _row[type_col] == 'spectral_peak_list':
                                _lp.data[_idx][count_col] = count
                                break

                    except (KeyError, ValueError):
                        pass

            if len(update_sf_name) > 0:
                for k, v in update_sf_name.items():
                    try:
                        sf = master_entry.get_saveframe_by_name(k)
                        sf.name = v
                    except KeyError:
                        pass

        count = len(master_entry.get_saveframes_by_category(sf_category))
        if count > 0:
            input_source = self.__reg.report.input_sources[0]
            input_source_dic = input_source.get()

            if isinstance(input_source_dic['content_subtype'], dict):
                input_source_dic['content_subtype']['spectral_peak'] = count

        # update sample and sample_condition of spectral peak list saveframe

        if sf_category in self.__reg.sf_category_list:
            for sf in master_entry.get_saveframes_by_category(sf_category):

                if any(True for t in sf.tags if t[0] == 'Chem_shift_reference_ID'):
                    chem_shift_ref_id = get_first_sf_tag(sf, 'Chem_shift_reference_ID')

                    if chem_shift_ref_id in EMPTY_VALUE:
                        _sf_category = 'chem_shift_reference'

                        if _sf_category in self.__reg.sf_category_list:
                            _sf = master_entry.get_saveframes_by_category(_sf_category)[0]

                            set_sf_tag(sf, 'Chem_shift_reference_ID', get_first_sf_tag(_sf, 'ID'))

                exp_id = get_first_sf_tag(sf, 'Experiment_ID')

                if exp_id not in EMPTY_VALUE:
                    exp_id = int(exp_id)

                    _sf_category = 'experiment_list'

                    if _sf_category in self.__reg.sf_category_list:
                        _sf = master_entry.get_saveframes_by_category(_sf_category)[0]

                        _lp_category = '_Experiment'

                        try:

                            _lp = _sf.get_loop(_lp_category)

                            _tags = ['ID', 'Sample_ID', 'Sample_label', 'Sample_condition_list_ID', 'Sample_condition_list_label']

                            _dat = _lp.get_tag(_tags)

                            for _row in _dat:
                                if _row[0] in EMPTY_VALUE or int(_row[0]) != exp_id:
                                    continue
                                if _row[1] not in EMPTY_VALUE:
                                    set_sf_tag(sf, 'Sample_ID', _row[1])
                                if _row[2] not in EMPTY_VALUE:
                                    set_sf_tag(sf, 'Sample_label', _row[2])
                                if _row[3] not in EMPTY_VALUE:
                                    set_sf_tag(sf, 'Sample_condition_list_ID', _row[3])
                                if _row[4] not in EMPTY_VALUE:
                                    set_sf_tag(sf, 'Sample_condition_list_label', _row[4])
                                break

                            continue

                        except KeyError:
                            pass

                sample_id = get_first_sf_tag(sf, 'Sample_ID')

                if sample_id in EMPTY_VALUE:
                    _sf_category = 'sample'

                    _sf_list = master_entry.get_saveframes_by_category(_sf_category)

                    if len(_sf_list) == 1:
                        set_sf_tag(sf, 'Sample_ID', get_first_sf_tag(_sf_list[0], 'ID'))
                        if len(get_first_sf_tag(_sf_list[0], 'Sf_framecode')) > 0:
                            set_sf_tag(sf, 'Sample_label', f"${get_first_sf_tag(_sf_list[0], 'Sf_framecode')}")

                sample_condition_list_id = get_first_sf_tag(sf, 'Sample_condition_list_ID')

                if sample_condition_list_id in EMPTY_VALUE:
                    _sf_category = 'sample_conditions'

                    _sf_list = master_entry.get_saveframes_by_category(_sf_category)

                    if len(_sf_list) == 1:
                        set_sf_tag(sf, 'Sample_condition_list_ID', get_first_sf_tag(_sf_list[0], 'ID'))

        # section 10: chemical shift reference

        sf_category = 'chem_shift_reference'

        if has_cs_lp:
            if sf_category in self.__reg.sf_category_list:
                for sf in master_entry.get_saveframes_by_category(sf_category):
                    try:
                        cs_ref_id = int(get_first_sf_tag(sf, 'ID'))
                        isotope_numbers = isotope_nums[cs_ref_id]

                        set_sf_tag(sf, 'Proton_shifts_flag',
                                   'yes' if any(n in isotope_numbers for n in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS['H']) else 'no')
                        set_sf_tag(sf, 'Carbon_shifts_flag',
                                   'yes' if any(n in isotope_numbers for n in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS['C']) else 'no')
                        set_sf_tag(sf, 'Nitrogen_shifts_flag',
                                   'yes' if any(n in isotope_numbers for n in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS['N']) else 'no')
                        set_sf_tag(sf, 'Phosphorus_shifts_flag',
                                   'yes' if any(n in isotope_numbers for n in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS['P']) else 'no')
                        set_sf_tag(sf, 'Other_shifts_flag',
                                   'yes' if any(n not in WELL_KNOWN_ISOTOPE_NUMBERS for n in isotope_numbers) else 'no')

                        lp_category = '_Chem_shift_ref'

                        try:

                            lp = sf.get_loop(lp_category)

                            dat = lp.get_tag(['Atom_isotope_number', 'Indirect_shift_ratio'])

                            tags = lp.tags
                            ind_shift_ratio_col = tags.index('Indirect_shift_ratio')

                            _isotope_numbers = set()

                            for idx, row in enumerate(dat):
                                if isinstance(row[0], int):
                                    _isotope_numbers.add(row[0])
                                else:
                                    if row[0] in EMPTY_VALUE or not row[0].isdigit():
                                        continue
                                    _isotope_numbers.add(int(row[0]))
                                if row[1] not in EMPTY_VALUE:
                                    try:
                                        n = int(row[0]) if isinstance(row[0], str) else row[0]
                                        ratio = float(row[1])
                                        if ratio <= 0.0 or ratio >= 1.0:
                                            if n in (1, 19, 31):
                                                lp.data[idx][ind_shift_ratio_col] = '1.0'
                                            elif n == 2:
                                                lp.data[idx][ind_shift_ratio_col] = '0.153506088'
                                            elif n == 13:
                                                lp.data[idx][ind_shift_ratio_col] = '0.251449530'
                                            elif n == 15:
                                                lp.data[idx][ind_shift_ratio_col] = '0.101329118'
                                            else:
                                                lp.data[idx][ind_shift_ratio_col] = '?'
                                        elif ratio < 0.5:
                                            if n in (1, 19, 31):
                                                lp.data[idx][ind_shift_ratio_col] = '1.0'
                                    except ValueError:
                                        if n in (1, 19, 31):
                                            lp.data[idx][ind_shift_ratio_col] = '1.0'
                                        elif n == 2:
                                            lp.data[idx][ind_shift_ratio_col] = '0.153506088'
                                        elif n == 13:
                                            lp.data[idx][ind_shift_ratio_col] = '0.251449530'
                                        elif n == 15:
                                            lp.data[idx][ind_shift_ratio_col] = '0.101329118'
                                        else:
                                            lp.data[idx][ind_shift_ratio_col] = '?'

                            if isotope_numbers == _isotope_numbers:
                                continue

                            isotope_number_not_in_lp = isotope_numbers - _isotope_numbers

                            if len(isotope_number_not_in_lp) > 0:
                                atom_type_col = tags.index('Atom_type') if 'Atom_type' in tags else -1
                                atom_iso_num_col = tags.index('Atom_isotope_number') if 'Atom_isotope_number' in tags else -1
                                ref_type_col = tags.index('Ref_type') if 'Ref_type' in tags else -1
                                mol_com_name_col = tags.index('Mol_common_name') if 'Mol_common_name' in tags else -1
                                atom_grp_col = tags.index('Atom_group') if 'Atom_group' in tags else -1

                                for n in isotope_number_not_in_lp:
                                    if n in ALLOWED_ISOTOPE_NUMBERS:
                                        row = copy.copy(lp.data[0])

                                        if atom_type_col != -1:
                                            row[atom_type_col] = next(k for k, v in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS.items()
                                                                      if n in v)
                                        if atom_iso_num_col != -1:
                                            row[atom_iso_num_col] = n
                                        if ref_type_col != -1:
                                            row[ref_type_col] = 'direct' if n in (1, 19, 31) else 'indirect'
                                        if ind_shift_ratio_col != -1:
                                            if n in (1, 19, 31):
                                                row[ind_shift_ratio_col] = '1.0'
                                            elif n == 2:
                                                row[ind_shift_ratio_col] = '0.153506088'
                                            elif n == 13:
                                                row[ind_shift_ratio_col] = '0.251449530'
                                            elif n == 15:
                                                row[ind_shift_ratio_col] = '0.101329118'
                                            else:
                                                row[ind_shift_ratio_col] = '?'
                                        if n in (19, 31):
                                            if mol_com_name_col != -1:
                                                row[mol_com_name_col] = 'CCl3F' if n == 19 else '(MeO)3PO'
                                            if atom_grp_col != -1:
                                                row[atom_grp_col] = 'fluorine' if n == 19 else 'phosphorus'

                                        lp.add_data(row)

                            isotope_number_exe_in_lp = _isotope_numbers - isotope_numbers

                            if len(isotope_number_exe_in_lp) > 0:

                                del_row_idx = []

                                for n in isotope_number_exe_in_lp:

                                    try:
                                        del_row_idx.append(next(idx for idx, row in enumerate(dat)
                                                                if isinstance(row[0], str) and row[0] not in EMPTY_VALUE
                                                                and row[0].isdigit() and int(row[0]) == n))
                                    except StopIteration:
                                        continue

                                if len(del_row_idx) > 0:
                                    for idx in reversed(del_row_idx):
                                        try:
                                            del lp.data[idx]
                                        except IndexError:
                                            pass

                        except KeyError:

                            items = ['Atom_type',
                                     'Atom_isotope_number',
                                     'Mol_common_name',
                                     'Atom_group',
                                     'Concentration_val',
                                     'Concentration_units',
                                     'Solvent',
                                     'Rank',
                                     'Chem_shift_units',
                                     'Chem_shift_val',
                                     'Ref_method',
                                     'Ref_type',
                                     'Indirect_shift_ratio',
                                     'External_ref_loc',
                                     'External_ref_sample_geometry',
                                     'External_ref_axis',
                                     'Indirect_shift_ratio_cit_ID',
                                     'Indirect_shift_ratio_cit_label',
                                     'Ref_correction_type',
                                     'Correction_val',
                                     'Correction_val_cit_ID',
                                     'Correction_val_cit_label',
                                     'Entry_ID',
                                     'Chem_shift_reference_ID']

                            lp = pynmrstar.Loop.from_scratch(lp_category)

                            tags = [f'{lp_category}.{item}' for item in items]

                            lp.add_tag(tags)

                            for n in isotope_numbers:
                                if n in ALLOWED_ISOTOPE_NUMBERS:
                                    row = ['.'] * len(items)
                                    row[0] = next(k for k, v in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS.items() if n in v)
                                    row[1] = n
                                    row[8], row[9], row[10] = 'ppm', '0.000', 'internal'
                                    if n in (1, 2, 13, 15, 19, 31):
                                        if n == 1:
                                            row[2], row[3], row[11], row[12] =\
                                                default_internal_reference, 'methyl protons', 'direct', '1.0'
                                        elif n == 2:
                                            row[2], row[3], row[11], row[12] =\
                                                default_internal_reference, 'methyl protons', 'indirect', '0.153506088'
                                        elif n == 13:
                                            row[2], row[3], row[11], row[12] =\
                                                default_internal_reference, 'methyl protons', 'indirect', '0.251449530'
                                        elif n == 15:
                                            row[2], row[3], row[11], row[12] =\
                                                default_internal_reference, 'methyl protons', 'indirect', '0.101329118'
                                        elif n == 19:
                                            row[2], row[3], row[11], row[12] = 'CCl3F', 'fluorine', 'direct', '1.0'
                                        else:
                                            row[2], row[3], row[11], row[12] = '(MeO)3PO', 'phosphorus', 'direct', '1.0'
                                    else:
                                        row[2], row[3], row[11], row[12] = '?', '?', '?', '?'

                                    row[-2], row[-1] = self.__reg.entry_id, cs_ref_id

                                    lp.add_data(row)

                            sf.add_loop(lp)

                    except ValueError:
                        continue

            elif len(cs_ref_sf_framecode) > 0:
                for cs_ref_id, sf_framecode in cs_ref_sf_framecode.items():
                    isotope_numbers = isotope_nums[cs_ref_id]

                    _sf = pynmrstar.Saveframe.from_scratch(sf_framecode, '_Chem_shift_reference')
                    _sf.add_tag('Sf_category', 'chem_shift_reference')
                    _sf.add_tag('Sf_framecode', sf_framecode)
                    _sf.add_tag('Entry_ID', self.__reg.entry_id)
                    _sf.add_tag('ID', cs_ref_id)
                    _sf.add_tag('Name', '.')
                    _sf.add_tag('Proton_shifts_flag',
                                'yes' if any(n in isotope_numbers for n in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS['H']) else 'no')
                    _sf.add_tag('Carbon_shifts_flag',
                                'yes' if any(n in isotope_numbers for n in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS['C']) else 'no')
                    _sf.add_tag('Nitrogen_shifts_flag',
                                'yes' if any(n in isotope_numbers for n in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS['N']) else 'no')
                    _sf.add_tag('Phosphorus_shifts_flag',
                                'yes' if any(n in isotope_numbers for n in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS['P']) else 'no')
                    _sf.add_tag('Other_shifts_flag',
                                'yes' if any(n not in WELL_KNOWN_ISOTOPE_NUMBERS for n in isotope_numbers) else 'no')
                    _sf.add_tag('Details', '.')

                    lp_category = '_Chem_shift_ref'

                    items = ['Atom_type',
                             'Atom_isotope_number',
                             'Mol_common_name',
                             'Atom_group',
                             'Concentration_val',
                             'Concentration_units',
                             'Solvent',
                             'Rank',
                             'Chem_shift_units',
                             'Chem_shift_val',
                             'Ref_method',
                             'Ref_type',
                             'Indirect_shift_ratio',
                             'External_ref_loc',
                             'External_ref_sample_geometry',
                             'External_ref_axis',
                             'Indirect_shift_ratio_cit_ID',
                             'Indirect_shift_ratio_cit_label',
                             'Ref_correction_type',
                             'Correction_val',
                             'Correction_val_cit_ID',
                             'Correction_val_cit_label',
                             'Entry_ID',
                             'Chem_shift_reference_ID']

                    lp = pynmrstar.Loop.from_scratch(lp_category)

                    tags = [f'{lp_category}.{item}' for item in items]

                    lp.add_tag(tags)

                    for n in isotope_numbers:
                        if n in ALLOWED_ISOTOPE_NUMBERS:
                            row = ['.'] * len(items)
                            row[0] = next(k for k, v in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS.items() if n in v)
                            row[1] = n
                            row[8], row[9], row[10] = 'ppm', '0.000', 'internal'
                            if n in (1, 2, 13, 15, 19, 31):
                                if n == 1:
                                    row[2], row[3], row[11], row[12] =\
                                        default_internal_reference, 'methyl protons', 'direct', '1.0'
                                elif n == 2:
                                    row[2], row[3], row[11], row[12] =\
                                        default_internal_reference, 'methyl protons', 'indirect', '0.153506088'
                                elif n == 13:
                                    row[2], row[3], row[11], row[12] =\
                                        default_internal_reference, 'methyl protons', 'indirect', '0.251449530'
                                elif n == 15:
                                    row[2], row[3], row[11], row[12] =\
                                        default_internal_reference, 'methyl protons', 'indirect', '0.101329118'
                                elif n == 19:
                                    row[2], row[3], row[11], row[12] = 'CCl3F', 'fluorine', 'direct', '1.0'
                                else:
                                    row[2], row[3], row[11], row[12] = '(MeO)3PO', 'phosphorus', 'direct', '1.0'
                            else:
                                row[2], row[3], row[11], row[12] = '?', '?', '?', '?'

                            row[-2], row[-1] = self.__reg.entry_id, cs_ref_id

                            lp.add_data(row)

                    _sf.add_loop(lp)

            sf_wo_lp = []

            for sf in master_entry.get_saveframes_by_category(sf_category):
                try:
                    lp = sf.get_loop('_Chem_shift_ref')
                    if len(lp) == 0:
                        sf_wo_lp.append(sf.name)
                except KeyError:
                    sf_wo_lp.append(sf.name)

            if len(sf_wo_lp) > 0:
                for sf_framecode in sf_wo_lp:
                    master_entry.remove_saveframe(sf_framecode)

        # check order of experiment id and synchronize them in the entry

        sf_category = 'experiment_list'

        if sf_category in self.__reg.sf_category_list:
            sf = master_entry.get_saveframes_by_category(sf_category)[0]

            lp_category = '_Experiment'

            try:

                lp = sf.get_loop(lp_category)

                dat = lp.get_tag(['ID'])

                exp_id_mapping = {}

                for idx, row in enumerate(dat, start=1):
                    if row not in EMPTY_VALUE and int(row) != idx:
                        exp_id_mapping[int(row)] = idx

                if len(exp_id_mapping) > 0:
                    id_col = lp.tags.index('ID')

                    for idx, row in enumerate(lp):
                        if row[id_col] not in EMPTY_VALUE and int(row[id_col]) in exp_id_mapping:
                            lp.data[idx][id_col] = exp_id_mapping[int(row[id_col])]

                    _sf_category = 'assigned_chemical_shifts'

                    if _sf_category in self.__reg.sf_category_list:
                        for _sf in master_entry.get_saveframes_by_category(_sf_category):

                            _lp_category = '_Chem_shift_experiment'

                            try:

                                _lp = _sf.get_loop(_lp_category)

                                exp_id_col = _lp.tags.index('Experiment_ID')

                                for _idx, _row in enumerate(_lp):
                                    if _row[exp_id_col] not in EMPTY_VALUE and int(_row[exp_id_col]) in exp_id_mapping:
                                        _lp.data[_idx][exp_id_col] = exp_id_mapping[int(_row[exp_id_col])]

                            except KeyError:
                                continue

                    _sf_category = 'spectral_peak_list'

                    if _sf_category in self.__reg.sf_category_list:
                        for _sf in master_entry.get_saveframes_by_category(_sf_category):
                            exp_id = get_first_sf_tag(_sf, 'Experiment_ID')

                            if exp_id not in EMPTY_VALUE and int(exp_id) in exp_id_mapping:
                                set_sf_tag(_sf, 'Experiment_ID', exp_id_mapping[int(exp_id)])

                tags = ['ID', 'Name']

                dat = lp.get_tag(tags)

                ref_row = lp.data[0]

                _sf_category = 'assigned_chemical_shifts'

                if _sf_category in self.__reg.sf_category_list:
                    for _sf in master_entry.get_saveframes_by_category(_sf_category):

                        _lp_category = '_Chem_shift_experiment'

                        try:

                            _lp = _sf.get_loop(_lp_category)

                            _tags = ['Experiment_ID', 'Experiment_name', 'Sample_ID', 'Sample_label', 'Sample_state']

                            _dat = _lp.get_tag(_tags)

                            for _idx, _row in enumerate(_dat):

                                if _row[1] in EMPTY_VALUE or any(True for row in dat if row[1] == _row[1]):
                                    continue

                                new_exp_id = str(max(int(row[0]) for row in dat) + 1)

                                tmp = copy.copy(ref_row)
                                tmp[lp.tags.index('ID')] = new_exp_id
                                tmp[lp.tags.index('Name')] = _row[1]
                                tmp[lp.tags.index('Sample_ID')] = _row[2]
                                tmp[lp.tags.index('Sample_label')] = _row[3]
                                tmp[lp.tags.index('Sample_state')] = _row[4]
                                tmp[lp.tags.index('Sample_condition_list_ID')] =\
                                    get_first_sf_tag(_sf, 'Sample_condition_list_ID',
                                                     ref_row[lp.tags.index('Sample_condition_list_ID')])
                                tmp[lp.tags.index('Sample_condition_list_label')] = '$'\
                                    + get_first_sf_tag(_sf, 'Sample_condition_list_label',
                                                       ref_row[lp.tags.index('Sample_condition_list_label')])

                                lp.add_data(tmp)
                                dat = lp.get_tag(tags)

                                _lp.data[_idx][_lp.tags.index('Experiment_ID')] = new_exp_id

                        except KeyError:
                            continue

                _sf_category = 'spectral_peak_list'

                if _sf_category in self.__reg.sf_category_list:
                    for _sf in master_entry.get_saveframes_by_category(_sf_category):
                        exp_id = get_first_sf_tag(_sf, 'Experiment_ID')
                        exp_name = get_first_sf_tag(_sf, 'Experiment_name')

                        if exp_name in EMPTY_VALUE or any(True for row in dat if row[1] == exp_name):
                            continue

                        new_exp_id = str(max(int(row[0]) for row in dat) + 1)

                        tmp = copy.copy(ref_row)
                        tmp[lp.tags.index('ID')] = new_exp_id
                        tmp[lp.tags.index('Name')] = exp_name
                        tmp[lp.tags.index('Sample_ID')] = get_first_sf_tag(_sf, 'Sample_ID')
                        tmp[lp.tags.index('Sample_label')] = get_first_sf_tag(_sf, 'Sample_label')
                        tmp[lp.tags.index('Sample_condition_list_ID')] =\
                            get_first_sf_tag(_sf, 'Sample_condition_list_ID',
                                             ref_row[lp.tags.index('Sample_condition_list_ID')])
                        tmp[lp.tags.index('Sample_condition_list_label')] = '$'\
                            + get_first_sf_tag(_sf, 'Sample_condition_list_label',
                                               ref_row[lp.tags.index('Sample_condition_list_label')])

                        lp.add_data(tmp)
                        dat = lp.get_tag(tags)

                        set_sf_tag(_sf, 'Experiment_ID', new_exp_id)

                exp_tags = ['ID', 'Name', 'NMR_spectrometer_ID']

                if set(exp_tags) & set(lp.tags) == set(exp_tags):
                    exp_list = lp.get_tag(exp_tags)

                _sf_category = 'spectral_peak_list'

                if _sf_category in self.__reg.sf_category_list:

                    for _sf in master_entry.get_saveframes_by_category(_sf_category):
                        exp_id = get_first_sf_tag(_sf, 'Experiment_ID')

                        if exp_id not in EMPTY_VALUE:

                            if isinstance(exp_id, int):
                                exp_id = str(int)

                            exp_row = next((row for row in exp_list
                                            if (row[0] == exp_id or (isinstance(row[0], int) and str(row[0]) == exp_id))), None)

                            if exp_row is not None:
                                set_sf_tag(_sf, 'Experiment_name', exp_row[1])

                                spectrometer_id = exp_row[2]

                                if spectrometer_id in EMPTY_VALUE:
                                    continue

                                if isinstance(spectrometer_id, str):

                                    if not spectrometer_id.isdigit():
                                        continue

                                    spectrometer_id = int(spectrometer_id)

                                if spectrometer_id not in spectrometer_dict:
                                    continue

                                field_strength = spectrometer_dict[spectrometer_id]['field_strength']

                                try:

                                    if isinstance(field_strength, str):
                                        field_strength = float(field_strength)

                                except ValueError:
                                    continue

                                cs_ref_ratio_map = {}

                                _lp_category = 'Chem_shift_ref'

                                try:

                                    for _lp in master_entry.get_loops_by_category(_lp_category):
                                        isotope_num_col = _lp.tags.index('Atom_isotope_number')
                                        ratio_col = _lp.tags.index('Indirect_shift_ratio')

                                        for _row in _lp:
                                            isotope_num = _row[isotope_num_col]
                                            ratio = float(_row[ratio_col])

                                            if isotope_num not in cs_ref_ratio_map:
                                                cs_ref_ratio_map[isotope_num] = ratio

                                except (KeyError, ValueError):
                                    continue

                                if len(cs_ref_ratio_map) == 0:
                                    continue

                                _lp_category = 'Spectral_dim'

                                try:

                                    _lp = _sf.get_loop(_lp_category)

                                    isotope_num_col = _lp.tags.index('Atom_isotope_number')
                                    spec_freq_col = _lp.tags.index('Spectrometer_frequency')

                                    for idx, _row in enumerate(_lp):
                                        isotope_num = _row[isotope_num_col]
                                        spec_freq = _row[spec_freq_col]

                                        if spec_freq in EMPTY_VALUE and isotope_num in cs_ref_ratio_map:
                                            _lp.data[idx][spec_freq_col] =\
                                                roundString(f'{field_strength * cs_ref_ratio_map[isotope_num]}',
                                                            field_strength_max_didits)

                                except (KeyError, ValueError):
                                    continue

                            else:
                                num_of_dim = get_first_sf_tag(_sf, 'Number_of_spectral_dimensions')
                                if num_of_dim in EMPTY_VALUE:
                                    set_sf_tag(_sf, 'Experiment_ID', None)
                                    continue
                                if isinstance(num_of_dim, str):
                                    if num_of_dim.isdigit():
                                        num_of_dim = int(num_of_dim)
                                    else:
                                        set_sf_tag(_sf, 'Experiment_ID', None)
                                        continue
                                exp_class = get_first_sf_tag(_sf, 'Experiment_class')
                                if exp_class in EMPTY_VALUE:
                                    set_sf_tag(_sf, 'Experiment_ID', None)
                                    continue
                                if 'through-space' in exp_class:
                                    exp_rows = [row for row in exp_list if f'{num_of_dim}D' in row[1] and 'NOE' in row[1]]
                                    if len(exp_rows) != 1:
                                        set_sf_tag(_sf, 'Experiment_ID', None)
                                        continue
                                    exp_row = exp_rows[0]
                                    set_sf_tag(_sf, 'Experiment_ID', exp_row[0])
                                    set_sf_tag(_sf, 'Experiment_name', exp_row[1])
                                if 'relayed' in exp_class:
                                    exp_rows = [row for row in exp_list if f'{num_of_dim}D' in row[1] and 'TOCSY' in row[1]]
                                    if len(exp_rows) != 1:
                                        set_sf_tag(_sf, 'Experiment_ID', None)
                                        continue
                                    exp_row = exp_rows[0]
                                    set_sf_tag(_sf, 'Experiment_ID', exp_row[0])
                                    set_sf_tag(_sf, 'Experiment_name', exp_row[1])
                                if 'jcoupling' in exp_class:
                                    exp_rows = [row for row in exp_list if f'{num_of_dim}D' in row[1] and 'COSY' in row[1]]
                                    if len(exp_rows) != 1:
                                        set_sf_tag(_sf, 'Experiment_ID', None)
                                        continue
                                    exp_row = exp_rows[0]
                                    set_sf_tag(_sf, 'Experiment_ID', exp_row[0])
                                    set_sf_tag(_sf, 'Experiment_name', exp_row[1])

                                if exp_row is None:
                                    continue

                                spectrometer_id = exp_row[2]

                                if spectrometer_id in EMPTY_VALUE:
                                    continue

                                if isinstance(spectrometer_id, str):

                                    if not spectrometer_id.isdigit():
                                        continue

                                    spectrometer_id = int(spectrometer_id)

                                if spectrometer_id not in spectrometer_dict:
                                    continue

                                field_strength = spectrometer_dict[spectrometer_id]['field_strength']

                                try:

                                    if isinstance(field_strength, str):
                                        field_strength = float(field_strength)

                                except ValueError:
                                    continue

                                cs_ref_ratio_map = {}

                                _lp_category = 'Chem_shift_ref'

                                try:

                                    for _lp in master_entry.get_loops_by_category(_lp_category):
                                        isotope_num_col = _lp.tags.index('Atom_isotope_number')
                                        ratio_col = _lp.tags.index('Indirect_shift_ratio')

                                        for _row in _lp:
                                            isotope_num = _row[isotope_num_col]
                                            ratio = float(_row[ratio_col])

                                            if isotope_num not in cs_ref_ratio_map:
                                                cs_ref_ratio_map[isotope_num] = ratio

                                except (KeyError, ValueError):
                                    continue

                                if len(cs_ref_ratio_map) == 0:
                                    continue

                                _lp_category = 'Spectral_dim'

                                try:

                                    _lp = _sf.get_loop(_lp_category)

                                    isotope_num_col = _lp.tags.index('Atom_isotope_number')
                                    spec_freq_col = _lp.tags.index('Spectrometer_frequency')

                                    for idx, _row in enumerate(_lp):
                                        isotope_num = _row[isotope_num_col]
                                        spec_freq = _row[spec_freq_col]

                                        if spec_freq in EMPTY_VALUE and isotope_num in cs_ref_ratio_map:
                                            _lp.data[idx][spec_freq_col] =\
                                                roundString(f'{field_strength * cs_ref_ratio_map[isotope_num]}',
                                                            field_strength_max_didits)

                                except (KeyError, ValueError):
                                    continue

            except KeyError:
                pass

        # section 5: polymer residues and ligands

        if int_sf is not None:
            if 'assigned_chemical_shifts' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Assigned_chem_shifts', 'yes')
            if 'general_distance_constraints' in self.__reg.sf_category_list\
               or 'torsion_angle_constraints' in self.__reg.sf_category_list\
               or 'RDC_constraints' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Constraints', 'yes')
            if 'coupling_constant' in self.__reg.sf_category_list\
               or 'J_three_bond_constraints' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Coupling_constants', 'yes')
            if 'chem_shift_isotope_effect' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Chem_shift_isotope_effect', 'yes')
            if 'chem_shift_perturbation' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Chem_shift_perturbation', 'yes')
            if 'auto_relaxation' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Auto_relaxation', 'yes')
            if 'tensor' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Tensor', 'yes')
            if 'general_distance_constraints' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Interatomic_distance', 'yes')
            if 'chem_shift_anisotropy' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Chem_shift_anisotropy', 'yes')
            if 'heteronucl_NOEs' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Heteronucl_NOEs', 'yes')
            if 'heteronucl_T1_relaxation' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Heteronucl_T1_relaxation', 'yes')
            if 'heteronucl_T2_relaxation' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Heteronucl_T2_relaxation', 'yes')
            if 'heteronucl_T1rho_relaxation' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Heteronucl_T1rho_relaxation', 'yes')
            if 'order_parameters' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Order_parameters', 'yes')
            if 'RDC_constraints' in self.__reg.sf_category_list\
               or 'RDCs' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Residual_dipolar_couplings', 'yes')
            if 'H_exch_rates' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'H_exchange_rate', 'yes')
            if 'H_exch_protection_factors' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'H_exchange_protection_factors', 'yes')
            if 'chemical_rates' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Chem_rate_constants', 'yes')
            if 'spectral_peak_list' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Spectral_peak_lists', 'yes')
            if 'dipolar_couplings' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Dipole_dipole_couplings', 'yes')
            if 'homonucl_NOEs' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Homonucl_NOEs', 'yes')
            if 'dipole_dipole_relaxation' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Dipole_dipole_relaxation', 'yes')
            if 'dipole_dipole_cross_correlations' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'DD_cross_correlation', 'yes')
            if 'dipole_CSA_cross_correlations' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Dipole_CSA_cross_correlation', 'yes')
            if 'binding_constants' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Binding_constants', 'yes')
            if 'pKa_value_data_set' in self.__reg.sf_category_list\
               or 'pH_NMR_param_list' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'PKa_value_data_set', 'yes')
            if 'D_H_fractionation_factors' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'D_H_fractionation_factors', 'yes')
            if 'theoretical_chem_shifts' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Theoretical_chem_shifts', 'yes')
            if 'theoretical_coupling_constants' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Theoretical_coupling_constants', 'yes')
            if 'theoretical_heteronucl_NOEs' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Theoretical_heteronucl_NOEs', 'yes')
            if 'theoretical_heteronucl_T1_relaxation' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Theoretical_T1_relaxation', 'yes')
            if 'theoretical_heteronucl_T2_relaxation' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Theoretical_T2_relaxation', 'yes')
            if 'theoretical_auto_relaxation' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Theoretical_auto_relaxation', 'yes')
            if 'theoretical_dipole_dipole_cross_correlations' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Theoretical_DD_cross_correlation', 'yes')
            if 'spectral_density_values' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Spectral_density_values', 'yes')
            if 'other_data_types' in self.__reg.sf_category_list:
                set_sf_tag(int_sf, 'Other_kind_of_data', 'yes')

        if has_non_polymer or has_nstd_monomer:

            if has_non_polymer:

                if int_sf is not None:
                    set_sf_tag(int_sf, 'Ligands', 'yes')

                # wo coordinates (bmrbdep)
                if self.__reg.report.getInputSourceIdOfCoord() < 0:

                    for sf in master_entry.get_saveframes_by_category(sf_category):
                        try:
                            entity_id = int(get_first_sf_tag(sf, 'ID'))

                            if get_first_sf_tag(sf, 'Type') == 'non-polymer':
                                comp_id = get_first_sf_tag(sf, 'BMRB_code')
                                if comp_id in EMPTY_VALUE:
                                    comp_id = get_first_sf_tag(sf, 'Nonpolymer_comp_ID')
                                if comp_id in EMPTY_VALUE:
                                    continue

                        except ValueError:
                            continue

            if has_nstd_monomer:

                if int_sf is not None:
                    set_sf_tag(int_sf, 'Non_standard_residues', 'yes')

        # section 4: biological polymers and ligands

        total_cys_in_assembly = 0

        if 'assembly' in self.__reg.sf_category_list:

            lp_category = '_Entity_assembly'

            try:

                lp = asm_sf.get_loop(lp_category)

                dat = lp.get_tag(['Entity_ID'])

                for row in dat:
                    if row not in EMPTY_VALUE:
                        entity_id = int(row)
                        if entity_id in entity_dict:
                            entity = entity_dict[entity_id]
                            sample_type = entity['sample_type']
                            if sample_type in ('protein', 'peptide'):
                                total_cys_in_assembly += entity['total_cys']

                annotate_assembly_experimental_data(lp, ent_asym_id_with_exptl_data)

            except KeyError:
                pass

        for entity_id, entity in entity_dict.items():
            sample_type = entity['sample_type']

            sf_category = 'entity'
            sf = next((sf for sf in master_entry.get_saveframes_by_category(sf_category)
                       if str(get_first_sf_tag(sf, 'ID')) == str(entity_id)), None)
            if sf is None:
                continue

            sf_framecode = get_first_sf_tag(sf, 'Sf_framecode')
            thiol_state = get_first_sf_tag(sf, 'Thiol_state')
            if sample_type in ('protein', 'peptide', 'DNA', 'RNA', 'DNA/RNA hybrid'):

                lp_category = '_Entity_comp_index'
                cys_total = 0

                try:

                    one_letter_code = ''

                    lp = sf.get_loop(lp_category)

                    dat = lp.get_tag(['Comp_ID'])

                    for idx, row in enumerate(dat, start=1):
                        one_letter_code += getOneLetterCodeCan(row)
                        if idx % 20 == 0:
                            one_letter_code += '\n'
                        if row not in EMPTY_VALUE:
                            if row in ('CYS', 'DCY'):
                                cys_total += 1

                    if len(one_letter_code) > 0 and len(dat) % 20 != 0:
                        one_letter_code += '\n'

                except KeyError:
                    pass

                if cys_total == 0:
                    if thiol_state in EMPTY_VALUE or thiol_state not in ALLOWED_THIOL_STATES:
                        set_sf_tag(sf, 'Thiol_state', 'not present')

                elif has_reducing_agent:
                    if thiol_state in EMPTY_VALUE or thiol_state not in ALLOWED_THIOL_STATES:
                        set_sf_tag(sf, 'Thiol_state', 'all free')

                elif total_cys_in_assembly == 1:
                    if thiol_state in EMPTY_VALUE or thiol_state not in ALLOWED_THIOL_STATES:
                        set_sf_tag(sf, 'Thiol_state', 'all free')

            else:
                thiol_state = get_first_sf_tag(sf, 'Thiol_state')
                if thiol_state in EMPTY_VALUE or thiol_state not in ALLOWED_THIOL_STATES:
                    set_sf_tag(sf, 'Thiol_state', 'not present')

        if total_cys_in_assembly > 1:
            target_entity_sf_framecodes, detailed_target_entity_sf_framecodes = [], []

            for entity_id, entity in entity_dict.items():
                sample_type = entity['sample_type']

                sf_category = 'entity'
                sf = next((sf for sf in master_entry.get_saveframes_by_category(sf_category)
                           if str(get_first_sf_tag(sf, 'ID')) == str(entity_id)), None)
                if sf is None:
                    continue

                sf_framecode = get_first_sf_tag(sf, 'Sf_framecode')
                thiol_state = get_first_sf_tag(sf, 'Thiol_state')
                if thiol_state in EMPTY_VALUE\
                   or thiol_state in ('not available', 'not present', 'not reported', 'unknown', 'all free'):
                    continue

                if sample_type in ('protein', 'peptide') and entity['total_cys'] > 0:
                    target_entity_sf_framecodes.append(f'save_{sf_framecode}')
                    detailed_target_entity_sf_framecodes.append(
                        f'_Entity.Thiol_state is {thiol_state} of save_{sf_framecode} frame')

            if len(target_entity_sf_framecodes) > 0:

                lp_category = '_Bond'

                try:

                    lp = asm_sf.get_loop(lp_category)

                    has_bond = len(lp) > 0

                except KeyError:
                    has_bond = False

                if has_bond:
                    _lp = pynmrstar.Loop.from_scratch(lp_category)

                    tags = ['ID', 'Type', 'Value_order',
                            'Entity_assembly_name_1', 'Comp_ID_1', 'Seq_ID_1', 'Atom_ID_1',
                            'Entity_assembly_name_2', 'Comp_ID_2', 'Seq_ID_2', 'Atom_ID_2']

                    _lp.add_tag(tags)

                    dat = lp.get_tag(tags)

                    for row in dat:
                        _lp.add_data(row)

                else:
                    _lp = pynmrstar.Loop.from_scratch(lp_category)

                    tags = ['ID', 'Type',
                            'Comp_ID_1', 'Seq_ID_1', 'Atom_ID_1',
                            'Comp_ID_2', 'Seq_ID_2', 'Atom_ID_2']

                    _lp.add_tag(tags)

                    _lp.add_data(['?'] * len(tags))

        for sf in master_entry.frame_list:
            sf_tag_prefix = sf.tag_prefix
            label_tags = [sf_tag.split('.')[1] for sf_tag in DEFAULT_SF_LABEL_TAGS
                          if sf_tag.startswith(f'{sf_tag_prefix}.') and sf_tag.endswith('_label')]
            if len(label_tags) == 0:
                continue
            for label_tag in label_tags:
                label = get_first_sf_tag(sf, label_tag)
                if label in EMPTY_VALUE:
                    parent_sf_tag_prefix = f'_{label_tag[:-6]}'
                    parent_list_id = get_first_sf_tag(sf, f'{parent_sf_tag_prefix[1:]}_ID')
                    if parent_list_id in EMPTY_VALUE:
                        continue
                    try:
                        parent_sf = master_entry.get_saveframes_by_tag_and_value(f'{parent_sf_tag_prefix}.ID', parent_list_id)[0]
                        parent_sf_framecode = get_first_sf_tag(parent_sf, 'Sf_framecode')
                        if len(parent_sf_framecode) > 0:
                            set_sf_tag(sf, label_tag, f'${parent_sf_framecode}')
                    except IndexError:
                        try:
                            parent_sf =\
                                master_entry.get_saveframes_by_tag_and_value(f'{parent_sf_tag_prefix}.ID',
                                                                             int(parent_list_id) if isinstance(parent_list_id, str)
                                                                             else str(parent_list_id))[0]
                            parent_sf_framecode = get_first_sf_tag(parent_sf, 'Sf_framecode')
                            if len(parent_sf_framecode) > 0:
                                set_sf_tag(sf, label_tag, f'${parent_sf_framecode}')
                        except IndexError:
                            set_sf_tag(sf, label_tag, '?')

        # cleanup
        if not self.__reg.annotation_mode:
            self.__reg.c2S.cleanup_str(master_entry)
            self.__reg.c2S.set_entry_id(master_entry, self.__reg.entry_id)
            self.__reg.c2S.normalize_str(master_entry)

        return True
