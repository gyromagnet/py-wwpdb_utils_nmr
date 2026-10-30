##
# File: PeakReaderRegressionTests.py
# Date: 29-Oct-2026
##
"""Regressions for independent peak readers and coordinates without experiment metadata."""
import io
import json
import sys
import unittest
from collections import defaultdict
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import MagicMock, Mock, patch

# Configure the test CCD before importing any scientific modules.
if __package__ is None or __package__ == "":
    from os import path

    sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
    from commonsetup import TESTOUTPUT  # noqa: F401 pylint: disable=import-error,unused-import
else:
    from .commonsetup import TESTOUTPUT  # noqa: F401 pylint: disable=relative-beyond-top-level

import pynmrstar

# PairwiseAlign is exported by the compiled alignment extension, invisible to static inference.
from wwpdb.utils.align.alignlib import PairwiseAlign  # pylint: disable=no-name-in-module
from wwpdb.utils.nmr.CifToNmrStar import CifToNmrStar
from wwpdb.utils.nmr.NmrDpUtility import NmrDpUtility
from wwpdb.utils.nmr.NmrVrptUtility import NmrVrptUtility
from wwpdb.utils.nmr.ann.BmrbAnnTasks import annotate_assembly_experimental_data
from wwpdb.utils.nmr.nef.NefTranslator import NefTranslator
from wwpdb.utils.nmr.pk.BasePKParserListener import BasePKParserListener
from wwpdb.utils.nmr.pk.XeasyPKReader import XeasyPKReader


class PeakReaderRegressionTests(unittest.TestCase):

    def test_normalization_stably_sorts_unknown_tags_without_mutating_schema(self):
        converter = CifToNmrStar(log=io.StringIO())
        original_order = list(converter.pynmrstar_schema.schema_order)
        results = []
        for reverse in (False, True):
            entry = pynmrstar.Entry.from_scratch('test')
            frame = pynmrstar.Saveframe.from_scratch('entry', '_Entry')
            frame.add_tag('Sf_category', 'entry_information')
            tags = [('User_extension_z', 'last'), ('User_extension_a', 'first')]
            for tag, value in reversed(tags) if reverse else tags:
                frame.add_tag(tag, value)
            loop = pynmrstar.Loop.from_scratch('_User_extension')
            tags = [('z', 'last'), ('a', 'first')]
            ordered = list(reversed(tags)) if reverse else tags
            loop.add_tag([tag for tag, _ in ordered])
            loop.add_data([value for _, value in ordered])
            frame.add_loop(loop)
            entry.add_saveframe(frame)
            converter.normalize_str(entry)
            self.assertEqual([tag[0] for tag in frame.tags][-2:], ['User_extension_a', 'User_extension_z'])
            self.assertEqual(loop.tags, ['a', 'z'])
            self.assertEqual(loop.data, [['first', 'last']])
            results.append(str(entry))
        self.assertEqual(results[0], results[1])
        self.assertEqual(converter.pynmrstar_schema.schema_order, original_order)

    def test_coordinate_merge_retains_grouped_author_atom_labels(self):
        coordinates = Path(__file__).resolve().parents[1] / 'tests-nmr/mock-data/2l9r.cif'
        entry = pynmrstar.Entry.from_file(str(coordinates.with_suffix('.str')))
        for frame in list(entry.frame_list):
            if frame.category not in ('entry_information', 'assembly', 'assigned_chemical_shifts'):
                entry.remove_saveframe(frame.name)
        # Use the existing matched structure and metadata. Remove optional original
        # naming so this exercises chain reconstruction from canonical identities.
        loop = entry.get_loops_by_category('Atom_chem_shift')[0]
        for tag in list(loop.tags):
            if tag.startswith('Original_PDB_'):
                loop.remove_tag(tag)
        atom_col = loop.tags.index('Atom_ID')
        seq_col = loop.tags.index('Comp_index_ID')
        auth_col = loop.tags.index('Auth_atom_ID')
        chain_col = loop.tags.index('Auth_asym_ID')
        for row in loop.data:
            row[chain_col] = '.'
            if row[seq_col] == '14' and row[atom_col] in ('HG21', 'HG22', 'HG23'):
                row[auth_col] = 'HG21'
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / 'shifts.str'
            entry.write_to_file(str(source))
            utility = NmrDpUtility()
            utility.setWorkspace(str(root / 'work'), str(root / 'cache'))
            utility.addInput(name='chem_shift_file_path_list',
                             value=[{'file_name': str(source), 'file_type': 'nm-shi'}], type='file_dict_list')
            utility.addInput(name='coordinate_file_path', value=str(coordinates), type='file')
            utility.addInput(name='remediation', value=True, type='param')
            utility.addInput(name='conversion_server', value=True, type='param')
            utility.addInput(name='nonblk_anomalous_cs', value=True, type='param')
            utility.setLog(str(root / 'report.json'))
            utility.setDestination(str(root / 'merged.str'))
            self.assertTrue(utility.op('nmr-cs-mr-merge'), (root / 'report.json').read_text(encoding='utf-8'))
            result = pynmrstar.Entry.from_file(str(root / 'merged.str'))
            rows = result.get_loops_by_category('Atom_chem_shift')[0].get_tag(
                ['Comp_index_ID', 'Atom_ID', 'Auth_atom_ID', 'Original_PDB_atom_name'])
            matched = [row for row in rows if row[0] == '14' and row[1] in ('HG21', 'HG22', 'HG23')]
            self.assertEqual(len(matched), 3)
            self.assertTrue(all(row[2:] == ['HG21', 'HG21'] for row in matched), matched)

    def test_shift_only_deposit_does_not_invent_metadata_references(self):
        source = Path(__file__).resolve().parent / 'data/explicit-shifts.str'
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            utility = NmrDpUtility()
            utility.setWorkspace(str(root / 'work'), str(root / 'cache'))
            utility.setSource(str(source))
            utility.setLog(str(root / 'report.json'))
            utility.setDestination(str(root / 'result.str'))
            utility.addInput(name='bmrb_only', value=True, type='param')
            utility.addOutput(name='entry_id', value='SHIFT_CONTROL', type='param')
            utility.addOutput(name='insert_entry_id_to_loops', value=True, type='param')
            self.assertTrue(utility.op('nmr-str2str-deposit'), (root / 'report.json').read_text(encoding='utf-8'))
            result = pynmrstar.Entry.from_file(str(root / 'result.str'))
            shifts = result.get_saveframes_by_category('assigned_chemical_shifts')[0]
            self.assertEqual(len(result.get_loops_by_category('Atom_chem_shift')[0].data), 6)
            self.assertFalse(result.get_saveframes_by_category('sample_conditions'))
            self.assertFalse(result.get_saveframes_by_category('chem_shift_reference'))
            for tag in ('Sample_condition_list_ID', 'Chem_shift_reference_ID'):
                self.assertTrue(all(value in (None, '.', '?') for value in shifts.get_tag(tag)))

    def test_nef_molecular_table_without_redundant_sequence_column(self):
        for optional in ({}, {'Seq_ID': '.'}, {'Seq_ID': '12'}):
            with self.subTest(optional=optional):
                utility = NmrVrptUtility(verbose=False, log=io.StringIO())
                reader = Mock()
                reader.getDataBlockNameList.return_value = ['entry']
                reader.hasCategory.side_effect = lambda category, _block: category == 'Chem_comp_assembly'
                reader.getDictList.return_value = [dict(Entity_assembly_ID='1', Comp_index_ID='1',
                                                        Comp_ID='ALA', Auth_seq_ID='12', **optional)]
                component = Mock()
                component.updateChemCompDict.return_value = False
                # Exercise the real extraction method with a bounded table reader.
                # Pylint cannot resolve mangled private fields on this external class.
                # pylint: disable=protected-access
                # __rR and __ccU are declared slots; Pylint misses their mangled spelling.
                # pylint: disable-next=assigning-non-slot
                utility._NmrVrptUtility__rR = reader
                # pylint: disable-next=assigning-non-slot
                utility._NmrVrptUtility__ccU = component
                self.assertTrue(utility._NmrVrptUtility__extractEntityInstancesFromNmr())
                key = '12' if optional.get('Seq_ID') == '12' else '1'
                self.assertIn((key, 'ALA'), utility._NmrVrptUtility__entityInstance['1'])
                # pylint: enable=protected-access

    def test_empty_required_shift_list_is_rejected_before_workflow(self):
        utility = NmrDpUtility()
        utility.addInput(name='chem_shift_file_path_list', value=[], type='file_dict_list')
        with self.assertRaisesRegex(ValueError, 'No input provided'):
            utility.op('nmr-cs-mr-merge')

    def test_unidentified_restraint_format_returns_diagnostics_without_crashing(self):
        coordinates = Path(__file__).resolve().parents[1] / 'tests-nmr/mock-data/2l9r.cif'
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / 'unknown.txt'
            source.write_text('not a recognized restraint format\n', encoding='utf-8')
            utility = NmrDpUtility()
            utility.setWorkspace(str(root / 'work'), str(root / 'cache'))
            utility.addInput(name='chem_shift_file_path_list',
                             value=[{'file_name': str(coordinates.with_suffix('.str')), 'file_type': 'nm-shi'}],
                             type='file_dict_list')
            utility.addInput(name='atypical_restraint_file_path_list',
                             value=[{'file_name': str(source), 'file_type': 'nm-res-bio'}], type='file_dict_list')
            utility.addInput(name='coordinate_file_path', value=str(coordinates), type='file')
            utility.addInput(name='remediation', value=True, type='param')
            utility.addInput(name='conversion_server', value=True, type='param')
            utility.setLog(str(root / 'report.json'))
            utility.setDestination(str(root / 'merged.str'))
            self.assertFalse(utility.op('nmr-cs-mr-merge'))
            self.assertTrue((root / 'report.json').is_file())
            report = json.loads((root / 'report.json').read_text(encoding='utf-8'))
            self.assertEqual(report['information']['input_sources'][1]['file_type'], 'nm-res-bio')
            self.assertTrue(report['error'])
            self.assertNotIn('internal_error', report['error'])

    def test_xeasy_creates_default_translator(self):
        reader = XeasyPKReader(verbose=False, log=io.StringIO())
        # The public API has no translator getter; verify the constructed dependency directly.
        # pylint: disable-next=protected-access
        self.assertIsInstance(reader._XeasyPKReader__nefT, NefTranslator)

    def test_coordinates_without_nmr_method_have_defined_default(self):
        cache = defaultdict(dict)
        cache.update(polymer_sequence=[], non_polymer=[], branched=[], entity_assembly=[])
        for metadata, expected in (([], '.'), ([{'method': 'X-RAY DIFFRACTION'}], '.'),
                                   ([{'method': 'SOLUTION NMR'}], 'SOLUTION NMR'),
                                   ([{'method': 'SOLID-STATE NMR'}], 'SOLID-STATE NMR')):
            with self.subTest(metadata=metadata):
                coordinates = Mock()
                coordinates.getDictList.return_value = metadata
                with patch('wwpdb.utils.nmr.pk.BasePKParserListener.coordAssemblyChecker', return_value=cache):
                    listener = BasePKParserListener(verbose=False, log=io.StringIO(), cR=coordinates, nefT=Mock())
                self.assertEqual(listener.exptlMethod, expected)

    def test_missing_list_identifiers_produce_diagnostic(self):
        for values in ([], ['.'], ['?', '.']):
            with self.subTest(values=values):
                entry = pynmrstar.Entry.from_scratch('test')
                frame = pynmrstar.Saveframe.from_scratch('shifts', '_Assigned_chem_shift_list')
                frame.add_tag('Sf_category', 'assigned_chemical_shifts')
                frame.add_tag('Sf_framecode', 'shifts')
                loop = pynmrstar.Loop.from_scratch('_Atom_chem_shift')
                loop.add_tag('Assigned_chem_shift_list_ID')
                for value in values:
                    loop.add_data([value])
                frame.add_loop(loop)
                entry.add_saveframe(frame)
                registry = MagicMock()
                registry.combined_mode = True
                registry.star_data = [entry]
                registry.report.input_sources[0].get.return_value = {
                    'file_name': 'missing-list.str', 'file_type': 'nmr-star', 'content_subtype': ['chem_shift']}
                registry.report.warning.getDescription.return_value = None
                utility = NmrDpUtility.__new__(NmrDpUtility)
                # __reg is a real mangled slot; Pylint does not resolve this external test assignment.
                # pylint: disable-next=protected-access,assigning-non-slot
                utility._NmrDpUtility__reg = registry
                self.assertFalse(utility._NmrDpUtility__appendParentSfTag())  # pylint: disable=protected-access
                args = registry.report.error.appendDescription.call_args.args
                self.assertEqual(args[0], 'missing_data')
                self.assertIn('Assigned_chem_shift_list_ID', args[1]['description'])
                self.assertEqual(frame.get_tag('ID'), [])

    def test_coordinate_alignment_uses_reference_sequence_without_auth_names(self):
        registry = MagicMock()
        model = [{'chain_id': 'A', 'seq_id': list(range(1, 32)), 'comp_id': ['ALA'] * 31},
                 {'chain_id': 'B', 'seq_id': [1, 2, 3], 'comp_id': ['GLY'] * 3}]
        nmr = [dict(model[0], chain_id='1'), dict(model[1], chain_id='2')]
        sources = [Mock(), Mock()]
        sources[0].get.return_value = {'polymer_sequence': nmr}
        sources[1].get.return_value = {'polymer_sequence': model}
        registry.report.input_sources = sources
        registry.report.getInputSourceIdOfCoord.return_value = 1
        registry.pA = PairwiseAlign()
        utility = NmrDpUtility.__new__(NmrDpUtility)
        # __reg is a real mangled slot; replace it to isolate alignment from file I/O.
        # pylint: disable-next=protected-access,assigning-non-slot
        utility._NmrDpUtility__reg = registry
        self.assertTrue(utility._NmrDpUtility__appendCoordPolymerSequenceAlignment())  # pylint: disable=protected-access
        reports = {call.args[0]: call.args[1] for call in registry.report.sequence_alignment.setItemValue.call_args_list}
        # Inspect actual alignment records, not a clamped percentage.
        records = reports['model_poly_seq_vs_nmr_poly_seq']
        self.assertEqual({(row['ref_chain_id'], row['test_chain_id']) for row in records}, {('A', '1'), ('B', '2')})
        self.assertTrue(all(row['sequence_coverage'] == 1.0 for row in records))

    def test_translation_keeps_setlog_path_for_the_next_phase(self):
        with TemporaryDirectory() as temporary:
            for method, output_key, suffix in (
                    ('_NmrDpUtility__initResourceForNef2Str', 'nmr-star_file_path', 'str'),
                    ('_NmrDpUtility__initResourceForStr2Nef', 'nef_file_path', 'nef')):
                for explicit in (False, True):
                    with self.subTest(method=method, explicit=explicit):
                        utility = NmrDpUtility(verbose=False, log=io.StringIO())
                        current = str(Path(temporary) / 'current.json')
                        override = str(Path(temporary) / 'override.json')
                        utility.setLog(current)
                        utility.addOutput(name=output_key, value=str(Path(temporary) / f'next.{suffix}'), type='file')
                        if explicit:
                            utility.addOutput(name='report_file_path', value=override, type='file')
                        self.assertTrue(getattr(utility, method)())
                        # Inspect the active workflow descriptor; no public getter exposes it.
                        # pylint: disable-next=protected-access
                        self.assertEqual(utility._NmrDpUtility__reg.inputParamDict['report_file_path'],
                                         override if explicit else current)

    def test_assembly_annotation_handles_missing_column_without_inventing_no(self):
        loop = pynmrstar.Loop.from_scratch('_Entity_assembly')
        loop.add_tag('ID')
        loop.add_data(['1'])
        loop.add_data(['2'])
        annotate_assembly_experimental_data(loop, {1})
        self.assertEqual(loop.get_tag('Experimental_data_reported'), ['yes', None])
        # Existing author assertions are retained even without matching extracted data.
        loop.data[1][1] = 'yes'
        annotate_assembly_experimental_data(loop, set())
        self.assertEqual(loop.get_tag('Experimental_data_reported'), ['yes', 'yes'])
