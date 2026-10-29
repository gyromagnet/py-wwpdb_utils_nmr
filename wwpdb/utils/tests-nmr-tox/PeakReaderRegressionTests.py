##
# File: PeakReaderRegressionTests.py
# Date: 29-Oct-2026
##
"""Regressions for independent peak readers and coordinates without experiment metadata."""
import io
import sys
import unittest
from collections import defaultdict
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
from wwpdb.utils.nmr.NmrDpUtility import NmrDpUtility
from wwpdb.utils.nmr.nef.NefTranslator import NefTranslator
from wwpdb.utils.nmr.pk.BasePKParserListener import BasePKParserListener
from wwpdb.utils.nmr.pk.XeasyPKReader import XeasyPKReader


class PeakReaderRegressionTests(unittest.TestCase):

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
