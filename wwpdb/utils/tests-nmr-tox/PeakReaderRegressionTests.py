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

from wwpdb.utils.nmr.NmrDpUtility import NmrDpUtility
from wwpdb.utils.nmr.nef.NefTranslator import NefTranslator
from wwpdb.utils.nmr.pk.BasePKParserListener import BasePKParserListener
from wwpdb.utils.nmr.pk.XeasyPKReader import XeasyPKReader


class PeakReaderRegressionTests(unittest.TestCase):

    def test_xeasy_creates_default_translator(self):
        reader = XeasyPKReader(verbose=False, log=io.StringIO())
        self.assertIsInstance(reader._XeasyPKReader__nefT, NefTranslator)  # pylint: disable=protected-access

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
                utility._NmrDpUtility__reg = registry  # pylint: disable=protected-access
                self.assertFalse(utility._NmrDpUtility__appendParentSfTag())  # pylint: disable=protected-access
                args = registry.report.error.appendDescription.call_args.args
                self.assertEqual(args[0], 'missing_data')
                self.assertIn('Assigned_chem_shift_list_ID', args[1]['description'])
                self.assertEqual(frame.get_tag('ID'), [])
