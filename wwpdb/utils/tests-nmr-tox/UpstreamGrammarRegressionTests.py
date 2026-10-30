##
# File: UpstreamGrammarRegressionTests.py
# Date: 30-Oct-2026
##
"""Rebased grammar changes must agree on native and Python parsing paths."""
import sys
import unittest

# Configure scientific imports before any NMR package is loaded.
if __package__ is None or __package__ == "":
    from os import path

    sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
    from commonsetup import TESTOUTPUT  # noqa: F401 pylint: disable=import-error,unused-import
else:
    from .commonsetup import TESTOUTPUT  # noqa: F401 pylint: disable=relative-beyond-top-level

from wwpdb.utils.nmr.AntlrParseUtil import parseAntlr, usingCppParser
from wwpdb.utils.nmr.mr import sa_ambermr, sa_cnsmr, sa_xplormr
from wwpdb.utils.nmr.mr.AmberMRLexer import AmberMRLexer
from wwpdb.utils.nmr.mr.AmberMRParser import AmberMRParser
from wwpdb.utils.nmr.mr.CnsMRLexer import CnsMRLexer
from wwpdb.utils.nmr.mr.CnsMRParser import CnsMRParser
from wwpdb.utils.nmr.mr.XplorMRLexer import XplorMRLexer
from wwpdb.utils.nmr.mr.XplorMRParser import XplorMRParser


class UpstreamGrammarRegressionTests(unittest.TestCase):

    def check_paths(self, lexer, parser, entry, shim, valid, invalid):
        # Pure-Python installs still exercise every syntax control. The local
        # qualification gate separately requires all native accelerators, so
        # both paths are exercised there without making optional binaries mandatory.
        implementations = (None, shim) if usingCppParser(shim) else (None,)
        trees = []
        for implementation in implementations:
            tree, parser_errors, lexer_errors = parseAntlr(lexer, parser, entry, valid, saModule=implementation)
            self.assertIsNone(parser_errors.getMessageList())
            self.assertIsNone(lexer_errors.getMessageList())
            trees.append(tree.getText())
            _, parser_errors, lexer_errors = parseAntlr(lexer, parser, entry, invalid, saModule=implementation)
            self.assertTrue(parser_errors.getMessageList() or lexer_errors.getMessageList())
        self.assertTrue(all(tree == trees[0] for tree in trees))

    def test_cns_keyword_class_name(self):
        self.check_paths(CnsMRLexer, CnsMRParser, 'cns_mr', sa_cnsmr,
                         'NOE CLASS BOND END\n', 'NOE CLASS BOND\n')

    def test_xplor_keyword_class_name(self):
        self.check_paths(XplorMRLexer, XplorMRParser, 'xplor_nih_mr', sa_xplormr,
                         'NOE CLASS BOND END\n', 'NOE CLASS BOND\n')

    def test_amber_alignment_parameters(self):
        self.check_paths(AmberMRLexer, AmberMRParser, 'amber_mr', sa_ambermr,
                         '&ALIGN FJ=1.0, DFAC_AL=1.0,2.0, /\n', '&ALIGN FJ=, /\n')

    def test_cns_averaging_class_remains_in_averaging_mode(self):
        for selection in ('*', 'NOE', 'DISTANCE_CLASS'):
            self.check_paths(CnsMRLexer, CnsMRParser, 'cns_mr', sa_cnsmr,
                             f'NOE AVERAGE {selection} R-6 END\n', f'NOE AVERAGE {selection} END\n')
