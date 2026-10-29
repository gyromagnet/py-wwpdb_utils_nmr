##
# XeasyPKReader.py
#
# Update:
##
""" A collection of classes for parsing XEASY PK files.
"""
__docformat__ = "restructuredtext en"
__author__ = "Masashi Yokochi"
__email__ = "yokochi@protein.osaka-u.ac.jp"
__license__ = "Apache License 2.0"
__version__ = "1.2.0"

import os
import sys
from typing import IO, List, Optional, Tuple

from antlr4 import ParseTreeWalker

try:
    from wwpdb.utils.nmr.NmrDpConstant import (MAX_ERROR_REPORT,
                                               REPRESENTATIVE_MODEL_ID,
                                               REPRESENTATIVE_ALT_ID)
    from wwpdb.utils.nmr.ChemCompUtil import ChemCompUtil
    from wwpdb.utils.nmr.BmrbChemShiftStat import BmrbChemShiftStat
    from wwpdb.utils.nmr.nef.NefTranslator import NefTranslator
    from wwpdb.utils.nmr.io.CifReader import CifReader
    from wwpdb.utils.nmr.AntlrParseUtil import parseAntlr
    from wwpdb.utils.nmr.pk import sa_xeasypk
    from wwpdb.utils.nmr.mr.LexerErrorListener import LexerErrorListener
    from wwpdb.utils.nmr.mr.ParserErrorListener import ParserErrorListener
    from wwpdb.utils.nmr.mr.ParserListenerUtil import coordAssemblyChecker
    from wwpdb.utils.nmr.pk.XeasyPKLexer import XeasyPKLexer
    from wwpdb.utils.nmr.pk.XeasyPKParser import XeasyPKParser
    from wwpdb.utils.nmr.pk.XeasyPKParserListener import XeasyPKParserListener
    from wwpdb.utils.nmr.pk.XeasyPROTReader import XeasyPROTReader
except ImportError:
    from nmr.NmrDpConstant import (MAX_ERROR_REPORT,
                                   REPRESENTATIVE_MODEL_ID,
                                   REPRESENTATIVE_ALT_ID)
    from nmr.ChemCompUtil import ChemCompUtil
    from nmr.BmrbChemShiftStat import BmrbChemShiftStat
    from nmr.nef.NefTranslator import NefTranslator
    from nmr.io.CifReader import CifReader
    from nmr.AntlrParseUtil import parseAntlr
    from nmr.pk import sa_xeasypk
    from nmr.mr.LexerErrorListener import LexerErrorListener
    from nmr.mr.ParserErrorListener import ParserErrorListener
    from nmr.mr.ParserListenerUtil import coordAssemblyChecker
    from nmr.pk.XeasyPKLexer import XeasyPKLexer
    from nmr.pk.XeasyPKParser import XeasyPKParser
    from nmr.pk.XeasyPKParserListener import XeasyPKParserListener
    from nmr.pk.XeasyPROTReader import XeasyPROTReader


class XeasyPKReader:
    """ Accessor methods for parsing XEASY PK files.
    """
    __slots__ = ('__class_name__',
                 '__version__',
                 '__verbose',
                 '__log',
                 '__debug',
                 '__enforcePeakRowFormat',
                 '__internal',
                 '__maxLexerErrorReport',
                 '__maxParserErrorReport',
                 '__representativeModelId',
                 '__representativeAltId',
                 '__mrAtomNameMapping',
                 '__ccU',
                 '__cR',
                 '__caC',
                 '__csStat',
                 '__nefT',
                 '__atomNumberDict',
                 '__reasons')

    def __init__(self, verbose: bool = True, log: IO = sys.stdout,
                 representativeModelId: int = REPRESENTATIVE_MODEL_ID,
                 representativeAltId: str = REPRESENTATIVE_ALT_ID,
                 mrAtomNameMapping: Optional[List[dict]] = None,
                 cR: Optional[CifReader] = None, caC: Optional[dict] = None, ccU: Optional[ChemCompUtil] = None,
                 csStat: Optional[BmrbChemShiftStat] = None, nefT: Optional[NefTranslator] = None,
                 atomNumberDict: Optional[dict] = None,
                 reasons: Optional[dict] = None) -> None:
        self.__class_name__ = self.__class__.__name__
        self.__version__ = __version__

        self.__verbose = verbose
        self.__log = log
        self.__debug = False
        self.__enforcePeakRowFormat = False
        self.__internal = False

        self.__maxLexerErrorReport = MAX_ERROR_REPORT
        self.__maxParserErrorReport = MAX_ERROR_REPORT

        self.__representativeModelId = representativeModelId
        self.__representativeAltId = representativeAltId
        self.__mrAtomNameMapping = mrAtomNameMapping

        # CCD accessing utility
        self.__ccU = ChemCompUtil(verbose, log) if ccU is None else ccU

        if cR is not None and caC is None:
            caC = coordAssemblyChecker(verbose, log, representativeModelId, representativeAltId,
                                       cR, self.__ccU, None, None, fullCheck=False)

        self.__cR = cR
        self.__caC = caC

        # BMRB chemical shift statistics
        self.__csStat = BmrbChemShiftStat(verbose, log, self.__ccU) if csStat is None else csStat

        # NefTranslator
        self.__nefT = NefTranslator(verbose, log, self.__ccU, self.__csStat) if nefT is None else nefT
        if nefT is None:
            self.__nefT.set_remediation_mode(True)

        # XeasyPROTParserListener.getAtomNumberDict()
        self.__atomNumberDict = atomNumberDict

        # reasons for re-parsing request from the previous trial
        self.__reasons = reasons

    def setDebugMode(self, debug: bool) -> None:
        """ Set debug mode.
        """

        self.__debug = debug

    def enforcePeakRowFormat(self, enforcePeakRowFormat: bool) -> None:
        """ Whether to enforce peak row format or not.
        """

        self.__enforcePeakRowFormat = enforcePeakRowFormat

    def setInternalMode(self, internal: bool) -> None:
        """ Set internal mode.
        """

        self.__internal = internal

    def setLexerMaxErrorReport(self, maxErrReport: int) -> None:
        """ Set the maximum number of lexer error messages to save.
        """

        self.__maxLexerErrorReport = maxErrReport

    def setParserMaxErrorReport(self, maxErrReport: int) -> None:
        """ Set the maximum number of parser error messages to save.
        """

        self.__maxParserErrorReport = maxErrReport

    def parse(self, pkFilePath: str, cifFilePath: Optional[str] = None, protFilePath: Optional[str] = None, isFilePath: bool = True,
              createSfDict: bool = False, originalFileName: Optional[str] = None, listIdCounter: Optional[dict] = None,
              reservedListIds: Optional[dict] = None, entryId: Optional[str] = None, csLoops: Optional[List[dict]] = None
              ) -> Tuple[Optional[XeasyPKParserListener], Optional[ParserErrorListener], Optional[LexerErrorListener]]:
        """ Parse XEASY PK file.
            @return: XeasyPKParserListener for success or None otherwise, ParserErrorListener, LexerErrorListener.
        """

        ifh = None

        try:

            if isFilePath:
                pkString = None

                if not os.access(pkFilePath, os.R_OK):
                    if self.__verbose:
                        self.__log.write(f"+{self.__class_name__}.parse() {pkFilePath} is not accessible.\n")
                    return None, None, None

                ifh = open(pkFilePath, 'r', encoding='utf-8', errors='ignore')  # pylint: disable=consider-using-with
                pkText = ifh.read()

            else:
                pkFilePath, pkString = None, pkFilePath

                if pkString is None or len(pkString) == 0:
                    if self.__verbose:
                        self.__log.write(f"+{self.__class_name__}.parse() Empty string.\n")
                    return None, None, None

                pkText = pkString

            if cifFilePath is not None:
                if not os.access(cifFilePath, os.R_OK):
                    if self.__verbose:
                        self.__log.write(f"+{self.__class_name__}.parse() {cifFilePath} is not accessible.\n")
                    return None, None, None

                if self.__cR is None:
                    self.__cR = CifReader(self.__verbose, self.__log)
                    if not self.__cR.parse(cifFilePath):
                        if self.__verbose:
                            self.__log.write(f"+{self.__class_name__}.parse() {cifFilePath} is not CIF file.\n")
                        return None, None, None

            if protFilePath is not None and self.__atomNumberDict is None:
                ptR = XeasyPROTReader(self.__verbose, self.__log,
                                      self.__representativeModelId,
                                      self.__representativeAltId,
                                      self.__mrAtomNameMapping,
                                      self.__cR, self.__caC,
                                      self.__ccU, self.__csStat, self.__nefT)
                protPL, _, _ = ptR.parse(protFilePath, cifFilePath)
                if protPL is not None:
                    self.__atomNumberDict = protPL.getAtomNumberDict()

            tree, parser_error_listener, lexer_error_listener =\
                parseAntlr(XeasyPKLexer, XeasyPKParser, 'xeasy_pk', pkText,
                           filePath=pkFilePath, saModule=sa_xeasypk,
                           maxLexerErrorReport=self.__maxLexerErrorReport,
                           maxParserErrorReport=self.__maxParserErrorReport,
                           ignoreCodicError=True)

            messageList = lexer_error_listener.getMessageList()

            if messageList is not None and self.__verbose:
                for description in messageList:
                    self.__log.write(f"[Syntax error] line {description['line_number']}:{description['column_position']} "
                                     f"{description['message']}\n")
                    if 'input' in description:
                        self.__log.write(f"{description['input']}\n")
                        self.__log.write(f"{description['marker']}\n")

            walker = ParseTreeWalker()
            listener = XeasyPKParserListener(self.__verbose, self.__log,
                                             self.__representativeModelId,
                                             self.__representativeAltId,
                                             self.__mrAtomNameMapping,
                                             self.__cR, self.__caC,
                                             self.__nefT,
                                             self.__atomNumberDict, self.__reasons)
            listener.debug = self.__debug
            listener.internal = self.__internal
            listener.createSfDict = createSfDict
            listener.enforcePeakRowFormat = self.__enforcePeakRowFormat
            if createSfDict:
                if originalFileName is not None:
                    listener.originalFileName = originalFileName
                if listIdCounter is not None:
                    listener.listIdCounter = listIdCounter
                if reservedListIds is not None:
                    listener.reservedListIds = reservedListIds
                if entryId is not None:
                    listener.entryId = entryId
                if csLoops is not None:
                    listener.csLoops = csLoops
            walker.walk(listener, tree)

            messageList = parser_error_listener.getMessageList()

            if messageList is not None and self.__verbose:
                for description in messageList:
                    self.__log.write(f"[Syntax error] line {description['line_number']}:{description['column_position']} "
                                     f"{description['message']}\n")
                    if 'input' in description:
                        self.__log.write(f"{description['input']}\n")
                        self.__log.write(f"{description['marker']}\n")

            if self.__verbose and self.__debug:
                if listener.warningMessage is not None and len(listener.warningMessage) > 0:
                    self.__log.write(f"+{self.__class_name__}.parse() ++ Info  -\n" + '\n'.join(listener.warningMessage) + '\n')
                if isFilePath:
                    self.__log.write(f"+{self.__class_name__}.parse() ++ Info  - {listener.getContentSubtype()}\n")

            return listener, parser_error_listener, lexer_error_listener

        except IOError as e:
            if self.__verbose:
                self.__log.write(f"+{self.__class_name__}.parse() ++ Error  - {str(e)}\n")
            return None, None, None
        finally:
            if isFilePath and ifh is not None:
                ifh.close()


if __name__ == "__main__":
    reader = XeasyPKReader(True)
    reader.setDebugMode(True)
    reader_listener, _, _ =\
        reader.parse('../../tests-nmr/mock-data-remediation/2khp/c13_noe.peaks',
                     '../../tests-nmr/mock-data-remediation/2khp/2khp.cif')

    reader = XeasyPKReader(True)
    reader.setDebugMode(True)
    reader_listener, _, _ =\
        reader.parse('../../tests-nmr/mock-data-remediation/5i2p/bmr30011/work/data/D_1000218151_nmr-peaks-upload_P2.dat.V1',
                     '../../tests-nmr/mock-data-remediation/5i2p/5i2p.cif')

    reader = XeasyPKReader(True)
    reader.setDebugMode(True)
    reader.setInternalMode(True)
    reader_listener, _, _ =\
        reader.parse('../../tests-nmr/mock-data-remediation/6lek/bmr36301/work/data/D_1300014649_nmr-peaks-upload_P1.dat.V1',
                     '../../tests-nmr/mock-data-remediation/6lek/6lek.cif')

    reader = XeasyPKReader(True)
    reader.setDebugMode(True)
    reader_listener, _, _ =\
        reader.parse('../../tests-nmr/mock-data-remediation/8eo9/bmr31044/work/data/D_1000268248_nmr-peaks-upload_P1.dat.V1',
                     '../../tests-nmr/mock-data-remediation/8eo9/8eo9.cif')

    reader = XeasyPKReader(True)
    reader.setDebugMode(True)
    reader_listener, _, _ =\
        reader.parse('../../tests-nmr/mock-data-remediation/7jyz/bmr30791/work/data/D_1000251588_nmr-peaks-upload_P2.dat.V1',
                     '../../tests-nmr/mock-data-remediation/7jyz/7jyz.cif')
    assert reader_listener.getReasonsForReparsing() is None

    reader = XeasyPKReader(True)
    reader.setDebugMode(True)
    reader.parse('../../tests-nmr/mock-data-remediation/2n6r/bmr25778/work/data/SF6_ali.peaks',
                 '../../tests-nmr/mock-data-remediation/2n6r/2n6r.cif')

    reader = XeasyPKReader(True)
    reader.setDebugMode(True)
    reader.parse('../../tests-nmr/mock-data-remediation/7cq1/bmr36368/work/data/D_1300017932_nmr-peaks-upload_P2.dat.V2',
                 '../../tests-nmr/mock-data-remediation/7cq1/7cq1.cif')

    reader = XeasyPKReader(True)
    reader.setDebugMode(True)
    reader.parse('../../tests-nmr/mock-data-remediation/7aaf/bmr34555/work/data/D_1292111069_nmr-peaks-upload_P1.dat.V2',
                 '../../tests-nmr/mock-data-remediation/7aaf/7aaf.cif',
                 '../../tests-nmr/mock-data-remediation/7aaf/bmr34555/work/data/D_1292111069_nmr-peaks-upload_P2.dat.V1')

    reader = XeasyPKReader(True)
    reader.setDebugMode(True)
    reader.parse('../../tests-nmr/mock-data-remediation/2kj5/2kj5-corrected-div_ext.mr',
                 '../../tests-nmr/mock-data-remediation/2kj5/2kj5.cif')

    reader = XeasyPKReader(True)
    reader.setDebugMode(True)
    reader.parse('../../tests-nmr/mock-data-remediation/7a2d/aro.peaks',  # -ignored-as-pea-any',
                 '../../tests-nmr/mock-data-remediation/7a2d/7a2d.cif',
                 '../../tests-nmr/mock-data-remediation/7a2d/all.prot')
