##
# File: test_NmrDpUtility.py
# Date:  26-Sep-2019  M. Yokochi
#
# Updates:
# 09-Oct-2019  M. Yokochi - add unit test for Xplor-NIH enabling 'check_mandatory_tag' option
# 28-Nov-2019  M. Yokochi - add unit test for NEF-Xplor-NIH-20191016-remediated.nef
# 26-Feb-2020  M. Yokochi - add unit test for CCPN_2mtv_docr.nef (DAOTHER-4785)
# 24-Mar-2020  M. Yokochi - add unit test for chemical shift reference (DAOTHER-1682)
# 14-May-2020  M. Yokochi - add unit tests for missing mandatory content (DAOTHER-5681 and 5682)
# 30-May-2020  M. Yokochi - add unit test for detection of total number of models (DAOTHER-5650)
# 08-Jul-2020  M. Yokochi - add unit test for combined NMR-STAR (DAOTHER-5926)
# 24-Aug-2021  M. Yokochi - add support for XPLOR-NIH planarity restraints (DAOTHER-7265)
# 27-Jan-2022  M. Yokochi - add restraint types described by XPLOR-NIH, CNS, CYANA, and AMBER systems (NMR restraint remediation)
# 04-Mar-2022  M. Yokochi - add coordinate geometry restraint (DAOTHER-7690, NMR restraint remediation)
# 30-Sep-2026  M. Yokochi - add standalone-mode deposit test for NmrVrptUtility's temporary files (DAOTHER-9785)
#
import glob
import json
import os
import shutil
import subprocess
import sys
import unittest

from testfixtures import LogCapture

# import commonsetup first: it mocks ConfigInfo, which ChemCompUtil reads when it is imported
if __package__ is None or __package__ == "":
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from commonsetup import TESTOUTPUT  # noqa:  F401 pylint: disable=import-error,unused-import
else:
    from .commonsetup import TESTOUTPUT  # noqa: F401 pylint: disable=relative-beyond-top-level

from wwpdb.utils.nmr.NmrDpUtility import NmrDpUtility

skipsome = True
if os.getenv("FULLTEST") is not None:
    skipsome = False


class TestNmrDpUtility(unittest.TestCase):
    def setUp(self):
        here = os.path.abspath(os.path.dirname(__file__))
        self.data_dir_path = os.path.join(here, os.pardir, "tests-nmr", "mock-data")
        self.utility = NmrDpUtility()

    def tearDown(self):
        pass

    # """
    # def test_init(self):
    #     nmr_content_subtypes = set(self.utility.nmr_content_subtypes)
    #
    #     self.assertEqual(nmr_content_subtypes, set(self.utility.sf_categories["nef"].keys()))
    #     self.assertEqual(nmr_content_subtypes, set(self.utility.sf_categories["nmr-star"].keys()))
    #
    #     self.assertEqual(nmr_content_subtypes, set(self.utility.lp_categories["nef"].keys()))
    #     self.assertEqual(nmr_content_subtypes, set(self.utility.lp_categories["nmr-star"].keys()))
    #
    #     # compare NMR content subtypes in NmrDpReportInputSource
    #     input_source = NmrDpReportInputSource()
    #     self.assertEqual(
    #         set(nmr_content_subtypes),
    #         set(input_source.content_subtypes)
    #         - {
    #             "plane_restraint",
    #             "adist_restraint",
    #             "rama_restraint",
    #             "radi_restraint",
    #             "diff_restraint",
    #             "nbase_restraint",
    #             "ang_restraint",
    #             "pre_restraint",
    #             "pcs_restraint",
    #             "prdc_restraint",
    #             "pang_restraint",
    #             "pccr_restraint",
    #             "hbond_restraint",
    #             "ssbond_restraint",
    #             "geo_restraint",
    #             "coordinate",
    #             "branched",
    #             "non_poly",
    #             "topology",
    #         },
    #     )
    #
    #     # data directory exists
    #     self.assertEqual(os.path.isdir(self.data_dir_path), True)
    # """

    def test_nmr_nef_consistency(self):
        # no input
        with LogCapture() as _logs:
            with self.assertRaises(ValueError):
                self.utility.op("nmr-nef-consistency-check")

        with LogCapture() as _logs:
            with self.assertRaises(IOError):
                self.utility.setSource("dummydummy")

        srcPath = os.path.join(self.data_dir_path, "2l9r.nef")
        self.utility.setSource(srcPath)
        self.utility.addInput(name="coordinate_file_path", value=srcPath, type="file")

        # invalid workflow operation
        with LogCapture() as _logs:  # noqa: F841
            with self.assertRaises(KeyError):
                self.utility.op("nmr")

        self.utility.setLog(os.path.join(TESTOUTPUT, "2l9r-nef-consistency-log.json"))

        self.utility.op("nmr-nef-consistency-check")

    def test_nmr_str_consistency(self):
        # no input
        with LogCapture() as _logs:
            with self.assertRaises(ValueError):
                self.utility.op("nmr-str-consistency-check")

        with LogCapture() as _logs:
            with self.assertRaises(IOError):
                self.utility.setSource("dummydummy")

        self.utility.setSource(os.path.join(self.data_dir_path, "2l9r.str"))
        self.utility.addInput(name="coordinate_file_path", value=os.path.join(self.data_dir_path, "2l9r.cif"), type="file")

        # invalid workflow operation
        with LogCapture() as _logs:  # noqa: F841
            with self.assertRaises(KeyError):
                self.utility.op("nmr")

        self.utility.setLog(os.path.join(TESTOUTPUT, "2l9r-str-consistency-log.json"))

        self.utility.op("nmr-str-consistency-check")

    def test_standalone_deposit(self):
        # Run the deposit operations the way the Docker image does, importing the package as bare nmr.* so that
        # every except ImportError branch is taken. Blocking wwpdb.utils.nmr is needed because CI installs the
        # package. The CIF output must survive the NmrVrptUtility run inside calculateOutputStats(), and no
        # temporary *.str2cif file may be left next to the outputs.
        utils_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
        work = os.path.join(TESTOUTPUT, "standalone-deposit")
        shutil.rmtree(work, ignore_errors=True)
        os.makedirs(work)
        for f in ("2l9r.str", "2l9r.cif"):
            shutil.copyfile(os.path.join(self.data_dir_path, f), os.path.join(work, f))

        script = """
import importlib.abc, os, sys
class BlockPackage(importlib.abc.MetaPathFinder):
    def find_spec(self, name, path, target=None):
        if name == 'wwpdb.utils.nmr' or name.startswith('wwpdb.utils.nmr.'):
            raise ImportError(name)
sys.meta_path.insert(0, BlockPackage())
utils_dir, work = sys.argv[1], sys.argv[2]
sys.path.insert(0, utils_dir)
from nmr.NmrDpUtility import NmrDpUtility
src, cif = os.path.join(work, '2l9r.str'), os.path.join(work, '2l9r.cif')
report = os.path.join(work, '2l9r-str-consistency-log.json')
u = NmrDpUtility()
u.setSource(src)
u.addInput(name='coordinate_file_path', value=cif, type='file')
u.setLog(report)
u.op('nmr-str-consistency-check')
for op, dst, nmr_cif in (('nmr-str2str-deposit', '2l9r_str2str.str', None),
                         ('nmr-str2cif-deposit', '2l9r_str2cif.str', '2l9r_str2cif.cif')):
    u = NmrDpUtility()
    u.setSource(src)
    u.addInput(name='coordinate_file_path', value=cif, type='file')
    u.addInput(name='report_file_path', value=report, type='file')
    u.setLog(os.path.join(work, op + '-log.json'))
    u.setDestination(os.path.join(work, dst))
    if nmr_cif is not None:
        u.addOutput(name='nmr_cif_file_path', value=os.path.join(work, nmr_cif), type='file')
    u.op(op)
"""
        proc = subprocess.run([sys.executable, "-c", script, utils_dir, work],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=1800, check=False)
        self.assertEqual(proc.returncode, 0, proc.stderr[-3000:])

        self.assertTrue(os.path.exists(os.path.join(work, "2l9r_str2str.str")))
        self.assertTrue(os.path.exists(os.path.join(work, "2l9r_str2cif.str")))
        self.assertTrue(os.path.exists(os.path.join(work, "2l9r_str2cif.cif")))
        self.assertEqual(glob.glob(os.path.join(work, "*.str2cif")), [])

    def test_release(self):
        # release() drops the parsed coordinates, and the next op() re-reads them with the same result;
        # nmr-cs*/nmr-str* operations release them on their own once they end
        cR = self.utility._NmrDpUtility__reg.cR  # pylint: disable=protected-access

        def check(log_name, src_name="2l9r.nef", op="nmr-nef-consistency-check"):
            log_path = os.path.join(TESTOUTPUT, log_name)
            self.utility.setSource(os.path.join(self.data_dir_path, src_name))
            self.utility.addInput(name="coordinate_file_path", value=os.path.join(self.data_dir_path, "2l9r.cif"), type="file")
            self.utility.setLog(log_path)
            self.utility.op(op)
            with open(log_path, "r", encoding="utf-8") as ifh:
                report = json.load(ifh)
            return report["information"]["status"], report["error"], report["warning"]

        check("2l9r-str-consistency-auto-release-log.json", "2l9r.str", "nmr-str-consistency-check")
        self.assertIsNone(cR.getDataBlock())

        before = check("2l9r-nef-consistency-before-release-log.json")
        self.assertIsNotNone(cR.getDataBlock())

        self.utility.release()
        self.assertIsNone(cR.getDataBlock())

        after = check("2l9r-nef-consistency-after-release-log.json")
        self.assertIsNotNone(cR.getDataBlock())
        self.assertEqual(before, after)

    def test_nmr_nef_consistency_check_non_std_residue(self):
        self.utility.setSource(os.path.join(self.data_dir_path, "2l9rnonstandard.nef"))
        self.utility.setLog(os.path.join(TESTOUTPUT, "2l9rnonstandard-nef-consistency-log.json"))

        self.utility.op("nmr-nef-consistency-check")

    @unittest.skipIf(skipsome is True, "Skip some tests")
    def test_nmr_nef_consistency_check_former_nef(self):
        self.utility.setSource(os.path.join(self.data_dir_path, "2l9rold.nef"))
        self.utility.setLog(os.path.join(TESTOUTPUT, "2l9rold-nef-consistency-log.json"))

        self.utility.op("nmr-nef-consistency-check")

    @unittest.skipIf(skipsome is True, "Skip some tests")
    def test_nmr_nef_consistency_check_cys(self):
        self.utility.setSource(os.path.join(self.data_dir_path, "2l9r-cys.nef"))
        self.utility.setLog(os.path.join(TESTOUTPUT, "2l9r-cys-nef-consistency-log.json"))

        self.utility.op("nmr-nef-consistency-check")

    @unittest.skipIf(skipsome is True, "Skip some tests")
    def test_nmr_nef_consistency_check_xplor_nih(self):
        self.utility.setSource(os.path.join(self.data_dir_path, "mth1743-test-20190919.nef"))
        self.utility.addInput(name="coordinate_file_path", value=os.path.join(self.data_dir_path, "1ryg.cif"), type="file")
        self.utility.addInput(name="resolve_conflict", value=True, type="param")
        self.utility.addInput(name="check_mandatory_tag", value=True, type="param")
        self.utility.setLog(os.path.join(TESTOUTPUT, "mth1743-test-20190919-nef-consistency-log.json"))

        self.utility.op("nmr-nef-consistency-check")

    @unittest.skipIf(skipsome is True, "Skip some tests")
    def test_nmr_nef_consistency_check_xplor_nih_remediated(self):
        self.utility.setSource(os.path.join(self.data_dir_path, "NEF-Xplor-NIH-20191016-remediated.nef"))
        self.utility.addInput(name="coordinate_file_path", value=os.path.join(self.data_dir_path, "1ryg.cif"), type="file")
        self.utility.addInput(name="resolve_conflict", value=True, type="param")
        self.utility.addInput(name="check_mandatory_tag", value=True, type="param")
        self.utility.setLog(os.path.join(TESTOUTPUT, "NEF-Xplor-NIH-20191016-remediated-nef-consistency-log.json"))

        self.utility.op("nmr-nef-consistency-check")

    @unittest.skipIf(skipsome is True, "Skip some tests")
    def test_nmr_nef_consistency_check_daother_4785(self):
        self.utility.setSource(os.path.join(self.data_dir_path, "CCPN_2mtv_docr.nef"))
        self.utility.addInput(name="coordinate_file_path", value=os.path.join(self.data_dir_path, "2mtv.cif"), type="file")
        self.utility.addInput(name="resolve_conflict", value=True, type="param")
        self.utility.addInput(name="check_mandatory_tag", value=True, type="param")
        self.utility.setLog(os.path.join(TESTOUTPUT, "ccpn_2mtv_docr-nef-consistency-log.json"))

        self.utility.op("nmr-nef-consistency-check")

    @unittest.skipIf(skipsome is True, "Skip some tests")
    def test_nmr_nef_consistency_check_daother_5681(self):
        self.utility.setSource(os.path.join(self.data_dir_path, "2k2e.nef-withoutrestraints"))
        self.utility.addInput(name="coordinate_file_path", value=os.path.join(self.data_dir_path, "2k2e.cif"), type="file")
        self.utility.addInput(name="resolve_conflict", value=True, type="param")
        self.utility.addInput(name="check_mandatory_tag", value=True, type="param")
        self.utility.setLog(os.path.join(TESTOUTPUT, "2k2e.nef-withoutrestraints-consistency-log.json"))
        self.utility.setVerbose(False)

        self.utility.op("nmr-nef-consistency-check")

    @unittest.skipIf(skipsome is True, "Skip some tests")
    def test_nmr_nef_consistency_check_daother_5682(self):
        self.utility.setSource(os.path.join(self.data_dir_path, "2k2e.nef-withoutchemicalshifts"))
        self.utility.addInput(name="coordinate_file_path", value=os.path.join(self.data_dir_path, "2k2e.cif"), type="file")
        self.utility.addInput(name="resolve_conflict", value=True, type="param")
        self.utility.addInput(name="check_mandatory_tag", value=True, type="param")
        self.utility.setLog(os.path.join(TESTOUTPUT, "2k2e.nef-withoutchemicalshifts-consistency-log.json"))
        self.utility.setVerbose(False)

        self.utility.op("nmr-nef-consistency-check")

    @unittest.skipIf(skipsome is True, "Skip some tests")
    def test_nmr_nef_consistency_check_daother_5682_modified(self):
        self.utility.setSource(os.path.join(self.data_dir_path, "2k2e.nef-withoutchemicalshifts-modified"))
        self.utility.addInput(name="coordinate_file_path", value=os.path.join(self.data_dir_path, "2k2e.cif"), type="file")
        self.utility.addInput(name="resolve_conflict", value=True, type="param")
        self.utility.addInput(name="check_mandatory_tag", value=True, type="param")
        self.utility.setLog(os.path.join(TESTOUTPUT, "2k2e.nef-withoutchemicalshifts-modified-consistency-log.json"))
        self.utility.setVerbose(False)

        self.utility.op("nmr-nef-consistency-check")

    @unittest.skipIf(skipsome is True, "Skip some tests")
    def test_nmr_nef_consistency_check_daother_5650(self):
        self.utility.setSource(os.path.join(self.data_dir_path, os.pardir, "NMR-VTF", "PDBStat_it2", "2k2e/2k2e.nef"))
        self.utility.addInput(name="coordinate_file_path",
                              value=os.path.join(self.data_dir_path, "2k2e.cif-only2models"),
                              type="file")
        self.utility.addInput(name="resolve_conflict", value=True, type="param")
        self.utility.addInput(name="check_mandatory_tag", value=True, type="param")
        self.utility.setLog(os.path.join(TESTOUTPUT, "2k2e.nef-only2models-consistency-log.json"))
        self.utility.setVerbose(False)

        self.utility.op("nmr-nef-consistency-check")

    @unittest.skipIf(skipsome is True, "Skip some tests")
    def test_nmr_str_consistency_check_bmrb_merged(self):
        self.utility.setSource(os.path.join(self.data_dir_path, "merged_30562_6nox.str"))
        self.utility.addInput(name="coordinate_file_path", value=os.path.join(self.data_dir_path, "6nox.cif"), type="file")
        self.utility.setLog(os.path.join(TESTOUTPUT, "merged_30562_6nox-str-consistency-log.json"))

        self.utility.op("nmr-str-consistency-check")

    def test_nmr_nef2str_deposit_str(self):
        nefSrcPath = os.path.join(self.data_dir_path, "2l9r.nef")
        modelSrcPath = os.path.join(self.data_dir_path, "2l9r.cif")
        consistJsonPath = os.path.join(TESTOUTPUT, "2l9r-nef2str-str-consistency-log.json")

        self.utility.setSource(nefSrcPath)
        self.utility.addInput(name="coordinate_file_path", value=modelSrcPath, type="file")
        self.utility.setLog(consistJsonPath)

        self.utility.op("nmr-nef-consistency-check")

        self.utility.setSource(nefSrcPath)
        self.utility.addInput(name="coordinate_file_path", value=modelSrcPath, type="file")
        self.utility.addInput(name="report_file_path", value=consistJsonPath, type="file")
        self.utility.setLog(os.path.join(TESTOUTPUT, "2l9r-nef2str-str-deposit-log.json"))
        self.utility.setDestination(os.path.join(TESTOUTPUT, "2l9r-nef2str-str.str"))
        self.utility.addOutput(name="nmr-star_file_path",
                               value=os.path.join(TESTOUTPUT, "2l9r-nef2str-str-nef2str.str"),
                               type="file")
        self.utility.addOutput(name="report_file_path",
                               value=os.path.join(TESTOUTPUT, "2l9r-nef2str-str-str-deposit-log.json"),
                               type="file")
        self.utility.addOutput(name="entry_id", value="2l9r", type="param")
        self.utility.addOutput(name="insert_entry_id_to_loops", value=True, type="param")
        self.utility.setVerbose(False)

        self.utility.op("nmr-nef2str-deposit")

    def test_nmr_str2str_deposit_nef(self):
        nefSrcPath = os.path.join(self.data_dir_path, "2l9r.nef")
        modelSrcPath = os.path.join(self.data_dir_path, "2l9r.cif")
        jsonLogPath = os.path.join(TESTOUTPUT, "2l9r-str2str-nef-consistency-log.json")

        self.utility.setSource(nefSrcPath)
        self.utility.addInput(name="coordinate_file_path", value=modelSrcPath, type="file")
        self.utility.setLog(jsonLogPath)

        self.utility.op("nmr-str-consistency-check")

        self.utility.setSource(nefSrcPath)
        self.utility.addInput(name="coordinate_file_path", value=modelSrcPath, type="file")
        self.utility.addInput(name="report_file_path", value=jsonLogPath, type="file")
        self.utility.setLog(os.path.join(TESTOUTPUT, "2l9r-str2str-nef-deposit-log.json"))
        self.utility.addOutput(name="entry_id", value="2l9r", type="param")
        self.utility.addOutput(name="insert_entry_id_to_loops", value=True, type="param")
        self.utility.setVerbose(False)

        self.utility.op("nmr-str2str-deposit")

    def test_nmr_nef2str_deposit(self):
        jsonLogPath = os.path.join(TESTOUTPUT, "2l9r-nef-consistency-log.json")
        if not os.access(jsonLogPath, os.F_OK):
            self.test_nmr_nef_consistency()

        self.utility.setSource(os.path.join(self.data_dir_path, "2l9r.nef"))
        self.utility.addInput(name="coordinate_file_path", value=os.path.join(self.data_dir_path, "2l9r.cif"), type="file")
        self.utility.addInput(name="report_file_path", value=jsonLogPath, type="file")
        self.utility.setLog(os.path.join(TESTOUTPUT, "2l9r-nef2str-deposit-log.json"))
        self.utility.setDestination(os.path.join(TESTOUTPUT, "2l9r-next.nef"))
        self.utility.addOutput(name="nmr-star_file_path", value=os.path.join(TESTOUTPUT, "2l9r-nef2str.str"), type="file")
        self.utility.addOutput(name="report_file_path",
                               value=os.path.join(TESTOUTPUT, "2l9r-nef2str-str-deposit-log.json"),
                               type="file")
        self.utility.addOutput(name="entry_id", value="NEED_ACC_NO", type="param")
        self.utility.setVerbose(False)

        self.utility.op("nmr-nef2str-deposit")

    @unittest.skipIf(skipsome is True, "Skip some tests")
    def test_nmr_nef2str_deposit_check_cys(self):
        jsonLogPath = os.path.join(TESTOUTPUT, "2l9r-cys-nef-consistency-log.json")
        if not os.access(jsonLogPath, os.F_OK):
            self.test_nmr_nef_consistency_check_cys()

        self.utility.setSource(os.path.join(self.data_dir_path, "2l9r-cys.nef"))
        self.utility.addInput(name="coordinate_file_path", value=os.path.join(self.data_dir_path, "2l9r-cys.cif"), type="file")
        self.utility.addInput(name="report_file_path", value=jsonLogPath, type="file")
        self.utility.setLog(os.path.join(TESTOUTPUT, "2l9r-cys-nef2str-deposit-log.json"))
        self.utility.setDestination(os.path.join(TESTOUTPUT, "2l9r-cys-next.nef"))
        self.utility.addOutput(name="nmr-star_file_path", value=os.path.join(TESTOUTPUT, "2l9r-cys-nef2str.str"), type="file")
        self.utility.addOutput(name="report_file_path",
                               value=os.path.join(TESTOUTPUT, "2l9r-cys-nef2str-str-deposit-log.json"),
                               type="file")
        self.utility.addOutput(name="entry_id", value="NEED_ACC_NO", type="param")
        self.utility.setVerbose(False)

        self.utility.op("nmr-nef2str-deposit")

    @unittest.skipIf(skipsome is True, "Skip some tests")
    def test_nmr_nef2str_deposit_check_xplor_nih(self):
        jsonLogPath = os.path.join(TESTOUTPUT, "mth1743-test-20190919-nef-consistency-log.json")
        if not os.access(jsonLogPath, os.F_OK):
            self.test_nmr_nef_consistency_check_xplor_nih()

        self.utility.setSource(os.path.join(self.data_dir_path, "mth1743-test-20190919.nef"))
        self.utility.addInput(name="coordinate_file_path", value=os.path.join(self.data_dir_path, "1ryg.cif"), type="file")
        self.utility.addInput(name="report_file_path", value=jsonLogPath, type="file")
        self.utility.addInput(name="resolve_conflict", value=True, type="param")
        self.utility.addInput(name="check_mandatory_tag", value=True, type="param")
        self.utility.setLog(os.path.join(TESTOUTPUT, "mth1743-test-20190919-nef2str-deposit-log.json"))
        self.utility.setDestination(os.path.join(TESTOUTPUT, "mth1743-test-20190919-next.nef"))
        self.utility.addOutput(name="nmr-star_file_path",
                               value=os.path.join(TESTOUTPUT, "mth1743-test-20190919-nef2str.str"),
                               type="file")
        self.utility.addOutput(name="report_file_path",
                               value=os.path.join(TESTOUTPUT, "mth1743-test-20190919-nef2str-str-deposit-log.json"),
                               type="file")
        self.utility.addOutput(name="entry_id", value="NEED_ACC_NO", type="param")
        self.utility.setVerbose(False)

        self.utility.op("nmr-nef2str-deposit")

    @unittest.skipIf(skipsome is True, "Skip some tests")
    def test_nmr_nef2str_deposit_check_daother_4785(self):
        jsonLogPath = os.path.join(TESTOUTPUT, "ccpn_2mtv_docr-nef-consistency-log.json")
        if not os.access(jsonLogPath, os.F_OK):
            self.test_nmr_nef_consistency_check_daother_4785()

        self.utility.setSource(os.path.join(self.data_dir_path, "CCPN_2mtv_docr.nef"))
        self.utility.addInput(name="coordinate_file_path", value=os.path.join(self.data_dir_path, "2mtv.cif"), type="file")
        self.utility.addInput(name="report_file_path", value=jsonLogPath, type="file")
        self.utility.addInput(name="resolve_conflict", value=True, type="param")
        self.utility.addInput(name="check_mandatory_tag", value=True, type="param")
        self.utility.setLog(os.path.join(TESTOUTPUT, "ccpn_2mtv_docr-nef2str-deposit-log.json"))
        self.utility.setDestination(os.path.join(TESTOUTPUT, "ccpn_2mtv_docr-next.nef"))
        self.utility.addOutput(name="nmr-star_file_path", value=os.path.join(TESTOUTPUT, "ccpn_2mtv_docr-nef2str.str"), type="file")
        self.utility.addOutput(name="report_file_path",
                               value=os.path.join(TESTOUTPUT, "ccpn_2mtv_docr-nef2str-str-deposit-log.json"),
                               type="file")
        self.utility.addOutput(name="entry_id", value="NEED_ACC_NO", type="param")
        self.utility.setVerbose(False)

        self.utility.op("nmr-nef2str-deposit")

    @unittest.skipIf(skipsome is True, "Skip some tests")
    def test_nmr_str2str_deposit(self):
        jsonLogPath = os.path.join(TESTOUTPUT, "2l9r-str-consistency-log.json")
        if not os.access(jsonLogPath, os.F_OK):
            self.test_nmr_str_consistency()

        self.utility.setSource(os.path.join(self.data_dir_path, "2l9r.str"))
        self.utility.addInput(name="coordinate_file_path", value=os.path.join(self.data_dir_path, "2l9r.cif"), type="file")
        self.utility.addInput(name="report_file_path", value=jsonLogPath, type="file")
        self.utility.setLog(os.path.join(TESTOUTPUT, "2l9r-str2str-deposit-log.json"))
        self.utility.setDestination(os.path.join(TESTOUTPUT, "2l9r-next.str"))
        self.utility.addOutput(name="entry_id", value="NEED_ACC_NO", type="param")
        self.utility.addOutput(name="insert_entry_id_to_loops", value=True, type="param")
        self.utility.setVerbose(False)

        self.utility.op("nmr-str2str-deposit")

    @unittest.skipIf(skipsome is True, "Skip some tests")
    def test_nmr_nef2str_deposit_check_non_std_residue(self):
        jsonLogPath = os.path.join(TESTOUTPUT, "2l9rnonstandard-nef-consistency-log.json")
        if not os.access(jsonLogPath, os.F_OK):
            self.test_nmr_nef_consistency_check_non_std_residue()

        self.utility.setSource(os.path.join(self.data_dir_path, "2l9rnonstandard.nef"))
        self.utility.addInput(name="coordinate_file_path", value=os.path.join(self.data_dir_path, "2l9r.cif"), type="file")
        self.utility.addInput(name="report_file_path", value=jsonLogPath, type="file")
        self.utility.setLog(os.path.join(TESTOUTPUT, "2l9rnonstandard-nef2str-deposit-log.json"))
        self.utility.setDestination(os.path.join(TESTOUTPUT, "2l9rnonstandard-next.nef"))
        self.utility.addOutput(name="nmr-star_file_path",
                               value=os.path.join(TESTOUTPUT, "2l9rnonstandard-nef2str.str"),
                               type="file")
        self.utility.addOutput(name="report_file_path",
                               value=os.path.join(TESTOUTPUT, "2l9rnonstandard-nef2str-str-deposit-log.json"),
                               type="file")
        self.utility.addOutput(name="entry_id", value="NEED_ACC_NO", type="param")
        self.utility.setVerbose(False)

        self.utility.op("nmr-nef2str-deposit")

    @unittest.skipIf(skipsome is True, "Skip some tests")
    def test_nmr_str2str_deposit_bmrb_merged(self):
        strSrcPath = os.path.join(self.data_dir_path, "merged_30562_6nox.str")
        modelSrcPath = os.path.join(self.data_dir_path, "6nox.cif")
        jsonLogPath = os.path.join(TESTOUTPUT, "merged_30562_6nox-str-consistency-log.json")

        self.utility.setSource(strSrcPath)
        self.utility.addInput(name="coordinate_file_path", value=modelSrcPath, type="file")
        self.utility.setLog(jsonLogPath)

        self.utility.op("nmr-str-consistency-check")

        self.utility.setSource(strSrcPath)
        self.utility.addInput(name="coordinate_file_path", value=modelSrcPath, type="file")
        self.utility.addInput(name="report_file_path", value=jsonLogPath, type="file")
        self.utility.addInput(name="nonblk_anomalous_cs", value=True, type="param")
        self.utility.addInput(name="nonblk_bad_nterm", value=True, type="param")
        self.utility.setLog(os.path.join(TESTOUTPUT, "merged_30562_6nox-str2str-deposit-log.json"))
        self.utility.setDestination(os.path.join(TESTOUTPUT, "merged_30562_6nox-next.str"))
        self.utility.addOutput(name="entry_id", value="6nox", type="param")
        self.utility.addOutput(name="insert_entry_id_to_loops", value=True, type="param")
        self.utility.setVerbose(False)

        self.utility.op("nmr-str2str-deposit")

    def test_nmr_str2str_deposit_cleaned(self):
        self.utility.setSource(os.path.join(self.data_dir_path, "24642_2mqq-clean.str"))
        self.utility.addInput(name="coordinate_file_path", value=os.path.join(self.data_dir_path, "2mqq.cif"), type="file")
        self.utility.setLog(os.path.join(TESTOUTPUT, "24642_2mqq-str2str-consistency-log.json"))

        self.assertFalse(self.utility.op("nmr-str-consistency-check"))

        self.utility.setSource(os.path.join(self.data_dir_path, "24642_2mqq-clean.str"))
        self.utility.addInput(name="coordinate_file_path", value=os.path.join(self.data_dir_path, "2mqq.cif"), type="file")
        self.utility.addInput(name="report_file_path",
                              value=os.path.join(TESTOUTPUT, "24642_2mqq-str2str-consistency-log.json"),
                              type="file")
        self.utility.setLog(os.path.join(TESTOUTPUT, "24642-str2str-deposit-log.json"))
        self.utility.setDestination(os.path.join(TESTOUTPUT, "24642_2mqq-clean-next.str"))
        self.utility.addOutput(name="entry_id", value="2mqq", type="param")
        self.utility.addOutput(name="insert_entry_id_to_loops", value=True, type="param")
        self.utility.setVerbose(False)

        self.assertFalse(self.utility.op("nmr-str2str-deposit"))
        # This historical "clean" input still contains anomalous assignments.
        # Exercise the conversion and require their diagnostics instead of skipping it
        # or confusing a generated file with a scientifically valid deposition.
        with open(os.path.join(TESTOUTPUT, "24642-str2str-deposit-log.json"), encoding="utf-8") as handle:
            report = json.load(handle)
        self.assertEqual(report['information']['status'], 'Error')
        self.assertGreater(len(report['error']['anomalous_data']), 0)
        self.assertFalse(report['error'].get('internal_error'))
        self.assertGreater(os.path.getsize(os.path.join(TESTOUTPUT, "24642_2mqq-clean-next.str")), 0)

    @unittest.skipIf(skipsome is True, "Skip some tests")
    def test_nmr_str2str_deposit_chem_shift_ref(self):
        self.utility.setSource(os.path.join(self.data_dir_path, "2la6-chem-shift-ref.str"))
        self.utility.addInput(name="coordinate_file_path", value=os.path.join(self.data_dir_path, "2la6.cif"), type="file")
        self.utility.setLog(os.path.join(TESTOUTPUT, "2la6-chem-shift-ref-consistency-log.json"))

        self.utility.op("nmr-str-consistency-check")

        self.utility.setSource(os.path.join(self.data_dir_path, "2la6-chem-shift-ref.str"))
        self.utility.addInput(name="coordinate_file_path", value=os.path.join(self.data_dir_path, "2la6.cif"), type="file")
        self.utility.addInput(name="report_file_path",
                              value=os.path.join(TESTOUTPUT, "2la6-chem-shift-ref-consistency-log.json"),
                              type="file")
        self.utility.setLog(os.path.join(TESTOUTPUT, "2la6-chem-shift-ref-deposit-log.json"))
        self.utility.setDestination(os.path.join(TESTOUTPUT, "2la6-chem-shift-ref-next.str"))
        self.utility.addOutput(name="entry_id", value="2la6", type="param")
        self.utility.addOutput(name="insert_entry_id_to_loops", value=True, type="param")
        self.utility.setVerbose(False)

        self.utility.op("nmr-str2str-deposit")

    @unittest.skipIf(skipsome is True, "Skip some tests")
    def test_nmr_str_consistency_check_review(self):
        self.utility.setSource(os.path.join(self.data_dir_path, "D_800107_nmr-data-str-review_P1.str.V20.rev"))
        self.utility.addInput(name="coordinate_file_path",
                              value=os.path.join(self.data_dir_path, "D_800107_model_P1.cif.V3"),
                              type="file")
        self.utility.setLog(os.path.join(TESTOUTPUT, "D_800107_nmr-data-str-consistency-log.json"))

        self.utility.op("nmr-str-consistency-check")

        # pynmrstar.Entry.from_file(self.data_dir_path + 'D_800107_nmr-data-str-review_P1.str.V20.rev')

    @unittest.skipIf(skipsome is True, "Skip some tests")
    def test_nmr_str_consistency_check_daother_5926(self):
        self.utility.setSource(os.path.join(self.data_dir_path, "swallow_NMR-Star_3-1.str"))
        self.utility.addInput(name="coordinate_file_path",
                              value=os.path.join(self.data_dir_path, "D_800365_model_P1.cif.V4"),
                              type="file")
        self.utility.addInput(name="nonblk_anomalous_cs", value=True, type="param")
        self.utility.addInput(name="nonblk_bad_nterm", value=True, type="param")
        self.utility.addInput(name="resolve_conflict", value=True, type="param")
        self.utility.addInput(name="check_mandatory_tag", value=True, type="param")
        self.utility.setLog(os.path.join(TESTOUTPUT, "daother-5926-str-consistency-log.json"))
        self.utility.setVerbose(False)

        self.utility.op("nmr-str-consistency-check")

    @unittest.skipIf(skipsome is True, "Skip some tests")
    def test_nmr_str2str_deposit_daother_5926(self):
        jsonLogPath = os.path.join(TESTOUTPUT, "daother-5926-str-consistency-log.json")
        if not os.access(jsonLogPath, os.F_OK):
            self.test_nmr_str_consistency_check_daother_5926()

        self.utility.setSource(os.path.join(self.data_dir_path, "swallow_NMR-Star_3-1.str"))
        self.utility.addInput(name="coordinate_file_path",
                              value=os.path.join(self.data_dir_path, "D_800365_model_P1.cif.V4"),
                              type="file")
        self.utility.addInput(name="report_file_path", value=jsonLogPath, type="file")
        self.utility.addInput(name="nonblk_anomalous_cs", value=True, type="param")
        self.utility.addInput(name="nonblk_bad_nterm", value=True, type="param")
        self.utility.addInput(name="resolve_conflict", value=True, type="param")
        self.utility.addInput(name="check_mandatory_tag", value=True, type="param")
        self.utility.setLog(os.path.join(TESTOUTPUT, "daother-5926-str2str-deposit-log.json"))
        self.utility.setDestination(os.path.join(TESTOUTPUT, "swallow_NMR-Star_3-1-next.str"))
        self.utility.setVerbose(False)

        self.utility.op("nmr-str2str-deposit")


if __name__ == "__main__":
    unittest.main()
