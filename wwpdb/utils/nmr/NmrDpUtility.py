##
# File: NmrDpUtility.py
# Date: 26-Sep-2019
#
# Updates:
# 10-Oct-2019  M. Yokochi - add 'check_mandatory_tag' option to detect missing mandatory tags as errors
# 15-Oct-2019  M. Yokochi - revise criteria on discrepancy in distance restraints using normalized value
# 01-Nov-2019  M. Yokochi - revise error message, instead of Python ValueError message
# 05-Nov-2019  M. Yokochi - revise error messages and detect empty sequence information
# 28-Nov-2019  M. Yokochi - fix saveframe name of nef_molecular_system and add 'nmr-str2nef-deposit' workflow operation
# 29-Nov-2019  M. Yokochi - relax allowable range of weight values in restraint data and support index pointer in auxiliary loops
# 11-Dec-2019  M. Yokochi - fix internal errors while processing NMR-VTF/PDBStat_examples and NMR-VTF/BMRB
# 24-Jan-2020  M. Yokochi - add histogram of distance restraints per residue and distance restraints on contact map
# 27-Jan-2020  M. Yokochi - add contact map for inter-chain distance restraints
# 28-Jan-2020  M. Yokochi - add struct_conf and struct_sheet_range data in dp report
# 29-Jan-2020  M. Yokochi - change plot type of dihedral angle and RDC restraints per residue
# 05-Feb-2020  M. Yokochi - add 'circular-shift' and 'void-zero' constraint for dihedral angle restraint
# 05-Feb-2020  M. Yokochi - move conflicted_data error to warning
# 07-Feb-2020  M. Yokochi - replace 'number_of_potential_types' by 'potential_type_of_constraints' in dp report
# 07-Feb-2020  M. Yokochi - allow multiple values in a data type on per residue plots
# 13-Feb-2020  M. Yokochi - add 'number_of_constraints_per_polymer_type' for apilayer.postModifyNMRRestraint
# 14-Feb-2020  M. Yokochi - add 'spectram_dim' for apilayer.postModifyNMRPeaks
# 21-Feb-2020  M. Yokochi - update content-type definitions and add release mode (nmr-str2nef-release workflow operation)
# 02-Mar-2020  M. Yokochi - add 'nmr-cs-nef-consistency-check' and 'nmr-cs-str-consistency-check' workflow operation (DAOTHER-4515)
# 05-Mar-2020  M. Yokochi - revise warning message (disordered_index) and enumerations (DAOTHER-5485)
# 06-Mar-2020  M. Yokochi - fix invalid ambiguity_code while parsing
# 13-Mar-2020  M. Yokochi - revise error/warning messages
# 17-Mar-2020  M. Yokochi - add 'undefined' value for potential_type (DAOTHER-5508)
# 17-Mar-2020  M. Yokochi - revise warning message about enumeration mismatch for potential_type and restraint_origin (DAOTHER-5508)
# 17-Mar-2020  M. Yokochi - check total number of models (DAOTHER-436)
# 17-Mar-2020  M. Yokochi - check consistency between saveframe name and sf_framecode value
# 18-Mar-2020  M. Yokochi - rename warning type from skipped_sf/lp_category to skipped_saveframe/loop_category
# 18-Mar-2020  M. Yokochi - support 'Saveframe' data type as conventional NMR data (DAOTHER-2737)
# 19-Mar-2020  M. Yokochi - atom nomenclature should not become a blocker (DAOTHER-5527)
# 24-Mar-2020  M. Yokochi - add support for chemical shift reference (DAOTHER-1682)
# 24-Mar-2020  M. Yokochi - revise chain assignment for identical dimer case (DAOTHER-3343)
# 30-Mar-2020  M. Yokochi - preserve original sf_framecode for nef_molecular_system (NEF) or assembly (NMR-STAR)
# 31-Mar-2020  M. Yokochi - enable processing without log file
# 03-Apr-2020  M. Yokochi - preserve case code of atom_name (NEF) and Auth_atom_ID/Original_PDB_atom_name (NMR-STAR)
# 06-Apr-2020  M. Yokochi - synchronize with coordinates' auth_asym_id and auth_seq_id for combined NMR-STAR deposition
# 10-Apr-2020  M. Yokochi - fix crash in case of format issue
# 14-Apr-2020  M. Yokochi - fix dependency on label_seq_id, instead of using auth_seq_id in case (DAOTHER-5584)
# 18-Apr-2020  M. Yokochi - fix no model error in coordinate and allow missing 'sf_framecode' in NMR conventional deposition
#                           (DAOTHER-5594)
# 19-Apr-2020  M. Yokochi - support concatenated CS data in NMR conventional deposition (DAOTHER-5594)
# 19-Apr-2020  M. Yokochi - report warning against not superimposed models (DAOTHER-4060)
# 22-Apr-2020  M. Yokochi - convert comp_id in capital letters (DAOTHER-5600)
# 22-Apr-2020  M. Yokochi - fix GLY:HA1/HA2 to GLY:HA2/HA3 (DAOTHER-5600)
# 22-Apr-2020  M. Yokochi - fix ambiguity code mismatch if possible (DAOTHER-5601)
# 22-Apr-2020  M. Yokochi - fix None type object is not iterable error (DAOTHER-5602)
# 23-Apr-2020  M. Yokochi - support conventional atom name for methyl group without wildcard character, e.g. ALA:HB (DAOTHER-5603)
# 23-Apr-2020  M. Yokochi - change missing ambiguity_set_id error to warning (DAOTHER-5609)
# 23-Apr-2020  M. Yokochi - make sure to parse chem_shift_ref saveframe tag (DAOTHER-5610)
# 23-Apr-2020  M. Yokochi - implement automatic format correction (DAOTHER-5603, 5610)
# 24-Apr-2020  M. Yokochi - separate format_issue error and missing_mandatory_content error (DAOTHER-5611)
# 24-Apr-2020  M. Yokochi - support 'QR' pseudoatom name (DAOTHER-5611)
# 24-Apr-2020  M. Yokochi - allow mandatory value is missing in NMR conventional deposition (DAOTHER-5611)
# 25-Apr-2020  M. Yokochi - implement automatic format correction for 6NZN, 6PQF, 6PSI entry (DAOTHE-5611)
# 25-Apr-2020  M. Yokochi - add 'entity' content subtype (DAOTHER-5611)
# 25-Apr-2020  M. Yokochi - add 'corrected_format_issue' warning type (DAOTHER-5611)
# 27-Apr-2020  M. Yokochi - add 'auth_atom_nomenclature_mismatch' warning type (DAOTHER-5611)
# 27-Apr-2020  M. Yokochi - implement recursive format corrections (DAOTHER-5602)
# 28-Apr-2020  M. Yokochi - copy the normalized CS/MR files if output file path list is set (DAOTHER-5611)
# 28-Apr-2020  M. Yokochi - catch 'range-float' error as 'unusual data' warning (DAOTHER-5611)
# 28-Apr-2020  M. Yokochi - extract sequence from CS/MR loop with gap (DAOTHER-5611)
# 29-Apr-2020  M. Yokochi - support diagnostic message of PyNMRSTAR v2.6.5.1 or later (DAOTHER-5611)
# 29-Apr-2020  M. Yokochi - implement more automatic format corrections with PyNMRSTAR v2.6.5.1 (DAOTHER-5611)
# 29-Apr-2020  M. Yokochi - fix different CS warning between NEF and NMR-STAR (DAOTHER-5621)
# 29-Apr-2020  M. Yokochi - add 'number_of_constraint_sets' of experiment data in report (DAOTHER-5622)
# 29-Apr-2020  M. Yokochi - sort 'conflicted_data' and 'inconsistent_data' warning items (DAOTHER-5622)
# 30-Apr-2020  M. Yokochi - allow NMR conventional atom naming scheme in NMR-STAR V3.2 (DAOTHER-5634)
# 01-May-2020  M. Yokochi - allow NMR conventional atom naming scheme in NMR-STAR V3.2 (DAOTHER-5634)
# 02-May-2020  M. Yokochi - additional support for format issue correction while STAR to NEF conversion (DAOTHER-5577)
# 02-May-2020  M. Yokochi - re-implement basic mathematical functions using Numpy library
# 07-May-2020  M. Yokochi - revise warning type (from 'insufficient_data' to 'encouragement')
#                           if total number of models is less than 8 (DAOTHER-5650)
# 07-May-2020  M. Yokochi - add preventive code for infinite loop while format issue correction
# 08-May-2020  M. Yokochi - sync update with wwpdb.utils.nmr.CifReader (DAOTHER-5654)
# 09-May-2020  M. Yokochi - add support for submitted coordinate file (allow missing of pdbx_poly_seq_scheme) (DAOTHER-5654)
# 12-May-2020  M. Yokochi - fix diselenide bond detection
# 14-May-2020  M. Yokochi - fix error detection for missing mandatory content (DAOTHER-5681 and 5682)
# 15-May-2020  M. Yokochi - add 'content_mismatch' error for NMR legacy deposition (DAOTHER-5687)
# 15-May-2020  M. Yokochi - revise encouragement message if total number of models is less than 5 (DAOTHER-5650)
# 16-May-2020  M. Yokochi - block NEF file upload in NMR legacy deposition (DAOTHER-5687)
# 30-May-2020  M. Yokochi - refer to atom_site to get total number of models (DAOTHER-5650)
# 01-Jun-2020  M. Yokochi - let RMSD cutoff value configurable (DAOTHER-4060)
# 05-Jun-2020  M. Yokochi - be compatible with wwpdb.utils.align.alignlib using Python 3 (DAOTHER-5766)
# 06-Jun-2020  M. Yokochi - be compatible with pynmrstar v3 (DAOTHER-5765)
# 12-Jun-2020  M. Yokochi - overall performance improvement by reusing cached data and code revision
# 19-Jun-2020  M. Yokochi - do not generate invalid restraints that include self atom
# 26-Jun-2020  M. Yokochi - add support for covalent bond information (_nef_covalent_links and _Bond categories)
# 30-Jun-2020  M. Yokochi - ignore third party loops and items gracefully (DAOTHER-5896)
# 30-Jun-2020  M. Yokochi - prevent pynmrstar's exception due to empty string (DAOTHER-5894)
# 08-Jul-2020  M. Yokochi - bug fix release for DAOTHER-5926
# 09-Jul-2020  M. Yokochi - add support for categories in NMR-STAR specific peak list (DAOTHER-5926)
# 09-Jul-2020  M. Yokochi - adjust arguments of pynmrstar write_to_file() to prevent data losses (v2.6.1, DAOTHER-5926)
# 17-Aug-2020  M. Yokochi - add support for residue variant (DAOTHER-5906)
# 20-Aug-2020  M. Yokochi - add 'leave_intl_note' output parameter decides whether to leave internal commentary note
#                           in processed NMR-STAR file,
#                           set False for OneDep environment (DAOTHER-6030)
# 10-Sep-2020  M. Yokochi - add 'transl_pseudo_name' input parameter decides whether to translate conventional pseudoatom
#                           nomenclature in combined NMR-STAR file (DAOTHER-6128)
# 16-Sep-2020  M. Yokochi - bug fix release based on internal test using BMRB NMR restraint archive of 6.3k entries (DAOTHER-6128)
# 18-Sep-2020  M. Yokochi - bug fix release for negative sequence numbers (DAOTHER-6128)
# 25-Sep-2020  M. Yokochi - add 'tolerant_seq_align' input parameter which enables tolerant sequence alignment for residue variant,
#                           set False for OneDep environment (DAOTHER-6128)
# 09-Oct-2020  M. Yokochi - support circular chain_id re-mapping with seq_id shifts in data loops if it is necessary,
#                           'tolerant_seq_align' input parameter is required (DAOTHER-6128)
# 22-Oct-2020  M. Yokochi - run diagnostic routine for case of sequence mismatch between defined polymer sequence
#                           and sequence in data loop (DAOTHER-6128)
# 11-Nov-2020  M. Yokochi - set NEF v1.1 as the default specification
# 12-Nov-2020  M. Yokochi - improve NMR warning messages (DAOTHER-6109, 6167)
# 18-Nov-2020  M. Yokochi - fix calculation of CS completeness, fix empty polymer_sequence_in_loop
#                           due to atom_site.pdbx_PDB_ins_code (DAOTHER-6128)
# 20-Nov-2020  M. Yokochi - rename 'remarkable_data' warning category to 'unusual/rare_data' (DAOTHER-6372)
# 26-Nov-2020  M. Yokochi - detect the nearest ferromagnetic atom, in addition to paramagnetic atom (DAOTHER-6366)
# 27-Nov-2020  M. Yokochi - add support for non-IUPAC atom names for standard amino acids, e.g. ARG:HB1/HB2 -> HB2/HB3
#                           (DAOTHER-6373)
# 17-Dec-2020  M. Yokochi - support 'atom_not_found' error with message revision (DAOTHER-6345)
# 25-Jan-2021  M. Yokochi - simplify code for Entity_assemble_ID and chain_code
# 25-Jan-2021  M. Yokochi - add CS validation code about rotameric state of ILE/LEU/VAL residue
# 03-Feb-2021  M. Yokochi - update polymer sequence which shares the same entity and missing
#                           in the molecular assembly information if necessary,
#                           e.g. double stranded DNA (DAOTHER-6128, BMRB entry: 16812, PDB ID: 6kae)
# 10-Mar-2021  M. Yokochi - block NEF deposition missing '_nef_sequence' category and turn off salvage routine for the case
#                           (DAOTHER-6694)
# 10-Mar-2021  M. Yokochi - add support for audit loop in NEF (DAOTHER-6327)
# 12-Mar-2021  M. Yokochi - add diagnostic routine to fix inconsistent sf_framecode of conventional CS file (DAOTHER-6693)
# 14-May-2021  M. Yokochi - add support for PyNMRSTAR v3.1.1 (DAOTHER-6693)
# 20-May-2021  M. Yokochi - fix duplicating pynmrstar data objects during format issue correction
#                           that leads to empty upload summary page (DAOTHER-6834)
# 24-May-2021  M. Yokochi - fix tautomer detection of coordinate (DAOTHER-6809)
# 17-Jun-2021  M. Yokochi - fix error in handling lower/upper linear limits (DAOTHER-6963)
# 17-Jun-2021  M. Yokochi - relax tolerance on chemical shift difference (DAOTHER-6963)
# 23-Jun-2021  M. Yokochi - send back the initial error message when format remediation fails (DAOTHER-6830)
# 25-Jun-2021  M. Yokochi - block restraint files that have no distance restraints (DAOTHER-6830)
# 28-Jun-2021  M. Yokochi - support cif-formatted CS file for reupload without changing CS data (DAOTHER-6830, 7097)
# 29-Jun-2021  M. Yokochi - include auth_asym_id in NMR data processing report (DAOTHER-7108)
# 29-Jun-2021  M. Yokochi - add support for PyNMRSTAR v3.2.0 (DAOTHER-7107)
# 02-Jul-2021  M. Yokochi - detect content type of AMBER restraint file and AMBER auxiliary file (DAOTHER-6830, 1901)
# 12-Jul-2021  M. Yokochi - add RCI validation code for graphical representation of NMR data
# 24-Aug-2021  M. Yokochi - detect content type of XPLOR-NIH planarity restraints (DAOTHER-7265)
# 10-Sep-2021  M. Yokochi - prevent system crash for an empty loop case of CS/MR data (D_1292117593)
# 13-Oct-2021  M. Yokochi - fix/adjust tolerances for spectral peak list (DAOTHER-7389, issue #1 and #2)
# 13-Oct-2021  M. Yokochi - code revision according to PEP8 using Pylint (DAOTHER-7389, issue #5)
# 14-Oct-2021  M. Yokochi - remove unassigned chemical shifts, clear incompletely assigned spectral peaks (DAOTHER-7389, issue #3)
# 27-Oct-2021  M. Yokochi - fix collection of unmapped sequences and utilize Auth_asym_ID* tag for chain_id
#                           if Entity_assembly_ID* is not available (DAOTHER-7421)
# 28-Oct-2021  M. Yokochi - resolve case-insensitive saveframe name collision for CIF (DAOTHER-7389, issue #4)
# 16-Nov-2021  M. Yokochi - fix sequence conflict in case that large sequence gap in CS implies multi chain complex (DAOTHER-7465)
# 16-Nov-2021  M. Yokochi - fix server crash with disulfide bond, which is not supported by chemical shifts (DAOTHER-7475)
# 16-Nov-2021  M. Yokochi - revised error message for malformed XPLOR-NIH RDC restraints (DAOTHER-7478)
# 18-Nov-2021  M. Yokochi - detect content type of XPLOR-NIH hydrogen bond geometry restraints (DAOTHER-7478)
# 18-Nov-2021  M. Yokochi - relax detection of distance restraints for nm-res-cya and nm-res-oth (DAOTHER-7491)
# 13-Dec-2021  M. Yokochi - append sequence spacer between large gap to prevent failure of sequence alignment
#                           (DAOTHER-7465, issue #2)
# 14-Dec-2021  M. Yokochi - report detailed warning message against not superimposed models and exactly overlaid models
#                           (DAOTHER-4060, 7544)
# 15-Dec-2021  M. Yokochi - fix server crash while uploading NMR restraint file in NMR-STAR format (DAOTHER-7545)
# 21-Dec-2021  M. Yokochi - fix wrong missing_mandatory_content error when uploading NMR restraint files in NMR-STAR format
#                           (DAOTHER-7545, issue #2)
# 14-Jan-2022  M. Yokochi - report exactly overlaid models in the coordinate file (DAOTHER-7544)
# 17-Feb-2022  M. Yokochi - aware of presence of _atom_site.pdbx_auth_atom_name for N-terminal protonation change
#                           while upload-conversion of the coordinate file (DAOTHER-7665)
# 17-Feb-2022  M. Yokochi - do report incompletely assigned chemical shifts for conventional deposition (DAOTHER-7662)
# 21-Feb-2022  M. Yokochi - verify 'onebond' coherence transfer type using CCD (DAOTHER-7681, issue #2)
# 21-Feb-2022  M. Yokochi - verify pseudoatom names in NMR restraints are in assigned chemical shifts (DAOTHER-7681, issue #1)
# 24-Mar-2022  M. Yokochi - utilize software specific MR parsers for sanity check of NMR restraint files (DAOTHER-7690)
# 20-Mar-2022  M. Yokochi - add support for _atom_site.label_alt_id (DAOTHER-4060, 7544, NMR restraint remediation)
# 06-Apr-2022  M. Yokochi - detect other possible MR format if the first parsing fails (DAOTHER-7690)
# 02-May-2022  M. Yokochi - implement recursive MR splitter guided by MR parsers (NMR restraint remediation)
# 17-May-2022  M. Yokochi - add support for BIOSYM MR format (DAOTHER-7825, NMR restraint remediation)
# 01-Jun-2022  M. Yokochi - add support for GROMACS PT/MR format (DAOTHER-7769, NMR restraint remediation)
# 17-Jun-2022  M. Yokochi - add support for DYNAMO/PALES/TALOS MR format (DAOTHER-7872, NMR restraint remediation)
# 06-Jul-2022  M. Yokochi - add support for SYBYL MR format (DAOTHER-7902, NMR restraint remediation)
# 05-Aug-2022  M. Yokochi - do not add a saveframe tag if there is already the tag (DAOTHER-7947)
# 31-Aug-2022  M. Yokochi - separate atom_not_found error and hydrogen_not_instantiated error (NMR restraint remediation)
# 06-Sep-2022  M. Yokochi - add support for branched entity (NMR restraint remediation)
# 13-Sep-2022  M. Yokochi - add 'nm-res-isd' file type for IDS (inference structure determination) restraint format
#                           (DAOTHER-8059, NMR restraint remediation)
# 22-Sep-2022  M. Yokochi - add 'nm-res-cha' file type for CHARMM restraint format (DAOTHER-8058, NMR restraint remediation)
# 20-Oct-2022  M. Yokochi - report recommendation message when there is no distance restraints for NMR deposition,
#                           instead of blocker (DAOTHER-8088 1.b, 8108)
# 24-Oct-2022  M. Yokochi - add support for floating chiral stereo assignments (NMR restraint remediation)
# 15-Dec-2022  M. Yokochi - merge CS and MR as a single NMR data file in CIF format with comprehensive molecular assembly
#                           information (DAOTHER-7407, NMR restraint remediation)
# 13-Jan-2023  M. Yokochi - add support for small angle X-ray scattering restraints (NMR restraint remediation)
# 24-Jan-2023  M. Yokochi - add support for heteronuclear relaxation data (NOE, T1, T2, T1rho, Order parameter)
#                           (NMR restraint remediation)
# 23-Feb-2023  M. Yokochi - combine spectral peak lists in any plain text format into single NMR-STAR until Phase 2 release
#                           (DAOTHER-7407)
# 24-Mar-2023  M. Yokochi - add 'nmr-nef2cif-deposit' and 'nmr-str2cif-deposit' workflow operations (DAOTHER-7407)
# 22-Jun-2023  M. Yokochi - convert model file when pdbx_poly_seq category is missing for reuploading nmr_data after unlock
#                           (DAOTHER-8580)
# 19-Jul-2023  M. Yokochi - fix not to merge restraint id (_Gen_dist_constraint.ID) if lower and upper bounds are different
#                           (DAOTHER-8705)
# 20-Jul-2023  M. Yokochi - throw 'format_issue' error when polymer sequence extraction fails (DAOTHER-8644)
# 09-Aug-2023  M. Yokochi - remediate combined nmr_data by default and improve robustness of sequence alignment (DAOTHER-8751)
# 13-Sep-2023  M. Yokochi - construct pseudo CCD from the coordinates (DAOTHER-8817)
# 29-Sep-2023  M. Yokochi - add 'nmr-str2cif-annotation' workflow operation (DAOTHER-8817, 8828)
# 02-Oct-2023  M. Yokochi - do not reorganize _Gen_dist_constraint.ID of native combined NMR data (DAOTHER-8855)
# 10-Nov-2023  M. Yokochi - raise a content mismatch error properly for spectral peak list when the file is irrelevant
#                           (DAOTHER-8949)
# 13-Dec-2023  M. Yokochi - add 'hydrogen_non_instantiated' warning (DAOTHER-8945)
# 11-Jan-2024  M. Yokochi - convert RTF to ASCII file if necessary (DAOTHER-9063)
# 12-Jan-2024  M. Yokochi - preserve the original sequence offset of CS loop of UNMAPPED residue (DAOTHER-9065)
# 12-Jan-2024  M. Yokochi - fix sequence merge of entity loop and CS loop (DAOTHER-9065)
# 16-Jan-2024  M. Yokochi - add 'nm-res-ari' file type for ARIA restraint format (DAOTHER-9079, NMR restraint remediation)
# 17-Jan-2024  M. Yokochi - detect coordinate issue (DAOTHER-9084, type_symbol mismatches label_atom_id)
# 24-Jan-2024  M. Yokochi - reconstruct polymer/non-polymer sequence based on pdb_mon_id, instead of auth_mon_id (D_1300043061)
# 21-Feb-2024  M. Yokochi - add support for discontinuous model_id (NMR restraint remediation, 2n6j)
# 07-Mar-2024  M. Yokochi - extract pdbx_poly_seq_scheme.auth_mon_id as alt_cmop_id to prevent sequence mismatch
#                           due to 5-letter CCD ID (DAOTHER-9158 vs D_1300043061)
# 22-Mar-2024  M. Yokochi - test tautomeric states of histidine-like residue across models (DAOTHER-9252)
# 01-May-2024  M. Yokochi - merge cs/mr sequence extensions containing unknown residues (e.g UNK, DN, N) if necessary
#                           (NMR restraint remediation, 6fw4)
# 22-May-2024  M. Yokochi - block deposition using a peak list file in any binary format and prevent 'nm-pea-any'
#                           occasionally matches with 'nm-res-cya' (DAOTHER-9425)
# 11-Jun-2024  M. Yokochi - add support for ligand remapping in annotation process (DAOTHER-9286)
# 25-Jun-2024  M. Yokochi - strip white spaces in a datablock name derived from the model file (DAOTHER-9511)
# 28-Jun-2024  M. Yokochi - ignore extraneous input value for numeric tags and replace statistics of chemical shifts
#                           using remediated loop (DAOTHER-9520)
# 28-Jun-2024  M. Yokochi - convert conventional NMR name of carbonyl carbon 'CO' to valid one 'C' (DAOTHER-9520, 2nd case)
# 19-Aug-2024  M. Yokochi - fill 'Data_file_name' saveframe tag by default (DAOTHER-9520, 4th case)
# 20-Aug-2024  M. Yokochi - support truncated loop sequence in the model (DAOTHER-9644)
# 14-Nov-2024  M. Yokochi - add support for CHARMM extended CRD (topology) file. file type: 'nm-aux-cha'
# 19-Nov-2024  M. Yokochi - add support for pH titration data (NMR restraint remediation)
# 22-Nov-2024  M. Yokochi - add support for CYANA NOA (NOE Assignment) file. file type: 'nm-res-noa'
#                           (DAOTHER-7829, 9785, NMR data remediation)
# 27-Nov-2024  M. Yokochi - implement atom name mapping history as requirement of standalone NMR data conversion service
# 28-Nov-2024  M. Yokochi - drop support for old pynmrstar versions less than 3.2
# 16-Dec-2024  M. Yokochi - combine spectral peak lists written in software native formats into single NMR-STAR file
#                           (DAOTHER-8905, NMR data remediation Phase 2)
# 19-Dec-2024  M. Yokochi - add 'nmr-if-merge-deposit' workflow operation (DAOTHER-8905, NMR data remediation Phase 2)
# 27-Dec-2024  M. Yokochi - extract NMRIF metadata from NMR-STAR (DAOTHER-1728, 9846)
# 09-Jan-2025  M. Yokochi - extract NMRIF metadata from NMR-STAR (as primary source) and model (as secondary source)
#                           (DAOTHER-1728, 9846)
# 31-Jan-2025  M. Yokochi - add 'coordinate_issue' and 'assigned_peak_atom_not_found' warnings
#                           used in NMR data remediation with peak list (DAOTHER-8905, 9785)
# 07-Feb-2025  M. Yokochi - add support for 'ignore_error' attribute in addInput() for test processing of spectral peak list files
#                           derived from legacy ADIT system (DAOTHER-8905)
# 13-Feb-2025  M. Yokochi - set _Spectral_dim_transfer.Type 'through-space?' temporarily and resolve it
#                           based on related experiment type (DAOTHER-8905, 1728, 9846)
# 17-Feb-2025  M. Yokochi - discard remediated spectral peak list in OneDep environment (DAOTHER-8905, 9785)
# 18-Feb-2025  M. Yokochi - add support for PONDEROSA spectral peak list (DAOTHER-8905, 9785)
# 26-Feb-2025  M. Yokochi - add support for CCPN tabular spectral peak list (DAOTHER-8905, 9785)
# 05-Mar-2025  M. Yokochi - add support for bare spectral peak list (DAOTHER-8905, 9785)
# 06-Mar-2025  M. Yokochi - add support for coupling constant data (NMR data remediation Phase 2)
# 28-Mar-2025  M. Yokochi - add support for SPARKY's 'save' (aka. ornament) peak list
#                           (DAOTHER-8905, 9785, NMR data remediation Phase 2)
# 09-Apr-2025  M. Yokochi - enable to convert chemical shifts in any software-native format
#                           in the standalone NMR data conversion service (v4.4.0, DAOTHER-9785)
# 23-Apr-2025  M. Yokochi - enable to inherit previous warnings/errors to report out failed restraint conversions
#                           and detailed messages (DAOTHER-9785)
# 29-May-2025  M. Yokochi - analyze spectral peak list files and provide warning/error messages to depositor
#                           (DAOTHER-8905, 8949, 10096, 10097, 10098, 10099, 10100, 10101)
# 11-Jun-2025  M. Yokochi - reconstruct atom name mapping from revision history and PDB Versioned Archive if possible
#                           (DAOTHER-7829, 8905)
# 25-Jul-2025  M. Yokochi - enable to configure whether to enforce to use _Peak_row_format loop for spectral peak list remediation
#                           (DAOTHER-8905, 9785)
# 06-Aug-2025  M. Yokochi - add support for SCHRODINGER/ASL MR format (DAOTHER-7902, 10172, NMR data remediation)
# 22-Aug-2025  M. Yokochi - add support for OLIVIA spectral peak list (DAOTHER-8905, 9785)
# 22-Aug-2025  M. Yokochi - add 'nm-shi-oli' file type for OLIVIA spectral peak list file (DAOTHER-9785)
# 19-Sep-2025  M. Yokochi - add 'nm-aux-pdb' file type for bare PDB file as AMBER/CHARMM/GROMACS topology
#                           (DAOTHER-7829, 9785, NMR data remediation)
# 02-Oct-2025  M. Yokochi - add support for ARIA NOE restraint (XML) file. file type: 'nm-res-arx'
#                           (DAOTHER-7829, 9785, NMR data remediation)
# 02-Oct-2025  M. Yokochi - add 'nm-res-bar' file type for Bare WSV/TSV/CSV restraint file
#                           (DAOTHER-7829, 9785, NMR data remediation)
# 20-Oct-2025  M. Yokochi - enable to parse concatenated notation of chain code and sequence code in ROSETTA restraints
#                           (DAOTHER-7829, 9785, NMR data remediation)
# 21-Oct-2025  M. Yokochi - enable to parse concatenated notation of chain code and sequence code in CYANA restraints
#                           (DAOTHER-7829, 9785, NMR data remediation)
# 19-Nov-2025  M. Yokochi - add 'nmr-str-replace-cs' workflow operation (DATAQUALITY-2178, NMR data remediation)
# 03-Dec-2025  M. Yokochi - split concatenation of auth_seq_id and ins_code in CS/MR loops (DAOTHER-10418, NMR data remediation)
# 07-Jan-2026  M. Yokochi - code refactoring: NmrDpConstant, NmrDpRegistry, NmrDpMrSplitter, NmrDpFirstAid, NmrDpValidation,
#                           and NmrDpRemediation classes, v5.0.0)
# 27-Jan-2026  M. Yokochi - raise error when entity exists and sequence inconsistency between the entity and loops,
#                           instead of warning, do not remediate CS loop in case of the sequence mismatch error (DAOTHER-10487)
# 16-Jun-2026  M. Yokochi - add setWorkspace() method to set current working directory and cache file directory (DAOTHER-9785)
# 30-Jun-2026  M. Yokochi - add support for chemical shift perturbation experiment by adding 'nm-csp-*' file types (DAOTHER-9785)
# 09-Jul-2026  M. Yokochi - implement BMRB's data provenance check in standalone NMR data conversion service (DAOTHER-9785)
# 13-Jul-2026  M. Yokochi - implement ensemble composition analysis including cluster analysis (DAOTHER-9785)
# 24-Jul-2026  M. Yokochi - explain exact copy of the multimer's chemical shifts using comma-separated Auth_asym_IDs (DAOTHER-10898)
# 07-Aug-2026  M. Yokochi - bridge ANTLR4 C++ lexer/parser into the Python listeners (speedy-antlr-tool, DAOTHER-9785)
# 17-Aug-2026  M. Yokochi - code refactoring: split NmrDpValidation and NmrDpRemediation classes into topic modules
#                           (NmrDpValidation{Base,Input,Coord,Loop,Nomencl,Cs,Mr,Pk,CoordChk,CsStats,MrStats,OutStats}
#                           and NmrDpRemediation{Base,Enum,PolySeq,CsLoop,Cs,Mr,Pk,LegacyCs,LegacyMr,LegacyPk,Stats,Merge}),
#                           no API change
# 24-Aug-2026  M. Yokochi - replace the deepcopy() serialization round-trip by minimal copies at the call sites
#                           (copyFactor/copyPolySeq/atomKey in ParserListenerUtil, list()/dict()/copy.copy(),
#                           or no copy at all where the retrieved object is immutable or never mutated),
#                           and compare atom selections of the 'not' clause via a hashed key set (DAOTHER-10315)
# 24-Sep-2026  M. Yokochi - add an optional C accelerator (cpp_src/c_listener_util.c, built by
#                           WWPDB_NMR_BUILD_C_ACCEL=1) for ParserListenerUtil.copyFactor/copyPolySeq/atomKey,
#                           replace the str(factor) key of the factor cache by ParserListenerUtil.factorKey(),
#                           and store the SCHRODINGER 'store' clause as a factor dictionary (DAOTHER-7829, 9785)
# 25-Sep-2026  M. Yokochi - clear the shared PairwiseAlign before each pairwise alignment, which otherwise re-aligns
#                           every test sequence ever added, and remove quadratic lookups in
#                           AlignUtil.alignPolymerSequence/assignPolymerSequence (DAOTHER-7829, 9785)
# 28-Sep-2026  M. Yokochi - return the glibc heap freed by the speedy-antlr C++ lexer/parser to the OS right after
#                           each parse (AntlrParseUtil.parseAntlr), which otherwise stays resident through the
#                           listener walk (DAOTHER-7829, 9785)
# 28-Sep-2026  M. Yokochi - collect the cyclic garbage of the previous NMR unified data before validating a large one
#                           (NmrDpValidationInput.validateInputSource), which otherwise sets the peak memory of
#                           'nmr-str2str-deposit' (DAOTHER-7829, 9785)
# 29-Sep-2026  M. Yokochi - keep only the run_gc() calls measured to free memory: a full collection at the registry
#                           reset of op(), trim_heap() alone after a C++ parse, and none at the end of
#                           NmrDpRemediationMerge.mergeLegacyData(), NmrVrptUtility.op() or CifReader.__calculateRmsd(),
#                           where nothing was left to collect (DAOTHER-7829, 9785)
# 29-Sep-2026  M. Yokochi - add release() to drop the parsed coordinates between workflow operations, e.g. before
#                           running NmrVrptUtility in the same process (DAOTHER-7829, 9785)
# 29-Sep-2026  M. Yokochi - collect the model numbers and the representative alt_id of atom_site with
#                           CifReader.getDistinctValues()/getFirstValue() in NmrDpUtility and NmrVrptUtility, instead of
#                           a dictionary per atom row (+263 MB for 1.36 M rows) (DAOTHER-7829, 9785)
# 29-Sep-2026  M. Yokochi - hand the NMR-STAR entry parsed by NefTranslator.validate_file() over to the following
#                           read_input_file() of the same, unchanged large file, instead of parsing it again (DAOTHER-7829, 9785)
# 30-Sep-2026  M. Yokochi - merge atom selections with a size-adaptive hashed index (appendUniqueAtoms() in
#                           BaseStackedMRParserListener) in exitSelection() of the XPLOR-NIH, CNS, CHARMM and SCHRODINGER
#                           listeners and in doIntersectionFactor_expressions(), instead of O(N*M) list scans (DAOTHER-7829, 9785)
# 01-Oct-2026  M. Yokochi - the speedy-antlr C++ accelerators lex a str stored one byte per character in place
#                           (cpp_src/latin1_input_stream.h), instead of a UTF-32 copy of the input (DAOTHER-7829, 9785)
# 02-Oct-2026  M. Yokochi - make the output independent of PYTHONHASHSEED: deduplicate atom selections in first-occurrence
#                           order, take the first atom as the reference of isAmbigAtomSelection(), keep the data item order
#                           in NefTranslator.check_data(), order coordinate chains and comp_ids, and NmrVrptUtility's
#                           chains and residues by first occurrence;
#                           reset distance_sub_type and bond_flag per restraint in NmrVrptUtility's dist_violation_seq,
#                           which carried over from the previous restraint when no atom matched the residue (DAOTHER-7829, 8905)
# 02-Oct-2026  M. Yokochi - hash the text of input files in chunks (NmrDpConstant.text_md5()) instead of holding
#                           the whole file twice, in CifReader, NmrDpMrSplitter and NmrDpValidationOutStats (DAOTHER-7829, 9785)
# 05-Oct-2026  M. Yokochi - pass chem_comp_type, the type_symbol of each atom of a pseudo CCD, to NefTranslator, and rebuild
#                           a cached coordinate assembly checker that lacks it (DAOTHER-8817)
# 07-Oct-2026  M. Yokochi - fold BaseStackedMRParserListener's validate{Rdc,T1T2,Csa,Pre,Pcs,Ccr}Range() into one helper,
#                           and move selectRealistic{Bond,Chi2Angle}Constraint(), duplicated in both MR listener bases,
#                           to NmrVrptUtility (DAOTHER-7829)
# 07-Oct-2026  M. Yokochi - delegate the copies of selectRealistic{Bond,Chi2Angle}Constraint() in Amber, Rosetta and
#                           Gromacs MR parser listeners to the same NmrVrptUtility functions (DAOTHER-7829)
# 08-Oct-2026  M. Yokochi - share getSfDict()/trimSfWoLp() and the local sequence offset searches of MR parser listeners in
#                           ParserListenerUtil, and fold repeated blocks of BaseStackedMRParserListener's exit() and
#                           doConsumeFactor_expressions() into nested functions (DAOTHER-7829)
# 08-Oct-2026  M. Yokochi - trim the last saveframe without any row of MR parser listeners (trimSfWoLp()), which never
#                           matched, so a later list took a list id that the next file reused; keep the last saveframes
#                           per listener instance (DAOTHER-7829)
# 09-Oct-2026  M. Yokochi - merge assignCoordPolymerSequenceWithChainIdWithoutCompId() into ...WithoutCompId(), share
#                           getSfDictOf(), and fold repeated code of BasePKParserListener (validatePeak{2,3,4}D(), J-coupling/
#                           relayed transfer remediation, extractPeakAssignment()) and BaseCSParserListener
#                           (extractChemShiftAssignment()) (DAOTHER-7829)
# 13-Oct-2026  M. Yokochi - performance: skip the NMR-STAR -> CIF -> CifReader round trip of NmrVrptUtility, remove
#                           per-residue scans in CifReader.getPolymerSequence() and coordAssemblyChecker(), and patch two
#                           hot spots of the ANTLR Python runtime in AntlrParseUtil (CPU time, 6x63: -74%, 6cc9: -22%,
#                           6pts: -16%, DAOTHER-10315: -6%) (DAOTHER-7829, 9785)
# 14-Oct-2026  M. Yokochi - admit a compiled accelerator only when it matches its Python code (AccelUtil.py), and build
#                           binary wheels with the accelerators next to the pure-Python package (cibuildwheel.toml) (DAOTHER-9785)
# 16-Oct-2026  M. Yokochi - add 'permit_missing_cs_file' input parameter (DAOTHER-7829, 8905, NMR data remediation Phase 3)
# 21-Oct-2026  M. Yokochi - fix system clash whole handling case of missing _Spectral_dim.Atom_isotope_number tag (DAOTHER-11129)
# 22-Oct-2026  M. Yokochi - enable in-place remediation in nmr-str-consistency-check workflow
#                           when sequence mismatch errors are raised for cs and mr loops (DAOTHER-11129)
##
""" Main class for NMR data processing.
    @author: Masashi Yokochi
"""
__docformat__ = "restructuredtext en"
__author__ = "Masashi Yokochi"
__email__ = "yokochi@protein.osaka-u.ac.jp"
__license__ = "Apache License 2.0"
__version__ = "5.4.4"

import collections
import copy
import gc
import hashlib
import itertools
import os
import re
import shutil
import sys
import time
from operator import itemgetter
from typing import Any, IO, List, Optional, Union

from mmcif.io.IoAdapterPy import IoAdapterPy

from munkres import Munkres

import pynmrstar

from wwpdb.utils.align.alignlib import PairwiseAlign  # pylint: disable=no-name-in-module

try:
    from wwpdb.utils.nmr.NmrDpConstant import (MODEL_FILE_PATH_KEY,
                                               ALT_MODEL_FILE_PATH_KEY,
                                               CS_FILE_PATH_LIST_KEY,
                                               MR_FILE_PATH_LIST_KEY,
                                               AR_FILE_PATH_LIST_KEY,
                                               AC_FILE_PATH_LIST_KEY,
                                               REPORT_FILE_PATH_KEY,
                                               NMR_CIF_FILE_PATH_KEY,
                                               NMRIF_FILE_PATH_KEY,
                                               NEXT_NEF_FILE_PATH_KEY,
                                               NEXT_STAR_FILE_PATH_KEY,
                                               DP_INPUT_PARAM_KEYS,
                                               DP_INPUT_FILE_KEYS,
                                               DP_INPUT_FILE_LIST_KEYS,
                                               DP_INPUT_FILE_DICT_KEYS,
                                               DP_OUTPUT_PARAM_KEYS,
                                               DP_OUTPUT_FILE_KEYS,
                                               DP_OUTPUT_FILE_LIST_KEYS,
                                               DP_WORKFLOW_OPS,
                                               SUB_DIR_NAME_FOR_CACHE,
                                               DEFAULT_ENTRY_ID,
                                               INITIAL_ENTRY_ID,
                                               NMR_CONTENT_SUBTYPES,
                                               MR_CONTENT_SUBTYPES,
                                               PK_CONTENT_SUBTYPES,
                                               CIF_CONTENT_SUBTYPES,
                                               READABLE_FILE_TYPE,
                                               CONTENT_TYPE,
                                               SF_CATEGORIES,
                                               LP_CATEGORIES,
                                               MR_MAX_SPACER_LINES,
                                               INDEX_TAGS,
                                               WEIGHT_TAGS,
                                               ANGLE_TYPE_TAGS,
                                               CONSIST_ID_TAGS,
                                               PK_KEY_ITEMS,
                                               DATA_ITEMS,
                                               NUM_DIM_ITEMS,
                                               ALLOWED_TAGS,
                                               DISALLOWED_PK_TAGS,
                                               ERR_TEMPLATE_FOR_MISSING_MANDATORY_LP_TAG,
                                               SF_TAG_PREFIXES,
                                               SF_TAG_ITEMS,
                                               MANDATORY_SF_TAG_ITEMS,
                                               SF_ALLOWED_TAGS,
                                               WARN_TEMPLATE_FOR_MISSING_MANDATORY_SF_TAG,
                                               INSTRUCTION_FOR_FULL_SEQUENCE,
                                               AUX_LP_CATEGORIES,
                                               LINKED_LP_CATEGORIES,
                                               AUX_ALLOWED_TAGS,
                                               ITEM_NAMES_IN_CS_LOOP,
                                               ITEM_NAMES_IN_PK_LOOP,
                                               ITEM_NAMES_IN_DIST_LOOP,
                                               ITEM_NAMES_IN_DIHED_LOOP,
                                               ITEM_NAMES_IN_RDC_LOOP,
                                               CS_LIST_SF_TAG_NAME,
                                               LOW_SEQ_COVERAGE,
                                               MIN_SEQ_COVERAGE_W_CONFLICT,
                                               LARGE_ASYM_ID,
                                               LEN_MAJOR_ASYM_ID,
                                               EMPTY_VALUE,
                                               TRUE_VALUE,
                                               STD_MON_DICT,
                                               PROTON_BEGIN_CODE,
                                               PARAMAGNETIC_ELEMENTS,
                                               FERROMAGNETIC_ELEMENTS,
                                               MAX_DIM_NUM_OF_SPECTRA,
                                               ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS,
                                               MAX_CONFLICT_ATTEMPT,
                                               PDB_ID_PAT,
                                               BMRB_ID_PAT,
                                               CNV_ID_PAT,
                                               WORK_MODEL_FILE_NAME_PAT,
                                               CHK_DESC_PAT,
                                               CHK_DESC_ONE_PAT,
                                               CHK_DESC_MAND_PAT,
                                               CHK_DESC_MAND_ONE_PAT,
                                               ARCHIVAL_MR_FILE_TYPES,
                                               PARSABLE_PK_FILE_TYPES,
                                               CS_RANGE_MIN,
                                               CS_RANGE_MAX,
                                               CS_UNCERT_MAX,
                                               REPRESENTATIVE_MODEL_ID,
                                               REPRESENTATIVE_ASYM_ID,
                                               REPRESENTATIVE_ALT_ID,
                                               SPECTRAL_DIM_TEMPLATE,
                                               DEFAULT_COORD_PROPERTIES)
    from wwpdb.utils.nmr.NmrDpRegistry import (NmrDpRegistry,
                                               get_next_path,
                                               test_path_with_suffix)
    from wwpdb.utils.nmr.NmrDpFirstAid import NmrDpFirstAid
    from wwpdb.utils.nmr.NmrDpMrSplitter import (NmrDpMrSplitter,
                                                 detect_bom,
                                                 convert_codec,
                                                 is_binary_file,
                                                 get_type_of_star_file,
                                                 get_peak_list_format)
    from wwpdb.utils.nmr.NmrDpValidation import (NmrDpValidation,
                                                 predict_redox_state_of_cysteine,
                                                 is_like_planality_boundary)
    from wwpdb.utils.nmr.NmrDpRemediation import (NmrDpRemediation,
                                                  get_chem_shift_format_from_string)
    from wwpdb.utils.nmr.NmrDpReport import (NmrDpReport,
                                             NmrDpReportInputSource)
    from wwpdb.utils.nmr.ChemCompUtil import ChemCompUtil
    from wwpdb.utils.nmr.BmrbChemShiftStat import BmrbChemShiftStat
    from wwpdb.utils.nmr.AlignUtil import (hasLargeInnerSeqGap,
                                           hasLargeSeqGap,
                                           fillInnerBlankCompId,
                                           fillBlankCompId,
                                           fillBlankCompIdWithOffset,
                                           beautifyPolySeq,
                                           getMiddleCode,
                                           getGaugeCode,
                                           getScoreOfSeqAlign,
                                           getOneLetterCodeCan,
                                           getOneLetterCodeCanSequence,
                                           alignPolymerSequence,
                                           alignPolymerSequenceWithConflicts,
                                           assignPolymerSequence,
                                           retrieveAtomNameMappingFromRevisions,
                                           retrieveAtomNameMappingFromInternal)
    from wwpdb.utils.nmr.CifToNmrStar import (CifToNmrStar,
                                              has_key_value,
                                              get_first_sf_tag,
                                              set_sf_tag)
    from wwpdb.utils.nmr.NmrVrptUtility import (uncompress_gzip_file,
                                                load_from_pickle,
                                                write_as_pickle,
                                                get_temp_path)
    from wwpdb.utils.nmr.nef.NefTranslator import NefTranslator
    from wwpdb.utils.nmr.io.CifReader import CifReader
    from wwpdb.utils.nmr.io.PdbxUtil import abandon_symbolic_labels
    from wwpdb.utils.nmr.mr.ParserListenerUtil import (coordAssemblyChecker,
                                                       isAmbigAtomSelection,
                                                       getTypeOfDihedralRestraint,
                                                       isLikeHis)
    from wwpdb.utils.nmr.pk.BasePKParserListener import guess_primary_dim_transfer_type
    from wwpdb.utils.nmr.ann.OneDepAnnTasks import OneDepAnnTasks
    from wwpdb.utils.nmr.ann.BmrbAnnTasks import BmrbAnnTasks
except ImportError:
    from nmr.NmrDpConstant import (MODEL_FILE_PATH_KEY,
                                   ALT_MODEL_FILE_PATH_KEY,
                                   CS_FILE_PATH_LIST_KEY,
                                   MR_FILE_PATH_LIST_KEY,
                                   AR_FILE_PATH_LIST_KEY,
                                   AC_FILE_PATH_LIST_KEY,
                                   REPORT_FILE_PATH_KEY,
                                   NMR_CIF_FILE_PATH_KEY,
                                   NMRIF_FILE_PATH_KEY,
                                   NEXT_NEF_FILE_PATH_KEY,
                                   NEXT_STAR_FILE_PATH_KEY,
                                   DP_INPUT_PARAM_KEYS,
                                   DP_INPUT_FILE_KEYS,
                                   DP_INPUT_FILE_LIST_KEYS,
                                   DP_INPUT_FILE_DICT_KEYS,
                                   DP_OUTPUT_PARAM_KEYS,
                                   DP_OUTPUT_FILE_KEYS,
                                   DP_OUTPUT_FILE_LIST_KEYS,
                                   DP_WORKFLOW_OPS,
                                   SUB_DIR_NAME_FOR_CACHE,
                                   DEFAULT_ENTRY_ID,
                                   INITIAL_ENTRY_ID,
                                   NMR_CONTENT_SUBTYPES,
                                   MR_CONTENT_SUBTYPES,
                                   PK_CONTENT_SUBTYPES,
                                   CIF_CONTENT_SUBTYPES,
                                   READABLE_FILE_TYPE,
                                   CONTENT_TYPE,
                                   SF_CATEGORIES,
                                   LP_CATEGORIES,
                                   MR_MAX_SPACER_LINES,
                                   INDEX_TAGS,
                                   WEIGHT_TAGS,
                                   ANGLE_TYPE_TAGS,
                                   CONSIST_ID_TAGS,
                                   PK_KEY_ITEMS,
                                   DATA_ITEMS,
                                   NUM_DIM_ITEMS,
                                   ALLOWED_TAGS,
                                   DISALLOWED_PK_TAGS,
                                   ERR_TEMPLATE_FOR_MISSING_MANDATORY_LP_TAG,
                                   SF_TAG_PREFIXES,
                                   SF_TAG_ITEMS,
                                   MANDATORY_SF_TAG_ITEMS,
                                   SF_ALLOWED_TAGS,
                                   WARN_TEMPLATE_FOR_MISSING_MANDATORY_SF_TAG,
                                   INSTRUCTION_FOR_FULL_SEQUENCE,
                                   AUX_LP_CATEGORIES,
                                   LINKED_LP_CATEGORIES,
                                   AUX_ALLOWED_TAGS,
                                   ITEM_NAMES_IN_CS_LOOP,
                                   ITEM_NAMES_IN_PK_LOOP,
                                   ITEM_NAMES_IN_DIST_LOOP,
                                   ITEM_NAMES_IN_DIHED_LOOP,
                                   ITEM_NAMES_IN_RDC_LOOP,
                                   CS_LIST_SF_TAG_NAME,
                                   LOW_SEQ_COVERAGE,
                                   MIN_SEQ_COVERAGE_W_CONFLICT,
                                   LARGE_ASYM_ID,
                                   LEN_MAJOR_ASYM_ID,
                                   EMPTY_VALUE,
                                   TRUE_VALUE,
                                   STD_MON_DICT,
                                   PROTON_BEGIN_CODE,
                                   PARAMAGNETIC_ELEMENTS,
                                   FERROMAGNETIC_ELEMENTS,
                                   MAX_DIM_NUM_OF_SPECTRA,
                                   ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS,
                                   MAX_CONFLICT_ATTEMPT,
                                   PDB_ID_PAT,
                                   BMRB_ID_PAT,
                                   CNV_ID_PAT,
                                   WORK_MODEL_FILE_NAME_PAT,
                                   CHK_DESC_PAT,
                                   CHK_DESC_ONE_PAT,
                                   CHK_DESC_MAND_PAT,
                                   CHK_DESC_MAND_ONE_PAT,
                                   ARCHIVAL_MR_FILE_TYPES,
                                   PARSABLE_PK_FILE_TYPES,
                                   CS_RANGE_MIN,
                                   CS_RANGE_MAX,
                                   CS_UNCERT_MAX,
                                   REPRESENTATIVE_MODEL_ID,
                                   REPRESENTATIVE_ASYM_ID,
                                   REPRESENTATIVE_ALT_ID,
                                   SPECTRAL_DIM_TEMPLATE,
                                   DEFAULT_COORD_PROPERTIES)
    from nmr.NmrDpRegistry import (NmrDpRegistry,
                                   get_next_path,
                                   test_path_with_suffix)
    from nmr.NmrDpFirstAid import NmrDpFirstAid
    from nmr.NmrDpMrSplitter import (NmrDpMrSplitter,
                                     detect_bom,
                                     convert_codec,
                                     is_binary_file,
                                     get_type_of_star_file,
                                     get_peak_list_format)
    from nmr.NmrDpValidation import (NmrDpValidation,
                                     predict_redox_state_of_cysteine,
                                     is_like_planality_boundary)
    from nmr.NmrDpRemediation import (NmrDpRemediation,
                                      get_chem_shift_format_from_string)
    from nmr.NmrDpReport import (NmrDpReport,
                                 NmrDpReportInputSource)
    from nmr.ChemCompUtil import ChemCompUtil
    from nmr.BmrbChemShiftStat import BmrbChemShiftStat
    from nmr.AlignUtil import (hasLargeInnerSeqGap,
                               hasLargeSeqGap,
                               fillInnerBlankCompId,
                               fillBlankCompId,
                               fillBlankCompIdWithOffset,
                               beautifyPolySeq,
                               getMiddleCode,
                               getGaugeCode,
                               getScoreOfSeqAlign,
                               getOneLetterCodeCan,
                               getOneLetterCodeCanSequence,
                               alignPolymerSequence,
                               alignPolymerSequenceWithConflicts,
                               assignPolymerSequence,
                               retrieveAtomNameMappingFromRevisions,
                               retrieveAtomNameMappingFromInternal)
    from nmr.CifToNmrStar import (CifToNmrStar,
                                  has_key_value,
                                  get_first_sf_tag,
                                  set_sf_tag)
    from nmr.NmrVrptUtility import (uncompress_gzip_file,
                                    load_from_pickle,
                                    write_as_pickle,
                                    get_temp_path)
    from nmr.nef.NefTranslator import NefTranslator
    from nmr.io.CifReader import CifReader
    from nmr.io.PdbxUtil import abandon_symbolic_labels
    from nmr.mr.ParserListenerUtil import (coordAssemblyChecker,
                                           isAmbigAtomSelection,
                                           getTypeOfDihedralRestraint,
                                           isLikeHis)
    from nmr.pk.BasePKParserListener import guess_primary_dim_transfer_type
    from nmr.ann.OneDepAnnTasks import OneDepAnnTasks
    from nmr.ann.BmrbAnnTasks import BmrbAnnTasks


class NmrDpUtility:
    """ Main class for NMR data processing.
    """
    __slots__ = ('__class_name__',
                 '__version__',
                 '__reg',
                 '__procTasksDict',
                 '__alt_chain',
                 '__valid_seq',
                 '__remediation_loop_count',
                 '__ensemble_composition',
                 '__dstPath__',
                 '__logPath',
                 '__tmpPath',
                 '__cifHashCode',
                 '__inputParamDict__',
                 '__authSeqMap',
                 '__nmrIfR')

    def __init__(self, verbose: bool = False, log: IO = sys.stderr) -> None:
        self.__class_name__ = self.__class__.__name__
        self.__version__ = __version__

        self.__reg = NmrDpRegistry()

        self.__reg.verbose = verbose
        self.__reg.log = log

        self.__reg.ccU = ChemCompUtil(verbose, log)

        self.__reg.csStat = BmrbChemShiftStat(verbose, log, self.__reg.ccU)

        self.__reg.c2S = CifToNmrStar(log)

        self.__reg.nefT = NefTranslator(verbose, log, self.__reg.ccU, self.__reg.csStat, self.__reg.c2S)
        self.__reg.nefT.permit_missing_dist_restraint(self.__reg.permit_missing_legacy_dist_restraint)

        self.__reg.pA = PairwiseAlign()
        self.__reg.pA.setVerbose(verbose)

        self.__reg.cR = CifReader(verbose, log, use_cache=True)

        mr_content_subtypes = MR_CONTENT_SUBTYPES
        nmr_rep_content_subtypes = ['chem_shift', 'spectral_peak']
        nmr_rep_content_subtypes.extend(mr_content_subtypes)

        def sf_key(content_subtype):
            return self.__reg.c2S.category_order.index(SF_TAG_PREFIXES['nmr-star'][content_subtype])

        mr_content_subtypes.sort(key=sf_key)
        nmr_rep_content_subtypes.sort(key=sf_key)

        self.__reg.mr_content_subtypes = mr_content_subtypes
        self.__reg.nmr_rep_content_subtypes = nmr_rep_content_subtypes

        self.__reg.dpS = NmrDpMrSplitter(self.__reg)
        self.__reg.dpA = NmrDpFirstAid(self.__reg)
        self.__reg.dpV = NmrDpValidation(self.__reg)
        self.__reg.dpR = NmrDpRemediation(self.__reg)

        # validation tasks for NMR data only
        __nmrCheckTasks = [self.__detectContentSubType,
                           self.__extractPublicMrFileIntoLegacyMr,
                           self.__detectContentSubTypeOfLegacyPk,
                           self.__detectContentSubTypeOfLegacyMr,
                           self.__extractPolymerSequence,
                           self.__extractPolymerSequenceInLoop,
                           self.__extractCommonPolymerSequence,
                           self.__extractNonStandardResidue,
                           self.__appendPolymerSequenceAlignment,
                           self.__testSequenceConsistency,
                           self.__validateAtomNomenclature,
                           self.__appendElemAndIsoNumOfNefCsLoop,
                           self.__validateAtomTypeOfCsLoop,
                           self.__validateAmbigCodeOfCsLoop,
                           self.__detectConflictDataInLoop,
                           self.__appendIndexTag,
                           self.__testIndexConsistency,
                           self.__appendWeightInLoop,
                           self.__appendDihedAngleType,
                           self.__remediateRdcLoop,
                           self.__testDataConsistencyInLoop,
                           self.__testDataConsistencyInAuxLoop,
                           self.__testNmrCovalentBond,
                           self.__appendSfTagItem,
                           self.__testSfTagConsistency,
                           self.__testCsPseudoAtomNameConsistencyInMrLoop,
                           self.__testCsValueConsistencyInPkLoop,
                           self.__testCsValueConsistencyInPkAltLoop
                           ]

        # validation tasks for coordinate file only
        __cifCheckTasks = [self.__validateCoordInputSource,
                           self.__detectCoordContentSubType,
                           self.__extractCoordPolymerSequence,
                           self.__extractCoordPolymerSequenceInLoop,
                           self.__extractCoordCommonPolymerSequence,
                           self.__extractCoordNonStandardResidue,
                           self.__appendCoordPolymerSequenceAlignment,
                           self.__testTautomerOfHistidinePerModel
                           ]

        # cross validation tasks
        __crossCheckTasks = [self.__assignCoordPolymerSequence,
                             self.__testCoordAtomIdConsistency,
                             self.__testCoordCovalentBond,
                             self.__testResidueVariant,
                             self.__validateCsValue,
                             self.__testRdcVector,
                             self.__extractCoordDisulfideBond,
                             self.__extractCoordOtherBond,
                             self.__validateStrMr,
                             self.__validateLegacyMr,
                             self.__validateLegacyPk,
                             self.__validateLegacyCs,
                             self.__validateLegacyCsp,
                             self.__validateSaxsMr,
                             self.__validateStrPk,
                             self.__updateConstraintStats,
                             self.__detectSimpleDistanceRestraint,
                             self.__calculateStatsOfExptlData,
                             self.__detectDimTransferTypeViaThroughSpace
                             ]

        # nmr-*-consistency-check tasks
        __checkTasks = [self.__initializeDpReport,
                        self.__validateInputSource
                        ]
        __checkTasks.extend(__nmrCheckTasks)
        __checkTasks.extend(__cifCheckTasks)
        __checkTasks.extend(__crossCheckTasks)

        # nmr-*-deposit tasks
        __depositTasks = [self.__retrieveDpReport,
                          self.__validateInputSource,
                          self.__parseCoordinate,
                          self.__detectCoordContentSubType,
                          self.__extractCoordPolymerSequence,
                          # resolve conflict
                          self.__resolveConflictsInLoop,
                          self.__resolveConflictsInAuxLoop,
                          # resolve minor issues
                          self.__validateAtomNomenclature,
                          self.__appendIndexTag,
                          self.__appendWeightInLoop,
                          self.__appendDihedAngleType,
                          self.__appendSfTagItem,
                          self.__deleteSkippedSf,
                          self.__deleteSkippedLoop,
                          self.__deleteUnparsedEntryLoop,
                          self.__updatePolymerSequence,
                          self.__remediateRawTextPk,
                          self.__updateAuthSequence,
                          self.__updateDihedralAngleType,
                          self.__fixDisorderedIndex,
                          self.__removeNonSenseZeroValue,
                          self.__fixNonSenseNegativeValue,
                          self.__fixEnumMismatch,
                          self.__fixEnumMismatchIgnorable,
                          self.__resetCapitalStringInLoop,
                          self.__resetBoolValueInLoop,
                          self.__resetBoolValueInAuxLoop,
                          self.__appendParentSfTag,
                          self.__addUnnamedEntryId,
                          self.__removeUnusedPdbInsCode,
                          self.__depositNmrData,
                          # re-setup for next
                          self.__initializeDpReportForNext,
                          self.__validateInputSourceForNext
                          ]

        __depositTasks.extend(__nmrCheckTasks)
        __depositTasks.extend(__cifCheckTasks)
        __depositTasks.extend(__crossCheckTasks)

        # additional nmr-nef2str/nef2cif tasks
        __nef2strTasks = [self.__translateNef2Str,
                          self.__dumpDpReport,
                          self.__initResourceForNef2Str
                          ]

        __nef2strTasks.extend(__checkTasks)
        __nef2strTasks.append(self.__dumpDpReport)
        __nef2strTasks.extend(__depositTasks)

        # additional nmr-str2nef tasks
        __str2nefTasks = [self.__translateStr2Nef,
                          self.__dumpDpReport,
                          self.__initResourceForStr2Nef
                          ]

        __str2nefTasks.extend(__checkTasks)
        __str2nefTasks.append(self.__dumpDpReport)
        __str2nefTasks.extend(__depositTasks)

        __mergeCsAndMrTasks = __checkTasks
        __mergeCsAndMrTasks.remove(self.__detectSimpleDistanceRestraint)
        __mergeCsAndMrTasks.remove(self.__calculateStatsOfExptlData)
        __mergeCsAndMrTasks.remove(self.__detectDimTransferTypeViaThroughSpace)
        __mergeCsAndMrTasks.append(self.__updatePolymerSequence)
        __mergeCsAndMrTasks.append(self.__mergeLegacyData)
        __mergeCsAndMrTasks.append(self.__detectSimpleDistanceRestraint)
        __mergeCsAndMrTasks.append(self.__remediateRawTextPk)
        __mergeCsAndMrTasks.append(self.__mergeCoordAsNmrIf)  # DAOTHER-8905: NMR data remediation Phase 2 (internal remediation)
        __mergeCsAndMrTasks.append(self.__performBmrbAnnTasks)
        __mergeCsAndMrTasks.append(self.__testDataConsistencyInPkLoop)  # refresh statistics of spectral peak list
        __mergeCsAndMrTasks.append(self.__testDataConsistencyInPkAuxLoop)  # refresh statistics of spectral peak list
        __mergeCsAndMrTasks.append(self.__calculateStatsOfExptlData)
        __mergeCsAndMrTasks.append(self.__detectDimTransferTypeViaThroughSpace)
        __mergeCsAndMrTasks.append(self.__discardPeakListRemediation)  # OneDep only

        __annotateTasks = [self.__initializeDpReport,
                           self.__validateInputSource,
                           self.__detectContentSubType,
                           self.__extractPolymerSequence,
                           self.__extractPolymerSequenceInLoop,
                           self.__extractCommonPolymerSequence,
                           self.__extractNonStandardResidue,
                           self.__appendPolymerSequenceAlignment
                           ]

        __annotateTasks.extend(__cifCheckTasks)
        __annotateTasks.extend(__crossCheckTasks)
        __annotateTasks.append(self.__updatePolymerSequence)
        __annotateTasks.append(self.__remediateRawTextPk)
        __annotateTasks.append(self.__performBmrbAnnTasks)
        __annotateTasks.append(self.__depositNmrData)
        __annotateTasks.extend(__depositTasks)
        __annotateTasks.append(self.__depositNmrData)

        __mergeNmrIfTasks = [self.__parseNmrIf,
                             self.__mergeNmrIf,
                             self.__performBmrbAnnTasks
                             ]

        __replaceCsTasks = copy.copy(__checkTasks)
        __replaceCsTasks.append(self.__replaceCsSf)
        __replaceCsTasks.append(self.__updatePolymerSequence)
        __replaceCsTasks.append(self.__performBmrbAnnTasks)
        __replaceCsTasks.append(self.__depositNmrData)

        # dictionary of processing tasks of each workflow operation
        self.__procTasksDict = {'consistency-check': __checkTasks,
                                'deposit': __depositTasks,
                                'nmr-nef2str-deposit': __nef2strTasks,
                                'nmr-nef2cif-deposit': __nef2strTasks,
                                'nmr-str2nef-release': __str2nefTasks,
                                'nmr-cs-nef-consistency-check': [self.__depositLegacyNmrData],
                                'nmr-cs-str-consistency-check': [self.__depositLegacyNmrData],
                                'nmr-cs-mr-merge': __mergeCsAndMrTasks,
                                'nmr-str2cif-annotate': __annotateTasks,
                                'nmr-if-merge-deposit': __mergeNmrIfTasks,
                                'nmr-str-replace-cs': __replaceCsTasks
                                }

        # internal statuses
        self.__alt_chain: bool = False
        self.__valid_seq: bool = False

        # loop count of remediation
        self.__remediation_loop_count: int = 0

        # ensemble composition
        self.__ensemble_composition: dict = None

        # copy of dstPath
        self.__dstPath__: str = None

        # log file path
        self.__logPath: str = None

        # temporary file path to be removed (release mode)
        self.__tmpPath: str = None

        # hash code of the coordinate file
        self.__cifHashCode = None

        # copy of inputParamDict to restart remediation
        self.__inputParamDict__: dict = None

        # temporary dictionaries used in mapping auth sequence scheme
        self.__authSeqMap = {}

        # NMRIF reader
        self.__nmrIfR = None

    def getNextPath(self, src_path: str, suffix: str = '~') -> str:
        """ Return candidate next file path.
        """
        return get_next_path(self.__reg, src_path, suffix)

    def testPathWithSuffix(self, src_path: str, suffix: str, defer_check: bool = False) -> str:
        """ Return basename(src_path) + suffix file path in either current workspace or default workspace if possible.
        """
        return test_path_with_suffix(self.__reg, src_path, suffix, defer_check)

    def setVerbose(self, verbose: bool) -> None:
        """ Set verbose mode.
        """

        self.__reg.verbose = verbose
        self.__reg.debug = verbose

    def setMrDebugMode(self, debug: bool) -> None:
        """ Set debug mode for MR splitter.
        """

        self.__reg.mr_debug = debug

    def setSource(self, fPath: str, originalName: str = None) -> None:
        """ Set primary source file path.
        """

        if os.access(fPath, os.F_OK):
            self.__reg.srcPath = os.path.abspath(fPath)
            if originalName is not None:
                self.__reg.srcName = originalName
            else:
                self.__reg.srcName = os.path.basename(self.__reg.srcPath)

        else:
            raise IOError(f"+{self.__class_name__}.setSource() ++ Error  - Could not access to file path {fPath}.")

    def setDestination(self, fPath: str) -> None:
        """ Set primary destination file path.
        """

        if fPath is not None:
            self.__reg.dstPath = os.path.abspath(fPath)
            self.__dstPath__ = copy.copy(self.__reg.dstPath)

    def setLog(self, fPath: str) -> None:
        """ Set a log file path for the primary input source.
        """

        if fPath is not None:
            self.__logPath = os.path.abspath(fPath)

    def setWorkspace(self, dirPath: str, cacheDirPath: str = None) -> None:
        """ Set working directory and cache directory (optional).
        """

        if dirPath is not None:
            if not os.path.isdir(dirPath):
                os.makedirs(dirPath)
            self.__reg.dirPath = dirPath

            # share current working directory with Cif2NmrStar
            self.__reg.c2S.dirPath = dirPath

        if cacheDirPath is not None:
            if not os.path.isdir(cacheDirPath):
                os.makedirs(cacheDirPath)
            self.__reg.cacheDirPath = cacheDirPath

            # share the cache directory with CifReader
            self.__reg.cR.cacheDirPath = cacheDirPath

    def addInput(self, name: Optional[str] = None, value: Any = None, type: str = 'file') -> None:  # noqa: E501, pylint: disable=redefined-builtin,line-too-long
        """ Add a named input and value to the dictionary of input parameters.
        """

        try:

            if type == 'param':
                if name not in DP_INPUT_PARAM_KEYS:
                    raise KeyError(f"+{self.__class_name__}.addInput() ++ Error  - Unknown input param {name!r}.")
                self.__reg.inputParamDict[name] = value
            elif type == 'file':
                if name not in DP_INPUT_FILE_KEYS:
                    raise KeyError(f"+{self.__class_name__}.addInput() ++ Error  - Unknown input file {name!r}.")
                self.__reg.inputParamDict[name] = os.path.abspath(value)
            elif type == 'file_list':
                if name not in DP_INPUT_FILE_LIST_KEYS:
                    raise KeyError(f"+{self.__class_name__}.addInput() ++ Error  - Unknown input file_list {name!r}.")
                self.__reg.inputParamDict[name] = [os.path.abspath(f) for f in value]
            elif type == 'file_dict_list':
                if name not in DP_INPUT_FILE_DICT_KEYS:
                    raise KeyError(f"+{self.__class_name__}.addInput() ++ Error  - Unknown input file_dict_list {name!r}.")
                if any(True for f in value if 'original_file_name' in f):
                    self.__reg.inputParamDict[name] = []
                    for f in value:
                        if 'original_file_name' in f:
                            self.__reg.inputParamDict[name].append({'file_name': os.path.abspath(f['file_name']),
                                                                    'file_type': f['file_type'],
                                                                    'original_file_name': f['original_file_name'],
                                                                    'ignore_error':
                                                                    False if 'ignore_error' not in f else f['ignore_error']})
                        else:
                            self.__reg.inputParamDict[name].append({'file_name': os.path.abspath(f['file_name']),
                                                                    'file_type': f['file_type'],
                                                                    'ignore_error':
                                                                    False if 'ignore_error' not in f else f['ignore_error']})
                else:
                    self.__reg.inputParamDict[name] = [{'file_name': os.path.abspath(f['file_name']),
                                                        'file_type': f['file_type'],
                                                        'ignore_error': False if 'ignore_error' not in f else f['ignore_error']}
                                                       for f in value]
            else:
                raise KeyError(f"+{self.__class_name__}.addInput() ++ Error  - Unknown input type {type!r}.")

        except Exception as e:  # pylint: disable=broad-exception-caught
            raise ValueError(f"+{self.__class_name__}.addInput() ++ Error  - " + str(e)) from e

    def addOutput(self, name: Optional[str] = None, value: Any = None, type: str = 'file') -> None:  # noqa: E501, pylint: disable=redefined-builtin,line-too-long
        """ Add a named input and value to the dictionary of output parameters.
        """

        try:

            if type == 'param':
                if name not in DP_OUTPUT_PARAM_KEYS:
                    raise KeyError(f"+{self.__class_name__}.addOutput() ++ Error  - Unknown output param {name!r}.")
                self.__reg.outputParamDict[name] = value
            elif type == 'file':
                if name not in DP_OUTPUT_FILE_KEYS:
                    raise KeyError(f"+{self.__class_name__}.addOutput() ++ Error  - Unknown output file {name!r}.")
                self.__reg.outputParamDict[name] = os.path.abspath(value)
            elif type == 'file_list':
                if name not in DP_OUTPUT_FILE_LIST_KEYS:
                    raise KeyError(f"+{self.__class_name__}.addOutput() ++ Error  - Unknown output file_list {name!r}.")
                self.__reg.outputParamDict[name] = [os.path.abspath(f) for f in value]
            else:
                raise KeyError(f"+{self.__class_name__}.addOutput() ++ Error  - Unknown output type {type!r}.")

        except Exception as e:  # pylint: disable=broad-exception-caught
            raise ValueError(f"+{self.__class_name__}.addOutput() ++ Error  - " + str(e)) from e

    def release(self) -> None:
        """ Release the memory held for the coordinate file between workflow operations, e.g. before running
            an NmrVrptUtility in the same process that is not given this utility's CifReader, and so parses the
            coordinate file on its own. The parsed coordinates otherwise stay in memory so that the next op() on
            the same file skips re-parsing it; after release(), the next op() re-reads it. On a 134 MB coordinate
            file (a single model of 1.36 M atoms) this frees 1.3 GB, and the peak of the following restraint
            validation drops from 3.5 GB to 2.2 GB (DAOTHER-7829, 9785).
        """

        self.__reg.cR.release()
        self.__reg.caC = None  # rebuilt by every op(), see __parseCoordFilePath()

        if not self.__reg.nefT.release():
            gc.collect(2)

    def op(self, op: str) -> bool:
        """ Perform a series of tasks for a given workflow operation.
        """

        if op not in DP_WORKFLOW_OPS:
            raise KeyError(f"+{self.__class_name__}.op() ++ Error  - Unknown workflow operation {op!r}.")

        self.__reg.rescue_mode = True

        if has_key_value(self.__reg.inputParamDict, 'validation_server'):
            if isinstance(self.__reg.inputParamDict['validation_server'], bool):
                self.__reg.validation_server = self.__reg.inputParamDict['validation_server']
            else:
                self.__reg.validation_server = self.__reg.inputParamDict['validation_server'] in TRUE_VALUE

        if has_key_value(self.__reg.inputParamDict, 'conversion_server'):
            if isinstance(self.__reg.inputParamDict['conversion_server'], bool):
                self.__reg.conversion_server = self.__reg.inputParamDict['conversion_server']
            else:
                self.__reg.conversion_server = self.__reg.inputParamDict['conversion_server'] in TRUE_VALUE

            if self.__reg.conversion_server:
                self.__reg.nefT.permit_missing_chem_shift(True)
                self.__reg.bmrb_only = True  # self.__reg.internal_mode = True

        if has_key_value(self.__reg.inputParamDict, 'bmrb_only'):
            if isinstance(self.__reg.inputParamDict['bmrb_only'], bool):
                self.__reg.bmrb_only = self.__reg.inputParamDict['bmrb_only']
            else:
                self.__reg.bmrb_only = self.__reg.inputParamDict['bmrb_only'] in TRUE_VALUE

        if has_key_value(self.__reg.inputParamDict, 'permit_missing_cs_file'):
            if isinstance(self.__reg.inputParamDict['permit_missing_cs_file'], bool):
                self.__reg.permit_missing_cs_file = self.__reg.inputParamDict['permit_missing_cs_file']
            else:
                self.__reg.permit_missing_cs_file = self.__reg.inputParamDict['permit_missing_cs_file'] in TRUE_VALUE

        if self.__reg.permit_missing_cs_file:
            self.__reg.nefT.permit_missing_chem_shift(True)
            self.__reg.bmrb_only = True

        self.__reg.combined_mode = 'cs' not in op or op == 'nmr-str-replace-cs'\
            or self.__reg.conversion_server and self.__reg.srcPath is not None

        if self.__reg.combined_mode:
            if self.__reg.srcPath is None:
                raise ValueError(f"+{self.__class_name__}.op() ++ Error  - No input provided for workflow operation {op}.")

            self.__reg.cs_file_path_list_len = 0
            self.__reg.file_path_list_len = 1

            if op == 'nmr-str-replace-cs':
                self.__reg.cs_file_path_list_len = len(self.__reg.inputParamDict[CS_FILE_PATH_LIST_KEY])
                self.__reg.file_path_list_len += self.__reg.cs_file_path_list_len

        else:
            if self.__reg.srcPath is not None and self.__reg.bmrb_only:
                self.__reg.combined_mode = True

                if CS_FILE_PATH_LIST_KEY in self.__reg.inputParamDict:
                    self.__reg.cs_file_path_list_len = len(self.__reg.inputParamDict[CS_FILE_PATH_LIST_KEY])
                    self.__reg.file_path_list_len += self.__reg.cs_file_path_list_len
                else:
                    self.__reg.cs_file_path_list_len = 0
                    self.__reg.file_path_list_len = 1

            if CS_FILE_PATH_LIST_KEY not in self.__reg.inputParamDict:
                if not self.__reg.permit_missing_cs_file:
                    raise ValueError(f"+{self.__class_name__}.op() ++ Error  - No input provided for workflow operation {op}.")

            else:
                self.__reg.cs_file_path_list_len = len(self.__reg.inputParamDict[CS_FILE_PATH_LIST_KEY])
                self.__reg.file_path_list_len = self.__reg.cs_file_path_list_len

            if MR_FILE_PATH_LIST_KEY in self.__reg.inputParamDict:
                self.__reg.file_path_list_len += len(self.__reg.inputParamDict[MR_FILE_PATH_LIST_KEY])

        self.__reg.cifPath = self.__cifHashCode = None
        self.__reg.cifChecked = False

        # incomplete assignments are edited by biocurators for conventional assigned chemical shifts (DAOTHER-7662)
        for key in self.__reg.key_items['nmr-star']['chem_shift']:
            if 'remove-bad-pattern' in key:
                key['remove-bad-pattern'] = self.__reg.combined_mode

        if has_key_value(self.__reg.inputParamDict, 'remediation'):
            if isinstance(self.__reg.inputParamDict['remediation'], bool):
                self.__reg.remediation_mode = self.__reg.inputParamDict['remediation']
            else:
                self.__reg.remediation_mode = self.__reg.inputParamDict['remediation'] in TRUE_VALUE

        if has_key_value(self.__reg.inputParamDict, 'internal'):
            if isinstance(self.__reg.inputParamDict['internal'], bool):
                self.__reg.internal_mode = self.__reg.inputParamDict['internal']
            else:
                self.__reg.internal_mode = self.__reg.inputParamDict['internal'] in TRUE_VALUE

        if op == 'nmr-cs-mr-merge':

            self.__reg.remediation_mode = True
            self.__reg.has_star_chem_shift = True

            if self.__inputParamDict__ is None:
                self.__inputParamDict__ = copy.deepcopy(self.__reg.inputParamDict)

            for v in self.__reg.key_items['nmr-star'].values():
                if v is None:
                    continue
                for d in v:
                    if d['name'].startswith('Entity_assembly_ID'):
                        d['type'] = 'str'
                        d['default'] = REPRESENTATIVE_ASYM_ID
                        if 'default-from' in d:
                            del d['default-from']

            for v in self.__reg.consist_key_items['nmr-star'].values():
                if v is None:
                    continue
                for d in v:
                    if d['name'].startswith('Entity_assembly_ID'):
                        d['type'] = 'str'
                        d['default'] = REPRESENTATIVE_ASYM_ID
                        if 'default-from' in d:
                            del d['default-from']

            for d in self.__reg.pk_data_items['nmr-star']:
                if d['name'].startswith('Entity_assembly_ID'):
                    d['type'] = 'str'
                    d['default'] = REPRESENTATIVE_ASYM_ID
                    if 'default-from' in d:
                        del d['default-from']
                    if 'enforce-non-zero' in d:
                        del d['enforce-non-zero']

            for v in self.__reg.aux_key_items['nmr-star'].values():
                if v is None:
                    continue
                for v2 in v.values():
                    for d in v2:
                        if d['name'].startswith('Entity_assembly_ID'):
                            d['type'] = 'str'
                            d['default'] = REPRESENTATIVE_ASYM_ID
                            if 'default-from' in d:
                                del d['default-from']

            for v in self.__reg.aux_data_items['nmr-star'].values():
                if v is None:
                    continue
                for v2 in v.values():
                    for d in v2:
                        if d['name'].startswith('Entity_assembly_ID'):
                            d['type'] = 'str'
                            d['default'] = REPRESENTATIVE_ASYM_ID
                            if 'default-from' in d:
                                del d['default-from']

        elif self.__reg.combined_mode and not self.__reg.remediation_mode:
            self.__reg.native_combined = True

        self.__remediation_loop_count = 0

        self.__reg.sll_pred_holder.clear()
        self.__reg.sll_pred_forced.clear()

        self.__reg.submission_mode = 'merge-deposit' in op
        self.__reg.annotation_mode = 'annotate' in op
        self.__reg.release_mode = 'release' in op or 'replace-cs' in op

        self.__reg.nefT.set_remediation_mode(self.__reg.remediation_mode)
        self.__reg.nefT.set_annotation_mode(self.__reg.annotation_mode)
        self.__reg.nefT.set_internal_mode(self.__reg.internal_mode or self.__reg.conversion_server)
        self.__reg.nefT.set_merge_rescue_mode(op in ('nmr-cs-mr-merge', 'nmr-str-replace-cs')  # DAOTHER-9927
                                              or (op == 'nmr-str2cif-annotate'
                                                  and self.__reg.remediation_mode))  # DAOTHER-10616
        self.__reg.nefT.cache_clear()

        if not self.__reg.permit_missing_legacy_dist_restraint and self.__reg.remediation_mode:
            self.__reg.nefT.permit_missing_dist_restraint(True)
            self.__reg.permit_missing_dist_restraint = self.__reg.permit_missing_legacy_dist_restraint = True

        if self.__reg.verbose:
            self.__reg.log.write(f"+{self.__class_name__}.op() starting op {op}\n")

        if 'cif' in op:
            if NMR_CIF_FILE_PATH_KEY not in self.__reg.outputParamDict:
                raise KeyError(f"+{self.__class_name__}.op() "
                               f"++ Error  - Could not find {NMR_CIF_FILE_PATH_KEY!r} output parameter.")
            if self.__reg.dstPath is None:
                self.__reg.dstPath = get_temp_path(self.__reg.outputParamDict[NMR_CIF_FILE_PATH_KEY])
                self.__dstPath__ = copy.copy(self.__reg.dstPath)
                self.__tmpPath = self.__dstPath__

        if self.__reg.release_mode and self.__reg.dstPath is None:
            self.__reg.dstPath = get_temp_path(self.__reg.srcPath)
            try:
                shutil.copyfile(self.__reg.srcPath, self.__reg.dstPath)
                self.__dstPath__ = copy.copy(self.__reg.dstPath)
                self.__tmpPath = self.__dstPath__
            except Exception:  # pylint: disable=broad-exception-caught
                self.__reg.dstPath = self.__reg.srcPath

        if has_key_value(self.__reg.outputParamDict, 'entry_id'):
            # DAOTHER-9511: replace white space in a datablock name to underscore
            self.__reg.entry_id = self.__reg.outputParamDict['entry_id'].strip().replace(' ', '_')

        if self.__reg.bmrb_only:
            self.__reg.cs_anomalous_error_scaled_by_sigma = 4.0
            self.__reg.cs_unusual_error_scaled_by_sigma = 3.5
            self.__reg.cs_diff_error_scaled_by_sigma = 5.0
            self.__reg.nefT.set_bmrb_only_mode(True)
            self.__reg.resolve_conflict = True

            if self.__reg.conversion_server:
                if self.__reg.entry_id is not None and CNV_ID_PAT.match(self.__reg.entry_id):
                    self.__reg.bmrb_id = self.__reg.entry_id
                # mitigation: in case bmrb_id has valid conversion_id, instead of entry_id
                elif has_key_value(self.__reg.inputParamDict, 'bmrb_id'):
                    if isinstance(self.__reg.inputParamDict['bmrb_id'], int):
                        if CNV_ID_PAT.match(f"C_{self.__reg.inputParamDict['bmrb_id']}"):
                            self.__reg.entry_id = self.__reg.bmrb_id = f"C_{self.__reg.inputParamDict['bmrb_id']}"
                    elif isinstance(self.__reg.inputParamDict['bmrb_id'], str):
                        if CNV_ID_PAT.match(self.__reg.inputParamDict['bmrb_id']):
                            self.__reg.entry_id = self.__reg.bmrb_id = self.__reg.inputParamDict['bmrb_id']

            elif has_key_value(self.__reg.inputParamDict, 'bmrb_id'):
                if isinstance(self.__reg.inputParamDict['bmrb_id'], int):
                    self.__reg.bmrb_id = str(self.__reg.inputParamDict['bmrb_id'])
                elif isinstance(self.__reg.inputParamDict['bmrb_id'], str):
                    self.__reg.bmrb_id = self.__reg.inputParamDict['bmrb_id']
                if self.__reg.bmrb_id is not None:
                    if BMRB_ID_PAT.match(self.__reg.bmrb_id):
                        if self.__reg.bmrb_id.startswith('bmr'):
                            self.__reg.bmrb_id = self.__reg.bmrb_id[3:]
                    else:
                        self.__reg.bmrb_id = None

                if self.__reg.bmrb_id is not None:
                    # DAOTHER-9511: replace white space in a datablock name to underscore
                    self.__reg.entry_id = self.__reg.bmrb_id.strip().replace(' ', '_')

        self.__reg.assembly_name = '?'

        entity_name_item = next(item for item in SF_TAG_ITEMS['nmr-star']['entity'] if item['name'] == 'Name')
        entity_name_item['mandatory'] = self.__reg.bmrb_only

        if has_key_value(self.__reg.inputParamDict, 'merge_any_pk_as_is'):
            if isinstance(self.__reg.inputParamDict['merge_any_pk_as_is'], bool):
                self.__reg.merge_any_pk_as_is = self.__reg.inputParamDict['merge_any_pk_as_is']
            else:
                self.__reg.merge_any_pk_as_is = self.__reg.inputParamDict['merge_any_pk_as_is'] in TRUE_VALUE

        if has_key_value(self.__reg.inputParamDict, 'enforce_peak_row_format'):
            if isinstance(self.__reg.inputParamDict['enforce_peak_row_format'], bool):
                self.__reg.enforce_peak_row_format = self.__reg.inputParamDict['enforce_peak_row_format']
            else:
                self.__reg.enforce_peak_row_format = self.__reg.inputParamDict['enforce_peak_row_format'] in TRUE_VALUE

        if has_key_value(self.__reg.inputParamDict, 'nonblk_anomalous_cs'):
            if isinstance(self.__reg.inputParamDict['nonblk_anomalous_cs'], bool):
                self.__reg.nonblk_anomalous_cs = self.__reg.inputParamDict['nonblk_anomalous_cs']
            else:
                self.__reg.nonblk_anomalous_cs = self.__reg.inputParamDict['nonblk_anomalous_cs'] in TRUE_VALUE

        if has_key_value(self.__reg.inputParamDict, 'nonblk_bad_nterm'):
            if isinstance(self.__reg.inputParamDict['nonblk_bad_nterm'], bool):
                self.__reg.nonblk_bad_nterm = self.__reg.inputParamDict['nonblk_bad_nterm']
            else:
                self.__reg.nonblk_bad_nterm = self.__reg.inputParamDict['nonblk_bad_nterm'] in TRUE_VALUE

        if has_key_value(self.__reg.inputParamDict, 'update_poly_seq'):
            if isinstance(self.__reg.inputParamDict['update_poly_seq'], bool):
                self.__reg.update_poly_seq = self.__reg.inputParamDict['update_poly_seq']
            else:
                self.__reg.update_poly_seq = self.__reg.inputParamDict['update_poly_seq'] in TRUE_VALUE

        if has_key_value(self.__reg.inputParamDict, 'resolve_conflict'):
            if isinstance(self.__reg.inputParamDict['resolve_conflict'], bool):
                self.__reg.resolve_conflict = self.__reg.inputParamDict['resolve_conflict']
            else:
                self.__reg.resolve_conflict = self.__reg.inputParamDict['resolve_conflict'] in TRUE_VALUE

        if has_key_value(self.__reg.inputParamDict, 'check_mandatory_tag'):
            if isinstance(self.__reg.inputParamDict['check_mandatory_tag'], bool):
                self.__reg.check_mandatory_tag = self.__reg.inputParamDict['check_mandatory_tag']
            else:
                self.__reg.check_mandatory_tag = self.__reg.inputParamDict['check_mandatory_tag'] in TRUE_VALUE

        if has_key_value(self.__reg.inputParamDict, 'check_auth_seq'):
            if isinstance(self.__reg.inputParamDict['check_auth_seq'], bool):
                self.__reg.check_auth_seq = self.__reg.inputParamDict['check_auth_seq']
            else:
                self.__reg.check_auth_seq = self.__reg.inputParamDict['check_auth_seq'] in TRUE_VALUE

        if has_key_value(self.__reg.inputParamDict, 'transl_pseudo_name'):
            if isinstance(self.__reg.inputParamDict['transl_pseudo_name'], bool):
                self.__reg.transl_pseudo_name = self.__reg.inputParamDict['transl_pseudo_name']
            else:
                self.__reg.transl_pseudo_name = self.__reg.inputParamDict['transl_pseudo_name'] in TRUE_VALUE
        elif op in ('nmr-str-consistency-check', 'nmr-str2str-deposit', 'nmr-str2cif-deposit',
                    'nmr-str2nef-release', 'nmr-str2cif-annotate'):
            self.__reg.transl_pseudo_name = True

        if has_key_value(self.__reg.inputParamDict, 'tolerant_seq_align'):
            if isinstance(self.__reg.inputParamDict['tolerant_seq_align'], bool):
                self.__reg.tolerant_seq_align = self.__reg.inputParamDict['tolerant_seq_align']
            else:
                self.__reg.tolerant_seq_align = self.__reg.inputParamDict['tolerant_seq_align'] in TRUE_VALUE

        if has_key_value(self.__reg.inputParamDict, 'fix_format_issue'):
            if isinstance(self.__reg.inputParamDict['fix_format_issue'], bool):
                self.__reg.fix_format_issue = self.__reg.inputParamDict['fix_format_issue']
            else:
                self.__reg.fix_format_issue = self.__reg.inputParamDict['fix_format_issue'] in TRUE_VALUE
        elif not self.__reg.combined_mode or self.__reg.release_mode:
            self.__reg.fix_format_issue = True

        if has_key_value(self.__reg.inputParamDict, 'excl_missing_data'):
            if isinstance(self.__reg.inputParamDict['excl_missing_data'], bool):
                self.__reg.excl_missing_data = self.__reg.inputParamDict['excl_missing_data']
            else:
                self.__reg.excl_missing_data = self.__reg.inputParamDict['excl_missing_data'] in TRUE_VALUE
        elif not self.__reg.combined_mode:
            self.__reg.excl_missing_data = True

        if has_key_value(self.__reg.inputParamDict, 'cmpl_missing_data'):
            if isinstance(self.__reg.inputParamDict['cmpl_missing_data'], bool):
                self.__reg.cmpl_missing_data = self.__reg.inputParamDict['cmpl_missing_data']
            else:
                self.__reg.cmpl_missing_data = self.__reg.inputParamDict['cmpl_missing_data'] in TRUE_VALUE
        elif not self.__reg.combined_mode:
            self.__reg.cmpl_missing_data = True

        if has_key_value(self.__reg.inputParamDict, 'trust_pdbx_nmr_ens'):
            if isinstance(self.__reg.inputParamDict['trust_pdbx_nmr_ens'], bool):
                self.__reg.trust_pdbx_nmr_ens = self.__reg.inputParamDict['trust_pdbx_nmr_ens']
            else:
                self.__reg.trust_pdbx_nmr_ens = self.__reg.inputParamDict['trust_pdbx_nmr_ens'] in TRUE_VALUE
        elif self.__reg.release_mode:
            self.__reg.trust_pdbx_nmr_ens = True

        if has_key_value(self.__reg.inputParamDict, 'rmsd_not_superimposed'):
            if isinstance(self.__reg.inputParamDict['rmsd_not_superimposed'], float):
                self.__reg.rmsd_not_superimposed = self.__reg.inputParamDict['rmsd_not_superimposed']

        if has_key_value(self.__reg.inputParamDict, 'rmsd_overlaid_exactly'):
            if isinstance(self.__reg.inputParamDict['rmsd_overlaid_exactly'], float):
                self.__reg.rmsd_overlaid_exactly = self.__reg.inputParamDict['rmsd_overlaid_exactly']

        if has_key_value(self.__reg.outputParamDict, 'retain_original'):
            if isinstance(self.__reg.outputParamDict['retain_original'], bool):
                self.__reg.retain_original = self.__reg.outputParamDict['retain_original']
            else:
                self.__reg.retain_original = self.__reg.outputParamDict['retain_original'] in TRUE_VALUE

        if has_key_value(self.__reg.outputParamDict, 'leave_intl_note'):
            if isinstance(self.__reg.outputParamDict['leave_intl_note'], bool):
                self.__reg.leave_intl_note = self.__reg.outputParamDict['leave_intl_note']
            else:
                self.__reg.leave_intl_note = self.__reg.outputParamDict['leave_intl_note'] in TRUE_VALUE

        if has_key_value(self.__reg.outputParamDict, 'reduced_atom_notation'):
            if isinstance(self.__reg.outputParamDict['reduced_atom_notation'], bool):
                self.__reg.reduced_atom_notation = self.__reg.outputParamDict['reduced_atom_notation']
            else:
                self.__reg.reduced_atom_notation = self.__reg.outputParamDict['reduced_atom_notation'] in TRUE_VALUE

        self.__reg.op = op

        if op.endswith('consistency-check'):

            for task in self.__procTasksDict['consistency-check']:

                if self.__reg.verbose:
                    self.__reg.log.write(f"+{self.__class_name__}.op() starting op {op} - task {task.__name__}\n")

                start_time = time.time()

                if not task():
                    pass

                if self.__reg.verbose:
                    end_time = time.time()
                    if end_time - start_time > 1.0:
                        self.__reg.log.write(f"op: {op}, task: {task.__name__}, elapsed time: {end_time - start_time:.1f} sec\n")

            # rerun workflow to check whether remediation task resolved existing sequence mismatch errors (DAOTHER-11129)

            if self.__reg.native_combined and self.__reg.dstPath is None and op == 'nmr-str-consistency-check'\
               and (self.__reg.report.error.hasSequenceMismatchErrorInCsLoop()
                    or self.__reg.report.error.hasSequenceMismatchErrorInMrLoop()):

                self.__reg.star_data[0].write_to_file(self.__reg.srcPath,
                                                      show_comments=(self.__reg.bmrb_only and self.__reg.internal_mode),
                                                      skip_empty_loops=True, skip_empty_tags=False)

                for task in self.__procTasksDict['consistency-check']:

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.op() starting op {op} - task {task.__name__}\n")

                    start_time = time.time()

                    if not task():
                        pass

                    if self.__reg.verbose:
                        end_time = time.time()
                        if end_time - start_time > 1.0:
                            self.__reg.log.write(f"op: {op}, task: {task.__name__}, "
                                                 f"elapsed time: {end_time - start_time:.1f} sec\n")

        elif op.endswith('deposit') or op.endswith('release'):

            for task in self.__procTasksDict['deposit']:

                if self.__reg.verbose:
                    self.__reg.log.write(f"+{self.__class_name__}.op() starting op {op} - task {task.__name__}\n")

                start_time = time.time()

                if not task():
                    pass

                if self.__reg.verbose:
                    end_time = time.time()
                    if end_time - start_time > 1.0:
                        self.__reg.log.write(f"op: {op}, task: {task.__name__}, elapsed time: {end_time - start_time:.1f} sec\n")

        # run workflow operation specific tasks
        if op in self.__procTasksDict:

            for task in self.__procTasksDict[op]:

                if self.__reg.verbose:
                    self.__reg.log.write(f"+{self.__class_name__}.op() starting op {op} - task {task.__name__}\n")

                start_time = time.time()

                if not task():
                    if task.__name__ in (self.__translateNef2Str.__name__, self.__translateStr2Nef.__name__):
                        break

                if self.__reg.verbose:
                    end_time = time.time()
                    if end_time - start_time > 1.0:
                        self.__reg.log.write(f"op: {op}, task: {task.__name__}, elapsed time: {end_time - start_time:.1f} sec\n")

        self.__dumpDpReport()

        try:

            input_source = self.__reg.report.input_sources[0]
            input_source_dic = input_source.get()

            if ((op == 'nmr-cs-mr-merge'
                 and self.__reg.report.error.getValueList('missing_mandatory_content',
                                                          input_source_dic['file_name'],
                                                          key='_Atom_chem_shift') is not None)
                or (op in ('nmr-str2str-deposit', 'nmr-str2cif-deposit', 'nmr-str2cif-annotate')
                    and self.__reg.remediation_mode))\
               and self.__reg.report.isError() and self.__reg.dstPath is not None:

                dir_path = os.path.dirname(self.__reg.dstPath) if self.__reg.dirPath is None else self.__reg.dirPath

                rem_dir = os.path.join(dir_path, 'remediation')

                if os.path.isdir(rem_dir):

                    for link_file in os.listdir(rem_dir):

                        link_path = os.path.join(rem_dir, link_file)

                        if os.path.islink(link_path):
                            os.remove(link_path)

                    os.removedirs(rem_dir)

                pk_dir = os.path.join(dir_path, 'nmr_peak_lists')

                if os.path.isdir(pk_dir):

                    for link_file in os.listdir(pk_dir):

                        link_path = os.path.join(pk_dir, link_file)

                        if os.path.islink(link_path):
                            os.remove(link_path)

                    os.removedirs(pk_dir)

            if (self.__reg.submission_mode or self.__reg.annotation_mode or self.__reg.release_mode) and self.__tmpPath is not None:
                os.remove(self.__tmpPath)
                self.__tmpPath = None

            return not self.__reg.report.isError()

        except OSError:
            return False

        finally:
            # NOTE: The specified workspace by setWorkspace() is updated to its default settings
            # after each workflow operation is completed.
            # Therefore, in order to continue applying the workspace settings, you need to configure the workspace again.
            self.__reg.dirPath = None
            self.__reg.cacheDirPath = None

            self.__reg.c2S.dirPath = None
            self.__reg.cR.cacheDirPath = None

            self.__reg.report = None
            self.__reg.report_prev = None

            self.__reg.inputParamDict.clear()
            self.__inputParamDict__ = None
            self.__reg.outputParamDict.clear()

            self.__reg.star_data_type.clear()
            self.__reg.star_data.clear()
            self.__reg.sf_name_corrections.clear()

            self.__reg.original_error_message.clear()
            self.__reg.divide_mr_error_message.clear()
            self.__reg.strip_mr_error_message.clear()

            self.__reg.sf_category_list.clear()
            self.__reg.lp_category_list.clear()

            self.__reg.suspended_errors_for_lazy_eval.clear()
            self.__reg.suspended_warnings_for_lazy_eval.clear()

            for v in self.__reg.lp_data.values():
                v.clear()

            for v in self.__reg.aux_data.values():
                v.clear()

            for v in self.__reg.sf_tag_data.values():
                v.clear()

            # release memory after executing the following tasks, which may involve processing large files
            if (op.startswith('nmr-cs') or op.startswith('nmr-str')) and self.__reg.cifChecked:
                self.release()

    def __dumpDpReport(self) -> bool:
        """ Dump current NMR data processing report.
        """

        if self.__reg.report_prev is not None:
            self.__reg.report.inheritFormatIssueErrors(self.__reg.report_prev)
            self.__reg.report.inheritCorrectedFormatIssueWarnings(self.__reg.report_prev)
            self.__reg.report.inheritCorrectedSaveframeNameWarnings(self.__reg.report_prev)

            if self.__reg.report_prev.error.get() is not None:
                self.__reg.report.setCorrectedError(self.__reg.report_prev)

            if self.__reg.report_prev.warning.get() is not None:
                self.__reg.report.setCorrectedWarning(self.__reg.report_prev)

        if self.__reg.output_statistics is not None:
            self.__reg.report.setOutputStatistics(self.__reg.output_statistics)

        self.__reg.report.error.sortFormatIssueError()
        self.__reg.report.warning.sortChemicalShiftValidation()
        self.__reg.report.warning.sortBySigma('conflicted_data')
        self.__reg.report.warning.sortBySigma('inconsistent_data')

        self.__reg.report.clean()

        if self.__logPath is None:
            return False

        return self.__reg.report.writeFile(self.__logPath)

    def __initializeDpReport(self, srcPath: str = None, calcOutputStats: bool = False) -> bool:
        """ Initialize NMR data processing report.
        """

        srcName = None
        if srcPath is None:
            srcPath = self.__reg.srcPath
            if self.__reg.srcName is not None:
                srcName = self.__reg.srcName

        self.__reg.report = NmrDpReport(self.__reg.verbose, self.__reg.log)

        if REPORT_FILE_PATH_KEY in self.__reg.inputParamDict:
            fPath = self.__reg.inputParamDict[REPORT_FILE_PATH_KEY]

            if os.access(fPath, os.F_OK) and os.path.getsize(fPath) > 0:
                self.__reg.report_prev = NmrDpReport(self.__reg.verbose, self.__reg.log)
                self.__reg.report_prev.loadFile(fPath)
                self.__reg.report.inheritFormatIssueErrors(self.__reg.report_prev)
                self.__reg.report.inheritPreviousErrors(self.__reg.report_prev)
                self.__reg.report.inheritPreviousWarnings(self.__reg.report_prev)

                if calcOutputStats and self.__reg.combined_mode and self.__reg.dstPath is not None:

                    if self.__reg.dstPath == self.__reg.srcPath and self.__reg.release_mode:
                        pass

                    elif not self.__reg.submission_mode and not self.__reg.annotation_mode\
                            or self.__reg.dstPath != self.__reg.srcPath:

                        if not self.__reg.op.endswith('consistency-check'):
                            self.__reg.dpV.calculateOutputStats()

        def proc_cs_file_path_list(offset):

            if not has_key_value(self.__reg.inputParamDict, CS_FILE_PATH_LIST_KEY) and self.__reg.permit_missing_cs_file:
                input_source = self.__reg.report.input_sources[0]
                input_source.setItemValue('file_type', 'nmr-star')
                input_source.setItemValue('content_type', 'nmr-chemical-shifts')
                input_source.setItemValue('ignore_error', False)

                return True

            for csListId, cs in enumerate(self.__reg.inputParamDict[CS_FILE_PATH_LIST_KEY], start=offset):

                if csListId > 0:
                    self.__reg.report.appendInputSource()

                input_source = self.__reg.report.input_sources[csListId]

                file_type = 'nmr-star'  # 'nef' in self.__reg.op else 'nmr-star' # DAOTHER-5673

                if isinstance(cs, str):

                    if cs.endswith('.gz'):

                        _cs = os.path.splitext(cs)[0]

                        if not os.path.exists(_cs):

                            try:

                                uncompress_gzip_file(cs, _cs)

                            except Exception as e:  # pylint: disable=broad-exception-caught

                                self.__reg.report.error.appendDescription('internal_error',
                                                                          f"+{self.__class_name__}.__initializeDpReport() "
                                                                          "++ Error  - " + str(e))

                                if self.__reg.verbose:
                                    self.__reg.log.write(f"+{self.__class_name__}.__initializeDpReport() "
                                                         f"++ Error  - {str(e)}\n")

                                return False

                        cs = _cs

                    if not os.path.basename(cs).startswith('bmr')\
                       and (self.__reg.op == 'nmr-cs-mr-merge'
                            or get_type_of_star_file(cs) == 'cif'
                            or self.__reg.nefT.read_input_file(cs)[1] == 'Saveframe'):

                        input_source.setItemValue('original_file_name', os.path.basename(cs))

                        _cs = self.getNextPath(cs, '.cif2str')
                        if not self.__reg.c2S.convert(cs, _cs):
                            _cs = cs

                        cs = _cs

                    input_source.setItemValue('file_name', re.sub(r'\.cif2str$', '', os.path.basename(cs)))
                    input_source.setItemValue('file_type', file_type)
                    input_source.setItemValue('content_type', 'nmr-chemical-shifts')
                    input_source.setItemValue('ignore_error', False)

                else:

                    if cs['file_name'].endswith('.gz'):

                        _cs = os.path.splitext(cs['file_name'])[0]

                        if not os.path.exists(_cs):

                            try:

                                uncompress_gzip_file(cs['file_name'], _cs)

                            except Exception as e:  # pylint: disable=broad-exception-caught

                                self.__reg.report.error.appendDescription('internal_error',
                                                                          f"+{self.__class_name__}.__initializeDpReport() "
                                                                          "++ Error  - " + str(e))

                                if self.__reg.verbose:
                                    self.__reg.log.write(f"+{self.__class_name__}.__initializeDpReport() "
                                                         f"++ Error  - {str(e)}\n")

                                return False

                        cs['file_name'] = _cs

                    if not os.path.basename(cs['file_name']).startswith('bmr')\
                       and (self.__reg.op == 'nmr-cs-mr-merge'
                            or get_type_of_star_file(cs['file_name']) == 'cif'
                            or self.__reg.nefT.read_input_file(cs['file_name'])[1] == 'Saveframe'):

                        if 'original_file_name' not in cs:
                            input_source.setItemValue('original_file_name', os.path.basename(cs['file_name']))

                        _cs = self.getNextPath(cs['file_name'], '.cif2str')
                        if not self.__reg.c2S.convert(cs['file_name'], _cs, originalFileName=cs.get('original_file_name')):
                            _cs = cs['file_name']

                        cs['file_name'] = _cs

                    input_source.setItemValue('file_name', re.sub(r'\.cif2str$', '', os.path.basename(cs['file_name'])))
                    input_source.setItemValue('file_type', file_type)
                    input_source.setItemValue('content_type', 'nmr-chemical-shifts')
                    if 'original_file_name' in cs:
                        input_source.setItemValue('original_file_name', cs['original_file_name'])
                    input_source.setItemValue('ignore_error', False)

            return True

        def proc_mr_file_path_list():

            if MR_FILE_PATH_LIST_KEY in self.__reg.inputParamDict:

                for mr in self.__reg.inputParamDict[MR_FILE_PATH_LIST_KEY]:

                    self.__reg.report.appendInputSource()

                    input_source = self.__reg.report.input_sources[-1]

                    file_type = 'nmr-star'  # 'nef' if 'nef' in self.__reg.op else 'nmr-star' # DAOTHER-5673

                    if isinstance(mr, str):

                        if get_type_of_star_file(mr) == 'cif'\
                           or self.__reg.nefT.read_input_file(mr)[1] == 'Saveframe':

                            input_source.setItemValue('original_file_name', os.path.basename(mr))

                            _mr = self.getNextPath(mr, '.cif2str')
                            if not self.__reg.c2S.convert(mr, _mr):
                                _mr = mr

                            mr = _mr

                        input_source.setItemValue('file_name', re.sub(r'\.cif2str$', '', os.path.basename(mr)))
                        input_source.setItemValue('file_type', file_type)
                        input_source.setItemValue('content_type', 'nmr-restraints')
                        input_source.setItemValue('ignore_error', False)

                    else:

                        if get_type_of_star_file(mr['file_name']) == 'cif'\
                           or self.__reg.nefT.read_input_file(mr['file_name'])[1] == 'Saveframe':

                            if 'original_file_name' not in mr:
                                input_source.setItemValue('original_file_name', os.path.basename(mr['file_name']))

                            _mr = self.getNextPath(mr['file_name'], '.cif2str')
                            if not self.__reg.c2S.convert(mr['file_name'], _mr, originalFileName=mr.get('original_file_name')):
                                _mr = mr['file_name']

                            mr['file_name'] = _mr

                        input_source.setItemValue('file_name', re.sub(r'\.cif2str$', '', os.path.basename(mr['file_name'])))
                        input_source.setItemValue('file_type', file_type)
                        input_source.setItemValue('content_type', 'nmr-restraints')
                        if 'original_file_name' in mr:
                            input_source.setItemValue('original_file_name', mr['original_file_name'])
                        input_source.setItemValue('ignore_error', False)

        def proc_ar_file_path_list():

            if AR_FILE_PATH_LIST_KEY in self.__reg.inputParamDict:

                for ar in self.__reg.inputParamDict[AR_FILE_PATH_LIST_KEY]:

                    self.__reg.report.appendInputSource()

                    input_source = self.__reg.report.input_sources[-1]

                    arPath = ar['file_name']

                    if arPath.endswith('.gz'):

                        _arPath = os.path.splitext(arPath)[0]

                        if not os.path.exists(_arPath):

                            try:

                                uncompress_gzip_file(arPath, _arPath)

                            except Exception as e:  # pylint: disable=broad-exception-caught

                                self.__reg.report.error.appendDescription('internal_error',
                                                                          f"+{self.__class_name__}.__initializeDpReport() "
                                                                          "++ Error  - " + str(e))

                                if self.__reg.verbose:
                                    self.__reg.log.write(f"+{self.__class_name__}.__initializeDpReport() "
                                                         f"++ Error  - {str(e)}\n")

                                return False

                        arPath = _arPath

                    if ar['file_type'] == 'nm-res-oth':

                        has_ext = in_atoms = False
                        atom_names = 0

                        three_letter_codes = STD_MON_DICT.keys()
                        atom_like_names_oth = self.__reg.csStat.getAtomLikeNameSet(1)

                        with open(arPath, 'r', encoding='utf-8', errors='ignore') as ifh:

                            for line in ifh:

                                if 'EXT' in line:
                                    l_split = line.split()
                                    _line = ' '.join(l_split)

                                    if len(_line) > 1 and _line[0].isdigit() and _line[1] == 'EXT':
                                        has_ext = in_atoms = True
                                        continue

                                elif in_atoms:
                                    l_split = line.split()
                                    _line = ' '.join(l_split)

                                    if len(_line) == 0 or _line.startswith('#') or _line.startswith('!') or _line.startswith(';'):
                                        continue

                                    if len(l_split) >= 10:
                                        try:
                                            atom_num = int(l_split[0])
                                            seq_id = int(l_split[8])
                                            comp_id = l_split[2]
                                            atom_id = l_split[3]
                                            if atom_num > 0 and seq_id > 0 and comp_id in three_letter_codes\
                                               and atom_id in atom_like_names_oth:
                                                atom_names += 1
                                        except ValueError:
                                            pass

                        if has_ext and atom_names > 0:
                            ar['file_type'] = 'nm-aux-cha'

                    if ar['file_type'] in ('nm-res-ari', 'nm-res-oth'):

                        with open(arPath, 'r', encoding='utf-8', errors='ignore') as ifh:

                            for idx, line in enumerate(ifh):
                                if '<!DOCTYPE noe_restraint_list SYSTEM' in line or '<noe_restraint_list>' in line:
                                    ar['file_type'] = 'nm-res-arx'
                                    break
                                if '<!DOCTYPE data_set SYSTEM' in line or '<date_set>' in line:
                                    ar['file_type'] = 'nm-res-arx'
                                    break
                                if idx >= MR_MAX_SPACER_LINES:
                                    break

                    if ar['file_type'] in ('nm-res-cya', 'nm-res-oth'):

                        atom_names = 0

                        three_letter_codes = STD_MON_DICT.keys()
                        atom_like_names_oth = self.__reg.csStat.getAtomLikeNameSet(1)

                        with open(arPath, 'r', encoding='utf-8', errors='ignore') as ifh:

                            for line in ifh:

                                if ' OK ' in line and ' + ' in line:
                                    l_split = line.split()
                                    _line = ' '.join(l_split)

                                    if len(_line) == 0 or _line.startswith('#') or _line.startswith('!') or _line.startswith(';'):
                                        continue

                                    if len(l_split) >= 12:
                                        try:
                                            if l_split[3] == '+' and l_split[7] == 'OK':
                                                seq_id = int(l_split[2])
                                                comp_id = l_split[1]
                                                atom_id = l_split[0]
                                                if seq_id > 0 and comp_id in three_letter_codes\
                                                   and atom_id in atom_like_names_oth:
                                                    seq_id = int(l_split[6])
                                                    comp_id = l_split[5]
                                                    atom_id = l_split[4]
                                                    if seq_id > 0 and comp_id in three_letter_codes\
                                                       and atom_id in atom_like_names_oth:
                                                        atom_names += 1
                                        except ValueError:
                                            pass

                        if atom_names > 0:
                            ar['file_type'] = 'nm-aux-noa'

                    if ar['file_type'] == 'nm-pea-any':

                        for test_file_type in ARCHIVAL_MR_FILE_TYPES:
                            if test_file_type == 'nmr-star':
                                continue
                            if os.path.exists(self.testPathWithSuffix(arPath, f'-selected-as-{test_file_type[-7:]}', True)):
                                ar['file_type'] = test_file_type
                                break

                        if ar['file_type'] == 'nm-pea-any':
                            file_type = get_peak_list_format(arPath, True)

                            if file_type is not None:
                                ar['file_type'] = file_type

                            for test_file_type in PARSABLE_PK_FILE_TYPES:
                                if os.path.exists(self.testPathWithSuffix(arPath, f'-selected-as-{test_file_type[-7:]}', True)):
                                    ar['file_type'] = test_file_type
                                    break

                        if not is_binary_file(arPath) or ar['file_type'] != 'nm-pea-any':

                            codec = detect_bom(arPath, 'utf-8')

                            if codec != 'utf-8':
                                _arPath = self.getNextPath(arPath)
                                convert_codec(arPath, _arPath, codec, 'utf-8')
                                arPath = _arPath

                    input_source.setItemValue('file_name', os.path.basename(arPath))
                    input_source.setItemValue('file_type', ar['file_type'])
                    input_source.setItemValue('content_type',
                                              'nmr-restraints' if not ar['file_type'].startswith('nm-pea') else 'nmr-peaks')
                    if 'original_file_name' in ar:
                        input_source.setItemValue('original_file_name', ar['original_file_name'])
                    input_source.setItemValue('ignore_error', False if 'ignore_error' not in ar else ar['ignore_error'])

            return True

        def proc_ac_file_path_list():

            if AC_FILE_PATH_LIST_KEY in self.__reg.inputParamDict and self.__reg.conversion_server:

                for acs in self.__reg.inputParamDict[AC_FILE_PATH_LIST_KEY]:

                    self.__reg.report.appendInputSource()

                    input_source = self.__reg.report.input_sources[-1]

                    acsPath = acs['file_name']

                    input_source.setItemValue('file_name', os.path.basename(acsPath))
                    input_source.setItemValue('file_type', acs['file_type'])
                    input_source.setItemValue('content_type', 'nmr-chemical-shifts')
                    if 'original_file_name' in acs:
                        input_source.setItemValue('original_file_name', acs['original_file_name'])
                    input_source.setItemValue('ignore_error', False if 'ignore_error' not in acs else acs['ignore_error'])

            if self.__reg.bmrb_only and (self.__reg.internal_mode or self.__reg.conversion_server)\
               and NMR_CIF_FILE_PATH_KEY in self.__reg.inputParamDict:

                nmr_cif = self.__reg.inputParamDict[NMR_CIF_FILE_PATH_KEY]

                _nmr_cif = self.getNextPath(nmr_cif, '.cif2str')
                if self.__reg.c2S.convert(nmr_cif, _nmr_cif):
                    self.__reg.srcNmrCifPath = _nmr_cif
                    self.__reg.native_combined = True  # DAOTHER-8855

        input_source = None

        if self.__reg.combined_mode:

            # set primary input source as NMR unified data
            input_source = self.__reg.report.input_sources[0]

            file_type = 'nef' if 'nef' in self.__reg.op and 'str2nef' not in self.__reg.op else 'nmr-star'
            content_type = CONTENT_TYPE[file_type]

            input_source.setItemValue('file_name', os.path.basename(srcPath))
            input_source.setItemValue('file_type', file_type)
            input_source.setItemValue('content_type', content_type)
            if srcName is not None:
                input_source.setItemValue('original_file_name', srcName)
            input_source.setItemValue('ignore_error', False)

            if self.__reg.op == 'nmr-str-replace-cs'\
               or (self.__reg.op == 'nmr-cs-mr-merge' and self.__reg.bmrb_only and self.__reg.internal_mode):  # DAOTHER-9785
                if self.__reg.conversion_server and self.__reg.combined_mode\
                   and CS_FILE_PATH_LIST_KEY not in self.__reg.inputParamDict:
                    pass
                elif not proc_cs_file_path_list(1):
                    return False

            # DAOTHER-9785
            if self.__reg.op == 'nmr-cs-mr-merge' and self.__reg.conversion_server:
                proc_mr_file_path_list()

                if not proc_ar_file_path_list():
                    return False

                proc_ac_file_path_list()

        else:

            if not proc_cs_file_path_list(0):
                return False

            proc_mr_file_path_list()

            if not proc_ar_file_path_list():
                return False

            proc_ac_file_path_list()

        self.__reg.star_data_type.clear()
        self.__reg.star_data.clear()
        self.__reg.sf_name_corrections.clear()

        self.__reg.original_error_message.clear()

        if self.__parseCoordFilePath():

            try:

                chem_comp = self.__reg.cR.getDictList('chem_comp')

                nstd_comp_ids = [item['id'] for item in chem_comp if item['mon_nstd_flag'] != 'y']

                for comp_id in nstd_comp_ids:

                    if self.__reg.ccU.updateChemCompDict(comp_id):  # matches with comp_id in CCD

                        ref_elems = set(a['type_symbol'] for a in self.__reg.ccU.lastAtomDictList
                                        if a['leaving_atom_flag'] != 'Y')

                        for elem in ref_elems:
                            if elem in PARAMAGNETIC_ELEMENTS or elem in FERROMAGNETIC_ELEMENTS:
                                self.__reg.report.setDiamagnetic(False)
                                break

            except Exception:  # pylint: disable=broad-exception-caught
                pass

        return input_source is not None

    def __validateInputSource(self, srcPath: str = None) -> bool:
        """ Validate NMR data as primary input source.
        """

        return self.__reg.dpV.validateInputSource(srcPath)

    def __detectContentSubType(self) -> bool:
        """ Detect content subtype of NMR data file in any STAR format.
        """

        return self.__reg.dpV.detectContentSubType()

    def __detectContentSubTypeOfLegacyMr(self) -> bool:
        """ Detect content subtype of legacy restraint files.
        """

        if self.__reg.combined_mode:
            return True

        if AR_FILE_PATH_LIST_KEY not in self.__reg.inputParamDict:
            return True

        corrected = self.__reg.dpS.detectContentSubTypeOfLegacyMr(self.__remediation_loop_count)

        # restart using format issue resolved input files

        if self.__reg.remediation_mode and corrected:

            self.__reg.report = None
            self.__reg.report_prev = None

            self.__reg.star_data_type.clear()
            self.__reg.star_data.clear()
            self.__reg.sf_name_corrections.clear()

            self.__reg.original_error_message.clear()

            self.__reg.sf_category_list.clear()
            self.__reg.lp_category_list.clear()

            self.__reg.suspended_errors_for_lazy_eval.clear()
            self.__reg.suspended_warnings_for_lazy_eval.clear()

            for v in self.__reg.lp_data.values():
                v.clear()

            for v in self.__reg.aux_data.values():
                v.clear()

            for v in self.__reg.sf_tag_data.values():
                v.clear()

            self.__reg.inputParamDict = copy.deepcopy(self.__inputParamDict__)

            self.__initializeDpReport()
            self.__validateInputSource()
            self.__detectContentSubType()
            self.__extractPublicMrFileIntoLegacyMr()
            self.__detectContentSubTypeOfLegacyPk()
            self.__detectContentSubTypeOfLegacyMr()

            self.__remediation_loop_count += 1

            self.__reg.sll_pred_holder.clear()

            if self.__reg.mr_debug:
                if self.__remediation_loop_count > 5:
                    self.__reg.log.write(f'repetiation of remediation: {self.__inputParamDict__}\n')

        return not self.__reg.report.isError()

    def __detectContentSubTypeOfLegacyPk(self) -> bool:
        """ Detect content subtype of legacy spectral peak files.
        """

        return self.__reg.dpS.detectContentSubTypeOfLegacyPk()

    def __extractPublicMrFileIntoLegacyMr(self) -> bool:
        """ Extract/split public MR file into legacy restraint files for NMR restraint remediation.
        """

        return self.__reg.dpS.extractPublicMrFileIntoLegacyMr()

    def __getPolymerSequence__(self, file_list_id: int, sf: Union[pynmrstar.Saveframe, pynmrstar.Loop], content_subtype: str
                               ) -> List[List[dict]]:
        """ Wrapper function to retrieve polymer sequence from loop of a specified saveframe and content subtype via NefTranslator.
        """

        input_source = self.__reg.report.input_sources[file_list_id]
        input_source_dic = input_source.get()

        file_type = input_source_dic['file_type']

        check_identity = content_subtype not in self.__reg.mr_content_subtypes

        if file_type == 'nef':  # DAOTHER-7389, issue #3, allow empty for 'chem_shift'
            return self.__reg.nefT.get_nef_seq(sf, lp_category=LP_CATEGORIES[file_type][content_subtype],
                                               allow_empty=(content_subtype in ('chem_shift', 'spectral_peak')),
                                               allow_gap=(content_subtype not in ('poly_seq', 'entity')),
                                               check_identity=check_identity)

        if content_subtype == 'spectral_peak_alt':
            return self.__reg.nefT.get_star_seq(sf, lp_category='_Assigned_peak_chem_shift',
                                                allow_empty=True,
                                                allow_gap=True,
                                                check_identity=check_identity)

        poly_seq = None

        if not self.__reg.bmrb_only or not self.__reg.internal_mode:
            if self.__reg.caC is None:
                self.__retrieveCoordAssemblyChecker__()

        elif self.__reg.caC is None:
            self.__extractPolymerSequenceInEntityAssembly__(file_list_id)

            input_source = self.__reg.report.input_sources[file_list_id]
            input_source_dic = input_source.get()

            if has_key_value(input_source_dic, 'polymer_sequence'):
                poly_seq = input_source_dic['polymer_sequence']

        # DAOTHER-7389, issue #3, allow empty for 'chem_shift'
        return self.__reg.nefT.get_star_seq(sf, lp_category=LP_CATEGORIES[file_type][content_subtype],
                                            allow_empty=(content_subtype in ('chem_shift', 'spectral_peak')),
                                            allow_gap=(content_subtype not in ('poly_seq', 'entity')),
                                            check_identity=check_identity,
                                            coord_assembly_checker=self.__reg.caC if self.__reg.native_combined
                                            or not self.__reg.combined_mode
                                            or self.__reg.op == 'nmr-str-replace-cs'
                                            or (self.__reg.op == 'nmr-str2cif-annotate'
                                                and self.__reg.remediation_mode)  # DAOTHER-10616
                                            else None,
                                            nmr_poly_seq=poly_seq)

    def __extractPolymerSequence(self) -> bool:
        """ Extract reference polymer sequence.
        """

        is_done = True

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            if input_source_dic['content_subtype'] is None:
                is_done = False
                continue

            content_subtype = 'poly_seq'

            if content_subtype not in input_source_dic['content_subtype']:
                continue

            sf_category = SF_CATEGORIES[file_type][content_subtype]
            lp_category = LP_CATEGORIES[file_type][content_subtype]

            for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                try:

                    poly_seq = self.__getPolymerSequence__(fileListId, sf, content_subtype)[0]

                    input_source.setItemValue('polymer_sequence', poly_seq)

                    if file_type == 'nmr-star':

                        auth_poly_seq = self.__reg.nefT.get_star_auth_seq(sf, lp_category)[0]

                        for ps in poly_seq:
                            chain_id, seq_ids, comp_ids = ps['chain_id'], ps['seq_id'], ps['comp_id']

                            for aps in auth_poly_seq:

                                if aps['chain_id'] != chain_id:
                                    continue

                                _seq_ids = aps['seq_id']

                                auth_asym_ids, auth_seq_ids, auth_comp_ids =\
                                    aps['auth_asym_id'], aps['auth_seq_id'], aps['auth_comp_id']

                                auth_asym_id_set = sorted(set(auth_asym_ids))

                                for auth_asym_id in auth_asym_id_set:

                                    offsets = []
                                    total = 0

                                    for _seq_id, _auth_asym_id, auth_seq_id, auth_comp_id in zip(_seq_ids, auth_asym_ids,
                                                                                                 auth_seq_ids, auth_comp_ids):

                                        if _auth_asym_id != auth_asym_id or auth_seq_id in EMPTY_VALUE:
                                            continue

                                        try:

                                            _auth_seq_id = int(auth_seq_id)

                                            offsets.append(_auth_seq_id - _seq_id)
                                            total += 1

                                        except ValueError:

                                            if self.__reg.check_auth_seq:
                                                warn = f"Auth_seq_ID {str(auth_seq_id)!r} "\
                                                    f"(Auth_asym_ID {auth_asym_id}, Auth_comp_ID {auth_comp_id}) "\
                                                    "should be an integer."

                                                self.__reg.report.warning.appendDescription('sequence_mismatch',
                                                                                            {'file_name': file_name,
                                                                                             'sf_framecode': sf_framecode,
                                                                                             'category': lp_category,
                                                                                             'description': warn})

                                                if self.__reg.verbose:
                                                    self.__reg.log.write(f"+{self.__class_name__}.__extractPolymerSequence() "
                                                                         f"++ Warning  - {warn}\n")

                                    if total > 1:

                                        offset = collections.Counter(offsets).most_common()[0][0]

                                        for _seq_id, _auth_asym_id, auth_seq_id, auth_comp_id in zip(_seq_ids, auth_asym_ids,
                                                                                                     auth_seq_ids, auth_comp_ids):

                                            if _auth_asym_id != auth_asym_id or auth_seq_id in EMPTY_VALUE:
                                                continue

                                            try:

                                                _auth_seq_id = int(auth_seq_id)

                                            except ValueError:
                                                continue

                                            if _auth_seq_id - _seq_id != offset:

                                                if self.__reg.check_auth_seq:
                                                    warn = f"Auth_seq_ID {str(auth_seq_id)!r} is inconsistent with "\
                                                        f"{str(_seq_id + offset)!r} "\
                                                        f"(Auth_asym_ID {auth_asym_id}, Auth_comp_ID {auth_comp_id})."

                                                    self.__reg.report.warning.appendDescription('sequence_mismatch',
                                                                                                {'file_name': file_name,
                                                                                                 'sf_framecode': sf_framecode,
                                                                                                 'category': lp_category,
                                                                                                 'description': warn})

                                                    if self.__reg.verbose:
                                                        self.__reg.log.write(f"+{self.__class_name__}.__extractPolymerSequence() "
                                                                             f"++ Warning  - {warn}\n")

                                for seq_id, comp_id in zip(seq_ids, comp_ids):

                                    for _seq_id, _auth_asym_id, auth_seq_id, auth_comp_id in zip(_seq_ids, auth_asym_ids,
                                                                                                 auth_seq_ids, auth_comp_ids):

                                        if _seq_id != seq_id:
                                            continue

                                        if comp_id == auth_comp_id:
                                            continue

                                        if self.__reg.check_auth_seq:
                                            warn = f"Auth_comp_ID {auth_comp_id!r} "\
                                                f"(Auth_asym_ID {_auth_asym_id}, Auth_seq_ID {auth_seq_id}) "\
                                                f"is inconsistent with {comp_id} (Entity_assembly_ID {chain_id}, Seq_ID {seq_id})."

                                            self.__reg.report.warning.appendDescription('sequence_mismatch',
                                                                                        {'file_name': file_name,
                                                                                         'sf_framecode': sf_framecode,
                                                                                         'category': lp_category,
                                                                                         'description': warn})

                                            if self.__reg.verbose:
                                                self.__reg.log.write(f"+{self.__class_name__}.__extractPolymerSequence() "
                                                                     f"++ Warning  - {warn}\n")

                                        break

                    continue

                except KeyError as e:

                    self.__reg.report.error.appendDescription('sequence_mismatch',
                                                              {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                               'category': lp_category, 'description': str(e).strip("'")})

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__extractPolymerSequence() "
                                             f"++ KeyError  - {str(e)}\n")

                except LookupError:
                    # """
                    # self.__reg.report.error.appendDescription('missing_mandatory_item',
                    #                                           {'file_name': file_name, 'sf_framecode': sf_framecode,
                    #                                            'category': lp_category, 'description': str(e).strip("'")})
                    #
                    # self.__reg.log.write(f"+{self.__class_name__}.__extractPolymerSequence() "
                    #                      f"++ LookupError  - {file_name} {sf_framecode} {lp_category} {str(e)}\n")
                    # """
                    pass
                except ValueError as e:

                    self.__reg.report.error.appendDescription('invalid_data',
                                                              {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                               'category': lp_category, 'description': str(e).strip("'")})

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__extractPolymerSequence() "
                                             f"++ ValueError  - {str(e)}\n")

                except UserWarning as e:

                    errs = str(e).strip("'").split('\n')

                    for err in errs:

                        if len(err) == 0:
                            continue

                        if err.startswith('[Invalid data]'):

                            p = err.index(']') + 2
                            err = err[p:]

                            self.__reg.report.error.appendDescription('format_issue',
                                                                      {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                       'category': lp_category, 'description': err})

                            if self.__reg.verbose:
                                self.__reg.log.write(f"+{self.__class_name__}.__extractPolymerSequence() "
                                                     f"++ ValueError  - {err}\n")

                        else:

                            self.__reg.report.error.appendDescription('internal_error',
                                                                      f"+{self.__class_name__}.__extractPolymerSequence() "
                                                                      "++ Error  - " + err)

                            if self.__reg.verbose:
                                self.__reg.log.write(f"+{self.__class_name__}.__extractPolymerSequence() "
                                                     f"++ Error  - {err}\n")

                except Exception as e:  # pylint: disable=broad-exception-caught

                    self.__reg.report.error.appendDescription('internal_error',
                                                              f"+{self.__class_name__}.__extractPolymerSequence() "
                                                              "++ Error  - " + str(e))

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__extractPolymerSequence() "
                                             f"++ Error  - {str(e)}\n")

                is_done = False

        return is_done

    def __extractPolymerSequenceInLoop(self) -> bool:
        """ Extract polymer sequence in interesting loops.
        """

        is_done = True

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            if input_source_dic['content_subtype'] is None:
                is_done = False
                continue

            poly_seq_list_set = {}

            for content_subtype in NMR_CONTENT_SUBTYPES:

                if content_subtype in ('entry_info', 'poly_seq', 'entity', 'ph_param_data')\
                   or (not has_key_value(input_source_dic['content_subtype'], content_subtype)):
                    continue

                poly_seq_list_set[content_subtype] = []

                sf_category = SF_CATEGORIES[file_type][content_subtype]
                lp_category = LP_CATEGORIES[file_type][content_subtype]

                if file_type == 'nmr-star' and content_subtype == 'spectral_peak_alt':
                    lp_category = '_Assigned_peak_chem_shift'

                has_poly_seq = False

                list_id = 1

                if self.__reg.star_data_type[fileListId] == 'Loop':
                    sf = self.__reg.star_data[fileListId]
                    sf_framecode = ''

                    has_poly_seq |= self.__extractPolymerSequenceInLoop__(fileListId, file_name, file_type, content_subtype, sf,
                                                                          list_id, sf_framecode, lp_category, poly_seq_list_set)

                elif self.__reg.star_data_type[fileListId] == 'Saveframe':
                    sf = self.__reg.star_data[fileListId]
                    sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                    has_poly_seq |= self.__extractPolymerSequenceInLoop__(fileListId, file_name, file_type, content_subtype, sf,
                                                                          list_id, sf_framecode, lp_category, poly_seq_list_set)

                else:

                    for list_id, sf in enumerate(self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category), start=1):
                        sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                        if not any(True for loop in sf.loops if loop.category == lp_category):
                            continue

                        has_poly_seq |= self.__extractPolymerSequenceInLoop__(fileListId, file_name, file_type, content_subtype, sf,
                                                                              list_id, sf_framecode, lp_category, poly_seq_list_set)

                if not has_poly_seq:
                    poly_seq_list_set.pop(content_subtype)

            if len(poly_seq_list_set) > 0:
                input_source.setItemValue('polymer_sequence_in_loop', poly_seq_list_set)

        return is_done

    def __extractPolymerSequenceInLoop__(self, file_list_id: int, file_name: str, file_type: str, content_subtype: str,
                                         sf: Union[pynmrstar.Saveframe, pynmrstar.Loop],
                                         list_id: int, sf_framecode: str, lp_category: str, poly_seq_list_set: dict
                                         ) -> bool:
        """ Extract polymer sequence in interesting loops.
        """

        has_poly_seq = False

        try:

            poly_seq = self.__getPolymerSequence__(file_list_id, sf, content_subtype)[0]

            if len(poly_seq) > 0:

                poly_seq_list_set[content_subtype].append({'list_id': list_id,
                                                           'sf_framecode': sf_framecode,
                                                           'polymer_sequence': poly_seq})

                has_poly_seq = True

                if file_type == 'nmr-star':

                    try:
                        auth_poly_seq = self.__reg.nefT.get_star_auth_seq(sf, lp_category)[0]
                    except LookupError:
                        return has_poly_seq

                    for ps in poly_seq:
                        chain_id, seq_ids, comp_ids = ps['chain_id'], ps['seq_id'], ps['comp_id']

                        for aps in auth_poly_seq:

                            if aps['chain_id'] != chain_id:
                                continue

                            _seq_ids = aps['seq_id']

                            auth_asym_ids, auth_seq_ids, auth_comp_ids =\
                                aps['auth_asym_id'], aps['auth_seq_id'], aps['auth_comp_id']

                            auth_asym_id_set = sorted(set(auth_asym_ids))

                            for auth_asym_id in auth_asym_id_set:

                                offsets = []
                                total = 0

                                for _seq_id, _auth_asym_id, auth_seq_id, auth_comp_id\
                                        in zip(_seq_ids, auth_asym_ids, auth_seq_ids, auth_comp_ids):

                                    if _auth_asym_id != auth_asym_id or auth_seq_id in EMPTY_VALUE:
                                        continue

                                    try:

                                        _auth_seq_id = int(auth_seq_id)

                                        offsets.append(_auth_seq_id - _seq_id)
                                        total += 1

                                    except ValueError:

                                        if self.__reg.check_auth_seq:
                                            warn = f"Auth_seq_ID {str(auth_seq_id)!r} (Auth_asym_ID {auth_asym_id}, "\
                                                f"Auth_comp_ID {auth_comp_id}) should be an integer."

                                            self.__reg.report.warning.appendDescription('sequence_mismatch',
                                                                                        {'file_name': file_name,
                                                                                         'sf_framecode': sf_framecode,
                                                                                         'category': lp_category,
                                                                                         'description': warn})

                                            if self.__reg.verbose:
                                                self.__reg.log.write(f"+{self.__class_name__}.__extractPolymerSequenceInLoop() "
                                                                     f"++ Warning  - {warn}\n")

                                if total > 1:

                                    offset = collections.Counter(offsets).most_common()[0][0]

                                    for _seq_id, _auth_asym_id, auth_seq_id, auth_comp_id in zip(_seq_ids, auth_asym_ids,
                                                                                                 auth_seq_ids, auth_comp_ids):

                                        if _auth_asym_id != auth_asym_id or auth_seq_id in EMPTY_VALUE:
                                            continue

                                        try:

                                            _auth_seq_id = int(auth_seq_id)

                                        except ValueError:
                                            continue

                                        if _auth_seq_id - _seq_id != offset:

                                            if self.__reg.check_auth_seq:
                                                warn = f"Auth_seq_ID {str(auth_seq_id)!r} is inconsistent with "\
                                                    f"{str(_seq_id + offset)!r} "\
                                                    f"(Auth_asym_ID {auth_asym_id}, Auth_comp_ID {auth_comp_id})."

                                                self.__reg.report.warning.appendDescription('sequence_mismatch',
                                                                                            {'file_name': file_name,
                                                                                             'sf_framecode': sf_framecode,
                                                                                             'category': lp_category,
                                                                                             'description': warn})

                                                if self.__reg.verbose:
                                                    self.__reg.log.write(f"+{self.__class_name__}.__extractPolymerSequenceInLoop() "
                                                                         f"++ Warning  - {warn}\n")

                            for seq_id, comp_id in zip(seq_ids, comp_ids):

                                for _seq_id, _auth_asym_id, auth_seq_id, auth_comp_id\
                                        in zip(_seq_ids, auth_asym_ids, auth_seq_ids, auth_comp_ids):

                                    if _seq_id != seq_id:
                                        continue

                                    if comp_id == auth_comp_id:
                                        continue

                                    if self.__reg.check_auth_seq:
                                        warn = f"Auth_comp_ID {auth_comp_id!r} "\
                                            f"(Auth_asym_ID {_auth_asym_id}, Auth_seq_ID {auth_seq_id}) "\
                                            f"is inconsistent with {comp_id!r} (Entity_assembly_ID {chain_id}, Seq_ID {seq_id})."

                                        self.__reg.report.warning.appendDescription('sequence_mismatch',
                                                                                    {'file_name': file_name,
                                                                                     'sf_framecode': sf_framecode,
                                                                                     'category': lp_category,
                                                                                     'description': warn})

                                        if self.__reg.verbose:
                                            self.__reg.log.write(f"+{self.__class_name__}.__extractPolymerSequenceInLoop() "
                                                                 f"++ Warning  - {warn}\n")

                                    break

        except KeyError as e:

            if 'Auth' not in str(e) or self.__reg.check_auth_seq:
                self.__reg.report.error.appendDescription('sequence_mismatch',
                                                          {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                           'category': lp_category, 'description': str(e).strip("'")})

                if self.__reg.verbose:
                    self.__reg.log.write(f"+{self.__class_name__}.__extractPolymerSequenceInLoop() "
                                         f"++ KeyError  - {str(e)}\n")

        except LookupError as e:

            self.__reg.report.error.appendDescription('missing_mandatory_item',
                                                      {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                       'category': lp_category, 'description': str(e).strip("'")})

            self.__reg.log.write(f"+{self.__class_name__}.__extractPolymerSequenceInLoop() "
                                 f"++ LookupError  - {file_name} {sf_framecode} {lp_category} {str(e)}\n")

        except ValueError as e:

            self.__reg.report.error.appendDescription('invalid_data',
                                                      {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                       'category': lp_category, 'description': str(e).strip("'")})

            if self.__reg.verbose:
                self.__reg.log.write(f"+{self.__class_name__}.__extractPolymerSequenceInLoop() "
                                     f"++ ValueError  - {str(e)}\n")

        except UserWarning as e:

            errs = str(e).strip("'").split('\n')

            for err in errs:

                if len(err) == 0:
                    continue

                if err.startswith('[Invalid data]'):

                    if not content_subtype.startswith('spectral_peak'):

                        _err = err.split('#')[0]

                        if 'Auth' not in _err or self.__reg.check_auth_seq:

                            p = err.index(']') + 2
                            err = err[p:]

                            self.__reg.report.error.appendDescription('format_issue',
                                                                      {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                       'category': lp_category, 'description': err})

                            if self.__reg.verbose:
                                self.__reg.log.write(f"+{self.__class_name__}.__extractPolymerSequenceInLoop() "
                                                     f"++ ValueError  - {err}\n")

                else:

                    self.__reg.report.error.appendDescription('internal_error',
                                                              f"+{self.__class_name__}.__extractPolymerSequenceInLoop() "
                                                              "++ Error  - " + err)

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__extractPolymerSequenceInLoop() "
                                             f"++ Error  - {err}\n")

        except Exception as e:  # pylint: disable=broad-exception-caught

            self.__reg.report.error.appendDescription('internal_error',
                                                      f"+{self.__class_name__}.__extractPolymerSequenceInLoop() "
                                                      "++ Error  - " + str(e))

            if self.__reg.verbose:
                self.__reg.log.write(f"+{self.__class_name__}.__extractPolymerSequenceInLoop() "
                                     f"++ Error  - {str(e)}\n")

        return has_poly_seq

    def __testSequenceConsistency(self) -> bool:
        """ Perform sequence consistency test among extracted polymer sequences.
        """

        if self.__valid_seq:
            return True

        update_poly_seq = False

        poly_seq0 = None

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            has_poly_seq = has_key_value(input_source_dic, 'polymer_sequence')
            has_poly_seq_in_lp = has_key_value(input_source_dic, 'polymer_sequence_in_loop')

            if (not has_poly_seq) or (not has_poly_seq_in_lp):
                continue

            poly_seq = input_source_dic['polymer_sequence']
            poly_seq_in_lp = input_source_dic['polymer_sequence_in_loop']

            subtype_with_poly_seq = ['poly_seq' if has_poly_seq else None]

            for subtype in poly_seq_in_lp.keys():
                subtype_with_poly_seq.append(subtype)

            if self.__reg.bmrb_only and self.__reg.internal_mode:

                to_entity_id = {}

                for sf in self.__reg.star_data[fileListId].get_saveframes_by_category('assembly'):

                    try:

                        loop = sf.get_loop('_Entity_assembly')

                        if loop is not None:
                            dat = loop.get_tag(['ID', 'Entity_ID'])

                            for row in dat:
                                to_entity_id[row[0]] = row[1] if isinstance(row[1], str) else str(row[1])

                    except KeyError:
                        pass

            for subtype_pair in itertools.combinations_with_replacement(subtype_with_poly_seq, 2):

                # poly_seq is reference sequence and suppress tests on combinations of two sequences in loop
                if has_poly_seq and ('poly_seq' not in subtype_pair or subtype_pair == ('poly_seq', 'poly_seq')):
                    continue

                subtype1 = subtype_pair[0]  # poly_seq will appear only on subtype1
                subtype2 = subtype_pair[1]

                if None in (subtype1, subtype2):
                    continue

                lp_category2 = LP_CATEGORIES[file_type][subtype2]

                if file_type == 'nmr-star':
                    if subtype2 == 'spectral_peak_alt':
                        lp_category2 = '_Assigned_peak_chem_shift'

                # reference polymer sequence exists
                if has_poly_seq and subtype1 == 'poly_seq':

                    if fileListId == 0 and poly_seq0 is None:
                        poly_seq0 = poly_seq

                    poly_seq1 = poly_seq0 if fileListId > 0 and poly_seq0 is not None else poly_seq

                    ref_chain_ids = {ps1['chain_id'] for ps1 in poly_seq1}

                    for _poly_seq_in_lp in poly_seq_in_lp[subtype2]:
                        poly_seq2 = _poly_seq_in_lp['polymer_sequence']
                        sf_framecode2 = _poly_seq_in_lp['sf_framecode']

                        for ps2 in poly_seq2:
                            chain_id = ps2['chain_id']

                            if self.__reg.bmrb_only and self.__reg.internal_mode\
                               and chain_id not in ref_chain_ids\
                               and not ('identical_chain_id' in ps2 and chain_id not in ps2['identical_chain_id']):

                                chain_id = to_entity_id.get(chain_id, chain_id)

                            if chain_id not in ref_chain_ids and not chain_id.isdigit()\
                               and self.__reg.combined_mode and self.__reg.caC is not None:
                                chain_id = next((str(item['entity_assembly_id']) for item in self.__reg.caC['entity_assembly']
                                                 if chain_id in item['auth_asym_id'].split(',')), chain_id)

                            if chain_id not in ref_chain_ids\
                               and not ('identical_chain_id' in ps2 and chain_id not in ps2['identical_chain_id']):

                                if self.__reg.op == 'nmr-str-replace-cs':
                                    seq_align, _ = alignPolymerSequence(self.__reg.pA, [poly_seq1[0]], [poly_seq2[0]])
                                    if len(seq_align) > 0:
                                        sa = seq_align[0]
                                        if sa['matched'] > 0 and sa['conflict'] == 0:
                                            continue

                                err = f"Invalid chain_id {chain_id!r} in a loop {lp_category2}."

                                single_poly = False
                                if self.__reg.caC is not None:
                                    cif_ps = self.__reg.caC['polymer_sequence']
                                    cif_br = self.__reg.caC['branched']
                                    cif_np = self.__reg.caC['non_polymer']
                                    single_poly = len(cif_ps) == 1 and cif_br is None and cif_np is None\
                                        and len(poly_seq1) == 1 and len(poly_seq2) == 1

                                if self.__reg.has_star_entity and single_poly:
                                    err = f"The sequence between the entity and {lp_category2} loop is inconsistent.\n"\
                                        f"entity: {poly_seq1[0]}\nloop: {poly_seq2[0]}"

                                if self.__reg.remediation_mode\
                                   and not (self.__reg.has_star_entity and single_poly):  # DAOTHER-10487

                                    self.__reg.report.warning.appendDescription('sequence_mismatch',
                                                                                {'file_name': file_name,
                                                                                 'sf_framecode': sf_framecode2,
                                                                                 'category': lp_category2,
                                                                                 'description': err})

                                    if self.__reg.verbose:
                                        self.__reg.log.write(f"+{self.__class_name__}.__testSequenceConsistency() "
                                                             f"++ Warning  - {err}\n")

                                else:

                                    self.__reg.report.error.appendDescription('sequence_mismatch',
                                                                              {'file_name': file_name,
                                                                               'sf_framecode': sf_framecode2,
                                                                               'category': lp_category2,
                                                                               'description': err})

                                    if self.__reg.verbose:
                                        self.__reg.log.write(f"+{self.__class_name__}.__testSequenceConsistency() "
                                                             f"++ Error  - {err}\n")

                            else:

                                for ps1 in poly_seq1:

                                    if ps1['chain_id'] != chain_id\
                                       and not ('identical_chain_id' in ps2 and ps1['chain_id'] in ps2['identical_chain_id']):
                                        continue

                                    if 'identical_chain_id' in ps2:
                                        _ps1 = next((_ps1 for _ps1 in poly_seq1
                                                     if _ps1['chain_id'] == chain_id), None)
                                        __ps1 = next((_ps1 for _ps1 in poly_seq1
                                                      if _ps1['chain_id'] in ps2['identical_chain_id']), None)
                                        if _ps1 is not None and len(ps1['seq_id']) != len(_ps1['seq_id']):
                                            continue
                                        if __ps1 is not None and len(ps1['seq_id']) != len(__ps1['seq_id']):
                                            continue

                                    has_gap_in_auth_seq = False
                                    if self.__reg.caC is not None:
                                        cif_ps = self.__reg.caC['polymer_sequence']
                                        nmr_ps = [ps1]

                                        seq_align, _ = alignPolymerSequence(self.__reg.pA, cif_ps, nmr_ps)
                                        if len(seq_align) > 0:
                                            chain_assign, _ = assignPolymerSequence(self.__reg.pA, self.__reg.ccU, 'nmr-star',
                                                                                    cif_ps, nmr_ps, seq_align)

                                            for ca in chain_assign:
                                                if ca['matched'] == 0 or ca['conflict'] > 0:
                                                    continue
                                                ref_chain_id = ca['ref_chain_id']
                                                ps = next(ps for ps in cif_ps if ps['auth_chain_id'] == ref_chain_id)
                                                has_gap_in_auth_seq = 'gap_in_auth_seq' in ps and ps['gap_in_auth_seq']
                                                break

                                    if any(True for seq_id, comp_id in zip(ps2['seq_id'], ps2['comp_id'])
                                           if seq_id in ps1['seq_id'] and comp_id not in EMPTY_VALUE
                                           and ps1['comp_id'][ps1['seq_id'].index(seq_id)] not in EMPTY_VALUE
                                           and comp_id != ps1['comp_id'][ps1['seq_id'].index(seq_id)]):
                                        seq_align, _ = alignPolymerSequence(self.__reg.pA, [ps1], [ps2])

                                        if len(seq_align) == 1:
                                            sa = seq_align[0]
                                            if sa['matched'] > 0 and sa['conflict'] == 0:
                                                continue
                                        else:
                                            has_alt_seq_align = False
                                            for ps0 in poly_seq1:
                                                if ps0['chain_id'] == ps1['chain_id']:
                                                    continue
                                                seq_align, _ = alignPolymerSequence(self.__reg.pA, [ps0], [ps2])
                                                if len(seq_align) == 1:
                                                    sa = seq_align[0]
                                                    if sa['matched'] > 0 and sa['conflict'] == 0:
                                                        has_alt_seq_align = True
                                                        break

                                            if has_alt_seq_align:
                                                continue

                                            for c in range(1, MAX_CONFLICT_ATTEMPT):
                                                seq_align, _ = alignPolymerSequenceWithConflicts(self.__reg.pA, [ps1], [ps2], c)

                                                if len(seq_align) == 1:
                                                    sa = seq_align[0]
                                                    if sa['matched'] > sa['conflict'] and sa['conflict'] <= c:
                                                        _ps1 = copy.copy(ps1)
                                                        _ps2 = copy.copy(ps2)
                                                        _ps1['seq_id'] = sa['test_seq_id']
                                                        _ps1['comp_id'] = sa['ref_comp_id']
                                                        _ps2['seq_id'] = sa['test_seq_id']
                                                        _ps2['comp_id'] = sa['test_comp_id']
                                                        ps1, ps2 = _ps1, _ps2
                                                        break

                                    for seq_id, comp_id in zip(ps2['seq_id'], ps2['comp_id']):

                                        if seq_id not in ps1['seq_id']:

                                            if comp_id != '.':

                                                if not self.__reg.remediation_mode\
                                                   and ((min(set(ps2['seq_id']) - set(ps1['seq_id'])) > 0 and seq_id > 0)
                                                        or not self.__reg.nonblk_bad_nterm):

                                                    err = f"Invalid seq_id {str(seq_id)!r} "\
                                                        f"(chain_id {chain_id}) in a loop {lp_category2}."

                                                    self.__reg.report.error.appendDescription('sequence_mismatch',
                                                                                              {'file_name': file_name,
                                                                                               'sf_framecode': sf_framecode2,
                                                                                               'category': lp_category2,
                                                                                               'description': err})

                                                    if self.__reg.verbose:
                                                        self.__reg.log.write(f"+{self.__class_name__}.__testSequenceConsistency() "
                                                                             f"++ Error  - {err}\n")

                                                else:

                                                    warn = f"Unmapped seq_id {str(seq_id)!r} (chain_id {chain_id}) "\
                                                        f"in a loop {lp_category2}. {INSTRUCTION_FOR_FULL_SEQUENCE}"

                                                    self.__reg.report.warning.appendDescription('sequence_mismatch',
                                                                                                {'file_name': file_name,
                                                                                                 'sf_framecode': sf_framecode2,
                                                                                                 'category': lp_category2,
                                                                                                 'description': warn})

                                                    if self.__reg.verbose:
                                                        self.__reg.log.write(f"+{self.__class_name__}.__testSequenceConsistency() "
                                                                             f"++ Warning  - {warn}\n")

                                        else:

                                            try:
                                                _comp_id = ps1['comp_id'][ps1['seq_id'].index(seq_id)]
                                            except IndexError:
                                                continue

                                            if comp_id not in EMPTY_VALUE and _comp_id not in EMPTY_VALUE and comp_id != _comp_id:

                                                err = f"Invalid comp_id {comp_id!r} vs {_comp_id!r} "\
                                                    f"(seq_id {seq_id}, chain_id {chain_id}) in a loop {lp_category2}."

                                                if (self.__reg.tolerant_seq_align
                                                    and self.__reg.dpV.equalsToRepCompId(comp_id, _comp_id))\
                                                   or (self.__reg.remediation_mode
                                                       and (self.__valid_seq or comp_id not in STD_MON_DICT
                                                            or has_gap_in_auth_seq)):
                                                    self.__reg.report.warning.appendDescription('sequence_mismatch',
                                                                                                {'file_name': file_name,
                                                                                                 'sf_framecode': sf_framecode2,
                                                                                                 'category': lp_category2,
                                                                                                 'description': err})

                                                    if self.__reg.verbose:
                                                        self.__reg.log.write(f"+{self.__class_name__}.__testSequenceConsistency() "
                                                                             f"++ Warning  - {err}\n")

                                                elif self.__reg.tolerant_seq_align\
                                                        and getOneLetterCodeCan(comp_id) == getOneLetterCodeCan(_comp_id):
                                                    self.__reg.report.warning.appendDescription('sequence_mismatch',
                                                                                                {'file_name': file_name,
                                                                                                 'sf_framecode': sf_framecode2,
                                                                                                 'category': lp_category2,
                                                                                                 'description': err})

                                                    if self.__reg.verbose:
                                                        self.__reg.log.write(f"+{self.__class_name__}.__testSequenceConsistency() "
                                                                             f"++ Warning  - {err}\n")

                                                    comp_id_conv_dict = {comp_id: _comp_id}

                                                    self.__reg.dpR.fixCompIdInLoop(fileListId, file_type, subtype2, sf_framecode2,
                                                                                   chain_id, seq_id, comp_id_conv_dict)

                                                    update_poly_seq = True

                                                else:
                                                    self.__reg.report.error.appendDescription('sequence_mismatch',
                                                                                              {'file_name': file_name,
                                                                                               'sf_framecode': sf_framecode2,
                                                                                               'category': lp_category2,
                                                                                               'description': err})

                                                    if self.__reg.verbose:
                                                        self.__reg.log.write(f"+{self.__class_name__}.__testSequenceConsistency() "
                                                                             f"++ Error  - {err}\n")

                # brute force check
                else:

                    for _poly_seq_in_lp in poly_seq_in_lp[subtype1]:
                        poly_seq1 = _poly_seq_in_lp['polymer_sequence']
                        sf_framecode1 = _poly_seq_in_lp['sf_framecode']

                        for _poly_seq_in_lp2 in poly_seq_in_lp[subtype2]:
                            poly_seq2 = _poly_seq_in_lp2['polymer_sequence']
                            sf_framecode2 = _poly_seq_in_lp2['sf_framecode']

                            # suppress redundant tests inside the same subtype
                            if subtype1 == subtype2 and _poly_seq_in_lp['list_id'] >= _poly_seq_in_lp2['list_id']:
                                continue

                            for ps2 in poly_seq2:
                                chain_id = ps2['chain_id']

                                for ps1 in poly_seq1:

                                    if chain_id != ps1['chain_id']:
                                        continue

                                    for seq_id, comp_id in zip(ps2['seq_id'], ps2['comp_id']):

                                        if seq_id in ps1['seq_id']:
                                            _comp_id = ps1['comp_id'][ps1['seq_id'].index(seq_id)]

                                            if comp_id not in EMPTY_VALUE and _comp_id not in EMPTY_VALUE and comp_id != _comp_id:

                                                err = f"Unmatched comp_id {comp_id!r} vs {_comp_id!r} "\
                                                    f"(seq_id {seq_id}, chain_id {chain_id}) exists "\
                                                    f"against {sf_framecode1!r} saveframe."

                                                if self.__reg.tolerant_seq_align\
                                                   and self.__reg.dpV.equalsToRepCompId(comp_id, _comp_id):
                                                    self.__reg.report.warning.appendDescription('sequence_mismatch',
                                                                                                {'file_name': file_name,
                                                                                                 'sf_framecode': sf_framecode2,
                                                                                                 'category': lp_category2,
                                                                                                 'description': err})

                                                    if self.__reg.verbose:
                                                        self.__reg.log.write(f"+{self.__class_name__}.__testSequenceConsistency() "
                                                                             f"++ Warning  - {err}\n")

                                                else:
                                                    self.__reg.report.error.appendDescription('sequence_mismatch',
                                                                                              {'file_name': file_name,
                                                                                               'sf_framecode': sf_framecode2,
                                                                                               'category': lp_category2,
                                                                                               'description': err})

                                                    if self.__reg.verbose:
                                                        self.__reg.log.write(f"+{self.__class_name__}.__testSequenceConsistency() "
                                                                             f"++ Error  - {err}\n")

                            # inverse check required for unverified sequences
                            for ps1 in poly_seq1:
                                chain_id = ps1['chain_id']

                                for ps2 in poly_seq2:

                                    if chain_id != ps2['chain_id']:
                                        continue

                                    for seq_id, comp_id in zip(ps1['seq_id'], ps1['comp_id']):

                                        if seq_id in ps2['seq_id']:
                                            j = ps2['seq_id'].index(seq_id)
                                            _comp_id = ps2['comp_id'][j]

                                            if comp_id not in EMPTY_VALUE and _comp_id not in EMPTY_VALUE and comp_id != _comp_id:

                                                err = f"Unmatched comp_id {comp_id!r} vs {_comp_id!r} "\
                                                    f"(seq_id {seq_id}, chain_id {chain_id}) exists "\
                                                    f"against {sf_framecode2!r} saveframe."

                                                if self.__reg.tolerant_seq_align\
                                                   and self.__reg.dpV.equalsToRepCompId(comp_id, _comp_id):
                                                    self.__reg.report.warning.appendDescription('sequence_mismatch',
                                                                                                {'file_name': file_name,
                                                                                                 'sf_framecode': sf_framecode2,
                                                                                                 'category': lp_category2,
                                                                                                 'description': err})

                                                    if self.__reg.verbose:
                                                        self.__reg.log.write(f"+{self.__class_name__}.__testSequenceConsistency() "
                                                                             f"++ Warning  - {err}\n")

                                                else:
                                                    self.__reg.report.error.appendDescription('sequence_mismatch',
                                                                                              {'file_name': file_name,
                                                                                               'sf_framecode': sf_framecode2,
                                                                                               'category': lp_category2,
                                                                                               'description': err})

                                                    if self.__reg.verbose:
                                                        self.__reg.log.write(f"+{self.__class_name__}.__testSequenceConsistency() "
                                                                             f"++ Error  - {err}\n")

        if update_poly_seq:
            self.__extractPolymerSequenceInLoop()
            self.__depositNmrData()

        return not self.__reg.report.isError()

    def __extractCommonPolymerSequence(self) -> bool:
        """ Extract common polymer sequence if required.
        """

        common_poly_seq = {}

        cs_has_alt_comp_id = False

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            content_type = input_source_dic['content_type']

            has_poly_seq = has_key_value(input_source_dic, 'polymer_sequence')
            has_poly_seq_in_lp = has_key_value(input_source_dic, 'polymer_sequence_in_loop')

            # pass if poly_seq exists
            if has_poly_seq or (not has_poly_seq_in_lp and not self.__reg.conversion_server):

                if not has_poly_seq or not has_poly_seq_in_lp:
                    continue

                self.__mergePolymerSequenceInCsLoop__(fileListId)

                continue

            if self.__extractPolymerSequenceInEntityAssembly__(fileListId):
                continue

            if not has_poly_seq_in_lp:
                continue

            poly_seq_in_lp = input_source_dic['polymer_sequence_in_loop']

            content_subtype = 'chem_shift'

            # DAOTHER-7545 NMR-STAR formatted MR has no chem_shift
            if content_subtype not in poly_seq_in_lp or content_type == 'nmr-restraints':

                if 'dist_restraint' in poly_seq_in_lp:
                    content_subtype = 'dist_restraint'
                elif 'dihed_restraint' in poly_seq_in_lp:
                    content_subtype = 'dihed_restraint'
                elif 'rdc_restraint' in poly_seq_in_lp:
                    content_subtype = 'rdc_restraint'
                else:
                    continue

            for _poly_seq_in_lp in poly_seq_in_lp[content_subtype]:
                poly_seq = _poly_seq_in_lp['polymer_sequence']

                for ps in poly_seq:
                    chain_id = ps['chain_id']

                    if chain_id not in common_poly_seq:
                        common_poly_seq[chain_id] = set()

            chain_ids = common_poly_seq.keys()
            offset_seq_ids = {c: 0 for c in chain_ids}

            for _poly_seq_in_lp in poly_seq_in_lp[content_subtype]:
                poly_seq = _poly_seq_in_lp['polymer_sequence']

                for ps in poly_seq:
                    chain_id = ps['chain_id']

                    min_seq_id = min(ps['seq_id'])
                    if min_seq_id < 0:
                        offset_seq_ids[chain_id] = min_seq_id * -1

                    if 'alt_comp_id' in ps:
                        for seq_id, comp_id, alt_comp_id in zip(ps['seq_id'], ps['comp_id'], ps['alt_comp_id']):
                            common_poly_seq[chain_id].add((seq_id + offset_seq_ids[chain_id], comp_id, alt_comp_id))
                            if (seq_id + offset_seq_ids[chain_id], comp_id) in common_poly_seq[chain_id]:
                                common_poly_seq[chain_id].remove((seq_id + offset_seq_ids[chain_id], comp_id))
                            cs_has_alt_comp_id = True
                    else:
                        for seq_id, comp_id in zip(ps['seq_id'], ps['comp_id']):
                            if not any(True for item in common_poly_seq[chain_id]
                                       if item[0] == seq_id + offset_seq_ids[chain_id] and item[1] == comp_id and len(item) == 3):
                                common_poly_seq[chain_id].add((seq_id + offset_seq_ids[chain_id], comp_id))

        def get_array_expr(c):
            seq_ids = sorted(set(item[0] - offset_seq_ids[c] for item in common_poly_seq[c]))
            comp_ids, alt_comp_ids = [], []

            for seq_id in seq_ids:
                _comp_ids = [item[1] for item in common_poly_seq[c] if item[0] - offset_seq_ids[c] == seq_id]
                if cs_has_alt_comp_id:
                    _alt_comp_ids = [item[len(item) - 1] for item in common_poly_seq[c]
                                     if item[0] - offset_seq_ids[c] == seq_id]
                if len(_comp_ids) == 1:
                    comp_ids.append(_comp_ids[0])
                    if cs_has_alt_comp_id:
                        alt_comp_ids.append(_alt_comp_ids[0])
                else:
                    comp_ids.append(next(comp_id for comp_id in _comp_ids if comp_id not in EMPTY_VALUE))
                    if cs_has_alt_comp_id:
                        alt_comp_ids.append(next(alt_comp_id for alt_comp_id in _alt_comp_ids
                                                 if alt_comp_id not in EMPTY_VALUE))
            return seq_ids, comp_ids, alt_comp_ids

        asm, failed_entity = [], []  # molecular assembly of a loop

        for chain_id in sorted(common_poly_seq.keys()):

            if len(common_poly_seq[chain_id]) > 0:
                seq_ids, comp_ids, alt_comp_ids = get_array_expr(chain_id)

                if self.__reg.combined_mode and self.__reg.has_star_entity:
                    ent = self.__extractPolymerSequenceInEntityLoopOfChain__(fileListId, chain_id)

                    if ent is not None:
                        if any(comp_id in comp_ids for comp_id in ent['comp_id']):
                            asm.append(ent)

                        else:
                            failed_entity.append(chain_id)

                    elif len(failed_entity) > 0:
                        hit = False
                        for _chain_id in copy.copy(failed_entity):
                            ent = self.__extractPolymerSequenceInEntityLoopOfChain__(fileListId, _chain_id)

                            if ent is not None and any(comp_id in comp_ids for comp_id in ent['comp_id']):
                                _seq_ids, _comp_ids, _alt_comp_ids = get_array_expr(_chain_id)
                                _ent = {'chain_id': _chain_id, 'seq_id': _seq_ids, 'comp_id': _comp_ids}
                                if cs_has_alt_comp_id:
                                    _ent['alt_comp_id'] = _alt_comp_ids

                                asm.append(_ent)

                                ent['chain_id'] = chain_id
                                asm.append(ent)

                                failed_entity.remove(_chain_id)

                                hit = True
                                break

                        if not hit:
                            failed_entity.append(chain_id)

                    continue

                ent = {'chain_id': chain_id, 'seq_id': seq_ids, 'comp_id': comp_ids}
                if cs_has_alt_comp_id:
                    ent['alt_comp_id'] = alt_comp_ids

                asm.append(ent)

        if len(asm) > 0:

            for fileListId in range(self.__reg.file_path_list_len):

                input_source = self.__reg.report.input_sources[fileListId]
                input_source_dic = input_source.get()

                has_poly_seq = has_key_value(input_source_dic, 'polymer_sequence')
                has_poly_seq_in_lp = has_key_value(input_source_dic, 'polymer_sequence_in_loop')

                # pass if poly_seq exists
                if has_poly_seq or (not has_poly_seq_in_lp):
                    continue

                if self.__extractPolymerSequenceInEntityAssembly__(fileListId):
                    continue

                input_source.setItemValue('polymer_sequence', asm)

        return True

    def __mergePolymerSequenceInCsLoop__(self, file_list_id: int) -> bool:
        """ Merge polymer sequence in CS loops.
        """

        input_source = self.__reg.report.input_sources[file_list_id]
        input_source_dic = input_source.get()

        has_poly_seq = has_key_value(input_source_dic, 'polymer_sequence')
        has_poly_seq_in_lp = has_key_value(input_source_dic, 'polymer_sequence_in_loop')

        if not has_poly_seq or not has_poly_seq_in_lp:
            return False

        poly_seq = input_source_dic['polymer_sequence']
        poly_seq_in_lp = input_source_dic['polymer_sequence_in_loop']

        content_subtype = 'chem_shift'

        if content_subtype not in poly_seq_in_lp:
            return False

        ext_seq_key_set = set()

        for poly_seq_in_cs_loop in poly_seq_in_lp[content_subtype]:
            for ps_in_cs_loop in poly_seq_in_cs_loop['polymer_sequence']:
                _chain_id = ps_in_cs_loop['chain_id']

                _peptide = _nucleotide = _carbohydrate = False
                for _comp_id in ps_in_cs_loop['comp_id']:
                    _peptide, _nucleotide, _carbohydrate = self.__reg.csStat.getTypeOfCompId(_comp_id)
                    if _peptide or _nucleotide or _carbohydrate:
                        break

                for ps in poly_seq:
                    chain_id = ps['chain_id']

                    if chain_id == _chain_id\
                       or 'identical_chain_id' in ps and _chain_id in ps['identical_chain_id']:

                        peptide = nucleotide = carbohydrate = False
                        for comp_id in ps['comp_id']:
                            peptide, nucleotide, carbohydrate = self.__reg.csStat.getTypeOfCompId(comp_id)
                            if peptide or nucleotide or carbohydrate:
                                break

                        if peptide != _peptide or nucleotide != _nucleotide or carbohydrate != _carbohydrate:
                            continue

                        if 'alt_comp_id' in ps_in_cs_loop:
                            for _seq_id, _comp_id, _alt_comp_id in zip(ps_in_cs_loop['seq_id'], ps_in_cs_loop['comp_id'],
                                                                       ps_in_cs_loop['alt_comp_id']):
                                if _seq_id not in ps['seq_id'] and _seq_id is not None:
                                    ext_seq_key_set.add((chain_id, _seq_id, _comp_id, _alt_comp_id))
                        else:
                            for _seq_id, _comp_id in zip(ps_in_cs_loop['seq_id'], ps_in_cs_loop['comp_id']):
                                if _seq_id not in ps['seq_id'] and _seq_id is not None:
                                    ext_seq_key_set.add((chain_id, _seq_id, _comp_id))

        if len(ext_seq_key_set) > 0:

            cs_has_alt_comp_id = False

            for ext_seq_key in ext_seq_key_set:
                ps = next(ps for ps in poly_seq if ps['chain_id'] == ext_seq_key[0])

                seq_id = ext_seq_key[1]
                comp_id = ext_seq_key[2]

                if len(ext_seq_key) > 3:
                    cs_has_alt_comp_id = True

                if seq_id in ps['seq_id']:
                    continue

                pos = None

                if ps['seq_id'][0] is not None and seq_id < ps['seq_id'][0]:
                    pos = 0
                elif ps['seq_id'][-1] is not None and seq_id > ps['seq_id'][-1]:
                    pos = len(ps['seq_id'])
                else:
                    for idx, _seq_id in enumerate(ps['seq_id']):
                        if _seq_id is None or seq_id > _seq_id:
                            continue
                        pos = idx
                        break

                if pos is not None:
                    ps['seq_id'].insert(pos, seq_id)
                    ps['comp_id'].insert(pos, comp_id)

            if cs_has_alt_comp_id:
                chain_ids = []
                for ext_seq_key in ext_seq_key_set:
                    if ext_seq_key[0] not in chain_ids:
                        chain_ids.append(ext_seq_key[0])

                for chain_id in chain_ids:
                    ps = next(ps for ps in poly_seq if ps['chain_id'] == chain_id)
                    if 'alt_comp_id' not in ps or len(ps['alt_comp_id']) != len(ps['comp_id']):
                        ps['alt_comp_id'] = list(ps['comp_id'])
                    for ext_seq_key in ext_seq_key_set:
                        if ext_seq_key[0] != chain_id:
                            continue
                        if ext_seq_key[1] in ps['seq_id']:
                            cs_has_alt_comp_id = len(ext_seq_key) > 3
                            pos = ps['seq_id'].index(ext_seq_key[1])
                            ps['alt_comp_id'][pos] = ext_seq_key[3 if cs_has_alt_comp_id else 2]

        if self.__reg.op == 'nmr-str-replace-cs':
            self.__valid_seq = False
            self.__testSequenceConsistency()

        return True

    def __extractPolymerSequenceInEntityAssembly__(self, file_list_id: int) -> bool:
        """ Extract polymer sequence in entity loops. (NMR combined deposition)
        """

        input_source = self.__reg.report.input_sources[file_list_id]
        input_source_dic = input_source.get()

        file_type = input_source_dic['file_type']

        if file_type != 'nmr-star' or not self.__reg.has_star_entity:
            return False

        if not self.__reg.combined_mode:
            return self.__extractPolymerSequenceInEntityLoop__(file_list_id)

        if file_list_id >= len(self.__reg.star_data):
            return False

        for sf in self.__reg.star_data[file_list_id].get_saveframes_by_category('assembly'):

            try:
                loop = sf.get_loop('_Entity_assembly')
            except KeyError:
                return False

            if loop is None:
                return False

            tags = ['ID', 'Entity_assembly_name', 'Entity_ID']

            if 'Entity_label' in loop.tags:
                tags.append('Entity_label')

            if set(tags) & set(loop.tags) != set(tags):
                return False

            dat = loop.get_tag(tags)

            asm = []  # molecular assembly of a loop

            chain_ids = set()
            entity_sfs = {}

            for c in dat:

                if c[0] in EMPTY_VALUE or c[1] in EMPTY_VALUE or c[2] in EMPTY_VALUE:
                    return False

                try:
                    chain_id = str(c[0])
                    # Entity_assemble_name or Entity_label
                    entity_sf = c[1] if len(c) < 4 else (c[3][1:] if c[3][0] == '$' else c[3])
                    entity_id = int(c[2])

                    if chain_id in chain_ids:
                        return False

                    chain_ids.add(chain_id)

                    for k, v in entity_sfs.items():
                        if (k != entity_sf and v == entity_id) or (k == entity_sf and v != entity_id):
                            return False

                    entity_sfs[entity_sf] = entity_id

                except ValueError:
                    return False

                _sf = self.__reg.dpA.getSaveframeByName(file_list_id, entity_sf)

                if _sf is None:
                    return False

                content_subtype = 'entity'

                try:
                    _loop = _sf.get_loop(LP_CATEGORIES[file_type][content_subtype])
                except KeyError:
                    return False

                if _loop is None:
                    return False

                _tags = ['ID', 'Comp_ID', 'Entity_ID']

                if set(_tags) & set(_loop.tags) != set(_tags):
                    return False

                _dat = _loop.get_tag(_tags)

                seq = set()

                for _row in _dat:

                    if _row[0] in EMPTY_VALUE or _row[1] in EMPTY_VALUE or _row[2] in EMPTY_VALUE:
                        return False

                    try:
                        seq_id = int(_row[0])
                        comp_id = _row[1].upper()
                        _entity_id = int(_row[2])
                    except ValueError:
                        return False

                    if entity_id != _entity_id:
                        return False

                    seq.add((seq_id, comp_id))

                sorted_seq = sorted(seq, key=itemgetter(0))

                asm.append({'chain_id': chain_id,
                            'seq_id': [x[0] for x in sorted_seq],
                            'comp_id': [x[1] for x in sorted_seq]})

            if len(asm) > 0:
                input_source.setItemValue('polymer_sequence', asm)
                return True

            break

        return False

    def __extractPolymerSequenceInEntityLoop__(self, file_list_id: int) -> bool:
        """ Extract polymer sequence in entity loops. (NMR conventional deposition)
        """

        input_source = self.__reg.report.input_sources[file_list_id]
        input_source_dic = input_source.get()

        file_type = input_source_dic['file_type']

        if file_type != 'nmr-star' or not self.__reg.has_star_entity:
            return False

        if self.__reg.combined_mode:
            return self.__extractPolymerSequenceInEntityAssembly__(file_list_id)

        star_data = self.__reg.star_data[file_list_id]

        content_subtype = 'entity'

        if self.__reg.star_data_type[file_list_id] != 'Loop':
            for sf in self.__reg.star_data[file_list_id].get_saveframes_by_category(content_subtype):
                entity_id = get_first_sf_tag(sf, 'ID')
                if len(entity_id) == 0 or entity_id in EMPTY_VALUE or not entity_id.isdigit():
                    entity_id = 1
                    set_sf_tag(sf, 'ID', entity_id)
                else:
                    entity_id = int(entity_id)
                self.__reg.c2S.set_local_sf_id(sf, str(entity_id))

        lp_category = LP_CATEGORIES[file_type][content_subtype]

        try:
            loops = star_data.get_loops_by_category(lp_category)
        except AttributeError:
            try:
                loops = [star_data.get_loop(lp_category)]
            except AttributeError:
                return False

        for sf in self.__reg.star_data[file_list_id].get_saveframes_by_category('assembly'):

            try:
                loop = sf.get_loop('_Entity_assembly')
            except KeyError:
                return False

            if loop is None:
                return False

            dat = loop.get_tag(['Entity_ID'])

            entity_id_set = set()

            for row in dat:
                if row not in EMPTY_VALUE:
                    entity_id_set.add(row)

            # DAOTHER-8800: make sure all _Entity saveframes have an _Entity Comp Index loop before relying on these loops
            if len(loops) != len(entity_id_set):
                return False

        asm = []  # molecular assembly of a loop

        chain_ids = set()
        seq = {}

        for loop in loops:

            if loop is None:
                continue

            tags = ['ID', 'Comp_ID', 'Entity_ID']
            tags_ = ['ID', 'Comp_ID']

            dat = []

            if set(tags) & set(loop.tags) == set(tags):
                dat = loop.get_tag(tags)
                for row in dat:
                    if row[2] in EMPTY_VALUE:
                        row[2] = '1'
            elif set(tags_) & set(loop.tags) == set(tags_):  # No Entity_ID tag case
                dat = loop.get_tag(tags_)
                for row in dat:
                    row.append('1')

            for row in dat:

                if row[0] in EMPTY_VALUE or row[1] in EMPTY_VALUE or row[2] in EMPTY_VALUE:
                    return False

                try:
                    c = str(row[2])

                    chain_ids.add(c)
                    if c not in seq:
                        seq[c] = set()
                    seq[c].add((int(row[0]), row[1].upper()))
                except (ValueError, TypeError):
                    return False

        for chain_id in chain_ids:

            sorted_seq = sorted(seq[chain_id], key=itemgetter(0))

            asm.append({'chain_id': chain_id,
                        'seq_id': [x[0] for x in sorted_seq],
                        'comp_id': [x[1] for x in sorted_seq]})

        if len(asm) > 0:
            input_source.setItemValue('polymer_sequence', asm)
            return True

        return False

    def __extractPolymerSequenceInEntityLoopOfChain__(self, file_list_id: int, chain_id: str
                                                      ) -> Optional[dict]:
        """ Extract polymer sequence in entity loops of a given chain id.
        """

        input_source = self.__reg.report.input_sources[file_list_id]
        input_source_dic = input_source.get()

        file_type = input_source_dic['file_type']

        if file_type != 'nmr-star' or not self.__reg.has_star_entity:
            return None

        star_data = self.__reg.star_data[file_list_id]

        content_subtype = 'entity'

        lp_category = LP_CATEGORIES[file_type][content_subtype]

        try:
            loops = star_data.get_loops_by_category(lp_category)
        except AttributeError:
            try:
                loops = [star_data.get_loop(lp_category)]
            except AttributeError:
                return None

        chain_ids = set()
        seq = {}

        for loop in loops:

            if loop is None:
                continue

            tags = ['ID', 'Comp_ID', 'Entity_ID']
            tags_ = ['ID', 'Comp_ID']

            dat = []

            if set(tags) & set(loop.tags) == set(tags):
                dat = loop.get_tag(tags)
                for row in dat:
                    if row[2] in EMPTY_VALUE:
                        row[2] = '1'
            elif set(tags_) & set(loop.tags) == set(tags_):  # No Entity_ID tag case
                dat = loop.get_tag(tags_)
                for row in dat:
                    row.append('1')

            for row in dat:

                if row[0] in EMPTY_VALUE or row[1] in EMPTY_VALUE or row[2] in EMPTY_VALUE:
                    continue

                try:
                    c = str(row[2])

                    chain_ids.add(c)
                    if c not in seq:
                        seq[c] = set()
                    seq[c].add((int(row[0]), row[1].upper()))
                except (ValueError, TypeError):
                    return None

        if chain_id in chain_ids:

            sorted_seq = sorted(seq[chain_id], key=itemgetter(0))

            return {'chain_id': chain_id,
                    'seq_id': [x[0] for x in sorted_seq],
                    'comp_id': [x[1] for x in sorted_seq]}

        return None

    def __extractNonStandardResidue(self) -> bool:
        """ Extract non-standard residue.
        """

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            content_subtype = 'poly_seq'

            has_poly_seq = has_key_value(input_source_dic, 'polymer_sequence')

            if not has_poly_seq:
                continue

            sf_category = SF_CATEGORIES[file_type][content_subtype]
            lp_category = LP_CATEGORIES[file_type][content_subtype]

            if self.__reg.star_data_type[fileListId] == 'Loop':
                sf = self.__reg.star_data[fileListId]
                sf_framecode = ''

                self.__extractNonStandardResidue__(file_name, sf_framecode, lp_category, input_source)

            elif self.__reg.star_data_type[fileListId] == 'Saveframe':
                sf = self.__reg.star_data[fileListId]
                sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                self.__extractNonStandardResidue__(file_name, sf_framecode, lp_category, input_source)

            else:

                for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                    sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                    if not any(True for loop in sf.loops if loop.category == lp_category):
                        continue

                    self.__extractNonStandardResidue__(file_name, sf_framecode, lp_category, input_source)

        return True

    def __extractNonStandardResidue__(self, file_name: str, sf_framecode: str, lp_category: str,
                                      input_source: NmrDpReportInputSource) -> None:
        """ Extract non-standard residue.
        """

        input_source_dic = input_source.get()

        poly_seq = input_source_dic['polymer_sequence']

        asm = []

        for ps in poly_seq:

            has_nstd_res = False

            ent = {'chain_id': ps['chain_id'], 'seq_id': [], 'comp_id': [], 'chem_comp_name': [], 'exptl_data': []}

            for seq_id, comp_id in zip(ps['seq_id'], ps['comp_id']):

                if comp_id not in STD_MON_DICT:
                    has_nstd_res = True

                    ent['seq_id'].append(seq_id)
                    ent['comp_id'].append(comp_id)

                    is_valid, cc_name, cc_rel_status = self.__reg.dpV.getChemCompNameAndStatusOf(comp_id)

                    if is_valid:  # matches with comp_id in CCD
                        if cc_rel_status == 'REL' or cc_name is not None:
                            ent['chem_comp_name'].append(cc_name)
                        else:
                            ent['chem_comp_name'].append(f"(Not available due to CCD status code {cc_rel_status})")

                    else:
                        ent['chem_comp_name'].append(cc_name)

                        if comp_id != '.':
                            warn = f"Non standard residue ({ps['chain_id']}:{seq_id}:{comp_id}) "\
                                "did not match with chemical component dictionary (CCD)."

                            self.__reg.report.warning.appendDescription('ccd_mismatch',
                                                                        {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                         'category': lp_category, 'description': warn})

                            if self.__reg.verbose:
                                self.__reg.log.write(f"+{self.__class_name__}.__extractNonStandardResidue() ++ Warning  - {warn}\n")

                        # DAOTHER-9065
                        else:
                            warn = f"Residue ({ps['chain_id']}:{seq_id}:{comp_id}) was not specified. "\
                                f"{INSTRUCTION_FOR_FULL_SEQUENCE}"

                            self.__reg.report.warning.appendDescription('sequence_mismatch',
                                                                        {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                         'category': lp_category, 'description': warn})

                            if self.__reg.verbose:
                                self.__reg.log.write(f"+{self.__class_name__}.__extractNonStandardResidue() ++ Warning  - {warn}\n")

                    ent['exptl_data'].append({'chem_shift': False, 'dist_restraint': False, 'dihed_restraint': False,
                                              'rdc_restraint': False, 'spectral_peak': False, 'coordinate': False})

            if has_nstd_res:
                asm.append(ent)

        if len(asm) > 0:
            input_source.setItemValue('non_standard_residue', asm)

    def __appendPolymerSequenceAlignment(self) -> bool:
        """ Append polymer sequence alignment of interesting loops.
        """

        is_done = True
        update_poly_seq = False

        self.__alt_chain = False
        self.__valid_seq = False

        if not self.__reg.tolerant_seq_align:
            self.__valid_seq = self.__reg.dpV.isConsistentSequence()

            if not self.__valid_seq:
                self.__reg.tolerant_seq_align = True

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_type = input_source_dic['file_type']

            has_poly_seq = has_key_value(input_source_dic, 'polymer_sequence')
            has_poly_seq_in_lp = has_key_value(input_source_dic, 'polymer_sequence_in_loop')

            if not has_poly_seq:
                is_done = False
                continue

            if not has_poly_seq_in_lp:
                continue

            poly_seq = input_source_dic['polymer_sequence']
            poly_seq_in_lp = input_source_dic['polymer_sequence_in_loop']

            for content_subtype in poly_seq_in_lp.keys():

                seq_align_set = []

                dst_chain_ids, ref_chain_ids, map_chain_ids, map_seq_ids, proc_chain_ids =\
                    {}, {}, {}, {}, {}

                for ps1 in poly_seq:
                    chain_id = ps1['chain_id']

                    for _poly_seq_in_lp in poly_seq_in_lp[content_subtype]:
                        poly_seq2 = _poly_seq_in_lp['polymer_sequence']
                        sf_framecode2 = _poly_seq_in_lp['sf_framecode']

                        for ps2 in poly_seq2:

                            if sf_framecode2 in ref_chain_ids and chain_id in ref_chain_ids[sf_framecode2]:
                                continue

                            chain_id2 = ps2['chain_id']

                            if chain_id != chain_id2:
                                continue

                            _ps2 = fillBlankCompIdWithOffset(ps2, 0)

                            if len(_ps2['seq_id']) > len(ps2['seq_id']) and len(_ps2['seq_id']) < len(ps1['seq_id']):
                                ps2 = _ps2

                            self.__reg.pA.clear()
                            self.__reg.pA.setReferenceSequence(ps1['comp_id'], f'REF{chain_id}')
                            self.__reg.pA.addTestSequence(ps2['comp_id'], chain_id)
                            self.__reg.pA.doAlign()

                            myAlign = self.__reg.pA.getAlignment(chain_id)

                            length = len(myAlign)

                            if length == 0:
                                continue

                            _matched, unmapped, conflict, offset_1, offset_2 = getScoreOfSeqAlign(myAlign)

                            alt_chain = False

                            if length == unmapped + conflict or _matched <= conflict + (1 if length > 1 else 0)\
                               or (len(poly_seq) > 1 and _matched < 4 and offset_1 > 0):

                                if self.__reg.tolerant_seq_align and _matched <= conflict + (1 if length > 1 else 0)\
                                   and len(poly_seq) > 1:

                                    __length = length
                                    __matched = _matched
                                    __unmapped = unmapped
                                    __conflict = conflict
                                    __chain_id = __ps1 = __offset_1 = __offset_2 = None

                                    for _ps1 in poly_seq:

                                        if _ps1 == ps1:
                                            continue

                                        chain_id_ = _ps1['chain_id']

                                        if sf_framecode2 in ref_chain_ids and chain_id_ in ref_chain_ids[sf_framecode2]:
                                            continue

                                        self.__reg.pA.clear()
                                        self.__reg.pA.setReferenceSequence(_ps1['comp_id'], f'REF{chain_id_}')
                                        self.__reg.pA.addTestSequence(ps2['comp_id'], chain_id_)
                                        self.__reg.pA.doAlign()

                                        myAlign = self.__reg.pA.getAlignment(chain_id_)

                                        length = len(myAlign)

                                        if length == 0:
                                            continue

                                        _matched, unmapped, conflict, offset_1, offset_2 = getScoreOfSeqAlign(myAlign)

                                        if length == unmapped + conflict or _matched <= conflict + (1 if length > 1 else 0):
                                            continue

                                        if _matched - conflict < __matched - __conflict\
                                           or unmapped + conflict > __unmapped + __conflict:
                                            continue

                                        __length = length
                                        __matched = _matched
                                        __unmapped = unmapped
                                        __conflict = conflict
                                        __chain_id = chain_id_
                                        __offset_1 = offset_1
                                        __offset_2 = offset_2
                                        __ps1 = copy.copy(_ps1)

                                        alt_chain = True

                                        break

                                if not alt_chain\
                                   or (sf_framecode2 in dst_chain_ids and __chain_id in dst_chain_ids[sf_framecode2])\
                                   or (sf_framecode2 in map_chain_ids and chain_id in map_chain_ids[sf_framecode2]):
                                    continue

                                if sf_framecode2 not in dst_chain_ids:
                                    dst_chain_ids[sf_framecode2] = set()

                                dst_chain_ids[sf_framecode2].add(__chain_id)

                                if sf_framecode2 not in map_chain_ids:
                                    map_chain_ids[sf_framecode2] = {}

                                map_chain_ids[sf_framecode2][chain_id] = __chain_id

                                length = __length
                                _matched = __matched
                                unmapped = __unmapped
                                conflict = __conflict
                                chain_id = __ps1['chain_id']
                                chain_id = __chain_id
                                offset_1 = __offset_1
                                offset_2 = __offset_2
                                ps1 = __ps1

                                update_poly_seq = True

                            if conflict == 0 and self.__alt_chain and not alt_chain and chain_id != ps2['chain_id']\
                               and (sf_framecode2 not in dst_chain_ids or chain_id not in dst_chain_ids[sf_framecode2])\
                               and (sf_framecode2 not in map_chain_ids or ps2['chain_id'] not in map_chain_ids[sf_framecode2])\
                               and unmapped != offset_1 + 1 and unmapped != offset_2 + 1\
                               and unmapped <= _matched + offset_1 and unmapped <= _matched + offset_2:

                                if sf_framecode2 not in dst_chain_ids:
                                    dst_chain_ids[sf_framecode2] = set()

                                dst_chain_ids[sf_framecode2].add(chain_id)

                                if sf_framecode2 not in map_chain_ids:
                                    map_chain_ids[sf_framecode2] = {}

                                map_chain_ids[sf_framecode2][ps2['chain_id']] = chain_id

                                alt_chain = True

                            _ps1 = ps1 if offset_1 == 0 else fillBlankCompIdWithOffset(ps1, offset_1)
                            _ps2 = ps2 if offset_2 == 0 else fillBlankCompIdWithOffset(ps2, offset_2)

                            if conflict > 0 and _ps1['seq_id'][0] < 0 and _ps2['seq_id'][0] < 0:  # noqa: E501, pylint: disable=chained-comparison,line-too-long
                                continue

                            if self.__reg.tolerant_seq_align:  # and not alt_chain:
                                seq_mismatch = any((__s1, __s2) for (__s1, __s2, __c1, __c2)
                                                   in zip(_ps1['seq_id'], _ps2['seq_id'], _ps1['comp_id'], _ps2['comp_id'])
                                                   if __s1 != '.' and __s2 != '.' and __s1 != __s2
                                                   and __c1 != '.' and __c2 != '.' and __c1 == __c2)
                                comp_mismatch = any((__s1, __s2) for (__s1, __s2, __c1, __c2)
                                                    in zip(_ps1['seq_id'], _ps2['seq_id'], _ps1['comp_id'], _ps2['comp_id'])
                                                    if __c1 != '.' and __c2 != '.' and __c1 != __c2)

                                if seq_mismatch or comp_mismatch:
                                    if all(__c1 == _ps2['comp_id'][_ps2['seq_id'].index(__s1)] for __s1, __c1
                                           in zip(_ps1['seq_id'], _ps1['comp_id'])
                                           if __s1 != '.' and __c1 != '.' and __s1 in _ps2['seq_id']
                                           and _ps2['comp_id'][_ps2['seq_id'].index(__s1)] != '.'):
                                        if all(__c2 == _ps1['comp_id'][_ps1['seq_id'].index(__s2)] for __s2, __c2
                                               in zip(_ps2['seq_id'], _ps2['comp_id'])
                                               if __s2 != '.' and __c2 != '.' and __s2 in _ps1['seq_id']
                                               and _ps1['comp_id'][_ps1['seq_id'].index(__s2)] != '.'):
                                            continue

                                if 0 < _matched < 4 and unmapped // _matched > 20 and seq_mismatch and len(poly_seq) > 1:
                                    continue

                            if not alt_chain:
                                if sf_framecode2 not in dst_chain_ids:
                                    dst_chain_ids[sf_framecode2] = set()

                                dst_chain_ids[sf_framecode2].add(chain_id)

                            if sf_framecode2 not in ref_chain_ids:
                                ref_chain_ids[sf_framecode2] = []

                            ref_chain_ids[sf_framecode2].append(chain_id)

                            ref_length = len(ps1['seq_id'])

                            ref_code = getOneLetterCodeCanSequence(_ps1['comp_id'])
                            test_code = getOneLetterCodeCanSequence(_ps2['comp_id'])
                            mid_code = getMiddleCode(ref_code, test_code)
                            ref_gauge_code = getGaugeCode(_ps1['seq_id'])
                            test_gauge_code = getGaugeCode(_ps2['seq_id'])

                            self.__alt_chain |= alt_chain

                            if self.__reg.tolerant_seq_align and (seq_mismatch or comp_mismatch):  # and not alt_chain:
                                if sf_framecode2 not in map_seq_ids:
                                    map_seq_ids[sf_framecode2] = set()
                                map_seq_ids[sf_framecode2].add(chain_id)
                                if _ps2['seq_id'] == list(range(_ps2['seq_id'][0], _ps2['seq_id'][-1] + 1)):
                                    seq_id_conv_dict = {str(__s2): str(__s1) for __s1, __s2
                                                        in zip(_ps1['seq_id'], _ps2['seq_id']) if __s2 != '.'}
                                    if comp_mismatch:
                                        _seq_align =\
                                            self.__getSeqAlignCodeWithChainRemap__(fileListId, file_type,
                                                                                   content_subtype, sf_framecode2,
                                                                                   chain_id, _ps1, _ps2, myAlign,
                                                                                   map_chain_ids.get(sf_framecode2),
                                                                                   ref_gauge_code, ref_code, mid_code,
                                                                                   test_code, test_gauge_code)
                                        _ps2['seq_id'] = _seq_align['test_seq_id']
                                        if _ps1['seq_id'][0] < 0 and _ps2['seq_id'][0] < 0:
                                            continue
                                        ref_gauge_code = _seq_align['ref_gauge_code']
                                        ref_code = _seq_align['ref_code']
                                        mid_code = _seq_align['mid_code']
                                        test_code = _seq_align['test_code']
                                        test_gauge_code = _seq_align['test_gauge_code']
                                    else:
                                        chain_id2 = chain_id
                                        if sf_framecode2 in map_chain_ids and chain_id in map_chain_ids[sf_framecode2].values():
                                            chain_id2 = next(k for k, v in map_chain_ids[sf_framecode2].items() if v == chain_id)

                                        if sf_framecode2 not in proc_chain_ids:
                                            proc_chain_ids[sf_framecode2] = set()

                                        if chain_id2 not in proc_chain_ids[sf_framecode2]:
                                            self.__reg.dpR.fixSeqIdInLoop(fileListId, file_type, content_subtype,
                                                                          sf_framecode2, chain_id2, seq_id_conv_dict)
                                            proc_chain_ids[sf_framecode2].add(chain_id2)

                                            if 'identical_chain_id' in ps2:
                                                for chain_id2_ in ps2['identical_chain_id']:
                                                    if chain_id2_ not in proc_chain_ids[sf_framecode2]:
                                                        self.__reg.dpR.fixSeqIdInLoop(fileListId, file_type, content_subtype,
                                                                                      sf_framecode2, chain_id2_, seq_id_conv_dict)
                                                        proc_chain_ids[sf_framecode2].add(chain_id2_)

                                        _ps2['seq_id'] = _ps1['seq_id']
                                        mid_code = getMiddleCode(ref_code, test_code)
                                        test_gauge_code = ref_gauge_code
                                else:
                                    if seq_mismatch:
                                        _seq_align =\
                                            self.__getSeqAlignCodeWithChainRemap__(fileListId, file_type,
                                                                                   content_subtype, sf_framecode2,
                                                                                   chain_id, _ps1, _ps2, myAlign,
                                                                                   map_chain_ids.get(sf_framecode2),
                                                                                   ref_gauge_code, ref_code, mid_code,
                                                                                   test_code, test_gauge_code)
                                        _ps2['seq_id'] = _seq_align['test_seq_id']
                                        if _ps1['seq_id'][0] < 0 and _ps2['seq_id'][0] < 0:
                                            continue
                                        ref_gauge_code = _seq_align['ref_gauge_code']
                                        ref_code = _seq_align['ref_code']
                                        mid_code = _seq_align['mid_code']
                                        test_code = _seq_align['test_code']
                                        test_gauge_code = _seq_align['test_gauge_code']
                                    else:
                                        _ps2 = fillBlankCompId(_ps1, _ps2)
                                        if _ps1['seq_id'][0] < 0 and _ps2['seq_id'][0] < 0:
                                            continue
                                        test_code = getOneLetterCodeCanSequence(_ps2['comp_id'])
                                        mid_code = getMiddleCode(ref_code, test_code)
                                        test_gauge_code = ref_gauge_code

                                update_poly_seq = True

                            matched = mid_code.count('|')

                            if self.__reg.tolerant_seq_align and len(poly_seq) > 1:  # and not alt_chain:
                                if 0 < matched < 4 and unmapped // matched > 20:
                                    continue

                            seq_align = {'list_id': _poly_seq_in_lp['list_id'], 'sf_framecode': sf_framecode2, 'chain_id': chain_id,
                                         'length': ref_length, 'matched': matched, 'conflict': conflict, 'unmapped': unmapped,
                                         'sequence_coverage': round(float(length - (unmapped + conflict)) / ref_length, 3),
                                         'ref_seq_id': _ps1['seq_id'], 'test_seq_id': _ps2['seq_id'],
                                         'ref_gauge_code': ref_gauge_code, 'ref_code': ref_code, 'mid_code': mid_code,
                                         'test_code': test_code, 'test_gauge_code': test_gauge_code}

                            if seq_align in seq_align_set:
                                continue

                            seq_align_set.append(seq_align)

                            if not self.__reg.combined_mode and input_source_dic['non_standard_residue'] is None:
                                has_nstd_res = False
                                for j, rc in enumerate(ref_code):
                                    if rc == 'X' and j < len(test_code) and test_code[j] == 'X':
                                        has_nstd_res = True
                                        break

                                if not has_nstd_res:
                                    continue

                                asm = []

                                for _ps in poly_seq:

                                    ent = {'chain_id': _ps['chain_id'], 'seq_id': [], 'comp_id': [],
                                           'chem_comp_name': [], 'exptl_data': []}

                                    for _seq_id, _comp_id in zip(_ps['seq_id'], _ps['comp_id']):

                                        if _comp_id not in STD_MON_DICT:

                                            ent['seq_id'].append(_seq_id)
                                            ent['comp_id'].append(_comp_id)

                                            is_valid, cc_name, cc_rel_status = self.__reg.dpV.getChemCompNameAndStatusOf(_comp_id)

                                            if is_valid:  # matches with comp_id in CCD
                                                if cc_rel_status == 'REL' or cc_name is not None:
                                                    ent['chem_comp_name'].append(cc_name)
                                                else:
                                                    ent['chem_comp_name'].append(f"(Not available due to CCD status code {cc_rel_status})")  # noqa: E501, pylint: disable=line-too-long

                                            else:
                                                ent['chem_comp_name'].append(cc_name)

                                            ent['exptl_data'].append({'coordinate': False})

                                    asm.append(ent)

                                input_source.setItemValue('non_standard_residue', asm)

                            for r_code, t_code, seq_id in zip(ref_code, test_code, ps1['seq_id']):
                                if r_code == 'X' and t_code == 'X':
                                    input_source.updateNonStandardResidueByExptlData(chain_id, seq_id, content_subtype)

                for ps1 in poly_seq:

                    chain_id = ps1['chain_id']

                    for _poly_seq_in_lp in poly_seq_in_lp[content_subtype]:
                        poly_seq2 = _poly_seq_in_lp['polymer_sequence']
                        sf_framecode2 = _poly_seq_in_lp['sf_framecode']

                        for ps2 in poly_seq2:

                            if sf_framecode2 in ref_chain_ids and chain_id in ref_chain_ids[sf_framecode2]:
                                continue

                            chain_id2 = ps2['chain_id']

                            if sf_framecode2 in dst_chain_ids and chain_id2 in dst_chain_ids[sf_framecode2]:
                                continue

                            if chain_id != chain_id2 and not self.__reg.tolerant_seq_align:
                                continue

                            _ps2 = fillBlankCompIdWithOffset(ps2, 0)

                            if len(_ps2['seq_id']) > len(ps2['seq_id']) and len(_ps2['seq_id']) < len(ps1['seq_id']):
                                ps2 = _ps2

                            self.__reg.pA.clear()
                            self.__reg.pA.setReferenceSequence(ps1['comp_id'], f'REF{chain_id}')
                            self.__reg.pA.addTestSequence(ps2['comp_id'], chain_id)
                            self.__reg.pA.doAlign()

                            myAlign = self.__reg.pA.getAlignment(chain_id)

                            length = len(myAlign)

                            if length == 0:
                                continue

                            _matched, unmapped, conflict, offset_1, offset_2 = getScoreOfSeqAlign(myAlign)

                            alt_chain = False

                            if length == unmapped + conflict or _matched <= conflict + (1 if length > 1 else 0)\
                               or (len(poly_seq) > 1 and _matched < 4 and offset_1 > 0):

                                if self.__reg.tolerant_seq_align and _matched <= conflict + (1 if length > 1 else 0)\
                                   and len(poly_seq) > 1:

                                    __length = length
                                    __matched = _matched
                                    __unmapped = unmapped
                                    __conflict = conflict
                                    __chain_id = __ps1 = __offset_1 = __offset_2 = None

                                    for _ps1 in poly_seq:

                                        if _ps1 == ps1:
                                            continue

                                        chain_id_ = _ps1['chain_id']

                                        if sf_framecode2 in ref_chain_ids and chain_id_ in ref_chain_ids[sf_framecode2]:
                                            continue

                                        self.__reg.pA.clear()
                                        self.__reg.pA.setReferenceSequence(_ps1['comp_id'], f'REF{chain_id_}')
                                        self.__reg.pA.addTestSequence(ps2['comp_id'], chain_id_)
                                        self.__reg.pA.doAlign()

                                        myAlign = self.__reg.pA.getAlignment(chain_id_)

                                        length = len(myAlign)

                                        if length == 0:
                                            continue

                                        _matched, unmapped, conflict, offset_1, offset_2 = getScoreOfSeqAlign(myAlign)

                                        if length == unmapped + conflict or _matched <= conflict + (1 if length > 1 else 0):
                                            continue

                                        if _matched - conflict < __matched - __conflict\
                                           or (unmapped + conflict > __unmapped + __conflict and __matched > 0):
                                            continue

                                        __length = length
                                        __matched = _matched
                                        __unmapped = unmapped
                                        __conflict = conflict
                                        __chain_id = chain_id_
                                        __offset_1 = offset_1
                                        __offset_2 = offset_2
                                        __ps1 = copy.copy(_ps1)

                                        alt_chain = True

                                        break

                                if not alt_chain\
                                   or (sf_framecode2 in dst_chain_ids and __chain_id in dst_chain_ids[sf_framecode2])\
                                   or (sf_framecode2 in map_chain_ids and chain_id in map_chain_ids[sf_framecode2]):
                                    continue

                                if sf_framecode2 not in dst_chain_ids:
                                    dst_chain_ids[sf_framecode2] = set()

                                dst_chain_ids[sf_framecode2].add(__chain_id)

                                if sf_framecode2 not in map_chain_ids:
                                    map_chain_ids[sf_framecode2] = {}

                                map_chain_ids[sf_framecode2][chain_id] = __chain_id

                                length = __length
                                _matched = __matched
                                unmapped = __unmapped
                                conflict = __conflict
                                chain_id = __ps1['chain_id']
                                chain_id = __chain_id
                                offset_1 = __offset_1
                                offset_2 = __offset_2
                                ps1 = __ps1

                                update_poly_seq = True

                            _ps1 = ps1 if offset_1 == 0 else fillBlankCompIdWithOffset(ps1, offset_1)
                            _ps2 = ps2 if offset_2 == 0 else fillBlankCompIdWithOffset(ps2, offset_2)

                            if conflict > 0 and _ps1['seq_id'][0] < 0 and _ps2['seq_id'][0] < 0:  # noqa: E501, pylint: disable=chained-comparison,line-too-long
                                continue

                            if self.__reg.tolerant_seq_align:  # and not alt_chain:
                                seq_mismatch = any((__s1, __s2) for (__s1, __s2, __c1, __c2)
                                                   in zip(_ps1['seq_id'], _ps2['seq_id'], _ps1['comp_id'], _ps2['comp_id'])
                                                   if __s1 != '.' and __s2 != '.' and __s1 != __s2
                                                   and __c1 != '.' and __c2 != '.' and __c1 == __c2)
                                comp_mismatch = any((__s1, __s2) for (__s1, __s2, __c1, __c2)
                                                    in zip(_ps1['seq_id'], _ps2['seq_id'], _ps1['comp_id'], _ps2['comp_id'])
                                                    if __c1 != '.' and __c2 != '.' and __c1 != __c2)

                                if seq_mismatch or comp_mismatch:
                                    if all(__c1 == _ps2['comp_id'][_ps2['seq_id'].index(__s1)] for __s1, __c1
                                           in zip(_ps1['seq_id'], _ps1['comp_id'])
                                           if __s1 != '.' and __c1 != '.' and __s1 in _ps2['seq_id']
                                           and _ps2['comp_id'][_ps2['seq_id'].index(__s1)] != '.'):
                                        if all(__c2 == _ps1['comp_id'][_ps1['seq_id'].index(__s2)] for __s2, __c2
                                               in zip(_ps2['seq_id'], _ps2['comp_id'])
                                               if __s2 != '.' and __c2 != '.' and __s2 in _ps1['seq_id']
                                               and _ps1['comp_id'][_ps1['seq_id'].index(__s2)] != '.'):
                                            continue

                                if 0 < _matched < 4 and unmapped // _matched > 20 and seq_mismatch and len(poly_seq) > 1:
                                    continue

                            if sf_framecode2 not in ref_chain_ids:
                                ref_chain_ids[sf_framecode2] = []

                            if sf_framecode2 not in map_chain_ids:
                                map_chain_ids[sf_framecode2] = {}

                            if sf_framecode2 not in ref_chain_ids or chain_id not in ref_chain_ids[sf_framecode2]:
                                map_chain_ids[sf_framecode2][ps2['chain_id']] = chain_id

                            ref_chain_ids[sf_framecode2].append(chain_id)

                            ref_length = len(ps1['seq_id'])

                            ref_code = getOneLetterCodeCanSequence(_ps1['comp_id'])
                            test_code = getOneLetterCodeCanSequence(_ps2['comp_id'])
                            mid_code = getMiddleCode(ref_code, test_code)
                            ref_gauge_code = getGaugeCode(_ps1['seq_id'])
                            test_gauge_code = getGaugeCode(_ps2['seq_id'])

                            self.__alt_chain |= not alt_chain

                            matched = mid_code.count('|')

                            if self.__reg.tolerant_seq_align and len(poly_seq) > 1:  # and not alt_chain:
                                if 0 < matched < 4 and unmapped // matched > 20:
                                    continue

                            seq_align = {'list_id': _poly_seq_in_lp['list_id'], 'sf_framecode': sf_framecode2, 'chain_id': chain_id,
                                         'length': ref_length, 'matched': matched, 'conflict': conflict, 'unmapped': unmapped,
                                         'sequence_coverage': round(float(length - (unmapped + conflict)) / ref_length, 3),
                                         'ref_seq_id': _ps1['seq_id'], 'test_seq_id': _ps2['seq_id'],
                                         'ref_gauge_code': ref_gauge_code, 'ref_code': ref_code, 'mid_code': mid_code,
                                         'test_code': test_code, 'test_gauge_code': test_gauge_code}

                            if seq_align in seq_align_set:
                                continue

                            seq_align_set.append(seq_align)

                            if not self.__reg.combined_mode and input_source_dic['non_standard_residue'] is None:
                                has_nstd_res = False
                                for j, rc in enumerate(ref_code):
                                    if rc == 'X' and j < len(test_code) and test_code[j] == 'X':
                                        has_nstd_res = True
                                        break

                                if not has_nstd_res:
                                    continue

                                asm = []

                                for _ps in poly_seq:

                                    ent = {'chain_id': _ps['chain_id'], 'seq_id': [], 'comp_id': [],
                                           'chem_comp_name': [], 'exptl_data': []}

                                    for _seq_id, _comp_id in zip(_ps['seq_id'], _ps['comp_id']):

                                        if _comp_id not in STD_MON_DICT:

                                            ent['seq_id'].append(_seq_id)
                                            ent['comp_id'].append(_comp_id)

                                            is_valid, cc_name, cc_rel_status = self.__reg.dpV.getChemCompNameAndStatusOf(_comp_id)

                                            if is_valid:  # matches with comp_id in CCD
                                                if cc_rel_status == 'REL' or cc_name is not None:
                                                    ent['chem_comp_name'].append(cc_name)
                                                else:
                                                    ent['chem_comp_name'].append(f"(Not available due to CCD status code {cc_rel_status})")  # noqa: E501, pylint: disable=line-too-long

                                            else:
                                                ent['chem_comp_name'].append(cc_name)

                                            ent['exptl_data'].append({'coordinate': False})

                                    asm.append(ent)

                                input_source.setItemValue('non_standard_residue', asm)

                            for r_code, t_code, seq_id in zip(ref_code, test_code, ps1['seq_id']):
                                if r_code == 'X' and t_code == 'X':
                                    input_source.updateNonStandardResidueByExptlData(chain_id, seq_id, content_subtype)

                if len(seq_align_set) > 0:
                    self.__reg.report.sequence_alignment.setItemValue(f'nmr_poly_seq_vs_{content_subtype}', seq_align_set)

                if self.__alt_chain:

                    for _poly_seq_in_lp in poly_seq_in_lp[content_subtype]:
                        poly_seq2 = _poly_seq_in_lp['polymer_sequence']
                        sf_framecode2 = _poly_seq_in_lp['sf_framecode']

                        if sf_framecode2 in map_chain_ids:
                            mapping = map_chain_ids[sf_framecode2]

                            total = set(mapping.keys()) | set(mapping.values())

                            k_rests = list(total - set(mapping.keys()))
                            v_rests = list(total - set(mapping.values()))

                            circular = cross = False

                            for k, v in mapping.items():
                                for _k, _v in mapping.items():
                                    if v == _k:
                                        circular = True
                                        break
                                if circular:
                                    break

                            if len(k_rests) == 1 and len(v_rests) == 1:

                                src_chain = k_rests[0]
                                dst_chain = v_rests[0]

                                if circular:
                                    mapping[src_chain] = dst_chain

                                else:

                                    for ps1 in poly_seq:
                                        chain_id = ps1['chain_id']

                                        if chain_id != dst_chain:
                                            continue

                                        for ps2 in poly_seq2:
                                            _ps2 = fillBlankCompIdWithOffset(ps2, 0)

                                            if len(_ps2['seq_id']) > len(ps2['seq_id'])\
                                               and len(_ps2['seq_id']) < len(ps1['seq_id']):
                                                ps2 = _ps2

                                            self.__reg.pA.clear()
                                            self.__reg.pA.setReferenceSequence(ps1['comp_id'], f'REF{chain_id}')
                                            self.__reg.pA.addTestSequence(ps2['comp_id'], chain_id)
                                            self.__reg.pA.doAlign()

                                            myAlign = self.__reg.pA.getAlignment(chain_id)

                                            length = len(myAlign)

                                            if length == 0:
                                                break

                                            _matched, unmapped, conflict, offset_1, offset_2 = getScoreOfSeqAlign(myAlign)

                                            if length == unmapped + conflict or _matched <= conflict + (1 if length > 1 else 0):
                                                break

                                            cross = True
                                            mapping[src_chain] = dst_chain

                                            break

                            for ps1 in poly_seq:
                                chain_id = ps1['chain_id']

                                for ps2 in poly_seq2:

                                    if chain_id != ps2['chain_id']:
                                        continue

                                    _ps2 = fillBlankCompIdWithOffset(ps2, 0)

                                    if len(_ps2['seq_id']) > len(ps2['seq_id']) and len(_ps2['seq_id']) < len(ps1['seq_id']):
                                        ps2 = _ps2

                                    self.__reg.pA.clear()
                                    self.__reg.pA.setReferenceSequence(ps1['comp_id'], f'REF{chain_id}')
                                    self.__reg.pA.addTestSequence(ps2['comp_id'], chain_id)
                                    self.__reg.pA.doAlign()

                                    myAlign = self.__reg.pA.getAlignment(chain_id)

                                    length = len(myAlign)

                                    if length == 0:
                                        continue

                                    _matched, unmapped, conflict, offset_1, offset_2 = getScoreOfSeqAlign(myAlign)

                                    if length == unmapped + conflict or _matched <= conflict + (1 if length > 1 else 0):
                                        continue

                                    _ps1 = ps1 if offset_1 == 0 else fillBlankCompIdWithOffset(ps1, offset_1)
                                    _ps2 = ps2 if offset_2 == 0 else fillBlankCompIdWithOffset(ps2, offset_2)

                                    if conflict > 0 and _ps1['seq_id'][0] < 0 and _ps2['seq_id'][0] < 0:  # noqa: E501, pylint: disable=chained-comparison,line-too-long
                                        continue

                                    if self.__reg.tolerant_seq_align:
                                        seq_mismatch = any((__s1, __s2) for (__s1, __s2, __c1, __c2)
                                                           in zip(_ps1['seq_id'], _ps2['seq_id'], _ps1['comp_id'], _ps2['comp_id'])
                                                           if __s1 != '.' and __s2 != '.' and __s1 != __s2
                                                           and __c1 != '.' and __c2 != '.' and __c1 == __c2)
                                        comp_mismatch = any((__s1, __s2) for (__s1, __s2, __c1, __c2)
                                                            in zip(_ps1['seq_id'], _ps2['seq_id'], _ps1['comp_id'], _ps2['comp_id'])
                                                            if __c1 != '.' and __c2 != '.' and __c1 != __c2)

                                        if seq_mismatch or comp_mismatch:
                                            if all(__c1 == _ps2['comp_id'][_ps2['seq_id'].index(__s1)] for __s1, __c1
                                                   in zip(_ps1['seq_id'], _ps1['comp_id'])
                                                   if __s1 != '.' and __c1 != '.' and __s1 in _ps2['seq_id']
                                                   and _ps2['comp_id'][_ps2['seq_id'].index(__s1)] != '.'):
                                                if all(__c2 == _ps1['comp_id'][_ps1['seq_id'].index(__s2)] for __s2, __c2
                                                       in zip(_ps2['seq_id'], _ps2['comp_id'])
                                                       if __s2 != '.' and __c2 != '.' and __s2 in _ps1['seq_id']
                                                       and _ps1['comp_id'][_ps1['seq_id'].index(__s2)] != '.'):
                                                    continue

                                        if 0 < _matched < 4 and unmapped // _matched > 20 and seq_mismatch and len(poly_seq) > 1:
                                            continue

                                    ref_length = len(ps1['seq_id'])

                                    ref_code = getOneLetterCodeCanSequence(_ps1['comp_id'])
                                    test_code = getOneLetterCodeCanSequence(_ps2['comp_id'])
                                    mid_code = getMiddleCode(ref_code, test_code)
                                    ref_gauge_code = getGaugeCode(_ps1['seq_id'])
                                    test_gauge_code = getGaugeCode(_ps2['seq_id'])

                                    if self.__reg.tolerant_seq_align and (seq_mismatch or comp_mismatch):
                                        if sf_framecode2 in map_seq_ids and chain_id in map_seq_ids[sf_framecode2]:
                                            continue
                                        if _ps2['seq_id'] == list(range(_ps2['seq_id'][0], _ps2['seq_id'][-1] + 1)):
                                            seq_id_conv_dict = {str(__s2): str(__s1) for __s1, __s2
                                                                in zip(_ps1['seq_id'], _ps2['seq_id']) if __s2 != '.'}
                                            if comp_mismatch:
                                                _seq_align =\
                                                    self.__getSeqAlignCodeWithChainRemap__(fileListId, file_type,
                                                                                           content_subtype, sf_framecode2,
                                                                                           chain_id, _ps1, _ps2, myAlign, mapping,
                                                                                           ref_gauge_code, ref_code, mid_code,
                                                                                           test_code, test_gauge_code)
                                                _ps2['seq_id'] = _seq_align['test_seq_id']
                                                if _ps1['seq_id'][0] < 0 and _ps2['seq_id'][0] < 0:
                                                    continue
                                                ref_gauge_code = _seq_align['ref_gauge_code']
                                                ref_code = _seq_align['ref_code']
                                                mid_code = _seq_align['mid_code']
                                                test_code = _seq_align['test_code']
                                                test_gauge_code = _seq_align['test_gauge_code']
                                            else:
                                                chain_id2 = chain_id
                                                if chain_id in mapping.values():
                                                    chain_id2 = next(k for k, v in mapping.items() if v == chain_id)

                                                if sf_framecode2 not in proc_chain_ids:
                                                    proc_chain_ids[sf_framecode2] = set()

                                                if chain_id2 not in proc_chain_ids[sf_framecode2]:
                                                    self.__reg.dpR.fixSeqIdInLoop(fileListId, file_type, content_subtype,
                                                                                  sf_framecode2, chain_id2,
                                                                                  seq_id_conv_dict)
                                                    proc_chain_ids[sf_framecode2].add(chain_id2)

                                                    if 'identical_chain_id' in ps2:
                                                        for chain_id2_ in ps2['identical_chain_id']:
                                                            if chain_id2_ not in proc_chain_ids[sf_framecode2]:
                                                                self.__reg.dpR.fixSeqIdInLoop(fileListId, file_type,
                                                                                              content_subtype, sf_framecode2,
                                                                                              chain_id2_, seq_id_conv_dict)
                                                                proc_chain_ids[sf_framecode2].add(chain_id2_)

                                                _ps2['seq_id'] = _ps1['seq_id']
                                                mid_code = getMiddleCode(ref_code, test_code)
                                                test_gauge_code = ref_gauge_code
                                        else:
                                            if seq_mismatch:
                                                _seq_align =\
                                                    self.__getSeqAlignCodeWithChainRemap__(fileListId, file_type,
                                                                                           content_subtype, sf_framecode2,
                                                                                           chain_id, _ps1, _ps2, myAlign, mapping,
                                                                                           ref_gauge_code, ref_code, mid_code,
                                                                                           test_code, test_gauge_code)
                                                _ps2['seq_id'] = _seq_align['test_seq_id']
                                                if _ps1['seq_id'][0] < 0 and _ps2['seq_id'][0] < 0:
                                                    continue
                                                ref_gauge_code = _seq_align['ref_gauge_code']
                                                ref_code = _seq_align['ref_code']
                                                mid_code = _seq_align['mid_code']
                                                test_code = _seq_align['test_code']
                                                test_gauge_code = _seq_align['test_gauge_code']
                                            else:
                                                _ps2 = fillBlankCompId(_ps1, _ps2)
                                                if _ps1['seq_id'][0] < 0 and _ps2['seq_id'][0] < 0:
                                                    continue
                                                test_code = getOneLetterCodeCanSequence(_ps2['comp_id'])
                                                mid_code = getMiddleCode(ref_code, test_code)
                                                test_gauge_code = ref_gauge_code

                                    matched = mid_code.count('|')

                                    _seq_align = next((_seq_align for _seq_align in seq_align_set
                                                       if _seq_align['list_id'] == _poly_seq_in_lp['list_id']
                                                       and _seq_align['sf_framecode'] == sf_framecode2
                                                       and _seq_align['chain_id'] == chain_id), None)

                                    if _seq_align is not None:
                                        seq_align_set.remove(_seq_align)

                                    seq_align = {'list_id': _poly_seq_in_lp['list_id'],
                                                 'sf_framecode': sf_framecode2, 'chain_id': chain_id,
                                                 'length': ref_length, 'matched': matched,
                                                 'conflict': conflict, 'unmapped': unmapped,
                                                 'sequence_coverage':
                                                 round(float(length - (unmapped + conflict)) / ref_length, 3),
                                                 'ref_seq_id': _ps1['seq_id'], 'test_seq_id': _ps2['seq_id'],
                                                 'ref_gauge_code': ref_gauge_code, 'ref_code': ref_code, 'mid_code': mid_code,
                                                 'test_code': test_code, 'test_gauge_code': test_gauge_code}

                                    if seq_align in seq_align_set:
                                        continue

                                    seq_align_set.append(seq_align)

                            if circular or cross:
                                for k, v in mapping.items():

                                    for ps2 in poly_seq2:

                                        if ps2['chain_id'] != k:
                                            continue

                                        ps2['chain_id'] = f'{v}_'

                                        break

                                    self.__reg.dpR.fixChainIdInLoop(fileListId, file_type, content_subtype, sf_framecode2,
                                                                    k, f'{v}_')

                                for v in mapping.values():

                                    for ps2 in poly_seq2:

                                        if ps2['chain_id'] != f'{v}_':
                                            continue

                                        ps2['chain_id'] = v

                                        break

                                    self.__reg.dpR.fixChainIdInLoop(fileListId, file_type, content_subtype, sf_framecode2,
                                                                    f'{v}_', v)

                            else:
                                for k, v in mapping.items():

                                    for ps2 in poly_seq2:

                                        if ps2['chain_id'] != k:
                                            continue

                                        ps2['chain_id'] = v

                                        break

                                    self.__reg.dpR.fixChainIdInLoop(fileListId, file_type, content_subtype, sf_framecode2,
                                                                    k, v)

        if update_poly_seq:
            self.__extractPolymerSequenceInLoop()
            self.__depositNmrData()

        return is_done

    def __getSeqAlignCodeWithChainRemap__(self, file_list_id: int, file_type: str, content_subtype: str, sf_framecode: str,
                                          chain_id: str, ps1: dict, ps2: dict, myAlign: list, mapping: dict,
                                          ref_gauge_code: str, ref_code: str, mid_code: str, test_code: str, test_gauge_code: str
                                          ) -> dict:
        """ Return human-readable sequence alignment codes with chain remapping if necessary.
        """

        len_ps1 = len(ps1['seq_id'])
        len_ps2 = len(ps2['seq_id'])

        length = len(myAlign)

        seq_id1, seq_id2, comp_id1, comp_id2 = [], [], [], []

        idx1 = idx2 = 0
        for i in range(length):
            myPr0, myPr1 = str(myAlign[i][0]), str(myAlign[i][1])

            if myPr0 != '.':
                while idx1 < len_ps1:
                    if ps1['comp_id'][idx1] == myPr0:
                        seq_id1.append(ps1['seq_id'][idx1])
                        comp_id1.append(myPr0)
                        idx1 += 1
                        break
                    idx1 += 1
            else:
                seq_id1.append(None)
                comp_id1.append('.')

            if myPr1 != '.':
                while idx2 < len_ps2:
                    if ps2['comp_id'][idx2] == myPr1:
                        seq_id2.append(ps2['seq_id'][idx2])
                        comp_id2.append(myPr1)
                        idx2 += 1
                        break
                    idx2 += 1
            else:
                seq_id2.append(None)
                comp_id2.append('.')

        seq_id_conv_dict = {str(_s2): str(_s1) for _s1, _s2
                            in zip(seq_id1, seq_id2) if _s1 is not None and _s2 is not None}
        if ps1['seq_id'] != list(range(ps1['seq_id'][0], ps1['seq_id'][-1] + 1))\
           and not any(True for k in seq_id_conv_dict.keys() if seq_id_conv_dict[k] != k):
            ps2['seq_id'] = ps1['seq_id']
            ref_code = test_code
            mid_code = getMiddleCode(ref_code, test_code)
            ref_gauge_code = test_gauge_code
        else:
            chain_id2 = chain_id
            if mapping is not None and chain_id in mapping.values():
                chain_id2 = next(k for k, v in mapping.items() if v == chain_id)
            self.__reg.dpR.fixSeqIdInLoop(file_list_id, file_type, content_subtype, sf_framecode, chain_id2, seq_id_conv_dict)
            ps2['seq_id'] = ps1['seq_id']
            ref_code = getOneLetterCodeCanSequence(comp_id1)
            test_code = getOneLetterCodeCanSequence(comp_id2)
            mid_code = getMiddleCode(ref_code, test_code)
            ref_gauge_code = getGaugeCode(seq_id1)
            test_gauge_code = ref_gauge_code
            if ' ' in ref_gauge_code:
                for p, g in enumerate(ref_gauge_code):
                    if g == ' ':
                        ref_code = f'{ref_code[0:p]}-{ref_code[p + 1:]}'
            if ' ' in test_gauge_code:
                for p, g in enumerate(test_gauge_code):
                    if g == ' ':
                        test_code = f'{test_code[0:p]}-{test_code[p + 1:]}'

        return {'ref_seq_id': ps1['seq_id'], 'test_seq_id': ps2['seq_id'],
                'ref_gauge_code': ref_gauge_code, 'ref_code': ref_code, 'mid_code': mid_code,
                'test_code': test_code, 'test_gauge_code': test_gauge_code}

    def __validateAtomNomenclature(self) -> bool:
        """ Validate atom nomenclature using NefTranslator and CCD.
        """

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            has_poly_seq_in_lp = has_key_value(input_source_dic, 'polymer_sequence_in_loop')

            if not has_poly_seq_in_lp:
                continue

            poly_seq_in_lp = input_source_dic['polymer_sequence_in_loop']

            content_subtypes = ['poly_seq']
            content_subtypes.extend(poly_seq_in_lp.keys())

            for content_subtype in content_subtypes:

                sf_category = SF_CATEGORIES[file_type][content_subtype]
                lp_category = LP_CATEGORIES[file_type][content_subtype]

                if content_subtype == 'poly_seq':
                    lp_category = AUX_LP_CATEGORIES[file_type][content_subtype][0]

                    if lp_category not in self.__reg.lp_category_list:
                        continue

                if file_type == 'nmr-star' and content_subtype == 'spectral_peak_alt':
                    lp_category = '_Assigned_peak_chem_shift'

                if self.__reg.star_data_type[fileListId] == 'Loop':
                    sf = self.__reg.star_data[fileListId]
                    sf_framecode = ''

                    self.__reg.dpV.validateAtomNomenclature(file_name, file_type, content_subtype,
                                                            sf, sf_framecode, lp_category)  # , first_comp_ids)

                elif self.__reg.star_data_type[fileListId] == 'Saveframe':
                    sf = self.__reg.star_data[fileListId]
                    sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                    self.__reg.dpV.validateAtomNomenclature(file_name, file_type, content_subtype,
                                                            sf, sf_framecode, lp_category)  # , first_comp_ids)

                else:

                    for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                        sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                        if not any(True for loop in sf.loops if loop.category == lp_category):
                            continue

                        self.__reg.dpV.validateAtomNomenclature(file_name, file_type, content_subtype,
                                                                sf, sf_framecode, lp_category)  # , first_comp_ids)

        return not self.__reg.report.isError()

    def __validateAtomTypeOfCsLoop(self) -> bool:
        """ Validate atom type, isotope number on assigned chemical shifts.
        """

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            if input_source_dic['content_subtype'] is None:
                continue

            content_subtype = 'chem_shift'

            if content_subtype not in input_source_dic['content_subtype']:
                continue

            sf_category = SF_CATEGORIES[file_type][content_subtype]
            lp_category = LP_CATEGORIES[file_type][content_subtype]

            if self.__reg.star_data_type[fileListId] == 'Loop':
                sf = self.__reg.star_data[fileListId]
                sf_framecode = ''

                self.__reg.dpV.validateAtomTypeOfCsLoop(file_name, file_type, sf, sf_framecode, lp_category)

            elif self.__reg.star_data_type[fileListId] == 'Saveframe':
                sf = self.__reg.star_data[fileListId]
                sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                self.__reg.dpV.validateAtomTypeOfCsLoop(file_name, file_type, sf, sf_framecode, lp_category)

            else:

                for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                    sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                    if not any(True for loop in sf.loops if loop.category == lp_category):
                        continue

                    self.__reg.dpV.validateAtomTypeOfCsLoop(file_name, file_type, sf, sf_framecode, lp_category)

        return not self.__reg.report.isError()

    def __validateAmbigCodeOfCsLoop(self) -> bool:
        """ Validate ambiguity code on assigned chemical shifts.
        """

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            if input_source_dic['content_subtype'] is None:
                continue

            content_subtype = 'chem_shift'

            if content_subtype not in input_source_dic['content_subtype']:
                continue

            # NEF file has no ambiguity code
            if file_type == 'nef':
                continue

            if not self.__reg.combined_mode or (self.__reg.op == 'nmr-str-replace-cs' and fileListId > 0):
                continue

            sf_category = SF_CATEGORIES[file_type][content_subtype]
            lp_category = LP_CATEGORIES[file_type][content_subtype]

            modified = False

            if self.__reg.star_data_type[fileListId] == 'Loop':
                sf = self.__reg.star_data[fileListId]
                sf_framecode = ''

                modified |= self.__reg.dpV.validateAmbigCodeOfCsLoop(file_name, sf, sf_framecode, lp_category)

            elif self.__reg.star_data_type[fileListId] == 'Saveframe':
                sf = self.__reg.star_data[fileListId]
                sf_framecode = get_first_sf_tag(sf, 'Sf_framecode')

                modified |= self.__reg.dpV.validateAmbigCodeOfCsLoop(file_name, sf, sf_framecode, lp_category)

            else:

                for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                    sf_framecode = get_first_sf_tag(sf, 'Sf_framecode')

                    if not any(True for loop in sf.loops if loop.category == lp_category):
                        continue

                    modified |= self.__reg.dpV.validateAmbigCodeOfCsLoop(file_name, sf, sf_framecode, lp_category)

            if modified:
                self.__depositNmrData()

        return not self.__reg.report.isError()

    def __testIndexConsistency(self) -> bool:
        """ Perform consistency test on index of interesting loops.
        """

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            if input_source_dic['content_type'] == 'nmr-restraints' or input_source_dic['content_subtype'] is None:
                continue

            for content_subtype in input_source_dic['content_subtype']:

                if content_subtype == 'entity':
                    continue

                sf_category = SF_CATEGORIES[file_type][content_subtype]
                lp_category = LP_CATEGORIES[file_type][content_subtype]

                index_tag = INDEX_TAGS[file_type][content_subtype]

                if index_tag is None:
                    continue

                if self.__reg.star_data_type[fileListId] == 'Loop':
                    sf = self.__reg.star_data[fileListId]
                    sf_framecode = ''

                    self.__reg.dpV.testIndexConsistency(file_name, sf, sf_framecode, lp_category, index_tag)

                elif self.__reg.star_data_type[fileListId] == 'Saveframe':
                    sf = self.__reg.star_data[fileListId]
                    sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                    self.__reg.dpV.testIndexConsistency(file_name, sf, sf_framecode, lp_category, index_tag)

                else:

                    for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                        sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                        if not any(True for loop in sf.loops if loop.category == lp_category):
                            continue

                        self.__reg.dpV.testIndexConsistency(file_name, sf, sf_framecode, lp_category, index_tag)

        return not self.__reg.report.isError()

    def __testDataConsistencyInLoop(self, target_content_subtype: str = None) -> bool:
        """ Perform consistency test on data of interesting loops.
        """

        __errors = self.__reg.report.getTotalErrors()

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            if input_source_dic['content_subtype'] is None:
                continue

            modified = False

            for content_subtype in input_source_dic['content_subtype']:

                if content_subtype in ('entry_info', 'entity'):
                    continue

                if target_content_subtype is not None and content_subtype != target_content_subtype:
                    continue

                sf_category = SF_CATEGORIES[file_type][content_subtype]
                lp_category = LP_CATEGORIES[file_type][content_subtype]

                if self.__reg.star_data_type[fileListId] == 'Loop':
                    sf = self.__reg.star_data[fileListId]
                    sf_framecode = ''

                    modified |= self.__reg.dpV.testDataConsistencyInLoop(fileListId, file_name, file_type, content_subtype,
                                                                         sf, sf_framecode, lp_category, 1)

                elif self.__reg.star_data_type[fileListId] == 'Saveframe':
                    sf = self.__reg.star_data[fileListId]
                    sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                    modified |= self.__reg.dpV.testDataConsistencyInLoop(fileListId, file_name, file_type, content_subtype,
                                                                         sf, sf_framecode, lp_category, 1)

                else:

                    parent_pointer = 0

                    for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                        sf_framecode = get_first_sf_tag(sf, 'sf_framecode')
                        parent_pointer += 1

                        if not any(True for loop in sf.loops if loop.category == lp_category):
                            continue

                        modified |= self.__reg.dpV.testDataConsistencyInLoop(fileListId, file_name, file_type, content_subtype,
                                                                             sf, sf_framecode, lp_category, parent_pointer)

            if modified and fileListId == 0:
                self.__depositNmrData()

        return self.__reg.report.getTotalErrors() == __errors

    def __testDataConsistencyInPkLoop(self) -> bool:
        """ Perform consistency test on data of interesting loops.
        """

        fileListId = 0

        input_source = self.__reg.report.input_sources[fileListId]
        input_source_dic = input_source.get()

        file_name = input_source_dic['file_name']
        file_type = input_source_dic['file_type']

        if file_type != 'nmr-star':
            return True

        if input_source_dic['content_subtype'] is None:
            return False

        content_subtype = 'spectral_peak'

        if content_subtype not in input_source_dic['content_subtype']:
            return False

        __errors = self.__reg.report.getTotalErrors()

        self.__reg.lp_data[content_subtype] = []

        sf_category = SF_CATEGORIES[file_type][content_subtype]
        lp_category = LP_CATEGORIES[file_type][content_subtype]

        if self.__reg.star_data_type[fileListId] != 'Entry':
            return False

        parent_pointer = 0

        for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
            sf_framecode = get_first_sf_tag(sf, 'sf_framecode')
            parent_pointer += 1

            if not any(True for loop in sf.loops if loop.category == lp_category):
                continue

            self.__reg.dpV.testDataConsistencyInLoop(fileListId, file_name, file_type, content_subtype,
                                                     sf, sf_framecode, lp_category, parent_pointer)

        return self.__reg.report.getTotalErrors() == __errors

    def __detectConflictDataInLoop(self) -> bool:
        """ Detect redundant/inconsistent data of interesting loops.
        """

        __errors = self.__reg.report.getTotalErrors()

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            if input_source_dic['content_subtype'] is None:
                continue

            for content_subtype in input_source_dic['content_subtype']:

                if content_subtype in ('dist_restraint', 'dihed_restraint', 'rdc_restraint'):

                    sf_category = SF_CATEGORIES[file_type][content_subtype]
                    lp_category = LP_CATEGORIES[file_type][content_subtype]

                    if self.__reg.star_data_type[fileListId] == 'Loop':
                        sf = self.__reg.star_data[fileListId]
                        sf_framecode = ''

                        self.__reg.dpV.detectConflictDataInLoop(file_name, file_type, content_subtype,
                                                                sf, sf_framecode, lp_category)

                    elif self.__reg.star_data_type[fileListId] == 'Saveframe':
                        sf = self.__reg.star_data[fileListId]
                        sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                        self.__reg.dpV.detectConflictDataInLoop(file_name, file_type, content_subtype,
                                                                sf, sf_framecode, lp_category)

                    else:

                        for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                            sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                            if not any(True for loop in sf.loops if loop.category == lp_category):
                                continue

                            self.__reg.dpV.detectConflictDataInLoop(file_name, file_type, content_subtype,
                                                                    sf, sf_framecode, lp_category)

        return self.__reg.report.getTotalErrors() == __errors

    def __testNmrCovalentBond(self) -> bool:
        """ Test the agreement between covalent bonds and experimental NMR data.
        """

        __errors = self.__reg.report.getTotalErrors()

        for fileListId in range(self.__reg.file_path_list_len):

            if fileListId >= len(self.__reg.star_data_type) or self.__reg.star_data_type[fileListId] != 'Entry':
                continue

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            if file_type not in ITEM_NAMES_IN_CS_LOOP:
                continue

            item_names = ITEM_NAMES_IN_CS_LOOP[file_type]
            chain_id_name = item_names['chain_id']
            seq_id_name = item_names['seq_id']
            comp_id_name = item_names['comp_id']
            atom_id_name = item_names['atom_id']
            value_name = item_names['value']

            item_names = ITEM_NAMES_IN_RDC_LOOP[file_type]
            chain_id_1_name = item_names['chain_id_1']
            chain_id_2_name = item_names['chain_id_2']
            seq_id_1_name = item_names['seq_id_1']
            seq_id_2_name = item_names['seq_id_2']
            comp_id_1_name = item_names['comp_id_1']
            comp_id_2_name = item_names['comp_id_2']
            atom_id_1_name = item_names['atom_id_1']
            atom_id_2_name = item_names['atom_id_2']

            auth_chain_id_1_name = 'Auth_asym_ID_1'
            auth_chain_id_2_name = 'Auth_asym_ID_2'
            auth_seq_id_1_name = 'Auth_seq_ID_1'
            auth_seq_id_2_name = 'Auth_seq_ID_2'
            auth_atom_id_1_name = 'Auth_atom_ID_1'
            auth_atom_id_2_name = 'Auth_atom_ID_2'

            # pylint: disable=cell-var-from-loop
            def ext_atom_names(row):
                return (row[chain_id_1_name], row[chain_id_2_name],
                        row[seq_id_1_name], row[seq_id_2_name],
                        row[comp_id_1_name], row[comp_id_2_name],
                        row[atom_id_1_name], row[atom_id_2_name])

            content_subtype = 'chem_shift'

            cs_sf_category = SF_CATEGORIES[file_type][content_subtype]
            cs_lp_category = LP_CATEGORIES[file_type][content_subtype]

            cs_lp_data = None

            for cs_sf in self.__reg.star_data[fileListId].get_saveframes_by_category(cs_sf_category):
                cs_sf_framecode = get_first_sf_tag(cs_sf, 'sf_framecode')

                if not any(True for loop in cs_sf.loops if loop.category == cs_lp_category):
                    continue

                cs_lp_data = next((lp['data'] for lp in self.__reg.lp_data[content_subtype]
                                   if lp['file_name'] == file_name and lp['sf_framecode'] == cs_sf_framecode), None)

                break

            content_subtype = 'poly_seq'

            sf_category = SF_CATEGORIES[file_type][content_subtype]
            lp_category = LP_CATEGORIES[file_type][content_subtype]

            parent_pointer = 0

            for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                sf_framecode = get_first_sf_tag(sf, 'sf_framecode')
                parent_pointer += 1

                for loop in sf.loops:

                    lp_category = loop.category

                    if lp_category is None:
                        continue

                    if (file_type == 'nef' and lp_category == '_nef_covalent_links')\
                       or (file_type == 'nmr-star' and lp_category == '_Bond'):

                        key_items = self.__reg.aux_key_items[file_type][content_subtype][lp_category]
                        data_items = self.__reg.aux_data_items[file_type][content_subtype][lp_category]
                        allowed_tags = AUX_ALLOWED_TAGS[file_type][content_subtype][lp_category]

                        try:

                            aux_data = self.__reg.nefT.check_data(sf, lp_category, key_items, data_items,
                                                                  allowed_tags, None, parent_pointer=parent_pointer,
                                                                  test_on_index=False, enforce_non_zero=False, enforce_sign=False,
                                                                  enforce_range=False, enforce_enum=False,
                                                                  enforce_allowed_tags=(file_type == 'nmr-star'),
                                                                  excl_missing_data=self.__reg.excl_missing_data)[0]

                            disulf_asm, other_asm = [], []

                            for row in aux_data:
                                chain_id_1, chain_id_2, seq_id_1, seq_id_2, \
                                    comp_id_1, comp_id_2, atom_id_1, atom_id_2 = ext_atom_names(row)

                                if atom_id_1 in EMPTY_VALUE or atom_id_2 in EMPTY_VALUE:
                                    continue

                                atom_id_1_ = atom_id_1[0]
                                atom_id_2_ = atom_id_2[0]

                                if atom_id_1_ == 'S' and atom_id_2_ == 'S'\
                                   and not atom_id_1.startswith('SE') and not atom_id_2.startswith('SE'):

                                    disulf = {}
                                    disulf['chain_id_1'] = chain_id_1
                                    disulf['seq_id_1'] = seq_id_1
                                    disulf['comp_id_1'] = comp_id_1
                                    disulf['atom_id_1'] = atom_id_1
                                    disulf['chain_id_2'] = chain_id_2
                                    disulf['seq_id_2'] = seq_id_2
                                    disulf['comp_id_2'] = comp_id_2
                                    disulf['atom_id_2'] = atom_id_2
                                    disulf['distance_value'] = None
                                    bond = self.__reg.dpV.getNmrBondLength(chain_id_1, seq_id_1, atom_id_1,
                                                                           chain_id_2, seq_id_2, atom_id_2)
                                    if bond is not None:
                                        disulf['distance_value'] =\
                                            next((b['distance'] for b in bond
                                                  if b['model_id'] == self.__reg.representative_model_id), None)
                                    if disulf['distance_value'] is None and file_type == 'nmr-star':
                                        cif_chain_id_1 = row[auth_chain_id_1_name]
                                        cif_chain_id_2 = row[auth_chain_id_2_name]
                                        cif_seq_id_1 = row[auth_seq_id_1_name]
                                        cif_seq_id_2 = row[auth_seq_id_2_name]
                                        cif_atom_id_1 = row[auth_atom_id_1_name]
                                        cif_atom_id_2 = row[auth_atom_id_2_name]
                                        bond = self.__reg.dpV.getCoordBondLength(cif_chain_id_1, cif_seq_id_1, cif_atom_id_1,
                                                                                 cif_chain_id_2, cif_seq_id_2, cif_atom_id_2,
                                                                                 label_scheme=False)
                                        if bond is not None:
                                            disulf['distance_value'] =\
                                                next((b['distance'] for b in bond
                                                      if b['model_id'] == self.__reg.representative_model_id), None)
                                    disulf['warning_description_1'] = None
                                    disulf['warning_description_2'] = None

                                    if 'Type' in row and row['Type'] not in EMPTY_VALUE:
                                        disulf['bond_type'] = row['Type']

                                    if cs_lp_data is not None:

                                        ca_chem_shift_1 = cb_chem_shift_1 = None

                                        for _row in cs_lp_data:

                                            atom_id = _row[atom_id_name]

                                            if _row[chain_id_name] == chain_id_1 and _row[seq_id_name] == seq_id_1\
                                               and _row[comp_id_name] in ('CYS', 'DCY'):
                                                if atom_id == 'CA':
                                                    ca_chem_shift_1 = _row[value_name]
                                                elif atom_id == 'CB':
                                                    cb_chem_shift_1 = _row[value_name]

                                            if None in (ca_chem_shift_1, cb_chem_shift_1):
                                                if _row[chain_id_name] == chain_id_1 and _row[seq_id_name] > seq_id_1:
                                                    break
                                            else:
                                                break

                                        disulf['ca_chem_shift_1'] = ca_chem_shift_1
                                        disulf['cb_chem_shift_1'] = cb_chem_shift_1

                                        ca_chem_shift_2 = cb_chem_shift_2 = None

                                        for _row in cs_lp_data:

                                            atom_id = _row[atom_id_name]

                                            if _row[chain_id_name] == chain_id_2 and _row[seq_id_name] == seq_id_2\
                                               and _row[comp_id_name] in ('CYS', 'DCY'):
                                                if atom_id == 'CA':
                                                    ca_chem_shift_2 = _row[value_name]
                                                elif atom_id == 'CB':
                                                    cb_chem_shift_2 = _row[value_name]

                                            if None in (ca_chem_shift_2, cb_chem_shift_2):
                                                if _row[chain_id_name] == chain_id_2 and _row[seq_id_name] > seq_id_2:
                                                    break
                                            else:
                                                break

                                        disulf['ca_chem_shift_2'] = ca_chem_shift_2
                                        disulf['cb_chem_shift_2'] = cb_chem_shift_2

                                        if cb_chem_shift_1 is not None:
                                            if cb_chem_shift_1 < 32.0:
                                                disulf['redox_state_pred_1'] = 'reduced'
                                            elif cb_chem_shift_1 > 35.0:
                                                disulf['redox_state_pred_1'] = 'oxidized'
                                            elif cb_chem_shift_2 is not None:
                                                if cb_chem_shift_2 < 32.0:
                                                    disulf['redox_state_pred_1'] = 'reduced'
                                                elif cb_chem_shift_2 > 35.0:
                                                    disulf['redox_state_pred_1'] = 'oxidized'
                                                else:
                                                    disulf['redox_state_pred_1'] = 'ambiguous'
                                            else:
                                                disulf['redox_state_pred_1'] = 'ambiguous'
                                        else:
                                            disulf['redox_state_pred_1'] = 'unknown'

                                        if cb_chem_shift_2 is not None:
                                            if cb_chem_shift_2 < 32.0:
                                                disulf['redox_state_pred_2'] = 'reduced'
                                            elif cb_chem_shift_2 > 35.0:
                                                disulf['redox_state_pred_2'] = 'oxidized'
                                            elif cb_chem_shift_1 is not None:
                                                if cb_chem_shift_1 < 32.0:
                                                    disulf['redox_state_pred_2'] = 'reduced'
                                                elif cb_chem_shift_1 > 35.0:
                                                    disulf['redox_state_pred_2'] = 'oxidized'
                                                else:
                                                    disulf['redox_state_pred_2'] = 'ambiguous'
                                            else:
                                                disulf['redox_state_pred_2'] = 'ambiguous'
                                        else:
                                            disulf['redox_state_pred_2'] = 'unknown'

                                        if disulf['redox_state_pred_1'] == 'ambiguous'\
                                           and ((ca_chem_shift_1 is not None) or (cb_chem_shift_1 is not None)):
                                            oxi, red = predict_redox_state_of_cysteine(ca_chem_shift_1, cb_chem_shift_1)
                                            disulf['redox_state_pred_1'] = f"oxidized {oxi:.1%}, reduced {red:.1%}"

                                        if disulf['redox_state_pred_2'] == 'ambiguous'\
                                           and ((ca_chem_shift_2 is not None) or (cb_chem_shift_2 is not None)):
                                            oxi, red = predict_redox_state_of_cysteine(ca_chem_shift_2, cb_chem_shift_2)
                                            disulf['redox_state_pred_2'] = f"oxidized {oxi:.1%}, reduced {red:.1%}"

                                        if disulf['redox_state_pred_1'] != 'oxidized' and disulf['redox_state_pred_1'] != 'unknown':

                                            warn = "Disulfide bond "\
                                                f"({chain_id_1}:{seq_id_1}:{comp_id_1} - {chain_id_2}:{seq_id_2}:{comp_id_2}) "\
                                                "can not be verified with the assigned chemical shift values "\
                                                f"({chain_id_1}:{seq_id_1}:{comp_id_1}:CA {ca_chem_shift_1} ppm, "\
                                                f"{chain_id_1}:{seq_id_1}:{comp_id_1}:CB {cb_chem_shift_1} ppm, "\
                                                f"redox_state_pred {disulf['redox_state_pred_1']})."

                                            item = 'anomalous_chemical_shift' if disulf['redox_state_pred_1'] == 'reduced'\
                                                else 'unusual_chemical_shift'

                                            disulf['warning_description_1'] = f'{item}: {warn}'

                                        if disulf['redox_state_pred_2'] != 'oxidized' and disulf['redox_state_pred_2'] != 'unknown':

                                            warn = "Disulfide bond "\
                                                f"({chain_id_1}:{seq_id_1}:{comp_id_1} - {chain_id_2}:{seq_id_2}:{comp_id_2}) "\
                                                "can not be verified with the assigned chemical shift values "\
                                                f"({chain_id_2}:{seq_id_2}:{comp_id_2}:CA {ca_chem_shift_2} ppm, "\
                                                f"{chain_id_2}:{seq_id_2}:{comp_id_2}:CB {cb_chem_shift_2} ppm, "\
                                                f"redox_state_pred {disulf['redox_state_pred_2']})."

                                            item = 'anomalous_chemical_shift' if disulf['redox_state_pred_2'] == 'reduced'\
                                                else 'unusual_chemical_shift'

                                            disulf['warning_description_2'] = f'{item}: {warn}'

                                    disulf_asm.append(disulf)

                                elif chain_id_1 == chain_id_2 and seq_id_1 == 1 and atom_id_1 == 'N'\
                                        and seq_id_2 > 1 and atom_id_2 == 'C':
                                    self.__reg.report.setCyclicPolymer(True)

                                elif atom_id_1 == 'SE' and atom_id_2 == 'SE':  # diselenide bond
                                    pass

                                # hydrogen bond begins

                                elif (atom_id_1_ == 'F' and atom_id_2_ in PROTON_BEGIN_CODE)\
                                        or (atom_id_2_ == 'F' and atom_id_1_ in PROTON_BEGIN_CODE):
                                    pass

                                elif (atom_id_1_ == 'F' and atom_id_2_ == 'F') or (atom_id_2_ == 'F' and atom_id_1_ == 'F'):
                                    pass

                                elif (atom_id_1_ == 'O' and atom_id_2_ in PROTON_BEGIN_CODE)\
                                        or (atom_id_2_ == 'O' and atom_id_1_ in PROTON_BEGIN_CODE):
                                    pass

                                elif (atom_id_1_ == 'O' and atom_id_2_ == 'N') or (atom_id_2_ == 'O' and atom_id_1_ == 'N'):
                                    pass

                                elif (atom_id_1_ == 'O' and atom_id_2_ == 'O') or (atom_id_2_ == 'O' and atom_id_1_ == 'O'):
                                    pass

                                elif (atom_id_1_ == 'N' and atom_id_2_ in PROTON_BEGIN_CODE)\
                                        or (atom_id_2_ == 'N' and atom_id_1_ in PROTON_BEGIN_CODE):
                                    pass

                                elif (atom_id_1_ == 'N' and atom_id_2_ == 'N') or (atom_id_2_ == 'N' and atom_id_1_ == 'N'):
                                    pass

                                # hydrogen bond ends

                                else:

                                    other = {}
                                    other['chain_id_1'] = chain_id_1
                                    other['seq_id_1'] = seq_id_1
                                    other['comp_id_1'] = comp_id_1
                                    other['atom_id_1'] = atom_id_1
                                    other['chain_id_2'] = chain_id_2
                                    other['seq_id_2'] = seq_id_2
                                    other['comp_id_2'] = comp_id_2
                                    other['atom_id_2'] = atom_id_2
                                    other['distance_value'] = None
                                    bond = self.__reg.dpV.getNmrBondLength(chain_id_1, seq_id_1, atom_id_1,
                                                                           chain_id_2, seq_id_2, atom_id_2)
                                    if bond is not None:
                                        other['distance_value'] =\
                                            next((b['distance'] for b in bond
                                                  if b['model_id'] == self.__reg.representative_model_id), None)
                                    if other['distance_value'] is None and file_type == 'nmr-star':
                                        cif_chain_id_1 = row[auth_chain_id_1_name]
                                        cif_chain_id_2 = row[auth_chain_id_2_name]
                                        cif_seq_id_1 = row[auth_seq_id_1_name]
                                        cif_seq_id_2 = row[auth_seq_id_2_name]
                                        cif_atom_id_1 = row[auth_atom_id_1_name]
                                        cif_atom_id_2 = row[auth_atom_id_2_name]
                                        bond = self.__reg.dpV.getCoordBondLength(cif_chain_id_1, cif_seq_id_1, cif_atom_id_1,
                                                                                 cif_chain_id_2, cif_seq_id_2, cif_atom_id_2,
                                                                                 label_scheme=False)
                                        if bond is not None:
                                            other['distance_value'] =\
                                                next((b['distance'] for b in bond
                                                      if b['model_id'] == self.__reg.representative_model_id), None)
                                    other['warning_description_1'] = None
                                    other['warning_description_2'] = None

                                    if 'Type' in row and row['Type'] not in EMPTY_VALUE:
                                        other['bond_type'] = row['Type']

                                    if cs_lp_data is not None:

                                        ca_chem_shift_1 = cb_chem_shift_1 = None

                                        for _row in cs_lp_data:

                                            atom_id = _row[atom_id_name]

                                            if _row[chain_id_name] == chain_id_1 and _row[seq_id_name] == seq_id_1\
                                               and _row[comp_id_name] in ('CYS', 'DCY'):
                                                if atom_id == 'CA':
                                                    ca_chem_shift_1 = _row[value_name]
                                                elif atom_id == 'CB':
                                                    cb_chem_shift_1 = _row[value_name]

                                            if None in (ca_chem_shift_1, cb_chem_shift_1):
                                                if _row[chain_id_name] == chain_id_1 and _row[seq_id_name] > seq_id_1:
                                                    break
                                            else:
                                                break

                                        other['ca_chem_shift_1'] = ca_chem_shift_1
                                        other['cb_chem_shift_1'] = cb_chem_shift_1

                                        ca_chem_shift_2 = cb_chem_shift_2 = None

                                        for _row in cs_lp_data:

                                            atom_id = _row[atom_id_name]

                                            if _row[chain_id_name] == chain_id_2 and _row[seq_id_name] == seq_id_2\
                                               and _row[comp_id_name] in ('CYS', 'DCY'):
                                                if atom_id == 'CA':
                                                    ca_chem_shift_2 = _row[value_name]
                                                elif atom_id == 'CB':
                                                    cb_chem_shift_2 = _row[value_name]

                                            if None in (ca_chem_shift_2, cb_chem_shift_2):
                                                if _row[chain_id_name] == chain_id_2 and _row[seq_id_name] > seq_id_2:
                                                    break
                                            else:
                                                break

                                        other['ca_chem_shift_2'] = ca_chem_shift_2
                                        other['cb_chem_shift_2'] = cb_chem_shift_2

                                        if cb_chem_shift_1 is not None:
                                            if cb_chem_shift_1 < 32.0:
                                                other['redox_state_pred_1'] = 'reduced'
                                            elif cb_chem_shift_1 > 35.0:
                                                other['redox_state_pred_1'] = 'oxidized'
                                            elif cb_chem_shift_2 is not None:
                                                if cb_chem_shift_2 < 32.0:
                                                    other['redox_state_pred_1'] = 'reduced'
                                                elif cb_chem_shift_2 > 35.0:
                                                    other['redox_state_pred_1'] = 'oxidized'
                                                else:
                                                    other['redox_state_pred_1'] = 'ambiguous'
                                            else:
                                                other['redox_state_pred_1'] = 'ambiguous'
                                        else:
                                            other['redox_state_pred_1'] = 'unknown'

                                        if cb_chem_shift_2 is not None:
                                            if cb_chem_shift_2 < 32.0:
                                                other['redox_state_pred_2'] = 'reduced'
                                            elif cb_chem_shift_2 > 35.0:
                                                other['redox_state_pred_2'] = 'oxidized'
                                            elif cb_chem_shift_1 is not None:
                                                if cb_chem_shift_1 < 32.0:
                                                    other['redox_state_pred_2'] = 'reduced'
                                                elif cb_chem_shift_1 > 35.0:
                                                    other['redox_state_pred_2'] = 'oxidized'
                                                else:
                                                    other['redox_state_pred_2'] = 'ambiguous'
                                            else:
                                                other['redox_state_pred_2'] = 'ambiguous'
                                        else:
                                            other['redox_state_pred_2'] = 'unknown'

                                        if other['redox_state_pred_1'] == 'ambiguous'\
                                           and ((ca_chem_shift_1 is not None) or (cb_chem_shift_1 is not None)):
                                            oxi, red = predict_redox_state_of_cysteine(ca_chem_shift_1, cb_chem_shift_1)
                                            other['redox_state_pred_1'] = f"oxidized {oxi:.1%}, reduced {red:.1%}"

                                        if other['redox_state_pred_2'] == 'ambiguous'\
                                           and ((ca_chem_shift_2 is not None) or (cb_chem_shift_2 is not None)):
                                            oxi, red = predict_redox_state_of_cysteine(ca_chem_shift_2, cb_chem_shift_2)
                                            other['redox_state_pred_2'] = f"oxidized {oxi:.1%}, reduced {red:.1%}"

                                        if other['redox_state_pred_1'] != 'oxidized' and other['redox_state_pred_1'] != 'unknown':

                                            warn = "Other bond "\
                                                f"({chain_id_1}:{seq_id_1}:{comp_id_1} - {chain_id_2}:{seq_id_2}:{comp_id_2}) "\
                                                "can not be verified with the assigned chemical shift values "\
                                                f"({chain_id_1}:{seq_id_1}:{comp_id_1}:CA {ca_chem_shift_1} ppm, "\
                                                f"{chain_id_1}:{seq_id_1}:{comp_id_1}:CB {cb_chem_shift_1} ppm, "\
                                                f"redox_state_pred {other['redox_state_pred_1']})."

                                            item = 'anomalous_chemical_shift' if other['redox_state_pred_1'] == 'reduced'\
                                                else 'unusual_chemical_shift'

                                            other['warning_description_1'] = f'{item}: {warn}'

                                        if other['redox_state_pred_2'] != 'oxidized' and other['redox_state_pred_2'] != 'unknown':

                                            warn = "Other bond "\
                                                f"({chain_id_1}:{seq_id_1}:{comp_id_1} - {chain_id_2}:{seq_id_2}:{comp_id_2}) "\
                                                "can not be verified with the assigned chemical shift values "\
                                                f"({chain_id_2}:{seq_id_2}:{comp_id_2}:CA {ca_chem_shift_2} ppm, "\
                                                f"{chain_id_2}:{seq_id_2}:{comp_id_2}:CB {cb_chem_shift_2} ppm, "\
                                                f"redox_state_pred {other['redox_state_pred_2']})."

                                            item = 'anomalous_chemical_shift' if other['redox_state_pred_2'] == 'reduced'\
                                                else 'unusual_chemical_shift'

                                            other['warning_description_2'] = f'{item}: {warn}'

                                    other_asm.append(other)

                            if len(disulf_asm) > 0:
                                input_source.setItemValue('disulfide_bond', disulf_asm)

                                self.__reg.report.setDisulfideBond(True)

                            if len(other_asm) > 0:
                                input_source.setItemValue('other_bond', other_asm)

                                self.__reg.report.setOtherBond(True)

                        except KeyError:
                            # """
                            # self.__reg.report.error.appendDescription('multiple_data',
                            #                                           {'file_name': file_name, 'sf_framecode': sf_framecode,
                            #                                            'category': lp_category, 'description': str(e).strip("'")})
                            #
                            # if self.__reg.verbose:
                            #     self.__reg.log.write(f"+{self.__class_name__}.__testNmrCovalentBond() "
                            #                          f"++ KeyError  - {str(e)}\n")
                            # """
                            pass
                        except LookupError as e:

                            item = 'format_issue' if 'Unauthorized' in str(e) else 'missing_mandatory_item'

                            self.__reg.report.error.appendDescription(item,
                                                                      {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                       'category': lp_category, 'description': str(e).strip("'")})

                            self.__reg.log.write(f"+{self.__class_name__}.__testNmrCovalentBond() "
                                                 f"++ LookupError  - {file_name} {sf_framecode} {lp_category} {str(e)}\n")

                        except ValueError as e:

                            self.__reg.report.error.appendDescription('invalid_data',
                                                                      {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                       'category': lp_category, 'description': str(e).strip("'")})

                            if self.__reg.verbose:
                                self.__reg.log.write(f"+{self.__class_name__}.__testNmrCovalentBond() "
                                                     f"++ ValueError  - {str(e)}\n")

                        except UserWarning:  # pylint: disable=broad-exception-caught
                            pass

                        except Exception as e:  # pylint: disable=broad-exception-caught

                            self.__reg.report.error.appendDescription('internal_error',
                                                                      f"+{self.__class_name__}.__testNmrCovalentBond() "
                                                                      "++ Error  - " + str(e))

                            if self.__reg.verbose:
                                self.__reg.log.write(f"+{self.__class_name__}.__testNmrCovalentBond() "
                                                     f"++ Error  - {str(e)}\n")

        return self.__reg.report.getTotalErrors() == __errors

    def __testDataConsistencyInAuxLoop(self) -> bool:
        """ Perform consistency test on data of auxiliary loops.
        """

        __errors = self.__reg.report.getTotalErrors()

        for fileListId in range(self.__reg.file_path_list_len):

            if fileListId >= len(self.__reg.star_data_type) or self.__reg.star_data_type[fileListId] != 'Entry':
                continue

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            if input_source_dic['content_subtype'] is None:
                continue

            modified = False

            for content_subtype in input_source_dic['content_subtype']:

                if content_subtype == 'entity':
                    continue

                if self.__reg.op == 'nmr-str-replace-cs' and fileListId > 0 and content_subtype == 'spectral_peak':
                    continue

                sf_category = SF_CATEGORIES[file_type][content_subtype]
                lp_category = LP_CATEGORIES[file_type][content_subtype]

                parent_pointer = 0

                for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                    sf_framecode = get_first_sf_tag(sf, 'sf_framecode')
                    parent_pointer += 1

                    if content_subtype.startswith('spectral_peak'):

                        try:

                            _num_dim = get_first_sf_tag(sf, NUM_DIM_ITEMS[file_type])
                            num_dim = int(_num_dim)

                            if num_dim not in range(1, MAX_DIM_NUM_OF_SPECTRA):
                                raise ValueError()

                        except ValueError:

                            err = f"{NUM_DIM_ITEMS[file_type]} {str(_num_dim)!r} "\
                                f"must be in {set(range(1, MAX_DIM_NUM_OF_SPECTRA))}."

                            self.__reg.report.error.appendDescription('invalid_data',
                                                                      {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                       'category': lp_category, 'description': err})

                            if self.__reg.verbose:
                                self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInAuxLoop() "
                                                     f"++ ValueError  - {err}\n")

                            continue

                    for loop in sf.loops:

                        lp_category = loop.category

                        if lp_category is None:
                            continue

                        # main content of loop has been processed in __testDataConsistencyInLoop()
                        if lp_category in LP_CATEGORIES[file_type][content_subtype]:
                            continue

                        if AUX_LP_CATEGORIES[file_type][content_subtype] is None:
                            continue

                        if lp_category in AUX_LP_CATEGORIES[file_type][content_subtype]:

                            key_items = self.__reg.aux_key_items[file_type][content_subtype][lp_category]
                            data_items = self.__reg.aux_data_items[file_type][content_subtype][lp_category]
                            allowed_tags = AUX_ALLOWED_TAGS[file_type][content_subtype][lp_category]

                            try:

                                aux_data = self.__reg.nefT.check_data(sf, lp_category, key_items, data_items,
                                                                      allowed_tags, None, parent_pointer=parent_pointer,
                                                                      test_on_index=True, enforce_non_zero=True, enforce_sign=True,
                                                                      enforce_range=True, enforce_enum=True,
                                                                      enforce_allowed_tags=(file_type == 'nmr-star'),
                                                                      excl_missing_data=self.__reg.excl_missing_data)[0]

                                self.__reg.aux_data[content_subtype].append({'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                             'category': lp_category, 'data': aux_data})

                                if content_subtype == 'spectral_peak':
                                    self.__reg.dpV.testDataConsistencyInAuxLoopOfSpectralPeak(file_name, file_type, sf_framecode,
                                                                                              num_dim, lp_category, aux_data)
                                if file_type == 'nmr-star' and content_subtype == 'spectral_peak_alt':
                                    self.__reg.dpV.testDataConsistencyInAuxLoopOfSpectralPeakAlt(file_name, file_type, sf_framecode,
                                                                                                 num_dim, lp_category, aux_data,
                                                                                                 sf, parent_pointer)

                            except KeyError as e:

                                self.__reg.report.error.appendDescription('multiple_data',
                                                                          {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                           'category': lp_category,
                                                                           'description': str(e).strip("'")})

                                if self.__reg.verbose:
                                    self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInAuxLoop() "
                                                         f"++ KeyError  - {str(e)}\n")

                            except LookupError as e:

                                item = 'format_issue' if 'Unauthorized' in str(e) else 'missing_mandatory_item'

                                self.__reg.report.error.appendDescription(item,
                                                                          {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                           'category': lp_category,
                                                                           'description': str(e).strip("'")})

                                self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInAuxLoop() "
                                                     f"++ LookupError  - {file_name} {sf_framecode} {lp_category} {str(e)}\n")

                            except ValueError as e:

                                self.__reg.report.error.appendDescription('invalid_data',
                                                                          {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                           'category': lp_category,
                                                                           'description': str(e).strip("'")})

                                if self.__reg.verbose:
                                    self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInAuxLoop() "
                                                         f"++ ValueError  - {str(e)}\n")

                            except UserWarning as e:

                                warns = str(e).strip("'").split('\n')

                                has_multiple_data = has_bad_pattern = False

                                for warn in warns:

                                    if len(warn) == 0:
                                        continue

                                    zero = warn.startswith('[Zero value error]')
                                    nega = warn.startswith('[Negative value error]')
                                    rang = warn.startswith('[Range value error]')
                                    enum = warn.startswith('[Enumeration error]')
                                    mult = warn.startswith('[Multiple data]')
                                    remo = warn.startswith('[Remove bad pattern]')
                                    clea = warn.startswith('[Clear bad pattern]')

                                    if zero or nega or rang or enum or mult or remo or clea:

                                        p = warn.index(']') + 2
                                        warn = warn[p:]

                                        if zero or nega or rang:
                                            item = 'unusual_data'
                                        elif enum:
                                            item = 'enum_mismatch'
                                        elif remo:
                                            if content_subtype == 'chem_shift':
                                                warn += ' Your unassigned chemical shifts have been removed.'
                                                item = 'incompletely_assigned_chemical_shift'
                                            else:
                                                item = 'insufficient_data'
                                            has_bad_pattern = True
                                        elif clea:
                                            if content_subtype.startswith('spectral_peak'):
                                                warn += ' Unassigned spectral peaks can be included in your peak list(s).'
                                                item = 'incompletely_assigned_spectral_peak'
                                            else:
                                                item = 'insufficient_data'
                                        elif self.__reg.resolve_conflict:
                                            # item = 'redundant_data'
                                            has_multiple_data = True
                                            continue
                                        else:
                                            item = 'multiple_data'

                                        if zero or nega or rang or enum or remo or clea or self.__reg.resolve_conflict:

                                            self.__reg.report.warning.appendDescription(item,
                                                                                        {'file_name': file_name,
                                                                                         'sf_framecode': sf_framecode,
                                                                                         'category': lp_category,
                                                                                         'description': warn})

                                            if self.__reg.verbose:
                                                self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInAuxLoop() "
                                                                     f"++ Warning  - {warn}\n")

                                        else:

                                            self.__reg.report.error.appendDescription(item,
                                                                                      {'file_name': file_name,
                                                                                       'sf_framecode': sf_framecode,
                                                                                       'category': lp_category,
                                                                                       'description': warn})

                                            if self.__reg.verbose:
                                                self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInAuxLoop() "
                                                                     f"++ KeyError  - {warn}\n")

                                    else:

                                        self.__reg.report.error.appendDescription('internal_error',
                                                                                  f"+{self.__class_name__}.__testDataConsistencyInAuxLoop() "  # noqa: E501, pylint: disable=line-too-long
                                                                                  "++ Error  - " + warn)

                                        if self.__reg.verbose:
                                            self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInAuxLoop() "
                                                                 f"++ Error  - {warn}\n")

                                # try to parse data without constraints
                                if has_multiple_data:
                                    conflict_id = self.__reg.nefT.get_conflict_id(sf, lp_category, key_items)[0]

                                    if len(conflict_id) > 0:
                                        _loop = sf.get_loop(lp_category)

                                        for lcid in conflict_id:
                                            del _loop.data[lcid]

                                        index_tag = INDEX_TAGS[file_type][content_subtype]
                                        if index_tag is not None:
                                            index_col = loop.tags.index(index_tag) if index_tag in loop.tags else -1
                                            if index_col != -1:
                                                for idx, row in enumerate(loop, start=1):
                                                    row[index_col] = idx

                                        modified = True

                                # try to parse data without bad patterns
                                if has_bad_pattern:
                                    conflict_id = self.__reg.nefT.get_bad_pattern_id(sf, lp_category, key_items, data_items)[0]

                                    if len(conflict_id) > 0:
                                        _loop = sf.get_loop(lp_category)

                                        for lcid in conflict_id:
                                            del _loop.data[lcid]

                                        modified = True

                                # collect unresolved redundant data
                                if has_multiple_data:

                                    try:

                                        aux_data = self.__reg.nefT.check_data(sf, lp_category, key_items, data_items,
                                                                              allowed_tags, None, parent_pointer=parent_pointer,
                                                                              test_on_index=True,  # important
                                                                              enforce_allowed_tags=(file_type == 'nmr-star'),
                                                                              excl_missing_data=self.__reg.excl_missing_data)[0]

                                    except UserWarning as e2:

                                        warns = str(e2).strip("'").split('\n')

                                        for warn in warns:

                                            if len(warn) == 0 or not warn.startswith('[Multiple data]'):
                                                continue

                                            p = warn.index(']') + 2
                                            warn = warn[p:]

                                            self.__reg.report.warning.appendDescription('redundant_data',
                                                                                        {'file_name': file_name,
                                                                                         'sf_framecode': sf_framecode,
                                                                                         'category': lp_category,
                                                                                         'description': warn})

                                            if self.__reg.verbose:
                                                self.__reg.log.write(f"+{self.__class_name__}.testDataConsistencyInLoop() "
                                                                     f"++ Warning  - {warn}\n")

                                    except Exception:  # pylint: disable=broad-exception-caught
                                        pass

                                try:

                                    aux_data = self.__reg.nefT.check_data(sf, lp_category, key_items, data_items,
                                                                          allowed_tags, None, parent_pointer=parent_pointer,
                                                                          enforce_allowed_tags=(file_type == 'nmr-star'),
                                                                          excl_missing_data=self.__reg.excl_missing_data)[0]

                                    self.__reg.aux_data[content_subtype].append({'file_name': file_name,
                                                                                 'sf_framecode': sf_framecode,
                                                                                 'category': lp_category, 'data': aux_data})

                                    if content_subtype == 'spectral_peak':
                                        self.__reg.dpV.testDataConsistencyInAuxLoopOfSpectralPeak(file_name, file_type,
                                                                                                  sf_framecode,
                                                                                                  num_dim, lp_category, aux_data)
                                    if file_type == 'nmr-star' and content_subtype == 'spectral_peak_alt':
                                        self.__reg.dpV.testDataConsistencyInAuxLoopOfSpectralPeakAlt(file_name, file_type,
                                                                                                     sf_framecode,
                                                                                                     num_dim, lp_category, aux_data,
                                                                                                     sf, parent_pointer)

                                except Exception:  # pylint: disable=broad-exception-caught
                                    pass

                            except Exception as e:  # pylint: disable=broad-exception-caught

                                self.__reg.report.error.appendDescription('internal_error',
                                                                          f"+{self.__class_name__}.__testDataConsistencyInAuxLoop() "  # noqa: E501, pylint: disable=line-too-long
                                                                          "++ Error  - " + str(e))

                                if self.__reg.verbose:
                                    self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInAuxLoop() "
                                                         f"++ Error  - {str(e)}\n")

                        elif lp_category in LINKED_LP_CATEGORIES[file_type][content_subtype]:

                            if not self.__reg.bmrb_only:

                                warn = f"Ignored {lp_category!r} loop in {sf_framecode!r} saveframe."

                                self.__reg.report.warning.appendDescription('skipped_loop_category',
                                                                            {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                             'category': lp_category, 'description': warn})

                                if self.__reg.verbose:
                                    self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInAuxLoop() "
                                                         f"++ Warning  - {warn}\n")

                        else:

                            if not self.__reg.bmrb_only:

                                if file_type == 'nef':
                                    warn = f"Ignored third party software's loop {lp_category!r} in {sf_framecode!r} saveframe."
                                else:
                                    warn = f"Ignored {lp_category!r} loop in {sf_framecode!r} saveframe."

                                self.__reg.report.warning.appendDescription('skipped_loop_category',
                                                                            {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                             'category': lp_category, 'description': warn})

                                if self.__reg.verbose:
                                    self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInAuxLoop() "
                                                         f"++ Warning  - {warn}\n")

            if modified and fileListId == 0:
                self.__depositNmrData()

        return self.__reg.report.getTotalErrors() == __errors

    def __testDataConsistencyInPkAuxLoop(self) -> bool:
        """ Perform consistency test on data of auxiliary loops.
        """

        fileListId = 0

        input_source = self.__reg.report.input_sources[fileListId]
        input_source_dic = input_source.get()

        file_name = input_source_dic['file_name']
        file_type = input_source_dic['file_type']

        if file_type != 'nmr-star':
            return True

        if input_source_dic['content_subtype'] is None:
            return False

        content_subtype = 'spectral_peak'

        if content_subtype not in input_source_dic['content_subtype']:
            return False

        # self.__reg.aux_data[content_subtype].clear()

        __errors = self.__reg.report.getTotalErrors()

        sf_category = SF_CATEGORIES[file_type][content_subtype]
        lp_category = LP_CATEGORIES[file_type][content_subtype]

        parent_pointer = 0

        modified = False

        for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
            sf_framecode = get_first_sf_tag(sf, 'sf_framecode')
            parent_pointer += 1

            if content_subtype.startswith('spectral_peak'):

                try:

                    _num_dim = get_first_sf_tag(sf, NUM_DIM_ITEMS[file_type])
                    num_dim = int(_num_dim)

                    if num_dim not in range(1, MAX_DIM_NUM_OF_SPECTRA):
                        raise ValueError()

                except ValueError:

                    err = f"{NUM_DIM_ITEMS[file_type]} {str(_num_dim)!r} must be in {set(range(1, MAX_DIM_NUM_OF_SPECTRA))}."

                    self.__reg.report.error.appendDescription('invalid_data',
                                                              {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                               'category': lp_category, 'description': err})

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInPkAuxLoop() "
                                             f"++ ValueError  - {err}\n")

                    continue

            for loop in sf.loops:

                lp_category = loop.category

                if lp_category is None:
                    continue

                # main content of loop has been processed in testDataConsistencyInLoop()
                if lp_category in LP_CATEGORIES[file_type][content_subtype]:
                    continue

                if AUX_LP_CATEGORIES[file_type][content_subtype] is None:
                    continue

                if lp_category in AUX_LP_CATEGORIES[file_type][content_subtype]:

                    key_items = self.__reg.aux_key_items[file_type][content_subtype][lp_category]
                    data_items = self.__reg.aux_data_items[file_type][content_subtype][lp_category]
                    allowed_tags = AUX_ALLOWED_TAGS[file_type][content_subtype][lp_category]

                    try:

                        aux_data = self.__reg.nefT.check_data(sf, lp_category, key_items, data_items,
                                                              allowed_tags, None, parent_pointer=parent_pointer,
                                                              test_on_index=True, enforce_non_zero=True, enforce_sign=True,
                                                              enforce_range=True, enforce_enum=True,
                                                              enforce_allowed_tags=(file_type == 'nmr-star'),
                                                              excl_missing_data=self.__reg.excl_missing_data)[0]

                        self.__reg.aux_data[content_subtype].append({'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                     'category': lp_category, 'data': aux_data})

                        if content_subtype == 'spectral_peak':
                            self.__reg.dpV.testDataConsistencyInAuxLoopOfSpectralPeak(file_name, file_type, sf_framecode,
                                                                                      num_dim, lp_category, aux_data)
                        if file_type == 'nmr-star' and content_subtype == 'spectral_peak_alt':
                            self.__reg.dpV.testDataConsistencyInAuxLoopOfSpectralPeakAlt(file_name, file_type, sf_framecode,
                                                                                         num_dim, lp_category, aux_data,
                                                                                         sf, parent_pointer)

                    except KeyError as e:

                        self.__reg.report.error.appendDescription('multiple_data',
                                                                  {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                   'category': lp_category, 'description': str(e).strip("'")})

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInPkAuxLoop() "
                                                 f"++ KeyError  - {str(e)}\n")

                    except LookupError as e:

                        item = 'format_issue' if 'Unauthorized' in str(e) else 'missing_mandatory_item'

                        self.__reg.report.error.appendDescription(item,
                                                                  {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                   'category': lp_category, 'description': str(e).strip("'")})

                        self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInPkAuxLoop() "
                                             f"++ LookupError  - {file_name} {sf_framecode} {lp_category} {str(e)}\n")

                    except ValueError as e:

                        self.__reg.report.error.appendDescription('invalid_data',
                                                                  {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                   'category': lp_category, 'description': str(e).strip("'")})

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInPkAuxLoop() "
                                                 f"++ ValueError  - {str(e)}\n")

                    except UserWarning as e:

                        warns = str(e).strip("'").split('\n')

                        has_multiple_data = has_bad_pattern = False

                        for warn in warns:

                            if len(warn) == 0:
                                continue

                            zero = warn.startswith('[Zero value error]')
                            nega = warn.startswith('[Negative value error]')
                            rang = warn.startswith('[Range value error]')
                            enum = warn.startswith('[Enumeration error]')
                            mult = warn.startswith('[Multiple data]')
                            remo = warn.startswith('[Remove bad pattern]')
                            clea = warn.startswith('[Clear bad pattern]')

                            if zero or nega or rang or enum or mult or remo or clea:

                                p = warn.index(']') + 2
                                warn = warn[p:]

                                if zero or nega or rang:
                                    item = 'unusual_data'
                                elif enum:
                                    item = 'enum_mismatch'
                                elif remo:
                                    if content_subtype == 'chem_shift':
                                        warn += ' Your unassigned chemical shifts have been removed.'
                                        item = 'incompletely_assigned_chemical_shift'
                                    else:
                                        item = 'insufficient_data'
                                    has_bad_pattern = True
                                elif clea:
                                    if content_subtype.startswith('spectral_peak'):
                                        warn += ' Unassigned spectral peaks can be included in your peak list(s).'
                                        item = 'incompletely_assigned_spectral_peak'
                                    else:
                                        item = 'insufficient_data'
                                elif self.__reg.resolve_conflict:
                                    # item = 'redundant_data'
                                    has_multiple_data = True
                                    continue
                                else:
                                    item = 'multiple_data'

                                if zero or nega or rang or enum or remo or clea or self.__reg.resolve_conflict:

                                    self.__reg.report.warning.appendDescription(item,
                                                                                {'file_name': file_name,
                                                                                 'sf_framecode': sf_framecode,
                                                                                 'category': lp_category, 'description': warn})

                                    if self.__reg.verbose:
                                        self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInPkAuxLoop() "
                                                             f"++ Warning  - {warn}\n")

                                else:

                                    self.__reg.report.error.appendDescription(item,
                                                                              {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                               'category': lp_category, 'description': warn})

                                    if self.__reg.verbose:
                                        self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInPkAuxLoop() "
                                                             f"++ KeyError  - {warn}\n")

                            else:

                                self.__reg.report.error.appendDescription('internal_error',
                                                                          f"+{self.__class_name__}.__testDataConsistencyInPkAuxLoop() "  # noqa: E501, pylint: disable=line-too-long
                                                                          "++ Error  - " + warn)

                                if self.__reg.verbose:
                                    self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInPkAuxLoop() "
                                                         f"++ Error  - {warn}\n")

                        # try to parse data without constraints
                        if has_multiple_data:
                            conflict_id = self.__reg.nefT.get_conflict_id(sf, lp_category, key_items)[0]

                            if len(conflict_id) > 0:
                                _loop = sf.get_loop(lp_category)

                                for lcid in conflict_id:
                                    del _loop.data[lcid]

                                index_tag = INDEX_TAGS[file_type][content_subtype]
                                if index_tag is not None:
                                    index_col = loop.tags.index(index_tag) if index_tag in loop.tags else -1
                                    if index_col != -1:
                                        for idx, row in enumerate(loop, start=1):
                                            row[index_col] = idx

                                modified = True

                        # try to parse data without bad patterns
                        if has_bad_pattern:
                            conflict_id = self.__reg.nefT.get_bad_pattern_id(sf, lp_category, key_items, data_items)[0]

                            if len(conflict_id) > 0:
                                _loop = sf.get_loop(lp_category)

                                for lcid in conflict_id:
                                    del _loop.data[lcid]

                                modified = True

                        # collect unresolved redundant data
                        if has_multiple_data:

                            try:

                                aux_data = self.__reg.nefT.check_data(sf, lp_category, key_items, data_items,
                                                                      allowed_tags, None, parent_pointer=parent_pointer,
                                                                      test_on_index=True,  # important
                                                                      enforce_allowed_tags=(file_type == 'nmr-star'),
                                                                      excl_missing_data=self.__reg.excl_missing_data)[0]

                            except UserWarning as e2:

                                warns = str(e2).strip("'").split('\n')

                                for warn in warns:

                                    if len(warn) == 0 or not warn.startswith('[Multiple data]'):
                                        continue

                                    p = warn.index(']') + 2
                                    warn = warn[p:]

                                    self.__reg.report.warning.appendDescription("redundant_data",
                                                                                {'file_name': file_name,
                                                                                 'sf_framecode': sf_framecode,
                                                                                 'category': lp_category, 'description': warn})

                                    if self.__reg.verbose:
                                        self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInPkAuxLoop() "
                                                             f"++ Warning  - {warn}\n")

                            except Exception:  # pylint: disable=broad-exception-caught
                                pass

                        try:

                            aux_data = self.__reg.nefT.check_data(sf, lp_category, key_items, data_items,
                                                                  allowed_tags, None, parent_pointer=parent_pointer,
                                                                  enforce_allowed_tags=(file_type == 'nmr-star'),
                                                                  excl_missing_data=self.__reg.excl_missing_data)[0]

                            self.__reg.aux_data[content_subtype].append({'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                         'category': lp_category, 'data': aux_data})

                            if content_subtype == 'spectral_peak':
                                self.__reg.dpV.testDataConsistencyInAuxLoopOfSpectralPeak(file_name, file_type, sf_framecode,
                                                                                          num_dim, lp_category, aux_data)
                            if file_type == 'nmr-star' and content_subtype == 'spectral_peak_alt':
                                self.__reg.dpV.testDataConsistencyInAuxLoopOfSpectralPeakAlt(file_name, file_type, sf_framecode,
                                                                                             num_dim, lp_category, aux_data,
                                                                                             sf, parent_pointer)

                        except Exception:  # pylint: disable=broad-exception-caught
                            pass

                    except Exception as e:  # pylint: disable=broad-exception-caught

                        self.__reg.report.error.appendDescription('internal_error',
                                                                  f"+{self.__class_name__}.__testDataConsistencyInPkAuxLoop() "
                                                                  "++ Error  - " + str(e))

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInPkAuxLoop() "
                                                 f"++ Error  - {str(e)}\n")

                elif lp_category in LINKED_LP_CATEGORIES[file_type][content_subtype]:

                    if not self.__reg.bmrb_only:

                        warn = f"Ignored {lp_category!r} loop in {sf_framecode!r} saveframe."

                        self.__reg.report.warning.appendDescription('skipped_loop_category',
                                                                    {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                     'category': lp_category, 'description': warn})

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInPkAuxLoop() "
                                                 f"++ Warning  - {warn}\n")

                else:

                    if not self.__reg.bmrb_only:

                        if file_type == 'nef':
                            warn = f"Ignored third party software's loop {lp_category!r} in {sf_framecode!r} saveframe."
                        else:
                            warn = f"Ignored {lp_category!r} loop in {sf_framecode!r} saveframe."

                        self.__reg.report.warning.appendDescription('skipped_loop_category',
                                                                    {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                     'category': lp_category, 'description': warn})

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__testDataConsistencyInPkAuxLoop() "
                                                 f"++ Warning  - {warn}\n")

        if modified:
            self.__depositNmrData()

        return self.__reg.report.getTotalErrors() == __errors

    def __testSfTagConsistency(self) -> bool:
        """ Perform consistency test on saveframe tags.
        """

        __errors = self.__reg.report.getTotalErrors()

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            if input_source_dic['content_subtype'] is None:
                continue

            if self.__reg.star_data_type[fileListId] != 'Entry':
                continue

            for content_subtype in input_source_dic['content_subtype']:

                if self.__reg.op == 'nmr-cs-mr-merge' or content_subtype == 'entity':
                    continue

                sf_category = SF_CATEGORIES[file_type][content_subtype]

                parent_keys = set()
                sf_framecode_dict = {}

                list_id = 1  # tentative parent key if not exists

                for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                    sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                    if self.__reg.combined_mode and sf.tag_prefix != SF_TAG_PREFIXES[file_type][content_subtype]:

                        err = f"Saveframe tag prefix {sf.tag_prefix!r} did not match with "\
                            f"{SF_TAG_PREFIXES[file_type][content_subtype]!r} in {sf_framecode!r} saveframe."

                        self.__reg.report.error.appendDescription('internal_error',
                                                                  f"+{self.__class_name__}.__testSfTagConsistency() "
                                                                  "++ Error  - " + err)

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__testSfTagConsistency() "
                                                 f"++ Error  - {err}\n")

                    try:

                        sf_tag_items = copy.copy(SF_TAG_ITEMS[file_type][content_subtype])

                        if not self.__reg.combined_mode:
                            for sf_tag_item in sf_tag_items:
                                if sf_tag_item['name'] == 'sf_framecode' if file_type == 'nef' else 'Sf_framecode':
                                    sf_tag_item['mandatory'] = False

                        sf_tag_data = self.__reg.nefT.check_sf_tag(sf, file_type, sf_category, sf_tag_items,
                                                                   SF_ALLOWED_TAGS[file_type][content_subtype],
                                                                   enforce_non_zero=True, enforce_sign=True,
                                                                   enforce_range=True, enforce_enum=True)

                        self.__reg.dpV.testParentChildRelation(file_name, file_type, content_subtype,
                                                               parent_keys, list_id, sf_framecode, sf_framecode_dict, sf_tag_data)

                        self.__reg.sf_tag_data[content_subtype].append({'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                        'data': sf_tag_data})

                    except LookupError as e:

                        self.__reg.report.error.appendDescription('missing_mandatory_item',
                                                                  {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                   'description': str(e).strip("'")})

                        self.__reg.log.write(f"+{self.__class_name__}.__testSfTagConsistency() "
                                             f"++ LookupError  - {file_name} {sf_framecode} {str(e)}\n")

                    except ValueError as e:

                        self.__reg.report.error.appendDescription('invalid_data',
                                                                  {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                   'description': str(e).strip("'")})

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__testSfTagConsistency() "
                                                 f"++ ValueError  - {str(e)}\n")

                    except UserWarning as e:

                        warns = str(e).strip("'").split('\n')

                        for warn in warns:

                            if len(warn) == 0:
                                continue

                            zero = warn.startswith('[Zero value error]')
                            nega = warn.startswith('[Negative value error]')
                            rang = warn.startswith('[Range value error]')
                            enum = warn.startswith('[Enumeration error]')

                            ignorable = False

                            if zero or nega or rang or enum:

                                p = warn.index(']') + 2
                                warn = warn[p:]

                                if zero or nega or rang:
                                    item = 'unusual_data'
                                else:  # enum

                                    if warn.startswith('The mandatory type'):
                                        try:
                                            g = CHK_DESC_MAND_PAT.search(warn).groups()
                                        except AttributeError:
                                            g = CHK_DESC_MAND_ONE_PAT.search(warn).groups()
                                    else:
                                        try:
                                            g = CHK_DESC_PAT.search(warn).groups()
                                        except AttributeError:
                                            g = CHK_DESC_ONE_PAT.search(warn).groups()

                                    if has_key_value(MANDATORY_SF_TAG_ITEMS[file_type], content_subtype):

                                        if any(True for item in MANDATORY_SF_TAG_ITEMS[file_type][content_subtype] if item == g[0]):
                                            if not self.__reg.nefT.is_mandatory_tag(f'_{sf_category}.{g[0]}', file_type):
                                                ignorable = True  # author provides the meta data through DepUI after upload

                                    item = 'enum_mismatch_ignorable' if ignorable else 'enum_mismatch'

                                self.__reg.report.warning.appendDescription(item,
                                                                            {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                             'description': warn})

                                if self.__reg.verbose:
                                    self.__reg.log.write(f"+{self.__class_name__}.__testSfTagConsistency() "
                                                         f"++ Warning  - {warn}\n")

                            else:

                                self.__reg.report.error.appendDescription('internal_error',
                                                                          f"+{self.__class_name__}.__testSfTagConsistency() "
                                                                          "++ Error  - " + warn)

                                if self.__reg.verbose:
                                    self.__reg.log.write(f"+{self.__class_name__}.__testSfTagConsistency() "
                                                         f"++ Error  - {warn}\n")

                        # try to parse data without constraints
                        try:

                            sf_tag_data = self.__reg.nefT.check_sf_tag(sf, file_type, sf_category, sf_tag_items,
                                                                       SF_ALLOWED_TAGS[file_type][content_subtype],
                                                                       enforce_non_zero=False, enforce_sign=False,
                                                                       enforce_range=False, enforce_enum=False)

                            self.__reg.dpV.testParentChildRelation(file_name, file_type, content_subtype,
                                                                   parent_keys, list_id, sf_framecode,
                                                                   sf_framecode_dict, sf_tag_data)

                            self.__reg.sf_tag_data[content_subtype].append({'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                            'data': sf_tag_data})

                        except Exception:  # pylint: disable=broad-exception-caught
                            pass

                    except Exception as e:  # pylint: disable=broad-exception-caught

                        self.__reg.report.error.appendDescription('internal_error',
                                                                  f"+{self.__class_name__}.__testSfTagConsistency() "
                                                                  "++ Error  - " + str(e))

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__testSfTagConsistency() "
                                                 f"++ Error  - {str(e)}\n")

                    parent_keys.add(list_id)
                    if str(list_id) not in sf_framecode_dict:
                        sf_framecode_dict = {list_id: sf_framecode}

                    list_id += 1

        return self.__reg.report.getTotalErrors() == __errors

    def __validateCsValue(self) -> bool:
        """ Validate assigned chemical shift value based on BMRB chemical shift statistics.
        """

        __errors = self.__reg.report.getTotalErrors()

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            if input_source_dic['content_subtype'] is None:
                continue

            content_subtype = 'chem_shift'

            if content_subtype not in input_source_dic['content_subtype']:
                continue

            sf_category = SF_CATEGORIES[file_type][content_subtype]
            lp_category = LP_CATEGORIES[file_type][content_subtype]

            modified = False

            if self.__reg.star_data_type[fileListId] == 'Loop':
                sf = self.__reg.star_data[fileListId]
                sf_framecode = ''

                modified |= self.__reg.dpV.validateCsValue(fileListId, file_name, file_type, content_subtype,
                                                           sf, sf_framecode, lp_category)

            elif self.__reg.star_data_type[fileListId] == 'Saveframe':
                sf = self.__reg.star_data[fileListId]
                sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                modified |= self.__reg.dpV.validateCsValue(fileListId, file_name, file_type, content_subtype,
                                                           sf, sf_framecode, lp_category)

            else:

                for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                    sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                    if not any(True for loop in sf.loops if loop.category == lp_category):
                        continue

                    modified |= self.__reg.dpV.validateCsValue(fileListId, file_name, file_type, content_subtype,
                                                               sf, sf_framecode, lp_category)

            if modified and fileListId == 0:
                self.__depositNmrData()

        return self.__reg.report.getTotalErrors() == __errors

    def __extractToNmrIf__(self) -> bool:
        """ Extract NMR metadata of NMR-STAR file (as primary source) and model file (as secondary source)
            to NMRIF file, if possible.
        """

        if NMRIF_FILE_PATH_KEY not in self.__reg.outputParamDict:
            return False

        if os.path.exists(self.__reg.outputParamDict[NMRIF_FILE_PATH_KEY]):
            return False

        if len(self.__reg.star_data) == 0 or not isinstance(self.__reg.star_data[0], pynmrstar.Entry):
            return self.__extractToNmrIfFromModel__()

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_type = input_source_dic['file_type']

        if file_type != 'nmr-star':
            return self.__extractToNmrIfFromModel__()

        master_entry = self.__reg.star_data[0]

        self.__reg.sf_category_list, self.__reg.lp_category_list = self.__reg.nefT.get_inventory_list(master_entry)

        ann = BmrbAnnTasks(self.__reg)

        if not self.__reg.internal_mode and self.__reg.report.getInputSourceIdOfCoord() >= 0\
           and self.__reg.cR.hasCategory('database_2'):

            sf_category = 'entry_information'

            database_code = self.__reg.cR.getDictListWithFilter('database_2',
                                                                [{'name': 'database_code', 'type': 'str'}],
                                                                [{'name': 'database_id', 'type': 'str', 'value': 'BMRB'}])

            if len(database_code) > 0:
                derived_entry_id = database_code[0]['database_code']
                derived_entry_title = None

                if sf_category in self.__reg.sf_category_list:
                    sf = master_entry.get_saveframes_by_category(sf_category)[0]
                    derived_entry_title = get_first_sf_tag(sf, 'Title', None)

                ann.setProvenanceInfo(derived_entry_id, derived_entry_title)

        ann.perform(master_entry)

        ann = OneDepAnnTasks(self.__reg)

        return ann.extract(master_entry, self.__reg.cR, self.__reg.outputParamDict[NMRIF_FILE_PATH_KEY])

    def __extractToNmrIfFromModel__(self) -> bool:
        """ Extract NMR metadata of the model file to NMRIF file, if possible.
        """

        if NMRIF_FILE_PATH_KEY not in self.__reg.outputParamDict:
            return False

        if os.path.exists(self.__reg.outputParamDict[NMRIF_FILE_PATH_KEY]):
            return False

        src_id = self.__reg.report.getInputSourceIdOfCoord()

        if src_id < 0:
            return False

        ann = OneDepAnnTasks(self.__reg)

        return ann.extract(None, self.__reg.cR, self.__reg.outputParamDict[NMRIF_FILE_PATH_KEY])

    def __remediateCsLoop__(self) -> bool:
        """ Remediate assigned chemical shift loop based on coordinates.
        """

        if not self.__reg.combined_mode and not self.__reg.bmrb_only and self.__reg.has_star_entity\
           and self.__reg.report.error.hasSequenceMismatchErrorInCsLoop():  # DAOTHER-10487
            return False

        __errors = self.__reg.report.getTotalErrors()

        content_subtype = 'chem_shift'

        if not self.__reg.native_combined:
            self.__reg.lp_data[content_subtype] = []

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_type = input_source_dic['file_type']
            original_file_name = input_source_dic['original_file_name']

            if input_source_dic['content_subtype'] is None:
                continue

            if content_subtype not in input_source_dic['content_subtype']:
                continue

            sf_category = SF_CATEGORIES[file_type][content_subtype]
            lp_category = LP_CATEGORIES[file_type][content_subtype]

            modified = False

            list_id = 1

            if self.__reg.star_data_type[fileListId] == 'Loop':
                sf = self.__reg.star_data[fileListId]
                sf_framecode = ''

                if file_type == 'nmr-star' and original_file_name not in EMPTY_VALUE:
                    set_sf_tag(sf, 'Data_file_name', original_file_name)

                modified |= self.__reg.dpR.remediateCsLoop(fileListId, file_type, content_subtype,
                                                           sf, list_id, sf_framecode, lp_category)

            elif self.__reg.star_data_type[fileListId] == 'Saveframe':
                sf = self.__reg.star_data[fileListId]
                sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                if file_type == 'nmr-star':
                    data_file_name = get_first_sf_tag(sf, 'Data_file_name')
                    len_data_file_name = len(data_file_name)
                    if len_data_file_name > 0:
                        data_file_name.strip("'").strip('"')
                        _len_data_file_name = len(data_file_name)
                        if _len_data_file_name > 0 and _len_data_file_name != len_data_file_name:
                            set_sf_tag(sf, 'Data_file_name', data_file_name)
                    elif original_file_name not in EMPTY_VALUE:
                        set_sf_tag(sf, 'Data_file_name', original_file_name)

                modified |= self.__reg.dpR.remediateCsLoop(fileListId, file_type, content_subtype,
                                                           sf, list_id, sf_framecode, lp_category)

            else:

                for list_id, sf in enumerate(self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category), start=1):
                    sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                    if not any(True for loop in sf.loops if loop.category == lp_category):
                        continue

                    if file_type == 'nmr-star':
                        data_file_name = get_first_sf_tag(sf, 'Data_file_name')
                        len_data_file_name = len(data_file_name)
                        if len_data_file_name > 0:
                            data_file_name.strip("'").strip('"')
                            _len_data_file_name = len(data_file_name)
                            if _len_data_file_name > 0 and _len_data_file_name != len_data_file_name:
                                set_sf_tag(sf, 'Data_file_name', data_file_name)
                        elif original_file_name not in EMPTY_VALUE:
                            set_sf_tag(sf, 'Data_file_name', original_file_name)

                    modified |= self.__reg.dpR.remediateCsLoop(fileListId, file_type, content_subtype,
                                                               sf, list_id, sf_framecode, lp_category)

            if modified:

                # update _Entity_assembly.Experimental_data_reported
                if file_type == 'nmr-star' and len(self.__reg.ent_asym_id_with_exptl_data) > 0:

                    _content_subtype = 'poly_seq'

                    _sf_category = SF_CATEGORIES[file_type][_content_subtype]

                    try:

                        _sf = self.__reg.star_data[fileListId].get_saveframes_by_category(_sf_category)[0]

                        try:
                            _loop = _sf.get_loop('_Entity_assembly')

                            if 'Experimental_data_reported' in _loop.tags:
                                id_col = _loop.tags.index('ID')
                                exptl_data_rep_col = _loop.tags.index('Experimental_data_reported')

                                for _row in _loop:
                                    if _row[id_col] in self.__reg.ent_asym_id_with_exptl_data:
                                        _row[exptl_data_rep_col] = 'yes'

                        except KeyError:
                            pass

                    except IndexError:
                        pass

                if fileListId == 0:
                    self.__depositNmrData()

        return self.__reg.report.getTotalErrors() == __errors

    def __removeUnusedPdbInsCode(self) -> bool:
        """ Remove unused PDB_ind_code tags from loops.
        """

        __errors = self.__reg.report.getTotalErrors()

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_type = input_source_dic['file_type']

            if file_type != 'nmr-star':
                continue

            if input_source_dic['content_subtype'] is None:
                continue

            modified = False

            for content_subtype in input_source_dic['content_subtype']:

                if content_subtype not in ('chem_shift', 'dist_restraint', 'dihed_restraint', 'rdc_restraint'):
                    continue

                sf_category = SF_CATEGORIES[file_type][content_subtype]
                lp_category = LP_CATEGORIES[file_type][content_subtype]

                if self.__reg.star_data_type[fileListId] == 'Loop':
                    sf = self.__reg.star_data[fileListId]

                    modified |= self.__reg.dpR.removeUnusedPdbInsCode(fileListId, content_subtype, sf, lp_category)

                elif self.__reg.star_data_type[fileListId] == 'Saveframe':
                    sf = self.__reg.star_data[fileListId]

                    modified |= self.__reg.dpR.removeUnusedPdbInsCode(fileListId, content_subtype, sf, lp_category)

                else:

                    for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):

                        if not any(True for loop in sf.loops if loop.category == lp_category):
                            continue

                        modified |= self.__reg.dpR.removeUnusedPdbInsCode(fileListId, content_subtype, sf, lp_category)

            if modified and fileListId == 0:
                self.__depositNmrData()

        return self.__reg.report.getTotalErrors() == __errors

    def __testCsPseudoAtomNameConsistencyInMrLoop(self) -> bool:
        """ Perform consistency test on pseudoatom names between assigned chemical shifts and restraints. (DAOTHER-7681, issue #1)
        """

        __errors = self.__reg.report.getTotalErrors()

        for fileListId in range(self.__reg.file_path_list_len):

            if fileListId >= len(self.__reg.star_data_type) or self.__reg.star_data_type[fileListId] != 'Entry':
                continue

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            # DAOTHER-9405
            if file_type == 'nef':
                continue

            if input_source_dic['content_subtype'] is None\
               or 'chem_shift' not in input_source_dic['content_subtype']:
                continue

            rescue_mode = self.__reg.cmpl_missing_data and input_source_dic['content_subtype']['chem_shift'] == 1

            cs_item_names = ITEM_NAMES_IN_CS_LOOP[file_type]
            cs_chain_id_name = cs_item_names['chain_id']
            cs_seq_id_name = cs_item_names['seq_id']
            cs_comp_id_name = cs_item_names['comp_id']
            cs_atom_id_name = cs_item_names['atom_id']
            cs_value_name = cs_item_names['value']

            missing_cs_atoms = []

            for content_subtype in input_source_dic['content_subtype']:

                if content_subtype in ('dist_restraint', 'dihed_restraint', 'rdc_restraint'):

                    sf_category = SF_CATEGORIES[file_type][content_subtype]
                    lp_category = LP_CATEGORIES[file_type][content_subtype]

                    for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                        sf_framecode = get_first_sf_tag(sf, 'Sf_framecode')

                        try:
                            cs_data, cs_list = next((lp['data'], lp['sf_framecode']) for lp in self.__reg.lp_data['chem_shift']
                                                    if lp['file_name'] == file_name)
                        except StopIteration:
                            continue

                        max_dim = 3 if content_subtype in ('dist_restraint', 'rdc_restraint') else 5

                        item_names = []
                        for j in range(1, max_dim):
                            _item_names = {}
                            for k, v in ITEM_NAMES_IN_PK_LOOP[file_type].items():
                                if '%s' in v:
                                    v = v % j
                                _item_names[k] = v
                            item_names.append(_item_names)

                        num_dim = max_dim - 1

                        chain_id_names, seq_id_names, comp_id_names, atom_id_names = [], [], [], []

                        for d in range(num_dim):

                            chain_id_names.append(item_names[d]['chain_id'])
                            seq_id_names.append(item_names[d]['seq_id'])
                            comp_id_names.append(item_names[d]['comp_id'])
                            atom_id_names.append(item_names[d]['atom_id'])

                        id_tag = CONSIST_ID_TAGS[file_type][content_subtype]

                        try:

                            lp_data = next((lp['data'] for lp in self.__reg.lp_data[content_subtype]
                                            if lp['file_name'] == file_name and lp['sf_framecode'] == sf_framecode), None)

                            if lp_data is not None:

                                for row in lp_data:
                                    for d in range(num_dim):
                                        chain_id = row.get(chain_id_names[d])
                                        seq_id = row.get(seq_id_names[d])
                                        comp_id = row.get(comp_id_names[d])
                                        atom_id = row.get(atom_id_names[d])

                                        if chain_id in EMPTY_VALUE or seq_id in EMPTY_VALUE\
                                           or comp_id in EMPTY_VALUE or atom_id in EMPTY_VALUE:
                                            continue

                                        _atom_ids = self.__reg.dpV.getAtomIdList(comp_id, atom_id)

                                        len_atom_id = len(_atom_ids)

                                        if len_atom_id == 0:
                                            atom_id_ = atom_id

                                        elif len_atom_id == 1 and atom_id == _atom_ids[0]:
                                            atom_id_ = atom_id

                                        else:  # representative atom id
                                            atom_id_ = _atom_ids[0]

                                        if self.__reg.csStat.getMaxAmbigCodeWoSetId(comp_id, atom_id_) < 2:
                                            continue

                                        _atom_id_ = atom_id_

                                        if file_type == 'nmr-star' and self.__reg.dpV.isNmrAtomName(comp_id, atom_id):
                                            pass
                                        else:
                                            atom_id_ = atom_id

                                        atom_ids_w_cs = [_row[cs_atom_id_name] for _row in cs_data
                                                         if _row[cs_chain_id_name] == chain_id
                                                         and _row[cs_seq_id_name] == seq_id
                                                         and _row[cs_comp_id_name] == comp_id]

                                        if atom_id_ in atom_ids_w_cs:
                                            continue

                                        has_chem_shift = False

                                        for atom_id_w_cs in atom_ids_w_cs:
                                            _atom_id_w_cs = self.__reg.dpV.getAtomIdList(comp_id, atom_id_w_cs)
                                            if any(True for _atom_id in _atom_ids if _atom_id in _atom_id_w_cs):
                                                has_chem_shift = True
                                                break

                                        if has_chem_shift:
                                            continue

                                        gem_atom_id = self.__reg.csStat.getGeminalAtom(comp_id, _atom_id_)

                                        if gem_atom_id is None:
                                            continue

                                        gem_atom_id_w_cs = None

                                        atom_ids_w_cs = [_row[cs_atom_id_name] for _row in cs_data
                                                         if _row[cs_chain_id_name] == chain_id
                                                         and _row[cs_seq_id_name] == seq_id
                                                         and _row[cs_comp_id_name] == comp_id]

                                        for atom_id_w_cs in atom_ids_w_cs:
                                            _atom_id_w_cs = self.__reg.dpV.getAtomIdList(comp_id, atom_id_w_cs)
                                            if gem_atom_id in _atom_id_w_cs:
                                                gem_atom_id_w_cs = atom_id_w_cs
                                                break

                                        if gem_atom_id_w_cs is None:
                                            continue

                                        if content_subtype == 'dist_restraint':
                                            subtype_name = "distance restraint"
                                        elif content_subtype == 'dihed_restraint':
                                            subtype_name = "dihedral angle restraint"
                                        else:
                                            subtype_name = "RDC restraint"

                                        if _atom_id_ in self.__reg.csStat.getMethylAtoms(comp_id)\
                                           and content_subtype == 'dist_restraint'\
                                           and not self.__reg.remediation_mode:

                                            cs_atom_id_map = {'chain_id': chain_id, 'seq_id': seq_id, 'comp_id': comp_id,
                                                              'src_atom_id': gem_atom_id_w_cs, 'dst_atom_id': atom_id,
                                                              'content_subtype_name': f'{subtype_name}s'}

                                            if cs_atom_id_map not in missing_cs_atoms:
                                                missing_cs_atoms.append(cs_atom_id_map)

                                            if rescue_mode:
                                                continue

                                            err = f"[Check row of {id_tag} {row[id_tag]}] Assignment of {subtype_name} "\
                                                + self.__reg.dpV.getReducedAtomNotation(chain_id_names[d], chain_id,
                                                                                        seq_id_names[d], seq_id,
                                                                                        comp_id_names[d], comp_id,
                                                                                        atom_id_names[d], atom_id)\
                                                + " was not found in assigned chemical shifts. In contrast, "\
                                                + self.__reg.dpV.getReducedAtomNotation(chain_id_names[d], chain_id,
                                                                                        seq_id_names[d], seq_id,
                                                                                        comp_id_names[d], comp_id,
                                                                                        atom_id_names[d], gem_atom_id_w_cs)\
                                                + f" is in the assigned chemical shifts of {cs_list!r} saveframe."

                                            self.__reg.report.error.appendDescription('invalid_data',
                                                                                      {'file_name': file_name,
                                                                                       'sf_framecode': sf_framecode,
                                                                                       'category': lp_category, 'description': err})

                                            if self.__reg.verbose:
                                                self.__reg.log.write(f"+{self.__class_name__}.__testCsPseudoAtomNameConsistencyInMrLoop() "  # noqa: E501, pylint: disable=line-too-long
                                                                     f"++ ValueError  - {err}\n")

                                        else:

                                            warn = f"[Check row of {id_tag} {row[id_tag]}] Assignment of {subtype_name} "\
                                                + self.__reg.dpV.getReducedAtomNotation(chain_id_names[d], chain_id,
                                                                                        seq_id_names[d], seq_id,
                                                                                        comp_id_names[d], comp_id,
                                                                                        atom_id_names[d], atom_id)\
                                                + " was not found in assigned chemical shifts. In contrast, "\
                                                + self.__reg.dpV.getReducedAtomNotation(chain_id_names[d], chain_id,
                                                                                        seq_id_names[d], seq_id,
                                                                                        comp_id_names[d], comp_id,
                                                                                        atom_id_names[d], gem_atom_id_w_cs)\
                                                + f" is in the assigned chemical shifts of {cs_list!r} saveframe."

                                            self.__reg.report.warning.appendDescription('missing_data',
                                                                                        {'file_name': file_name,
                                                                                         'sf_framecode': sf_framecode,
                                                                                         'category': lp_category,
                                                                                         'description': warn})

                                            if self.__reg.verbose:
                                                self.__reg.log.write(f"+{self.__class_name__}.__testCsPseudoAtomNameConsistencyInMrLoop() "  # noqa: E501, pylint: disable=line-too-long
                                                                     f"++ Warning  - {warn}\n")

                        except Exception as e:  # pylint: disable=broad-exception-caught

                            self.__reg.report.error.appendDescription('internal_error',
                                                                      f"+{self.__class_name__}.__testCsPseudoAtomNameConsistencyInMrLoop() "  # noqa: E501, pylint: disable=line-too-long
                                                                      "++ Error  - " + str(e))

                            if self.__reg.verbose:
                                self.__reg.log.write(f"+{self.__class_name__}.__testCsPseudoAtomNameConsistencyInMrLoop() "
                                                     f"++ Error  - {str(e)}\n")

            if rescue_mode and len(missing_cs_atoms) > 0:

                content_subtype = 'chem_shift'

                sf_category = SF_CATEGORIES[file_type][content_subtype]
                lp_category = LP_CATEGORIES[file_type][content_subtype]

                star_data = copy.copy(self.__reg.star_data[fileListId])

                for sf in star_data.get_saveframes_by_category(sf_category):
                    sf_framecode = get_first_sf_tag(sf, 'Sf_framecode')

                    loop = sf.get_loop(lp_category)

                    lp = pynmrstar.Loop.from_scratch(lp_category)

                    lp.add_tag(loop.tags)

                    chain_id_col = loop.tags.index(cs_chain_id_name)
                    seq_id_col = loop.tags.index(cs_seq_id_name)
                    comp_id_col = loop.tags.index(cs_comp_id_name)
                    atom_id_col = loop.tags.index(cs_atom_id_name)
                    value_col = loop.tags.index(cs_value_name)

                    for row in loop:
                        lp.add_data(row)

                        chain_id = row[chain_id_col]
                        try:
                            seq_id = int(row[seq_id_col])
                        except (ValueError, TypeError):
                            continue
                        comp_id = row[comp_id_col]
                        atom_id = row[atom_id_col]
                        value = row[value_col]

                        _missing_cs_atoms = [missing_cs_atom for missing_cs_atom in missing_cs_atoms
                                             if missing_cs_atom['chain_id'] == chain_id
                                             and missing_cs_atom['seq_id'] == seq_id
                                             and missing_cs_atom['comp_id'] == comp_id
                                             and missing_cs_atom['src_atom_id'] == atom_id]

                        if len(_missing_cs_atoms) == 0:
                            continue

                        _subtype_name = ' and '.join([missing_cs_atom['content_subtype_name']
                                                      for missing_cs_atom in _missing_cs_atoms])

                        missing_cs_atom = _missing_cs_atoms[0]

                        _row = copy.copy(row)
                        _row[atom_id_col] = missing_cs_atom['dst_atom_id']
                        lp.data.append(_row)

                        warn = "The unbound resonance assignment "\
                            + self.__reg.dpV.getReducedAtomNotation(cs_chain_id_name, chain_id, cs_seq_id_name, seq_id,
                                                                    cs_comp_id_name, comp_id,
                                                                    cs_atom_id_name, missing_cs_atom['dst_atom_id'])\
                            + f" in {_subtype_name} has been added to the assigned chemical shifts by referring to "\
                            + self.__reg.dpV.getReducedAtomNotation(cs_chain_id_name, chain_id, cs_seq_id_name, seq_id,
                                                                    cs_comp_id_name, comp_id,
                                                                    cs_atom_id_name, missing_cs_atom['src_atom_id'])\
                            + f", {value} ppm."

                        self.__reg.report.warning.appendDescription('complemented_chemical_shift',
                                                                    {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                     'category': lp_category, 'description': warn})

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__testCsPseudoAtomNameConsistencyInMrLoop() "
                                                 f"++ Warning  - {warn}\n")

                    del sf[loop]

                    sf.add_loop(lp)

                    parent_pointer = 1
                    for idx, lp_data in enumerate(self.__reg.lp_data[content_subtype]):
                        if lp_data['file_name'] == file_name and lp_data['sf_framecode'] == sf_framecode:
                            del self.__reg.lp_data[content_subtype][idx]
                            parent_pointer = idx + 1
                            break

                    key_items = self.__reg.key_items[file_type][content_subtype]
                    data_items = DATA_ITEMS[file_type][content_subtype]
                    allowed_tags = ALLOWED_TAGS[file_type][content_subtype]
                    disallowed_tags = None

                    try:

                        lp_data = self.__reg.nefT.check_data(sf, lp_category, key_items, data_items,
                                                             allowed_tags, disallowed_tags, parent_pointer=parent_pointer,
                                                             test_on_index=True, enforce_non_zero=True, enforce_sign=True,
                                                             enforce_range=True, enforce_enum=True,
                                                             enforce_allowed_tags=(file_type == 'nmr-star'),
                                                             excl_missing_data=self.__reg.excl_missing_data)[0]

                        self.__reg.lp_data[content_subtype].append({'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                    'data': lp_data})

                    except Exception as e:  # pylint: disable=broad-exception-caught

                        self.__reg.report.error.appendDescription('internal_error',
                                                                  f"+{self.__class_name__}.__testCsPseudoAtomNameConsistencyInMrLoop() "  # noqa: E501, pylint: disable=line-too-long
                                                                  "++ Error  - " + str(e))

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__testCsPseudoAtomNameConsistencyInMrLoop() "
                                                 f"++ Error  - {str(e)}\n")

                # self.__depositNmrData()

        return self.__reg.report.getTotalErrors() == __errors

    def __testCsValueConsistencyInPkLoop(self) -> bool:
        """ Perform consistency test on peak position and assignment of spectral peaks.
        """

        __errors = self.__reg.report.getTotalErrors()

        for fileListId in range(self.__reg.file_path_list_len):

            if fileListId >= len(self.__reg.star_data_type) or self.__reg.star_data_type[fileListId] != 'Entry':
                continue

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            if input_source_dic['content_subtype'] is None:
                continue

            content_subtype = 'spectral_peak'

            if content_subtype not in input_source_dic['content_subtype']:
                continue

            sf_category = SF_CATEGORIES[file_type][content_subtype]
            lp_category = LP_CATEGORIES[file_type][content_subtype]

            cs_item_names = ITEM_NAMES_IN_CS_LOOP[file_type]
            cs_chain_id_name = cs_item_names['chain_id']
            cs_seq_id_name = cs_item_names['seq_id']
            cs_comp_id_name = cs_item_names['comp_id']
            cs_atom_id_name = cs_item_names['atom_id']
            cs_value_name = cs_item_names['value']
            cs_error_name = cs_item_names['error']
            cs_atom_type = cs_item_names['atom_type']
            cs_iso_number = cs_item_names['isotope_number']

            for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                try:

                    cs_list = get_first_sf_tag(sf, CS_LIST_SF_TAG_NAME[file_type])

                except Exception:  # pylint: disable=broad-exception-caught
                    continue

                cs_data = None

                try:

                    cs_data = next(lp['data'] for lp in self.__reg.lp_data['chem_shift']
                                   if lp['file_name'] == file_name
                                   and (lp['sf_framecode'] == cs_list
                                        or (lp['sf_framecode'].startswith('nef_')  # DAOTHER-9623, 10781
                                            and lp['sf_framecode'].endswith(cs_list))))

                except StopIteration:

                    if cs_list not in EMPTY_VALUE:

                        if fileListId == 0:

                            err = "Assigned chemical shifts are required to verify the consistency of assigned peak list. "\
                                f"Referred {cs_list!r} saveframe containing the assigned chemical shift does not exist."

                            self.__reg.report.error.appendDescription('missing_mandatory_content',
                                                                      {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                       'description': err})

                            if self.__reg.verbose:
                                self.__reg.log.write(f"+{self.__class_name__}.__testCsValueConsistencyInPkLoop() "
                                                     f"++ Error  - {err}\n")

                            continue

                        cs_input_source = self.__reg.report.input_sources[0]
                        cs_input_source_dic = cs_input_source.get()
                        cs_file_name = cs_input_source_dic['file_name']

                        try:
                            cs_data = next(lp['data'] for lp in self.__reg.lp_data['chem_shift'] if lp['file_name'] == cs_file_name)
                        except StopIteration:
                            continue

                try:

                    _num_dim = get_first_sf_tag(sf, NUM_DIM_ITEMS[file_type])
                    num_dim = int(_num_dim)

                    if num_dim not in range(1, MAX_DIM_NUM_OF_SPECTRA):
                        raise ValueError()

                except ValueError:  # raised error already at __testIndexConsistency()
                    return False

                max_dim = num_dim + 1

                aux_data = next((lp['data'] for lp in self.__reg.aux_data[content_subtype]
                                 if lp['file_name'] == file_name and lp['sf_framecode'] == sf_framecode
                                 and lp['category'] == AUX_LP_CATEGORIES[file_type][content_subtype][0]), None)

                axis_codes, alt_axis_codes, abs_pk_pos, sp_widths = [], [], [], []

                if aux_data is not None and len(aux_data) > 0:
                    for i in range(1, max_dim):
                        for sp_dim in aux_data:
                            if file_type == 'nef':
                                if sp_dim['dimension_id'] != i:
                                    continue
                                axis_codes.append(sp_dim['axis_code'])
                                alt_axis_codes.append(sp_dim['axis_code'])
                                abs_pk_pos.append(False if 'absolute_peak_poistions' not in sp_dim
                                                  else sp_dim['absolute_peak_positions'])
                                sp_width = None if 'axis_unit' not in sp_dim else sp_dim.get('spectral_width')
                                if 'axis_unit' in sp_dim and sp_dim['axis_unit'] == 'Hz'\
                                   and 'spectrometer_frequency' in sp_dim and sp_width is not None:
                                    sp_freq = sp_dim['spectrometer_frequency']
                                    if sp_freq not in EMPTY_VALUE:
                                        sp_width /= sp_freq
                                sp_widths.append(sp_width)
                            else:
                                if sp_dim['ID'] != i:
                                    continue
                                axis_codes.append(sp_dim['Axis_code'])
                                if 'Atom_isotope_number' in sp_dim and 'Atom_type' in sp_dim and sp_dim['Atom_type'] is not None:
                                    alt_axis_codes.append(str(sp_dim['Atom_isotope_number']) + sp_dim['Atom_type'])
                                else:
                                    alt_axis_codes.append(sp_dim['Axis_code'])
                                abs_pk_pos.append(False if 'Absolute_peak_positions' not in sp_dim
                                                  else sp_dim['Absolute_peak_positions'])
                                sp_width = None if 'Sweep_width_units' not in sp_dim else sp_dim.get('Sweep_width')
                                if 'Sweep_width_units' in sp_dim and sp_dim['Sweep_width_units'] == 'Hz'\
                                   and 'Spectrometer_frequency' in sp_dim and sp_width is not None:
                                    sp_freq = sp_dim['Spectrometer_frequency']
                                    if sp_freq not in EMPTY_VALUE:
                                        sp_width /= sp_freq
                                sp_widths.append(sp_width)
                            break
                else:
                    for i in range(num_dim):
                        axis_codes.append(None)
                        alt_axis_codes.append(None)
                        abs_pk_pos.append(False)
                        sp_widths.append(None)

                aux_data = next((lp['data'] for lp in self.__reg.aux_data[content_subtype]
                                 if lp['file_name'] == file_name and lp['sf_framecode'] == sf_framecode
                                 and lp['category'] == AUX_LP_CATEGORIES[file_type][content_subtype][1]), None)

                onebond = [[False] * num_dim for i in range(num_dim)]
                if aux_data is not None:
                    for sp_dim_trans in aux_data:
                        if file_type == 'nef':
                            if sp_dim_trans['transfer_type'] == 'onebond':
                                dim_1 = sp_dim_trans['dimension_1']
                                dim_2 = sp_dim_trans['dimension_2']
                                onebond[dim_1 - 1][dim_2 - 1] = True
                                onebond[dim_2 - 1][dim_1 - 1] = True
                        else:
                            if sp_dim_trans['Type'] == 'onebond':
                                dim_1 = sp_dim_trans['Spectral_dim_ID_1']
                                dim_2 = sp_dim_trans['Spectral_dim_ID_2']
                                onebond[dim_1 - 1][dim_2 - 1] = True
                                onebond[dim_2 - 1][dim_1 - 1] = True

                jcoupling = [[False] * num_dim for i in range(num_dim)]
                if aux_data is not None:
                    for sp_dim_trans in aux_data:
                        if file_type == 'nef':
                            if sp_dim_trans['transfer_type'] is not None and sp_dim_trans['transfer_type'].startswith('j'):
                                dim_1 = sp_dim_trans['dimension_1']
                                dim_2 = sp_dim_trans['dimension_2']
                                jcoupling[dim_1 - 1][dim_2 - 1] = True
                                jcoupling[dim_2 - 1][dim_1 - 1] = True
                        else:
                            if sp_dim_trans['Type'] is not None and sp_dim_trans['Type'].startswith('j'):
                                dim_1 = sp_dim_trans['Spectral_dim_ID_1']
                                dim_2 = sp_dim_trans['Spectral_dim_ID_2']
                                jcoupling[dim_1 - 1][dim_2 - 1] = True
                                jcoupling[dim_2 - 1][dim_1 - 1] = True

                relayed = [[False] * num_dim for i in range(num_dim)]
                if aux_data is not None:
                    for sp_dim_trans in aux_data:
                        if file_type == 'nef':
                            if sp_dim_trans['transfer_type'] is not None and sp_dim_trans['transfer_type'].startswith('relayed'):
                                dim_1 = sp_dim_trans['dimension_1']
                                dim_2 = sp_dim_trans['dimension_2']
                                relayed[dim_1 - 1][dim_2 - 1] = True
                                relayed[dim_2 - 1][dim_1 - 1] = True
                        else:
                            if sp_dim_trans['Type'] is not None and sp_dim_trans['Type'].startswith('relayed'):
                                dim_1 = sp_dim_trans['Spectral_dim_ID_1']
                                dim_2 = sp_dim_trans['Spectral_dim_ID_2']
                                relayed[dim_1 - 1][dim_2 - 1] = True
                                relayed[dim_2 - 1][dim_1 - 1] = True

                item_names = []
                for dim in range(1, max_dim):
                    _d = {}
                    for k, v in ITEM_NAMES_IN_PK_LOOP[file_type].items():
                        if '%s' in v:
                            v = v % dim
                        _d[k] = v
                    item_names.append(_d)

                chain_id_names, seq_id_names, comp_id_names, atom_id_names, position_names = [], [], [], [], []

                for d in range(num_dim):
                    chain_id_names.append(item_names[d]['chain_id'])
                    seq_id_names.append(item_names[d]['seq_id'])
                    comp_id_names.append(item_names[d]['comp_id'])
                    atom_id_names.append(item_names[d]['atom_id'])
                    position_names.append(item_names[d]['position'])

                id_tag = CONSIST_ID_TAGS[file_type][content_subtype]

                try:

                    lp_data = next((lp['data'] for lp in self.__reg.lp_data[content_subtype]
                                    if lp['file_name'] == file_name and lp['sf_framecode'] == sf_framecode), None)

                    if lp_data is not None and cs_data is not None:

                        for row in lp_data:
                            for d in range(num_dim):

                                if not (chain_id_names[d] in row and seq_id_names[d] in row
                                        and comp_id_names[d] in row and atom_id_names[d] in row):
                                    continue

                                chain_id = row.get(chain_id_names[d])
                                seq_id = row.get(seq_id_names[d])
                                comp_id = row.get(comp_id_names[d])
                                atom_id = _atom_id = row.get(atom_id_names[d])

                                if chain_id in EMPTY_VALUE or seq_id in EMPTY_VALUE\
                                   or comp_id in EMPTY_VALUE or atom_id in EMPTY_VALUE:
                                    continue

                                position = row[position_names[d]]

                                _atom_ids = self.__reg.dpV.getAtomIdList(comp_id, atom_id)

                                len_atom_id = len(_atom_ids)

                                if len_atom_id == 0:
                                    atom_id_ = atom_id

                                elif len_atom_id == 1 and atom_id == _atom_ids[0]:
                                    atom_id_ = atom_id

                                else:  # representative atom id
                                    atom_id_ = _atom_ids[0]

                                cs_idx = -1

                                if file_type == 'nmr-star' and self.__reg.dpV.isNmrAtomName(comp_id, atom_id):
                                    pass
                                else:
                                    atom_id_ = atom_id

                                atom_ids_w_cs = [_row[cs_atom_id_name] for _row in cs_data
                                                 if _row[cs_chain_id_name] == chain_id
                                                 and _row[cs_seq_id_name] == seq_id
                                                 and _row[cs_comp_id_name] == comp_id]

                                if atom_id_ in atom_ids_w_cs:
                                    cs_idx = atom_ids_w_cs.index(atom_id_)

                                else:
                                    for atom_id_w_cs in atom_ids_w_cs:
                                        _atom_id_w_cs = self.__reg.dpV.getAtomIdList(comp_id, atom_id_w_cs)
                                        if any(True for _atom_id_ in _atom_ids if _atom_id_ in _atom_id_w_cs):
                                            cs_idx = atom_ids_w_cs.index(atom_id_w_cs)
                                            break

                                if cs_idx == -1:

                                    err = f"[Check row of {id_tag} {row[id_tag]}] Assignment of spectral peak "\
                                        + self.__reg.dpV.getReducedAtomNotation(chain_id_names[d], chain_id,
                                                                                seq_id_names[d], seq_id,
                                                                                comp_id_names[d], comp_id,
                                                                                atom_id_names[d], atom_id)\
                                        + f" was not found in assigned chemical shifts of {cs_list!r} saveframe."

                                    self.__reg.report.warning.appendDescription('insufficient_data',
                                                                                {'file_name': file_name,
                                                                                 'sf_framecode': sf_framecode,
                                                                                 'category': lp_category, 'description': err})

                                    if self.__reg.verbose:
                                        self.__reg.log.write(f"+{self.__class_name__}.__testCsValueConsistencyInPkLoop() "
                                                             f"++ Warning  - {err}\n")

                                else:

                                    cs_intra = [_row for _row in cs_data
                                                if _row[cs_chain_id_name] == chain_id
                                                and _row[cs_seq_id_name] == seq_id
                                                and _row[cs_comp_id_name] == comp_id]

                                    cs = cs_intra[cs_idx]

                                    value = cs[cs_value_name]
                                    error = cs[cs_error_name]

                                    if value in EMPTY_VALUE:
                                        continue

                                    if error is None or error < 1.0e-3\
                                       or error * self.__reg.cs_diff_error_scaled_by_sigma > CS_UNCERT_MAX:
                                        error = CS_UNCERT_MAX
                                    else:
                                        error *= self.__reg.cs_diff_error_scaled_by_sigma

                                    if abs(position - value) > error:

                                        if not abs_pk_pos[d] and sp_widths[d] is not None:
                                            if position < value:
                                                while position < value:
                                                    position += sp_widths[d]
                                            elif position > value:
                                                while position > value:
                                                    position -= sp_widths[d]

                                        if abs(position - value) > error and sp_widths[d] is not None:

                                            if CS_RANGE_MIN < sp_widths[d] < CS_RANGE_MAX:

                                                err = f"[Check row of {id_tag} {row[id_tag]}] "\
                                                    f"Peak position of spectral peak {position_names[d]} {position} ("\
                                                    + self.__reg.dpV.getReducedAtomNotation(chain_id_names[d], chain_id,
                                                                                            seq_id_names[d], seq_id,
                                                                                            comp_id_names[d], comp_id,
                                                                                            atom_id_names[d], atom_id)\
                                                    + f") in {sf_framecode!r} saveframe is inconsistent "\
                                                    f"with the assigned chemical shift value {value} "\
                                                    f"(difference {position - value:.3f}, tolerance {error}) "\
                                                    f"in {cs_list!r} saveframe."

                                                if error >= CS_UNCERT_MAX and not self.__reg.remediation_mode:

                                                    self.__reg.report.error.appendDescription('invalid_data',
                                                                                              {'file_name': file_name,
                                                                                               'sf_framecode': sf_framecode,
                                                                                               'category': lp_category,
                                                                                               'description': err})

                                                    if self.__reg.verbose:
                                                        self.__reg.log.write(f"+{self.__class_name__}.__testCsValueConsistencyInPkLoop() "  # noqa: E501, pylint: disable=line-too-long
                                                                             f"++ ValueError  - {err}\n")

                                                else:

                                                    self.__reg.report.warning.appendDescription('unusual_chemical_shift',
                                                                                                {'file_name': file_name,
                                                                                                 'sf_framecode': sf_framecode,
                                                                                                 'category': lp_category,
                                                                                                 'description': err})

                                                    if self.__reg.verbose:
                                                        self.__reg.log.write(f"+{self.__class_name__}.__testCsValueConsistencyInPkLoop() "  # noqa: E501, pylint: disable=line-too-long
                                                                             f"++ Warning  - {err}\n")

                                    axis_code = str(cs[cs_iso_number]) + cs[cs_atom_type]

                                    if axis_codes[d] is not None and d < num_dim\
                                       and axis_code not in (axis_codes[d], alt_axis_codes[d])\
                                       and cs[cs_atom_type] != axis_codes[d]:

                                        err = f"[Check row of {id_tag} {row[id_tag]}] Assignment of spectral peak "\
                                            + self.__reg.dpV.getReducedAtomNotation(chain_id_names[d], chain_id,
                                                                                    seq_id_names[d], seq_id,
                                                                                    comp_id_names[d], comp_id,
                                                                                    atom_id_names[d], atom_id)\
                                            + f" is inconsistent with axis code {axis_code} vs {axis_codes[d]}."

                                        self.__reg.report.error.appendDescription('invalid_data',
                                                                                  {'file_name': file_name,
                                                                                   'sf_framecode': sf_framecode,
                                                                                   'category': lp_category, 'description': err})

                                        if self.__reg.verbose:
                                            self.__reg.log.write(f"+{self.__class_name__}.__testCsValueConsistencyInPkLoop() "
                                                                 f"++ ValueError  - {err}\n")

                                if True in onebond[d]:
                                    for d2 in range(num_dim):
                                        if onebond[d][d2]:
                                            chain_id2 = row[chain_id_names[d2]]
                                            seq_id2 = row[seq_id_names[d2]]
                                            comp_id2 = row[comp_id_names[d2]]
                                            atom_id2 = _atom_id2 = row[atom_id_names[d2]]

                                            if atom_id2 is not None:
                                                diff = len(atom_id) != len(atom_id2)
                                                _atom_id = '_' + (atom_id[1:-1] if atom_id.startswith('H') and diff
                                                                  else atom_id[1:])
                                                _atom_id2 = '_' + (atom_id2[1:-1] if atom_id2.startswith('H') and diff
                                                                   else atom_id2[1:])

                                            if chain_id2 in EMPTY_VALUE or seq_id2 in EMPTY_VALUE\
                                               or comp_id2 in EMPTY_VALUE or atom_id2 in EMPTY_VALUE\
                                               or (d < d2 and (chain_id2 != chain_id or seq_id2 != seq_id
                                                               or comp_id2 != comp_id or _atom_id2 != _atom_id)):

                                                # DAOTHER-7681, issue #2
                                                if d < d2 and chain_id2 == chain_id and seq_id2 == seq_id\
                                                   and comp_id2 == comp_id and _atom_id2 != _atom_id\
                                                   and self.__reg.ccU.updateChemCompDict(comp_id):
                                                    _atom_ids = self.__reg.dpV.getAtomIdList(comp_id, atom_id)
                                                    _atom_ids2 = self.__reg.dpV.getAtomIdList(comp_id, atom_id2)
                                                    if any(True for b in self.__reg.ccU.lastBondDictList
                                                           if ((b['atom_id_1'] in _atom_ids
                                                                and b['atom_id_2'] in _atom_ids2)
                                                               or (b['atom_id_1'] in _atom_ids2
                                                                   and b['atom_id_2'] in _atom_ids))):
                                                        continue

                                                err = f"[Check row of {id_tag} {row[id_tag]}] Coherence transfer type is onebond. "\
                                                    "However, assignment of spectral peak is inconsistent with the type, ("\
                                                    + self.__reg.dpV.getReducedAtomNotation(chain_id_names[d], chain_id,
                                                                                            seq_id_names[d], seq_id,
                                                                                            comp_id_names[d], comp_id,
                                                                                            atom_id_names[d], atom_id)\
                                                    + ") vs ("\
                                                    + self.__reg.dpV.getReducedAtomNotation(chain_id_names[d2], chain_id2,
                                                                                            seq_id_names[d2], seq_id2,
                                                                                            comp_id_names[d2], comp_id2,
                                                                                            atom_id_names[d2], atom_id2)\
                                                    + ")."

                                                self.__reg.report.error.appendDescription('invalid_data',
                                                                                          {'file_name': file_name,
                                                                                           'sf_framecode': sf_framecode,
                                                                                           'category': lp_category,
                                                                                           'description': err})

                                                if self.__reg.verbose:
                                                    self.__reg.log.write(f"+{self.__class_name__}.__testCsValueConsistencyInPkLoop() "  # noqa: E501, pylint: disable=line-too-long
                                                                         f"++ ValueError  - {err}\n")

                                if True in jcoupling[d]:
                                    for d2 in range(num_dim):
                                        if jcoupling[d][d2]:
                                            chain_id2 = row[chain_id_names[d2]]
                                            seq_id2 = row[seq_id_names[d2]]
                                            comp_id2 = row[comp_id_names[d2]]
                                            atom_id2 = row[atom_id_names[d2]]

                                            # DAOTHER-7389, issue #2
                                            if chain_id2 in EMPTY_VALUE or seq_id2 in EMPTY_VALUE\
                                               or comp_id2 in EMPTY_VALUE or atom_id2 in EMPTY_VALUE\
                                               or (d < d2 and (chain_id2 != chain_id or seq_id2 != seq_id
                                                               or comp_id2 != comp_id)):

                                                err = f"[Check row of {id_tag} {row[id_tag]}] "\
                                                    "Coherence transfer type is jcoupling. "\
                                                    "However, assignment of spectral peak is inconsistent with the type, ("\
                                                    + self.__reg.dpV.getReducedAtomNotation(chain_id_names[d], chain_id,
                                                                                            seq_id_names[d], seq_id,
                                                                                            comp_id_names[d], comp_id,
                                                                                            atom_id_names[d], atom_id)\
                                                    + ") vs ("\
                                                    + self.__reg.dpV.getReducedAtomNotation(chain_id_names[d2], chain_id2,
                                                                                            seq_id_names[d2], seq_id2,
                                                                                            comp_id_names[d2], comp_id2,
                                                                                            atom_id_names[d2], atom_id2)\
                                                    + ")."

                                                self.__reg.report.error.appendDescription('invalid_data',
                                                                                          {'file_name': file_name,
                                                                                           'sf_framecode': sf_framecode,
                                                                                           'category': lp_category,
                                                                                           'description': err})

                                                if self.__reg.verbose:
                                                    self.__reg.log.write(f"+{self.__class_name__}.__testCsValueConsistencyInPkLoop() "  # noqa: E501, pylint: disable=line-too-long
                                                                         f"++ ValueError  - {err}\n")

                                if True in relayed[d]:
                                    for d2 in range(num_dim):
                                        if relayed[d][d2]:
                                            chain_id2 = row[chain_id_names[d2]]
                                            seq_id2 = row[seq_id_names[d2]]
                                            comp_id2 = row[comp_id_names[d2]]
                                            atom_id2 = row[atom_id_names[d2]]

                                            if chain_id2 in EMPTY_VALUE or seq_id2 in EMPTY_VALUE\
                                               or comp_id2 in EMPTY_VALUE or atom_id2 in EMPTY_VALUE\
                                               or (d < d2 and (chain_id2 != chain_id or abs(seq_id2 - seq_id) > 1)):
                                                # DAOTHER-7389, issue #2

                                                err = f"[Check row of {id_tag} {row[id_tag]}] Coherence transfer type is relayed. "\
                                                    "However, assignment of spectral peak is inconsistent with the type, ("\
                                                    + self.__reg.dpV.getReducedAtomNotation(chain_id_names[d], chain_id,
                                                                                            seq_id_names[d], seq_id,
                                                                                            comp_id_names[d], comp_id,
                                                                                            atom_id_names[d], atom_id)\
                                                    + ") vs ("\
                                                    + self.__reg.dpV.getReducedAtomNotation(chain_id_names[d2], chain_id2,
                                                                                            seq_id_names[d2], seq_id2,
                                                                                            comp_id_names[d2], comp_id2,
                                                                                            atom_id_names[d2], atom_id2)\
                                                    + ")."

                                                self.__reg.report.error.appendDescription('invalid_data',
                                                                                          {'file_name': file_name,
                                                                                           'sf_framecode': sf_framecode,
                                                                                           'category': lp_category,
                                                                                           'description': err})

                                                if self.__reg.verbose:
                                                    self.__reg.log.write(f"+{self.__class_name__}.__testCsValueConsistencyInPkLoop() "  # noqa: E501, pylint: disable=line-too-long
                                                                         f"++ ValueError  - {err}\n")

                except Exception as e:  # pylint: disable=broad-exception-caught

                    self.__reg.report.error.appendDescription('internal_error',
                                                              f"+{self.__class_name__}.__testCsValueConsistencyInPkLoop() "
                                                              "++ Error  - " + str(e))

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__testCsValueConsistencyInPkLoop() "
                                             f"++ Error  - {str(e)}\n")

        return self.__reg.report.getTotalErrors() == __errors

    def __testCsValueConsistencyInPkAltLoop(self) -> bool:
        """ Perform consistency test on peak position and assignment of spectral peaks.
        """

        __errors = self.__reg.report.getTotalErrors()

        for fileListId in range(self.__reg.file_path_list_len):

            if fileListId >= len(self.__reg.star_data_type) or self.__reg.star_data_type[fileListId] != 'Entry':
                continue

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            if file_type == 'nef' or input_source_dic['content_subtype'] is None:
                continue

            content_subtype = 'spectral_peak_alt'

            if content_subtype not in input_source_dic['content_subtype']:
                continue

            sf_category = SF_CATEGORIES[file_type][content_subtype]
            lp_category = '_Assigned_peak_chem_shift'

            cs_item_names = ITEM_NAMES_IN_CS_LOOP[file_type]
            cs_chain_id_name = cs_item_names['chain_id']
            cs_seq_id_name = cs_item_names['seq_id']
            cs_comp_id_name = cs_item_names['comp_id']
            cs_atom_id_name = cs_item_names['atom_id']
            cs_value_name = cs_item_names['value']
            cs_error_name = cs_item_names['error']
            cs_atom_type = cs_item_names['atom_type']
            cs_iso_number = cs_item_names['isotope_number']

            for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                sf_framecode = get_first_sf_tag(sf, 'Sf_framecode')

                cs_data = None

                cs_file_name = file_name
                csFileListId = fileListId

                try:

                    cs_list = get_first_sf_tag(sf, CS_LIST_SF_TAG_NAME[file_type])
                    _cs_list_id = get_first_sf_tag(sf, 'ID')

                    try:

                        cs_data = next(lp['data'] for lp in self.__reg.lp_data['chem_shift']
                                       if lp['file_name'] == file_name and lp['sf_framecode'] == cs_list)

                    except StopIteration:

                        if cs_list not in EMPTY_VALUE:

                            if fileListId == 0:

                                err = "Assigned chemical shifts are required to verify the consistency of assigned peak lists. "\
                                    f"Referred {cs_list!r} saveframe containing the assigned chemical shift does not exist."

                                self.__reg.report.error.appendDescription('missing_mandatory_content',
                                                                          {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                           'description': err})

                                if self.__reg.verbose:
                                    self.__reg.log.write(f"+{self.__class_name__}.__testCsValueConsistencyInPkAltLoop() "
                                                         f"++ Error  - {err}\n")

                                continue

                            cs_input_source = self.__reg.report.input_sources[0]
                            cs_input_source_dic = cs_input_source.get()
                            cs_file_name = cs_input_source_dic['file_name']
                            csFileListId = 0

                except Exception:  # pylint: disable=broad-exception-caught
                    pass

                if cs_data is None:

                    try:

                        _cs_data = next(lp for lp in self.__reg.lp_data['chem_shift'] if lp['file_name'] == cs_file_name)

                    except StopIteration:
                        continue

                    cs_data = _cs_data['data']
                    cs_list = _cs_data['sf_framecode']

                    cs_sf = self.__reg.dpA.getSaveframeByName(csFileListId, cs_list)

                    if cs_sf is None:
                        continue

                    _cs_list_id = get_first_sf_tag(sf, 'ID')

                try:

                    _num_dim = get_first_sf_tag(sf, NUM_DIM_ITEMS[file_type])
                    num_dim = int(_num_dim)

                    if num_dim not in range(1, MAX_DIM_NUM_OF_SPECTRA):
                        raise ValueError()

                except ValueError:  # raised error already at __testIndexConsistency()
                    return False

                max_dim = num_dim + 1

                item_names = []
                for dim in range(1, max_dim):
                    _d = {}
                    for k, v in ITEM_NAMES_IN_PK_LOOP[file_type].items():
                        if '%s' in v:
                            v = v % dim
                        _d[k] = v
                    item_names.append(_d)

                chain_id_names, seq_id_names, comp_id_names, atom_id_names = [], [], [], []

                for i in range(num_dim):
                    chain_id_names.append(item_names[i]['chain_id'])
                    seq_id_names.append(item_names[i]['seq_id'])
                    comp_id_names.append(item_names[i]['comp_id'])
                    atom_id_names.append(item_names[i]['atom_id'])

                aux_data = next((lp['data'] for lp in self.__reg.aux_data[content_subtype]
                                 if lp['file_name'] == file_name and lp['sf_framecode'] == sf_framecode
                                 and lp['category'] == AUX_LP_CATEGORIES[file_type][content_subtype][0]), None)

                axis_codes, alt_axis_codes, abs_pk_pos, sp_widths = [], [], [], []

                if aux_data is not None and len(aux_data) > 0:
                    for i in range(1, max_dim):
                        for sp_dim in aux_data:
                            if sp_dim['ID'] != i:
                                continue
                            axis_codes.append(sp_dim['Axis_code'])
                            if 'Atom_isotope_number' in sp_dim and 'Atom_type' in sp_dim and sp_dim['Atom_type'] is not None:
                                alt_axis_codes.append(str(sp_dim['Atom_isotope_number']) + sp_dim['Atom_type'])
                            else:
                                alt_axis_codes.append(sp_dim['Axis_code'])
                            abs_pk_pos.append(False if 'Absolute_peak_positions' not in sp_dim
                                              else sp_dim['Absolute_peak_positions'])
                            sp_width = None if 'Sweep_width_units' not in sp_dim else sp_dim.get('Sweep_width')
                            if 'Sweep_width_units' in sp_dim and sp_dim['Sweep_width_units'] == 'Hz'\
                               and 'Spectrometer_frequency' in sp_dim and sp_width is not None:
                                sp_freq = sp_dim['Spectrometer_frequency']
                                if sp_freq not in EMPTY_VALUE:
                                    sp_width /= sp_freq
                            sp_widths.append(sp_width)
                            break
                else:
                    for i in range(num_dim):
                        axis_codes.append(None)
                        alt_axis_codes.append(None)
                        abs_pk_pos.append(False)
                        sp_widths.append(None)

                aux_data = next((lp['data'] for lp in self.__reg.aux_data[content_subtype]
                                 if lp['file_name'] == file_name and lp['sf_framecode'] == sf_framecode
                                 and lp['category'] == AUX_LP_CATEGORIES[file_type][content_subtype][1]), None)

                onebond = [[False] * num_dim for i in range(num_dim)]
                if aux_data is not None:
                    for sp_dim_trans in aux_data:
                        if sp_dim_trans['Type'] == 'onebond':
                            dim_1 = sp_dim_trans['Spectral_dim_ID_1']
                            dim_2 = sp_dim_trans['Spectral_dim_ID_2']
                            onebond[dim_1 - 1][dim_2 - 1] = True
                            onebond[dim_2 - 1][dim_1 - 1] = True

                jcoupling = [[False] * num_dim for i in range(num_dim)]
                if aux_data is not None:
                    for sp_dim_trans in aux_data:
                        if sp_dim_trans['Type'].startswith('j'):
                            dim_1 = sp_dim_trans['Spectral_dim_ID_1']
                            dim_2 = sp_dim_trans['Spectral_dim_ID_2']
                            jcoupling[dim_1 - 1][dim_2 - 1] = True
                            jcoupling[dim_2 - 1][dim_1 - 1] = True

                relayed = [[False] * num_dim for i in range(num_dim)]
                if aux_data is not None:
                    for sp_dim_trans in aux_data:
                        if sp_dim_trans['Type'].startswith('relayed'):
                            dim_1 = sp_dim_trans['Spectral_dim_ID_1']
                            dim_2 = sp_dim_trans['Spectral_dim_ID_2']
                            relayed[dim_1 - 1][dim_2 - 1] = True
                            relayed[dim_2 - 1][dim_1 - 1] = True

                pk_id_name = 'Peak_ID'
                dim_id_name = 'Spectral_dim_ID'
                set_id_name = 'Set_ID'

                try:

                    lp_data = next((lp['data'] for lp in self.__reg.aux_data[content_subtype]
                                    if lp['file_name'] == file_name and lp['sf_framecode'] == sf_framecode
                                    and lp['category'] == lp_category), None)

                    if lp_data is not None:

                        for row in lp_data:

                            if not (cs_chain_id_name in row and cs_seq_id_name in row
                                    and cs_comp_id_name in row and cs_atom_id_name in row):
                                continue

                            chain_id = row[cs_chain_id_name]
                            if chain_id in EMPTY_VALUE:
                                continue

                            seq_id = row[cs_seq_id_name]
                            if seq_id in EMPTY_VALUE:
                                continue

                            comp_id = row[cs_comp_id_name]
                            if comp_id in EMPTY_VALUE:
                                continue

                            atom_id = _atom_id = row[cs_atom_id_name]
                            if atom_id in EMPTY_VALUE:
                                continue

                            cs_list_id = row['Assigned_chem_shift_list_ID']

                            if cs_list_id != _cs_list_id:

                                for _cs_data in self.__reg.lp_data['chem_shift']:

                                    if _cs_data['file_name'] == file_name:

                                        cs_data = _cs_data['data']
                                        cs_list = _cs_data['sf_framecode']

                                        cs_sf = self.__reg.dpA.getSaveframeByName(csFileListId, cs_list)

                                        if cs_sf is None:
                                            continue

                                        _cs_list_id = get_first_sf_tag(sf, 'ID')

                                        if cs_list_id == _cs_list_id:
                                            break

                            pk_id = row[pk_id_name]
                            d = row[dim_id_name] - 1
                            set_id = row[set_id_name]

                            position = row[cs_value_name]

                            _atom_ids = self.__reg.dpV.getAtomIdList(comp_id, atom_id)

                            len_atom_id = len(_atom_ids)

                            if len_atom_id == 0:
                                atom_id_ = atom_id

                            elif len_atom_id == 1 and atom_id == _atom_ids[0]:
                                atom_id_ = atom_id

                            else:  # representative atom id
                                atom_id_ = _atom_ids[0]

                            cs_idx = -1

                            if file_type == 'nmr-star' and self.__reg.dpV.isNmrAtomName(comp_id, atom_id):
                                pass
                            else:
                                atom_id_ = atom_id

                            atom_ids_w_cs = [_row[cs_atom_id_name] for _row in cs_data
                                             if _row[cs_chain_id_name] == chain_id
                                             and _row[cs_seq_id_name] == seq_id
                                             and _row[cs_comp_id_name] == comp_id]

                            if atom_id_ in atom_ids_w_cs:
                                cs_idx = atom_ids_w_cs.index(atom_id_)

                            else:
                                for atom_id_w_cs in atom_ids_w_cs:
                                    _atom_id_w_cs = self.__reg.dpV.getAtomIdList(comp_id, atom_id_w_cs)
                                    if any(True for _atom_id_ in _atom_ids if _atom_id_ in _atom_id_w_cs):
                                        cs_idx = atom_ids_w_cs.index(atom_id_w_cs)
                                        break

                            if cs_idx == -1:

                                err = f"[Check row of {pk_id_name} {row[pk_id_name]}] Assignment of spectral peak "\
                                    + self.__reg.dpV.getReducedAtomNotation(chain_id_names[d], chain_id, seq_id_names[d], seq_id,
                                                                            comp_id_names[d], comp_id, atom_id_names[d], atom_id)\
                                    + f" was not found in assigned chemical shifts of {cs_list!r} saveframe."

                                self.__reg.report.warning.appendDescription('insufficient_data',
                                                                            {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                             'category': lp_category, 'description': err})

                                if self.__reg.verbose:
                                    self.__reg.log.write(f"+{self.__class_name__}.__testCsValueConsistencyInPkAltLoop() "
                                                         f"++ Warning  - {err}\n")

                            else:

                                cs_intra = [_row for _row in cs_data
                                            if _row[cs_chain_id_name] == chain_id
                                            and _row[cs_seq_id_name] == seq_id
                                            and _row[cs_comp_id_name] == comp_id]

                                cs = cs_intra[cs_idx]

                                value = cs[cs_value_name]
                                error = cs[cs_error_name]

                                if value in EMPTY_VALUE:
                                    continue

                                if error is None or error < 1.0e-3\
                                   or error * self.__reg.cs_diff_error_scaled_by_sigma > CS_UNCERT_MAX:
                                    error = CS_UNCERT_MAX
                                else:
                                    error *= self.__reg.cs_diff_error_scaled_by_sigma

                                if abs(position - value) > error:

                                    if d < num_dim and not abs_pk_pos[d] and sp_widths[d] is not None:
                                        if position < value:
                                            while position < value:
                                                position += sp_widths[d]
                                        elif position > value:
                                            while position > value:
                                                position -= sp_widths[d]

                                    if abs(position - value) > error and sp_widths[d] is not None:

                                        if CS_RANGE_MIN < sp_widths[d] < CS_RANGE_MAX:

                                            err = f"[Check row of {pk_id_name} {row[pk_id_name]}] Peak position of "\
                                                f"spectral peak {cs_value_name} {position} ("\
                                                + self.__reg.dpV.getReducedAtomNotation(cs_chain_id_name, chain_id,
                                                                                        cs_seq_id_name, seq_id,
                                                                                        cs_comp_id_name, comp_id,
                                                                                        cs_atom_id_name, atom_id)\
                                                + f") in {sf_framecode!r} saveframe is inconsistent "\
                                                f"with the assigned chemical shift value {value} "\
                                                f"(difference {position - value:.3f}, tolerance {error}) in {cs_list!r} saveframe."

                                            if error >= CS_UNCERT_MAX and not self.__reg.remediation_mode:

                                                self.__reg.report.error.appendDescription('invalid_data',
                                                                                          {'file_name': file_name,
                                                                                           'sf_framecode': sf_framecode,
                                                                                           'category': lp_category,
                                                                                           'description': err})

                                                if self.__reg.verbose:
                                                    self.__reg.log.write(f"+{self.__class_name__}.__testCsValueConsistencyInPkAltLoop() "  # noqa: E501, pylint: disable=line-too-long
                                                                         f"++ ValueError  - {err}\n")

                                            else:

                                                self.__reg.report.warning.appendDescription('unusual_chemical_shift',
                                                                                            {'file_name': file_name,
                                                                                             'sf_framecode': sf_framecode,
                                                                                             'category': lp_category,
                                                                                             'description': err})

                                                if self.__reg.verbose:
                                                    self.__reg.log.write(f"+{self.__class_name__}.__testCsValueConsistencyInPkAltLoop() "  # noqa: E501, pylint: disable=line-too-long
                                                                         f"++ Warning  - {err}\n")

                                axis_code = str(cs[cs_iso_number]) + cs[cs_atom_type]

                                if aux_data is not None and d < num_dim\
                                   and axis_code not in (axis_codes[d], alt_axis_codes[d]) and cs[cs_atom_type] != axis_codes[d]:

                                    err = f"[Check row of {pk_id_name} {row[pk_id_name]}] Assignment of spectral peak "\
                                        + self.__reg.dpV.getReducedAtomNotation(cs_chain_id_name, chain_id, cs_seq_id_name, seq_id,
                                                                                cs_comp_id_name, comp_id, cs_atom_id_name, atom_id)\
                                        + f" is inconsistent with axis code {axis_code} vs {axis_codes[d]}."

                                    self.__reg.report.error.appendDescription('invalid_data',
                                                                              {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                               'category': lp_category, 'description': err})

                                    if self.__reg.verbose:
                                        self.__reg.log.write(f"+{self.__class_name__}.__testCsValueConsistencyInPkAltLoop() "
                                                             f"++ ValueError  - {err}\n")

                            if d < num_dim and True in onebond[d]:
                                for d2 in range(num_dim):
                                    if onebond[d][d2]:

                                        try:
                                            _row = next(_row for _row in lp_data
                                                        if _row[pk_id_name] == pk_id
                                                        and _row[dim_id_name] - 1 == d2
                                                        and _row[set_id_name] is set_id)
                                        except StopIteration:
                                            continue

                                        chain_id2 = _row[cs_chain_id_name]
                                        seq_id2 = _row[cs_seq_id_name]
                                        comp_id2 = _row[cs_comp_id_name]
                                        atom_id2 = _atom_id2 = _row[cs_atom_id_name]

                                        if atom_id2 is not None:
                                            diff = len(atom_id) != len(atom_id2)
                                            _atom_id = '_' + (atom_id[1:-1] if atom_id.startswith('H') and diff
                                                              else atom_id[1:])
                                            _atom_id2 = '_' + (atom_id2[1:-1] if atom_id2.startswith('H') and diff
                                                               else atom_id2[1:])

                                        if chain_id2 in EMPTY_VALUE or seq_id2 in EMPTY_VALUE\
                                           or comp_id2 in EMPTY_VALUE or atom_id2 in EMPTY_VALUE\
                                           or (d < d2 and (chain_id2 != chain_id or seq_id2 != seq_id
                                                           or comp_id2 != comp_id or _atom_id2 != _atom_id)):

                                            # DAOTHER-7681, issue #2
                                            if d < d2 and chain_id2 == chain_id and seq_id2 == seq_id\
                                               and comp_id2 == comp_id and _atom_id2 != _atom_id\
                                               and self.__reg.ccU.updateChemCompDict(comp_id):
                                                _atom_ids = self.__reg.dpV.getAtomIdList(comp_id, atom_id)
                                                _atom_ids2 = self.__reg.dpV.getAtomIdList(comp_id, atom_id2)
                                                if any(True for b in self.__reg.ccU.lastBondDictList
                                                       if ((b['atom_id_1'] in _atom_ids
                                                            and b['atom_id_2'] in _atom_ids2)
                                                           or (b['atom_id_1'] in _atom_ids2
                                                               and b['atom_id_2'] in _atom_ids))):
                                                    continue

                                            err = f"[Check row of {pk_id_name} {row[pk_id_name]}] "\
                                                "Coherence transfer type is onebond. "\
                                                "However, assignment of spectral peak is inconsistent with the type, ("\
                                                + self.__reg.dpV.getReducedAtomNotation(chain_id_names[d], chain_id,
                                                                                        seq_id_names[d], seq_id,
                                                                                        comp_id_names[d], comp_id,
                                                                                        atom_id_names[d], atom_id)\
                                                + ") vs ("\
                                                + self.__reg.dpV.getReducedAtomNotation(chain_id_names[d2], chain_id2,
                                                                                        seq_id_names[d2], seq_id2,
                                                                                        comp_id_names[d2], comp_id2,
                                                                                        atom_id_names[d2], atom_id2)\
                                                + ")."

                                            self.__reg.report.error.appendDescription('invalid_data',
                                                                                      {'file_name': file_name,
                                                                                       'sf_framecode': sf_framecode,
                                                                                       'category': lp_category, 'description': err})

                                            if self.__reg.verbose:
                                                self.__reg.log.write(f"+{self.__class_name__}.__testCsValueConsistencyInPkAltLoop() "  # noqa: E501, pylint: disable=line-too-long
                                                                     f"++ ValueError  - {err}\n")

                            if d < num_dim and True in jcoupling[d]:
                                for d2 in range(num_dim):
                                    if jcoupling[d][d2]:

                                        try:
                                            _row = next(_row for _row in lp_data
                                                        if _row[pk_id_name] == pk_id
                                                        and _row[dim_id_name] - 1 == d2
                                                        and _row[set_id_name] is set_id)
                                        except StopIteration:
                                            continue

                                        chain_id2 = _row[cs_chain_id_name]
                                        seq_id2 = _row[cs_seq_id_name]
                                        comp_id2 = _row[cs_comp_id_name]
                                        atom_id2 = _row[cs_atom_id_name]

                                        if chain_id2 in EMPTY_VALUE or seq_id2 in EMPTY_VALUE\
                                           or comp_id2 in EMPTY_VALUE or atom_id2 in EMPTY_VALUE\
                                           or (d < d2 and (chain_id2 != chain_id or seq_id2 != seq_id
                                                           or comp_id2 != comp_id)):  # DAOTHER-7389, issue #2

                                            err = f"[Check row of {pk_id_name} {row[pk_id_name]}] "\
                                                "Coherence transfer type is jcoupling. "\
                                                "However, assignment of spectral peak is inconsistent with the type, ("\
                                                + self.__reg.dpV.getReducedAtomNotation(chain_id_names[d], chain_id,
                                                                                        seq_id_names[d], seq_id,
                                                                                        comp_id_names[d], comp_id,
                                                                                        atom_id_names[d], atom_id)\
                                                + ") vs ("\
                                                + self.__reg.dpV.getReducedAtomNotation(chain_id_names[d2], chain_id2,
                                                                                        seq_id_names[d2], seq_id2,
                                                                                        comp_id_names[d2], comp_id2,
                                                                                        atom_id_names[d2], atom_id2)\
                                                + ")."

                                            self.__reg.report.error.appendDescription('invalid_data',
                                                                                      {'file_name': file_name,
                                                                                       'sf_framecode': sf_framecode,
                                                                                       'category': lp_category, 'description': err})

                                            if self.__reg.verbose:
                                                self.__reg.log.write(f"+{self.__class_name__}.__testCsValueConsistencyInPkAltLoop() "  # noqa: E501, pylint: disable=line-too-long
                                                                     f"++ ValueError  - {err}\n")

                            if d < num_dim and True in relayed[d]:
                                for d2 in range(num_dim):
                                    if relayed[d][d2]:

                                        try:
                                            _row = next(_row for _row in lp_data
                                                        if _row[pk_id_name] == pk_id
                                                        and _row[dim_id_name] - 1 == d2
                                                        and _row[set_id_name] is set_id)
                                        except StopIteration:
                                            continue

                                        chain_id2 = _row[cs_chain_id_name]
                                        seq_id2 = _row[cs_seq_id_name]
                                        comp_id2 = _row[cs_comp_id_name]
                                        atom_id2 = _row[cs_atom_id_name]

                                        # DAOTHER-7389, issue #2
                                        if chain_id2 in EMPTY_VALUE or seq_id2 in EMPTY_VALUE\
                                           or comp_id2 in EMPTY_VALUE or atom_id2 in EMPTY_VALUE\
                                           or (d < d2 and (chain_id2 != chain_id or abs(seq_id2 - seq_id) > 1)):

                                            err = f"[Check row of {pk_id_name} {row[pk_id_name]}] "\
                                                "Coherence transfer type is relayed. "\
                                                "However, assignment of spectral peak is inconsistent with the type, ("\
                                                + self.__reg.dpV.getReducedAtomNotation(chain_id_names[d], chain_id,
                                                                                        seq_id_names[d], seq_id,
                                                                                        comp_id_names[d], comp_id,
                                                                                        atom_id_names[d], atom_id)\
                                                + ") vs ("\
                                                + self.__reg.dpV.getReducedAtomNotation(chain_id_names[d2], chain_id2,
                                                                                        seq_id_names[d2], seq_id2,
                                                                                        comp_id_names[d2], comp_id2,
                                                                                        atom_id_names[d2], atom_id2)\
                                                + ")."

                                            self.__reg.report.error.appendDescription('invalid_data',
                                                                                      {'file_name': file_name,
                                                                                       'sf_framecode': sf_framecode,
                                                                                       'category': lp_category, 'description': err})

                                            if self.__reg.verbose:
                                                self.__reg.log.write(f"+{self.__class_name__}.__testCsValueConsistencyInPkAltLoop() "  # noqa: E501, pylint: disable=line-too-long
                                                                     f"++ ValueError  - {err}\n")

                except Exception as e:  # pylint: disable=broad-exception-caught

                    self.__reg.report.error.appendDescription('internal_error',
                                                              f"+{self.__class_name__}.__testCsValueConsistencyInPkAltLoop() "
                                                              "++ Error  - " + str(e))

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__testCsValueConsistencyInPkAltLoop() "
                                             f"++ Error  - {str(e)}\n")

        return self.__reg.report.getTotalErrors() == __errors

    def __testRdcVector(self) -> bool:
        """ Perform consistency test on RDC bond vectors.
        """

        __errors = self.__reg.report.getTotalErrors()

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            if input_source_dic['content_subtype'] is None:
                continue

            content_subtype = 'rdc_restraint'

            if content_subtype not in input_source_dic['content_subtype']:
                continue

            sf_category = SF_CATEGORIES[file_type][content_subtype]
            lp_category = LP_CATEGORIES[file_type][content_subtype]

            if self.__reg.star_data_type[fileListId] == 'Loop':
                sf = self.__reg.star_data[fileListId]
                sf_framecode = ''

                self.__reg.dpV.testRdcVector(file_name, file_type, content_subtype, sf_framecode, lp_category)

            elif self.__reg.star_data_type[fileListId] == 'Saveframe':
                sf = self.__reg.star_data[fileListId]
                sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                self.__reg.dpV.testRdcVector(file_name, file_type, content_subtype, sf_framecode, lp_category)

            else:

                for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                    sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                    if not any(True for loop in sf.loops if loop.category == lp_category):
                        continue

                    self.__reg.dpV.testRdcVector(file_name, file_type, content_subtype, sf_framecode, lp_category)

        return self.__reg.report.getTotalErrors() == __errors

    def __testCoordCovalentBond(self) -> bool:
        """ Perform consistency test on covalent bonds.
        """

        __errors = self.__reg.report.getTotalErrors()

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            if input_source_dic['content_subtype'] is None:
                continue

            content_subtype = 'poly_seq'

            if content_subtype not in input_source_dic['content_subtype']:
                continue

            sf_category = SF_CATEGORIES[file_type][content_subtype]
            lp_category = AUX_LP_CATEGORIES[file_type][content_subtype][0]

            if self.__reg.star_data_type[fileListId] == 'Loop':
                sf = self.__reg.star_data[fileListId]
                sf_framecode = ''

                self.__reg.dpV.testCoordCovalentBond(file_name, file_type, content_subtype, sf_framecode, lp_category)

            elif self.__reg.star_data_type[fileListId] == 'Saveframe':
                sf = self.__reg.star_data[fileListId]
                sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                self.__reg.dpV.testCoordCovalentBond(file_name, file_type, content_subtype, sf_framecode, lp_category)

            else:

                for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                    sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                    if not any(True for loop in sf.loops if loop.category == lp_category):
                        continue

                    self.__reg.dpV.testCoordCovalentBond(file_name, file_type, content_subtype, sf_framecode, lp_category)

        return self.__reg.report.getTotalErrors() == __errors

    def __testResidueVariant(self) -> bool:
        """ Perform consistency test on residue variants.
        """

        src_id = self.__reg.report.getInputSourceIdOfCoord()

        if src_id < 0:
            return False

        cif_input_source = self.__reg.report.input_sources[src_id]
        cif_input_source_dic = cif_input_source.get()

        has_poly_seq = has_key_value(cif_input_source_dic, 'polymer_sequence')

        if not has_poly_seq:
            return False

        cif_poly_seq = cif_input_source_dic['polymer_sequence']

        __errors = self.__reg.report.getTotalErrors()

        for fileListId in range(self.__reg.file_path_list_len):

            nmr_input_source = self.__reg.report.input_sources[fileListId]
            nmr_input_source_dic = nmr_input_source.get()

            file_name = nmr_input_source_dic['file_name']
            file_type = nmr_input_source_dic['file_type']

            if nmr_input_source_dic['content_subtype'] is None:
                continue

            content_subtype = 'poly_seq'

            if content_subtype not in nmr_input_source_dic['content_subtype']:
                continue

            sf_category = SF_CATEGORIES[file_type][content_subtype]
            lp_category = AUX_LP_CATEGORIES[file_type][content_subtype][1]  # nef: _nef_sequence, nmr-star: _Entity_deleted_atom

            if lp_category not in self.__reg.lp_category_list:
                continue

            seq_align_dic = self.__reg.report.sequence_alignment.get()
            chain_assign_dic = self.__reg.report.chain_assignment.get()

            if 'nmr_poly_seq_vs_model_poly_seq' not in chain_assign_dic:

                err = "Chain assignment does not exist, __assignCoordPolymerSequence() should be invoked."

                self.__reg.report.error.appendDescription('internal_error',
                                                          f"+{self.__class_name__}.__testCoordResidueVariant() ++ Error  - " + err)

                if self.__reg.verbose:
                    self.__reg.log.write(f"+{self.__class_name__}.__testCoordResidueVariant() ++ Error  - {err}\n")

                continue

            if not has_key_value(seq_align_dic, 'nmr_poly_seq_vs_model_poly_seq'):
                continue

            if not has_key_value(chain_assign_dic, 'nmr_poly_seq_vs_model_poly_seq'):
                continue

            nmr2ca = {}

            for ca in chain_assign_dic['nmr_poly_seq_vs_model_poly_seq']:

                ref_chain_id = ca['ref_chain_id']
                test_chain_id = ca['test_chain_id']

                result = next((seq_align for seq_align in seq_align_dic['nmr_poly_seq_vs_model_poly_seq']
                               if seq_align['ref_chain_id'] == ref_chain_id and seq_align['test_chain_id'] == test_chain_id), None)

                if ref_chain_id not in nmr2ca:
                    nmr2ca[ref_chain_id] = []

                sa = {'seq_align': result}  # DAOTHER-7465

                if 'unmapped_sequence' in ca:
                    sa['seq_unmap'] = [unmapped['ref_seq_id'] for unmapped in ca['unmapped_sequence']]

                nmr2ca[ref_chain_id].append(sa)

            if self.__reg.star_data_type[fileListId] == 'Loop':
                sf = self.__reg.star_data[fileListId]
                sf_framecode = ''

                self.__reg.dpV.testResidueVariant(file_name, file_type, content_subtype,
                                                  sf, sf_framecode, lp_category, cif_poly_seq, nmr2ca)

            elif self.__reg.star_data_type[fileListId] == 'Saveframe':
                sf = self.__reg.star_data[fileListId]
                sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                self.__reg.dpV.testResidueVariant(file_name, file_type, content_subtype,
                                                  sf, sf_framecode, lp_category, cif_poly_seq, nmr2ca)

            else:

                for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                    sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                    if not any(True for loop in sf.loops if loop.category == lp_category):
                        continue

                    self.__reg.dpV.testResidueVariant(file_name, file_type, content_subtype,
                                                      sf, sf_framecode, lp_category, cif_poly_seq, nmr2ca)

        return self.__reg.report.getTotalErrors() == __errors

    def __retrieveCoordAssemblyChecker__(self) -> None:
        """ Wrapper function for ParserListenerUtil.coordAssemblyChecker.
        """

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        has_poly_seq = has_key_value(input_source_dic, 'polymer_sequence')

        nmrPolySeq = input_source_dic['polymer_sequence'] if has_poly_seq and self.__reg.bmrb_only and self.__reg.internal_mode\
            else None

        hash_code_ext = ''
        if nmrPolySeq is not None:
            hash_code_ext = f'_{hashlib.md5(str(nmrPolySeq).encode()).hexdigest()[:4]}'

        self.__reg.asmChkCachePath = None

        if self.__cifHashCode is not None:

            self.__reg.asmChkCachePath =\
                os.path.join(self.__reg.cacheDirPath, f"{self.__cifHashCode}{hash_code_ext}_asm_chk.pkl")
            self.__reg.coordPropCachePath =\
                os.path.join(self.__reg.cacheDirPath, f"{self.__cifHashCode}{hash_code_ext}_coord_prop.pkl")

            self.__reg.caC = load_from_pickle(self.__reg.asmChkCachePath)
            self.__reg.cpC = load_from_pickle(self.__reg.coordPropCachePath, default=copy.copy(DEFAULT_COORD_PROPERTIES))
            self.__reg.cpcHashCode = hash(str(self.__reg.cpC))

            # DAOTHER-8817
            if self.__reg.caC is not None and 'chem_comp_atom' in self.__reg.caC\
               and 'chem_comp_type' in self.__reg.caC\
               and 'auth_atom_name_to_id' in self.__reg.caC\
               and 'auth_atom_name_to_id_ext' in self.__reg.caC\
               and 'auth_to_star_seq_ann' in self.__reg.caC\
               and 'mod_residue' in self.__reg.caC\
               and 'split_ligand' in self.__reg.caC:
                self.__reg.nefT.set_chem_comp_dict(self.__reg.caC['chem_comp_atom'],
                                                   self.__reg.caC['chem_comp_bond'],
                                                   self.__reg.caC['chem_comp_topo'],
                                                   self.__reg.caC['auth_atom_name_to_id'],
                                                   self.__reg.caC['chem_comp_type'])
                return

        self.__parseCoordinate()  # need to set representative_model/alt_id values

        self.__reg.caC = coordAssemblyChecker(self.__reg.verbose, self.__reg.log,
                                              self.__reg.representative_model_id,
                                              self.__reg.representative_alt_id,
                                              self.__reg.cR, self.__reg.ccU, None, nmrPolySeq)

        if None not in (self.__reg.caC, self.__reg.asmChkCachePath):
            write_as_pickle(self.__reg.caC, self.__reg.asmChkCachePath)

        # DAOTHER-8817
        self.__reg.nefT.set_chem_comp_dict(self.__reg.caC['chem_comp_atom'],
                                           self.__reg.caC['chem_comp_bond'],
                                           self.__reg.caC['chem_comp_topo'],
                                           self.__reg.caC['auth_atom_name_to_id'],
                                           self.__reg.caC['chem_comp_type'])

    def __validateStrMr(self) -> bool:
        """ Validate restraints of NMR-STAR restraint files.
        """

        if self.__reg.release_mode:
            return True

        src_id = self.__reg.report.getInputSourceIdOfCoord()

        if src_id < 0:
            return False

        cif_input_source = self.__reg.report.input_sources[src_id]
        cif_input_source_dic = cif_input_source.get()

        has_poly_seq = has_key_value(cif_input_source_dic, 'polymer_sequence')

        if not has_poly_seq:
            return False

        if self.__reg.list_id_counter is None:
            self.__reg.list_id_counter = {}
        if self.__reg.mr_sf_dict_holder is None:
            self.__reg.mr_sf_dict_holder = {}

        # DAOTHER-8751, 8817 (D_1300043061)
        if self.__reg.combined_mode\
           and (not self.__reg.remediation_mode or self.__reg.annotation_mode or self.__reg.native_combined):

            if len(self.__reg.star_data) == 0 or self.__reg.star_data[0] is None or self.__reg.star_data_type[0] != 'Entry':
                return True

            master_entry = self.__reg.star_data[0]

            fileListId = 0

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_type = input_source_dic['file_type']
            content_subtype = input_source_dic['content_subtype']

            if file_type != 'nmr-star':
                return True

            file_name = input_source_dic['file_name']

            if input_source_dic['content_subtype'] is None:
                return True

            for content_subtype in self.__reg.mr_content_subtypes:

                if content_subtype not in input_source_dic['content_subtype']:
                    continue

                if content_subtype not in self.__reg.mr_sf_dict_holder:
                    self.__reg.mr_sf_dict_holder[content_subtype] = []

                sf_category = SF_CATEGORIES[file_type][content_subtype]
                lp_category = LP_CATEGORIES[file_type][content_subtype]

                for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                    sf_framecode = get_first_sf_tag(sf, 'Sf_framecode')
                    original_file_name = get_first_sf_tag(sf, 'Data_file_name')
                    if len(original_file_name) == 0:
                        original_file_name = file_name.replace('-corrected', '')
                        if 'original_file_name' in input_source_dic:
                            if input_source_dic['original_file_name'] is not None:
                                original_file_name = os.path.basename(input_source_dic['original_file_name'])

                    if self.__reg.dpV.validateStrMr(fileListId, file_type, original_file_name, content_subtype,
                                                    sf, sf_framecode, lp_category):
                        del master_entry[sf]

                        _sf = self.__reg.mr_sf_dict_holder[content_subtype][-1]['saveframe']

                        master_entry.add_saveframe(_sf)

            if self.__reg.dstPath is not None:

                # __validateStrPk() will do the same task later
                if not any(content_subtype in PK_CONTENT_SUBTYPES for content_type in input_source_dic['content_subtype']):

                    if not self.__reg.annotation_mode:

                        master_entry.write_to_file(self.__reg.dstPath,
                                                   show_comments=(self.__reg.bmrb_only and self.__reg.internal_mode),
                                                   skip_empty_loops=True, skip_empty_tags=False)

                        if NMR_CIF_FILE_PATH_KEY in self.__reg.outputParamDict:

                            try:

                                myIo = IoAdapterPy(False, sys.stderr)
                                containerList = myIo.readFile(self.__reg.dstPath)

                                if containerList is not None and len(containerList) > 1:

                                    if self.__reg.verbose:
                                        self.__reg.log.write("Input container list is "
                                                             f"{[(c.getName(), c.getType()) for c in containerList]!r}\n")

                                    eff_block_id = 1
                                    # if len(containerList[0].getObjNameList()) == 0 and not self.__reg.internal_mode else 0
                                    abandon_symbolic_labels(containerList)
                                    myIo.writeFile(self.__reg.outputParamDict[NMR_CIF_FILE_PATH_KEY],
                                                   containerList=containerList[eff_block_id:])

                            except Exception as e:  # pylint: disable=broad-exception-caught
                                self.__reg.log.write(f"+{self.__class_name__}.__validateStrMr() ++ Error  - {str(e)}\n")

            return True

        if self.__reg.combined_mode and self.__reg.remediation_mode:

            if len(self.__reg.star_data) == 0 or self.__reg.star_data[0] is None or self.__reg.star_data_type[0] != 'Entry':
                return True

            master_entry = self.__reg.star_data[0]

            fileListId = 0

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_type = input_source_dic['file_type']
            content_subtype = input_source_dic['content_subtype']

            if file_type != 'nmr-star':
                return True

            file_name = input_source_dic['file_name']

            if input_source_dic['content_subtype'] is None:
                return True

            content_subtype = 'dihed_restraint'

            if content_subtype in input_source_dic['content_subtype']:

                sf_category = SF_CATEGORIES[file_type][content_subtype]
                lp_category = LP_CATEGORIES[file_type][content_subtype]

                modified = False

                for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):

                    if not any(True for loop in sf.loops if loop.category == lp_category):
                        continue

                    modified |= self.__reg.dpR.remediateDihedLoop(file_type, sf.get_loop(lp_category))

                if modified and fileListId == 0:
                    self.__reg.lp_data[content_subtype].clear()
                    self.__testDataConsistencyInLoop(content_subtype)
                    self.__depositNmrData()

        if MR_FILE_PATH_LIST_KEY not in self.__reg.inputParamDict:
            return True

        for fileListId in range(self.__reg.cs_file_path_list_len, self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_type = input_source_dic['file_type']
            content_subtype = input_source_dic['content_subtype']

            if file_type != 'nmr-star':
                continue

            file_name = input_source_dic['file_name']

            original_file_name = file_name.replace('-corrected', '')
            if 'original_file_name' in input_source_dic:
                if input_source_dic['original_file_name'] is not None:
                    original_file_name = os.path.basename(input_source_dic['original_file_name'])

            if input_source_dic['content_subtype'] is None:
                continue

            for content_subtype in self.__reg.mr_content_subtypes:

                if content_subtype not in input_source_dic['content_subtype']:
                    continue

                if content_subtype not in self.__reg.mr_sf_dict_holder:
                    self.__reg.mr_sf_dict_holder[content_subtype] = []

                sf_category = SF_CATEGORIES[file_type][content_subtype]
                lp_category = LP_CATEGORIES[file_type][content_subtype]

                if self.__reg.star_data_type[fileListId] == 'Loop':
                    sf = self.__reg.star_data[fileListId]
                    sf_framecode = ''

                    self.__reg.dpV.validateStrMr(fileListId, file_type, original_file_name, content_subtype,
                                                 sf, sf_framecode, lp_category)

                elif self.__reg.star_data_type[fileListId] == 'Saveframe':
                    sf = self.__reg.star_data[fileListId]
                    sf_framecode = get_first_sf_tag(sf, 'Sf_framecode')

                    self.__reg.dpV.validateStrMr(fileListId, file_type, original_file_name, content_subtype,
                                                 sf, sf_framecode, lp_category)

                else:

                    for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                        sf_framecode = get_first_sf_tag(sf, 'Sf_framecode')

                        self.__reg.dpV.validateStrMr(fileListId, file_type, original_file_name, content_subtype,
                                                     sf, sf_framecode, lp_category)

        return True

    def __validateLegacyMr(self) -> bool:
        """ Validate data content of legacy restraint files.
        """

        return self.__reg.dpR.validateLegacyMr()

    def __validateLegacyPk(self) -> bool:
        """ Validate data content of legacy spectral peak files and merge them if possible.
        """

        return self.__reg.dpR.validateLegacyPk()

    def __validateLegacyCs(self) -> bool:
        """ Validate data content of legacy NMR chemical shift files and merge them if possible.
        """

        return self.__reg.dpR.validateLegacyCs()

    def __validateLegacyCsp(self) -> bool:
        """ Validate data content of legacy NMR chemical shift perturbation files and merge them if possible.
        """

        return self.__reg.dpR.validateLegacyCsp()

    def __validateSaxsMr(self) -> bool:
        """ Validate SAXS restraint files.
        """

        return self.__reg.dpR.validateSaxsMr()

    def __validateStrPk(self) -> bool:
        """ Validate spectral peak lists in NMR-STAR restraint files.
        """

        if self.__reg.release_mode:
            return True

        src_id = self.__reg.report.getInputSourceIdOfCoord()

        if src_id < 0:
            return False

        cif_input_source = self.__reg.report.input_sources[src_id]
        cif_input_source_dic = cif_input_source.get()

        has_poly_seq = has_key_value(cif_input_source_dic, 'polymer_sequence')

        if not has_poly_seq:
            return False

        list_id = 1

        # DAOTHER-8751, 8817 (D_1300043061)
        if self.__reg.combined_mode\
           and (not self.__reg.remediation_mode or self.__reg.annotation_mode or self.__reg.native_combined):

            if len(self.__reg.star_data) == 0 or not isinstance(self.__reg.star_data[0], pynmrstar.Entry):
                return True

            master_entry = self.__reg.star_data[0]

            fileListId = 0

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_type = input_source_dic['file_type']
            content_subtype = input_source_dic['content_subtype']

            if file_type != 'nmr-star':
                return True

            if input_source_dic['content_subtype'] is None:
                return True

            for content_subtype in PK_CONTENT_SUBTYPES:

                if content_subtype not in input_source_dic['content_subtype']:
                    continue

                sf_category = SF_CATEGORIES[file_type][content_subtype]
                lp_category = LP_CATEGORIES[file_type][content_subtype]

                for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                    sf_framecode = get_first_sf_tag(sf, 'Sf_framecode')

                    self.__reg.dpV.validateStrPk(fileListId, file_type, content_subtype, list_id, sf, sf_framecode, lp_category)

                    list_id += 1

            if list_id > 1 and self.__reg.dstPath is not None:

                if not self.__reg.annotation_mode:

                    master_entry.write_to_file(self.__reg.dstPath,
                                               show_comments=(self.__reg.bmrb_only and self.__reg.internal_mode),
                                               skip_empty_loops=True, skip_empty_tags=False)

                    if NMR_CIF_FILE_PATH_KEY in self.__reg.outputParamDict:

                        try:

                            myIo = IoAdapterPy(False, sys.stderr)
                            containerList = myIo.readFile(self.__reg.dstPath)

                            if containerList is not None and len(containerList) > 1:

                                if self.__reg.verbose:
                                    self.__reg.log.write("Input container list is "
                                                         f"{[(c.getName(), c.getType()) for c in containerList]!r}\n")

                                eff_block_id = 1
                                abandon_symbolic_labels(containerList)
                                myIo.writeFile(self.__reg.outputParamDict[NMR_CIF_FILE_PATH_KEY],
                                               containerList=containerList[eff_block_id:])

                        except Exception as e:  # pylint: disable=broad-exception-caught
                            self.__reg.log.write(f"+{self.__class_name__}.__validateStrPk() ++ Error  - {str(e)}\n")

            return True

        if MR_FILE_PATH_LIST_KEY not in self.__reg.inputParamDict:
            return True

        for fileListId in range(self.__reg.cs_file_path_list_len, self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_type = input_source_dic['file_type']
            content_subtype = input_source_dic['content_subtype']

            if file_type != 'nmr-star':
                continue

            file_name = input_source_dic['file_name']

            if input_source_dic['content_subtype'] is None:
                continue

            for content_subtype in PK_CONTENT_SUBTYPES:

                if content_subtype not in input_source_dic['content_subtype']:
                    continue

                sf_category = SF_CATEGORIES[file_type][content_subtype]
                lp_category = LP_CATEGORIES[file_type][content_subtype]

                if self.__reg.star_data_type[fileListId] == 'Loop':

                    err = f"Mandatory loops with categories {AUX_LP_CATEGORIES[file_type][content_subtype]} are missing. "\
                        f"Please re-upload the {file_type.upper()} file."

                    self.__reg.report.error.appendDescription('missing_data',
                                                              {'file_name': file_name, 'category': lp_category,
                                                               'description': err})

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__validateStrPk() ++ Error  - {err}\n")

                elif self.__reg.star_data_type[fileListId] == 'Saveframe':
                    sf = self.__reg.star_data[fileListId]
                    sf_framecode = get_first_sf_tag(sf, 'Sf_framecode')

                    self.__reg.dpV.validateStrPk(fileListId, file_type, content_subtype, list_id, sf, sf_framecode, lp_category)

                    list_id += 1

                else:

                    for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                        sf_framecode = get_first_sf_tag(sf, 'Sf_framecode')

                        self.__reg.dpV.validateStrPk(fileListId, file_type, content_subtype, list_id, sf, sf_framecode, lp_category)

                        list_id += 1

        return True

    def __calculateStatsOfExptlData(self) -> bool:
        """ Calculate statistics of experimental data.
        """

        __errors = self.__reg.report.getTotalErrors()

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            if input_source_dic['content_subtype'] is None:
                continue

            seq_align_dic = self.__reg.report.sequence_alignment.get()

            stats = {}

            for content_subtype in input_source_dic['content_subtype']:

                if content_subtype in ('entry_info', 'entity', 'ph_param_data'):
                    continue

                if self.__reg.report_prev is not None and content_subtype != 'chem_shift':
                    prev_stats = self.__reg.report_prev.getNmrLegacyStatsOfExptlData(fileListId, content_subtype)
                    if prev_stats is not None:
                        stats[content_subtype] = prev_stats
                        continue

                sf_category = SF_CATEGORIES[file_type][content_subtype]
                lp_category = LP_CATEGORIES[file_type][content_subtype]

                if content_subtype == 'poly_seq':
                    lp_category = AUX_LP_CATEGORIES[file_type][content_subtype][0]

                asm = []

                list_id = 1

                if self.__reg.star_data_type[fileListId] == 'Loop':
                    sf = self.__reg.star_data[fileListId]
                    sf_framecode = ''

                    self.__reg.dpV.calculateStatsOfExptlData(fileListId, file_name, file_type, content_subtype,
                                                             sf, list_id, sf_framecode, lp_category, seq_align_dic, asm)

                elif self.__reg.star_data_type[fileListId] == 'Saveframe':
                    sf = self.__reg.star_data[fileListId]
                    sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                    self.__reg.dpV.calculateStatsOfExptlData(fileListId, file_name, file_type, content_subtype,
                                                             sf, list_id, sf_framecode, lp_category, seq_align_dic, asm)

                else:

                    for list_id, sf in enumerate(self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category), start=1):
                        sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                        if not any(True for loop in sf.loops if loop.category == lp_category):
                            if content_subtype != 'spectral_peak':
                                continue

                        self.__reg.dpV.calculateStatsOfExptlData(fileListId, file_name, file_type, content_subtype,
                                                                 sf, list_id, sf_framecode, lp_category, seq_align_dic, asm)

                if len(asm) > 0:
                    stats[content_subtype] = asm

            input_source.setItemValue('stats_of_exptl_data', stats)

            if self.__reg.op == 'nmr-str-consistency-check' and self.__reg.combined_mode\
               and REPORT_FILE_PATH_KEY in self.__reg.inputParamDict:
                fPath = self.__reg.inputParamDict[REPORT_FILE_PATH_KEY]

                if os.access(fPath, os.F_OK) and os.path.getsize(fPath) > 0:
                    self.__reg.dpV.calculateOutputStats()

        return self.__reg.report.getTotalErrors() == __errors

    def __detectDimTransferTypeViaThroughSpace(self) -> bool:
        """ Detect spectral peak transfer type via through-space.
        """

        if len(self.__reg.star_data) == 0 or not isinstance(self.__reg.star_data[0], pynmrstar.Entry):
            return False

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_name = input_source_dic['file_name']
        file_type = input_source_dic['file_type']

        content_subtype = 'spectral_peak'

        sf_category = SF_CATEGORIES[file_type][content_subtype]
        lp_category = LP_CATEGORIES[file_type][content_subtype]
        pk_char_category = '_Peak_char'

        sf_list = self.__reg.star_data[0].get_saveframes_by_category(sf_category)

        if len(sf_list) == 0:
            return False  # raised proper recommendation already

        solid_state_nmr = self.__reg.exptl_method == 'SOLID-STATE NMR'

        exp_classes = []

        for sf in sf_list:

            data_file_name = ''

            # natively combined nmr_data file needs to be investigated
            if not self.__reg.native_combined or self.__reg.internal_mode:

                if file_type == 'nmr-star':
                    exp_class = get_first_sf_tag(sf, 'Experiment_class')

                    if exp_class.endswith('through-space') or exp_class.endswith('through-space?'):
                        return True

                    data_file_name = get_first_sf_tag(sf, 'Data_file_name')

                # '_nef_spectrum_dimension_transfer' or '_Spectral_dim_transfer'
                aux_lp_category = AUX_LP_CATEGORIES[file_type][content_subtype][1]

                try:

                    aux_loop = sf.get_loop(aux_lp_category)

                    dat = aux_loop.get_tag(['Type' if file_type == 'nmr-star' else 'transfer_type'])

                    for row in dat:
                        if row in ('through-space', 'through-space?'):
                            return True

                except KeyError:
                    pass

            try:

                _num_dim = get_first_sf_tag(sf, NUM_DIM_ITEMS[file_type])
                num_dim = int(_num_dim)

                if num_dim not in range(1, MAX_DIM_NUM_OF_SPECTRA):
                    raise ValueError()

            except ValueError:
                continue

            if num_dim > 5:
                continue

            max_dim = num_dim + 1

            cur_spectral_dim = {}

            try:

                loop = sf.get_loop(lp_category)

                key_items = []
                for dim in range(1, max_dim):
                    for k in PK_KEY_ITEMS[file_type]:
                        if k['type'] == 'float':  # position
                            _k = copy.copy(k)
                            if '%s' in k['name']:
                                _k['name'] = k['name'] % dim
                            key_items.append(_k['name'])
                    cur_spectral_dim[dim] = copy.copy(SPECTRAL_DIM_TEMPLATE)
                    cur_spectral_dim[dim]['freq_hint'] = []

                dat = loop.get_tag(key_items)

                for row in dat:

                    for dim in range(1, max_dim):
                        freq = row[dim - 1]
                        if freq in EMPTY_VALUE:
                            continue
                        if isinstance(freq, str):
                            freq = float(freq)
                        cur_spectral_dim[dim]['freq_hint'].append(freq)

                exp_class = guess_primary_dim_transfer_type(solid_state_nmr, data_file_name, num_dim, cur_spectral_dim)

                if exp_class in ('through-space', 'through-space?'):
                    return True

                if exp_class not in EMPTY_VALUE:
                    exp_class = f'{exp_class!r}'

                    if exp_class not in exp_classes:
                        exp_classes.append(exp_class)

            except KeyError:

                if file_type == 'nef':
                    continue

                try:

                    loop = sf.get_loop(pk_char_category)

                    dat = loop.get_tag(['Spectral_dim_ID', 'Chem_shift_val'])

                    for row in dat:
                        dim, freq = row
                        if isinstance(dim, str):
                            try:
                                dim = int(dim)
                            except ValueError:
                                continue
                        if freq in EMPTY_VALUE:
                            continue
                        if isinstance(freq, str):
                            freq = float(freq)
                        cur_spectral_dim[dim]['freq_hint'].append(freq)

                    exp_class = guess_primary_dim_transfer_type(solid_state_nmr, data_file_name, num_dim, cur_spectral_dim)

                    if exp_class in ('through-space', 'through-space?'):
                        return True

                    if exp_class not in EMPTY_VALUE:
                        exp_class = f'{exp_class!r}'

                        if exp_class not in exp_classes:
                            exp_classes.append(exp_class)

                except KeyError:
                    continue

        primary_spectra_for_structure_determination =\
            'NOESY or ROESY' if not solid_state_nmr else 'DARR, REDOR, TEDOR or RFDR'

        hint = f" classified into {', '.join(exp_classes)}, respectively" if len(exp_classes) > 0 else ''

        warn = f"You have uploaded spectral peak list(s){hint}. "\
            "However, none of spectral peak list(s) appears to be derived from expected spectra "\
            f"such as the {primary_spectra_for_structure_determination}. "\
            "The wwPDB NMR Validation Task Force strongly encourages the submission of spectral peak lists, "\
            f"in particular those generated from the {primary_spectra_for_structure_determination} spectra."

        self.__reg.report.warning.appendDescription('encouragement',
                                                    {'file_name': file_name, 'description': warn})

        if self.__reg.verbose:
            self.__reg.log.write(f"+{self.__class_name__}.__detectDimTransferTypeViaThroughSpace() ++ Warning  - {warn}\n")

        return False

    def __validateCoordInputSource(self) -> bool:
        """ Validate coordinate file as secondary input resource.
        """

        file_type = 'pdbx'

        content_type = CONTENT_TYPE[file_type]

        if self.__parseCoordinate():

            self.__reg.report.appendInputSource()

            input_source = self.__reg.report.input_sources[-1]

            input_source.setItemValue('file_name', os.path.basename(self.__reg.cifPath))
            input_source.setItemValue('file_type', file_type)
            input_source.setItemValue('content_type', content_type)
            input_source.setItemValue('ignore_error', False)

            return True

        if self.__reg.entry_id == INITIAL_ENTRY_ID:
            self.__reg.entry_id = DEFAULT_ENTRY_ID

        return False

    def __parseCoordinate(self) -> bool:
        """ Parse coordinate file.
        """

        if self.__reg.cifChecked:
            return True

        file_type = 'pdbx'

        if not self.__parseCoordFilePath():

            if MODEL_FILE_PATH_KEY in self.__reg.inputParamDict:

                err = f"No such {self.__reg.inputParamDict[MODEL_FILE_PATH_KEY]!r} file."

                self.__reg.report.error.appendDescription('internal_error',
                                                          f"+{self.__class_name__}.__parseCoordinate() ++ Error  - " + err)

                if self.__reg.verbose:
                    self.__reg.log.write(f"+{self.__class_name__}.__parseCoordinate() ++ Error  - {err}\n")

            elif not self.__reg.bmrb_only:

                err = f"{READABLE_FILE_TYPE[file_type]} formatted coordinate file is mandatory."

                self.__reg.report.error.appendDescription('internal_error',
                                                          f"+{self.__class_name__}.__parseCoordinate() ++ Error  - " + err)

                if self.__reg.verbose:
                    self.__reg.log.write(f"+{self.__class_name__}.__parseCoordinate() ++ Error  - {err}\n")

            return False

        file_name = os.path.basename(self.__reg.cifPath)

        try:

            if self.__reg.cifPath is None:

                err = f"{file_name!r} is invalid {READABLE_FILE_TYPE[file_type]} file."

                self.__reg.report.error.appendDescription('internal_error',
                                                          f"+{self.__class_name__}.__parseCoordinate() ++ Error  - " + err)

                if self.__reg.verbose:
                    self.__reg.log.write(f"+{self.__class_name__}.__parseCoordinate() ++ Error  - {err}\n")

                return False

            if self.__reg.entry_id == INITIAL_ENTRY_ID:
                entry = self.__reg.cR.getDictList('entry')

                if len(entry) == 0 or ('id' not in entry[0]):
                    self.__reg.entry_id = DEFAULT_ENTRY_ID
                else:
                    # DAOTHER-9511: replace white space in a datablock name to underscore
                    self.__reg.entry_id = entry[0]['id'].strip().replace(' ', '_')

            exptl = self.__reg.cR.getDictList('exptl')

            if len(exptl) > 0 and 'method' in exptl[0]:
                self.__reg.exptl_method = exptl[0]['method']

            self.__reg.total_models = 0
            self.__reg.eff_model_ids.clear()
            self.__reg.representative_model_id = REPRESENTATIVE_MODEL_ID
            self.__reg.representative_alt_id = REPRESENTATIVE_ALT_ID

            ensemble = self.__reg.cR.getDictList('pdbx_nmr_ensemble')

            if self.__reg.trust_pdbx_nmr_ens and len(ensemble) > 0 and 'conformers_submitted_total_number' in ensemble[0]:

                try:
                    self.__reg.total_models = int(ensemble[0]['conformers_submitted_total_number'])
                except ValueError:
                    pass

            if len(ensemble) == 0 or self.__reg.total_models == 0:

                ensemble = self.__reg.cR.getDictList('rcsb_nmr_ensemble')

                if self.__reg.trust_pdbx_nmr_ens and len(ensemble) > 0 and 'conformers_submitted_total_number' in ensemble[0]:

                    try:
                        self.__reg.total_models = int(ensemble[0]['conformers_submitted_total_number'])
                    except ValueError:
                        pass

                else:

                    try:

                        model_num_name = 'pdbx_PDB_model_num' if 'pdbx_PDB_model_num' in self.__reg.coord_atom_site_tags\
                            else 'ndb_model'

                        model_ids = self.__reg.cR.getDistinctValues('atom_site', model_num_name, 'int')

                        if len(model_ids) > 0:
                            self.__reg.representative_model_id = min(model_ids)
                            self.__reg.total_models = len(model_ids)
                            self.__reg.eff_model_ids = sorted(model_ids)

                    except Exception as e:  # pylint: disable=broad-exception-caught

                        self.__reg.report.error.appendDescription('internal_error',
                                                                  f"+{self.__class_name__}.__parseCoordinate() "
                                                                  "++ Error  - " + str(e))

            if self.__reg.trust_pdbx_nmr_ens and len(ensemble) > 0 and 'representative_conformer' in ensemble[0]:

                try:

                    rep_model_id = int(ensemble[0]['representative_conformer'])

                    if 1 <= rep_model_id <= self.__reg.total_models:
                        self.__reg.representative_model_id = rep_model_id

                except ValueError:
                    pass

            if len(self.__reg.eff_model_ids) == 0:

                try:

                    model_num_name = 'pdbx_PDB_model_num' if 'pdbx_PDB_model_num' in self.__reg.coord_atom_site_tags\
                        else 'ndb_model'

                    model_ids = self.__reg.cR.getDistinctValues('atom_site', model_num_name, 'int')

                    if len(model_ids) > 0:
                        self.__reg.total_models = len(model_ids)
                        self.__reg.eff_model_ids = sorted(model_ids)

                except Exception as e:  # pylint: disable=broad-exception-caught

                    self.__reg.report.error.appendDescription('internal_error',
                                                              f"+{self.__class_name__}.__parseCoordinate() ++ Error  - " + str(e))

            if self.__reg.total_models < 2:

                if not self.__reg.remediation_mode:

                    models = 'no' if self.__reg.total_models == 0\
                        else ('only one' if self.__reg.total_models == 1 else self.__reg.total_models)

                    err = f"Coordinate file has {models} model(s). "\
                        "Deposition of minimized average structure must be accompanied with ensemble "\
                        "and must be homogeneous with the ensemble."

                    self.__reg.report.error.appendDescription('missing_mandatory_content',
                                                              {'file_name': file_name, 'description': err})

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__parseCoordinate() ++ Error  - {err}\n")

            elif self.__reg.total_models < 5:

                if not self.__reg.remediation_mode:

                    warn = f"Coordinate file has {self.__reg.total_models} models. "\
                        "We encourage you to deposit a sufficient number of models in the ensemble."

                    self.__reg.report.warning.appendDescription('encouragement',
                                                                {'file_name': file_name, 'description': warn})

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__parseCoordinate() ++ Warning  - {warn}\n")

            if self.__reg.cR.hasItem('atom_site', 'label_alt_id'):
                alt_id = self.__reg.cR.getFirstValue('atom_site', 'label_alt_id')

                if alt_id is not None:
                    self.__reg.representative_alt_id = alt_id

            self.__ensemble_composition = {'total_models': self.__reg.total_models,
                                           'eff_model_ids': self.__reg.eff_model_ids,
                                           'representative_model_id': self.__reg.representative_model_id,
                                           'selection_criteria': 'unknown'
                                           }

            nmr_representative = self.__reg.cR.getDictList('pdbx_nmr_representative')

            if len(nmr_representative) > 0:

                try:
                    self.__reg.representative_model_id = int(nmr_representative[0]['conformer_id'])
                    self.__ensemble_composition['representative_model_id'] = self.__reg.representative_model_id
                    if 'selection_criteria' in nmr_representative[0]\
                       and nmr_representative[0]['selection_criteria'] not in EMPTY_VALUE:
                        self.__ensemble_composition['selection_criteria'] = nmr_representative[0]['selection_criteria']
                except ValueError:
                    pass

            self.__reg.recvd_nmr_constraints = False
            if self.__reg.cR.hasItem('pdbx_database_status', 'recvd_nmr_constraints'):
                pdbx_database_status = self.__reg.cR.getDictList('pdbx_database_status')
                self.__reg.recvd_nmr_constraints = pdbx_database_status[0]['recvd_nmr_constraints'] == 'Y'
                if not self.__reg.recvd_nmr_constraints and self.__reg.cR.hasItem('pdbx_database_status', 'date_nmr_constraints'):
                    date_nmr_constraints = pdbx_database_status[0]['date_nmr_constraints']
                    if date_nmr_constraints not in EMPTY_VALUE:
                        self.__reg.recvd_nmr_constraints = True
            if self.__reg.recvd_nmr_constraints:
                self.__reg.remediation_mode = True

            self.__reg.recvd_nmr_data = False
            if self.__reg.cR.hasItem('pdbx_database_status', 'recvd_nmr_data'):
                pdbx_database_status = self.__reg.cR.getDictList('pdbx_database_status')
                self.__reg.recvd_nmr_data = pdbx_database_status[0]['recvd_nmr_data'] == 'Y'

            self.__reg.versioned_atom_name_mapping = None
            self.__reg.internal_atom_name_mapping.clear()

            maxitpath = os.getenv('MAXITPATH')
            if maxitpath is None:
                package_dir = os.getenv('PACKAGE_DIR')
                maxitpath = os.path.join(package_dir, 'maxit', 'bin', 'maxit') if package_dir is not None else 'maxit'

            def has_coordinates(file_name):
                with open(file_name, 'r', encoding='utf-8', errors='ignore') as ifh:
                    for line in ifh:
                        if line.startswith('ATOM ') and line.count('.') >= 3:
                            return True
                return False

            def has_cif_coordinates(file_name):
                has_atom_site = False
                with open(file_name, 'r', encoding='utf-8', errors='ignore') as ifh:
                    for line in ifh:
                        if not has_atom_site:
                            if line.startswith('_atom_site.'):
                                has_atom_site = True
                            continue
                        if line.startswith('ATOM ') and line.count('.') >= 3:
                            return True
                return False

            if self.__reg.internal_mode and self.__reg.cR.hasCategory('database_2')\
               and self.__reg.cR.hasCategory('pdbx_audit_revision_history'):
                extended_pdb_id = None
                if self.__reg.cR.hasItem('database_2', 'pdbx_database_accession'):
                    database_2 = self.__reg.cR.getDictListWithFilter('database_2',
                                                                     [{'name': 'pdbx_database_accession', 'type': 'str'}],
                                                                     [{'name': 'database_id', 'type': 'str', 'value': 'PDB'}])
                    if len(database_2) > 0:
                        extended_pdb_id = database_2[0]['pdbx_database_accession']
                        if extended_pdb_id is not None and not PDB_ID_PAT.match(extended_pdb_id):
                            extended_pdb_id = None
                if extended_pdb_id is None:
                    if self.__reg.cR.hasCategory('entry'):
                        entry = self.__reg.cR.getDictList('entry')
                        if len(entry) > 0 and 'id' in entry[0]:
                            entry_id = entry[0]['id']
                            if len(entry_id) == 4 and PDB_ID_PAT.match(entry_id):
                                extended_pdb_id = f"pdb_0000{entry[0]['id'].lower()}"

                revision_history = {}
                for r in self.__reg.cR.getDictList('pdbx_audit_revision_history'):
                    major_revision = int(r['major_revision'])
                    minor_revision = int(r['minor_revision'])
                    if major_revision not in revision_history or minor_revision > revision_history[major_revision]:
                        revision_history[major_revision] = minor_revision

                if extended_pdb_id is not None and self.__reg.cR.hasCategory('chem_comp')\
                   and any(d['id'] not in EMPTY_VALUE and d['id'] not in STD_MON_DICT
                           for d in self.__reg.cR.getDictList('chem_comp')):

                    if len(revision_history) > 1:
                        self.__reg.versioned_atom_name_mapping =\
                            retrieveAtomNameMappingFromRevisions(self.__reg.cR, self.__reg.cacheDirPath,
                                                                 extended_pdb_id, revision_history,
                                                                 self.__reg.representative_model_id,
                                                                 self.__reg.representative_alt_id,
                                                                 self.__reg.csStat)

                    internal_cif_file =\
                        os.path.join(self.__reg.cR.getDirPath(), f'{extended_pdb_id[-4:]}_model-upload_P1.cif.V1')
                    internal_cif_file0 =\
                        os.path.join(self.__reg.cR.getDirPath(), f'{extended_pdb_id[-4:]}_model-upload_P1.cif.V0')
                    internal_cif_file_prefix =\
                        os.path.join(self.__reg.cR.getDirPath(), f'{extended_pdb_id[-4:]}_model-upload_P1.cif.V')
                    other_internal_cif_file =\
                        os.path.join(self.__reg.cR.getDirPath(), f'{extended_pdb_id[-4:]}_model-upload_P1.cif.V-1')

                    max_ver = 1
                    if not os.path.exists(internal_cif_file):

                        try:

                            import subprocess  # pylint: disable=import-outside-toplevel

                            database_2 = self.__reg.cR.getDictListWithFilter('database_2',
                                                                             [{'name': 'database_code', 'type': 'str'}],
                                                                             [{'name': 'database_id', 'type': 'str',
                                                                               'value': 'BMRB'}])
                            if len(database_2) > 0:

                                bmrb_id = database_2[0]['database_code']
                                if bmrb_id is not None and bmrb_id.isdigit():
                                    ret_code = -1
                                    intnl_upload_dir = os.path.join(self.__reg.cR.getDirPath(), f'bmr{bmrb_id}', 'work', 'upload')
                                    if os.path.isdir(intnl_upload_dir):
                                        for file_name in os.listdir(intnl_upload_dir):
                                            file_path = os.path.join(intnl_upload_dir, file_name)
                                            if not os.path.isfile(file_path):
                                                continue
                                            if file_name.endswith('.pdb'):
                                                if not has_coordinates(file_path):
                                                    continue
                                                com = [maxitpath, '-input', file_path,
                                                       '-output', f'{internal_cif_file_prefix}{max_ver}',
                                                       '-o', '1']
                                                result = subprocess.run(com, check=False)
                                                ret_code = result.returncode
                                                if ret_code == 0:
                                                    max_ver += 1
                                                    # break
                                            elif file_name.endswith('.cif'):
                                                if not has_cif_coordinates(file_path):
                                                    continue
                                                com = [maxitpath, '-input', file_path,
                                                       '-output', f'{internal_cif_file_prefix}{max_ver}',
                                                       '-o', '8']
                                                result = subprocess.run(com, check=False)
                                                ret_code = result.returncode
                                                if ret_code == 0:
                                                    max_ver += 1
                                                    # break
                                    if ret_code != 0:
                                        intnl_upload_dir = os.path.join(self.__reg.cR.getDirPath(), f'bmr{bmrb_id}', 'work')
                                        if os.path.isdir(intnl_upload_dir):
                                            for file_name in os.listdir(intnl_upload_dir):
                                                file_path = os.path.join(intnl_upload_dir, file_name)
                                                if not os.path.isfile(file_path):
                                                    continue
                                                if file_name.endswith('.pdb'):
                                                    if not has_coordinates(file_path):
                                                        continue
                                                    com = [maxitpath, '-input', file_path, '-output', internal_cif_file, '-o', '1']
                                                    result = subprocess.run(com, check=False)
                                                    ret_code = result.returncode
                                                    if ret_code == 0:
                                                        break
                                                elif file_name.endswith('.cif'):
                                                    if not has_cif_coordinates(file_path):
                                                        continue
                                                    com = [maxitpath, '-input', file_path, '-output', internal_cif_file, '-o', '8']
                                                    result = subprocess.run(com, check=False)
                                                    ret_code = result.returncode
                                                    if ret_code == 0:
                                                        break
                                    if ret_code != 0:
                                        intnl_upload_dir = os.path.join(self.__reg.cR.getDirPath(), f'bmr{bmrb_id}', 'work', 'data')
                                        if os.path.isdir(intnl_upload_dir):
                                            for file_name in os.listdir(intnl_upload_dir):
                                                if file_name.endswith('model-upload_P1.pdb.V1'):
                                                    file_path = os.path.join(intnl_upload_dir, file_name)
                                                    com = [maxitpath, '-input', file_path, '-output', internal_cif_file, '-o', '1']
                                                    result = subprocess.run(com, check=False)
                                                    ret_code = result.returncode
                                                    if ret_code == 0:
                                                        break
                                            if ret_code != 0:
                                                for file_name in os.listdir(intnl_upload_dir):
                                                    if file_name.endswith('model-upload_P1.cif.V1'):
                                                        file_path = os.path.join(intnl_upload_dir, file_name)
                                                        com = [maxitpath, '-input', file_path, '-output', internal_cif_file, '-o', '8']  # noqa: E501, pylint: disable=line-too-long
                                                        result = subprocess.run(com, check=False)
                                                        ret_code = result.returncode
                                                        if ret_code == 0:
                                                            break

                        except ImportError:
                            pass
                        except Exception:  # pylint: disable=broad-exception-caught
                            pass

                    if not os.path.exists(internal_cif_file0):

                        try:

                            database_2 = self.__reg.cR.getDictListWithFilter('database_2',
                                                                             [{'name': 'database_code', 'type': 'str'}],
                                                                             [{'name': 'database_id', 'type': 'str',
                                                                               'value': 'BMRB'}])
                            if len(database_2) > 0:

                                bmrb_id = database_2[0]['database_code']
                                if bmrb_id is not None and bmrb_id.isdigit():
                                    ret_code = -1
                                    intnl_upload_dir = os.path.join(self.__reg.cR.getDirPath(), f'bmr{bmrb_id}', 'work', 'upload')
                                    if os.path.isdir(intnl_upload_dir):
                                        for file_name in os.listdir(intnl_upload_dir):
                                            file_path = os.path.join(intnl_upload_dir, file_name)
                                            if not os.path.isfile(file_path):
                                                continue
                                            if file_name.endswith('.cif'):
                                                if not has_cif_coordinates(file_path):
                                                    continue
                                                ret_path = shutil.copyfile(file_path, internal_cif_file0)
                                                if os.path.exists(ret_path):
                                                    ret_code = 0
                                                    break
                                    if ret_code != 0:
                                        intnl_upload_dir = os.path.join(self.__reg.cR.getDirPath(), f'bmr{bmrb_id}', 'work')
                                        if os.path.isdir(intnl_upload_dir):
                                            for file_name in os.listdir(intnl_upload_dir):
                                                file_path = os.path.join(intnl_upload_dir, file_name)
                                                if not os.path.isfile(file_path):
                                                    continue
                                                if file_name.endswith('.cif'):
                                                    if not has_cif_coordinates(file_path):
                                                        continue
                                                    ret_path = shutil.copyfile(file_path, internal_cif_file0)
                                                    if os.path.exists(ret_path):
                                                        ret_code = 0
                                    if ret_code != 0:
                                        intnl_upload_dir = os.path.join(self.__reg.cR.getDirPath(), f'bmr{bmrb_id}', 'work', 'data')
                                        if os.path.isdir(intnl_upload_dir):
                                            for file_name in os.listdir(intnl_upload_dir):
                                                if file_name.endswith('model-upload_P1.cif.V1'):
                                                    file_path = os.path.join(intnl_upload_dir, file_name)
                                                    ret_path = shutil.copyfile(file_path, internal_cif_file0)
                                                    if os.path.exists(ret_path):
                                                        ret_code = 0
                                                        break

                        except Exception:  # pylint: disable=broad-exception-caught
                            pass

                    if os.path.exists(internal_cif_file):
                        for ver in range(1, 10):
                            if not os.path.exists(f'{internal_cif_file_prefix}{ver}'):
                                break
                            max_ver = ver + 1
                        for ver in range(1, max_ver):
                            major = max(revision_history)
                            minor = revision_history[major] + ver - 1
                            revision_history[major] = minor
                            self.__reg.internal_atom_name_mapping[ver] =\
                                retrieveAtomNameMappingFromInternal(self.__reg.cR, self.__reg.cacheDirPath,
                                                                    revision_history, f'{internal_cif_file_prefix}{ver}',
                                                                    self.__reg.representative_model_id,
                                                                    self.__reg.representative_alt_id,
                                                                    self.__reg.csStat)

                    if os.path.exists(internal_cif_file0):
                        major = max(revision_history)
                        minor = revision_history[major] - 1
                        revision_history[major] = minor
                        self.__reg.internal_atom_name_mapping[0] =\
                            retrieveAtomNameMappingFromInternal(self.__reg.cR, self.__reg.cacheDirPath,
                                                                revision_history, internal_cif_file0,
                                                                self.__reg.representative_model_id,
                                                                self.__reg.representative_alt_id,
                                                                self.__reg.csStat)

                        try:

                            import subprocess  # pylint: disable=import-outside-toplevel

                            database_2 = self.__reg.cR.getDictListWithFilter('database_2',
                                                                             [{'name': 'database_code', 'type': 'str'}],
                                                                             [{'name': 'database_id', 'type': 'str',
                                                                               'value': 'BMRB'}])
                            if len(database_2) > 0:

                                bmrb_id = database_2[0]['database_code']
                                if bmrb_id is not None and bmrb_id.isdigit() and not os.path.exists(other_internal_cif_file):
                                    for dir_name in os.listdir(self.__reg.cR.getDirPath()):
                                        if dir_name.startswith('bmr') and dir_name != f'bmr{bmrb_id}':
                                            dir_path = os.path.join(self.__reg.cR.getDirPath(), dir_name)
                                            if os.path.isdir(dir_path):
                                                ret_code = -1
                                                other_intnl_upload_dir = os.path.join(dir_path, 'work')
                                                if os.path.isdir(other_intnl_upload_dir):
                                                    for file_name in os.listdir(other_intnl_upload_dir):
                                                        file_path = os.path.join(other_intnl_upload_dir, file_name)
                                                        if not os.path.isfile(file_path):
                                                            continue
                                                        if file_name.endswith('.pdb'):
                                                            if not has_coordinates(file_path):
                                                                continue
                                                            com = [maxitpath, '-input', file_path,
                                                                   '-output', other_internal_cif_file,
                                                                   '-o', '1']
                                                            result = subprocess.run(com, check=False)
                                                            ret_code = result.returncode
                                                            if ret_code == 0:
                                                                break
                                                    if ret_code != 0:
                                                        cif_path_list = []
                                                        for file_name in os.listdir(other_intnl_upload_dir):
                                                            file_path = os.path.join(other_intnl_upload_dir, file_name)
                                                            if not os.path.isfile(file_path):
                                                                continue
                                                            if file_name.endswith('.cif'):
                                                                if not has_cif_coordinates(file_path):
                                                                    continue
                                                                cif_path_list.append((file_path, os.path.getsize(file_path)))
                                                        if len(cif_path_list) > 0:
                                                            cif_path_list = sorted(cif_path_list, key=lambda item: item[1])
                                                            file_path = cif_path_list[0][0]
                                                            shutil.copyfile(file_path, other_internal_cif_file)
                                                            ret_code = 0
                                                            break
                                                if ret_code != 0:
                                                    other_intnl_upload_dir = os.path.join(dir_path, 'work', 'upload')
                                                    if os.path.isdir(other_intnl_upload_dir):
                                                        for file_name in os.listdir(other_intnl_upload_dir):
                                                            file_path = os.path.join(other_intnl_upload_dir, file_name)
                                                            if not os.path.isfile(file_path):
                                                                continue
                                                            if file_name.endswith('.pdb'):
                                                                if not has_coordinates(file_path):
                                                                    continue
                                                                com = [maxitpath, '-input', file_path,
                                                                       '-output', other_internal_cif_file,
                                                                       '-o', '1']
                                                                result = subprocess.run(com, check=False)
                                                                ret_code = result.returncode
                                                                if ret_code == 0:
                                                                    break
                                                            if ret_code != 0:
                                                                cif_path_list = []
                                                                for file_name in os.listdir(other_intnl_upload_dir):
                                                                    file_path = os.path.join(other_intnl_upload_dir, file_name)
                                                                    if not os.path.isfile(file_path):
                                                                        continue
                                                                    if file_name.endswith('.cif'):
                                                                        if not has_cif_coordinates(file_path):
                                                                            continue
                                                                        cif_path_list.append((file_path,
                                                                                              os.path.getsize(file_path)))
                                                                if len(cif_path_list) > 0:
                                                                    cif_path_list = sorted(cif_path_list, key=lambda item: item[1])
                                                                    file_path = cif_path_list[0][0]
                                                                    shutil.copyfile(file_path, other_internal_cif_file)
                                                                    ret_code = 0
                                                                    break

                            if os.path.exists(internal_cif_file0):
                                major = max(revision_history)
                                minor = revision_history[major] - 2
                                revision_history[major] = minor
                                self.__reg.internal_atom_name_mapping[-1] =\
                                    retrieveAtomNameMappingFromInternal(self.__reg.cR, self.__reg.cacheDirPath,
                                                                        revision_history, other_internal_cif_file,
                                                                        self.__reg.representative_model_id,
                                                                        self.__reg.representative_alt_id,
                                                                        self.__reg.csStat)

                        except Exception:  # pylint: disable=broad-exception-caught
                            pass

            if len(self.__reg.internal_atom_name_mapping) == 0:
                cif_file_name = os.path.basename(self.__reg.cifPath)

                if WORK_MODEL_FILE_NAME_PAT.match(cif_file_name):

                    try:

                        import subprocess  # pylint: disable=import-outside-toplevel

                        cur_dir_path = self.__reg.cR.getDirPath()
                        dep_id = WORK_MODEL_FILE_NAME_PAT.search(cif_file_name).groups()[0]
                        internal_cif_file = os.path.join(cur_dir_path, self.__reg.cacheDirPath, f'{dep_id}_model-upload_P1.cif.V1')

                        if not os.path.exists(internal_cif_file):
                            ret_code = -1
                            for file_name in os.listdir(cur_dir_path):
                                if file_name == f'{dep_id}_model-upload_P1.pdb.V1':
                                    file_path = os.path.join(cur_dir_path, file_name)
                                    com = [maxitpath, '-input', file_path, '-output', internal_cif_file, '-o', '1']
                                    result = subprocess.run(com, check=False)
                                    ret_code = result.returncode
                                    if ret_code == 0:
                                        break
                            if ret_code != 0:
                                for file_name in os.listdir(cur_dir_path):
                                    if file_name == f'{dep_id}_model_P1.pdb.V1':
                                        file_path = os.path.join(cur_dir_path, file_name)
                                        com = [maxitpath, '-input', file_path, '-output', internal_cif_file, '-o', '1']
                                        result = subprocess.run(com, check=False)
                                        ret_code = result.returncode
                                        if ret_code == 0:
                                            break
                            if ret_code != 0:
                                for file_name in os.listdir(cur_dir_path):
                                    if file_name == f'{dep_id}_model-upload_P1.cif.V1':
                                        file_path = os.path.join(cur_dir_path, file_name)
                                        com = [maxitpath, '-input', file_path, '-output', internal_cif_file, '-o', '8']
                                        result = subprocess.run(com, check=False)
                                        ret_code = result.returncode
                                        if ret_code == 0:
                                            break
                            if ret_code != 0:
                                for file_name in os.listdir(cur_dir_path):
                                    if file_name == f'{dep_id}_model_P1.cif.V1':
                                        file_path = os.path.join(cur_dir_path, file_name)
                                        com = [maxitpath, '-input', file_path, '-output', internal_cif_file, '-o', '8']
                                        result = subprocess.run(com, check=False)
                                        ret_code = result.returncode
                                        if ret_code == 0:
                                            break

                        if os.path.exists(internal_cif_file):
                            self.__reg.internal_atom_name_mapping[1] =\
                                retrieveAtomNameMappingFromInternal(self.__reg.cR, self.__reg.cacheDirPath, {0: 0},
                                                                    internal_cif_file,
                                                                    self.__reg.representative_model_id,
                                                                    self.__reg.representative_alt_id,
                                                                    self.__reg.csStat)

                    except Exception:  # pylint: disable=broad-exception-caught
                        pass

            # DAOTHER-8580: convert working model file if pdbx_poly_seq_scheme category is missing
            # @see: wwpdb.utils.wf.plugins.FormatUtils.pdb2pdbxDepositOp

            if self.__reg.remediation_mode and self.__reg.recvd_nmr_data\
               and not self.__reg.cR.hasCategory('pdbx_poly_seq_scheme') and not self.__reg.cifPath.endswith('~'):

                srcCifPath = self.__reg.cifPath
                dstCifPath = self.getNextPath(self.__reg.cifPath)

                done = False

                if maxitpath is not None:

                    try:

                        import subprocess  # pylint: disable=import-outside-toplevel

                        com = [maxitpath, '-input', srcCifPath, '-output', dstCifPath, '-o', '8']
                        result = subprocess.run(com, check=False)

                        done = result.returncode == 0

                    except ImportError:
                        pass

                if not done:

                    try:
                        from wwpdb.utils.config.ConfigInfo import ConfigInfo  # pylint: disable=import-outside-toplevel
                        from wwpdb.utils.dp.RcsbDpUtility import RcsbDpUtility  # pylint: disable=import-outside-toplevel
                    except ImportError:
                        return False

                    try:

                        if self.__reg.cacheDirPath is None:
                            self.__reg.cacheDirPath = os.path.join(self.__reg.dirPath, SUB_DIR_NAME_FOR_CACHE)

                        if not os.path.isdir(self.__reg.cacheDirPath):
                            os.makedirs(self.__reg.cacheDirPath)

                        tmp_path = os.path.join(self.__reg.cacheDirPath, 'cif2cif')
                        if not os.path.isdir(tmp_path):
                            os.makedirs(tmp_path)

                        cI = ConfigInfo()
                        siteId = cI.get('SITE_PREFIX')
                        rdU = RcsbDpUtility(tmpPath=tmp_path, siteId=siteId, verbose=self.__reg.verbose, log=self.__reg.log)
                        rdU.imp(srcCifPath)
                        rdU.op('annot-cif2cif-dep')
                        rdU.exp(dstCifPath)
                        rdU.cleanup()
                        os.rmdir(tmp_path)

                    except Exception as e:  # pylint: disable=broad-exception-caught

                        self.__reg.report.error.appendDescription('internal_error',
                                                                  f"+{self.__class_name__}.__parseCoordinate() "
                                                                  "++ Error  - " + str(e))

                        return False

                self.__reg.inputParamDict[MODEL_FILE_PATH_KEY] = dstCifPath
                self.__reg.cifPath = None
                self.__reg.cifChecked = False

                return self.__parseCoordinate()

            self.__reg.cifChecked = True

            if self.__reg.caC is None:
                self.__retrieveCoordAssemblyChecker__()

            return True

        except Exception as e:  # pylint: disable=broad-exception-caught

            self.__reg.report.error.appendDescription('internal_error',
                                                      f"+{self.__class_name__}.__parseCoordinate() ++ Error  - " + str(e))

            return False

    def __parseCoordFilePath(self) -> bool:
        """ Parse effective coordinate file path.
        """

        if self.__reg.cifPath is not None:
            return True

        self.__cifHashCode = None

        if MODEL_FILE_PATH_KEY in self.__reg.inputParamDict:

            fPath = self.__reg.inputParamDict[MODEL_FILE_PATH_KEY]

            if fPath.endswith('.gz'):

                _fPath = os.path.splitext(fPath)[0]

                if not os.path.exists(_fPath):

                    try:

                        uncompress_gzip_file(fPath, _fPath)

                    except Exception as e:  # pylint: disable=broad-exception-caught

                        self.__reg.report.error.appendDescription('internal_error',
                                                                  f"+{self.__class_name__}.__parseCoordFilePath() "
                                                                  "++ Error  - " + str(e))

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__parseCoordFilePath() ++ Error  - {str(e)}\n")

                        return False

                fPath = _fPath

            try:

                self.__reg.spareDirPath = os.path.dirname(fPath)
                if self.__reg.dirPath is None:
                    self.__reg.dirPath = self.__reg.spareDirPath

                if self.__reg.cacheDirPath is None:
                    self.__reg.cacheDirPath = os.path.join(self.__reg.dirPath, SUB_DIR_NAME_FOR_CACHE)

                # share the cache directory with CifReader
                self.__reg.cR.cacheDirPath = self.__reg.cacheDirPath

                self.__reg.cifPath = fPath

                if self.__reg.cR.parse(fPath):
                    return True

                # try deposit storage if possible
                if ALT_MODEL_FILE_PATH_KEY in self.__reg.inputParamDict:

                    fPath = self.__reg.inputParamDict[ALT_MODEL_FILE_PATH_KEY]

                    self.__reg.cifPath = fPath

                    if self.__reg.cR.parse(fPath):
                        return True

            except IOError:
                pass

            finally:
                self.__reg.caC = None
                self.__reg.cpC = copy.deepcopy(DEFAULT_COORD_PROPERTIES)
                self.__reg.exptl_method = None
                self.__reg.symmetric = None
                self.__reg.is_cyclic_polymer.clear()
                self.__reg.nmr_struct_conf.clear()
                self.__reg.ent_asym_id_with_exptl_data.clear()
                self.__reg.label_asym_id_with_exptl_data.clear()
                self.__reg.auth_asym_ids_with_chem_exch.clear()
                self.__reg.auth_seq_ids_with_chem_exch.clear()
                self.__reg.chain_id_map_for_remediation.clear()
                self.__reg.seq_id_map_for_remediation.clear()

                self.__cifHashCode = self.__reg.cR.getHashCode()

                self.__reg.coord_atom_site_tags = self.__reg.cR.getAttributeList('atom_site')

        return False

    def __detectCoordContentSubType(self) -> bool:
        """ Detect content subtype of coordinate file.
        """

        src_id = self.__reg.report.getInputSourceIdOfCoord()

        if src_id < 0:
            return False

        cif_input_source = self.__reg.report.input_sources[src_id]
        cif_input_source_dic = cif_input_source.get()

        if has_key_value(cif_input_source_dic, 'content_subtype'):
            return False

        file_type = cif_input_source_dic['file_type']

        lp_counts = {t: 0 for t in CIF_CONTENT_SUBTYPES}

        for content_subtype in CIF_CONTENT_SUBTYPES:

            lp_category = LP_CATEGORIES[file_type][content_subtype]

            if self.__reg.cR.hasCategory(lp_category):
                lp_counts[content_subtype] = 1

            elif content_subtype != 'branched':

                if content_subtype != 'non_poly':

                    if content_subtype == 'poly_seq'\
                       and self.__reg.cR.hasCategory(LP_CATEGORIES[file_type][f'{content_subtype}_alias']):
                        lp_counts[content_subtype] = 1

                elif self.__reg.cR.hasCategory(LP_CATEGORIES[file_type][f'{content_subtype}_alias']):
                    lp_counts[content_subtype] = 1

        content_subtypes = {k: lp_counts[k] for k in lp_counts if lp_counts[k] > 0}

        cif_input_source.setItemValue('content_subtype', content_subtypes)

        return True

    def __extractCoordPolymerSequence(self) -> bool:
        """ Extract reference polymer sequence of coordinate file.
        """

        src_id = self.__reg.report.getInputSourceIdOfCoord()

        if src_id < 0:

            if (not self.__reg.combined_mode and self.__reg.permit_missing_legacy_dist_restraint)\
               or self.__reg.conversion_server:  # no exception

                if len(self.__reg.suspended_errors_for_lazy_eval) > 0:
                    for msg in self.__reg.suspended_errors_for_lazy_eval:
                        for k, v in msg.items():
                            self.__reg.report.error.appendDescription(k, v)

                            if k == 'missing_mandatory_content'\
                               and 'Deposition of assigned chemical shifts is mandatory' in v['description']\
                               and self.__reg.remediation_mode:
                                dir_path = os.path.dirname(self.__reg.dstPath)\
                                    if self.__reg.dirPath is None else self.__reg.dirPath

                                touch_file = os.path.join(dir_path, '.entry_without_cs')
                                if not os.path.exists(touch_file):
                                    with open(touch_file, 'w', encoding='utf-8') as ofh:
                                        ofh.write('')

                    self.__reg.suspended_errors_for_lazy_eval.clear()

                if len(self.__reg.suspended_warnings_for_lazy_eval) > 0:
                    for msg in self.__reg.suspended_warnings_for_lazy_eval:
                        for k, v in msg.items():
                            self.__reg.report.warning.appendDescription(k, v)
                    self.__reg.suspended_warnings_for_lazy_eval.clear()

            return False

        cif_input_source = self.__reg.report.input_sources[src_id]
        cif_input_source_dic = cif_input_source.get()

        file_name = cif_input_source_dic['file_name']
        file_type = cif_input_source_dic['file_type']

        if cif_input_source_dic['content_subtype'] is None:
            return False

        content_subtype = 'poly_seq'

        if content_subtype not in cif_input_source_dic['content_subtype']:
            return False

        if has_key_value(cif_input_source_dic, 'polymer_sequence'):
            return False

        alias = False
        lp_category = LP_CATEGORIES[file_type][content_subtype]
        key_items = self.__reg.key_items[file_type][content_subtype]

        if self.__reg.cR.hasItem(lp_category, 'pdb_mon_id'):
            _key_items = copy.copy(key_items)
            _key_items.append({'name': 'pdb_mon_id', 'type': 'str', 'alt_name': 'auth_comp_id', 'default-from': 'mon_id'})
            key_items = _key_items

        if self.__reg.cR.hasItem(lp_category, 'auth_mon_id'):
            _key_items = copy.copy(key_items)
            _key_items.append({'name': 'auth_mon_id', 'type': 'str', 'alt_name': 'alt_comp_id', 'default-from': 'mon_id'})
            key_items = _key_items

        if not self.__reg.cR.hasCategory(lp_category):
            alias = True
            lp_category = LP_CATEGORIES[file_type][f'{content_subtype}_alias']
            key_items = self.__reg.key_items[file_type][f'{content_subtype}_alias']

        try:

            poly_seq = cluster = poly_seq_cache_path = cluster_cache_path = None

            if self.__cifHashCode is not None:
                poly_seq_cache_path = os.path.join(self.__reg.cacheDirPath, f"{self.__cifHashCode}_poly_seq_full.pkl")
                poly_seq = load_from_pickle(poly_seq_cache_path)
                cluster_cache_path = os.path.join(self.__reg.cacheDirPath, f"{self.__cifHashCode}_cluster.pkl")
                cluster = load_from_pickle(cluster_cache_path)

            if poly_seq is None:

                try:
                    poly_seq, cluster = self.__reg.cR.getPolymerSequence(lp_category, key_items,
                                                                         withStructConf=True, withRmsd=True, alias=alias,
                                                                         totalModels=self.__reg.total_models,
                                                                         effModelIds=self.__reg.eff_model_ids,
                                                                         repAltId=self.__reg.representative_alt_id)
                except KeyError:  # pdbx_PDB_ins_code throws KeyError
                    if content_subtype + ('_ins_alias' if alias else '_ins') in self.__reg.key_items[file_type]:
                        key_items = self.__reg.key_items[file_type][content_subtype + ('_ins_alias' if alias else '_ins')]
                        poly_seq, cluster = self.__reg.cR.getPolymerSequence(lp_category, key_items,
                                                                             withStructConf=True, withRmsd=True, alias=alias,
                                                                             totalModels=self.__reg.total_models,
                                                                             effModelIds=self.__reg.eff_model_ids,
                                                                             repAltId=self.__reg.representative_alt_id)
                    else:
                        poly_seq = []

                if len(poly_seq) == 0:
                    return False

                content_subtype = 'branched'

                lp_category = LP_CATEGORIES[file_type][content_subtype]

                if self.__reg.cR.hasCategory(lp_category):

                    key_items = self.__reg.key_items[file_type][content_subtype]

                    if self.__reg.cR.hasItem(lp_category, 'pdb_mon_id'):
                        _key_items = copy.copy(key_items)
                        _key_items.append({'name': 'pdb_mon_id', 'type': 'str', 'alt_name': 'auth_comp_id',
                                           'default-from': 'mon_id'})
                        key_items = _key_items

                    if self.__reg.cR.hasItem(lp_category, 'auth_mon_id'):
                        _key_items = copy.copy(key_items)
                        _key_items.append({'name': 'auth_mon_id', 'type': 'str', 'alt_name': 'alt_comp_id',
                                           'default-from': 'mon_id'})
                        key_items = _key_items

                    try:

                        branched_seq, _ = self.__reg.cR.getPolymerSequence(lp_category, key_items,
                                                                           withStructConf=False, withRmsd=False, alias=False,
                                                                           totalModels=self.__reg.total_models,
                                                                           effModelIds=self.__reg.eff_model_ids,
                                                                           repAltId=self.__reg.representative_alt_id)
                        if len(branched_seq) > 0:
                            poly_seq.extend(branched_seq)

                    except Exception:  # pylint: disable=broad-exception-caught
                        pass

                content_subtype = 'non_poly'

                lp_category = LP_CATEGORIES[file_type][content_subtype]

                if self.__reg.cR.hasCategory(lp_category):

                    key_items = self.__reg.key_items[file_type][content_subtype]

                    if self.__reg.cR.hasItem(lp_category, 'pdb_mon_id'):
                        _key_items = copy.copy(key_items)
                        _key_items.append({'name': 'pdb_mon_id', 'type': 'str', 'alt_name': 'auth_comp_id',
                                           'default-from': 'mon_id'})
                        key_items = _key_items

                    if self.__reg.cR.hasItem(lp_category, 'auth_mon_id'):
                        _key_items = copy.copy(key_items)
                        _key_items.append({'name': 'auth_mon_id', 'type': 'str', 'alt_name': 'alt_comp_id',
                                           'default-from': 'mon_id'})
                        key_items = _key_items

                    try:

                        non_poly, _ = self.__reg.cR.getPolymerSequence(lp_category, key_items,
                                                                       withStructConf=False, withRmsd=False, alias=False,
                                                                       totalModels=self.__reg.total_models,
                                                                       effModelIds=self.__reg.eff_model_ids,
                                                                       repAltId=self.__reg.representative_alt_id)

                        if len(non_poly) > 0:
                            poly_seq.extend(non_poly)

                    except Exception:  # pylint: disable=broad-exception-caught
                        pass

                if len(poly_seq) > 0 and poly_seq_cache_path is not None:
                    write_as_pickle(poly_seq, poly_seq_cache_path)

                if cluster is not None and len(cluster) > 0 and cluster_cache_path is not None:
                    write_as_pickle(cluster, cluster_cache_path)

            cif_input_source.setItemValue('polymer_sequence', poly_seq)

            not_superimposed_ensemble, exactly_overlaid_ensemble, exactly_overlaid_models = {}, {}, {}

            if not self.__reg.combined_mode and self.__reg.permit_missing_legacy_dist_restraint:  # no exception

                if len(self.__reg.suspended_errors_for_lazy_eval) > 0:
                    for msg in self.__reg.suspended_errors_for_lazy_eval:
                        for k, v in msg.items():
                            self.__reg.report.error.appendDescription(k, v)

                            if k == 'missing_mandatory_content'\
                               and 'Deposition of assigned chemical shifts is mandatory' in v['description']\
                               and self.__reg.remediation_mode:
                                dir_path = os.path.dirname(self.__reg.dstPath)\
                                    if self.__reg.dirPath is None else self.__reg.dirPath

                                touch_file = os.path.join(dir_path, '.entry_without_cs')
                                if not os.path.exists(touch_file):
                                    with open(touch_file, 'w', encoding='utf-8') as ofh:
                                        ofh.write('')

                    self.__reg.suspended_errors_for_lazy_eval.clear()

                if len(self.__reg.suspended_warnings_for_lazy_eval) > 0:
                    for msg in self.__reg.suspended_warnings_for_lazy_eval:
                        for k, v in msg.items():
                            self.__reg.report.warning.appendDescription(k, v)
                    self.__reg.suspended_warnings_for_lazy_eval.clear()

            domain_ids = set()

            for ps in poly_seq:

                if 'type' in ps:

                    poly_type = ps['type']

                    if 'polypeptide' in poly_type:
                        rmsd_label = 'ca_rmsd'

                        if not self.__reg.combined_mode:

                            if len(self.__reg.suspended_errors_for_lazy_eval) > 0:
                                for msg in self.__reg.suspended_errors_for_lazy_eval:
                                    for k, v in msg.items():
                                        self.__reg.report.error.appendDescription(k, v)
                                self.__reg.suspended_errors_for_lazy_eval.clear()

                            if len(self.__reg.suspended_warnings_for_lazy_eval) > 0:
                                for msg in self.__reg.suspended_warnings_for_lazy_eval:
                                    for k, v in msg.items():
                                        self.__reg.report.warning.appendDescription(k, v)
                                self.__reg.suspended_warnings_for_lazy_eval.clear()

                    elif 'ribonucleotide' in poly_type:
                        rmsd_label = 'p_rmsd'
                    else:
                        continue

                    chain_id = ps['chain_id']

                    if rmsd_label in ps and 'well_defined_region' in ps:
                        rmsd = ps[rmsd_label]
                        region = ps['well_defined_region']

                        for r in region:
                            domain_ids.add(r['domain_id'])

                        for r in rmsd:
                            model_id = r['model_id']

                            if 'raw_rmsd_in_well_defined_region' in r and 'rmsd_in_well_defined_region' in r:

                                if r['raw_rmsd_in_well_defined_region'] - r['rmsd_in_well_defined_region'] >\
                                        self.__reg.rmsd_not_superimposed:
                                    rmsd_item = {'model_id': model_id,
                                                 'raw_rmsd': r['raw_rmsd_in_well_defined_region'],
                                                 'rmsd': r['rmsd_in_well_defined_region']}
                                    domain_id = r['domain_id']
                                    domain = next((r for r in region if r['domain_id'] == domain_id), None)
                                    if domain is not None:
                                        rmsd_item['monomers'] = domain['number_of_monomers']
                                        rmsd_item['gaps'] = domain['number_of_gaps']
                                        rmsd_item['core'] = domain['percent_of_core']
                                        rmsd_item['range'] = domain['range_of_seq_id']
                                        if chain_id not in not_superimposed_ensemble:
                                            not_superimposed_ensemble[chain_id] = []
                                        not_superimposed_ensemble[chain_id].append(rmsd_item)

                                if r['rmsd_in_well_defined_region'] < self.__reg.rmsd_overlaid_exactly:
                                    domain_id = r['domain_id']
                                    domain = next((r for r in region if r['domain_id'] == domain_id), None)
                                    if domain is not None and domain['mean_rmsd'] < self.__reg.rmsd_overlaid_exactly:
                                        region_item = {'monomers': domain['number_of_monomers'],
                                                       'gaps': domain['number_of_gaps'],
                                                       'core': domain['percent_of_core'],
                                                       'mean_rmsd': domain['mean_rmsd'],
                                                       'range': domain['range_of_seq_id']}
                                        exactly_overlaid_ensemble[chain_id] = region_item

                                elif 'exactly_overlaid_model' in r:
                                    domain_id = r['domain_id']
                                    domain = next((r for r in region if r['domain_id'] == domain_id), None)
                                    if domain is not None:
                                        for m in r['exactly_overlaid_model']:
                                            rmsd_item = {'model_id_1': m['ref_model_id'],
                                                         'model_id_2': m['test_model_id'],
                                                         'rmsd': m['rmsd_in_well_defined_region']}
                                            rmsd_item['monomers'] = domain['number_of_monomers']
                                            rmsd_item['gaps'] = domain['number_of_gaps']
                                            rmsd_item['core'] = domain['percent_of_core']
                                            rmsd_item['range'] = domain['range_of_seq_id']
                                            if chain_id not in exactly_overlaid_models:
                                                exactly_overlaid_models[chain_id] = []
                                            exactly_overlaid_models[chain_id].append(rmsd_item)

            if len(domain_ids) > 0:
                single_seq_id_pat = re.compile(r'(-?\d+)')
                range_seq_id_pat = re.compile(r'(-?\d+)-(-?\d+)')
                domain_ids = sorted(list(domain_ids))

                well_defined_region = []

                for domain_id in domain_ids:
                    _medoid_model_id = -1
                    _mean_rmsd = _medoid_rmsd = -1.0
                    _monomers = _total_seq = 0
                    _range_of_seq_id = []

                    for ps in poly_seq:

                        if 'type' in ps:

                            poly_type = ps['type']

                            if 'polypeptide' in poly_type:
                                rmsd_label = 'ca_rmsd'

                                if not self.__reg.combined_mode:

                                    if len(self.__reg.suspended_errors_for_lazy_eval) > 0:
                                        for msg in self.__reg.suspended_errors_for_lazy_eval:
                                            for k, v in msg.items():
                                                self.__reg.report.error.appendDescription(k, v)
                                        self.__reg.suspended_errors_for_lazy_eval.clear()

                                    if len(self.__reg.suspended_warnings_for_lazy_eval) > 0:
                                        for msg in self.__reg.suspended_warnings_for_lazy_eval:
                                            for k, v in msg.items():
                                                self.__reg.report.warning.appendDescription(k, v)
                                        self.__reg.suspended_warnings_for_lazy_eval.clear()

                            elif 'ribonucleotide' in poly_type:
                                rmsd_label = 'p_rmsd'

                            else:
                                continue

                            auth_chain_id = ps['auth_chain_id']
                            _total_seq += len(ps['seq_id'])

                            if rmsd_label in ps and 'well_defined_region' in ps:

                                for region in ps['well_defined_region']:

                                    if region['domain_id'] != domain_id:
                                        continue

                                    if _medoid_model_id < 0:
                                        _medoid_model_id = region['medoid_model_id']
                                        _mean_rmsd = region['mean_rmsd']
                                        _medoid_rmsd = region['medoid_rmsd']

                                    _monomers += region['number_of_monomers']

                                    ranges = re.sub(r'\]', '', re.sub(r'\[', '', region['range_of_seq_id'])).split(',')
                                    _ranges = []
                                    for r in ranges:
                                        if range_seq_id_pat.match(r):
                                            g = range_seq_id_pat.search(r).groups()
                                            _ranges.append(f'{auth_chain_id}:{g[0]}-{auth_chain_id}:{g[1]}')
                                        elif single_seq_id_pat.match(r):
                                            g = single_seq_id_pat.search(r).groups()
                                            _ranges.append(f'{auth_chain_id}:{g[0]}')

                                    _range_of_seq_id.append(f"{','.join(_ranges)}")

                    well_defined_region.append({'domain_id': domain_id,
                                                'medoid_model_id': _medoid_model_id,
                                                'number_of_monomers': _monomers,
                                                'percent_of_core': round(float(_monomers) / _total_seq * 100.0, 1),
                                                'mean_rmsd': _mean_rmsd,
                                                'medoid_rmsd': _medoid_rmsd,
                                                'range_of_seq_id': f"{', '.join(_range_of_seq_id)}"})

                if len(well_defined_region) > 0:
                    self.__ensemble_composition['well_defined_region'] = well_defined_region
                    if cluster is not None:
                        self.__ensemble_composition['cluster_analysis'] = cluster

            cif_input_source.setItemValue('ensemble_composition', self.__ensemble_composition)

            if len(not_superimposed_ensemble) > 0:

                for chain_id, rmsd in not_superimposed_ensemble.items():

                    conformer_id = 1

                    nmr_representative = self.__reg.cR.getDictList('pdbx_nmr_representative')

                    if len(nmr_representative) > 0:

                        try:
                            conformer_id = int(nmr_representative[0]['conformer_id'])
                        except ValueError:
                            conformer_id = 1

                    r = next((r for r in rmsd if r['model_id'] == conformer_id), rmsd[0])

                    warn = f"The coordinates (chain_id {chain_id}) are not superimposed. "\
                        f"RMSD ({r['raw_rmsd']}Å) for a well-defined region "\
                        f"(Sequence ranges {r['range']}) is greater than predicted value ({r['rmsd']}Å). "\
                        "Please superimpose the coordinates and re-upload the model file."

                    self.__reg.report.warning.appendDescription('not_superimposed_model',
                                                                {'file_name': file_name, 'category': 'atom_site',
                                                                 'description': warn})

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__extractCoordPolymerSequence() ++ Warning  - {warn}\n")

            elif len(exactly_overlaid_ensemble) > 0:

                for chain_id, r in exactly_overlaid_ensemble.items():

                    warn = f"The coordinates (chain_id {chain_id}) are overlaid exactly. "\
                        "Please check there has not been an error during the creation of your model file. "\
                        "You are receiving this message because the mean RMSD for a well-defined region "\
                        f"(Sequence ranges {r['range']}) is {r['mean_rmsd']}Å. "\
                        "We require you to deposit an appropriate ensemble of coordinate models."

                    self.__reg.report.warning.appendDescription('exactly_overlaid_model',
                                                                {'file_name': file_name, 'category': 'atom_site',
                                                                 'description': warn})

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__extractCoordPolymerSequence() ++ Warning  - {warn}\n")

            elif len(exactly_overlaid_models) > 0:

                for chain_id, rs in exactly_overlaid_models.items():

                    for r in rs:

                        warn = f"Two models in the coordinate file (chain_id {chain_id}) are overlaid exactly. "\
                            "Please check there has not been an error during the creation of your model file. "\
                            "You are receiving this message because the RMSD for a well-defined region "\
                            f"(Sequence ranges {r['range']}) between model {r['model_id_1']!r} and model {r['model_id_2']!r} "\
                            f"is {r['rmsd']}Å. "\
                            "We require you to deposit an appropriate ensemble of coordinate models."

                        self.__reg.report.warning.appendDescription('exactly_overlaid_model',
                                                                    {'file_name': file_name, 'category': 'atom_site',
                                                                     'description': warn})

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__extractCoordPolymerSequence() ++ Warning  - {warn}\n")

            return True

        except KeyError as e:

            self.__reg.report.error.appendDescription('sequence_mismatch',
                                                      {'file_name': file_name, 'category': lp_category,
                                                       'description': str(e).strip("'")})

            if self.__reg.verbose:
                self.__reg.log.write(f"+{self.__class_name__}.__extractCoordPolymerSequence() ++ KeyError  - {str(e)}\n")

        except LookupError as e:

            self.__reg.report.error.appendDescription('missing_mandatory_item',
                                                      {'file_name': file_name, 'category': lp_category,
                                                       'description': str(e).strip("'")})

            self.__reg.log.write(f"+{self.__class_name__}.__extractCoordPolymerSequence() "
                                 f"++ LookupError  - {file_name} {lp_category} {str(e)}\n")

        except ValueError as e:

            self.__reg.report.error.appendDescription('invalid_data',
                                                      {'file_name': file_name, 'category': lp_category,
                                                       'description': str(e).strip("'")})

            if self.__reg.verbose:
                self.__reg.log.write(f"+{self.__class_name__}.__extractCoordPolymerSequence() ++ ValueError  - {str(e)}\n")

        except Exception as e:  # pylint: disable=broad-exception-caught

            self.__reg.report.error.appendDescription('internal_error',
                                                      f"+{self.__class_name__}.__extractCoordPolymerSequence() "
                                                      "++ Error  - " + str(e))

            if self.__reg.verbose:
                self.__reg.log.write(f"+{self.__class_name__}.__extractCoordPolymerSequence() ++ Error  - {str(e)}\n")

        return False

    def __extractCoordPolymerSequenceInLoop(self) -> bool:
        """ Extract polymer sequence in interesting loops of coordinate file.
        """

        src_id = self.__reg.report.getInputSourceIdOfCoord()

        if src_id < 0:
            return False

        __errors = self.__reg.report.getTotalErrors()

        cif_input_source = self.__reg.report.input_sources[src_id]
        cif_input_source_dic = cif_input_source.get()

        file_name = cif_input_source_dic['file_name']
        file_type = cif_input_source_dic['file_type']

        poly_seq_list_set = {}

        for content_subtype in CIF_CONTENT_SUBTYPES:

            if content_subtype in ('entry_info', 'poly_seq', 'branched', 'non_poly')\
               or (not has_key_value(cif_input_source_dic['content_subtype'], content_subtype)):
                continue

            poly_seq_list_set[content_subtype] = []

            alias = False
            lp_category = LP_CATEGORIES[file_type][content_subtype]
            key_items = self.__reg.key_items[file_type][content_subtype]

            if not self.__reg.cR.hasCategory(lp_category):
                alias = True
                lp_category = LP_CATEGORIES[file_type][f'{content_subtype}_alias']
                key_items = self.__reg.key_items[file_type][f'{content_subtype}_alias']

            elif content_subtype == 'coordinate' and 'pdbx_PDB_model_num' not in self.__reg.coord_atom_site_tags:
                alias = True
                key_items = self.__reg.key_items[file_type][f'{content_subtype}_alias']

            has_poly_seq = False

            list_id = 1

            try:

                poly_seq = poly_seq_cache_path = None

                if self.__cifHashCode is not None:
                    poly_seq_cache_path = os.path.join(self.__reg.cacheDirPath, f"{self.__cifHashCode}_poly_seq.pkl")
                    poly_seq = load_from_pickle(poly_seq_cache_path)

                if poly_seq is None:

                    try:
                        poly_seq, _ = self.__reg.cR.getPolymerSequence(lp_category, key_items)
                    except KeyError:  # pdbx_PDB_ins_code throws KeyError
                        if content_subtype + ('_ins_alias' if alias else '_ins') in self.__reg.key_items[file_type]:
                            key_items = self.__reg.key_items[file_type][content_subtype + ('_ins_alias' if alias else '_ins')]
                            poly_seq, _ = self.__reg.cR.getPolymerSequence(lp_category, key_items)
                        else:
                            poly_seq = []

                    if len(poly_seq) > 0 and poly_seq_cache_path is not None:
                        write_as_pickle(poly_seq, poly_seq_cache_path)

                if len(poly_seq) > 0:
                    poly_seq_list_set[content_subtype].append({'list_id': list_id, 'polymer_sequence': poly_seq})

                    has_poly_seq = True

            except KeyError as e:

                self.__reg.report.error.appendDescription('sequence_mismatch',
                                                          {'file_name': file_name, 'category': lp_category,
                                                           'description': str(e).strip("'")})

                if self.__reg.verbose:
                    self.__reg.log.write(f"+{self.__class_name__}.__extractCoordPolymerSequenceInLoop() "
                                         f"++ KeyError  - {str(e)}\n")

            except LookupError as e:

                self.__reg.report.error.appendDescription('missing_mandatory_item',
                                                          {'file_name': file_name, 'category': lp_category,
                                                           'description': str(e).strip("'")})

                self.__reg.log.write(f"+{self.__class_name__}.__extractCoordPolymerSequenceInLoop() "
                                     f"++ LookupError  - {file_name} {lp_category} {str(e)}\n")

            except ValueError as e:

                if not (content_subtype == 'non_poly' and alias):
                    self.__reg.report.error.appendDescription('invalid_data',
                                                              {'file_name': file_name, 'category': lp_category,
                                                               'description': str(e).strip("'")})

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__extractCoordPolymerSequenceInLoop() "
                                             f"++ ValueError  - {str(e)}\n")

            except Exception as e:  # pylint: disable=broad-exception-caught

                self.__reg.report.error.appendDescription('internal_error',
                                                          f"+{self.__class_name__}.__extractCoordPolymerSequenceInLoop() "
                                                          "++ Error  - " + str(e))

                if self.__reg.verbose:
                    self.__reg.log.write(f"+{self.__class_name__}.__extractCoordPolymerSequenceInLoop() "
                                         f"++ Error  - {str(e)}\n")

            list_id += 1

            if not has_poly_seq:
                poly_seq_list_set.pop(content_subtype)

        if self.__reg.report.getTotalErrors() > __errors:
            return False

        if len(poly_seq_list_set) > 0:
            cif_input_source.setItemValue('polymer_sequence_in_loop', poly_seq_list_set)

        return True

    def __extractCoordCommonPolymerSequence(self) -> bool:
        """ Extract common polymer sequence of coordinate file if required.
        """

        common_poly_seq = {}

        src_id = self.__reg.report.getInputSourceIdOfCoord()

        if src_id < 0:
            return False

        cif_input_source = self.__reg.report.input_sources[src_id]
        cif_input_source_dic = cif_input_source.get()

        has_poly_seq = has_key_value(cif_input_source_dic, 'polymer_sequence')
        has_poly_seq_in_lp = has_key_value(cif_input_source_dic, 'polymer_sequence_in_loop')

        if has_poly_seq or (not has_poly_seq_in_lp):
            return False

        poly_seq_in_lp = cif_input_source_dic['polymer_sequence_in_loop']

        chain_ids = set()

        for content_subtype in poly_seq_in_lp.keys():

            for _poly_seq_in_lp in poly_seq_in_lp[content_subtype]:
                poly_seq = _poly_seq_in_lp['polymer_sequence']

                for ps in poly_seq:
                    chain_id = ps['chain_id']
                    chain_ids.add(chain_id)

                    if chain_id not in common_poly_seq:
                        common_poly_seq[chain_id] = set()

        _offset_seq_ids = {c: 0 for c in chain_ids}

        for content_subtype in poly_seq_in_lp.keys():

            for _poly_seq_in_lp in poly_seq_in_lp[content_subtype]:
                poly_seq = _poly_seq_in_lp['polymer_sequence']

                for ps in poly_seq:
                    chain_id = ps['chain_id']

                    min_seq_id = min(ps['seq_id'])
                    if min_seq_id < _offset_seq_ids[chain_id]:
                        _offset_seq_ids[chain_id] = min_seq_id

        offset_seq_ids = {k: (0 if v >= 0 else -v) for k, v in _offset_seq_ids.items()}

        for content_subtype in poly_seq_in_lp.keys():

            for _poly_seq_in_lp in poly_seq_in_lp[content_subtype]:
                poly_seq = _poly_seq_in_lp['polymer_sequence']

                for ps in poly_seq:
                    chain_id = ps['chain_id']

                    for seq_id, comp_id in zip(ps['seq_id'], ps['comp_id']):
                        common_poly_seq[chain_id].add((seq_id + offset_seq_ids[chain_id], comp_id))

        asm = []  # molecular assembly of a loop

        for chain_id in sorted(common_poly_seq.keys()):

            if len(common_poly_seq[chain_id]) > 0:
                seq_ids = sorted(set(item[0] - offset_seq_ids[chain_id] for item in common_poly_seq[chain_id]))
                comp_ids = []

                for seq_id in seq_ids:
                    _comp_ids = [item[1] for item in common_poly_seq[chain_id]
                                 if item[0] - offset_seq_ids[chain_id] == seq_id]
                    if len(_comp_ids) == 1:
                        comp_ids.append(_comp_ids[0])
                    else:
                        comp_ids.append(next(comp_id for comp_id in _comp_ids if comp_id not in EMPTY_VALUE))

                asm.append({'chain_id': chain_id, 'seq_id': seq_ids, 'comp_id': comp_ids})

        if len(asm) > 0:
            cif_input_source.setItemValue('polymer_sequence', asm)

        return True

    def __extractCoordNonStandardResidue(self) -> bool:
        """ Extract non-standard residue of coordinate file.
        """

        src_id = self.__reg.report.getInputSourceIdOfCoord()

        if src_id < 0:
            return False

        cif_input_source = self.__reg.report.input_sources[src_id]
        cif_input_source_dic = cif_input_source.get()

        has_poly_seq = has_key_value(cif_input_source_dic, 'polymer_sequence')

        if not has_poly_seq:
            return False

        poly_seq = cif_input_source_dic['polymer_sequence']

        asm = []

        for ps in poly_seq:

            has_nstd_res = False

            ent = {'chain_id': ps['chain_id'], 'seq_id': [], 'comp_id': [], 'chem_comp_name': [], 'exptl_data': []}

            for seq_id, comp_id in zip(ps['seq_id'], ps['comp_id']):

                if comp_id not in STD_MON_DICT:
                    has_nstd_res = True

                    ent['seq_id'].append(seq_id)
                    ent['comp_id'].append(comp_id)

                    is_valid, cc_name, cc_rel_status = self.__reg.dpV.getChemCompNameAndStatusOf(comp_id)

                    if is_valid:  # matches with comp_id in CCD
                        if cc_rel_status == 'REL' or cc_name is not None:
                            ent['chem_comp_name'].append(cc_name)
                        else:
                            ent['chem_comp_name'].append(f"(Not available due to CCD status code {cc_rel_status})")

                    else:
                        ent['chem_comp_name'].append(cc_name)

                    ent['exptl_data'].append({'coordinate': False})

            if has_nstd_res:
                asm.append(ent)

        if len(asm) > 0:
            cif_input_source.setItemValue('non_standard_residue', asm)

        return True

    def __appendCoordPolymerSequenceAlignment(self) -> bool:
        """ Append polymer sequence alignment between coordinate and NMR data.
        """

        src_id = self.__reg.report.getInputSourceIdOfCoord()

        if src_id < 0:
            return False

        # sequence alignment inside coordinate file

        cif_input_source = self.__reg.report.input_sources[src_id]
        cif_input_source_dic = cif_input_source.get()

        has_poly_seq = has_key_value(cif_input_source_dic, 'polymer_sequence')
        has_poly_seq_in_lp = has_key_value(cif_input_source_dic, 'polymer_sequence_in_loop')

        if not has_poly_seq:
            return False

        cif_poly_seq = cif_input_source_dic['polymer_sequence']

        if has_poly_seq_in_lp:

            poly_seq_in_lp = cif_input_source_dic['polymer_sequence_in_loop']

            for content_subtype in poly_seq_in_lp.keys():

                if content_subtype in ('non_poly', 'branched'):
                    continue

                seq_align_set = []

                for i1, ps1 in enumerate(cif_poly_seq):
                    chain_id = ps1['chain_id']

                    if i1 >= LEN_MAJOR_ASYM_ID:  # to process large assembly avoiding forced timeout
                        continue

                    for _poly_seq_in_lp in poly_seq_in_lp[content_subtype]:
                        poly_seq2 = _poly_seq_in_lp['polymer_sequence']

                        for ps2 in poly_seq2:

                            if chain_id != ps2['chain_id']:
                                continue

                            self.__reg.pA.clear()
                            self.__reg.pA.setReferenceSequence(ps1['comp_id'], f'REF{chain_id}')
                            self.__reg.pA.addTestSequence(ps2['comp_id'], chain_id)
                            self.__reg.pA.doAlign()

                            myAlign = self.__reg.pA.getAlignment(chain_id)

                            length = len(myAlign)

                            if length == 0:
                                continue

                            _matched, unmapped, conflict, offset_1, offset_2 = getScoreOfSeqAlign(myAlign)

                            if length == unmapped + conflict or _matched <= conflict + (1 if length > 1 else 0):
                                continue

                            _ps1 = ps1 if offset_1 == 0 else fillBlankCompIdWithOffset(ps1, offset_1)
                            _ps2 = ps2 if offset_2 == 0 else fillBlankCompIdWithOffset(ps2, offset_2)

                            seq_align = self.__getSeqAlignCodeWithBlankFill__(length, len(ps1['seq_id']), _matched,
                                                                              conflict, unmapped,
                                                                              offset_1, offset_2, _ps1, _ps2, myAlign)
                            seq_align['list_id'] = _poly_seq_in_lp['list_id']
                            seq_align['chain_id'] = chain_id

                            seq_align_set.append(seq_align)

                if len(seq_align_set) > 0:
                    self.__reg.report.sequence_alignment.setItemValue(f'model_poly_seq_vs_{content_subtype}', seq_align_set)

        # sequence alignment between model and NMR data

        nmr_input_source = self.__reg.report.input_sources[0]
        nmr_input_source_dic = nmr_input_source.get()

        has_nmr_poly_seq = has_key_value(nmr_input_source_dic, 'polymer_sequence')

        if not has_nmr_poly_seq:
            return False

        nmr_poly_seq = nmr_input_source_dic['polymer_sequence']

        seq_align_set = []

        for i1, ps1 in enumerate(cif_poly_seq):
            chain_id = ps1['chain_id']

            if i1 >= LEN_MAJOR_ASYM_ID / 2:  # to process large assembly avoiding forced timeout
                continue

            for i2, ps2 in enumerate(nmr_poly_seq):
                chain_id2 = ps2['chain_id']

                if i2 >= LEN_MAJOR_ASYM_ID / 2:  # to process large assembly avoiding forced timeout
                    continue

                self.__reg.pA.clear()
                self.__reg.pA.setReferenceSequence(ps1['auth_comp_id'] if 'auth_comp_id' in ps1 else ps2['comp_id'],
                                                   f'REF{chain_id}')
                self.__reg.pA.addTestSequence(ps2['comp_id'], chain_id)
                self.__reg.pA.doAlign()

                myAlign = self.__reg.pA.getAlignment(chain_id)

                length = len(myAlign)

                if length == 0:
                    continue

                _matched, unmapped, conflict, offset_1, offset_2 = getScoreOfSeqAlign(myAlign)

                if length == unmapped + conflict:
                    if len(ps1['seq_id']) == 1 and 'alt_comp_id' in ps1 and ps1['alt_comp_id'][0] in ps2['comp_id']:
                        self.__reg.pA.clear()
                        self.__reg.pA.setReferenceSequence(ps1['alt_comp_id'], f'REF{chain_id}')
                        self.__reg.pA.addTestSequence(ps2['comp_id'], chain_id)
                        self.__reg.pA.doAlign()

                        myAlign = self.__reg.pA.getAlignment(chain_id)

                        length = len(myAlign)

                        if length == 0:
                            continue

                        _matched, unmapped, conflict, offset_1, offset_2 = getScoreOfSeqAlign(myAlign)

                        if length == unmapped + conflict:
                            continue

                        self.__reg.native_combined = True  # DAOTHER-8817

                    else:
                        continue

                _ps1 = ps1 if offset_1 == 0 else fillBlankCompIdWithOffset(ps1, offset_1)
                _ps2 = ps2 if offset_2 == 0 else fillBlankCompIdWithOffset(ps2, offset_2)

                if conflict > 0 and (hasLargeSeqGap(_ps1, _ps2)
                                     or (not hasLargeInnerSeqGap(ps1) and hasLargeInnerSeqGap(ps2))):  # DAOTHER-7465
                    _ps2 = self.__compensateLadderHistidinTag__(chain_id, _ps1, _ps2)
                    __ps1, __ps2 = beautifyPolySeq(_ps1, _ps2)
                    _ps1_ = __ps1
                    _ps2_ = __ps2

                    self.__reg.pA.clear()
                    self.__reg.pA.setReferenceSequence(_ps1_['comp_id'], f'REF{chain_id}')
                    self.__reg.pA.addTestSequence(_ps2_['comp_id'], chain_id)
                    self.__reg.pA.doAlign()

                    myAlign = self.__reg.pA.getAlignment(chain_id)

                    length = len(myAlign)

                    __matched, _unmapped, _conflict, _offset_1, _offset_2 = getScoreOfSeqAlign(myAlign)

                    if __matched > 0 and _conflict == 0:  # and len(__ps2['comp_id']) - len(ps2['comp_id']) == conflict:
                        _matched, unmapped = __matched, _unmapped
                        conflict = 0
                        offset_1 = _offset_1
                        offset_2 = _offset_2
                        _ps1 = __ps1
                        _ps2 = __ps2

                if _matched <= conflict + (1 if length > 1 else 0):
                    continue

                seq_align = self.__getSeqAlignCodeWithBlankFill__(length, len(ps1['seq_id']), _matched, conflict, unmapped,
                                                                  offset_1, offset_2, _ps1, _ps2, myAlign)
                seq_align['ref_chain_id'] = chain_id
                seq_align['test_chain_id'] = chain_id2

                if 'auth_seq_id' in _ps1:
                    seq_align['ref_auth_seq_id'] = _ps1['auth_seq_id']

                seq_align_set.append(seq_align)

        if len(seq_align_set) > 0:
            self.__reg.report.sequence_alignment.setItemValue('model_poly_seq_vs_nmr_poly_seq', seq_align_set)

        seq_align_set = []

        for i1, ps1 in enumerate(nmr_poly_seq):
            chain_id = ps1['chain_id']

            if i1 >= LEN_MAJOR_ASYM_ID / 2:  # to process large assembly avoiding forced timeout
                continue

            for i2, ps2 in enumerate(cif_poly_seq):
                chain_id2 = ps2['chain_id']

                if i2 >= LEN_MAJOR_ASYM_ID / 2:  # to process large assembly avoiding forced timeout
                    continue

                self.__reg.pA.clear()
                self.__reg.pA.setReferenceSequence(ps1['comp_id'], f'REF{chain_id}')
                self.__reg.pA.addTestSequence(ps2['auth_comp_id'] if 'auth_comp_id' in ps2 else ps2['comp_id'], chain_id)
                self.__reg.pA.doAlign()

                myAlign = self.__reg.pA.getAlignment(chain_id)

                length = len(myAlign)

                if length == 0:
                    continue

                _matched, unmapped, conflict, offset_1, offset_2 = getScoreOfSeqAlign(myAlign)

                if length == unmapped + conflict:
                    if len(ps2['seq_id']) == 1 and 'alt_comp_id' in ps2 and ps2['alt_comp_id'][0] in ps1['comp_id']:
                        self.__reg.pA.clear()
                        self.__reg.pA.setReferenceSequence(ps1['comp_id'], f'REF{chain_id}')
                        self.__reg.pA.addTestSequence(ps2['alt_comp_id'], chain_id)
                        self.__reg.pA.doAlign()

                        myAlign = self.__reg.pA.getAlignment(chain_id)

                        length = len(myAlign)

                        if length == 0:
                            continue

                        _matched, unmapped, conflict, offset_1, offset_2 = getScoreOfSeqAlign(myAlign)

                        if length == unmapped + conflict:
                            continue

                        self.__reg.native_combined = True  # DAOTHER-8817

                    else:
                        continue

                _ps1 = ps1 if offset_1 == 0 else fillBlankCompIdWithOffset(ps1, offset_1)
                _ps2 = ps2 if offset_2 == 0 else fillBlankCompIdWithOffset(ps2, offset_2)

                if conflict > 0 and (hasLargeSeqGap(_ps1, _ps2)
                                     or (hasLargeInnerSeqGap(ps1) and not hasLargeInnerSeqGap(ps2))):  # DAOTHER-7465
                    _ps1 = self.__compensateLadderHistidinTag__(chain_id, _ps2, _ps1)
                    __ps1, __ps2 = beautifyPolySeq(_ps1, _ps2)
                    _ps1_ = __ps1
                    _ps2_ = __ps2

                    self.__reg.pA.clear()
                    self.__reg.pA.setReferenceSequence(_ps1_['comp_id'], f'REF{chain_id}')
                    self.__reg.pA.addTestSequence(_ps2_['comp_id'], chain_id)
                    self.__reg.pA.doAlign()

                    myAlign = self.__reg.pA.getAlignment(chain_id)

                    length = len(myAlign)

                    __matched, _unmapped, _conflict, _offset_1, _offset_2 = getScoreOfSeqAlign(myAlign)

                    if __matched > 0 and _conflict == 0:  # and len(__ps1['comp_id']) - len(ps1['comp_id']) == conflict:
                        _matched, unmapped = __matched, _unmapped
                        conflict = 0
                        offset_1 = _offset_1
                        offset_2 = _offset_2
                        _ps1 = __ps1
                        _ps2 = __ps2

                if _matched <= conflict + (1 if length > 1 else 0):
                    continue

                seq_align = self.__getSeqAlignCodeWithBlankFill__(length, len(ps1['seq_id']), _matched, conflict, unmapped,
                                                                  offset_1, offset_2, _ps1, _ps2, myAlign)
                seq_align['ref_chain_id'] = chain_id
                seq_align['test_chain_id'] = chain_id2

                if 'auth_seq_id' in _ps2:
                    seq_align['test_auth_seq_id'] = _ps2['auth_seq_id']

                seq_align_set.append(seq_align)

        if len(seq_align_set) > 0:
            self.__reg.report.sequence_alignment.setItemValue('nmr_poly_seq_vs_model_poly_seq', seq_align_set)

        return True

    def __compensateLadderHistidinTag__(self, chain_id: str, ps1: dict, ps2: dict) -> dict:
        """ Compensate ladder-like Histidine tag in polymer sequence 2.
        """

        self.__reg.pA.clear()
        self.__reg.pA.setReferenceSequence(ps1['comp_id'], f'REF{chain_id}')
        self.__reg.pA.addTestSequence(ps2['comp_id'], chain_id)
        self.__reg.pA.doAlign()

        _ps2 = copy.copy(ps2)

        len_ps2 = len(ps2['comp_id'])

        myAlign = self.__reg.pA.getAlignment(chain_id)

        length = len(myAlign)

        _myPr0 = '.'

        idx2 = 0
        for p in range(length):
            myPr0, myPr1 = str(myAlign[p][0]), str(myAlign[p][1])

            if myPr0 == myPr1:
                pass

            elif myPr0 == 'HIS' and myPr1 == '.' and _myPr0 == 'HIS':
                if idx2 < len_ps2:
                    _ps2['comp_id'][idx2] = 'HIS'
                    idx2 += 1

            _myPr0 = myPr0

            if myPr1 != '.':
                while idx2 < len_ps2:
                    if ps2['comp_id'][idx2] == myPr1:
                        idx2 += 1
                        break
                    idx2 += 1

        return _ps2

    def __getSeqAlignCodeWithBlankFill__(self, length: int, ref_length: int, matched: int,  # pylint: disable=no-self-use
                                         conflict: int, unmapped: int,
                                         offset_1: int, offset_2: int,
                                         ps1: dict, ps2: dict, myAlign: list) -> dict:
        """ Return human-readable sequence alignment codes with blank filling if necessary.
        """

        ref_code = getOneLetterCodeCanSequence(ps1['comp_id'])
        test_code = getOneLetterCodeCanSequence(ps2['comp_id'])
        mid_code = getMiddleCode(ref_code, test_code)
        ref_gauge_code = getGaugeCode(ps1['seq_id'])
        test_gauge_code = getGaugeCode(ps2['seq_id'])

        if any((__s1, __s2) for (__s1, __s2, __c1, __c2)
               in zip(ps1['seq_id'], ps2['seq_id'], ps1['comp_id'], ps2['comp_id'])
               if __c1 != '.' and __c2 != '.' and __c1 != __c2):
            len_ps1 = len(ps1['seq_id'])
            len_ps2 = len(ps2['seq_id'])

            seq_id1, seq_id2, comp_id1, comp_id2 = [], [], [], []

            idx1 = idx2 = 0
            for i in range(length):
                myPr0, myPr1 = str(myAlign[i][0]), str(myAlign[i][1])

                if myPr0 != '.':
                    while idx1 < len_ps1:
                        if ps1['comp_id'][idx1] == myPr0:
                            seq_id1.append(ps1['seq_id'][idx1])
                            comp_id1.append(myPr0)
                            idx1 += 1
                            break
                        idx1 += 1
                else:
                    seq_id1.append(None)
                    comp_id1.append('.')

                if myPr1 != '.':
                    while idx2 < len_ps2:
                        if ps2['comp_id'][idx2] == myPr1:
                            seq_id2.append(ps2['seq_id'][idx2])
                            comp_id2.append(myPr1)
                            idx2 += 1
                            break
                        idx2 += 1
                else:
                    seq_id2.append(None)
                    comp_id2.append('.')

            ref_code = getOneLetterCodeCanSequence(comp_id1)
            test_code = getOneLetterCodeCanSequence(comp_id2)
            mid_code = getMiddleCode(ref_code, test_code)
            ref_gauge_code = getGaugeCode(seq_id1, offset_1)
            test_gauge_code = getGaugeCode(seq_id2, offset_2)
            if ' ' in ref_gauge_code:
                for p, g in enumerate(ref_gauge_code):
                    if g == ' ':
                        ref_code = f'{ref_code[0:p]}-{ref_code[p + 1:]}'
            if ' ' in test_gauge_code:
                for p, g in enumerate(test_gauge_code):
                    if g == ' ':
                        test_code = f'{test_code[0:p]}-{test_code[p + 1:]}'

        matched = mid_code.count('|')

        return {'length': ref_length, 'matched': matched, 'conflict': conflict, 'unmapped': unmapped,
                'sequence_coverage': round(float(length - (unmapped + conflict)) / ref_length, 3),
                'ref_seq_id': ps1['seq_id'], 'test_seq_id': ps2['seq_id'],
                'ref_gauge_code': ref_gauge_code, 'ref_code': ref_code, 'mid_code': mid_code,
                'test_code': test_code, 'test_gauge_code': test_gauge_code}

    def __assignCoordPolymerSequence(self) -> bool:
        """ Assign polymer sequences of coordinate file.
        """

        src_id = self.__reg.report.getInputSourceIdOfCoord()

        if src_id < 0:
            return False

        cif_input_source = self.__reg.report.input_sources[src_id]
        cif_input_source_dic = cif_input_source.get()

        cif_file_name = cif_input_source_dic['file_name']

        has_cif_poly_seq = has_key_value(cif_input_source_dic, 'polymer_sequence')

        if not has_cif_poly_seq:
            return False

        __errors = self.__reg.report.getTotalErrors()

        ent_asm_id_map = {}
        ign_chain_ids = []

        if len(self.__reg.star_data) > 0 and isinstance(self.__reg.star_data[0], pynmrstar.Entry):
            input_source = self.__reg.report.input_sources[0]
            input_source_dic = input_source.get()

            file_type = input_source_dic['file_type']

            if file_type == 'nmr-star':

                sf_category = 'assembly'

                try:

                    sf = self.__reg.star_data[0].get_saveframes_by_category(sf_category)[0]

                    lp_category = '_Entity_assembly'

                    lp = sf.get_loop(lp_category)

                    tags = ['ID', 'PDB_chain_ID', 'Experimental_data_reported']

                    if set(tags) & set(lp.tags) == set(tags):

                        dat = lp.get_tag(tags)

                        for row in dat:

                            if row[2] != 'yes':
                                ign_chain_ids.append(row[0])

                            if ',' in row[1]:
                                continue

                            ent_asm_id_map[row[0]] = row[1]

                        if len(ent_asm_id_map) + len(ign_chain_ids) != len(dat):
                            ent_asm_id_map = {}

                except (IndexError, KeyError):
                    pass

        for fileListId in range(self.__reg.file_path_list_len):

            nmr_input_source = self.__reg.report.input_sources[fileListId]
            nmr_input_source_dic = nmr_input_source.get()

            nmr_file_name = nmr_input_source_dic['file_name']

            has_nmr_poly_seq = has_key_value(nmr_input_source_dic, 'polymer_sequence')

            if not has_nmr_poly_seq:
                continue

            seq_align_dic = self.__reg.report.sequence_alignment.get()

            if has_key_value(seq_align_dic, 'model_poly_seq_vs_nmr_poly_seq')\
                    and has_key_value(seq_align_dic, 'nmr_poly_seq_vs_model_poly_seq'):

                cif_poly_seq = cif_input_source_dic['polymer_sequence']
                nmr_poly_seq = nmr_input_source_dic['polymer_sequence']

                if nmr_poly_seq is None:
                    continue

                cif_chains = len(cif_poly_seq)
                nmr_chains = len(nmr_poly_seq)

                # map polymer sequences between coordinate and NMR data using Hungarian algorithm

                m = Munkres()

                # from model to nmr (first trial, never raise a warning or an error)

                mat, indices, lig_chain_ids, lig_chain_ids2 = [], [], [], []

                valid_ent_asm_id_map = fileListId == 0 and len(ent_asm_id_map) > 0\
                    and all('auth_chain_id' in ps1 and ps1['auth_chain_id'] in ent_asm_id_map.values() for ps1 in cif_poly_seq)\
                    and all(ps2['chain_id'] in ent_asm_id_map or ps2['chain_id'] in ign_chain_ids for ps2 in nmr_poly_seq)

                if valid_ent_asm_id_map:
                    for k, v in ent_asm_id_map.items():
                        ps1 = next((ps for ps in cif_poly_seq if ps['auth_chain_id'] == v), None)
                        ps2 = next((ps for ps in nmr_poly_seq if ps['chain_id'] == k), None)
                        if ps1 is None or ps2 is None:
                            valid_ent_asm_id_map = False
                            break

                        result = next((seq_align for seq_align in seq_align_dic['model_poly_seq_vs_nmr_poly_seq']
                                       if seq_align['ref_chain_id'] == v
                                       and seq_align['test_chain_id'] == k), None)
                        if result is None:
                            valid_ent_asm_id_map = False
                            break

                        if result['unmapped'] + result['conflict'] - result['length'] >= 0:
                            valid_ent_asm_id_map = False
                            break

                for ps1 in cif_poly_seq:
                    chain_id = ps1['chain_id']

                    cost = [0 for i in range(nmr_chains)]

                    for ps2 in nmr_poly_seq:
                        chain_id2 = ps2['chain_id']

                        if valid_ent_asm_id_map and (chain_id2 in ign_chain_ids or ent_asm_id_map[chain_id2] != chain_id):
                            continue

                        result = next((seq_align for seq_align in seq_align_dic['model_poly_seq_vs_nmr_poly_seq']
                                       if seq_align['ref_chain_id'] == chain_id
                                       and seq_align['test_chain_id'] == chain_id2), None)

                        if result is not None:
                            cost[nmr_poly_seq.index(ps2)] = result['unmapped'] + result['conflict'] - result['length']
                            if not self.__reg.native_combined and result['length'] >= len(ps1['seq_id']) - result['unmapped']:
                                indices.append((cif_poly_seq.index(ps1), nmr_poly_seq.index(ps2)))

                    mat.append(cost)

                if self.__reg.native_combined:
                    indices = m.compute(mat)

                concatenated_nmr_chain = {}

                for row, col in indices:

                    if mat[row][col] >= 0:

                        if self.__reg.native_combined:
                            continue

                        # DAOTHER-8751
                        has_row = has_col = False
                        for _row, _col in indices:
                            if mat[_row][_col] < 0:
                                if _row == row:
                                    has_row = True
                                if _col == col:
                                    has_col = True

                        if has_row and has_col:
                            continue

                        _cif_chain_ids = [cif_poly_seq[_row]['chain_id'] for _row, _col in indices if col == _col]

                        if len(_cif_chain_ids) > 1:
                            chain_id2 = nmr_poly_seq[col]['chain_id']
                            concatenated_nmr_chain[chain_id2] = _cif_chain_ids

                    chain_id = cif_poly_seq[row]['chain_id']
                    chain_id2 = nmr_poly_seq[col]['chain_id']

                    result = next(seq_align for seq_align in seq_align_dic['model_poly_seq_vs_nmr_poly_seq']
                                  if seq_align['ref_chain_id'] == chain_id and seq_align['test_chain_id'] == chain_id2)
                    _result = next((seq_align for seq_align in seq_align_dic['nmr_poly_seq_vs_model_poly_seq']
                                    if seq_align['ref_chain_id'] == chain_id2 and seq_align['test_chain_id'] == chain_id), None)

                    if result['matched'] == 0\
                       or (result['conflict'] > 0
                           and result['sequence_coverage'] < LOW_SEQ_COVERAGE < float(result['conflict']) / float(result['matched'])):  # noqa: E501, pylint: disable=line-too-long
                        continue

                    ca = {'ref_chain_id': chain_id, 'test_chain_id': chain_id2, 'length': result['length'],
                          'matched': result['matched'], 'conflict': result['conflict'], 'unmapped': result['unmapped'],
                          'sequence_coverage': result['sequence_coverage']}

                    # DAOTHER-8751
                    low_evid_chain_mapping = result['sequence_coverage'] < LOW_SEQ_COVERAGE
                    if low_evid_chain_mapping:
                        low_evid_chain_mapping = False
                        for _row, _col in indices:
                            if mat[_row][_col] >= 0:
                                if _row == row or _col == col:
                                    low_evid_chain_mapping = True

                    auth_chain_id = chain_id
                    if 'auth_chain_id' in cif_poly_seq[row]:
                        auth_chain_id = cif_poly_seq[row]['auth_chain_id']
                        ca['ref_auth_chain_id'] = auth_chain_id

                    ps1 = next(ps for ps in cif_poly_seq if ps['chain_id'] == chain_id)
                    ps2 = next(ps for ps in nmr_poly_seq if ps['chain_id'] == chain_id2)

                    self.__reg.pA.clear()
                    self.__reg.pA.setReferenceSequence(ps1['comp_id'], f'REF{chain_id}')
                    self.__reg.pA.addTestSequence(ps2['comp_id'], chain_id)
                    self.__reg.pA.doAlign()

                    myAlign = self.__reg.pA.getAlignment(chain_id)

                    length = len(myAlign)

                    _matched, unmapped, conflict, offset_1, offset_2 = getScoreOfSeqAlign(myAlign)

                    _ps1 = ps1 if offset_1 == 0 else fillBlankCompIdWithOffset(ps1, offset_1)
                    _ps2 = ps2 if offset_2 == 0 else fillBlankCompIdWithOffset(ps2, offset_2)

                    if conflict == 0:
                        has_inner_gap_1 = hasLargeInnerSeqGap(_ps1)
                        has_inner_gap_2 = hasLargeInnerSeqGap(_ps2)

                        if has_inner_gap_2 and not has_inner_gap_1:
                            _ps2 = fillInnerBlankCompId(_ps2)
                        elif has_inner_gap_1 and not has_inner_gap_2:
                            _ps1 = fillInnerBlankCompId(_ps1)

                    if conflict > 0 and (hasLargeSeqGap(_ps1, _ps2)
                                         or (not hasLargeInnerSeqGap(ps1) and hasLargeInnerSeqGap(ps2))):  # DAOTHER-7465
                        __ps1, __ps2 = beautifyPolySeq(_ps1, _ps2)
                        _ps1 = __ps1
                        _ps2 = __ps2

                        self.__reg.pA.clear()
                        self.__reg.pA.setReferenceSequence(_ps1['comp_id'], f'REF{chain_id}')
                        self.__reg.pA.addTestSequence(_ps2['comp_id'], chain_id)
                        self.__reg.pA.doAlign()

                        myAlign = self.__reg.pA.getAlignment(chain_id)

                        length = len(myAlign)

                        _matched, unmapped, _conflict, _, _ = getScoreOfSeqAlign(myAlign)

                        if _conflict == 0 and len(__ps2['comp_id']) - len(ps2['comp_id']) == conflict:
                            result['conflict'] = 0
                            ps2 = __ps2

                    # update residue name in CS loop to follow CCD replacement (6vu1)
                    if conflict > 0 and ('alt_comp_id' in ps1 or 'alt_comp_id' in ps2) and len(ps1['seq_id']) == len(ps2['seq_id']):
                        for (k1, k2) in zip(['alt_comp_id', 'comp_id'], ['alt_comp_id', 'comp_id']):

                            if k1 == k2 == 'comp_id' or k1 not in ps1 or k2 not in ps2:
                                continue

                            self.__reg.pA.clear()
                            self.__reg.pA.setReferenceSequence(ps1[k1], f'REF{chain_id}')
                            self.__reg.pA.addTestSequence(ps2[k2], chain_id)
                            self.__reg.pA.doAlign()

                            myAlign = self.__reg.pA.getAlignment(chain_id)

                            length = len(myAlign)

                            _matched, unmapped, _conflict, _offset_1, _offset_2 = getScoreOfSeqAlign(myAlign)

                            if _conflict == 0:
                                self.__reg.dpR.updateCompIdInCsLoop(fileListId, ps1, ps2)

                                result['conflict'] = _result['conflict'] = 0
                                ps2['comp_id'] = ps1['comp_id']

                                break

                    if conflict == 0:
                        # resolve unmapped author sequence in CS loop based on sequence alignment (2kxc)
                        self.__reg.dpR.resolveUnmappedAuthSequenceInCsLoop(fileListId, ps1, ps2)

                    ref_code = getOneLetterCodeCanSequence(ps1['comp_id'])
                    test_code = getOneLetterCodeCanSequence(ps2['comp_id'])

                    for r_code, t_code, seq_id, seq_id2 in zip(ref_code, test_code, ps1['seq_id'], ps2['seq_id']):
                        if r_code == 'X' and t_code == 'X':
                            nmr_input_source.updateNonStandardResidueByExptlData(chain_id2, seq_id2, 'coordinate')
                            cif_input_source.updateNonStandardResidueByExptlData(chain_id, seq_id, 'coordinate')

                    if length == _matched == 1:  # to avoid redundant ligand mapping
                        if chain_id in lig_chain_ids or chain_id2 in lig_chain_ids2:
                            continue
                        lig_chain_ids.append(chain_id)
                        lig_chain_ids2.append(chain_id2)

                    if result['unmapped'] > 0 or result['conflict'] > 0:

                        aligned = [True] * length
                        seq_id1, seq_id2 = [], []

                        j = 0
                        for i in range(length):
                            if j < len(ps1['seq_id']) and (str(myAlign[i][0]) != '.' or ps1['comp_id'][j] == '.'):
                                seq_id1.append(ps1['seq_id'][j])
                                j += 1
                            else:
                                seq_id1.append(None)

                        j = 0
                        for i in range(length):
                            if j < len(ps2['seq_id']) and (str(myAlign[i][1]) != '.' or ps2['comp_id'][j] == '.'):
                                seq_id2.append(ps2['seq_id'][j])
                                j += 1
                            else:
                                seq_id2.append(None)

                        for i in range(length):
                            myPr0, myPr1 = str(myAlign[i][0]), str(myAlign[i][1])

                            if myPr0 == '.' or myPr1 == '.':
                                aligned[i] = False
                            elif myPr0 != myPr1:
                                pass
                            else:
                                break

                        for i in reversed(range(length)):
                            myPr0, myPr1 = str(myAlign[i][0]), str(myAlign[i][1])

                            if myPr0 == '.' or myPr1 == '.':
                                aligned[i] = False
                            elif myPr0 != myPr1:
                                pass
                            else:
                                break

                        if not self.__reg.native_combined:

                            _conflicts = 0

                            for i in range(length):
                                cif_comp_id, nmr_comp_id = str(myAlign[i][0]), str(myAlign[i][1])

                                if cif_comp_id == nmr_comp_id:
                                    continue

                                if nmr_comp_id == '.' and cif_comp_id != '.':
                                    pass

                                elif nmr_comp_id != cif_comp_id and aligned[i]:
                                    _conflicts += 1

                            if _conflicts > ca['unmapped'] and ca['sequence_coverage'] < MIN_SEQ_COVERAGE_W_CONFLICT:
                                continue

                            # DAOTHER-7825 (2lyw)
                            if _conflicts + offset_1 > _matched and ca['sequence_coverage'] < LOW_SEQ_COVERAGE:
                                if not low_evid_chain_mapping:  # DAOTHER-8751
                                    continue

                        unmapped, conflict = [], []

                        for i in range(length):
                            cif_comp_id, nmr_comp_id = str(myAlign[i][0]), str(myAlign[i][1])

                            if cif_comp_id == nmr_comp_id:
                                continue

                            if nmr_comp_id == '.' and cif_comp_id != '.':
                                unmapped.append({'ref_seq_id': seq_id1[i], 'ref_comp_id': cif_comp_id})

                            elif nmr_comp_id != cif_comp_id and aligned[i]:
                                conflict.append({'ref_seq_id': seq_id1[i], 'ref_comp_id': cif_comp_id,
                                                 'test_seq_id': seq_id2[i], 'test_comp_id': nmr_comp_id})

                        if len(unmapped) > 0:
                            ca['unmapped_sequence'] = unmapped

                        if len(conflict) > 0:
                            ca['conflict_sequence'] = conflict
                            ca['conflict'] = len(conflict)
                            ca['unmapped'] = ca['unmapped'] - len(conflict)
                            if ca['unmapped'] < 0:
                                ca['conflict'] -= ca['unmapped']
                                ca['unmapped'] = 0

                            result['conflict'] = ca['conflict']
                            result['unmapped'] = ca['unmapped']

                            if _result is not None:
                                _result['conflict'] = ca['conflict']
                                _result['unmapped'] = ca['unmapped']

                # from nmr to model

                ca_idx = 0

                mat, indices, lig_chain_ids, lig_chain_ids2 = [], [], [], []

                for ps1 in nmr_poly_seq:
                    chain_id = ps1['chain_id']

                    cost = [0 for i in range(cif_chains)]

                    for ps2 in cif_poly_seq:
                        chain_id2 = ps2['chain_id']

                        if valid_ent_asm_id_map and (chain_id in ign_chain_ids or ent_asm_id_map[chain_id] != chain_id2):
                            continue

                        result = next((seq_align for seq_align in seq_align_dic['nmr_poly_seq_vs_model_poly_seq']
                                       if seq_align['ref_chain_id'] == chain_id and seq_align['test_chain_id'] == chain_id2), None)

                        if result is not None:
                            cost[cif_poly_seq.index(ps2)] = result['unmapped'] + result['conflict'] - result['length']
                            if not self.__reg.native_combined and result['length'] >= len(ps2['seq_id']) - result['unmapped']:
                                indices.append((nmr_poly_seq.index(ps1), cif_poly_seq.index(ps2)))

                    mat.append(cost)

                if self.__reg.native_combined:
                    indices = m.compute(mat)

                chain_assign = []

                for row, col in indices:

                    if mat[row][col] >= 0:

                        if self.__reg.native_combined:
                            continue

                        # DAOTHER-8751
                        has_row = has_col = False
                        for _row, _col in indices:
                            if mat[_row][_col] < 0:
                                if _row == row:
                                    has_row = True
                                if _col == col:
                                    has_col = True

                        if has_row and has_col:
                            continue

                    chain_id = nmr_poly_seq[row]['chain_id']
                    chain_id2 = cif_poly_seq[col]['chain_id']

                    result = next(seq_align for seq_align in seq_align_dic['nmr_poly_seq_vs_model_poly_seq']
                                  if seq_align['ref_chain_id'] == chain_id and seq_align['test_chain_id'] == chain_id2)
                    _result = next((seq_align for seq_align in seq_align_dic['model_poly_seq_vs_nmr_poly_seq']
                                    if seq_align['ref_chain_id'] == chain_id2 and seq_align['test_chain_id'] == chain_id), None)

                    if result['matched'] == 0\
                       or (result['conflict'] > 0
                           and result['sequence_coverage'] < LOW_SEQ_COVERAGE < float(result['conflict']) / float(result['matched'])):  # noqa: E501, pylint: disable=line-too-long
                        continue

                    ca = {'ref_chain_id': chain_id, 'test_chain_id': chain_id2, 'length': result['length'],
                          'matched': result['matched'], 'conflict': result['conflict'], 'unmapped': result['unmapped'],
                          'sequence_coverage': result['sequence_coverage']}
                    if 'auth_chain_id' in cif_poly_seq[col]:
                        ca['test_auth_chain_id'] = cif_poly_seq[col]['auth_chain_id']
                    else:
                        ca['test_auth_chain_id'] = chain_id2

                    # DAOTHER-8751
                    low_evid_chain_mapping = result['sequence_coverage'] < LOW_SEQ_COVERAGE
                    if low_evid_chain_mapping:
                        low_evid_chain_mapping = False
                        for _row, _col in indices:
                            if mat[_row][_col] >= 0:
                                if _row == row or _col == col:
                                    low_evid_chain_mapping = True

                    auth_chain_id2 = chain_id2
                    if 'auth_chain_id' in cif_poly_seq[col]:
                        auth_chain_id2 = cif_poly_seq[col]['auth_chain_id']
                        ca['test_auth_chain_id'] = auth_chain_id2

                    ps1 = next(ps for ps in nmr_poly_seq if ps['chain_id'] == chain_id)
                    ps2 = next(ps for ps in cif_poly_seq if ps['chain_id'] == chain_id2)

                    self.__reg.pA.clear()
                    self.__reg.pA.setReferenceSequence(ps1['comp_id'], f'REF{chain_id}')
                    self.__reg.pA.addTestSequence(ps2['comp_id'], chain_id)
                    self.__reg.pA.doAlign()

                    myAlign = self.__reg.pA.getAlignment(chain_id)

                    length = len(myAlign)

                    _matched, unmapped, conflict, offset_1, offset_2 = getScoreOfSeqAlign(myAlign)

                    if conflict > 0 and any(True for c in ps2['comp_id'] if len(c) > 3) and 'alt_comp_id' in ps2:
                        self.__reg.pA.addTestSequence(ps2['alt_comp_id'], chain_id)
                        self.__reg.pA.doAlign()

                        myAlign = self.__reg.pA.getAlignment(chain_id)

                        length = len(myAlign)

                        _matched, unmapped, conflict, offset_1, offset_2 = getScoreOfSeqAlign(myAlign)

                        if conflict > 0:
                            self.__reg.pA.addTestSequence(ps2['comp_id'], chain_id)
                            self.__reg.pA.doAlign()

                            myAlign = self.__reg.pA.getAlignment(chain_id)

                            length = len(myAlign)

                            _matched, unmapped, conflict, offset_1, offset_2 = getScoreOfSeqAlign(myAlign)

                    _ps1 = ps1 if offset_1 == 0 else fillBlankCompIdWithOffset(ps1, offset_1)
                    _ps2 = ps2 if offset_2 == 0 else fillBlankCompIdWithOffset(ps2, offset_2)

                    if conflict == 0:
                        has_inner_gap_1 = hasLargeInnerSeqGap(_ps1)
                        has_inner_gap_2 = hasLargeInnerSeqGap(_ps2)

                        if has_inner_gap_2 and not has_inner_gap_1:
                            _ps2 = fillInnerBlankCompId(_ps2)
                        elif has_inner_gap_1 and not has_inner_gap_2:
                            _ps1 = fillInnerBlankCompId(_ps1)

                    if conflict > 0 and (hasLargeSeqGap(_ps1, _ps2)
                                         or (hasLargeInnerSeqGap(ps1) and not hasLargeInnerSeqGap(ps2))):  # DAOTHER-7465
                        __ps1, __ps2 = beautifyPolySeq(_ps1, _ps2)
                        _ps1 = __ps1
                        _ps2 = __ps2

                        self.__reg.pA.clear()
                        self.__reg.pA.setReferenceSequence(_ps1['comp_id'], f'REF{chain_id}')
                        self.__reg.pA.addTestSequence(_ps2['comp_id'], chain_id)
                        self.__reg.pA.doAlign()

                        myAlign = self.__reg.pA.getAlignment(chain_id)

                        length = len(myAlign)

                        _matched, unmapped, _conflict, _, _ = getScoreOfSeqAlign(myAlign)

                        if _conflict == 0 and len(__ps1['comp_id']) - len(ps1['comp_id']) == conflict:
                            result['conflict'] = 0
                            ps1 = __ps1

                    if conflict > 0 and 'gap_in_auth_seq' in _ps2 and _ps2['gap_in_auth_seq'] and 'auth_seq_id' in _ps2:
                        __ps1 = copy.deepcopy(_ps1)
                        for p in range(len(_ps2['auth_seq_id']) - 1):
                            s_p = _ps2['auth_seq_id'][p]
                            s_q = _ps2['auth_seq_id'][p + 1]
                            if None in (s_p, s_q) or s_p + 1 == s_q:
                                continue
                            for s_o in range(s_p + 1, s_q):
                                if s_o in __ps1['seq_id']:
                                    idx = __ps1['seq_id'].index(s_o)
                                    if __ps1['comp_id'][idx] in EMPTY_VALUE:
                                        __ps1['seq_id'].pop(idx)
                                        __ps1['comp_id'].pop(idx)

                        if len(_ps1['seq_id']) != len(__ps1['seq_id']):
                            __ps1, __ps2 = beautifyPolySeq(__ps1, _ps2)
                            _ps1 = __ps1
                            _ps2 = __ps2

                            self.__reg.pA.clear()
                            self.__reg.pA.setReferenceSequence(_ps1['comp_id'], f'REF{chain_id}')
                            self.__reg.pA.addTestSequence(_ps2['comp_id'], chain_id)
                            self.__reg.pA.doAlign()

                            myAlign = self.__reg.pA.getAlignment(chain_id)

                            length = len(myAlign)

                            _matched, unmapped, _conflict, _, _ = getScoreOfSeqAlign(myAlign)

                            if _conflict == 0:
                                result['conflict'] = 0
                                ps1 = __ps1

                    if length == _matched == 1:  # to avoid redundant ligand mapping
                        if chain_id in lig_chain_ids or chain_id2 in lig_chain_ids2:
                            continue
                        lig_chain_ids.append(chain_id)
                        lig_chain_ids2.append(chain_id2)

                    if result['unmapped'] > 0 or result['conflict'] > 0:

                        aligned = [True] * length
                        seq_id1, seq_id2 = [], []

                        j = 0
                        for i in range(length):
                            if j < len(ps1['seq_id']) and (str(myAlign[i][0]) != '.' or ps1['comp_id'][j] == '.'):  # DAOTHER-7421
                                seq_id1.append(ps1['seq_id'][j])
                                j += 1
                            else:
                                seq_id1.append(None)

                        j = 0
                        for i in range(length):
                            if j < len(ps2['seq_id']) and (str(myAlign[i][1]) != '.' or ps2['comp_id'][j] == '.'):  # DAOTHER-7421
                                seq_id2.append(ps2['seq_id'][j])
                                j += 1
                            else:
                                seq_id2.append(None)

                        for i in range(length):
                            myPr0, myPr1 = str(myAlign[i][0]), str(myAlign[i][1])

                            if myPr0 == '.' or myPr1 == '.':
                                aligned[i] = False
                            elif myPr0 != myPr1:
                                pass
                            else:
                                break

                        for i in reversed(range(length)):
                            myPr0, myPr1 = str(myAlign[i][0]), str(myAlign[i][1])

                            if myPr0 == '.' or myPr1 == '.':
                                aligned[i] = False
                            elif myPr0 != myPr1:
                                pass
                            else:
                                break

                        for i in range(length):
                            myPr = myAlign[i]
                            if aligned[i]:
                                if str(myPr[0]) == '.':
                                    if (seq_id2[i] is not None)\
                                       and ((i > 0 and seq_id2[i - 1] is not None
                                             and seq_id2[i - 1] + 1 == seq_id2[i])
                                            or (i + 1 < len(seq_id2) and seq_id2[i + 1] is not None
                                                and seq_id2[i + 1] - 1 == seq_id2[i])):
                                        aligned[i] = False
                                if str(myPr[1]) == '.':
                                    if (seq_id1[i] is not None)\
                                       and ((i > 0 and seq_id1[i - 1] is not None
                                             and seq_id1[i - 1] + 1 == seq_id1[i])
                                            or (i + 1 < len(seq_id1) and seq_id1[i + 1] is not None
                                                and seq_id1[i + 1] - 1 == seq_id1[i])):
                                        aligned[i] = False

                        if not self.__reg.native_combined:

                            _conflicts = 0

                            for i in range(length):
                                nmr_comp_id, cif_comp_id = str(myAlign[i][0]), str(myAlign[i][1])

                                if nmr_comp_id == cif_comp_id:
                                    continue

                                if cif_comp_id == '.' and nmr_comp_id != '.':
                                    pass

                                elif cif_comp_id != nmr_comp_id and aligned[i]:
                                    _conflicts += 1

                            if _conflicts > ca['unmapped'] and ca['sequence_coverage'] < MIN_SEQ_COVERAGE_W_CONFLICT:
                                continue

                            # DAOTHER-7825 (2lyw)
                            if _conflicts + offset_1 > _matched and ca['sequence_coverage'] < LOW_SEQ_COVERAGE:
                                if not low_evid_chain_mapping:  # DAOTHER-8751
                                    continue

                        unmapped, conflict = [], []

                        for i in range(length):
                            nmr_comp_id, cif_comp_id = str(myAlign[i][0]), str(myAlign[i][1])

                            if nmr_comp_id == cif_comp_id:
                                continue

                            if cif_comp_id == '.' and nmr_comp_id != '.':

                                _seq_id1 = seq_id1[i] - offset_1 if seq_id1[i] is not None else None

                                unmapped.append({'ref_seq_id': _seq_id1, 'ref_comp_id': nmr_comp_id})

                                if not aligned[i]:

                                    if self.__reg.native_combined or chain_id not in concatenated_nmr_chain\
                                       or chain_id2 not in concatenated_nmr_chain[chain_id]:

                                        warn = f"{chain_id}:{_seq_id1}:{nmr_comp_id} is not present "\
                                            f"in the coordinates (chain_id {chain_id2}). {INSTRUCTION_FOR_FULL_SEQUENCE}"

                                        self.__reg.suspended_warnings_for_lazy_eval.append({'sequence_mismatch':
                                                                                            {'ca_idx': ca_idx,
                                                                                             'file_name': nmr_file_name,
                                                                                             'description': warn}})

                                        if self.__reg.verbose:
                                            self.__reg.log.write(f"+{self.__class_name__}.__assignCoordPolymerSequence() "
                                                                 f"++ Warning  - {warn}\n")

                            elif cif_comp_id != nmr_comp_id and aligned[i]:

                                _seq_id1 = seq_id1[i] - offset_1 if seq_id1[i] is not None else None
                                _seq_id2 = seq_id2[i] - offset_2 if seq_id2[i] is not None else None

                                if _seq_id1 is None and _seq_id2 is None:
                                    continue

                                if _seq_id1 == _seq_id2 and len(nmr_comp_id) + 1 == len(cif_comp_id)\
                                   and cif_comp_id.startswith(nmr_comp_id) and cif_comp_id[-1] in ('5', '3'):
                                    continue

                                conflict.append({'ref_seq_id': _seq_id1, 'ref_comp_id': nmr_comp_id,
                                                 'test_seq_id': _seq_id2, 'test_comp_id': cif_comp_id})

                                try:
                                    label_seq_id = _seq_id2
                                    auth_seq_id = ps2['auth_seq_id'][ps2['seq_id'].index(_seq_id2)]
                                except (KeyError, IndexError, ValueError):
                                    label_seq_id = _seq_id2
                                    auth_seq_id = label_seq_id
                                cif_seq_code = f"{chain_id2}:{label_seq_id}:{cif_comp_id}"
                                if cif_comp_id == '.':
                                    cif_seq_code += ', insertion error'
                                    # DAOTHER-9644: skip insertion error due to truncated loop
                                    cif_ps = next(cif_ps for cif_ps in cif_poly_seq if cif_ps['chain_id'] == chain_id2)
                                    if label_seq_id is None:
                                        if self.__reg.caC is not None and 'missing_polymer_linkage' in self.__reg.caC\
                                           and any(True for mis in self.__reg.caC['missing_polymer_linkage']
                                                   if mis['auth_chain_id'] == (cif_ps['auth_chain_id'] if 'auth_chain_id' in cif_ps
                                                                               else cif_ps['chain_id'])):
                                            continue
                                    else:
                                        if self.__reg.caC is not None and 'missing_polymer_linkage' in self.__reg.caC\
                                           and any(True for mis in self.__reg.caC['missing_polymer_linkage']
                                                   if mis['auth_chain_id'] == (cif_ps['auth_chain_id'] if 'auth_chain_id' in cif_ps
                                                                               else cif_ps['chain_id'])
                                                   and mis['auth_seq_id_1'] < auth_seq_id < mis['auth_seq_id_2']):
                                            continue
                                nmr_seq_code = f"{chain_id}:{_seq_id1}:{nmr_comp_id}"
                                if nmr_comp_id == '.':
                                    nmr_seq_code += ', insertion error'

                                if cif_comp_id != '.':
                                    if chain_id2 != auth_chain_id2 or auth_seq_id != label_seq_id:
                                        cif_seq_code += f", or {auth_chain_id2}:{auth_seq_id}:{cif_comp_id} "\
                                            "in author sequence scheme"

                                err = f"Sequence alignment error between the NMR data ({nmr_seq_code}) "\
                                    f"and the coordinate ({cif_seq_code}). "\
                                    "Please verify the two sequences and re-upload the correct file(s)."

                                if self.__reg.tolerant_seq_align and self.__reg.dpV.equalsToRepCompId(cif_comp_id, nmr_comp_id):
                                    self.__reg.suspended_warnings_for_lazy_eval.append({'sequence_mismatch':
                                                                                        {'ca_idx': ca_idx,
                                                                                         'file_name': nmr_file_name,
                                                                                         'description': err}})

                                    if self.__reg.verbose:
                                        self.__reg.log.write(f"+{self.__class_name__}.__assignCoordPolymerSequence() "
                                                             f"++ Warning  - {err}\n")

                                elif not self.__reg.annotation_mode:
                                    self.__reg.suspended_errors_for_lazy_eval.append({'sequence_mismatch':
                                                                                      {'ca_idx': ca_idx,
                                                                                       'file_name': nmr_file_name,
                                                                                       'description': err}})

                                    if self.__reg.verbose:
                                        self.__reg.log.write(f"+{self.__class_name__}.__assignCoordPolymerSequence() "
                                                             f"++ Error  - {err}\n")

                        if len(unmapped) > 0:
                            ca['unmapped_sequence'] = unmapped

                        if len(conflict) > 0:
                            ca['conflict_sequence'] = conflict
                            ca['conflict'] = len(conflict)
                            ca['unmapped'] = ca['unmapped'] - len(conflict)
                            if ca['unmapped'] < 0:
                                ca['conflict'] -= ca['unmapped']
                                ca['unmapped'] = 0

                            result['conflict'] = ca['conflict']
                            result['unmapped'] = ca['unmapped']

                            if _result is not None:
                                _result['conflict'] = ca['conflict']
                                _result['unmapped'] = ca['unmapped']

                    chain_assign.append(ca)
                    ca_idx += 1

                if len(chain_assign) > 0 and fileListId == 0:

                    if len(cif_poly_seq) > 1:

                        if len(self.__reg.suspended_errors_for_lazy_eval) + len(self.__reg.suspended_warnings_for_lazy_eval) > 0:

                            _del_ca_idx = []

                            for ca_idx, ca in enumerate(chain_assign):

                                if ca['conflict'] == 0:
                                    continue

                                ref_chain_id = ca['ref_chain_id']
                                test_chain_id = ca['test_chain_id']

                                if any(True for _ca in chain_assign
                                       if ((_ca['ref_chain_id'] == ref_chain_id and _ca['test_chain_id'] != test_chain_id)
                                           or (_ca['ref_chain_id'] != ref_chain_id and _ca['test_chain_id'] == test_chain_id))
                                       and _ca['conflict'] == 0):
                                    _del_ca_idx.append(ca_idx)

                            if len(_del_ca_idx) > 0:
                                for ca_idx in reversed(_del_ca_idx):
                                    del chain_assign[ca_idx]
                                if len(self.__reg.suspended_errors_for_lazy_eval) > 0:
                                    _del_msg_idx = set()
                                    for msg_idx, msg in enumerate(self.__reg.suspended_errors_for_lazy_eval):
                                        for k, v in msg.items():
                                            if v['ca_idx'] in _del_ca_idx:
                                                _del_msg_idx.add(msg_idx)
                                    if len(_del_msg_idx) > 0:
                                        for msg_idx in reversed(list(_del_msg_idx)):
                                            if msg_idx < len(self.__reg.suspended_errors_for_lazy_eval):
                                                del self.__reg.suspended_errors_for_lazy_eval[msg_idx]
                                if len(self.__reg.suspended_warnings_for_lazy_eval) > 0:
                                    _del_msg_idx = set()
                                    for msg_idx, msg in enumerate(self.__reg.suspended_warnings_for_lazy_eval):
                                        for k, v in msg.items():
                                            if v['ca_idx'] in _del_ca_idx:
                                                _del_msg_idx.add(msg_idx)
                                    if len(_del_msg_idx) > 0:
                                        for msg_idx in reversed(list(_del_msg_idx)):
                                            if msg_idx < len(self.__reg.suspended_warnings_for_lazy_eval):
                                                del self.__reg.suspended_warnings_for_lazy_eval[msg_idx]

                        if any(True for ps in cif_poly_seq if 'identical_chain_id' in ps):

                            _chain_assign = chain_assign.copy()

                            for ca in _chain_assign:

                                if ca['conflict'] > 0:
                                    continue

                                _chain_id = ca['test_chain_id']
                                _auth_chain_id = ca.get('test_auth_chain_id')

                                if _chain_id in lig_chain_ids2:
                                    continue

                                try:
                                    identity = next(ps['identical_chain_id'] for ps in cif_poly_seq
                                                    if ps['chain_id'] == _chain_id and 'identical_chain_id' in ps)

                                    for _chain_id in identity:

                                        if not any(True for _ca in chain_assign if _ca['test_chain_id'] == _chain_id):
                                            _ca = ca.copy()
                                            _ca['test_chain_id'] = _chain_id
                                            if _auth_chain_id is not None:
                                                _ca['test_auth_chain_id'] = _auth_chain_id

                                            if _ca['ref_chain_id'] in lig_chain_ids:
                                                continue

                                            chain_assign.append(_ca)

                                except StopIteration:
                                    pass

                            if self.__reg.combined_mode:  # DAOTHER-9660, issue #2
                                nmr_ident_chain_id = set()
                                cif_ident_chain_id = set()

                                for ca in chain_assign:

                                    if ca['conflict'] > 0:
                                        continue

                                    nmr_chain_id = ca['ref_chain_id']
                                    cif_chain_id = ca['test_chain_id']

                                    try:
                                        identity = next(ps['identical_chain_id'] for ps in nmr_poly_seq
                                                        if ps['chain_id'] == nmr_chain_id and 'identical_chain_id' in ps)
                                        nmr_ident_chain_id.add(nmr_chain_id)
                                    except StopIteration:
                                        pass

                                    try:
                                        identity = next(ps['identical_chain_id'] for ps in cif_poly_seq
                                                        if ps['chain_id'] == cif_chain_id and 'identical_chain_id' in ps)
                                        cif_ident_chain_id.add(cif_chain_id)
                                    except StopIteration:
                                        pass

                                if len(cif_ident_chain_id) == len(nmr_ident_chain_id) and len(cif_ident_chain_id) > 1:
                                    _chain_assign = chain_assign.copy()

                                    nmr_mapped_chain_id, cif_mapped_chain_id = [], []

                                    for ca in _chain_assign:

                                        if ca['conflict'] > 0:
                                            continue

                                        nmr_chain_id = ca['ref_chain_id']
                                        cif_chain_id = ca['test_chain_id']

                                        if nmr_chain_id not in nmr_ident_chain_id or cif_chain_id not in cif_ident_chain_id:
                                            continue

                                        if nmr_chain_id not in nmr_mapped_chain_id and cif_chain_id not in cif_mapped_chain_id:
                                            nmr_mapped_chain_id.append(nmr_chain_id)
                                            cif_mapped_chain_id.append(cif_chain_id)

                                        else:
                                            chain_assign.remove(ca)

                    self.__reg.report.chain_assignment.setItemValue('nmr_poly_seq_vs_model_poly_seq', chain_assign)

                    if len(self.__reg.suspended_errors_for_lazy_eval) > 0:
                        for msg in self.__reg.suspended_errors_for_lazy_eval:
                            for k, v in msg.items():
                                if 'ca_idx' in v:
                                    del v['ca_idx']
                                self.__reg.report.error.appendDescription(k, v)
                        self.__reg.suspended_errors_for_lazy_eval.clear()

                    if len(self.__reg.suspended_warnings_for_lazy_eval) > 0:
                        for msg in self.__reg.suspended_warnings_for_lazy_eval:
                            for k, v in msg.items():
                                if 'ca_idx' in v:
                                    del v['ca_idx']
                                self.__reg.report.warning.appendDescription(k, v)
                        self.__reg.suspended_warnings_for_lazy_eval.clear()

                # from model to nmr (final)

                ca_idx = 0

                mat, indices = [], []

                for ps1 in cif_poly_seq:
                    chain_id = ps1['chain_id']

                    cost = [0 for i in range(nmr_chains)]

                    for ps2 in nmr_poly_seq:
                        chain_id2 = ps2['chain_id']

                        if valid_ent_asm_id_map and (chain_id2 in ign_chain_ids or ent_asm_id_map[chain_id2] != chain_id):
                            continue

                        result = next((seq_align for seq_align in seq_align_dic['model_poly_seq_vs_nmr_poly_seq']
                                       if seq_align['ref_chain_id'] == chain_id
                                       and seq_align['test_chain_id'] == chain_id2), None)

                        if result is not None:
                            cost[nmr_poly_seq.index(ps2)] = result['unmapped'] + result['conflict'] - result['length']
                            if not self.__reg.native_combined and result['length'] >= len(ps1['seq_id']) - result['unmapped']:
                                indices.append((cif_poly_seq.index(ps1), nmr_poly_seq.index(ps2)))

                    mat.append(cost)

                if self.__reg.native_combined:
                    indices = m.compute(mat)

                chain_assign = []

                concatenated_nmr_chain = {}

                for row, col in indices:

                    if mat[row][col] >= 0:

                        if self.__reg.native_combined:
                            continue

                        # DAOTHER-8751
                        has_row = has_col = False
                        for _row, _col in indices:
                            if mat[_row][_col] < 0:
                                if _row == row:
                                    has_row = True
                                if _col == col:
                                    has_col = True

                        if has_row and has_col:
                            continue

                        _cif_chain_ids = [cif_poly_seq[_row]['chain_id'] for _row, _col in indices if col == _col]

                        if len(_cif_chain_ids) > 1:
                            chain_id2 = nmr_poly_seq[col]['chain_id']
                            concatenated_nmr_chain[chain_id2] = _cif_chain_ids

                            warn = f"The chain ID {chain_id2!r} of the sequences in the NMR data "\
                                f"will be re-assigned to the chain IDs {_cif_chain_ids} in the coordinates during biocuration."

                            self.__reg.report.warning.appendDescription('concatenated_sequence',
                                                                        {'file_name': nmr_file_name, 'description': warn})

                            if self.__reg.verbose:
                                self.__reg.log.write(f"+{self.__class_name__}.__assignCoordPolymerSequence() "
                                                     f"++ Warning  - {warn}\n")

                    chain_id = cif_poly_seq[row]['chain_id']
                    chain_id2 = nmr_poly_seq[col]['chain_id']

                    result = next(seq_align for seq_align in seq_align_dic['model_poly_seq_vs_nmr_poly_seq']
                                  if seq_align['ref_chain_id'] == chain_id and seq_align['test_chain_id'] == chain_id2)
                    _result = next((seq_align for seq_align in seq_align_dic['nmr_poly_seq_vs_model_poly_seq']
                                    if seq_align['ref_chain_id'] == chain_id2 and seq_align['test_chain_id'] == chain_id), None)

                    if result['matched'] == 0\
                       or (result['conflict'] > 0
                           and result['sequence_coverage'] < LOW_SEQ_COVERAGE < float(result['conflict']) / float(result['matched'])):  # noqa: E501, pylint: disable=line-too-long
                        continue

                    ca = {'ref_chain_id': chain_id, 'test_chain_id': chain_id2, 'length': result['length'],
                          'matched': result['matched'], 'conflict': result['conflict'], 'unmapped': result['unmapped'],
                          'sequence_coverage': result['sequence_coverage']}
                    if 'auth_chain_id' in cif_poly_seq[row]:
                        ca['ref_auth_chain_id'] = cif_poly_seq[row]['auth_chain_id']
                    else:
                        ca['ref_auth_chain_id'] = chain_id

                    # DAOTHER-8751
                    low_evid_chain_mapping = result['sequence_coverage'] < LOW_SEQ_COVERAGE
                    if low_evid_chain_mapping:
                        low_evid_chain_mapping = False
                        for _row, _col in indices:
                            if mat[_row][_col] >= 0:
                                if _row == row or _col == col:
                                    low_evid_chain_mapping = True

                    auth_chain_id = chain_id
                    if 'auth_chain_id' in cif_poly_seq[row]:
                        auth_chain_id = cif_poly_seq[row]['auth_chain_id']
                        ca['ref_auth_chain_id'] = auth_chain_id

                    ps1 = next(ps for ps in cif_poly_seq if ps['chain_id'] == chain_id)
                    ps2 = next(ps for ps in nmr_poly_seq if ps['chain_id'] == chain_id2)

                    self.__reg.pA.clear()
                    self.__reg.pA.setReferenceSequence(ps1['comp_id'], f'REF{chain_id}')
                    self.__reg.pA.addTestSequence(ps2['comp_id'], chain_id)
                    self.__reg.pA.doAlign()

                    myAlign = self.__reg.pA.getAlignment(chain_id)

                    length = len(myAlign)

                    _matched, unmapped, conflict, offset_1, offset_2 = getScoreOfSeqAlign(myAlign)

                    if conflict > 0 and any(True for c in ps1['comp_id'] if len(c) > 3) and 'alt_comp_id' in ps1:
                        self.__reg.pA.setReferenceSequence(ps1['alt_comp_id'], f'REF{chain_id}')
                        self.__reg.pA.doAlign()

                        myAlign = self.__reg.pA.getAlignment(chain_id)

                        length = len(myAlign)

                        _matched, unmapped, conflict, offset_1, offset_2 = getScoreOfSeqAlign(myAlign)

                        if conflict > 0:
                            self.__reg.pA.setReferenceSequence(ps1['comp_id'], f'REF{chain_id}')
                            self.__reg.pA.doAlign()

                            myAlign = self.__reg.pA.getAlignment(chain_id)

                            length = len(myAlign)

                            _matched, unmapped, conflict, offset_1, offset_2 = getScoreOfSeqAlign(myAlign)

                    _ps1 = ps1 if offset_1 == 0 else fillBlankCompIdWithOffset(ps1, offset_1)
                    _ps2 = ps2 if offset_2 == 0 else fillBlankCompIdWithOffset(ps2, offset_2)

                    if conflict == 0:
                        has_inner_gap_1 = hasLargeInnerSeqGap(_ps1)
                        has_inner_gap_2 = hasLargeInnerSeqGap(_ps2)

                        if has_inner_gap_2 and not has_inner_gap_1:
                            _ps2 = fillInnerBlankCompId(_ps2)
                        elif has_inner_gap_1 and not has_inner_gap_2:
                            _ps1 = fillInnerBlankCompId(_ps1)

                    if conflict > 0 and (hasLargeSeqGap(_ps1, _ps2)
                                         or (not hasLargeInnerSeqGap(ps1) and hasLargeInnerSeqGap(ps2))):  # DAOTHER-7465
                        __ps1, __ps2 = beautifyPolySeq(_ps1, _ps2)
                        _ps1 = __ps1
                        _ps2 = __ps2

                        self.__reg.pA.clear()
                        self.__reg.pA.setReferenceSequence(_ps1['comp_id'], f'REF{chain_id}')
                        self.__reg.pA.addTestSequence(_ps2['comp_id'], chain_id)
                        self.__reg.pA.doAlign()

                        myAlign = self.__reg.pA.getAlignment(chain_id)

                        length = len(myAlign)

                        _matched, unmapped, _conflict, _, _ = getScoreOfSeqAlign(myAlign)

                        if _conflict == 0 and len(__ps2['comp_id']) - len(ps2['comp_id']) == conflict:
                            result['conflict'] = 0
                            ps2 = __ps2

                    if conflict > 0 and 'gap_in_auth_seq' in _ps1 and _ps1['gap_in_auth_seq'] and 'auth_seq_id' in _ps1:
                        __ps2 = copy.deepcopy(_ps2)
                        for p in range(len(_ps1['auth_seq_id']) - 1):
                            s_p = _ps1['auth_seq_id'][p]
                            s_q = _ps1['auth_seq_id'][p + 1]
                            if None in (s_p, s_q) or s_p + 1 == s_q:
                                continue
                            for s_o in range(s_p + 1, s_q):
                                if s_o in __ps2['seq_id']:
                                    idx = __ps2['seq_id'].index(s_o)
                                    if __ps2['comp_id'][idx] in EMPTY_VALUE:
                                        __ps2['seq_id'].pop(idx)
                                        __ps2['comp_id'].pop(idx)

                        if len(_ps2['seq_id']) != len(__ps2['seq_id']):
                            __ps1, __ps2 = beautifyPolySeq(_ps1, __ps2)
                            _ps1 = __ps1
                            _ps2 = __ps2

                            self.__reg.pA.clear()
                            self.__reg.pA.setReferenceSequence(_ps1['comp_id'], f'REF{chain_id}')
                            self.__reg.pA.addTestSequence(_ps2['comp_id'], chain_id)
                            self.__reg.pA.doAlign()

                            myAlign = self.__reg.pA.getAlignment(chain_id)

                            length = len(myAlign)

                            _matched, unmapped, _conflict, _, _ = getScoreOfSeqAlign(myAlign)

                            if _conflict == 0:
                                result['conflict'] = 0
                                ps2 = __ps2

                    ref_code = getOneLetterCodeCanSequence(ps1['comp_id'])
                    test_code = getOneLetterCodeCanSequence(ps2['comp_id'])

                    for r_code, t_code, seq_id, seq_id2 in zip(ref_code, test_code, ps1['seq_id'], ps2['seq_id']):
                        if r_code == 'X' and t_code == 'X':
                            nmr_input_source.updateNonStandardResidueByExptlData(chain_id2, seq_id2, 'coordinate')
                            cif_input_source.updateNonStandardResidueByExptlData(chain_id, seq_id, 'coordinate')

                    if length == _matched == 1:  # to avoid redundant ligand mapping
                        if chain_id in lig_chain_ids or chain_id2 in lig_chain_ids2:
                            continue
                        lig_chain_ids.append(chain_id)
                        lig_chain_ids2.append(chain_id2)

                    if result['unmapped'] > 0 or result['conflict'] > 0:

                        aligned = [True] * length
                        seq_id1, seq_id2 = [], []

                        j = 0
                        for i in range(length):
                            if j < len(ps1['seq_id']) and (str(myAlign[i][0]) != '.' or ps1['comp_id'][j] == '.'):
                                seq_id1.append(ps1['seq_id'][j])
                                j += 1
                            else:
                                seq_id1.append(None)

                        j = 0
                        for i in range(length):
                            if j < len(ps2['seq_id']) and (str(myAlign[i][1]) != '.' or ps2['comp_id'][j] == '.'):
                                seq_id2.append(ps2['seq_id'][j])
                                j += 1
                            else:
                                seq_id2.append(None)

                        for i in range(length):
                            myPr0, myPr1 = str(myAlign[i][0]), str(myAlign[i][1])

                            if myPr0 == '.' or myPr1 == '.':
                                aligned[i] = False
                            elif myPr0 != myPr1:
                                pass
                            else:
                                break

                        for i in reversed(range(length)):
                            myPr0, myPr1 = str(myAlign[i][0]), str(myAlign[i][1])

                            if myPr0 == '.' or myPr1 == '.':
                                aligned[i] = False
                            elif myPr0 != myPr1:
                                pass
                            else:
                                break

                        if not self.__reg.native_combined:

                            _conflicts = 0

                            for i in range(length):
                                cif_comp_id, nmr_comp_id = str(myAlign[i][0]), str(myAlign[i][1])

                                if cif_comp_id == nmr_comp_id:
                                    continue

                                if nmr_comp_id == '.' and cif_comp_id != '.':
                                    pass

                                elif nmr_comp_id != cif_comp_id and aligned[i]:
                                    _conflicts += 1

                            if _conflicts > ca['unmapped'] and ca['sequence_coverage'] < MIN_SEQ_COVERAGE_W_CONFLICT:
                                continue

                            # DAOTHER-7825 (2lyw)
                            if _conflicts + offset_1 > _matched and ca['sequence_coverage'] < LOW_SEQ_COVERAGE:
                                if not low_evid_chain_mapping:  # DAOTHER-8751
                                    continue

                        unmapped, conflict = [], []

                        for i in range(length):
                            cif_comp_id, nmr_comp_id = str(myAlign[i][0]), str(myAlign[i][1])

                            if cif_comp_id == nmr_comp_id:
                                continue

                            if nmr_comp_id == '.' and cif_comp_id != '.':

                                unmapped.append({'ref_seq_id': seq_id1[i], 'ref_comp_id': cif_comp_id})

                                try:
                                    label_seq_id = seq_id1[i]
                                    auth_seq_id = ps1['auth_seq_id'][ps1['seq_id'].index(label_seq_id)]
                                except (KeyError, IndexError, ValueError):
                                    label_seq_id = seq_id1[i]
                                    auth_seq_id = label_seq_id

                                if not aligned[i]:
                                    cif_seq_code = f"{chain_id}:{label_seq_id}:{cif_comp_id}"
                                    if chain_id != auth_chain_id or label_seq_id != auth_seq_id:
                                        cif_seq_code += f" ({auth_chain_id}:{auth_seq_id}:{cif_comp_id} in author sequence scheme)"

                                    warn = f"{cif_seq_code} is not present in the NMR data (chain_id {chain_id2})."

                                    self.__reg.suspended_warnings_for_lazy_eval.append({'sequence_mismatch':
                                                                                        {'ca_idx': ca_idx,
                                                                                         'file_name': cif_file_name,
                                                                                         'description': warn}})

                                    if self.__reg.verbose:
                                        self.__reg.log.write(f"+{self.__class_name__}.__assignCoordPolymerSequence() "
                                                             f"++ Warning  - {warn}\n")

                            elif nmr_comp_id != cif_comp_id and aligned[i]:

                                _seq_id1 = seq_id1[i]
                                _seq_id2 = seq_id2[i]

                                if _seq_id1 is None and _seq_id2 is None:
                                    continue

                                if _seq_id1 == _seq_id2 and len(nmr_comp_id) + 1 == len(cif_comp_id)\
                                   and cif_comp_id.startswith(nmr_comp_id) and cif_comp_id[-1] in ('5', '3'):
                                    continue

                                conflict.append({'ref_seq_id': _seq_id1, 'ref_comp_id': cif_comp_id,
                                                 'test_seq_id': _seq_id2, 'test_comp_id': nmr_comp_id})

                                try:
                                    label_seq_id = _seq_id1
                                    auth_seq_id = ps1['auth_seq_id'][ps1['seq_id'].index(label_seq_id)]
                                except (KeyError, IndexError, ValueError):
                                    label_seq_id = _seq_id1
                                    auth_seq_id = label_seq_id

                                cif_seq_code = f"{chain_id}:{label_seq_id}:{cif_comp_id}"
                                if cif_comp_id == '.':
                                    cif_seq_code += ', insertion error'
                                    # DAOTHER-9644: skip insertion error due to truncated loop
                                    cif_ps = next(cif_ps for cif_ps in cif_poly_seq if cif_ps['chain_id'] == chain_id)
                                    if label_seq_id is None:
                                        if self.__reg.caC is not None and 'missing_polymer_linkage' in self.__reg.caC\
                                           and any(True for mis in self.__reg.caC['missing_polymer_linkage']
                                                   if mis['auth_chain_id'] == (cif_ps['auth_chain_id'] if 'auth_chain_id' in cif_ps
                                                                               else cif_ps['chain_id'])):
                                            continue
                                    else:
                                        if self.__reg.caC is not None and 'missing_polymer_linkage' in self.__reg.caC\
                                           and any(True for mis in self.__reg.caC['missing_polymer_linkage']
                                                   if mis['auth_chain_id'] == (cif_ps['auth_chain_id'] if 'auth_chain_id' in cif_ps
                                                                               else cif_ps['chain_id'])
                                                   and mis['auth_seq_id_1'] < auth_seq_id < mis['auth_seq_id_2']):
                                            continue
                                nmr_seq_code = f"{chain_id2}:{_seq_id2}:{nmr_comp_id}"
                                if nmr_comp_id == '.':
                                    nmr_seq_code += ', insertion error'
                                if cif_comp_id != '.':
                                    if chain_id != auth_chain_id or label_seq_id != auth_seq_id:
                                        cif_seq_code += f", or {auth_chain_id}:{auth_seq_id}:{cif_comp_id} "\
                                            "in author sequence scheme"

                                err = f"Sequence alignment error between the coordinate ({cif_seq_code}) "\
                                    f"and the NMR data ({nmr_seq_code}). "\
                                    "Please verify the two sequences and re-upload the correct file(s)."

                                if self.__reg.tolerant_seq_align and self.__reg.dpV.equalsToRepCompId(nmr_comp_id, cif_comp_id):
                                    self.__reg.suspended_warnings_for_lazy_eval.append({'sequence_mismatch':
                                                                                        {'ca_idx': ca_idx,
                                                                                         'file_name': cif_file_name,
                                                                                         'description': err}})

                                    if self.__reg.verbose:
                                        self.__reg.log.write(f"+{self.__class_name__}.__assignCoordPolymerSequence() "
                                                             f"++ Warning  - {err}\n")

                                elif not self.__reg.annotation_mode:
                                    self.__reg.suspended_errors_for_lazy_eval.append({'sequence_mismatch':
                                                                                      {'ca_idx': ca_idx,
                                                                                       'file_name': cif_file_name,
                                                                                       'description': err}})

                                    if self.__reg.verbose:
                                        self.__reg.log.write(f"+{self.__class_name__}.__assignCoordPolymerSequence() "
                                                             f"++ Error  - {err}\n")

                        if len(unmapped) > 0:
                            ca['unmapped_sequence'] = unmapped

                        if len(conflict) > 0:
                            ca['conflict_sequence'] = conflict
                            ca['conflict'] = len(conflict)
                            ca['unmapped'] = ca['unmapped'] - len(conflict)
                            if ca['unmapped'] < 0:
                                ca['conflict'] -= ca['unmapped']
                                ca['unmapped'] = 0

                            result['conflict'] = ca['conflict']
                            result['unmapped'] = ca['unmapped']

                            if _result is not None:
                                _result['conflict'] = ca['conflict']
                                _result['unmapped'] = ca['unmapped']

                    chain_assign.append(ca)
                    ca_idx += 1

                if len(chain_assign) > 0 and fileListId == 0:

                    if len(cif_poly_seq) > 1:

                        if len(self.__reg.suspended_errors_for_lazy_eval) + len(self.__reg.suspended_warnings_for_lazy_eval) > 0:

                            _del_ca_idx = []

                            for ca_idx, ca in enumerate(chain_assign):

                                if ca['conflict'] == 0:
                                    continue

                                ref_chain_id = ca['ref_chain_id']
                                test_chain_id = ca['test_chain_id']

                                if any(True for _ca in chain_assign
                                       if ((_ca['ref_chain_id'] == ref_chain_id and _ca['test_chain_id'] != test_chain_id)
                                           or (_ca['ref_chain_id'] != ref_chain_id and _ca['test_chain_id'] == test_chain_id))
                                       and _ca['conflict'] == 0):
                                    _del_ca_idx.append(ca_idx)

                            if len(_del_ca_idx) > 0:
                                for ca_idx in reversed(_del_ca_idx):
                                    del chain_assign[ca_idx]
                                if len(self.__reg.suspended_errors_for_lazy_eval) > 0:
                                    _del_msg_idx = set()
                                    for msg_idx, msg in enumerate(self.__reg.suspended_errors_for_lazy_eval):
                                        for k, v in msg.items():
                                            if v['ca_idx'] in _del_ca_idx:
                                                _del_msg_idx.add(msg_idx)
                                    if len(_del_msg_idx) > 0:
                                        for msg_idx in reversed(list(_del_msg_idx)):
                                            if msg_idx < len(self.__reg.suspended_errors_for_lazy_eval):
                                                del self.__reg.suspended_errors_for_lazy_eval[msg_idx]
                                if len(self.__reg.suspended_warnings_for_lazy_eval) > 0:
                                    _del_msg_idx = set()
                                    for msg_idx, msg in enumerate(self.__reg.suspended_warnings_for_lazy_eval):
                                        for k, v in msg.items():
                                            if v['ca_idx'] in _del_ca_idx:
                                                _del_msg_idx.add(msg_idx)
                                    if len(_del_msg_idx) > 0:
                                        for msg_idx in reversed(list(_del_msg_idx)):
                                            if msg_idx < len(self.__reg.suspended_warnings_for_lazy_eval):
                                                del self.__reg.suspended_warnings_for_lazy_eval[msg_idx]

                        if any(True for ps in cif_poly_seq if 'identical_chain_id' in ps):

                            _chain_assign = chain_assign.copy()

                            for ca in _chain_assign:

                                if ca['conflict'] > 0:
                                    continue

                                chain_id = ca['ref_chain_id']
                                auth_chain_id = ca.get('ref_auth_chain_id')

                                if chain_id in lig_chain_ids:
                                    continue

                                try:
                                    identity = next(ps['identical_chain_id'] for ps in cif_poly_seq
                                                    if ps['chain_id'] == chain_id and 'identical_chain_id' in ps)

                                    for chain_id in identity:

                                        if not any(True for _ca in chain_assign if _ca['ref_chain_id'] == chain_id):
                                            _ca = ca.copy()
                                            _ca['ref_chain_id'] = chain_id
                                            if auth_chain_id is not None:
                                                _ca['ref_auth_chain_id'] = auth_chain_id

                                            if _ca['test_chain_id'] in lig_chain_ids2:
                                                continue

                                            chain_assign.append(_ca)

                                except StopIteration:
                                    pass

                            if self.__reg.combined_mode:  # DAOTHER-9660, issue #2
                                nmr_ident_chain_id = set()
                                cif_ident_chain_id = set()

                                for ca in chain_assign:

                                    if ca['conflict'] > 0:
                                        continue

                                    nmr_chain_id = ca['test_chain_id']
                                    cif_chain_id = ca['ref_chain_id']

                                    try:
                                        identity = next(ps['identical_chain_id'] for ps in nmr_poly_seq
                                                        if ps['chain_id'] == nmr_chain_id and 'identical_chain_id' in ps)
                                        nmr_ident_chain_id.add(nmr_chain_id)
                                    except StopIteration:
                                        pass

                                    try:
                                        identity = next(ps['identical_chain_id'] for ps in cif_poly_seq
                                                        if ps['chain_id'] == cif_chain_id and 'identical_chain_id' in ps)
                                        cif_ident_chain_id.add(cif_chain_id)
                                    except StopIteration:
                                        pass

                                if len(cif_ident_chain_id) == len(nmr_ident_chain_id) and len(cif_ident_chain_id) > 1:
                                    _chain_assign = chain_assign.copy()

                                    nmr_mapped_chain_id, cif_mapped_chain_id = [], []

                                    for ca in _chain_assign:

                                        if ca['conflict'] > 0:
                                            continue

                                        nmr_chain_id = ca['test_chain_id']
                                        cif_chain_id = ca['ref_chain_id']

                                        if nmr_chain_id not in nmr_ident_chain_id or cif_chain_id not in cif_ident_chain_id:
                                            continue

                                        if nmr_chain_id not in nmr_mapped_chain_id and cif_chain_id not in cif_mapped_chain_id:
                                            nmr_mapped_chain_id.append(nmr_chain_id)
                                            cif_mapped_chain_id.append(cif_chain_id)

                                        else:
                                            chain_assign.remove(ca)

                    self.__reg.report.chain_assignment.setItemValue('model_poly_seq_vs_nmr_poly_seq', chain_assign)

                    if len(self.__reg.suspended_errors_for_lazy_eval) > 0:
                        for msg in self.__reg.suspended_errors_for_lazy_eval:
                            for k, v in msg.items():
                                if 'ca_idx' in v:
                                    del v['ca_idx']
                                self.__reg.report.error.appendDescription(k, v)
                        self.__reg.suspended_errors_for_lazy_eval.clear()

                    if len(self.__reg.suspended_warnings_for_lazy_eval) > 0:
                        for msg in self.__reg.suspended_warnings_for_lazy_eval:
                            for k, v in msg.items():
                                if 'ca_idx' in v:
                                    del v['ca_idx']
                                self.__reg.report.warning.appendDescription(k, v)
                        self.__reg.suspended_warnings_for_lazy_eval.clear()

                chain_assign_dic = self.__reg.report.chain_assignment.get()

                if has_key_value(chain_assign_dic, 'nmr_poly_seq_vs_model_poly_seq'):
                    ref_chain_ids = {}
                    for ca in chain_assign_dic['nmr_poly_seq_vs_model_poly_seq']:
                        if ca['ref_chain_id'] in ref_chain_ids or (ca['length'] == 1 and 1 in ref_chain_ids.values()):
                            continue
                        self.__reg.label_asym_id_with_exptl_data.add(ca['test_chain_id'])
                        ref_chain_ids[ca['ref_chain_id']] = ca['length']

            elif not self.__reg.annotation_mode:

                poly_seq = input_source_dic['polymer_sequence']

                single_ligand = poly_seq is not None and len(poly_seq) == 1 and len(poly_seq[0]['seq_id']) == 1\
                    and poly_seq[0]['comp_id'][0] not in STD_MON_DICT

                if single_ligand:

                    err = "Failed in sequence alignment because of single ligand."

                    self.__reg.report.warning.appendDescription('sequence_mismatch',
                                                                {'file_name': cif_file_name, 'description': err})

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__assignCoordPolymerSequence() "
                                             f"++ Warning  - {err}\n")

                else:

                    err = "No polymer sequence alignments found."

                    self.__reg.report.error.appendDescription('sequence_mismatch',
                                                              {'file_name': cif_file_name, 'description': err})

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__assignCoordPolymerSequence() "
                                             f"++ Error  - {err}\n")

                return False

        return self.__reg.report.getTotalErrors() == __errors

    def __testCoordAtomIdConsistency(self) -> bool:
        """ Perform consistency test on atom names of coordinate file.
        """

        src_id = self.__reg.report.getInputSourceIdOfCoord()

        if src_id < 0:
            return False

        cif_input_source = self.__reg.report.input_sources[src_id]
        cif_input_source_dic = cif_input_source.get()

        cif_poly_seq = cif_input_source_dic['polymer_sequence']

        __errors = self.__reg.report.getTotalErrors()

        for fileListId in range(self.__reg.file_path_list_len):

            nmr_input_source = self.__reg.report.input_sources[fileListId]
            nmr_input_source_dic = nmr_input_source.get()

            file_name = nmr_input_source_dic['file_name']
            file_type = nmr_input_source_dic['file_type']

            seq_align_dic = self.__reg.report.sequence_alignment.get()
            chain_assign_dic = self.__reg.report.chain_assignment.get()

            if 'nmr_poly_seq_vs_model_poly_seq' not in chain_assign_dic:

                err = "Chain assignment does not exist, __assignCoordPolymerSequence() should be invoked."

                self.__reg.report.error.appendDescription('internal_error',
                                                          f"+{self.__class_name__}.__testCoordAtomIdConsistency() "
                                                          "++ Error  - " + err)

                if self.__reg.verbose:
                    self.__reg.log.write(f"+{self.__class_name__}.__testCoordAtomIdConsistency() ++ Error  - {err}\n")

                continue

            if not has_key_value(seq_align_dic, 'nmr_poly_seq_vs_model_poly_seq'):
                continue

            if not has_key_value(chain_assign_dic, 'nmr_poly_seq_vs_model_poly_seq'):
                continue

            nmr2ca = {}

            _ref_chain_id = None

            for ca in chain_assign_dic['nmr_poly_seq_vs_model_poly_seq']:

                ref_chain_id = ca['ref_chain_id']
                test_chain_id = ca['test_chain_id']

                if _ref_chain_id is None:
                    _ref_chain_id = ref_chain_id

                result = next((seq_align for seq_align in seq_align_dic['nmr_poly_seq_vs_model_poly_seq']
                               if seq_align['ref_chain_id'] == ref_chain_id and seq_align['test_chain_id'] == test_chain_id), None)

                if ref_chain_id not in nmr2ca:
                    nmr2ca[ref_chain_id] = []

                sa = {'seq_align': result}  # DAOTHER-7465

                if 'unmapped_sequence' in ca:
                    sa['seq_unmap'] = [unmapped['ref_seq_id'] for unmapped in ca['unmapped_sequence']]

                nmr2ca[ref_chain_id].append(sa)

            if nmr_input_source_dic['content_subtype'] is None:
                continue

            ref_chain_id = _ref_chain_id

            modified = False

            for content_subtype in nmr_input_source_dic['content_subtype']:

                if content_subtype in ('entry_info', 'entity'):
                    continue

                sf_category = SF_CATEGORIES[file_type][content_subtype]
                lp_category = LP_CATEGORIES[file_type][content_subtype]

                if content_subtype == 'poly_seq':
                    lp_category = AUX_LP_CATEGORIES[file_type][content_subtype][0]

                list_id = 1

                if self.__reg.star_data_type[fileListId] == 'Loop':
                    sf = self.__reg.star_data[fileListId]
                    sf_framecode = ''

                    modified |= self.__reg.dpV.testCoordAtomIdConsistency(fileListId, file_name, file_type, content_subtype, sf,
                                                                          list_id, sf_framecode, lp_category, cif_poly_seq,
                                                                          seq_align_dic, nmr2ca, ref_chain_id)

                elif self.__reg.star_data_type[fileListId] == 'Saveframe':
                    sf = self.__reg.star_data[fileListId]
                    sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                    modified |= self.__reg.dpV.testCoordAtomIdConsistency(fileListId, file_name, file_type, content_subtype, sf,
                                                                          list_id, sf_framecode, lp_category, cif_poly_seq,
                                                                          seq_align_dic, nmr2ca, ref_chain_id)

                else:

                    for list_id, sf in enumerate(self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category), start=1):
                        sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                        if not any(True for loop in sf.loops if loop.category == lp_category):
                            continue

                        modified |= self.__reg.dpV.testCoordAtomIdConsistency(fileListId, file_name, file_type, content_subtype,
                                                                              sf, list_id, sf_framecode, lp_category, cif_poly_seq,
                                                                              seq_align_dic, nmr2ca, ref_chain_id)

            if modified and fileListId == 0:
                self.__depositNmrData()

        return self.__reg.report.getTotalErrors() == __errors

    def __retrieveDpReport(self) -> bool:
        """ Retrieve NMR data processing report from JSON file.
        """

        if not self.__reg.combined_mode:
            return True

        # retrieve sf_category_list which is required to resolve minor issues
        if len(self.__reg.sf_category_list) == 0:

            _, star_data_type, star_data = self.__reg.nefT.read_input_file(self.__reg.srcPath)

            self.__reg.sf_category_list, self.__reg.lp_category_list = self.__reg.nefT.get_inventory_list(star_data)

            if len(self.__reg.star_data_type) == 0:
                self.__reg.star_data_type.append(star_data_type)
            else:
                self.__reg.star_data_type[0] = star_data_type

            if len(self.__reg.star_data) == 0:
                self.__reg.star_data.append(star_data)
            else:
                self.__reg.star_data[0] = star_data

            corrections = self.__reg.nefT.resolve_sf_names_for_cif(star_data)  # DAOTHER-7389, issue #4

            if len(self.__reg.sf_name_corrections) == 0:
                self.__reg.sf_name_corrections.append(corrections)
            else:
                self.__reg.sf_name_corrections[0] = corrections

        if REPORT_FILE_PATH_KEY not in self.__reg.inputParamDict or self.__reg.annotation_mode:
            self.__initializeDpReport()
            self.__reg.dstPath = self.__dstPath__ if self.__reg.cifPath is None else self.__reg.srcPath

            return False

        fPath = self.__reg.inputParamDict[REPORT_FILE_PATH_KEY]

        if not os.access(fPath, os.F_OK):
            raise IOError(f"+{self.__class_name__}.__retrieveDpReport() "
                          f"++ Error  - Could not access to file path {fPath}.")

        if os.path.getsize(fPath) == 0:
            raise IOError(f"+{self.__class_name__}.__retrieveDpReport() "
                          f"++ Error  - Could not find any content in file path {fPath}.")

        self.__reg.report = NmrDpReport(self.__reg.verbose, self.__reg.log)
        self.__reg.report.loadFile(fPath)

        self.__reg.report_prev = NmrDpReport(self.__reg.verbose, self.__reg.log)
        self.__reg.report_prev.loadFile(fPath)

        return True

    def __resolveConflictsInLoop(self) -> bool:
        """ Resolve conflicted rows in loops.
        """

        if not self.__reg.combined_mode:
            return True

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_name = input_source_dic['file_name']
        file_type = input_source_dic['file_type']

        if input_source_dic['content_subtype'] is None:
            return False

        if not self.__reg.resolve_conflict:
            return True

        for content_subtype in input_source_dic['content_subtype']:

            if content_subtype in ('entry_info', 'entity'):
                continue

            sf_category = SF_CATEGORIES[file_type][content_subtype]
            lp_category = LP_CATEGORIES[file_type][content_subtype]

            for sf in self.__reg.star_data[0].get_saveframes_by_category(sf_category):

                if not any(True for loop in sf.loops if loop.category == lp_category):
                    continue

                sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                if content_subtype == 'spectral_peak':

                    try:

                        _num_dim = get_first_sf_tag(sf, NUM_DIM_ITEMS[file_type])
                        num_dim = int(_num_dim)

                        if num_dim not in range(1, MAX_DIM_NUM_OF_SPECTRA):
                            raise ValueError()

                    except ValueError:  # raised error already at __testIndexConsistency()
                        continue

                    max_dim = num_dim + 1

                    key_items = []
                    for dim in range(1, max_dim):
                        for k in PK_KEY_ITEMS[file_type]:
                            if k['type'] == 'float':  # position
                                _k = copy.copy(k)
                                if '%s' in k['name']:
                                    _k['name'] = k['name'] % dim
                                key_items.append(_k)
                    for k in PK_KEY_ITEMS[file_type]:
                        if k['type'] == 'positive-int':  # peak_id
                            key_items.append(k)

                    data_items = []
                    for d in DATA_ITEMS[file_type][content_subtype]:
                        data_items.append(d)
                    for dim in range(1, max_dim):
                        for d in self.__reg.pk_data_items[file_type]:
                            _d = copy.copy(d)
                            if '%s' in d['name']:
                                _d['name'] = d['name'] % dim
                            if 'default-from' in d and '%s' in d['default-from']:  # DAOTHER-7421
                                _d['default-from'] = d['default-from'] % dim
                            data_items.append(_d)

                else:

                    key_items = self.__reg.key_items[file_type][content_subtype]
                    data_items = DATA_ITEMS[file_type][content_subtype]

                    if len(key_items) == 0:
                        continue

                modified = False

                if content_subtype in ('dist_restraint', 'dihed_restraint', 'rdc_restraint'):

                    try:

                        conflict_id = self.__reg.nefT.get_conflict_atom_id(sf, file_type, lp_category, key_items)[0]

                        if len(conflict_id) > 0:
                            loop = sf.get_loop(lp_category)

                            for _id in conflict_id:
                                del loop.data[_id]

                            modified = True

                        conflict_id = self.__reg.nefT.get_bad_pattern_id(sf, lp_category, key_items, data_items)[0]

                        if len(conflict_id) > 0:
                            loop = sf.get_loop(lp_category)

                            for _id in conflict_id:
                                del loop.data[_id]

                            modified = True

                        if modified:

                            lp = next((lp for lp in self.__reg.lp_data[content_subtype]
                                       if lp['file_name'] == file_name and lp['sf_framecode'] == sf_framecode), None)

                            if lp is not None:
                                lp['data'] = self.__reg.nefT.check_data(sf, lp_category, key_items, data_items, None, None, None,
                                                                        enforce_allowed_tags=(file_type == 'nmr-star'),
                                                                        excl_missing_data=self.__reg.excl_missing_data)[0]

                    except Exception:  # pylint: disable=broad-exception-caught
                        pass

        return True

    def __resolveConflictsInAuxLoop(self) -> bool:
        """ Resolve conflicted rows in auxiliary loops.
        """

        if not self.__reg.combined_mode:
            return True

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_type = input_source_dic['file_type']

        if input_source_dic['content_subtype'] is None:
            return False

        if not self.__reg.resolve_conflict:
            return True

        for content_subtype in input_source_dic['content_subtype']:

            if content_subtype == 'entity':
                continue

            sf_category = SF_CATEGORIES[file_type][content_subtype]
            lp_category = LP_CATEGORIES[file_type][content_subtype]

            for sf in self.__reg.star_data[0].get_saveframes_by_category(sf_category):

                if content_subtype.startswith('spectral_peak'):

                    try:

                        _num_dim = get_first_sf_tag(sf, NUM_DIM_ITEMS[file_type])
                        num_dim = int(_num_dim)

                        if num_dim not in range(1, MAX_DIM_NUM_OF_SPECTRA):
                            raise ValueError()

                    except ValueError:  # raised error already at __testIndexConsistency()
                        pass

                for loop in sf.loops:

                    lp_category = loop.category

                    if lp_category is None:
                        continue

                    # main content of loop has been processed in testDataConsistencyInLoop()
                    if lp_category in LP_CATEGORIES[file_type][content_subtype]:
                        continue

                    if AUX_LP_CATEGORIES[file_type][content_subtype] is None:
                        continue

                    if lp_category in AUX_LP_CATEGORIES[file_type][content_subtype]:

                        key_items = self.__reg.aux_key_items[file_type][content_subtype][lp_category]

                        if len(key_items) == 0:
                            continue

                        try:

                            conflict_id = self.__reg.nefT.get_conflict_id(sf, lp_category, key_items)[0]

                            if len(conflict_id) > 0:
                                _loop = sf.get_loop(lp_category)

                                for _id in conflict_id:
                                    del _loop.data[_id]

                        except Exception:  # pylint: disable=broad-exception-caught
                            pass

        return True

    def __appendIndexTag(self) -> bool:
        """ Append index tag if required.
        """

        if not self.__reg.combined_mode:
            return True

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            if input_source_dic['content_subtype'] is None:
                continue

            for content_subtype in input_source_dic['content_subtype']:

                if content_subtype == 'entity':
                    continue

                sf_category = SF_CATEGORIES[file_type][content_subtype]
                lp_category = LP_CATEGORIES[file_type][content_subtype]

                index_tag = INDEX_TAGS[file_type][content_subtype]

                if index_tag is None:
                    continue

                for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                    sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                    try:
                        loop = sf.get_loop(lp_category)
                    except KeyError:
                        continue

                    if index_tag in loop.tags:
                        continue

                    lp_tag = f'{lp_category}.{index_tag}'
                    err = ERR_TEMPLATE_FOR_MISSING_MANDATORY_LP_TAG % (lp_tag, file_type.upper())

                    if self.__reg.check_mandatory_tag and self.__reg.nefT.is_mandatory_tag(lp_tag, file_type):

                        if self.__reg.rescue_mode:
                            self.__reg.report.error.appendDescription('missing_mandatory_item',
                                                                      {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                       'category': lp_category, 'description': err})

                            if self.__reg.verbose:
                                self.__reg.log.write(f"+{self.__class_name__}.__appendIndexTag() "
                                                     f"++ LookupError  - {file_name} {sf_framecode} {lp_category} {err}\n")

                    lp = pynmrstar.Loop.from_scratch(lp_category)

                    lp.add_tag(lp_tag)
                    lp.add_tag(loop.tags)

                    for idx, row in enumerate(loop, start=1):
                        lp.add_data([idx] + row)

                    del sf[loop]

                    sf.add_loop(lp)

        return True

    def __deleteSkippedSf(self) -> bool:
        """ Delete skipped saveframes.
        """

        if not self.__reg.combined_mode or self.__reg.submission_mode or self.__reg.annotation_mode or self.__reg.release_mode:
            return True

        if len(self.__reg.star_data) == 0:
            return False

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_name = input_source_dic['file_name']

        warnings = self.__reg.report.warning.getValueList('skipped_saveframe_category', file_name)

        if warnings is None:
            return True

        if self.__reg.retain_original:
            return True

        for w in warnings:

            if self.__reg.star_data_type[0] == 'Entry' or self.__reg.star_data_type[0] == 'Saveframe':

                if 'sf_category' not in w:

                    err = "Could not specify 'sf_category' in NMR data processing report."

                    self.__reg.report.error.appendDescription('internal_error',
                                                              f"+{self.__class_name__}.__deleteSkippedSf() ++ Error  - " + err)

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__deleteSkippedSf() ++ Error  - {err}\n")

                else:

                    sf_list = self.__reg.star_data[0].get_saveframes_by_category(w['sf_category'])

                    if sf_list is None:

                        err = f"Could not specify sf_category {w['sf_category']} unexpectedly."

                        self.__reg.report.error.appendDescription('internal_error',
                                                                  f"+{self.__class_name__}.__deleteSkippedSf() ++ Error  - " + err)

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__deleteSkippedSf() ++ Error  - {err}\n")

                    else:

                        for sf in reversed(sf_list):
                            del self.__reg.star_data[0][sf]

            else:

                err = f"Unexpected PyNMRSTAR object type {self.__reg.star_data_type[0]} found about {file_name!r} file."

                self.__reg.report.error.appendDescription('internal_error',
                                                          f"+{self.__class_name__}.__deleteSkippedSf() ++ Error  - " + err)

                if self.__reg.verbose:
                    self.__reg.log.write(f"+{self.__class_name__}.__deleteSkippedSf() ++ Error  - {err}\n")

        return True

    def __deleteSkippedLoop(self) -> bool:
        """ Delete skipped loops.
        """

        if not self.__reg.combined_mode:
            return True

        if len(self.__reg.star_data) == 0:
            return False

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_name = input_source_dic['file_name']

        warnings = self.__reg.report.warning.getValueList('skipped_loop_category', file_name)

        if warnings is None:
            return True

        if self.__reg.retain_original:
            return True

        for w in warnings:

            if self.__reg.star_data_type[0] == 'Entry' or self.__reg.star_data_type[0] == 'Saveframe':

                if 'sf_framecode' not in w:

                    err = "Could not specify 'sf_framecode' in NMR data processing report."

                    self.__reg.report.error.appendDescription('internal_error',
                                                              f"+{self.__class_name__}.__deleteSkippedLoop() ++ Error  - " + err)

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__deleteSkippedLoop() ++ Error  - {err}\n")

                else:

                    sf = self.__reg.dpA.getSaveframeByName(0, w['sf_framecode'])

                    if sf is None:

                        err = f"Could not specify {w['sf_framecode']!r} saveframe unexpectedly in {file_name!r} file."

                        self.__reg.report.error.appendDescription('internal_error',
                                                                  f"+{self.__class_name__}.__deleteSkippedLoop() "
                                                                  "++ Error  - " + err)

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__deleteSkippedLoop() ++ Error  - {err}\n")

                        continue

                    if 'category' not in w:

                        err = "Could not specify 'category' in NMR data processing report."

                        self.__reg.report.error.appendDescription('internal_error',
                                                                  f"+{self.__class_name__}.__deleteSkippedLoop() "
                                                                  "++ Error  - " + err)

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__deleteSkippedLoop() ++ Error  - {err}\n")

                    else:

                        del sf[sf.get_loop(w['category'])]

            else:

                err = f"Unexpected PyNMRSTAR object type {self.__reg.star_data_type[0]} found about {file_name!r} file."

                self.__reg.report.error.appendDescription('internal_error',
                                                          f"+{self.__class_name__}.__deleteSkippedLoop() ++ Error  - " + err)

                if self.__reg.verbose:
                    self.__reg.log.write(f"+{self.__class_name__}.__deleteSkippedLoop() ++ Error  - {err}\n")

        return True

    def __deleteUnparsedEntryLoop(self) -> bool:
        """ Delete unparsed entry loops.
        """

        if not self.__reg.combined_mode:
            return True

        if len(self.__reg.star_data) == 0:
            return False

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_name = input_source_dic['file_name']
        file_type = input_source_dic['file_type']

        if self.__reg.retain_original:
            return True

        content_subtype = 'entry_info'

        sf_category = SF_CATEGORIES[file_type][content_subtype]
        lp_category = LP_CATEGORIES[file_type][content_subtype]

        if self.__reg.star_data_type[0] == 'Entry' or self.__reg.star_data_type[0] == 'Saveframe':

            if sf_category in self.__reg.sf_category_list:

                for sf in self.__reg.star_data[0].get_saveframes_by_category(sf_category):

                    loops = []

                    for loop in sf.loops:

                        if loop.category == lp_category:
                            continue

                        loops.append(loop)

                    for loop in reversed(loops):
                        del sf[loop]

        else:

            err = f"Unexpected PyNMRSTAR object type {self.__reg.star_data_type[0]} found about {file_name!r} file."

            self.__reg.report.error.appendDescription('internal_error',
                                                      f"+{self.__class_name__}.__deleteUnparsedEntryLoop() ++ Error  - " + err)

            if self.__reg.verbose:
                self.__reg.log.write(f"+{self.__class_name__}.__deleteUnparsedEntryLoop() ++ Error  - {err}\n")

        return True

    def __updatePolymerSequence(self) -> bool:
        """ Update polymer sequence.
        """

        if len(self.__reg.star_data) == 0 or self.__reg.star_data[0] is None or self.__reg.star_data_type[0] != 'Entry':
            return False

        self.__reg.nefT.set_merge_rescue_mode(False)  # D_1300021766 vs DAOTHER-9927

        self.__extractPolymerSequence()
        self.__extractPolymerSequenceInLoop()

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_type = input_source_dic['file_type']
        orig_poly_seq = input_source_dic['polymer_sequence']

        has_poly_seq = has_key_value(input_source_dic, 'polymer_sequence')

        if not has_poly_seq:
            return False

        if not self.__reg.submission_mode and not self.__reg.annotation_mode and not self.__reg.release_mode:
            self.__extractToNmrIf__()
            self.__reg.dpR.cleanUpSf()

        if not self.__reg.cifChecked:
            return self.__remediateCsLoop__()

        if self.__reg.caC is None:
            self.__retrieveCoordAssemblyChecker__()

        is_done, asm_sf = self.__reg.dpR.updatePolymerSequence()

        if not is_done:
            return False

        content_subtype = 'poly_seq'

        try:
            poly_seq = self.__getPolymerSequence__(0, asm_sf, content_subtype)[0]
        except KeyError:
            return False
        except UserWarning:
            poly_seq = []

        identical = True

        # to process large assembly avoiding forced timeout
        if len(poly_seq) < LEN_MAJOR_ASYM_ID or len(poly_seq) != len(orig_poly_seq):
            seq_align, _ = alignPolymerSequence(self.__reg.pA, poly_seq, orig_poly_seq, conservative=False)
            chain_assign, _ = assignPolymerSequence(self.__reg.pA, self.__reg.ccU, file_type, poly_seq, orig_poly_seq, seq_align)

            # trim cross chain assignment (DAOTHER-10898)
            if len(chain_assign) > len(poly_seq) > 1:

                _del_ca_idx = []

                for ca_idx, ca in enumerate(chain_assign):

                    if ca['conflict'] != 0:
                        continue

                    ref_chain_id = ca['ref_chain_id']
                    test_chain_id = ca['test_chain_id']

                    if any(True for _ca in chain_assign
                           if ((_ca['ref_chain_id'] == ref_chain_id and _ca['test_chain_id'] != test_chain_id)
                               or (_ca['ref_chain_id'] != ref_chain_id and _ca['test_chain_id'] == test_chain_id))
                           and _ca['conflict'] == 0):
                        _del_ca_idx.append(ca_idx)

                if len(_del_ca_idx) > 0:
                    for ca_idx in reversed(_del_ca_idx):
                        del chain_assign[ca_idx]

            self.__reg.chain_id_map_for_remediation.clear()
            self.__reg.seq_id_map_for_remediation.clear()

            for ps in orig_poly_seq:
                chain_id = ps['chain_id']
                for seq_id, comp_id in zip(ps['seq_id'], ps['comp_id']):
                    seq_key = (chain_id, seq_id)

                    _chain_id = _seq_id = _comp_id = None

                    if chain_assign is not None:
                        _chain_id = next((ca['ref_chain_id'] for ca in chain_assign if ca['test_chain_id'] == chain_id), None)
                        if _chain_id is not None:
                            sa = next((sa for sa in seq_align
                                       if sa['ref_chain_id'] == _chain_id and sa['test_chain_id'] == chain_id
                                       and seq_id in sa['test_seq_id']), None)
                            if sa is not None:
                                _seq_id = next((ref_seq_id for ref_seq_id, test_seq_id in zip(sa['ref_seq_id'], sa['test_seq_id'])
                                                if test_seq_id == seq_id), None)
                                if _seq_id is not None:
                                    _ps = next(_ps for _ps in poly_seq if _ps['chain_id'] == _chain_id)
                                    if _seq_id in _ps['seq_id']:
                                        try:
                                            _comp_id = _ps['comp_id'][_ps['seq_id'].index(_seq_id)]
                                        except IndexError:
                                            pass

                    if _chain_id is not None and _seq_id is not None and comp_id == _comp_id:
                        if chain_id not in self.__reg.chain_id_map_for_remediation:
                            self.__reg.chain_id_map_for_remediation[chain_id] = _chain_id
                        self.__reg.seq_id_map_for_remediation[seq_key] = (_chain_id, _seq_id)
                        if chain_id != _chain_id or seq_id != _seq_id:
                            identical = False
                    else:
                        _ps = next((_ps for _ps in poly_seq if _ps['chain_id'] == _chain_id and _ps['seq_id'] == _seq_id), None)
                        if _ps is not None:
                            try:
                                _comp_id = _ps['comp_id'][_ps['seq_id'].index(_seq_id)]
                                if comp_id == _comp_id:
                                    if chain_id not in self.__reg.chain_id_map_for_remediation:
                                        self.__reg.chain_id_map_for_remediation[chain_id] = _chain_id
                                    self.__reg.seq_id_map_for_remediation[seq_key] = (_chain_id, _seq_id)
                            except IndexError:
                                pass

        self.__mergePolymerSequenceInCsLoop__(0)

        self.__remediateCsLoop__()

        if not identical:
            self.__reg.dpR.syncMrLoop()

        self.__removeUnusedPdbInsCode()

        if file_type == 'nef':
            return True

        return self.__reg.dpR.updateEntitySaveframe()

    def __remediateRdcLoop(self) -> bool:
        """ Remediate rdc target value due to the known OneDep bug, if required.
        """

        try:

            for fileListId in range(self.__reg.file_path_list_len):

                modified = False

                input_source = self.__reg.report.input_sources[fileListId]
                input_source_dic = input_source.get()

                file_type = input_source_dic['file_type']

                if file_type != 'nmr-star':
                    continue

                content_subtype = 'rdc_restraint'

                if not has_key_value(input_source_dic['content_subtype'], content_subtype):
                    continue

                sf_category = SF_CATEGORIES[file_type][content_subtype]
                lp_category = LP_CATEGORIES[file_type][content_subtype]

                if self.__reg.star_data_type[fileListId] == 'Loop':
                    if self.__reg.star_data[fileListId].category == lp_category:
                        modified |= self.__reg.dpR.remediateRdcLoop(file_type, self.__reg.star_data[fileListId])

                elif self.__reg.star_data_type[fileListId] == 'Saveframe':
                    try:
                        modified |= self.__reg.dpR.remediateRdcLoop(file_type,
                                                                    self.__reg.star_data[fileListId].get_loop(lp_category))
                    except KeyError:
                        continue

                else:
                    for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                        try:
                            modified |= self.__reg.dpR.remediateRdcLoop(file_type, sf.get_loop(lp_category))
                        except KeyError:
                            continue

                if modified:
                    self.__reg.lp_data[content_subtype].clear()
                    self.__testDataConsistencyInLoop(content_subtype)
                    self.__depositNmrData()

            return True

        except ValueError as e:

            self.__reg.report.error.appendDescription('internal_error',
                                                      f"+{self.__class_name__}.__remediateRdcLoop() ++ Error  - " + str(e))

            if self.__reg.verbose:
                self.__reg.log.write(f"+{self.__class_name__}.__remediateRdcLoop() ++ Error  - {str(e)}\n")

            return False

    def __remediateRawTextPk(self) -> bool:
        """ Remediate raw text data in saveframe of spectral peak list (for NMR data remediation upgrade to Phase 2).
        """

        if len(self.__reg.star_data) == 0 or self.__reg.star_data[0] is None or self.__reg.star_data_type[0] != 'Entry':
            return False

        # This ridiculous reverse implementation is for OneDep only
        if self.__reg.op != 'nmr-cs-mr-merge' and not self.__reg.internal_mode:

            if self.__reg.bmrb_only:
                input_source = self.__reg.report.input_sources[0]
                input_source_dic = input_source.get()

                file_type = input_source_dic['file_type']

                if file_type == 'nef':
                    return True

                master_entry = self.__reg.star_data[0]

                self.__reg.dpR.remediateSpectralPeakListSaveframe(master_entry)

                if self.__reg.srcPath is not None:
                    self.__reg.dpR.performBmrbJAnnTasks(True)

                    self.__reg.c2S.set_entry_id(master_entry, self.__reg.bmrb_id)
                    self.__reg.c2S.normalize_str(master_entry)

                    master_entry.write_to_file(self.__reg.srcPath, show_comments=True, skip_empty_loops=True, skip_empty_tags=False)

            return True

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_type = input_source_dic['file_type']

        if file_type == 'nef':
            return True

        has_poly_seq = has_key_value(input_source_dic, 'polymer_sequence')

        if not has_poly_seq:
            return False

        if self.__reg.srcPath == self.__reg.dstPath:
            return True

        content_subtype = 'spectral_peak'

        sf_category = SF_CATEGORIES[file_type][content_subtype]
        lp_category = LP_CATEGORIES[file_type][content_subtype]

        def get_reserved_list_ids():
            rlist_ids_dict = {content_subtype: []}
            for idx, sf in enumerate(self.__reg.star_data[0].get_saveframes_by_category(sf_category), start=1):
                list_id = get_first_sf_tag(sf, 'ID')
                rlist_ids_dict[content_subtype].append(int(list_id) if list_id not in EMPTY_VALUE else idx)
            return rlist_ids_dict

        reserved_list_ids = get_reserved_list_ids()
        multiple_pk_sf_idx = []

        for idx, sf in enumerate(self.__reg.star_data[0].get_saveframes_by_category(sf_category), start=1):

            if any(True for loop in sf.loops if loop.category == lp_category):
                continue

            text_data = get_first_sf_tag(sf, 'Text_data')

            if text_data in EMPTY_VALUE:
                continue

            if get_chem_shift_format_from_string(text_data) is not None:
                self.__reg.star_data[0].remove_saveframe(sf.name)
                continue

            data_format = get_first_sf_tag(sf, 'Text_data_format')

            if data_format == 'json':
                continue

            data_file_name = get_first_sf_tag(sf, 'Data_file_name')

            file_type = content_subtype_dict = None

            if data_format == 'ARIA':
                file_type = 'nm-pea-ari'
            elif data_format == 'CCPN':
                file_type = 'nm-pea-ccp'
            elif data_format == 'Olivia':
                file_type = 'nm-pea-oli'
            elif data_format == 'NMRPipe':
                file_type = 'nm-pea-pip'
            elif data_format == 'PONDEROSA':
                file_type = 'nm-pea-pon'
            elif data_format == 'Sparky':
                file_type = 'nm-pea-spa'
            elif data_format == 'TopSpin':
                file_type = 'nm-pea-top'
            elif data_format == 'NMRView':
                file_type = 'nm-pea-vie'
            elif data_format == 'VNMR':
                file_type = 'nm-pea-vnm'
            elif data_format == 'XwinNMR':
                file_type = 'nm-pea-xwi'
            else:

                for _file_type in PARSABLE_PK_FILE_TYPES:

                    if _file_type.startswith('nm-aux'):
                        continue

                    reader = self.__reg.dpS.getSimpleFileReader(_file_type, False)

                    listener, parser_err_listener, lexer_err_listener = reader.parse(text_data, None, isFilePath=False)

                    has_lexer_error = lexer_err_listener is not None and lexer_err_listener.getMessageList() is not None
                    has_parser_error = parser_err_listener is not None and parser_err_listener.getMessageList() is not None
                    content_subtype_dict = listener.getContentSubtype() if listener is not None else None
                    if not has_lexer_error and not has_parser_error and content_subtype_dict is not None\
                       and len(content_subtype_dict) > 0:
                        file_type = _file_type
                        break

            if file_type is None:
                continue

            if content_subtype_dict is None:
                reader = self.__reg.dpS.getSimpleFileReader(file_type, False)

                listener, parser_err_listener, lexer_err_listener = reader.parse(text_data, None, isFilePath=False)

                has_lexer_error = lexer_err_listener is not None and lexer_err_listener.getMessageList() is not None
                has_parser_error = parser_err_listener is not None and parser_err_listener.getMessageList() is not None
                content_subtype_dict = listener.getContentSubtype() if listener is not None else None

            if content_subtype_dict is None\
               or content_subtype not in content_subtype_dict\
               or content_subtype_dict[content_subtype] == 0:
                continue

            content_subtype_len = content_subtype_dict[content_subtype]

            if content_subtype_len == 1:  # enable to process text data in place

                _reserved_list_ids = {k: list(v) for k, v in reserved_list_ids.items()}
                list_id = get_first_sf_tag(sf, 'ID')
                _reserved_list_ids[content_subtype].remove(int(list_id) if list_id not in EMPTY_VALUE else idx)

                self.__reg.dpR.remediateRawTextPk(sf, file_type, data_file_name, text_data, _reserved_list_ids)

            else:
                multiple_pk_sf_idx.append(idx)

        if len(multiple_pk_sf_idx) > 0:  # need to split into multiple peak lists

            for idx, sf in enumerate(self.__reg.star_data[0].get_saveframes_by_category(sf_category), start=1):

                if idx not in multiple_pk_sf_idx:
                    continue

                reserved_list_ids = get_reserved_list_ids()
                list_id = get_first_sf_tag(sf, 'ID')
                reserved_list_ids[content_subtype].remove(int(list_id) if list_id not in EMPTY_VALUE else idx)

                text_data = get_first_sf_tag(sf, 'Text_data')
                data_format = get_first_sf_tag(sf, 'Text_data_format')
                data_file_name = get_first_sf_tag(sf, 'Data_file_name')

                file_type = content_subtype_dict = None

                if data_format == 'ARIA':
                    file_type = 'nm-pea-ari'
                elif data_format == 'CCPN':
                    file_type = 'nm-pea-ccp'
                elif data_format == 'Olivia':
                    file_type = 'nm-pea-oli'
                elif data_format == 'NMRPipe':
                    file_type = 'nm-pea-pip'
                elif data_format == 'PONDEROSA':
                    file_type = 'nm-pea-pon'
                elif data_format == 'Sparky':
                    file_type = 'nm-pea-spa'
                elif data_format == 'TopSpin':
                    file_type = 'nm-pea-top'
                elif data_format == 'NMRView':
                    file_type = 'nm-pea-vie'
                elif data_format == 'VNMR':
                    file_type = 'nm-pea-vnm'
                elif data_format == 'XwinNMR':
                    file_type = 'nm-pea-xwi'
                else:

                    for _file_type in PARSABLE_PK_FILE_TYPES:

                        if _file_type.startswith('nm-aux'):
                            continue

                        reader = self.__reg.dpS.getSimpleFileReader(_file_type, False)

                        listener, parser_err_listener, lexer_err_listener = reader.parse(text_data, None, isFilePath=False)

                        has_lexer_error = lexer_err_listener is not None and lexer_err_listener.getMessageList() is not None
                        has_parser_error = parser_err_listener is not None and parser_err_listener.getMessageList() is not None
                        content_subtype_dict = listener.getContentSubtype() if listener is not None else None
                        if not has_lexer_error and not has_parser_error\
                           and content_subtype_dict is not None and len(content_subtype_dict) > 0:
                            file_type = _file_type
                            break

                if file_type is None:
                    continue

                self.__reg.dpR.remediateRawTextPk(sf, file_type, data_file_name, text_data, reserved_list_ids)

        return True

    def __updateAuthSequence(self) -> bool:
        """ Update auth sequence in NMR-STAR.
        """

        if not self.__reg.combined_mode:
            return True

        if len(self.__reg.star_data) == 0:
            return False

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_type = input_source_dic['file_type']

        if file_type == 'nef':
            return True

        has_poly_seq = has_key_value(input_source_dic, 'polymer_sequence')

        if not has_poly_seq:
            return False

        if self.__reg.srcPath == self.__reg.dstPath:
            return True

        chain_assign_dic = self.__reg.report.chain_assignment.get()

        if 'nmr_poly_seq_vs_model_poly_seq' not in chain_assign_dic:
            return False

        if not has_key_value(chain_assign_dic, 'nmr_poly_seq_vs_model_poly_seq'):
            return False

        if input_source_dic['content_subtype'] is None:
            return False

        seqAlignMap = {}

        poly_seq = input_source_dic['polymer_sequence']

        for ps in poly_seq:
            chain_id = ps['chain_id']
            seqAlignMap[chain_id] = self.__reg.report.getSequenceAlignmentWithNmrChainId(chain_id)

        if len(seqAlignMap) == 0:
            return False

        tags = ['Entity_assembly_ID', 'Comp_index_ID', 'Auth_asym_ID', 'Auth_seq_ID']

        self.__authSeqMap.clear()

        content_subtype = 'poly_seq'

        sf_category = SF_CATEGORIES[file_type][content_subtype]
        lp_category = LP_CATEGORIES[file_type][content_subtype]

        for sf in self.__reg.star_data[0].get_saveframes_by_category(sf_category):

            try:
                loop = sf.get_loop(lp_category)
            except KeyError:
                continue

            star_chain_index = loop.tags.index(tags[0])
            star_seq_index = loop.tags.index(tags[1])

            for row in loop:
                star_chain = row[star_chain_index]
                star_seq = row[star_seq_index]

                if star_chain in seqAlignMap:
                    seq_align = seqAlignMap[star_chain]

                    if seq_align is None:
                        continue

                    try:
                        auth_seq = seq_align['test_seq_id'][seq_align['ref_seq_id'].index(star_seq)]
                        self.__authSeqMap[(star_chain, star_seq)] = (seq_align['test_chain_id'], auth_seq)
                    except (IndexError, ValueError):
                        pass

        if len(self.__authSeqMap) == 0:
            return False

        def update_auth_seq(lp, _tags_):
            # Entity_assembly_ID*
            star_chain_index = lp.tags.index(_tags_[0])
            # Comp_index_ID*
            star_seq_index = lp.tags.index(_tags_[1])
            # Auth_asym_ID*
            auth_chain_index = lp.tags.index(_tags_[2])
            # Auth_seq_ID*
            auth_seq_index = lp.tags.index(_tags_[3])

            for row in lp:
                star_chain = row[star_chain_index]
                star_seq = row[star_seq_index]

                if star_chain in EMPTY_VALUE or star_seq in EMPTY_VALUE:
                    continue

                seq_key = (star_chain, star_seq)

                if seq_key in self.__authSeqMap:
                    row[auth_chain_index], row[auth_seq_index] = self.__authSeqMap[seq_key]

        for content_subtype in input_source_dic['content_subtype']:

            if content_subtype in ('entry_info', 'entity'):
                continue

            sf_category = SF_CATEGORIES[file_type][content_subtype]
            lp_category = LP_CATEGORIES[file_type][content_subtype]

            for sf in self.__reg.star_data[0].get_saveframes_by_category(sf_category):

                try:
                    loop = sf.get_loop(lp_category)
                except KeyError:
                    continue

                if set(tags) & set(loop.tags) == set(tags):
                    update_auth_seq(loop, tags)

                else:
                    for i in range(1, MAX_DIM_NUM_OF_SPECTRA):
                        _tags = [f'{t}_{i}' for t in tags]

                        if set(_tags) & set(loop.tags) == set(_tags):
                            update_auth_seq(loop, _tags)
                        else:
                            break

        return True

    def __testTautomerOfHistidinePerModel(self) -> bool:
        """ Check tautomeric state of a given histidine per model. (DAOTHER-9252)
        """

        src_id = self.__reg.report.getInputSourceIdOfCoord()

        if src_id < 0:
            return False

        cif_input_source = self.__reg.report.input_sources[src_id]
        cif_input_source_dic = cif_input_source.get()

        has_poly_seq = has_key_value(cif_input_source_dic, 'polymer_sequence')

        if not has_poly_seq:
            return False

        file_name = cif_input_source_dic['file_name']
        cif_poly_seq = cif_input_source_dic['polymer_sequence']

        if len(self.__reg.cpC['tautomer_per_model']) > 0:

            for inst in self.__reg.cpC['tautomer_per_model']:
                tautomer_per_model = inst['tautomer_per_model']

                try:
                    rep_tautomer = tautomer_per_model[self.__reg.representative_model_id]
                except KeyError:
                    try:
                        rep_tautomer = tautomer_per_model[self.__reg.eff_model_ids[0]]
                    except KeyError:
                        continue

                if any(tautomer != rep_tautomer for tautomer in tautomer_per_model.values()):
                    chain_id, auth_chain_id = inst['chain_id'], inst['auth_chain_id']
                    seq_id, auth_seq_id = inst['seq_id'], inst['auth_seq_id']
                    comp_id = inst['comp_id']
                    cif_seq_code = f"{chain_id}:{seq_id}:{comp_id}"
                    if chain_id != auth_chain_id or seq_id != auth_seq_id:
                        cif_seq_code += f" ({auth_chain_id}:{auth_seq_id}:{comp_id} in author sequence scheme)"

                    err = f"{cif_seq_code} has been instantiated with different tautomeric states across models, "\
                        f"{tautomer_per_model}. Please re-upload the model file."

                    if self.__reg.internal_mode:  # and not self.__reg.conversion_server:

                        self.__reg.report.warning.appendDescription('coordinate_issue',
                                                                    {'file_name': file_name, 'category': 'atom_site',
                                                                     'description': err})

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__testTautomerOfHistidinePerModel() "
                                                 f"++ Warning  - {err}\n")

                    else:

                        self.__reg.report.error.appendDescription('coordinate_issue',
                                                                  {'file_name': file_name, 'category': 'atom_site',
                                                                   'description': err})

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__testTautomerOfHistidinePerModel() ++ Error  - {err}\n")

            return True

        model_num_name = 'pdbx_PDB_model_num' if 'pdbx_PDB_model_num' in self.__reg.coord_atom_site_tags else 'ndb_model'

        for ps in cif_poly_seq:
            chain_id = ps['chain_id']

            auth_chain_id = chain_id
            if 'auth_chain_id' in ps:
                auth_chain_id = ps['auth_chain_id']

            if len(cif_poly_seq) >= LEN_MAJOR_ASYM_ID:
                if auth_chain_id not in LARGE_ASYM_ID:
                    continue

            for seq_id, comp_id in zip(ps['seq_id'], ps['comp_id']):

                if not isLikeHis(comp_id, self.__reg.ccU):
                    continue

                if comp_id == 'HIS':
                    hd1_name = 'HD1'
                    he2_name = 'HE2'
                else:
                    _hd1_name = self.__reg.ccU.getBondedAtoms(comp_id, 'ND1', onlyProton=True)
                    _he2_name = self.__reg.ccU.getBondedAtoms(comp_id, 'NE2', onlyProton=True)
                    if len(_hd1_name) != 1 or len(_he2_name) != 1:
                        continue
                    hd1_name = _hd1_name[0]
                    he2_name = _he2_name[0]

                try:
                    auth_seq_id = ps['auth_seq_id'][ps['seq_id'].index(seq_id)]
                except (KeyError, IndexError, ValueError):
                    auth_seq_id = seq_id

                try:

                    protons = self.__reg.cR.getDictListWithFilter('atom_site',
                                                                  [{'name': 'label_atom_id', 'type': 'starts-with-alnum',
                                                                    'alt_name': 'atom_id'},
                                                                   {'name': model_num_name, 'type': 'int', 'alt_name': 'model_id'},
                                                                   ],
                                                                  [{'name': 'label_asym_id', 'type': 'str', 'value': chain_id},
                                                                   {'name': 'label_seq_id', 'type': 'int', 'value': seq_id},
                                                                   {'name': 'label_comp_id', 'type': 'str', 'value': comp_id},
                                                                   {'name': 'type_symbol', 'type': 'str', 'value': 'H'},
                                                                   {'name': 'label_alt_id', 'type': 'enum',
                                                                    'enum': (self.__reg.representative_alt_id,)}
                                                                   ])

                except Exception as e:  # pylint: disable=broad-exception-caught

                    self.__reg.report.error.appendDescription('internal_error',
                                                              f"+{self.__class_name__}.__testTautomerOfHistidinePerModel() "
                                                              "++ Error  - " + str(e))

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__testTautomerOfHistidinePerModel() ++ Error  - {str(e)}\n")

                    return False

                if len(protons) > 0:

                    tautomer_per_model = {}

                    for model_id in self.__reg.eff_model_ids:

                        _protons = [h for h in protons if h['model_id'] == model_id]

                        has_hd1 = has_he2 = False

                        for h in _protons:
                            if h['atom_id'] == hd1_name:
                                has_hd1 = True
                            elif h['atom_id'] == he2_name:
                                has_he2 = True

                        if has_hd1 and has_he2:
                            tautomer_per_model[model_id] = 'biprotonated'

                        elif has_hd1:
                            tautomer_per_model[model_id] = 'pi-tautomer'

                        elif has_he2:
                            tautomer_per_model[model_id] = 'tau-tautomer'

                        else:
                            tautomer_per_model[model_id] = 'unknown'

                    try:
                        rep_tautomer = tautomer_per_model[self.__reg.representative_model_id]
                    except KeyError:
                        try:
                            rep_tautomer = tautomer_per_model[self.__reg.eff_model_ids[0]]
                        except KeyError:
                            continue

                    self.__reg.cpC['tautomer_per_model'].append({'chain_id': chain_id, 'seq_id': seq_id, 'comp_id': comp_id,
                                                                 'auth_chain_id': auth_chain_id, 'auth_seq_id': auth_seq_id,
                                                                 'tautomer_per_model': tautomer_per_model})

                    if any(tautomer != rep_tautomer for tautomer in tautomer_per_model.values()):
                        cif_seq_code = f"{chain_id}:{seq_id}:{comp_id}"
                        if chain_id != auth_chain_id or seq_id != auth_seq_id:
                            cif_seq_code += f" ({auth_chain_id}:{auth_seq_id}:{comp_id} in author sequence scheme)"

                        err = f"{cif_seq_code} has been instantiated with different tautomeric states across models, "\
                            f"{tautomer_per_model}. Please re-upload the model file."

                        if self.__reg.internal_mode:  # and not self.__reg.conversion_server:

                            self.__reg.report.warning.appendDescription('coordinate_issue',
                                                                        {'file_name': file_name, 'category': 'atom_site',
                                                                         'description': err})

                            if self.__reg.verbose:
                                self.__reg.log.write(f"+{self.__class_name__}.__testTautomerOfHistidinePerModel() "
                                                     f"++ Warning  - {err}\n")

                        else:

                            self.__reg.report.error.appendDescription('coordinate_issue',
                                                                      {'file_name': file_name, 'category': 'atom_site',
                                                                       'description': err})

                            if self.__reg.verbose:
                                self.__reg.log.write(f"+{self.__class_name__}.__testTautomerOfHistidinePerModel() "
                                                     f"++ Error  - {err}\n")

        if self.__reg.coordPropCachePath is not None:
            hash_value = hash(str(self.__reg.cpC))
            if hash_value != self.__reg.cpcHashCode:
                write_as_pickle(self.__reg.cpC, self.__reg.coordPropCachePath)
                self.__reg.cpcHashCode = hash_value

        return True

    def __extractCoordDisulfideBond(self) -> bool:
        """ Extract disulfide bond of coordinate file.
        """

        src_id = self.__reg.report.getInputSourceIdOfCoord()

        if src_id < 0:
            return False

        cif_input_source = self.__reg.report.input_sources[src_id]

        chain_assign_dic = self.__reg.report.chain_assignment.get()

        if 'model_poly_seq_vs_nmr_poly_seq' not in chain_assign_dic:

            err = "Chain assignment does not exist, __assignCoordPolymerSequence() should be invoked."

            self.__reg.report.error.appendDescription('internal_error',
                                                      f"+{self.__class_name__}.__extractCoordDisulfideBond() ++ Error  - " + err)

            if self.__reg.verbose:
                self.__reg.log.write(f"+{self.__class_name__}.__extractCoordDisulfideBond() ++ Error  - {err}\n")

            return False

        if not has_key_value(chain_assign_dic, 'model_poly_seq_vs_nmr_poly_seq'):
            return False

        try:

            if self.__reg.cR.hasCategory('struct_conn'):
                struct_conn = self.__reg.cR.getDictListWithFilter('struct_conn',
                                                                  [{'name': 'conn_type_id', 'type': 'str'},
                                                                   {'name': 'ptnr1_label_asym_id', 'type': 'str'},
                                                                   {'name': 'ptnr1_label_seq_id', 'type': 'int'},
                                                                   {'name': 'ptnr1_label_comp_id', 'type': 'str'},
                                                                   {'name': 'ptnr1_label_atom_id', 'type': 'str'},
                                                                   {'name': 'ptnr2_label_asym_id', 'type': 'str'},
                                                                   {'name': 'ptnr2_label_seq_id', 'type': 'int'},
                                                                   {'name': 'ptnr2_label_comp_id', 'type': 'str'},
                                                                   {'name': 'ptnr2_label_atom_id', 'type': 'str'},
                                                                   {'name': 'pdbx_dist_value', 'type': 'float'}
                                                                   ])

            else:
                struct_conn = []

        except Exception as e:  # pylint: disable=broad-exception-caught

            self.__reg.report.error.appendDescription('internal_error',
                                                      f"+{self.__class_name__}.__extractCoordDisulfideBond() ++ Error  - " + str(e))

            if self.__reg.verbose:
                self.__reg.log.write(f"+{self.__class_name__}.__extractCoordDisulfideBond() ++ Error  - {str(e)}\n")

            return False

        if len(struct_conn) > 0:

            asm = []

            for sc in struct_conn:

                if sc['conn_type_id'] != 'disulf':
                    continue

                disulf = {}
                disulf['chain_id_1'] = sc['ptnr1_label_asym_id']
                disulf['seq_id_1'] = sc['ptnr1_label_seq_id']
                disulf['comp_id_1'] = sc['ptnr1_label_comp_id']
                disulf['atom_id_1'] = sc['ptnr1_label_atom_id']
                disulf['chain_id_2'] = sc['ptnr2_label_asym_id']
                disulf['seq_id_2'] = sc['ptnr2_label_seq_id']
                disulf['comp_id_2'] = sc['ptnr2_label_comp_id']
                disulf['atom_id_2'] = sc['ptnr2_label_atom_id']
                disulf['distance_value'] = sc['pdbx_dist_value']
                # DAOTHER-7475
                disulf['warning_description_1'] = None
                disulf['warning_description_2'] = None
                asm.append(disulf)

            if len(asm) > 0:
                cif_input_source.setItemValue('disulfide_bond', asm)

                self.__reg.report.setDisulfideBond(True)

                return self.__reg.dpV.mapCoordDisulfideBond2Nmr(asm)

        return True

    def __extractCoordOtherBond(self) -> bool:
        """ Extract other bond (neither disulfide nor covalent bond) of coordinate file.
        """

        src_id = self.__reg.report.getInputSourceIdOfCoord()

        if src_id < 0:
            return False

        cif_input_source = self.__reg.report.input_sources[src_id]

        chain_assign_dic = self.__reg.report.chain_assignment.get()

        if 'model_poly_seq_vs_nmr_poly_seq' not in chain_assign_dic:

            err = "Chain assignment does not exist, __assignCoordPolymerSequence() should be invoked."

            self.__reg.report.error.appendDescription('internal_error',
                                                      f"+{self.__class_name__}.__extractCoordOtherBond() ++ Error  - " + err)

            if self.__reg.verbose:
                self.__reg.log.write(f"+{self.__class_name__}.__extractCoordOtherBond() ++ Error  - {err}\n")

            return False

        if not has_key_value(chain_assign_dic, 'model_poly_seq_vs_nmr_poly_seq'):
            return False

        try:

            if self.__reg.cR.hasCategory('struct_conn'):
                struct_conn = self.__reg.cR.getDictListWithFilter('struct_conn',
                                                                  [{'name': 'conn_type_id', 'type': 'str'},
                                                                   {'name': 'ptnr1_label_asym_id', 'type': 'str'},
                                                                   {'name': 'ptnr1_label_seq_id', 'type': 'int'},
                                                                   {'name': 'ptnr1_label_comp_id', 'type': 'str'},
                                                                   {'name': 'ptnr1_label_atom_id', 'type': 'str'},
                                                                   {'name': 'ptnr2_label_asym_id', 'type': 'str'},
                                                                   {'name': 'ptnr2_label_seq_id', 'type': 'int'},
                                                                   {'name': 'ptnr2_label_comp_id', 'type': 'str'},
                                                                   {'name': 'ptnr2_label_atom_id', 'type': 'str'},
                                                                   {'name': 'pdbx_dist_value', 'type': 'float'}
                                                                   ])

            else:
                struct_conn = []

        except Exception as e:  # pylint: disable=broad-exception-caught

            self.__reg.report.error.appendDescription('internal_error',
                                                      f"+{self.__class_name__}.__extractCoordOtherBond() ++ Error  - " + str(e))

            if self.__reg.verbose:
                self.__reg.log.write(f"+{self.__class_name__}.__extractCoordOtherBond() ++ Error  - {str(e)}\n")

            return False

        if len(struct_conn) > 0:

            asm = []

            for sc in struct_conn:

                if sc['conn_type_id'] in ('disulf', 'hydrog') or sc['conn_type_id'].startswith('covale'):
                    continue

                other = {}
                other['bond_type'] = sc['conn_type_id']
                other['chain_id_1'] = sc['ptnr1_label_asym_id']
                other['seq_id_1'] = sc['ptnr1_label_seq_id']
                other['comp_id_1'] = sc['ptnr1_label_comp_id']
                other['atom_id_1'] = sc['ptnr1_label_atom_id']
                other['chain_id_2'] = sc['ptnr2_label_asym_id']
                other['seq_id_2'] = sc['ptnr2_label_seq_id']
                other['comp_id_2'] = sc['ptnr2_label_comp_id']
                other['atom_id_2'] = sc['ptnr2_label_atom_id']
                other['distance_value'] = sc['pdbx_dist_value']
                # DAOTHER-7475
                other['warning_description_1'] = None
                other['warning_description_2'] = None
                asm.append(other)

            if len(asm) > 0:
                cif_input_source.setItemValue('other_bond', asm)

                self.__reg.report.setOtherBond(True)

                return self.__reg.dpV.mapCoordOtherBond2Nmr(asm)

        return True

    def __appendElemAndIsoNumOfNefCsLoop(self) -> bool:
        """ Append element and isotope_number columns in NEF CS loop if required.
        """

        if not self.__reg.combined_mode:
            return True

        try:

            for fileListId in range(self.__reg.file_path_list_len):

                input_source = self.__reg.report.input_sources[fileListId]
                input_source_dic = input_source.get()

                file_type = input_source_dic['file_type']

                if file_type != 'nef':
                    continue

                content_subtype = 'chem_shift'

                if not has_key_value(input_source_dic['content_subtype'], content_subtype):
                    continue

                sf_category = SF_CATEGORIES[file_type][content_subtype]
                lp_category = LP_CATEGORIES[file_type][content_subtype]

                cs_item_names = ITEM_NAMES_IN_CS_LOOP[file_type]
                cs_atom_id_name = cs_item_names['atom_id']
                cs_atom_type = cs_item_names['atom_type']
                cs_iso_number = cs_item_names['isotope_number']

                for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):

                    try:
                        loop = sf.get_loop(lp_category)
                    except KeyError:
                        continue

                    has_atom_type = cs_atom_type in loop.tags
                    has_iso_number = cs_iso_number in loop.tags

                    atomIdCol = loop.tags.index(cs_atom_id_name)

                    if has_atom_type and has_iso_number:

                        atomTypeCol = loop.tags.index(cs_atom_type)
                        isoNumCol = loop.tags.index(cs_iso_number)

                        for row in loop:

                            atom_id = row[atomIdCol]

                            if row[atomTypeCol] in EMPTY_VALUE:
                                row[atomTypeCol] = atom_id[0]

                            if row[isoNumCol] in EMPTY_VALUE:

                                try:
                                    row[isoNumCol] = ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS[atom_id[0]][0]
                                except KeyError:
                                    pass

                    elif has_atom_type:

                        atomTypeCol = loop.tags.index(cs_atom_type)

                        for row in loop:

                            atom_id = row[atomIdCol]

                            if row[atomTypeCol] in EMPTY_VALUE:
                                row[atomTypeCol] = atom_id[0]

                            try:
                                iso_num = ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS[atom_id[0]][0]
                                row.append(iso_num)
                            except KeyError:
                                row.append('.')

                        loop.add_tag(cs_iso_number)

                    elif has_iso_number:

                        isoNumCol = loop.tags.index(cs_iso_number)

                        for row in loop:

                            atom_id = row[atomIdCol]

                            if row[isoNumCol] in EMPTY_VALUE:

                                try:
                                    row[isoNumCol] = ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS[atom_id[0]][0]
                                except KeyError:
                                    pass

                            row.append(atom_id[0] if atom_id[0] in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS else '.')

                        loop.add_tag(cs_atom_type)

                    else:

                        for row in loop:

                            atom_id = row[atomIdCol]

                            row.append(atom_id[0] if atom_id[0] in ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS else '.')

                            try:
                                iso_num = ISOTOPE_NUMBERS_OF_NMR_OBS_NUCS[atom_id[0]][0]
                                row.append(iso_num)
                            except KeyError:
                                row.append('.')

                        loop.add_tag(cs_atom_type)
                        loop.add_tag(cs_iso_number)

            return True

        except ValueError as e:

            self.__reg.report.error.appendDescription('internal_error',
                                                      f"+{self.__class_name__}.__appendElemAndIsoNumOfNefCsLoop() "
                                                      "++ Error  - " + str(e))

            if self.__reg.verbose:
                self.__reg.log.write(f"+{self.__class_name__}.__appendElemAndIsoNumOfNefCsLoop() ++ Error  - {str(e)}\n")

            return False

    def __appendWeightInLoop(self) -> bool:
        """ Append weight column in interesting loops, if required.
        """

        if not self.__reg.combined_mode:
            return True

        try:

            is_done = True

            for fileListId in range(self.__reg.file_path_list_len):

                input_source = self.__reg.report.input_sources[fileListId]
                input_source_dic = input_source.get()

                file_name = input_source_dic['file_name']
                file_type = input_source_dic['file_type']

                if input_source_dic['content_subtype'] is None:
                    is_done = False
                    continue

                for content_subtype in input_source_dic['content_subtype']:

                    if content_subtype == 'entity':
                        continue

                    sf_category = SF_CATEGORIES[file_type][content_subtype]
                    lp_category = LP_CATEGORIES[file_type][content_subtype]

                    weight_tag = WEIGHT_TAGS[file_type][content_subtype]

                    if weight_tag is None:
                        continue

                    for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                        sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                        try:
                            loop = sf.get_loop(lp_category)
                        except KeyError:
                            continue

                        if weight_tag in loop.tags:
                            continue

                        lp_tag = f'{lp_category}.{weight_tag}'
                        err = ERR_TEMPLATE_FOR_MISSING_MANDATORY_LP_TAG % (lp_tag, file_type.upper())

                        if self.__reg.check_mandatory_tag and self.__reg.nefT.is_mandatory_tag(lp_tag, file_type):

                            if self.__reg.rescue_mode:
                                self.__reg.report.error.appendDescription('missing_mandatory_item',
                                                                          {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                           'category': lp_category, 'description': err})

                                if self.__reg.verbose:
                                    self.__reg.log.write(f"+{self.__class_name__}.__appendWeightInLoop() "
                                                         f"++ LookupError  - {file_name} {sf_framecode} {lp_category} {err}\n")

                        for row in loop:
                            row.append('1.0')

                        loop.add_tag(weight_tag)

            return is_done

        except ValueError as e:

            self.__reg.report.error.appendDescription('internal_error',
                                                      f"+{self.__class_name__}.__appendWeightInLoop() ++ Error  - " + str(e))

            if self.__reg.verbose:
                self.__reg.log.write(f"+{self.__class_name__}.__appendWeightInLoop() ++ Error  - {str(e)}\n")

            return False

    def __appendDihedAngleType(self) -> bool:
        """ Append dihedral angle type column, if required.
        """

        if not self.__reg.combined_mode:
            return True

        try:

            for fileListId in range(self.__reg.file_path_list_len):

                input_source = self.__reg.report.input_sources[fileListId]
                input_source_dic = input_source.get()

                file_type = input_source_dic['file_type']

                content_subtype = 'dihed_restraint'

                if not has_key_value(input_source_dic['content_subtype'], content_subtype):
                    continue

                sf_category = SF_CATEGORIES[file_type][content_subtype]
                lp_category = LP_CATEGORIES[file_type][content_subtype]

                angle_type_tag = ANGLE_TYPE_TAGS[file_type]

                for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):

                    try:
                        loop = sf.get_loop(lp_category)
                    except KeyError:
                        continue

                    if angle_type_tag in loop.tags:
                        continue

                    loop.add_tag(angle_type_tag, update_data=True)

            return True

        except ValueError as e:

            self.__reg.report.error.appendDescription('internal_error',
                                                      f"+{self.__class_name__}.__appendDihedAngleType() ++ Error  - " + str(e))

            if self.__reg.verbose:
                self.__reg.log.write(f"+{self.__class_name__}.__appendDihedAngleType() ++ Error  - {str(e)}\n")

            return False

    def __appendSfTagItem(self) -> bool:
        """ Append saveframe tag items, if required.
        """

        if not self.__reg.combined_mode:
            return True

        try:

            for fileListId in range(self.__reg.file_path_list_len):

                input_source = self.__reg.report.input_sources[fileListId]
                input_source_dic = input_source.get()

                file_name = input_source_dic['file_name']
                file_type = input_source_dic['file_type']

                if input_source_dic['content_subtype'] is None:
                    continue

                for content_subtype in input_source_dic['content_subtype']:

                    if content_subtype in ('entry_info', 'poly_seq', 'entity'):
                        continue

                    sf_category = SF_CATEGORIES[file_type][content_subtype]

                    has_data_file_name = file_type == 'nmr-star' and 'Data_file_name' in SF_ALLOWED_TAGS[file_type][content_subtype]
                    tag_items = MANDATORY_SF_TAG_ITEMS[file_type][content_subtype]

                    if not has_data_file_name and tag_items is None:
                        continue

                    for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                        sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                        tagNames = [t[0] for t in sf.tags]

                        # DAOTHER-9520: fill data_file_name by default
                        if has_data_file_name:
                            original_file_name = input_source_dic['file_name']
                            if 'original_file_name' in input_source_dic and input_source_dic['original_file_name'] is not None:
                                original_file_name = os.path.basename(input_source_dic['original_file_name'])

                            if 'Data_file_name' in tagNames:
                                data_file_name = get_first_sf_tag(sf, 'Data_file_name').strip()
                                if len(data_file_name) == 0 or data_file_name in EMPTY_VALUE:
                                    set_sf_tag(sf, 'Data_file_name', original_file_name)
                            else:
                                set_sf_tag(sf, 'Data_file_name', original_file_name)

                        if tag_items is None:
                            continue

                        for tag_item in tag_items:

                            if tag_item in tagNames:
                                continue

                            sf_tag = f'_{sf_category}.{tag_item}'
                            warn = WARN_TEMPLATE_FOR_MISSING_MANDATORY_SF_TAG % (sf_tag, file_type.upper())

                            if self.__reg.check_mandatory_tag and self.__reg.nefT.is_mandatory_tag(sf_tag, file_type):

                                if self.__reg.rescue_mode:
                                    self.__reg.report.warning.appendDescription('missing_data',
                                                                                {'file_name': file_name,
                                                                                 'sf_framecode': sf_framecode,
                                                                                 'category': sf_category, 'description': warn})

                                    if self.__reg.verbose:
                                        self.__reg.log.write(f"+{self.__class_name__}.__appendSfTagItem() ++ Warning  - {warn}\n")

                            sf.add_tag(tag_item, '.')

            return True

        except ValueError as e:

            self.__reg.report.error.appendDescription('internal_error',
                                                      f"+{self.__class_name__}.__appendSfTagItem() ++ Error  - " + str(e))

            if self.__reg.verbose:
                self.__reg.log.write(f"+{self.__class_name__}.__appendSfTagItem() ++ Error  - {str(e)}\n")

            return False

    def __updateDihedralAngleType(self) -> bool:
        """ Update dihedral angle types if possible.
        """

        if not self.__reg.combined_mode:
            return True

        for fileListId in range(self.__reg.file_path_list_len):

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_name = input_source_dic['file_name']
            file_type = input_source_dic['file_type']

            content_subtype = 'dihed_restraint'

            if not has_key_value(input_source_dic['content_subtype'], content_subtype):
                continue

            item_names = ITEM_NAMES_IN_DIHED_LOOP[file_type]
            index_id_name = INDEX_TAGS[file_type][content_subtype]
            chain_id_1_name = item_names['chain_id_1']
            chain_id_2_name = item_names['chain_id_2']
            chain_id_3_name = item_names['chain_id_3']
            chain_id_4_name = item_names['chain_id_4']
            seq_id_1_name = item_names['seq_id_1']
            seq_id_2_name = item_names['seq_id_2']
            seq_id_3_name = item_names['seq_id_3']
            seq_id_4_name = item_names['seq_id_4']
            comp_id_1_name = item_names['comp_id_1']
            comp_id_2_name = item_names['comp_id_2']
            comp_id_3_name = item_names['comp_id_3']
            comp_id_4_name = item_names['comp_id_4']
            atom_id_1_name = item_names['atom_id_1']
            atom_id_2_name = item_names['atom_id_2']
            atom_id_3_name = item_names['atom_id_3']
            atom_id_4_name = item_names['atom_id_4']
            angle_type_name = item_names['angle_type']
            lower_limit_name = item_names['lower_limit']
            upper_limit_name = item_names['upper_limit']

            # pylint: disable=cell-var-from-loop
            def ext_atoms(row):
                return ({'chain_id': row[chain_id_1_name], 'seq_id': row[seq_id_1_name],
                         'comp_id': row[comp_id_1_name], 'atom_id': row[atom_id_1_name]},
                        {'chain_id': row[chain_id_2_name], 'seq_id': row[seq_id_2_name],
                         'comp_id': row[comp_id_2_name], 'atom_id': row[atom_id_2_name]},
                        {'chain_id': row[chain_id_3_name], 'seq_id': row[seq_id_3_name],
                         'comp_id': row[comp_id_3_name], 'atom_id': row[atom_id_3_name]},
                        {'chain_id': row[chain_id_4_name], 'seq_id': row[seq_id_4_name],
                         'comp_id': row[comp_id_4_name], 'atom_id': row[atom_id_4_name]})

            sf_category = SF_CATEGORIES[file_type][content_subtype]
            lp_category = LP_CATEGORIES[file_type][content_subtype]

            for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                sf_framecode = get_first_sf_tag(sf, 'sf_framecode')

                lp_data = next((lp['data'] for lp in self.__reg.lp_data[content_subtype]
                                if lp['file_name'] == file_name and lp['sf_framecode'] == sf_framecode), None)

                if lp_data is None:

                    key_items = self.__reg.key_items[file_type][content_subtype]
                    data_items = DATA_ITEMS[file_type][content_subtype]

                    try:

                        lp_data = self.__reg.nefT.check_data(sf, lp_category, key_items, data_items, None, None, None,
                                                             enforce_allowed_tags=(file_type == 'nmr-star'),
                                                             excl_missing_data=self.__reg.excl_missing_data)[0]

                        self.__reg.lp_data[content_subtype].append({'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                    'data': lp_data})

                    except Exception:  # pylint: disable=broad-exception-caught
                        pass

                if lp_data is not None:

                    update = False
                    update_index = {}

                    try:

                        for row in lp_data:
                            index_id = row[index_id_name]

                            atom1, atom2, atom3, atom4 = ext_atoms(row)

                            angle_type = row[angle_type_name]

                            if angle_type not in EMPTY_VALUE:
                                continue

                            peptide, nucleotide, carbohydrate = self.__reg.csStat.getTypeOfCompId(atom2['comp_id'])
                            plane_like = is_like_planality_boundary(row, lower_limit_name, upper_limit_name)

                            data_type = getTypeOfDihedralRestraint(peptide, nucleotide, carbohydrate,
                                                                   [atom1, atom2, atom3, atom4], plane_like)

                            if data_type in EMPTY_VALUE or data_type.startswith('pseudo'):
                                continue

                            update = True

                            if data_type not in update_index:
                                update_index[data_type] = []

                            update_index[data_type].append(index_id)

                    except Exception as e:  # pylint: disable=broad-exception-caught

                        self.__reg.report.error.appendDescription('internal_error',
                                                                  f"+{self.__class_name__}.__updateDihedralAngleType() "
                                                                  "++ Error  - " + str(e))

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__updateDihedralAngleType() "
                                                 f"++ Error  - {str(e)}\n")

                        continue

                    if update:

                        try:
                            loop = sf.get_loop(lp_category)
                        except KeyError:
                            continue

                        idxCol = loop.tags.index(index_id_name)
                        aglCol = loop.tags.index(angle_type_name)

                        for row in loop:

                            index_id = int(row[idxCol])

                            for k, v in update_index.items():
                                if index_id in v:
                                    row[aglCol] = k

        return True

    def __fixDisorderedIndex(self) -> bool:
        """ Fix disordered indices.
        """

        if not self.__reg.combined_mode:
            return True

        if len(self.__reg.star_data) == 0:
            return False

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_name = input_source_dic['file_name']
        file_type = input_source_dic['file_type']

        warnings = self.__reg.report.warning.getValueList('disordered_index', file_name)

        if warnings is None:
            return True

        for w in warnings:

            if self.__reg.star_data_type[0] == 'Entry' or self.__reg.star_data_type[0] == 'Saveframe':

                if 'sf_framecode' not in w:

                    err = "Could not specify 'sf_framecode' in NMR data processing report."

                    self.__reg.report.error.appendDescription('internal_error',
                                                              f"+{self.__class_name__}.__fixDisorderedIndex() ++ Error  - " + err)

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__fixDisorderedIndex() ++ Error  - {err}\n")

                else:

                    sf = self.__reg.dpA.getSaveframeByName(0, w['sf_framecode'])

                    if sf is None:

                        err = f"Could not specify {w['sf_framecode']!r} saveframe unexpectedly in {file_name!r} file."

                        self.__reg.report.error.appendDescription('internal_error',
                                                                  f"+{self.__class_name__}.__fixDisorderedIndex() "
                                                                  "++ Error  - " + err)

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__fixDisorderedIndex() ++ Error  - {err}\n")

                        continue

                    if 'category' not in w:

                        err = "Could not specify 'category' in NMR data processing report."

                        self.__reg.report.error.appendDescription('internal_error',
                                                                  f"+{self.__class_name__}.__fixDisorderedIndex() "
                                                                  "++ Error  - " + err)

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__fixDisorderedIndex() ++ Error  - {err}\n")

                    else:

                        try:

                            category = w['category'] if w['category'].startswith('_') else f"_{w['category']}"  # pynmrstar v2.6.5.1

                            content_subtype = next(c for c in input_source_dic['content_subtype']
                                                   if LP_CATEGORIES[file_type][c] == category
                                                   and INDEX_TAGS[file_type][c] is not None)

                            loop = sf.get_loop(w['category'])

                            loop.renumber_rows(INDEX_TAGS[file_type][content_subtype])

                        except StopIteration:

                            err = "Could not specify content_subtype in NMR data processing report."

                            self.__reg.report.error.appendDescription('internal_error',
                                                                      f"+{self.__class_name__}.__fixDisorderedIndex() "
                                                                      "++ Error  - " + err)

                            if self.__reg.verbose:
                                self.__reg.log.write(f"+{self.__class_name__}.__fixDisorderedIndex() ++ Error  - {err}\n")

            else:

                err = f"Unexpected PyNMRSTAR object type {self.__reg.star_data_type[0]} found about {file_name!r} file."

                self.__reg.report.error.appendDescription('internal_error',
                                                          f"+{self.__class_name__}.__fixDisorderedIndex() ++ Error  - " + err)

                if self.__reg.verbose:
                    self.__reg.log.write(f"+{self.__class_name__}.__fixDisorderedIndex() ++ Error  - {err}\n")

        return True

    def __removeNonSenseZeroValue(self) -> bool:
        """ Remove non-sense zero values.
        """

        if not self.__reg.combined_mode:
            return True

        if len(self.__reg.star_data) == 0:
            return False

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_name = input_source_dic['file_name']

        warnings = self.__reg.report.warning.getValueList('missing_data', file_name)

        if warnings is None:
            return True

        for w in warnings:

            if "should not have zero value" not in w['description']:
                continue

            if self.__reg.star_data_type[0] == 'Entry' or self.__reg.star_data_type[0] == 'Saveframe':

                if 'sf_framecode' not in w:

                    err = "Could not specify 'sf_framecode' in NMR data processing report."

                    self.__reg.report.error.appendDescription('internal_error',
                                                              f"+{self.__class_name__}.__removeNonSenseZeroValue() "
                                                              "++ Error  - " + err)

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__removeNonSenseZeroValue() "
                                             f"++ Error  - {err}\n")

                else:

                    sf = self.__reg.dpA.getSaveframeByName(0, w['sf_framecode'])

                    if sf is None:

                        err = f"Could not specify {w['sf_framecode']!r} saveframe unexpectedly in {file_name!r} file."

                        self.__reg.report.error.appendDescription('internal_error',
                                                                  f"+{self.__class_name__}.__removeNonSenseZeroValue() "
                                                                  "++ Error  - " + err)

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__removeNonSenseZeroValue() "
                                                 f"++ Error  - {err}\n")

                        continue

                    if 'category' not in w:

                        err = "Could not specify 'category' in NMR data processing report."

                        self.__reg.report.error.appendDescription('internal_error',
                                                                  f"+{self.__class_name__}.__removeNonSenseZeroValue() "
                                                                  "++ Error  - " + err)

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__removeNonSenseZeroValue() "
                                                 f"++ Error  - {err}\n")

                    else:

                        itName = w['description'].split(' ')[0]

                        loop = sf.get_loop(w['category'])

                        if itName not in loop.tags:

                            err = f"Could not find loop tag {itName} in {w['category']} category, "\
                                f"{w['sf_framecode']!r} saveframe, {file_name!r} file."

                            self.__reg.report.error.appendDescription('internal_error',
                                                                      f"+{self.__class_name__}.__removeNonSenseZeroValue() "
                                                                      "++ Error  - " + err)

                            if self.__reg.verbose:
                                self.__reg.log.write(f"+{self.__class_name__}.__removeNonSenseZeroValue() "
                                                     f"++ Error  - {err}\n")

                        else:

                            itCol = loop.tags.index(itName)

                            for row in loop:

                                val = row[itCol]

                                if val in EMPTY_VALUE:
                                    continue

                                try:
                                    if float(val) == 0:
                                        row[itCol] = '.'
                                except ValueError:
                                    row[itCol] = '.'

            else:

                err = f"Unexpected PyNMRSTAR object type {self.__reg.star_data_type[0]} found about {file_name!r} file."

                self.__reg.report.error.appendDescription('internal_error',
                                                          f"+{self.__class_name__}.__removeNonSenseZeroValue() "
                                                          "++ Error  - " + err)

                if self.__reg.verbose:
                    self.__reg.log.write(f"+{self.__class_name__}.__removeNonSenseZeroValue() "
                                         f"++ Error  - {err}\n")

        return True

    def __fixNonSenseNegativeValue(self) -> bool:
        """ Fix non-sense negative values.
        """

        if not self.__reg.combined_mode:
            return True

        if len(self.__reg.star_data) == 0:
            return False

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_name = input_source_dic['file_name']

        warnings = self.__reg.report.warning.getValueList('unusual_data', file_name)

        if warnings is None:
            return True

        for w in warnings:

            if "should not have negative value" not in w['description']:
                continue

            if self.__reg.star_data_type[0] == 'Entry' or self.__reg.star_data_type[0] == 'Saveframe':

                if 'sf_framecode' not in w:

                    err = "Could not specify 'sf_framecode' in NMR data processing report."

                    self.__reg.report.error.appendDescription('internal_error',
                                                              f"+{self.__class_name__}.__fixNonSenseNegativeValue() "
                                                              "++ Error  - " + err)

                    if self.__reg.verbose:
                        self.__reg.log.write(f"+{self.__class_name__}.__fixNonSenseNegativeValue() "
                                             f"++ Error  - {err}\n")

                else:

                    sf = self.__reg.dpA.getSaveframeByName(0, w['sf_framecode'])

                    if sf is None:

                        err = f"Could not specify {w['sf_framecode']!r} saveframe unexpectedly in {file_name!r} file."

                        self.__reg.report.error.appendDescription('internal_error',
                                                                  f"+{self.__class_name__}.__fixNonSenseNegativeValue() "
                                                                  "++ Error  - " + err)

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__fixNonSenseNegativeValue() "
                                                 f"++ Error  - {err}\n")

                        continue

                    if 'category' not in w:

                        err = "Could not specify 'category' in NMR data processing report."

                        self.__reg.report.error.appendDescription('internal_error',
                                                                  f"+{self.__class_name__}.__fixNonSenseNegativeValue() "
                                                                  "++ Error  - " + err)

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__fixNonSenseNegativeValue() "
                                                 f"++ Error  - {err}\n")

                    else:

                        itName = w['description'].split(' ')[0]

                        loop = sf.get_loop(w['category'])

                        if itName not in loop.tags:

                            err = f"Could not find loop tag {itName} in {w['category']} category, "\
                                f"{w['sf_framecode']!r} saveframe, {file_name!r} file."

                            self.__reg.report.error.appendDescription('internal_error',
                                                                      f"+{self.__class_name__}.__fixNonSenseNegativeValue() "
                                                                      "++ Error  - " + err)

                            if self.__reg.verbose:
                                self.__reg.log.write(f"+{self.__class_name__}.__fixNonSenseNegativeValue() "
                                                     f"++ Error  - {err}\n")

                        else:

                            itCol = loop.tags.index(itName)

                            for row in loop:

                                val = row[itCol]

                                if val in EMPTY_VALUE:
                                    continue

                                try:
                                    if float(val) < 0.0:
                                        row[itCol] = abs(float(val))
                                except ValueError:
                                    row[itCol] = '.'

            else:

                err = f"Unexpected PyNMRSTAR object type {self.__reg.star_data_type[0]} found about {file_name!r} file."

                self.__reg.report.error.appendDescription('internal_error',
                                                          f"+{self.__class_name__}.__fixNonSenseNegativeValue() "
                                                          "++ Error  - " + err)

                if self.__reg.verbose:
                    self.__reg.log.write(f"+{self.__class_name__}.__fixNonSenseNegativeValue() "
                                         f"++ Error  - {err}\n")

        return True

    def __fixEnumMismatch(self) -> bool:
        """ Fix enumeration mismatches if possible.
        """

        if not self.__reg.combined_mode:
            return True

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_name = input_source_dic['file_name']

        warnings = self.__reg.report.warning.getValueList('enum_mismatch', file_name)

        if warnings is None:
            return True

        return self.__reg.dpR.fixEnumerationFailure(warnings)

    def __fixEnumMismatchIgnorable(self) -> bool:
        """ Fix enumeration mismatches (ignorable) if possible.
        """

        if not self.__reg.combined_mode:
            return True

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_name = input_source_dic['file_name']

        warnings = self.__reg.report.warning.getValueList('enum_mismatch_ignorable', file_name)

        if warnings is None:
            return True

        return self.__reg.dpR.fixEnumerationFailure(warnings)

    def __resetCapitalStringInLoop(self) -> bool:
        """ Reset capital string values (chain_id, comp_id, atom_id) in loops depending on file type.
        """

        if not self.__reg.combined_mode:
            return True

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_type = input_source_dic['file_type']

        if input_source_dic['content_subtype'] is None:
            return False

        for content_subtype in input_source_dic['content_subtype']:

            if content_subtype in ('entry_info', 'entity'):
                continue

            sf_category = SF_CATEGORIES[file_type][content_subtype]
            lp_category = LP_CATEGORIES[file_type][content_subtype]

            for sf in self.__reg.star_data[0].get_saveframes_by_category(sf_category):

                if not any(True for loop in sf.loops if loop.category == lp_category):
                    continue

                if content_subtype == 'spectral_peak':

                    try:

                        _num_dim = get_first_sf_tag(sf, NUM_DIM_ITEMS[file_type])
                        num_dim = int(_num_dim)

                        if num_dim not in range(1, MAX_DIM_NUM_OF_SPECTRA):
                            raise ValueError()

                    except ValueError:  # raised error already at __testIndexConsistency()
                        continue

                    max_dim = num_dim + 1

                    key_items = []
                    for dim in range(1, max_dim):
                        for k in PK_KEY_ITEMS[file_type]:
                            if k['type'] == 'float':  # position
                                _k = copy.copy(k)
                                if '%s' in k['name']:
                                    _k['name'] = k['name'] % dim
                                key_items.append(_k)
                    for k in PK_KEY_ITEMS[file_type]:
                        if k['type'] == 'positive-int':  # peak_id
                            key_items.append(k)

                    data_items = []
                    for d in DATA_ITEMS[file_type][content_subtype]:
                        data_items.append(d)
                    for dim in range(1, max_dim):
                        for d in self.__reg.pk_data_items[file_type]:
                            _d = copy.copy(d)
                            if '%s' in d['name']:
                                _d['name'] = d['name'] % dim
                            if 'default-from' in d and '%s' in d['default-from']:  # DAOTHER-7421
                                _d['default-from'] = d['default-from'] % dim
                            data_items.append(_d)

                    if max_dim < MAX_DIM_NUM_OF_SPECTRA:
                        disallowed_tags = []
                        for dim in range(max_dim, MAX_DIM_NUM_OF_SPECTRA):
                            for t in DISALLOWED_PK_TAGS[file_type]:
                                if '%s' in t:
                                    t = t % dim
                                disallowed_tags.append(t)

                else:

                    key_items = self.__reg.key_items[file_type][content_subtype]
                    data_items = DATA_ITEMS[file_type][content_subtype]

                loop = sf.get_loop(lp_category)

                if file_type == 'nef':
                    key_names = [k['name'] for k in key_items
                                 if k['name'].startswith('chain_code') or k['name'].startswith('residue_name')
                                 or k['name'].startswith('atom_name') or k['name'] == 'element']
                else:
                    key_names = [k['name'] for k in key_items
                                 if k['name'].startswith('Comp_ID') or k['name'].startswith('Atom_ID') or k['name'] == 'Atom_type']

                for itName in key_names:

                    if itName in loop.tags:

                        itCol = loop.tags.index(itName)

                        for row in loop:

                            val = row[itCol]

                            if val in EMPTY_VALUE:
                                continue

                            if (file_type == 'nef' and itName.startswith('atom_name'))\
                               or (file_type == 'nmr-star'
                                   and (itName.startswith('Auth_atom_ID') or itName == 'Original_PDB_atom_name')):
                                continue

                            row[itCol] = val.upper()

                if file_type == 'nef':
                    data_names = [d['name'] for d in data_items
                                  if d['name'].startswith('chain_code') or d['name'].startswith('residue_name')
                                  or d['name'].startswith('atom_name') or d['name'] == 'element']
                else:
                    data_names = [d['name'] for d in data_items
                                  if d['name'].startswith('Comp_ID') or d['name'].startswith('Atom_ID') or d['name'] == 'Atom_type']

                for itName in data_names:

                    if itName in loop.tags:

                        itCol = loop.tags.index(itName)

                        for row in loop:

                            val = row[itCol]

                            if val in EMPTY_VALUE:
                                continue

                            if (file_type == 'nef' and itName.startswith('atom_name'))\
                               or (file_type == 'nmr-star'
                                   and (itName.startswith('Auth_atom_ID') or itName == 'Original_PDB_atom_name')):
                                continue

                            row[itCol] = val.upper()

        return True

    def __resetBoolValueInLoop(self) -> bool:
        """ Reset bool values in loops depending on file type.
        """

        if not self.__reg.combined_mode:
            return True

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_type = input_source_dic['file_type']

        yes_value = 'true' if file_type == 'nef' else 'yes'
        no_value = 'false' if file_type == 'nef' else 'no'

        if input_source_dic['content_subtype'] is None:
            return False

        for content_subtype in input_source_dic['content_subtype']:

            if content_subtype in ('entry_info', 'entity'):
                continue

            sf_category = SF_CATEGORIES[file_type][content_subtype]
            lp_category = LP_CATEGORIES[file_type][content_subtype]

            for sf in self.__reg.star_data[0].get_saveframes_by_category(sf_category):

                if not any(True for loop in sf.loops if loop.category == lp_category):
                    continue

                if content_subtype == 'spectral_peak':

                    try:

                        _num_dim = get_first_sf_tag(sf, NUM_DIM_ITEMS[file_type])
                        num_dim = int(_num_dim)

                        if num_dim not in range(1, MAX_DIM_NUM_OF_SPECTRA):
                            raise ValueError()

                    except ValueError:  # raised error already at __testIndexConsistency()
                        continue

                    max_dim = num_dim + 1

                    key_items = []
                    for dim in range(1, max_dim):
                        for k in PK_KEY_ITEMS[file_type]:
                            if k['type'] == 'float':  # position
                                _k = copy.copy(k)
                                if '%s' in k['name']:
                                    _k['name'] = k['name'] % dim
                                key_items.append(_k)
                    for k in PK_KEY_ITEMS[file_type]:
                        if k['type'] == 'positive-int':  # peak_id
                            key_items.append(k)

                    data_items = []
                    for d in DATA_ITEMS[file_type][content_subtype]:
                        data_items.append(d)
                    for dim in range(1, max_dim):
                        for d in self.__reg.pk_data_items[file_type]:
                            _d = copy.copy(d)
                            if '%s' in d['name']:
                                _d['name'] = d['name'] % dim
                            if 'default-from' in d and '%s' in d['default-from']:  # DAOTHER-7421
                                _d['default-from'] = d['default-from'] % dim
                            data_items.append(_d)

                    if max_dim < MAX_DIM_NUM_OF_SPECTRA:
                        disallowed_tags = []
                        for dim in range(max_dim, MAX_DIM_NUM_OF_SPECTRA):
                            for t in DISALLOWED_PK_TAGS[file_type]:
                                if '%s' in t:
                                    t = t % dim
                                disallowed_tags.append(t)

                else:

                    key_items = self.__reg.key_items[file_type][content_subtype]
                    data_items = DATA_ITEMS[file_type][content_subtype]

                has_bool_key = has_bool_data = False

                if key_items is not None:
                    has_bool_key = next((k['type'] == 'bool' for k in key_items if k['type'] == 'bool'), False)

                if data_items is not None:
                    has_bool_data = next((d['type'] == 'bool' for d in data_items if d['type'] == 'bool'), False)

                if has_bool_key or has_bool_data:

                    loop = sf.get_loop(lp_category)

                    if has_bool_key:

                        for itName in [k['name'] for k in key_items if k['type'] == 'bool']:

                            if itName in loop.tags:

                                itCol = loop.tags.index(itName)

                                for row in loop:

                                    val = row[itCol]

                                    if val in EMPTY_VALUE:
                                        continue

                                    if val.lower() in TRUE_VALUE:
                                        row[itCol] = yes_value
                                    else:
                                        row[itCol] = no_value

                    if has_bool_data:

                        for itName in [d['name'] for d in data_items if d['type'] == 'bool']:

                            if itName in loop.tags:

                                itCol = loop.tags.index(itName)

                                for row in loop:

                                    val = row[itCol]

                                    if val in EMPTY_VALUE:
                                        continue

                                    if val.lower() in TRUE_VALUE:
                                        row[itCol] = yes_value
                                    else:
                                        row[itCol] = no_value

        return True

    def __resetBoolValueInAuxLoop(self) -> bool:
        """ Reset bool values in auxiliary loops depending on file type.
        """

        if not self.__reg.combined_mode:
            return True

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_type = input_source_dic['file_type']

        yes_value = 'true' if file_type == 'nef' else 'yes'
        no_value = 'false' if file_type == 'nef' else 'no'

        if input_source_dic['content_subtype'] is None:
            return False

        for content_subtype in input_source_dic['content_subtype']:

            if content_subtype in ('entry_info', 'entity'):
                continue

            sf_category = SF_CATEGORIES[file_type][content_subtype]

            for sf in self.__reg.star_data[0].get_saveframes_by_category(sf_category):

                for loop in sf.loops:

                    lp_category = loop.category

                    if lp_category is None:
                        continue

                    # main content of loop has been processed in __resetBoolValueInLoop()
                    if lp_category in LP_CATEGORIES[file_type][content_subtype]:
                        continue

                    if AUX_LP_CATEGORIES[file_type][content_subtype] is None:
                        continue

                    if lp_category in AUX_LP_CATEGORIES[file_type][content_subtype]:

                        key_items = self.__reg.aux_key_items[file_type][content_subtype][lp_category]
                        data_items = self.__reg.aux_data_items[file_type][content_subtype][lp_category]

                        has_bool_key = has_bool_data = False

                        if key_items is not None:
                            has_bool_key = next((k['type'] == 'bool' for k in key_items if k['type'] == 'bool'), False)

                        if data_items is not None:
                            has_bool_data = next((d['type'] == 'bool' for d in data_items if d['type'] == 'bool'), False)

                        if has_bool_key or has_bool_data:

                            _loop = sf.get_loop(lp_category)

                            if has_bool_key:

                                for itName in [k['name'] for k in key_items if k['type'] == 'bool']:

                                    if itName in _loop.tags:

                                        itCol = _loop.tags.index(itName)

                                        for row in _loop:

                                            val = row[itCol]

                                            if val in EMPTY_VALUE:
                                                continue

                                            if val.lower() in TRUE_VALUE:
                                                row[itCol] = yes_value
                                            else:
                                                row[itCol] = no_value

                            if has_bool_data:

                                for itName in [d['name'] for d in data_items if d['type'] == 'bool']:

                                    if itName in _loop.tags:

                                        itCol = _loop.tags.index(itName)

                                        for row in _loop:

                                            val = row[itCol]

                                            if val in EMPTY_VALUE:
                                                continue

                                            if val.lower() in TRUE_VALUE:
                                                row[itCol] = yes_value
                                            else:
                                                row[itCol] = no_value

        return True

    def __appendParentSfTag(self) -> bool:
        """ Append parent tag of saveframe if not exists.
        """

        if not self.__reg.combined_mode:
            return True

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_name = input_source_dic['file_name']
        file_type = input_source_dic['file_type']

        if file_type == 'nef':
            return True

        if input_source_dic['content_subtype'] is None:
            return False

        try:

            for content_subtype in input_source_dic['content_subtype']:

                if content_subtype in ('entry_info', 'entity'):
                    continue

                sf_category = SF_CATEGORIES[file_type][content_subtype]
                lp_category = LP_CATEGORIES[file_type][content_subtype]

                data_items = DATA_ITEMS[file_type][content_subtype]

                list_id_tag_in_lp = None

                if data_items is not None:
                    list_id_tag_in_lp = next((d for d in data_items if d['type'] == 'pointer-index'), None)

                if list_id_tag_in_lp is not None:

                    for sf in self.__reg.star_data[0].get_saveframes_by_category(sf_category):
                        sf_framecode = get_first_sf_tag(sf, 'Sf_framecode')

                        if not any(True for loop in sf.loops if loop.category == lp_category):
                            continue

                        warn_desc = self.__reg.report.warning.getDescription('duplicated_index', file_name, sf_framecode)

                        if (warn_desc is not None)\
                           and warn_desc.split(' ')[0] == f"{SF_TAG_PREFIXES[file_type][content_subtype].lstrip('_')}.ID":
                            continue

                        loop = sf.get_loop(lp_category)

                        itName = list_id_tag_in_lp['name']

                        if itName in loop.tags:

                            itCol = loop.tags.index(itName)

                            list_ids = []

                            for row in loop:

                                val = row[itCol]

                                if val in EMPTY_VALUE:
                                    continue

                                list_ids.append(val)

                            if not list_ids:
                                # No parent identity can be inferred from an empty/all-missing column.
                                # Report the source defect instead of indexing an empty Counter.
                                self.__reg.report.error.appendDescription('missing_data',
                                                                          {'file_name': file_name,
                                                                           'sf_framecode': sf_framecode,
                                                                           'category': lp_category,
                                                                           'description': f"No defined {itName} values; "
                                                                                          "cannot infer saveframe ID."})
                                return False

                            list_id = collections.Counter(list_ids).most_common()[0][0]

                            tagNames = [t[0] for t in sf.tags]

                            if 'ID' in tagNames:

                                itCol = tagNames.index('ID')

                                sf.tags[itCol][1] = list_id

                            else:

                                sf_tag = f'_{sf_category}.ID'
                                warn = WARN_TEMPLATE_FOR_MISSING_MANDATORY_SF_TAG % (sf_tag, file_type.upper())

                                if self.__reg.check_mandatory_tag and self.__reg.nefT.is_mandatory_tag(sf_tag, file_type):

                                    if self.__reg.rescue_mode:
                                        self.__reg.report.warning.appendDescription('missing_data',
                                                                                    {'file_name': file_name,
                                                                                     'sf_framecode': sf_framecode,
                                                                                     'category': sf_category, 'description': warn})

                                        if self.__reg.verbose:
                                            self.__reg.log.write(f"+{self.__class_name__}.__appendParentSfTag() "
                                                                 f"++ Warning  - {warn}\n")

                                sf.add_tag('ID', list_id)

            return True

        except ValueError as e:

            self.__reg.report.error.appendDescription('internal_error',
                                                      f"+{self.__class_name__}.__appendParentSfTag() ++ Error  - " + str(e))

            if self.__reg.verbose:
                self.__reg.log.write(f"+{self.__class_name__}.__appendParentSfTag() ++ Error  - {str(e)}\n")

            return False

    def __addUnnamedEntryId(self) -> bool:
        """ Add UNNAMED entry id.
        """

        if not self.__reg.combined_mode:
            return True

        if len(self.__reg.star_data) == 0 or self.__reg.star_data[0] is None:
            return False

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_type = input_source_dic['file_type']

        if self.__reg.star_data_type[0] == 'Entry':
            if self.__reg.bmrb_only and self.__reg.internal_mode and self.__reg.bmrb_id is not None:
                self.__reg.star_data[0].entry_id = self.__reg.bmrb_id
            else:
                self.__reg.star_data[0].entry_id = f'nef_{self.__reg.entry_id.lower()}'

        if file_type == 'nmr-star':
            self.__reg.c2S.set_entry_id(self.__reg.star_data[0], self.__reg.entry_id)

        self.__reg.dpR.sortCsLoop()

        if file_type == 'nef':
            return True

        # if self.__reg.c2S.set_entry_id(self.__reg.star_data[0], self.__reg.entry_id):
        #     self.__depositNmrData()

        return True

    def __depositNmrData(self) -> bool:
        """ Deposit next NMR unified data file.
        """

        if not self.__reg.combined_mode:

            if (self.__reg.bmrb_only or self.__reg.internal_mode) and self.__reg.dstPath is not None:
                master_entry = self.__reg.c2S.normalize(self.__reg.star_data[0])

                self.__reg.dpR.remediateSpectralPeakListSaveframe(master_entry)

                master_entry.write_to_file(self.__reg.dstPath,
                                           show_comments=(self.__reg.bmrb_only and self.__reg.internal_mode),
                                           skip_empty_loops=True, skip_empty_tags=False)

            return True

        if self.__reg.dstPath is None:

            if not self.__reg.op.endswith('consistency-check'):

                err = "Not found destination file path."

                self.__reg.report.error.appendDescription('internal_error',
                                                          f"+{self.__class_name__}.__depositNmrData() ++ Error  - " + err)

                if self.__reg.verbose:
                    self.__reg.log.write(f"+{self.__class_name__}.__depositNmrData() ++ Error  - {err}\n")

            return False

        if self.__reg.dstPath == self.__reg.srcPath and self.__reg.release_mode:
            return True

        if len(self.__reg.star_data) == 0 or self.__reg.star_data[0] is None:
            return False

        master_entry = self.__reg.c2S.normalize(self.__reg.star_data[0])

        if not self.__reg.submission_mode and not self.__reg.annotation_mode or self.__reg.dstPath != self.__reg.srcPath:
            master_entry.write_to_file(self.__reg.dstPath, show_comments=(self.__reg.bmrb_only and self.__reg.internal_mode),
                                       skip_empty_loops=True, skip_empty_tags=False)

        if self.__reg.op in ('nmr-str2str-deposit', 'nmr-str2cif-deposit', 'nmr-str2cif-annotate') and self.__reg.remediation_mode:

            dir_path = os.path.dirname(self.__reg.dstPath) if self.__reg.dirPath is None else self.__reg.dirPath

            rem_dir = os.path.join(dir_path, 'remediation')

            try:

                if not os.path.isdir(rem_dir):
                    os.makedirs(rem_dir)

                nmr_file_name = os.path.basename(self.__reg.dstPath)

                if nmr_file_name.endswith('_nmr_data.str'):
                    nmr_file_link = os.path.join(rem_dir, nmr_file_name)

                    if os.path.exists(nmr_file_link):
                        os.remove(nmr_file_link)

                    os.symlink(self.__reg.dstPath, nmr_file_link)

            except OSError:
                pass

        if 'nef' not in self.__reg.op and ('deposit' in self.__reg.op or 'annotate' in self.__reg.op
                                           or 'replace-cs' in self.__reg.op)\
           and NMR_CIF_FILE_PATH_KEY in self.__reg.outputParamDict:

            if self.__reg.cifPath is None or self.__reg.submission_mode:

                if self.__dstPath__ is None:
                    self.__dstPath__ = self.__reg.outputParamDict[NMR_CIF_FILE_PATH_KEY]

                self.__reg.dpR.remediateSpectralPeakListSaveframe(master_entry)

                master_entry.write_to_file(self.__dstPath__,
                                           show_comments=(self.__reg.bmrb_only and self.__reg.internal_mode),
                                           skip_empty_loops=True, skip_empty_tags=False)

            try:

                myIo = IoAdapterPy(False, sys.stderr)
                containerList = myIo.readFile(self.__dstPath__)

                if containerList is not None and len(containerList) > 1:

                    if self.__reg.verbose:
                        self.__reg.log.write(f"Input container list is {[(c.getName(), c.getType()) for c in containerList]!r}\n")

                    eff_block_id = 1  # if len(containerList[0].getObjNameList()) == 0 and not self.__reg.internal_mode else 0
                    abandon_symbolic_labels(containerList)
                    myIo.writeFile(self.__reg.outputParamDict[NMR_CIF_FILE_PATH_KEY], containerList=containerList[eff_block_id:])

            except Exception as e:  # pylint: disable=broad-exception-caught
                self.__reg.log.write(f"+{self.__class_name__}.__depositNmrData() ++ Error  - {str(e)}\n")

        return not self.__reg.report.isError()

    def __depositLegacyNmrData(self) -> bool:
        """ Deposit next NMR legacy data files.
        """

        if self.__reg.combined_mode or self.__reg.dstPath is None:
            return True

        if len(self.__reg.star_data) == 0 or self.__reg.star_data[0] is None or self.__reg.star_data_type[0] != 'Entry':
            return False

        master_entry = self.__reg.star_data[0]

        master_entry.entry_id = f'cs_{self.__reg.entry_id.lower()}'

        self.__reg.c2S.set_entry_id(master_entry, self.__reg.entry_id)

        self.__reg.c2S.normalize_str(master_entry)

        try:

            master_entry.write_to_file(self.__reg.dstPath,
                                       show_comments=(self.__reg.bmrb_only and self.__reg.internal_mode),
                                       skip_empty_loops=True, skip_empty_tags=False)

        except Exception:  # pylint: disable=broad-exception-caught
            return False

        if NMR_CIF_FILE_PATH_KEY in self.__reg.outputParamDict:

            try:

                myIo = IoAdapterPy(False, sys.stderr)
                containerList = myIo.readFile(self.__reg.dstPath)

                if containerList is not None and len(containerList) > 1:

                    if self.__reg.verbose:
                        self.__reg.log.write(f"Input container list is {[(c.getName(), c.getType()) for c in containerList]!r}\n")

                    eff_block_id = 1  # if len(containerList[0].getObjNameList()) == 0 and not self.__reg.internal_mode else 0
                    abandon_symbolic_labels(containerList)
                    myIo.writeFile(self.__reg.outputParamDict[NMR_CIF_FILE_PATH_KEY], containerList=containerList[eff_block_id:])

                    return True

            except Exception as e:  # pylint: disable=broad-exception-caught
                self.__reg.log.write(f"+{self.__class_name__}.__depositLegacyNmrData() ++ Error  - {str(e)}\n")

        return False

    def __mergeLegacyData(self) -> bool:
        """ Merge CS+MR+PK into next NMR unified data files.
        """

        if not self.__reg.dpR.mergeLegacyData():
            return False

        # reset cache dictionaries

        for v in self.__reg.lp_data.values():
            v.clear()

        for v in self.__reg.aux_data.values():
            v.clear()

        for v in self.__reg.sf_tag_data.values():
            v.clear()

        self.__testDataConsistencyInLoop()
        self.__testDataConsistencyInAuxLoop()
        self.__testSfTagConsistency()

        return True

    def __updateConstraintStats(self) -> bool:
        """ Update _Constraint_stat_list saveframe.
        """

        return self.__reg.dpR.updateConstraintStats()

    def __detectSimpleDistanceRestraint(self) -> bool:
        """ Detect simple distance restraints.
        """

        if self.__reg.dstPath is None:
            return True

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        file_type = input_source_dic['file_type']

        if file_type == 'nef':
            return True

        if len(self.__reg.star_data) == 0 or not isinstance(self.__reg.star_data[0], pynmrstar.Entry):
            return False

        master_entry = self.__reg.star_data[0]

        if not isinstance(master_entry, pynmrstar.Entry):
            return False

        sf_category = 'constraint_statistics'
        lp_category = '_Constraint_file'

        try:

            sf = master_entry.get_saveframes_by_category(sf_category)[0]

            data_file_name = get_first_sf_tag(sf, 'Data_file_name')
            if len(data_file_name) == 0:
                data_file_name = self.__reg.srcName

            try:
                lp = sf.get_loop(lp_category)
            except KeyError:
                return False

            try:
                block_id_col = lp.tags.index('Block_ID')
            except ValueError:
                return False
            try:
                file_name_col = lp.tags.index('Constraint_filename')
            except ValueError:
                return False

            constraint_type_col = lp.tags.index('Constraint_type')
            constraint_subtype_col = lp.tags.index('Constraint_subtype')
            constraint_subsubtype_col = lp.tags.index('Constraint_subsubtype')

            dist_rows = [row for row in lp if row[constraint_type_col] == 'distance']

            subtypes_not_derived_from_noes = ('paramagnetic relaxation',
                                              'photo cidnp',
                                              'chemical shift perturbation',
                                              'mutation',
                                              'symmetry',
                                              'metal coordination',
                                              'diselenide bond',
                                              'disulfide bond',
                                              'hydrogen bond',
                                              'covalent bond')

            if len(dist_rows) == 0\
               or any(True for row in dist_rows
                      if row[constraint_subtype_col] not in subtypes_not_derived_from_noes
                      and row[constraint_subsubtype_col] == 'simple'):
                return True

            content_subtype = 'dist_restraint'

            sf_category = SF_CATEGORIES[file_type][content_subtype]
            lp_category = LP_CATEGORIES[file_type][content_subtype]

            if not any(True for row in dist_rows
                       if row[constraint_subtype_col] not in subtypes_not_derived_from_noes):

                subtypes = ','.join([row[constraint_subtype_col] for row in dist_rows])

                warn = f"There is no unique distance restraints derived from NOE/ROE experiment, except for {subtypes}. "\
                       "The wwPDB NMR Validation Task Force highly recommends the submission of unambiguous distance restraints "\
                       "used for the structure determination."

                self.__reg.report.warning.appendDescription('missing_content',
                                                            {'file_name': data_file_name, 'category': lp_category,
                                                             'description': warn})

                if self.__reg.verbose:
                    self.__reg.log.write(f"+{self.__class_name__}.__detectSimpleDistanceRestraint() ++ Warning  - {warn}\n")

                return False

            block_ids = {row[block_id_col]: row[file_name_col] for row in dist_rows
                         if row[constraint_subtype_col] not in subtypes_not_derived_from_noes}

            for block_id in block_ids:
                for sf in master_entry.get_saveframes_by_category(sf_category):
                    if get_first_sf_tag(sf, 'Block_ID') == block_id:

                        try:
                            lp = sf.get_loop(lp_category)
                        except KeyError:
                            continue

                        item_names = ITEM_NAMES_IN_DIST_LOOP[file_type]
                        id_col = lp.tags.index('ID')
                        member_logic_code_col = lp.tags.index('Member_logic_code') if 'Member_logic_code' in lp.tags else -1
                        auth_asym_id_1_col = lp.tags.index('Auth_asym_ID_1')
                        auth_seq_id_1_col = lp.tags.index('Auth_seq_ID_1')
                        auth_asym_id_2_col = lp.tags.index('Auth_asym_ID_2')
                        auth_seq_id_2_col = lp.tags.index('Auth_seq_ID_2')
                        comp_id_1_col = lp.tags.index(item_names['comp_id_1'])
                        comp_id_2_col = lp.tags.index(item_names['comp_id_2'])
                        atom_id_1_col = lp.tags.index(item_names['atom_id_1'])
                        atom_id_2_col = lp.tags.index(item_names['atom_id_2'])

                        for row in lp:
                            if member_logic_code_col != -1 and row[member_logic_code_col] != 'OR':
                                return True

                        prev_id = -1
                        for row in lp:
                            if member_logic_code_col != -1 and row[member_logic_code_col] == 'OR':
                                _id = int(row[id_col])
                                if _id != prev_id:
                                    _atom1 = {'chain_id': row[auth_asym_id_1_col],
                                              'seq_id':
                                              int(row[auth_seq_id_1_col]) if row[auth_seq_id_1_col] not in EMPTY_VALUE else None,
                                              'comp_id': row[comp_id_1_col],
                                              'atom_id': row[atom_id_1_col]}
                                    _atom2 = {'chain_id': row[auth_asym_id_2_col],
                                              'seq_id':
                                              int(row[auth_seq_id_2_col]) if row[auth_seq_id_2_col] not in EMPTY_VALUE else None,
                                              'comp_id': row[comp_id_2_col],
                                              'atom_id': row[atom_id_2_col]}
                                    prev_id = _id
                                    continue
                                atom1 = {'chain_id': row[auth_asym_id_1_col],
                                         'seq_id':
                                         int(row[auth_seq_id_1_col]) if row[auth_seq_id_1_col] not in EMPTY_VALUE else None,
                                         'comp_id': row[comp_id_1_col],
                                         'atom_id': row[atom_id_1_col]}
                                atom2 = {'chain_id': row[auth_asym_id_2_col],
                                         'seq_id':
                                         int(row[auth_seq_id_2_col]) if row[auth_seq_id_2_col] not in EMPTY_VALUE else None,
                                         'comp_id': row[comp_id_2_col],
                                         'atom_id': row[atom_id_2_col]}
                                if not isAmbigAtomSelection([_atom1, atom1], self.__reg.csStat)\
                                   and not isAmbigAtomSelection([_atom2, atom2], self.__reg.csStat):
                                    return True
                                _atom1, _atom2 = atom1, atom2

            for block_id, file_name in block_ids.items():
                for sf in master_entry.get_saveframes_by_category(sf_category):
                    if block_id == get_first_sf_tag(sf, 'Block_ID'):
                        sf_framecode = get_first_sf_tag(sf, 'Sf_framecode')

                        warn = "There is no unique distance restraints derived from NOE/ROE experiment "\
                            "in the set of uploaded restraint file(s). "\
                            "The wwPDB NMR Validation Task Force highly recommends "\
                            "the submission of unambiguous distance restraints used for the structure determination."

                        self.__reg.report.warning.appendDescription('encouragement',
                                                                    {'file_name': file_name, 'sf_framecode': sf_framecode,
                                                                     'category': lp_category, 'description': warn})

                        if self.__reg.verbose:
                            self.__reg.log.write(f"+{self.__class_name__}.__detectSimpleDistanceRestraint() ++ Warning  - {warn}\n")

            has_upper_bound = False
            lower_bound_sf_names = []
            for block_id in block_ids:
                for sf in master_entry.get_saveframes_by_category(sf_category):
                    if block_id == get_first_sf_tag(sf, 'Block_ID'):
                        potential_type = get_first_sf_tag(sf, 'Potential_type')
                        if potential_type not in EMPTY_VALUE:
                            if 'upper-bound' in potential_type\
                               or 'square-well' in potential_type\
                               or 'log-harmonic' in potential_type\
                               or potential_type == 'parabolic':
                                has_upper_bound = True
                            if 'lower-bound' in potential_type:
                                lower_bound_sf_names.append(f'{sf.name!r}')

            if not has_upper_bound and len(lower_bound_sf_names) > 0:
                if len(lower_bound_sf_names) == 0:
                    lower_bound_sf_names = lower_bound_sf_names[0]

                err = "There are no attractive distance restraints expressed by upper limit potential, "\
                      f"except for repulsive distance restraints, {lower_bound_sf_names}. "\
                      "It means that your restraints do not have sufficient power to determine the structure. "\
                      "Please re-upload all restraints used for the structure determination."

                self.__reg.report.error.appendDescription('missing_mandatory_content',
                                                          {'file_name': data_file_name, 'category': lp_category,
                                                           'description': err})

                if self.__reg.verbose:
                    self.__reg.log.write(f"+{self.__class_name__}.__detectSimpleDistanceRestraint() ++ Error  - {err}\n")

            return False

        except IndexError:
            return True

    def __initializeDpReportForNext(self) -> bool:
        """ Initialize NMR data processing report using the next version of NMR unified data.
        """

        return self.__initializeDpReport(srcPath=self.__reg.dstPath, calcOutputStats=True)

    def __validateInputSourceForNext(self) -> bool:
        """ Validate the next version of NMR unified data as primary input source.
        """

        return self.__validateInputSource(srcPath=self.__reg.dstPath)

    def __translateNef2Str(self) -> bool:
        """ Translate NEF to NMR-STAR V3.2 file.
        """

        if not self.__reg.combined_mode:
            return True

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        if self.__reg.dstPath is None:
            raise KeyError(f"+{self.__class_name__}.__translateNef2Str() "
                           "++ Error  - Could not find destination path as input NEF file for NefTranslator.")

        file_name = os.path.basename(self.__reg.dstPath)
        file_type = input_source_dic['file_type']

        if NEXT_STAR_FILE_PATH_KEY not in self.__reg.outputParamDict:
            raise KeyError(f"+{self.__class_name__}.__translateNef2Str() "
                           f"++ Error  - Could not find {NEXT_STAR_FILE_PATH_KEY!r} output parameter.")

        fPath = self.__reg.outputParamDict[NEXT_STAR_FILE_PATH_KEY]

        try:

            is_valid, message = self.__reg.nefT.nef_to_nmrstar(self.__reg.dstPath, fPath,
                                                               report=self.__reg.report,
                                                               leave_unmatched=self.__reg.leave_intl_note)

            if self.__reg.release_mode and self.__tmpPath is not None:
                os.remove(self.__tmpPath)
                self.__tmpPath = None

        except Exception as e:  # pylint: disable=broad-exception-caught

            err = f"{file_name!r} is not compliant with the {READABLE_FILE_TYPE[file_type]} dictionary."

            if 'No such file or directory' not in str(e):
                err += ' ' + re.sub('not in list', 'unknown item.', str(e))

            if not self.__reg.report.isError():
                self.__reg.report.error.appendDescription('format_issue',
                                                          {'file_name': file_name, 'description': err})

            self.__reg.log.write(f"+{self.__class_name__}.__translateNef2Str() ++ Error  - "
                                 f"{file_name} {err}\n")

            if os.path.exists(fPath):
                os.remove(fPath)

            return False

        if is_valid:

            if 'deposit' in self.__reg.op and NMR_CIF_FILE_PATH_KEY in self.__reg.outputParamDict:

                try:

                    myIo = IoAdapterPy(False, sys.stderr)
                    containerList = myIo.readFile(fPath)

                    if containerList is not None and len(containerList) > 1:

                        if self.__reg.verbose:
                            self.__reg.log.write("Input container list is "
                                                 f"{[(c.getName(), c.getType()) for c in containerList]!r}\n")

                        eff_block_id = 1  # if len(containerList[0].getObjNameList()) == 0 and not self.__reg.internal_mode else 0
                        abandon_symbolic_labels(containerList)
                        myIo.writeFile(self.__reg.outputParamDict[NMR_CIF_FILE_PATH_KEY],
                                       containerList=containerList[eff_block_id:])

                except Exception as e:  # pylint: disable=broad-exception-caught
                    self.__reg.log.write(f"+{self.__class_name__}.__translateNef2Str() ++ Error  - {str(e)}\n")

            return True

        err = f"{file_name!r} is not compliant with the {READABLE_FILE_TYPE[file_type]} dictionary."

        if len(message['error']) > 0:
            for err_message in message['error']:
                if 'No such file or directory' not in err_message:
                    err += ' ' + re.sub('not in list', 'unknown item.', err_message)

        if not self.__reg.report.isError():
            self.__reg.report.error.appendDescription('format_issue',
                                                      {'file_name': file_name, 'description': err})

        self.__reg.log.write(f"+{self.__class_name__}.__translateNef2Str() ++ Error  - "
                             f"{file_name} {err}\n")

        if os.path.exists(fPath):
            os.remove(fPath)

        return False

    def __initResourceForNef2Str(self) -> bool:
        """ Initialize resources for the translated NMR-STAR V3.2 file.
        """

        self.__reg.rescue_mode = False

        self.__reg.report_prev = None

        try:

            self.__reg.srcPath = self.__reg.outputParamDict[NEXT_STAR_FILE_PATH_KEY]
            self.__reg.dstPath = self.__reg.srcPath
            self.__logPath = self.__reg.outputParamDict.get(REPORT_FILE_PATH_KEY)
            if self.__logPath is not None:
                self.addInput(REPORT_FILE_PATH_KEY, self.__logPath, type='file')
            self.__reg.op = 'nmr-str-consistency-check'

            # reset cache dictionaries

            for v in self.__reg.lp_data.values():
                v.clear()

            for v in self.__reg.aux_data.values():
                v.clear()

            for v in self.__reg.sf_tag_data.values():
                v.clear()

            return True

        except Exception as e:
            raise KeyError(f"+{self.__class_name__}.__initReousrceForNef2Str() ++ Error  - "
                           f"Could not find {NEXT_STAR_FILE_PATH_KEY!r} or {REPORT_FILE_PATH_KEY!r} output parameter.") from e

        return False

    def __translateStr2Nef(self) -> bool:
        """ Translate NMR-STAR V3.2 to NEF file.
        """

        if not self.__reg.combined_mode:
            return True

        input_source = self.__reg.report.input_sources[0]
        input_source_dic = input_source.get()

        if self.__reg.dstPath is None:
            raise KeyError(f"+{self.__class_name__}.__translateStr2Nef() "
                           "++ Error  - Could not find destination path as input NMR-STAR file for NefTranslator.")

        file_name = os.path.basename(self.__reg.dstPath)
        file_type = input_source_dic['file_type']

        if NEXT_NEF_FILE_PATH_KEY not in self.__reg.outputParamDict:
            raise KeyError(f"+{self.__class_name__}.__translateStr2Nef() "
                           f"++ Error  - Could not find {NEXT_NEF_FILE_PATH_KEY!r} output parameter.")

        fPath = self.__reg.outputParamDict[NEXT_NEF_FILE_PATH_KEY]

        try:

            is_valid, message = self.__reg.nefT.nmrstar_to_nef(self.__reg.dstPath, fPath, report=self.__reg.report)

            if self.__reg.release_mode and self.__tmpPath is not None:
                os.remove(self.__tmpPath)
                self.__tmpPath = None

        except Exception as e:  # pylint: disable=broad-exception-caught

            err = f"{file_name!r} is not compliant with the {READABLE_FILE_TYPE[file_type]} dictionary."

            if 'No such file or directory' not in str(e):
                err += ' ' + re.sub('not in list', 'unknown item.', str(e))

            if not self.__reg.report.isError():
                self.__reg.report.error.appendDescription('format_issue',
                                                          {'file_name': file_name, 'description': err})

            self.__reg.log.write(f"+{self.__class_name__}.__translateStr2Nef() ++ Error  - "
                                 f"{file_name} {err}\n")

            if os.path.exists(fPath):
                os.remove(fPath)

            return False

        if is_valid:
            return True

        err = f"{file_name!r} is not compliant with the {READABLE_FILE_TYPE[file_type]} dictionary."

        if len(message['error']) > 0:
            for err_message in message['error']:
                if 'No such file or directory' not in err_message:
                    err += ' ' + re.sub('not in list', 'unknown item.', err_message)

        if not self.__reg.report.isError():
            self.__reg.report.error.appendDescription('format_issue',
                                                      {'file_name': file_name, 'description': err})

        self.__reg.log.write(f"+{self.__class_name__}.__translateStr2Nef() ++ Error  - "
                             f"{file_name} {err}\n")

        if os.path.exists(fPath):
            os.remove(fPath)

        return False

    def __initResourceForStr2Nef(self) -> bool:
        """ Initialize resources for the translated NEF file.
        """

        self.__reg.rescue_mode = False

        self.__reg.report_prev = None

        try:

            self.__reg.srcPath = self.__reg.outputParamDict[NEXT_NEF_FILE_PATH_KEY]
            self.__reg.dstPath = self.__reg.srcPath
            self.__logPath = self.__reg.outputParamDict.get(REPORT_FILE_PATH_KEY)
            if self.__logPath is not None:
                self.addInput(REPORT_FILE_PATH_KEY, self.__logPath, type='file')
            self.__reg.op = 'nmr-nef-consistency-check'

            # reset cache dictionaries

            for v in self.__reg.lp_data.values():
                v.clear()

            for v in self.__reg.aux_data.values():
                v.clear()

            for v in self.__reg.sf_tag_data.values():
                v.clear()

            return True

        except Exception as e:
            raise KeyError(f"+{self.__class_name__}.__initReousrceForStr2Nef() ++ Error  - "
                           f"Could not find {NEXT_NEF_FILE_PATH_KEY!r} or {REPORT_FILE_PATH_KEY!r} output parameter.") from e

        return False

    def __parseNmrIf(self) -> bool:
        """ Parse NMRIF file.
        """

        if NMRIF_FILE_PATH_KEY not in self.__reg.inputParamDict:

            err = f"No such {self.__reg.inputParamDict[NMRIF_FILE_PATH_KEY]!r} file."

            self.__reg.report.error.appendDescription('internal_error',
                                                      f"+{self.__class_name__}.__parseNmrIf() ++ Error  - " + err)

            if self.__reg.verbose:
                self.__reg.log.write(f"+{self.__class_name__}.__parseNmrIf() ++ Error  - {err}\n")

            return False

        if self.__nmrIfR is not None:
            return True

        fPath = self.__reg.inputParamDict[NMRIF_FILE_PATH_KEY]

        self.__nmrIfR = CifReader(self.__reg.verbose, self.__reg.log)

        try:

            if self.__nmrIfR.parse(fPath):
                return True

        except IOError:
            pass

        return False

    def __mergeCoordAsNmrIf(self) -> bool:
        """ Merge NMRIF metadata of the coordinates (DAOTHER-8905: NMR data remediation Phase 2).
        """

        if len(self.__reg.star_data) == 0 or self.__reg.star_data[0] is None or self.__reg.star_data_type[0] != 'Entry':
            return False

        if not self.__reg.internal_mode or not self.__reg.remediation_mode or self.__reg.report.getInputSourceIdOfCoord() < 0:
            return False

        master_entry = self.__reg.star_data[0]

        self.__reg.sf_category_list, self.__reg.lp_category_list = self.__reg.nefT.get_inventory_list(master_entry)

        if not self.__reg.bmrb_only and self.__reg.cR.hasCategory('entry'):
            entry = self.__reg.cR.getDictList('entry')

            if len(entry) > 0 and 'id' in entry[0]:
                self.__reg.entry_id = entry[0]['id'].strip().replace(' ', '_')

        ann = OneDepAnnTasks(self.__reg)

        return ann.merge(master_entry, self.__reg.cR, self.__reg.bmrb_only and self.__reg.internal_mode)

    def __mergeNmrIf(self) -> bool:
        """ Merge NMRIF metadata.
        """

        if len(self.__reg.star_data) == 0 or self.__reg.star_data[0] is None or self.__reg.star_data_type[0] != 'Entry':
            return False

        if self.__nmrIfR is None or not self.__reg.submission_mode:
            return False

        master_entry = self.__reg.star_data[0]

        self.__reg.sf_category_list, self.__reg.lp_category_list = self.__reg.nefT.get_inventory_list(master_entry)

        if not self.__reg.bmrb_only and self.__nmrIfR.hasCategory('entry'):
            entry = self.__nmrIfR.getDictList('entry')

            if len(entry) > 0 and 'id' in entry[0]:
                self.__reg.entry_id = entry[0]['id'].strip().replace(' ', '_')

        ann = OneDepAnnTasks(self.__reg)

        return ann.merge(master_entry, self.__nmrIfR, self.__reg.bmrb_only and self.__reg.internal_mode)

    def __performBmrbAnnTasks(self) -> bool:
        """ Perform a series of standalone BMRB annotation tasks.
        """

        if len(self.__reg.star_data) == 0 or self.__reg.star_data[0] is None or self.__reg.star_data_type[0] != 'Entry':
            return False

        # if not self.__reg.submission_mode and not self.__reg.internal_mode:
        #     return False

        master_entry = self.__reg.star_data[0]

        self.__reg.sf_category_list, self.__reg.lp_category_list = self.__reg.nefT.get_inventory_list(master_entry)

        if not self.__reg.bmrb_only and self.__nmrIfR is not None and self.__nmrIfR.hasCategory('entry'):
            entry = self.__nmrIfR.getDictList('entry')

            if len(entry) > 0 and 'id' in entry[0]:
                self.__reg.entry_id = entry[0]['id'].strip().replace(' ', '_')

        ann = BmrbAnnTasks(self.__reg)

        if self.__reg.report.getInputSourceIdOfCoord() >= 0 and self.__reg.cR.hasCategory('database_2'):

            sf_category = 'entry_information'

            database_code = self.__reg.cR.getDictListWithFilter('database_2',
                                                                [{'name': 'database_code', 'type': 'str'}],
                                                                [{'name': 'database_id', 'type': 'str', 'value': 'BMRB'}])

            if len(database_code) > 0:
                derived_entry_id = database_code[0]['database_code']

                if not self.__reg.internal_mode:
                    derived_entry_title = None

                    if sf_category in self.__reg.sf_category_list:
                        sf = master_entry.get_saveframes_by_category(sf_category)[0]
                        derived_entry_title = get_first_sf_tag(sf, 'Title', None)

                    ann.setProvenanceInfo(derived_entry_id, derived_entry_title)

                elif not self.__reg.bmrb_only:

                    if sf_category in self.__reg.sf_category_list:
                        sf = master_entry.get_saveframes_by_category(sf_category)[0]
                        derived_entry_title = get_first_sf_tag(sf, 'Title', None)

                        if derived_entry_title in EMPTY_VALUE:
                            dir_path = self.__reg.cR.getDirPath()

                            if os.path.exists(os.path.join(dir_path, f'bmr{derived_entry_id}_3.str')):

                                is_done, _, star_data =\
                                    self.__reg.nefT.read_input_file(os.path.join(dir_path, f'bmr{derived_entry_id}_3.str'))

                                if is_done:

                                    try:

                                        _sf = star_data.get_saveframes_by_category(sf_category)[0]
                                        derived_entry_title = get_first_sf_tag(_sf, 'Title', None)

                                    except IndexError:
                                        pass

                            if derived_entry_title not in EMPTY_VALUE:
                                set_sf_tag(sf, 'Title', derived_entry_title)

                        lp_category = '_Related_entries'

                        try:

                            lp = sf.get_loop(lp_category)

                            tags = ['Database_name', 'Database_accession_code']

                            if set(tags) & set(lp.tags) == set(tags):
                                related_entries = lp.get_tag(tags)

                                has_provenance = False
                                for idx, row in enumerate(related_entries):
                                    if row == ['BMRB', derived_entry_id]:
                                        has_provenance = True
                                        if derived_entry_title not in EMPTY_VALUE:
                                            lp.data[idx][lp.tags.index('Relationship')] = derived_entry_title
                                        break

                                if not has_provenance:
                                    row = [None] * len(lp.tags)
                                    row[lp.tags.index('Database_name')] = 'BMRB'
                                    row[lp.tags.index('Database_accession_code')] = derived_entry_id
                                    if derived_entry_title not in EMPTY_VALUE:
                                        row[lp.tags.index('Relationship')] = derived_entry_title
                                    row[lp.tags.index('Entry_ID')] = self.__reg.entry_id

                                    lp.add_data(row)

                                    lp.sort_rows(['Database_name', 'Database_accession_code'])

                        except KeyError:
                            items = ['Database_name', 'Database_accession_code', 'Relationship', 'Entry_ID']

                            lp = pynmrstar.Loop.from_scratch(lp_category)

                            tags = [f'{lp_category}.{item}' for item in items]

                            lp.add_tag(tags)

                            lp.add_data(['BMRB', derived_entry_id, derived_entry_title, self.__reg.entry_id])

                            database_code = self.__reg.cR.getDictListWithFilter('database_2',
                                                                                [{'name': 'database_code', 'type': 'str'},
                                                                                 {'name': 'database_id', 'type': 'str'}])

                            if len(database_code) > 0:
                                for d in database_code:
                                    if d['database_code'] in EMPTY_VALUE or d['database_id'] in EMPTY_VALUE\
                                       or d['database_id'] in ('PDB', 'BMRB', 'WWPDB'):
                                        continue
                                    lp.add_data([d['database_id'], d['database_code'], None, self.__reg.entry_id])

                            if len(lp) > 1:
                                lp.sort_rows(['Database_name', 'Database_accession_code'])

                            sf.add_loop(lp)

        if self.__reg.op == 'nmr-cs-mr-merge' and self.__reg.bmrb_only\
                and CS_FILE_PATH_LIST_KEY in self.__reg.inputParamDict:
            cs = self.__reg.inputParamDict[CS_FILE_PATH_LIST_KEY][0]

            if isinstance(cs, str):
                cs_path = cs
            else:
                cs_path = cs['file_name']

            self.__reg.c2S.set_entry_id(master_entry, self.__reg.bmrb_id)
            self.__reg.c2S.normalize_str(master_entry)

            master_entry.write_to_file(cs_path, show_comments=True, skip_empty_loops=True, skip_empty_tags=False)

        is_done = ann.perform(master_entry)

        if is_done:

            if self.__reg.bmrb_only:
                input_source = self.__reg.report.input_sources[0]
                input_source_dic = input_source.get()

                file_type = input_source_dic['file_type']

                self.__reg.sf_category_list, self.__reg.lp_category_list = self.__reg.nefT.get_inventory_list(master_entry)

                lp_counts = {t: 0 for t in NMR_CONTENT_SUBTYPES}
                for lp_category in self.__reg.lp_category_list:
                    if lp_category in LP_CATEGORIES[file_type].values():
                        lp_counts[[k for k, v in LP_CATEGORIES[file_type].items() if v == lp_category][0]] += 1

                content_subtypes = {k: lp_counts[k] for k in lp_counts if lp_counts[k] > 0}

                input_source.setItemValue('content_subtype', content_subtypes)

            self.__depositNmrData()

        return is_done

    def __discardPeakListRemediation(self) -> bool:
        """ Discard remediated spectral peak list (NMR data remediation, Phase 2) in OneDep environment
            @note: This ridiculous reverse implementation is for OneDep only
        """

        if self.__reg.native_combined or not self.__reg.merge_any_pk_as_is:
            return True

        if len(self.__reg.star_data) == 0 or self.__reg.star_data[0] is None or self.__reg.star_data_type[0] != 'Entry':
            return False

        if AR_FILE_PATH_LIST_KEY not in self.__reg.inputParamDict:
            return True

        src_id = self.__reg.report.getInputSourceIdOfCoord()

        if src_id < 0:
            return False

        cif_input_source = self.__reg.report.input_sources[src_id]
        cif_input_source_dic = cif_input_source.get()

        has_poly_seq = has_key_value(cif_input_source_dic, 'polymer_sequence')

        if not has_poly_seq:
            return False

        file_type = 'nmr-star'
        content_subtype = 'spectral_peak'

        sf_category = SF_CATEGORIES[file_type][content_subtype]

        if sf_category not in self.__reg.sf_category_list:
            return True

        modified = False

        for sf in self.__reg.star_data[0].get_saveframes_by_category(sf_category):

            text_data = get_first_sf_tag(sf, 'Text_data')

            if text_data not in EMPTY_VALUE:
                continue

            file_name = get_first_sf_tag(sf, 'Data_file_name')

            if file_name in EMPTY_VALUE:
                continue

            for ar in self.__reg.inputParamDict[AR_FILE_PATH_LIST_KEY]:
                file_type = ar['file_type']

                if file_type is None or not file_type.startswith('nm-pea'):
                    continue

                file_path = ar['file_name']
                original_file_name = ar['original_file_name'] if 'original_file_name' in ar else None

                if file_name in (os.path.basename(file_path), original_file_name):
                    data_format = get_peak_list_format(file_path)

                    if data_format is not None:

                        for lp in sf.loops:
                            if lp.category in ('_Peak_row_format', '_Peak', '_Peak_general_char',
                                               '_Peak_char', '_Assigned_peak_chem_shift'):
                                del sf[lp]  # What a waste!

                        set_sf_tag(sf, 'Text_data_format', data_format)

                        with open(file_path, 'r', encoding='ascii', errors='ignore') as ifh:
                            set_sf_tag(sf, 'Text_data', ifh.read())

                        modified = True

                    break

        if modified:
            master_entry = self.__reg.c2S.normalize(self.__reg.star_data[0])

            master_entry.write_to_file(self.__reg.dstPath,
                                       show_comments=(self.__reg.bmrb_only and self.__reg.internal_mode),
                                       skip_empty_loops=True, skip_empty_tags=False)

        return True

    def __replaceCsSf(self) -> bool:
        """ Replace assigned chemical shift saveframe(s).
        """

        if len(self.__reg.star_data) == 0 or self.__reg.star_data[0] is None or self.__reg.star_data_type[0] != 'Entry':
            return False

        master_entry = self.__reg.star_data[0]

        has_chem_shift = False

        file_type = 'nmr-star'
        content_subtype = 'chem_shift'

        cs_files = []

        for fileListId in range(self.__reg.cs_file_path_list_len):

            fileListId += 1

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_type = input_source_dic['file_type']
            file_name = input_source_dic['file_name']
            cs_files.append(f'{file_name!r}')

            if file_type != 'nmr-star':
                continue

            if input_source_dic['content_subtype'] is None:
                continue

            if content_subtype not in input_source_dic['content_subtype']:
                continue

            if self.__reg.star_data_type[fileListId] not in ('Saveframe', 'Entry'):
                continue

            has_chem_shift = True

        if not has_chem_shift:
            err = f"There is no input assigned chemical shifts in {', '.join(cs_files)} file(s)."

            self.__reg.report.error.appendDescription('internal_error',
                                                      f"+{self.__class_name__}.__replaceCsSf() ++ Error  - " + err)

            if self.__reg.verbose:
                self.__reg.log.write(f"+{self.__class_name__}.__replaceCsSf() ++ Error  - {err}\n")

            return False

        sf_category = SF_CATEGORIES[file_type][content_subtype]
        lp_category = LP_CATEGORIES[file_type][content_subtype]

        removed_sf = []
        removed_sf_framecode = []
        for sf in master_entry.get_saveframes_by_category(sf_category):
            removed_sf.append(sf)
            removed_sf_framecode.append(sf.name)
            master_entry.remove_saveframe(sf.name)

        added_sf_framecode = []

        list_id = 1

        for fileListId in range(self.__reg.cs_file_path_list_len):

            fileListId += 1

            input_source = self.__reg.report.input_sources[fileListId]
            input_source_dic = input_source.get()

            file_type = input_source_dic['file_type']

            if file_type != 'nmr-star':
                continue

            if input_source_dic['content_subtype'] is None:
                continue

            if content_subtype not in input_source_dic['content_subtype']:
                continue

            if self.__reg.star_data_type[fileListId] not in ('Saveframe', 'Entry'):
                continue

            _list_id = list_id

            if isinstance(self.__reg.star_data[fileListId], pynmrstar.Saveframe):
                sf = self.__reg.star_data[fileListId]
                self.__reg.c2S.set_local_sf_id(sf, list_id)

                if list_id - 1 < len(removed_sf):
                    for src_cs_lp in removed_sf[list_id - 1]:
                        if not any(True for lp in sf if lp.category == src_cs_lp.category):
                            sf.add_loop(src_cs_lp)

                added_sf_framecode.append(sf.name)
                master_entry.add_saveframe(sf)

                list_id += 1

            else:

                for sf in self.__reg.star_data[fileListId].get_saveframes_by_category(sf_category):
                    self.__reg.c2S.set_local_sf_id(sf, list_id)

                    if list_id - 1 < len(removed_sf):
                        for src_cs_lp in removed_sf[list_id - 1]:
                            if not any(True for lp in sf if lp.category == src_cs_lp.category):
                                sf.add_loop(src_cs_lp)

                    added_sf_framecode.append(sf.name)
                    master_entry.add_saveframe(sf)

                    list_id += 1

            if list_id == _list_id:

                for lp in self.__reg.star_data[fileListId].get_loops_by_category(lp_category):
                    sf = pynmrstar.Saveframe.from_scratch(f'assigned_chemical_shift_{list_id}', '_Assigned_chem_shift_list')
                    sf.add_tag('Sf_category', sf_category)
                    sf.add_tag('Sf_framecode', f'assigned_chemical_shift_{list_id}')
                    sf.add_tag('Entry_ID', self.__reg.entry_id)
                    sf.add_tag('ID', list_id)
                    sf.add_tag('Data_file_name', input_source_dic['file_name'])
                    sf.add_loop(lp)

                    added_sf_framecode.append(sf.name)
                    master_entry.add_saveframe(sf)

                    list_id += 1

        if list_id > 1:
            content_subtype = 'spectral_peak'

            sf_category = SF_CATEGORIES[file_type][content_subtype]

            for sf in master_entry.get_saveframes_by_category(sf_category):
                cs_framecode = get_first_sf_tag(sf, 'Chemical_shift_list')
                if cs_framecode in removed_sf_framecode:
                    idx = removed_sf_framecode.index(cs_framecode)
                    set_sf_tag(sf, 'Chemical_shift_list', added_sf_framecode[idx] if idx < len(added_sf_framecode) else None)
                else:
                    set_sf_tag(sf, 'Chemical_shift_list', None)

        return list_id > 1
