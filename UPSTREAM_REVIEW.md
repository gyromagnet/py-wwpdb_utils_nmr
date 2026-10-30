# Upstream review and rebase — 2026-10-30

This review compares the previously qualified fork base `4ccd0e2d` with all
advertised branches from both upstream repositories. The previous local head
`4f1d3bea` is preserved under
`refs/archive/pre-upstream-rebase-20261030/fix-peak-parser-exptl-method`;
the old main is also archived. A recovery bundle is retained outside the repo.

## Reviewed branches

| Repository | Branch | Reviewed head | Disposition |
| --- | --- | --- | --- |
| wwPDB | master | `64410c1e8a4004e9ca2e1c6fe21fa014fffce398` | Historical v0.53 release line; not the active BMRB Extract baseline. |
| wwPDB | develop | `c24316fa6e8566880b9ccbc05be201b022a9161f` | Divergent release/integration line; reviewed separately, not stacked wholesale. |
| wwPDB | bmrb-extract_devel | `3f40da35f03e964598967d29b5039b4e8f8fa47c` | Already an ancestor of the selected target. |
| yokochi47 | devel | `9fcbfde31c4e3241b5883ec5e4718be6163e16bb` | Already an ancestor of the selected target. |
| yokochi47 | main | `a450780df1c0c445b5617edfb9ef2ac98af2adb7` | Selected rebase target; integrates both active feature lines. |

The selected target adds 36 commits (19 non-merge changes) after the old base,
changing 51 files. Ten local commits replayed without conflicts;
`git range-diff` marked all ten patches equivalent. No local repair was dropped.
Both upstreams are configured as read-only HTTPS fetch remotes (`wwpdb`,
`yokochi`); the fork remains the push remote.

## Changes affecting the packet work

- AMBER RDC/alignment grammar now accepts `FJ` and `DFAC_AL`; atom-name lexical
  handling also changed. This affects reader support, not the still-missing
  independently verified topology/atom-index mappings in the official corpus.
- CNS/XPLOR grammar changes cover classification names, optional assignment
  punctuation, flags and annotations. Their Python and generated C++ parsers
  changed together. These affect native restraint packets and format detection.
- CCPN, Sparky and VNMR comment lexers no longer require a trailing newline.
  Sparky-derived readers share lexer code, so installed native/Python grammar
  fingerprints must be checked rather than assuming an old binary remains valid.
- `BaseStackedMRParserListener` changes duplex chain/sequence remapping and
  polymer/non-polymer classification. `AlignUtil` adjusts identical-chain lookup.
  Mapping-sensitive cases require fresh evidence after the rebase.
- `NmrDpReport` and `NmrVrptUtility` distinguish unmodeled chemical shifts when
  coordinate unobserved-residue metadata is unavailable. This affects validation
  classification; it does not prove a packet's author-numbering correspondence.
- `BmrbAnnTasks` initializes the chemical-shift-perturbation loop before using it.
  This is separate from the fork's nullable assembly-annotation repair.
- First-aid and input validation include pynmrstar diagnostic text. Exact wording
  can therefore change while the underlying rejection remains correct.

Representative upstream changes:
[grammar/mapping work](https://github.com/yokochi47/py-wwpdb_utils_nmr/commit/0e523766),
[unmodeled shifts](https://github.com/yokochi47/py-wwpdb_utils_nmr/commit/5d2f050b),
[annotation initialization](https://github.com/yokochi47/py-wwpdb_utils_nmr/commit/90daa0e4),
[parser diagnostic text](https://github.com/wwPDB/py-wwpdb_utils_nmr/commit/3f40da35).

## Divergent release branches

Ancestry alone overstates differences: several `wwPDB/develop` cherry-picks have
counterparts on the active feature line. Direct comparison finds identical
`BmrbChemShiftStat.py`, chemical-shift-statistics data and `setup.py` between
`develop` and the selected target, including the DATAQUALITY-2320/TUV work and
accelerator build configuration. The release line's Python 3.11 Azure changes
are packaging/CI differences, not additional scientific parser fixes.

Its remaining source differences include older selector/unknown-type handling,
output-statistics guards, coordinate ensemble processing and grammar versions.
Rebasing the current BMRB Extract fork onto that line, or merging it wholesale,
would replace newer active-line behavior. No additional release-line source
change was identified that should override the selected target for this work.

## Validation boundary

The v8 image and all historical packet snapshots remain unchanged. Qualification
v9 rebuilds native accelerators against the rebased source while retaining the
same pinned science-resource base. `just accelerators` requires every native
parser to be compiled and admitted, rather than silently falling back to Python.
New CNS/XPLOR/AMBER regressions compare native and Python parse trees and reject
malformed controls. Compiler and test tools are local qualification dependencies.

Rebase success and parser agreement do not resolve the outstanding scientific
source/producer/mapping qualification gaps. Fresh checks and packet results are
recorded in `QUALIFICATION.md`; historical results are not promoted automatically.


## Confirmed upstream CNS regression repaired locally

Commit [7969f224](https://github.com/yokochi47/py-wwpdb_utils_nmr/commit/7969f224)
placed the new `CLAS_MODE` declaration before the existing `Class_name_AM` and
`SPACE_AM` lexer rules. In the active 4.10 grammar those averaging rules thereby
became members of the classification mode, leaving averaging unable to tokenize
its class selector and whitespace correctly.

`NOE AVERAGE * R-6 END` has no lexer/parser errors on the v8 Python parser, but
on the newly rebased upstream parser produces a lexer error at the selector and
`missing Class_name_AM at 'R-6'`. This is a syntax control, not a complete
scientific deposition. The repair moves the classification-mode block after the
averaging rules, retaining both existing averaging and the new reserved-word
classification-name support. XPLOR's averaging mode and the archived 4.9 CNS
mode did not have this misplaced declaration and are unchanged.

`just regenerate-cns` rebuilds the Python lexer/parser and native C++ source with
ANTLR 4.13.2 and speedy-antlr-tool 1.4.4, preserving the hand-written listener.
Generation uses grammar basenames so checkout paths do not leak into generated
headers. The new regression checks wildcard and named averaging classes, plus
malformed controls, on both Python and compiled paths.

The initial native build was deliberately cancelled after the regression was
confirmed. Its incomplete build is not test evidence. The repaired build uses
bounded compilation concurrency, preserves normal compiler object ordering and
failure propagation, omits temporary debug information and removes build objects
before the builder layer is finalized. It still rebuilds all native parsers.


## Completed checks

The final image admits all 43 accelerators and matches installed/source Python
bytes. The full suite passes 133 tests without failures or skips. Four grammar
regressions pass in the native-enabled profile and again with acceleration
explicitly disabled. Ten hash-equivalent selected packet workflows retain the
same native outcomes as v8 (three OK, five warnings, two expected negative
errors), and nine separate shift sources retain the bounded accounting profile.
The archived 2MTX bundle completes with zero native errors; its qualification
limits remain explicit. See `QUALIFICATION.md` for the image ID and evidence scope.

Generated ANTLR output retains the generator's canonical formatting, consistent
with the repository's existing generated-code lint exclusions. Hand-written
changes pass the repository's lint gates; the hand-written CNS listener is
unchanged.
