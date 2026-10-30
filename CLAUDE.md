# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`NmrDpUtility` is the backend NMR data-processing tool for the wwPDB OneDep and
BMRB `bmrb_extract` systems. It takes a coordinate file plus assorted NMR data
files (chemical shifts, restraints, peak lists in ~50 vendor formats) and emits a
combined NEF or NMR-STAR file, reporting status through a JSON file validated
against `wwpdb/utils/tests-nmr/json-schema/`.

`wwpdb/utils/nmr/README.md` documents the public API (`setSource`, `setLog`,
`addInput`, `op`), the 13 workflow operations, the report schema, standalone mode
and the ANTLR C++ accelerators. Read it before changing anything user-facing.

## Commands

CI (`azure-pipelines.yml`, py3.9) runs four tox environments. Reproduce them exactly:

```bash
tox -e format_pep8-py39     # flake8 --max-line-length=132
tox -e lint_pylint-py39     # pylint --disable=R,C --rcfile=pylintrc
tox -e py39                 # unit tests
tox -e test_coverage-py39   # coverage, --fail-under=20
```

Directly, without tox:

```bash
python -m unittest discover -v -s wwpdb/utils/tests-nmr-tox -p "*Tests.py"
FULLTEST=1 python -m unittest discover -v -s wwpdb/utils/tests-nmr-tox -p "*Tests.py"
python -m unittest discover -v -s wwpdb/utils/tests-nmr-tox -p "NmrDpUtilityTests.py" -k <test_name>
```

- **`FULLTEST=1` matters.** The suite now has 137 tests. Heavy workflow tests
  are skipped without it. The local qualification harness sets it automatically
  and runs all tests without the FULLTEST-dependent skips. Tox may additionally skip the three optional
  `c_listener_util` tests when that accelerator is not built. See
  `QUALIFICATION.md` for the exact tested interpreter and image scope.
- `format_black` is in `tox.ini`'s envlist but **not** wired into CI, and the
  codebase is not black-formatted. Do not run black.
- Two test directories: `tests-nmr-tox/` (`*Tests.py`) is what CI runs;
  `tests-nmr/` (`test_*.py`, `fulltest_*.py`) is only linted, never run by CI,
  yet holds the real per-ticket coverage and all the mock data.
- Test runs leave artifacts *next to their input data* (`utils_nmr/` cache dirs
  and `*-log.json` under `tests-nmr/mock-data*`). These are gitignored now, but
  **never `git add -A`** here — check `git diff --name-only` before committing.

### Known failures

The complete local v8 qualification passed 128 tests with no failures or skips.
The formerly skipped `test_nmr_str2str_deposit_cleaned` is an active diagnostic
regression: its historical input contains anomalous shifts and is expected to
retain an Error result. Whole-repository Pylint has additional diagnostics under
the pinned modern tool version; focused checks on new regressions pass. See
`QUALIFICATION.md` for validation limits.

`test_get_nef_atom` used to fail. The cause was test import order:
`commonsetup` mocks `wwpdb.utils.config.ConfigInfo` with the test CCD
(`tests-nmr-tox/data/components`), and `ChemCompUtil` reads `CC_CVS_PATH` from
`ConfigInfo` when it is imported. A test module that imported `wwpdb.utils.nmr`
before `commonsetup` therefore bound the real CCD path, and `HEM` became
`"Unknown non-standard residue HEM found."`. Every module in `tests-nmr-tox/`
now imports `commonsetup` first. Keep it that way in new test modules,
or results start depending on which module unittest happens to load first.

Flake8 and the focused regression checks pass at the local harness invocations.

## Architecture

### Shared state: NmrDpRegistry

`NmrDpUtility` owns one `NmrDpRegistry` dataclass (`NmrDpRegistry.py`) holding
*all* mutable processing state — parsed data, report, coordinate caches, ~80
mode flags. Every collaborator receives it as `registry` and reaches everything
through `self._reg`. Collaborators find each other through the registry rather
than by import:

| field | class | role |
|---|---|---|
| `dpS` | `NmrDpMrSplitter` | split public PDB-MR files |
| `dpA` | `NmrDpFirstAid` | rescue malformed input |
| `dpV` | `NmrDpValidation` | validation and statistics |
| `dpR` | `NmrDpRemediation` | remediation and legacy merge |

Constants live in `NmrDpConstant.py` (~9k lines, no logic). `NmrDpReport.py` is
the JSON report model.

`NmrDpUtility.op(op)` looks a workflow op up in `__procTasksDict` and runs an
ordered list of bound `self.__*` task methods. Task order is a contract — e.g.
`updatePolymerSequence` must precede `updateEntitySaveframe`, which reads state
the former sets.

### NmrDpValidation / NmrDpRemediation are assembled from mixins

Both were single 20k-line classes. Each is now a thin facade over 12 modules:

```
NmrDpValidation.py        facade: class NmrDpValidation(NmrDpValidationInput, ..., NmrDpValidationOutStats)
NmrDpValidationBase.py    __init__, __slots__, and every helper with >1 caller
NmrDpValidation<Topic>.py one mixin each: Input, Coord, Loop, Nomencl, Cs, Mr, Pk,
                          CoordChk, CsStats, MrStats, OutStats
```

Rules that keep this working — breaking them fails at import or, worse, silently:

- **Only the base declares non-empty `__slots__`.** Two non-empty-slots bases
  raise `TypeError: multiple bases have instance lay-out conflict`. Mixins and
  the facade use `__slots__ = ()`.
- **Private attributes use one underscore** (`self._reg`, not `self.__reg`).
  Name mangling binds to the *defining* class, so `self.__reg` in a mixin would
  resolve to `_NmrDpValidationCs__reg` and fail.
- **A helper called from more than one mixin belongs on the base**, otherwise
  pylint raises `no-member` (CI gates on E and W).
- **The concrete class must keep its name.** `self.__class_name__` is embedded in
  `report.error.appendDescription('internal_error', ...)`, so it reaches the
  persisted JSON report. Registry-wired delegate objects would silently change
  report text; mixins preserve `self.__class__.__name__` for free.
- The facades re-export module-level helpers so existing
  `from ...NmrDpValidation import is_like_planality_boundary` imports keep
  working. Those need `# noqa: F401 pylint: disable=unused-import` — the noqa
  alone does not satisfy pylint.

### Format readers: mr/, pk/, cs/

~50 input formats, each with a generated ANTLR4 `<Fmt>Lexer.py` / `<Fmt>Parser.py`
(excluded from lint), a hand-written `<Fmt>ParserListener.py`, and a `<Fmt>Reader.py`
entry point. Readers share a two-pass protocol: `parse()` returns
`(listener, parser_err_listener, lexer_err_listener)`; if
`listener.getReasonsForReparsing()` is not None the caller rebuilds the reader
with those reasons and parses again.

The per-format dispatch in `NmrDpRemediationPk`, `NmrDpRemediationLegacy{Cs,Mr,Pk}`
is table-driven (`PK_READERS`, `LEGACY_CS_READERS`, `LEGACY_MR_READERS`, ...).
Formats whose control flow genuinely differs keep explicit branches — SPARKY and
NmrView fall back to alternate readers, AMBER sources reasons from
`reader.getReasons()`, CHARMM/CNS/XPLOR/SCHRODINGER configure extra reader modes.
**Add a new format to the table; do not add a branch** unless it truly cannot be
described by the table's keys.

`sa_*.py` modules gate optional speedy-antlr C++ accelerators behind
`USE_CPP_IMPLEMENTATION`, falling back to the ANTLR Python runtime. Their sources
live in `cpp_src/` and are compiled only for the container image.

### Row layouts are addressed by named column constants

`NmrDpRemediationCsLoop` builds NMR-STAR `_Atom_chem_shift` rows as plain lists.
Positions are named (`CS_AUTH_ASYM`, `CS_ORIG_SEQ`, ...) and derived from the same
`CS_TAGS` tuple that builds the loop tags, so indices and column order have one
definition and a transposed multi-column assignment is visible. `CS_TAGS`/`CS_*`
and the NEF-layout `NEF_*` set are in `NmrDpRemediationCsLoop.py`; `CCA_TAGS`/
`CCA_*` for `registry.chem_comp_asm_dat` rows are in `NmrDpRemediationBase.py`.
Keep extending these; do not reintroduce a bare index in a row assignment.

## The dual-import idiom

Every module in the package must do:

```python
try:
    from wwpdb.utils.nmr.X import a, b     # OneDep environment
except ImportError:
    from nmr.X import a, b                 # standalone mode
```

Standalone mode sets `PYTHONPATH=<repo>/wwpdb/utils` so the package resolves as
bare `nmr.*`; this is how the Docker image runs. **CI never exercises the
fallback branch**, so a name present in only one of the two lists raises
`NameError` in production. There is already a fix commit for exactly that. When
touching imports, verify both branches bind the same names, and smoke-test:

```bash
PYTHONPATH=$PWD/wwpdb/utils python -c "from nmr.NmrDpUtility import NmrDpUtility; NmrDpUtility()"
```

## Conventions

- Module header is `##` / `# File:` / `# Date:` / `# Updates:` / `##`, then the
  docstring and `__docformat__`/`__author__`/`__email__`/`__license__`/`__version__`.
  Changelog entries are `DD-Mon-YYYY  M. Yokochi - <description> (DAOTHER-NNNNN)`.
- **`__version__` is bumped repo-wide in lockstep**, in its own "Bump versions"
  commit, including `wwpdb/utils/nmr/__init__.py` (which drives the package
  version). Never bump one module alone. Record cross-module work as a single
  `# Updates:` line in `NmrDpUtility.py`.
- `setup.py` uses `find_packages()`, so a new `.py` beside its siblings needs no
  registration. `MANIFEST.in` concerns data files only.
- `flake8-import-order` applies to `wwpdb/utils/nmr` (waived only for the test
  dirs). Mirror neighbouring files' import order rather than sorting.
- Line length 132 for both flake8 and pylint. Existing over-length lines carry
  `# noqa: E501`.
- `ChemCompUtil` accessors guard with `if not self.updateChemCompDict(compId):
  return ...`. Do **not** reintroduce a `compId != self.lastCompId and`
  short-circuit in front of it: `updateChemCompDict` sets `lastCompId` even when
  the lookup fails, so the short-circuit lets a failed comp_id skip the guard and
  silently inherit the previously loaded compound's atom and bond lists.

## Refactoring these modules

Hard-won lessons; the modules are large and the test suite reaches few branches,
so static checks carry the weight.

- **"The moved code is AST-identical" proves nothing about its environment.**
  Deleting a `for offset in ...` loop while its body still referenced `offset`
  shipped a `NameError` past an AST-equivalence gate. It does not look like an
  undefined name either, because the name was bound elsewhere in the enclosing
  function, so Python resolves it as a free variable and fails only when that
  branch runs. **Check for bindings too:** compare the local/free classification
  of every name in every scope against a baseline revision and flag
  `local -> free`, which is the signature of a deleted binding.
- **Verify before rewriting.** For per-format branch families, canonicalise them
  (hide the class name, the label, and continuation-line indentation, which
  varies with name length) and count *distinct shapes*; that number, not the line
  count, decides how much can collapse. Then prove a parameterised template
  reproduces each original branch byte-for-byte.
- **A clone detector's "reclaimable lines" runs 3-4x optimistic here.** Most
  matches are fragments straddling a block boundary (a "24-line clone" ending
  mid-`if`), or contain `continue`/`break`, or rebind 20+ names. Per candidate,
  check: does it `ast.parse` as a standalone suite, what does it rebind, does it
  contain `continue`/`break`/`return`.
- Extract into a **nested closure with `nonlocal` for every name the block
  rebinds**; that is exactly equivalent to inlining, so no reasoning is needed
  about which values outlive the block. Compute the `nonlocal` set with a
  **scope-aware** walk — a plain `ast.walk` descends into nested functions and
  wrongly reports names as shared.
- When collapsing a `for offset` loop whose forward/backward halves are
  duplicated, replace the header with a **generator** yielding the rows in order.
  A nested `for` over `(+offset, -offset)` silently retargets every `break`.
- Renaming private names shortens identifiers, which breaks visually-aligned
  continuation lines (~1400 E127/E124 in one pass). Re-align by iterating flake8
  `--select=E12x`; use **flake8, not pycodestyle**, because pycodestyle treats a
  bare `# noqa` as suppressing every code and so hides E127 on lines carrying
  `# noqa: E501`.
- Slice assignment does **not** enforce arity: CPython normalises
  `slice(16, 20, 1)` to a plain slice, so a wrong-length right-hand side silently
  resizes the row, where tuple unpacking raises.
