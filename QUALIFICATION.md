# Local converter qualification

Use the local qualification harness to test this fork with fixed science resources:

```sh
just build-docker
just fix
just test
```

`Dockerfile.qualification` overlays the entire current NMR package onto the
previously installed resource image pinned by digest. Its isolated builder rebuilds
all native parser accelerators against the current generated source;
`just accelerators` rejects missing, disabled or fingerprint-mismatched binaries.
The builder installs a compiler and pinned setuptools, which are build-only.
The runtime does not run CCD,
chemical-shift-statistics or dictionary updaters. Test-tool versions are pinned
in `qualification-requirements.txt`. The image must already exist locally; this
is not a public, independently downloadable release build.

The active qualification image is `bmrb-extract/utils-nmr:qualification-20261030-v10`.
The v8 image is retained for historical evidence. `WWPDB_NMR_QUALIFICATION_IMAGE`
can select an explicit existing image. `just activate` promotes an already tested
image to `local-fork` after source/admission checks; `build-docker` builds then
activates. Historical versioned tags are retained. Before testing, `image-source-check` compares
all shipped Python files with the source checkout. Obsolete modules are excluded
as declared in `.dockerignore`. Rebuild after source changes. Tests are mounted
from this checkout and use the installed image modules and accelerators.

On 2026-10-29, the full v8 suite ran **128 tests with no failures or skips**,
using Python 3.14.10 on Linux amd64 through Docker (717.590 seconds). This
includes the formerly skipped cleaned-input case, now asserting its real
anomaly/rejection diagnostics, and the new input-boundary regressions.
`just fix` passes flake8, eight peak/workflow regressions, the unknown-content
schema regression and focused Pylint checks on the new regression modules.
The second upstream discovery root (`just test-library`) contains no tests.

The eight repaired defects are:

- Undefined experiment method when coordinates omit NMR experiment metadata.
- Tuple used instead of the standalone XEASY translator dependency.
- Empty parent list-ID Counter causing an unhandled IndexError.
- Coordinate/NMR alignment using the test sequence as its own reference.
- Unidentified content (`null`) violating the diagnostic report schema.
- NEF/STAR translation dropping the report path configured through `setLog()`.
- Assembly annotation indexing an absent nullable Experimental_data_reported column.
- Empty required chemical-shift lists bypassing the input precondition and causing
  later source/type failures; unidentified input types are also guarded in legacy
  restraint and peak validation.

The nullable-content schema repair must also be applied to clients that vendor
`nmr-data-processing-report-schema-v5.json`. The local BMRB Extract backend has
been synchronized. Unknown content remains distinct from zero observations;
non-object values and invalid count values still fail validation.

The full-repository Pylint 4.1.3 invocation has additional mixin-inference and
optional-return diagnostics. Focused Pylint on the added regression modules
passes; this qualification does not certify the entire repository Pylint-clean.
No broad diagnostic suppressions or formatter rewrite are included.

Unit tests and local native conversion do not certify the scientific validity
of every fixture or a hosted deployment. The evaluation repository retains
hash-bound packet inputs, actual outputs, native reports, stage events, source
checks and documented unknowns. Its common synthetic geometry remains an atom
inventory; unsupported producer dialects and topology mappings remain explicit
qualification gaps. See the evaluation `docs/FIXTURE_QUALIFICATION.md`.


## Final-output diagnostic schema checkpoint

Analyzing the final conversion phase exposed report-schema incompatibilities
that the initial phase did not exercise. Output statistics can legitimately
retain unknown entry titles and assembly names/component counts, or unavailable
chemical-shift/restraint summaries as null. The v5 schema now represents those
specific unknown informational values explicitly. Present title strings and
positive component counts retain their constraints; wrong types and invalid
counts remain rejected. This does not waive required scientific input checks.
The bundled backend schema is synchronized.

`just fix` passes flake8, eight existing peak/workflow regressions, both schema
regressions and focused Pylint. Executable converter source and the pinned v8
image are unchanged, so the previous 128-test full executable-suite result is
historical evidence for that image; the two current schema regressions are the
focused validation of this checkpoint, not a newly claimed full suite run.


## Upstream rebase checkpoint

See [UPSTREAM_REVIEW.md](UPSTREAM_REVIEW.md) for all five reviewed branch heads,
changes affecting packets and the release-line disposition. The selected base
is `yokochi47/main` at `a450780df1c0c445b5617edfb9ef2ac98af2adb7`,
which contains the reviewed `wwPDB/bmrb-extract_devel` and `yokochi47/devel`
heads. All ten local commits replayed with identical patches; the original head
is preserved under an archive ref and in an external recovery bundle.

The v8 image and historical campaign outputs remain unchanged. New runtime tests
must use the v9 image and its recorded immutable ID. A rebase or an unchanged
5.4.4 software version does not establish behavioral or hosted equivalence.


### Completed v9 validation

The repaired v9 image has immutable ID
`sha256:204fc73559e846349e3cb70175fa29127a77dfb7829ed212b3147463c8c03832`.
All 42 parser accelerators plus `c_listener_util` are compiled and admitted;
`image-source-check` finds no installed/source Python differences. The complete
`FULLTEST=1` suite passed **133 tests with no failures or skips** (261.066 seconds).
The second upstream discovery root still contains zero tests and is not counted
as additional coverage. `just fix` passes flake8, peak/workflow and schema
regressions, focused Pylint, accelerator admission, and native/Python grammar
controls. The same four grammar tests also pass with acceleration explicitly
disabled; ordinary pure-Python installations retain their supported test path.

A fresh selected packet snapshot at evaluation
`analysis/qualification-20261030-upstream-v9/` executed ten workflows on v9:
three native OK, five warnings and two expected negative errors. All ten packet
manifests match the prior v8 campaign byte-for-byte; no native status or error
count changed, and there were no runtime or report-analysis failures. All nine
separately supplied shift inputs pass the bounded 13-field accounting profile.
The two intentional negative errors make the packet runner exit nonzero as
expected; they are not silently counted as accepted positive depositions.

The matched 2MTX source bundle also completed on v9 with zero native errors and
eight visible phase-warning occurrences (four underlying rare assignments).
Its bounded inventory still retains all 134 shifts and 322 coordinate atoms.
The strict expanded-source shift profile remains inconclusive, including the
24 grouped author-name normalizations. Full restraint/physical qualification
remains open. These selected checks do not replace a fresh full corpus campaign
or establish hosted or scientific equivalence.

The source ledger, recovery bundle, generation/build logs and full test output
are retained privately under `bmrb-upstream-review-20261030/` beside the checkouts.
The prior v8 image and historical evidence remain intact.


## Explicit shift-core qualification (v10)

The v10 image is
`sha256:99801e53d00e031de908a3fbb427ea25563be170f8cfb55b941939554255a6dd`.
It repairs two additional defects exposed by independent shift-only controls:

- Annotation generated `Sample_condition_list_ID=1` and
  `Chem_shift_reference_ID=1` without corresponding frames. Missing metadata
  now remains missing instead of creating dangling foreign keys.
- Validation of a NEF-derived molecular table accessed absent
  `Chem_comp_assembly.Seq_ID`. When that redundant column is absent or missing,
  extraction now uses canonical `Comp_index_ID`; explicit Seq_ID behavior is
  retained. A reader-level regression exercises missing, null and present IDs.

The full suite passed **135 tests with no failures or skips** in 239.523 seconds.
`just fix` passes, including ten peak/workflow regressions, both schema tests,
focused Pylint, admission of all 43 accelerators, and native/Python grammar
regressions. This retains the existing Linux amd64 / Python 3.14 scope.

For this Python-only change, `just qualification-python-build` verifies the
exact v9 base image ID, compares every native C/C++ source/header byte before
reusing its compiled accelerators, overlays the complete current NMR package,
and preserves science resources. A native-source difference stops that build;
use the full `qualification-build` with a new tag when changing native code.
`just activate` promotes a tested image only after source/admission checks.

The evaluation repository's `analysis/shift-core-v1-qualified/` records two
separately scoped controls, each run twice on this exact image. Both have native
OK, zero warnings/errors, no internal diagnostics and byte-identical repeated
outputs. The independent explicit-shifts-v1 profile verifies all admitted
observations (six STAR, three NEF), identifiers, relationships, uncertainties,
precision, missingness, metadata and derived counts. Qualification includes
source/output/checker/reference hashes and actual attempt evidence.

These are constructed interchange-conversion controls. Missing experimental
sample-condition/calibration metadata remains unknown, and the production
NMR-STAR archival requiredness audit still reports those two unavailable
references. They are not full dictionary/deposition certificates, experimental
producer runs, or replacements for the historical 66-packet campaign. No
coordinate/restraint baseline has been qualified by this shift-only work.
