# Local converter qualification

Use the local qualification harness to test this fork with fixed science resources:

```sh
just build-docker
just fix
just test
```

`Dockerfile.qualification` overlays the entire current NMR package onto the
previously installed resource image pinned by digest. It does not run CCD,
chemical-shift-statistics or dictionary updaters. Test-tool versions are pinned
in `qualification-requirements.txt`. The image must already exist locally; this
is not a public, independently downloadable release build.

The image is tagged `bmrb-extract/utils-nmr:qualification-20261029-v8` and
`bmrb-extract/utils-nmr:local-fork`. Before testing, `image-source-check` compares
all shipped Python files with the source checkout. Obsolete modules are excluded
as declared in `.dockerignore`. Rebuild after source changes. Tests are mounted
from this checkout and use the installed image modules and accelerators.

On 2026-10-29, the v7 full suite ran **126 tests with no failures or skips**, using
Python 3.14.10 on Linux amd64 through Docker. This includes the formerly skipped
cleaned-input case, now asserting its real anomaly/rejection diagnostics.
The current suite includes two additional input-boundary regressions (128 total);
its full acceptance run is pending. `just fix` runs flake8, the eight peak/workflow regressions, the unknown-content
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
