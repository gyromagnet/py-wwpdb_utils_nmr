set shell := ["bash", "-euo", "pipefail", "-c"]

QUALIFICATION_IMAGE := env("WWPDB_NMR_QUALIFICATION_IMAGE", "bmrb-extract/utils-nmr:qualification-20261030-v16")

help:
    @just --list

# Install all current fork fixes, preserving the already installed science resources.
build-docker: qualification-build activate

qualification-build:
    docker build --platform linux/amd64 -f Dockerfile.qualification -t {{QUALIFICATION_IMAGE}} .

# Build current source first; tests use its image with fixed science resources.
test *args: image-source-check
    docker run --rm --pull=never --platform linux/amd64 --user "$(id -u):$(id -g)" -e FULLTEST=1 -e HOME=/tmp -e PYTHONPATH=/opt/py-packages:/opt/py-wwpdb_utils_nmr -v "$(pwd):/workspace" -w /tmp {{QUALIFICATION_IMAGE}} python -m unittest discover -v -s /workspace/wwpdb/utils/tests-nmr-tox -p '*Tests.py' {{args}}

lint:
    docker run --rm --pull=never --platform linux/amd64 -e PYTHONPATH=/workspace:/opt/py-packages:/opt/py-wwpdb_utils_nmr -v "$(pwd):/workspace:ro" -w /workspace {{QUALIFICATION_IMAGE}} flake8 --max-line-length=132 --config tox.ini wwpdb/utils/nmr wwpdb/utils/tests-nmr wwpdb/utils/tests-nmr-tox setup.py tools/build_qualification_accelerators.py tools/gen_speedy_antlr.py

# Pylint follows the repository's canonical tox invocation.
pylint:
    docker run --rm --pull=never --platform linux/amd64 -e PYTHONPATH=/workspace:/opt/py-packages:/opt/py-wwpdb_utils_nmr -v "$(pwd):/workspace:ro" -w /workspace {{QUALIFICATION_IMAGE}} pylint --disable=R,C --reports=n --rcfile=pylintrc --ignore=obsolete --ignore-paths='.*(Parser|Lexer)[.]py$,.*sa_.*[.]py$' wwpdb/utils/nmr wwpdb/utils/tests-nmr-tox

peak-regressions: image-source-check
    docker run --rm --pull=never --platform linux/amd64 --user "$(id -u):$(id -g)" -e PYTHONPATH=/opt/py-packages:/opt/py-wwpdb_utils_nmr -v "$(pwd):/workspace" -w /tmp {{QUALIFICATION_IMAGE}} python -m unittest discover -v -s /workspace/wwpdb/utils/tests-nmr-tox -p 'PeakReaderRegressionTests.py'

# The second discovery root is part of the upstream tox test contract.
test-library:
    docker run --rm --pull=never --platform linux/amd64 -e FULLTEST=1 -e PYTHONPATH=/opt/py-packages:/opt/py-wwpdb_utils_nmr {{QUALIFICATION_IMAGE}} python -c "import unittest; suite = unittest.defaultTestLoader.discover('/opt/py-wwpdb_utils_nmr/wwpdb/utils/nmr', pattern='*Tests.py'); result = unittest.TextTestRunner(verbosity=2).run(suite); raise SystemExit(not result.wasSuccessful())"

fix: lint peak-regressions schema-regressions lint-regressions accelerators upstream-regressions upstream-python-regressions

# Refuse to test stale installed source after editing the fork.
# Obsolete modules are excluded by .dockerignore and are not shipped in this image.
image-source-check:
    docker run --rm --pull=never --platform linux/amd64 -e PYTHONPATH=/workspace:/opt/py-packages:/opt/py-wwpdb_utils_nmr -v "$(pwd):/workspace:ro" {{QUALIFICATION_IMAGE}} python -c "from pathlib import Path; import hashlib; root = Path('/workspace/wwpdb/utils/nmr'); installed = Path('/opt/py-wwpdb_utils_nmr/wwpdb/utils/nmr'); changed = [str(p.relative_to(root)) for p in sorted(root.rglob('*.py')) if 'obsolete' not in p.relative_to(root).parts if not (installed / p.relative_to(root)).is_file() or hashlib.sha256(p.read_bytes()).digest() != hashlib.sha256((installed / p.relative_to(root)).read_bytes()).digest()]; print('Source/image differences:', changed); raise SystemExit(bool(changed))"

lint-regressions:
    docker run --rm --pull=never --platform linux/amd64 -e PYTHONPATH=/workspace:/opt/py-packages:/opt/py-wwpdb_utils_nmr -v "$(pwd):/workspace:ro" -w /workspace {{QUALIFICATION_IMAGE}} pylint --disable=R,C --reports=n --rcfile=pylintrc wwpdb/utils/tests-nmr-tox/PeakReaderRegressionTests.py wwpdb/utils/tests-nmr-tox/ReportSchemaRegressionTests.py wwpdb/utils/tests-nmr-tox/UpstreamGrammarRegressionTests.py tools/build_qualification_accelerators.py

schema-regressions:
    just test -k ReportSchemaRegressionTests

# A changed grammar must not silently fall back to Python or use a stale native parser.
accelerators: image-source-check
    docker run --rm --pull=never --platform linux/amd64 -e PYTHONPATH=/opt/py-packages:/opt/py-wwpdb_utils_nmr -v "$(pwd):/workspace:ro" {{QUALIFICATION_IMAGE}} python /workspace/tools/check_accelerators.py --require-compiled

upstream-regressions:
    just test -k UpstreamGrammarRegressionTests

# Keep the hand-written Python listener intact; copy only generated lexer/parser files.
regenerate-cns:
    #!/usr/bin/env bash
    set -euo pipefail
    grammar_root="$(pwd)/wwpdb/utils/tests-nmr/antlr-grammars-v4.10"
    generated_root="$(mktemp -d)"
    trap 'rm -rf "$generated_root"' EXIT
    (cd "$grammar_root" && uv run --no-project --with antlr4-tools==0.2.2 --with antlr4-python3-runtime==4.13.2 python -c 'from antlr4_tool_runner import tool; tool()' -v 4.13.2 -Dlanguage=Python3 -no-visitor -Xexact-output-dir -o "$generated_root" CnsMRLexer.g4 CnsMRParser.g4)
    cp "$generated_root/CnsMRLexer.py" "$generated_root/CnsMRParser.py" wwpdb/utils/nmr/mr/
    uv run --no-project --with antlr4-tools==0.2.2 --with speedy-antlr-tool==1.4.4 --with antlr4-python3-runtime==4.13.2 python tools/gen_speedy_antlr.py CnsMR

upstream-python-regressions: image-source-check
    docker run --rm --pull=never --platform linux/amd64 --user "$(id -u):$(id -g)" -e WWPDB_NMR_DISABLE_ACCEL=1 -e PYTHONPATH=/opt/py-packages:/opt/py-wwpdb_utils_nmr -v "$(pwd):/workspace" -w /tmp {{QUALIFICATION_IMAGE}} python -m unittest discover -v -s /workspace/wwpdb/utils/tests-nmr-tox -p 'UpstreamGrammarRegressionTests.py'

# Promote an already tested image without rebuilding or changing historical tags.
activate: image-source-check accelerators
    docker tag {{QUALIFICATION_IMAGE}} bmrb-extract/utils-nmr:local-fork

# Reuse v9 accelerators only after proving exact base and native-source identity.
qualification-python-build:
    docker image inspect bmrb-extract/utils-nmr:qualification-20261030-v9 | python3 -c 'import json, sys; assert json.load(sys.stdin)[0]["Id"] == "sha256:204fc73559e846349e3cb70175fa29127a77dfb7829ed212b3147463c8c03832", "Unexpected v9 base"'
    docker build --pull=false --platform linux/amd64 -f Dockerfile.qualification-python -t {{QUALIFICATION_IMAGE}} .
