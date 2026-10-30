set shell := ["bash", "-euo", "pipefail", "-c"]

help:
    @just --list

# Install all current fork fixes, preserving the already installed science resources.
build-docker: qualification-build
    docker tag bmrb-extract/utils-nmr:qualification-20261029-v8 bmrb-extract/utils-nmr:local-fork

qualification-build:
    docker build --platform linux/amd64 -f Dockerfile.qualification -t bmrb-extract/utils-nmr:qualification-20261029-v8 .

# Build current source first; tests use its image with fixed science resources.
test *args: image-source-check
    docker run --rm --pull=never --platform linux/amd64 --user "$(id -u):$(id -g)" -e FULLTEST=1 -e HOME=/tmp -e PYTHONPATH=/opt/py-packages:/opt/py-wwpdb_utils_nmr -v "$(pwd):/workspace" -w /tmp bmrb-extract/utils-nmr:qualification-20261029-v8 python -m unittest discover -v -s /workspace/wwpdb/utils/tests-nmr-tox -p '*Tests.py' {{args}}

lint:
    docker run --rm --pull=never --platform linux/amd64 -e PYTHONPATH=/workspace:/opt/py-packages:/opt/py-wwpdb_utils_nmr -v "$(pwd):/workspace:ro" -w /workspace bmrb-extract/utils-nmr:qualification-20261029-v8 flake8 --max-line-length=132 --config tox.ini wwpdb/utils/nmr wwpdb/utils/tests-nmr wwpdb/utils/tests-nmr-tox setup.py

# Pylint follows the repository's canonical tox invocation.
pylint:
    docker run --rm --pull=never --platform linux/amd64 -e PYTHONPATH=/workspace:/opt/py-packages:/opt/py-wwpdb_utils_nmr -v "$(pwd):/workspace:ro" -w /workspace bmrb-extract/utils-nmr:qualification-20261029-v8 pylint --disable=R,C --reports=n --rcfile=pylintrc --ignore=obsolete --ignore-paths='.*(Parser|Lexer)[.]py$,.*sa_.*[.]py$' wwpdb/utils/nmr wwpdb/utils/tests-nmr-tox

peak-regressions: image-source-check
    docker run --rm --pull=never --platform linux/amd64 --user "$(id -u):$(id -g)" -e PYTHONPATH=/opt/py-packages:/opt/py-wwpdb_utils_nmr -v "$(pwd):/workspace" -w /tmp bmrb-extract/utils-nmr:qualification-20261029-v8 python -m unittest discover -v -s /workspace/wwpdb/utils/tests-nmr-tox -p 'PeakReaderRegressionTests.py'

# The second discovery root is part of the upstream tox test contract.
test-library:
    docker run --rm --pull=never --platform linux/amd64 -e FULLTEST=1 -e PYTHONPATH=/opt/py-packages:/opt/py-wwpdb_utils_nmr bmrb-extract/utils-nmr:qualification-20261029-v8 python -c "import unittest; suite = unittest.defaultTestLoader.discover('/opt/py-wwpdb_utils_nmr/wwpdb/utils/nmr', pattern='*Tests.py'); result = unittest.TextTestRunner(verbosity=2).run(suite); raise SystemExit(not result.wasSuccessful())"

fix: lint peak-regressions schema-regressions lint-regressions

# Refuse to test stale installed source after editing the fork.
# Obsolete modules are excluded by .dockerignore and are not shipped in this image.
image-source-check:
    docker run --rm --pull=never --platform linux/amd64 -e PYTHONPATH=/workspace:/opt/py-packages:/opt/py-wwpdb_utils_nmr -v "$(pwd):/workspace:ro" bmrb-extract/utils-nmr:qualification-20261029-v8 python -c "from pathlib import Path; import hashlib; root = Path('/workspace/wwpdb/utils/nmr'); installed = Path('/opt/py-wwpdb_utils_nmr/wwpdb/utils/nmr'); changed = [str(p.relative_to(root)) for p in sorted(root.rglob('*.py')) if 'obsolete' not in p.relative_to(root).parts if not (installed / p.relative_to(root)).is_file() or hashlib.sha256(p.read_bytes()).digest() != hashlib.sha256((installed / p.relative_to(root)).read_bytes()).digest()]; print('Source/image differences:', changed); raise SystemExit(bool(changed))"

lint-regressions:
    docker run --rm --pull=never --platform linux/amd64 -e PYTHONPATH=/workspace:/opt/py-packages:/opt/py-wwpdb_utils_nmr -v "$(pwd):/workspace:ro" -w /workspace bmrb-extract/utils-nmr:qualification-20261029-v8 pylint --disable=R,C --reports=n --rcfile=pylintrc wwpdb/utils/tests-nmr-tox/PeakReaderRegressionTests.py wwpdb/utils/tests-nmr-tox/ReportSchemaRegressionTests.py

schema-regressions:
    just test -k ReportSchemaRegressionTests
