#!/usr/bin/env python3
##
# File: build_qualification_accelerators.py
# Date: 30-Oct-2026
##
"""Bound native compilation concurrency for the pinned qualification builder.

setuptools build_clib otherwise compiles the unchanged ANTLR runtime serially
under amd64 emulation. Invoke the original public compiler method per source;
retain source/object ordering and propagate every compilation failure.
"""
import runpy
from concurrent.futures import ThreadPoolExecutor

# The builder pins setuptools 84.2.1, which supplies this non-stdlib namespace.
# Pylint runs in the resource image without the build-only setuptools .pth shim;
# real native compilation validates the import. Its generic stdlib deprecation
# warning does not apply to this vendored compiler implementation.
from distutils.compilers.C.base import Compiler  # pylint: disable=import-error,deprecated-module

ORIGINAL_COMPILE = Compiler.compile
WORKERS = 4


def parallel_compile(compiler, sources, *args, **kwargs):
    """Compile sources concurrently while retaining the normal compiler contract."""
    def compile_source(source):
        return ORIGINAL_COMPILE(compiler, [source], *args, **kwargs)

    with ThreadPoolExecutor(max_workers=WORKERS) as executor:
        results = executor.map(compile_source, sources)
        return [obj for objects in results for obj in objects]


if __name__ == '__main__':
    Compiler.compile = parallel_compile
    runpy.run_path('setup.py', run_name='__main__')
