#!/usr/bin/env python3
##
# File: tools/gen_speedy_antlr.py
# Date: 7-Aug-2026
#
# Updates:
# 01-Oct-2026  M. Yokochi - lex a one-byte-per-character str in place (latin1_input_stream.h) (DAOTHER-7829, 9785)
# 08-Oct-2026  M. Yokochi - require speedy-antlr-tool 1.4.4, which ships the parse-tree reference leak fix
#                           (upstream PR #19); drop the local support-library patch
# 14-Oct-2026  M. Yokochi - export GRAMMAR_FINGERPRINT from each accelerator and admit it in the shim only when it
#                           matches the Python parser (AccelUtil.py); add --repatch (DAOTHER-9785)
##
""" Generate the speedy-antlr-tool C++ accelerators for this package's grammars.

    Grammars are DISCOVERED from the *Reader.py classes rather than listed here:
    each Reader names its lexer/parser in its imports, its grammar entry rule in
    the parse call, and whether it wants the SLL prediction mode. Discovery keeps
    the accelerators from drifting away from the Python path, and means a new
    grammar needs no edit to this script.

    For each grammar this
      1. runs the ANTLR tool with the C++ target over the .g4 pair,
      2. runs speedy_antlr_tool.generate() against the *already committed* Python
         parser to emit the sa_<grammar>.py shim and the C++ translator,
      3. patches do_parse() and the shim so the prediction mode is a per-call
         argument (see PREDICTION_MODE_NOTE below), and do_parse() so that a str
         stored one byte per character is lexed in place by
         cpp_src/latin1_input_stream.h instead of a UTF-32 copy,
      4. patches the module init to export GRAMMAR_FINGERPRINT, and the shim to
         use the accelerator only when that matches its Python parser
         (AccelUtil.py),
      5. writes cpp_src/speedy_antlr_manifest.json, which setup.py reads to
         declare the extensions and to compute the fingerprints.

    --repatch re-applies steps 3 to 5 to the committed sources without
    regenerating them, after a patch is added here; it needs neither Java nor
    speedy-antlr-tool.

    The Python target is deliberately NOT regenerated. The committed
    <Grammar>Lexer.py/<Grammar>Parser.py are byte-identical to a fresh
    'antlr4 -Dlanguage=Python3 -no-visitor' run with ANTLR 4.13.2, and
    regenerating in place would overwrite the hand-written
    <Grammar>ParserListener.py, which shadows ANTLR's generated listener base of
    the same name. (Note the -no-listener flag that speedy-antlr-tool's own
    example uses would strip the enterRule/exitRule hooks ParseTreeWalker needs.)

    Requires the 'antlr4-tools' and 'speedy-antlr-tool' (>= 1.4.4) pip packages,
    and Java (antlr4-tools provisions a JRE on first use). 1.4.3's shared
    speedy_antlr.cpp leaked the whole translated parse tree on every parse; the
    fix is upstream since 1.4.4 (amykyta3/speedy-antlr-tool#19).

    Usage:
        python3 tools/gen_speedy_antlr.py --list
        python3 tools/gen_speedy_antlr.py XplorMR NmrPipeCS
        python3 tools/gen_speedy_antlr.py --all
        python3 tools/gen_speedy_antlr.py --repatch --all
    @author: Masashi Yokochi
"""
import argparse
import collections
import glob
import io
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.request
import zipfile

ANTLR_VERSION = '4.13.2'

# The first release with the convert_ctx() reference leak fix (upstream PR #19),
# and the do_parse() template CPP_PATCHES are written against.
MIN_SPEEDY_ANTLR_VERSION = (1, 4, 4)

CPP_RUNTIME_URL = f'https://www.antlr.org/download/antlr4-cpp-runtime-{ANTLR_VERSION}-source.zip'

REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NMR_DIR = os.path.join(REPO_DIR, 'wwpdb', 'utils', 'nmr')
GRAMMAR_DIR = os.path.join(REPO_DIR, 'wwpdb', 'utils', 'tests-nmr', 'antlr-grammars-v4.10')
CPP_SRC_DIR = os.path.join(NMR_DIR, 'cpp_src')
CPP_RUNTIME_DIR = os.path.join(CPP_SRC_DIR, 'antlr4-cpp-runtime')
MANIFEST_PATH = os.path.join(CPP_SRC_DIR, 'speedy_antlr_manifest.json')

SUB_PACKAGES = ('mr', 'pk', 'cs')

# --- Reader introspection -----------------------------------------------------
# Imports are indented (they sit inside the try/except ImportError pair), hence
# the leading \s*.
LEXER_IMPORT = re.compile(r'^\s*from (?:wwpdb\.utils\.)?nmr\.(mr|pk|cs)\.(\w+Lexer) import', re.M)
PARSER_IMPORT = re.compile(r'^\s*from (?:wwpdb\.utils\.)?nmr\.(mr|pk|cs)\.(\w+Parser) import', re.M)
# The entry rule, in either the original inline block or the converted call:
#   inline:    tree = parser.<rule>()
#   converted: parseAntlr(<Lexer>, <Parser>, '<rule>', ...)
ENTRY_RULE_PATTERNS = (re.compile(r'^\s*tree\s*=\s*parser\.(\w+)\(\)', re.M),
                       re.compile(r'parseAntlr\(\s*\w+\s*,\s*\w+\s*,\s*[\'"](\w+)[\'"]', re.M))
# An *uncommented* SLL prediction mode, in either form. 26 of the 47 call sites
# have it commented out, so an 'in source' test would be wrong.
SLL_PATTERNS = (re.compile(r'^\s*parser\._interp\.predictionMode\s*=\s*PredictionMode\.SLL', re.M),
                re.compile(r'predictionModeSll\s*=\s*(?!False\b)', re.M))

GrammarSpec = collections.namedtuple(
    'GrammarSpec', 'name subPackage parserName lexerSubPackage lexerName entryRule readers useSll')

# --- Generated-code patches ---------------------------------------------------
PREDICTION_MODE_NOTE = """
The prediction mode must be settable per call rather than baked in at compile
time. Six mr/ Readers (XplorMR, CnsMR, CyanaMR, CharmmMR, SchrodingerMR,
CyanaNOA) choose it at runtime via 'if not isFilePath or self.__sll_pred', and
NmrDpMrSplitter escalates it on retry. The XML grammar makes it unavoidable:
AriaMRXReader drives it with SLL while AriaCSReader, AriaPKReader and
TopSpinPKReader drive the same accelerator with LL. speedy-antlr-tool offers no
hook, so do_parse() and the shim's parse() gain a flag.
"""

# Each entry is (description, [(old, new), ...]) with %(parser)s substituted.
# Every replacement must match exactly once, so a speedy-antlr-tool upgrade
# cannot silently drop a patch.
CPP_PATCHES = [
    ('do_parse() accepts a prediction-mode flag', [
        ('        PyObject *sa_err_listener = NULL;\n',
         '        PyObject *sa_err_listener = NULL;\n'
         '        int use_sll = 0;  // wwPDB: PredictionMode.SLL when true\n'),
        ('            "OOsO:do_parse",\n'
         '            &parser_cls, &stream, &entry_rule_name, &sa_err_listener\n',
         '            "OOsOp:do_parse",\n'
         '            &parser_cls, &stream, &entry_rule_name, &sa_err_listener, &use_sll\n'),
        ('        antlr4::tree::ParseTree *parse_tree;\n',
         '        // wwPDB: mirror the Python driver\'s prediction mode, or the two\n'
         '        // paths accept different inputs and report different errors.\n'
         '        if (use_sll) {\n'
         '            parser.getInterpreter<antlr4::atn::ParserATNSimulator>()\n'
         '                ->setPredictionMode(antlr4::atn::PredictionMode::SLL);\n'
         '        }\n'
         '        antlr4::tree::ParseTree *parse_tree;\n'),
    ]),
    # A str stored one byte per character (code points <= U+00FF, which covers every
    # ASCII input) is lexed in place by wwpdb::Latin1InputStream (latin1_input_stream.h)
    # with the code points and indices ANTLRInputStream would see. Otherwise the str
    # takes the UTF-8 copy, which ANTLRInputStream decodes into a UTF-32 copy of its
    # own, 4 bytes per character (DAOTHER-7829, 9785).
    ('do_parse() lexes a one-byte-per-character str in place', [
        ('#include "speedy_antlr.h"\n',
         '#include "speedy_antlr.h"\n'
         '#include "latin1_input_stream.h"  // wwPDB\n'),
        ('        // Extract input stream\'s string\n'
         '        char *cstrdata;\n'
         '        Py_ssize_t bufsize;\n'
         '        strdata = PyObject_GetAttrString(stream, "strdata");\n'
         '        if(!strdata) throw speedy_antlr::PythonException();\n'
         '\n'
         '        // PyUnicode_AsUTF8AndSize is not part of the stable ABI until python3.10\n'
         '        // To maximize backwards compatibility, Working around by converting to\n'
         '        // bytes, then to char instead\n'
         '        strdata_as_bytes = PyCodec_Encode(strdata, "utf-8", NULL);\n'
         '        if(!strdata_as_bytes) throw speedy_antlr::PythonException();\n'
         '        PyBytes_AsStringAndSize(strdata_as_bytes, &cstrdata, &bufsize);\n'
         '        if(!cstrdata) throw speedy_antlr::PythonException();\n'
         '\n'
         '        // Create an antlr InputStream object\n'
         '        antlr4::ANTLRInputStream cpp_stream(cstrdata, bufsize);\n',
         '        // Extract input stream\'s string\n'
         '        strdata = PyObject_GetAttrString(stream, "strdata");\n'
         '        if(!strdata) throw speedy_antlr::PythonException();\n'
         '        if(!PyUnicode_Check(strdata)) {\n'
         '            PyErr_SetString(PyExc_TypeError, "InputStream.strdata must be a str");\n'
         '            throw speedy_antlr::PythonException();\n'
         '        }\n'
         '\n'
         '        // Create an antlr InputStream object. wwPDB: a str stored one byte per\n'
         '        // character is lexed in place (strdata is held until the end of the parse);\n'
         '        // any other goes through the UTF-32 copy of ANTLRInputStream as before.\n'
         '        std::unique_ptr<antlr4::CharStream> cpp_stream_holder;\n'
         '        if(PyUnicode_KIND(strdata) == PyUnicode_1BYTE_KIND) {\n'
         '            cpp_stream_holder.reset(new wwpdb::Latin1InputStream(\n'
         '                reinterpret_cast<const unsigned char *>(PyUnicode_1BYTE_DATA(strdata)),\n'
         '                static_cast<size_t>(PyUnicode_GET_LENGTH(strdata))));\n'
         '        } else {\n'
         '            char *cstrdata;\n'
         '            Py_ssize_t bufsize;\n'
         '            strdata_as_bytes = PyCodec_Encode(strdata, "utf-8", NULL);\n'
         '            if(!strdata_as_bytes) throw speedy_antlr::PythonException();\n'
         '            if(PyBytes_AsStringAndSize(strdata_as_bytes, &cstrdata, &bufsize) < 0)\n'
         '                throw speedy_antlr::PythonException();\n'
         '            cpp_stream_holder.reset(new antlr4::ANTLRInputStream(cstrdata, bufsize));\n'
         '        }\n'
         '        antlr4::CharStream &cpp_stream = *cpp_stream_holder;\n'),
    ]),
    # setup.py defines WWPDB_GRAMMAR_FINGERPRINT from the serialized ATNs of the
    # generated C++ lexer and parser (AccelUtil.atnFingerprint). A build that does
    # not go through setup.py exports an empty one, which the shim rejects.
    ('the module exports the fingerprint of the grammar it was compiled from', [
        ('PyMODINIT_FUNC\n',
         '// wwPDB: see AccelUtil.py\n'
         '#ifndef WWPDB_GRAMMAR_FINGERPRINT\n'
         '#define WWPDB_GRAMMAR_FINGERPRINT ""\n'
         '#endif\n'
         '\n'
         'PyMODINIT_FUNC\n'),
        ('    PyObject *m = PyModule_Create(&module);\n'
         '    return m;\n',
         '    PyObject *m = PyModule_Create(&module);\n'
         '    if (m != NULL && PyModule_AddStringConstant(m, "GRAMMAR_FINGERPRINT", WWPDB_GRAMMAR_FINGERPRINT) < 0) {\n'
         '        Py_DECREF(m);\n'
         '        return NULL;\n'
         '    }\n'
         '    return m;\n'),
    ]),
]

# speedy-antlr-tool names the lexer after the parser ('<X>Parser' -> '<X>Lexer'),
# which breaks the three grammars that borrow another grammar's lexer through
# tokenVocab: NmrViewNPKParser uses NmrViewPKLexer, and SparkyNPKParser and
# SparkyRPKParser use SparkyPKLexer. Without this the generated C++ fails to
# compile on a missing '<X>Lexer.h'.
BORROWED_LEXER_CPP_PATCHES = [
    ('the generated C++ references the borrowed lexer', [
        ('#include "%(assumedLexer)s.h"', '#include "%(lexer)s.h"'),
        ('        %(assumedLexer)s lexer(&cpp_stream);', '        %(lexer)s lexer(&cpp_stream);'),
    ]),
]

# The shim imports the assumed lexer at module scope, so without this the whole
# sa_<grammar> module raises ModuleNotFoundError - which would take the Reader
# down rather than degrade to the Python parser.
BORROWED_LEXER_SHIM_PATCHES = [
    ('the shim imports the borrowed lexer', [
        ('from .%(assumedLexer)s import %(assumedLexer)s', 'from .%(lexer)s import %(lexer)s'),
        ('    lexer = %(assumedLexer)s(stream)', '    lexer = %(lexer)s(stream)'),
    ]),
]

SHIM_PATCHES = [
    ('parse()/_cpp_parse() forward the prediction-mode flag', [
        ('def parse(stream:InputStream, entry_rule_name:str, sa_err_listener:SA_ErrorListener=None) -> ParseTree:',
         'def parse(stream:InputStream, entry_rule_name:str, sa_err_listener:SA_ErrorListener=None,\n'
         '          prediction_mode_sll:bool=False) -> ParseTree:'),
        ('        return _cpp_parse(stream, entry_rule_name, sa_err_listener)\n'
         '    else:\n'
         '        return _py_parse(stream, entry_rule_name, sa_err_listener)\n',
         '        return _cpp_parse(stream, entry_rule_name, sa_err_listener, prediction_mode_sll)\n'
         '    else:\n'
         '        return _py_parse(stream, entry_rule_name, sa_err_listener, prediction_mode_sll)\n'),
        ('def _cpp_parse(stream:InputStream, entry_rule_name:str, sa_err_listener:SA_ErrorListener=None) -> ParseTree:',
         'def _cpp_parse(stream:InputStream, entry_rule_name:str, sa_err_listener:SA_ErrorListener=None,\n'
         '               prediction_mode_sll:bool=False) -> ParseTree:'),
        ('def _py_parse(stream:InputStream, entry_rule_name:str, sa_err_listener:SA_ErrorListener=None) -> ParseTree:',
         'def _py_parse(stream:InputStream, entry_rule_name:str, sa_err_listener:SA_ErrorListener=None,\n'
         '              prediction_mode_sll:bool=False) -> ParseTree:'),
        ('    parser = %(parser)s(token_stream)\n',
         '    parser = %(parser)s(token_stream)\n'
         '    if prediction_mode_sll:\n'
         '        from antlr4 import PredictionMode\n'
         '        parser._interp.predictionMode = PredictionMode.SLL\n'),
        ('do_parse(%(parser)s, stream, entry_rule_name, sa_err_listener)',
         'do_parse(%(parser)s, stream, entry_rule_name, sa_err_listener, prediction_mode_sll)'),
    ]),
    # An extension that imports is not necessarily one built from this grammar:
    # a stale in-place build after a grammar change would parse differently
    # without raising anything.
    ('the accelerator is used only when its grammar fingerprint matches', [
        ('try:\n'
         '    from . import %(shim)s_cpp_parser\n'
         'except ImportError:\n'
         '    USE_CPP_IMPLEMENTATION = False\n',
         'try:\n'
         '    from . import %(shim)s_cpp_parser\n'
         'except ImportError:\n'
         '    %(shim)s_cpp_parser = None\n'
         '\n'
         '# wwPDB: use the C++ parser only when it was compiled from the same grammar as the\n'
         '# Python parser, and record which path this grammar takes (see AccelUtil.py).\n'
         'from ..AccelUtil import acceptAccelerator\n'
         'USE_CPP_IMPLEMENTATION = acceptAccelerator(__name__, %(shim)s_cpp_parser, %(lexer)s, %(parser)s)\n'),
    ]),
]


def readerFacts(path: str) -> tuple:
    """ Extract (lexer, parser, entry rule, uses SLL) from one *Reader.py.
        @return: None when the file does not drive an ANTLR parse
    """

    with open(path, 'r', encoding='utf-8') as ifh:
        source = ifh.read()

    lexers = sorted(set(LEXER_IMPORT.findall(source)))
    parsers = sorted(set(PARSER_IMPORT.findall(source)))
    entryRules = sorted({rule for pattern in ENTRY_RULE_PATTERNS for rule in pattern.findall(source)})

    if not (lexers or parsers or entryRules):
        return None  # e.g. io/CifReader - not an ANTLR reader

    if not (len(lexers) == len(parsers) == len(entryRules) == 1):
        raise RuntimeError(f'{path}: expected exactly one lexer/parser/entry rule, got '
                           f'{lexers} / {parsers} / {entryRules}.')

    return lexers[0], parsers[0], entryRules[0], any(p.search(source) for p in SLL_PATTERNS)


def discoverGrammars() -> dict:
    """ Map grammar name -> GrammarSpec, by introspecting every *Reader.py.
    """

    found = {}

    for subPackage in SUB_PACKAGES:
        for path in sorted(glob.glob(os.path.join(NMR_DIR, subPackage, '*Reader.py'))):
            facts = readerFacts(path)
            if facts is None:
                continue
            (lexerSub, lexerName), (parserSub, parserName), entryRule, useSll = facts
            readerModule = os.path.basename(path)[:-3]
            name = parserName[:-len('Parser')]

            spec = found.get(name)
            if spec is None:
                found[name] = GrammarSpec(name, parserSub, parserName, lexerSub, lexerName,
                                          entryRule, [(readerModule, subPackage, useSll)], useSll)
                continue

            # Several Readers can drive one grammar (XMLParser has four,
            # XeasyPROTParser two). They must agree on the entry rule; they need
            # NOT agree on the prediction mode, which is why it is a per-call
            # argument rather than a compile-time constant.
            if (spec.lexerName, spec.entryRule) != (lexerName, entryRule):
                raise RuntimeError(f'{readerModule}: drives {parserName} with '
                                   f'({lexerName}, {entryRule}) but another Reader uses '
                                   f'({spec.lexerName}, {spec.entryRule}).')
            spec.readers.append((readerModule, subPackage, useSll))
            found[name] = spec._replace(useSll=spec.useSll or useSll)

    return found


def checkSpeedyAntlrVersion() -> None:
    """ Refuse an older speedy-antlr-tool, whose support library leaks the parse tree.
    """

    from importlib.metadata import version

    installed = version('speedy-antlr-tool')
    if tuple(int(part) for part in re.findall(r'\d+', installed)[:3]) < MIN_SPEEDY_ANTLR_VERSION:
        raise RuntimeError(f'speedy-antlr-tool {installed} is too old; '
                           f'{".".join(map(str, MIN_SPEEDY_ANTLR_VERSION))} or later is required.')


def fetchCppRuntime() -> None:
    """ Download and unpack the ANTLR4 C++ runtime source, if not already present.
    """

    if os.path.isfile(os.path.join(CPP_RUNTIME_DIR, 'antlr4-runtime.h')):
        return

    print(f'Downloading {CPP_RUNTIME_URL}')
    with urllib.request.urlopen(CPP_RUNTIME_URL) as response:  # nosec - fixed https URL
        payload = response.read()

    os.makedirs(CPP_RUNTIME_DIR, exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(payload)) as zfh:
        for member in zfh.namelist():
            if member.endswith('/'):
                continue
            if member == 'LICENSE.txt':
                target = os.path.join(CPP_RUNTIME_DIR, 'LICENSE.txt')
            elif member.startswith('runtime/src/'):
                target = os.path.join(CPP_RUNTIME_DIR, os.path.relpath(member, 'runtime/src'))
            else:
                continue
            os.makedirs(os.path.dirname(target), exist_ok=True)
            with zfh.open(member) as src, open(target, 'wb') as dst:
                shutil.copyfileobj(src, dst)

    print(f'ANTLR4 C++ runtime {ANTLR_VERSION} unpacked into {CPP_RUNTIME_DIR}')


def runAntlrCpp(spec: GrammarSpec) -> None:
    """ Generate the C++ target for a grammar's lexer/parser pair.
    """

    # 'python -m antlr4_tool_runner' would launch ANTLR's tree interpreter;
    # tool() is the entry point behind the 'antlr4' console script.
    command = [sys.executable, '-c', 'from antlr4_tool_runner import tool; tool()',
               '-v', ANTLR_VERSION,
               '-Dlanguage=Cpp', '-visitor', '-no-listener',
               '-o', CPP_SRC_DIR,
               spec.lexerName + '.g4', spec.parserName + '.g4']

    # Basenames keep generated headers independent of the checkout's absolute path.
    subprocess.run(command, cwd=GRAMMAR_DIR, check=True)


def applyPatches(path: str, patches: list, substitutions: dict) -> None:
    """ Apply exact-match patches to a generated file, in place.
    """

    with open(path, 'r', encoding='utf-8') as ifh:
        source = ifh.read()

    for description, replacements in patches:
        for old, new in replacements:
            old, new = old % substitutions, new % substitutions
            if new in source:
                continue  # already patched
            if source.count(old) != 1:
                raise RuntimeError(f'{os.path.basename(path)}: cannot apply "{description}" - '
                                   f'expected exactly one occurrence of\n{old!r}\nbut found '
                                   f'{source.count(old)}. speedy-antlr-tool output has changed.')
            source = source.replace(old, new)

    with open(path, 'w', encoding='utf-8') as ofh:
        ofh.write(source)


def generatedSources(spec: GrammarSpec) -> list:
    """ The C++ sources belonging to one accelerator, relative to cpp_src.

        Enumerated explicitly rather than globbed by prefix: NmrViewNPKParser,
        SparkyNPKParser and SparkyRPKParser borrow another grammar's lexer, so a
        '<parser prefix>*.cpp' glob would omit the lexer they need.
    """

    shim = f'sa_{spec.name.lower()}'
    sources = [f'{spec.lexerName}.cpp',
               f'{spec.parserName}.cpp',
               f'{spec.parserName}BaseVisitor.cpp',
               f'{spec.parserName}Visitor.cpp',
               f'{shim}_translator.cpp',
               f'{shim}_cpp_parser.cpp']

    missing = [name for name in sources if not os.path.isfile(os.path.join(CPP_SRC_DIR, name))]
    if missing:
        raise RuntimeError(f'{spec.name}: generated sources missing: {missing}')

    return sources


def generate(spec: GrammarSpec) -> dict:
    """ Generate the C++ accelerator for one grammar.
        @return: its manifest entry
    """

    import speedy_antlr_tool

    shim = f'sa_{spec.name.lower()}'
    readers = ', '.join(f'{module}[{"SLL" if sll else "LL"}]' for module, _, sll in spec.readers)
    print(f'{spec.name} -> {spec.subPackage}/{shim}.py  entry={spec.entryRule!r}  readers: {readers}')

    os.makedirs(CPP_SRC_DIR, exist_ok=True)
    runAntlrCpp(spec)

    # Read the committed Python parser; the shim is written alongside it.
    speedy_antlr_tool.generate(
        py_parser_path=os.path.join(NMR_DIR, spec.subPackage, spec.parserName + '.py'),
        cpp_output_dir=CPP_SRC_DIR,
        entry_rule_names=[spec.entryRule],
    )

    return patch(spec)


def patch(spec: GrammarSpec) -> dict:
    """ Apply the wwPDB patches to one grammar's generated C++ entry and shim.
        Every patch is idempotent, so this is safe on already patched sources.
        @return: its manifest entry
    """

    shim = f'sa_{spec.name.lower()}'
    assumedLexer = spec.name + 'Lexer'
    substitutions = {'parser': spec.parserName,
                     'lexer': spec.lexerName,
                     'assumedLexer': assumedLexer,
                     'shim': shim}

    cppEntry = os.path.join(CPP_SRC_DIR, f'{shim}_cpp_parser.cpp')
    shimPath = os.path.join(NMR_DIR, spec.subPackage, f'{shim}.py')

    applyPatches(cppEntry, CPP_PATCHES, substitutions)
    applyPatches(shimPath, SHIM_PATCHES, substitutions)

    if spec.lexerName != assumedLexer:
        applyPatches(cppEntry, BORROWED_LEXER_CPP_PATCHES, substitutions)
        applyPatches(shimPath, BORROWED_LEXER_SHIM_PATCHES, substitutions)
        print(f'  borrowed lexer: {assumedLexer} -> {spec.lexerName}')

    return {'module': f'{shim}_cpp_parser',
            'subpackage': spec.subPackage,
            'shim': shim,
            'grammar': spec.name,
            'lexer': spec.lexerName,
            'parser': spec.parserName,
            'entry_rule': spec.entryRule,
            'readers': [module for module, _, _ in spec.readers],
            'sources': generatedSources(spec)}


def writeManifest(entries: list) -> None:
    """ Record the built accelerators for setup.py, which must not depend on this
        script (the container deletes tools/ after building).
    """

    existing = {}
    if os.path.isfile(MANIFEST_PATH):
        with open(MANIFEST_PATH, 'r', encoding='utf-8') as ifh:
            existing = {entry['module']: entry for entry in json.load(ifh).get('extensions', [])}

    for entry in entries:
        existing[entry['module']] = entry

    manifest = {'antlr_version': ANTLR_VERSION,
                'generated_by': 'tools/gen_speedy_antlr.py',
                'shared_sources': ['speedy_antlr.cpp'],
                'extensions': [existing[key] for key in sorted(existing)]}

    with open(MANIFEST_PATH, 'w', encoding='utf-8') as ofh:
        json.dump(manifest, ofh, indent=2)
        ofh.write('\n')

    print(f'\n{os.path.relpath(MANIFEST_PATH, REPO_DIR)}: {len(manifest["extensions"])} accelerator(s)')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('grammars', nargs='*', help='grammar names, e.g. XplorMR')
    parser.add_argument('--all', action='store_true', help='generate every discovered grammar')
    parser.add_argument('--list', action='store_true', help='list discovered grammars and exit')
    parser.add_argument('--repatch', action='store_true',
                        help='re-apply the patches to the committed sources without regenerating them')
    args = parser.parse_args()

    grammars = discoverGrammars()

    if args.list:
        print(f'{"grammar":16s} {"pkg":4s} {"lexer":22s} {"entry rule":16s} {"mode":5s} readers')
        for name in sorted(grammars):
            spec = grammars[name]
            modes = {sll for _, _, sll in spec.readers}
            mode = 'mixed' if len(modes) > 1 else ('SLL' if spec.useSll else 'LL')
            print(f'{name:16s} {spec.subPackage:4s} '
                  f'{spec.lexerSubPackage + "/" + spec.lexerName:22s} {spec.entryRule:16s} '
                  f'{mode:5s} '
                  + ', '.join(f'{m}[{"SLL" if s else "LL"}]' for m, _, s in spec.readers))
        print(f'\n{len(grammars)} grammars driven by '
              f'{sum(len(s.readers) for s in grammars.values())} readers')
        return 0

    names = sorted(grammars) if args.all else args.grammars
    if not names:
        parser.error('name at least one grammar, or pass --all')

    unknown = [name for name in names if name not in grammars]
    if unknown:
        parser.error(f'unknown grammar(s): {", ".join(unknown)}; '
                     f'known: {", ".join(sorted(grammars))}')

    if args.repatch:
        entries = [patch(grammars[name]) for name in names]
        writeManifest(entries)
        return 0

    checkSpeedyAntlrVersion()
    fetchCppRuntime()
    entries = [generate(grammars[name]) for name in names]
    writeManifest(entries)

    print('\nNow build the accelerators:\n'
          '  WWPDB_NMR_BUILD_SPEEDY_ANTLR=1 python setup.py build_clib build_ext --inplace -j $(nproc)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
