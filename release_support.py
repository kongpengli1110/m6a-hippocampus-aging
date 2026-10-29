"""Packaging-only safety and old-provenance path resolution; no statistics."""
from pathlib import Path
import json, os, sys
OLD_ROOT = '/workspace/projects/m6a-aging-paper'
def release_root():
    return Path(__file__).resolve().parent

def resolve_source(value):
    """Resolve original absolute evidence paths WITHOUT editing frozen JSON bytes.
    Explicit map entries are old source-code/QA attestations, not fresh results.
    Newly generated absolute paths are returned as-is.
    """
    text=str(value)
    if text.startswith(OLD_ROOT + '/'):
        rel=text[len(OLD_ROOT)+1:]
        mapping=json.loads((release_root()/'provenance/frozen_path_map.json').read_text())
        return release_root()/mapping.get(rel,rel)
    return Path(value)

def gate(filename, module_name):
    # The H2 runner has function definitions only and is imported by the toy template test.
    if module_name != '__main__':
        if Path(filename).name == 'run_h2_v1.py': return
        raise RuntimeError('Analysis entrypoints may not be imported. Use run.py or AST extraction.')
    if '--help' in sys.argv or '-h' in sys.argv:
        print('Release safety adapter: --help / --dry-run never import scientific packages or run analysis.')
        print('Use: python run.py <stage> [--dry-run | --execute]. See stages.json and README.')
        print('Direct execution requires --execute; H2 also needs --validate-lock-only or --execute-locked-h2.')
        raise SystemExit(0)
    if '--dry-run' in sys.argv:
        print('DRY RUN ONLY:', Path(filename).name, '; no analysis imports or writes.')
        raise SystemExit(0)
    if '--execute' not in sys.argv:
        raise SystemExit('Not executed. Explicit --execute is required; prefer run.py for input/overwrite checks.')
    sys.argv.remove('--execute')
