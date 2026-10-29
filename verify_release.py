#!/usr/bin/env python3
"""Read-only payload SHA256, AST/compile and recorded source-copy audit. No analysis imports."""
from pathlib import Path
import argparse,ast,hashlib,json

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--compile',action='store_true',help='also parse/compile all Python files in memory, without importing them');a=p.parse_args()
    human=Path(__file__).resolve().parent;root=human.parent
    manifest=root/'SHA256SUMS'
    if not manifest.is_file():p.error('Run against the final extracted release, which contains SHA256SUMS.')
    n=0
    for line in manifest.read_text().splitlines():
        expected,rel=line.split('  ',1);file=root/rel
        if not file.is_file() or digest(file)!=expected:raise RuntimeError('Payload mismatch: '+rel)
        n+=1
    source=json.loads((human/'provenance/source_map.json').read_text())
    for r in source:
        if digest(root/r['release_relative'])!=r['release_sha256']:raise RuntimeError('Source-copy map mismatch: '+r['release_relative'])
    compiled=0
    if a.compile:
        for file in root.rglob('*.py'):
            tree=ast.parse(file.read_text(),filename=str(file));compile(tree,str(file),'exec');compiled+=1
    print(json.dumps({'status':'PASS','payload_files_checked':n,'python_files_compiled_without_import':compiled,'source_copy_records':len(source),'cohort_models_or_nulls_executed':False},indent=2))
if __name__=='__main__':main()
