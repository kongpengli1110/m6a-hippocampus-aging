#!/usr/bin/env python3
"""Small explicit stage launcher. Default is read-only input/overwrite inspection."""
import argparse,json,subprocess,sys
from pathlib import Path
from release_support import release_root,resolve_source

def main():
    p=argparse.ArgumentParser(description=__doc__)
    stages=json.loads((release_root()/'stages.json').read_text())
    p.add_argument('stage',nargs='?',choices=list(stages));p.add_argument('--list',action='store_true')
    mode=p.add_mutually_exclusive_group();mode.add_argument('--dry-run',action='store_true');mode.add_argument('--execute',action='store_true')
    a=p.parse_args()
    if a.list or not a.stage:
        for k,s in stages.items():print(k,s['kind'],s['script'])
        return
    s=stages[a.stage];root=release_root()
    inputs=[resolve_source('/workspace/projects/m6a-aging-paper/'+r) if a.stage in ('h2','h2-lock') else root/r for r in s['inputs']]
    missing=[str(x.relative_to(root)) for x in inputs if not x.is_file()]
    existing=[x for x in s['outputs'] if (root/x).exists()]
    command=[sys.executable,str(root/s['script']),'--execute',*s['args']]
    print(json.dumps(dict(stage=a.stage,cwd=str(root),command=command,missing_inputs=missing,existing_outputs=existing,mode='execute' if a.execute else 'dry-run',note=s['note']),indent=2))
    if not a.execute:return
    if missing or existing:p.error('Execution blocked: missing inputs or existing outputs. Preserve prior outputs; use a fresh scratch copy. No automatic overwrite.')
    for d in ['human_analysis/gse48350_v1/logs','human_analysis/gse48350_v1/results','human_analysis/gse48350_v1/figures','human_analysis/sea_hip_v1/logs','human_analysis/sea_hip_v1/derived','human_analysis/sea_hip_v1/results','human_analysis/sea_hip_v1/figures','human_analysis/h2_v1/readiness','human_data/sea_hip/logs','human_data/sea_hip/provenance']:(root/d).mkdir(parents=True,exist_ok=True)
    subprocess.run(command,cwd=root,check=True)
if __name__=='__main__':main()
