#!/usr/bin/env python3
"""Inspect exact raw inputs by default. Network requires --download --execute."""
import argparse,hashlib,json,os,urllib.request
from pathlib import Path

def main():
    root=Path(__file__).resolve().parent
    rows=json.loads((root/'input_manifest.json').read_text())['external_and_raw']
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--id',choices=[r['id'] for r in rows]);p.add_argument('--download',action='store_true');p.add_argument('--execute',action='store_true');p.add_argument('--check',action='store_true')
    a=p.parse_args();selected=[r for r in rows if not a.id or r['id']==a.id]
    if a.download and (not a.execute or not a.id):p.error('Download requires both --id and --execute; no bulk downloads.')
    for r in selected:
        dest=root/r['path'];print(json.dumps({**r,'local_exists':dest.exists()},indent=2))
        if a.check:
            if not dest.is_file():print('MISSING (expected for omitted raw input)');continue
            if dest.stat().st_size!=r['bytes']:raise RuntimeError('Size mismatch: '+r['id'])
            h=hashlib.sha256()
            with dest.open('rb') as f:
                for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
            if h.hexdigest()!=r['sha256']:raise RuntimeError('Checksum mismatch: '+r['id'])
            print('SHA256 PASS')
        if not a.download:continue
        if r['delivery_status'] not in ['external_raw_omitted','external_annotation_omitted']:p.error('This input is bundled, derived or not needed. No download performed.')
        if dest.exists():p.error('Refuse to overwrite existing data.')
        dest.parent.mkdir(parents=True,exist_ok=True);part=dest.with_name(dest.name+'.part')
        if part.exists():p.error('Prior partial file exists; preserve and review it. No automatic overwrite/resume.')
        h=hashlib.sha256();n=0
        with urllib.request.urlopen(r['url'],timeout=120) as response,part.open('xb') as out:
            while b:=response.read(4*1024*1024):
                n+=len(b)
                if n>r['bytes']:raise RuntimeError('Response exceeds exact frozen size; partial retained')
                out.write(b);h.update(b)
        if n!=r['bytes'] or h.hexdigest()!=r['sha256']:raise RuntimeError('Frozen identity mismatch; partial retained; do not use')
        part.rename(dest);print('Downloaded exact verified object; NO analysis executed.')
if __name__=='__main__':main()
