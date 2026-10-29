"""Packaging-only minimal report after successful actual aggregation.
Not the historical full report; does not claim old downloader/integration tests ran.
"""
from pathlib import Path
import hashlib,json

def main():
    root=Path(__file__).resolve().parents[1]
    q=json.loads((root/'QA.json').read_text())
    if q['status']!='PASS':raise RuntimeError('Aggregation QA did not pass')
    report=root/'QA_report.md'
    if report.exists():raise RuntimeError('Refuse to overwrite report')
    report.write_text('# Aggregation QA (release adapter)\n\nGenerated solely from this aggregation QA.json. No historical download-resume or integration test is claimed.\n\n```json\n'+json.dumps(q,indent=2)+'\n```\n')
    p=root/'output_manifest.json';m=json.loads(p.read_text());m['QA_report.md']={'bytes':report.stat().st_size,'sha256':hashlib.sha256(report.read_bytes()).hexdigest()};p.write_text(json.dumps(m,indent=2)+'\n')
if __name__=='__main__':
    import sys
    if '--execute' not in sys.argv or '--help' in sys.argv or '--dry-run' in sys.argv:print(__doc__+'\nUse --execute only after completed aggregation QA.')
    else:main()
