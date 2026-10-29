"""Synthetic figure/table/report smoke test only, no real H2 scores or matrices."""

# Release-only safety/path adapter; see provenance/adaptation_diffs.
from pathlib import Path as _ReleasePath
import sys as _release_sys
_release_home = next(p for p in _ReleasePath(__file__).resolve().parents if (p / 'release_support.py').is_file())
_release_sys.path.insert(0, str(_release_home))
from release_support import gate, release_root, resolve_source
gate(__file__, __name__)
from pathlib import Path
import sys,json,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
import pandas as pd
from run_h2_v1 import plot_results,write_report,SET_ORDER
from h2_core import empirical_summary,bh_fixed
O=Path(__file__).resolve().parents[1]/'readiness/synthetic_output_smoke'
O.mkdir(exist_ok=True)
rows=[]
for family,prefix in [('H2_BASE','H2'),('H2_LENGTH','H2L')]:
    for i,name in enumerate(SET_ORDER):
        # Deterministic artificial reference values, not draws from any real gene set.
        null=np.linspace(.15+.02*i,.5+.02*i,21)
        row={'test_id':f'{prefix}_{i+1:02}','family':family,'set_name':name,'status':'EXECUTED'}
        row.update(empirical_summary(.3+.04*i,null));rows.append(row)
r=pd.DataFrame(rows);r['q_BH']=np.nan
for family in ['H2_BASE','H2_LENGTH']:
    mask=r.family==family;r.loc[mask,'q_BH']=bh_fixed(r.loc[mask,'raw_p'],8)
r.to_csv(O/'SYNTHETIC_results.tsv',sep='\t',index=False)
plot_results(r,O/'SYNTHETIC_figure',synthetic=True)
write_report(r,{'status':'SYNTHETIC_TEMPLATE_QA_ONLY_NOT_HUMAN_RESULTS'},O/'SYNTHETIC_report_template.md')
report=(O/'SYNTHETIC_report_template.md').read_text()
# Guard synthetic output from being mistaken for study results.
(O/'SYNTHETIC_report_template.md').write_text('# SYNTHETIC TEMPLATE SMOKE TEST — NOT HUMAN RESULTS\n\n'+report)
assert (O/'SYNTHETIC_figure.pdf').stat().st_size>1000
assert (O/'SYNTHETIC_figure.png').stat().st_size>1000
assert 'NOT confidence intervals' in report or 'NOT effect' in report
assert len(pd.read_csv(O/'SYNTHETIC_results.tsv',sep='\t'))==16
print('PASS: 16-row synthetic table/report and labelled English PDF/PNG; no real data read.')
