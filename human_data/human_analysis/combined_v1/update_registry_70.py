"""Registry-only excerpt of original stage1_delivery/build_delivery.py. No inference.
Unrelated Word/Excel reports and external genetics retry dependencies are not executed.
Original H4 gate prose is historical; no H4/cognition result is created.
"""
from pathlib import Path
import sys,csv
_home=next(p for p in Path(__file__).resolve().parents if (p/'release_support.py').exists())
sys.path.insert(0,str(_home))
from release_support import gate
gate(__file__,__name__)
P=_home;D=P/'human_analysis/combined_v1/integrated_registry';C=P/'human_analysis/combined_v1';H=P/'human_analysis/h2_v1'
D.mkdir(parents=True,exist_ok=True)
def tsv(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(name,rows):
    cols=list(dict.fromkeys(k for r in rows for k in r))
    with (D/name).open('x') as f:
        w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def n(r,k):return float(r[k])
oldreg=tsv(C/'results/registry_tracking_62.tsv');family=tsv(C/'results/family_summary.tsv')
h2=tsv(H/'h2_results_long.tsv');base=[r for r in h2 if r['family']=='H2_BASE'];length=[r for r in h2 if r['family']=='H2_LENGTH']
tracking=[]
for old in oldreg:
    r=dict(old);r['original_family']=r['family'];r['original_status']=r['status'];r['delivery_gate']='source_status_retained'
    if r['test_id'].startswith('H2_'):
        x=next(z for z in base if z['test_id']==r['test_id'])
        r.update(family='H2_BASE',family_m='8',status=x['status'],cohort='SEAAD_HIP_same40',target='Glut_vs_three_class_expression_reference',gene_set=x['set_name'],reason='',raw_P=x['raw_p'],adjusted_P=x['q_BH'],p_for_adjustment_internal=x['p_for_family_adjustment'],adjustment_status='complete_fixed_family',source_result='h2_v1/h2_results_long.tsv',position_resolution='tested',delivery_gate='coverage_and_matching_PASS',exploratory_discovery_status=('positive_matched_background_enrichment' if n(x,'Delta')>0 and n(x,'q_BH')<.05 else 'does_not_meet_enrichment_rule'))
    elif r['test_id']=='H4_01':
        r['delivery_gate']='conditional_input_retry; one Pericytes chr21 summary integrity PASS, analysis inputs NOT READY; no genetics/MR result'
        r['delivery_note']='H4 pending, not executed. Pericytes chr21 summary integrity PASS; companion annotation prefix only65536B/HTTP206 shows EA/OA,hg38/hg19,MAF (no full CRC). Author A1 effect, pooled combined_7 MAF and normal-approx derived-SE code known; no joining/harmonization or SE computation. External Pericytes150 is not verified per-pair runtime N. Final N/df, variance/frequency strategy and H4 plan unfrozen. coloc.abf itself does not require LD.'
    elif 'pending' in r['status']: r['delivery_gate']='dictionary_and_time_precision_BLOCKED'
    tracking.append(r)
for x in length:
    tracking.append({'test_id':x['test_id'],'hypothesis':'H2_length_sensitivity','cohort':'SEAAD_HIP_same40','target':'Glut_vs_three_class_expression_reference','gene_set':x['set_name'],'family':'H2_LENGTH','family_m':'8','status':x['status'],'reason':x['reason'],'source_result':'h2_v1/h2_results_long.tsv','position_resolution':'tested','raw_P':x['raw_p'],'p_for_adjustment_internal':x['p_for_family_adjustment'],'adjusted_P':x['q_BH'],'adjustment_status':'complete_fixed_family','primary_support_status':'not_applicable','exploratory_discovery_status':('positive_matched_background_enrichment' if n(x,'Delta')>0 and n(x,'q_BH')<.05 else 'does_not_meet_enrichment_rule'),'original_family':'new_separate_H2_LENGTH','original_status':'registered_in_H2_appendix','delivery_gate':'coverage_length_and_matching_PASS'})
assert len(tracking)==70 and len({r['test_id'] for r in tracking})==70
pending=[r for r in tracking if r['position_resolution']=='unknown_pending']; terminal=[r for r in tracking if r['position_resolution']=='terminal_NA_preserve_slot']
assert len(pending)==13 and len(terminal)==5
assert sum(bool(r['raw_P']) for r in tracking)==52
write_tsv('registry_tracking_70.tsv',tracking)
newfamilies=[dict(r) for r in family if r['family']!='H2_FUTURE']
for f in ['H2_BASE','H2_LENGTH']:
    newfamilies.append({'family':f,'fixed_m':'8','registered_positions':'8','tested_n':'8','terminal_NA_positions_n':'0','unknown_pending_n':'0','adjustment_status':'complete_fixed_family','method':'fdr_bh','unknown_P_not_imputed':'True','terminal_NA_internal_P1':'0'})
write_tsv('family_summary_current.tsv',newfamilies)

