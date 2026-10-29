"""Join frozen registry + completed GSE results + new SEA/S16 results, no refits."""

# Release-only safety/path adapter; see provenance/adaptation_diffs.
from pathlib import Path as _ReleasePath
import sys as _release_sys
_release_home = next(p for p in _ReleasePath(__file__).resolve().parents if (p / 'release_support.py').is_file())
_release_sys.path.insert(0, str(_release_home))
from release_support import gate, release_root, resolve_source
gate(__file__, __name__)
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,sys
import numpy as np,pandas as pd,statsmodels
from statsmodels.stats.multitest import multipletests
ROOT=Path(__file__).resolve().parents[1];AN=ROOT.parent;G=AN/'gse48350_v1';C=AN/'combined_v1';R=C/'results';L=C/'logs'
for p in [R,L,C/'figures']:p.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,p):p.write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))
def manual_adjust(p,method):
 p=np.asarray(p,float);m=len(p);order=np.argsort(p,kind='stable');s=p[order]
 if method=='holm':a=np.minimum(1,np.maximum.accumulate((m-np.arange(m))*s))
 else:a=np.minimum(1,np.minimum.accumulate((m*s/(np.arange(m)+1))[::-1])[::-1])
 out=np.empty(m);out[order]=a;return out
# Fixed-slot including terminal NA, vs shrinking m, unit test before real adjustment.
for method in ['holm','fdr_bh']:
 p=np.array([.01,.04,1.]);assert np.allclose(manual_adjust(p,method),multipletests(p,method=method)[1])
assert np.allclose(manual_adjust([.01,.04,1.],'holm'),[.03,.08,1.])
reg=pd.read_csv(AN/'prespecified_tests_v1_frozen.tsv',sep='\t',dtype=str,keep_default_na=False);old=pd.read_csv(G/'results/registered_results_long.tsv',sep='\t');old=old[old.test_id.isin(['P01','E01','D01','S07','S08'])].copy();new=pd.read_csv(ROOT/'results/registered_results_long.tsv',sep='\t');assert len(old)==5 and len(new)==36
assert json.loads((ROOT/'logs/execution_PASS.json').read_text())['status']=='PASS';assert sha(ROOT/'results/registered_results_long.tsv')==json.loads((ROOT/'logs/execution_PASS.json').read_text())['results_SHA256']
src={**{r.test_id:('GSE_original',r._asdict()) for r in old.itertuples(index=False)},**{r.test_id:('SEA_v1_or_S16',r._asdict()) for r in new.itertuples(index=False)}};assert len(src)==41
rows=[]
for _,rr in reg.iterrows():
 tid=rr.test_id;row=rr.to_dict()
 if tid in src:
  origin,result=src[tid]
  # Verify registry identity fields were not changed by either execution.
  for k in ['family','cohort','target','gene_set','term_tested','model']:
   assert str(result[k])==str(rr[k]),(tid,k,result[k],rr[k])
  row.update(result);row['source_result']=origin
 else:
  row.update(status='pending_dictionary_and_time_QC' if rr.family=='E_COG' or tid in ['S19','S20','S21','S22'] else ('registered_not_executable' if rr.family=='H2_FUTURE' else 'conditional_later_pending'),reason=rr.stage,source_result='not_executed',raw_P=np.nan)
 status=row['status'];p=row.get('raw_P',np.nan)
 if status.startswith('tested'):
  assert np.isfinite(p) and 0<=p<=1;row['p_for_adjustment_internal']=p;row['position_resolution']='tested'
 elif status.startswith('not_testable') or status=='identical_prespecified_subset':
  assert pd.isna(p);row['p_for_adjustment_internal']=1.;row['position_resolution']='terminal_NA_preserve_slot'
 else:
  assert pd.isna(p);row['p_for_adjustment_internal']=np.nan;row['position_resolution']='unknown_pending'
 row['adjusted_P']=np.nan;row['adjustment_status']='pending-family';row['family_method']='';row['primary_support_status']='not_applicable';row['exploratory_discovery_status']='not_applicable';rows.append(row)
out=pd.DataFrame(rows);assert len(out)==62 and out.test_id.is_unique
summary=[]
for family,g in out.groupby('family',sort=False):
 fixed=reg.loc[reg.family.eq(family),'family_m'].unique();assert len(fixed)==1;fixed_m=int(fixed[0]) if fixed[0].isdigit() else None
 if fixed_m is not None:assert len(g)==fixed_m
 unknown=int(g.position_resolution.eq('unknown_pending').sum());terminal=int(g.position_resolution.eq('terminal_NA_preserve_slot').sum());tested=int(g.position_resolution.eq('tested').sum());method='holm' if family=='PRIMARY' else 'fdr_bh'
 complete=family in ['PRIMARY','E_EXPR','E_MODULE','E_AD'];assert not complete or unknown==0
 if complete:
  p=g.p_for_adjustment_internal.to_numpy();a=manual_adjust(p,method);assert np.allclose(a,multipletests(p,method=method)[1],atol=1e-14)
  for idx,adj in zip(g.index,a):
   out.loc[idx,'adjustment_status']='complete_fixed_family';out.loc[idx,'family_method']=method
   # Public NA for an untested terminal slot, even if internal adjusted1.
   if out.loc[idx,'position_resolution']=='tested':out.loc[idx,'adjusted_P']=adj
   if family=='PRIMARY':out.loc[idx,'primary_support_status']='meets_negative_Holm_rule' if out.loc[idx,'position_resolution']=='tested' and out.loc[idx,'beta']<0 and adj<.05 else 'does_not_meet_primary_support_rule'
   else:out.loc[idx,'exploratory_discovery_status']='BH_q_below_0.05_exploratory_only' if out.loc[idx,'position_resolution']=='tested' and adj<.05 else 'does_not_meet_exploratory_FDR_rule'
 summary.append({'family':family,'fixed_m':fixed_m,'registered_positions':len(g),'tested_n':tested,'terminal_NA_positions_n':terminal,'unknown_pending_n':unknown,'adjustment_status':'complete_fixed_family' if complete else 'pending-family','method':method if complete else 'not_applied','unknown_P_not_imputed':True,'terminal_NA_internal_P1':terminal})
out.to_csv(R/'combined_results_long.tsv',sep='\t',index=False);out[out.source_result.ne('not_executed')].to_csv(R/'completed_positions_long.tsv',sep='\t',index=False);out[['test_id','hypothesis','cohort','target','gene_set','family','family_m','status','reason','source_result','position_resolution','raw_P','p_for_adjustment_internal','adjusted_P','adjustment_status','primary_support_status','exploratory_discovery_status']].to_csv(R/'registry_tracking_62.tsv',sep='\t',index=False);pd.DataFrame(summary).to_csv(R/'family_summary.tsv',sep='\t',index=False)
for family in ['PRIMARY','E_EXPR','E_MODULE','E_AD','E_SENS']:out[out.family.eq(family)].to_csv(R/(family+'_results.tsv'),sep='\t',index=False)
# Original GSE numerical inference remains exactly unchanged, not recomputed.
for _,r in old.iterrows():
 dest=out.set_index('test_id').loc[r.test_id]
 for col in ['n','beta','SE_HC3','df','CI95_low','CI95_high','raw_P']:assert np.isclose(float(dest[col]),float(r[col]),rtol=0,atol=0),(r.test_id,col)
assert out[out.family.eq('E_SENS')].adjusted_P.isna().all();assert out[out.test_id.isin(['S19','S20','S21','S22'])].p_for_adjustment_internal.isna().all();assert out.position_resolution.value_counts().to_dict()=={'tested':36,'unknown_pending':21,'terminal_NA_preserve_slot':5}
input_paths=[AN/'human_analysis_plan_v1.md',AN/'prespecified_tests_v1_frozen.tsv',AN/'freeze_manifest.json',G/'results/registered_results_long.tsv',G/'logs/final_QA.json',ROOT/'results/registered_results_long.tsv',ROOT/'logs/execution_PASS.json'];dump({'files':{str(p):sha(p) for p in input_paths},'no_model_refits':True},C/'input_manifest.json');dump({'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'statsmodels':statsmodels.__version__,'assembly_script':str(Path(__file__)),'assembly_script_SHA256':sha(Path(__file__))},C/'software_versions.json')
dump({'status':'PASS','completed_at_utc':datetime.now(timezone.utc).isoformat(),'registry_n':62,'resolved_positions':41,'unknown_pending':21,'terminal_NA':5,'GSE_original_numeric_results_unchanged':True,'fixed_m_manual_vs_statsmodels':'PASS','E_SENS_not_adjusted':True,'outputs':{p.name:sha(p) for p in R.glob('*')}},L/'assembly_PASS.json');print(pd.DataFrame(summary).to_string(index=False));print(out[out.family.isin(['PRIMARY','E_EXPR','E_MODULE'])][['test_id','family','target','gene_set','beta','CI95_low','CI95_high','raw_P','adjusted_P','primary_support_status','exploratory_discovery_status']].to_string(index=False))
