"""Execute only authorized GSE registered tests: P01 E01 D01 S07 S08.
S16 remains pending actual SEA gene-A intersection. No network, gene scans,
SEA models, mouse reanalysis, cognition or adaptive/model/probe selection.
"""

# Release-only safety/path adapter; see provenance/adaptation_diffs.
from pathlib import Path as _ReleasePath
import sys as _release_sys
_release_home = next(p for p in _ReleasePath(__file__).resolve().parents if (p / 'release_support.py').is_file())
_release_sys.path.insert(0, str(_release_home))
from release_support import gate, release_root, resolve_source
gate(__file__, __name__)
from pathlib import Path
from datetime import datetime,timezone
import gzip,json,hashlib,re,sys,platform,math
import numpy as np,pandas as pd,scipy,statsmodels,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from fixed_pipeline import *

ROOT=Path(__file__).resolve().parents[1];PROJECT=ROOT.parents[1];A=ROOT.parent;PRE=PROJECT/'human_preflight';F=PROJECT/'handoff/frozen_v1';RES=ROOT/'results';LOG=ROOT/'logs';FIG=ROOT/'figures'
for p in [RES,LOG,FIG]:p.mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def safe(o):
 if isinstance(o,np.generic):return o.item()
 if isinstance(o,Path):return str(o)
 raise TypeError(type(o).__name__)
def dump(obj,p):p.write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=safe))
def tab(df,name):df.to_csv(RES/name,sep='\t',index=False)
def now():return datetime.now(timezone.utc).isoformat()
if (LOG/'run_complete.json').exists():raise SystemExit('Completed run exists: do not repeat successful analysis or overwrite results.')
started=now();freeze=json.loads((A/'freeze_manifest.json').read_text());assert freeze['status']=='INTERNALLY_PROSPECTIVELY_LOCKED';assert datetime.fromisoformat(started)>datetime.fromisoformat(freeze['lock_recorded_at'])
for r in freeze['frozen_files']:assert sha(resolve_source(r['path']))==r['sha256'],r['path']
assert hashlib.md5(resolve_source(freeze['candidate_xlsx']['path']).read_bytes()).hexdigest()==freeze['candidate_xlsx']['md5']
for line in (F/'FROZEN_md5.txt').read_text().splitlines():
 if re.match(r'^[a-f0-9]{32}\s',line):
  md5,n=line.split(maxsplit=1);p=F/Path(n).name;assert hashlib.md5(p.read_bytes()).hexdigest()==md5
synthetic=json.loads((LOG/'synthetic_validation.json').read_text());assert synthetic['status']=='PASS';assert datetime.fromisoformat(synthetic['completed_at_utc'])<datetime.fromisoformat(started)
source_paths=[A/'freeze_manifest.json',*[resolve_source(r['path']) for r in freeze['frozen_files']],*sorted(F.glob('*')),PRE/'raw/GSE48350_series_matrix.txt.gz',PRE/'raw/download_manifest.json',PRE/'human_candidate_crosswalk.tsv',PRE/'derived/GPL570_probe_hgnc_crosswalk.tsv',PRE/'derived/GSE48350_sample_manifest.tsv',PRE/'derived/GSE48350_matrix_scale.tsv',PRE/'logs/final_validation.json']
source_paths=[p for p in source_paths if p.is_file()];source_hashes={str(p):sha(p) for p in source_paths}
manifest=json.loads((PRE/'raw/download_manifest.json').read_text());expected=next(r for r in manifest if r['file']=='GSE48350_series_matrix.txt.gz');assert source_hashes[str(PRE/'raw/GSE48350_series_matrix.txt.gz')]==expected['sha256']
pv=json.loads((PRE/'logs/final_validation.json').read_text());assert source_hashes[str(PRE/'human_candidate_crosswalk.tsv')]==pv['outputs']['human_candidate_crosswalk.tsv']['sha256']
params={'version':'GSE48350_v1','started_at_utc':started,'prospective_lock':freeze['lock_recorded_at'],'authorized_test_ids':['P01','E01','D01','S07','S08'],'pending_GSE_test_ids':['S16'],'primary_scale':'log2_ratio','primary_log2_count':1,'pseudocount':0,'quantile_normalization':False,'probe_aggregation':'median_after_scale_transform_all_single_identity_probes','reference':'normal43','gene_sd_ddof':1,'score_sd_ddof':1,'SD_epsilon':1e-8,'coverage_fraction':.8,'minimum_genes':10,'weighting':'equal','age10':'(age-80)/10','sex':{'female':0,'male':1},'AD_age_range_inclusive':[60,94],'OLS_covariance':'HC3','inference':'two-sided Student t df=n-rank(X); pointwise 95% CI','design_n_gate':'max(20,3p,p+15)','age_unique_min':5,'age_range_min':10,'sex_min_each':3,'VIF_age_flag':10,'hat_flag':'2p/n','Cook_flag':'4/n','LOO':'fixed-score beta only, no P/no selection','multiplicity':{'PRIMARY':'pending-family fixed Holm m3','E_EXPR':'pending-family fixed BH m11','E_SENS':'pending-family fixed BH m23','E_AD':'BH m1 complete'},'actual_expression_GSE':'not_assessable_from_ratio_matrix'}
versions={'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'statsmodels':statsmodels.__version__,'matplotlib':matplotlib.__version__,'scripts_sha256':{p.name:sha(p) for p in sorted((ROOT/'scripts').glob('*.py'))}}
dump(params,ROOT/'run_parameters.json');dump(versions,ROOT/'software_versions.json');dump({'checked_at_utc':now(),'files':source_hashes},ROOT/'input_manifest.json')
registry=pd.read_csv(A/'prespecified_tests_v1_frozen.tsv',sep='\t',keep_default_na=False);reg=registry[registry.cohort.eq('GSE48350')].set_index('test_id');assert set(reg.index)=={'P01','E01','D01','S07','S08','S16'}
meta=pd.read_csv(PRE/'derived/GSE48350_sample_manifest.tsv',sep='\t');hip=meta[meta.is_hippocampus].copy();normal=hip[hip.primary_age_manifest].copy();ad=hip[hip.diagnosis.eq('AD')]
assert len(hip)==62 and len(normal)==43 and len(ad)==19;assert normal.diagnosis.eq('normal').all();assert normal.sex.value_counts().to_dict()=={'male':23,'female':20}
# Explicit diagnosis namespace follows GEO's distinct indiv versus AD labels.
hip['donor_id']=hip.apply(lambda r: 'normal_indiv_'+re.search(r'_indiv(.+)$',r.title).group(1) if r.diagnosis=='normal' else 'AD_'+re.search(r'_AD_(.+)$',r.title).group(1),axis=1)
assert hip.donor_id.is_unique and hip['sample'].is_unique
for _,r in hip.iterrows():
 if r.diagnosis=='normal':assert r.individual.split(',')[0].strip()==r.donor_id.removeprefix('normal_indiv_')
hip['sex_label']=hip.sex;hip['sex']=hip.sex.map({'female':0,'male':1});assert hip.sex.notna().all();hip['age10']=(hip.age-80)/10;hip['AD']=hip.diagnosis.eq('AD').astype(int);hip['D01_population']=hip.age.between(60,94)
assert max(normal.age.min(),ad.age.min())==60 and min(normal.age.max(),ad.age.max())==94
hip=hip.set_index('sample');ref=hip.index[hip.primary_age_manifest].tolist();normal_ids=ref;tab(hip.reset_index(),'sample_manifest.tsv')
# Mapping is REUSED, never re-inferred. Validate listed single-identity probes.
x=pd.read_csv(PRE/'human_candidate_crosswalk.tsv',sep='\t',keep_default_na=False);probe_map=pd.read_csv(PRE/'derived/GPL570_probe_hgnc_crosswalk.tsv',sep='\t',keep_default_na=False).set_index('probe_id')
sets={'core23_1to1':x[x.set_kind.eq('core')].copy(),'METTL3cKO_HIP_down':x[x.set_name.eq('METTL3cKO_HIP_down')].copy()};maps={};allprobes=set();genemeta={}
for name,g in sets.items():
 g0=g[g.human_mapping_unambiguous];assert g0.current_HGNC_IDs.is_unique;assert len(g0)==(22 if name=='core23_1to1' else 120)
 tech=g0[g0.GPL570_primary_feature_available];gm={r.current_HGNC_IDs:r.GPL570_unique_probe_ids.split(';') for r in tech.itertuples()};maps[name]=gm;genemeta[name]=g.set_index('current_HGNC_IDs',drop=False)
 for hid,probes in gm.items():
  assert len(probes)==len(set(probes)) and probes
  for p in probes:
   pr=probe_map.loc[p];assert bool(pr.unambiguous) and pr.current_HGNC_IDs==hid and not pr.stable_symbol_conflict and not pr.has_unresolved_entrez
  allprobes.update(probes)
# Read the already downloaded matrix once, only HIP columns. No re-download.
with gzip.open(PRE/'raw/GSE48350_series_matrix.txt.gz','rt') as f:
 for line in f:
  if line.startswith('!series_matrix_table_begin'):break
 matrix=pd.read_csv(f,sep='\t',comment='!',index_col=0,usecols=['ID_REF',*hip.index])
assert matrix.shape==(54675,62) and matrix.index.is_unique
m=matrix.to_numpy(dtype=float);assert np.isfinite(m).all() and (m>0).all()
scale=pd.read_csv(PRE/'derived/GSE48350_matrix_scale.tsv',sep='\t').set_index('scope').loc['hippocampus_all'];assert np.isclose(m.min(),scale['min'],rtol=0,atol=1e-12) and np.isclose(m.max(),scale['max'],rtol=0,atol=1e-10)
source_ratio=matrix.loc[sorted(allprobes),hip.index].copy();source_ratio.to_csv(RES/'eligible_probe_source_ratios.tsv.gz',sep='\t',index=True);original_ratios=source_ratio.copy(deep=True);del matrix,m
scores={};gene_selection=[];probe_rows=[];gene_refs=[];score_refs=[];coverage=[];score_records=[];stage_values={}
for name,g in sets.items():
 gm=maps[name];g0=g[g.human_mapping_unambiguous];samples=list(hip.index) if name=='core23_1to1' else ref;main_A=None
 for input_scale in ['log2_ratio','raw_ratio']:
  transformed,logcount=transform_probes(source_ratio.loc[:,samples],scale=input_scale,input_log2_count=0);assert logcount==(1 if input_scale=='log2_ratio' else 0)
  gv=aggregate_probes(transformed,gm);out=standardize_fixed_A(gv,ref,fixed_A=main_A if input_scale=='raw_ratio' else None)
  if input_scale=='log2_ratio':main_A=out['A']
  else:assert out['A']==main_A
  pass_cov,min_genes=coverage_gate(len(out['A']),len(g0));out['coverage_pass']=pass_cov;out['minimum_genes']=min_genes
  scores[(name,input_scale)]=out;stage_values[(name,input_scale)]=gv
  coverage.append({'gene_set':name,'scale':input_scale,'source_member_n':len(g),'G0_one_to_one_n':len(g0),'technical_available_n':len(gm),'actual_expression_n':np.nan,'actual_expression_status':'not_assessable_from_ratio_matrix','reference_variance_pass_n':int(out['variance_pass'].sum()),'A_n':len(out['A']),'number_coverage':len(out['A'])/len(g0),'weight_coverage':len(out['A'])/len(g0),'minimum_genes':min_genes,'coverage_gate_pass':pass_cov,'score_status':out['status'],'reference_n':len(ref),'gene_A_sha256':hashlib.sha256('\n'.join(out['A']).encode()).hexdigest(),'source_log2_count':0,'applied_log2_count':logcount,'reference_samples_sha256':hashlib.sha256('\n'.join(ref).encode()).hexdigest()})
  for _,rr in g.iterrows():
   hid=rr.current_HGNC_IDs;tech=hid in gm;keep=hid in out['A'];reason='' if keep else ('not_one_to_one' if not rr.human_mapping_unambiguous else ('no_single_identity_probe' if not tech else 'nonfinite_or_SD_le_1e-8'))
   gene_selection.append({'gene_set':name,'scale':input_scale,'source_row_id':rr.source_row_id,'mouse_symbol':rr.mouse_symbol,'human_symbol':rr.current_human_symbols,'HGNC_ID':hid,'Entrez':rr.current_human_Entrez,'orthology':rr.ortholog_type,'G0_one_to_one':rr.human_mapping_unambiguous,'technical_available':tech,'included_in_A':keep,'exclusion_reason':reason,'qualified_probe_ids':';'.join(gm.get(hid,[])),'actual_expression_status':'not_assessable_from_ratio_matrix'})
  for hid in out['A']:
   rr=genemeta[name].loc[hid];gene_refs.append({'gene_set':name,'scale':input_scale,'HGNC_ID':hid,'human_symbol':rr.current_human_symbols,'reference_n':43,'mu':out['mu'][hid],'sample_SD':out['sd'][hid],'ddof':1,'fixed_weight':1.0,'fixed_mean_denominator':len(out['A'])})
   if input_scale=='log2_ratio':
    for p in gm[hid]:probe_rows.append({'gene_set':name,'HGNC_ID':hid,'human_symbol':rr.current_human_symbols,'probe_id':p,'mapping':'preflight_single_identity','applied_log2_count_main':1})
  if out['status']=='ready':
   score_refs.append({'gene_set':name,'scale':input_scale,'reference_n':43,'S_mean_reference':out['S_mean_reference'],'S_sample_SD_reference':out['S_sd_reference'],'gene_n':len(out['A']),'score_SD_ddof':1,'reference_samples':';'.join(ref)})
   assert abs(out['Y'].loc[ref].mean())<1e-12 and np.isclose(out['Y'].loc[ref].std(ddof=1),1,atol=1e-12)
   for sample in samples:
    score_records.append({'sample':sample,'donor_id':hip.loc[sample,'donor_id'],'gene_set':name,'scale':input_scale,'in_normal43_reference':sample in ref,'S_mean_gene_z':out['S'][sample],'Y_score_reference_SD':out['Y'][sample],'A_n':len(out['A']),'fixed_mean_denominator':len(out['A']),'source_log2_count':0,'applied_log2_count':logcount})
  gv.reset_index(names='HGNC_ID').to_csv(RES/f'gene_aggregates_{name}_{input_scale}.tsv.gz',sep='\t',index=False)
assert source_ratio.equals(original_ratios)
for frame,fn in [(pd.DataFrame(gene_selection),'gene_selection.tsv'),(pd.DataFrame(probe_rows),'qualified_probe_lists.tsv'),(pd.DataFrame(gene_refs),'gene_reference_parameters.tsv'),(pd.DataFrame(score_refs),'score_reference_parameters.tsv'),(pd.DataFrame(coverage),'set_coverage_and_gates.tsv'),(pd.DataFrame(score_records),'donor_scores.tsv')]:tab(frame,fn)
# Technical score manifests above are saved before any association models.
dump({'technical_scoring_completed_at_utc':now(),'human_models_started':False,'reference_n':43,'genes_and_scoring_fixed':True,'scoring_outputs_sha256':{n:sha(RES/n) for n in ['gene_selection.tsv','qualified_probe_lists.tsv','gene_reference_parameters.tsv','score_reference_parameters.tsv','set_coverage_and_gates.tsv','donor_scores.tsv']}},LOG/'pre_model_scoring_lock.json')
results=[];flows=[];diagrows=[];loorows=[];manual_checks=[];coefficients=[]
for tid in ['P01','E01','D01','S07','S08']:
 rr=reg.loc[tid].to_dict();name=rr['gene_set'];input_scale='raw_ratio' if tid in ['S07','S08'] else 'log2_ratio';out=scores[(name,input_scale)];cr=next(c for c in coverage if c['gene_set']==name and c['scale']==input_scale)
 population=hip.D01_population if tid=='D01' else hip.primary_age_manifest;d=hip.loc[population].copy();d['Y_score']=out.get('Y',pd.Series(dtype=float)).reindex(d.index);terms=['age10','sex','AD'] if tid=='D01' else ['age10','sex'];target=rr['term_tested'];model,X,gate=design_gate(d,terms)
 row={'test_id':tid,**rr,'input_scale':input_scale,'n_population':len(d),'n':len(model),'n_normal':int(model.diagnosis.eq('normal').sum()),'n_AD':int(model.diagnosis.eq('AD').sum()),'reference_n':43,'G0_n':cr['G0_one_to_one_n'],'technical_n':cr['technical_available_n'],'actual_expression_status':cr['actual_expression_status'],'A_n':cr['A_n'],'number_coverage':cr['number_coverage'],'weight_coverage':cr['weight_coverage'],'minimum_genes':cr['minimum_genes'],'coverage_gate_pass':cr['coverage_gate_pass'],'score_status':out['status'],**gate,'adjusted_P':np.nan,'adjustment_status':'pending-family','CI_type':'pointwise_95pct_t_HC3','beta_unit':'normal43_score_SD_per_AD_vs_normal' if tid=='D01' else 'normal43_score_SD_per_10_years','model_score_reference':'full_normal43_no_model_specific_rescaling'}
 for sample,mr in hip.iterrows():
  inside=sample in d.index;included=sample in model.index;flows.append({'test_id':tid,'sample':sample,'donor_id':mr.donor_id,'diagnosis':mr.diagnosis,'age':mr.age,'sex_code':mr.sex,'population_eligible':inside,'complete_case':included,'exclusion_reason':'' if included else ('outside_registered_population' if not inside else 'technical_NA_or_covariate_missing'),'Y_score':out.get('Y',pd.Series(dtype=float)).get(sample,np.nan)})
 if not cr['coverage_gate_pass']:row.update(status='not_testable_coverage',reason='below_fixed_gene_coverage')
 elif out['status']!='ready':row.update(status=out['status'],reason='fixed_A_or_score_variance')
 elif gate['status']=='ready':
  fitted,diag,loo=fit_hc3(X,model.Y_score,target);row.update(fitted)
  if fitted['status']=='tested':
   # Independent normal-equation HC3 spot check for EVERY target coefficient.
   xx=X.to_numpy();yy=model.Y_score.to_numpy();bread=np.linalg.inv(xx.T@xx);b=bread@xx.T@yy;e=yy-xx@b;h=np.einsum('ij,jk,ik->i',xx,bread,xx);cov=bread@(xx.T@np.diag(e**2/(1-h)**2)@xx)@bread;j=list(X.columns).index(target);manual_se=math.sqrt(cov[j,j]);df=len(yy)-np.linalg.matrix_rank(xx);manual_p=float(2*stats.t.sf(abs(b[j]/manual_se),df));manual_ci=b[j]+np.array([-1,1])*stats.t.ppf(.975,df)*manual_se
   assert np.allclose([fitted['beta'],fitted['SE_HC3'],fitted['raw_P']],[b[j],manual_se,manual_p],rtol=1e-10,atol=1e-12);assert np.allclose([fitted['CI95_low'],fitted['CI95_high']],manual_ci,rtol=1e-10,atol=1e-12)
   manual_checks.append({'test_id':tid,'manual_beta':float(b[j]),'manual_SE_HC3':manual_se,'manual_df':int(df),'manual_raw_P':manual_p,'manual_CI95_low':manual_ci[0],'manual_CI95_high':manual_ci[1],'match':True})
   # All model coefficients saved, but only registered target term is tested/reported.
   for k,t in enumerate(X.columns):coefficients.append({'test_id':tid,'term':t,'coefficient':b[k],'registered_tested_term':t==target})
   diag.index.name='sample';diag=diag.reset_index();diag['test_id']=tid;diag['donor_id']=diag['sample'].map(hip.donor_id);diagrows.append(diag);loo['test_id']=tid;loo['omitted_donor']=loo.omitted_sample.map(hip.donor_id);loorows.append(loo)
   tab(pd.concat([model[['donor_id','Y_score','age','diagnosis']],X],axis=1).reset_index(),'model_design_'+tid+'.tsv')
   fig,ax=plt.subplots(1,3,figsize=(12,3.3));ax[0].scatter(diag.fitted,diag.residual,s=22,alpha=.8);ax[0].axhline(0,color='black',lw=.8);ax[0].set(xlabel='Fitted score (reference SD)',ylabel='OLS residual',title=tid+' residuals')
   ax[1].scatter(np.arange(len(diag)),diag.hat,s=18,label='Hat');ax[1].axhline(2*X.shape[1]/len(X),ls='--',color='red',lw=.8);ax[1].set(xlabel='Donor order',ylabel='Leverage',title='Fixed-design leverage')
   ax[2].scatter(np.arange(len(loo)),loo.beta,s=18);ax[2].axhline(fitted['beta'],ls='--',color='black',lw=.8);ax[2].axhline(0,color='grey',lw=.6);ax[2].set(xlabel='Omitted donor order',ylabel='Registered-term beta',title='Leave-one-donor-out (no P)');fig.tight_layout();fig.savefig(FIG/f'{tid}_diagnostics.png',dpi=170);fig.savefig(FIG/f'{tid}_diagnostics.pdf');plt.close(fig)
 if tid=='D01':
  row['adjustment_status']='complete_E_AD_BH_m1';row['adjusted_P']=row.get('raw_P',np.nan)
 results.append(row)
# Pending registered GSE intersection: no A invented, no P=1 substituted in public output.
rr=reg.loc['S16'].to_dict();results.append({'test_id':'S16',**rr,'status':'pending_SEA_actual_A_intersection','reason':'Await verified main actual-expression/variance gene-A for SEA CA3 and SUB; do not precompute intersection','adjusted_P':np.nan,'adjustment_status':'pending-family','raw_P':np.nan})
result=pd.DataFrame(results);tab(result,'registered_results_long.tsv');tab(pd.DataFrame(flows),'model_sample_manifest.tsv');tab(pd.concat(diagrows,ignore_index=True),'model_diagnostics.tsv');tab(pd.concat(loorows,ignore_index=True),'leave_one_donor_out_beta_only.tsv');tab(pd.DataFrame(coefficients),'model_coefficients_descriptive.tsv');tab(pd.DataFrame(manual_checks),'manual_HC3_spotchecks.tsv')
# Manual score spotcheck using the saved source values, no second normalizer.
name='core23_1to1';sample=ref[0];hid=next(iter(maps[name]));ps=maps[name][hid];observed=source_ratio.loc[ps,sample].to_numpy();manual_gene=float(np.median([math.log2(float(v)) for v in observed]));computed=float(stage_values[(name,'log2_ratio')].loc[hid,sample]);assert np.isclose(manual_gene,computed,atol=1e-12)
o=scores[(name,'log2_ratio')];manual_z=[(float(stage_values[(name,'log2_ratio')].loc[h,sample])-float(o['mu'][h]))/float(o['sd'][h]) for h in o['A']];manual_S=sum(manual_z)/len(manual_z);manual_Y=(manual_S-o['S_mean_reference'])/o['S_sd_reference'];assert np.isclose(manual_Y,o['Y'][sample],atol=1e-12)
dump({'sample':sample,'gene':hid,'probe_ids':ps,'original_ratios':observed.tolist(),'log2_values':[math.log2(float(v)) for v in observed],'median_log2_manual':manual_gene,'median_log2_pipeline':computed,'all_gene_z_manual':manual_z,'S_manual':manual_S,'S_pipeline':float(o['S'][sample]),'Y_manual':manual_Y,'Y_pipeline':float(o['Y'][sample]),'PASS':True},LOG/'manual_score_spotcheck.json')
# Forest display includes raw sensitivity as separate registered tests, never picks a winner.
r=result[result.status.eq('tested')].copy();fig,ax=plt.subplots(figsize=(8,4));order=list(r.test_id)
for k,rr in enumerate(r.itertuples()):ax.errorbar(rr.beta,k,xerr=[[rr.beta-rr.CI95_low],[rr.CI95_high-rr.beta]],fmt='o',capsize=3,color='#1f4e79' if rr.input_scale=='log2_ratio' else '#8c564b')
ax.axvline(0,color='grey',lw=1);ax.set_yticks(range(len(order)),[f'{t}: '+('AD extension' if t=='D01' else ('raw-ratio sensitivity' if t.startswith('S') else 'log2 primary scale')) for t in order]);ax.invert_yaxis();ax.set_xlabel('Normal43 score SD / 10 years (D01: AD vs normal)');ax.set_title('Registered GSE effects: pointwise 95% HC3 t intervals\nFamily adjustments pending except E_AD (m=1)');fig.tight_layout();fig.savefig(FIG/'registered_effects.png',dpi=170);fig.savefig(FIG/'registered_effects.pdf');plt.close(fig)
after={str(p):sha(p) for p in source_paths};assert source_hashes==after
# Write completion marker only after all outputs verified.
assert len(result)==6 and set(result.loc[result.status.eq('tested'),'test_id'])=={'P01','E01','D01','S07','S08'}
assert result.set_index('test_id').loc['D01','n']==40
assert result.set_index('test_id').loc['D01','n_normal']==21 and result.set_index('test_id').loc['D01','n_AD']==19
assert result.set_index('test_id').loc[['P01','E01','S07','S08'],'n'].eq(43).all()
assert result.loc[result.family.ne('E_AD'),'adjusted_P'].isna().all()
checks={'status':'PASS','completed_at_utc':now(),'human_models_first_started_after':json.loads((LOG/'pre_model_scoring_lock.json').read_text())['technical_scoring_completed_at_utc'],'source_files_unchanged':True,'source_file_n':len(source_paths),'freeze_SHA_and_five_candidate_MD5':'PASS','source_ratio_unchanged':True,'double_log_assertion_test':'PASS','synthetic_validation':'PASS','manual_score_spotcheck':'PASS','manual_HC3_all_five_models':'PASS','normal_n':43,'normal_female':20,'normal_male':23,'unique_hip_donors':62,'AD_overlap_total_n':40,'AD_overlap_normal_n':21,'AD_overlap_AD_n':19,'models_executed':['P01','E01','D01','S07','S08'],'pending':['S16'],'family_adjustment_pending':['PRIMARY_m3','E_EXPR_m11','E_SENS_m23'],'E_AD_m1_complete':True,'software_versions':versions,'result_SHA256':sha(RES/'registered_results_long.tsv')}
dump(checks,LOG/'run_complete.json');print(result[['test_id','family','status','n','A_n','number_coverage','beta','SE_HC3','df','CI95_low','CI95_high','raw_P','adjusted_P','adjustment_status','LOO_beta_min','LOO_beta_max','LOO_direction_flip_n','unstable_flag']].to_string(index=False));print('ALL INPUT HASHES UNCHANGED; SYNTHETIC, SCORE AND HC3 MANUAL VALIDATION PASS')
