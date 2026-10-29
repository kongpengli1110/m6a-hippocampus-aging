"""Execute only frozen SEA expression positions and newly released GSE S16.
No cognition, H2/H4, single-gene associations, model search or old GSE refits.
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
import json,hashlib,sys,platform,resource
import numpy as np,pandas as pd,scipy,statsmodels
from scipy import stats
from sea_helpers import *
ROOT=Path(__file__).resolve().parents[1];AN=ROOT.parent;PROJ=AN.parent;PRE=PROJ/'human_preflight';GSE=AN/'gse48350_v1';D=ROOT/'derived';R=ROOT/'results';L=ROOT/'logs'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,p):p.write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))
def tab(df,n):df.to_csv(R/n,sep='\t',index=False)
def now():return datetime.now(timezone.utc).isoformat()
if (L/'execution_PASS.json').exists():raise SystemExit('Completed SEA run exists; do not repeat successful analyses.')
prep=json.loads((L/'preparation_PASS.json').read_text());assert prep['status']=='PASS';assert json.loads((L/'synthetic_PASS.json').read_text())['status']=='PASS'
inputs=json.loads((ROOT/'input_manifest.json').read_text())['files'];assert all(sha(p)==h for p,h in inputs.items())
for n,h in prep['derived_files'].items():assert sha(D/n)==h
registry=pd.read_csv(AN/'prespecified_tests_v1_frozen.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('test_id');authorized=registry[(registry.cohort=='SEAAD_HIP')&registry.family.isin(['PRIMARY','E_EXPR','E_MODULE','E_SENS'])&~registry.index.isin(['S19','S20','S21','S22'])].copy();assert len(authorized)==35
x=pd.read_csv(PRE/'human_candidate_crosswalk.tsv',sep='\t',dtype=str,keep_default_na=False);x['human_mapping_unambiguous']=strict_bool(x.human_mapping_unambiguous);x['SEAAD_primary_feature_available']=strict_bool(x.SEAAD_primary_feature_available);x.loc[x.set_kind.eq('core'),'set_name']='core23_1to1'
clinical=pd.read_csv(D/'clinical_covariates_verified.tsv',sep='\t').set_index('donor');tm=pd.read_csv(D/'target_donor_manifest.tsv',sep='\t');
for c in ['eligible_ge10','eligible_ge20','eligible_ge30','has_observed_cells']:tm[c]=strict_bool(tm[c])
audit=pd.read_csv(D/'candidate_feature_identity_audit.tsv',sep='\t',dtype=str,keep_default_na=False);gene_features={hid:g.gene_index.astype(int).to_numpy() for hid,g in audit.groupby('HGNC_IDs',sort=True)}
params={'frozen_plan_SHA256':sha(AN/'human_analysis_plan_v1.md'),'registry_SHA256':sha(AN/'prespecified_tests_v1_frozen.tsv'),'started_at_utc':now(),'source_layer':'UMIs','all_feature_library_denominator':36601,'expression_rule':'CPM>=1 in >=max(5,ceil(0.2*nR)) reference donors','transform':'log2(CPM+1) once from raw counts','primary_nuclei_min':20,'sensitivity_nuclei_min':[10,30],'reference':'all technical eligible >=20 per target; no phenotype stratification','minimum_gene_rule':'max(10,ceil(.8*G0))','module_weight_coverage_min':.8,'gene_and_score_sd_ddof':1,'sd_epsilon':1e-8,'age10':'(age-80)/10','sex':'Female0 Male1','PMI_RIN_center':'same primary target reference','Braak':'strip explicit Braak prefix; 0,I,II,III,IV,V,VI to 0..6','OLS':'HC3 t(df=n-rank); two-sided pointwise95CI','feature_duplicates':'official-HGNC-identity-audited preflight features summed once; no extra weight','blocked_tests':['C01-C08','S19-S22','H2_01-H2_08','H4_01'],'completed_GSE_models_refit':False}
dump(params,ROOT/'run_parameters.json');dump({'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'statsmodels':statsmodels.__version__,'scripts_SHA256':{p.name:sha(p) for p in (ROOT/'scripts').glob('*.py')},'reused_OLS_helper':str(GSE/'scripts/fixed_pipeline.py'),'reused_OLS_helper_SHA256':sha(GSE/'scripts/fixed_pipeline.py')},ROOT/'software_versions.json')
data={};gates=[];covrefs=[]
for target in tm.target.unique():
 loaded=np.load(D/(target.replace('+','_')+'_all_feature_counts.npz'));raw=loaded['counts'];donors=loaded['donors'].tolist();m=tm[tm.target.eq(target)].set_index('donor').loc[donors].copy();assert donors==clinical.index.tolist();fullL=raw.sum(axis=1);assert np.array_equal(fullL,m.totalUMI_all36601)
 gc=np.stack([raw[:,ii].sum(axis=1) for ii in gene_features.values()],axis=1);assert gc.dtype==np.int64;cp=full_background_cpm(gc,fullL);vals=np.log2(cp+1);ids=list(gene_features);df=pd.DataFrame(vals.T,index=ids,columns=donors);cpm=pd.DataFrame(cp.T,index=ids,columns=donors);ref=m.index[m.eligible_ge20].tolist();minimum=max(5,int(np.ceil(.2*len(ref))));nexpressed=(cpm.loc[:,ref]>=1).sum(axis=1);sds=df.loc[:,ref].std(axis=1,ddof=1,skipna=False);finite=np.isfinite(df.loc[:,ref]).all(axis=1);eg=nexpressed>=minimum;vg=finite & sds.gt(1e-8);eligible=eg&vg
 gate=pd.DataFrame({'target':target,'HGNC_ID':ids,'reference_n':len(ref),'expression_min_donors':minimum,'CPM_ge1_reference_donors':nexpressed,'actual_expression_pass':eg,'reference_finite':finite,'reference_sample_SD':sds,'variance_pass':vg,'final_global_gene_gate':eligible}).set_index('HGNC_ID',drop=False);gates.append(gate.reset_index(drop=True))
 meta=clinical.loc[donors].join(m.drop(columns='target'));meta['PMI_c']=meta.PMI-meta.loc[ref,'PMI'].mean();meta['RIN_c']=meta.RIN-meta.loc[ref,'RIN'].mean();meta['target']=target;data[target]={'x':df,'cpm':cpm,'counts':pd.DataFrame(gc.T,index=ids,columns=donors),'gate':gate,'R':ref,'meta':meta,'L':pd.Series(fullL,index=donors)};covrefs.append({'target':target,'reference_n':len(ref),'reference_donors':';'.join(ref),'PMI_mean_R':meta.loc[ref,'PMI'].mean(),'RIN_mean_R':meta.loc[ref,'RIN'].mean(),'expression_min_donors':minimum,'all_gene_feature_n':36601})
 df.to_csv(D/(target.replace('+','_')+'_candidate_log2CPM1.tsv.gz'),sep='\t',index_label='HGNC_ID');data[target]['counts'].to_csv(D/(target.replace('+','_')+'_candidate_counts.tsv.gz'),sep='\t',index_label='HGNC_ID');del raw,gc,cp,vals,loaded
# Per-score tables include each original source member and stage-specific exclusion.
score_map={};score_defs={};coverage_rows=[];selection_rows=[];gene_refrows=[];score_refrows=[];donor_scores=[]
def create_score(key,target,name,variant='main',A_override=None,basekey=None):
 d=data[target];src=x[x.set_name.eq(name)].copy()
 if variant=='kME_ge0.6':src=src[pd.to_numeric(src.kME,errors='coerce')>=.6]
 g0=src[src.human_mapping_unambiguous];assert g0.current_HGNC_IDs.is_unique;G0=g0.current_HGNC_IDs.tolist();module=src.set_kind.eq('modules').all();weights=pd.Series(pd.to_numeric(g0.kME).to_numpy() if module else np.ones(len(g0)),index=G0,dtype=float)
 tech=[g for g in G0 if g in d['x'].index];eg=[g for g in tech if d['gate'].loc[g,'actual_expression_pass']];vg=[g for g in eg if d['gate'].loc[g,'variance_pass']];A=[g for g in G0 if g in set(vg)] if A_override is None else [g for g in G0 if g in set(A_override)];assert set(A)<=set(vg)
 cr={'score_id':key,'target':target,'gene_set':name,'variant':variant,'source_member_n':len(src),'one_to_one_n':len(G0),'technical_n':len(tech),'actual_expression_pass_n':len(eg),'variance_pass_after_expression_n':len(vg),'reference_n':len(d['R']),'expression_min_donors':max(5,int(np.ceil(.2*len(d['R'])))),**coverage(A,G0,weights)}
 out=weighted_fixed_score(d['x'],d['R'],A,weights);cr['score_status']=out['status'];cr['A_SHA256']=hashlib.sha256('\n'.join(A).encode()).hexdigest();cr['base_score_id']=basekey or '';out.update(coverage=cr,target=target,gene_set=name,variant=variant);score_map[key]=out;score_defs[key]={'G0':G0,'A':A,'weights':weights};coverage_rows.append(cr)
 for _,r in src.iterrows():
  hid=r.current_HGNC_IDs;reason='' if hid in A else ('not_one_to_one' if not r.human_mapping_unambiguous else ('no_unambiguous_platform_feature' if hid not in tech else ('below_actual_expression_gate' if hid not in eg else ('nonfinite_or_low_reference_SD' if hid not in vg else 'outside_prespecified_subset'))))
  selection_rows.append({'score_id':key,'target':target,'gene_set':name,'variant':variant,'source_row_id':r.source_row_id,'mouse_symbol':r.mouse_symbol,'human_symbol':r.current_human_symbols,'HGNC_ID':hid,'Entrez':r.current_human_Entrez,'source_feature_ids':r.SEAAD_unique_feature_ids,'G0_one_to_one':r.human_mapping_unambiguous,'technical_pass':hid in tech,'actual_expression_pass':hid in eg,'variance_pass':hid in vg,'included_A':hid in A,'weight':weights.get(hid,np.nan),'exclusion_reason':reason})
 if out['status']=='ready':
  assert np.isclose(out['Y'][d['R']].mean(),0,atol=1e-12) and np.isclose(out['Y'][d['R']].std(ddof=1),1,atol=1e-12)
  for hid in A:gene_refrows.append({'score_id':key,'target':target,'gene_set':name,'HGNC_ID':hid,'mu':out['mu'][hid],'sample_SD':out['sd'][hid],'ddof':1,'weight':out['weights'][hid],'fixed_weight_denominator':out['weight_denominator']})
  score_refrows.append({'score_id':key,'target':target,'gene_set':name,'variant':variant,'reference_n':len(d['R']),'reference_donors':';'.join(d['R']),'S_mean':out['S_mean'],'S_sample_SD':out['S_SD'],'weight_denominator':out['weight_denominator'],'A_n':len(A)})
  for donor in d['meta'].index:donor_scores.append({'score_id':key,'target':target,'gene_set':name,'variant':variant,'donor':donor,'in_primary20_reference':donor in d['R'],'S':out['S'][donor],'Y':out['Y'][donor],'A_n':len(A),'fixed_weight_denominator':out['weight_denominator'],'source':'raw_UMIs','transform':'log2(CPM+1)'})
 return out
base_ids=authorized[~authorized.index.str.startswith('S')]
for tid,r in base_ids.iterrows():
 key=r.target+'::'+r.gene_set
 if key not in score_map:create_score(key,r.target,r.gene_set)
assert len(score_map)==19
# High-kME module subsets use existing main A, never rerun expression screening.
for tid in [f'S{i:02d}' for i in range(9,16)]:
 r=registry.loc[tid];base=r.target+'::'+r.gene_set;g0=x[(x.set_name==r.gene_set)&x.human_mapping_unambiguous];hi=set(g0.loc[pd.to_numeric(g0.kME)>=.6,'current_HGNC_IDs']);A=set(score_map[base]['A'])&hi;create_score(tid,r.target,r.gene_set,'kME_ge0.6',A_override=A,basekey=base)
# Formal cross-platform intersection uses actual SEA main A, not preflight var coverage.
gsel=pd.read_csv(GSE/'results/gene_selection.tsv',sep='\t',dtype=str,keep_default_na=False);gsel['included_in_A']=strict_bool(gsel.included_in_A);intersections={}
for name in ['core23_1to1','METTL3cKO_HIP_down']:
 ga=set(gsel.loc[gsel.gene_set.eq(name)&gsel.scale.eq('log2_ratio')&gsel.included_in_A,'HGNC_ID']);ca=set(score_map['CA3_1::'+name]['A']);su=set(score_map['Sub_1+Sub_2::'+name]['A']);intersections[name]=sorted(ga&ca&su)
 dump({'gene_set':name,'GSE_main_A':sorted(ga),'CA3_actual_main_A':sorted(ca),'SUB_actual_main_A':sorted(su),'intersection':intersections[name],'no_new_expression_gate':True},D/(name+'_crossplatform_intersection.json'))
for tid in ['S17','S18']:
 r=registry.loc[tid];create_score(tid,r.target,r.gene_set,'crossplatform_intersection',A_override=intersections[r.gene_set],basekey=r.target+'::'+r.gene_set)
# S16: cached already-log2 GSE per-gene data, no matrix re-read/logging/base refit.
gmeta=pd.read_csv(GSE/'results/sample_manifest.tsv',sep='\t').set_index('sample');gmeta['primary_age_manifest']=strict_bool(gmeta.primary_age_manifest);gref=gmeta.index[gmeta.primary_age_manifest].tolist();gv=pd.read_csv(GSE/'results/gene_aggregates_core23_1to1_log2_ratio.tsv.gz',sep='\t',index_col=0);Aint=intersections['core23_1to1'];gg0=x[x.set_name.eq('core23_1to1')&x.human_mapping_unambiguous].current_HGNC_IDs.tolist();gw=pd.Series(1.,index=gg0);s16=weighted_fixed_score(gv,gref,Aint,gw);s16.update(target='bulk_HIP',gene_set='core23_1to1',variant='crossplatform_intersection',coverage={'score_id':'S16','target':'bulk_HIP','gene_set':'core23_1to1','variant':'crossplatform_intersection','reference_n':43,'technical_n':21,'actual_expression_pass_n':np.nan,'actual_expression_status':'GSE_not_assessable;SEA_actual_A_used_for_intersection',**coverage(Aint,gg0,gw)});score_map['S16']=s16;coverage_rows.append(s16['coverage']);data['bulk_HIP']={'meta':gmeta,'R':gref,'x':gv}
if s16['status']=='ready':
 score_refrows.append({'score_id':'S16','target':'bulk_HIP','gene_set':'core23_1to1','variant':'crossplatform_intersection','reference_n':43,'reference_donors':';'.join(gref),'S_mean':s16['S_mean'],'S_sample_SD':s16['S_SD'],'weight_denominator':s16['weight_denominator'],'A_n':len(Aint)})
 for hid in Aint:gene_refrows.append({'score_id':'S16','target':'bulk_HIP','gene_set':'core23_1to1','HGNC_ID':hid,'mu':s16['mu'][hid],'sample_SD':s16['sd'][hid],'ddof':1,'weight':1,'fixed_weight_denominator':len(Aint)})
 for donor in gv.columns:donor_scores.append({'score_id':'S16','target':'bulk_HIP','gene_set':'core23_1to1','variant':'crossplatform_intersection','donor':donor,'in_primary20_reference':donor in gref,'S':s16['S'][donor],'Y':s16['Y'][donor],'A_n':len(Aint),'fixed_weight_denominator':len(Aint),'source':'cached_GSE_log2_gene_aggregates','transform':'NO_additional_log'})
for df,name in [(pd.concat(gates),'actual_expression_and_variance_gates.tsv'),(pd.DataFrame(covrefs),'covariate_reference_parameters.tsv'),(pd.DataFrame(coverage_rows),'score_coverage.tsv'),(pd.DataFrame(selection_rows),'gene_selection.tsv'),(pd.DataFrame(gene_refrows),'gene_reference_parameters.tsv'),(pd.DataFrame(score_refrows),'score_reference_parameters.tsv'),(pd.DataFrame(donor_scores),'donor_scores.tsv')]:tab(df,name)
dump({'status':'PASS','completed_at_utc':now(),'human_model_fitting_started':False,'source_expression_transform_count':1,'scoring_blind_to_age_diagnosis_outcome':True,'outputs':{n:sha(R/n) for n in ['actual_expression_and_variance_gates.tsv','score_coverage.tsv','gene_selection.tsv','gene_reference_parameters.tsv','score_reference_parameters.tsv','donor_scores.tsv']}},L/'pre_model_score_lock.json')
# Registered models only. NaNs retain test slots and explicit reasons.
res=[];res_by={};flow=[];diags=[];loos=[];manual=[];design_ranges=[];allids=list(authorized.index)+['S16'];gseres=pd.read_csv(GSE/'results/registered_results_long.tsv',sep='\t').set_index('test_id')
for tid in allids:
 r=registry.loc[tid];key=tid if tid in [f'S{i:02d}' for i in range(9,19)] else r.target+'::'+r.gene_set;o=score_map[key];cr=o['coverage'];m=data[r.target]['meta'].copy();threshold=int(r.nuclei_min) if r.cohort=='SEAAD_HIP' else None;population=m['eligible_ge'+str(threshold)] if threshold else m.primary_age_manifest;m['Y_score']=o.get('Y',pd.Series(dtype=float)).reindex(m.index);terms=['age10','sex'] if r.cohort=='GSE48350' else ['age10','sex','PMI_c','RIN_c']
 if tid in ['S05','S06']:terms+=['Braak_numeric']
 if tid=='S23':terms+=['Sub_2_nuclei_fraction']
 md,X,gate=design_gate(m.loc[population],terms);row={'test_id':tid,**r.to_dict(),'score_id':key,'input_scale':'cached_log2_ratio' if tid=='S16' else 'log2_CPM_plus1','n_population':int(population.sum()),'n':len(md),**cr,**gate,'beta':np.nan,'SE_HC3':np.nan,'CI95_low':np.nan,'CI95_high':np.nan,'raw_P':np.nan,'adjusted_P':np.nan,'adjustment_status':'pending_combined_family' if r.family!='E_SENS' else 'pending-family','beta_unit':'fixed_reference_score_SD_per_10_years','CI_type':'pointwise_HC3_t_95pct','model_execution':'not_fitted'}
 for donor,mr in m.iterrows():flow.append({'test_id':tid,'score_id':key,'target':r.target,'donor':donor,'n_nuclei':mr.get('n_nuclei',np.nan),'totalUMI_all36601':mr.get('totalUMI_all36601',np.nan),'registered_population_eligible':bool(population.loc[donor]),'complete_case':donor in md.index,'Y_score':mr.Y_score,'age':mr.age,'sex':mr.sex,'PMI_c':mr.get('PMI_c',np.nan),'RIN_c':mr.get('RIN_c',np.nan),'Braak_numeric':mr.get('Braak_numeric',np.nan),'exclusion_reason':'' if donor in md.index else ('outside_registered_population' if not population.loc[donor] else 'score_or_covariate_missing')})
 for term in terms:design_ranges.append({'test_id':tid,'covariate':term,'complete_case_min':md[term].min(),'complete_case_max':md[term].max(),'missing_in_population_n':int(m.loc[population,term].isna().sum())})
 if not cr['coverage_gate_pass']:row.update(status='not_testable_coverage',reason='below_fixed_count_or_weight_coverage')
 elif o['status']!='ready':row.update(status='not_testable_score',reason='nonfinite_or_low_score_SD')
 elif tid in [f'S{i:02d}' for i in range(9,16)] and o['A']==score_map[r.target+'::'+r.gene_set]['A']:
  parent=res_by[r.parent_test];row.update(status='identical_prespecified_subset',reason='same_actual_A_and_weights_as_parent;no_new_P',model_execution='not_refit_identical');
  for col in ['beta','SE_HC3','df','CI95_low','CI95_high']:row[col]=parent.get(col,np.nan)
  row['effect_source_test_id']=r.parent_test;row['raw_P']=np.nan
 elif gate['status']=='ready':
  reused=False
  if tid=='S16':
   prior=pd.read_csv(GSE/'results/donor_scores.tsv',sep='\t');prior=prior[(prior.gene_set=='core23_1to1')&(prior.scale=='log2_ratio')].set_index('sample');gA=set(gsel.loc[gsel.gene_set.eq('core23_1to1')&gsel.scale.eq('log2_ratio')&gsel.included_in_A,'HGNC_ID'])
   if set(Aint)==gA:
    assert np.allclose(o['Y'].loc[gref],prior.loc[gref,'Y_score_reference_SD'],rtol=1e-12,atol=1e-12);p=gseres.loc['P01'];keys=['beta','SE_HC3','df','t_HC3','CI95_low','CI95_high','raw_P','max_hat','max_Cook','age_VIF','hat_threshold','Cook_threshold','n_hat_flag','n_Cook_flag','age_VIF_flag','unstable_flag','LOO_beta_min','LOO_beta_max','LOO_direction_flip_n','LOO_not_estimable_n'];row.update({k:p[k] for k in keys if k in p});row.update(status='tested_reused_identical_GSE_P01',reason='intersection equals GSE main A; scores exactly match; reuse completed P01 without refit',model_execution='reused_completed_GSE_P01');reused=True
  if not reused:
   fit,diag,loo=fit_hc3(X,md.Y_score,'age10');row.update(fit);row['model_execution']='fitted_registered_position_once'
   if fit['status']=='tested':
    xx=X.to_numpy();yy=md.Y_score.to_numpy();B=np.linalg.inv(xx.T@xx);b=B@xx.T@yy;e=yy-xx@b;h=np.einsum('ij,jk,ik->i',xx,B,xx);V=B@(xx.T@np.diag(e**2/(1-h)**2)@xx)@B;j=list(X.columns).index('age10');se=float(np.sqrt(V[j,j]));pv=float(2*stats.t.sf(abs(b[j]/se),len(yy)-np.linalg.matrix_rank(xx)));assert np.allclose([fit['beta'],fit['SE_HC3'],fit['raw_P']],[b[j],se,pv],rtol=1e-9,atol=1e-11);manual.append({'test_id':tid,'manual_beta':b[j],'manual_SE_HC3':se,'manual_df':len(yy)-X.shape[1],'manual_raw_P':pv,'match':True})
    diag.index.name='donor';diag=diag.reset_index();diag['test_id']=tid;diag['target']=r.target;diags.append(diag);loo['test_id']=tid;loos.append(loo);tab(pd.concat([md[['Y_score','age']],X],axis=1).reset_index(),'design_'+tid+'.tsv')
    # Braak-complete-case comparison is needed only if its covariate caused loss.
    if tid in ['S05','S06'] and len(md)<len(m.loc[population]):
     baseX=X.drop(columns='Braak_numeric');bb=np.linalg.lstsq(baseX.to_numpy(),md.Y_score.to_numpy(),rcond=None)[0];row['paired_complete_case_main_beta_descriptive']=bb[list(baseX.columns).index('age10')]
 res.append(row);res_by[tid]=row
out=pd.DataFrame(res);assert len(out)==36 and out.test_id.is_unique;tab(out,'registered_results_long.tsv');tab(pd.DataFrame(flow),'model_sample_manifest.tsv');tab(pd.DataFrame(design_ranges),'design_covariate_ranges.tsv');tab(pd.concat(diags,ignore_index=True) if diags else pd.DataFrame(),'model_diagnostics.tsv');tab(pd.concat(loos,ignore_index=True) if loos else pd.DataFrame(),'leave_one_donor_out_beta_only.tsv');tab(pd.DataFrame(manual),'manual_HC3_spotchecks.tsv')
# Scalar primary score check; independent arithmetic uses full-library UMI denominator.
target='CA3_1';name='METTL3cKO_HIP_down';o=score_map[target+'::'+name];donor=data[target]['R'][0]
if o['status']=='ready':
 hid=o['A'][0];count=int(data[target]['counts'].loc[hid,donor]);library=int(data[target]['L'][donor]);cpm=count*1e6/library;xx=float(np.log2(cpm+1));assert np.isclose(xx,data[target]['x'].loc[hid,donor]);zs=[(float(data[target]['x'].loc[g,donor])-o['mu'][g])/o['sd'][g] for g in o['A']];s=sum(zs)/len(zs);y=(s-o['S_mean'])/o['S_SD'];assert np.isclose(y,o['Y'][donor],atol=1e-12);dump({'PASS':True,'target':target,'donor':donor,'HGNC_ID':hid,'raw_gene_count':count,'full_all36601_library':library,'manual_CPM':cpm,'manual_log2_CPM1':xx,'manual_S':s,'pipeline_S':o['S'][donor],'manual_Y':y,'pipeline_Y':o['Y'][donor],'A_n':len(o['A'])},L/'scalar_primary_score_spotcheck.json')
assert all(sha(p)==h for p,h in inputs.items());passed={'status':'PASS','completed_at_utc':now(),'source_files_unchanged':True,'GSE_five_models_not_refit':True,'registered_positions':36,'status_counts':out.status.value_counts().to_dict(),'fit_n':int(out.model_execution.eq('fitted_registered_position_once').sum()),'manual_HC3_matches':len(manual),'reference_n':{t:len(d['R']) for t,d in data.items()},'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'results_SHA256':sha(R/'registered_results_long.tsv'),'all_outputs':{str(p.relative_to(ROOT)):sha(p) for p in [*R.glob('*'),*D.glob('*')]}};dump(passed,L/'execution_PASS.json');print(json.dumps({k:v for k,v in passed.items() if k!='all_outputs'},indent=2));print(out[['test_id','family','target','gene_set','status','n','A_n','G0_n','number_coverage','weight_coverage','beta','CI95_low','CI95_high','raw_P']].to_string(index=False))
