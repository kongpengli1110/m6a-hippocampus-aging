"""Offline verified-source 4-target sums; never re-aggregate cells or alter source."""

# Release-only safety/path adapter; see provenance/adaptation_diffs.
from pathlib import Path as _ReleasePath
import sys as _release_sys
_release_home = next(p for p in _ReleasePath(__file__).resolve().parents if (p / 'release_support.py').is_file())
_release_sys.path.insert(0, str(_release_home))
from release_support import gate, release_root, resolve_source
gate(__file__, __name__)
from pathlib import Path
import json,hashlib,re,sys,resource
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy import sparse
# Read-only fallback for the already installed XLSX reader; no environment mutation.
# Release: install openpyxl normally; no private site-packages fallback.
ROOT=Path(__file__).resolve().parents[1];PROJ=ROOT.parents[1];SRC=PROJ/'human_data/sea_hip';PRE=PROJ/'human_preflight';AN=ROOT.parent;D=ROOT/'derived';L=ROOT/'logs'
for p in [D,L,ROOT/'results',ROOT/'figures']:p.mkdir(parents=True,exist_ok=True)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,p):p.write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))
def boolean(s):
 a=s.astype(str);assert set(a.unique())<= {'True','False'},set(a.unique());return a.map({'True':True,'False':False})
if (L/'preparation_PASS.json').exists():raise SystemExit('Preparation succeeded already; reuse it.')
f=json.loads((AN/'freeze_manifest.json').read_text());
for r in f['frozen_files']:assert sha(resolve_source(r['path']))==r['sha256']
for line in (PROJ/'handoff/frozen_v1/FROZEN_md5.txt').read_text().splitlines():
 if re.match(r'^[a-f0-9]{32}\s',line):
  md5,n=line.split(maxsplit=1);assert hashlib.md5((PROJ/'handoff/frozen_v1'/Path(n).name).read_bytes()).hexdigest()==md5
qa=json.loads((SRC/'QA.json').read_text());assert qa['status']=='PASS' and qa['counts_source']=='/layers/UMIs' and qa['raw_total_umi']==6196354669
manifest=json.loads((SRC/'output_manifest.json').read_text());primary=['counts.npz','features.tsv','pseudobulk_samples.tsv','QA.json','QA_report.md','qa_gene_totals.tsv','qa_donor_totals.tsv'];
for n in primary:assert sha(SRC/n)==manifest[n]['sha256'],n
paths=[SRC/n for n in primary]+[SRC/'output_manifest.json',AN/'freeze_manifest.json',*[resolve_source(r['path']) for r in f['frozen_files']],*sorted((PROJ/'handoff/frozen_v1').glob('*')),PRE/'human_candidate_crosswalk.tsv',PRE/'derived/SEAAD_HIP_feature_hgnc_crosswalk.tsv',PRE/'raw/SEAAD_HIP_var.tsv',PRE/'raw/hgnc_complete_set.txt',PRE/'derived/SEAAD_HIP_donor_clinical_qc.tsv',PRE/'raw/SEAAD_donors_analysis_columns.tsv',PRE/'donor_supertype_qc.tsv']
GSE=AN/'gse48350_v1';gqa=json.loads((GSE/'logs/final_QA.json').read_text())
for n in ['results/registered_results_long.tsv','results/gene_selection.tsv','results/gene_aggregates_core23_1to1_log2_ratio.tsv.gz','results/gene_aggregates_METTL3cKO_HIP_down_log2_ratio.tsv.gz','results/sample_manifest.tsv','results/gene_reference_parameters.tsv','results/score_reference_parameters.tsv','results/donor_scores.tsv']:
 p=GSE/n;assert sha(p)==gqa['outputs'][n]['sha256'];paths.append(p)
paths += [GSE/'logs/final_QA.json',GSE/'scripts/fixed_pipeline.py'];inputs={str(p):sha(p) for p in paths if p.is_file()};dump({'files':inputs,'verified_at_utc':datetime.now(timezone.utc).isoformat(),'frozen_lock':f['lock_recorded_at']},ROOT/'input_manifest.json')
features=pd.read_csv(SRC/'features.tsv',sep='\t',dtype=str,keep_default_na=False);var=pd.read_csv(PRE/'raw/SEAAD_HIP_var.tsv',sep='\t',dtype=str,keep_default_na=False);assert features.gene_index.astype(int).tolist()==list(range(36601));pd.testing.assert_frame_equal(features[var.columns],var)
fm=pd.read_csv(PRE/'derived/SEAAD_HIP_feature_hgnc_crosswalk.tsv',sep='\t',dtype=str,keep_default_na=False);fm['unambiguous']=boolean(fm.unambiguous);fm['stable_symbol_conflict']=boolean(fm.stable_symbol_conflict);assert fm.ensembl_id.tolist()==features.gene_ids.tolist();assert fm.feature_symbol.tolist()==features['index'].tolist();fm['gene_index']=features.gene_index.astype(int)
x=pd.read_csv(PRE/'human_candidate_crosswalk.tsv',sep='\t',dtype=str,keep_default_na=False);x['human_mapping_unambiguous']=boolean(x.human_mapping_unambiguous);x['SEAAD_primary_feature_available']=boolean(x.SEAAD_primary_feature_available)
ids=set(x.loc[x.human_mapping_unambiguous,'current_HGNC_IDs']);cand=fm[fm.unambiguous & ~fm.stable_symbol_conflict & fm.HGNC_IDs.isin(ids)].copy();assert cand.gene_index.is_unique
# Audit, don't guess: every multi-feature candidate is explicitly the current symbol,
# previous symbol or alias of the SAME official HGNC identity in the frozen crosswalk.
h=pd.read_csv(PRE/'raw/hgnc_complete_set.txt',sep='\t',dtype=str,keep_default_na=False).set_index('hgnc_id');audits=[]
for hid,g in cand.groupby('HGNC_IDs'):
 for _,r in g.iterrows():
  off=h.loc[hid];allowed={off['symbol'],*off.prev_symbol.split('|'),*off.alias_symbol.split('|')};assert r.feature_symbol in allowed or r.mapping_route=='ensembl'
  audits.append({**r.to_dict(),'current_human_symbol':off['symbol'],'HGNC_official_Ensembl':off.ensembl_gene_id,'HGNC_previous_symbols':off.prev_symbol,'HGNC_alias_symbols':off.alias_symbol,'n_source_features_for_HGNC':len(g),'count_policy':'sum_each_distinct_preflight_feature_once;one_gene_one_weight','identity_review':'official_HGNC_and_exact_preflight_match'})
 for _,rr in x[x.current_HGNC_IDs.eq(hid)&x.human_mapping_unambiguous].iterrows():assert set(rr.SEAAD_unique_feature_ids.split(';'))==set(g.ensembl_id)
pd.DataFrame(audits).to_csv(D/'candidate_feature_identity_audit.tsv',sep='\t',index=False)
# Four multi-feature IDs: source Ensembl rows remain distinct; sum counts once into
# one current HGNC identity; do NOT add gene weight or remove features from library L.
cols=['Donor ID','Age at Death','Sex','Braak','RIN','PMI'];clinical=pd.read_csv(PRE/'derived/SEAAD_HIP_donor_clinical_qc.tsv',sep='\t',usecols=cols);orig=pd.read_csv(PRE/'raw/SEAAD_donors_analysis_columns.tsv',sep='\t',usecols=cols).set_index('Donor ID').loc[clinical['Donor ID']].reset_index();pd.testing.assert_frame_equal(clinical[cols].sort_values('Donor ID').reset_index(drop=True),orig[cols].sort_values('Donor ID').reset_index(drop=True),check_dtype=False)
assert clinical['Donor ID'].is_unique and len(clinical)==40;clinical=clinical.set_index('Donor ID').sort_index();clinical.index.name='donor';clinical['age']=pd.to_numeric(clinical['Age at Death'],errors='raise');clinical['sex']=clinical.Sex.map({'Female':0,'Male':1});assert clinical.sex.notna().all();clinical['age10']=(clinical.age-80)/10
braak={'0':0,'I':1,'II':2,'III':3,'IV':4,'V':5,'VI':6};b=clinical.Braak.str.replace(r'^Braak\s+','',regex=True);assert set(b.dropna())<=set(braak);clinical['Braak_numeric']=b.map(braak)
for c in ['PMI','RIN']:clinical[c]=pd.to_numeric(clinical[c],errors='raise');assert np.isfinite(clinical[c]).all()
clinical.to_csv(D/'clinical_covariates_verified.tsv',sep='\t')
samples=pd.read_csv(SRC/'pseudobulk_samples.tsv',sep='\t');assert samples.pseudobulk_index.tolist()==list(range(3402));assert not samples.duplicated(['donor','supertype']).any();counts=sparse.load_npz(SRC/'counts.npz');assert sparse.isspmatrix_csr(counts) and counts.shape==(3402,36601) and counts.dtype==np.int64;assert counts.nnz==41420887 and counts.data.min()>=0;assert int(counts.sum())==6196354669
assert np.array_equal(np.asarray(counts.sum(axis=1)).ravel(),samples.totalUMI.to_numpy())
prior=pd.read_csv(PRE/'donor_supertype_qc.tsv',sep='\t');targets={'CA3_1':['CA3_1'],'Sub_1+Sub_2':['Sub_1','Sub_2'],'CA1_1':['CA1_1'],'DG_1+DG_2':['DG_1','DG_2']};donors=clinical.index.tolist();dindex={d:i for i,d in enumerate(donors)};target_rows=[];spot=[]
for target,parts in targets.items():
 select=samples[samples.supertype.isin(parts)];rows=np.array([dindex[d] for d in select.donor]);idx=select.pseudobulk_index.to_numpy();grouper=sparse.csr_matrix((np.ones(len(idx),dtype=np.int64),(rows,idx)),shape=(40,3402));combined=(grouper@counts).toarray();assert combined.dtype==np.int64 and combined.shape==(40,36601)
 for d in donors:
  i=dindex[d];s=select[select.donor.eq(d)];umi=int(combined[i].sum());n=int(s.n_nuclei.sum());assert umi==int(s.totalUMI.sum());pn=prior[(prior.target==target)&prior['Donor ID'].eq(d)].n_nuclei.iloc[0];assert n==pn
  libs=sorted(set(';'.join(s.libraries).split(';'))-{''});sub2=int(s.loc[s.supertype.eq('Sub_2'),'n_nuclei'].sum());r={'target':target,'donor':d,'n_nuclei':n,'totalUMI_all36601':umi,'n_libraries':len(libs),'libraries':';'.join(libs),'source_pseudobulk_indices':';'.join(map(str,s.pseudobulk_index)),'Sub_2_nuclei_fraction':sub2/n if n>0 and target=='Sub_1+Sub_2' else np.nan,'has_observed_cells':n>0}
  for threshold in [10,20,30]:r['eligible_ge'+str(threshold)]=n>=threshold and umi>0
  target_rows.append(r)
  if d in [donors[0],'H20.33.046']:
   manual=np.zeros(36601,dtype=np.int64)
   for j in s.pseudobulk_index:manual+=counts.getrow(j).toarray().ravel()
   assert np.array_equal(manual,combined[i]);spot.append({'target':target,'donor':d,'n_component_rows':len(s),'totalUMI':umi,'all36601_manual_sum_match':True,'libraries_combined_upstream_only':True})
 np.savez_compressed(D/(target.replace('+','_')+'_all_feature_counts.npz'),counts=combined,donors=np.array(donors),gene_ids=features.gene_ids.to_numpy(dtype=str));del combined
pd.DataFrame(target_rows).to_csv(D/'target_donor_manifest.tsv',sep='\t',index=False);pd.DataFrame(spot).to_csv(D/'component_sum_spotchecks.tsv',sep='\t',index=False);del counts
assert all(sha(Path(p))==hsh for p,hsh in inputs.items());result={'status':'PASS','completed_at_utc':datetime.now(timezone.utc).isoformat(),'source_QA_accepted':True,'npz_allgene_sum':6196354669,'source_shapes':[3402,36601],'target_rows':160,'all_feature_n':36601,'no_source_H5AD_reread':True,'multi_feature_candidate_HGNC_IDs':sorted(cand.loc[cand.HGNC_IDs.duplicated(False),'HGNC_IDs'].unique()),'multi_feature_rule':'identity-audited sum of frozen preflight mapped source features; gene weight once; full feature library denominator unchanged','source_files_unchanged':True,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'derived_files':{p.name:sha(p) for p in D.glob('*')}};dump(result,L/'preparation_PASS.json');print(json.dumps(result,indent=2));print(pd.DataFrame(target_rows).groupby('target')[['eligible_ge10','eligible_ge20','eligible_ge30']].sum().to_string())
