"""All-gene donor x supertype pseudobulk from backed CSR /layers/UMIs only.

Bounded 512-cell CSR chunks; int64 summation; dense only for ~3400 x 36601
pseudobulk accumulator (~0.93 GiB), NEVER cells x genes. Final output CSR.
No normalization, filtering, inferential statistics, gene sets or mouse data.
"""

# Release-only safety/path adapter; see provenance/adaptation_diffs.
from pathlib import Path as _ReleasePath
import sys as _release_sys
_release_home = next(p for p in _ReleasePath(__file__).resolve().parents if (p / 'release_support.py').is_file())
_release_sys.path.insert(0, str(_release_home))
from release_support import gate, release_root, resolve_source
gate(__file__, __name__)
from pathlib import Path
import sys,os,json,time,hashlib,datetime,resource,platform,importlib.metadata,fcntl,gc
import h5py,numpy as np,pandas as pd,scipy,scipy.sparse as sp
ROOT=Path(__file__).resolve().parents[1]
PREFLIGHT=ROOT.parents[1]/'human_preflight'
RAW=ROOT/'raw/SEAAD_HIP_RNAseq_final-nuclei.2026-06-22.h5ad'
SEED=20261012;CHUNK=512
utc=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
def digest(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  while b:=f.read(8*1024**2):h.update(b)
 return h.hexdigest()
def write_json(path,obj):Path(path).write_text(json.dumps(obj,indent=2,ensure_ascii=False))
def field(a):
 if isinstance(a,h5py.Group):
  assert set(a)>=set(['categories','codes'])
  cats=a['categories'].asstr()[:] if a['categories'].dtype.kind in 'OS' else a['categories'][:]
  return pd.Categorical.from_codes(a['codes'][:],cats)
 return a.asstr()[:] if a.dtype.kind in 'OS' else a[:]
def jsonable(x):
 if isinstance(x,np.generic):return x.item()
 if isinstance(x,np.ndarray):return x.tolist()
 return str(x)
def fingerprint(indices,values):
 return hashlib.sha256(np.asarray(indices,dtype='<i4').tobytes()+np.asarray(values,dtype='<i8').tobytes()).hexdigest()
def accumulate_chunk(block,sample_codes,target):
 """Sparse grouping of a bounded cell block; dense target is pseudobulk only."""
 ug,inverse=np.unique(sample_codes,return_inverse=True)
 selector=sp.csr_matrix((np.ones(block.shape[0],dtype=np.int64),(inverse,np.arange(block.shape[0]))),shape=(len(ug),block.shape[0]))
 summed=selector@block
 summed.sum_duplicates();summed.sort_indices()
 assert summed.has_canonical_format
 reduced=summed.tocoo(copy=False)
 # Distinct group x gene coordinates permit advanced-indexed += without loss.
 target[ug[reduced.row],reduced.col]+=reduced.data

def checked_equal(got,expected,cols,label):
 a=got[cols].copy();b=expected[cols].copy()
 for c in cols:
  if c not in ['n_nuclei','n_libraries']:a[c]=a[c].astype(str);b[c]=b[c].astype(str)
 a=a.sort_values(cols).reset_index(drop=True);b=b.sort_values(cols).reset_index(drop=True)
 pd.testing.assert_frame_equal(a,b,check_dtype=False,check_exact=True)
 checks[label]={'status':'PASS','rows':len(a),'columns':cols}
lock=open(ROOT/'aggregate.lock','a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
assert not (ROOT/'QA.json').exists(),'Completed QA exists: do not overwrite without review.'
start=time.monotonic();checks={};qa={'status':'RUNNING','started_utc':utc(),'seed':SEED,'cell_chunk_rows':CHUNK,'counts_source':'/layers/UMIs','checks':checks}
write_json(ROOT/'logs/aggregate_status.json',qa)
def record_uncaught_exception(kind,error,tb):
 qa.update(status='FAILED',failed_utc=utc(),error_type=kind.__name__,error=str(error))
 try:write_json(ROOT/'logs/aggregate_status.json',qa)
 except Exception:pass
 sys.__excepthook__(kind,error,tb)
sys.excepthook=record_uncaught_exception
manifest=json.loads((ROOT/'download_manifest.json').read_text())
assert manifest['status']=='complete' and RAW.stat().st_size==manifest['source_identity']['content_length']==4976204965
# Reverify raw bytes independently of download process before using them.
assert digest(RAW)==manifest['sha256'];checks['source_sha256']='PASS'
inspect=json.loads((ROOT/'logs/source_inspection.json').read_text())['http_headers']
assert manifest['source_identity']=={'content_length':int(inspect['Content-Length']),'etag':inspect['ETag'],'last_modified':inspect['Last-Modified'],'version_id':inspect.get('x-amz-version-id')}
checks['same_object_identity_as_layer_inspection']='PASS'
inputs={
 'var':PREFLIGHT/'raw/SEAAD_HIP_var.tsv',
 'donor_libraries':PREFLIGHT/'derived/SEAAD_HIP_donor_libraries.tsv',
 'donor_supertype_all':PREFLIGHT/'derived/SEAAD_HIP_donor_supertype_all.tsv',
 'target_qc':PREFLIGHT/'donor_supertype_qc.tsv'}
input_hashes={k:{'path':str(v),'sha256':digest(v)} for k,v in inputs.items()}
write_json(ROOT/'provenance/prior_input_manifest.json',input_hashes)
env={'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'hdf5':h5py.version.hdf5_version,'packages':{d.metadata['Name']:d.version for d in importlib.metadata.distributions()},'raw_sha256':manifest['sha256']}
write_json(ROOT/'environment.json',env)
(ROOT/'requirements.lock.txt').write_text(''.join(f'{k}=={v}\n' for k,v in sorted(env['packages'].items())))
with h5py.File(RAW,'r',rdcc_nbytes=16*1024**2) as h:
 assert 'UMIs' in h['layers'];a=h['layers/UMIs']
 assert a.attrs['encoding-type']=='csr_matrix'
 ncell,ngene=map(int,a.attrs['shape']);assert (ncell,ngene)==(208843,36601)
 assert h['uns/X_normalization'].asstr()[()]=='ln(UP10K+1)'
 v=pd.DataFrame({k:field(h['var'][k]) for k in h['var']})
 oldv=pd.read_csv(inputs['var'],sep='\t',keep_default_na=False)
 pd.testing.assert_frame_equal(v.astype(str),oldv.astype(str),check_exact=True)
 assert len(v)==ngene and v['gene_ids'].notna().all() and v['gene_ids'].is_unique
 assert v['index'].notna().all() and v['index'].is_unique
 checks['gene_order_all_36601_and_metadata']='PASS'
 features=v.astype(str).copy();features.insert(0,'gene_index',np.arange(ngene,dtype=np.int64))
 features.to_csv(ROOT/'features.tsv',sep='\t',index=False)
 obsindex=h['obs'].attrs['_index']
 names=[obsindex,'sample_id','Donor ID','library_prep','Supertype','Subclass','method','Number of UMIs','Genes detected']
 o=pd.DataFrame({k:field(h['obs'][k]) for k in names})
 assert len(o)==ncell and o.notna().all().all(),'Missing IDs or UMI/gene metadata'
 for k in names[:7]:o[k]=o[k].astype(str);assert (o[k].str.len()>0).all()
 assert o[obsindex].is_unique
 assert o['Donor ID'].nunique()==40 and o['library_prep'].nunique()==41
 assert (o.groupby('library_prep')['Donor ID'].nunique()==1).all()
 assert (o.groupby('sample_id')['Donor ID'].nunique()==1).all()
 assert (o.groupby('Supertype')['Subclass'].nunique()==1).all()
 obsum=pd.to_numeric(o['Number of UMIs']).to_numpy()
 obsgenes=pd.to_numeric(o['Genes detected']).to_numpy()
 for vals in [obsum,obsgenes]:assert np.isfinite(vals).all() and (vals>=0).all() and (vals==np.floor(vals)).all()
 obsum=obsum.astype(np.int64);obsgenes=obsgenes.astype(np.int64)
 lib=o.groupby(['Donor ID','library_prep','method'],observed=True).size().rename('n_nuclei').reset_index()
 oldlib=pd.read_csv(inputs['donor_libraries'],sep='\t')
 checked_equal(lib,oldlib,['Donor ID','library_prep','method','n_nuclei'],'prior_libraries_every_row')
 groups=o.groupby(['Donor ID','Supertype','Subclass'],observed=True,sort=True).agg(n_nuclei=('library_prep','size'),n_libraries=('library_prep','nunique'),libraries=('library_prep',lambda x:';'.join(sorted(set(x))))).reset_index()
 oldgroups=pd.read_csv(inputs['donor_supertype_all'],sep='\t')
 checked_equal(groups,oldgroups,['Donor ID','Supertype','Subclass','n_nuclei','n_libraries'],'prior_donor_supertype_every_row')
 keys=list(zip(groups['Donor ID'],groups['Supertype']));assert len(set(keys))==len(keys)
 lookup={k:i for i,k in enumerate(keys)}
 codes=np.array([lookup[k] for k in zip(o['Donor ID'],o['Supertype'])],dtype=np.int32)
 nsamp=len(groups);assert np.array_equal(np.bincount(codes,minlength=nsamp),groups.n_nuclei.to_numpy())
 target=pd.read_csv(inputs['target_qc'],sep='\t')
 for _,r in target.iterrows():
  count=groups.loc[(groups['Donor ID']==r['Donor ID']) & groups.Supertype.isin(r.included_supertypes.split(';')),'n_nuclei'].sum()
  assert int(count)==int(r.n_nuclei)
 checks['prior_target_qc_including_combined_counts']={'status':'PASS','rows':len(target),'note':'QA nuclei counts only; combined expression rows NOT created'}
 checks['nuclei_and_donor_library_ids']='PASS'
 meta=o[[obsindex,'sample_id','Donor ID','library_prep','Supertype','Subclass','method']].copy()
 meta.insert(0,'cell_index',np.arange(ncell));meta['pseudobulk_index']=codes
 meta.to_csv(ROOT/'cell_to_pseudobulk.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
 del meta
 ptr=a['indptr'][:].astype(np.int64)
 nnz=int(a['data'].shape[0]);assert a['indices'].shape==(nnz,)
 assert ptr.shape==(ncell+1,) and ptr[0]==0 and ptr[-1]==nnz and (np.diff(ptr)>=0).all()
 qa.update(source_shape=[ncell,ngene],source_stored_entries=nnz,source_dtype=str(a['data'].dtype),donors=40,libraries=41,observed_supertypes=int(o.Supertype.nunique()),pseudobulk_shape=[nsamp,ngene],pseudobulk_accumulator_bytes=nsamp*ngene*8,unobserved_donor_supertype_combinations=40*int(o.Supertype.nunique())-nsamp)
 rng=np.random.default_rng(SEED)
 cell_spots=sorted(map(int,rng.choice(ncell,size=24,replace=False)))
 random_groups=sorted(map(int,rng.choice(nsamp,size=6,replace=False)))
 double_candidates=groups[(groups['Donor ID']=='H20.33.046') & (groups.n_libraries==2)]
 assert len(double_candidates)>0
 # Include a two-library group without replacing six unbiased random selections.
 double_id=int(double_candidates.sort_values('n_nuclei').index[len(double_candidates)//2])
 group_spots=sorted(set(random_groups+[double_id]));cell_fingerprints={}
 pb=np.zeros((nsamp,ngene),dtype=np.int64)
 raw_gene=np.zeros(ngene,dtype=np.int64);cell_umi=np.zeros(ncell,dtype=np.int64);cell_ngene=np.zeros(ncell,dtype=np.int64)
 raw_total=0;maximum=0;minimum=np.inf;zero_stored=0;noncanonical_chunks=0
 for lo in range(0,ncell,CHUNK):
  hi=min(lo+CHUNK,ncell);left,right=int(ptr[lo]),int(ptr[hi])
  values=a['data'][left:right];idx=a['indices'][left:right]
  assert np.isfinite(values).all(),f'Nonfinite counts at row block {lo}'
  assert (values>=0).all(),f'Negative counts at row block {lo}'
  assert (values==np.floor(values)).all(),f'Noninteger counts at row block {lo}'
  assert values.size==right-left and idx.size==values.size
  assert (idx>=0).all() and (idx<ngene).all()
  if values.size:
   maximum=max(maximum,int(values.max()));minimum=min(minimum,int(values.min()))
   assert maximum<2**24,'float32 UMI precision requires review at >=2^24'
  zero_stored+=int((values==0).sum())
  vint=values.astype(np.int64);del values
  raw_total+=int(vint.sum(dtype=np.int64));assert raw_total<2**63
  block=sp.csr_matrix((vint,idx,ptr[lo:hi+1]-left),shape=(hi-lo,ngene),copy=False)
  if not block.has_canonical_format:noncanonical_chunks+=1;block.sum_duplicates();block.sort_indices()
  block.eliminate_zeros()
  cell_umi[lo:hi]=np.asarray(block.sum(axis=1)).ravel()
  cell_ngene[lo:hi]=np.diff(block.indptr)
  assert np.array_equal(cell_umi[lo:hi],obsum[lo:hi]),f'Raw layer UMI sums differ from Number of UMIs at rows {lo}:{hi}; stop and review.'
  assert np.array_equal(cell_ngene[lo:hi],obsgenes[lo:hi]),f'Raw layer detected genes differ from metadata at rows {lo}:{hi}'
  raw_gene+=np.asarray(block.sum(axis=0)).ravel()
  accumulate_chunk(block,codes[lo:hi],pb)
  for ci in cell_spots:
   if lo<=ci<hi:
    rr=ci-lo;s,e=block.indptr[rr:rr+2]
    cell_fingerprints[ci]=fingerprint(block.indices[s:e],block.data[s:e])
  if lo==0 or hi==ncell or (lo//CHUNK)%20==0:
   qa.update(cells_processed=hi,raw_total_umi_so_far=raw_total,elapsed_seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
   write_json(ROOT/'logs/aggregate_status.json',qa)
   print(f'{utc()} cells={hi}/{ncell} raw_UMI={raw_total} elapsed_s={time.monotonic()-start:.1f} peak_RSS_MiB={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.1f}',flush=True)
  # ug/inverse are local to accumulate_chunk(), not this chunk-loop scope.
  del block,vint,idx
 # Exact integer conservation on total, each gene, each sample.
 assert int(pb.sum(dtype=np.int64))==raw_total==int(raw_gene.sum())==int(cell_umi.sum())==int(obsum.sum())
 assert np.array_equal(pb.sum(axis=0,dtype=np.int64),raw_gene)
 expected_sample=np.zeros(nsamp,dtype=np.int64);np.add.at(expected_sample,codes,cell_umi)
 assert np.array_equal(pb.sum(axis=1,dtype=np.int64),expected_sample)
 checks.update(all_input_counts_finite_nonnegative_integer='PASS',source_csr_structure_and_gene_bounds='PASS',all_cell_umi_metadata_exact='PASS',all_cell_detected_genes_metadata_exact='PASS',all_total_umi_conservation='PASS',every_gene_total_conservation='PASS',every_pseudobulk_total_conservation='PASS')
 qa.update(raw_total_umi=raw_total,raw_min_stored_count=int(minimum),raw_max_stored_count=maximum,explicit_zero_stored_entries=zero_stored,noncanonical_source_chunks=noncanonical_chunks)
 # Independent, direct per-cell reread for 7 groups: no chunk selector multiplication.
 spotgroups=[]
 for gi in group_spots:
  cells=np.flatnonzero(codes==gi);reference=np.zeros(ngene,dtype=np.int64)
  for ci in cells:
   s,e=map(int,ptr[ci:ci+2]);ix=a['indices'][s:e];val=a['data'][s:e].astype(np.int64)
   np.add.at(reference,ix,val)
  assert np.array_equal(reference,pb[gi]),f'Independent full-gene spotcheck failed for {keys[gi]}'
  spotgroups.append({'pseudobulk_index':gi,'donor':keys[gi][0],'supertype':keys[gi][1],'n_nuclei':len(cells),'n_libraries':int(groups.iloc[gi].n_libraries),'totalUMI':int(reference.sum()),'all_genes_equal':True,'selection':'fixed_seed_random' if gi in random_groups else 'forced_two_library_donor'})
 spotcells=[]
 for ci in cell_spots:
  s,e=map(int,ptr[ci:ci+2]);ix=a['indices'][s:e];val=a['data'][s:e].astype(np.int64)
  row=sp.csr_matrix((val,ix,np.array([0,len(ix)])),shape=(1,ngene));row.sum_duplicates();row.sort_indices();row.eliminate_zeros()
  fp=fingerprint(row.indices,row.data)
  assert fp==cell_fingerprints[ci] and int(row.sum())==cell_umi[ci] and row.nnz==cell_ngene[ci]
  spotcells.append({'cell_index':ci,'cell_id':str(o.iloc[ci][obsindex]),'donor':str(o.iloc[ci]['Donor ID']),'supertype':str(o.iloc[ci].Supertype),'totalUMI':int(row.sum()),'genes_detected':row.nnz,'counts_fingerprint_sha256':fp,'all_genes_equal':True})
 pd.DataFrame(spotgroups).to_csv(ROOT/'qa_group_spotchecks.tsv',sep='\t',index=False)
 pd.DataFrame(spotcells).to_csv(ROOT/'qa_cell_spotchecks.tsv',sep='\t',index=False)
 checks['independent_full_gene_group_spotchecks']={'status':'PASS','n':len(spotgroups),'includes_two_library_donor':True,'seed':SEED}
 checks['independent_cell_reread_spotchecks']={'status':'PASS','n':len(spotcells),'seed':SEED}
 counts=sp.csr_matrix(pb);del pb;gc.collect()
 assert counts.dtype==np.int64 and counts.has_canonical_format and counts.nnz>0 and counts.data.min()>0
 # Save unconditionally only after every count-level check passes.
 sp.save_npz(ROOT/'counts.npz',counts,compressed=True)
 samples=groups.rename(columns={'Donor ID':'donor','Supertype':'supertype','Subclass':'subclass'}).copy()
 samples.insert(0,'pseudobulk_id',[f'{d}__{s}' for d,s in keys]);samples.insert(0,'pseudobulk_index',np.arange(nsamp,dtype=np.int64))
 assert samples.pseudobulk_id.is_unique
 samples['totalUMI']=expected_sample
 samples.to_csv(ROOT/'pseudobulk_samples.tsv',sep='\t',index=False)
 genes=pd.DataFrame({'gene_index':np.arange(ngene),'gene_ids':v.gene_ids,'index':v['index'],'raw_totalUMI':raw_gene})
 genes.to_csv(ROOT/'qa_gene_totals.tsv',sep='\t',index=False)
 donor_summary=samples.groupby('donor').agg(n_nuclei=('n_nuclei','sum'),n_observed_supertypes=('supertype','nunique'),totalUMI=('totalUMI','sum')).reset_index()
 donor_summary['n_libraries']=donor_summary.donor.map(o.groupby('Donor ID').library_prep.nunique())
 donor_summary.to_csv(ROOT/'qa_donor_totals.tsv',sep='\t',index=False)
 super_summary=samples.groupby(['supertype','subclass']).agg(n_nuclei=('n_nuclei','sum'),n_observed_donors=('donor','nunique'),totalUMI=('totalUMI','sum')).reset_index()
 super_summary.to_csv(ROOT/'qa_supertype_totals.tsv',sep='\t',index=False)
 assert donor_summary.loc[donor_summary.donor=='H20.33.046','n_nuclei'].iloc[0]==4501
 assert donor_summary.loc[donor_summary.donor=='H20.33.046','n_libraries'].iloc[0]==2
 qa.update(output_nnz=int(counts.nnz),output_dtype=str(counts.dtype),output_format='csr',output_nonzero_min=int(counts.data.min()),output_nonzero_max=int(counts.data.max()),genes_with_zero_total=int((raw_gene==0).sum()),pseudobulk_zero_umi_rows=int((expected_sample==0).sum()))
 # Save a compact structure/data fingerprint, release memory, then reload deliverable.
 original_fp={k:hashlib.sha256(np.asarray(getattr(counts,k)).tobytes()).hexdigest() for k in ['data','indices','indptr']}
 del counts;gc.collect()
 loaded=sp.load_npz(ROOT/'counts.npz')
 assert loaded.shape==(nsamp,ngene) and loaded.dtype==np.int64 and loaded.has_canonical_format
 assert all(hashlib.sha256(np.asarray(getattr(loaded,k)).tobytes()).hexdigest()==v for k,v in original_fp.items())
 assert int(loaded.sum())==raw_total
 assert np.array_equal(np.asarray(loaded.sum(axis=0)).ravel(),raw_gene)
 assert np.array_equal(np.asarray(loaded.sum(axis=1)).ravel(),expected_sample)
 pd.testing.assert_frame_equal(pd.read_csv(ROOT/'pseudobulk_samples.tsv',sep='\t'),samples,check_dtype=False)
 pd.testing.assert_frame_equal(pd.read_csv(ROOT/'features.tsv',sep='\t',keep_default_na=False),features,check_dtype=False)
 checks['saved_npz_and_tsv_reloaded_exact']='PASS'
 del loaded
# Ensure every read-only input is unchanged.
for k,p in inputs.items():assert digest(p)==input_hashes[k]['sha256']
checks['prior_inputs_unmodified']='PASS'
qa.update(status='PASS',completed_utc=utc(),elapsed_seconds=time.monotonic()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,raw_sha256=manifest['sha256'],no_inference_or_gene_set_scoring=True,no_mouse_work=True,matrix_orientation='rows=observed donor x supertype; columns=all source genes in original order',zero_cell_group_policy='absent combinations omitted, not fabricated zero-expression samples',combined_type_policy='Sub and DG component supertypes retained separately; no duplicated combined count rows')
write_json(ROOT/'QA.json',qa);write_json(ROOT/'logs/aggregate_status.json',qa)
outputs=['counts.npz','features.tsv','pseudobulk_samples.tsv','cell_to_pseudobulk.tsv.gz','QA.json','environment.json','requirements.lock.txt','qa_gene_totals.tsv','qa_donor_totals.tsv','qa_supertype_totals.tsv','qa_group_spotchecks.tsv','qa_cell_spotchecks.tsv']
write_json(ROOT/'output_manifest.json',{n:{'bytes':(ROOT/n).stat().st_size,'sha256':digest(ROOT/n)} for n in outputs})
print('QA PASS',json.dumps(qa,ensure_ascii=False),flush=True)
# Finish human-readable documentation automatically in the same successful pipeline.
import subprocess
subprocess.run([sys.executable,str(ROOT/'scripts/finalize_release_manifest.py'),'--execute'],check=True)
