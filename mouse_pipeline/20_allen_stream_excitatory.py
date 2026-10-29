"""Stream ALL unbiased-sampled cells of 016 CA1-ProS, 017 CA3, 037 DG, 023 SUB-ProS Glut from Allen HPF ROIs (HIP + PAR-POST-PRE-SUB-ProS)
for co-expression module analysis. Raw UMI counts -> data/allen/excit_HPF_raw.h5ad"""
import fsspec,h5py,numpy as np,pandas as pd,scipy.sparse as sp,anndata as ad,time,sys
D='./data/allen/'
URL='https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/expression_matrices/Zeng-Aging-Mouse-10Xv3/20241130/Zeng-Aging-Mouse-10Xv3-raw.h5ad'
SUB=['016 CA1-ProS Glut','017 CA3 Glut','037 DG Glut','023 SUB-ProS Glut']
cl=np.load(D+'h5ad_cell_label.npy',allow_pickle=True); ip=np.load(D+'h5ad_indptr.npy')
m=pd.read_csv(D+'cell_metadata.csv',usecols=['cell_label','region_of_interest_label','donor_label','donor_sex','donor_age','donor_age_category','population_sampling','cluster_alias','library_label','umi_count','gene_count','doublet_score'])
p=pd.read_csv(D+'cluster_to_taxonomy.csv',index_col=0)
m=m.set_index('cell_label').reindex(cl).join(p,on='cluster_alias'); m['row']=np.arange(len(m))
m=m[m.region_of_interest_label.isin(['HPF - HIP','HPF - PAR-POST-PRE-SUB-ProS'])&(m.population_sampling=='unbiased')&m.subclass.isin(SUB)].sort_values('row')
print(m.shape, m.groupby(['subclass','donor_age_category']).size().to_dict(),flush=True)
f=fsspec.open(URL,'rb',block_size=2**24,cache_type='readahead').open(); X=h5py.File(f,'r')['X']
rows=m.row.values; segs=np.split(rows,np.where(np.diff(rows)>1)[0]+1); mats=[]; t=time.time()
# merge segments that are close (gap<200 rows) to reduce requests
blocks=[];cur=[segs[0]]
for s in segs[1:]:
    if s[0]-cur[-1][-1]<300 and sum(len(c) for c in cur)<4000: cur.append(s)
    else: blocks.append(cur); cur=[s]
blocks.append(cur)
for i,b in enumerate(blocks):
    lo,hi=b[0][0],b[-1][-1]; a,z=ip[lo],ip[hi+1]
    d=X['data'][a:z].astype(np.float32); ix=X['indices'][a:z].astype(np.int32)
    M=sp.csr_matrix((d,ix,ip[lo:hi+2]-a),shape=(hi-lo+1,32285))
    want=np.concatenate(b)-lo; mats.append(M[want])
    if i%20==0: print(i,len(blocks),round(time.time()-t),flush=True)
A=ad.AnnData(sp.vstack(mats).tocsr(),obs=m.drop(columns=['row']).assign(row=m.row.values),var=pd.read_csv(D+'genes.csv').set_index('gene_identifier'))
A.obs.index=A.obs.index.astype(str); A.write_h5ad(D+'excit_HPF_raw.h5ad',compression='gzip'); print('done',A.shape,round(time.time()-t))
