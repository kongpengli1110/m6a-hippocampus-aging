"""Stream only HPF-HIP cells from Allen Zeng-Aging-Mouse-10Xv3 raw h5ad (S3 public, 13.8 GB) via HTTP range reads.
Outputs: pseudobulk raw counts per donor x population_sampling x subclass; stratified cell subsample (raw counts) as h5ad."""
import fsspec,h5py,time,numpy as np,pandas as pd,scipy.sparse as sp,sys,multiprocessing as mp,anndata as ad
D='./data/allen/'
URL='https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/expression_matrices/Zeng-Aging-Mouse-10Xv3/20241130/Zeng-Aging-Mouse-10Xv3-raw.h5ad'
REGION=sys.argv[1] if len(sys.argv)>1 else 'HPF - HIP'
TAG=REGION.split(' - ')[-1].split('-')[0]
NG=32285
def load_meta():
    cl=np.load(D+'h5ad_cell_label.npy',allow_pickle=True)
    m=pd.read_csv(D+'cell_metadata.csv',usecols=['cell_label','region_of_interest_label','donor_label','donor_sex','donor_age','donor_age_category','population_sampling','cluster_alias','library_label','umi_count','gene_count','doublet_score'])
    p=pd.read_csv(D+'cluster_to_taxonomy.csv',index_col=0)
    m=m.set_index('cell_label').reindex(cl).join(p,on='cluster_alias')
    m['row']=np.arange(len(m))
    return m[m.region_of_interest_label==REGION].copy()
def work(args):
    rows,gkeys,keepmask,wid=args
    f=fsspec.open(URL,'rb',block_size=2**24,cache_type='readahead').open(); h=h5py.File(f,'r'); X=h['X']
    ip=np.load(D+'h5ad_indptr.npy')
    ng=NGRP; pb=np.zeros((ng,NG)); keep=[]; t=time.time()
    # chunk over contiguous runs of rows
    CH=1500
    for s in range(0,len(rows),CH):
        r=rows[s:s+CH]; 
        # split into contiguous segments
        segs=np.split(np.arange(len(r)),np.where(np.diff(r)>1)[0]+1)
        for sg in segs:
            rr=r[sg]; a,b=ip[rr[0]],ip[rr[-1]+1]
            d=X['data'][a:b].astype(np.float32); ix=X['indices'][a:b].astype(np.int32)
            lip=(ip[rr[0]:rr[-1]+2]-a)
            M=sp.csr_matrix((d,ix,lip),shape=(len(rr),NG))
            G=sp.csr_matrix((np.ones(len(rr)),(gkeys[s+sg],np.arange(len(rr)))),shape=(ng,len(rr)))
            pb+= (G@M).toarray()
            km=keepmask[s+sg]
            if km.any(): keep.append((rr[km],M[km]))
        print(wid,s+len(r),'/',len(rows),round(time.time()-t),flush=True)
    return pb,keep
if __name__=='__main__':
    m=load_meta(); print(REGION,m.shape,flush=True)
    m=m.sort_values('row')
    m['gkey_str']=m.donor_label+'|'+m.population_sampling+'|'+m.subclass
    cats=pd.Categorical(m.gkey_str); m['gkey']=cats.codes
    global NGRP; NGRP=len(cats.categories)
    rng=np.random.default_rng(1)
    # stratified subsample: up to 600 cells per subclass x age (unbiased sampling only)
    m['keep']=False
    for _,g in m[m.population_sampling=='unbiased'].groupby(['subclass','donor_age_category']):
        k=rng.choice(g.index.values,size=min(600,len(g)),replace=False); m.loc[k,'keep']=True
    print('subsample',m.keep.sum(),flush=True)
    rows=m.row.values; half=len(rows)//2
    jobs=[(rows[:half],m.gkey.values[:half],m.keep.values[:half],0),(rows[half:],m.gkey.values[half:],m.keep.values[half:],1)]
    import subprocess,pickle
    if len(sys.argv)>2:  # worker mode
        w=int(sys.argv[2]); r=work(jobs[w]); pickle.dump(r,open(D+f'_tmp_{TAG}_{w}.pkl','wb'),protocol=4); sys.exit(0)
    ps=[subprocess.Popen([sys.executable,__file__,REGION,str(w)]) for w in (0,1)]
    assert all(p.wait()==0 for p in ps)
    res=[pickle.load(open(D+f'_tmp_{TAG}_{w}.pkl','rb')) for w in (0,1)]
    pb=res[0][0]+res[1][0]
    genes=pd.read_csv(D+'genes.csv')
    pbdf=pd.DataFrame(pb.astype(np.int64),index=cats.categories,columns=genes.gene_identifier)
    pbdf.to_parquet(D+f'pseudobulk_{TAG}_donor_sampling_subclass.parquet')
    gm=m.groupby('gkey_str').agg(donor_label=('donor_label','first'),population_sampling=('population_sampling','first'),subclass=('subclass','first'),cls=('class','first'),donor_sex=('donor_sex','first'),donor_age=('donor_age','first'),age=('donor_age_category','first'),n_cells=('row','size')).reindex(cats.categories)
    gm.to_csv(D+f'pseudobulk_{TAG}_meta.csv')
    keep=[k for r in res for k in r[1]]
    krows=np.concatenate([k[0] for k in keep]); KM=sp.vstack([k[1] for k in keep]).tocsr()
    obs=m.set_index('row').loc[krows].drop(columns=['keep','gkey'])
    obs.index=obs.index.astype(str)
    A=ad.AnnData(KM,obs=obs,var=genes.set_index('gene_identifier'))
    A.write_h5ad(D+f'subsample_{TAG}_raw.h5ad',compression='gzip')
    print('done',A.shape,pbdf.shape,flush=True)
    import os; [os.remove(D+f'_tmp_{TAG}_{w}.pkl') for w in (0,1)]
