"""Replication sc sets -> per-sample x cell-type pseudobulk (raw counts) + per-type mean CPM (for EWCE).
B2 GSE161340 (Ogrodnik 2021): hippocampus, 4 vs 24 mo F; sc (3'v2) and sn (5') libraries, each library = pool of 2 mice -> 2 vs 2 per modality.
B3 GSE233363 (DG scRNA): 3 / 9-11 / 16-21(x2) mo -> 1 library per group (no replication: descriptive only, pooled old vs young).
B4 GSE129788 (Ximerakis 2019, whole brain, 8 young 2-3 mo vs 8 old 21-22 mo): GEO matrices are ln(1+CPM*1e4)-normalised; raw counts
recovered exactly as (exp(x)-1)*nUMI/1e4 (verified: integers, sum==nUMI). Author cell-type labels used.
Annotation for B2/B3 (no author labels on GEO): per-cell marker scoring (scanpy score_genes) with Allen-HPF top markers, assign argmax; neurons split CA1/CA3/DG/SUB by Allen markers."""
import scanpy as sc,anndata as ad,pandas as pd,numpy as np,scipy.sparse as sp,scipy.io as sio,gzip,sys
D='./'; R=D+'data/'
sc.settings.verbosity=0
spec=pd.read_csv(D+'analysis/tables/allen_HPF_specificity_EWCE.csv',index_col=0); mean=pd.read_csv(D+'analysis/tables/allen_HPF_meanCPM_by_type.csv',index_col=0)
TYPES=['CA1','CA3','DG','SUB','GABA','Astro','Oligo','OPC','Microglia','Endo','Mural','VLMC','BAM','DG_IMN']
def markers(t,n=40):
    if t not in spec: return []
    s=spec[t][(mean[t]>=20)]; return list(s.sort_values(ascending=False).head(n).index)
MK={t:markers(t) for t in TYPES if t in spec}
def annotate(A):
    A.layers['counts']=A.X.copy(); sc.pp.normalize_total(A,1e4); sc.pp.log1p(A)
    for t,g in MK.items():
        g=[x for x in g if x in A.var_names]
        if len(g)>=10: sc.tl.score_genes(A,g,score_name='sc_'+t,ctrl_size=100,random_state=0)
    S=A.obs[[c for c in A.obs if c.startswith('sc_')]]
    A.obs['celltype']=S.idxmax(axis=1).str[3:]; srt=np.sort(S.values,1); A.obs['margin']=srt[:,-1]-srt[:,-2]
    A.obs.loc[A.obs.margin<0.1,'celltype']='ambiguous'
    A.X=A.layers['counts']; return A
def read10x(d,gz,genes_file):
    m=sio.mmread(gzip.open(d+'/matrix.mtx.gz') if gz else d+'/matrix.mtx').T.tocsr()
    b=pd.read_csv(d+('/barcodes.tsv.gz' if gz else '/barcodes.tsv'),header=None)[0]
    g=pd.read_csv(d+'/'+genes_file,sep='\t',header=None)
    A=ad.AnnData(m.astype(np.float32),obs=pd.DataFrame(index=b.values),var=pd.DataFrame({'ensembl':g[0].values},index=g[1].values)); A.var_names_make_unique(); return A
def qc(A,mito='mt-'):
    A.var['mt']=A.var_names.str.startswith(mito); sc.pp.calculate_qc_metrics(A,qc_vars=['mt'],inplace=True,percent_top=None,log1p=False)
    return A[(A.obs.n_genes_by_counts>=500)&(A.obs.pct_counts_mt<10)].copy()
def pseudobulk(A,key,typ='celltype'):
    A.obs['pb']=A.obs[key].astype(str)+'|'+A.obs[typ].astype(str)
    G=pd.get_dummies(A.obs.pb); M=sp.csr_matrix(G.values.T.astype(np.float32))@A.X
    P=pd.DataFrame(M.toarray() if sp.issparse(M) else M,index=G.columns,columns=A.var_names)
    meta=A.obs.groupby('pb').agg(sample=(key,'first'),celltype=(typ,'first'),n_cells=(key,'size')).loc[P.index]; return P,meta
out={}
which=sys.argv[1:] or ['B2','B3','B4']
if 'B2' in which:
    for mod in ['sc','sn']:
        As=[]
        for lib,age in [('OBrain1','old'),('OBrain3','old'),('YBrain1','young'),('YBrain3','young')]:
            d=R+f'B2_GSE161340/raw/{mod}_{lib}_mm10'+('_pre_mRNA' if mod=='sn' else '')
            A=read10x(d,mod=='sn','features.tsv.gz' if mod=='sn' else 'genes.tsv'); A.obs['sample']=f'{mod}_{lib}'; A.obs['age']=age; As.append(qc(A))
        A=ad.concat(As,join='inner',index_unique='-'); A=annotate(A)
        P,M=pseudobulk(A,'sample'); M['age']=M['sample'].map(dict(zip(A.obs['sample'],A.obs.age)))
        P.to_parquet(R+f'B2_GSE161340/pb_{mod}.parquet'); M.to_csv(R+f'B2_GSE161340/pb_{mod}_meta.csv')
        print('B2',mod,A.shape,A.obs.groupby(['sample','celltype']).size().unstack(fill_value=0).to_string(),flush=True)
if 'B3' in which:
    As=[]
    for gsm,s,age in [('GSM7426316_P23360_1001_S1','DG_young','young'),('GSM7426317_P23360_1002_S2','DG_middle','middle'),('GSM7426318_P23360_1003_S3','DG_old1','old'),('GSM7426319_P23360_1004_S4','DG_old2','old')]:
        p=R+'B3_GSE233363/raw/'+gsm
        m=sio.mmread(gzip.open(p+'_matrix.mtx.gz')).T.tocsr(); b=pd.read_csv(p+'_barcodes.tsv.gz',header=None)[0]; g=pd.read_csv(p+'_features.tsv.gz',sep='\t',header=None)
        A=ad.AnnData(m.astype(np.float32),obs=pd.DataFrame(index=b.values),var=pd.DataFrame({'ensembl':g[0].values},index=g[1].values)); A.var_names_make_unique()
        A.obs['sample']=s; A.obs['age']=age; As.append(qc(A))
    A=ad.concat(As,join='inner',index_unique='-'); A=annotate(A)
    P,M=pseudobulk(A,'sample'); M['age']=M['sample'].map(dict(zip(A.obs['sample'],A.obs.age)))
    P.to_parquet(R+'B3_GSE233363/pb.parquet'); M.to_csv(R+'B3_GSE233363/pb_meta.csv')
    print('B3',A.shape,A.obs.groupby(['sample','celltype']).size().unstack(fill_value=0).to_string(),flush=True)
if 'B4' in which:
    meta=pd.read_csv(R+'B4_GSE129788/GSE129788_Supplementary_meta_data_Cell_Types_Etc.txt.gz',sep='\t',skiprows=[1],index_col=0)
    meta.index=meta.index.str.replace('Aging_mouse_brain_portal_data_','',regex=False)
    import glob; Ps=[];Ms=[]
    for f in sorted(glob.glob(R+'B4_GSE129788/raw/*.txt.gz')):
        x=pd.read_csv(f,sep='\t',index_col=0); mm=meta.reindex(x.columns)
        ok=mm.nUMI.notna(); x=x.loc[:,ok.values]; mm=mm[ok]
        cnt=np.round((np.expm1(x.values))*mm.nUMI.values/1e4)
        assert np.allclose(cnt.sum(0),mm.nUMI.values,rtol=1e-3)
        lab=f.split('/')[-1].split('_')[1]
        typ=mm.cell_type_age.str.rsplit('_',n=1).str[0] if 'cell_type_age' in mm else mm.cluster
        typ=mm['cluster'].astype(str)
        P=pd.DataFrame(cnt,index=x.index,columns=x.columns).T.groupby(typ.values).sum()
        P.index=[f'{lab}|{t}' for t in P.index]; Ps.append(P)
        Ms.append(pd.DataFrame({'sample':lab,'celltype':P.index.str.split('|').str[1],'n_cells':typ.value_counts().reindex(P.index.str.split('|').str[1]).values,'age':'old' if lab.startswith('O') else 'young'},index=P.index))
        print(lab,x.shape,flush=True)
    P=pd.concat(Ps).fillna(0); M=pd.concat(Ms)
    P.to_parquet(R+'B4_GSE129788/pb.parquet'); M.to_csv(R+'B4_GSE129788/pb_meta.csv'); print('B4',P.shape,M.groupby(['celltype','age']).size().unstack().to_string())
