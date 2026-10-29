"""Allen Zeng-Aging-Mouse-10Xv3, HPF (HIP+SUB/PAR+ENT ROIs), unbiased-sampled libraries only.
Pseudobulk = sum of raw UMI per donor x major cell type (one region per donor). DESeq2 (pydeseq2) ~ region + sex + age per cell type.
Also: EWCE-style specificity matrix from mean CPM across donors."""
import sys; sys.path.insert(0,'./scripts')
from common import *
import warnings; warnings.filterwarnings('ignore')
import pydeseq2_patch
from pydeseq2.dds import DeseqDataSet; from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
X,M,genes=load_pb()
sym=dict(zip(genes.gene_identifier,genes.gene_symbol))
M.to_csv(D+'analysis/tables/allen_HPF_pseudobulk_samples.csv')
print(M.groupby(['major','age']).size().unstack(fill_value=0))
# --- specificity (EWCE): mean CPM per type over donors with >=20 cells
cpm=X.div(X.sum(1),axis=0)*1e6
ok=M.n_cells>=20
types=[t for t in M.major.unique() if t not in ('Other','OtherGlut')]
mean=pd.DataFrame({t:cpm[(M.major==t)&ok].mean() for t in types if ((M.major==t)&ok).sum()>=3})
mean.index=mean.index.map(sym); mean=mean.groupby(level=0).sum()
spec=mean.div(mean.sum(1),axis=0).fillna(0)
mean.to_csv(D+'analysis/tables/allen_HPF_meanCPM_by_type.csv'); spec.to_csv(D+'analysis/tables/allen_HPF_specificity_EWCE.csv')
# --- DE per cell type
inf=DefaultInference(n_cpus=1)
res=[]
for t in ['CA1','CA3','DG','SUB','GABA','Astro','Oligo','OPC','Microglia','Endo','Mural','VLMC','ENT_Glut','CR','BAM']:
    sel=(M.major==t)&(M.n_cells>=20)
    m=M[sel]; 
    if m.groupby('age').size().min()<2: print('skip',t,m.groupby('age').size().to_dict()); continue
    c=X[sel.values]; c=c.loc[:,(c>=10).sum()>=max(3,int(0.5*len(c)))]
    des='~sex+age' if m.region.nunique()==1 else '~region+sex+age'
    if m.region.nunique()>1 and (m.groupby('region').age.nunique()<2).any():
        # keep design estimable
        des='~sex+age'
    meta=m[['region','sex','age']].copy()
    dds=DeseqDataSet(counts=c.astype(int),metadata=meta,design=des,refit_cooks=True,inference=inf,quiet=True)
    dds.deseq2()
    st=DeseqStats(dds,contrast=['age','aged','adult'],inference=inf,quiet=True); st.summary()
    r=st.results_df.copy(); r['gene']=r.index.map(sym); r['celltype']=t; r['n_adult']=(m.age=='adult').sum(); r['n_aged']=(m.age=='aged').sum(); r['design']=des
    res.append(r); print(t,des,len(m),'genes',c.shape[1],'padj<0.05',(r.padj<0.05).sum(),flush=True)
R=pd.concat(res); R.index.name='ensembl'; R.to_csv(D+'analysis/tables/allen_HPF_pseudobulk_DESeq2_age.csv')
