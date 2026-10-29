"""2b. Cell-type attribution of bulk aging diff-methylated genes.
EWCE-like test: mean Allen-HPF specificity of a gene set per cell type vs 10,000 random sets drawn from the
proper background (genes methylated in that MeRIP region AND expressed in Allen HPF), matched on overall HPF expression decile.
Plus overlap with Tegowski cortex cell-type aging DM-RNAs (orthogonal, conditional; background = methylated in both)."""
import pandas as pd,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, seaborn as sns
from scipy.stats import hypergeom; from statsmodels.stats.multitest import multipletests
D='./'; T=D+'analysis/tables/'; F=D+'analysis/figures/'
rng=np.random.default_rng(42); NPERM=10000
spec=pd.read_csv(T+'allen_HPF_specificity_EWCE.csv',index_col=0); mean=pd.read_csv(T+'allen_HPF_meanCPM_by_type.csv',index_col=0)
types=['CA1','CA3','DG','SUB','ENT_Glut','GABA','Astro','Oligo','OPC','Microglia','Endo','Mural','VLMC','BAM']
spec=spec[types]; tot=mean[types].sum(1)
expressed=set(tot[tot>=1].index)   # >=1 CPM summed over types
gs=pd.read_csv(T+'genesets_m6a_aging.csv'); G={k:set(v) for k,v in gs.groupby('set').gene}
bg=pd.read_csv(T+'merip_backgrounds.csv'); B={k:set(v) for k,v in bg.groupby('region').gene}
B['A2_HIP']=set.union(*[B[f'A1_{r}'] for r in ['CA1','CA3','DG']])  # A2 provides no peak-level background; use union of A1 HIP methylated genes (approximation, stated)
def ewce(genes,background):
    bgl=sorted((background&expressed)&set(spec.index)); g=sorted(set(genes)&set(bgl))
    dec=pd.qcut(np.log1p(tot[bgl]),10,labels=False,duplicates='drop'); dec=pd.Series(dec,index=bgl)
    obs=spec.loc[g].mean().values
    # matched sampling
    pools={d:np.array(dec.index[dec==d]) for d in dec.unique()}; need=dec[g].value_counts()
    S=spec.loc[bgl].values; idx={x:i for i,x in enumerate(bgl)}
    null=np.zeros((NPERM,len(types)))
    pidx={d:np.array([idx[x] for x in pools[d]]) for d in pools}
    for i in range(NPERM):
        sel=np.concatenate([rng.choice(pidx[d],n,replace=False) for d,n in need.items()])
        null[i]=S[sel].mean(0)
    p=(1+(null>=obs).sum(0))/(NPERM+1); z=(obs-null.mean(0))/null.std(0); pl=(1+(null<=obs).sum(0))/(NPERM+1)
    return pd.DataFrame({'celltype':types,'n_genes':len(g),'n_background':len(bgl),'obs_mean_spec':obs,'null_mean':null.mean(0),'z':z,'p_enrich':p,'p_deplete':pl})
res=[]
tests=[(f'A1_{r}_{d}',f'A1_{r}') for r in ['CA1','CA3','DG','ACC'] for d in ['hypo','hyper']]+[('A2_HIP_hyper','A2_HIP'),('A2_HIP_hypo','A2_HIP')]
for s,b in tests:
    r=ewce(G[s],B[b]); r.insert(0,'geneset',s); res.append(r); print(s,len(G[s]),r.loc[r.z.idxmax(),['celltype','z','p_enrich']].tolist(),flush=True)
R=pd.concat(res); R['q_enrich']=multipletests(R.p_enrich,method='fdr_bh')[1]; R['q_deplete']=multipletests(R.p_deplete,method='fdr_bh')[1]
R.to_csv(T+'2b_EWCE_celltype_enrichment.csv',index=False)
# sensitivity: unmatched background = all Allen-expressed genes (shows the bias the proper background removes)
res2=[]
for s,b in tests[:6]:
    r=ewce(G[s],expressed); r.insert(0,'geneset',s); res2.append(r)
pd.concat(res2).to_csv(T+'2b_EWCE_sensitivity_allExpressedBackground.csv',index=False)
# heatmap
Z=R.pivot(index='geneset',columns='celltype',values='z')[types].loc[[t[0] for t in tests]]
Q=R.pivot(index='geneset',columns='celltype',values='q_enrich')[types].loc[Z.index]
ann=Q.map(lambda q:'***' if q<0.001 else '**' if q<0.01 else '*' if q<0.05 else '')
fig,ax=plt.subplots(figsize=(11,5.2)); sns.heatmap(Z,cmap='RdBu_r',center=0,annot=ann,fmt='',ax=ax,cbar_kws={'label':'EWCE z (vs expression-matched methylated background)'})
ax.set_title('Cell-type specificity of aging differentially m6A-methylated genes (Allen HPF 2 vs 18 mo reference)\n* BH q<0.05 (enrichment), 10,000 permutations'); ax.set_xlabel(''); ax.set_ylabel('')
plt.tight_layout(); plt.savefig(F+'Fig2b_EWCE_celltype_attribution.png',dpi=300)
# --- orthogonal: overlap with Tegowski cortex cell-type DM-RNAs
dart_bg=set(pd.read_csv(T+'A3_scDART_cortex_methylated_RNAs_background.csv').gene)
a3=[k for k in G if k.startswith('A3_ctx_') and 'Sites' not in k and len(G[k])>=15]
ov=[]
for s,b in tests:
    U=B[b]&dart_bg
    for t in a3:
        A=G[s]&U; Bb=G[t]&U; k=len(A&Bb); p=hypergeom.sf(k-1,len(U),len(A),len(Bb))
        ov.append(dict(bulk_set=s,scDART_set=t,background=len(U),n_bulk=len(A),n_sc=len(Bb),overlap=k,expected=len(A)*len(Bb)/max(len(U),1),fold=k/max(len(A)*len(Bb)/max(len(U),1),1e-9),p=p,genes=';'.join(sorted(A&Bb))))
OV=pd.DataFrame(ov); OV['q']=multipletests(OV.p,method='fdr_bh')[1]; OV.to_csv(T+'2b_overlap_bulkMeRIP_vs_scDART_celltype.csv',index=False)
print(OV.sort_values('p').head(12)[['bulk_set','scDART_set','background','n_bulk','n_sc','overlap','expected','fold','p','q']].to_string())
FO=OV.pivot(index='bulk_set',columns='scDART_set',values='fold').loc[[t[0] for t in tests]]
QO=OV.pivot(index='bulk_set',columns='scDART_set',values='q').loc[FO.index]; KO=OV.pivot(index='bulk_set',columns='scDART_set',values='overlap').loc[FO.index]
an=KO.astype(int).astype(str)+QO.map(lambda q:'\n***' if q<0.001 else '\n**' if q<0.01 else '\n*' if q<0.05 else '')
fig,ax=plt.subplots(figsize=(10,5.5)); sns.heatmap(np.log2(FO+0.05),cmap='RdBu_r',center=0,annot=an,fmt='',ax=ax,annot_kws={'fontsize':7},cbar_kws={'label':'log2 fold over expected'})
ax.set_title('Overlap: bulk hippocampal aging DM genes vs Tegowski 2024 cortex scDART cell-type aging DM-RNAs\n(conditional orthogonal support; n=2 bio reps/age, age confounded with kit)'); ax.set_xlabel(''); ax.set_ylabel('')
plt.tight_layout(); plt.savefig(F+'Fig2b2_overlap_scDART.png',dpi=300)
