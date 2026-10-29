"""(also run with LEN=1: + gene-length tertile matching, GENCODE vM25 span)
Sensitivity for 2b: EWCE null matched jointly on HPF expression (5 bins) x number of MeTPeak peaks per gene
in the same region (young+old union, 1/2/3-4/5+), to control for detection power (long, multi-peak neuronal transcripts)."""
import pandas as pd,numpy as np
from statsmodels.stats.multitest import multipletests
import sys; LEN=len(sys.argv)>1
gl=pd.read_csv('./data/ref/mm10_gencodeM25_genes.tsv',sep='\t',header=None); GL=(gl[4]-gl[3]).groupby(gl[1]).max()
D='./'; T=D+'analysis/tables/'
rng=np.random.default_rng(7); NPERM=10000
spec=pd.read_csv(T+'allen_HPF_specificity_EWCE.csv',index_col=0); mean=pd.read_csv(T+'allen_HPF_meanCPM_by_type.csv',index_col=0)
types=['CA1','CA3','DG','SUB','ENT_Glut','GABA','Astro','Oligo','OPC','Microglia','Endo','Mural','VLMC','BAM']
spec=spec[types]; tot=mean[types].sum(axis=1)
genes=pd.read_csv(D+'data/allen/genes.csv'); e2s=dict(zip(genes.gene_identifier,genes.gene_symbol))
gs=pd.read_csv(T+'genesets_m6a_aging.csv'); G={k:set(v) for k,v in gs.groupby('set').gene}
out=[]
for r in ['CA1','CA3','DG','ACC']:
    d=pd.concat([pd.read_csv(D+f'data/A1_GSE198526/GSE198526_{r}{t}_peak.bed.gz',sep='\t') for t in ['young','old']])
    npk=d.iloc[:,3].map(e2s).value_counts()
    bgl=[g for g in npk.index if g in spec.index and tot.get(g,0)>=1 and g in GL.index]
    eb=pd.qcut(np.log1p(tot[bgl]),5,labels=False); pb=pd.cut(npk[bgl],[0,1,2,4,1e9],labels=False)
    strat=pd.Series(eb.values*10+pb.values,index=bgl)
    if LEN: strat=strat*10+pd.qcut(np.log10(GL[bgl]),3,labels=False).values
    S=spec.loc[bgl].values; idx={x:i for i,x in enumerate(bgl)}
    pools={s:np.array([idx[x] for x in strat.index[strat==s]]) for s in strat.unique()}
    for dname in ['hypo','hyper']:
        g=[x for x in G[f'A1_{r}_{dname}'] if x in idx]
        need=strat[g].value_counts(); obs=S[[idx[x] for x in g]].mean(0)
        null=np.array([S[np.concatenate([rng.choice(pools[s],min(n,len(pools[s])),replace=n>len(pools[s])) for s,n in need.items()])].mean(0) for _ in range(NPERM)])
        z=(obs-null.mean(0))/null.std(0); p=(1+(null>=obs).sum(0))/(NPERM+1)
        out.append(pd.DataFrame({'geneset':f'A1_{r}_{dname}','celltype':types,'n_genes':len(g),'n_background':len(bgl),'z':z,'p_enrich':p}))
        print(r,dname,len(g),dict(zip(types,np.round(z,1))),flush=True)
R=pd.concat(out); R['q_enrich']=multipletests(R.p_enrich,method='fdr_bh')[1]; R.to_csv(T+('2b_EWCE_sensitivity_peak_length_matched.csv' if LEN else '2b_EWCE_sensitivity_peakcount_matched.csv'),index=False)
