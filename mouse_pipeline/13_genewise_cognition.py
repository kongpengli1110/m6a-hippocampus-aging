"""Gene-wise within-old (24+28 mo) association with Barnes MedianRankCombinedScore in GSE302188; adj. age group. Global signal check."""
import pandas as pd,numpy as np,numpy.linalg as la
from scipy import stats; from statsmodels.stats.multitest import multipletests
D='./'; T=D+'analysis/tables/'
cnt=pd.read_csv(D+'data/B11_GSE302188/GSE302188_brain_RNA_gene_counts.csv.gz',index_col=0)
e2s=pd.read_csv(D+'data/B11_GSE302188/zenodo/GRCm38.102.EnsemblToGeneSymbol.tsv.gz',sep='\t',header=None,index_col=0)[1]
cp=cnt/cnt.sum()*1e6; L=np.log2(cp[(cp>=1).mean(axis=1)>=0.5]+1)
S=pd.read_csv(T+'2e_GSE302188_scores.csv',index_col=0); o=S[S.age.isin(['24M','28M'])].dropna(subset=['cog'])
X=L[o.index]
A=np.column_stack([np.ones(len(o)),(o.age=='28M').astype(float),(o.cog-o.cog.mean())/o.cog.std()])
inv=la.inv(A.T@A); B=X.values@A@inv; res=X.values-B@A.T; df=len(o)-3; s2=(res**2).sum(1)/df; t=B[:,2]/np.sqrt(s2*inv[2,2])
p=2*stats.t.sf(abs(t),df); q=multipletests(p,method='fdr_bh')[1]
R=pd.DataFrame({'ensembl':X.index,'gene':X.index.map(e2s),'t':t,'p':p,'q':q}).sort_values('p'); R.to_csv(T+'2e_GSE302188_genewise_cognition_withinOld.csv',index=False)
print('genes',len(p),'p<0.001',(p<0.001).sum(),'(expected',round(len(p)*0.001),') q<0.1',(q<0.1).sum(),'pi0-ish frac p<0.05',(p<0.05).mean().round(3))
print(R.head(15).to_string())
