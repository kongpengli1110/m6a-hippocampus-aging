"""B10 GSE212336 (Hahn 2023 Cell): hippocampus (anterior+posterior) bulk RNA-seq 3-28 mo; trajectories of m6A program scores and regulators."""
import pandas as pd,numpy as np,gzip,re,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import sys; sys.path.insert(0,'./scripts'); from common import REG2CLS
D='./'; T=D+'analysis/tables/'; F=D+'analysis/figures/'
L=[l.rstrip('\n').split('\t') for l in gzip.open(D+'data/B10_GSE212336/GSE212336_series_matrix.txt.gz','rt') if l.startswith(('!Sample_title','!Sample_geo_accession'))]
tit=[x.strip('"') for x in L[0][1:]]; desc=[l for l in gzip.open(D+'data/B10_GSE212336/GSE212336_series_matrix.txt.gz','rt') if l.startswith('!Sample_description')][2].rstrip('\n').split('\t')[1:]
meta=pd.DataFrame({'title':tit,'col':[d.strip('"') for d in desc]})
pt=meta.title.str.extract(r'^(.*) (Male|Female) (\d+)mo (\d+)ID$'); meta[['tissue','sex','age','mouse']]=pt; meta['age']=meta.age.astype(int)
meta=meta[meta.tissue.str.startswith('Hippocampus')].set_index('col')
cnt=pd.read_csv(D+'data/B10_GSE212336/GSE212336_BulkSeq_Counttable.txt.gz',sep='\t',index_col=0)
cnt=cnt[[c for c in meta.index if c in cnt.columns]]; meta=meta.loc[cnt.columns]
print(meta.groupby(['tissue','sex','age']).size().unstack(fill_value=0))
cpm=cnt/cnt.sum()*1e6; lg=np.log2(cpm[(cpm>=1).mean(axis=1)>=0.5]+1)
Z=lg.sub(lg.mean(axis=1),axis=0).div(lg.std(axis=1),axis=0)
sets=pd.read_csv(T+'2e_score_genesets.csv'); SS={k:[g for g in v if g in Z.index] for k,v in sets.groupby('set').gene}
S=pd.DataFrame({k:Z.loc[v].mean() for k,v in SS.items() if len(v)>=10}); S=S.join(meta)
for g in REG2CLS:
    if g in lg.index: S['gene_'+g]=lg.loc[g]
S.to_csv(T+'B10_GSE212336_HIP_scores.csv')
rows=[]
# one value per mouse (mean of anterior/posterior) to avoid pseudo-replication
sc=[c for c in S.columns if c.startswith(('MK_','M6A','REG_','gene_'))]
m=S[S.sex=='Male'].groupby('mouse').agg({**{c:'mean' for c in sc},'age':'first'}); print('male mice',len(m), m.age.value_counts().sort_index().to_dict())
for c in sc:
    d=m.assign(y=(m[c]-m[c].mean())/m[c].std())
    f=smf.ols('y~age',d).fit(); f2=smf.ols('y~age+I(age**2)',d).fit()
    rows.append(dict(score=c,n=len(d),beta_SD_per_month=f.params['age'],p_linear=f.pvalues['age'],p_quadratic=f2.pvalues['I(age ** 2)']))
R=pd.DataFrame(rows); R['q_linear']=multipletests(R.p_linear,method='fdr_bh')[1]; R.to_csv(T+'B10_trajectory_models_male.csv',index=False)
print(R.sort_values('p_linear').head(25).to_string())
show=['M6A_CA1','M6A_DG','M6A_core87','MK_CA1','MK_DG','MK_Microglia','MK_Oligo','REG_METTL3cKO_HIP_down','gene_Mettl3','gene_Fto','gene_Ythdf1','gene_Alkbh5']
fig,axs=plt.subplots(3,4,figsize=(14,9),sharex=True); 
for ax,c in zip(axs.ravel(),show):
    if c not in S: continue
    for (tis,sx),d in S.groupby(['tissue','sex']):
        g=d.groupby('age')[c].agg(['mean','sem']); ax.errorbar(g.index,g['mean'],yerr=g['sem'],marker='o',ms=3,label=f"{tis.split('(')[1][:-1]} {sx[0]}",ls='-' if sx=='Male' else ':')
    r=R.set_index('score').loc[c]; ax.set_title(f"{c}\nmale q={r.q_linear:.2g}",fontsize=9)
axs[0,0].legend(fontsize=6); [a.set_xlabel('age (mo)') for a in axs[-1]]
plt.suptitle('GSE212336 hippocampus bulk RNA-seq: trajectories (z-score means)'); plt.tight_layout(); plt.savefig(F+'FigS_B10_trajectories.png',dpi=300)
