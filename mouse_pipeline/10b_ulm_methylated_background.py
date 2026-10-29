"""Sensitivity for 2d(i): restrict the competitive ULM universe to genes methylated in A1 hippocampal MeRIP (any region)
so that 'aging diff-methylated' sets are compared against other methylated genes, not against all genes."""
import pandas as pd,numpy as np,decoupler as dc,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt,seaborn as sns
from statsmodels.stats.multitest import multipletests
D='./'; T=D+'analysis/tables/'; F=D+'analysis/figures/'
R=pd.read_csv(T+'allen_HPF_pseudobulk_DESeq2_age.csv')
types=['CA1','CA3','DG','SUB','ENT_Glut','GABA','Astro','Oligo','OPC','Microglia','Endo','Mural','VLMC','CR','BAM']
bg=pd.read_csv(T+'merip_backgrounds.csv'); meth=set(bg[bg.region.isin(['A1_CA1','A1_CA3','A1_DG'])].gene)
St=R.dropna(subset=['stat']).pivot_table(index='celltype',columns='gene',values='stat',aggfunc='mean').loc[types].fillna(0)
Sm=St[[g for g in St.columns if g in meth]]
gsA=pd.read_csv(T+'genesets_m6a_aging.csv'); gsR=pd.read_csv(T+'genesets_regulator_targets.csv')
keepA=['A1_CA1_hypo','A1_CA3_hypo','A1_DG_hypo','A1_CA3_hyper','A1_DG_hyper','A2_HIP_hyper','A2_HIP_hypo','A3_ctx_Glut_hyper','A3_ctx_Glut_hypo','A3_ctx_Astro_hyper','A3_ctx_Oligo_hyper','A3_ctx_Endo_hyper','A3_ctx_GABA_hyper']
net=pd.concat([gsA[gsA.set.isin(keepA)][['set','gene']],gsR[['set','gene']]]).rename(columns={'set':'source','gene':'target'}).drop_duplicates(); net['weight']=1.0
out={}
for lab,M in [('allgenes',St),('methylated_only',Sm)]:
    es,pv=dc.mt.ulm(data=M,net=net,tmin=10); q=pd.DataFrame(multipletests(pv.values.ravel(),method='fdr_bh')[1].reshape(pv.shape),index=pv.index,columns=pv.columns)
    es.to_csv(T+f'2d_ULM_{lab}_score.csv'); q.to_csv(T+f'2d_ULM_{lab}_q.csv'); out[lab]=(es,q)
    print(lab,M.shape[1],'genes'); print(es[[c for c in keepA if c in es]].round(1).to_string())
# mean DE stat of methylated vs non-methylated genes per type (global shift)
g=pd.DataFrame({t:[St.loc[t,[c for c in St.columns if c in meth]].mean(),St.loc[t,[c for c in St.columns if c not in meth]].mean()] for t in types},index=['methylated','not_methylated']).T
g.to_csv(T+'2d_mean_DEstat_methylated_vs_not.csv'); print(g.round(2))
es,q=out['methylated_only']; order=[s for s in keepA if s in es.columns]+sorted([s for s in es.columns if s not in keepA])
fig,ax=plt.subplots(figsize=(15,7)); an=q[order].map(lambda v:'***' if v<0.001 else '**' if v<0.01 else '*' if v<0.05 else '')
sns.heatmap(es[order].clip(-8,8),cmap='RdBu_r',center=0,annot=an,fmt='',ax=ax,cbar_kws={'label':'ULM t (universe = MeRIP-methylated genes)'})
ax.set_title('Age-related shift of m6A programs per HPF cell type, competitive within methylated genes (BH q * <0.05)'); plt.xticks(fontsize=7); ax.set_ylabel('')
plt.tight_layout(); plt.savefig(F+'Fig2d_target_program_activity_methylatedUniverse.png',dpi=300)
