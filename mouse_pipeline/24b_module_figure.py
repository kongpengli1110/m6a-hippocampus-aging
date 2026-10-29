"""Figure for modules: enrichment (-log10 q) heatmap for all modules x m6A sets + donor-level age beta; hub-gene lists for m6A-enriched modules."""
import pandas as pd,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt,seaborn as sns
T='./analysis/tables/'; F='./analysis/figures/'
en=pd.read_csv(T+'mod_enrichment.csv'); ag=pd.read_csv(T+'mod_age_effect_donor.csv'); mm=pd.read_csv(T+'mod_gene_membership_kME.csv')
sets=['A1_HIP_hypo_union','A1_DG_hyper','A3_ctx_Glut_hyper','METTL3cKO_HIP_down','METTL3cKO_HIP_up','METTL3KO_Cb_m6Aloss','YTHDF1_CLIP_HIP']
Q=en.pivot(index='module',columns='geneset',values='q')[sets]; Fo=en.pivot(index='module',columns='geneset',values='fold')[sets]
mods=ag.sort_values(['celltype','module']).module.tolist()
H=(-np.log10(Q.loc[mods])).clip(0,15)*np.sign(np.log(Fo.loc[mods]+1e-9))
fig,ax=plt.subplots(1,2,figsize=(13,max(8,len(mods)*0.25)),gridspec_kw={'width_ratios':[3,1]})
an=Q.loc[mods].map(lambda q:'***' if q<0.001 else '**' if q<0.01 else '*' if q<0.05 else '')
sns.heatmap(H,cmap='RdBu_r',center=0,annot=an,fmt='',ax=ax[0],cbar_kws={'label':'signed -log10 q (hypergeometric)'})
lab=[f"{m} (n={int(ag.set_index('module').loc[m,'size'])})" for m in mods]; ax[0].set_yticklabels(lab,fontsize=7); ax[0].set_title('Co-expression modules in Allen HPF excitatory neurons (metacells)')
a=ag.set_index('module').loc[mods]; y=np.arange(len(mods))+.5
ax[1].barh(y,a.beta_aged,color=['#d1495b' if p<0.05 else '#bbbbbb' for p in a.p]); ax[1].set_ylim(len(mods),0); ax[1].set_yticks([]); ax[1].axvline(0,color='k',lw=.5)
ax[1].set_xlabel('ME aged-vs-adult (SD; donor-level OLS)'); ax[1].set_title('Module age effect\n(red: nominal p<0.05)')
plt.tight_layout(); plt.savefig(F+'P2b_Fig4_excitatory_modules.png',dpi=300)
hub=[]
for m in en[(en.q<0.01)&en.geneset.isin(['A1_HIP_hypo_union','METTL3cKO_HIP_down','METTL3KO_Cb_m6Aloss'])].module.unique():
    g=mm[mm.module==m].sort_values('kME_own',ascending=False)
    hub.append(dict(module=m,size=len(g),enriched_for=';'.join(en[(en.module==m)&(en.q<0.05)].geneset),age_beta=ag.set_index('module').loc[m,'beta_aged'],age_p=ag.set_index('module').loc[m,'p'],top15_hubs=';'.join(g.gene.head(15))))
pd.DataFrame(hub).to_csv(T+'mod_m6A_enriched_hubs.csv',index=False); print(pd.DataFrame(hub).to_string())
