"""Figures/summary for replication: forest plots of length-adjusted gene-set age effect per cell type x dataset; EWCE z by dataset; regulator consistency."""
import pandas as pd,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt,seaborn as sns
T='./analysis/tables/'; F='./analysis/figures/'
SE=pd.read_csv(T+'rep_geneset_age_effect_lenadj.csv'); EW=pd.read_csv(T+'rep_EWCE_by_dataset.csv'); RE=pd.read_csv(T+'rep_regulator_DE_by_dataset.csv')
DSO=['Allen','GSE161340_sn','GSE161340_sc','GSE233363_DG','GSE129788']; COL=dict(zip(DSO,['k','#1f77b4','#9ecae1','#2ca02c','#ff7f0e']))
CT=['CA1','CA3','DG','SUB','Neuron(all)','GABA','Astro','Oligo','OPC','Microglia','Endo','Mural']
sets=['A1_HIP_hypo_union','METTL3cKO_HIP_down','METTL3KO_Cb_m6Aloss','YTHDF1_CLIP_HIP']
fig,axs=plt.subplots(1,len(sets),figsize=(5*len(sets),8),sharey=True)
for ax,s in zip(axs,sets):
    d=SE[SE.geneset==s]; y=0; ticks=[]
    for ct in CT:
        dd=d[d.celltype==ct]
        if dd.empty: continue
        for j,ds in enumerate(DSO):
            r=dd[dd.dataset==ds]
            if r.empty: continue
            r=r.iloc[0]; ax.errorbar(r.coef_lenadj,y-j*0.13,xerr=1.96*r.se_lenadj,fmt='o',ms=3,color=COL[ds],label=ds if (ct==CT[0] or y==0) else None)
        ticks.append((y-0.25,ct)); y-=1
    ax.axvline(0,color='grey',ls='--'); ax.set_yticks([t[0] for t in ticks]); ax.set_yticklabels([t[1] for t in ticks]); ax.set_title(s+'\n(length/expression-adjusted, methylated universe)',fontsize=9)
    ax.set_xlabel('set coefficient on aged-vs-young log2FC')
from matplotlib.lines import Line2D
axs[-1].legend([Line2D([0],[0],marker='o',color=COL[d],ls='') for d in DSO],[d+(' (8v15 donors)' if d=='Allen' else ' (2v2 pools)' if 'GSE161340' in d else ' (1v2 libs, descriptive)' if 'GSE233363' in d else ' (8v8 mice)') for d in DSO],fontsize=7,loc='lower right')
plt.tight_layout(); plt.savefig(F+'P2b_Fig1_replication_forest_programs.png',dpi=300)
# EWCE heat
E=EW.pivot_table(index=['geneset','dataset'],columns='celltype',values='z'); cols=[c for c in CT+['VLMC','BAM','DG_IMN','ENT_Glut','CR'] if c in E.columns]
E=E[cols]; Q=EW.pivot_table(index=['geneset','dataset'],columns='celltype',values='q')[cols]
fig,ax=plt.subplots(figsize=(13,9)); sns.heatmap(E.clip(-10,10),cmap='RdBu_r',center=0,annot=Q.map(lambda q:'' if pd.isna(q) else '*' if q<0.05 else ''),fmt='',ax=ax,cbar_kws={'label':'EWCE z (own young-sample reference)'})
ax.set_title('Cell-type attribution of A1 aging-DM genes using each dataset\'s own cell-type specificity (* BH q<0.05)'); plt.tight_layout(); plt.savefig(F+'P2b_Fig1b_replication_EWCE.png',dpi=300)
# consistency summary: sign agreement of lenadj coef with Allen per set x celltype
rows=[]
for s in SE.geneset.unique():
    for ct in SE.celltype.unique():
        d=SE[(SE.geneset==s)&(SE.celltype==ct)].set_index('dataset')
        if 'Allen' not in d.index or len(d)<2: continue
        a=np.sign(d.loc['Allen','coef_lenadj']); o=d.drop('Allen')
        rows.append(dict(geneset=s,celltype=ct,allen_coef=d.loc['Allen','coef_lenadj'],allen_q=d.loc['Allen','q_lenadj'],n_rep=len(o),n_same_sign=int((np.sign(o.coef_lenadj)==a).sum()),n_same_sign_p05=int(((np.sign(o.coef_lenadj)==a)&(o.p_lenadj<0.05)).sum()),n_opposite_p05=int(((np.sign(o.coef_lenadj)!=a)&(o.p_lenadj<0.05)).sum()),
                         inv_var_meta=(d.coef_lenadj/d.se_lenadj**2).sum()/(1/d.se_lenadj**2).sum(),meta_se=np.sqrt(1/(1/d.se_lenadj**2).sum())))
C=pd.DataFrame(rows); C['meta_z']=C.inv_var_meta/C.meta_se; C.to_csv(T+'rep_consistency_summary.csv',index=False)
print(C[C.geneset.isin(sets)].round(3).to_string())
# regulators: concordance with Allen
pv=RE.pivot_table(index=['celltype','gene'],columns='dataset',values='log2FC')
print(pv.corr(min_periods=20).round(2))
