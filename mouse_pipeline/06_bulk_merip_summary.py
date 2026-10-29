"""2a. Bulk MeRIP aging diff-methylated genes per hippocampal subregion (authors' processed outputs only).
A1 GSE198526 (PNAS 2023, 3 vs 16 mo, CA1/CA3/DG/ACC; Dataset S17 FC>=1.2 & padj<=0.05) ; A2 GSE217882 (Aging Cell 2023, whole HIP 3 vs 20 mo, FDR<0.05)."""
import pandas as pd,numpy as np,itertools,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, seaborn as sns
from scipy.stats import hypergeom
D='./'; T=D+'analysis/tables/'; F=D+'analysis/figures/'
gs=pd.read_csv(T+'genesets_m6a_aging.csv'); bg=pd.read_csv(T+'merip_backgrounds.csv')
G={k:set(v) for k,v in gs.groupby('set').gene}; B={k:set(v) for k,v in bg.groupby('region').gene}
# --- counts
rows=[]
for r in ['CA1','CA3','DG','ACC']:
    rows.append(dict(region=r,hypo=len(G[f'A1_{r}_hypo']),hyper=len(G[f'A1_{r}_hyper']),mixed=len(G[f'A1_{r}_mixed']),methylated_background=len(B[f'A1_{r}'])))
rows.append(dict(region='wholeHIP (A2)',hypo=len(G['A2_HIP_hypo']),hyper=len(G['A2_HIP_hyper']),mixed=np.nan,methylated_background=np.nan))
C=pd.DataFrame(rows); C['pct_hypo']=C.hypo/(C.hypo+C.hyper+C.mixed.fillna(0))*100; C.to_csv(T+'2a_merip_counts.csv',index=False); print(C)
# GEO diff-peak BED (ExomePeak) peak-level summary
pk=[]
for r in ['CA1','CA3','DG','ACC']:
    d=pd.read_csv(D+f'data/A1_GSE198526/GSE198526_{r}_sig_diff_peak.bed.gz',sep='\t'); pk.append(dict(region=r,diff_peaks_BED=len(d),peaks_p_le_005=(d.iloc[:,4]<=0.05).sum(),genes=d.iloc[:,3].nunique()))
pd.DataFrame(pk).to_csv(T+'2a_GEO_diffpeak_bed_summary.csv',index=False); print(pd.DataFrame(pk))
# --- overlaps among hypo sets (hippocampal subregions + ACC) and with A2 whole HIP
sets=['A1_CA1_hypo','A1_CA3_hypo','A1_DG_hypo','A1_ACC_hypo','A1_CA1_hyper','A1_CA3_hyper','A1_DG_hyper','A2_HIP_hypo','A2_HIP_hyper']
U=set.union(*B.values())  # background: genes with MeRIP peak in any A1 region
ov=[]
for a,b in itertools.combinations(sets,2):
    A=G[a]&U; Bb=G[b]&U; k=len(A&Bb)
    p=hypergeom.sf(k-1,len(U),len(A),len(Bb)); ov.append(dict(set1=a,set2=b,n1=len(A),n2=len(Bb),overlap=k,expected=len(A)*len(Bb)/len(U),OR_like=k/(len(A)*len(Bb)/len(U)+1e-9),p=p,jaccard=k/len(A|Bb)))
OV=pd.DataFrame(ov); OV.to_csv(T+'2a_merip_set_overlaps.csv',index=False)
core=G['A1_CA1_hypo']&G['A1_CA3_hypo']&G['A1_DG_hypo']; pd.Series(sorted(core)).to_csv(T+'2a_core_hypo_CA1_CA3_DG.csv',index=False,header=['gene'])
print('core hypo in all 3 HIP subregions:',len(core))
# --- figure
fig,ax=plt.subplots(1,3,figsize=(15,4.6),gridspec_kw={'width_ratios':[1.1,1.3,1]})
c=C.set_index('region')[['hypo','hyper','mixed']].fillna(0); c.plot.bar(stacked=True,ax=ax[0],color=['#3b75af','#d1495b','#aaaaaa'])
ax[0].set_ylabel('Differentially methylated genes'); ax[0].set_title('Aging m6A changes (authors\' calls)'); ax[0].tick_params(axis='x',rotation=30)
for i,(n,v) in enumerate(C.pct_hypo.items()): ax[0].text(i,c.sum(1).iloc[i]+30,f'{v:.0f}% hypo',ha='center',fontsize=8)
lab=['A1_CA1_hypo','A1_CA3_hypo','A1_DG_hypo','A1_ACC_hypo','A1_DG_hyper','A2_HIP_hypo','A2_HIP_hyper']
M=pd.DataFrame(np.nan,index=lab,columns=lab)
for _,r in OV.iterrows():
    if r.set1 in lab and r.set2 in lab: M.loc[r.set1,r.set2]=M.loc[r.set2,r.set1]=np.log2(r.OR_like+1e-3)
P=pd.DataFrame('',index=lab,columns=lab)
for _,r in OV.iterrows():
    if r.set1 in lab and r.set2 in lab: P.loc[r.set1,r.set2]=P.loc[r.set2,r.set1]=f"{int(r.overlap)}\n{'***' if r.p<1e-3 else '**' if r.p<1e-2 else '*' if r.p<0.05 else ''}"
sns.heatmap(M,annot=P,fmt='',cmap='RdBu_r',center=0,ax=ax[1],cbar_kws={'label':'log2(obs/exp)'},annot_kws={'fontsize':7})
ax[1].set_title('Overlap (background: MeRIP-methylated genes)')
# A1 vs A2 direction concordance
a2h=G['A2_HIP_hyper']; a2o=G['A2_HIP_hypo']
tab=[]
for r in ['CA1','CA3','DG']:
    tab.append([len(G[f'A1_{r}_hypo']&a2o),len(G[f'A1_{r}_hypo']&a2h),len(G[f'A1_{r}_hyper']&a2o),len(G[f'A1_{r}_hyper']&a2h)])
TB=pd.DataFrame(tab,index=['CA1','CA3','DG'],columns=['A1hypo&A2hypo','A1hypo&A2hyper','A1hyper&A2hypo','A1hyper&A2hyper']); TB.to_csv(T+'2a_A1_vs_A2_direction.csv')
TB.plot.bar(ax=ax[2],color=['#3b75af','#9b59b6','#e67e22','#d1495b']); ax[2].set_title('A1 subregion vs A2 whole-HIP direction'); ax[2].tick_params(axis='x',rotation=0); ax[2].legend(fontsize=7)
plt.tight_layout(); plt.savefig(F+'Fig2a_bulk_MeRIP_aging.png',dpi=300); print(TB)
