"""2e. Cognition association of cell-type m6A program scores.
Discovery: GSE302188 brain (left hippocampus) bulk RNA-seq counts (GEO) x Zenodo 15849341 Barnes-maze MedianRankCombinedScore,
sample IDs matched by expression correlation (analysis/tables/B11_sample_matching.csv; 82/82 reciprocal best matches).
Replication: GSE13799 (Affymetrix 430 2.0, MAS5; aged 23 mo impaired n=6 vs unimpaired n=6; young n=6).
Pre-specified scores (no selection on outcome):
  MK_<type>   : composition proxy = top-100 Allen-HPF specific markers of that type
  M6A_<type>  : aging-hypomethylated genes (A1 union CA1/CA3/DG hypo) within top-decile specificity for that type
  M6Ahyper_<type>: same for hypermethylated (A1 union hyper + A2 hyper)
  REG_<set>   : regulator-perturbation target sets (scripts/08)
Score = mean row-z of log2CPM (or log2 MAS5) over genes present. Models (OLS):
  old-only (24+28 mo): score ~ C(agegroup) + cognition ; all ages: score ~ C(agegroup) + cognition ;
  composition-adjusted: M6A_t ~ C(agegroup) + cognition + MK_t """
import pandas as pd,numpy as np,statsmodels.formula.api as smf,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, seaborn as sns
from statsmodels.stats.multitest import multipletests; from scipy import stats
D='./'; T=D+'analysis/tables/'; F=D+'analysis/figures/'
spec=pd.read_csv(T+'allen_HPF_specificity_EWCE.csv',index_col=0); mean=pd.read_csv(T+'allen_HPF_meanCPM_by_type.csv',index_col=0)
types=['CA1','CA3','DG','SUB','GABA','Astro','Oligo','OPC','Microglia','Endo','Mural']
gs=pd.read_csv(T+'genesets_m6a_aging.csv'); G={k:set(v) for k,v in gs.groupby('set').gene}
rg=pd.read_csv(T+'genesets_regulator_targets.csv'); RG={k:set(v) for k,v in rg.groupby('set').gene}
expressed=mean.index[mean[types].max(axis=1)>=5]
hypo=G['A1_CA1_hypo']|G['A1_CA3_hypo']|G['A1_DG_hypo']; hyper=G['A1_CA1_hyper']|G['A1_CA3_hyper']|G['A1_DG_hyper']|G['A2_HIP_hyper']
SETS={}
for t in types:
    s=spec.loc[expressed,t]; mk=s[(mean.loc[expressed,t]>=10)].sort_values(ascending=False).head(100).index
    SETS[f'MK_{t}']=set(mk)
    top=set(s[s>=s.quantile(0.9)].index)
    SETS[f'M6A_{t}']=top&hypo; SETS[f'M6Ahyper_{t}']=top&hyper
for k,v in RG.items(): SETS[f'REG_{k}']=v
SETS['M6A_allHIPhypo']=hypo&set(expressed); SETS['M6A_core87']=G['A1_CA1_hypo']&G['A1_CA3_hypo']&G['A1_DG_hypo']
pd.DataFrame([(k,g) for k,v in SETS.items() for g in sorted(v)],columns=['set','gene']).to_csv(T+'2e_score_genesets.csv',index=False)
def score(E):  # E genes x samples (log scale)
    Z=E.sub(E.mean(axis=1),axis=0).div(E.std(axis=1)+1e-9,axis=0)
    out={};n={}
    for k,v in SETS.items():
        g=[x for x in v if x in Z.index]; n[k]=len(g)
        if len(g)>=10: out[k]=Z.loc[g].mean()
    return pd.DataFrame(out),pd.Series(n)
# ---------- GSE302188
cnt=pd.read_csv(D+'data/B11_GSE302188/GSE302188_brain_RNA_gene_counts.csv.gz',index_col=0)
e2s=pd.read_csv(D+'data/B11_GSE302188/zenodo/GRCm38.102.EnsemblToGeneSymbol.tsv.gz',sep='\t',header=None,index_col=0)[1]
cpm=cnt/cnt.sum()*1e6; keep=(cpm>=1).sum(axis=1)>=0.5*cnt.shape[1]; L=np.log2(cpm[keep]+1)
L.index=L.index.map(e2s); L=L[L.index.notna()]; L=L.groupby(level=0).mean()
mt=pd.read_csv(T+'B11_sample_matching.csv').set_index('geo_id')
cog=pd.read_csv(D+'data/B11_GSE302188/zenodo/df_allCognitionData_2025-02-17.tsv',sep='\t',index_col=0)
S,N=score(L); S['zen']=mt.loc[S.index,'zenodo_best']; S['age']=S.index.str.split('.').str[0]
S=S.join(cog,on='zen'); S['agegroup']=S.age.map({'3M':'m03','9M':'m09','15M':'m15','24M':'m24','28M':'m28'})
S['cog']=S.MedianRankCombinedScore
S.to_csv(T+'2e_GSE302188_scores.csv'); N.to_csv(T+'2e_geneset_sizes_GSE302188.csv')
print('n samples',len(S),'with cognition',S.cog.notna().sum(), S.groupby('age').cog.count().to_dict())
cols=[c for c in S.columns if c.startswith(('MK_','M6A','REG_'))]
rows=[]
for c in cols:
    for lab,d in [('old24_28',S[S.age.isin(['24M','28M'])]),('all_ages',S),('old_only28',S[S.age=='28M']),('old_only24',S[S.age=='24M'])]:
        d=d.dropna(subset=['cog']).copy(); d['y']=(d[c]-d[c].mean())/d[c].std(); d['x']=(d.cog-d.cog.mean())/d.cog.std()
        f='y~x' if d.agegroup.nunique()==1 else 'y~C(agegroup)+x'
        m=smf.ols(f,d).fit(); rows.append(dict(score=c,cohort=lab,n=len(d),beta_std=m.params['x'],se=m.bse['x'],p=m.pvalues['x'],model=f))
        if c.startswith('M6A_') and c.replace('M6A_','MK_') in S:
            mk=c.replace('M6A_','MK_'); d['mk']=(d[mk]-d[mk].mean())/d[mk].std(); m2=smf.ols(f+'+mk',d).fit()
            rows.append(dict(score=c+'|adjMK',cohort=lab,n=len(d),beta_std=m2.params['x'],se=m2.bse['x'],p=m2.pvalues['x'],model=f+'+MK'))
    # age effect (all ages, continuous months) for reference
    a=S.copy(); a['mo']=a.age.str[:-1].astype(int); m=smf.ols(f'Q("{c}")~mo',a).fit()
    rows.append(dict(score=c,cohort='age_trend_months',n=len(a),beta_std=m.params['mo']/a[c].std()*12,se=m.bse['mo'],p=m.pvalues['mo'],model='score~months (beta=SD per year)'))
R=pd.DataFrame(rows)
for coh in R.cohort.unique():
    i=R.cohort==coh; R.loc[i,'q']=multipletests(R.loc[i,'p'],method='fdr_bh')[1]
R.to_csv(T+'2e_GSE302188_cognition_models.csv',index=False)
print(R[R.cohort=='old24_28'].sort_values('p').head(25).to_string())
# ---------- GSE13799
ann=pd.read_csv(D+'data/ref/GPL1261.annot.gz',sep='\t',skiprows=27,usecols=['ID','Gene symbol'],low_memory=False).dropna()
ann=ann[~ann['Gene symbol'].str.contains('///')]
mx=pd.read_csv(D+'data/B12_GSE13799/GSE13799_series_matrix.txt.gz',sep='\t',comment='!',index_col=0)
lg=np.log2(mx.clip(lower=1)); lg=lg[(lg>np.log2(100)).mean(axis=1)>=0.5]
lg=lg.join(ann.set_index('ID')); lg['mean']=lg.drop(columns='Gene symbol').mean(axis=1)
lg=lg.sort_values('mean',ascending=False).drop_duplicates('Gene symbol').set_index('Gene symbol').drop(columns='mean')
grp=['impaired']*6+['unimpaired']*6+['young']*6; g=pd.Series(grp,index=mx.columns)
S2,N2=score(lg); S2['group']=g; S2.to_csv(T+'2e_GSE13799_scores.csv'); N2.to_csv(T+'2e_geneset_sizes_GSE13799.csv')
r2=[]
for c in [c for c in S2.columns if c!='group']:
    a=S2.loc[S2.group=='unimpaired',c]; b=S2.loc[S2.group=='impaired',c]; y=S2.loc[S2.group=='young',c]
    t,p=stats.ttest_ind(a,b); pooled=np.sqrt((a.var()+b.var())/2)
    ta,pa=stats.ttest_ind(pd.concat([a,b]),y)
    r2.append(dict(score=c,d_unimp_vs_imp=(a.mean()-b.mean())/pooled,t=t,p=p,p_mw=stats.mannwhitneyu(a,b).pvalue,aged_vs_young_t=ta,aged_vs_young_p=pa))
R2=pd.DataFrame(r2); R2['q']=multipletests(R2.p,method='fdr_bh')[1]; R2.to_csv(T+'2e_GSE13799_impaired_vs_unimpaired.csv',index=False)
# ---------- replication table
o=R[R.cohort=='old24_28'].set_index('score'); rep=o[['beta_std','p','q']].join(R2.set_index('score')[['d_unimp_vs_imp','p','q']],rsuffix='_GSE13799',how='inner')
rep['same_direction']=np.sign(rep.beta_std)==np.sign(rep.d_unimp_vs_imp); rep.to_csv(T+'2e_replication_GSE302188_vs_GSE13799.csv')
print(rep.sort_values('p').head(25).to_string())
print('direction agreement (all scores):',rep.same_direction.mean(), ' among p<0.05 in discovery:',rep[rep.p<0.05].same_direction.mean(), (rep.p<0.05).sum())
# ---------- figure
show=[c for c in cols if c.startswith(('M6A_','MK_'))]
P=R[(R.cohort=='old24_28')&R.score.isin(show)].set_index('score').loc[show]
fig,ax=plt.subplots(1,2,figsize=(13,max(6,len(show)*0.22)))
y=np.arange(len(show)); ax[0].errorbar(P.beta_std,y,xerr=1.96*P.se,fmt='o',color='k',ms=3); ax[0].axvline(0,color='grey',ls='--')
ax[0].set_yticks(y); ax[0].set_yticklabels([f"{s} (q={q:.2g})" for s,q in zip(P.index,P.q)],fontsize=7); ax[0].set_xlabel('std. beta for Barnes cognition (24+28 mo, adj. age group)'); ax[0].set_title(f'GSE302188 old mice (n={int(P.n.iloc[0])})')
Q=R2.set_index('score').reindex(show); ax[1].barh(y,Q.d_unimp_vs_imp,color=['#3b75af' if v>0 else '#d1495b' for v in Q.d_unimp_vs_imp]); ax[1].set_yticks(y); ax[1].set_yticklabels([f"p={p:.2g}" for p in Q.p],fontsize=7); ax[1].axvline(0,color='grey')
ax[1].set_xlabel("Cohen's d (unimpaired - impaired)"); ax[1].set_title('GSE13799 aged 23 mo (6 vs 6)')
plt.tight_layout(); plt.savefig(F+'Fig2e_cognition_association.png',dpi=300)
