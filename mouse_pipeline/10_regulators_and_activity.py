"""2c regulator expression by cell type x age (pseudobulk DESeq2 results) ; 2d target-set activity per cell type.
2d(i) competitive: decoupler ULM of each gene set on the per-cell-type DESeq2 Wald statistic (aged vs adult) -> which programs shift with age in which cell type.
2d(ii) donor-level: per cell type, ULM score on donor pseudobulk log-CPM (z within type); OLS score ~ sex + age (+region) -> donor-level p."""
import sys; sys.path.insert(0,'./scripts')
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, seaborn as sns, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import decoupler as dc
T=D+'analysis/tables/'; F=D+'analysis/figures/'
R=pd.read_csv(T+'allen_HPF_pseudobulk_DESeq2_age.csv')
types=[t for t in ['CA1','CA3','DG','SUB','ENT_Glut','GABA','Astro','Oligo','OPC','Microglia','Endo','Mural','VLMC','CR','BAM'] if t in R.celltype.unique()]
# ---- 2c
X,M,genes=load_pb(); sym=dict(zip(genes.gene_identifier,genes.gene_symbol))
cpm=X.div(X.sum(1),axis=0)*1e6; cpm.columns=cpm.columns.map(sym); cpm=cpm.T.groupby(level=0).sum().T
regs=[g for g in REG2CLS if g in cpm.columns]
ok=M.n_cells>=20
E=[]
for t in types:
    for a in ['adult','aged']:
        s=(M.major==t)&(M.age==a)&ok
        E.append(pd.Series(np.log2(cpm.loc[s,regs]+1).mean(),name=(t,a)))
E=pd.DataFrame(E); E.index=pd.MultiIndex.from_tuples(E.index,names=['celltype','age'])
E.to_csv(T+'2c_regulator_log2CPM_by_type_age.csv')
RR=R[R.gene.isin(regs)][['celltype','gene','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj','n_adult','n_aged','design']].copy()
RR['class']=RR.gene.map(REG2CLS); RR.to_csv(T+'2c_regulator_pseudobulk_DE_age.csv',index=False)
sig=RR[RR.padj<0.05].sort_values('padj'); print('regulator DE padj<0.05:',len(sig)); print(sig[['celltype','gene','log2FoldChange','padj']].to_string())
L=RR.pivot(index='gene',columns='celltype',values='log2FoldChange').reindex(index=regs,columns=types)
Pq=RR.pivot(index='gene',columns='celltype',values='padj').reindex(index=regs,columns=types)
fig,ax=plt.subplots(1,2,figsize=(17,8),gridspec_kw={'width_ratios':[1,1]})
ad=E.xs('adult',level='age').T.reindex(index=regs,columns=types)
sns.heatmap(ad.sub(ad.mean(1),axis=0),cmap='PuOr_r',center=0,ax=ax[0],cbar_kws={'label':'log2CPM (row-centred), 2 mo'})
ax[0].set_title('m6A regulator expression by HPF cell type (young adult)'); ax[0].set_yticklabels([f'{g} ({REG2CLS[g][0]})' for g in regs],fontsize=8)
an=Pq.map(lambda q:'' if pd.isna(q) else '***' if q<0.001 else '**' if q<0.01 else '*' if q<0.05 else '')
sns.heatmap(L.clip(-1,1),cmap='RdBu_r',center=0,annot=an,fmt='',ax=ax[1],cbar_kws={'label':'log2FC aged(18mo)/adult(2mo), pseudobulk DESeq2'})
ax[1].set_title('Aging change (donor-level pseudobulk; * padj<0.05)'); ax[1].set_yticklabels(regs,fontsize=8)
plt.tight_layout(); plt.savefig(F+'Fig2c_m6A_regulators_celltype_age.png',dpi=300)
# ---- 2d gene sets
gsA=pd.read_csv(T+'genesets_m6a_aging.csv'); gsR=pd.read_csv(T+'genesets_regulator_targets.csv')
keepA=['A1_CA1_hypo','A1_CA3_hypo','A1_DG_hypo','A1_CA3_hyper','A1_DG_hyper','A2_HIP_hyper','A2_HIP_hypo','A3_ctx_Glut_hyper','A3_ctx_Glut_hypo','A3_ctx_Astro_hyper','A3_ctx_Oligo_hyper','A3_ctx_Endo_hyper','A3_ctx_GABA_hyper']
net=pd.concat([gsA[gsA.set.isin(keepA)][['set','gene']],gsR[['set','gene']]]).rename(columns={'set':'source','gene':'target'}).drop_duplicates(); net['weight']=1.0
# (i) competitive on DE stat
St=R.dropna(subset=['stat']).pivot_table(index='celltype',columns='gene',values='stat',aggfunc='mean').loc[types]
St=St.fillna(0)
es,pv=dc.mt.ulm(data=St,net=net,tmin=10)
es.to_csv(T+'2d_ULM_on_DEstat_score.csv'); pv.to_csv(T+'2d_ULM_on_DEstat_pval.csv')
q=pd.DataFrame(multipletests(pv.values.ravel(),method='fdr_bh')[1].reshape(pv.shape),index=pv.index,columns=pv.columns); q.to_csv(T+'2d_ULM_on_DEstat_qval.csv')
order=[s for s in keepA if s in es.columns]+sorted([s for s in es.columns if s not in keepA])
fig,ax=plt.subplots(figsize=(15,7)); an=q[order].map(lambda v:'***' if v<0.001 else '**' if v<0.01 else '*' if v<0.05 else '')
sns.heatmap(es[order].clip(-8,8),cmap='RdBu_r',center=0,annot=an,fmt='',ax=ax,cbar_kws={'label':'ULM t (gene set on aged-vs-adult DESeq2 stat)'})
ax.set_title('Age-related shift of m6A target programs within each HPF cell type (Allen, donor pseudobulk; BH q: * <0.05)'); ax.set_ylabel(''); plt.xticks(fontsize=7)
plt.tight_layout(); plt.savefig(F+'Fig2d_target_program_activity_by_celltype.png',dpi=300)
# (ii) donor-level
rows=[]
for t in types:
    s=(M.major==t)&ok; m=M[s]
    lg=np.log2(cpm[s]+1); lg=lg.loc[:,(lg>1).mean()>=0.5]; z=(lg-lg.mean())/(lg.std()+1e-9)
    e2,p2=dc.mt.ulm(data=z,net=net,tmin=10)
    for gset in e2.columns:
        d=m.assign(y=e2[gset].values)
        f='y~sex+age' if (m.region.nunique()==1 or (m.groupby('region').age.nunique()<2).any()) else 'y~region+sex+age'
        fit=smf.ols(f,d).fit(); rows.append(dict(celltype=t,geneset=gset,n_adult=(m.age=='adult').sum(),n_aged=(m.age=='aged').sum(),coef_aged=fit.params['age[T.aged]'],p=fit.pvalues['age[T.aged]'],model=f))
DL=pd.DataFrame(rows); DL['q']=multipletests(DL.p,method='fdr_bh')[1]; DL.to_csv(T+'2d_donorlevel_ULM_age_OLS.csv',index=False)
print(DL.sort_values('p').head(30).to_string())
print('competitive q<0.05 counts per celltype:'); print((q<0.05).sum(1))
