"""Gene-length confound checks. (1) Do aging-hypomethylated (A1) genes simply = long genes? (2) Is the age-related down-shift of
hypomethylated sets in Allen pseudobulk explained by the known age-associated long-gene decline? -> regress DE stat on log(gene length)
(GENCODE vM25 genomic span) within each cell type, re-run ULM on residuals (universe=methylated genes).
(3) 2b EWCE re-run with null matched on expression x peak count x length tertile."""
import pandas as pd,numpy as np,decoupler as dc,statsmodels.api as sm,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt,seaborn as sns
from statsmodels.stats.multitest import multipletests; from scipy.stats import mannwhitneyu
D='./'; T=D+'analysis/tables/'; F=D+'analysis/figures/'
gl=pd.read_csv(D+'data/ref/mm10_gencodeM25_genes.tsv',sep='\t',header=None,names=['ens','sym','chr','start','end','type'])
gl['len']=gl.end-gl.start+1; L=gl.drop_duplicates('sym').set_index('sym')['len']; Le=gl.set_index('ens')['len']
gs=pd.read_csv(T+'genesets_m6a_aging.csv'); G={k:set(v) for k,v in gs.groupby('set').gene}
bg=pd.read_csv(T+'merip_backgrounds.csv'); meth=set(bg[bg.region.isin(['A1_CA1','A1_CA3','A1_DG'])].gene)
rows=[]
for s in ['A1_CA1_hypo','A1_CA3_hypo','A1_DG_hypo','A1_DG_hyper','A2_HIP_hyper','A3_ctx_Glut_hyper']:
    a=L.reindex(list(G[s]&meth)).dropna(); b=L.reindex(list(meth-G[s])).dropna()
    rows.append(dict(set=s,median_len_kb=a.median()/1e3,median_len_other_methylated_kb=b.median()/1e3,p_mw=mannwhitneyu(a,b).pvalue))
print(pd.DataFrame(rows).to_string()); pd.DataFrame(rows).to_csv(T+'2x_length_of_sets.csv',index=False)
R=pd.read_csv(T+'allen_HPF_pseudobulk_DESeq2_age.csv'); R['len']=R.ensembl.map(Le)
types=['CA1','CA3','DG','SUB','ENT_Glut','GABA','Astro','Oligo','OPC','Microglia','Endo','Mural','VLMC','CR','BAM']
res={}; lr=[]
for t in types:
    d=R[(R.celltype==t)&R.stat.notna()&R.len.notna()]
    X=sm.add_constant(np.log10(d.len)); f=sm.OLS(d.stat,X).fit(); lr.append(dict(celltype=t,slope_stat_per_log10len=f.params.iloc[1],p=f.pvalues.iloc[1],r2=f.rsquared))
    res[t]=pd.Series(f.resid.values,index=d.gene.values).groupby(level=0).mean()
LR=pd.DataFrame(lr); LR.to_csv(T+'2x_DEstat_vs_genelength.csv',index=False); print(LR.round(4).to_string())
Rs=pd.DataFrame(res).T.fillna(0); Rs=Rs[[g for g in Rs.columns if g in meth]]
gsR=pd.read_csv(T+'genesets_regulator_targets.csv')
keepA=['A1_CA1_hypo','A1_CA3_hypo','A1_DG_hypo','A1_CA3_hyper','A1_DG_hyper','A2_HIP_hyper','A2_HIP_hypo','A3_ctx_Glut_hyper','A3_ctx_Glut_hypo','A3_ctx_Astro_hyper','A3_ctx_Oligo_hyper','A3_ctx_Endo_hyper','A3_ctx_GABA_hyper']
net=pd.concat([gs[gs.set.isin(keepA)][['set','gene']],gsR[['set','gene']]]).rename(columns={'set':'source','gene':'target'}).drop_duplicates(); net['weight']=1.0
es,pv=dc.mt.ulm(data=Rs,net=net,tmin=10); q=pd.DataFrame(multipletests(pv.values.ravel(),method='fdr_bh')[1].reshape(pv.shape),index=pv.index,columns=pv.columns)
es.to_csv(T+'2d_ULM_lengthResid_methylated_score.csv'); q.to_csv(T+'2d_ULM_lengthResid_methylated_q.csv')
print(es[[c for c in keepA if c in es]].round(1).to_string())
order=[s for s in keepA if s in es.columns]+sorted([s for s in es.columns if s not in keepA])
fig,ax=plt.subplots(figsize=(15,7)); an=q[order].map(lambda v:'***' if v<0.001 else '**' if v<0.01 else '*' if v<0.05 else '')
sns.heatmap(es[order].clip(-6,6),cmap='RdBu_r',center=0,annot=an,fmt='',ax=ax,cbar_kws={'label':'ULM t (gene-length-residualised DE stat; methylated universe)'})
ax.set_title('Aging shift of m6A programs per HPF cell type after removing gene-length trend (Allen 2 vs 18 mo; BH q * <0.05)'); plt.xticks(fontsize=7); ax.set_ylabel('')
plt.tight_layout(); plt.savefig(F+'Fig2d_target_program_activity_lengthAdjusted.png',dpi=300)
