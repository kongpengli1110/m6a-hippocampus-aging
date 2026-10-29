"""Task 1: replication in GSE161340 (sn & sc; 2 vs 2 pooled libraries), GSE233363 (DG; 1 young vs 2 old libraries, descriptive),
GSE129788 (whole brain; 8 vs 8 animals) + Allen HPF (8 vs 15 donors) for reference.
(a) EWCE of A1 hippocampal aging-DM genes using each dataset's own young-sample cell-type specificity.
(b) m6A regulator pseudobulk DE (PyDESeq2 ~age; B3 descriptive log2FC only).
(c) Gene-set age effects per cell type, competitive & length-adjusted: OLS over methylated genes
    lfc_g ~ inSet + log10(gene length) + mean log2CPM  -> coef (log2 units) +/- SE  (forest plot);
    and sample-level score (mean z of set genes) aged-vs-young."""
import sys; sys.path.insert(0,'./scripts')
import pydeseq2_patch
from common import *
import statsmodels.api as sm, warnings; warnings.filterwarnings('ignore')
from statsmodels.stats.multitest import multipletests
from pydeseq2.dds import DeseqDataSet; from pydeseq2.ds import DeseqStats; from pydeseq2.default_inference import DefaultInference
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, seaborn as sns
T=D+'analysis/tables/'; F=D+'analysis/figures/'; inf=DefaultInference(n_cpus=1); rng=np.random.default_rng(3)
gl=pd.read_csv(D+'data/ref/mm10_gencodeM25_genes.tsv',sep='\t',header=None); GL=(gl[4]-gl[3]).groupby(gl[1]).max()
gs=pd.read_csv(T+'genesets_m6a_aging.csv'); G={k:set(v) for k,v in gs.groupby('set').gene}
rg=pd.read_csv(T+'genesets_regulator_targets.csv'); RG={k:set(v) for k,v in rg.groupby('set').gene}
bg=pd.read_csv(T+'merip_backgrounds.csv'); METH=set(bg[bg.region.isin(['A1_CA1','A1_CA3','A1_DG'])].gene)
HYPO=G['A1_CA1_hypo']|G['A1_CA3_hypo']|G['A1_DG_hypo']
SETS={'A1_HIP_hypo_union':HYPO,'A1_core87_hypo':G['A1_CA1_hypo']&G['A1_CA3_hypo']&G['A1_DG_hypo'],'A1_DG_hyper':G['A1_DG_hyper'],
      'METTL3cKO_HIP_down':RG['METTL3cKO_HIP_down'],'METTL3cKO_HIP_up':RG['METTL3cKO_HIP_up'],'METTL3KO_Cb_m6Aloss':RG['METTL3KO_Cb_m6Aloss'],'YTHDF1_CLIP_HIP':RG['YTHDF1_CLIP_HIP'],'FTOcKO_HIP_up':RG['FTOcKO_HIP_up']}
# harmonise cell-type names
MAPB4={'mNEUR':'Neuron(all)','ASC':'Astro','OLG':'Oligo','OPC':'OPC','MG':'Microglia','EC':'Endo','PC':'Mural','VSMC':'Mural','VLMC':'VLMC','ImmN':'DG_IMN','MAC':'BAM'}
def load(ds):
    if ds=='Allen':
        X,M,genes=load_pb(); sym=dict(zip(genes.gene_identifier,genes.gene_symbol)); X=X.T.groupby(X.columns.map(sym)).sum().T
        M=M.rename(columns={'major':'celltype','donor':'sample'}); M['age']=M.age.map({'adult':'young','aged':'old'}); M=M[M.n_cells>=20]; return X.loc[M.index],M
    if ds=='GSE129788':
        X=pd.read_parquet(D+'data/B4_GSE129788/pb.parquet'); M=pd.read_csv(D+'data/B4_GSE129788/pb_meta.csv',index_col=0)
        M['celltype']=M.celltype.map(MAPB4); M=M.dropna(subset=['celltype']); M['key']=M['sample']+'|'+M.celltype
        X=X.loc[M.index].groupby(M.key.values).sum(); M=M.groupby('key').agg(sample=('sample','first'),celltype=('celltype','first'),age=('age','first'),n_cells=('n_cells','sum')).loc[X.index]
    else:
        f={'GSE161340_sn':'B2_GSE161340/pb_sn','GSE161340_sc':'B2_GSE161340/pb_sc','GSE233363_DG':'B3_GSE233363/pb'}[ds]
        X=pd.read_parquet(D+f'data/{f}.parquet'); M=pd.read_csv(D+f'data/{f}_meta.csv',index_col=0)
        M=M[M.celltype!='ambiguous']
        if ds=='GSE233363_DG': M=M[M.age!='middle']
    M=M[M.n_cells>=20]; return X.loc[M.index],M
DS=['Allen','GSE161340_sn','GSE161340_sc','GSE233363_DG','GSE129788']
ewce_rows=[];reg_rows=[];set_rows=[];score_rows=[]
for ds in DS:
    X,M=load(ds); X=X.loc[:,~X.columns.duplicated()]
    cpm=X.div(X.sum(1),axis=0)*1e6; L=np.log2(cpm+1)
    types=[t for t,g in M.groupby('celltype') if g.age.nunique()==2 and g.groupby('age').size().min()>=1 and t!='Other']
    print(ds,X.shape,{t:M[M.celltype==t].groupby('age').size().to_dict() for t in types},flush=True)
    # (a) EWCE using young samples' mean CPM per type
    yt=[t for t in types if ((M.celltype==t)&(M.age=='young')).sum()>=1]
    mean=pd.DataFrame({t:cpm[(M.celltype==t)&(M.age=='young')].mean() for t in yt}); spec=mean.div(mean.sum(1),axis=0).fillna(0); tot=mean.sum(1)
    bgl=[g for g in spec.index if g in METH and tot[g]>=1]; dec=pd.qcut(np.log1p(tot[bgl]),10,labels=False,duplicates='drop')
    S=spec.loc[bgl].values; idx={g:i for i,g in enumerate(bgl)}; pools={d:np.array([idx[g] for g in dec.index[dec==d]]) for d in dec.unique()}
    for sname,gset in [('A1_CA1_hypo',G['A1_CA1_hypo']),('A1_CA3_hypo',G['A1_CA3_hypo']),('A1_DG_hypo',G['A1_DG_hypo']),('A1_DG_hyper',G['A1_DG_hyper']),('A1_CA3_hyper',G['A1_CA3_hyper'])]:
        g=[x for x in gset if x in idx]; need=dec[g].value_counts(); obs=S[[idx[x] for x in g]].mean(0)
        null=np.array([S[np.concatenate([rng.choice(pools[d],n,replace=False) for d,n in need.items()])].mean(0) for _ in range(5000)])
        z=(obs-null.mean(0))/null.std(0); p=(1+(null>=obs).sum(0))/5001
        for t,zz,pp in zip(yt,z,p): ewce_rows.append(dict(dataset=ds,geneset=sname,celltype=t,n_genes=len(g),z=zz,p=pp))
    # (b)+(c)
    for t in types:
        sel=(M.celltype==t); m=M[sel]; c=X[sel.values]
        ny,no=(m.age=='young').sum(),(m.age=='old').sum()
        keep=(cpm[sel.values]>=5).mean()>=0.5; c=c.loc[:,keep]; l=L[sel.values].loc[:,keep]
        lfc=l[(m.age=='old').values].mean()-l[(m.age=='young').values].mean(); ave=l.mean()
        # regulator
        regs=[g for g in REG2CLS if g in c.columns]
        if ny>=2 and no>=2:
            meta=m[['age']].copy(); des='~age'
            if ds=='Allen': meta=m[['age','sex','region']].copy(); des='~region+sex+age' if (m.groupby('region').age.nunique()==2).all() else '~sex+age'
            try:
                dds=DeseqDataSet(counts=c.astype(int).round(),metadata=meta,design=des,inference=inf,quiet=True); dds.deseq2()
                st=DeseqStats(dds,contrast=['age','old','young'],inference=inf,quiet=True); st.summary(); r=st.results_df
                for g in regs: reg_rows.append(dict(dataset=ds,celltype=t,gene=g,n_young=ny,n_old=no,log2FC=r.loc[g,'log2FoldChange'],lfcSE=r.loc[g,'lfcSE'],pvalue=r.loc[g,'pvalue'],padj=r.loc[g,'padj'],method='DESeq2'))
            except Exception as e: print('deseq fail',ds,t,e)
        else:
            for g in regs: reg_rows.append(dict(dataset=ds,celltype=t,gene=g,n_young=ny,n_old=no,log2FC=lfc[g],lfcSE=np.nan,pvalue=np.nan,padj=np.nan,method='descriptive_meanlog2CPM'))
        # gene-set competitive, length-adjusted
        U=[g for g in lfc.index if g in METH and g in GL.index]
        dfu=pd.DataFrame({'lfc':lfc[U],'len':np.log10(GL[U].values),'ave':ave[U]})
        for sname,gset in SETS.items():
            dfu['in']=dfu.index.isin(gset).astype(float)
            if dfu['in'].sum()<10: continue
            fit=sm.OLS(dfu.lfc,sm.add_constant(dfu[['in','len','ave']])).fit(); fit0=sm.OLS(dfu.lfc,sm.add_constant(dfu[['in']])).fit()
            set_rows.append(dict(dataset=ds,celltype=t,geneset=sname,n_young=ny,n_old=no,n_set=int(dfu['in'].sum()),n_universe=len(dfu),coef_lenadj=fit.params['in'],se_lenadj=fit.bse['in'],p_lenadj=fit.pvalues['in'],coef_raw=fit0.params['in'],se_raw=fit0.bse['in'],p_raw=fit0.pvalues['in']))
            # sample-level score
            z=l.sub(l.mean(),axis=1).div(l.std()+1e-9,axis=1); gg=[x for x in gset if x in z.columns]
            sc_=z[gg].mean(1); score_rows+= [dict(dataset=ds,celltype=t,geneset=sname,sample=s,age=a,score=v) for s,a,v in zip(m['sample'],m.age,sc_)]
EW=pd.DataFrame(ewce_rows); EW['q']=EW.groupby('dataset').p.transform(lambda p:multipletests(p,method='fdr_bh')[1]); EW.to_csv(T+'rep_EWCE_by_dataset.csv',index=False)
RE=pd.DataFrame(reg_rows); RE.to_csv(T+'rep_regulator_DE_by_dataset.csv',index=False)
SE=pd.DataFrame(set_rows); SE['q_lenadj']=SE.groupby(['dataset','geneset']).p_lenadj.transform(lambda p:multipletests(p,method='fdr_bh')[1]); SE.to_csv(T+'rep_geneset_age_effect_lenadj.csv',index=False)
SC=pd.DataFrame(score_rows); SC.to_csv(T+'rep_sample_scores.csv',index=False)
print(EW.pivot_table(index=['dataset','geneset'],columns='celltype',values='z').round(1).to_string())
print(SE[SE.geneset.isin(['A1_HIP_hypo_union','METTL3cKO_HIP_down'])].pivot_table(index=['geneset','dataset'],columns='celltype',values='coef_lenadj').round(3).to_string())
