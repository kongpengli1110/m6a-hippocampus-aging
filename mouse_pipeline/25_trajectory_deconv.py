"""Task 4: multi-age trajectories of METTL3 program / hypo set / excitatory modules in GSE212336 (male, 3-28 mo, per-mouse mean of ant+post HIP)
and GSE302188 (3-28 mo, n=82), change-point by segmented (broken-stick) regression with grid search over breakpoints + bootstrap CI;
NNLS deconvolution of GSE302188 & GSE212336 with Allen HPF young reference (mean CPM over donors, marker genes) ->
test whether program age effect persists after adjusting for estimated composition."""
import pandas as pd,numpy as np,statsmodels.formula.api as smf,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt,gzip
from scipy.optimize import nnls; from statsmodels.stats.multitest import multipletests
D='./'; T=D+'analysis/tables/'; F=D+'analysis/figures/'
rng=np.random.default_rng(5)
gs=pd.read_csv(T+'genesets_m6a_aging.csv'); G={k:set(v) for k,v in gs.groupby('set').gene}
rg=pd.read_csv(T+'genesets_regulator_targets.csv'); RG={k:set(v) for k,v in rg.groupby('set').gene}
SETS={'METTL3cKO_HIP_down':RG['METTL3cKO_HIP_down'],'METTL3KO_Cb_m6Aloss':RG['METTL3KO_Cb_m6Aloss'],'A1_HIP_hypo_union':G['A1_CA1_hypo']|G['A1_CA3_hypo']|G['A1_DG_hypo'],'A1_core87_hypo':G['A1_CA1_hypo']&G['A1_CA3_hypo']&G['A1_DG_hypo']}
mm=pd.read_csv(T+'mod_gene_membership_kME.csv'); ae=pd.read_csv(T+'mod_age_effect_donor.csv')
en=pd.read_csv(T+'mod_enrichment.csv')
selm=en[(en.q<0.01)&en.geneset.isin(['A1_HIP_hypo_union','METTL3cKO_HIP_down','METTL3KO_Cb_m6Aloss'])&(en['size']<=300)].module.unique()  # pre-specified by m6A-set enrichment (not by age)
for m in selm:
    g=mm[(mm.module==m)].sort_values('kME_own',ascending=False).head(50).gene; SETS['MOD_'+m]=set(g)
# ---- Allen reference for deconvolution
mean=pd.read_csv(T+'allen_HPF_meanCPM_by_type.csv',index_col=0); spec=pd.read_csv(T+'allen_HPF_specificity_EWCE.csv',index_col=0)
# collapse excitatory subtypes (collinear in NNLS) into one Glut reference (unweighted mean of CA1/CA3/DG/SUB profiles)
mean['Glut']=mean[['CA1','CA3','DG','SUB']].mean(axis=1); spec=mean.div(mean[['Glut','GABA','Astro','Oligo','OPC','Microglia','Endo','Mural','VLMC']].sum(axis=1),axis=0)
REFT=['Glut','GABA','Astro','Oligo','OPC','Microglia','Endo','Mural','VLMC']
mk=set()
for t in REFT: s=spec[t][mean[t]>=20].sort_values(ascending=False); mk|=set(s.head(50).index)
def deconv(L):   # L: genes x samples linear CPM
    g=[x for x in mk if x in L.index]; Rf=mean.loc[g,REFT].values; out={}
    for s in L.columns:
        w,_=nnls(Rf,L.loc[g,s].values); out[s]=w/w.sum()
    return pd.DataFrame(out,index=REFT).T
def scores(L):  # log2CPM genes x samples
    Z=L.sub(L.mean(axis=1),axis=0).div(L.std(axis=1)+1e-9,axis=0)
    return pd.DataFrame({k:Z.loc[[g for g in v if g in Z.index]].mean() for k,v in SETS.items() if len([g for g in v if g in Z.index])>=10})
def brokenstick(age,y,grid):
    best=None
    for b in grid:
        X=np.column_stack([np.ones_like(age),age,np.clip(age-b,0,None)]); beta,res,*_=np.linalg.lstsq(X,y,rcond=None); rss=((y-X@beta)**2).sum()
        if best is None or rss<best[1]: best=(b,rss,beta)
    lin=np.polyfit(age,y,1); rss0=((y-np.polyval(lin,age))**2).sum(); n=len(y)
    F=((rss0-best[1])/1)/(best[1]/(n-3)); return best[0],best[2],F,rss0,best[1]
def run(name,L,meta,agecol='mo'):
    S=scores(np.log2(L+1)); P=deconv(L); S=S.join(meta); P.to_csv(T+f'traj_{name}_deconv_props.csv')
    rows=[]
    for c in [c for c in S.columns if c in SETS]:
        d=S.dropna(subset=[c]).copy(); age=d[agecol].values.astype(float); y=((d[c]-d[c].mean())/d[c].std()).values
        grid=np.unique(age)[1:-1]; b,beta,Fst,_,_=brokenstick(age,y,grid)
        bs=[]
        for _ in range(500):
            i=rng.integers(0,len(y),len(y))
            if len(np.unique(age[i]))<4: continue
            bs.append(brokenstick(age[i],y[i],grid)[0])
        m0=smf.ols('y~a',pd.DataFrame({'y':y,'a':age})).fit()
        dd=pd.DataFrame({'y':y,'a':age}).join(P.loc[d.index].reset_index(drop=True)[['Glut','Oligo','Microglia','Astro']].rename(columns=lambda x:'p_'+x))
        m1=smf.ols('y~a+p_Glut+p_Oligo+p_Microglia+p_Astro',dd).fit()
        rows.append(dict(dataset=name,score=c,n=len(y),slope_SD_per_month=m0.params['a'],p_linear=m0.pvalues['a'],slope_adj_composition=m1.params['a'],p_adj_composition=m1.pvalues['a'],
                         breakpoint_mo=b,bp_CI_lo=np.percentile(bs,2.5),bp_CI_hi=np.percentile(bs,97.5),slope_before=beta[1],slope_after=beta[1]+beta[2],F_vs_linear=Fst))
    R=pd.DataFrame(rows); R['q_linear']=multipletests(R.p_linear,method='fdr_bh')[1]; R['q_adj']=multipletests(R.p_adj_composition,method='fdr_bh')[1]
    pa=[]
    for t in REFT:
        f=smf.ols('p~a',pd.DataFrame({'p':P[t].values,'a':S.loc[P.index,agecol].values.astype(float)})).fit(); pa.append(dict(dataset=name,celltype=t,mean_prop=P[t].mean(),slope_per_month=f.params['a'],p=f.pvalues['a']))
    return S,R,P,pd.DataFrame(pa)
# GSE212336 male HIP, one value per mouse
L=[l.rstrip('\n').split('\t') for l in gzip.open(D+'data/B10_GSE212336/GSE212336_series_matrix.txt.gz','rt') if l.startswith(('!Sample_title','!Sample_description'))]
tit=[x.strip('"') for x in L[0][1:]]; col=[x.strip('"') for x in L[3][1:]]
meta=pd.DataFrame({'title':tit,'col':col}); meta[['tissue','sex','age','mouse']]=meta.title.str.extract(r'^(.*) (Male|Female) (\d+)mo (\d+)ID$')
meta=meta[meta.tissue.str.startswith('Hippocampus')&(meta.sex=='Male')]; meta['mo']=meta.age.astype(int)
cnt=pd.read_csv(D+'data/B10_GSE212336/GSE212336_BulkSeq_Counttable.txt.gz',sep='\t',index_col=0)[meta.col]
cnt=cnt.T.groupby(meta.mouse.values).sum().T; mmeta=meta.groupby('mouse').agg(mo=('mo','first'))
cpm=cnt/cnt.sum()*1e6; cpm=cpm[(cpm>=1).mean(axis=1)>=0.5]
S1,R1,P1,PA1=run('GSE212336',cpm,mmeta)
# GSE302188
c2=pd.read_csv(D+'data/B11_GSE302188/GSE302188_brain_RNA_gene_counts.csv.gz',index_col=0)
e2s=pd.read_csv(D+'data/B11_GSE302188/zenodo/GRCm38.102.EnsemblToGeneSymbol.tsv.gz',sep='\t',header=None,index_col=0)[1]
cp2=c2/c2.sum()*1e6; cp2=cp2[(cp2>=1).mean(axis=1)>=0.5]; cp2.index=cp2.index.map(e2s); cp2=cp2[cp2.index.notna()].groupby(level=0).sum()
m2=pd.DataFrame({'mo':[int(c.split('M')[0]) for c in cp2.columns]},index=cp2.columns)
S2,R2,P2,PA2=run('GSE302188',cp2,m2)
R=pd.concat([R1,R2]); R.to_csv(T+'traj_breakpoints_composition.csv',index=False); PA=pd.concat([PA1,PA2]); PA.to_csv(T+'traj_deconv_age_slopes.csv',index=False)
print(R.round(4).to_string()); print(PA.round(4).to_string())
show=[c for c in SETS if c in S1.columns][:8]
fig,axs=plt.subplots(2,len(show),figsize=(3.2*len(show),6.5),squeeze=False)
for j,c in enumerate(show):
    for i,(nm,S,Rr) in enumerate([('GSE212336 (male, per mouse)',S1,R1),('GSE302188',S2,R2)]):
        ax=axs[i,j]; z=(S[c]-S[c].mean())/S[c].std(); g=pd.DataFrame({'z':z,'mo':S.mo}).groupby('mo').z.agg(['mean','sem'])
        ax.scatter(S.mo+rng.normal(0,.3,len(S)),z,s=6,alpha=.4,color='grey'); ax.errorbar(g.index,g['mean'],yerr=g['sem'],color='k',marker='o',ms=3)
        r=Rr.set_index('score').loc[c]; ax.axvline(r.breakpoint_mo,color='r',ls='--',lw=.8)
        ax.set_title(f"{c}\n{nm.split()[0]} q={r.q_linear:.2g}; adj q={r.q_adj:.2g}\nbp={r.breakpoint_mo:.0f} [{r.bp_CI_lo:.0f}-{r.bp_CI_hi:.0f}]",fontsize=7)
        if j==0: ax.set_ylabel('score (z)')
        if i==1: ax.set_xlabel('age (mo)')
plt.tight_layout(); plt.savefig(F+'P2b_Fig5_trajectory_breakpoints.png',dpi=300)
fig,ax=plt.subplots(1,2,figsize=(12,4))
for a,(P,S,nm) in zip(ax,[(P1,S1,'GSE212336'),(P2,S2,'GSE302188')]):
    P.assign(mo=S.mo).groupby('mo').mean().plot.bar(stacked=True,ax=a,colormap='tab20',legend=(nm=='GSE302188')); a.set_title(nm+' NNLS composition (Allen HPF ref)'); a.set_ylabel('proportion')
ax[1].legend(fontsize=6,bbox_to_anchor=(1,1)); plt.tight_layout(); plt.savefig(F+'P2b_FigS_deconvolution.png',dpi=300)
