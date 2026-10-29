"""Task 3: hdWGCNA-like co-expression modules within Allen HPF excitatory neurons (CA1, DG; CA3/SUB pooled as 'CA3SUB' if small).
Steps per cell type: metacells = k-NN aggregation (k=25 cells, within donor, PCA space, max shared 10%) of raw counts -> log2CPM;
genes = top 3000 HVG (by metacell variance) U {methylated & METTL3/hypo-set genes expressed} capped at 4000;
signed network: adjacency = ((1+cor)/2)^beta, beta = smallest power with scale-free R^2>=0.8 (signed); TOM; average-linkage on 1-TOM;
dynamic-tree-like cut (scipy fcluster on distance with min module size 30) -> modules; eigengene = PC1 of scaled module genes;
kME = cor(gene, ME). Module age effect: ME averaged per donor (pseudobulk by animal) -> OLS ME ~ age + sex (+ROI).
Enrichment: hypergeometric of METTL3cKO_HIP_down / METTL3KO_Cb_m6Aloss / A1 hypo union / YTHDF1 CLIP vs network gene universe."""
import scanpy as sc,anndata as ad,numpy as np,pandas as pd,scipy.sparse as sp,sys,warnings; warnings.filterwarnings('ignore')
from scipy.cluster.hierarchy import linkage,fcluster; from scipy.spatial.distance import squareform
from scipy.stats import hypergeom; import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
D='./'; T=D+'analysis/tables/'
A=ad.read_h5ad(D+'data/allen/excit_HPF_raw.h5ad'); A.var_names=A.var.gene_symbol.astype(str); A.var_names_make_unique()
A.obs['ct']=A.obs.subclass.map({'016 CA1-ProS Glut':'CA1','037 DG Glut':'DG','017 CA3 Glut':'CA3','023 SUB-ProS Glut':'SUB'})
gs=pd.read_csv(T+'genesets_m6a_aging.csv'); G={k:set(v) for k,v in gs.groupby('set').gene}
rg=pd.read_csv(T+'genesets_regulator_targets.csv'); RG={k:set(v) for k,v in rg.groupby('set').gene}
HYPO=G['A1_CA1_hypo']|G['A1_CA3_hypo']|G['A1_DG_hypo']
TEST={'A1_HIP_hypo_union':HYPO,'A1_DG_hyper':G['A1_DG_hyper'],'METTL3cKO_HIP_down':RG['METTL3cKO_HIP_down'],'METTL3cKO_HIP_up':RG['METTL3cKO_HIP_up'],'METTL3KO_Cb_m6Aloss':RG['METTL3KO_Cb_m6Aloss'],'YTHDF1_CLIP_HIP':RG['YTHDF1_CLIP_HIP'],'A3_ctx_Glut_hyper':G['A3_ctx_Glut_hyper']}
rng=np.random.default_rng(11)
def metacells(B,k=25):
    sc.pp.normalize_total(B,1e4); sc.pp.log1p(B); sc.pp.highly_variable_genes(B,n_top_genes=2000); sc.tl.pca(B,n_comps=30,mask_var='highly_variable')
    from sklearn.neighbors import NearestNeighbors
    rows=[];meta=[]
    for d,idx in B.obs.groupby('donor_label').indices.items():
        if len(idx)<k: continue
        Xp=B.obsm['X_pca'][idx]; nn=NearestNeighbors(n_neighbors=k).fit(Xp); _,nb=nn.kneighbors(Xp)
        used=np.zeros(len(idx)); order=rng.permutation(len(idx)); chosen=[]
        for i in order:
            s=nb[i]
            if (used[s]>0).mean()<=0.1 or len(chosen)==0: chosen.append(s); used[s]+=1
            if len(chosen)>=max(5,len(idx)//k*2): break
        for s in chosen: rows.append(idx[s]); meta.append(dict(donor=d,age=B.obs.donor_age_category.iloc[idx[0]],sex=B.obs.donor_sex.iloc[idx[0]],roi=B.obs.region_of_interest_label.iloc[idx[0]]))
    C=B.layers['counts']; MC=sp.vstack([sp.csr_matrix(C[r].sum(0)) for r in rows]).toarray()
    return pd.DataFrame(MC,columns=B.var_names),pd.DataFrame(meta)
def soft_threshold(Cr):
    k_res=[]
    for beta in [4,6,8,10,12,14,16,18,20]:
        Aj=((1+Cr)/2)**beta; np.fill_diagonal(Aj,0); k=Aj.sum(0)
        h,e=np.histogram(k,bins=20); x=(e[:-1]+e[1:])/2; ok=h>0
        sl,ic=np.polyfit(np.log10(x[ok]),np.log10(h[ok]/h.sum()),1); r2=np.corrcoef(np.log10(x[ok]),np.log10(h[ok]/h.sum()))[0,1]**2*np.sign(-sl)
        k_res.append((beta,r2,k.mean()))
        if r2>=0.8: return beta,k_res
    return max(k_res,key=lambda x:x[1])[0],k_res
out_mod=[];out_gene=[];out_enr=[];out_age=[];sft=[]
for ct in (sys.argv[1:] or ['CA1','DG','CA3','SUB']):
    B=A[A.obs.ct==ct].copy(); B.layers['counts']=B.X.copy()
    MC,meta=metacells(B); cpm=MC.div(MC.sum(1),axis=0)*1e6; L=np.log2(cpm+1)
    expd=L.columns[(cpm>=10).mean()>=0.5]
    var=L[expd].var().sort_values(ascending=False); genes=list(var.index[:3000])
    extra=[g for g in (TEST['METTL3cKO_HIP_down']|TEST['A1_HIP_hypo_union']) if g in expd and g not in genes]; genes=(genes+extra)[:4000]
    X=L[genes]; Xs=(X-X.mean())/X.std(); Cr=np.corrcoef(Xs.T.values)
    beta,res=soft_threshold(Cr); sft.append(dict(celltype=ct,beta=beta,n_metacells=len(MC),n_genes=len(genes),fits=str([(b,round(r,2),round(k,1)) for b,r,k in res])))
    Aj=((1+Cr)/2)**beta; np.fill_diagonal(Aj,0); k=Aj.sum(0)
    # TOM (signed)
    num=Aj@Aj+Aj; den=np.minimum.outer(k,k)+1-Aj; TOM=num/den; np.fill_diagonal(TOM,1)
    dist=1-TOM; np.fill_diagonal(dist,0); Z=linkage(squareform(dist,checks=False),'average')
    best=None
    for h in np.linspace(0.90,0.999,40):
        lab=fcluster(Z,t=h,criterion='distance'); sizes=pd.Series(lab).value_counts(); nmod=(sizes>=30).sum()
        if best is None or nmod>best[1]: best=(h,nmod,lab)
    lab=best[2]; sizes=pd.Series(lab).value_counts(); big=sizes[sizes>=30].index
    mods=pd.Series(np.where(pd.Series(lab).isin(big),lab,0),index=genes)
    ren={m:f'{ct}-M{i+1}' for i,m in enumerate(sizes[sizes>=30].index)}; ren[0]=f'{ct}-grey'; mods=mods.map(ren)
    MEs={}
    for m in [x for x in mods.unique() if not x.endswith('grey')]:
        g=mods.index[mods==m]; u,s,vt=np.linalg.svd(Xs[g].values-Xs[g].values.mean(0),full_matrices=False); me=u[:,0]*s[0]
        if np.corrcoef(me,Xs[g].mean(1))[0,1]<0: me=-me
        MEs[m]=me
    ME=pd.DataFrame(MEs)
    kME=pd.DataFrame({m:Xs.corrwith(ME[m]) for m in ME.columns})
    for g in genes: out_gene.append(dict(celltype=ct,gene=g,module=mods[g],kME_own=kME.loc[g,mods[g]] if mods[g] in kME else np.nan,connectivity=k[genes.index(g)]))
    U=set(genes)
    for m in ME.columns:
        mg=set(mods.index[mods==m])
        for sname,sset in TEST.items():
            a=len(mg&sset); K=len(U&sset); p=hypergeom.sf(a-1,len(U),K,len(mg)); out_enr.append(dict(celltype=ct,module=m,size=len(mg),geneset=sname,overlap=a,expected=K*len(mg)/len(U),fold=a/max(K*len(mg)/len(U),1e-9),p=p))
        d=meta.assign(me=ME[m].values).groupby('donor').agg(me=('me','mean'),age=('age','first'),sex=('sex','first'),roi=('roi','first'))
        f='me~C(age)+sex'+('+roi' if d.roi.nunique()>1 and (d.groupby('roi').age.nunique()==2).all() else '')
        fit=smf.ols(f,d).fit(); key=[x for x in fit.params.index if 'age' in x][0]
        out_age.append(dict(celltype=ct,module=m,size=len(mg),n_donors_young=(d.age=='adult').sum(),n_donors_old=(d.age=='aged').sum(),beta_aged=fit.params[key]/d.me.std(),p=fit.pvalues[key],model=f,
                            hubs=';'.join(kME[m].loc[list(mg)].sort_values(ascending=False).head(10).index)))
    print(ct,'beta',beta,'modules',len(ME.columns),'grey',(mods.str.endswith('grey')).sum(),flush=True)
R=pd.DataFrame(out_enr); R['q']=multipletests(R.p,method='fdr_bh')[1]; R.to_csv(T+'mod_enrichment.csv',index=False)
AG=pd.DataFrame(out_age); AG['q']=multipletests(AG.p,method='fdr_bh')[1]; AG.to_csv(T+'mod_age_effect_donor.csv',index=False)
pd.DataFrame(out_gene).to_csv(T+'mod_gene_membership_kME.csv',index=False); pd.DataFrame(sft).to_csv(T+'mod_soft_threshold.csv',index=False)
print(AG.sort_values('p').head(20).to_string()); print(R[R.q<0.05].sort_values('p').head(30).to_string())
