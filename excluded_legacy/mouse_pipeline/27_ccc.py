"""Phase 3.1: cell-cell communication, Allen HPF (unbiased libraries), animal-level.
Resource: LIANA 'mouseconsensus' (liana-py 1.x bundled; OmniPath-derived consensus of CellPhoneDB/CellChatDB/ICELLNET/connectomeDB/CellTalkDB ... orthology-mapped).
Score per donor d, sender s, receiver r, pair (L,R): geometric mean sqrt(L_s * R_r) of donor pseudobulk CPM; complexes = min over subunits.
Filters: sender & receiver each >=20 cells in donor; L mean CPM>=10 in sender and R>=10 in receiver (mean over donors of either age); present in >=3 donors per age.
Test: per interaction OLS log2(score+1) ~ age + sex (+ROI) ; empirical p from 2000 permutations of age labels within ROI; BH over all tested.
Highlight universe (pre-specified): L or R in neuron-attributed aging-DM genes (A1 CA1/CA3/DG hypo|hyper with Allen excitatory share>=0.5) or handoff core candidates.
Replication of direction: GSE161340 snRNA (2 vs 2 pools) same score, sign agreement only."""
import sys; sys.path.insert(0,'./scripts')
from common import *
import statsmodels.api as sm, itertools, warnings; warnings.filterwarnings('ignore')
from statsmodels.stats.multitest import multipletests
T=D+'analysis/tables/'; rng=np.random.default_rng(17)
res=pd.read_csv(D+'data/ref/liana_mouseconsensus.csv')
X,M,genes=load_pb(); sym=dict(zip(genes.gene_identifier,genes.gene_symbol))
X=X.T.groupby(X.columns.map(sym)).sum().T; cpm=X.div(X.sum(1),axis=0)*1e6
TYPES=['CA1','CA3','DG','SUB','GABA','Astro','Oligo','OPC','Microglia','Endo','Mural']
M=M[(M.n_cells>=20)&M.major.isin(TYPES)]
spec=pd.read_csv(T+'allen_HPF_specificity_EWCE.csv',index_col=0); exc=spec[['CA1','CA3','DG','SUB']].sum(axis=1)
gs=pd.read_csv(T+'genesets_m6a_aging.csv'); G={k:set(v) for k,v in gs.groupby('set').gene}
dm=set().union(*[G[f'A1_{r}_{d}'] for r in ['CA1','CA3','DG'] for d in ['hypo','hyper']])
NEU_DM={g for g in dm if exc.get(g,0)>=0.5}
core=set(pd.read_csv(D+'handoff/mouse_candidates_v1_core.csv').mouse_symbol)
HL=NEU_DM|core
def units(x): return x.split('_')
def expr(df,g):  # df donors x genes ; complex -> min
    u=units(g); 
    if not all(x in df.columns for x in u): return None
    return df[u].min(axis=1)
donors=M.donor.unique(); dmeta=M.groupby('donor').agg(age=('age','first'),sex=('sex','first'),region=('region','first')).loc[donors]
E={t:cpm.loc[M.index[M.major==t]].set_axis(M.donor[M.major==t].values) for t in TYPES}
rows=[]
for s,r in itertools.product(TYPES,TYPES):
    Es,Er=E[s],E[r]; common=Es.index.intersection(Er.index)
    if dmeta.loc[common].groupby('age').size().min()<3 if len(common) else True: continue
    for L,R in zip(res.ligand,res.receptor):
        l=expr(Es.loc[common],L); rr=expr(Er.loc[common],R)
        if l is None or rr is None: continue
        if l.mean()<10 or rr.mean()<10: continue
        sc=np.log2(np.sqrt(l*rr)+1)
        rows.append(dict(sender=s,receiver=r,ligand=L,receptor=R,donors=';'.join(common),**{f'd_{d}':v for d,v in sc.items()}))
W=pd.DataFrame(rows); print('interactions tested',len(W),flush=True)
dcols=[c for c in W.columns if c.startswith('d_')]
S=W[dcols].values; dn=[c[2:] for c in dcols]; meta=dmeta.loc[dn]
age=(meta.age=='aged').astype(float).values; sex=(meta.sex=='M').astype(float).values; reg=pd.get_dummies(meta.region,drop_first=True).astype(float).values
def ols_t(Y,a):   # Y: n_int x n_donor (nan allowed) ; returns coef,t for age
    co=np.full(len(Y),np.nan); tt=np.full(len(Y),np.nan)
    for i,y in enumerate(Y):
        ok=~np.isnan(y); Xd=np.column_stack([np.ones(ok.sum()),a[ok],sex[ok],reg[ok]])
        Xd=Xd[:,np.linalg.matrix_rank(Xd)>0] 
        try:
            b,res_,rk,_=np.linalg.lstsq(Xd,y[ok],rcond=None); dof=ok.sum()-rk
            if dof<2: continue
            e=y[ok]-Xd@b; s2=(e@e)/dof; cov=s2*np.linalg.pinv(Xd.T@Xd); co[i]=b[1]; tt[i]=b[1]/np.sqrt(cov[1,1])
        except Exception: pass
    return co,tt
co,tt=ols_t(S,age)
# permutation within ROI
NP=2000; exceed=np.zeros(len(S)); regl=meta.region.values
for k in range(NP):
    ap=age.copy()
    for g in np.unique(regl):
        i=np.where(regl==g)[0]; ap[i]=rng.permutation(ap[i])
    _,tp=ols_t(S,ap); exceed+=(np.abs(tp)>=np.abs(tt))
    if k%200==0: print('perm',k,flush=True)
W['log2FC_aged']=co; W['t']=tt; W['p_perm']=(1+exceed)/(NP+1); W['q_perm']=multipletests(W.p_perm,method='fdr_bh')[1]
W['highlight']=W.ligand.apply(lambda x:any(u in HL for u in units(x)))|W.receptor.apply(lambda x:any(u in HL for u in units(x)))
W['q_perm_within_highlight']=np.nan; h=W.highlight; W.loc[h,'q_perm_within_highlight']=multipletests(W.loc[h,'p_perm'],method='fdr_bh')[1]
W['m6A_gene_role']=[('L' if any(u in HL for u in units(l)) else '')+('R' if any(u in HL for u in units(r)) else '') for l,r in zip(W.ligand,W.receptor)]
# replication direction: GSE161340 sn
Xb=pd.read_parquet(D+'data/B2_GSE161340/pb_sn.parquet'); Mb=pd.read_csv(D+'data/B2_GSE161340/pb_sn_meta.csv',index_col=0); Mb=Mb[(Mb.n_cells>=20)&Mb.celltype.isin(TYPES)]
cb=Xb.div(Xb.sum(1),axis=0)*1e6; Eb={t:cb.loc[Mb.index[Mb.celltype==t]].set_axis(Mb['sample'][Mb.celltype==t].values) for t in TYPES if (Mb.celltype==t).any()}
agb=Mb.groupby('sample').age.first()
rep=[]
for i,w in W.iterrows():
    if w.sender not in Eb or w.receiver not in Eb: rep.append(np.nan); continue
    cm=Eb[w.sender].index.intersection(Eb[w.receiver].index)
    l=expr(Eb[w.sender].loc[cm],w.ligand); rr=expr(Eb[w.receiver].loc[cm],w.receptor)
    if l is None or rr is None or agb.loc[cm].nunique()<2: rep.append(np.nan); continue
    sc=np.log2(np.sqrt(l*rr)+1); rep.append(sc[agb.loc[cm]=='old'].mean()-sc[agb.loc[cm]=='young'].mean())
W['GSE161340sn_log2FC']=rep; W['rep_same_sign']=np.sign(W.GSE161340sn_log2FC)==np.sign(W.log2FC_aged)
W.drop(columns=dcols).to_csv(T+'ccc_allen_animal_level.csv',index=False)
W[['sender','receiver','ligand','receptor']+dcols].to_csv(T+'ccc_allen_donor_scores.csv',index=False)
print('q_perm<0.05 all:',(W.q_perm<0.05).sum(),' highlight tested:',h.sum(),' highlight q<0.05 (within):',(W.q_perm_within_highlight<0.05).sum(),' q<0.1:',(W.q_perm_within_highlight<0.1).sum())
print(W[W.highlight].sort_values('p_perm').head(30)[['sender','receiver','ligand','receptor','m6A_gene_role','log2FC_aged','p_perm','q_perm','q_perm_within_highlight','GSE161340sn_log2FC']].round(4).to_string())
