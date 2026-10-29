"""Sensitivity for task 4: (i) deconvolution markers EXCLUDING all tested-set genes (avoid circularity);
(ii) adjust only for glial proportions (Microglia, Astro, Oligo), not for the neuronal fraction that the program itself tracks."""
import pandas as pd,numpy as np,statsmodels.formula.api as smf,gzip
from scipy.optimize import nnls; from statsmodels.stats.multitest import multipletests
exec(open('./scripts/25_trajectory_deconv.py').read().split('# GSE212336 male HIP')[0].replace("def run(","def _unused("))
allset=set().union(*SETS.values()); mk2=mk-allset
def deconv2(L):
    g=[x for x in mk2 if x in L.index]; Rf=mean.loc[g,REFT].values
    return pd.DataFrame({s:(lambda w:w/w.sum())(nnls(Rf,L.loc[g,s].values)[0]) for s in L.columns},index=REFT).T
src=open('./scripts/25_trajectory_deconv.py').read()
blk=src.split('# GSE212336 male HIP')[1].split('S1,R1,P1,PA1')[0].split('\n',1)[1]; exec(blk)
c1=cpm
blk2=src.split('\n# GSE302188')[1].split('S2,R2,P2,PA2')[0]; exec(blk2)
rows=[]
for name,L,meta in [('GSE212336',c1,mmeta),('GSE302188',cp2,m2)]:
    S=scores(np.log2(L+1)).join(meta); P=deconv2(L)
    for c in [c for c in S.columns if c in SETS]:
        d=pd.DataFrame({'y':(S[c]-S[c].mean())/S[c].std(),'a':S.mo.astype(float)}).join(P.rename(columns=lambda x:'p_'+x))
        for lab,f in [('none','y~a'),('glia_only','y~a+p_Microglia+p_Astro+p_Oligo'),('all_incl_Glut','y~a+p_Glut+p_Microglia+p_Astro+p_Oligo')]:
            m=smf.ols(f,d).fit(); rows.append(dict(dataset=name,score=c,adjust=lab,slope_SD_per_month=m.params['a'],p=m.pvalues['a']))
    print(name,'markers used',len([x for x in mk2 if x in L.index]),'Glut prop slope p',smf.ols('g~a',pd.DataFrame({'g':P.Glut,'a':S.mo.astype(float)})).fit().pvalues['a'])
R=pd.DataFrame(rows); R['q']=R.groupby(['dataset','adjust']).p.transform(lambda p:multipletests(p,method='fdr_bh')[1]); R.to_csv(T+'traj_composition_sensitivity.csv',index=False)
print(R.pivot_table(index=['dataset','score'],columns='adjust',values=['slope_SD_per_month','q']).round(3).to_string())
