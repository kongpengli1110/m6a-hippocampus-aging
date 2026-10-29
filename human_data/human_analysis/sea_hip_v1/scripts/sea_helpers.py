from pathlib import Path
import sys
import numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'gse48350_v1/scripts'))
from fixed_pipeline import design_gate,fit_hc3

def strict_bool(s):
 a=s.astype(str);assert set(a.unique()) <= {'True','False'};return a.map({'True':True,'False':False})

def weighted_fixed_score(x,R,A,w):
 x=x.loc[A];mu=x.loc[:,R].mean(axis=1,skipna=False);sd=x.loc[:,R].std(axis=1,ddof=1,skipna=False)
 if not A or not np.isfinite(x.loc[:,R]).all().all() or not (sd>1e-8).all():return {'status':'not_testable_score','A':A,'mu':mu,'sd':sd}
 weights=w.loc[A].astype(float);assert (weights>0).all();z=x.sub(mu,axis=0).div(sd,axis=0);S=pd.Series(weights.to_numpy()@z.to_numpy()/weights.sum(),index=x.columns);sm=S.loc[R].mean(skipna=False);ss=S.loc[R].std(ddof=1,skipna=False)
 if not np.isfinite(ss) or ss<=1e-8:return {'status':'not_testable_score','A':A,'mu':mu,'sd':sd}
 Y=(S-sm)/ss
 return {'status':'ready','A':A,'mu':mu,'sd':sd,'weights':weights,'S':S,'Y':Y,'S_mean':float(sm),'S_SD':float(ss),'weight_denominator':float(weights.sum())}

def full_background_cpm(counts,L):
 a=np.asarray(counts);L=np.asarray(L);assert (a>=0).all() and (L>=0).all();return np.divide(a*1e6,L[:,None],out=np.full(a.shape,np.nan,dtype=float),where=L[:,None]>0)

def coverage(A,G0,weights):
 minimum=max(10,int(np.ceil(.8*len(G0))));frac=len(A)/len(G0);wf=float(weights.loc[A].sum()/weights.loc[G0].sum())
 return {'A_n':len(A),'G0_n':len(G0),'minimum_genes':minimum,'number_coverage':frac,'weight_coverage':wf,'coverage_gate_pass':len(A)>=minimum and wf>=.8}
