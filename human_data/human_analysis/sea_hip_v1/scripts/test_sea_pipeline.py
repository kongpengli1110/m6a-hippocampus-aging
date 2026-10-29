
# Release-only safety/path adapter; see provenance/adaptation_diffs.
from pathlib import Path as _ReleasePath
import sys as _release_sys
_release_home = next(p for p in _ReleasePath(__file__).resolve().parents if (p / 'release_support.py').is_file())
_release_sys.path.insert(0, str(_release_home))
from release_support import gate, release_root, resolve_source
gate(__file__, __name__)
from pathlib import Path
import json
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy import sparse,stats
from sea_helpers import *
ROOT=Path(__file__).resolve().parents[1]
# Component + two-library test. Donor A's two libraries were summed upstream;
# combine Sub components once, never multiply by n_libraries.
a=np.array([[1,2,97],[3,4,193],[5,6,289],[7,8,385]],dtype=np.int64)
# First two rows=A Sub1 libraries; third=A Sub2; fourth=B Sub1.
upstream=np.stack([a[0]+a[1],a[2],a[3]])
g=sparse.csr_matrix(np.array([[1,1,0],[0,0,1]]));combined=g@upstream
assert np.array_equal(combined[0],a[:3].sum(axis=0));assert combined[0].sum()==600 and combined[1].sum()==400
cpm=full_background_cpm(combined[:,:2],combined.sum(axis=1));assert np.isclose(cpm[0,0],9/600*1e6);assert not np.isclose(cpm[0,0],9/21*1e6)
assert strict_bool(pd.Series(['True','False'])).tolist()==[True,False]
# Fixed >=20 reference; additional >=10 donor must not alter A/mu/sd/score scaler.
x=pd.DataFrame([[1.,2.,4.,10.],[2.,3.,5.,15.],[4.,3.,1.,9.]],index=['g1','g2','g3'],columns=['d1','d2','d3','new10']);R=['d1','d2','d3'];w=pd.Series({'g1':.9,'g2':.6,'g3':.4});o=weighted_fixed_score(x,R,list(x.index),w);sub=weighted_fixed_score(x[R],R,list(x.index),w)
assert np.allclose(o['Y'][R],sub['Y']);assert np.isclose(o['Y'][R].std(ddof=1),1);assert np.isclose(o['Y'][R].mean(),0);manual=sum(w[g]*(x.loc[g,'new10']-o['mu'][g])/o['sd'][g] for g in x.index)/w.sum();assert np.isclose(o['Y']['new10'],(manual-o['S_mean'])/o['S_SD'])
x.loc['g1','new10']=np.nan;o2=weighted_fixed_score(x,R,list(x.index),w);assert pd.isna(o2['Y']['new10'])
# Known five-parameter SEA design and manual HC3/t/CI, before human results.
n=30;i=np.arange(n);age=np.linspace(69,99,n);X=pd.DataFrame({'Intercept':1.,'age10':(age-80)/10,'sex':i%2,'PMI_c':np.sin(i*.73),'RIN_c':np.cos(i*.41)});xx=X.to_numpy();truth=np.array([.4,-.6,.2,.15,-.3]);B=np.linalg.inv(xx.T@xx);e=np.sin(i*1.13)*(1+i/n);e=e-xx@(B@xx.T@e);y=xx@truth+e;r,diag,loo=fit_hc3(X,y,'age10');bb=B@xx.T@y;err=y-xx@bb;h=np.einsum('ij,jk,ik->i',xx,B,xx);cov=B@(xx.T@np.diag(err**2/(1-h)**2)@xx)@B;se=np.sqrt(cov[1,1]);p=2*stats.t.sf(abs(bb[1]/se),25)
assert np.allclose(bb,truth,atol=1e-12);assert np.allclose([r['beta'],r['SE_HC3'],r['raw_P']],[truth[1],se,p],rtol=1e-10,atol=1e-12)
d=X.assign(Y_score=y,age=age);_,_,gate=design_gate(d,['age10','sex','PMI_c','RIN_c']);assert gate['status']=='ready' and gate['df']==25
_,_,bad=design_gate(d.iloc[:19],['age10','sex','PMI_c','RIN_c']);assert bad['status']=='not_testable_design'
result={'status':'PASS','completed_at_utc':datetime.now(timezone.utc).isoformat(),'component_and_two_library_counts':True,'full_background_CPM_not_candidate_denominator':True,'False_string_not_truthy':True,'fixed20_reference_scores_for_10_new_donor':True,'weighted_mean_fixed_denominator':True,'missing_new_donor_no_partial_mean':True,'known_SEA_5_parameter_OLS_HC3_t_manual':{'true_beta':-.6,'beta':r['beta'],'SE':se,'df':25,'p':p},'design_gate_under20_rejected':True}
(ROOT/'logs/synthetic_PASS.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
