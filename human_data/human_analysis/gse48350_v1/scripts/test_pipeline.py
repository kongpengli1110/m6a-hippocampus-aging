"""Run before human models: deterministic known OLS, independent HC3, score/gate tests."""

# Release-only safety/path adapter; see provenance/adaptation_diffs.
from pathlib import Path as _ReleasePath
import sys as _release_sys
_release_home = next(p for p in _ReleasePath(__file__).resolve().parents if (p / 'release_support.py').is_file())
_release_sys.path.insert(0, str(_release_home))
from release_support import gate, release_root, resolve_source
gate(__file__, __name__)
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy import stats
from fixed_pipeline import transform_probes,aggregate_probes,standardize_fixed_A,coverage_gate,design_gate,fit_hc3
ROOT=Path(__file__).resolve().parents[1]
def run_tests():
    checks={}
    n=40;age=np.linspace(20,99,n);sex=np.arange(n)%2;X=pd.DataFrame({'Intercept':np.ones(n),'age10':(age-80)/10,'sex':sex});truth=np.array([.3,-.5,.7]);xx=X.to_numpy();bread=np.linalg.inv(xx.T@xx)
    noise=np.sin(np.arange(n)*1.31)*(.4+np.arange(n)/n);error=noise-xx@(bread@xx.T@noise);y=xx@truth+error
    r,diag,loo=fit_hc3(X,y,'age10');beta=bread@xx.T@y;e=y-xx@beta;hat=np.einsum('ij,jk,ik->i',xx,bread,xx);hc3=bread@(xx.T@np.diag(e**2/(1-hat)**2)@xx)@bread;se=np.sqrt(hc3[1,1]);t=beta[1]/se;pv=2*stats.t.sf(abs(t),37);ci=beta[1]+np.array([-1,1])*stats.t.ppf(.975,37)*se
    assert np.allclose(beta,truth,atol=1e-12);assert np.allclose([r['beta'],r['SE_HC3'],r['raw_P']],[truth[1],se,pv],rtol=1e-11,atol=1e-13);assert np.allclose([r['CI95_low'],r['CI95_high']],ci,atol=1e-12);assert r['df']==37
    checks['known_synthetic_OLS_HC3_t_df']={'PASS':True,'true_age_beta':-.5,'estimated_age_beta':r['beta'],'SE_HC3':se,'df':37,'manual_P':pv,'manual_CI':ci.tolist()}
    probes=pd.DataFrame([[1,4,9],[9,16,25],[2,8,18]],index=['p1','p2','p3'],columns=['r1','r2','new'],dtype=float)
    logged,flag=transform_probes(probes,scale='log2_ratio',input_log2_count=0);assert flag==1
    try:transform_probes(logged,scale='log2_ratio',input_log2_count=flag)
    except AssertionError:checks['double_log_rejected']=True
    else:raise AssertionError('Double log accepted')
    for bad in [0,-1,np.nan,np.inf]:
        z=probes.copy();z.iloc[0,0]=bad
        try:transform_probes(z,scale='log2_ratio',input_log2_count=0)
        except AssertionError:pass
        else:raise AssertionError('Invalid ratio accepted')
    g=aggregate_probes(logged,{'g1':['p1','p2'],'g2':['p3']});assert np.isclose(g.loc['g1','r1'],np.log2(3));assert not np.isclose(g.loc['g1','r1'],np.log2(5))
    z=standardize_fixed_A(g,['r1','r2']);assert np.isclose(z['Y'][['r1','r2']].mean(),0,atol=1e-12);assert np.isclose(z['Y'][['r1','r2']].std(ddof=1),1,atol=1e-12)
    g.loc['g1','new']=np.nan;zz=standardize_fixed_A(g,['r1','r2'],fixed_A=z['A']);assert pd.isna(zz['Y']['new']);assert zz['A']==z['A']
    raw,rf=transform_probes(probes,scale='raw_ratio',input_log2_count=0);assert rf==0;assert raw.equals(probes)
    assert coverage_gate(18,22)==(True,18);assert coverage_gate(17,22)==(False,18);assert coverage_gate(95,120)==(False,96);assert coverage_gate(96,120)==(True,96)
    d=pd.DataFrame({'Y_score':y,'age':age,'age10':X.age10,'sex':sex});_,_,gate=design_gate(d,['age10','sex']);assert gate['status']=='ready'
    for label,bad in [('small_n',d.iloc[:10]),('one_sex',d.assign(sex=0)),('narrow_age',d.assign(age=80+np.linspace(0,2,n),age10=np.linspace(0,.2,n))),('rank_deficient',d.assign(sex=d.age10))]:
        _,_,gg=design_gate(bad,['age10','sex']);assert gg['status']=='not_testable_design',label
    checks.update(log_then_probe_median=True,normal_reference_z_and_score_SD=True,missing_nonreference_no_variable_denominator=True,raw_ratio_unchanged=True,coverage_boundaries=True,design_gates=True,LOO_beta_only=len(loo)==40)
    from datetime import datetime,timezone
    out={'status':'PASS','completed_at_utc':datetime.now(timezone.utc).isoformat(),'checks':checks}
    (ROOT/'logs/synthetic_validation.json').write_text(json.dumps(out,indent=2));return out
if __name__=='__main__':print(json.dumps(run_tests(),indent=2))
