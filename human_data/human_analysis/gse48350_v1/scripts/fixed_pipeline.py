"""Fixed helpers for registered GSE48350 v1; no model selection or gene scans."""
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

SD_EPS=1e-8

def transform_probes(frame, *, scale, input_log2_count):
    if input_log2_count != 0:
        raise AssertionError('Input already logged: second log2 prohibited')
    a=frame.to_numpy(dtype=float)
    if not np.isfinite(a).all() or not (a>0).all():
        raise AssertionError('Input must be finite strictly positive source ratios')
    if scale=='log2_ratio':return pd.DataFrame(np.log2(a),index=frame.index,columns=frame.columns),1
    if scale=='raw_ratio':return frame.copy(),0
    raise ValueError(scale)

def aggregate_probes(probes, gene_to_probes):
    assert probes.index.is_unique
    # Transform is upstream: even-number probe medians are medians on that scale.
    return pd.DataFrame({g:probes.loc[ids].median(axis=0,skipna=False) for g,ids in gene_to_probes.items()}).T

def standardize_fixed_A(gene_values, reference_samples, fixed_A=None):
    ref=gene_values.loc[:,reference_samples]
    mu=ref.mean(axis=1,skipna=False);sd=ref.std(axis=1,ddof=1,skipna=False)
    good=np.isfinite(ref.to_numpy()).all(axis=1)&np.isfinite(sd.to_numpy())&(sd.to_numpy()>SD_EPS)
    A=list(gene_values.index[good]) if fixed_A is None else list(fixed_A)
    if fixed_A is not None and not set(A).issubset(set(gene_values.index[good])):
        return {'status':'not_testable_fixed_A_variance','A':A,'mu':mu,'sd':sd,'variance_pass':pd.Series(good,index=gene_values.index)}
    if not A:return {'status':'not_testable_score','A':[],'mu':mu,'sd':sd,'variance_pass':pd.Series(good,index=gene_values.index)}
    z=gene_values.loc[A].sub(mu.loc[A],axis=0).div(sd.loc[A],axis=0)
    S=z.mean(axis=0,skipna=False)  # Never sample-specific denominator.
    S_mu=float(S.loc[reference_samples].mean(skipna=False));S_sd=float(S.loc[reference_samples].std(ddof=1,skipna=False))
    ok=np.isfinite(S_sd) and S_sd>SD_EPS
    Y=(S-S_mu)/S_sd if ok else S*np.nan
    return {'status':'ready' if ok else 'not_testable_score','A':A,'mu':mu,'sd':sd,'variance_pass':pd.Series(good,index=gene_values.index),'z':z,'S':S,'Y':Y,'S_mean_reference':S_mu,'S_sd_reference':S_sd}

def coverage_gate(A_n,G0_n):
    minimum=max(10,int(np.ceil(.8*G0_n)))
    return A_n>=minimum,minimum

def design_gate(frame,terms):
    cols=['Y_score','age','sex',*terms]
    cols=list(dict.fromkeys(cols));complete=np.isfinite(frame[cols].to_numpy(dtype=float)).all(axis=1)
    d=frame.loc[complete].copy();X=pd.DataFrame({'Intercept':1.0},index=d.index)
    for t in terms:X[t]=d[t].astype(float)
    n,p=X.shape;rank=int(np.linalg.matrix_rank(X.to_numpy())) if n else 0
    counts=d.sex.value_counts().to_dict();age_unique=d.age.nunique();age_range=float(d.age.max()-d.age.min()) if n else np.nan
    reasons=[]
    if n<max(20,3*p,p+15):reasons.append('minimum_n')
    if rank!=p:reasons.append('rank_deficient')
    if age_unique<5 or not np.isfinite(age_range) or age_range<10:reasons.append('age_support')
    if counts.get(0,0)<3 or counts.get(1,0)<3:reasons.append('sex_support')
    return d,X,{'status':'not_testable_design' if reasons else 'ready','reason':';'.join(reasons),'n_complete':n,'n_missing':int((~complete).sum()),'p':p,'rank':rank,'df':n-rank,'minimum_n':max(20,3*p,p+15),'age_unique':age_unique,'age_min':float(d.age.min()) if n else np.nan,'age_max':float(d.age.max()) if n else np.nan,'age_range':age_range,'female_n':int(counts.get(0,0)),'male_n':int(counts.get(1,0))}

def fit_hc3(X,y,term):
    base=sm.OLS(np.asarray(y,float),X,hasconst=True).fit()
    inf=base.get_influence();hat=inf.hat_matrix_diag
    if not np.isfinite(hat).all() or (1-hat<=0).any():return {'status':'not_testable_HC3','reason':'undefined_leverage_denominator'},None,None
    fit=base.get_robustcov_results(cov_type='HC3',use_t=True)
    j=list(X.columns).index(term);n,p=X.shape;rank=np.linalg.matrix_rank(X.to_numpy());df=int(n-rank)
    beta=float(fit.params[j]);se=float(fit.bse[j]);critical=float(stats.t.ppf(.975,df))
    if not np.isfinite(se) or se<=0:return {'status':'not_testable_HC3','reason':'nonpositive_or_nonfinite_SE'},None,None
    statistic=beta/se;praw=float(2*stats.t.sf(abs(statistic),df))
    assert np.isclose(praw,fit.pvalues[j],rtol=1e-12,atol=1e-14)
    assert np.allclose([beta-critical*se,beta+critical*se],fit.conf_int(alpha=.05)[j],rtol=1e-12,atol=1e-14)
    age=X['age10'].to_numpy();other=X.drop(columns='age10').to_numpy();ar=age-other@np.linalg.lstsq(other,age,rcond=None)[0];sst=((age-age.mean())**2).sum();vif=float(sst/(ar@ar))
    cook=inf.cooks_distance[0]
    r={'status':'tested','reason':'','beta':beta,'SE_HC3':se,'df':df,'t_HC3':statistic,'CI95_low':beta-critical*se,'CI95_high':beta+critical*se,'raw_P':praw,'max_hat':float(hat.max()),'max_Cook':float(cook.max()),'age_VIF':vif,'hat_threshold':2*p/n,'Cook_threshold':4/n,'n_hat_flag':int((hat>2*p/n).sum()),'n_Cook_flag':int((cook>4/n).sum()),'age_VIF_flag':vif>=10}
    r['unstable_flag']=bool(r['n_hat_flag'] or r['n_Cook_flag'] or r['age_VIF_flag'])
    diag=pd.DataFrame({'fitted':base.fittedvalues,'residual':base.resid,'hat':hat,'Cook_distance':cook,'hat_flag':hat>2*p/n,'Cook_flag':cook>4/n},index=X.index)
    loo=[]
    for i in range(n):
        keep=np.arange(n)!=i;xx=X.to_numpy()[keep];yy=np.asarray(y,float)[keep]
        if np.linalg.matrix_rank(xx)!=p:loo.append((X.index[i],np.nan,'rank_deficient'));continue
        bb=np.linalg.lstsq(xx,yy,rcond=None)[0][j];loo.append((X.index[i],float(bb),'computed_beta_only'))
    loo=pd.DataFrame(loo,columns=['omitted_sample','beta','status']);finite=loo.beta.dropna()
    r.update(LOO_beta_min=float(finite.min()),LOO_beta_max=float(finite.max()),LOO_direction_flip_n=int(((finite*beta)<0).sum()),LOO_not_estimable_n=int(loo.beta.isna().sum()))
    return r,diag,loo
