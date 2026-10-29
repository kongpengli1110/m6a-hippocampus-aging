"""Reusable locked-H2 mathematics. Importing this module does not read real data."""
from dataclasses import dataclass
from collections import Counter
import math
import numpy as np
import pandas as pd
from scipy import sparse

class NotTestable(RuntimeError):
    """A technical/statistical gate failed; never silently switch the matching rule."""


def aggregate_reference(counts, samples, class_order, min_nuclei=20, min_donors=20):
    """Rows are existing donor-supertype pseudobulks; sum them, never mean their means.
    Library totals include ALL input features, not just tested/background columns.
    """
    counts=sparse.csr_matrix(counts)
    required={'donor','Class','n_nuclei'}
    if not required<=set(samples): raise ValueError('Missing sample fields')
    if counts.shape[0]!=len(samples): raise ValueError('Sample/matrix shape mismatch')
    if not np.isfinite(counts.data).all() or (counts.data<0).any(): raise ValueError('Invalid counts')
    if not samples.Class.isin(class_order).all(): raise ValueError('Unknown official class')
    keys=samples[['donor','Class']].drop_duplicates().sort_values(['Class','donor']).reset_index(drop=True)
    keymap={(r.donor,r.Class):i for i,r in keys.iterrows()}
    ix=np.asarray([keymap[(r.donor,r.Class)] for r in samples.itertuples()])
    selector=sparse.csr_matrix((np.ones(len(samples),dtype=np.int64),(ix,np.arange(len(samples)))),shape=(len(keys),len(samples)))
    grouped=selector@counts
    lib=np.asarray(grouped.sum(axis=1)).ravel()
    nuc=np.asarray(selector@samples.n_nuclei.to_numpy()).ravel()
    keys['n_nuclei']=nuc;keys['totalUMI']=lib;keys['eligible']=(nuc>=min_nuclei)&(lib>0)
    means=[];detect=[];n_by_class={};donors={};cpms={}
    for cls in class_order:
        mask=((keys.Class==cls)&keys.eligible).to_numpy()
        n=int(mask.sum())
        if n<min_donors: raise NotTestable(f'class_insufficient_donors:{cls}:{n}')
        cp=grouped[mask].astype(np.float64).multiply((1e6/lib[mask])[:,None]).tocsr()
        means.append(np.asarray(cp.mean(axis=0)).ravel())
        detect.append(np.asarray((cp>=1).sum(axis=0)).ravel()/n)
        n_by_class[cls]=n;donors[cls]=keys.loc[mask,'donor'].tolist();cpms[cls]=cp
    return {'means':np.stack(means),'detection':np.stack(detect),'manifest':keys,'n_by_class':n_by_class,'donors':donors,'grouped_counts':grouped,'cpm':cpms}


def specificity_shares(class_means):
    x=np.asarray(class_means,dtype=float)
    if x.ndim!=2 or not np.isfinite(x).all() or (x<0).any(): raise ValueError('Invalid class means')
    total=x.sum(axis=0);valid=total>0
    shares=np.full(x.shape,np.nan)
    shares[:,valid]=x[:,valid]/total[valid]
    if not np.allclose(shares[:,valid].sum(axis=0),1,rtol=0,atol=2e-15): raise AssertionError('Shares do not sum to one')
    return shares,valid


@dataclass
class SamplingPlan:
    """Fixed slot order: lexicographic bin tuple, then HGNC ID within bin.
    Pools contain only background noncandidate genes. Gene positions reference bg order.
    """
    bin_columns: tuple
    candidate_hgnc: tuple
    weights: np.ndarray
    slot_bins: tuple
    candidate_positions: np.ndarray
    bins: tuple
    background_hgnc: tuple
    background_ensembl: tuple
    g0_count: int
    g0_weight: float

    @property
    def weight_sum(self): return float(self.weights.sum())
    @property
    def n_slots(self): return len(self.weights)


def build_sampling_plan(bg, members, bin_columns, candidate_union, min_count=None):
    if bg.HGNC_ID.duplicated().any() or members.HGNC_ID.duplicated().any(): raise ValueError('Duplicate gene identity')
    if not np.isfinite(members.weight).all() or (members.weight<=0).any(): raise ValueError('Invalid weights')
    m=members.copy();g0_count=len(m);g0_weight=float(m.weight.sum())
    bg=bg.reset_index(drop=True).copy();pos={g:i for i,g in enumerate(bg.HGNC_ID)}
    retained=m[m.HGNC_ID.isin(pos)].copy()
    if min_count is None: min_count=max(10,math.ceil(.8*g0_count))
    if len(retained)<min_count or retained.weight.sum()/g0_weight<.8: raise NotTestable('coverage_gate_failed')
    actual_controls=~bg.HGNC_ID.isin(candidate_union)
    if not np.array_equal(actual_controls.to_numpy(),bg.noncandidate_pool.to_numpy(dtype=bool)):
        raise ValueError('Frozen noncandidate pool inconsistent with fixed all-set union')
    for col in bin_columns:
        if bg[col].isna().any(): raise ValueError('Missing frozen bin')
        retained[col]=[int(bg.loc[pos[g],col]) for g in retained.HGNC_ID]
    retained=retained.sort_values(list(bin_columns)+['HGNC_ID'],kind='stable')
    weights=retained.weight.to_numpy(dtype=float);weights.flags.writeable=False
    candidate_positions=np.asarray([pos[g] for g in retained.HGNC_ID],dtype=int)
    slot_bins=tuple(tuple(int(row[c]) for c in bin_columns) for _,row in retained.iterrows())
    bins=[]
    for key in sorted(set(slot_bins)):
        slots=np.asarray([i for i,k in enumerate(slot_bins) if k==key],dtype=int)
        mask=actual_controls.to_numpy().copy()
        for col,val in zip(bin_columns,key):mask&=(bg[col].to_numpy()==val)
        pool=np.flatnonzero(mask)
        if len(pool)<len(slots):
            raise NotTestable(f'insufficient_pool:bin={key}:pool={len(pool)}:slots={len(slots)}')
        if set(bg.loc[pool,'HGNC_ID'])&candidate_union:raise AssertionError('Candidate in control pool')
        bins.append((key,slots,pool))
    return SamplingPlan(tuple(bin_columns),tuple(retained.HGNC_ID),weights,slot_bins,candidate_positions,tuple(bins),tuple(bg.HGNC_ID),tuple(bg.ensembl_id),g0_count,g0_weight)


def draw_positions(plan,rng):
    """One null. All occupied bins sample without replacement, then random slot assignment."""
    assigned=np.empty(plan.n_slots,dtype=int)
    for key,slots,pool in plan.bins:
        sampled=rng.choice(pool,size=len(slots),replace=False)
        assigned[slots]=sampled[rng.permutation(len(slots))]
    if len(np.unique(assigned))!=plan.n_slots:raise AssertionError('Duplicate gene in null')
    return assigned


def weighted_score(values,weights):
    values=np.asarray(values,dtype=float);weights=np.asarray(weights,dtype=float)
    if values.shape!=weights.shape or len(values)==0 or not np.isfinite(values).all():raise ValueError('Invalid score input')
    if not np.isfinite(weights).all() or (weights<=0).any():raise ValueError('Invalid score weights')
    return float(np.dot(values,weights)/weights.sum())


def generate_null(plan,background_shares,rng,n_draws=10000,audit_draws=(0,1,17,9999)):
    """No gene selection/bin reconstruction occurs here. Seed is supplied by fixed child stream."""
    x=np.asarray(background_shares,dtype=float)
    if len(x)!=len(plan.background_hgnc) or not np.isfinite(x).all():raise ValueError('Invalid background shares')
    out=np.empty(n_draws,dtype=float);audits=[];selected_audits=set(audit_draws)
    original_weights=plan.weights.copy()
    for b in range(n_draws):
        ix=draw_positions(plan,rng)
        out[b]=weighted_score(x[ix],plan.weights)
        if b in selected_audits:
            independent=math.fsum(float(w)*float(s) for w,s in zip(plan.weights,x[ix]))/math.fsum(map(float,plan.weights))
            if not math.isclose(out[b],independent,rel_tol=1e-13,abs_tol=1e-15):raise AssertionError('Independent null score check failed')
            for slot,idx in enumerate(ix):
                audits.append({'draw_0based':b,'slot':slot,'candidate_slot_HGNC':plan.candidate_hgnc[slot],
                 'slot_weight':float(plan.weights[slot]),'bin':':'.join(map(str,plan.slot_bins[slot])),
                 'control_HGNC':plan.background_hgnc[idx],'control_Ensembl':plan.background_ensembl[idx],
                 'control_position':int(idx),'control_Glut_share':float(x[idx]),'null_score':float(out[b])})
    if not np.array_equal(original_weights,plan.weights):raise AssertionError('Slot weights changed')
    return out,audits


def empirical_summary(observed,null):
    x=np.asarray(null,dtype=float)
    if len(x)==0 or not np.isfinite(x).all() or not np.isfinite(observed):raise ValueError('Invalid empirical-test input')
    n=len(x);lower=(1+np.count_nonzero(x<=observed))/(n+1);upper=(1+np.count_nonzero(x>=observed))/(n+1)
    p=min(1.,2*min(lower,upper));center=float(np.median(x))
    lo,hi=np.quantile(x,[.025,.975],method='linear')
    # Independent ordered-search implementation, including ties.
    ordered=np.sort(x)
    lower2=(1+np.searchsorted(ordered,observed,side='right'))/(n+1)
    upper2=(1+n-np.searchsorted(ordered,observed,side='left'))/(n+1)
    if lower!=lower2 or upper!=upper2:raise AssertionError('Independent empirical tails differ')
    return dict(Tobs=float(observed),null_median=center,Delta=float(observed-center),
                null_reference_lower_025=float(lo),null_reference_upper_975=float(hi),
                p_lower=float(lower),p_upper=float(upper),raw_p=float(p),
                null_draws=n,tail_probability_resolution=1/(n+1),minimum_two_sided_p=2/(n+1),
                null_reference_interval_type='central_95_percent_matched_null_reference_NOT_effect_CI')


def bh_fixed(pvalues,m=8):
    """Unavailable tests retain their positions with internal p=1; public raw p stays NA."""
    p=np.asarray(pvalues,dtype=float)
    if len(p)!=m:raise ValueError(f'Family must retain exactly {m} registered positions')
    if np.any(np.isfinite(p)&((p<0)|(p>1))):raise ValueError('Invalid p')
    work=np.where(np.isfinite(p),p,1.0);order=np.argsort(work,kind='stable')
    ranked=work[order]*m/np.arange(1,m+1)
    adj=np.minimum.accumulate(ranked[::-1])[::-1].clip(0,1)
    out=np.empty(m);out[order]=adj
    return out
