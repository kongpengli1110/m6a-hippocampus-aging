#!/usr/bin/env python3
"""Single-command locked H2 runner. NOT run during code preparation.
Explicit --execute-locked-h2 is required. No age/cognition/bootstrap analyses.
"""

# Release-only safety/path adapter; see provenance/adaptation_diffs.
from pathlib import Path as _ReleasePath
import sys as _release_sys
_release_home = next(p for p in _ReleasePath(__file__).resolve().parents if (p / 'release_support.py').is_file())
_release_sys.path.insert(0, str(_release_home))
from release_support import gate, release_root, resolve_source
gate(__file__, __name__)
from pathlib import Path
from datetime import datetime,timezone
import argparse,csv,fcntl,hashlib,importlib.metadata,io,json,math,os,platform,shutil,subprocess,sys,traceback
import numpy as np
import pandas as pd
from scipy import sparse
from statsmodels.stats.multitest import multipletests
from h2_core import aggregate_reference,specificity_shares,build_sampling_plan,generate_null,weighted_score,empirical_summary,bh_fixed,NotTestable

ROOT=release_root()
DESIGN=ROOT/'H2_design'; D=DESIGN/'derived'; SEA=ROOT/'human_data/sea_hip'
OUT=ROOT/'human_analysis/h2_v1'
FROZEN_MANIFEST_SHA256='e16f9662914da57f99e63293fd5605b8e2f442ab462c5effa7ab0dc069faa3b6'
B=10000;BASE_SEED=20261012
CLASS_ORDER=['Neuronal: Glutamatergic','Neuronal: GABAergic','Non-neuronal and Non-neural']
SET_ORDER=['core23','CA3-M2','DG-M3','DG-M7','CA3-M4','CA1-M9','SUB-M4','SUB-M6']
AUDIT_DRAWS=(0,1,17,9999)

def now():return datetime.now(timezone.utc).isoformat()
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()
def json_write(path,obj):
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n');temp.replace(path)
def read_tsv(path):return pd.read_csv(path,sep='\t')

def verify_locked_contract():
    mf=DESIGN/'freeze_manifest.json'
    if sha(mf)!=FROZEN_MANIFEST_SHA256:raise RuntimeError('Frozen H2 manifest SHA changed; stop, do not adapt rules')
    m=json.loads(mf.read_text())
    if m['fixed_set_order']!=SET_ORDER or m['reference_group_labels']!=CLASS_ORDER:raise RuntimeError('Frozen order mismatch')
    if m['null_generation_planned_not_run']['B']!=B or m['null_generation_planned_not_run']['seed']!=BASE_SEED:raise RuntimeError('Frozen RNG specification mismatch')
    checked={}
    for rec in m['locked_and_supporting_files']+m['upstream_input_hashes']:
        path=resolve_source(rec['path']);expected=rec['sha256']
        if str(path) in checked:
            if checked[str(path)]!=expected:raise RuntimeError('Inconsistent duplicate frozen checksum')
            continue
        actual=sha(path)
        if actual!=expected:raise RuntimeError(f'Frozen input altered: {path}')
        checked[str(path)]=actual
    cand=resolve_source(m['candidate_xlsx']['path'])
    if hashlib.md5(cand.read_bytes()).hexdigest()!=m['candidate_xlsx']['md5']:raise RuntimeError('Candidate freeze changed')
    bg=read_tsv(D/'background_matching_manifest.tsv');bl=read_tsv(D/'length_background_matching_manifest.tsv')
    member=read_tsv(D/'candidate_technical_membership.tsv');ml=read_tsv(D/'length_candidate_technical_membership.tsv')
    if len(bg)!=18589 or len(bl)!=18589 or bg.noncandidate_pool.sum()!=17471 or bl.noncandidate_pool.sum()!=17471:raise RuntimeError('Frozen background totals changed')
    if len(member)!=1276 or len(ml)!=1276:raise RuntimeError('Frozen memberships changed')
    if not np.array_equal(bg.HGNC_ID,bl.HGNC_ID) or not np.array_equal(bg.gene_index,bl.gene_index):raise RuntimeError('Length identity/order mismatch')
    union=set(member.HGNC_ID)
    for table in [bg,bl]:
        if not np.array_equal(table.noncandidate_pool,~table.HGNC_ID.isin(union)):raise RuntimeError('Control pool includes candidate union')
    for ids,family in [([f'H2_{i:02}' for i in range(1,9)],'H2_BASE'),([f'H2L_{i:02}' for i in range(1,9)],'H2_LENGTH')]:
        actual=[r for r in m['registered_families'] if r['family']==family]
        if len(actual)!=1 or actual[0]['m']!=8 or actual[0]['ids']!=ids:raise RuntimeError('Registered family changed')
    # Check exact frozen capacities and memberships without any scores/RNG.
    plans={}
    for length,table,capfile in [(False,bg,'matching_bin_capacity.tsv'),(True,bl,'length_matching_bin_capacity.tsv')]:
        cols=('expression_bin','detection_bin','length_bin') if length else ('expression_bin','detection_bin')
        caps=read_tsv(D/capfile)
        if len(caps)!=(2000 if length else 400) or not caps.capacity_pass.all():raise RuntimeError('Frozen capacity gates not passed')
        for name in SET_ORDER:
            mm=member[member.set_name==name].copy()
            plan=build_sampling_plan(table,mm,cols,union)
            eligible_flag=mm.background_eligible
            if set(mm.loc[eligible_flag,'HGNC_ID'])!=set(plan.candidate_hgnc):raise RuntimeError('Runtime candidate membership differs from freeze')
            lm=ml[ml.set_name==name].set_index('HGNC_ID')
            for h,w in zip(plan.candidate_hgnc,plan.weights):
                if not bool(lm.loc[h,'length_background_eligible']) or float(lm.loc[h,'weight'])!=float(w):raise RuntimeError('Length membership/weight mismatch')
            for _,row in caps[caps.set_name==name].iterrows():
                key=tuple(int(row[c]) for c in cols)
                pool=table.noncandidate_pool.to_numpy().copy()
                for col,value in zip(cols,key):pool&=table[col].to_numpy()==value
                slots=sum(k==key for k in plan.slot_bins)
                if int(pool.sum())!=int(row.noncandidate_pool_n) or slots!=int(row.required_candidate_slots):raise RuntimeError('Frozen bin capacity mismatch')
            plans[('H2_LENGTH' if length else 'H2_BASE',name)]=plan
    return m,checked,bg,bl,member,ml,plans


def plot_results(results,prefix,synthetic=False):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    fig,axes=plt.subplots(1,2,figsize=(13.5,6.7),sharey=True)
    for ax,family,title in zip(axes,['H2_BASE','H2_LENGTH'],['Expression + detection matching','Expression + detection + exon-union length']):
        rr=results[results.family==family].set_index('set_name').reindex(SET_ORDER)
        for i,(name,row) in enumerate(rr.iterrows()):
            if row.status!='EXECUTED':
                ax.text(.03,i,'Not testable',va='center',color='gray');continue
            ax.plot([row.null_reference_lower_025,row.null_reference_upper_975],[i,i],color='#9aa3ad',lw=5,alpha=.7,solid_capstyle='butt')
            ax.scatter(row.null_median,i,marker='|',s=110,color='#333333',zorder=3)
            ax.scatter(row.Tobs,i,s=42,color='#1b6b93',zorder=4)
            ax.text(1.02,i,f"p={row.raw_p:.3g}; q={row.q_BH:.3g}",transform=ax.get_yaxis_transform(),va='center',fontsize=8.5)
        ax.set_xlim(0,1);ax.set_ylim(7.5,-.6);ax.set_yticks(range(8));ax.set_yticklabels(SET_ORDER)
        ax.set_xlabel('Weighted mean glutamatergic expression share')
        ax.set_title(title,fontsize=10.5);ax.grid(axis='x',alpha=.2)
        ax.spines[['top','right']].set_visible(False)
    fig.suptitle('SYNTHETIC QA ONLY — NOT HUMAN RESULTS' if synthetic else 'Human expression-reference enrichment of fixed mouse-derived gene sets',fontsize=12)
    legend=[Line2D([0],[0],marker='o',color='none',markerfacecolor='#1b6b93',markeredgecolor='#1b6b93',label='Observed set score'),Line2D([0],[0],marker='|',color='none',markeredgecolor='#333333',markersize=11,label='Matched-null median'),Line2D([0],[0],color='#9aa3ad',lw=5,label='Central 95% null reference band (NOT effect CI)')]
    fig.legend(handles=legend,loc='lower center',ncol=3,fontsize=9,bbox_to_anchor=(.5,.065))
    fig.text(.5,.025,'Three official broad classes; 40 equally weighted donors/class; 10,000 null draws/test. Separate BH families (m=8 each).\nExpression evidence only: not a measurement of m6A, enzyme activity, causality, or rescue of the age results.',ha='center',fontsize=8.5)
    fig.subplots_adjust(left=.08,right=.86,bottom=.19,top=.89,wspace=.72)
    fig.savefig(str(prefix)+'.pdf',bbox_inches='tight');fig.savefig(str(prefix)+'.png',dpi=220,bbox_inches='tight');plt.close(fig)


def write_report(results,qa,path):
    lines=['# Locked H2 results — expression-reference enrichment','',
      'H1/H3 results were already known when this H2 appendix was locked: all three primary Holm tests did not meet the threshold. H2 was not yet computed at that lock. These analyses do not rescue or replace the age tests.','',
      'Official Glut/GABA/non-neuronal-and-non-neural classes; 40 donors per class, equal donor means. Fixed stable single-feature background: 18,589 genes; noncandidate pool: 17,471. All 117 multifeature HGNC identities remain excluded from H2 projection, but retained in full-feature library denominators.','',
      '| Family | Set | Status | Tobs | Null median | Delta | Raw two-sided p | BH q (m=8) |','|---|---|---|---:|---:|---:|---:|---:|']
    def val(x):return 'NA' if pd.isna(x) else f'{x:.6g}'
    for r in results.itertuples():lines.append(f'| {r.family} | {r.set_name} | {r.status} | {val(r.Tobs)} | {val(r.null_median)} | {val(r.Delta)} | {val(r.raw_p)} | {val(r.q_BH)} |')
    lines+=['','## Interpretation and reporting boundaries']
    for fam in ['H2_BASE','H2_LENGTH']:
        r=results[(results.family==fam)&(results.status=='EXECUTED')]
        pos=r[(r.Delta>0)&(r.q_BH<.05)].set_name.tolist();neg=r[(r.Delta<0)&(r.q_BH<.05)].set_name.tolist()
        lines.append(f'- {fam}: positive matched-background enrichment meeting q<0.05: '+(', '.join(pos) if pos else 'none')+'. Negative deviation meeting q<0.05: '+(', '.join(neg) if neg else 'none')+'.')
    lines+=['- A significant result supports a relative glutamatergic expression share compared with the fixed expression/detection-matched gene background (and, separately, the registered length sensitivity). It does not establish single-cell m6A changes, METTL3 enzyme activity, a causal effect, or aging replication.',
      '- These three broad classes do not locate the signal specifically to CA3 or SUB. A significant set and a nonsignificant set are not a formal between-set or between-cell-type comparison. The same SEA donors are used; this is not an independent cohort replication.',
      '- Tobs/null median/Delta and both null reference limits are reported in h2_results_long.tsv for every registered test. Null reference bands are NOT confidence intervals for an observed enrichment effect. No bootstrap or observed-effect CI was run.',
      '- Exactly 10,000 null draws per executed test. Each tail uses inclusive ties and a +1 correction; tail grid is 1/10001, and the minimum two-sided empirical p is 2/10001 ≈ 0.000199980. No extra draws were added to pursue smaller p. These are finite-Monte-Carlo matched-background comparisons, not exact causal randomization tests.',
      '- Two separately registered families use BH with m=8 each. Unavailable positions, if any, retain internal p=1; their public raw p and q remain NA. The length family is not appended to the original E_SENS23 family.',
      '- Cognition and H4 remain pending; no age, cognitive or finer-cell-type tests were performed. No existing combined_v1 or age results were overwritten. No combined tracking copy is produced by this minimal runner.',
      '', '## QA and audit',f"- Final QA status: {qa['status']}; details in QA.json and QA_report.md.",'- Full 16×10000 null matrix and per-test distributions are saved; complete control gene IDs, bins and fixed slot weights for draws 0, 1, 17 and 9999 are in random_draw_audit.tsv.','- Input/freeze hashes, software versions, source code hashes, full memberships/background expression audit and output hashes are recorded.']
    path.write_text('\n'.join(lines)+'\n')


def execute():
    OUT.mkdir(exist_ok=True)
    lockfile=open(OUT/'execution.lock','a+')
    try:fcntl.flock(lockfile.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise RuntimeError('Another H2 runner is active; no duplicate worker')
    statusfile=OUT/'run_status.json'
    if statusfile.exists():raise RuntimeError('A prior H2 run exists; preserve it and do not overwrite/re-run silently')
    started=now();json_write(statusfile,{'status':'RUNNING','started_utc':started,'completed_registered_tests':0})
    for d in ['qa','figures','null_distributions','frozen_input_copies']: (OUT/d).mkdir(exist_ok=True)
    try:
        with (OUT/'qa/synthetic_tests_at_execution.log').open('w') as log:
            subprocess.run([sys.executable,str(OUT/'tests/test_h2_core.py')],stdout=log,stderr=subprocess.STDOUT,check=True)
        print('Synthetic tests passed; verifying frozen hashes and exact bins.',flush=True)
        m,checked,bg,bl,member,ml,plans=verify_locked_contract()
        software={'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'packages':{k:importlib.metadata.version(k) for k in ['numpy','scipy','pandas','matplotlib','statsmodels']}}
        json_write(OUT/'software.json',software)
        for fn in ['background_matching_manifest.tsv','length_background_matching_manifest.tsv','candidate_technical_membership.tsv','length_candidate_technical_membership.tsv','matching_bin_capacity.tsv','length_matching_bin_capacity.tsv','official_taxonomy_map.tsv','donor_class_technical_manifest.tsv','multifeature_identity_pending.tsv','source_gene_exon_union_lengths.tsv']:
            shutil.copy2(D/fn,OUT/'frozen_input_copies'/fn)
        json_write(OUT/'execution_manifest.json',{'started_utc':started,'authorization':'Owner authorized locked H2 execution; parent initiates this one explicit --execute-locked-h2 command','h2_freeze_manifest_sha256':FROZEN_MANIFEST_SHA256,'input_hashes':checked,'source_code_hashes':{str(p):sha(p) for p in (OUT/'scripts').glob('*.py')},'test_code_hashes':{str(p):sha(p) for p in (OUT/'tests').glob('*.py')},'B':B,'seed':BASE_SEED,'class_order':CLASS_ORDER,'set_order':SET_ORDER,'families':['H2_BASE','H2_LENGTH'],'family_m':[8,8],'not_run':['bootstrap','age','cognition','H4','subclass enrichment']})
        samples=read_tsv(SEA/'pseudobulk_samples.tsv');tax=read_tsv(D/'official_taxonomy_map.tsv');features=read_tsv(SEA/'features.tsv')
        samples=samples.merge(tax,left_on='supertype',right_on='Supertype',how='left',validate='many_to_one')
        if not (samples.subclass==samples.Subclass).all() or not np.array_equal(samples.pseudobulk_index,np.arange(3402)):raise RuntimeError('Sample taxonomy/order mismatch')
        counts=sparse.load_npz(SEA/'counts.npz')
        if counts.shape!=(3402,36601) or not np.array_equal(np.asarray(counts.sum(axis=1)).ravel(),samples.totalUMI):raise RuntimeError('Input count rows/totals mismatch')
        reference=aggregate_reference(counts,samples,CLASS_ORDER,20,20)
        if any(n!=40 for n in reference['n_by_class'].values()) or any(set(reference['donors'][c])!=set(reference['donors'][CLASS_ORDER[0]]) for c in CLASS_ORDER):raise RuntimeError('Expected common 40 donors not retained')
        expected=read_tsv(D/'donor_class_technical_manifest.tsv')
        pd.testing.assert_frame_equal(reference['manifest'].reset_index(drop=True),expected.reset_index(drop=True),check_dtype=False)
        reference['manifest'].to_csv(OUT/'donor_class_audit.tsv',sep='\t',index=False)
        ix=bg.gene_index.to_numpy(dtype=int)
        if not np.array_equal(features.loc[ix,'gene_ids'],bg.ensembl_id):raise RuntimeError('Background gene-index identity mismatch')
        means=reference['means'];det=reference['detection']
        np.testing.assert_allclose(means[:,ix].mean(axis=0),bg.mean_CPM_class_balanced,rtol=2e-12,atol=1e-10)
        np.testing.assert_allclose(det[:,ix].mean(axis=0),bg.detection_class_balanced,rtol=0,atol=5e-16)
        if not (np.max(det[:,ix],axis=0)>=.2).all():raise RuntimeError('Frozen background expression gate no longer valid')
        shares,valid=specificity_shares(means)
        if not valid[ix].all():raise RuntimeError('A fixed background gene has zero share denominator; do not silently drop')
        # Independent direct donor loops for deterministic source gene subset (includes out-of-background features).
        qa_positions=np.unique(np.concatenate([np.linspace(0,36600,31,dtype=int),ix[:5],ix[-5:]]))
        direct=[]
        for ci,cls in enumerate(CLASS_ORDER):
            donor_vectors=[]
            for donor in sorted(reference['donors'][cls]):
                rows=np.flatnonzero(((samples.donor==donor)&(samples.Class==cls)).to_numpy())
                cnt=np.asarray(counts[rows].sum(axis=0)).ravel()
                denominator=float(cnt.sum())
                donor_vectors.append(cnt[qa_positions]*1e6/denominator)
            independent=np.stack(donor_vectors).mean(axis=0)
            np.testing.assert_allclose(independent,means[ci,qa_positions],rtol=2e-12,atol=1e-10)
            direct.append({'class':cls,'donor_count':len(donor_vectors),'source_genes_independently_checked':len(qa_positions),'max_abs_difference':float(np.max(np.abs(independent-means[ci,qa_positions])))})
        audit_bg=bg.copy()
        for ci,short in enumerate(['Glut','GABA','NN']):
            audit_bg[f'mean_CPM_{short}']=means[ci,ix];audit_bg[f'share_{short}']=shares[ci,ix];audit_bg[f'detection_{short}']=det[ci,ix]
        audit_bg.to_csv(OUT/'background_expression_audit.tsv',sep='\t',index=False)
        # Same fixed gene identities/order for base and length; keep explicit lookup rather than assuming it.
        x_base=shares[0,bg.gene_index.to_numpy(dtype=int)];x_length=shares[0,bl.gene_index.to_numpy(dtype=int)]
        codebook=audit_bg.set_index('HGNC_ID')
        membership_audit=member.copy()
        membership_audit['Glut_share']=[float(codebook.loc[h,'share_Glut']) if h in codebook.index else np.nan for h in member.HGNC_ID]
        membership_audit['weighted_Glut_share']=membership_audit.weight*membership_audit.Glut_share
        membership_audit.to_csv(OUT/'candidate_membership_audit.tsv',sep='\t',index=False)
        del counts,reference
        children=np.random.SeedSequence(BASE_SEED).spawn(16)
        rows=[];all_audits=[];matrix=np.full((16,B),np.nan);null_ids=[]
        for family_i,(family,background_x,background) in enumerate([('H2_BASE',x_base,bg),('H2_LENGTH',x_length,bl)]):
            for set_i,name in enumerate(SET_ORDER):
                stream=8*family_i+set_i;test_id=f'H2L_{set_i+1:02}' if family_i else f'H2_{set_i+1:02}';null_ids.append(test_id)
                plan=plans[(family,name)]
                row={'test_id':test_id,'family':family,'family_m':8,'set_name':name,'status':'EXECUTED','reason':'','set_G0_1to1_n':plan.g0_count,'projected_member_n':plan.n_slots,'number_coverage':plan.n_slots/plan.g0_count,'weight_coverage':plan.weight_sum/plan.g0_weight,'projected_weight_sum':plan.weight_sum,'background_n':len(background),'noncandidate_pool_n':int(background.noncandidate_pool.sum()),'donors_per_class':40,'seed_entropy':BASE_SEED,'seed_child_index':stream,'seed_spawn_key':':'.join(map(str,children[stream].spawn_key))}
                try:
                    observed=weighted_score(background_x[plan.candidate_positions],plan.weights)
                    independent=math.fsum(float(w)*float(background_x[j]) for w,j in zip(plan.weights,plan.candidate_positions))/math.fsum(map(float,plan.weights))
                    if not math.isclose(observed,independent,rel_tol=1e-13,abs_tol=1e-15):raise AssertionError('Independent observed score failed')
                    # Audit drawing is part of the same original stream; no redraw/seed reset.
                    null,audits=generate_null(plan,background_x,np.random.default_rng(children[stream]),B,AUDIT_DRAWS)
                    matrix[stream]=null;np.save(OUT/'null_distributions'/f'{test_id}.npy',null,allow_pickle=False)
                    for a in audits:a.update(test_id=test_id,family=family,set_name=name)
                    all_audits.extend(audits);row.update(empirical_summary(observed,null))
                    # Deterministic independent audit against frozen full identities/bins/weights.
                    for draw in AUDIT_DRAWS:
                        aa=[a for a in audits if a['draw_0based']==draw]
                        if len(aa)!=plan.n_slots or len({a['control_HGNC'] for a in aa})!=plan.n_slots:raise AssertionError('Audit draw size/uniqueness failed')
                        if set(a['control_HGNC'] for a in aa)&set(member.HGNC_ID):raise AssertionError('Candidate appears as control')
                        for a in aa:
                            j=a['control_position'];slot=a['slot'];actual_bin=':'.join(str(int(background.iloc[j][c])) for c in plan.bin_columns)
                            if actual_bin!=a['bin'] or a['slot_weight']!=float(plan.weights[slot]):raise AssertionError('Audit bin/weight mismatch')
                except NotTestable as exc:
                    row.update(status='NOT_TESTABLE',reason=str(exc),Tobs=np.nan,null_median=np.nan,Delta=np.nan,null_reference_lower_025=np.nan,null_reference_upper_975=np.nan,raw_p=np.nan,p_lower=np.nan,p_upper=np.nan,null_draws=0,tail_probability_resolution=1/(B+1),minimum_two_sided_p=2/(B+1),null_reference_interval_type='not_computed')
                rows.append(row)
                pd.DataFrame(rows).to_csv(OUT/'qa/h2_results_progress.tsv',sep='\t',index=False)
                json_write(statusfile,{'status':'RUNNING','started_utc':started,'completed_registered_tests':len(rows),'last_test_id':test_id,'updated_utc':now()})
                print(f'{test_id} {name}: {row["status"]}; {len(rows)}/16 registered positions completed',flush=True)
        results=pd.DataFrame(rows)
        results['p_for_family_adjustment']=results.raw_p.fillna(1.)
        results['q_BH']=np.nan
        for family in ['H2_BASE','H2_LENGTH']:
            mask=results.family==family;pp=results.loc[mask,'raw_p'].to_numpy()
            qq=bh_fixed(pp,8);independent=multipletests(np.where(np.isfinite(pp),pp,1),method='fdr_bh')[1]
            np.testing.assert_allclose(qq,independent,rtol=0,atol=1e-15)
            results.loc[mask,'q_BH']=np.where(np.isfinite(pp),qq,np.nan)
        results.to_csv(OUT/'h2_results_long.tsv',sep='\t',index=False)
        np.savez_compressed(OUT/'null_distributions_16x10000.npz',test_ids=np.asarray(null_ids,dtype='U10'),values=matrix)
        pd.DataFrame(all_audits).to_csv(OUT/'random_draw_audit.tsv',sep='\t',index=False)
        # Reload every distribution and verify final table formula/shape alignment independently of serialization.
        saved=np.load(OUT/'null_distributions_16x10000.npz',allow_pickle=False)
        if saved['values'].shape!=(16,10000) or saved['test_ids'].tolist()!=null_ids:raise AssertionError('Null matrix output mismatch')
        for i,r in results.iterrows():
            if r.status=='EXECUTED':
                x=saved['values'][i];np.testing.assert_array_equal(x,np.load(OUT/'null_distributions'/f'{r.test_id}.npy',allow_pickle=False))
                again=empirical_summary(r.Tobs,x)
                for k in ['raw_p','null_median','Delta','null_reference_lower_025','null_reference_upper_975']:
                    if not math.isclose(float(r[k]),again[k],rel_tol=1e-14,abs_tol=1e-15):raise AssertionError('Reloaded formula failed')
        # Recheck protected inputs at end, including the freeze manifest itself.
        if sha(DESIGN/'freeze_manifest.json')!=FROZEN_MANIFEST_SHA256:raise AssertionError('Frozen manifest changed during execution')
        for path,digest in checked.items():
            if sha(Path(path))!=digest:raise AssertionError(f'Protected input changed during execution: {path}')
        qa={'status':'PASS' if (results.status=='EXECUTED').all() else 'COMPLETED_WITH_NOT_TESTABLE_POSITIONS','registered_tests':16,'executed_tests':int((results.status=='EXECUTED').sum()),'base_family_m':8,'length_family_m':8,'frozen_inputs_start_and_end_hash_pass':True,'sample_and_all_feature_library_totals_pass':True,'three_shares_sum_one_pass':True,'zero_denominator_source_genes_excluded_from_shares':int((~valid).sum()),'no_zero_denominator_in_fixed_background':True,'independent_donor_CPM_checks':direct,'null_no_duplicates_all_draws_checked':True,'candidate_union_absent_from_pools':True,'frozen_bin_capacities_and_membership_exact':True,'all_slot_weights_preserved':True,'random_full_gene_audit_draws_0based':list(AUDIT_DRAWS),'audited_draw_rows':len(all_audits),'independent_observed_and_audited_null_weighted_score_checks':True,'inclusive_tails_independent_ordered_search_check':True,'BH_independent_statsmodels_check':True,'saved_all_16_null_distributions_reloaded':True,'draws_per_test':B,'tail_grid':1/(B+1),'minimum_two_sided_empirical_p':2/(B+1),'no_bootstrap_effect_CI_age_cognition_H4_or_subclass_tests':True,'existing_age_and_combined_v1_results_not_modified':True,'H1_H3_known_before_H2_lock':True,'completed_utc':now()}
        json_write(OUT/'QA.json',qa)
        (OUT/'QA_report.md').write_text('# H2 execution QA\n\n'+ '\n'.join(f'- {k}: {json.dumps(v,ensure_ascii=False)}' for k,v in qa.items())+'\n\nThe null reference interval is not an observed-effect confidence interval. Full original and length families are reported regardless of outcome.\n')
        plot_results(results,OUT/'figures/h2_expression_reference_enrichment')
        write_report(results,qa,OUT/'concise_result.md')
        artifact_paths=[p for p in OUT.rglob('*') if p.is_file() and '__pycache__' not in str(p) and 'readiness' not in p.relative_to(OUT).parts and p.name not in ['execution.lock','output_manifest.json','run_status.json','execution.log']]
        json_write(OUT/'output_manifest.json',{'status':qa['status'],'created_utc':now(),'files':[{'path':str(p.relative_to(OUT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(artifact_paths)]})
        json_write(statusfile,{'status':'COMPLETED','qa_status':qa['status'],'started_utc':started,'completed_utc':now(),'completed_registered_tests':16,'no_more_null_draws_to_be_added':True})
        print('Completed. All registered results retained. See concise_result.md, h2_results_long.tsv and QA.json.',flush=True)
    except Exception as exc:
        json_write(statusfile,{'status':'NOT_EXECUTABLE_TECHNICAL_GATE' if isinstance(exc,NotTestable) else 'FAILED_PRESERVED_ARTIFACTS','started_utc':started,'failed_utc':now(),'error':str(exc),'traceback':traceback.format_exc(),'retry_policy':'Do not start another worker or overwrite results automatically; review failure and existing artifacts first.'})
        raise
    finally:
        fcntl.flock(lockfile.fileno(),fcntl.LOCK_UN);lockfile.close()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--validate-lock-only',action='store_true',help='Hashes/schema/capacity only, no real expression or null computation')
    modes.add_argument('--execute-locked-h2',action='store_true',help='Only the parent launches this actual 16-test analysis once')
    args=parser.parse_args()
    if args.validate_lock_only:
        m,checked,bg,bl,member,ml,plans=verify_locked_contract()
        report={'status':'READINESS_CONTRACT_PASS','real_expression_matrix_loaded':False,'actual_H2_score_or_null_or_P_computed':False,'verified_file_count':len(checked),'manifest_sha256':FROZEN_MANIFEST_SHA256,'background_n':len(bg),'control_pool_n':int(bg.noncandidate_pool.sum()),'registered_sampling_plans':len(plans),'fixed_families':m['registered_families'],'checked_at_utc':now()}
        (OUT/'readiness').mkdir(exist_ok=True);json_write(OUT/'readiness/locked_contract_check.json',report)
        print(json.dumps(report,ensure_ascii=False,indent=2))
    else:execute()

if __name__=='__main__':main()
