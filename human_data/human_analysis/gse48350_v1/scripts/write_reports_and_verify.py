"""Read completed registered outputs; no model refits or parameter changes."""

# Release-only safety/path adapter; see provenance/adaptation_diffs.
from pathlib import Path as _ReleasePath
import sys as _release_sys
_release_home = next(p for p in _ReleasePath(__file__).resolve().parents if (p / 'release_support.py').is_file())
_release_sys.path.insert(0, str(_release_home))
from release_support import gate, release_root, resolve_source
gate(__file__, __name__)
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib
import numpy as np,pandas as pd
from PIL import Image
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'results';L=ROOT/'logs'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def md(d):
 def text(x):
  if pd.isna(x):return 'NA'
  if isinstance(x,float):return f'{x:.6g}'
  return str(x).replace('|',' / ')
 return '| '+' | '.join(map(str,d.columns))+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(text(x) for x in row)+' |' for row in d.itertuples(index=False,name=None))
run=json.loads((L/'run_complete.json').read_text());assert run['status']=='PASS';res=pd.read_csv(R/'registered_results_long.tsv',sep='\t');scores=pd.read_csv(R/'donor_scores.tsv',sep='\t');coverage=pd.read_csv(R/'set_coverage_and_gates.tsv',sep='\t');genes=pd.read_csv(R/'gene_selection.tsv',sep='\t');probes=pd.read_csv(R/'qualified_probe_lists.tsv',sep='\t');sample=pd.read_csv(R/'sample_manifest.tsv',sep='\t');flow=pd.read_csv(R/'model_sample_manifest.tsv',sep='\t');loo=pd.read_csv(R/'leave_one_donor_out_beta_only.tsv',sep='\t');diag=pd.read_csv(R/'model_diagnostics.tsv',sep='\t');params=json.loads((ROOT/'run_parameters.json').read_text());versions=json.loads((ROOT/'software_versions.json').read_text());source=json.loads((ROOT/'input_manifest.json').read_text())
assert sha(R/'registered_results_long.tsv')==run['result_SHA256']
for path,hsh in source['files'].items():assert sha(Path(path))==hsh,path
for n,hsh in versions['scripts_sha256'].items():assert sha(ROOT/'scripts'/n)==hsh,n
assert set(res.test_id)=={'P01','E01','D01','S07','S08','S16'};assert res.test_id.is_unique
assert res.set_index('test_id').loc['S16','status']=='pending_SEA_actual_A_intersection';assert res.set_index('test_id').loc['S16',['raw_P','adjusted_P']].isna().all()
assert res[res.family.ne('E_AD')].adjusted_P.isna().all();assert res[res.family.ne('E_AD')].adjustment_status.eq('pending-family').all()
assert np.isclose(res.set_index('test_id').loc['D01','raw_P'],res.set_index('test_id').loc['D01','adjusted_P'])
for (name,scale),g in scores.groupby(['gene_set','scale']):
 ref=g[g.in_normal43_reference];assert len(ref)==43 and ref['sample'].is_unique;assert np.isclose(ref.Y_score_reference_SD.mean(),0,atol=1e-12);assert np.isclose(ref.Y_score_reference_SD.std(ddof=1),1,atol=1e-12);assert g.A_n.nunique()==1 and g.fixed_mean_denominator.nunique()==1
 assert g.applied_log2_count.eq(1 if scale=='log2_ratio' else 0).all();assert g.source_log2_count.eq(0).all()
for name,g in genes.groupby('gene_set'):
 log=set(g[(g.scale=='log2_ratio')&g.included_in_A].HGNC_ID);raw=set(g[(g.scale=='raw_ratio')&g.included_in_A].HGNC_ID);assert log==raw
for tid,g in flow.groupby('test_id'):
 rr=res.set_index('test_id').loc[tid];assert len(g)==62 and g['sample'].is_unique;assert g.complete_case.sum()==rr['n'];assert g.population_eligible.sum()==rr.n_population
for tid,g in loo.groupby('test_id'):
 rr=res.set_index('test_id').loc[tid];assert len(g)==rr['n'];assert np.isclose(g.beta.min(),rr.LOO_beta_min);assert np.isclose(g.beta.max(),rr.LOO_beta_max)
assert not any('P' in c or 'p_value' in c for c in loo.columns) # no LOO hypothesis testing
assert len(scores)==210 and len(flow)==310 and len(genes)==308
for p in (ROOT/'figures').glob('*.png'):
 with Image.open(p) as im:im.verify()
for p in (ROOT/'figures').glob('*.pdf'):assert p.read_bytes().startswith(b'%PDF-')
# Report exactly the fixed target coefficients, not covariate tests.
r=res[res.status.eq('tested')].copy();labels={'P01':'core年龄，主log2尺度','E01':'METTL3程序年龄，辅助log2尺度','D01':'core AD扩展，log2尺度','S07':'core年龄，raw-ratio敏感性','S08':'METTL3程序年龄，raw-ratio敏感性'};r['analysis']=r.test_id.map(labels)
show=r[['test_id','analysis','n','A_n','beta','SE_HC3','df','CI95_low','CI95_high','raw_P','adjusted_P']].rename(columns={'analysis':'分析','n':'n','A_n':'gene A','beta':'β','SE_HC3':'SE HC3','CI95_low':'95%CI下限','CI95_high':'95%CI上限','raw_P':'原始P','adjusted_P':'家族调整P'})
desc=r[['test_id','age_VIF','max_hat','n_hat_flag','max_Cook','n_Cook_flag','LOO_beta_min','LOO_beta_max','LOO_direction_flip_n']]
flagged=diag[diag.hat_flag|diag.Cook_flag][['test_id','sample','donor_id','hat','Cook_distance','hat_flag','Cook_flag']]
refs=pd.read_csv(R/'score_reference_parameters.tsv',sep='\t').drop(columns='reference_samples')
concise=f'''# GSE48350冻结v1：有界执行结果

内部前瞻锁定：2026-10-12 21:14:35（UTC+8）；本次只执行P01、E01、D01、S07、S08。计划/登记表SHA256、候选五文件MD5与执行后输入SHA256全部通过。未重下载、重映射、重跑小鼠；未做单基因扫描或SEA/认知/H2/遗传学分析。

## 固定方法与覆盖

主尺度严格一次`log2(positive ratio)`，不加1；每基因取全部合格单身份探针的变换后中位数。以正常海马43人为固定参照：逐gene样本SD标准化→集合等权mean-z→按参照score样本SD再缩放。core为21/22（95.45%，门槛18），METTL3cKO_HIP_down为113/120（94.17%，门槛96）；所有技术可用基因均通过参照有限值/方差门槛。Lynx1多同源与FBXO41歧义探针未纳，未依据结果救回或换候选。缺乏独立微阵列实际检出门控，注释覆盖不等于表达检出。

年龄模型固定`Y ~ 1 + (age-80)/10 + sex`（女0/男1），正常43、女20/男23、20–99岁。D01固定共同年龄60–94岁，实际**40人=正常21+AD19，女20/男20**，拟合`Y ~ 1 + age10 + sex + AD`并仅检验AD系数；score仍用全正常43参照，AD不进入参照。OLS HC3，双侧t，df=n−rank(X)，逐项95%CI。

## 登记结果

{md(show)}

β单位：年龄项为**正常43参照score SD/10岁**；D01为相同score SD的**AD相对正常差值**，不是年龄效应。CI为逐项区间，不是多重同时区间。

- **P01**：core年龄点估计为负但CI跨0；不能解释为已证明无效，亦不声称验证成功。
- **E01**：主METTL3表达程序呈负向年龄估计（−0.1504 SD/10岁；逐项CI −0.2746至−0.02625；raw P=0.01884）。这是E_EXPR中的辅助结果；**家族未齐，未宣布FDR通过或确认性复制**。
- **D01**：AD扩展β=−0.7779，CI跨0；唯一E_AD族m=1，因此本族调整P等于raw P=0.06151。它不是年龄验证，不作因果/疾病排除解释。
- **S07/S08**：两项已登记raw-ratio敏感性独立列出，基因/探针与normal43参照群不变，在raw尺度重算均值/SD。**没有因S08 raw P更小而把raw-ratio替换主log2尺度**。
- **S16**：`pending_SEA_actual_A_intersection`；未计算交集得分/关联，不能以SEA var存在替代实际表达门控后的A。

## 多重校正与稳定性

PRIMARY固定Holm m=3、E_EXPR固定BH m=11、E_SENS固定BH m=23均尚有未运行项目，所有相应`adjusted_P=NA`、`adjustment_status=pending-family`。未知P没有代入1来假装家族完成，也没有被当阴性结论。完整family汇总留给负责人后续授权。

5个模型均过覆盖/设计门槛；无缺失导致的模型样本流失。所有模型均触发至少一个Cook影响诊断标志，D01另有1人hat超阈值，不能省略该不稳定提示；未删除任何供体，未改变模型。固定score的LOO效应方向翻转均0，但这不消除小样本影响性、混杂或宽CI。

{md(desc)}

## 文件

- `results/registered_results_long.tsv`：5已检验+1pending登记位置，完整效应/CI/raw P/family/status/覆盖/诊断。
- `results/gene_selection.tsv`、`qualified_probe_lists.tsv`、`set_coverage_and_gates.tsv`。
- `results/gene_reference_parameters.tsv`、`score_reference_parameters.tsv`、`donor_scores.tsv`、`sample_manifest.tsv`、`model_sample_manifest.tsv`及每项model design。
- `results/model_diagnostics.tsv`、`leave_one_donor_out_beta_only.tsv`、`manual_HC3_spotchecks.tsv`。
- `figures/registered_effects.pdf/png`及5项diagnostics PDF/PNG（English labels）。
- `QA.md`、`run_parameters.json`、`input_manifest.json`、`software_versions.json`、`logs/run_complete.json`、`logs/final_QA.json`。

解释限于成人正常bulk海马横断面表达集合关联；不是单细胞m6A、METTL3酶活、衰老因果或认知改善。H2尚未执行，不能宣称人类细胞富集环节完成。
'''
(ROOT/'concise_result.md').write_text(concise)
qa=f'''# GSE48350 v1 — execution QA

## 1. Scope, frozen authorization, provenance

- Registered executed IDs: P01/E01/D01/S07/S08; S16 remains pending. No additional models, univariate gene scans, SEA/mouse/cognition/H2/MR analyses, or result-dependent exclusions were run.
- Frozen plan SHA256: `c8ecd963efff5ea4b078c223fddfc2a3ccbd9e2a6b1ff987aa1f74e17d5bc31a`.
- Frozen registry SHA256: `d20d4220ec0b5c1cffa563c02b9b510ab33f95657fdfb132302be2f1ac7a43a1`.
- Candidate XLSX MD5: `7fd7ac5097f9dad9977eea804edeb36e`; all five frozen file MD5s checked before execution.
- Internal prospective lock recorded `{params['prospective_lock']}`. Execution started `{params['started_at_utc']}`. Synthetic checks and fixed scoring manifests were completed before human model fitting. This is an internal lock, not an external registry claim.
- `{run['source_file_n']}` source files have SHA256 snapshots in `input_manifest.json`; every hash was rechecked after execution and again during report generation. Preflight/raw/frozen/plan files unchanged. Crosswalk reused, no remapping/download.

## 2. Implementation and software

Python `{versions['python'].split()[0]}`; NumPy `{versions['numpy']}`; pandas `{versions['pandas']}`; SciPy `{versions['scipy']}`; statsmodels `{versions['statsmodels']}`; matplotlib `{versions['matplotlib']}`. Isolated `.venv`; software versions and execution-script SHA256 in `software_versions.json`, exact distributions in `requirements.lock.txt`.

```bash
cd "{ROOT.parents[1]}"
# Read-only inspection; preserve completed models.
python run.py gse --dry-run
python run.py gse-report --dry-run
```

`run_gse_v1.py` refuses to overwrite a completed run (`logs/run_complete.json`). For authorized reproduction use an isolated output directory rather than deleting the successful outputs. Reports re-read completed tables only and do not refit models.

## 3. Scale, aggregation, missingness and reference checks

- Input source is the checksum-verified positive GeneSpring ratio matrix. Original HIP data 54,675×62, strictly finite positive, min0.01/max1438.7626 matching preflight. Source matrix and stored source-ratio data never overwritten with logged values.
- Main `applied_log2_count=1`, raw sensitivity=0; all `source_log2_count=0`. Transform rejects prelogged flag, nonpositive or nonfinite input. No +1, double log, CPM, quantile normalization or per-chip normalization.
- Median is calculated **after** probe transformation, including even probe counts. Gene→probe lists are fixed single-identity mappings from preflight; no variance/intensity/phenotype-based probe selection.
- Normal43 fixed reference, ddof=1 gene and score SD, SD>1e-8 gate. A is established on main log2 scale, then held identical for raw sensitivities; no AD-specific or model-specific gene screening.
- Missing/nonfinite reference gene would be wholly excluded with coverage recalculated; non-reference sample with any A-gene missing receives score NA, never a reduced-denominator mean. These behaviors have synthetic tests. Actual data had no such missingness/zero variance.
- Every set/scale has normal-reference Y mean≈0 and sample SD≈1 after saved-table roundtrip; fixed denominator and logflags checked on every score row.

{md(coverage[['gene_set','scale','source_member_n','G0_one_to_one_n','technical_available_n','reference_variance_pass_n','A_n','number_coverage','minimum_genes','coverage_gate_pass']])}

Independent actual-expression counts are deliberately NA (`not_assessable_from_ratio_matrix`), not filled with annotation availability. Both sets exceed the prelocked 80% and minimum10 rule. Equal-weight coverage equals member-count coverage. No modules or kME variants were run here.

Saved S reference constants:

{md(refs)}

## 4. Sample identity and registered population

- All62 HIP source sample IDs are unique. Donor IDs are parsed using the GEO title namespaces `normal_indiv_*` versus `AD_*`; normal title IDs match the source `individual` metadata field. No repeated HIP donor record within either source namespace. Clinical names/identities were not invented; uniqueness is based on available GEO donor labeling, not an independent genotyping proof.
- Normal43: 20F/23M;31 distinct ages; range20–99. No missing age/sex/score, all43 enter each age model. No age cutpoint search.
- Source AD19 range60–94. Prespecified common support [60,94] contains21 normal+19 AD=40,20F/20M;18 distinct ages. Normal reference still includes all43, not only common-support21. AD never affects reference mu/sd or gene A. No AD-only age model/interaction was run.
- `sample_manifest.tsv`62 rows; `model_sample_manifest.tsv`310 rows, one for each HIP sample×five registered attempted models with explicit population/complete-case status; `donor_scores.tsv`210 rows (core62×2 scales, program normal43×2 scales).

## 5. OLS, inference, gates and independent tests

For each model the design gate checks n≥max(20,3p,p+15), full rank, ≥5 ages spanning≥10y, sex each≥3. Actual models all pass: age models n43,p3,df40; AD n40,p4,df36. No gate was relaxed or covariate dropped. HC3 uses residual/(1-hat) squared; target t tests and 95%CI use df=n−rank(X), not normal z inference. CI is pointwise, not multiplicity-simultaneous.

- **Known-synthetic pipeline test before human fitting:** deterministic full-rank40-person design with true age β−0.5 and constructed orthogonal heteroscedastic residual. OLS recovered −0.4999999999999998; manual HC3 SE0.05008733647, df37, raw P4.811562967e−12; statsmodels and independent matrix-algebra HC3/t/CI agree. See `logs/synthetic_validation.json`.
- Unit tests also cover forbidden double log, invalid positive-ratio input, median-of-logs order, fixed reference scoring, constant denominator with non-reference NA, coverage boundary18/22 and96/120, minimum n/rank/age/sex failures, and LOO beta-only output.
- **Manual human score spotcheck:** GSM300168/CACNA1E,6 frozen probes. Explicit scalar log2 values yield probe median−0.007728276184377075; explicit mean of21 gene z values gives S0.30894362453332785 and Y0.44912156211342147, matching the vector pipeline. Full source numbers in `logs/manual_score_spotcheck.json`.
- **Independent manual HC3 spotchecks for all five registered target coefficients:** invert X'X, calculate OLS beta/residuals/hats, HC3 sandwich, Student t P and t critical CI, compared to statsmodels at rtol1e−10/atol1e−12. `results/manual_HC3_spotchecks.tsv` all match.
- Neither unit tests nor QA generate extra human hypotheses; non-target covariate coefficients are stored descriptively without P-value scans.

## 6. Diagnostics, LOO, and no deletion

Age VIF, leverage and ordinary OLS Cook distances are descriptive design/influence diagnostics; Cook is not represented as a new hypothesis test or a heteroscedasticity-corrected test. Fixed thresholds: VIF≥10,hat>2p/n,Cook>4/n. LOO refits only beta with fixed original scores/reference/A; no LOO P values or best-case result selection.

{md(desc)}

All five models have Cook flags, so `unstable_flag=True`; D01 also has one high-leverage observation. No age-VIF≥10. No donor removed. All LOO target coefficients remain negative (zero sign flips), but this is only a sensitivity description, not proof of robustness or causality. Flagged observations, retained in their registered models:

{md(flagged)}

Residual/fitted, leverage and LOO plots exist for each test; summary forest plot separates main log2 scale, disease extension and raw sensitivities. Six PDF/PNG pairs passed format integrity; summary and AD diagnostic PNGs visually inspected for clipping/label correctness. English labels only.

## 7. Family handling and unresolved tasks

- PRIMARY(P01–P03,m3),E_EXPR(m11),E_SENS(m23): raw P saved for executed positions; adjusted P stays NA and status `pending-family`. Unrun unknown P values were not assigned1 or interpreted as negative. Parent retains responsibility for complete frozen-family assembly.
- E_AD is a complete single-position family: D01 adjusted P=raw P0.0615089438 (full precision in TSV).
- S16 has public raw/adjusted P NA, status pending_SEA_actual_A_intersection. No intersection inferred from platform availability; must await SEA main actual-expression+variance A for CA3 and SUB.
- No H2 completion/cell enrichment claim, no m6A or METTL3 enzyme-activity claim, no individual longitudinal aging or cognition-improvement claim. Cross-sectional bulk effects only. Higher raw-ratio sensitivity significance never changes the prelocked log2 main analysis.

## 8. Final roundtrip and integrity

Results, scoring/selection/ref/model manifests and LOO ranges reread; counts, unique IDs, all normal-score means/SDs, A identity across scales, source/transform flags, model population and complete-case counts, family pending fields, source hashes and original execution-script hashes all asserted. PDF magic bytes and PNG integrity verified. `logs/final_QA.json` records final PASS and output SHA256.
'''
# Avoid manually rounded P in a machine-generated QA description.
qa=qa.replace('0.0615089438',f"{res.set_index('test_id').loc['D01','raw_P']:.12g}")
(ROOT/'QA.md').write_text(qa)
outputs=[ROOT/'QA.md',ROOT/'concise_result.md',ROOT/'run_parameters.json',ROOT/'software_versions.json',ROOT/'input_manifest.json',*sorted(R.glob('*')),*sorted((ROOT/'figures').glob('*'))]
final={'status':'PASS','completed_at_utc':datetime.now(timezone.utc).isoformat(),'source_files_unchanged':True,'execution_script_hashes_match':True,'results_reread':True,'score_mean0_SD1_roundtrip':True,'same_A_raw_vs_log':True,'model_populations_verified':True,'family_pending_NA_verified':True,'LOO_beta_only_verified':True,'PNG_PDF_integrity':'PASS','PNG_visual_QA':'registered_effects and D01_diagnostics inspected','executed_test_ids':['P01','E01','D01','S07','S08'],'pending_test_ids':['S16'],'scores_rows':len(scores),'model_manifest_rows':len(flow),'gene_selection_rows':len(genes),'qualified_probe_set_rows':len(probes),'qualified_probe_unique_n':probes.probe_id.nunique(),'report_script_sha256':sha(Path(__file__)),'outputs':{str(p.relative_to(ROOT)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in outputs}}
(L/'final_QA.json').write_text(json.dumps(final,ensure_ascii=False,indent=2));print(json.dumps({k:v for k,v in final.items() if k!='outputs'},indent=2));print(md(show))
