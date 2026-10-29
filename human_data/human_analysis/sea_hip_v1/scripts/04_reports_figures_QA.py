"""Final output verification, plots and reporting; no model refits."""

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
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image
ROOT=Path(__file__).resolve().parents[1];AN=ROOT.parent;G=AN/'gse48350_v1';C=AN/'combined_v1';R=ROOT/'results';L=ROOT/'logs';D=ROOT/'derived'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,p):p.write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))
def md(d):
 def cell(v):
  if pd.isna(v):return 'NA'
  if isinstance(v,float):return f'{v:.6g}'
  return str(v).replace('|',' / ').replace('\n',' ')
 return '| '+' | '.join(map(str,d.columns))+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(cell(x) for x in row)+' |' for row in d.itertuples(index=False,name=None))
execution=json.loads((L/'execution_PASS.json').read_text());assembly=json.loads((C/'logs/assembly_PASS.json').read_text());assert execution['status']==assembly['status']=='PASS'
for name,hsh in execution['all_outputs'].items():assert sha(ROOT/name)==hsh
for name,hsh in assembly['outputs'].items():assert sha(C/'results'/name)==hsh
for path,hsh in json.loads((ROOT/'input_manifest.json').read_text())['files'].items():assert sha(Path(path))==hsh
versions=json.loads((ROOT/'software_versions.json').read_text())
for n,hsh in versions['scripts_SHA256'].items():assert sha(ROOT/'scripts'/n)==hsh
assert sha(G/'scripts/fixed_pipeline.py')==versions['reused_OLS_helper_SHA256']
res=pd.read_csv(R/'registered_results_long.tsv',sep='\t');combined=pd.read_csv(C/'results/combined_results_long.tsv',sep='\t');fam=pd.read_csv(C/'results/family_summary.tsv',sep='\t');cov=pd.read_csv(R/'score_coverage.tsv',sep='\t');sel=pd.read_csv(R/'gene_selection.tsv',sep='\t');gre=pd.read_csv(R/'gene_reference_parameters.tsv',sep='\t');sre=pd.read_csv(R/'score_reference_parameters.tsv',sep='\t');scores=pd.read_csv(R/'donor_scores.tsv',sep='\t');flow=pd.read_csv(R/'model_sample_manifest.tsv',sep='\t');diag=pd.read_csv(R/'model_diagnostics.tsv',sep='\t');loo=pd.read_csv(R/'leave_one_donor_out_beta_only.tsv',sep='\t');tm=pd.read_csv(D/'target_donor_manifest.tsv',sep='\t');gates=pd.read_csv(R/'actual_expression_and_variance_gates.tsv',sep='\t');refcov=pd.read_csv(R/'covariate_reference_parameters.tsv',sep='\t')
# After serialization: all reference scores 0 mean and SD1, weights/genes fixed.
for sid,g in scores.groupby('score_id'):
 rr=g[g.in_primary20_reference];assert np.isclose(rr.Y.mean(),0,atol=1e-12) and np.isclose(rr.Y.std(ddof=1),1,atol=1e-12);assert g.A_n.nunique()==1 and g.fixed_weight_denominator.nunique()==1
 refs=gre[gre.score_id.eq(sid)];assert refs.HGNC_ID.is_unique;assert refs.shape[0]==int(g.A_n.iloc[0]);assert np.isclose(refs.weight.sum(),g.fixed_weight_denominator.iloc[0]);assert (refs.sample_SD>1e-8).all()
for _,r in cov.iterrows():assert r.coverage_gate_pass==(r.A_n>=max(10,int(np.ceil(.8*r.G0_n))) and r.weight_coverage>=.8)
for tid,g in flow.groupby('test_id'):
 r=res.set_index('test_id').loc[tid];assert g.donor.is_unique;assert g.registered_population_eligible.sum()==r.n_population;assert g.complete_case.sum()==r['n']
for tid,g in loo.groupby('test_id'):
 r=res.set_index('test_id').loc[tid];assert len(g)==r['n'];assert np.isclose(g.beta.min(),r.LOO_beta_min);assert np.isclose(g.beta.max(),r.LOO_beta_max)
assert not any('P' in c or 'p_value' in c for c in loo.columns)
assert res.set_index('test_id').loc[['S17','S18'],'status'].eq('not_testable_coverage').all();assert res.set_index('test_id').loc[['S17','S18'],'raw_P'].isna().all();assert res.set_index('test_id').loc[['S17','S18'],'A_n'].eq(92).all();assert res.set_index('test_id').loc[['S17','S18'],'minimum_genes'].eq(96).all()
assert res.set_index('test_id').loc[['S10','S11','S12'],'status'].eq('identical_prespecified_subset').all();assert res.set_index('test_id').loc[['S10','S11','S12'],'raw_P'].isna().all();assert res.set_index('test_id').loc['S16','status']=='tested_reused_identical_GSE_P01'
assert combined[combined.family.eq('E_SENS')].adjusted_P.isna().all();assert combined[combined.position_resolution.eq('unknown_pending')].p_for_adjustment_internal.isna().all();assert len(combined)==62
# Explicit S16 source-row view, referencing already frozen/logged GSE mapping.
gsel=pd.read_csv(G/'results/gene_selection.tsv',sep='\t');s16sel=gsel[(gsel.gene_set=='core23_1to1')&(gsel.scale=='log2_ratio')].copy();common=set(json.loads((D/'core23_1to1_crossplatform_intersection.json').read_text())['intersection']);s16sel['score_id']='S16';s16sel['intersection_included_A']=s16sel.HGNC_ID.isin(common);s16sel['subset_reason']=np.where(s16sel.intersection_included_A,'included_same_GSE_A','original_GSE_exclusion_or_crossplatform_intersection');s16sel.to_csv(R/'S16_gene_selection_source_rows.tsv',sep='\t',index=False)
flow['model_status']=flow.test_id.map(res.set_index('test_id').status);flow['fitted_this_run']=flow.complete_case & flow.test_id.map(res.set_index('test_id').model_execution).eq('fitted_registered_position_once');flow.to_csv(R/'model_flow_with_gate_status.tsv',sep='\t',index=False)
# Synthetic/manual checks were completed before any current-model fit.
assert json.loads((L/'synthetic_PASS.json').read_text())['status']=='PASS';assert pd.read_csv(R/'manual_HC3_spotchecks.tsv',sep='\t').match.all();assert json.loads((L/'scalar_primary_score_spotcheck.json').read_text())['PASS']
# Independent weighted module score reconstruction from SAVED scalar parameters.
sid='CA3_1::CA3-M2';rr=sre[sre.score_id.eq(sid)].iloc[0];donor=rr.reference_donors.split(';')[0];g=gre[gre.score_id.eq(sid)];xx=pd.read_csv(D/'CA3_1_candidate_log2CPM1.tsv.gz',sep='\t',index_col=0);manual_S=sum(float(row.weight)*(float(xx.loc[row.HGNC_ID,donor])-float(row.mu))/float(row.sample_SD) for row in g.itertuples())/float(g.weight.sum());manual_Y=(manual_S-rr.S_mean)/rr.S_sample_SD;actual=scores[(scores.score_id==sid)&(scores.donor==donor)].iloc[0];assert np.isclose(manual_S,actual.S,atol=1e-12) and np.isclose(manual_Y,actual.Y,atol=1e-12);dump({'PASS':True,'score_id':sid,'donor':donor,'gene_A_n':len(g),'fixed_kME_denominator':g.weight.sum(),'manual_S':manual_S,'saved_S':actual.S,'manual_Y':manual_Y,'saved_Y':actual.Y},L/'manual_weighted_module_spotcheck.json')
# Simple diagnostic pages for every fitted model; no tests calculated here.
with PdfPages(ROOT/'figures/all_fitted_model_diagnostics.pdf') as pdf:
 for tid,g in diag.groupby('test_id',sort=False):
  ll=loo[loo.test_id.eq(tid)];r=res.set_index('test_id').loc[tid];fig,ax=plt.subplots(1,4,figsize=(13,3.3));ax[0].scatter(g.fitted,g.residual,s=15);ax[0].axhline(0,color='grey',lw=.7);ax[0].set(xlabel='Fitted score',ylabel='OLS residual',title=tid+' residuals')
  ax[1].scatter(range(len(g)),g.hat,s=13);ax[1].axhline(r.hat_threshold,color='red',ls='--',lw=.7);ax[1].set(xlabel='Donor order',ylabel='Hat',title='Leverage')
  ax[2].scatter(range(len(g)),g.Cook_distance,s=13);ax[2].axhline(r.Cook_threshold,color='red',ls='--',lw=.7);ax[2].set(xlabel='Donor order',ylabel='Cook distance',title='Influence (no deletion)')
  ax[3].scatter(range(len(ll)),ll.beta,s=13);ax[3].axhline(r.beta,color='black',ls='--',lw=.7);ax[3].axhline(0,color='grey',lw=.7);ax[3].set(xlabel='Omitted donor',ylabel='Age beta',title='LOO beta only (no P)');fig.tight_layout();pdf.savefig(fig)
  if tid in ['P02','P03','M05']:fig.savefig(ROOT/'figures'/f'{tid}_diagnostics.png',dpi=160)
  plt.close(fig)
def forest(frame,title,name,base=C):
 z=frame[frame.status.str.startswith('tested')].copy();fig,ax=plt.subplots(figsize=(10,max(3.2,.38*len(z)+1.6)))
 for i,r in enumerate(z.itertuples()):ax.errorbar(r.beta,i,xerr=[[r.beta-r.CI95_low],[r.CI95_high-r.beta]],fmt='o',capsize=3,color='#24476b')
 labels=[f'{r.test_id}  {r.target}  {r.gene_set}' for r in z.itertuples()];ax.set_yticks(range(len(z)),labels,fontsize=8);ax.invert_yaxis();ax.axvline(0,color='grey',lw=.8);ax.set(xlabel='Fixed-reference score SD per 10 years',title=title+'\nPointwise 95% HC3 t intervals');fig.tight_layout();fig.savefig(base/'figures'/f'{name}.png',dpi=170);fig.savefig(base/'figures'/f'{name}.pdf');plt.close(fig)
forest(combined[combined.family=='PRIMARY'],'Primary family: fixed Holm m=3; no item passes 0.05','primary_effects')
forest(combined[combined.family=='E_EXPR'],'Exploratory expression family: BH m=11; no q<0.05','expression_effects')
forest(combined[combined.family=='E_MODULE'],'Exploratory modules: BH m=7; M05 q<0.05 (not a primary test)','module_effects')
forest(combined[combined.family=='E_SENS'],'Executed sensitivity positions: family adjustment PENDING','sensitivity_effects')
# SEA-only primary effect figure also lives with its own run.
forest(res[res.test_id.isin(['P02','P03'])],'SEA primary effects; consult combined Holm family','SEA_primary_effects',ROOT)
primary=combined[combined.family=='PRIMARY'][['test_id','cohort','target','gene_set','n','A_n','beta','SE_HC3','df','CI95_low','CI95_high','raw_P','adjusted_P']]
expl=combined[combined.family=='E_EXPR'][['test_id','target','gene_set','n','A_n','beta','CI95_low','CI95_high','raw_P','adjusted_P']];mods=combined[combined.family=='E_MODULE'][['test_id','target','gene_set','n','A_n','beta','CI95_low','CI95_high','raw_P','adjusted_P']];sens=combined[combined.family=='E_SENS'][['test_id','target','gene_set','status','n','A_n','beta','CI95_low','CI95_high','raw_P']]
maincoverage=cov[cov.variant.eq('main')][['target','gene_set','reference_n','expression_min_donors','technical_n','actual_expression_pass_n','A_n','G0_n','number_coverage','weight_coverage']]
pdiag=res[res.test_id.isin(['P02','P03','S05','S06','S23'])][['test_id','n','p','df','age_min','age_max','female_n','male_n','age_VIF','max_hat','max_Cook','n_hat_flag','n_Cook_flag','LOO_beta_min','LOO_beta_max','LOO_direction_flip_n']]
merges=pd.read_csv(D/'candidate_feature_identity_audit.tsv',sep='\t');merges=merges[merges.n_source_features_for_HGNC.gt(1)][['HGNC_IDs','feature_symbol','ensembl_id','current_human_symbol','HGNC_previous_symbols']]
combined_report=f'''# 人类表达v1组合结果（固定家族已组装）

2026-10-12。依据21:14:35（UTC+8）内部前瞻锁定计划，GSE已完成五项原结果原样复用；本次新增SEA登记表达位置及S16。**无单gene扫描、无小鼠重跑、无认知/H2/H4分析、无结果后阈值/协变量调整。**

## 1. 唯一主要家族：未达到预设Holm支持门槛

{md(primary)}

β为各自固定参照score SD/10岁；HC3、双侧t，CI是逐项95%区间，不是多重同时区间。**P01–P03均未通过固定Holm m=3的0.05门槛。** CA3与SUB的主程序点估计均为负；SUB P03逐项CI不含0、raw P≈0.03032，但Holm P≈0.09097，不能表述为确认性验证成功，更不能写“CA3与SUB均验证”。CA3 CI较宽；未过门槛不等于证明没有效应。GSE正常成人年龄跨度与SEA高龄病理谱估计对象不同，不合并效应、不作meta分析。

## 2. 辅助表达家族 E_EXPR：固定BH m=11，无q<0.05

{md(expl)}

GSE主程序E01的raw P≈0.01884对应q≈0.10359；CA1主程序E06 raw P≈0.006816对应q≈0.07498。它们不作为主要结果未达门槛后的“救援”或把CA1替换预设CA3/SUB。

## 3. 模块家族 E_MODULE：仅M05达到探索性BH门槛

{md(mods)}

**CA1-M9/M05为探索性负向年龄关联：β≈−0.5784，CI −0.9099至−0.2469，BH q≈0.008277。** 仍是独立登记的探索模块结果，不升级为主H3、不证明CA1比CA3/SUB更显著或更特异，也不声称复制了小鼠方向（冻结CA1-M9小鼠beta为正且原P不显著）。不同target各自标准化，未检验target×age。

## 4. SEA实际可评分基因与样本

原层UMIs准备QA为PASS；counts.npz为3402×36601 CSR int64；全部6,196,354,669 UMI守恒。只组合预设target，不使用X。全36,601 feature的UMI和作为CPM分母，`log2(CPM+1)`一次。主参考供体n：CA3=29、SUB=31、CA1=38、DG=39；actual-expression门槛依次6、7、8、8供体达到CPM≥1，再通过参照有限值与样本SD>1e−8。

{md(maincoverage)}

主METTL3程序技术覆盖120不等于actual A：CA3实际100/120（83.33%）、SUB106/120（88.33%），均通过≥96的固定门槛。所有模块按原正kME加权，数量和kME权重覆盖分别检查，未重估人类eigengene。core四target均22/22。每个主A/基因mu与sd/score缩放及PMI-RIN中心值均固定于本target≥20核参照；未按年龄、诊断或模型完整病例重新筛基因。

## 5. 有限敏感性与pending位置

{md(sens)}

- S17/S18的GSE主程序A∩SEA CA3实际A∩SUB实际A只有**92/120=76.67%**，低于预设96/120。标记`not_testable_coverage`，**没有拟合关联、未降低门槛**。
- S16 core交集为21，恰与原GSE主A相同；用已保存log2 gene aggregates构建一次交集score并确认与原score一致后，**直接引用P01效应而未重拟合任何GSE基础模型**。登记位置保留且注明结果完全相同，不视作新独立证据。
- S10/S11/S12高kME子集与相应主A和权重完全相同；引用父模型效应，公开P=NA，未额外拟合；这些位置只保留固定家族槽位。
- S03（SUB≥10核）与P03实际供体完全相同，结果相同，不能当独立重复。其余10/30阈值与+Braak、SUB2组成都保留主A、主参考、协变量编码；Braak无缺失，因此无需另作完整病例对照拟合。SUB组成项改变的是预设敏感性模型，不替代主模型。
- **E_SENS固定23位置：14项有P、5项terminal NA、4项认知敏感性仍未知。整个族调整P仍NA/pending-family，未用未知P=1冒充完成。** S19–S22未运行。
- E_COG8、H2未来8、H4未来1均pending/registered-not-executable，未读取认知-表达结果。E_AD唯一项沿原GSE D01：调整P=raw P≈0.06151，不是年龄验证。

## 6. 诊断与解释边界

{md(pdiag)}

P02/P03均有影响点/杠杆标志，供体全部保留；LOO年龄beta始终负，但不据此忽略宽CI或多重校正。全模型诊断和beta-only LOO均已保存，无LOO P值、无择优删除。SEA年龄使用来源提供数值；没有推断超过来源精度的真实年龄，不把高龄病理谱当健康lifespan aging。PMI小时、RIN及Braak罗马级别已与原donor表精确核对，Braak只作线性0–6；没有教育/认知/APOE追加协变量。

**结论限于表达集合关联。RNA并非m6A或METTL3酶活实测；本轮H2未完成，不能宣称人类细胞富集链路完成；不声称认知改善、因果或纵向个体下降。**

## 7. 交付索引

- `results/combined_results_long.tsv`、`registry_tracking_62.tsv`：全部62个登记位置，未知与不可检验明确区分。
- `results/family_summary.tsv`及各族`*_results.tsv`；PRIMARY/E_EXPR/E_MODULE/E_AD校正完成，E_SENS/E_COG/H2/H4 pending。
- `figures/primary_effects.*`、`expression_effects.*`、`module_effects.*`、`sensitivity_effects.*`（English）。
- 本目录`QA.md`、`input_manifest.json`、`software_versions.json`、`logs/assembly_PASS.json`、`logs/final_QA.json`。
- SEA逐源行/feature身份审计、gene A、参考参数、每供体score、模型manifest与完整诊断在`../sea_hip_v1/`；原GSE五项仍位于`../gse48350_v1/`，未被改写。
'''
(C/'concise_result.md').write_text(combined_report)
sea_report=f'''# SEA-HIP v1 有界执行摘要

数据准备与输入hash接受通过；依据冻结v1执行36个新放行登记位置（35个SEA表达位置+S16）。实际新拟合**30**项；高kME identical3项不重拟合；交集低覆盖2项不拟合；S16复用完全相同的原P01、未重跑GSE五项。原始counts、preflight、frozen计划与GSE结果均保持hash不变。

## 主结果（最终多重校正见combined_v1）

{md(primary[primary.test_id.isin(['P02','P03'])])}

CA3实际gene A=100/120，SUB=106/120；主n29/31、df24/26。均通过覆盖/设计门槛。两项估计负向，但最终Holm后均未达0.05；不称验证成功。详见`../combined_v1/concise_result.md`，其中还保留GSE主检验和全部辅助族，不能只截取SUB raw P。

## 技术与实际覆盖

{md(maincoverage)}

严格全gene库量CPM→log2(CPM+1)，参照固定≥20核R，actual门控不按表型。库在来源文件中已按donor×supertype合并；本轮只把Sub1/2与DG1/2各自按donor相加，不能按n_libraries再乘一次。四target独立，不把同供体跨target当独立n。

## 明确NA与引用

S17/S18交集92/120不足96，未拟合；S10/11/12 identical引用父模型效应且公开P NA；S16交集21与原GSE完全相同，引用P01而无新GSE拟合。S19–S22及全部认知/H2/H4仍未执行。

所有基因/feature身份与原源行可追溯于`derived/candidate_feature_identity_audit.tsv`及`results/gene_selection.tsv`；S16独立源行视图在`results/S16_gene_selection_source_rows.tsv`。参考参数/供体score/模型manifest/实际表达门控/手工HC3核验均在results。诊断见`figures/all_fitted_model_diagnostics.pdf`（30页）及主模型PNG；无选择性删人或LOO P值。

本轮唯一多feature处理为对原preflight已映射同一当前HGNC的4个身份做正式HGNC历史名称核验后各source feature只相加一次、gene权重只计一次；未改映射或另择feature，详QA。临床认知字段没有用于评分或模型。
'''
(ROOT/'concise_result.md').write_text(sea_report)
qa=f'''# SEA-HIP v1执行QA与组合家族核验

## A. 冻结/输入与范围
- 原计划SHA256 `c8ecd963efff5ea4b078c223fddfc2a3ccbd9e2a6b1ff987aa1f74e17d5bc31a`；登记表SHA256 `d20d4220ec0b5c1cffa563c02b9b510ab33f95657fdfb132302be2f1ac7a43a1`。候选五文件MD5通过，xlsx `7fd7ac5097f9dad9977eea804edeb36e`。内部锁定时间早于GSE/SEA关联运行，并非外部registry声称。
- 来源准备QA报告/JSON及主NPZ/TSV的SHA256均与output_manifest一致；NPZ int64/CSR3402×36601、41,420,887 nnz、每行allgene和与metadata一致，总UMI6,196,354,669。接受原准备任务对1,038,862,619原层数值及cell/group指纹的PASS，不重读5GB H5AD或重跑cell aggregation。
- 原UMIs有效；X为ln(UP10K+1)，本任务从未用X或反变换X。feature顺序与前次36601行var每字段完全一致。源counts和所有冻结/GSE输入前后hash一致。
- 实际执行30个新OLS模型；另外3 identical、2 not-testable coverage、1已存在GSE效应引用。未执行认知/H2/H4、单基因关联、额外target/age阈值/归一化，不删供体。

## B. 精确身份/多feature审计：不悄悄加权或择一
- TSV布尔字段通过严格字符串True/False映射，非法值会拒绝；synthetic明确断言字符串False→False。
- 每冻结集合1:1 HGNC唯一；每source feature只属于一个选定当前HGNC，输出score每HGNC只出现一次。
- 原preflight确有4个当前HGNC各对应2个不同来源Ensembl feature。已在任何SEA关联前对HGNC approved/current/previous名称逐项核验，所有8个feature正是同一官方当前身份或其明确previous symbol（不是大写猜测）。**保持冻结preflight映射，先把这两个source count相加至该HGNC，再以一个gene权重参与score。** 36601-feature库分母不变、每source count只计一次；不丢掉旧feature、也不赋双份kME。没有同源/候选选择规则变更。

{md(merges)}

这些多feature身份只涉及DG-M3、DG-M7、CA3-M4；core与主METTL3程序不存在该多feature情形。完整映射/旧名/当前Ensembl/source index可审计。不是把不同HGNC身份相加；若审核发现不一致代码将停止。

## C. 供体/target与元数据
- counts已有两库供体H20.33.046的库内合并；本轮component求和从donor×supertype行开始，不按库数重复。4target×40 donor可追溯视图保留无核行的技术absence标志，UMI0时CPM/score为NA，不伪造可分析零表达供体。
- Sub=Sub1+Sub2、DG=DG1+DG2，纯CA3=CA3_1、CA1=CA1_1。每个target分别矩阵乘法合并；全部target供体n与前次400行QC逐项一致。
- 两个固定donor（含两库供体）×4target共8次手工全36601-feature求和与矩阵方法完全一致。count整数/int64非负；all-feature total与组件totalUMI严格一致。
- Clinical只取donor ID、age、sex、Braak、RIN、PMI，与原donor XLSX相应字段精确一致。读取Excel时只读这些列；没有读取/拟合cognitive endpoint。Braak显式去`Braak `前缀后0/I/II/III/IV/V/VI映射0..6，未识别值将拒绝；实际无缺失。PMI单位由已有来源字典确认为小时，RIN为1–10指标。
- 直接使用来源数值Age at Death（总体68–100）；现有字段字典未说明额外top-code规则，未恢复或猜测隐藏年龄，不将100自行改为100+，不宣称个体年龄精度超出来源。

## D. Full-background CPM、actual A、权重、参照锁定
- 分母是每donor×target所有36601 feature的count总和，含未映射/非候选feature；非候选从未从L剔除。候选多feature合并只改变候选分子，不改变L。
- CPM=1e6*gene_count/L；log2(CPM+1)一次。主R人数29/31/38/39，表达门槛CPM≥1至少6/7/8/8供体；actual gate没有年龄/诊断分组、DE或P。
- 逐gene参照SD(ddof1)>1e−8，技术NA不填0；生物count0保留。源member→1:1→技术→actual→variance→A及排除原因逐行保留。
- gene mu/sd、等权或固定kME分母、S mean/SD，以及PMI/RIN中心值均固定于主R；10/30、Braak、组成模型不重筛/重缩放。高kME和intersection只在主A上取预设子集、同R重算新score尺度。
- 每score保存后重读，R内Y均值≈0/样本SD≈1，每gene身份唯一、权重总和与固定分母一致。S17/S18虽可形成92gene技术score，因覆盖低于96未拟合模型，未因已有score偷跑P。

## E. Synthetic和独立手工验证
- before-model synthetic PASS：两库合并+Sub组件count求和、全背景CPM与错误候选分母区分、False字符串、固定≥20参照用于新增≥10供体、kME加权固定分母、非参照NA不做变分母平均。
- 五参数SEA设计known truth ageβ−0.6由独立构造误差恢复−0.6000000000000005；手工HC3 sandwich、t df25与statsmodels一致，design n<20正确拒绝。
- 30个实际新模型全部使用独立normal-equation beta、manual HC3 sandwich和t双侧P逐项核验（rtol1e−9,atol1e−11），无失败。统计helper沿已synthetic/手工核验的GSE固定OLS helper只读复用；没有运行GSE旧模型脚本。
- 主程序CA3某供体的raw count→fullL CPM→log2(CPM+1)→逐gene z→S→Y手工spotcheck一致（`logs/scalar_primary_score_spotcheck.json`）。另从已落盘参数对CA3-M2做独立逐gene kME加权重算，S与Y一致（`logs/manual_weighted_module_spotcheck.json`）。
- 所有complete-case/模型n/母群/年龄范围/sex/秩/残差df留表；未按结果改变协变量。Braak完整，因此没有必要增加同子集对照回归。

## F. 影响性与图
{md(pdiag)}
- 各模型记录VIF、hat>2p/n、Cook>4/n以及LOO beta范围与方向翻转数；没有LOO P或best-case selection。原scores/基因A/参照在每次LOO固定。不因标志删除人。
- `all_fitted_model_diagnostics.pdf`含30个实际模型的残差/杠杆/Cook/LOO四面板。主P02/P03与M05另存PNG；combined有主/表达/模块/敏感性效应PDF/PNG。点wise CI不伪装同时区间。

## G. 组合家族
{md(fam)}
- 只读取既有GSE五项和新结果；原GSE五项n/beta/SE/df/CI/rawP逐值精确一致，未重估。
- PRIMARY固定Holm3；E_EXPR固定BH11；E_MODULE固定BH7；E_AD固定BH1。手写排序/累积调整与statsmodels multipletests逐项一致；synthetic [0.01,0.04,1]验证Holm [0.03,0.08,1]，不能通过删除terminal槽位缩小m。
- 已知terminal NA（S17/S18覆盖失败、S10/11/12 identical）公开raw/adjusted NA，内部待用p1保留槽位；未知pending21项内部仍NA，没有当阴性。
- E_SENS19个位置已处理（14有P+5terminalNA），认知4个未知使整族adjusted NA/pending-family；不擅称全部完成。E_COG/H2/H4仍pending。
- 目前只有探索模块M05 BH q<0.05；主要3项均未达Holm支持，E_EXPR无BH q<0.05。不得用辅助/敏感性替主、升格target或推断m6A/酶活/认知改善。

## H. 复现、资源与文件保护
- 同一已存在GSE Python环境只读调用，NumPy {versions['numpy']}、pandas {versions['pandas']}、SciPy {versions['scipy']}、statsmodels {versions['statsmodels']}；Excel reader来自已有cohort环境只读fallback，没有改动旧环境或重下载数据。精确版本/输入/脚本hash见JSON。
- source NPZ始终稀疏，仅4个40×36601目标小矩阵临时dense；prepare峰值RSS见preparation_PASS，计分拟合峰值{execution['max_rss_kib']/1024:.1f} MiB。没有3402×36601多份dense或全细胞dense。
- preparation与analysis脚本各有成功标志拒绝重复运行。初次准备仅因日志目录/缺少Excel reader本地环境问题停止，修复读路径后才进入源矩阵处理；没有原始文件损坏或重新下载。
- 源SHA256再次核验；所有结果/参数/score/模型/LOO表已重读断言；PDF/PNG完整性核验。最终凭据：SEA与combined各自`logs/final_QA.json`。无Word。

命令（已成功执行，不要重跑模型）：
```bash
cd "{ROOT.parents[1]}"
# Read-only inspection; preserve completed models.
python run.py sea-prepare --dry-run
python run.py sea-models --dry-run
python run.py combined --dry-run
python run.py sea-report --dry-run
```
'''
(ROOT/'QA.md').write_text(qa);(C/'QA.md').write_text('# Combined human v1 QA\n\n组合家族的详细方法/QA见`../sea_hip_v1/QA.md`。本目录不重拟合任何模型，只读取两次执行结果并按冻结registry固定m组装。\n\n'+md(fam)+'\n\n主3项Holm、E_EXPR11/E_MODULE7/E_AD1 BH手写实现与statsmodels相同；原GSE五项数值逐值不变。NA原因及未知pending与internal-P1槽位明确分开。E_SENS不校正，认知/H2/H4未执行。全部62位置见results/registry_tracking_62.tsv。\n\nSource及result hash见input_manifest.json、logs/assembly_PASS.json、logs/final_QA.json。参数/基因/参考/样本身份追溯分别保留在GSE和SEA原执行目录。\n')
# Final file format + saved values + source protection; immutable scientific outputs.
for base in [ROOT,C]:
 for p in (base/'figures').glob('*.png'):
  with Image.open(p) as im:im.verify()
 for p in (base/'figures').glob('*.pdf'):assert p.read_bytes().startswith(b'%PDF-')
for path,hsh in json.loads((ROOT/'input_manifest.json').read_text())['files'].items():assert sha(Path(path))==hsh
assert sha(ROOT/'results/registered_results_long.tsv')==execution['results_SHA256']
for base in [ROOT,C]:
 files=[base/'QA.md',base/'concise_result.md',*sorted((base/'results').glob('*')),*sorted((base/'figures').glob('*')),base/'input_manifest.json',base/'software_versions.json'];record={'status':'PASS','completed_at_utc':datetime.now(timezone.utc).isoformat(),'source_files_unchanged':True,'GSE_five_results_unchanged_no_refit':True,'gene_reference_score_roundtrip':'PASS','model_flow_LOO_roundtrip':'PASS','actual_coverage_gates':'PASS','synthetic_and_manual_QA':'PASS','manual_weighted_module_spotcheck':'PASS','fixed_family_adjustment':'PASS','unknown_pending_not_imputed':True,'PNG_PDF_integrity':'PASS','report_script_SHA256':sha(Path(__file__)),'files':{str(p.relative_to(base)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in files}}
 dump(record,base/'logs/final_QA.json')
print('SEA and combined final QA PASS');print(md(primary));print(md(fam))
