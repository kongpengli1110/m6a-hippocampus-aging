# CANDIDATE_SPEC — mouse_candidates_v1（冻结版，只读）

冻结时间：2026-10-12 20:45 UTC+8。生成脚本：scripts/26_handoff_candidates.py（脚本副本在 02_scripts_env/scripts/）。
本目录的 6 个文件是 `handoff/` 的逐字节拷贝。打包时用 md5sum 实测，5 个数据文件全部与 `FROZEN_md5.txt` 一致；FROZEN_md5.txt 的第 6 行是说明文字，不是 md5 行。
**本说明只描述冻结文件里已有的内容。冻结文件和脚本里都没有的字段，一律标"缺项：未冻结"；没有决定过的事项标"未定"。v1 不做任何修改。**

## 1. 内容清单
| 文件 / sheet | 行数 | 内容 |
|---|---|---|
| core_candidates（`_core.csv`） | 23 | 核心基因：Tier1 7 个（Cacna1e、Slitrk1、Dlgap2、Adgrb2、Chst1、Rnf165/Ark2c、Grin2a），Tier2 16 个（Dgkg、Galnt16、Trank1、Tmem150c、Miat、Gria2、Lynx1、Epha5、Fbxo41、Cx3cl1、Ddn、Prrt3、Frrs1l、Zbtb18、Gabra1、Sv2b） |
| all_E1E2_candidates（只在 xlsx 中） | 71 | 同时满足 E1 和 E2 的全部基因。evidence_score 分布：2 分 17 个、3 分 31 个、4 分 16 个、5 分 4 个、6 分 3 个（Tier3 共 48 个） |
| module_genesets（`_modules.csv`） | 1356 行，1221 个不同基因（123 个基因同时属于 >1 个模块） | 7 个 m6A 富集模块的**全部成员**及其 kME：CA1-M9 32、CA3-M2 737、CA3-M4 131、DG-M3 110、DG-M7 36、SUB-M4 195、SUB-M6 115 |
| program_genesets（`_programs.csv`） | 718 | 3 个程序集：METTL3cKO_HIP_down 131、METTL3KO_Cb_m6Aloss 500、A1_core87_hypo 87 |
| README（`_README.csv`） | 18 项 | 规则说明，以及预设的人类假设 H1–H4 |

## 2. 选择规则（scripts/26 + README.csv）
- **E1（MeRIP）**：在 GSE198526（PNAS 2023 Dataset S17：FC≥1.2、padj≤0.05、作者跨重复一致的调用）的 CA1/CA3/DG 中，至少 2 个区随龄（3 vs 16 月）低甲基化。这是组织层面的证据，不是细胞特异的。
- **E2（细胞类型）**：Allen HPF（unbiased 文库，8 只成年 vs 15 只老年供体）中，CA1+CA3+DG+SUB 谷氨酸能类型的特异性份额（平均 CPM 占比）之和 ≥0.5。
- **evidence_score** = n_regions_hypo（2 或 3）+ 1[in_m6A_module] + 1[METTL3 靶（cKO_HIP_down 或 Cb_m6Aloss）或 YTHDF1 CLIP] + 1[Tegowski 皮层 Glut DM-RNA，任一方向] + 1[Allen 兴奋性 pseudobulk 年龄效应在 CA1/CA3/DG/SUB 任一中 padj<0.1]，满分 7。
- **Tier**：≥5 为 Tier1；=4 为 Tier2；其余为 Tier3。**Core** = score≥4 的基因按 score 降序、再按 max_kME 降序排列，取前 50 个。实际只有 23 个达到 score≥4，所以上限 50 没有起作用。
- **in_m6A_module** 的定义：`max_kME ≥ 0.6`，max_kME 取该基因在 7 个 m6A 富集模块中 kME_own 的最大值；满足则得 1 分。m6A 富集模块指 mod_enrichment 中，对 A1_HIP_hypo_union / METTL3cKO_HIP_down / METTL3KO_Cb_m6Aloss 任一的 BH q<0.01 的模块，即 mod_m6A_enriched_hubs.csv 中的 7 个。
  - 注意 Gria2：max_kME=0.573，在模块中但 kME 低于 0.6，所以 in_m6A_module=False，不得分。

## 3. 必读澄清
### 3.1 "kME≥0.6 counts" 的含义
- README 的 "member with kME≥0.6 counts" **只**是 core 证据计分中模块项的门槛（scripts/26：`inmod = max_kME>=0.6` → +1 分）。
- `_modules.csv` 导出的是 7 个模块的**完整成员**，没有按 kME 过滤（scripts/26 对 `mm[mm.module==m]` 的每一行都写出）。
- **实测**（handoff_pkg/_work/chk1.py 复用）：共 **1356 行，其中 kME<0.6 的 144 行**，kME 缺失 0 行。与接手方报告的两个数字一致。
  - 各模块 kME<0.6 的行数：CA3-M2 110、SUB-M4 24、SUB-M6 7、CA1-M9 3，其余为 0。kME≥0.6 共 1212 行。
- **人类投射时用全部成员还是只用 kME≥0.6 子集：从未定过 → 未定。** 由接手方在看任何人类关联结果之前前瞻冻结（建议写入 v1 之外的分析计划并记录时间戳），v1 不改。

### 3.2 "95.7%" 的含义
- 只指 **core 中 22/23 = 95.7% 为 MGI 1:1 同源**（ortholog_type=='1:1'）。实测：1:1 有 22 个，1:2 有 1 个（Lynx1）。
- 这个比例**不**适用于模块和程序集。模块行中 1:1 占 1254/1356，NaN（无 MGI 人类同源）83 行，1:2 16 行，1:3/1:5/1:7 各 1 行。程序集：METTL3cKO_HIP_down 1:1 120/131（NaN 9），Cb_m6Aloss 467/500（NaN 19），core87 85/87（NaN 1）。
- **Lynx1**：MGI 同源组 DB Class Key 52333509，人类对应 LYNX1（Entrez 66004）和 LYNX1-SLURP2（Entrez 111188157，HGNC:52291）两个，属一对多。README 写明"1:1 recommended for tests; others flagged"。**主分析不纳入 Lynx1**（这是 README 规则的直接应用；Lynx1 仍保留在 v1 文件里，未删除）。
- **Rnf165 → Ark2c / ARK2C**：A1 数据中用的是旧符号 Rnf165。scripts/26 通过 MGI MRK_List2.rpt 同义词表映射到现行符号 Ark2c（MGI:2444521；同义词包括 LOC225743|2900024M11Rik|G630064H08Rik|Ark2c|Rnf165），再映射到人类 ARK2C（HGNC:31696，Entrez 494470，1:1，DB Class Key 52323930）。以上都已在 data/ref 文件中实查确认。人类数据里请用 **HGNC:31696 / Entrez 494470** 追溯，不要按符号 "RNF165" 检索。
- 模块中另有 34 个、程序集中另有 3 个旧符号经同义词映射（例如 March1→Marchf1、Sept9→Septin9、5830416P10Rik→Mirt1），对照见 `mgi_current_symbol` 列。

### 3.3 模块 age_p 与程序证据
- 7 个模块的 module_age_p（供体层面 OLS，`me~C(age)+sex(+roi)`）：DG-M3 0.172（最小）、CA3-M4 0.235、CA3-M2 0.243、SUB-M6 0.252、DG-M7 0.361、SUB-M4 0.484、CA1-M9 0.539。**全部 >0.17**（实测最小值 0.1715）；37 个模块中最小 q=0.71。
- 单个模块的年龄效应不显著，**不能用模块来替代程序层面的证据**（METTL3cKO_HIP_down 在 CA3/SUB 的基因集检验，见 05b）。但两者**并不矛盾**：两个检验的单位和功效完全不同。模块检验以供体为单位，CA3 只有 5 vs 4、SUB 4 vs 7、DG 5 vs 5；程序检验以基因为单位，做的是竞争性检验。模块组水平的弱证据见 07 检验卡 C3。

## 4. 逐项字段说明
| 项目 | v1 中的状态 |
|---|---|
| 基因清单 | 已冻结：core 23 / E1E2 71 / 模块 1356 行 / 程序 718 行 |
| 模块清单 | 已冻结：7 个模块及成员（含 kME） |
| **权重** | core：给出 evidence_score（整数 2–7）和 tier，但**没有规定**人类评分时是否按分数加权 → **缺项：未冻结**（人类评分用等权还是加权）。模块：只给出 kME，**没有定义**是否以 kME 加权 → **缺项：未冻结**。程序集：无权重。README H1 写的是 "core-list score (mean z)"，即隐含等权的均值 z，但没有写成明确的加权规则 |
| **方向** | core：`predicted_direction_human_aging_expression` = "down in aged hippocampus/cortex excitatory neurons (or bulk), relative to young"；`predicted_direction_m6A_if_measurable` = "lower m6A with age (hypomethylation), conditional evidence only"。H1–H3 为单侧，只允许按所述方向（H1 core 评分随龄降低；H2 富集于兴奋性神经元特异基因；H3 METTL3cKO_HIP_down 同源集随龄活性降低）；H4 遗传学部分无方向先验，属探索性。**模块方向**：模块文件只有 module_age_beta_SD（小鼠观测值），**没有给人类预测方向** → 缺项：未冻结。program 中 Cb_m6Aloss / core87 的人类方向 → 缺项：未冻结（H3 只覆盖 METTL3cKO_HIP_down） |
| MGI symbol | 已给出：mouse_symbol（原始）+ mgi_current_symbol（现行） |
| **MGI ID** | **缺项：未冻结**（没有导出 MGI:xxxx 列；可以从 HOM 文件按符号回查，但 v1 里没有） |
| **Ensembl mouse ID** | **缺项：未冻结**（core/模块/程序文件都没有 Ensembl 列） |
| HGNC | 已给出：human_symbol + HGNC_ID（一对多时用分号连接） |
| Entrez（人） | 已给出：human_Entrez（分号连接）。**小鼠 Entrez → 缺项：未冻结** |
| **Ensembl human ID** | **缺项：未冻结** |
| 细胞类型 | core：top_celltype（Allen 特异性最高的兴奋性类型）+ excitatory_specificity；模块：celltype（CA1/CA3/DG/SUB）。**人类细胞类型对应规则**（Allen 小鼠亚类 → SEA-AD/人类注释）→ **缺项：未冻结** |
| 同源映射规则及版本 | 已冻结：MGI HOM_MouseHumanSequence.rpt（下载于 2026-10-12，md5 89dc473ae85b11204fec5318bbf9ed81）；旧符号经 MRK_List2.rpt（md5 30e1088b64f45d2fd45c6befb3867e36）同义词映射到现行 MGI 符号。ortholog_type = 该 DB Class Key 下小鼠基因数 : 人类基因数。空白 = MGI 中没有人类同源。打包时实测两个文件的 md5 与 README 一致。MGI 报告文件内部没有版本号字段，只能靠下载日期和 md5 标识 |
| 一对多同源的处理 | README："1:1 recommended for tests; others flagged"。一对多时人类集合是否展开为全部人类基因（如 LYNX1;LYNX1-SLURP2）→ **未定**（只规定了推荐 1:1） |
| 模块用全部成员还是 kME≥0.6 子集 | **未定**（见 3.1） |
| 无同源基因（NaN）| 在人类分析中自然丢失；是否报告丢失比例 → 没有规定（建议报告） |
| 背景集（人类富集/EWCE 用） | H2 写的是 "expression-matched background"，具体背景基因集 → **缺项：未冻结** |
| 统计检验细节（人类） | H1–H3 只给了方向和单侧；具体模型、协变量、多重校正 → 按 step8_human_feasibility.md 的"分析口径"由接手方前瞻定，v1 中**未冻结** |

## 5. 不得声称
README 原文：不得把任何人类表达结果描述为"验证 m6A 变化"。小鼠 sc/bulk RNA 只支持表达或细胞状态层面的结论；bulk MeRIP 是组织层面的证据。
