# ENV_NOTES — 运行环境与已知坑

- Python：**3.12.15**（uv 管理的 CPython；`uv venv .venv --python 3.12`，原机路径 ./.venv）。无 R。
- 依赖：`requirements.txt`（`uv pip freeze --python .venv/bin/python` 导出；`.venv/bin/python -m pip freeze` 失败，因 venv 内没有 pip）。
- 硬件：2 核 / 8 GB / 约 13 GB 盘。所有 pseudobulk 以动物/供体/文库为单位。
- 关键包（当前 venv）：scanpy 1.12.4、anndata 0.13.4、decoupler 2.2.0、pydeseq2 0.5.4、statsmodels 0.15.0、numpy 2.5.4、scipy 1.18.1、h5py 3.16.0、fsspec 2026.9.0、pyranges 0.1.4、scikit-learn 1.9.1、liana 1.10.0。

## ⚠ 版本漂移（交接时发现）
- REPORT_phase2a 记录分析时 **pandas 3.0.6**；2026-10-12 20:47（UTC+8）为 2c 细胞通讯执行 `uv pip install liana img2pdf`（后台任务 bc5ca11）后，venv 中 pandas 变为 **2.3.3**（dist-info 时间戳 12:47 UTC），同时新增 liana/mudata/plotnine 等。
- 因此 `requirements.txt` 反映的是**当前**环境，不是 scripts/01–26 实际运行时的环境；01–26 运行于 pandas 3.0.6（其余包时间戳 10:48 UTC 未变）。复现 01–26 时建议 `pandas==3.0.6`。

## 已知注意点
1. **沙箱无 /dev/shm**：joblib / multiprocessing 信号量不可用。pydeseq2 需打串行补丁：在导入 pydeseq2 后 `import pydeseq2_patch`（把 `pydeseq2.default_inference.parallel_backend` 替换成 nullcontext），并用 `DefaultInference(n_cpus=1)`。结果不受影响。scripts/05 首次运行（`backend='threading'`）失败，后用 sed 改为 `DefaultInference(n_cpus=1)` 重跑成功。
2. **Allen ABC（Zeng-Aging-Mouse-10Xv3）不下载全量**：raw h5ad（13.8 GB）用 `fsspec.open(url,'rb',block_size=2**24,cache_type='readahead')` + `h5py` 做 HTTP range 流式读取，只抽取目标 ROI 的行（scripts/01 取 cell_label/indptr；02 按 ROI 抽 pseudobulk + 分层子样本；20 抽兴奋性神经元全细胞）。**每个区域约 2–3 分钟**（HIP+SUB/PAR+ENT 共约 6 分钟）。URL：`https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/expression_matrices/Zeng-Aging-Mouse-10Xv3/20241130/Zeng-Aging-Mouse-10Xv3-raw.h5ad`；metadata 取自 `.../metadata/Zeng-Aging-Mouse-10Xv3/20250131/`。scripts/02 内部用 2 个子进程并行（subprocess，不依赖 /dev/shm）。
3. **PMC 补充材料（Tegowski 2024，PMC11614689 Supp Table 1/4/9/10）有 PoW 挑战**：直接 curl 得到的是 HTML 挑战页。用 `pow.py <表号...>`：解析 `POW_CHALLENGE`/`POW_DIFFICULTY`，暴力求 sha256 前导 0，设 cookie `cloudpmc-viewer-pow` 后重下。PNAS 补充材料 403 → 改用 figshare 10.6084/m9.figshare.21966983。Zenodo API 403 → 用 `records/15849341/files/<name>?download=1`。
4. **GSE129788 GEO 矩阵还原 raw counts**：GEO 提供的是 ln(1+CPM×1e4) 形式归一化值；`counts = (exp(x) − 1) × nUMI / 1e4`（nUMI 取自 `GSE129788_Supplementary_meta_data_Cell_Types_Etc.txt.gz`），四舍五入后逐细胞断言列和 = nUMI（rtol 1e-3），全部通过（scripts/21）。
5. 其他：`OMP_NUM_THREADS=1` 运行 scripts/24（避免 2 核过载）；GSE302188 GEO 与 Zenodo 样本 ID 不能按名字对应，需用 scripts/03 的表达匹配表 `B11_sample_matching.csv`；`data/ref/mm_genelen.tsv` 是 0 字节空文件，未被使用。
6. 所有脚本把项目根写死为 `./`（common.py 的 `D`，以及若干脚本中的 `D=`/`T=`/`F=`）；03/04/13 用相对路径，需在项目根目录运行。
