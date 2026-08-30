# 合并流程文档：01–05 → `Rversion/00test/merged_pipeline_01_05.R`

> 本文件记录 01–05 五个分析流程**合并为单脚本**的逻辑、流程、文件结构，以及如何把改动**同步回**原有的 01–05 各脚本。
> 同步回原脚本**本次不执行**，仅在本文档化，供后续维护参考。

---

## 一、合并逻辑

| 项 | 说明 |
|---|---|
| 合并来源 | `Rversion/01_initiation/`、`02_quality_control/`、`03_extract_cd45/`、`04_integration_clustering/`、`05_celltype_annotation/` 五个目录下的主脚本 |
| 合并产物 | 单一脚本 `merged_pipeline_01_05.R`（逻辑原样保留，仅重写路径/加段分隔/加断点开关） |
| 样本范围 | **仅 GSE203115 子集**（data/ 下 `GSE203115_1/2/3` 三个样本目录），见 `SAMPLE_PREFIX <- "GSE203115"` |
| 处理逻辑 | **仅主流程**，不含 exclude（剔除 5 样本）分支 |
| 运行形态 | 一个脚本按顺序跑完 01→05；段间 rds **走磁盘传递**（与原有模块化设计一致） |
| 关键断点 | SECTION 03 末尾 `saveRDS(03_CD45_positive.rds)`，作为整条流程的存盘点 |
| 续跑开关 | 顶部 `RESUME_FROM`：可设为 `""/"01"/"02"/"03"/"04"/"05"`，跳过前段（前段 rds 需已存在） |
| 落点 | 全部生成文件统一在 `Rversion/00test/` |

**为何这样合并**：原 01–05 各自独立、靠"上一步 rds"串联。合并成单脚本后可一次性复现整条分析链，且因 rds 仍落盘，可随时从任意段续跑；GSE203115 子集 + `00test` 沙盒定位便于快速验证流程结构。

---

## 二、流程说明

```
data/GSE203115_1|2|3 + data/样本分组信息.xlsx
  │  SECTION 01  (Read10X → CreateSeuratObject → merge)
  ▼
01_seurat_combined.rds                → output/01_initiation/
  │  SECTION 02  (QC 指标 / scDblFinder / 固定阈值过滤)
  ▼
02_seurat_qc.rds                     → output/02_quality_control/
  │  SECTION 03  (PTPRC counts>0 判 CD45+ → subset)
  ▼
03_CD45_positive.rds   ★关键断点      → output/03_extract_cd45/
  │  SECTION 04  (Normalize → HVG → Scale → PCA → Harmony → UMAP/tSNE → 多分辨率聚类)
  ▼
04_CD45_integrated.rds               → output/04_integration_clustering/
  │  SECTION 05  (marker 打分 → 注释/生成模板 → 输出图与 05_annotated.rds)
  ▼
05_annotated.rds                     → output/05_celltype_annotation/
```

- **段间数据流**：每段开头的 `INPUT_RDS` 指向**上一段落盘的 rds**（例如 SECTION 03 读 `RDS_02`，SECTION 04 读 `RDS_03`）。因此即使单进程内内存中 `obj` 被覆盖，也始终以磁盘 rds 为权威输入，天然支持 `RESUME_FROM`。
- **03 断点语义**：无论从头还是从 `RESUME_FROM="04"` 续跑，SECTION 04 都读取已落盘的 `03_CD45_positive.rds`，该 rds 即整条流程的"中段检查点"。
- **SECTION 05 双分支**：存在 `celltype_map.csv` → 应用注释；不存在 → 生成 `celltype_map_template.csv` 供人工填写后重跑（与原 05 行为一致）。
- **范围提示**：GSE203115 子集仅 3 个样本（ypN0×2 + ypN+×1），下游 04/05 仅供流程结构验证，**不可作生物学结论**。

---

## 三、文件结构

```
Rversion/00test/
├── merged_pipeline_01_05.R     # 合并后的单脚本（本说明的对应代码）
├── README_merged.md            # 本文件
├── manual_markers.R            # 从 05 复制的经典 marker 字典（已内联进脚本，留作参照）
├── celltype_map.csv            # 从 05 复制的注释映射表（21 行，celltype 列待核定）
└── output/
    ├── 01_initiation/          # 01_seurat_combined.rds
    ├── 02_quality_control/     # 02_seurat_qc.rds + csv/pdf
    ├── 03_extract_cd45/        # 03_CD45_positive.rds ★断点 + csv/pdf
    ├── 04_integration_clustering/  # 04_CD45_integrated.rds + csv/pdf
    └── 05_celltype_annotation/ # 05_annotated.rds + csv/pdf
```

### 输入/输出 rds 映射表

| 段 | 输入 rds（来自上一段落盘） | 输出 rds |
|---|---|---|
| 01 | data/ （GSE203115 样本，经 xlsx 分组） | `output/01_initiation/01_seurat_combined.rds` |
| 02 | `output/01_initiation/01_seurat_combined.rds` | `output/02_quality_control/02_seurat_qc.rds` |
| 03 | `output/02_quality_control/02_seurat_qc.rds` | `output/03_extract_cd45/03_CD45_positive.rds` ★ |
| 04 | `output/03_extract_cd45/03_CD45_positive.rds` | `output/04_integration_clustering/04_CD45_integrated.rds` |
| 05 | `output/04_integration_clustering/04_CD45_integrated.rds` | `output/05_celltype_annotation/05_annotated.rds` |

> 各段输出目录独立，可避免 02/03 同名 pdf（如 `01_QC_violin_before_filter.pdf`）相互覆盖。

---

## 四、同步回 01–05 的映射（供后续维护，本次不执行）

合并脚本按"统一 CONFIG + 5 段"重构，改动若需回写原 01–05，按下表定位对应位置。

| 合并段 | 对应原文件 | 原文件关键行 | 合并版差异点 |
|---|---|---|---|
| 全局 CONFIG | （各原文件顶部 `## 0. 配置`） | 01:34-39 / 02:56-86 / 03:55-68 / 04:58-103 / 05:87-232 | 路径改为 `OUT_0X`/`RDS_0X`；新增 `RESUME_FROM`、`SAMPLE_PREFIX`；`DATA_DIR` 修正为两级 `..`；`ensure_pkg`/`save_pdf`/`run_from` 提为全局 |
| SECTION 01 | `01_initiation.R` | 21-171 | `DATA_DIR` 修正两级 `..`；样本发现后追加 `startsWith(...,"GSE203115")` 过滤；`OUT_DIR/OUT_RDS` 指向 `OUT_01/RDS_01` |
| SECTION 02 | `02_quality_control.R` | 28-505 | `INPUT_RDS`→`RDS_01`、`OUT_DIR`→`OUT_02`；`check_meta_consistency` 等函数原样保留 |
| SECTION 03 | `03_extract_cd45.R` | 26-341 | `INPUT_RDS`→`RDS_02`、`OUT_DIR`→`OUT_03`；`saveRDS(03_CD45_positive.rds)` 在绘图前（第 243 行原位置），即关键断点 |
| SECTION 04 | `04_integration_clustering.R` | 29-568 | `INPUT_RDS`→`RDS_03`、`OUT_DIR`→`OUT_04`；**保留三处硬编码**：`DIMS=1:30`（原 283）、harmony `theta=5/max.iter.harmony=50/sigma=0.1`（原 295-304）、`DEFAULT_RESOLUTION<-0.4`（原 401） |
| SECTION 05 | `05_celltype_annotation.R` | 57-772 | `INPUT_RDS`→`RDS_04`、`OUT_DIR`→`OUT_05`；原 `source(manual_markers.R)`（原 226-232）改为**内联 `MANUAL_MARKERS`**；13 个辅助函数原样保留 |

**回写建议**：在合并脚本中修改某段逻辑后，把该段整体（含段横幅注释）对照上表复制到对应原文件，并将路径变量还原为原文件的 `file.path(script_dir, "..", ...)` 形式（`DATA_DIR` 在原 01 中沿用其原有写法即可——原 01 的 bug 仅影响从子目录运行时，原文件历史产物已存在故保持原样更稳妥；若希望原 01 也可重跑，则一并修正其 `DATA_DIR` 为两级 `..`）。

---

## 五、已知修正与保留点

### ✅ 已修正
- **01 `DATA_DIR` 路径 bug**：原 `01_initiation.R:34` 写 `file.path(script_dir, "..", "data")`，少一级 `..`（脚本被移入子目录后失效，指向不存在的 `Rversion/data`）。合并版统一为 `file.path(script_dir, "..", "..", "data")`，正确指向 `A:\Workbuddy\singlecell\data`。

### ⚠️ 刻意保留（不可"上提 CONFIG"，否则结果与原始流程不一致）
- **04 三处硬编码**：
  - `DIMS = 1:30`（无条件覆盖自动选 PC 逻辑，原 04:283）
  - harmony `theta = 5 / max.iter.harmony = 50 / sigma = 0.1`（覆盖 `HARMONY_THETA = NULL` 分支，原 04:295-304）
  - `DEFAULT_RESOLUTION <- 0.4`（注释写 0.5，实为 0.4，原 04:401）
- **ypN0/ypN+ 配色** `#4DBBD5`/`#E64B35`（03:297、04:452、05 的 `HM_COLOR`），反引号 `` `ypN+` `` 写法保留。
- **各段 JoinLayers 兜底**（02:192 / 03:138 / 04:172）保留，保证 `LayerData(...,"counts")` 取数正确。

### ✅ 新增
- 样本发现限定 GSE203115（`startsWith(sample_dirs, "GSE203115")`）。
- `RESUME_FROM` 断点续跑开关 + `run_from()` 判定函数。
- `manual_markers.R` 内联为 `MANUAL_MARKERS`，去除对外部脚本的 `source` 依赖。

---

## 六、运行方式

> R 不在系统 PATH，必须用能绝对路径的 `Rscript.exe` 调用；不要用 `Rscript -e "..."` 内联（本沙箱会 segfault）。

```bash
# 从头跑完整 01->05（默认 RESUME_FROM=""）
"C:/Program Files/R/R-4.4.3/bin/Rscript.exe" "A:/Workbuddy/singlecell/Rversion/00test/merged_pipeline_01_05.R"

# 仅从 03 断点之后续跑（需 03_CD45_positive.rds 已存在）
#   修改脚本顶部 RESUME_FROM <- "03" 后运行同上命令
```

**资源提示**：机器 28GB 内存，单跑 04/05 约占用 8–13GB（加载 ~2.4GB rds 时峰值高），**不要并发跑两个 R**。

**产物自查**：跑完后以 `output/03_extract_cd45/03_CD45_positive.rds` 等实际落盘文件为准，勿仅凭日志判断。
