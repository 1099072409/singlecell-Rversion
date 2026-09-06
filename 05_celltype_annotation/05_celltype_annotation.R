# =====================================================================
# 05_celltype_annotation.R
# 单细胞数据分析流程 - 第五步（CD45+ 免疫细胞类型注释）【手动注释版】
#
# 本脚本的功能：
#   1) 读取第四步整合聚类后的对象（04_CD45_integrated.rds）
#   2) 【唯一注释来源】使用脚本顶部硬编码的手动注释表 MANUAL_CELLTYPE_MAP
#      （cluster -> 英文细胞类型），逐簇写入 celltype 列
#   3) 一致性校验：对象里的 cluster 是否全部被手动表覆盖、手动表是否有重复/多余簇
#   4) 输出：注释后 UMAP（整体 + 按分组）、细胞类型组成统计与组成图、
#      marker 验证图（UMAP 表达网格 + 手动经典 marker 点图）、注释后的 rds
#
# ---------------------------------------------------------------------
# 本版主要改动（相对上一版）：
#   * 【删除】AddModuleScore 免疫 marker 打分（Score_* 系列列全部不再生成）
#   * 【删除】01_celltype_score_by_cluster.csv 打分矩阵、02_score_heatmap.pdf 打分热图
#   * 【删除】top_hit 自动注释建议、celltype_map_template.csv 模板生成、
#     celltype_map.csv 的读取逻辑（不再需要"先看结果再填表"，注释由脚本内手写表决定）
#   * 【删除】derive_celltype_detailed()：不再按 CD8T/CD4T 打分把 T 细胞自动拆成
#     CD8T/CD4T（那是自动注释，本版只保留你手写的类型名），celltype_detailed 列已移除
#   * 【删除】FindAllMarkers 自动 topN 差异基因，以及 06_marker_dotplot_auto /
#     _combined 两份点图（自动找 marker 属于自动注释信息）
#   * 【新增】MANUAL_CELLTYPE_MAP：10 类手写注释（T / gdT / NK / B / Plasma /
#     Myeloid / Mast / Proliferating / Endothelial / Fibroblast），覆盖 cluster 0~22
#   * 【保留+改名】注释 UMAP、按分组 UMAP、组成统计、组成图；marker 验证仅保留
#     手动经典 marker 点图（manual_markers.R，已按新 10 类重写）
#   * 所有可调参数集中到 ## 0. 配置，并对每个默认值注释"为什么这样设"
#   * 全脚本逐行中文注释
#
# 关于类型名的写法：一律使用 ASCII 字符串（如 T / NK / B / Plasma / Myeloid / Mast /
# T_proliferating / B_proliferating / Epithelial / Fibroblast / Endothelial）。
# 原因：Windows 中文环境的本机编码是 GBK（CP936），不含希腊字形；若类型名含非 ASCII 字符
# （如希腊字母），即便写成 Unicode 转义，R 在 GBK locale 下也会解析为 "<U+XXXX>" 形式，
# 图例和 CSV 里都会显示成这种转义文本，既不专业又影响下游读取。故统一使用 ASCII 名称。
#
# 手动注释表（cluster -> celltype，2026-09-02 由用户确认；对象为 20 簇，编号 0~19，已全覆盖）：
#     9                       -> NK
#     3                       -> B
#     14                      -> Plasma
#     7, 8, 11, 18            -> Myeloid
#     12                      -> Mast
#     10                      -> T_proliferating
#     16                      -> B_proliferating
#     19                      -> Epithelial    （上皮细胞）
#     15                      -> Fibroblast    （成纤维细胞）
#     17                      -> Endothelial   （内皮细胞）
#     0, 1, 2, 4, 5, 6, 13    -> T
#
# 用法：
#   在 RStudio 中打开本文件并点击 Source，或在命令行执行：
#     "C:/Program Files/R/R-4.4.3/bin/Rscript.exe" 05_celltype_annotation.R
#   说明：所有路径按脚本自身所在目录自动推算，可从任意工作目录运行。
#   若要修改注释结果：直接改 ## 0.2 节的 MANUAL_CELLTYPE_MAP，重跑即可，
#   不需要（也不会）读取任何外部 csv。
#
# 输出文件清单（位于 output/ 目录）：
#   01_celltype_map_applied.csv   - 实际应用的注释表（cluster | celltype | n_cells）
#   02_annotation_UMAP.pdf        - 注释后的 UMAP（按细胞类型着色 + 标签）
#   （03_annotation_UMAP_by_group.pdf 已删除：与 13_* 组级别堆叠图语义重复）
#   （原 04_celltype_composition.pdf 已于 2026-09 删除：与 13_* 完全重复，此处不再补号）
#   05_celltype_summary.csv       - 总体细胞类型计数与占比
#   06_celltype_by_group.csv      - 分组 × 细胞类型计数
#   07_marker_UMAP.pdf            - QC 用：核心 marker 的 UMAP 表达网格
#   08_marker_dotplot_manual.pdf  - 验证用：手写经典 marker 点图（按新 10 类）
#   09_stacked_composition_percent.pdf/.png - 纵向堆叠条形图（比例，样本在 x 轴，标注所属分组）
#   10_stacked_composition_count.pdf/.png   - 纵向堆叠条形图（绝对细胞数，标注所属分组）
#   11_proportion_boxplot_by_group.pdf/.png  - 分组箱式图（每样本细胞类型比例 + Wilcoxon 标注）
#   12_ROE_observed_counts.csv    - RO/E 观察计数矩阵（行=样本，列=细胞类型）
#   12_ROE_ratio_matrix.csv       - RO/E 比值矩阵
#   12_ROE_heatmap.pdf/.png       - RO/E 富集热图（标注数值）
#   05_annotated.rds              - 注释后的 Seurat 对象（meta 含 celltype 列）
#
# 注意：output/ 下 2026-08-27 的旧产物（01_celltype_score_by_cluster.csv、
#   02_score_heatmap.pdf、06_marker_dotplot_*.pdf、celltype_map_template.csv 等）
#   来自上一版脚本且基于旧的 21 簇对象，本版不再更新它们，可自行归档删除。
# =====================================================================


## ---- 0. 配置（所有可调参数集中于此，逐条注释默认值设定原因）----

# -------------------------------------------------------------------
# 0.1 路径与输入/输出
# -------------------------------------------------------------------

# 自动定位脚本所在目录（兼容 Rscript 与 RStudio 两种运行方式）：
#   - RStudio 中优先通过 rstudioapi 获取当前文档路径
#   - 命令行 Rscript 时通过 sys.frames()/--file= 参数获取脚本路径
#   - 全部失败则退回当前工作目录
# 原因：脚本需"可独立运行"，不能依赖用户手动 setwd()，故用脚本自身位置推算相对路径。
script_dir <- tryCatch({
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
    dirname(rstudioapi::getActiveDocumentContext()$path)            # 方法1：RStudio 当前文档路径
  } else {
    f <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
    if (is.null(f) || !nzchar(f)) {
      arg <- commandArgs(trailingOnly = FALSE)
      f <- sub("^--file=", "", grep("^--file=", arg, value = TRUE)[1]) # 方法2：命令行 --file=
    }
    if (is.null(f) || !nzchar(f) || is.na(f)) {
      getwd()                                                         # 方法3：退回工作目录
    } else {
      dirname(normalizePath(f, mustWork = FALSE))
    }
  }
}, error = function(e) getwd())

cat(sprintf("脚本目录: %s\n", script_dir))

# 第四步整合聚类对象（CD45+ 免疫细胞）的路径。
# 默认：脚本目录上一级的 04_integration_clustering/output/04_CD45_integrated.rds
# 原因：项目目录约定 05 步依赖 04 步产物；相对路径保证可移植。
INPUT_RDS <- file.path(script_dir, "..", "04_integration_clustering", "output", "04_CD45_integrated.rds")

# 【追加 · 2026-09】兼容 04 步产物被归档到子目录的情况。
# 原因：04_integration_clustering/output/ 下的产物在 2026-09-02 被整理进了
#       output/包含HRS/ 子目录，上面的默认路径随之失效，脚本会直接
#       "找不到输入 rds" 而中断。这里按候选顺序探测，两种目录结构都能跑。
# 做法：只做"存在性探测 + 取第一个命中的文件"，不改动后续任何读取逻辑。
# 注意：若将来 04 步又换了别的子目录，把新路径追加进下面的候选向量即可。
# 做法：先按默认路径找；找不到就在 04 步的 output/ 下递归搜同名文件。
# 关键：候选路径由 list.files() 从文件系统回读，而不是在代码里硬编码中文目录名。
# 原因：Windows 下 R 以 GBK 解析源码里的 UTF-8 字面量，硬编码的中文路径
#       （如 04 步整理出的"包含HRS"子目录）在 file.exists() 里恒为 FALSE，
#       只能靠文件系统回读拿到编码正确的路径。
rds_root    <- file.path(script_dir, "..", "04_integration_clustering", "output")
rds_scanned <- list.files(rds_root, pattern = "^04_CD45_integrated\\.rds$",
                          recursive = TRUE, full.names = TRUE)
rds_scanned <- rds_scanned[order(nchar(rds_scanned))]        # 层级浅的优先
rds_candidates <- unique(c(INPUT_RDS, rds_scanned))          # 原默认路径排在最先
rds_hit <- rds_candidates[file.exists(rds_candidates)]
if (length(rds_hit) == 0) {
  stop("找不到输入 rds，已尝试以下路径:\n  ", paste(rds_candidates, collapse = "\n  "))
}
INPUT_RDS <- rds_hit[1]
if (!identical(INPUT_RDS, file.path(rds_root, "04_CD45_integrated.rds"))) {
  cat(sprintf("提示: 输入 rds 不在默认路径，已自动定位到: %s\n", INPUT_RDS))
}
rm(rds_root, rds_scanned, rds_candidates, rds_hit)

# 本步输出目录（自动创建）。原因：所有中间结果与图表统一落盘，便于复查。
OUT_DIR   <- file.path(script_dir, "output")

# -------------------------------------------------------------------
# 0.2 手动注释表（本脚本唯一的注释来源）
# -------------------------------------------------------------------

# 用哪一列作为聚类分组（04 步的主 cluster 列）。
# 原因：项目约定 04 步产物主聚类列名为 "cluster"。
# 注意：簇数随 04 步的聚类参数而变（2026-09-02 17:39 重跑 04 后由 23 簇变为 20 簇，编号 0~19）。
#   重跑 04 后务必核对下面的 MANUAL_CELLTYPE_MAP 是否仍覆盖全部 cluster，
#   未覆盖的 cluster 会被自动标为 "Unknown"（运行日志会给出提示）。
CLUSTER_COL <- "cluster"

# 手动注释：每个细胞类型对应哪些 cluster。
# 写法：类型名 = c(cluster 编号, ...)，类型名即最终出现在 UMAP/统计表中的英文标签。
# 原因：注释完全由人工判定（结合 marker 表达与文献），脚本不做任何自动推断，
#       故不再读取 celltype_map.csv，也不生成 top_hit 建议，改这里即改结果。
# 重要：一个 cluster 只能出现在一个类型里（脚本启动时会自动查重并报错）。
MANUAL_CELLTYPE_MAP <- list(
  "T"               = c(0, 1, 2, 4, 5, 6, 13),  # T 细胞
  "NK"              = c(9),                      # NK 细胞
  "B"               = c(3),                      # B 细胞
  "Plasma"          = c(14),                     # 浆细胞
  "Myeloid"         = c(7, 8, 11, 18),           # 髓系细胞（单核/巨噬/DC/中性粒等合并）
  "Mast"            = c(12),                     # 肥大细胞
  "T_proliferating" = c(10),                     # 增殖型 T 细胞（增殖状态，非独立谱系）
  "B_proliferating" = c(16),                     # 增殖型 B 细胞
  "Epithelial"      = c(19),                     # 上皮细胞
  "Fibroblast"      = c(15),                     # 成纤维细胞
  "Endothelial"     = c(17)                      # 内皮细胞
)
# 注意：当前对象共 20 个 cluster（0~19），上面的映射已完整覆盖，因此不会产生 "Unknown"。
#   若 04 步重跑后簇数变化，请同步更新本表；未覆盖的 cluster 会被自动标为 "Unknown"
#   （运行日志会打印"未注释（Unknown）"提示，可据此发现遗漏）。

# 细胞类型在图例/坐标轴中的排列顺序（按免疫谱系：T → NK → B → 浆 → 髓系 → 肥大 →
# 增殖（T/B）→ 基质（上皮/内皮/成纤维））。
# 原因：factor 的 levels 即绘图顺序，固定顺序保证多次运行配色与排列一致。
# 注意：未被 MANUAL_CELLTYPE_MAP 覆盖的 cluster 会被标为 "Unknown" 并追加到末尾。
TYPE_ORDER <- c("T", "NK", "B", "Plasma", "Myeloid", "Mast",
                "T_proliferating", "B_proliferating",
                "Epithelial", "Fibroblast", "Endothelial", "Unknown")

# 样本列名与分组列名（04 步产物中已有 sample / group / sample_id 三列）。
# 原因：堆叠图、箱式图、RO/E 都按"样本"统计，这里显式声明列名，换数据集只改这一行。
SAMPLE_COL   <- "sample"     # 样本标识（如 GSE203115_1 / HRR1430835）
GROUP_COL    <- "group"      # 分组（ypN0 / ypN+）

# 样本在图中的排列顺序。
# 原因：ClusterDistrBar/ggplot 默认会把样本按字母排序，会打乱整合时的原始样本顺序。
#   NULL   = 沿用样本在 Seurat 对象中首次出现的顺序（即 04 步合并顺序，推荐，最"原有"）
#   "sort" = 按字母排序
#   也可直接写字符向量，如 c("GSE203115_1", "HRR1430835", ...) 指定任意顺序
SAMPLE_ORDER <- NULL

# 09/10 堆叠图的 x 轴样本名是否追加所属的分组（ypN0 / ypN+）。
# 原因：只写样本编号（如 HRR1430835）看不出它属于哪一组，读图时必须另查样本表；
#       TRUE 时标签变成两行："样本名" 换行 "(ypN0)" / "(ypN+)"，一眼可辨分组。
# 注意：这只改绘图层级的因子标签（写进 meta.data 的副本），不改变任何统计结果；
#       也只作用于 09/10 两张按样本堆叠的图，13/14（按组合并）不受影响。
#       设为 FALSE 可恢复成只显示样本名。
SAMPLE_LABEL_GROUP <- TRUE

# -------------------------------------------------------------------
# 0.3 QC / 验证用 marker 基因
# -------------------------------------------------------------------

# QC 用：画在 UMAP 上的核心 marker（覆盖本版全部 10 类，用来肉眼核对注释是否站得住）。
# 原因：不同于上一版的"打分热图自动推断"，这里只做展示、不参与任何自动判定；
#       每类挑 1 个最经典的基因，共 10 个，网格为 4 列 × 3 行。
QC_GENES <- c("CD3D", "NKG7", "MS4A1",              # T / NK / B
              "MZB1", "LYZ", "TPSAB1", "MKI67",     # 浆 / 髓系 / 肥大 / 增殖
              "PECAM1", "COL1A1", "EPCAM")          # 内皮 / 成纤维 / 上皮

# -------------------------------------------------------------------
# 0.4 SeuratExtend 绘图参数
# -------------------------------------------------------------------

# FeaturePlot3.grid 的点大小（pt.size）。
# 原因：CD45+ 数据约 8 万细胞，点太小看不清表达分布、太大互相遮挡；
# 0.5 在密度与可读性之间取得平衡。
FP_PT_SIZE  <- 0.5

# FeaturePlot3.grid 的配色方案（color）。
# 原因："ryb"（红-黄-蓝）三通道渐变对"低-中-高"表达层次区分直观，是 SeuratExtend 推荐默认。
FP_COLOR    <- "ryb"

# FeaturePlot3.grid 每行展示的基因数（ncol）。
# 原因：4 列网格在常用 PDF 宽度下每幅图大小合适、不易过密。
FP_NCOL     <- 4

# DotPlot2 点图配色方案（color_scheme）。
# 原因：白→红 Sequential 更符合点图习惯（与 pythonver 的 Reds 风格一致）。
DOT_COLOR_SCHEME <- "Reds"

# DotPlot2 是否绘制点边框（border）。原因：FALSE 去除边框更整洁，点多时边框显杂乱。
DOT_BORDER       <- FALSE

# DotPlot2 是否显示网格线（show_grid）。原因：FALSE 减少视觉干扰，便于聚焦点与色。
DOT_SHOW_GRID    <- FALSE

# DotPlot2 是否翻转坐标（flip）。原因：TRUE 使基因名位于横轴、水平可读。
DOT_FLIP         <- TRUE

# 点图尺寸参数（每个基因/每行占位英寸）。原因：随基因数与类型数自适应画布，避免挤压。
DOT_W_PER_GENE <- 0.45
DOT_H_PER_ROW  <- 0.5
DOT_BASE_H     <- 2

# 细胞类型组成图是否堆叠（stack）。原因：TRUE 时按组（ypN0/ypN+）画比例堆叠条，
# 直观对比两组免疫组成差异；FALSE 则画并排绝对计数。
COMP_STACK  <- TRUE

# PDF 出图设备。原因：类型名含希腊字母 gamma-delta，Windows 默认 pdf() 设备的标准字体
# 不含希腊字形，容易出方块/问号；cairo_pdf 走系统字体、对 UTF-8 更友好。
# 若你的环境 cairo 不可用，把它改成 "pdf" 即可（gdT 可能显示为乱码）。
PDF_DEVICE  <- "cairo_pdf"

# -------------------------------------------------------------------
# 0.7 新增分析（堆叠图 / 箱式图 / RO/E 热图）统一尺寸与分辨率
# -------------------------------------------------------------------
# 三张新图统一以 PDF + PNG 双格式落盘，且尺寸/分辨率保持一致，便于并排比较与论文插入。
# 原因：用户要求"图片统一保存为 PDF/PNG，尺寸与分辨率一致"，故集中定义常量，避免各处写死。
IMG_W    <- 10     # PDF/PNG 统一宽度（英寸）
IMG_H    <- 7      # PDF/PNG 统一高度（英寸）
IMG_DPI  <- 300    # PNG 统一分辨率（dpi）；PDF 为矢量图，分辨率概念不适用

# x 轴刻度文本倾角（度）与对齐方式。
# 原因：样本名（如 HRR1430835）较长，水平排列会互相重叠，故统一倾斜 45 度并右对齐（hjust=1）。
#      设为 0 可恢复水平。所有新增图（09~15）统一套用，保证风格一致。
X_TEXT_ANGLE <- 45
X_TEXT_HJUST <- 1

# 分组（ypN0 / ypN+）顺序。原因：组级别堆叠图与 RO/E 分析需把所有样本合并到这两组，
#   固定顺序可保证两图与统计表中 ypN0 始终在左、ypN+ 始终在右，跨图对比不串位。
GROUP_LEVELS <- c("ypN0", "ypN+")

# -------------------------------------------------------------------
# 0.5 SeuratExtend 安装参数
# -------------------------------------------------------------------

# SeuratExtend 的 GitHub 仓库。原因：该包仅发布于 GitHub（huayc09/SeuratExtend），
# CRAN 无此包，必须指定仓库地址安装。
SE_REPO     <- "huayc09/SeuratExtend"

# 安装时是否升级依赖（upgrade）。原因："never" 避免安装过程因升级大量依赖而失败/变慢。
SE_UPGRADE  <- "never"

# -------------------------------------------------------------------
# 0.6 手动 marker 字典（用于 08_marker_dotplot_manual.pdf 验证图）
# -------------------------------------------------------------------

# 载入手动整理的经典 marker 字典（已按本版 10 类重写，见 manual_markers.R）。
# 原因：点图的基因列表来自人工整理，不使用 FindAllMarkers 自动差异基因。
#
# 查找顺序：脚本目录 -> output/ 目录 -> 脚本目录上一级。
# 原因：2026-09-01 实测根目录的 manual_markers.R 被误删，脚本静默跳过点图
# （只打印一句警告，很容易漏看），故加多路径兜底 + 下方键名校验。
marker_candidates <- c(
  file.path(script_dir, "manual_markers.R"),
  file.path(OUT_DIR, "manual_markers.R"),
  file.path(dirname(script_dir), "manual_markers.R")
)
MANUAL_MARKERS_FILE <- marker_candidates[file.exists(marker_candidates)][1]

if (!is.na(MANUAL_MARKERS_FILE)) {
  # 注意：不要加 encoding="UTF-8"！GBK 环境下 source(encoding="UTF-8") 会把 UTF-8 文件
  # 强行转成本机编码，遇到无法转换的字符（如希腊字母）直接报错中断；本文件键名已全部为
  # ASCII（gdT），中文注释按本机编码读取虽可能误读但不影响运行，故直接 source。
  source(MANUAL_MARKERS_FILE)
  cat(sprintf("已载入 marker 字典: %s（%d 类）\n",
              MANUAL_MARKERS_FILE, length(MANUAL_MARKERS)))
} else {
  warning("未找到 manual_markers.R（已查找: ", paste(marker_candidates, collapse = " / "),
          "），将跳过 08_marker_dotplot_manual.pdf。")
  MANUAL_MARKERS <- list()
}

# 键名校验：marker 字典的键必须与 MANUAL_CELLTYPE_MAP 的类型名一一对应。
# 原因：两者靠字符串匹配关联，键名对不上（例如字典还是旧的 CD8T/CD4T 版本）时
#       点图会整片空白却没有任何报错，属于最容易漏掉的静默失败。
if (length(MANUAL_MARKERS) > 0) {
  mk_missing <- setdiff(names(MANUAL_CELLTYPE_MAP), names(MANUAL_MARKERS))
  mk_extra   <- setdiff(names(MANUAL_MARKERS), names(MANUAL_CELLTYPE_MAP))
  if (length(mk_missing) > 0) {
    warning("manual_markers.R 缺少这些类型的 marker，点图中不会有它们: ",
            paste(mk_missing, collapse = ", "))
  }
  if (length(mk_extra) > 0) {
    cat(sprintf("提示: manual_markers.R 中有 %d 个类型不在本次注释中，已忽略: %s\n",
                length(mk_extra), paste(mk_extra, collapse = ", ")))
  }
}


## ---- 1. 依赖检查与加载 ----

# 通用依赖安装函数：ensure_pkg("包名")，未安装时自动从 CRAN 安装（安静模式）。
# 原因：让脚本在"干净环境"也能一键跑通，无需用户手动装包。
ensure_pkg <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) return(TRUE)             # 已装则跳过
  cat(sprintf("  正在安装依赖包: %s ...\n", pkg))                    # 提示开始安装
  install.packages(pkg, repos = "https://cloud.r-project.org")        # 从官方镜像安装
  requireNamespace(pkg, quietly = TRUE)                               # 返回是否安装成功
}

# SeuratExtend 专用安装函数：优先尝试从 GitHub 安装（源码包，需 Rtools 编译）。
# 原因：SeuratExtend 不在 CRAN，且安装可能失败（无网络/缺 Rtools），故单独处理。
ensure_seuratextend <- function(repo = SE_REPO, upgrade = SE_UPGRADE) {
  if (requireNamespace("SeuratExtend", quietly = TRUE)) return(TRUE)  # 已装则跳过
  cat(sprintf("  正在安装 SeuratExtend（GitHub 源码包 %s，需数分钟与 Rtools）...\n", repo))
  if (!requireNamespace("remotes", quietly = TRUE)) {                 # 确保 remotes 可用
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  ok <- tryCatch({
    remotes::install_github(repo, upgrade = upgrade, quiet = TRUE)    # 从 GitHub 安装
    requireNamespace("SeuratExtend", quietly = TRUE)                  # 校验是否装上
  }, error = function(e) {
    cat("  SeuratExtend 安装失败:", conditionMessage(e), "\n")        # 失败仅提示
    FALSE
  })
  ok
}

# 依次确保基础依赖可用（dplyr 数据整理 / ggplot2 绘图）。
# 说明：上一版的 pheatmap 仅服务于打分热图，本版已删除打分，故不再依赖。
for (p in c("Seurat", "dplyr", "ggplot2", "tidyr", "remotes")) ensure_pkg(p)

# 安装并加载 SeuratExtend（本脚本主要绘图方式）。失败则明确报错并退出。
if (!ensure_seuratextend()) {
  stop("SeuratExtend 安装/加载失败，无法继续（请确认网络可访问 GitHub 且已安装 Rtools）。")
}

# 加载核心包。SeuratExtend 提供 DimPlot2 / DotPlot2 / FeaturePlot3.grid / ClusterDistrBar。
library(Seurat)        # 单细胞分析主包（readRDS / 元数据操作）
library(SeuratExtend)  # 主要绘图方式
library(dplyr)         # 数据整理（group_by / summarise / mutate 等）
library(ggplot2)       # 组成图 ggplot 兜底
library(tidyr)        # RO/E 宽表转换（pivot_wider / pivot_longer）

# 创建输出目录（已存在则不报错）
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 打印运行环境信息
cat(sprintf("工作目录: %s\n", getwd()))
cat(sprintf("输出目录: %s\n", normalizePath(OUT_DIR)))

# 注：原"wmic 查可用物理内存"的提示已移除——wmic 在本环境被安全策略列入黑名单，
# 调用会触发进程异常。本步读 2.1GB 对象请自行确保可用内存 >=8GB。

# Seurat v5 全局对象大小限制（Future 框架）。
# 原因：默认上限较小，读取大对象时会报 "serialized object too large"，提高到 8GB。
options(future.globals.maxSize = 8 * 1024^3)


## ---- 1b. 辅助函数（每个任务一个函数，便于单独定位与修改）----

# 1b-1 统一的 PDF 落盘函数（支持 cairo_pdf / pdf 两种设备）
# 原因：类型名含 gamma-delta，需 UTF-8 友好的设备；同时统一收口宽高与 limitsize 设置。
save_pdf <- function(plot_obj, filename, width, height) {
  dev_fun <- tryCatch(get(PDF_DEVICE), error = function(e) grDevices::pdf)
  ggsave(filename = filename, plot = plot_obj, device = dev_fun,
         width = width, height = height, limitsize = FALSE)
  cat("  已保存图表:", basename(filename), "\n")
}

# 1b-1b 统一的 PNG 落盘函数（与 save_pdf 同尺寸/同分辨率，用于论文或网页插入）
# 原因：类型名含 gamma-delta，PNG 默认设备字体也可能乱码，故用 type="cairo" 走系统字体（与 cairo_pdf 同源）。
save_png <- function(plot_obj, filename, width, height, dpi = IMG_DPI) {
  ggsave(filename = filename, plot = plot_obj,
         device = "png", type = "cairo",
         width = width, height = height, dpi = dpi, limitsize = FALSE)
  cat("  已保存图表:", basename(filename), "\n")
}

# 1b-1c 统一的 x 轴文本倾斜主题（默认 45 度、右对齐）
# 原因：样本名（HRR1430835 等）较长，水平排列必然重叠，需倾斜；但每张图都手抄
#      element_text(angle=..., hjust=...) 易遗漏且难统一调整，故抽成函数集中收口。
# 用法：在任意 ggplot 对象上 `+ theme_x45()`；传 angle=0 可单独取消倾斜。
theme_x45 <- function(angle = X_TEXT_ANGLE, hjust = X_TEXT_HJUST) {
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = angle, hjust = hjust))
}

# 1b-1d 写 UTF-8（带 BOM）CSV 的辅助函数
# 原因：① 直接 write.csv(..., fileEncoding="UTF-8-BOM") 在本机 R 4.4.3 会报
#   "unsupported conversion from 'UTF-8-BOM' to ''"；② write.csv(..., append=TRUE)
#   在本机也会被忽略，无法先写 BOM 再追加内容。故改为：先写到临时文件（UTF-8），
#   再读取全部字节并在最前面插入 BOM。这样可保证带 BOM，Excel 与 pandas 都能正确识别。
# 注意：在 GBK locale 的 R 中，希腊字母（γδT）字符串字面量常被解析为 <U+03B3><U+03B4>T，
#   本函数无法把 <U+XXXX> 转义还原为希腊字母；如需图中显示希腊字母，需把类型名改为"gdT"
#   或在 UTF-8 locale 下运行 R。
write_utf8_csv <- function(df, path, row.names = FALSE) {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))                               # 退出时自动清理临时文件
  write.csv(df, tmp, row.names = row.names,
            fileEncoding = "UTF-8")                 # 临时文件为 UTF-8（无 BOM）
  body <- readBin(tmp, "raw", n = file.info(tmp)$size)
  con <- file(path, open = "wb")
  writeBin(as.raw(c(0xef, 0xbb, 0xbf)), con)         # UTF-8 BOM
  writeBin(body, con)                                 # 追加 CSV 内容
  close(con)
}

# 1b-2 注释 UMAP：仅画整体（按细胞类型着色 + 标签）
#   说明：原 03_annotation_UMAP_by_group.pdf（按 ypN0/ypN+ 分面的 UMAP）已删除，
#   原因是它与 13/14 组级别堆叠图语义重复，且单维 UMAP 不如堆叠条形图直观。
plot_annotation_umap <- function(obj, OUT_DIR) {
  p_anno <- DimPlot2(obj,
    features = "celltype",   # 用 celltype 作为着色变量
    label    = TRUE,         # 显示类型标签
    box      = TRUE,         # 标签加框更清晰
    repel    = TRUE,         # 标签自动避让避免重叠
    theme    = NoLegend())   # 右侧图例占位去除（标签已直接标在图上）
  save_pdf(p_anno, file.path(OUT_DIR, "02_annotation_UMAP.pdf"), width = 9, height = 7)
}

# 1b-3 细胞类型组成统计表：总体 + 按分组
write_composition_stats <- function(obj, OUT_DIR) {
  ct_summary <- obj@meta.data %>%
    dplyr::group_by(celltype) %>%
    dplyr::summarise(n_cells = dplyr::n(),
                     pct     = 100 * dplyr::n() / ncol(obj),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_cells))
  # 说明：用 write_utf8_csv（带 BOM 的 UTF-8），避免 GBK 环境下 gamma-delta 名称被转义成 <U+03B3><U+03B4>T
  write_utf8_csv(ct_summary, file.path(OUT_DIR, "05_celltype_summary.csv"),
                row.names = FALSE)
  cat("  已保存: 05_celltype_summary.csv\n")

  ct_group <- obj@meta.data %>%
    dplyr::group_by(celltype, group) %>%
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop")
  write_utf8_csv(ct_group, file.path(OUT_DIR, "06_celltype_by_group.csv"),
                row.names = FALSE)
  cat("  已保存: 06_celltype_by_group.csv\n")
}

# 1b-4 细胞类型组成图（按"样本 × 细胞类型"堆叠 / 横向版）：ClusterDistrBar 为主，ggplot 兜底
#   历史原因：早期版本此函数输出"按分组堆叠"的 04_celltype_composition.pdf，
#   后续在 4i 加入了更详尽的 13_ / 14_（按分组、按比例 / 绝对计数双份），
#   两张图完全重复，2026-09 删除 04_* 产物与本函数。如有特殊需求请改用 13_stacked_composition_percent_by_group。
#   保留此函数位置以保证行号稳定，调用方已不再触发。
plot_composition <- function(obj, OUT_DIR) {
  # 故意不产生任何产物；保留函数定义仅为兼容可能的外部调用。
  invisible(NULL)
}

# 1b-5 QC：核心 marker 的 UMAP 表达网格（只展示，不参与任何自动判定）
plot_marker_umap <- function(obj, feat_genes, OUT_DIR) {
  p_feat <- tryCatch(
    FeaturePlot3.grid(obj,
      features = feat_genes,   # 要展示的基因
      pt.size  = FP_PT_SIZE,   # 点大小（见配置 FP_PT_SIZE 注释）
      color    = FP_COLOR,     # 配色（见配置 FP_COLOR 注释）
      ncol     = FP_NCOL),     # 每行基因数（见配置 FP_NCOL 注释）
    error = function(e) {
      cat("  FeaturePlot3.grid 不支持 color/ncol 参数，回退最小参数集:", conditionMessage(e), "\n")
      FeaturePlot3.grid(obj, features = feat_genes, pt.size = FP_PT_SIZE)
    })
  save_pdf(p_feat, file.path(OUT_DIR, "07_marker_UMAP.pdf"),
           width  = 5.5 * FP_NCOL,                               # 画布宽随列数自适应
           height = 5.5 * ceiling(length(feat_genes) / FP_NCOL)) # 高随行数自适应
}

# 1b-6 验证：手动经典 marker 点图（基因来自 manual_markers.R，非自动差异基因）
plot_manual_dotplot <- function(obj, OUT_DIR) {
  if (!exists("MANUAL_MARKERS") || length(MANUAL_MARKERS) == 0) {
    cat("  注意: manual_markers.R 为空，跳过 marker 点图\n")
    return(invisible(NULL))
  }
  # 只保留对象里真实存在、且属于本次注释结果的类型；缺失基因一并过滤
  present_types <- intersect(TYPE_ORDER, as.character(unique(obj$celltype)))
  mk <- MANUAL_MARKERS[intersect(present_types, names(MANUAL_MARKERS))]
  mk <- lapply(mk, function(g) g[g %in% rownames(obj)])
  mk <- mk[lengths(mk) > 0]
  if (length(mk) == 0) {
    cat("  注意: 没有可用 marker，跳过点图\n")
    return(invisible(NULL))
  }
  p <- DotPlot2(obj,
    features     = mk,
    group.by     = "celltype",
    color_scheme = DOT_COLOR_SCHEME,
    border       = DOT_BORDER,
    show_grid    = DOT_SHOW_GRID,
    flip         = DOT_FLIP)
  n_cols <- sum(lengths(mk))
  n_rows <- length(mk)
  save_pdf(p, file.path(OUT_DIR, "08_marker_dotplot_manual.pdf"),
           width  = max(10, n_cols * DOT_W_PER_GENE),
           height = max(6,  n_rows * DOT_H_PER_ROW + DOT_BASE_H))
}

# 1b-7 纵向堆叠条形图（样本置于 x 轴，各细胞类型沿 y 轴堆叠）
#   输出两张：① 比例堆叠（每样本各细胞类型占比之和=1）；② 绝对细胞数堆叠
#   配色用 ClusterDistrBar(cols="light")，与 13/14 组级别堆叠图同源，跨图可比；
#   x 轴样本名按 SAMPLE_LABEL_GROUP 追加一行所属分组（ypN0 / ypN+），避免读图时另查样本表；
#   样本顺序由 SAMPLE_ORDER 控制（NULL=对象首次出现顺序，最"原有"）。
plot_stacked_bars <- function(obj, OUT_DIR) {
  # 解析样本顺序（与 13/14 组级别图一致逻辑）
  if (is.null(SAMPLE_ORDER)) {
    sample_levels <- unique(as.character(obj@meta.data[[SAMPLE_COL]]))   # 首次出现顺序
  } else if (identical(SAMPLE_ORDER, "sort")) {
    sample_levels <- sort(unique(as.character(obj@meta.data[[SAMPLE_COL]])))  # 字母序
  } else {
    sample_levels <- SAMPLE_ORDER                                        # 显式指定
  }
  # 把 sample 设为因子并固定 levels，ClusterDistrBar 会沿 levels 排 x 轴顺序
  if (isTRUE(SAMPLE_LABEL_GROUP)) {
    # 每个样本只需取它第一个细胞对应的分组即可（同一样本的 group 必然唯一）
    grp_per_sample <- vapply(sample_levels, function(s) {
      as.character(obj@meta.data[[GROUP_COL]][match(s, obj@meta.data[[SAMPLE_COL]])])
    }, character(1), USE.NAMES = FALSE)
    # levels 也要同步加后缀，否则因子值与 levels 对不上会被整体置为 NA
    sample_levels_lbl <- sprintf("%s\n(%s)", sample_levels, grp_per_sample)
    obj[[SAMPLE_COL]] <- factor(sprintf("%s\n(%s)",
                                        as.character(obj@meta.data[[SAMPLE_COL]]),
                                        as.character(obj@meta.data[[GROUP_COL]])),
                                levels = sample_levels_lbl)
  } else {
    obj[[SAMPLE_COL]] <- factor(as.character(obj@meta.data[[SAMPLE_COL]]), levels = sample_levels)
  }

  # ① 比例堆叠（flip=FALSE 纵向；percent=TRUE 表示每样本占比）
  #   注意：obj$sample / obj$celltype 经 "$" 取为向量（避开 obj[["x"]] 返回 data.frame 的坑）
  p_prop <- ClusterDistrBar(origin  = obj$sample,
                            cluster = obj$celltype,
                            flip    = FALSE,   # 纵向：样本在 x 轴，细胞类型沿 y 轴堆叠
                            stack   = TRUE,
                            percent = TRUE,    # 比例（每样本合计=100%）
                            cols    = "light") +  # 与 13/14 组级别图同源配色
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = X_TEXT_ANGLE, hjust = X_TEXT_HJUST,
                                                       lineheight = 0.9))
  save_pdf(p_prop, file.path(OUT_DIR, "09_stacked_composition_percent.pdf"),
           width = IMG_W, height = IMG_H)
  save_png(p_prop, file.path(OUT_DIR, "09_stacked_composition_percent.png"),
           width = IMG_W, height = IMG_H)

  # ② 绝对细胞数堆叠（percent=FALSE）
  p_count <- ClusterDistrBar(origin  = obj$sample,
                              cluster = obj$celltype,
                              flip    = FALSE,
                              stack   = TRUE,
                              percent = FALSE,   # 绝对计数
                              cols    = "light") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = X_TEXT_ANGLE, hjust = X_TEXT_HJUST,
                                                       lineheight = 0.9))
  save_pdf(p_count, file.path(OUT_DIR, "10_stacked_composition_count.pdf"),
           width = IMG_W, height = IMG_H)
  save_png(p_count, file.path(OUT_DIR, "10_stacked_composition_count.png"),
           width = IMG_W, height = IMG_H)
}

# 1b-8 分组箱式图：每样本各细胞类型比例，按 ypN0/ypN+ 分组，每个细胞类型一个分面
#   添加 Wilcoxon 秩和检验（两组比较）的显著性标注（ns / * p<0.05 / ** p<0.01 / *** p<0.001）
plot_proportion_boxplot <- function(obj, OUT_DIR) {
  # 每样本 × 细胞类型 计数，并换算成"每样本内该细胞类型的占比"
  cnt <- obj@meta.data %>%
    dplyr::count(.data[[SAMPLE_COL]], celltype) %>%
    dplyr::rename(sample = .data[[SAMPLE_COL]]) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(prop = n / sum(n)) %>%          # 每样本内占比（分组基础）
    dplyr::ungroup()
  # 关联分组（ypN0/ypN+）；同一 sample 的 group 必然一致，故去重
  grp <- obj@meta.data %>%
    dplyr::distinct(.data[[SAMPLE_COL]], .data[[GROUP_COL]]) %>%
    dplyr::rename(sample = .data[[SAMPLE_COL]], group = .data[[GROUP_COL]])
  prop_df <- cnt %>% dplyr::left_join(grp, by = "sample")

  # 逐细胞类型做 Wilcoxon 检验（prop ~ group），将 p 值映射为星号，并保留具体数值
  wil <- prop_df %>%
    dplyr::group_by(celltype) %>%
    dplyr::summarise(
      p = wilcox.test(prop ~ group, data = cur_data(), exact = FALSE)$p.value,  # 非精确（含并列）避免警告
      .groups = "drop")
  wil$stars <- cut(wil$p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                   labels = c("***", "**", "*", "ns"), right = FALSE)
  wil$ptxt  <- ifelse(wil$p < 0.001, sprintf("%.1e", wil$p),
                      sprintf("%.3f", wil$p))                       # 具体 p 值文本（论文引用用）
  cat("  Wilcoxon 检验结果（细胞类型 × ypN 组）:\n")
  print(wil)

  # 箱式图：x=分组，y=每样本占比，按细胞类型分面；右上角标注显著性
  p_box <- ggplot2::ggplot(prop_df, aes(x = group, y = prop, fill = group)) +
    ggplot2::geom_boxplot(width = 0.6, outlier.size = 1) +
    ggplot2::scale_fill_brewer(palette = "Set2") +                 # 两组两色，与细胞类型配色区分
    ggplot2::facet_wrap(~ celltype, scales = "free_y") +           # 每细胞类型一个分面，y 轴独立缩放
    ggplot2::labs(x = "Group", y = "Cell-type proportion (per sample)",
                  title = "Per-sample cell-type proportion by ypN group") +
    ggplot2::theme_minimal() +
    ggplot2::theme(strip.text = element_text(size = 9)) +
    theme_x45() +                                                  # x 轴（ypN0/ypN+）文本倾斜，与其余新图风格统一
    # 每个分面右上角叠加 Wilcoxon 结论（星号 + p 值）；celltype 列使各分面只显示自身结果
    ggplot2::geom_text(
      data = wil,
      aes(x = Inf, y = Inf, label = paste0(stars, "\n", "p=", ptxt)),
      inherit.aes = FALSE, hjust = 1.1, vjust = 1.5, size = 3, color = "black")
  save_pdf(p_box, file.path(OUT_DIR, "11_proportion_boxplot_by_group.pdf"),
           width = IMG_W, height = IMG_H)
  save_png(p_box, file.path(OUT_DIR, "11_proportion_boxplot_by_group.png"),
           width = IMG_W, height = IMG_H)
}

# 1b-9 RO/E（观察/期望）富集分析
#   列联表：行=样本，列=细胞类型
#     - 观察数 O_{ij} = 样本 i 中细胞类型 j 的细胞数
#     - 期望数 E_{ij} = (样本 i 的行合计) × (细胞类型 j 的列合计) / 总计
#       这正是卡方独立性检验的"期望频数"标准公式（行合计 × 列合计 / 总计），
#       也是单细胞组成分析中"相对富集"的标准做法（等价于 该样本内该类型占比 / 总体该类型占比）。
#     - RO/E_{ij} = O_{ij} / E_{ij}：>1 表示过代表（富集），<1 表示欠代表。
#   输出：① 观察计数矩阵 CSV；② RO/E 比值矩阵 CSV（保留原始数值以便复查/补分析）；
#        ③ 带数值标注的热图（PDF + PNG）。
#   关于"是否为标准做法"的确认：用户所给定义即标准卡方期望频数，无需改动，已按此实现。
roe_analysis <- function(obj, OUT_DIR) {
  # ① 观察计数矩阵（行=样本，列=细胞类型）
  cnt_wide <- obj@meta.data %>%
    dplyr::count(.data[[SAMPLE_COL]], celltype) %>%
    dplyr::rename(sample = .data[[SAMPLE_COL]]) %>%
    tidyr::pivot_wider(id_cols = sample, names_from = celltype,
                       values_from = n, values_fill = 0)
  # 按样本顺序（SAMPLE_ORDER）与细胞类型顺序（绘图 levels）排好矩阵行列
  if (is.null(SAMPLE_ORDER)) {
    sample_levels <- unique(as.character(obj@meta.data[[SAMPLE_COL]]))
  } else if (identical(SAMPLE_ORDER, "sort")) {
    sample_levels <- sort(unique(as.character(obj@meta.data[[SAMPLE_COL]])))
  } else {
    sample_levels <- SAMPLE_ORDER
  }
  ct_levels <- levels(obj$celltype)
  cnt_wide <- cnt_wide %>% dplyr::arrange(factor(sample, levels = sample_levels))
  obs_mat <- as.matrix(cnt_wide[, ct_levels, drop = FALSE])
  rownames(obs_mat) <- cnt_wide$sample

  # ② 期望计数与 RO/E（标准卡方期望：行合计 × 列合计 / 总计）
  row_tot <- rowSums(obs_mat)                 # 各样本总细胞数（行合计）
  col_tot <- colSums(obs_mat)                 # 各细胞类型总细胞数（列合计）
  grand  <- sum(obs_mat)                       # 总计
  exp_mat <- outer(row_tot, col_tot) / grand  # 期望频数矩阵 E_{ij}
  roe_mat <- obs_mat / exp_mat                # RO/E 富集矩阵

  # 落盘计数矩阵与比值矩阵（保留原始数值，便于复查/补分析）
  # 说明：用 write_utf8_csv（带 BOM 的 UTF-8），避免 GBK 环境下 gamma-delta 名称被转义失真。
  write_utf8_csv(as.data.frame(obs_mat), file.path(OUT_DIR, "12_ROE_observed_counts.csv"),
                 row.names = TRUE)
  write_utf8_csv(as.data.frame(roe_mat), file.path(OUT_DIR, "12_ROE_ratio_matrix.csv"),
                 row.names = TRUE)
  cat("  已保存: 12_ROE_observed_counts.csv / 12_ROE_ratio_matrix.csv\n")

  # ③ 热图：SeuratExtend::Heatmap 绘制 RO/E 矩阵（蓝-红双色，值越大越红=越富集），
  #    并叠加 geom_text 标注每个格子的比值（保留 2 位小数）。
  #    说明：Heatmap 返回 ggplot，内部数据列为 id(行)/variable(列)/value(比值)，可直接标注。
  p_roe <- SeuratExtend::Heatmap(roe_mat, color_scheme = "BuRd") +
    ggplot2::geom_text(aes(label = sprintf("%.2f", value)),
                       color = "black", size = 3) +
    theme_x45()   # x 轴细胞类型名倾斜（统一由 X_TEXT_ANGLE 控制）
  save_pdf(p_roe, file.path(OUT_DIR, "12_ROE_heatmap.pdf"),
           width = IMG_W, height = IMG_H)
  save_png(p_roe, file.path(OUT_DIR, "12_ROE_heatmap.png"),
           width = IMG_W, height = IMG_H)
}

# 1b-10 组级别（样本合并）纵向堆叠条形图：x 轴 = ypN0 / ypN+（所有样本按分组归并）
#   与 1b-7 完全同源（同一 ClusterDistrBar 调用 / 同一配色 / 同一尺寸），
#   差别仅 origin 用 group 列而非 sample 列，用于看"两组的总体构成差异"。
#   输出：13_ 比例堆叠（每组合计=100%）、14_ 绝对计数堆叠。
plot_stacked_bars_group <- function(obj, OUT_DIR) {
  # 分组顺序以 GROUP_LEVELS 为准，并把实际出现但未列出的分组追加到末尾，避免丢失/NA。
  present_grps <- unique(as.character(obj@meta.data[[GROUP_COL]]))
  grp_levels   <- c(intersect(GROUP_LEVELS, present_grps),
                    setdiff(present_grps, GROUP_LEVELS))
  grp <- factor(as.character(obj@meta.data[[GROUP_COL]]), levels = grp_levels)

  # ① 比例堆叠（每组合计=100%）
  p_prop <- ClusterDistrBar(origin  = grp,
                            cluster = obj$celltype,
                            flip    = FALSE,
                            stack   = TRUE,
                            percent = TRUE,
                            cols    = "light") +
    theme_x45()   # x 轴分组名倾斜（统一风格）
  save_pdf(p_prop, file.path(OUT_DIR, "13_stacked_composition_percent_by_group.pdf"),
           width = IMG_W, height = IMG_H)
  save_png(p_prop, file.path(OUT_DIR, "13_stacked_composition_percent_by_group.png"),
           width = IMG_W, height = IMG_H)

  # ② 绝对计数堆叠
  p_count <- ClusterDistrBar(origin  = grp,
                             cluster = obj$celltype,
                             flip    = FALSE,
                             stack   = TRUE,
                             percent = FALSE,
                             cols    = "light") +
    theme_x45()
  save_pdf(p_count, file.path(OUT_DIR, "14_stacked_composition_count_by_group.pdf"),
           width = IMG_W, height = IMG_H)
  save_png(p_count, file.path(OUT_DIR, "14_stacked_composition_count_by_group.png"),
           width = IMG_W, height = IMG_H)
}

# 1b-11 组级别 RO/E（观察/期望）富集分析：行=ypN0/ypN+ 两组，列=细胞类型
#   与 1b-9 完全同源（同一卡方期望公式、同一热图、同一尺寸），差别仅把样本合并到分组。
#   额外输出：卡方独立性检验（两组总体构成是否有显著差异）+ 组级别汇总表（计数/占比/ROE）。
#   输出：15_ROE_observed_counts_by_group.csv、15_ROE_ratio_matrix_by_group.csv、
#        15_ROE_group_summary.csv、15_ROE_heatmap_by_group.pdf/.png。
roe_analysis_group <- function(obj, OUT_DIR) {
  # ① 观察计数矩阵（行=分组，列=细胞类型）
  cnt_wide <- obj@meta.data %>%
    dplyr::count(.data[[GROUP_COL]], celltype) %>%
    dplyr::rename(group = .data[[GROUP_COL]]) %>%
    tidyr::pivot_wider(id_cols = group, names_from = celltype,
                       values_from = n, values_fill = 0)
  present_grps <- unique(as.character(obj@meta.data[[GROUP_COL]]))
  grp_levels   <- c(intersect(GROUP_LEVELS, present_grps),
                    setdiff(present_grps, GROUP_LEVELS))
  ct_levels    <- levels(obj$celltype)
  cnt_wide <- cnt_wide %>% dplyr::arrange(factor(group, levels = grp_levels))
  obs_mat <- as.matrix(cnt_wide[, ct_levels, drop = FALSE])
  rownames(obs_mat) <- cnt_wide$group

  # ② 期望计数与 RO/E（标准卡方期望：行合计 × 列合计 / 总计）
  row_tot <- rowSums(obs_mat)
  col_tot <- colSums(obs_mat)
  grand  <- sum(obs_mat)
  exp_mat <- outer(row_tot, col_tot) / grand
  roe_mat <- obs_mat / exp_mat

  # ③ 卡方独立性检验：两组总体细胞构成是否不同（2 行仍有效）
  cs <- suppressWarnings(chisq.test(obs_mat))
  cat(sprintf("  组级别卡方独立性检验: chi2=%.2f, df=%d, p=%.3e\n",
              cs$statistic, cs$parameter, cs$p.value))

  # ④ 落盘计数矩阵与比值矩阵（write_utf8_csv：带 BOM 的 UTF-8，避免 gamma-delta 名称被转义）
  write_utf8_csv(as.data.frame(obs_mat), file.path(OUT_DIR, "15_ROE_observed_counts_by_group.csv"),
                 row.names = TRUE)
  write_utf8_csv(as.data.frame(roe_mat), file.path(OUT_DIR, "15_ROE_ratio_matrix_by_group.csv"),
                 row.names = TRUE)

  # ④b 组级别汇总表：计数 + 组内占比 + RO/E（写论文/汇报最常用的一张表）
  prop_mat <- sweep(obs_mat, 1, row_tot, "/")   # 每组内各类型占比
  sum_tbl  <- data.frame(celltype = ct_levels, stringsAsFactors = FALSE)
  for (g in rownames(obs_mat)) {
    sum_tbl[[paste0(g, "_n")]]    <- obs_mat[g, ]
    sum_tbl[[paste0(g, "_prop")]] <- round(prop_mat[g, ], 4)
    sum_tbl[[paste0(g, "_ROE")]]  <- round(roe_mat[g, ], 3)
  }
  write_utf8_csv(sum_tbl, file.path(OUT_DIR, "15_ROE_group_summary.csv"),
                 row.names = FALSE)
  cat("  已保存: 15_ROE_observed_counts_by_group.csv / 15_ROE_ratio_matrix_by_group.csv / 15_ROE_group_summary.csv\n")

  # ⑤ 热图（2 行 × 细胞类型 列；蓝-红双色，叠加数值标注）
  p_roe <- SeuratExtend::Heatmap(roe_mat, color_scheme = "BuRd") +
    ggplot2::geom_text(aes(label = sprintf("%.2f", value)),
                       color = "black", size = 3) +
    theme_x45()
  save_pdf(p_roe, file.path(OUT_DIR, "15_ROE_heatmap_by_group.pdf"),
           width = IMG_W, height = IMG_H)
  save_png(p_roe, file.path(OUT_DIR, "15_ROE_heatmap_by_group.png"),
           width = IMG_W, height = IMG_H)
}

## ---- 2. 读入对象与基础检查 ----

# 2.1 检查输入文件是否存在，缺失则终止并给出明确路径。
if (!file.exists(INPUT_RDS)) stop("找不到输入 rds: ", INPUT_RDS)

# 2.2 读取第四步的整合聚类对象（约 8 万细胞，读取需 1~3 分钟）。
cat("正在读取 04_CD45_integrated.rds（约需 1~3 分钟）...\n")
obj <- readRDS(INPUT_RDS)
cat(sprintf("读入完成：%d 个基因 x %d 个细胞\n", nrow(obj), ncol(obj)))

# 2.3 检查主 cluster 列是否存在（04 步产物应包含）。
if (!CLUSTER_COL %in% colnames(obj@meta.data)) {
  stop("对象中找不到主 cluster 列 '", CLUSTER_COL, "'，请先运行 04_integration_clustering.R")
}

# 【重要】取 cluster 必须用 obj@meta.data[[CLUSTER_COL]]（返回向量），
# 不能用 obj[[CLUSTER_COL]]：Seurat v5 的 [[ 对 meta.data 返回的是"单列 data.frame"，
# 会导致 unique()/sort()/setdiff() 作用在 data.frame 上——不去重、按行比较，
# 结果把 8 万个细胞的逐细胞值当成一个元素（旧版日志里"簇数: 1"就是这个坑）。
# 本脚本统一用 cl_vec 这个字符向量做 cluster 相关运算。
cl_vec <- as.character(obj@meta.data[[CLUSTER_COL]])

# 防御性检查：cluster 向量长度必须等于细胞数，且簇数应在合理范围内。
# 原因：一旦取值方式退化成 data.frame（见上方说明），长度会对不上，
#       这里提前报错，避免后面静默地把注释打乱。
if (length(cl_vec) != ncol(obj)) {
  stop(sprintf("cluster 向量长度(%d)与细胞数(%d)不一致，请检查 CLUSTER_COL 的取值方式",
               length(cl_vec), ncol(obj)))
}
cat(sprintf("当前聚类（%s）簇数: %d\n", CLUSTER_COL, length(unique(cl_vec))))

# 若对象含 group 列（ypN0/ypN+）则用之；否则构造占位 "All"，保证下游组成图不报错。
if (!"group" %in% colnames(obj@meta.data)) {
  cat("提示: 对象中无 'group' 列，组成图将退化为单一分组 'All'。\n")
  obj$group <- "All"
}


## ---- 3. 应用手动注释（本步为分析，无图）----

cat("\n==== 3. 应用手动注释 ====\n")

# 3.1 把 MANUAL_CELLTYPE_MAP（类型 -> 多个 cluster）展开成 cluster -> 类型 的对照表。
# 原因：脚本里按"类型为单位"写更易读易改，但注释是逐 cluster 落的，需要反转。
map_df <- data.frame(
  cluster  = unlist(MANUAL_CELLTYPE_MAP, use.names = FALSE),
  celltype = rep(names(MANUAL_CELLTYPE_MAP), lengths(MANUAL_CELLTYPE_MAP)),
  stringsAsFactors = FALSE
)
map_df$cluster <- as.character(map_df$cluster)

# 3.2 查重：同一个 cluster 被写进两个类型时立即报错（手写表最易犯的错）。
# 原因：静默覆盖会让注释结果无声出错，必须在运行前挡住。
dup_clusters <- unique(map_df$cluster[duplicated(map_df$cluster)])
if (length(dup_clusters) > 0) {
  stop("MANUAL_CELLTYPE_MAP 中存在重复 cluster: ", paste(dup_clusters, collapse = ", "),
       " —— 请检查 ## 0.2 节，每个 cluster 只能归入一个细胞类型")
}

# 3.3 覆盖性校验：对象里出现了、但手动表没写的 cluster（标为 Unknown 并告警）。
# 说明：cl_vec 是 2.3 节定义的"逐细胞 cluster 字符向量"（从 meta.data 取，非 [[）。
obj_clusters <- sort(unique(cl_vec))
missing_clusters <- setdiff(obj_clusters, map_df$cluster)
if (length(missing_clusters) > 0) {
  cat(sprintf("警告: 有 %d 个 cluster 未在 MANUAL_CELLTYPE_MAP 中注释，将被标为 Unknown: %s\n",
              length(missing_clusters), paste(missing_clusters, collapse = ", ")))
}

# 3.4 手动表里写了、但对象里没有的 cluster（通常是换分辨率后编号变了，提示用户校对）。
extra_clusters <- setdiff(map_df$cluster, obj_clusters)
if (length(extra_clusters) > 0) {
  cat(sprintf("警告: 手动表中有 %d 个 cluster 在当前对象中不存在，已忽略: %s\n",
              length(extra_clusters), paste(extra_clusters, collapse = ", ")))
}

# 3.5 逐 cluster 写入 celltype；未覆盖的置为 "Unknown"。
# 原因：用 match 一次性映射，避免逐行循环；结果按 TYPE_ORDER 设为 factor 固定绘图顺序。
obj$celltype <- map_df$celltype[match(cl_vec, map_df$cluster)]
obj$celltype[is.na(obj$celltype)] <- "Unknown"

# 仅保留真实出现的类型作为 levels，顺序优先按 TYPE_ORDER。
# 注意：必须用 intersect + setdiff 合并，而不能只写 intersect(TYPE_ORDER, unique(...))。
# 原因：当类型标签含非 ASCII 字符（如希腊字母）时，R 内部字符串编码标记可能不一致，
# intersect 会把它判为"不相等"而剔除出 levels；factor() 遇到不在 levels 里的值会置为 NA。
# 本脚本已统一使用 ASCII 的 "gdT" 标签，问题本身不再触发；保留 intersect + setdiff
# 写法可作为防御性代码，确保任何情况下实际出现的类型都不会因编码标记差异而丢失。
# 用 setdiff 把"实际出现但被 intersect 漏掉"的类型补回 levels，可彻底避免该问题。
present_types <- unique(obj$celltype)
ordered_types <- c(intersect(TYPE_ORDER, present_types),
                  setdiff(present_types, TYPE_ORDER))   # 其余实际出现的类型追加在末尾
obj$celltype  <- factor(obj$celltype, levels = ordered_types)

cat(sprintf("注释完成：%d 个细胞、%d 种细胞类型\n", ncol(obj), nlevels(obj$celltype)))
print(table(obj$celltype, useNA = "ifany"))
if (any(obj$celltype == "Unknown", na.rm = TRUE)) {
  cat(sprintf("注意: %d 个细胞未注释（Unknown），请在 ## 0.2 MANUAL_CELLTYPE_MAP 中补全后重跑\n",
              sum(obj$celltype == "Unknown")))
}

# 3.6 落盘实际使用的注释表（cluster | celltype | n_cells），便于复查与写论文时引用。
# 原因：把"实际生效的注释"固化下来，避免脚本改动后无法追溯当时用的哪一版。
applied <- obj@meta.data %>%
  dplyr::group_by(cluster = .data[[CLUSTER_COL]], celltype) %>%
  dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(as.numeric(as.character(cluster)))
write_utf8_csv(applied, file.path(OUT_DIR, "01_celltype_map_applied.csv"),
               row.names = FALSE)
cat("  已保存: 01_celltype_map_applied.csv\n")


## ---- 4. 注释后可视化与统计 ----

cat("\n==== 4. 注释后可视化与统计 ====\n")
plot_annotation_umap(obj, OUT_DIR)          # 4a: 注释 UMAP（整体）
write_composition_stats(obj, OUT_DIR)       # 4b: 组成统计表（总体 + 按分组）
# 4c: 按"分组"汇总的组成图 —— 已与 4i (13_/14_) 完全重复，2026-09 删除以避免冗余产物。

# 4d: QC —— 核心 marker 的 UMAP 表达网格（只用于肉眼核对注释，不参与判定）
feat_genes <- QC_GENES[QC_GENES %in% rownames(obj)]                 # 过滤数据中缺失的基因
if (length(QC_GENES) > length(feat_genes)) {
  cat(sprintf("  提示: QC 基因中有 %d 个在数据中不存在，已跳过: %s\n",
              length(QC_GENES) - length(feat_genes),
              paste(setdiff(QC_GENES, feat_genes), collapse = "/")))
}
if (length(feat_genes) > 0) {
  cat(sprintf("  使用 %d 个 QC marker 基因作图\n", length(feat_genes)))
  plot_marker_umap(obj, feat_genes, OUT_DIR)
}

# 4e: 验证 —— 手动经典 marker 点图
plot_manual_dotplot(obj, OUT_DIR)

# 4f: 纵向堆叠条形图（比例 + 绝对细胞数，样本在 x 轴并标注所属分组，配色与 13/14 同源）
plot_stacked_bars(obj, OUT_DIR)

# 4g: 分组箱式图（每样本细胞类型比例，按 ypN0/ypN+ 分组，每细胞类型一分面 + Wilcoxon 标注）
plot_proportion_boxplot(obj, OUT_DIR)

# 4h: RO/E 富集分析（观察/期望比值热图 + 计数矩阵/比值矩阵落盘）→ 12_*
roe_analysis(obj, OUT_DIR)

# 4i: 组级别（ypN0/ypN+）合并堆叠图（比例 + 绝对计数）→ 13_* / 14_*
plot_stacked_bars_group(obj, OUT_DIR)

# 4j: 组级别 RO/E 富集分析（两组比较，含卡方检验与汇总表）→ 15_*
roe_analysis_group(obj, OUT_DIR)


## ---- 5. 保存注释对象与完成 ----

# 5.1 保存注释后的对象（meta 含 celltype 列；本版不再生成 Score_* 打分列）。
saveRDS(obj, file.path(OUT_DIR, "05_annotated.rds"))
cat(sprintf("  已保存: %s\n", file.path(OUT_DIR, "05_annotated.rds")))

# 5.2 打印最终结果摘要
cat("\n==== 05 细胞注释完成 ====\n")
cat(sprintf("已注释 %d 个细胞、%d 种细胞类型\n", ncol(obj), nlevels(obj$celltype)))
cat("\n输出文件清单（位于 output/ 目录）：\n")
cat("  01_celltype_map_applied.csv      - 实际应用的注释表（cluster | celltype | n_cells）\n")
cat("  02_annotation_UMAP.pdf           - 注释后的 UMAP（按细胞类型着色 + 标签）\n")
cat("  （03_annotation_UMAP_by_group.pdf 已删除：与 13_* 组级别堆叠图语义重复）\n")
cat("  （原 04_celltype_composition.pdf 已删除：与 13_* 完全重复）\n")
cat("  05_celltype_summary.csv          - 总体细胞类型计数与占比\n")
cat("  06_celltype_by_group.csv         - 分组 × 细胞类型计数\n")
cat("  07_marker_UMAP.pdf               - QC：核心 marker UMAP 表达网格\n")
cat("  08_marker_dotplot_manual.pdf     - 验证：手写经典 marker 点图\n")
cat("  09_stacked_composition_percent.pdf/.png - 纵向堆叠条形图（比例，样本在 x 轴，标注分组）\n")
cat("  10_stacked_composition_count.pdf/.png   - 纵向堆叠条形图（绝对细胞数，标注分组）\n")
cat("  11_proportion_boxplot_by_group.pdf/.png  - 分组箱式图（每样本细胞类型比例 + Wilcoxon）\n")
cat("  12_ROE_observed_counts.csv       - RO/E 观察计数矩阵（行=样本，列=细胞类型）\n")
cat("  12_ROE_ratio_matrix.csv          - RO/E 比值矩阵\n")
cat("  12_ROE_heatmap.pdf/.png          - RO/E 富集热图（标注数值）\n")
cat("  05_annotated.rds                 - 注释后的 Seurat 对象（meta 含 celltype 列）\n")
cat("\n05_celltype_annotation 完成。\n")
