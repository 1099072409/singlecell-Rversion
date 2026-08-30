# =====================================================================
# 05_celltype_annotation_exclude.R
# 单细胞数据分析流程 - 第五步（CD45+ 免疫细胞类型注释）【exclude 变体】
#   本脚本位于 Rversion/05_celltype_annotation/exclude/ 下，与
#   05_celltype_annotation.R 逻辑、参数、输出完全一致，差异如下：
#   a) 输入改为 Rversion/04_integration_clustering/exclude/04_CD45_integrated.rds
#      （exclude 分析产物：剔除 ypN07/011/012/013/014 后 18 个样本的 CD45+ 整合对象）
#   b) 脚本内含"按指定条件剔除样本"的防御性步骤（EXCLUDE_SAMPLES），
#      即使输入为未剔除数据也能正确执行；输入已剔除时自动确认并跳过
#   c) 输出写到脚本自身所在目录（05_celltype_annotation/exclude 文件夹）
#   细胞注释信息（CELLTYPE_MARKERS / MARKER_SETS）与绘图方式（SeuratExtend：
#   DimPlot2 / FeaturePlot3.grid / DotPlot2 / Heatmap / ClusterDistrBar）与
#   主脚本一致；详见 05_celltype_annotation.R 顶部说明。
#
# 功能（与主脚本一致）：
#   1) 读取第四步整合聚类后的对象（04_CD45_integrated.rds，exclude 版）
#   2) 按指定条件剔除样本（EXCLUDE_SAMPLES，防御性步骤；已剔除则确认跳过）
#   3) 用 AddModuleScore 对每个 cluster 做免疫 marker 打分
#   4) 可视化：打分热图 + marker UMAP 表达图
#   5) 注释流程：读取 celltype_map.csv（与主脚本同结构，cluster → celltype）
#      应用到每个 cluster；首次运行则生成模板
#   6) 输出：注释后对象、打分矩阵、组成统计、注释 UMAP、以及拆分三份的
#      marker 点图（自动 top5 / 手动经典 / 自动+经典 10 列）
#
# 用法：
#   在 RStudio 中打开本文件并点击 Source，或在命令行执行：
#     "C:/Program Files/R/R-4.4.3/bin/Rscript.exe" 05_celltype_annotation_exclude.R
#   说明：所有路径按脚本自身所在目录自动推算，可从任意工作目录运行。
#
# 输出文件清单（位于本 exclude/ 目录）：
#   05_annotated.rds                   - 注释后的 Seurat 对象（含打分与 celltype 列）
#   01_celltype_score_by_cluster.csv   - cluster × 细胞类型平均得分矩阵
#   02_score_heatmap.pdf               - 打分热图（注释主要依据）
#   03_marker_UMAP.pdf                 - 关键 marker 的 UMAP 表达图（FeaturePlot3.grid）
#   04_annotation_UMAP.pdf             - 注释后的 UMAP（应用注释后生成）
#   05_celltype_composition.pdf        - 细胞类型组成图（应用注释后生成）
#   06_marker_dotplot_auto.pdf        - Part 1：自动聚类 top5 marker 点图（第一行 CD8T）
#   06_marker_dotplot_manual.pdf      - Part 2：手动指定经典 marker 点图
#   06_marker_dotplot_combined.pdf    - Part 3：自动+经典 marker 组合点图（10 列/类型）
#   02_celltype_summary.csv / 03_celltype_by_group.csv - 组成统计表
#   celltype_map_template.csv          - 注释映射模板（cluster → celltype 待填）
#
# 参考：Rversion/04_integration_clustering/exclude/04_integration_clustering_exclude.R
#       （同名 exclude 变体的设计范式：输入/输出重定向 + 防御性剔除）
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

# 第四步整合聚类对象（CD45+ 免疫细胞）的【exclude 变体】路径。
# 默认：脚本目录上一级的 04_integration_clustering/exclude/04_CD45_integrated.rds
# 原因：本 exclude 脚本对应 04 的 exclude 变体，读取剔除 5 个样本后的 18 样本整合对象；
#       相对路径保证可移植，且两阶段 exclude 流程无缝衔接（输入即上一步 exclude 产物）。
INPUT_RDS <- file.path(script_dir, "..", "..", "04_integration_clustering", "exclude", "04_CD45_integrated.rds")

# 本步输出目录（自动创建）。原因：所有中间结果与图表统一落盘于脚本自身所在目录
# （05_celltype_annotation/exclude 文件夹），与 04_integration_clustering_exclude.R 一致，
# 不走 output/ 子目录，便于区分"全样本"与"剔除样本"两套结果。
OUT_DIR   <- script_dir

# 用户编辑的注释映射表（cluster -> celltype）。
# 不存在时脚本生成 celltype_map_template.csv 供填写，保存为同名文件后重跑即应用。
# 原因：手动注释需人机协同，模板/正式表用同一文件名，逻辑清晰；且本脚本已随附
#       主脚本产出的 celltype_map.csv（cluster 编号与 04 exclude 对象一致），可直接应用。
MAP_PATH  <- file.path(script_dir, "celltype_map.csv")

# -------------------------------------------------------------------
# 0.2 分组/注释相关
# -------------------------------------------------------------------

# 用哪一列作为聚类分组（04 步的主 cluster 列，默认 0.4 分辨率）。
# 原因：项目约定 04 步产物主聚类列名为 "cluster"；如需其它分辨率列在此改。
CLUSTER_COL <- "cluster"

# 物种。原因：marker 命名有物种差异，本项目为人类样本，故默认 "human"；
# 若为小鼠数据改为 "mouse"，并相应调整下方 CELLTYPE_MARKERS 与 MARKER_SETS。
SPECIES <- "human"

# -------------------------------------------------------------------
# 0.3 随机种子与可复现性
# -------------------------------------------------------------------

# 随机种子。原因：AddModuleScore 内部有随机打乱步骤、首遍运行模拟 sample 时也用随机；
# 固定种子保证两次运行打分一致、结果可复现。
SEED <- 123

# 打分使用的 assay。默认 NULL = 使用 DefaultAssay(obj)。
# 原因：AddModuleScore 基于标准化后的 data 层；04 步已 NormalizeData，用默认 assay 即可，
# 留空避免在多 assay 场景误指定。
SCORE_ASSAY <- NULL

# -------------------------------------------------------------------
# 0.4 差异基因（FindAllMarkers）参数 —— 用于每细胞类型 top5 marker
# -------------------------------------------------------------------

# 只保留在某类中"正向"高表达的基因（only.pos）。
# 原因：注释用 marker 关心"该类特异性高表达"，负向基因对身份判定帮助小，故默认 TRUE。
FM_ONLY_POS   <- TRUE

# 至少在 25% 的该类细胞中检测到表达（min.pct）。
# 原因：与 yijian.R 保持一致（min.pct = 0.25），过滤掉仅在极少数细胞表达的低普适基因，
# 降低噪声、聚焦稳定 marker。
FM_MIN_PCT    <- 0.25

# 最小 log2 倍率阈值（logfc.threshold）。
# 原因：与 yijian.R 保持一致（0.25），兼顾灵敏度与特异性，避免把轻微差异误判为 marker。
FM_LOGFC      <- 0.25

# 每个细胞类型取 top N 个差异基因用于点图。
# 原因：结合"注释文件 marker ≤5 个"，5+5≈10 个/类型，正好满足"每细胞约10个marker"；
# 若想更多/更少在此调整。
FM_TOP_N      <- 5

# 注释文件中每种细胞类型最多取几个 marker 放入点图（先于 top5）。
# 原因：注释文件里每种类型常列 5~7 个 marker，全部放入会超过"约10个/类型"的目标；
# 截断到 5 个后再追加 5 个 top5，合计约 10，版面紧凑且信息充分。
ANNOT_MARKER_CAP <- 5

# -------------------------------------------------------------------
# 0.5 SeuratExtend 绘图参数
# -------------------------------------------------------------------

# FeaturePlot3.grid 的点大小（pt.size）。
# 原因：CD45+ 数据约 8 万细胞，点太小看不清表达分布、太大互相遮挡；
# 0.5 在密度与可读性之间取得平衡（与 Lesson-3 教程示例一致）。
FP_PT_SIZE  <- 0.5

# FeaturePlot3.grid 的配色方案（color）。
# 原因："ryb"（红-黄-蓝）三通道渐变对"低-中-高"表达层次区分直观，是 SeuratExtend 推荐默认。
FP_COLOR    <- "ryb"

# FeaturePlot3.grid 每行展示的基因数（ncol）。
# 原因：4 列网格在 A4/常用 PDF 宽度下每幅图大小合适、不易过密。
FP_NCOL     <- 4

# DotPlot2 点图配色方案（color_scheme）。
# 原因：参考 pythonver/08_visualization.py 的 scanpy Reds 风格，白→红 Sequential 更符合点图习惯。
DOT_COLOR_SCHEME <- "Reds"

# DotPlot2 是否绘制点边框（border）。
# 原因：FALSE 去除边框更整洁；细胞数多时点密，边框会显得杂乱。
DOT_BORDER       <- FALSE

# DotPlot2 是否显示网格线（show_grid）。
# 原因：FALSE 去除网格减少视觉干扰，便于聚焦点与色。
DOT_SHOW_GRID    <- FALSE

# DotPlot2 是否翻转坐标（flip）。
# 原因：TRUE 使基因名位于横轴、水平可读，避免竖排基因名旋转挤压。
DOT_FLIP         <- TRUE

# 细胞类型在点图中的排列顺序（用于梯形向下分布）。
# 原因：按免疫谱系把 CD8T 放第一行，其余按 T→NK→B→髓系顺序排列。
TYPE_ORDER <- c("CD8T", "CD4T", "T cell", "T Regulatory", "T Exhausted",
                "NK cell", "B cell", "Plasma cell", "Monocyte",
                "Macrophage", "DC", "Mast cell", "Neutrophil")

# 点图尺寸参数（每个基因/每行占位英寸）。
DOT_W_PER_GENE <- 0.45
DOT_H_PER_ROW  <- 0.5
DOT_BASE_H     <- 2

# 打分热图的行标准化方式（scale）。
# 原因："row" 按 cluster 行标准化，便于在同一 cluster 内比较各细胞类型相对高低，
# 比全局标准化更能凸显"哪个类型在该 cluster 得分最高"。
HM_SCALE    <- "row"

# 打分热图配色（蓝-白-红）。原因：蓝低红高、白居中，符合打分热图常规语义。
HM_COLOR    <- c("#3C5488", "white", "#E64B35")

# 细胞类型组成图是否堆叠（stack）。原因：TRUE 时按组（ypN0/ypN+）画比例堆叠条，
# 直观对比两组免疫组成差异；FALSE 则画并排绝对计数。
COMP_STACK  <- TRUE

# -------------------------------------------------------------------
# 0.6 SeuratExtend 安装参数
# -------------------------------------------------------------------

# SeuratExtend 的 GitHub 仓库。原因：该包仅发布于 GitHub（huayc09/SeuratExtend），
# CRAN 无此包，必须指定仓库地址安装。
SE_REPO     <- "huayc09/SeuratExtend"

# 安装时是否升级依赖（upgrade）。原因："never" 避免安装过程因升级大量依赖而失败/变慢，
# 在已装 Seurat/ComplexHeatmap 等的机器上最稳妥。
SE_UPGRADE  <- "never"

# 是否打印详细安装/运行信息（verbose）。原因：FALSE 保持终端输出简洁；
# 排错时改为 TRUE 可见更多过程信息。
VERBOSE     <- FALSE

# 载入手动整理的经典 marker 字典（来自 细胞注释.md）。
# 原因：避免运行时解析 markdown，提高稳定性；MANUAL_MARKERS 供 Part 2/3 点图使用。
MANUAL_MARKERS_FILE <- file.path(script_dir, "manual_markers.R")
if (file.exists(MANUAL_MARKERS_FILE)) {
  source(MANUAL_MARKERS_FILE)
} else {
  cat("警告: 未找到 manual_markers.R，手动点图将使用 CELLTYPE_MARKERS 兜底。\n")
  MANUAL_MARKERS <- list()
}

# -------------------------------------------------------------------
# 0.7 样本剔除条件（按指定条件剔除样本）
# -------------------------------------------------------------------
# 需剔除的样本标签（与 03_extract_cd45/exclude、04_integration_clustering/exclude 剔除口径一致）：
#   ypN07、ypN011、ypN012、ypN013、ypN014（均为 ypN0 组样本）
# 本脚本读取的输入已是剔除后的对象，此步骤为防御性执行：
#   - 若对象中仍含这些样本 → 按条件剔除
#   - 若对象中已无这些样本 → 打印确认信息并继续（不重复误删）
EXCLUDE_SAMPLES <- c("ypN07", "ypN011", "ypN012", "ypN013", "ypN014")
EXCLUDE_BY_COL  <- "sample_id"  # 按 meta.data 中哪一列进行样本剔除（04 产物中为 sample_id）


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
# 返回 TRUE/FALSE，供后续决定是否 stop。
# 原因：SeuratExtend 不在 CRAN，且安装可能失败（无网络/缺 Rtools），故单独处理并兜底。
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

# 依次确保基础依赖可用（dplyr 数据整理 / ggplot2 绘图 / patchwork 拼图 /
# pheatmap 作为 SeuratExtend::Heatmap 的兜底热图）。
for (p in c("Seurat", "dplyr", "ggplot2", "patchwork", "pheatmap", "remotes")) ensure_pkg(p)

# 安装并加载 SeuratExtend（本脚本主要绘图方式）。失败则明确报错并退出。
if (!ensure_seuratextend()) {
  stop("SeuratExtend 安装/加载失败，无法继续（请确认网络可访问 GitHub 且已安装 Rtools）。")
}

# 加载核心包。注意：SeuratExtend 提供 DimPlot2/DotPlot2/FeaturePlot3/Heatmap/ClusterDistrBar
# 等增强绘图函数，须在 Seurat 之后加载以使用其增强版本。
library(Seurat)        # 单细胞分析主包（AddModuleScore / FindAllMarkers / readRDS 等）
library(SeuratExtend)  # 主要绘图方式：DimPlot2 / DotPlot2 / FeaturePlot3.grid / Heatmap / ClusterDistrBar
library(dplyr)         # 数据整理（group_by / summarise / mutate 等）
library(ggplot2)        # 兜底绘图（组成图 ggplot 兜底）
library(patchwork)      # 图组合（FeaturePlot3.grid 内部已用，显式加载保险）
library(pheatmap)       # SeuratExtend::Heatmap 失败时的兜底热图

# 创建输出目录（已存在则不报错）
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 打印运行环境信息
cat(sprintf("工作目录: %s\n", getwd()))
cat(sprintf("输出目录: %s\n", normalizePath(OUT_DIR)))

# 检查可用物理内存并给出提示（Windows 下通过 wmic 查询；查询失败时自动跳过）。
# 原因：本步需读取整合对象（约 1~2GB）并做多步计算，内存不足易卡死，提前提示。
tryCatch({
  w  <- system("wmic OS get FreePhysicalMemory /value", intern = TRUE)   # 可用内存(KB)
  fr <- as.numeric(sub(".*=", "", w[grepl("FreePhysicalMemory", w)]))
  if (length(fr) == 1 && !is.na(fr)) {
    free_gb <- fr / 1024^2
    cat(sprintf("当前可用物理内存: %.1f GB\n", free_gb))
    if (free_gb < 8) cat("提示: 本步建议可用内存 >=8GB。\n")
  }
}, error = function(e) invisible(NULL))

# Seurat v5 并行/全局对象大小限制（Future 框架）。
# 原因：默认上限较小，读取大对象或并行计算时会报 "serialized object too large"，
# 提高到 8GB 以避免该问题。
options(future.globals.maxSize = 8 * 1024^3)


## ---- 1b. 辅助函数（把长 if/tryCatch 块拆成独立小块，便于单独修改）----

# 说明：原版把"打分热图 / marker UMAP / 注释应用 / 注释后各图 / 首遍参考图"等
#       长逻辑写在一串 if/tryCatch 里（最长 100+ 行），难以定位与修改。
#       现将每个独立任务拆成一个函数：函数名即任务名，找到即可单独改；
#       所有可调参数仍集中在顶部 ## 0. 配置 中定义（函数内部直接引用全局配置）。
#       本节的每个函数在行为上与改动前完全一致，输出文件与图不变。

# 1b-1 打分热图：SeuratExtend::Heatmap 为主，pheatmap 兜底（对应原 4.5 节）
plot_score_heatmap <- function(mat_hm, OUT_DIR) {
  pdf(file.path(OUT_DIR, "02_score_heatmap.pdf"),
      width = 9, height = max(5, nrow(mat_hm) * 0.35 + 2))
  tryCatch({
    # SeuratExtend 增强热图（主绘图方式）
    SeuratExtend::Heatmap(mat_hm,
      cluster_rows = FALSE,        # 行保持 cluster 顺序（不打乱，便于对照聚类编号）
      cluster_columns = TRUE,      # 列按细胞类型聚类，凸显相似类型
      show_row_names = TRUE, show_column_names = TRUE)
  }, error = function(e) {
    # 兜底：pheatmap 画同口径热图，保证一定出图
    cat("  SeuratExtend::Heatmap 失败，回退 pheatmap:", conditionMessage(e), "\n")
    pheatmap::pheatmap(mat_hm,
      scale = HM_SCALE,
      cluster_cols = TRUE, cluster_rows = FALSE,
      color = colorRampPalette(HM_COLOR)(100),
      main = "Average marker score by cluster",
      fontsize = 10)
  })
  dev.off()
  cat("  已保存图表: 02_score_heatmap.pdf\n")
}

# 1b-2 marker UMAP 表达网格：FeaturePlot3.grid，带 color/ncol 参数回退（对应原 5.1 节）
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
  ggsave(file.path(OUT_DIR, "03_marker_UMAP.pdf"),
         p_feat,
         width  = 5.5 * FP_NCOL,                              # 画布宽随列数自适应
         height = 5.5 * ceiling(length(feat_genes) / FP_NCOL), # 高随行数自适应
         limitsize = FALSE)
  cat("  已保存图表: 03_marker_UMAP.pdf\n")
}

# 1b-3 应用注释映射表：读 celltype_map.csv，逐 cluster 写入 celltype 列（对应原 6 节 if 分支）
apply_celltype_map <- function(obj, MAP_PATH, CLUSTER_COL) {
  cat("检测到注释映射表 celltype_map.csv，正在应用...\n")
  map <- read.csv(MAP_PATH, stringsAsFactors = FALSE)
  if (!all(c("cluster", "celltype") %in% colnames(map))) {
    stop("celltype_map.csv 需要包含 cluster 和 celltype 两列")
  }
  map$celltype <- trimws(as.character(map$celltype))
  # 若某 cluster 的 celltype 为空，自动回退用 top_hit 作为初始注释，保证即使未手动填写也能先出一版。
  if ("top_hit" %in% colnames(map)) {
    empty_idx <- is.na(map$celltype) | map$celltype == ""
    if (sum(empty_idx) > 0) {
      map$celltype[empty_idx] <- map$top_hit[empty_idx]
      cat(sprintf("提示: %d 个 cluster 的 celltype 为空，已自动采用 top_hit 建议作为初始注释。\n", sum(empty_idx)))
      cat("      如需精调，请编辑 celltype_map.csv 中的 celltype 列后重新运行本脚本。\n")
    }
  }
  # 逐 cluster 应用注释（仍为空的 cluster 标为 Unknown）。
  obj$celltype <- "Unknown"
  for (i in seq_len(nrow(map))) {
    if (is.na(map$celltype[i]) || map$celltype[i] == "") next
    obj$celltype[obj[[CLUSTER_COL]] == as.character(map$cluster[i])] <- map$celltype[i]
  }
  cat("注释结果（细胞类型分布）：\n")
  print(table(obj$celltype, useNA = "ifany"))
  unknown_n <- sum(obj$celltype == "Unknown")
  if (unknown_n > 0) cat(sprintf("注意: %d 个细胞未注释（Unknown），请在 celltype_map.csv 中补全后重跑\n", unknown_n))
  obj   # 返回更新后的对象
}

# 1b-4 生成注释模板表（首次运行，对应原 6 节 else 分支）
write_map_template <- function(map_template, OUT_DIR, script_dir) {
  cat("尚未找到 celltype_map.csv（首次运行）。\n")
  cat("请先查看 02_score_heatmap.pdf 与 03_marker_UMAP.pdf，\n")
  cat("确认每个 cluster 的身份后，编辑 celltype_map_template.csv 中的 celltype 列，\n")
  cat("保存为 celltype_map.csv（与本脚本同目录），再重新运行本脚本即可应用注释。\n")
  write.csv(map_template, file.path(OUT_DIR, "celltype_map_template.csv"), row.names = FALSE, na = "")
  write.csv(map_template, file.path(script_dir, "celltype_map.csv"), row.names = FALSE, na = "")
  cat("  已生成注释模板: celltype_map.csv（脚本目录）与 celltype_map_template.csv（本目录）\n")
}

# 1b-5 注释 UMAP：整体 + 按分组两幅（对应原 7a / 7a-2 节）
plot_annotation_umap <- function(obj, OUT_DIR) {
  p_anno <- DimPlot2(obj,
    features = "celltype",   # 用 celltype 作为着色变量
    label    = TRUE,         # 显示类型标签
    box      = TRUE,         # 标签加框更清晰
    repel    = TRUE,         # 标签自动避让避免重叠
    theme    = NoLegend())   # 右侧图例占位去除（标签已直接标在图上）
  ggsave(file.path(OUT_DIR, "04_annotation_UMAP.pdf"), p_anno, width = 9, height = 7)
  cat("  已保存图表: 04_annotation_UMAP.pdf\n")
  p_anno_grp <- DimPlot2(obj,
    features  = "celltype",
    split.by  = "group",     # 按 group 拆成多面板
    label     = TRUE, box = TRUE, repel = TRUE)
  ggsave(file.path(OUT_DIR, "04_annotation_UMAP_by_group.pdf"), p_anno_grp, width = 16, height = 7)
  cat("  已保存图表: 04_annotation_UMAP_by_group.pdf\n")
}

# 1b-6 细胞类型组成统计表：总体 + 按分组（对应原 7b / 7b-2 节）
write_composition_stats <- function(obj, OUT_DIR) {
  ct_summary <- obj@meta.data %>%
    dplyr::group_by(celltype) %>%
    dplyr::summarise(n_cells = dplyr::n(),
                     pct     = 100 * dplyr::n() / ncol(obj),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_cells))
  write.csv(ct_summary, file.path(OUT_DIR, "02_celltype_summary.csv"), row.names = FALSE)
  cat("  已保存: 02_celltype_summary.csv\n")
  ct_group <- obj@meta.data %>%
    dplyr::group_by(celltype, group) %>%
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop")
  write.csv(ct_group, file.path(OUT_DIR, "03_celltype_by_group.csv"), row.names = FALSE)
  cat("  已保存: 03_celltype_by_group.csv\n")
}

# 1b-7 细胞类型组成图：ClusterDistrBar 为主，ggplot 兜底（对应原 7c 节）
plot_composition <- function(obj, OUT_DIR) {
  p_comp <- tryCatch({
    ClusterDistrBar(origin   = obj$group,     # 分组来源（ypN0/ypN+ 或 All）
                    cluster  = obj$celltype,  # 细胞类型
                    flip     = TRUE,          # 横向条，类型名更易读
                    stack    = COMP_STACK)    # 堆叠比例（见配置 COMP_STACK 注释）
  }, error = function(e) {
    cat("  ClusterDistrBar 失败，回退 ggplot:", conditionMessage(e), "\n")
    comp_df <- obj@meta.data %>%
      dplyr::count(group, celltype) %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(prop = n / sum(n))
    ggplot(comp_df, aes(x = group, y = prop, fill = celltype)) +
      geom_col(position = "fill") +
      scale_fill_brewer(palette = "Set3") +
      theme_minimal() +
      labs(title = "Celltype composition by group", x = "Group", y = "Proportion", fill = "Celltype")
  })
  ggsave(file.path(OUT_DIR, "05_celltype_composition.pdf"), p_comp, width = 8, height = 6)
  cat("  已保存图表: 05_celltype_composition.pdf\n")
}

# 1b-8 构建 celltype_detailed：归一化别名 + 对 T 细胞拆 CD8T/CD4T
# 原因：需求要求"第一行必须是 CD8T"；原 top_hit 只有 "Tcell"，需用 CD8T/CD4T 模块打分逐细胞拆分。
derive_celltype_detailed <- function(obj, CLUSTER_COL) {
  ct <- as.character(obj$celltype)
  # 别名归一化到 TYPE_ORDER 中的规范名
  alias <- list(
    Tcell = "T cell", `T cell` = "T cell", T = "T cell",
    B = "B cell", `B cell` = "B cell",
    Plasma = "Plasma cell", `Plasma cell` = "Plasma cell",
    Mast = "Mast cell", `Mast cell` = "Mast cell",
    NK = "NK cell", `NK cell` = "NK cell",
    Neutro = "Neutrophil", Neutrophil = "Neutrophil",
    Macro = "Macrophage", Macrophage = "Macrophage",
    Mono = "Monocyte", Monocyte = "Monocyte", DC = "DC"
  )
  ct <- unname(sapply(ct, function(x) if (x %in% names(alias)) alias[[x]] else x))
  # 对 T 细胞按 CD8 vs CD4 打分逐细胞拆分
  if (all(c("Score_CD8T1", "Score_CD4T1") %in% colnames(obj@meta.data))) {
    s8 <- obj$Score_CD8T1; s4 <- obj$Score_CD4T1
    isT <- ct %in% c("T cell", "Tcell", "T")
    # 仅当 CD8/CD4 中至少一个为正时才拆分，避免双低 T 细胞被强制归类
    ct[isT] <- ifelse(s8[isT] >= s4[isT] & pmax(s8[isT], s4[isT]) > 0,
                      "CD8T", ifelse(pmax(s8[isT], s4[isT]) > 0, "CD4T", "T cell"))
  }
  present <- intersect(TYPE_ORDER, unique(ct))
  obj$celltype_detailed <- factor(ct, levels = present)
  obj
}

# 1b-9 自动 topN 差异基因（FindAllMarkers），命名 list 按 TYPE_ORDER 排序
get_auto_markers <- function(obj, ident_col, top_n = FM_TOP_N) {
  cat(sprintf("  正在计算 %s 的 top%d 差异基因（FindAllMarkers）...\n", ident_col, top_n))
  Idents(obj) <- ident_col
  fm <- FindAllMarkers(obj,
    only.pos        = FM_ONLY_POS,
    min.pct         = FM_MIN_PCT,
    logfc.threshold = FM_LOGFC)
  top_df <- fm %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = top_n, with_ties = FALSE)
  lst <- split(top_df$gene, top_df$cluster)
  lst <- lst[intersect(TYPE_ORDER, names(lst))]
  lapply(lst, function(g) g[g %in% rownames(obj)])
}

# 1b-10 由 MANUAL_MARKERS 构建分组 list（按 TYPE_ORDER，过滤缺失基因）
get_manual_markers <- function(obj) {
  if (!exists("MANUAL_MARKERS") || length(MANUAL_MARKERS) == 0) return(list())
  lst <- MANUAL_MARKERS[intersect(TYPE_ORDER, names(MANUAL_MARKERS))]
  lapply(lst, function(g) g[g %in% rownames(obj)])
}

# 1b-11 组合 auto + manual（Part3：10 列/类型）
combine_auto_manual <- function(auto_lst, manual_lst) {
  common <- intersect(names(auto_lst), names(manual_lst))
  out <- list()
  for (ct in common) out[[ct]] <- c(auto_lst[[ct]], manual_lst[[ct]])
  out[intersect(TYPE_ORDER, names(out))]
}

# 1b-12 通用点图绘制（梯形向下布局：行与列块均按 TYPE_ORDER）
plot_dot_generic <- function(obj, grouped_features, fname, OUT_DIR) {
  if (length(grouped_features) == 0) {
    cat("  注意: 未整理出可用于", fname, "的 marker，跳过\n")
    return(invisible(NULL))
  }
  p <- DotPlot2(obj,
    features     = grouped_features,
    group.by     = "celltype_detailed",
    color_scheme = DOT_COLOR_SCHEME,
    border       = DOT_BORDER,
    show_grid    = DOT_SHOW_GRID,
    flip         = DOT_FLIP)
  n_cols <- sum(lengths(grouped_features))
  n_rows <- length(grouped_features)
  ggsave(file.path(OUT_DIR, fname), p,
         width  = max(10, n_cols * DOT_W_PER_GENE),
         height = max(6,  n_rows * DOT_H_PER_ROW + DOT_BASE_H),
         limitsize = FALSE)
  cat("  已保存图表:", fname, "\n")
}

# 1b-13 三个点图一次计算并落盘（避免 FindAllMarkers 重复计算）
run_three_dotplots <- function(obj, OUT_DIR) {
  auto_lst   <- get_auto_markers(obj, "celltype_detailed")
  manual_lst <- get_manual_markers(obj)
  plot_dot_generic(obj, auto_lst,                       "06_marker_dotplot_auto.pdf",     OUT_DIR)
  plot_dot_generic(obj, manual_lst,                     "06_marker_dotplot_manual.pdf",   OUT_DIR)
  plot_dot_generic(obj, combine_auto_manual(auto_lst, manual_lst),
                                                "06_marker_dotplot_combined.pdf", OUT_DIR)
}


## ---- 2. 读入对象与基础检查（仅检查，本节约不作图）----

# 2.1 检查输入文件是否存在，缺失则终止并给出明确路径。
if (!file.exists(INPUT_RDS)) stop("找不到输入 rds: ", INPUT_RDS)

# 2.2 读取第四步的整合聚类对象（约 8 万细胞，读取需 1~3 分钟）。
cat("正在读取 04_CD45_integrated.rds（约需 1~3 分钟）...\n")
obj <- readRDS(INPUT_RDS)
cat(sprintf("读入完成：%d 个基因 x %d 个细胞\n", nrow(obj), ncol(obj)))

# 2.3 检查主 cluster 列是否存在（04 步产物应包含）。
if (!CLUSTER_COL %in% colnames(obj@meta.data)) {
  stop("对象中找不到主 cluster 列 '", CLUSTER_COL, "'，请先运行 04_integration_clustering_exclude.R")
}
cat(sprintf("当前聚类（%s）簇数: %d\n", CLUSTER_COL, length(unique(obj[[CLUSTER_COL]]))))

# 2.4 确保有 data 层（AddModuleScore 基于标准化数据；04 步已 NormalizeData）。
if (!"data" %in% Layers(obj, assay = DefaultAssay(obj))) {
  cat("对象缺少 data 层，执行 NormalizeData（LogNormalize）...\n")
  obj <- NormalizeData(obj, verbose = FALSE)
}

# 若对象含 group 列（ypN0/ypN+）则用之；否则构造占位 "All"，保证下游组成图不报错。
if (!"group" %in% colnames(obj@meta.data)) {
  cat("提示: 对象中无 'group' 列，组成图将退化为单一分组 'All'。\n")
  obj$group <- "All"
}

## ---- 2.6 按指定条件剔除样本（防御性，用户需求专项） ----
# 输入 04_integration_clustering/exclude/04_CD45_integrated.rds 已是剔除 EXCLUDE_SAMPLES
# （ypN07/ypN011/ypN012/ypN013/ypN014）后的对象。本步骤确保脚本自包含：
#   1) 若对象中仍含这些样本 → 按 EXCLUDE_BY_COL 列剔除，并打印剔除前后细胞数；
#   2) 若对象中已无这些样本 → 打印确认信息后继续（不重复误删、不报错中止）；
#   3) 若剔除列缺失 → 打印警告并跳过（不中断主流程）。
cat("\n==== 2.6 按指定条件剔除样本 ====\n")
cat(sprintf("  排除样本列表: %s\n", paste(EXCLUDE_SAMPLES, collapse = ", ")))
if (!EXCLUDE_BY_COL %in% colnames(obj@meta.data)) {          # 剔除列不存在时
  cat(sprintf("  警告: meta.data 中无 %s 列，跳过样本剔除检查（继续运行）。\n", EXCLUDE_BY_COL))
} else {                                                     # 剔除列存在时
  all_ids <- unique(as.character(obj@meta.data[[EXCLUDE_BY_COL]])) # 对象中全部样本标签
  present <- intersect(EXCLUDE_SAMPLES, all_ids)             # 需要剔除且实际存在的样本
  missing <- setdiff(EXCLUDE_SAMPLES, all_ids)               # 列表中但不在对象中的样本
  if (length(present) > 0) {                                 # 仍有需剔除的样本
    cat(sprintf("  剔除前: %d 个样本 / %d 个细胞\n", length(all_ids), ncol(obj)))
    cat(sprintf("  剔除样本: %s\n", paste(present, collapse = ", ")))
    obj <- subset(obj, subset = !(get(EXCLUDE_BY_COL) %in% EXCLUDE_SAMPLES)) # 执行剔除
    n_after <- length(unique(as.character(obj@meta.data[[EXCLUDE_BY_COL]])))  # 剔除后样本数
    cat(sprintf("  剔除后: %d 个样本 / %d 个细胞\n", n_after, ncol(obj)))
  } else {                                                   # 已无待剔除样本
    cat(sprintf("  确认输入数据已不含排除样本（%s），无需剔除。\n",
                paste(EXCLUDE_SAMPLES, collapse = ", ")))
  }
  if (length(missing) > 0 && length(present) == 0) {         # 全部不在对象中
    cat("  （以上样本在 04 exclude 步骤已剔除，符合预期）\n")
  } else if (length(missing) > 0) {                          # 部分存在部分缺失
    cat(sprintf("  注意: 以下样本不在数据中（可能已在 04 剔除）: %s\n",
                paste(missing, collapse = ", ")))
  }
}


## ---- 3. 细胞注释信息整理（依据 细胞注释.md + Lesson 5/6）----

# 本节的 CELLTYPE_MARKERS 与 MARKER_SETS 即"按参考文件定义整理出的细胞注释信息"，
# 后续打分与作图均以此为准。
#
# 数据来源说明（用户确认参考 细胞注释.md + Lesson-5-Core-scRNA.md + Lesson-6-Advanced-scRNA.md）：
#   * 细胞注释.md 提供了 T/B/NK/髓系/Mast/Plasma 等各类型的 marker 基因表；
#   * Lesson-6 CellChat 分组（B cell / CD4 T / CD8 T / DC / Mono CD14 / Mono FCGR3A /
#     NK cell / Platelet）印证本项目免疫细胞的命名口径；
#   * 本项目实际注释结果（celltype_map.csv）共 7 种类型，故此处按这 7 种组织。

# 3.1 各细胞类型的"特征 marker"（取自 细胞注释.md 对应表，拼写已校正，如
#     FCG3RA→FCGR3A、TYOBP→TYROBP）。仅作注释/点图展示用，不要求全部存在。
CELLTYPE_MARKERS <- list(
  # T 细胞通用 + 代表性 CD4/CD8/记忆标志（细胞注释.md: T Cells=CD3D,CD3E; CD4+T=CD3D,CD4; CD8+T=CD3D,CD8A,GZMK）
  "T cell"      = c("CD3D", "CD3E", "CD4", "CD8A", "GZMK", "IL7R"),
  # B 细胞（细胞注释.md: B cells=CD19,CD79A,CD79B,MS4A1[CD20]）
  "B cell"      = c("CD19", "CD79A", "CD79B", "MS4A1"),
  # 浆细胞（细胞注释.md: Plasma cells=IGHG1,MZB1,SDC1,CD79A; Plasma=SDC1,IGHG1,MZB1,TNFRSF17）
  "Plasma cell" = c("IGHG1", "MZB1", "SDC1", "CD79A", "PRDM1", "JCHAIN"),
  # NK 细胞（细胞注释.md: NK=KLRD1,GNLY; NK cells=FGFBP2,FCGR3A,CX3CR1,GNLY,NKG7,TYROBP,PRF1）
  "NK cell"     = c("KLRD1", "GNLY", "NKG7", "CD247", "FCGR3A", "PRF1"),
  # 巨噬细胞（细胞注释.md: Monocytes and macrophages=CST3,LYZ,CD68,CD163,CD14; Myeloid=LYZ,CST3,CD14,CD68）
  "Macrophage"  = c("LYZ", "CST3", "CD14", "CD68", "CD163", "C1QA"),
  # 中性粒细胞（细胞注释.md: Neutrophils=CST3,LYZ,FCGR3B,CSF3R; Myeloid 子表 Neutrophil=MPO,AZU1,ELANE,SRGN,ENO1）
  "Neutrophil"  = c("CST3", "LYZ", "FCGR3B", "CSF3R", "MPO", "S100A8", "S100A9"),
  # 肥大细胞（细胞注释.md: mast cell=CST3,KIT,CPA3; Mast=TPSAB1,TPSB2）
  "Mast cell"   = c("TPSAB1", "TPSB2", "CPA3", "KIT", "CST3")
)
cat(sprintf("已整理 %d 种细胞类型的注释 marker 信息。\n", length(CELLTYPE_MARKERS)))

# 3.2 免疫细胞 marker 基因集（人类，CD45+ 免疫细胞常见亚群）——用于 AddModuleScore 打分。
# 名称用作打分列名（避免空格），基因须在数据中存在，缺失的会自动跳过。
MARKER_SETS <- list(
  Tcell    = c("CD3D", "CD3E", "CD3G"),                        # T 细胞通用
  CD8T     = c("CD8A", "CD8B"),                                # CD8 T 细胞
  CD4T     = c("CD4", "IL7R"),                                 # CD4 T 细胞
  NK       = c("NKG7", "GNLY", "KLRD1", "NCAM1"),              # NK 细胞
  B        = c("MS4A1", "CD79A", "CD79B"),                     # B 细胞
  Plasma   = c("MZB1", "SDC1", "JCHAIN"),                      # 浆细胞
  Mono     = c("LYZ", "FCGR3A", "CSF1R"),                      # 单核细胞
  Macro    = c("CD68", "C1QA", "C1QB", "CSF1R"),               # 巨噬细胞
  DC       = c("CLEC9A", "ITGAX", "FCER1A", "BATF3"),          # 树突状细胞
  Mast     = c("TPSAB1", "TPSB2", "CPA3"),                     # 肥大细胞
  Neutro   = c("FCGR3B", "S100A8", "S100A9")                   # 中性粒细胞
)

# 3.3 关键 marker 用于 UMAP 表达图（从 MARKER_SETS 与 CELLTYPE_MARKERS 中挑选代表性基因）。
# 原因：精选出覆盖各主要免疫亚群的核心基因，既能确认聚类身份又不过于冗长。
FEATURE_GENES <- c("CD3D", "CD8A", "CD4", "NKG7", "MS4A1", "MZB1",
                   "LYZ", "CD68", "CLEC9A", "TPSAB1", "FCGR3B")


## ---- 4. 免疫 marker 打分（AddModuleScore）+ 打分热图（SeuratExtend::Heatmap）----

# 4.1 逐细胞类型打分：对每个 marker 基因集，用 AddModuleScore 计算
#     该细胞类型的"特征得分"（高得分 = 更像该细胞类型）。
cat("\n==== 4. 免疫 marker 打分 ====\n")
used_sets <- list()                                  # 记录实际打分成功的基因集
for (ct in names(MARKER_SETS)) {                     # 遍历各细胞类型
  genes <- MARKER_SETS[[ct]]                         # 该类型 marker 基因
  present <- genes[genes %in% rownames(obj)]         # 仅保留数据中存在的基因
  if (length(present) == 0) {                        # 全部缺失则跳过
    cat(sprintf("  警告: %s 的 marker 基因全部缺失，跳过该类型打分\n", ct))
    next
  }
  if (length(present) < length(genes)) {             # 部分缺失时提示
    cat(sprintf("  %s: 缺失 %d 个基因（%s），使用剩余 %d 个\n",
                ct, length(genes) - length(present),
                paste(setdiff(genes, present), collapse = "/"), length(present)))
  }
  obj <- AddModuleScore(obj,                          # 执行特征打分
    features = list(present),                        # 基因集（列表形式）
    name = paste0("Score_", ct),                     # 打分列名前缀 -> Score_Tcell 等
    seed = SEED,                                     # 随机种子（可复现）
    assay = if (is.null(SCORE_ASSAY)) DefaultAssay(obj) else SCORE_ASSAY)
  used_sets[[ct]] <- present                         # 记录成功打分的基因集
  cat(sprintf("  %s 打分完成（%d 个基因）\n", ct, length(present)))
}

# 4.2 收集所有打分列名（AddModuleScore 会在 name 后追加 "1"，如 Score_Tcell1）。
score_cols <- paste0("Score_", names(used_sets), "1")
cat("打分列:", paste(score_cols, collapse = ", "), "\n")

# 4.3 计算每个 cluster 的各细胞类型平均得分（打分矩阵，注释的核心依据）。
score_mat <- obj@meta.data %>%
  dplyr::group_by(.data[[CLUSTER_COL]]) %>%          # 按 cluster 分组
  dplyr::summarise(
    dplyr::across(dplyr::all_of(score_cols), mean),   # 各打分列取均值
    n_cells = dplyr::n(),                            # 每簇细胞数
    .groups = "drop"
  )
colnames(score_mat)[2:(length(score_cols) + 1)] <- names(used_sets)  # 打分列改名为细胞类型名

# 4.4 保存打分矩阵 CSV。
write.csv(score_mat, file.path(OUT_DIR, "01_celltype_score_by_cluster.csv"), row.names = FALSE)
cat("  已保存: 01_celltype_score_by_cluster.csv（cluster × 细胞类型平均得分）\n")

# 4.5 打分热图（本节的"对应作图"）：行为 cluster、列为细胞类型，颜色深浅表示得分高低。
#     每个 cluster 得分最高的细胞类型即其最可能的身份（先看结果再定注释）。
#     已拆分为独立函数 plot_score_heatmap()（见 1b 节；SeuratExtend::Heatmap 为主、pheatmap 兜底）。
mat_hm <- as.matrix(score_mat[, names(used_sets), drop = FALSE])
rownames(mat_hm) <- score_mat[[CLUSTER_COL]]
plot_score_heatmap(mat_hm, OUT_DIR)

# 4.6 自动推断每个 cluster 的最可能细胞类型（得分最高者），作为注释建议。
top_hit   <- apply(mat_hm, 1, function(x) names(used_sets)[which.max(x)])
top_score <- apply(mat_hm, 1, max)


## ---- 5. marker UMAP 表达图（FeaturePlot3.grid，SeuratExtend 主要绘图）----

# 5.1 从 FEATURE_GENES 中挑选数据中存在的基因，用 SeuratExtend 的 FeaturePlot3.grid
#     一次性画出多基因 UMAP 表达网格，直观查看 marker 在 UMAP 上的表达分布。
#     已拆分为独立函数 plot_marker_umap()（见 1b 节；含 color/ncol 参数回退）。
cat(sprintf("\n==== 5. marker UMAP 表达图 ====\n"))
feat_genes <- FEATURE_GENES[FEATURE_GENES %in% rownames(obj)]   # 过滤缺失基因
cat(sprintf("  使用 %d 个代表性 marker 基因作图\n", length(feat_genes)))
if (length(feat_genes) > 0) {
  plot_marker_umap(obj, feat_genes, OUT_DIR)
}


## ---- 6. 应用注释（先看结果，再填表确认；本步为分析，无图）----

# 6.1 生成注释建议表：cluster | n_cells | top_hit（打分最高的类型建议）| celltype（待填）。
map_template <- data.frame(
  cluster  = score_mat[[CLUSTER_COL]],    # cluster 编号
  n_cells  = score_mat$n_cells,           # 簇细胞数
  top_hit  = top_hit,                     # 自动建议的细胞类型（供参考）
  celltype = NA_character_,               # 待用户填写的最终注释
  stringsAsFactors = FALSE
)
cat("\n==== 6. 注释流程 ====\n")
if (file.exists(MAP_PATH)) {              # 用户已填写映射表时 -> 应用注释
  obj <- apply_celltype_map(obj, MAP_PATH, CLUSTER_COL)   # 已拆分为函数（见 1b 节）
} else {                                  # 首次运行 -> 生成模板并提示
  write_map_template(map_template, OUT_DIR, script_dir)   # 已拆分为函数（见 1b 节）
}

# 6.2 确保 obj$celltype 存在（首遍运行时回退到 top_hit），然后构建 celltype_detailed。
# 原因：Part 1~3 点图统一按 celltype_detailed 分组，必须不依赖用户先填表。
if (!"celltype" %in% colnames(obj@meta.data)) {
  obj$celltype <- top_hit[as.character(obj[[CLUSTER_COL]])]
  obj$celltype[is.na(obj$celltype)] <- "Unknown"
}
obj <- derive_celltype_detailed(obj, CLUSTER_COL)
cat("celltype_detailed 分布:\n")
print(table(obj$celltype_detailed, useNA = "ifany"))


## ---- 7. 注释后可视化与统计（分析 + 对应作图同节）----

# celltype_detailed 已在 section 6 构建完成，本节统一产出注释后图表。
cat("\n==== 7. 注释后可视化与统计 ====\n")
plot_annotation_umap(obj, OUT_DIR)          # 7a: 注释 UMAP（整体 + 按分组）
write_composition_stats(obj, OUT_DIR)       # 7b: 组成统计表（总体 + 按分组）
plot_composition(obj, OUT_DIR)              # 7c: 组成图（ClusterDistrBar 为主，ggplot 兜底）

cat("\n==== 7d. 拆分三份 marker 点图（自动 / 手动 / 自动+经典）====\n")
run_three_dotplots(obj, OUT_DIR)


## ---- 8. 保存注释对象与完成 ----

# 8.1 保存注释后的对象（无论是否已应用注释，均保存；meta 含各 Score_* 打分列与 celltype 列）。
saveRDS(obj, file.path(OUT_DIR, "05_annotated.rds"))
cat(sprintf("  已保存: %s\n", file.path(OUT_DIR, "05_annotated.rds")))

# 打印最终结果摘要
annotated <- "celltype" %in% colnames(obj@meta.data) && sum(obj$celltype != "Unknown", na.rm = TRUE) > 0
cat("\n==== 05 细胞注释（exclude 变体）完成 ====\n")
if (annotated) {
  cat(sprintf("已注释 %d 个细胞、%d 种细胞类型\n",
              ncol(obj), length(unique(obj$celltype[obj$celltype != "Unknown"]))))
} else {
  cat("当前为首次运行：请查看打分热图与 marker UMAP 后填写 celltype_map.csv 并重跑。\n")
}
cat("\n输出文件清单（位于本 exclude/ 目录）：\n")
cat("  05_annotated.rds                   - 注释后的 Seurat 对象（含打分与 celltype 列）\n")
cat("  01_celltype_score_by_cluster.csv   - cluster × 细胞类型平均得分矩阵\n")
cat("  02_score_heatmap.pdf               - 打分热图（注释主要依据）\n")
cat("  03_marker_UMAP.pdf                 - 关键 marker 的 UMAP 表达图（FeaturePlot3.grid）\n")
cat("  04_annotation_UMAP.pdf             - 注释后的 UMAP（应用注释后生成）\n")
cat("  04_annotation_UMAP_by_group.pdf    - 按分组的注释 UMAP（应用注释后生成）\n")
cat("  05_celltype_composition.pdf        - 细胞类型组成图（应用注释后生成）\n")
cat("  06_marker_dotplot_auto.pdf        - Part 1：自动聚类 top5 marker 点图（第一行 CD8T）\n")
cat("  06_marker_dotplot_manual.pdf      - Part 2：手动指定经典 marker 点图\n")
cat("  06_marker_dotplot_combined.pdf    - Part 3：自动+经典 marker 组合点图（10 列/类型）\n")
cat("  02_celltype_summary.csv / 03_celltype_by_group.csv - 组成统计表（应用注释后生成）\n")
cat("  celltype_map_template.csv          - 注释映射模板（cluster → celltype 待填）\n")
cat("\n05_celltype_annotation_exclude 完成。\n")
