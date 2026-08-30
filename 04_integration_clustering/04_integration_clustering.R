# =====================================================================
# 04_integration_clustering.R
# 单细胞数据分析流程 - 第四步（标准化 + 去批次 + 降维 + 聚类）
#
# 功能：
#   1) 读取第三步提取的 CD45+ 细胞对象（03_CD45_positive.rds）
#   2) 标准化：NormalizeData（LogNormalize）+ FindVariableFeatures（vst，3000 高变基因）
#      + ScaleData（可选回归线粒体比例，对齐 yijian.R）
#   3) 去批次：RunPCA 后 RunHarmony（以 sample 为批次变量），消除样本间技术差异
#   3b) PC 选择诊断：ElbowPlot + DimHeatmap（看每个 PC 的 top 基因），然后
#       根据累积方差/标准差自动确定 DIMS（也支持手动覆盖）
#   4) 降维可视化：RunUMAP（基于 harmony）、RunTSNE（可选）
#   5) 聚类：FindNeighbors + FindClusters，支持多分辨率（默认 0.2 / 0.5 / 1.0），
#      各分辨率的聚类结果分别保存在 meta.data 的 cluster_0.2 / cluster_0.5 /
#      cluster_1.0 列中，并指定一个默认分辨率（0.5）作为主 cluster 列
#   6) 输出：整合聚类后的 Seurat 对象（04_CD45_integrated.rds）、聚类统计表
#      （01_clustering_summary.csv / 02_cluster_by_sample.csv）、
#      以及 PCA 诊断图（Elbow / DimHeatmap）、UMAP/聚类图（先多分辨率拼图，再默认分辨率图）
#
# 用法：
#   在 RStudio 中打开本文件并点击 Source，或在命令行执行：
#     "C:/Program Files/R/R-4.4.3/bin/Rscript.exe" 04_integration_clustering.R
#   说明：所有路径按脚本自身所在目录自动推算，可从任意工作目录运行。
#
# 常用修改点：见下方 ## ---- 0. 配置 ---- 中的 CONFIG 参数。
# 参考：Rversion/yijian.R 的 run_integration()（807-835 行，标准化/整合/聚类范式）
# =====================================================================

## ---- 0. 配置（按需修改） ----

# 自动定位脚本所在目录（兼容 Rscript 与 RStudio 两种运行方式）：
#   - RStudio 中优先通过 rstudioapi 获取当前文档路径
#   - 命令行 Rscript 时通过 sys.frames()/--file= 参数获取脚本路径
#   - 全部失败则退回当前工作目录
script_dir <- tryCatch({
  # 方法1：RStudio 中获取
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
    dirname(rstudioapi::getActiveDocumentContext()$path)
  } else {
    # 方法2：通过 sys.frames 获取（Rscript 或 source）
    f <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
    if (is.null(f) || !nzchar(f)) {
      # 方法3：从命令行参数获取
      arg <- commandArgs(trailingOnly = FALSE)
      f <- sub("^--file=", "", grep("^--file=", arg, value = TRUE)[1])
    }
    if (is.null(f) || !nzchar(f) || is.na(f)) {
      # 方法4：都失败则使用当前工作目录
      getwd()
    } else {
      dirname(normalizePath(f, mustWork = FALSE))
    }
  }
}, error = function(e) getwd())

cat(sprintf("脚本目录: %s\n", script_dir))

# ---- 输入/输出路径 ----
INPUT_RDS <- file.path(script_dir, "..", "03_extract_cd45", "output", "03_CD45_positive.rds") # 第三步 CD45+ 细胞对象
OUT_DIR   <- file.path(script_dir, "output")                                                  # 本步输出目录（自动创建）

# ---- 标准化参数 ----
NFEATURES     <- 3000      # 高变基因（HVG）数量（FindVariableFeatures 的 nfeatures）
REGRESS_MT    <- TRUE      # ScaleData 是否回归 percent.mt（TRUE 对齐 yijian.R；FALSE 保留线粒体差异）

# ---- 降维参数 ----
NPCS          <- 50        # PCA 主成分数量（跑足一些，供 ElbowPlot 与方差判断选择维度）
DIMS          <- NULL      # 下游分析（Harmony/UMAP/聚类）使用的主成分范围。
                           #   NULL = 自动选择：先跑 PCA 并输出 ElbowPlot + 方差表，
                           #   再根据累积方差/标准差阈值自动确定（"参数选择先有结果"）。
                           #   也可手动指定，例如 DIMS <- 1:30。
# DIMS 为 NULL 时自动选择的判定阈值：
CUMVAR_THRESHOLD <- 0.85   # 累积方差解释比例阈值（达到该比例所需的最小 PC 数）
STDEV_THRESHOLD  <- 1      # PC 标准差阈值（第一个标准差 < 该值的 PC 视为噪音，其前的 PC 保留）
MIN_PCS          <- 5      # 自动选择的最少 PC 数（下限）
MAX_PCS          <- 40     # 自动选择的最多 PC 数（上限，避免纳入纯噪音维度）

# ---- UMAP / tSNE 参数（扩展，参考 Seurat 最佳实践与 Lesson 教程）----
UMAP_METHOD     <- "uwot"  # UMAP 实现（R 原生 uwot；也可 "umap-learn" 需 Python）
UMAP_METRIC     <- "cosine"# 距离度量（单细胞数据常用 cosine；可选 "euclidean"）
UMAP_NEIGHBORS  <- 25      # n.neighbors：局部邻域大小（越小越强调局部结构）
UMAP_MIN_DIST   <- 0.2     # min.dist：嵌入点间最小距离（越大越分散；0.1/0.3/0.5 可先试图再选）
UMAP_SPREAD     <- 1       # spread：嵌入的有效尺度（配合 min.dist 调节）
TSNE_PERPLEXITY <- 30      # tSNE 困惑度（通常 5~50，越大越全局）

# ---- 去批次（Harmony）参数 ----
RUN_HARMONY   <- TRUE      # 是否运行 Harmony 去批次（TRUE 为默认推荐）
HARMONY_BATCH <- "sample"  # 批次变量（以样本为批次；如需按其他列可改）
HARMONY_THETA <- NULL      # Harmony 批次校正强度 theta；NULL=用默认值(2)。
                           #   如去批次不足（样本仍各自成团），可调大如 5 增强校正。
PLOT_HARMONY_CONVERGENCE <- TRUE # 是否绘制 Harmony 收敛诊断图（图 10_Harmony_convergence.pdf）；
                           #   TRUE 时 RunHarmony 会额外保存收敛过程数据（不影响结果，仅略增内存）。
                           #   若画图失败（不同版本 harmony 槽结构差异）脚本会自动跳过，不影响主流程。

# ---- 聚类参数 ----
RESOLUTIONS        <- c(0.1,0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8,0.9,1.0)  # 要尝试的聚类分辨率（可增删）
CLUSTER_ALGORITHM  <- 1               # 聚类算法：1=Louvain, 2=Leiden(需 leiden 包), 3=SLM, 4=Leiden(需 leiden 包)
CLUSTER_METHOD     <- "igraph"        # 图构建方法（igraph 为推荐默认）
GROUP_SINGLETONS   <- TRUE            # 是否将孤立单细胞簇并入最近簇（避免单细胞噪声簇）

# ---- 其他 ----
RUN_TSNE <- TRUE        # 是否运行 tSNE 降维（较慢，可关闭）
SEED     <- 123         # 随机种子（保证 UMAP/聚类结果可复现）

## ---- 1. 依赖检查与加载 ----

# 通用依赖安装函数：ensure_pkg("包名")，未安装时自动从 CRAN 安装
ensure_pkg <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) return(TRUE)             # 已装则跳过
  cat(sprintf("  正在安装依赖包: %s ...\n", pkg))                    # 提示开始安装
  install.packages(pkg, repos = "https://cloud.r-project.org")        # 从官方镜像安装
  requireNamespace(pkg, quietly = TRUE)                               # 返回是否安装成功
}

# 依次确保需要的包可用（Seurat/dplyr/ggplot2/patchwork/harmony 均已安装）
for (p in c("Seurat", "dplyr", "ggplot2", "patchwork", "harmony")) ensure_pkg(p)

# 加载核心包
library(Seurat)    # 单细胞分析主包
library(dplyr)     # 数据整理（group_by/summarise 等）
library(ggplot2)   # 绘图
library(patchwork) # 图组合
library(harmony)   # 批次整合（RunHarmony）

# 设置随机种子（保证可复现）
set.seed(SEED)

# 创建输出目录（已存在则不报错）
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 打印运行环境信息，便于排查问题
cat(sprintf("工作目录: %s\n", getwd()))                       # 当前工作目录
cat(sprintf("输出目录: %s\n", normalizePath(OUT_DIR)))        # 输出目录

# 检查可用物理内存并给出提示（Windows 下通过 wmic 查询；查询失败时自动跳过，不影响运行）。
# 本步涉及 PCA/Harmony/UMAP 等计算，内存需求较高，提前提醒用户关闭其他程序。
tryCatch({
  w  <- system("wmic OS get FreePhysicalMemory /value", intern = TRUE) # 查询可用内存(KB)
  fr <- as.numeric(sub(".*=", "", w[grepl("FreePhysicalMemory", w)]))  # 提取数值
  if (length(fr) == 1 && !is.na(fr)) {                                 # 查询成功时
    free_gb <- fr / 1024^2                                             # 转 GB
    cat(sprintf("当前可用物理内存: %.1f GB\n", free_gb))                # 打印可用内存
    if (free_gb < 10) {                                                # 低于 10GB 时提示
      cat("提示: 本步计算密集，建议可用内存 >=10GB，并关闭其他程序后再运行。\n")
    }
  }
}, error = function(e) invisible(NULL))                                 # 查询失败静默跳过

# Seurat v5 并行/全局对象大小限制（Harmony/UMAP 等大对象操作需要）
options(future.globals.maxSize = 8 * 1024^3)

## ---- 2. 读入对象并做基础检查 ----

# 2.1 检查输入文件是否存在
if (!file.exists(INPUT_RDS)) stop("找不到输入 rds: ", INPUT_RDS)

# 2.2 读取第三步的 CD45+ 细胞对象（约 8 万细胞，读取需 1~3 分钟）
cat("正在读取 03_CD45_positive.rds（约需 1~3 分钟）...\n")
obj <- readRDS(INPUT_RDS)                          # 读入 Seurat 对象
cat(sprintf("读入完成：%d 个基因 x %d 个细胞\n", nrow(obj), ncol(obj))) # 打印维度
cat("分组分布：\n")                                 # 打印分组分布
print(table(obj$group, useNA = "ifany"))           # ypN0 / ypN+ 细胞数

# 2.3 确保默认 assay 为 RNA（后续标准化基于它）
if (DefaultAssay(obj) != "RNA") {                  # 若默认 assay 不是 RNA
  warning("默认 assay 不是 RNA，已切换为 RNA")      # 提示后切换
  DefaultAssay(obj) <- "RNA"
}

# 2.4 统一 Layer 结构：若存在多个 layer（merge 产物），JoinLayers 合并为单一 counts
cat(sprintf("当前 RNA layer: %s\n", paste(Layers(obj, assay = "RNA"), collapse = ", ")))
if (length(Layers(obj, assay = "RNA")) > 1) {      # 若存在多个 layer
  cat("检测到多个 layer，执行 JoinLayers 合并...\n")
  obj <- JoinLayers(obj, assay = "RNA")            # 合并所有 layer
}

# 2.5 统一样本标识列（03 产物已含 sample 列，此处兜底，缺失时用 orig.ident 生成）
if (!"sample" %in% colnames(obj@meta.data)) {      # 输入对象没有 sample 列时
  obj$sample <- as.character(obj$orig.ident)       # 用 orig.ident（样本文件夹名）生成
}
cat(sprintf("样本数: %d\n", length(unique(obj$sample)))) # 打印样本数（应为 23）

## ---- 3. 标准化（LogNormalize）----

# 3.1 对数归一化：每个细胞的 counts 除以总计数并乘以 10000，再 log1p 变换
cat("\n==== 3. 标准化 ====\n")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", verbose = FALSE)

# 3.2 挑选高变基因（vst 法，3000 个）：这些基因最能反映细胞间差异，用于后续 PCA
obj <- FindVariableFeatures(obj, selection.method = "vst",
                            nfeatures = NFEATURES, verbose = FALSE)
cat(sprintf("高变基因数: %d\n", length(VariableFeatures(obj)))) # 打印 HVG 数

# 3.3 数据标准化（Z-score）：使基因表达均值为 0、方差为 1，PCA 才能正确计算
#     可选回归 percent.mt（REGRESS_MT=TRUE 时移除线粒体比例带来的变异，对齐 yijian.R）
if (REGRESS_MT) {                                  # 需要回归线粒体比例时
  cat("ScaleData 回归 percent.mt...\n")             # 打印提示
  obj <- ScaleData(obj, vars.to.regress = "percent.mt", verbose = FALSE) # 回归后缩放
} else {                                           # 不回归时
  obj <- ScaleData(obj, verbose = FALSE)           # 直接 Z-score 缩放
}

## ---- 4. PCA 降维与 PC 数选择 ----

# 通用 PDF 保存函数：把 ggplot 对象 p 输出到 OUT_DIR 下的 file，尺寸 w×h 英寸
save_pdf <- function(p, file, w, h) {              # 通用 PDF 保存函数（提前定义）
  pdf(file.path(OUT_DIR, file), width = w, height = h) # 打开 PDF 设备
  print(p)                                         # 打印图形
  dev.off()                                        # 关闭设备（必须，否则文件不完整）
  cat(sprintf("  已保存图表: %s\n", file))          # 提示保存成功
}

# 4.1 运行 PCA：先跑足 NPCS 个主成分（50），供后续 ElbowPlot/方差表判断保留多少个 PC
cat("\n==== 4. PCA 降维与 PC 数选择 ====\n")
obj <- RunPCA(obj, npcs = NPCS, verbose = FALSE)   # 运行 PCA（先不限定维度，跑足供判断）

# 4.2 输出 PC 方差表：每个 PC 的标准差、方差占比、累积方差解释率，
#     这是"选择 PC 数"的结果依据（配合 ElbowPlot 使用）
pca_stdev <- obj@reductions$pca@stdev              # 各 PC 的标准差（sqrt(特征值)）
var_ratio <- pca_stdev^2 / sum(pca_stdev^2)        # 各 PC 的方差解释比例
cumvar    <- cumsum(var_ratio)                     # 累积方差解释比例
pc_tbl <- data.frame(
  PC       = 1:NPCS,                               # PC 编号
  stdev    = pca_stdev,                            # 标准差
  var_ratio = var_ratio,                           # 单 PC 方差占比
  cumvar   = cumvar                                # 累积方差占比
)
write.csv(pc_tbl, file.path(OUT_DIR, "01_PC_variance_table.csv"), row.names = FALSE) # 写方差表
cat("  已保存: 01_PC_variance_table.csv（各 PC 方差/累积方差，供选择 PC 数）\n")       # 提示

# 4.3 绘制 Elbow 图：方差贡献拐点图，辅助判断保留多少个 PC（拐点之后的 PC 贡献趋于平缓）
p_elbow <- ElbowPlot(obj, ndims = NPCS) +          # Elbow 图
  geom_vline(xintercept = if (!is.null(DIMS)) max(DIMS) else NA, # 若手动指定则标出截断线
             linetype = "dashed", color = "red") + # 红色虚线标记选中的 PC 数
  theme_minimal() +                                # 简洁主题
  labs(title = "PCA Elbow Plot")                   # 图标题
save_pdf(p_elbow, "01_PCA_elbow.pdf", 8, 5)        # 保存 Elbow 图

# 4.3b 绘制 DimHeatmap（按 Lesson-4 教程规范）：
#       展示每个 PC 的 top 正/负相关基因在 top/bottom 500 细胞上的表达热图。
#       这是 PC 选择的另一维度诊断：每个 PC 是否捕获了真实的生物学变异
#       （看 top 基因功能），与 ElbowPlot 互补——Elbow 决定"砍到第几维"，
#       DimHeatmap 决定"前几个 PC 是否有效信号"。两图都属"计算 PC"阶段。
n_heat_pcs <- min(NPCS, 12)                          # 默认绘制前 12 个 PC 的热图（够判读且不过长）
heat_max   <- if (!is.null(DIMS)) max(DIMS) else n_heat_pcs
n_heat_pcs <- min(n_heat_pcs, heat_max)              # 与 ElbowPlot/DIMS 范围保持一致
cat(sprintf("绘制 PC_1 ~ PC_%d 的 DimHeatmap...\n", n_heat_pcs))
pdf(file.path(OUT_DIR, "01_PCA_dimheatmap.pdf"), width = 8.5, height = 6) # 多页 PDF，每 PC 一页
for (i in seq_len(n_heat_pcs)) {                     # 遍历 PC_1 ~ PC_n
  p_h <- DimHeatmap(                                 # 该 PC 的 dim heatmap
    object      = obj,                               # Seurat 对象
    dims        = i,                                 # 当前 PC
    cells       = 500,                               # 取 top/bottom 各 250 共 500 细胞
    balanced    = TRUE,                              # 按表达值正负平衡采样，避免一边倒
    nfeatures   = 30,                                # 取前 30 个高负荷基因
    reduction   = "pca",                             # 基于 PCA
    fast        = FALSE                              # fast=FALSE 显示完整基因标签（与截图一致）
  )
  print(p_h + ggtitle(paste0("PC_", i)))             # 打印到当前 PDF 设备，附 PC 编号标题
}
dev.off()                                            # 关闭 PDF 设备，文件落盘
cat(sprintf("  已保存: 01_PCA_dimheatmap.pdf（PC_1 ~ PC_%d，每个 PC 一页）\n", n_heat_pcs))

# 4.4 确定 DIMS（PC 数选择）：
#     遵循"参数选择先有结果"的规范——若 CONFIG 中 DIMS=NULL（默认），则根据刚算出的
#     方差表自动选择 PC 数，并打印选择依据；用户可查看 Elbow 图/方差表后手动改 DIMS。
if (is.null(DIMS)) {                               # 自动选择模式
  # 方法1：累积方差解释达到 CUMVAR_THRESHOLD（85%）所需的最小 PC 数
  k_cumvar <- which(cumvar >= CUMVAR_THRESHOLD)[1] # 首次达到阈值的位置
  if (is.na(k_cumvar)) k_cumvar <- NPCS            # 未达到则取全部
  # 方法2：第一个标准差 < STDEV_THRESHOLD（1）的 PC 之前的 PC 数（噪音截断）
  k_stdev <- which(pca_stdev < STDEV_THRESHOLD)[1] - 1 # 噪音 PC 的前一个
  if (is.na(k_stdev) || k_stdev < 1) k_stdev <- NPCS    # 无噪音 PC 则取全部
  # 取两种方法的较小者（保守选择显著 PC），并限制在 [MIN_PCS, MAX_PCS] 区间
  k <- max(MIN_PCS, min(k_cumvar, k_stdev, MAX_PCS))
  DIMS <- 1:k                                      # 生成 PC 范围
  cat(sprintf("自动选择 PC 数 = %d（累积方差法=%d，stdev<%.1f法=%d）\n",
              k, k_cumvar, STDEV_THRESHOLD, k_stdev)) # 打印选择依据
  cat(sprintf("DIMS = 1:%d；如想手动调整，请查看 01_PCA_elbow.pdf 后修改 CONFIG 中的 DIMS。\n", k))
} else {                                           # 手动指定模式
  cat(sprintf("使用手动指定 DIMS = %s\n", paste(range(DIMS), collapse = ":"))) # 打印范围
}
DIMS = 1:30
# 4.5 打印前 5 个主成分的方差解释率（快速概览）
cat("前 5 个主成分的方差解释率(%):\n")               # 打印提示
print(round(100 * head(pc_tbl$var_ratio, 5), 2))   # 前 5 PC 方差占比

## ---- 5. 去批次（Harmony 整合）----

# Harmony 以 sample 为批次变量，在 PCA 空间迭代校正样本间的技术差异，
# 使相同细胞类型跨样本对齐（同时保留生物学差异）。整合后的嵌入保存在
# reduction "harmony" 中，供下游 UMAP/聚类使用。
if (RUN_HARMONY) {                                 # 开启去批次时
  cat("\n==== 5. Harmony 去批次（批次变量: sample）====\n")
  harmony_args <- list(                            # 用列表统一管理参数（便于条件追加）
    object = obj,                                  # Seurat 对象
    group.by.vars = HARMONY_BATCH,                 # 指定批次变量（样本）
    reduction.use = "pca",                         # 基于 PCA 结果整合
    reduction.save = "harmony",                    # 保存为新降维 "harmony"
    theta = 5,           # 默认 2，增加到 3-5
    max.iter.harmony = 50,# 增加迭代次数
    sigma = 0.1,                  # 聚类带宽（默认0.1，可尝试0.05-0.2）
    verbose = F                                # 静默模式
  )
  if (!is.null(HARMONY_THETA)) {                   # 手动指定了校正强度 theta 时
    harmony_args$theta <- HARMONY_THETA            # 追加 theta 参数（如 5 增强去批次）
  }
  if (PLOT_HARMONY_CONVERGENCE) {                  # 需要收敛诊断图时
    harmony_args$plot_convergence <- TRUE          # 让 RunHarmony 返回收敛过程数据
  }
  # 执行 Harmony。注意：plot_convergence=TRUE 时 RunHarmony 返回一个列表
  # list(obj=整合后对象, final=Harmony 模型)，而非直接返回对象，需兼容两种返回类型
  harmony_res <- do.call(RunHarmony, harmony_args) # 执行 Harmony（参数完整传入）
  if (is.list(harmony_res) && !inherits(harmony_res, "Seurat")) { # 返回了列表（带收敛数据）
    obj <- harmony_res$obj                         # 取整合后的 Seurat 对象
    harmony_final <- harmony_res$final             # 取 Harmony 模型（含收敛过程，供诊断图使用）
  } else {                                         # 返回的是对象本身（未开启收敛）
    obj <- harmony_res                             # 直接用
    harmony_final <- NULL                          # 无收敛数据
  }
  reduction_use <- "harmony"                       # 下游使用 harmony 嵌入
} else {                                           # 关闭去批次时
  cat("\n==== 5. 跳过 Harmony（RUN_HARMONY=FALSE），使用 PCA 嵌入 ====\n")
  reduction_use <- "pca"                           # 下游使用 PCA 嵌入
  harmony_final <- NULL                            # 无 Harmony 模型（收敛图自动跳过）
}
cat(sprintf("下游降维将使用: %s（共 %d 维中的前 %d 维）\n",
            reduction_use, NPCS, max(DIMS)))       # 打印实际使用的降维与维度数

## ---- 6. UMAP / tSNE 可视化降维 ----

# 6.1 UMAP：非线性降维，用于二维可视化（基于 harmony/pca 嵌入的前 DIMS 个维度）。
#     参数按 Seurat 最佳实践扩展：uwot 实现 + cosine 距离 + 邻域/最小距离等可调
cat("\n==== 6. UMAP / tSNE ====\n")
cat("运行 UMAP（约需数分钟）...\n")
obj <- RunUMAP(obj,
  reduction = reduction_use,                       # 输入降维（harmony/pca）
  dims = DIMS,                                     # 使用的主成分范围（已按 Elbow 确定）
  umap.method = UMAP_METHOD,                       # UMAP 实现（R 原生 uwot）
  metric = UMAP_METRIC,                            # 距离度量（cosine 适合单细胞）
  n.neighbors = UMAP_NEIGHBORS,                    # 局部邻域大小
  min.dist = UMAP_MIN_DIST,                        # 嵌入点最小距离（控制聚集程度）
  spread = UMAP_SPREAD,                            # 嵌入有效尺度
  seed.use = SEED,                                 # 随机种子（可复现）
  verbose = FALSE)                                 # 静默模式
cat("UMAP 完成。\n")

# 6.2 tSNE：另一种非线性降维（可选，较慢；RUN_TSNE=FALSE 时跳过）
if (RUN_TSNE) {                                    # 开启 tSNE 时
  cat("运行 tSNE（较慢，约需 5~15 分钟）...\n")      # 提示耗时
  obj <- RunTSNE(obj,
    reduction = reduction_use,                     # 输入降维（harmony/pca）
    dims = DIMS,                                   # 使用的主成分范围
    perplexity = TSNE_PERPLEXITY,                  # 困惑度（平衡局部与全局结构）
    check_duplicates = FALSE,                      # 跳过重复点检查（单细胞无重复条码）
    verbose = FALSE)                               # 静默模式
  cat("tSNE 完成。\n")
}

## ---- 7. 聚类（FindNeighbors + FindClusters，多分辨率）----

# 7.1 构建 KNN 图与 SNN 图（基于 harmony/pca 嵌入，前 DIMS 维），聚类的基础
cat("\n==== 7. 聚类（多分辨率）====\n")
obj <- FindNeighbors(obj, reduction = reduction_use, dims = DIMS, verbose = FALSE)

# 7.2 依次跑多个分辨率：每个分辨率重跑一次 FindClusters（复用同一个图），
#     结果分别存为 cluster_0.2 / cluster_0.3 等列，供后续自行选择
for (res in RESOLUTIONS) {                         # 遍历分辨率向量
  cat(sprintf("  聚类分辨率 res=%s ...\n", res))    # 提示当前分辨率
  obj <- FindClusters(obj,
    resolution = res,                              # 分辨率（越大簇越多）
    algorithm = CLUSTER_ALGORITHM,                 # 聚类算法（1=Louvain）
    method = CLUSTER_METHOD,                       # 图划分方法（igraph）
    group.singletons = GROUP_SINGLETONS,           # 孤立单细胞并入最近簇
    verbose = FALSE)                               # 静默模式
  obj[[paste0("cluster_", res)]] <- as.character(obj$seurat_clusters) # 保存该分辨率结果
  cat(sprintf("    得到 %d 个 cluster\n", length(unique(obj$seurat_clusters)))) # 打印簇数
}
cat("可视化顺序：先多分辨率聚类拼图（辅助选 res），后默认分辨率下的常规图。\n")

# ==== 9.1 多分辨率 UMAP 拼图（先出，便于选 res） ====
# 将各分辨率的聚类结果并排展示，直观对比不同分辨率下的簇划分，
# 辅助选择合适的分辨率。后续 9.2~9.7 全部基于默认分辨率（DEFAULT_RESOLUTION）。
if (length(RESOLUTIONS) > 1) {                     # 多于 1 个分辨率时
  n_col <- min(3, length(RESOLUTIONS))             # 每行最多 3 幅
  plist <- lapply(RESOLUTIONS, function(r) {       # 逐分辨率生成 UMAP 图
    DimPlot(obj, reduction = "umap",               # 用 UMAP 嵌入
            group.by = paste0("cluster_", r),      # 按该分辨率聚类列着色
            label = TRUE, repel = TRUE) +          # 标注簇号
      NoLegend() +                                 # 隐藏图例（拼图省空间）
      ggtitle(paste0("res = ", r))                 # 标题注明分辨率
  })
  p_multi <- patchwork::wrap_plots(plist, ncol = n_col) # 拼接为多面板图
  save_pdf(p_multi, "02_UMAP_multi_resolution.pdf",
           7 * n_col, 6 * ceiling(length(RESOLUTIONS) / n_col)) # 保存拼图
}


## ---- 8. 聚类统计表输出 ----
# 7.3 指定默认分辨率作为主 cluster 列（默认 0.5），并设为 Seurat 的当前 idents
DEFAULT_RESOLUTION <- 0.4             # 默认分辨率：其聚类结果写入主 cluster 列
obj$cluster <- as.character(obj@meta.data[[paste0("cluster_", DEFAULT_RESOLUTION)]]) # 主 cluster 列
Idents(obj) <- obj$cluster                         # 设为当前 idents（后续绘图按此着色）

# 7.4 打印默认分辨率的聚类统计
cat(sprintf("\n默认分辨率（res=%s）聚类结果：\n", DEFAULT_RESOLUTION)) # 提示
print(table(obj$cluster))                           # 各 cluster 细胞数

# 8.1 默认分辨率的聚类汇总：cluster / 细胞数 / 占比
clust_summary <- obj@meta.data %>%
  dplyr::group_by(cluster) %>%                     # 按主 cluster 分组
  dplyr::summarise(
    n_cells = dplyr::n(),                          # 细胞数
    pct     = 100 * dplyr::n() / ncol(obj),        # 占比(%)
    .groups = "drop"
  ) %>%
  dplyr::arrange(cluster)                          # 按 cluster 排序
write.csv(clust_summary, file.path(OUT_DIR, "01_clustering_summary.csv"), row.names = FALSE) # 写汇总
cat("  已保存: 01_clustering_summary.csv\n")        # 提示

# 8.2 默认分辨率的 cluster × 样本交叉表（含分组信息），便于检查各簇的样本构成

clust_by_sample <- obj@meta.data %>%
  dplyr::group_by(cluster, sample, group, sample_id) %>% # 按簇+样本+分组+编号分组
  dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") # 细胞数
write.csv(clust_by_sample, file.path(OUT_DIR, "02_cluster_by_sample.csv"), row.names = FALSE) # 写交叉表
cat("  已保存: 02_cluster_by_sample.csv\n")        # 提示

# 8.3 各分辨率的簇数汇总（长格式，方便对比选择分辨率）
res_summary <- data.frame(
  resolution = RESOLUTIONS,                        # 分辨率
  n_clusters = sapply(RESOLUTIONS, function(r) {   # 对应簇数
    length(unique(obj@meta.data[[paste0("cluster_", r)]]))
  }),
  row.names = NULL
)
write.csv(res_summary, file.path(OUT_DIR, "03_resolution_summary.csv"), row.names = FALSE) # 写汇总
cat("  已保存: 03_resolution_summary.csv\n")        # 提示

## ---- 9. 可视化 ----

cat("\n==== 9. 可视化 ====\n")
# ==== 9.2 UMAP 按样本着色（默认分辨率下）：检查 Harmony 去批次效果 ====
p_smp <- DimPlot(obj, reduction = "umap", group.by = "sample",
                 label = FALSE) +                  # 按样本着色
  theme_minimal() +                                # 简洁主题
  labs(title = "UMAP by sample (after Harmony)")   # 图标题
save_pdf(p_smp, "03_UMAP_by_sample.pdf", 9, 7)     # 保存

# ==== 9.3 UMAP 按分组着色（默认分辨率下）：观察 ypN0 / ypN+ 分布差异 ====
p_grp <- DimPlot(obj, reduction = "umap", group.by = "group",
                 cols = c(ypN0 = "#4DBBD5", `ypN+` = "#E64B35")) + # ypN0=蓝, ypN+=红
  theme_minimal() +                                # 简洁主题
  labs(title = "UMAP by group")                    # 图标题
save_pdf(p_grp, "04_UMAP_by_group.pdf", 8, 6)      # 保存

# ==== 9.4 UMAP 按主 cluster 着色（默认分辨率下）：查看聚类结构 ====
p_clu <- DimPlot(obj, reduction = "umap", group.by = "cluster",
                 label = TRUE, repel = TRUE) +     # 标注簇号
  theme_minimal() +                                # 简洁主题
  labs(title = paste0("UMAP by cluster (res=", DEFAULT_RESOLUTION, ")")) # 图标题
save_pdf(p_clu, "05_UMAP_by_cluster.pdf", 9, 7)    # 保存

# ==== 9.5 UMAP 按样本分面（默认分辨率下）：检查各样本去批次效果 ====
p_split <- DimPlot(obj, reduction = "umap", group.by = "sample",
                   split.by = "sample", ncol = 6) + # 每样本一面板
  theme_minimal() +                                # 简洁主题
  theme(legend.position = "none",                  # 隐藏图例
        axis.text = element_blank(), axis.ticks = element_blank()) + # 简化坐标轴
  labs(title = "UMAP by sample (split)")           # 图标题
save_pdf(p_split, "06_UMAP_by_sample_split.pdf", 18, 12) # 保存（大画布容纳 23 面板）

# ==== 9.6 主 cluster 细胞数柱状图（默认分辨率下） ====
p_bar <- ggplot(clust_summary, aes(x = cluster, y = n_cells)) + # 柱状图
  geom_col(fill = "#3C5488") +                     # 深蓝色柱子
  geom_text(aes(label = n_cells), vjust = -0.3, size = 3) + # 柱顶标注细胞数
  theme_minimal() +                                # 简洁主题
  labs(title = paste0("Cell count per cluster (res=", DEFAULT_RESOLUTION, ")"), # 图标题
       x = "Cluster", y = "Cell count")            # 坐标轴标签
save_pdf(p_bar, "07_cluster_cell_count.pdf", 9, 5) # 保存

# ==== 9.7 tSNE 按主 cluster 着色（默认分辨率下，可选） ====
if (RUN_TSNE && "tsne" %in% names(obj@reductions)) { # tSNE 存在时
  p_tsne <- DimPlot(obj, reduction = "tsne", group.by = "cluster",
                    label = TRUE, repel = TRUE) +  # 标注簇号
    theme_minimal() +                              # 简洁主题
    labs(title = paste0("tSNE by cluster (res=", DEFAULT_RESOLUTION, ")")) # 图标题
  save_pdf(p_tsne, "08_TSNE_by_cluster.pdf", 9, 7) # 保存
}

# ==== 9.8 Harmony 去批次效果对比图（核心诊断图，默认分辨率下） ====
# 左 = 未整合（基于 PCA 的 UMAP），右 = 整合后（基于 Harmony 的 UMAP），均按样本着色。
# 判读标准：左图样本各自成团（批次效应明显）→ 右图样本混合均匀即去批次有效；
# 若右图仍存在整块样本聚团 → 去批次不足，可调大 HARMONY_THETA（如 5）后重跑；
# 右图出现单一类型跨样本对齐但与其他类型分离，属正常（生物学差异保留）。
if (RUN_HARMONY && "harmony" %in% names(obj@reductions)) { # 开启了 Harmony 时
  cat("  绘制 Harmony 去批次前后对比图（需额外运行一次 PCA-UMAP，约数分钟）...\n")
  # 基于 PCA 再跑一个独立 UMAP，存为 umap.pca（不覆盖整合后的 umap 降维）
  obj <- RunUMAP(obj,
    reduction = "pca", dims = DIMS,                # 输入未整合的 PCA 嵌入
    reduction.name = "umap.pca", reduction.key = "PCUMAP_", # 独立命名保存，避免冲突
    umap.method = UMAP_METHOD, metric = UMAP_METRIC, # 与主 UMAP 相同的参数（保证可比）
    n.neighbors = UMAP_NEIGHBORS, min.dist = UMAP_MIN_DIST,
    spread = UMAP_SPREAD, seed.use = SEED, verbose = FALSE) # 完整参数
  # 左图：整合前（PCA-UMAP）按样本着色 —— 预期各样本各自成团（批次效应）
  p_before <- DimPlot(obj, reduction = "umap.pca", group.by = "sample",
                      label = FALSE) +             # 按样本着色
    theme_minimal() +                              # 简洁主题
    labs(title = "Before Harmony (PCA-UMAP)")      # 标题注明未整合
  # 右图：整合后（Harmony-UMAP）按样本着色 —— 预期样本混合均匀（去批次效果）
  p_after <- DimPlot(obj, reduction = "umap", group.by = "sample",
                     label = FALSE) +              # 按样本着色
    theme_minimal() +                              # 简洁主题
    labs(title = "After Harmony (Harmony-UMAP)")   # 标题注明已整合
  p_ba <- patchwork::wrap_plots(p_before, p_after, ncol = 2) # 左右并排拼图
  save_pdf(p_ba, "09_Harmony_before_after_UMAP.pdf", 16, 7) # 保存对比图
}

# ==== 9.9 Harmony 收敛诊断图（可选） ====
# 展示各轮迭代批次校正目标（theta 散度）的下降情况。
# 曲线随迭代趋于平稳 → Harmony 已收敛（默认迭代轮数足够）；
# 若尚未平稳，说明需要更多迭代，可在 RunHarmony 中增大 max.iter.harmony。
# 注：harmony 不同版本的收敛数据槽结构有差异，绘制失败会自动跳过，不影响主流程。
if (PLOT_HARMONY_CONVERGENCE && !is.null(harmony_final)) { # 有收敛数据时
  p_conv <- tryCatch(                              # 版本差异兜底：失败返回 NULL
    harmony::plot_convergence(harmony_final),      # harmony 包自带的收敛图函数
    error = function(e) NULL)                      # 出错则跳过
  if (!is.null(p_conv)) {                          # 绘制成功时
    p_conv <- p_conv + theme_minimal() +           # 简洁主题
      labs(title = "Harmony convergence")          # 图标题
    save_pdf(p_conv, "10_Harmony_convergence.pdf", 7, 5) # 保存收敛图
  } else {                                         # 绘制失败时
    cat("  注意: Harmony 收敛图绘制失败（harmony 版本差异），已跳过，不影响主流程。\n")
  }
}

## ---- 10. 保存整合对象与完成 ----

# 保存整合聚类后的 Seurat 对象（含 counts/data/scale.data、pca/harmony/umap/tsne 降维、
# 多分辨率聚类列 cluster_0.2/0.3/0.4/0.5、主 cluster 列、以及全部 meta.data 信息，
# 供后续细胞注释/差异分析直接使用）
saveRDS(obj, file.path(OUT_DIR, "04_CD45_integrated.rds"))
cat(sprintf("  已保存: %s\n", file.path(OUT_DIR, "04_CD45_integrated.rds")))

# 打印最终结果摘要
cat("\n==== 04 整合聚类完成 ====\n")
cat(sprintf("细胞数: %d\n", ncol(obj)))                        # 总细胞数
cat(sprintf("默认分辨率 res=%s 的簇数: %d\n",                  # 簇数
            DEFAULT_RESOLUTION, length(unique(obj$cluster))))
cat(sprintf("本次使用的 PC 数（DIMS）: %d\n", max(DIMS)))      # 实际使用的 PC 数
cat("\n输出文件清单（位于 output/ 目录）：\n")
cat("  04_CD45_integrated.rds           - 整合聚类后的 Seurat 对象\n")
cat("  01_PC_variance_table.csv         - 各 PC 方差/累积方差表（PC 选择依据）\n")
cat("  01_clustering_summary.csv        - 主 cluster 细胞数与占比\n")
cat("  02_cluster_by_sample.csv         - cluster × 样本交叉表\n")
cat("  03_resolution_summary.csv        - 各分辨率簇数汇总\n")
cat("  01_PCA_elbow.pdf                 - PCA Elbow 图（PC 选择依据）\n")
cat("  01_PCA_dimheatmap.pdf            - PC_1~PC_n 的 DimHeatmap（每个 PC 一页，含 top 正负基因）\n")
cat("  02_UMAP_multi_resolution.pdf     - 多分辨率 UMAP 拼图（先出，辅助选 res）\n")
cat("  03_UMAP_by_sample.pdf            - UMAP 按样本着色（去批次效果，默认分辨率下）\n")
cat("  04_UMAP_by_group.pdf             - UMAP 按分组着色（默认分辨率下）\n")
cat("  05_UMAP_by_cluster.pdf           - UMAP 按主 cluster 着色（默认分辨率下）\n")
cat("  06_UMAP_by_sample_split.pdf      - UMAP 按样本分面（默认分辨率下）\n")
cat("  07_cluster_cell_count.pdf        - 各 cluster 细胞数柱状图（默认分辨率下）\n")
cat("  08_TSNE_by_cluster.pdf           - tSNE 按主 cluster 着色（可选，默认分辨率下）\n")
cat("  09_Harmony_before_after_UMAP.pdf - Harmony 去批次前后 UMAP 对比（核心诊断图）\n")
cat("  10_Harmony_convergence.pdf       - Harmony 收敛诊断图（可选，版本兼容）\n")
cat("\n04_integration_clustering 完成。\n")

