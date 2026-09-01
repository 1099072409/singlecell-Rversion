# =====================================================================
# 02_quality_control.R
# 单细胞数据分析流程 - 第二步（质量控制 QC）
#
# 功能：
#   1) 读取第一步合并好的 Seurat 对象（01_seurat_combined.rds）
#   2) 重新读取 data/样本分组信息.xlsx，与 rds 中的 meta.data 交叉校验，
#      确保"分组 Group + 样本编号 SampleID"的对应关系不丢失、不串位
#   3) 计算质控指标：线粒体基因比例 percent.mt、核糖体基因比例 percent.ribo、
#      血红蛋白比例 percent.hb（连同 Seurat 自带的 nCount_RNA / nFeature_RNA）
#   4) 过滤低质量细胞（两段式，fixed 模式）：
#        阶段1 基础过滤（先砍下界）：nFeature>200 且 nCount>600
#        阶段2 双细胞检测：scDblFinder（可关闭；按样本汇总输出 doublet_summary_by_sample.csv）
#        阶段3 精细过滤（砍上界）：singlet 且 nFeature<6000、nCount<30000、percent.mt<25
#        （"mad"：参考 yijian.R 的按样本 3-MAD 自适应离群检测，保留可选）
#   5) 输出：过滤后 Seurat 对象（02_seurat_qc.rds）、按样本汇总表
#      （01_QC_summary_by_sample.csv / cell_count_by_sample.csv / doublet_summary_by_sample.csv）、
#      QC 图（过滤前后小提琴图、散点图、细胞数柱状图，PDF 格式）
#
# 用法：
#   在 RStudio 中打开本文件并点击 Source，或在命令行执行：
#     "C:/Program Files/R/R-4.4.3/bin/Rscript.exe" 02_quality_control.R
#   说明：所有路径按脚本自身所在目录自动推算，可从任意工作目录运行。
#
# 常用修改点：见下方 ## ---- 0. 配置 ---- 中的 CONFIG 参数（阈值、模式等）。
# 参考：Rversion/yijian.R（MAD 质控范式）、Rversion/01_initiation/01_initiation.R（脚本风格）
# =====================================================================

## ---- 0. 配置（按需修改） ----

# 自动定位脚本所在目录（兼容 Rscript 与 RStudio 两种运行方式）：
#   - RStudio 中 Source 时，通过 sys.frames()[[1]]$ofile 拿到脚本路径
#   - 命令行 Rscript 时，通过 --file= 参数拿到脚本路径
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
INPUT_RDS <- file.path(script_dir, "..", "01_initiation", "output", "01_seurat_combined.rds") # 第一步合并好的 rds
XLSX_PATH <- file.path(script_dir, "..", "..", "data", "样本分组信息.xlsx")                     # 分组表（含 Group/Sample/样本编号）
OUT_DIR   <- file.path(script_dir, "output")                                                  # 本步输出目录（自动创建）

# ---- 质控指标（人类基因命名规则，切换物种时修改） ----
MT_PATTERN   <- "^MT-"        # 线粒体基因前缀（human）
RIBO_PATTERN <- "^RP[SL]"     # 核糖体蛋白基因前缀（RPS/RPL）
HB_PATTERN   <- "^HB[ABDEGQZ]"# 血红蛋白基因前缀

# ---- 质控模式与阈值 ----
QC_MODE          <- "fixed"   # 质控模式："fixed"=固定阈值（两段式）；"mad"=按样本 3-MAD 自适应（参考 yijian.R）
# 提示：阈值设定遵循"参数选择先有结果"的规范——首次运行时先用默认阈值跑一遍，
#       查看输出的 output/01_QC_violin_before_filter.pdf 与 output/00_QC_quantiles_by_sample.csv
#       （各样本 QC 指标 5%/50%/95% 分位数），再回来调整下面的阈值。
# fixed 模式参数（QC_MODE="fixed" 时生效）——两段式过滤：
FIXED_MIN_FEAT   <- 200       # 阶段1 下界：细胞最少检测到的基因数
FIXED_MAX_FEAT   <- 6000      # 阶段3 上界：细胞最多检测到的基因数（剔除极可能为双细胞/高复杂度细胞）
FIXED_MIN_COUNTS <- 600       # 阶段1 下界：细胞最少总 UMI 计数
FIXED_MAX_COUNTS <- 30000     # 阶段3 上界：细胞最多总 UMI 计数
FIXED_MAX_PCT_MT <- 20        # 阶段3 上界：细胞线粒体基因比例上限（%）
# mad 模式参数（QC_MODE="mad" 时生效）：
MAD_NMADS        <- 3         # MAD 离群判定倍数（中位数 ± nmads × MAD 之外视为离群）
MAD_MIN_FEAT     <- 200       # mad 模式下 nFeature_RNA 的硬性下限（与 yijian.R 一致）

# ---- 双细胞检测 ----
RUN_DOUBLET      <- TRUE      # 是否运行 scDblFinder 双细胞检测（需要安装 Bioconductor 包）
DOUBLET_THREADS  <- 2         # 双细胞检测线程数（Windows 下自动降级为单线程，见第 6 节说明）

# ---- 其他 ----
STOP_ON_MISMATCH <- FALSE     # 元数据校验不一致时：TRUE=报错终止；FALSE=警告后继续
SEED             <- 123       # 随机种子（保证 scDblFinder 结果可复现）

## ---- 1. 依赖检查与加载 ----

# 通用依赖安装函数：ensure_pkg("包名", bioc=是否 Bioconductor 包)
#   - 已安装则直接返回 TRUE
#   - 未安装则自动安装（CRAN 走 install.packages，Bioc 走 BiocManager::install）
ensure_pkg <- function(pkg, bioc = FALSE) {
  if (requireNamespace(pkg, quietly = TRUE)) return(TRUE)             # 已装则跳过
  cat(sprintf("  正在安装依赖包: %s ...\n", pkg))                    # 提示开始安装
  if (bioc) {                                                         # Bioconductor 包分支
    if (!requireNamespace("BiocManager", quietly = TRUE)) {           # BiocManager 本身未装则先装它
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, update = FALSE, ask = FALSE)            # 通过 BiocManager 安装
  } else {                                                            # CRAN 包分支
    install.packages(pkg, repos = "https://cloud.r-project.org")      # 通过官方镜像安装
  }
  requireNamespace(pkg, quietly = TRUE)                               # 返回是否安装成功
}

# 依次确保需要的包都可用（Seurat/readxl/dplyr/ggplot2 已装；patchwork 为绘图组合依赖；其余为双细胞检测所需）
for (p in c("Seurat", "readxl", "dplyr", "ggplot2", "patchwork"))  ensure_pkg(p, bioc = FALSE)
for (p in c("BiocParallel", "SingleCellExperiment", "SummarizedExperiment", "scDblFinder")) ensure_pkg(p, bioc = TRUE)

# 加载核心包（Bioc 包通过 :: 显式调用，避免命名空间冲突）
library(Seurat)    # 单细胞分析主包
library(readxl)    # 读取 xlsx 分组表
library(dplyr)     # 数据整理（group_by/summarise 等）
library(ggplot2)   # 绘图

# 创建输出目录（已存在则不报错）
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 打印运行环境信息，便于排查问题
cat(sprintf("脚本目录: %s\n", script_dir))        # 脚本所在目录
cat(sprintf("工作目录: %s\n", getwd()))           # 当前工作目录
cat(sprintf("输出目录: %s\n", normalizePath(OUT_DIR))) # 输出目录
cat(sprintf("质控模式: %s\n", QC_MODE))           # 当前质控模式

# Seurat v5 并行/全局对象大小限制（为后续标准化、Harmony 等大对象操作预留空间）
options(future.globals.maxSize = 8 * 1024^3)

# 检查可用物理内存并给出提示（Windows 下通过 wmic 查询；查询失败时自动跳过，不影响运行）。
# 合并对象约需 5~8GB 内存，若可用内存过少，提前提醒用户关闭其他程序后再运行。
tryCatch({
  w  <- system("wmic OS get FreePhysicalMemory /value", intern = TRUE) # 查询可用内存(KB)
  fr <- as.numeric(sub(".*=", "", w[grepl("FreePhysicalMemory", w)]))  # 提取数值
  if (length(fr) == 1 && !is.na(fr)) {                                 # 查询成功时
    free_gb <- fr / 1024^2                                             # 转 GB
    cat(sprintf("当前可用物理内存: %.1f GB\n", free_gb))                # 打印可用内存
    if (free_gb < 8) {                                                 # 低于 8GB 时提示
      cat("提示: 可用内存偏少，本步大对象约需 5~8GB，建议关闭其他程序后再运行。\n")
    }
  }
}, error = function(e) invisible(NULL))                                 # 查询失败静默跳过

## ---- 2. 读取分组表（data/样本分组信息.xlsx） ----

# 检查分组表是否存在
if (!file.exists(XLSX_PATH)) stop("找不到分组表: ", XLSX_PATH)
meta <- read_xlsx(XLSX_PATH)                       # 读取 xlsx 全部列

# 兼容中/英文列名的自动匹配函数：从候选名中返回第一个命中的真实列名
find_col <- function(candidates, df) {
  for (c in candidates) {                          # 逐个尝试候选名
    idx <- grep(c, names(df), ignore.case = TRUE)  # 在表头中模糊匹配（忽略大小写）
    if (length(idx) >= 1) return(names(df)[idx[1]]) # 命中即返回真实列名
  }
  return(NA_character_)                            # 全部未命中返回 NA
}

# 分别匹配三列：样本名 / 分组 / 样本编号（支持中英文表头）
sample_col <- find_col(c("Sample", "样本", "样本名"), meta)      # 样本名列（对应样本文件夹名）
group_col  <- find_col(c("Group", "分组", "组别"), meta)         # 分组列（ypN0/ypN+）
id_col     <- find_col(c("样本编号", "编号", "SampleID", "ID"), meta) # 样本编号列（ypN01...）
cat(sprintf("分组表列映射 -> Sample: %s | Group: %s | 编号: %s\n",
            sample_col, group_col, id_col))        # 打印列映射，便于确认

# 只保留三列并统一重命名为 Sample / Group / SampleID
meta <- meta[, c(sample_col, group_col, id_col)]   # 按匹配到的列名取子集
colnames(meta) <- c("Sample", "Group", "SampleID") # 统一英文列名
meta$Sample   <- as.character(meta$Sample)         # 转字符型（样本名）
meta$Group    <- as.character(meta$Group)          # 转字符型（分组）
meta$SampleID <- as.character(meta$SampleID)       # 转字符型（样本编号）
cat(sprintf("分组表共 %d 个样本。\n", nrow(meta)))  # 打印样本数（应为 23）

## ---- 3. 读入合并对象并校验元数据 ----

# 3.1 读取第一步保存的 Seurat 对象（538MB 左右，读取需 1~3 分钟）
if (!file.exists(INPUT_RDS)) stop("找不到输入 rds: ", INPUT_RDS)
cat("正在读取 01_seurat_combined.rds（约需 1~3 分钟）...\n")
obj <- readRDS(INPUT_RDS)                          # 读入 Seurat 对象
cat(sprintf("读入完成：%d 个基因 x %d 个细胞\n", nrow(obj), ncol(obj))) # 打印维度

# 3.2 确保默认 assay 为 RNA（后续百分比计算都基于它）
if (DefaultAssay(obj) != "RNA") {                  # 若默认 assay 不是 RNA
  warning("默认 assay 不是 RNA，已切换为 RNA")      # 提示后切换
  DefaultAssay(obj) <- "RNA"
}

# 3.3 统一 Layer 结构（Seurat v5 关键步骤）：
#     01_initiation.R 用 merge() 合并 23 个样本，v5 的 merge 会把每个样本的
#     计数矩阵拆成多个 layer（counts.1 / counts.2 / ...），必须 JoinLayers
#     合并为单一的 counts 层，后续 PercentageFeatureSet / LayerData 才能正确工作。
cat(sprintf("当前 RNA layer: %s\n", paste(Layers(obj, assay = "RNA"), collapse = ", ")))
if (length(Layers(obj, assay = "RNA")) > 1) {      # 若存在多个 layer
  cat("检测到多个 layer（merge 产物），执行 JoinLayers 合并...\n")
  obj <- JoinLayers(obj, assay = "RNA")            # 合并所有 layer 为单一 counts
}
cat(sprintf("合并后 RNA layer: %s\n", paste(Layers(obj, assay = "RNA"), collapse = ", ")))

# 3.4 统一样本标识列：obj$sample = 样本名（与 yijian.R 的 sample 语义一致，供按样本分组统计/MAD/双细胞使用）
obj$sample <- as.character(obj$orig.ident)         # orig.ident 即样本文件夹名（如 HRR1430840）

# 3.5 元数据一致性校验函数：对比 rds 中的 meta.data 与 xlsx 分组表
check_meta_consistency <- function(obj, meta) {
  cat("\n==== 3.5 元数据一致性校验（rds meta.data vs 样本分组信息.xlsx）====\n")
  md <- obj@meta.data                              # 取 Seurat 对象的 meta.data

  # 按样本（orig.ident）汇总 rds 侧的分组与编号（取每样本唯一值）
  rds_tbl <- md %>%
    dplyr::group_by(orig.ident) %>%                # 按样本名分组
    dplyr::summarise(
      rds_group     = paste(unique(stats::na.omit(group)), collapse = "/"),     # rds 中该样本的分组
      rds_sample_id = paste(unique(stats::na.omit(sample_id)), collapse = "/"), # rds 中该样本的编号
      n_cells       = dplyr::n(),                  # 该样本的细胞数
      .groups       = "drop"                       # 取消分组
    )
  names(rds_tbl)[names(rds_tbl) == "orig.ident"] <- "Sample"  # 统一列名为 Sample

  # 与 xlsx 表按 Sample 左连接，得到逐样本对照表
  cmp <- dplyr::left_join(rds_tbl, meta, by = "Sample")
  cat("逐样本对照表（列为空表示仅另一侧存在）：\n")
  print(as.data.frame(cmp))                        # 打印对照表

  # 检测各类不一致
  only_rds       <- setdiff(rds_tbl$Sample, meta$Sample)     # 只出现在 rds 中的样本
  only_xlsx      <- setdiff(meta$Sample, rds_tbl$Sample)     # 只出现在 xlsx 中的样本
  na_cells       <- sum(is.na(md$group) | is.na(md$sample_id)) # 分组或编号为 NA 的细胞数
  mismatch_grp   <- sum(!is.na(cmp$Group) & !is.na(cmp$rds_group) &
                          as.character(cmp$Group) != as.character(cmp$rds_group))    # 分组不一致的样本数
  mismatch_id    <- sum(!is.na(cmp$SampleID) & !is.na(cmp$rds_sample_id) &
                          as.character(cmp$SampleID) != as.character(cmp$rds_sample_id)) # 编号不一致的样本数

  # 汇总判断：有任何不一致时按 STOP_ON_MISMATCH 决定终止或警告
  if (length(only_rds) > 0 || length(only_xlsx) > 0 || na_cells > 0 ||
      mismatch_grp > 0 || mismatch_id > 0) {
    msg <- sprintf(paste0("元数据不一致：仅rds有样本=[%s]；仅xlsx有样本=[%s]；",
                          "分组/编号为NA的细胞=%d；分组不一致样本=%d；编号不一致样本=%d"),
                   paste(only_rds, collapse = ","), paste(only_xlsx, collapse = ","),
                   na_cells, mismatch_grp, mismatch_id)
    if (STOP_ON_MISMATCH) stop(msg) else warning(msg)  # 按配置终止或警告
  } else {
    cat("校验通过：rds 与 xlsx 的样本集合、分组、样本编号完全一致。\n")
  }
  invisible(cmp)                                   # 返回对照表（不打印）
}

# 执行元数据一致性校验（保证后续过滤与分组统计建立在正确的分组/编号对应关系上）
meta_check <- check_meta_consistency(obj, meta)

## ---- 4. 计算质控指标 ----

# 4.1 线粒体基因比例 percent.mt（线粒体基因高比例是细胞受损/死亡的标志，需过滤）
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = MT_PATTERN)

# 4.2 核糖体蛋白基因比例 percent.ribo（参考 yijian.R，辅助判断细胞状态）
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = RIBO_PATTERN)

# 4.3 血红蛋白基因比例 percent.hb（红细胞污染标志，参考 yijian.R）
obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = HB_PATTERN)

# 打印各指标的分布概况，快速了解数据质量
cat("\n==== QC 指标分布（过滤前）====\n")
print(summary(obj@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.hb")]))

# 4.4 按样本输出各 QC 指标的分位数表（供确定过滤阈值参考）：
#     遵循"参数选择先有结果"的规范——先查看本表与 01_QC_violin_before_filter.pdf
#     中各样本的指标分布（5%/50%/95% 分位数），再回头调整第 0 节 CONFIG 中的
#     FIXED_MIN_FEAT / FIXED_MAX_FEAT / FIXED_MIN_COUNTS / FIXED_MAX_PCT_MT 等阈值。
qc_q <- obj@meta.data %>%
  dplyr::group_by(sample, group, sample_id) %>%              # 按样本/分组/编号分组
  dplyr::summarise(
    nCount_q05      = stats::quantile(nCount_RNA, 0.05),     # nCount 5% 分位数
    nCount_median   = stats::median(nCount_RNA),             # nCount 中位数
    nCount_q95      = stats::quantile(nCount_RNA, 0.95),     # nCount 95% 分位数
    nFeature_q05    = stats::quantile(nFeature_RNA, 0.05),   # nFeature 5% 分位数
    nFeature_median = stats::median(nFeature_RNA),           # nFeature 中位数
    nFeature_q95    = stats::quantile(nFeature_RNA, 0.95),   # nFeature 95% 分位数
    mt_q05          = stats::quantile(percent.mt, 0.05),     # percent.mt 5% 分位数
    mt_median       = stats::median(percent.mt),             # percent.mt 中位数
    mt_q95          = stats::quantile(percent.mt, 0.95),     # percent.mt 95% 分位数
    .groups         = "drop"
  )
write.csv(qc_q, file.path(OUT_DIR, "00_QC_quantiles_by_sample.csv"), row.names = FALSE) # 写分位数表
cat("  已保存: 00_QC_quantiles_by_sample.csv（各样本 QC 指标分位数，供设定过滤阈值参考）\n") # 提示

## ---- 5. 基础过滤（先砍下界） ----

# 阶段1：先砍下界，剔除空液滴/低质量细胞（基因数过少或总 UMI 过少）。
# 只按下界过滤，不碰上界——上界（双细胞、高复杂度、高线粒体）留到 scDblFinder 之后的精细过滤再处理。
cat(sprintf("\n==== 5. 基础过滤（下界: nFeature>%d 且 nCount>%d）====\n",
            FIXED_MIN_FEAT, FIXED_MIN_COUNTS))
merged <- subset(obj, subset = nFeature_RNA > FIXED_MIN_FEAT &
                                nCount_RNA   > FIXED_MIN_COUNTS)
cat(sprintf("基础过滤：%d -> %d 个细胞（移除 %d，%.2f%%）\n",
            ncol(obj), ncol(merged), ncol(obj) - ncol(merged),
            100 * (ncol(obj) - ncol(merged)) / ncol(obj)))

## ---- 6. 双细胞检测（scDblFinder，可选） ----

# 双细胞（doublet）是两个细胞被一个液滴捕获，会干扰下游聚类，通常建议剔除。
# scDblFinder 是 Bioconductor 的常用双细胞检测包；未安装成功时自动降级（统一视为 singlet，不阻塞流程）。
# 注意：在阶段1 基础过滤后的对象上运行（先剔除空液滴，检测更准确、更省内存）。
run_doublet_detection <- function(obj) {
  cat("\n==== 6. 双细胞检测（scDblFinder，约需 30~90 分钟）====\n")
  if (!requireNamespace("scDblFinder", quietly = TRUE)) {          # 包不可用时降级
    warning("scDblFinder 未安装，跳过双细胞检测（所有细胞统一视为 singlet，不做双细胞过滤）")
    obj$scDblFinder.class <- "singlet"                              # 未检测时统一视为 singlet，避免精细过滤误删全部细胞
    obj$discard_doublet   <- FALSE                                  # 不因双细胞过滤任何细胞
    return(obj)
  }

  # 从 Seurat 对象提取 counts 矩阵（v5 语法：LayerData，不能再用旧的 slot= 参数）
  counts <- LayerData(obj, assay = "RNA", layer = "counts")
  cat(sprintf("counts 矩阵维度: %d x %d\n", nrow(counts), ncol(counts)))

  # 构造 SingleCellExperiment 对象（scDblFinder 的输入格式）
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
  sce$sample <- obj$sample                                          # 传入样本分组，scDblFinder 按样本独立检测

  # Windows 下 BiocParallel 多进程会复制大矩阵导致内存翻倍，故默认单线程（SerialParam）；
  # 非 Windows 平台可用 MulticoreParam 并行（线程数=DOUBLET_THREADS）。
  if (.Platform$OS.type == "windows") {
    bp <- BiocParallel::SerialParam(progressbar = TRUE)             # Windows：单线程最稳
  } else {
    bp <- BiocParallel::MulticoreParam(workers = DOUBLET_THREADS, progressbar = TRUE) # 其他平台可并行
  }

  # 运行 scDblFinder：samples 参数指定按样本分组检测，每个样本内部自适应阈值
  sce <- tryCatch(
    scDblFinder::scDblFinder(sce, samples = sce$sample, BPPARAM = bp), # 核心检测调用
    error = function(e) {                                            # 出错时降级
      warning("scDblFinder 运行失败: ", conditionMessage(e), "，跳过双细胞检测")
      return(NULL)
    }
  )
  if (is.null(sce)) {                                                # 降级分支
    obj$scDblFinder.class <- "singlet"                               # 检测失败时统一视为 singlet，避免精细过滤误删全部细胞
    obj$discard_doublet   <- FALSE                                   # 不过滤
    return(obj)
  }

  # 把检测结果写回 Seurat 对象的 meta.data
  cls <- as.character(SummarizedExperiment::colData(sce)$scDblFinder.class) # 提取 singlet/doublet 分类
  obj$scDblFinder.class <- cls                                       # 写入双细胞分类（singlet/doublet）
  obj$discard_doublet   <- cls %in% c("doublet", "Doublet")          # 是否判定为双细胞
  cat(sprintf("双细胞比例: %.2f%%\n", 100 * mean(obj$discard_doublet))) # 打印双细胞占比
  return(obj)
}

# 按配置决定是否运行双细胞检测（在基础过滤后的 merged 对象上运行）
if (RUN_DOUBLET) {                                                   # 开关打开
  merged <- run_doublet_detection(merged)                            # 执行检测
} else {                                                             # 开关关闭
  cat("\n==== 6. 双细胞检测已关闭（RUN_DOUBLET=FALSE）====\n")
  merged$scDblFinder.class <- "singlet"                              # 未检测时统一视为 singlet，避免精细过滤误删全部细胞
  merged$discard_doublet   <- FALSE                                  # 不过滤
}

# 按样本汇总双细胞检测结果（新增输出 doublet_summary_by_sample.csv）
doublet_summary <- merged@meta.data %>%
  dplyr::group_by(sample, group, sample_id) %>%
  dplyr::summarise(
    cells_tested = dplyr::n(),                                       # 参与检测的细胞数（基础过滤后）
    singlet      = sum(scDblFinder.class == "singlet"),              # 判定为单细胞的个数
    doublet      = sum(scDblFinder.class %in% c("doublet", "Doublet")), # 判定为双细胞的个数
    doublet_pct  = 100 * doublet / cells_tested,                     # 双细胞比例（%）
    .groups      = "drop"
  )
write.csv(doublet_summary, file.path(OUT_DIR, "doublet_summary_by_sample.csv"), row.names = FALSE)
cat("  已保存: doublet_summary_by_sample.csv（按样本双细胞检测汇总）\n")
print(as.data.frame(doublet_summary))                                # 控制台打印汇总表

## ---- 7. 精细过滤（砍上界） ----

# 7.1 MAD 离群检测函数（移植自 yijian.R 481-496 行）：
#     按 batch（样本）分组，计算中位数 ± nmads×MAD 区间，区间外判为离群。
#     log=TRUE 时先做 log10(x+1) 转换，更符合计数数据分布。
is_outlier_mad <- function(x, nmads = 3, type = c("both", "lower", "higher"),
                           log = FALSE, batch = NULL) {
  type <- match.arg(type)                                            # 校验方向参数
  xx <- if (log) log10(x + 1) else x                                 # 按需对数转换
  flag <- rep(FALSE, length(x))                                      # 初始化离群标记向量
  grp <- if (is.null(batch)) rep("all", length(x)) else as.character(batch) # 分组向量（无 batch 时视为一组）
  for (g in unique(grp)) {                                           # 逐组计算
    idx <- grp == g                                                  # 该组的细胞下标
    med <- stats::median(xx[idx], na.rm = TRUE)                      # 组内中位数
    md  <- stats::mad(xx[idx], na.rm = TRUE)                         # 组内 MAD
    if (!is.finite(md) || md == 0) md <- 1e-8                        # MAD 为 0 或非有限时兜底，避免除零
    if (type %in% c("both", "lower"))  flag[idx] <- flag[idx] | (xx[idx] <  med - nmads * md) # 下侧离群
    if (type %in% c("both", "higher")) flag[idx] <- flag[idx] | (xx[idx] >  med + nmads * md) # 上侧离群
  }
  flag[is.na(flag)] <- TRUE                                          # NA 值视为离群（保守处理）
  flag                                                                 # 返回逻辑向量
}

# 7.2 阶段3：连同上界一起砍。在基础过滤后的 merged 上执行，按 QC_MODE 选择固定阈值或 MAD。
cat(sprintf("\n==== 7. 精细过滤（模式: %s，砍上界 + 剔除双细胞）====\n", QC_MODE))
if (QC_MODE == "fixed") {
  # ---- 方案 A：固定阈值两段式（阶段1 已砍下界，此处只砍上界 + singlet）----
  cat(sprintf("使用固定阈值: nFeature<%d、nCount<%d、percent.mt<%d、且 singlet\n",
              FIXED_MAX_FEAT, FIXED_MAX_COUNTS, FIXED_MAX_PCT_MT))
  obj_qc <- subset(merged, subset = scDblFinder.class == "singlet" &
                      nFeature_RNA < FIXED_MAX_FEAT &
                      nCount_RNA   < FIXED_MAX_COUNTS &
                      percent.mt   < FIXED_MAX_PCT_MT)
} else if (QC_MODE == "mad") {
  # ---- 方案 B：按样本 3-MAD 自适应（完全对齐 yijian.R run_qc 的判定逻辑，作用于基础过滤后的 merged）----
  cat(sprintf("使用按样本 3-MAD 自适应（nmads=%d, 基因数硬下限=%d）\n", MAD_NMADS, MAD_MIN_FEAT))
  merged$discard_low_features <- is_outlier_mad(merged$nFeature_RNA, MAD_NMADS, "lower", log = TRUE, batch = merged$sample) |
                                 merged$nFeature_RNA <= MAD_MIN_FEAT                                        # 基因数显著低于同批样本中位数 或 低于硬下限
  merged$discard_low_counts   <- is_outlier_mad(merged$nCount_RNA, MAD_NMADS, "lower", log = TRUE, batch = merged$sample) # 总计数显著低于同批样本
  merged$discard_high_mt      <- is_outlier_mad(merged$percent.mt, MAD_NMADS, "higher", batch = merged$sample) # 线粒体比例显著高于同批样本
  obj_qc <- subset(merged, subset = scDblFinder.class == "singlet" &
                      !discard_low_features & !discard_low_counts & !discard_high_mt)
} else {
  stop("QC_MODE 必须是 'fixed' 或 'mad'，当前为: ", QC_MODE)        # 非法模式直接报错
}

# 7.3 打印精细过滤的移除统计
cat(sprintf("精细过滤：%d -> %d 个细胞（移除 %d，%.2f%%）\n",
            ncol(merged), ncol(obj_qc), ncol(merged) - ncol(obj_qc),
            100 * (ncol(merged) - ncol(obj_qc)) / ncol(merged)))

## ---- 8. 绘制 QC 图（过滤前后对比） ----

# 通用保存 PDF 函数：把 ggplot 对象 p 输出到 OUT_DIR 下的 file，尺寸 w×h 英寸
save_pdf <- function(p, file, w, h) {
  pdf(file.path(OUT_DIR, file), width = w, height = h)              # 打开 PDF 设备
  print(p)                                                           # 打印图形
  dev.off()                                                          # 关闭设备（必须，否则文件不完整）
  cat(sprintf("  已保存图表: %s\n", file))                           # 提示保存成功
}

# 8.1 过滤前小提琴图：展示 5 个指标在各样本中的分布（pt.size=0 表示不画散点，加快绘制）
cat("\n==== 8. 绘制 QC 图 ====\n")
p_vln_before <- VlnPlot(obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), # 5 个指标
  group.by = "sample",                          # 按样本分面/分组
  pt.size = 0,                                  # 不画点，只看分布
  ncol = 5) +                                   # 5 列排布
  theme(axis.text.x = element_text(angle = 45, hjust = 1))          # x 轴标签旋转 45° 防重叠
save_pdf(p_vln_before, "01_QC_violin_before_filter.pdf", 18, 5)

# 8.2 过滤前散点图：总计数 vs 基因数，直观显示低质量细胞分布（左下方聚集）
p_scatter <- FeatureScatter(obj,
  feature1 = "nCount_RNA",                      # x 轴：总 UMI 计数
  feature2 = "nFeature_RNA",                    # y 轴：检测基因数
  group.by = "sample")                          # 按样本着色
save_pdf(p_scatter, "02_QC_scatter_before_filter.pdf", 8, 6)

## ---- 9. 汇总统计并保存 ----

# 9.1 按样本（sample）+ 分组（group）+ 样本编号（sample_id）汇总质控统计量。
#     两段式过滤后，分别统计三个阶段的细胞数：
#       cells_before       = 原始对象 obj 中的细胞数
#       cells_after_basic  = 基础过滤（砍下界）后的细胞数（merged）
#       doublet/doublet_pct = 双细胞检测汇总（来自 doublet_summary）
#       kept               = 精细过滤（砍上界+singlet）后的细胞数（obj_qc）
qc_summary <- obj@meta.data %>%                  # 基于原始 meta.data 统计（含过滤前中位数）
  dplyr::group_by(sample, group, sample_id) %>%  # 按样本/分组/编号三列分组
  dplyr::summarise(
    cells_before    = dplyr::n(),                # 过滤前细胞数
    median_features = stats::median(nFeature_RNA), # 过滤前中位基因数
    median_counts   = stats::median(nCount_RNA),   # 过滤前中位总计数
    median_mt       = stats::median(percent.mt),   # 过滤前中位线粒体比例(%)
    .groups         = "drop"                     # 取消分组
  )

# 基础过滤后各样本细胞数（merged）
basic_cnt <- merged@meta.data %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(cells_after_basic = dplyr::n(), .groups = "drop")

# 精细过滤后各样本细胞数（obj_qc）
kept_cnt <- obj_qc@meta.data %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(kept = dplyr::n(), .groups = "drop")

# 左连接合并三阶段计数与双细胞汇总
qc_summary <- qc_summary %>%
  dplyr::left_join(basic_cnt, by = "sample") %>%
  dplyr::left_join(kept_cnt, by = "sample") %>%
  dplyr::left_join(dplyr::select(doublet_summary, sample, doublet, doublet_pct), by = "sample")

# 兜底：若某样本在某个阶段被完全过滤掉（计数为 NA），补 0，避免后续统计报错
qc_summary$cells_after_basic[is.na(qc_summary$cells_after_basic)] <- 0
qc_summary$kept[is.na(qc_summary$kept)] <- 0
qc_summary$doublet[is.na(qc_summary$doublet)] <- 0
qc_summary$doublet_pct[is.na(qc_summary$doublet_pct)] <- 0

# 丢弃细胞数 = 过滤前 - 保留
qc_summary$discarded <- qc_summary$cells_before - qc_summary$kept

write.csv(qc_summary, file.path(OUT_DIR, "01_QC_summary_by_sample.csv"), row.names = FALSE) # 写汇总表
cat("  已保存: 01_QC_summary_by_sample.csv\n")

# 9.2 生成长格式的细胞数统计表 cell_count_by_sample.csv（便于绘图/后续分析直接读取）
cell_count_long <- rbind(
  data.frame(sample_id = qc_summary$sample_id, sample = qc_summary$sample,
             group = qc_summary$group, status = "before", count = qc_summary$cells_before,
             stringsAsFactors = FALSE),          # 过滤前计数
  data.frame(sample_id = qc_summary$sample_id, sample = qc_summary$sample,
             group = qc_summary$group, status = "after", count = qc_summary$kept,
             stringsAsFactors = FALSE)           # 过滤后计数
)
write.csv(cell_count_long, file.path(OUT_DIR, "cell_count_by_sample.csv"), row.names = FALSE) # 写长格式表
cat("  已保存: cell_count_by_sample.csv\n")

# 9.3 过滤后小提琴图：仅展示 3 个核心指标，确认过滤后分布更合理
p_vln_after <- VlnPlot(obj_qc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), # 3 个核心指标
  group.by = "sample",                          # 按样本分组
  pt.size = 0,                                  # 不画点
  ncol = 3) +                                   # 3 列排布
  theme(axis.text.x = element_text(angle = 45, hjust = 1))      # 旋转标签
save_pdf(p_vln_after, "03_QC_violin_after_filter.pdf", 14, 5)

# 9.4 过滤前后细胞数柱状图：按样本编号展示 before/after 对比，直观看出各样本过滤比例
sample_id_levels <- meta$SampleID               # 用 xlsx 中的原始样本编号顺序作为 x 轴顺序
cell_count_long$sample_id <- factor(cell_count_long$sample_id, levels = sample_id_levels) # 转因子固定顺序
p_bar <- ggplot(cell_count_long, aes(x = sample_id, y = count, fill = status)) + # 柱状图
  geom_col(position = "dodge", width = 0.7) +    # 分组柱状图
  scale_fill_manual(values = c(before = "grey60", after = "#3C5488")) + # 灰色=过滤前，深蓝=过滤后
  theme_minimal() +                              # 简洁主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +        # 旋转 x 轴标签
  labs(title = "Cell count before/after QC by sample",              # 图标题（英文避免字体问题）
       x = "Sample ID", y = "Cell count", fill = "Status")          # 坐标轴与图例标签
save_pdf(p_bar, "04_QC_cell_count_by_sample.pdf", 14, 6)

# 9.5 保存过滤后的 Seurat 对象（meta.data 完整保留 sample/group/sample_id/percent.*/scDblFinder.class 等全部信息，
#     分组与样本编号的对应关系原样保留，供后续标准化/整合/聚类步骤直接使用）
saveRDS(obj_qc, file.path(OUT_DIR, "02_seurat_qc.rds"))
cat(sprintf("  已保存过滤后对象: %s\n", file.path(OUT_DIR, "02_seurat_qc.rds")))

## ---- 10. 完成 ----

# 打印最终结果摘要
cat("\n==== 质控完成 ====\n")
cat(sprintf("过滤前细胞数: %d\n", ncol(obj)))                        # 过滤前
cat(sprintf("过滤后细胞数: %d\n", ncol(obj_qc)))                     # 过滤后
cat(sprintf("保留比例: %.2f%%\n", 100 * ncol(obj_qc) / ncol(obj)))   # 保留比例
cat("\n输出文件清单（位于 output/ 目录）：\n")
cat("  02_seurat_qc.rds                  - 过滤后的 Seurat 对象\n")
cat("  01_QC_summary_by_sample.csv       - 按样本/分组/编号的质控汇总表\n")
cat("  cell_count_by_sample.csv          - 过滤前后细胞数长格式表\n")
cat("  doublet_summary_by_sample.csv     - 按样本双细胞检测汇总表\n")
cat("  01_QC_violin_before_filter.pdf    - 过滤前 QC 小提琴图\n")
cat("  02_QC_scatter_before_filter.pdf   - 过滤前计数-基因数散点图\n")
cat("  03_QC_violin_after_filter.pdf     - 过滤后 QC 小提琴图\n")
cat("  04_QC_cell_count_by_sample.pdf    - 过滤前后细胞数柱状图\n")
cat("\n02_quality_control 完成。\n")

