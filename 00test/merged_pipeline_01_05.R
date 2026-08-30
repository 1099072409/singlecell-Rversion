# =====================================================================
# merged_pipeline_01_05.R
# 单细胞分析流程 - 01~05 合并版（单脚本）
#
# 合并来源（Rversion/ 下原有 5 个独立脚本，逻辑原样保留）：
#   01_initiation.R            -> SECTION 01
#   02_quality_control.R       -> SECTION 02
#   03_extract_cd45.R          -> SECTION 03
#   04_integration_clustering.R-> SECTION 04
#   05_celltype_annotation.R   -> SECTION 05
#   06_parameter_sweep.R        -> SECTION 06 (DIMS x RESOLUTIONS x UMAP 参数扫描)
#
# 本合并版的关键约定（详见同目录 README_merged.md）：
#   1) 样本范围：仅处理 data/ 下 GSE203115_* 三个样本目录（SAMPLE_PREFIX）。
#   2) 逻辑：仅主流程，不含 exclude 剔除分支。
#   3) 单脚本一次跑完 01->05；段间 rds 走磁盘传递（与原有模块化设计一致）。
#   4) 断点：SECTION 03 末尾 saveRDS(03_CD45_positive.rds) 为关键存盘点。
#   5) 续跑：RESUME_FROM 开关可跳过前段（如前段 rds 已存在）。
#   6) 修正：原 01 的 DATA_DIR 少一级 ".." 已修正为两级 ".."（见全局 CONFIG）。
#   7) 保留：04 三处硬编码（DIMS=1:30 / theta=5 / DEFAULT_RESOLUTION=0.4）
#      与 ypN0/ypN+ 配色均原样保留，保证结果可复现。
#
# 运行方式（R 不在 PATH，须用绝对路径）：
#   "C:/Program Files/R/R-4.4.3/bin/Rscript.exe" "A:/Workbuddy/singlecell/Rversion/00test/merged_pipeline_01_05.R"
# =====================================================================


## ====================================================================
## 全局配置（CONFIG）—— 所有路径以脚本自身所在目录为基准推算
## ====================================================================

# 自动定位脚本所在目录（兼容 Rscript 与 RStudio 两种运行方式）
script_dir <- tryCatch({
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
    dirname(rstudioapi::getActiveDocumentContext()$path)
  } else {
    f <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
    if (is.null(f) || !nzchar(f)) {
      arg <- commandArgs(trailingOnly = FALSE)
      f <- sub("^--file=", "", grep("^--file=", arg, value = TRUE)[1])
    }
    if (is.null(f) || !nzchar(f) || is.na(f)) {
      getwd()
    } else {
      dirname(normalizePath(f, mustWork = FALSE))
    }
  }
}, error = function(e) getwd())

cat(sprintf("脚本目录: %s\n", script_dir))

# ---- 输入/输出根路径 ----
# 修正说明：原 01_initiation.R:34 写的是 file.path(script_dir, "..", "data")，
# 少一级 ".."（脚本被移入子目录后失效）。本合并版统一用两级 ".."。
DATA_DIR  <- file.path(script_dir, "..", "..", "data")        # 项目根 data/

# 分组主表：不在源码中硬编码中文路径（R 在 GBK 区域下读取 UTF-8 源码会把中文文件名
# 解析为乱码，导致 file.exists 报 "name too long"）。改用 list.files 按扩展名发现
# （返回的是系统原生编码，可靠），并支持 XLSX_PATH 环境变量 / 命令行覆盖。
XLSX_PATH <- Sys.getenv("XLSX_PATH", unset = "")
if (!nzchar(XLSX_PATH) || !file.exists(XLSX_PATH)) {
  cand <- list.files(DATA_DIR, pattern = "\\.xlsx$", full.names = TRUE, ignore.case = TRUE)
  cand <- cand[!grepl("backup", cand, ignore.case = TRUE)]   # 排除 backup 目录
  if (length(cand) >= 1) XLSX_PATH <- cand[1]
}

# 各段输出目录（各自独立子目录，避免 02/03 同名 pdf 冲突）
OUT_01 <- file.path(script_dir, "output", "01_initiation")
OUT_02 <- file.path(script_dir, "output", "02_quality_control")
OUT_03 <- file.path(script_dir, "output", "03_extract_cd45")
OUT_04 <- file.path(script_dir, "output", "04_integration_clustering")
OUT_05 <- file.path(script_dir, "output", "05_celltype_annotation")

# 各段产物 rds 路径（段间走磁盘传递，与"03 存盘后"断点语义一致）
RDS_01 <- file.path(OUT_01, "01_seurat_combined.rds")
RDS_02 <- file.path(OUT_02, "02_seurat_qc.rds")
RDS_03 <- file.path(OUT_03, "03_CD45_positive.rds")   # ★ SECTION 03 末尾存盘（关键断点）
RDS_04 <- file.path(OUT_04, "04_CD45_integrated.rds")
RDS_05 <- file.path(OUT_05, "05_annotated.rds")

# ---- 断点续跑开关 ----
# ""      = 从头跑（默认）
# "01".."05" = 从该段开始（其前各段的 rds 需已存在）
RESUME_FROM <- ""

# ---- 提前停止开关（用于只跑前段生成中间 rds）----
# ""      = 不提前停止（重头运行：01 -> 06 全程跑完）
# "03"    = 跑到 SECTION 03 末尾保存 03_CD45_positive.rds 后停止
STOP_AFTER <- ""

# ---- 样本范围（重头运行：仅纳入 data/ 下以 "GSE" 开头的样本，即 GSE203115_*）----
SAMPLE_PREFIX <- "GSE"

# 允许命令行覆盖开关：Rscript merged.R STOP_AFTER=03   /   RESUME_FROM=03
args <- commandArgs(trailingOnly = TRUE)
for (a in args) {
  if (grepl("^STOP_AFTER=", a)) STOP_AFTER <- sub("^STOP_AFTER=", "", a)
  if (grepl("^RESUME_FROM=", a)) RESUME_FROM <- sub("^RESUME_FROM=", "", a)
}

# ---- 创建所有输出目录 ----
for (d in c(OUT_01, OUT_02, OUT_03, OUT_04, OUT_05)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- 通用依赖安装函数（全局，供各段复用）----
ensure_pkg <- function(pkg, bioc = FALSE) {
  if (requireNamespace(pkg, quietly = TRUE)) return(TRUE)
  cat(sprintf("  正在安装依赖包: %s ...\n", pkg))
  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  requireNamespace(pkg, quietly = TRUE)
}

# SeuratExtend 专用安装函数（GitHub 源码包，可能失败则兜底）
ensure_seuratextend <- function(repo, upgrade) {
  if (requireNamespace("SeuratExtend", quietly = TRUE)) return(TRUE)
  cat(sprintf("  正在安装 SeuratExtend（GitHub 源码包 %s，需数分钟与 Rtools）...\n", repo))
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  ok <- tryCatch({
    remotes::install_github(repo, upgrade = upgrade, quiet = TRUE)
    requireNamespace("SeuratExtend", quietly = TRUE)
  }, error = function(e) {
    cat("  SeuratExtend 安装失败:", conditionMessage(e), "\n")
    FALSE
  })
  ok
}

# 通用 PDF 保存函数（使用当前段设置的全局 OUT_DIR）
save_pdf <- function(p, file, w, h) {
  pdf(file.path(OUT_DIR, file), width = w, height = h)
  print(p)
  dev.off()
  cat(sprintf("  已保存图表: %s\n", file))
}

# 判断是否运行某段（依据 RESUME_FROM）
run_from <- function(stage) {
  order <- c("01", "02", "03", "04", "05", "06")
  if (RESUME_FROM == "") return(TRUE)
  idx_run   <- which(order == RESUME_FROM)
  idx_stage <- which(order == stage)
  if (length(idx_run) == 0) return(TRUE)   # RESUME_FROM 非法时视为从头
  idx_stage >= idx_run
}

cat(sprintf("样本范围前缀: %s | 续跑起点: %s\n", SAMPLE_PREFIX,
            ifelse(RESUME_FROM == "", "从头(01)", RESUME_FROM)))


## ====================================================================
## SECTION 01 — INITIATION（数据整理 + 读入 Seurat）
##   读 data/ 下 GSE203115_* 样本 -> Read10X -> CreateSeuratObject
##   -> merge -> 保存 01_seurat_combined.rds
## ====================================================================
if (run_from("01")) {
  cat(sprintf("\n########## SECTION 01 — INITIATION ##########\n"))

  ## ---- 0. 配置 ----
  OUT_DIR   <- OUT_01
  OUT_RDS   <- RDS_01
  MIN_CELLS    <- 3        # CreateSeuratObject: 至少在这么多细胞中出现的基因才保留
  MIN_FEATURES <- 200      # CreateSeuratObject: 至少检测到这么多基因的细胞才保留

  ## ---- 1. 依赖与目录 ----
  required_pkgs <- c("Seurat", "readxl", "R.utils")
  for (p in required_pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
  library(Seurat)
  library(readxl)
  library(R.utils)

  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  cat(sprintf("脚本目录: %s\n", script_dir))
  cat(sprintf("工作目录: %s\n", getwd()))
  cat(sprintf("数据目录: %s\n", normalizePath(DATA_DIR)))

  ## ---- 2. 读取分组表 ----
  if (!file.exists(XLSX_PATH)) stop("找不到分组表: ", XLSX_PATH)
  meta <- read_xlsx(XLSX_PATH)

  # 兼容中/英文列名，自动匹配 Sample / Group / 编号 三列
  find_col <- function(candidates) {
    for (c in candidates) {
      idx <- grep(c, names(meta), ignore.case = TRUE)
      if (length(idx) >= 1) return(names(meta)[idx[1]])
    }
    return(NA_character_)
  }
  sample_col <- find_col(c("Sample", "样本", "样本名"))
  group_col  <- find_col(c("Group", "分组", "组别"))
  id_col     <- find_col(c("样本编号", "编号", "SampleID", "ID"))
  cat(sprintf("分组表列映射 -> Sample: %s | Group: %s | 编号: %s\n",
              sample_col, group_col, id_col))

  meta <- meta[, c(sample_col, group_col, id_col)]
  colnames(meta) <- c("Sample", "Group", "SampleID")
  meta$Sample   <- as.character(meta$Sample)
  meta$Group    <- as.character(meta$Group)
  meta$SampleID <- as.character(meta$SampleID)

  ## ---- 3. 格式统一：压缩 HRS（genes.tsv -> features.tsv） ----
  # 列出 data/ 下所有包含 barcodes 的样本文件夹
  all_dirs <- list.dirs(DATA_DIR, full.names = FALSE, recursive = FALSE)
  sample_dirs <- all_dirs[
    file.exists(file.path(DATA_DIR, all_dirs, "barcodes.tsv")) |
    file.exists(file.path(DATA_DIR, all_dirs, "barcodes.tsv.gz"))
  ]
  # 【范围限定】本次合并仅处理 GSE203115 开头的样本
  sample_dirs <- sample_dirs[startsWith(sample_dirs, SAMPLE_PREFIX)]
  cat(sprintf("发现 %d 个 GSE203115 样本文件夹。\n", length(sample_dirs)))

  gzipped_log <- c()
  for (s in sample_dirs) {
    d <- file.path(DATA_DIR, s)

    # HRS 用 genes.tsv，统一改名为 features.tsv
    if (file.exists(file.path(d, "genes.tsv"))) {
      file.rename(file.path(d, "genes.tsv"), file.path(d, "features.tsv"))
    }

    # 对三个核心文件，仅当 .gz 不存在时才原地压缩（幂等）
    for (base in c("barcodes.tsv", "features.tsv", "matrix.mtx")) {
      f  <- file.path(d, base)
      fz <- file.path(d, paste0(base, ".gz"))
      if (file.exists(f) && !file.exists(fz)) {
        R.utils::gzip(f, overwrite = TRUE)   # 压缩并删除原文件
        gzipped_log <- c(gzipped_log, s)
      }
    }
  }
  if (length(gzipped_log) > 0) {
    cat("本次新压缩的样本:", paste(unique(gzipped_log), collapse = ", "), "\n")
  } else {
    cat("没有需要新压缩的样本（GSE203115 样本已是 .gz，或此前已压缩）。\n")
  }

  ## ---- 4. 交叉校验与读入 ----
  obj_list <- list()
  skipped  <- c()
  for (s in sample_dirs) {
    d <- file.path(DATA_DIR, s)

    # 确认三个 .gz 都齐全
    needed <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
    if (!all(file.exists(file.path(d, needed)))) {
      cat(sprintf("警告: 样本 %s 缺少必要文件，跳过。\n", s))
      skipped <- c(skipped, s)
      next
    }

    # 取分组信息（匹配不到则填 NA，保留容错）
    mrow <- meta[meta$Sample == s, ]
    if (nrow(mrow) == 0) {
      cat(sprintf("警告: 样本 %s 不在分组表中，分组信息填 NA。\n", s))
      grp <- NA_character_; sid <- NA_character_
    } else {
      grp <- mrow$Group; sid <- mrow$SampleID
    }

    mat <- Read10X(data.dir = d)
    obj <- CreateSeuratObject(mat, project = s,
                              min.cells = MIN_CELLS, min.features = MIN_FEATURES)
    obj$group     <- grp
    obj$sample_id <- sid

    # 用样本名做 cell 前缀，避免跨样本 barcode 冲突
    colnames(obj) <- paste0(s, "_", colnames(obj))
    obj_list[[s]] <- obj
    cat(sprintf("  %-12s cells=%-7d features=%-6d group=%s\n",
                s, ncol(obj), nrow(obj), grp))
  }

  # 分组表中有但无文件夹的样本（本次仅 GSE203115 子集，其余为预期内）
  missing <- setdiff(meta$Sample, sample_dirs)
  if (length(missing) > 0) {
    cat("注意: 本次仅处理 GSE203115 子集，分组表其余样本未读入:",
        paste(missing, collapse = ", "), "\n")
  }

  ## ---- 5. 合并与保存 ----
  if (length(obj_list) == 0) stop("没有任何样本成功读入，已终止。")

  combined <- obj_list[[1]]
  if (length(obj_list) > 1) {
    combined <- merge(combined, y = obj_list[-1])
  }
  cat(sprintf("合并后维度: %d features x %d cells\n", nrow(combined), ncol(combined)))
  cat("分组分布:\n")
  print(table(combined$group, useNA = "ifany"))

  saveRDS(combined, OUT_RDS)
  cat(sprintf("已保存 Seurat 对象: %s\n", OUT_RDS))
  cat("01 initiation 完成。\n")
}


## ====================================================================
## SECTION 02 — QUALITY CONTROL（质量控制）
##   读取 01_seurat_combined.rds -> 计算 QC 指标 -> 双细胞检测 ->
##   过滤 -> 保存 02_seurat_qc.rds + csv/pdf
## ====================================================================
if (run_from("02")) {
  cat(sprintf("\n########## SECTION 02 — QUALITY CONTROL ##########\n"))

  ## ---- 0. 配置（按需修改） ----
  INPUT_RDS <- RDS_01    # 第一步合并好的 rds
  OUT_DIR   <- OUT_02    # 本步输出目录

  # ---- 质控指标（人类基因命名规则，切换物种时修改） ----
  MT_PATTERN   <- "^MT-"        # 线粒体基因前缀（human）
  RIBO_PATTERN <- "^RP[SL]"     # 核糖体蛋白基因前缀（RPS/RPL）
  HB_PATTERN   <- "^HB[ABDEGQZ]"# 血红蛋白基因前缀

  # ---- 质控模式与阈值 ----
  QC_MODE          <- "fixed"   # 质控模式："fixed"=固定阈值；"mad"=按样本 3-MAD 自适应
  FIXED_MIN_FEAT   <- 200       # 细胞最少检测到的基因数（下限）
  FIXED_MAX_FEAT   <- 6000      # 细胞最多检测到的基因数（上限）
  FIXED_MIN_COUNTS <- 500       # 细胞最少总 UMI 计数
  FIXED_MAX_PCT_MT <- 20        # 细胞线粒体基因比例上限（%）
  MAD_NMADS        <- 3         # MAD 离群判定倍数
  MAD_MIN_FEAT     <- 200       # mad 模式下 nFeature_RNA 的硬性下限

  # ---- 双细胞检测 ----
  RUN_DOUBLET      <- TRUE      # 是否运行 scDblFinder 双细胞检测
  DOUBLET_THREADS  <- 2         # 双细胞检测线程数（Windows 下自动降级为单线程）

  # ---- 其他 ----
  STOP_ON_MISMATCH <- FALSE     # 元数据校验不一致时：TRUE=报错终止；FALSE=警告后继续
  SEED             <- 123       # 随机种子

  ## ---- 1. 依赖检查与加载 ----
  for (p in c("Seurat", "readxl", "dplyr", "ggplot2", "patchwork"))  ensure_pkg(p, bioc = FALSE)
  for (p in c("BiocParallel", "SingleCellExperiment", "SummarizedExperiment", "scDblFinder")) ensure_pkg(p, bioc = TRUE)

  library(Seurat)
  library(readxl)
  library(dplyr)
  library(ggplot2)

  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  cat(sprintf("脚本目录: %s\n", script_dir))
  cat(sprintf("工作目录: %s\n", getwd()))
  cat(sprintf("输出目录: %s\n", normalizePath(OUT_DIR)))
  cat(sprintf("质控模式: %s\n", QC_MODE))

  options(future.globals.maxSize = 8 * 1024^3)

  tryCatch({
    w  <- system("wmic OS get FreePhysicalMemory /value", intern = TRUE)
    fr <- as.numeric(sub(".*=", "", w[grepl("FreePhysicalMemory", w)]))
    if (length(fr) == 1 && !is.na(fr)) {
      free_gb <- fr / 1024^2
      cat(sprintf("当前可用物理内存: %.1f GB\n", free_gb))
      if (free_gb < 8) {
        cat("提示: 可用内存偏少，本步大对象约需 5~8GB，建议关闭其他程序后再运行。\n")
      }
    }
  }, error = function(e) invisible(NULL))

  ## ---- 2. 读取分组表（data/样本分组信息.xlsx） ----
  if (!file.exists(XLSX_PATH)) stop("找不到分组表: ", XLSX_PATH)
  meta <- read_xlsx(XLSX_PATH)

  find_col <- function(candidates, df) {
    for (c in candidates) {
      idx <- grep(c, names(df), ignore.case = TRUE)
      if (length(idx) >= 1) return(names(df)[idx[1]])
    }
    return(NA_character_)
  }
  sample_col <- find_col(c("Sample", "样本", "样本名"), meta)
  group_col  <- find_col(c("Group", "分组", "组别"), meta)
  id_col     <- find_col(c("样本编号", "编号", "SampleID", "ID"), meta)
  cat(sprintf("分组表列映射 -> Sample: %s | Group: %s | 编号: %s\n",
              sample_col, group_col, id_col))
  meta <- meta[, c(sample_col, group_col, id_col)]
  colnames(meta) <- c("Sample", "Group", "SampleID")
  meta$Sample   <- as.character(meta$Sample)
  meta$Group    <- as.character(meta$Group)
  meta$SampleID <- as.character(meta$SampleID)
  cat(sprintf("分组表共 %d 个样本。\n", nrow(meta)))

  ## ---- 3. 读入合并对象并校验元数据 ----
  if (!file.exists(INPUT_RDS)) stop("找不到输入 rds: ", INPUT_RDS)
  cat("正在读取 01_seurat_combined.rds（约需 1~3 分钟）...\n")
  obj <- readRDS(INPUT_RDS)
  cat(sprintf("读入完成：%d 个基因 x %d 个细胞\n", nrow(obj), ncol(obj)))

  if (DefaultAssay(obj) != "RNA") {
    warning("默认 assay 不是 RNA，已切换为 RNA")
    DefaultAssay(obj) <- "RNA"
  }

  cat(sprintf("当前 RNA layer: %s\n", paste(Layers(obj, assay = "RNA"), collapse = ", ")))
  if (length(Layers(obj, assay = "RNA")) > 1) {
    cat("检测到多个 layer（merge 产物），执行 JoinLayers 合并...\n")
    obj <- JoinLayers(obj, assay = "RNA")
  }
  cat(sprintf("合并后 RNA layer: %s\n", paste(Layers(obj, assay = "RNA"), collapse = ", ")))

  obj$sample <- as.character(obj$orig.ident)

  check_meta_consistency <- function(obj, meta) {
    cat("\n==== 3.5 元数据一致性校验（rds meta.data vs 样本分组信息.xlsx）====\n")
    md <- obj@meta.data
    rds_tbl <- md %>%
      dplyr::group_by(orig.ident) %>%
      dplyr::summarise(
        rds_group     = paste(unique(stats::na.omit(group)), collapse = "/"),
        rds_sample_id = paste(unique(stats::na.omit(sample_id)), collapse = "/"),
        n_cells       = dplyr::n(),
        .groups       = "drop"
      )
    names(rds_tbl)[names(rds_tbl) == "orig.ident"] <- "Sample"
    cmp <- dplyr::left_join(rds_tbl, meta, by = "Sample")
    cat("逐样本对照表（列为空表示仅另一侧存在）：\n")
    print(as.data.frame(cmp))
    only_rds       <- setdiff(rds_tbl$Sample, meta$Sample)
    only_xlsx      <- setdiff(meta$Sample, rds_tbl$Sample)
    na_cells       <- sum(is.na(md$group) | is.na(md$sample_id))
    mismatch_grp   <- sum(!is.na(cmp$Group) & !is.na(cmp$rds_group) &
                            as.character(cmp$Group) != as.character(cmp$rds_group))
    mismatch_id    <- sum(!is.na(cmp$SampleID) & !is.na(cmp$rds_sample_id) &
                            as.character(cmp$SampleID) != as.character(cmp$rds_sample_id))
    if (length(only_rds) > 0 || length(only_xlsx) > 0 || na_cells > 0 ||
        mismatch_grp > 0 || mismatch_id > 0) {
      msg <- sprintf(paste0("元数据不一致：仅rds有样本=[%s]；仅xlsx有样本=[%s]；",
                            "分组/编号为NA的细胞=%d；分组不一致样本=%d；编号不一致样本=%d"),
                     paste(only_rds, collapse = ","), paste(only_xlsx, collapse = ","),
                     na_cells, mismatch_grp, mismatch_id)
      if (STOP_ON_MISMATCH) stop(msg) else warning(msg)
    } else {
      cat("校验通过：rds 与 xlsx 的样本集合、分组、样本编号完全一致。\n")
    }
    invisible(cmp)
  }
  meta_check <- check_meta_consistency(obj, meta)

  ## ---- 4. 计算质控指标 ----
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = MT_PATTERN)
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = RIBO_PATTERN)
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = HB_PATTERN)

  cat("\n==== QC 指标分布（过滤前）====\n")
  print(summary(obj@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.hb")]))

  qc_q <- obj@meta.data %>%
    dplyr::group_by(sample, group, sample_id) %>%
    dplyr::summarise(
      nCount_q05      = stats::quantile(nCount_RNA, 0.05),
      nCount_median   = stats::median(nCount_RNA),
      nCount_q95      = stats::quantile(nCount_RNA, 0.95),
      nFeature_q05    = stats::quantile(nFeature_RNA, 0.05),
      nFeature_median = stats::median(nFeature_RNA),
      nFeature_q95    = stats::quantile(nFeature_RNA, 0.95),
      mt_q05          = stats::quantile(percent.mt, 0.05),
      mt_median       = stats::median(percent.mt),
      mt_q95          = stats::quantile(percent.mt, 0.95),
      .groups         = "drop"
    )
  write.csv(qc_q, file.path(OUT_DIR, "00_QC_quantiles_by_sample.csv"), row.names = FALSE)
  cat("  已保存: 00_QC_quantiles_by_sample.csv（各样本 QC 指标分位数，供设定过滤阈值参考）\n")

  ## ---- 5. 双细胞检测（scDblFinder，可选） ----
  run_doublet_detection <- function(obj) {
    cat("\n==== 5. 双细胞检测（scDblFinder，约需 30~90 分钟）====\n")
    if (!requireNamespace("scDblFinder", quietly = TRUE)) {
      warning("scDblFinder 未安装，跳过双细胞检测（doublet_call 标记为 not_tested）")
      obj$doublet_call    <- "not_tested"
      obj$discard_doublet <- FALSE
      return(obj)
    }
    counts <- LayerData(obj, assay = "RNA", layer = "counts")
    cat(sprintf("counts 矩阵维度: %d x %d\n", nrow(counts), ncol(counts)))
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))
    sce$sample <- obj$sample
    if (.Platform$OS.type == "windows") {
      bp <- BiocParallel::SerialParam(progressbar = TRUE)
    } else {
      bp <- BiocParallel::MulticoreParam(workers = DOUBLET_THREADS, progressbar = TRUE)
    }
    sce <- tryCatch(
      scDblFinder::scDblFinder(sce, samples = sce$sample, BPPARAM = bp),
      error = function(e) {
        warning("scDblFinder 运行失败: ", conditionMessage(e), "，跳过双细胞检测")
        return(NULL)
      }
    )
    if (is.null(sce)) {
      obj$doublet_call    <- "not_tested"
      obj$discard_doublet <- FALSE
      return(obj)
    }
    cls <- as.character(SummarizedExperiment::colData(sce)$scDblFinder.class)
    obj$doublet_call    <- cls
    obj$discard_doublet <- cls %in% c("doublet", "Doublet")
    cat(sprintf("双细胞比例: %.2f%%\n", 100 * mean(obj$discard_doublet)))
    return(obj)
  }
  if (RUN_DOUBLET) {
    obj <- run_doublet_detection(obj)
  } else {
    cat("\n==== 5. 双细胞检测已关闭（RUN_DOUBLET=FALSE）====\n")
    obj$doublet_call    <- "not_tested"
    obj$discard_doublet <- FALSE
  }

  ## ---- 6. 过滤判定（低质量细胞） ----
  is_outlier_mad <- function(x, nmads = 3, type = c("both", "lower", "higher"),
                             log = FALSE, batch = NULL) {
    type <- match.arg(type)
    xx <- if (log) log10(x + 1) else x
    flag <- rep(FALSE, length(x))
    grp <- if (is.null(batch)) rep("all", length(x)) else as.character(batch)
    for (g in unique(grp)) {
      idx <- grp == g
      med <- stats::median(xx[idx], na.rm = TRUE)
      md  <- stats::mad(xx[idx], na.rm = TRUE)
      if (!is.finite(md) || md == 0) md <- 1e-8
      if (type %in% c("both", "lower"))  flag[idx] <- flag[idx] | (xx[idx] <  med - nmads * md)
      if (type %in% c("both", "higher")) flag[idx] <- flag[idx] | (xx[idx] >  med + nmads * md)
    }
    flag[is.na(flag)] <- TRUE
    flag
  }

  cat(sprintf("\n==== 6. 过滤判定（模式: %s）====\n", QC_MODE))
  if (QC_MODE == "fixed") {
    cat("使用固定阈值: nFeature∈[200,6000]、nCount>500、percent.mt<20%\n")
    obj$discard_low_features <- obj$nFeature_RNA < FIXED_MIN_FEAT | obj$nFeature_RNA > FIXED_MAX_FEAT
    obj$discard_low_counts   <- obj$nCount_RNA   < FIXED_MIN_COUNTS
    obj$discard_high_mt      <- obj$percent.mt   > FIXED_MAX_PCT_MT
  } else if (QC_MODE == "mad") {
    cat(sprintf("使用按样本 3-MAD 自适应（nmads=%d, 基因数硬下限=%d）\n", MAD_NMADS, MAD_MIN_FEAT))
    obj$discard_low_features <- is_outlier_mad(obj$nFeature_RNA, MAD_NMADS, "lower", log = TRUE, batch = obj$sample) |
                                obj$nFeature_RNA <= MAD_MIN_FEAT
    obj$discard_low_counts   <- is_outlier_mad(obj$nCount_RNA, MAD_NMADS, "lower", log = TRUE, batch = obj$sample)
    obj$discard_high_mt      <- is_outlier_mad(obj$percent.mt, MAD_NMADS, "higher", batch = obj$sample)
  } else {
    stop("QC_MODE 必须是 'fixed' 或 'mad'，当前为: ", QC_MODE)
  }

  obj$discard <- obj$discard_low_features | obj$discard_low_counts |
                 obj$discard_high_mt | obj$discard_doublet

  cat("各过滤原因丢弃的细胞数：\n")
  print(data.frame(
    reason      = c("low_features", "low_counts", "high_mt", "doublet", "total"),
    n_discarded = c(sum(obj$discard_low_features), sum(obj$discard_low_counts),
                    sum(obj$discard_high_mt), sum(obj$discard_doublet), sum(obj$discard)),
    row.names   = NULL
  ))
  cat(sprintf("总体：过滤前 %d 个细胞，将丢弃 %d 个（%.2f%%），保留 %d 个。\n",
              ncol(obj), sum(obj$discard), 100 * mean(obj$discard), sum(!obj$discard)))

  ## ---- 7. 绘制 QC 图（过滤前后对比） ----
  cat("\n==== 7. 绘制 QC 图 ====\n")
  p_vln_before <- VlnPlot(obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"),
    group.by = "sample", pt.size = 0, ncol = 5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_pdf(p_vln_before, "01_QC_violin_before_filter.pdf", 18, 5)

  p_scatter <- FeatureScatter(obj,
    feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample")
  save_pdf(p_scatter, "02_QC_scatter_before_filter.pdf", 8, 6)

  ## ---- 8. 执行过滤、汇总统计并保存 ----
  obj_qc <- subset(obj, subset = discard == FALSE)
  cat(sprintf("\n过滤完成：%d 个细胞 -> %d 个细胞（保留 %.2f%%）\n",
              ncol(obj), ncol(obj_qc), 100 * ncol(obj_qc) / ncol(obj)))

  qc_summary <- obj@meta.data %>%
    dplyr::group_by(sample, group, sample_id) %>%
    dplyr::summarise(
      cells_before    = dplyr::n(),
      discarded       = sum(discard),
      kept            = sum(!discard),
      median_features = stats::median(nFeature_RNA),
      median_counts   = stats::median(nCount_RNA),
      median_mt       = stats::median(percent.mt),
      .groups         = "drop"
    )
  write.csv(qc_summary, file.path(OUT_DIR, "01_QC_summary_by_sample.csv"), row.names = FALSE)
  cat("  已保存: 01_QC_summary_by_sample.csv\n")

  cell_count_long <- rbind(
    data.frame(sample_id = qc_summary$sample_id, sample = qc_summary$sample,
               group = qc_summary$group, status = "before", count = qc_summary$cells_before,
               stringsAsFactors = FALSE),
    data.frame(sample_id = qc_summary$sample_id, sample = qc_summary$sample,
               group = qc_summary$group, status = "after", count = qc_summary$kept,
               stringsAsFactors = FALSE)
  )
  write.csv(cell_count_long, file.path(OUT_DIR, "cell_count_by_sample.csv"), row.names = FALSE)
  cat("  已保存: cell_count_by_sample.csv\n")

  p_vln_after <- VlnPlot(obj_qc,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = "sample", pt.size = 0, ncol = 3) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_pdf(p_vln_after, "03_QC_violin_after_filter.pdf", 14, 5)

  sample_id_levels <- meta$SampleID
  cell_count_long$sample_id <- factor(cell_count_long$sample_id, levels = sample_id_levels)
  p_bar <- ggplot(cell_count_long, aes(x = sample_id, y = count, fill = status)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = c(before = "grey60", after = "#3C5488")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Cell count before/after QC by sample",
         x = "Sample ID", y = "Cell count", fill = "Status")
  save_pdf(p_bar, "04_QC_cell_count_by_sample.pdf", 14, 6)

  saveRDS(obj_qc, file.path(OUT_DIR, "02_seurat_qc.rds"))
  cat(sprintf("  已保存过滤后对象: %s\n", file.path(OUT_DIR, "02_seurat_qc.rds")))

  ## ---- 9. 完成 ----
  cat("\n==== 质控完成 ====\n")
  cat(sprintf("过滤前细胞数: %d\n", ncol(obj)))
  cat(sprintf("过滤后细胞数: %d\n", ncol(obj_qc)))
  cat(sprintf("保留比例: %.2f%%\n", 100 * ncol(obj_qc) / ncol(obj)))
  cat("\n输出文件清单（位于 output/02_quality_control/ 目录）：\n")
  cat("  02_seurat_qc.rds                  - 过滤后的 Seurat 对象\n")
  cat("  01_QC_summary_by_sample.csv       - 按样本/分组/编号的质控汇总表\n")
  cat("  cell_count_by_sample.csv          - 过滤前后细胞数长格式表\n")
  cat("  01_QC_violin_before_filter.pdf    - 过滤前 QC 小提琴图\n")
  cat("  02_QC_scatter_before_filter.pdf   - 过滤前计数-基因数散点图\n")
  cat("  03_QC_violin_after_filter.pdf     - 过滤后 QC 小提琴图\n")
  cat("  04_QC_cell_count_by_sample.pdf    - 过滤前后细胞数柱状图\n")
  cat("\n02 quality_control 完成。\n")
}


## ====================================================================
## SECTION 03 — EXTRACT CD45+（提取 CD45+ 免疫细胞）
##   读取 02_seurat_qc.rds -> 基于 PTPRC 表达判定 CD45+ ->
##   提取并 ★保存 03_CD45_positive.rds（关键断点）-> 输出 QC 图/统计
## ====================================================================
if (run_from("03")) {
  cat(sprintf("\n########## SECTION 03 — EXTRACT CD45+ ##########\n"))

  ## ---- 0. 配置（按需修改） ----
  INPUT_RDS <- RDS_02    # 第二步质控后的 rds
  OUT_DIR   <- OUT_03    # 本步输出目录

  # ---- CD45+ 判定参数 ----
  CD45_GENE      <- "PTPRC"   # CD45 标志基因名（人类；小鼠为 Ptprc）
  CD45_MIN_COUNTS <- 0        # PTPRC counts 大于该值即判定为 CD45+

  SEED <- 123

  ## ---- 1. 依赖检查与加载 ----
  for (p in c("Seurat", "dplyr", "ggplot2", "patchwork")) ensure_pkg(p)

  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)

  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  cat(sprintf("脚本目录: %s\n", script_dir))
  cat(sprintf("工作目录: %s\n", getwd()))
  cat(sprintf("输出目录: %s\n", normalizePath(OUT_DIR)))

  tryCatch({
    w  <- system("wmic OS get FreePhysicalMemory /value", intern = TRUE)
    fr <- as.numeric(sub(".*=", "", w[grepl("FreePhysicalMemory", w)]))
    if (length(fr) == 1 && !is.na(fr)) {
      free_gb <- fr / 1024^2
      cat(sprintf("当前可用物理内存: %.1f GB\n", free_gb))
      if (free_gb < 8) {
        cat("提示: 可用内存偏少，本步大对象约需 5~8GB，建议关闭其他程序后再运行。\n")
      }
    }
  }, error = function(e) invisible(NULL))

  options(future.globals.maxSize = 8 * 1024^3)

  ## ---- 2. 读入对象并检查 CD45 基因 ----
  if (!file.exists(INPUT_RDS)) stop("找不到输入 rds: ", INPUT_RDS)
  cat("正在读取 02_seurat_qc.rds（约需 1~3 分钟）...\n")
  obj <- readRDS(INPUT_RDS)
  cat(sprintf("读入完成：%d 个基因 x %d 个细胞\n", nrow(obj), ncol(obj)))

  if (DefaultAssay(obj) != "RNA") {
    warning("默认 assay 不是 RNA，已切换为 RNA")
    DefaultAssay(obj) <- "RNA"
  }

  if (!"sample" %in% colnames(obj@meta.data)) {
    obj$sample <- as.character(obj$orig.ident)
  }

  cat(sprintf("当前 RNA layer: %s\n", paste(Layers(obj, assay = "RNA"), collapse = ", ")))
  if (length(Layers(obj, assay = "RNA")) > 1) {
    cat("检测到多个 layer，执行 JoinLayers 合并...\n")
    obj <- JoinLayers(obj, assay = "RNA")
  }

  if (!CD45_GENE %in% rownames(obj)) {
    hit <- grep("PTPRC", rownames(obj), value = TRUE, ignore.case = TRUE)
    if (length(hit) > 0) {
      CD45_GENE <- hit[1]
      cat(sprintf("提示：基因名大小写与配置不同，已改用 '%s'\n", CD45_GENE))
    } else {
      stop("数据中找不到 CD45 标志基因 PTPRC，无法提取 CD45+ 细胞")
    }
  }
  cat(sprintf("CD45 标志基因: %s\n", CD45_GENE))

  ## ---- 3. 判定 CD45+ 细胞 ----
  counts_mat <- LayerData(obj, assay = "RNA", layer = "counts")
  ptprc_counts <- as.numeric(counts_mat[CD45_GENE, ])
  rm(counts_mat); gc()

  obj$CD45_counts   <- ptprc_counts
  obj$CD45_positive <- ptprc_counts > CD45_MIN_COUNTS

  cat("\n==== CD45+ 判定结果 ====\n")
  cat(sprintf("全部细胞: %d\n", ncol(obj)))
  cat(sprintf("CD45+ 细胞: %d（%.2f%%）\n",
              sum(obj$CD45_positive), 100 * mean(obj$CD45_positive)))
  cat(sprintf("CD45- 细胞: %d\n", sum(!obj$CD45_positive)))

  cat("PTPRC 表达量分布（仅对 CD45+ 细胞）：\n")
  cat(sprintf("  counts 分位数 [5%%, 25%%, 50%%, 75%%, 95%%]: %s\n",
              paste(quantile(ptprc_counts[obj$CD45_positive], c(0.05, 0.25, 0.5, 0.75, 0.95)), collapse = ", ")))
  cat(sprintf("  counts==1 的 CD45+ 细胞占比: %.2f%%\n",
              100 * mean(ptprc_counts[obj$CD45_positive] == 1)))

  grp_summary <- obj@meta.data %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      cells_total   = dplyr::n(),
      cells_CD45pos = sum(CD45_positive),
      pct_CD45pos   = 100 * mean(CD45_positive),
      .groups       = "drop"
    )
  cat("按分组汇总：\n")
  print(as.data.frame(grp_summary))

  ## ---- 4. 按样本统计 CD45+ 细胞数并写 CSV ----
  cd45_summary <- obj@meta.data %>%
    dplyr::group_by(sample, group, sample_id) %>%
    dplyr::summarise(
      cells_total   = dplyr::n(),
      cells_CD45pos = sum(CD45_positive),
      pct_CD45pos   = 100 * mean(CD45_positive),
      .groups       = "drop"
    )
  write.csv(cd45_summary, file.path(OUT_DIR, "CD45_summary_by_sample.csv"), row.names = FALSE)
  cat("  已保存: CD45_summary_by_sample.csv\n")

  ## ---- 5. QC 图（提取前：全部细胞） ----
  cat("\n==== 5. 绘制 QC 图（提取前：全部细胞）====\n")
  p_vln_before <- VlnPlot(obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"),
    group.by = "sample", pt.size = 0, ncol = 5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_pdf(p_vln_before, "01_QC_violin_before_filter.pdf", 18, 5)

  p_scatter <- FeatureScatter(obj,
    feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample")
  save_pdf(p_scatter, "02_QC_scatter_before_filter.pdf", 8, 6)

  ## ---- 6. 提取 CD45+ 细胞并保存 ★（关键断点）----
  # 6.1 按 CD45_positive 标记提取 CD45+ 细胞
  obj_cd45 <- subset(obj, subset = CD45_positive == TRUE)
  cat(sprintf("\n提取完成：%d 个细胞 -> %d 个 CD45+ 细胞（占比 %.2f%%）\n",
              ncol(obj), ncol(obj_cd45), 100 * ncol(obj_cd45) / ncol(obj)))

  # 6.2 保存 CD45+ 细胞对象（meta.data 完整保留 sample/group/sample_id/CD45_counts/CD45_positive）
  #     ★ 此为整条合并流程的关键存盘点：后续 SECTION 04 无论从头还是 RESUME_FROM="04"
  #       都会读取本 rds，可作为分段运行的中间产物。
  saveRDS(obj_cd45, file.path(OUT_DIR, "03_CD45_positive.rds"))
  cat(sprintf("  已保存: %s\n", file.path(OUT_DIR, "03_CD45_positive.rds")))

  ## ---- 7. QC 图（提取后：CD45+ 细胞） ----
  cat("\n==== 7. 绘制 QC 图（提取后：CD45+ 细胞）====\n")
  p_vln_after <- VlnPlot(obj_cd45,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = "sample", pt.size = 0, ncol = 3) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_pdf(p_vln_after, "03_QC_violin_after_filter.pdf", 14, 5)

  sample_id_levels <- unique(cd45_summary$sample_id)
  cell_count_long <- rbind(
    data.frame(sample_id = cd45_summary$sample_id, group = cd45_summary$group,
               status = "before", count = cd45_summary$cells_total, stringsAsFactors = FALSE),
    data.frame(sample_id = cd45_summary$sample_id, group = cd45_summary$group,
               status = "after",  count = cd45_summary$cells_CD45pos, stringsAsFactors = FALSE)
  )
  cell_count_long$sample_id <- factor(cell_count_long$sample_id, levels = sample_id_levels)
  p_bar <- ggplot(cell_count_long, aes(x = sample_id, y = count, fill = status)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = c(before = "grey60", after = "#E64B35")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Cell count before/after CD45+ extraction by sample",
         x = "Sample ID", y = "Cell count", fill = "Status")
  save_pdf(p_bar, "04_QC_cell_count_by_sample.pdf", 14, 6)

  ## ---- 8. 补充图 ----
  cat("\n==== 8. 绘制补充图 ====\n")
  p_ptprc <- ggplot(obj@meta.data, aes(x = log10(CD45_counts + 1))) +
    geom_histogram(bins = 60, fill = "steelblue", color = "white") +
    geom_vline(xintercept = log10(CD45_MIN_COUNTS + 1),
               linetype = "dashed", color = "red", linewidth = 0.8) +
    annotate("text", x = log10(CD45_MIN_COUNTS + 1) + 0.15,
             y = Inf, label = paste0("CD45+ threshold (counts>", CD45_MIN_COUNTS, ")"),
             vjust = 1.5, hjust = 0, color = "red", size = 3.5) +
    theme_minimal() +
    labs(title = "PTPRC (CD45) expression distribution",
         x = "log10(PTPRC counts + 1)", y = "Cell count")
  save_pdf(p_ptprc, "05_PTPRC_expression_distribution.pdf", 9, 5)

  p_prop <- ggplot(cd45_summary, aes(x = sample_id, y = pct_CD45pos, fill = group)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c(ypN0 = "#4DBBD5", `ypN+` = "#E64B35")) +  # ★ ypN0=蓝, ypN+=红（保留原配色）
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "CD45+ cell proportion by sample",
         x = "Sample ID", y = "CD45+ proportion (%)", fill = "Group")
  save_pdf(p_prop, "06_CD45_proportion_by_sample.pdf", 14, 6)

  cd45_stack_long <- rbind(
    data.frame(sample_id = cd45_summary$sample_id, group = cd45_summary$group,
               status = "CD45-", count = cd45_summary$cells_total - cd45_summary$cells_CD45pos,
               stringsAsFactors = FALSE),
    data.frame(sample_id = cd45_summary$sample_id, group = cd45_summary$group,
               status = "CD45+", count = cd45_summary$cells_CD45pos,
               stringsAsFactors = FALSE)
  )
  cd45_stack_long$sample_id <- factor(cd45_stack_long$sample_id, levels = sample_id_levels)
  p_stack <- ggplot(cd45_stack_long, aes(x = sample_id, y = count, fill = status)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("CD45-" = "#90A4AE", "CD45+" = "#E64B35")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "CD45+ / CD45- cell count by sample (stacked)",
         x = "Sample ID", y = "Cell count", fill = "CD45 status")
  save_pdf(p_stack, "07_CD45_pos_neg_stack_by_sample.pdf", 14, 6)

  ## ---- 9. 完成 ----
  cat("\n==== CD45+ 提取完成 ====\n")
  cat(sprintf("全部细胞: %d\n", ncol(obj)))
  cat(sprintf("CD45+ 细胞: %d（%.2f%%）\n", ncol(obj_cd45), 100 * ncol(obj_cd45) / ncol(obj)))
  cat("\n输出文件清单（位于 output/03_extract_cd45/ 目录）：\n")
  cat("  03_CD45_positive.rds                  - CD45+ 细胞 Seurat 对象（关键断点产物）\n")
  cat("  CD45_summary_by_sample.csv            - 按样本/分组/编号的 CD45+ 统计表\n")
  cat("  01_QC_violin_before_filter.pdf        - 提取前（全部细胞）QC 小提琴图\n")
  cat("  02_QC_scatter_before_filter.pdf       - 提取前计数-基因数散点图\n")
  cat("  03_QC_violin_after_filter.pdf         - 提取后（CD45+ 细胞）QC 小提琴图\n")
  cat("  04_QC_cell_count_by_sample.pdf        - 提取前后细胞数柱状图\n")
  cat("  05_PTPRC_expression_distribution.pdf  - PTPRC 表达分布直方图（判定阈值依据）\n")
  cat("  06_CD45_proportion_by_sample.pdf      - 各样本 CD45+ 比例柱状图\n")
  cat("  07_CD45_pos_neg_stack_by_sample.pdf    - 各样本 CD45+ / CD45- 堆叠柱状图\n")
  cat("\n03 extract_cd45 完成。\n")

  if (STOP_AFTER %in% c("03")) {
    cat(sprintf("\n*** STOP_AFTER=%s：已在 SECTION 03 后停止（已保存 03_CD45_positive.rds）。***\n", STOP_AFTER))
    quit(save = "no")
  }
}


## ====================================================================
## SECTION 04 — INTEGRATION & CLUSTERING（标准化 + 去批次 + 降维 + 聚类）
##   读取 03_CD45_positive.rds -> Normalize/Scale/PCA/Harmony/UMAP/tSNE ->
##   多分辨率聚类 -> 保存 04_CD45_integrated.rds + 图/CSV
## ====================================================================
if (run_from("04")) {
  cat(sprintf("\n########## SECTION 04 — INTEGRATION & CLUSTERING ##########\n"))

  ## ---- 0. 配置（按需修改） ----
  INPUT_RDS <- RDS_03    # 第三步 CD45+ 细胞对象
  OUT_DIR   <- OUT_04    # 本步输出目录

  # ---- 标准化参数 ----
  NFEATURES     <- 3000
  REGRESS_MT    <- TRUE

  # ---- 降维参数 ----
  NPCS          <- 50
  DIMS          <- NULL      # NULL=自动选择；也可手动指定如 1:30
  CUMVAR_THRESHOLD <- 0.85
  STDEV_THRESHOLD  <- 1
  MIN_PCS          <- 5
  MAX_PCS          <- 40

  # ---- UMAP / tSNE 参数 ----
  UMAP_METHOD     <- "uwot"
  UMAP_METRIC     <- "cosine"
  UMAP_NEIGHBORS  <- 25
  UMAP_MIN_DIST   <- 0.2
  UMAP_SPREAD     <- 1
  TSNE_PERPLEXITY <- 30

  # ---- 去批次（Harmony）参数 ----
  RUN_HARMONY   <- TRUE
  HARMONY_BATCH <- "sample"
  HARMONY_THETA <- NULL
  PLOT_HARMONY_CONVERGENCE <- TRUE

  # ---- 聚类参数 ----
  RESOLUTIONS        <- c(0.1,0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8,0.9,1.0)
  CLUSTER_ALGORITHM  <- 1
  CLUSTER_METHOD     <- "igraph"
  GROUP_SINGLETONS   <- TRUE

  # ---- 其他 ----
  RUN_TSNE <- TRUE
  SEED     <- 123

  ## ---- 1. 依赖检查与加载 ----
  for (p in c("Seurat", "dplyr", "ggplot2", "patchwork", "harmony")) ensure_pkg(p)

  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(harmony)

  set.seed(SEED)
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  cat(sprintf("工作目录: %s\n", getwd()))
  cat(sprintf("输出目录: %s\n", normalizePath(OUT_DIR)))

  tryCatch({
    w  <- system("wmic OS get FreePhysicalMemory /value", intern = TRUE)
    fr <- as.numeric(sub(".*=", "", w[grepl("FreePhysicalMemory", w)]))
    if (length(fr) == 1 && !is.na(fr)) {
      free_gb <- fr / 1024^2
      cat(sprintf("当前可用物理内存: %.1f GB\n", free_gb))
      if (free_gb < 10) {
        cat("提示: 本步计算密集，建议可用内存 >=10GB，并关闭其他程序后再运行。\n")
      }
    }
  }, error = function(e) invisible(NULL))

  options(future.globals.maxSize = 8 * 1024^3)

  ## ---- 2. 读入对象并做基础检查 ----
  if (!file.exists(INPUT_RDS)) stop("找不到输入 rds: ", INPUT_RDS)
  cat("正在读取 03_CD45_positive.rds（约需 1~3 分钟）...\n")
  obj <- readRDS(INPUT_RDS)
  cat(sprintf("读入完成：%d 个基因 x %d 个细胞\n", nrow(obj), ncol(obj)))
  cat("分组分布：\n")
  print(table(obj$group, useNA = "ifany"))

  if (DefaultAssay(obj) != "RNA") {
    warning("默认 assay 不是 RNA，已切换为 RNA")
    DefaultAssay(obj) <- "RNA"
  }

  cat(sprintf("当前 RNA layer: %s\n", paste(Layers(obj, assay = "RNA"), collapse = ", ")))
  if (length(Layers(obj, assay = "RNA")) > 1) {
    cat("检测到多个 layer，执行 JoinLayers 合并...\n")
    obj <- JoinLayers(obj, assay = "RNA")
  }

  if (!"sample" %in% colnames(obj@meta.data)) {
    obj$sample <- as.character(obj$orig.ident)
  }
  cat(sprintf("样本数: %d\n", length(unique(obj$sample))))

  ## ---- 3. 标准化（LogNormalize） ----
  cat("\n==== 3. 标准化 ====\n")
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = NFEATURES, verbose = FALSE)
  cat(sprintf("高变基因数: %d\n", length(VariableFeatures(obj))))
  if (REGRESS_MT) {
    cat("ScaleData 回归 percent.mt...\n")
    obj <- ScaleData(obj, vars.to.regress = "percent.mt", verbose = FALSE)
  } else {
    obj <- ScaleData(obj, verbose = FALSE)
  }

  ## ---- 4. PCA 降维与 PC 数选择 ----
  cat("\n==== 4. PCA 降维与 PC 数选择 ====\n")
  obj <- RunPCA(obj, npcs = NPCS, verbose = FALSE)

  pca_stdev <- obj@reductions$pca@stdev
  var_ratio <- pca_stdev^2 / sum(pca_stdev^2)
  cumvar    <- cumsum(var_ratio)
  pc_tbl <- data.frame(
    PC       = 1:NPCS,
    stdev    = pca_stdev,
    var_ratio = var_ratio,
    cumvar   = cumvar
  )
  write.csv(pc_tbl, file.path(OUT_DIR, "01_PC_variance_table.csv"), row.names = FALSE)
  cat("  已保存: 01_PC_variance_table.csv（各 PC 方差/累积方差，供选择 PC 数）\n")

  p_elbow <- ElbowPlot(obj, ndims = NPCS) +
    geom_vline(xintercept = if (!is.null(DIMS)) max(DIMS) else NA,
               linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(title = "PCA Elbow Plot")
  save_pdf(p_elbow, "01_PCA_elbow.pdf", 8, 5)

  n_heat_pcs <- min(NPCS, 12)
  heat_max   <- if (!is.null(DIMS)) max(DIMS) else n_heat_pcs
  n_heat_pcs <- min(n_heat_pcs, heat_max)
  cat(sprintf("绘制 PC_1 ~ PC_%d 的 DimHeatmap...\n", n_heat_pcs))
  pdf(file.path(OUT_DIR, "01_PCA_dimheatmap.pdf"), width = 8.5, height = 6)
  for (i in seq_len(n_heat_pcs)) {
    p_h <- DimHeatmap(
      object      = obj,
      dims        = i,
      cells       = 500,
      balanced    = TRUE,
      nfeatures   = 30,
      reduction   = "pca",
      fast        = FALSE
    )
    print(p_h + ggtitle(paste0("PC_", i)))
  }
  dev.off()
  cat(sprintf("  已保存: 01_PCA_dimheatmap.pdf（PC_1 ~ PC_%d，每个 PC 一页）\n", n_heat_pcs))

  if (is.null(DIMS)) {
    k_cumvar <- which(cumvar >= CUMVAR_THRESHOLD)[1]
    if (is.na(k_cumvar)) k_cumvar <- NPCS
    k_stdev <- which(pca_stdev < STDEV_THRESHOLD)[1] - 1
    if (is.na(k_stdev) || k_stdev < 1) k_stdev <- NPCS
    k <- max(MIN_PCS, min(k_cumvar, k_stdev, MAX_PCS))
    DIMS <- 1:k
    cat(sprintf("自动选择 PC 数 = %d（累积方差法=%d，stdev<%.1f法=%d）\n",
                k, k_cumvar, STDEV_THRESHOLD, k_stdev))
    cat(sprintf("DIMS = 1:%d；如想手动调整，请查看 01_PCA_elbow.pdf 后修改 CONFIG 中的 DIMS。\n", k))
  } else {
    cat(sprintf("使用手动指定 DIMS = %s\n", paste(range(DIMS), collapse = ":")))
  }
  # 注意：以下一行保留原 04 脚本的硬编码覆盖（原 04_integration_clustering.R:283），
  # 会无条件将 DIMS 设为 1:30。为保证与原始流程结果一致，此处予以保留，不并入上方自动选择。
  DIMS = 1:30
  cat("前 5 个主成分的方差解释率(%):\n")
  print(round(100 * head(pc_tbl$var_ratio, 5), 2))

  ## ---- 5. 去批次（Harmony 整合） ----
  if (RUN_HARMONY) {
    cat("\n==== 5. Harmony 去批次（批次变量: sample）====\n")
    harmony_args <- list(
      object = obj,
      group.by.vars = HARMONY_BATCH,
      reduction.use = "pca",
      reduction.save = "harmony",
      # 以下 theta/sigma/max.iter.harmony 为原 04 脚本硬编码（:295-304），覆盖 HARMONY_THETA=NULL 分支，保留以保结果一致：
      theta = 5,           # 默认 2，增加到 3-5
      max.iter.harmony = 50,# 增加迭代次数
      sigma = 0.1,         # 聚类带宽
      verbose = F
    )
    if (!is.null(HARMONY_THETA)) {
      harmony_args$theta <- HARMONY_THETA
    }
    if (PLOT_HARMONY_CONVERGENCE) {
      harmony_args$plot_convergence <- TRUE
    }
    harmony_res <- do.call(RunHarmony, harmony_args)
    if (is.list(harmony_res) && !inherits(harmony_res, "Seurat")) {
      obj <- harmony_res$obj
      harmony_final <- harmony_res$final
    } else {
      obj <- harmony_res
      harmony_final <- NULL
    }
    reduction_use <- "harmony"
  } else {
    cat("\n==== 5. 跳过 Harmony（RUN_HARMONY=FALSE），使用 PCA 嵌入 ====\n")
    reduction_use <- "pca"
    harmony_final <- NULL
  }
  cat(sprintf("下游降维将使用: %s（共 %d 维中的前 %d 维）\n",
              reduction_use, NPCS, max(DIMS)))

  ## ---- 6. UMAP / tSNE 可视化降维 ----
  cat("\n==== 6. UMAP / tSNE ====\n")
  cat("运行 UMAP（约需数分钟）...\n")
  obj <- RunUMAP(obj,
    reduction = reduction_use,
    dims = DIMS,
    umap.method = UMAP_METHOD,
    metric = UMAP_METRIC,
    n.neighbors = UMAP_NEIGHBORS,
    min.dist = UMAP_MIN_DIST,
    spread = UMAP_SPREAD,
    seed.use = SEED,
    verbose = FALSE)
  cat("UMAP 完成。\n")

  if (RUN_TSNE) {
    cat("运行 tSNE（较慢，约需 5~15 分钟）...\n")
    obj <- RunTSNE(obj,
      reduction = reduction_use,
      dims = DIMS,
      perplexity = TSNE_PERPLEXITY,
      check_duplicates = FALSE,
      verbose = FALSE)
    cat("tSNE 完成。\n")
  }

  ## ---- 7. 聚类（FindNeighbors + FindClusters，多分辨率） ----
  cat("\n==== 7. 聚类（多分辨率）====\n")
  obj <- FindNeighbors(obj, reduction = reduction_use, dims = DIMS, verbose = FALSE)
  for (res in RESOLUTIONS) {
    cat(sprintf("  聚类分辨率 res=%s ...\n", res))
    obj <- FindClusters(obj,
      resolution = res,
      algorithm = CLUSTER_ALGORITHM,
      method = CLUSTER_METHOD,
      group.singletons = GROUP_SINGLETONS,
      verbose = FALSE)
    obj[[paste0("cluster_", res)]] <- as.character(obj$seurat_clusters)
    cat(sprintf("    得到 %d 个 cluster\n", length(unique(obj$seurat_clusters))))
  }
  cat("可视化顺序：先多分辨率聚类拼图（辅助选 res），后默认分辨率下的常规图。\n")

  if (length(RESOLUTIONS) > 1) {
    n_col <- min(3, length(RESOLUTIONS))
    plist <- lapply(RESOLUTIONS, function(r) {
      DimPlot(obj, reduction = "umap",
              group.by = paste0("cluster_", r),
              label = TRUE, repel = TRUE) +
        NoLegend() +
        ggtitle(paste0("res = ", r))
    })
    p_multi <- patchwork::wrap_plots(plist, ncol = n_col)
    save_pdf(p_multi, "02_UMAP_multi_resolution.pdf",
             7 * n_col, 6 * ceiling(length(RESOLUTIONS) / n_col))
  }

  ## ---- 8. 聚类统计表输出 ----
  # 原 04 脚本硬编码 DEFAULT_RESOLUTION <- 0.4（:401，注释写 0.5 实为 0.4），保留以保结果一致。
  DEFAULT_RESOLUTION <- 0.4
  obj$cluster <- as.character(obj@meta.data[[paste0("cluster_", DEFAULT_RESOLUTION)]])
  Idents(obj) <- obj$cluster

  cat(sprintf("\n默认分辨率（res=%s）聚类结果：\n", DEFAULT_RESOLUTION))
  print(table(obj$cluster))

  clust_summary <- obj@meta.data %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      pct     = 100 * dplyr::n() / ncol(obj),
      .groups = "drop"
    ) %>%
    dplyr::arrange(cluster)
  write.csv(clust_summary, file.path(OUT_DIR, "01_clustering_summary.csv"), row.names = FALSE)
  cat("  已保存: 01_clustering_summary.csv\n")

  clust_by_sample <- obj@meta.data %>%
    dplyr::group_by(cluster, sample, group, sample_id) %>%
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop")
  write.csv(clust_by_sample, file.path(OUT_DIR, "02_cluster_by_sample.csv"), row.names = FALSE)
  cat("  已保存: 02_cluster_by_sample.csv\n")

  res_summary <- data.frame(
    resolution = RESOLUTIONS,
    n_clusters = sapply(RESOLUTIONS, function(r) {
      length(unique(obj@meta.data[[paste0("cluster_", r)]]))
    }),
    row.names = NULL
  )
  write.csv(res_summary, file.path(OUT_DIR, "03_resolution_summary.csv"), row.names = FALSE)
  cat("  已保存: 03_resolution_summary.csv\n")

  ## ---- 9. 可视化 ----
  cat("\n==== 9. 可视化 ====\n")
  p_smp <- DimPlot(obj, reduction = "umap", group.by = "sample", label = FALSE) +
    theme_minimal() +
    labs(title = "UMAP by sample (after Harmony)")
  save_pdf(p_smp, "03_UMAP_by_sample.pdf", 9, 7)

  p_grp <- DimPlot(obj, reduction = "umap", group.by = "group",
                   cols = c(ypN0 = "#4DBBD5", `ypN+` = "#E64B35")) +  # ★ ypN0=蓝, ypN+=红（保留原配色）
    theme_minimal() +
    labs(title = "UMAP by group")
  save_pdf(p_grp, "04_UMAP_by_group.pdf", 8, 6)

  p_clu <- DimPlot(obj, reduction = "umap", group.by = "cluster",
                   label = TRUE, repel = TRUE) +
    theme_minimal() +
    labs(title = paste0("UMAP by cluster (res=", DEFAULT_RESOLUTION, ")"))
  save_pdf(p_clu, "05_UMAP_by_cluster.pdf", 9, 7)

  p_split <- DimPlot(obj, reduction = "umap", group.by = "sample",
                     split.by = "sample", ncol = 6) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text = element_blank(), axis.ticks = element_blank()) +
    labs(title = "UMAP by sample (split)")
  save_pdf(p_split, "06_UMAP_by_sample_split.pdf", 18, 12)

  p_bar <- ggplot(clust_summary, aes(x = cluster, y = n_cells)) +
    geom_col(fill = "#3C5488") +
    geom_text(aes(label = n_cells), vjust = -0.3, size = 3) +
    theme_minimal() +
    labs(title = paste0("Cell count per cluster (res=", DEFAULT_RESOLUTION, ")"),
         x = "Cluster", y = "Cell count")
  save_pdf(p_bar, "07_cluster_cell_count.pdf", 9, 5)

  if (RUN_TSNE && "tsne" %in% names(obj@reductions)) {
    p_tsne <- DimPlot(obj, reduction = "tsne", group.by = "cluster",
                      label = TRUE, repel = TRUE) +
      theme_minimal() +
      labs(title = paste0("tSNE by cluster (res=", DEFAULT_RESOLUTION, ")"))
    save_pdf(p_tsne, "08_TSNE_by_cluster.pdf", 9, 7)
  }

  if (RUN_HARMONY && "harmony" %in% names(obj@reductions)) {
    cat("  绘制 Harmony 去批次前后对比图（需额外运行一次 PCA-UMAP，约数分钟）...\n")
    obj <- RunUMAP(obj,
      reduction = "pca", dims = DIMS,
      reduction.name = "umap.pca", reduction.key = "PCUMAP_",
      umap.method = UMAP_METHOD, metric = UMAP_METRIC,
      n.neighbors = UMAP_NEIGHBORS, min.dist = UMAP_MIN_DIST,
      spread = UMAP_SPREAD, seed.use = SEED, verbose = FALSE)
    p_before <- DimPlot(obj, reduction = "umap.pca", group.by = "sample", label = FALSE) +
      theme_minimal() +
      labs(title = "Before Harmony (PCA-UMAP)")
    p_after <- DimPlot(obj, reduction = "umap", group.by = "sample", label = FALSE) +
      theme_minimal() +
      labs(title = "After Harmony (Harmony-UMAP)")
    p_ba <- patchwork::wrap_plots(p_before, p_after, ncol = 2)
    save_pdf(p_ba, "09_Harmony_before_after_UMAP.pdf", 16, 7)
  }

  if (PLOT_HARMONY_CONVERGENCE && !is.null(harmony_final)) {
    p_conv <- tryCatch(
      harmony::plot_convergence(harmony_final),
      error = function(e) NULL)
    if (!is.null(p_conv)) {
      p_conv <- p_conv + theme_minimal() +
        labs(title = "Harmony convergence")
      save_pdf(p_conv, "10_Harmony_convergence.pdf", 7, 5)
    } else {
      cat("  注意: Harmony 收敛图绘制失败（harmony 版本差异），已跳过，不影响主流程。\n")
    }
  }

  ## ---- 10. 保存整合对象与完成 ----
  saveRDS(obj, file.path(OUT_DIR, "04_CD45_integrated.rds"))
  cat(sprintf("  已保存: %s\n", file.path(OUT_DIR, "04_CD45_integrated.rds")))

  cat("\n==== 04 整合聚类完成 ====\n")
  cat(sprintf("细胞数: %d\n", ncol(obj)))
  cat(sprintf("默认分辨率 res=%s 的簇数: %d\n",
              DEFAULT_RESOLUTION, length(unique(obj$cluster))))
  cat(sprintf("本次使用的 PC 数（DIMS）: %d\n", max(DIMS)))
  cat("\n输出文件清单（位于 output/04_integration_clustering/ 目录）：\n")
  cat("  04_CD45_integrated.rds           - 整合聚类后的 Seurat 对象\n")
  cat("  01_PC_variance_table.csv         - 各 PC 方差/累积方差表（PC 选择依据）\n")
  cat("  01_clustering_summary.csv        - 主 cluster 细胞数与占比\n")
  cat("  02_cluster_by_sample.csv         - cluster × 样本交叉表\n")
  cat("  03_resolution_summary.csv        - 各分辨率簇数汇总\n")
  cat("  01_PCA_elbow.pdf                 - PCA Elbow 图（PC 选择依据）\n")
  cat("  01_PCA_dimheatmap.pdf            - PC_1~PC_n 的 DimHeatmap（每个 PC 一页）\n")
  cat("  02_UMAP_multi_resolution.pdf     - 多分辨率 UMAP 拼图（先出，辅助选 res）\n")
  cat("  03_UMAP_by_sample.pdf            - UMAP 按样本着色（去批次效果）\n")
  cat("  04_UMAP_by_group.pdf             - UMAP 按分组着色\n")
  cat("  05_UMAP_by_cluster.pdf           - UMAP 按主 cluster 着色\n")
  cat("  06_UMAP_by_sample_split.pdf      - UMAP 按样本分面\n")
  cat("  07_cluster_cell_count.pdf        - 各 cluster 细胞数柱状图\n")
  cat("  08_TSNE_by_cluster.pdf           - tSNE 按主 cluster 着色（可选）\n")
  cat("  09_Harmony_before_after_UMAP.pdf - Harmony 去批次前后 UMAP 对比\n")
  cat("  10_Harmony_convergence.pdf       - Harmony 收敛诊断图（可选）\n")
  cat("\n04 integration_clustering 完成。\n")
}


## ====================================================================
## SECTION 05 — CELL TYPE ANNOTATION（CD45+ 免疫细胞类型注释）
##   读取 04_CD45_integrated.rds -> marker 打分 -> 注释（首遍生成模板/
##   再遍应用 celltype_map.csv）-> 保存 05_annotated.rds + 图/CSV
## ====================================================================
if (run_from("05")) {
  cat(sprintf("\n########## SECTION 05 — CELL TYPE ANNOTATION ##########\n"))

  ## ---- 0. 配置 ----
  INPUT_RDS <- RDS_04    # 第四步整合聚类对象
  OUT_DIR   <- OUT_05    # 本步输出目录
  MAP_PATH  <- file.path(script_dir, "celltype_map.csv")  # 用户编辑的注释映射表（已复制到本目录）

  CLUSTER_COL <- "cluster"
  SPECIES <- "human"
  SEED <- 123
  SCORE_ASSAY <- NULL

  # 差异基因（FindAllMarkers）参数
  FM_ONLY_POS   <- TRUE
  FM_MIN_PCT    <- 0.25
  FM_LOGFC      <- 0.25
  FM_TOP_N      <- 5
  ANNOT_MARKER_CAP <- 5

  # SeuratExtend 绘图参数
  FP_PT_SIZE  <- 0.5
  FP_COLOR    <- "ryb"
  FP_NCOL     <- 4
  DOT_COLOR_SCHEME <- "Reds"
  DOT_BORDER       <- FALSE
  DOT_SHOW_GRID    <- FALSE
  DOT_FLIP         <- TRUE
  TYPE_ORDER <- c("CD8T", "CD4T", "T cell", "T Regulatory", "T Exhausted",
                 "NK cell", "B cell", "Plasma cell", "Monocyte",
                 "Macrophage", "DC", "Mast cell", "Neutrophil")
  DOT_W_PER_GENE <- 0.45
  DOT_H_PER_ROW  <- 0.5
  DOT_BASE_H     <- 2
  HM_SCALE    <- "row"
  HM_COLOR    <- c("#3C5488", "white", "#E64B35")   # ★ 蓝-白-红（保留原配色）
  COMP_STACK  <- TRUE

  SE_REPO     <- "huayc09/SeuratExtend"
  SE_UPGRADE  <- "never"
  VERBOSE     <- FALSE

  # 内联经典 marker 字典（原 manual_markers.R 内容，供 Part 2 点图使用）
  MANUAL_MARKERS <- list(
    "CD8T"         = c("CD8A", "CD8B", "GZMK", "GZMA", "CCL5"),
    "CD4T"         = c("CD4", "IL7R", "CD40LG", "CCR7", "ICOS"),
    "T cell"       = c("CD3D", "CD3E", "CD3G", "TRAC", "CD2"),
    "T Regulatory" = c("FOXP3", "IL2RA", "IKZF2", "CTLA4", "TIGIT"),
    "T Exhausted"  = c("PDCD1", "TIGIT", "LAG3", "HAVCR2", "CXCL13"),
    "NK cell"      = c("KLRD1", "GNLY", "NKG7", "PRF1", "FCGR3A"),
    "B cell"       = c("CD19", "MS4A1", "CD79A", "CD79B", "IGHD"),
    "Plasma cell"  = c("MZB1", "SDC1", "JCHAIN", "PRDM1", "IGHG1"),
    "Monocyte"     = c("CD14", "LYZ", "FCN1", "CSF1R", "FCGR3A"),
    "Macrophage"   = c("CD68", "CD163", "C1QA", "C1QB", "CST3"),
    "DC"           = c("CLEC9A", "CADM1", "FCER1A", "CD1C", "CLEC10A"),
    "Mast cell"    = c("TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"),
    "Neutrophil"   = c("S100A8", "S100A9", "FCGR3B", "MPO", "ELANE")
  )

  ## ---- 1. 依赖检查与加载 ----
  for (p in c("Seurat", "dplyr", "ggplot2", "patchwork", "pheatmap", "remotes")) ensure_pkg(p)
  if (!ensure_seuratextend(SE_REPO, SE_UPGRADE)) {
    stop("SeuratExtend 安装/加载失败，无法继续（请确认网络可访问 GitHub 且已安装 Rtools）。")
  }
  library(Seurat)
  library(SeuratExtend)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)

  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  cat(sprintf("工作目录: %s\n", getwd()))
  cat(sprintf("输出目录: %s\n", normalizePath(OUT_DIR)))

  tryCatch({
    w  <- system("wmic OS get FreePhysicalMemory /value", intern = TRUE)
    fr <- as.numeric(sub(".*=", "", w[grepl("FreePhysicalMemory", w)]))
    if (length(fr) == 1 && !is.na(fr)) {
      free_gb <- fr / 1024^2
      cat(sprintf("当前可用物理内存: %.1f GB\n", free_gb))
      if (free_gb < 8) cat("提示: 本步建议可用内存 >=8GB。\n")
    }
  }, error = function(e) invisible(NULL))

  options(future.globals.maxSize = 8 * 1024^3)

  ## ---- 1b. 辅助函数 ----
  plot_score_heatmap <- function(mat_hm, OUT_DIR) {
    pdf(file.path(OUT_DIR, "02_score_heatmap.pdf"),
        width = 9, height = max(5, nrow(mat_hm) * 0.35 + 2))
    tryCatch({
      SeuratExtend::Heatmap(mat_hm,
        cluster_rows = FALSE, cluster_columns = TRUE,
        show_row_names = TRUE, show_column_names = TRUE)
    }, error = function(e) {
      cat("  SeuratExtend::Heatmap 失败，回退 pheatmap:", conditionMessage(e), "\n")
      pheatmap::pheatmap(mat_hm,
        scale = HM_SCALE, cluster_cols = TRUE, cluster_rows = FALSE,
        color = colorRampPalette(HM_COLOR)(100),
        main = "Average marker score by cluster", fontsize = 10)
    })
    dev.off()
    cat("  已保存图表: 02_score_heatmap.pdf\n")
  }

  plot_marker_umap <- function(obj, feat_genes, OUT_DIR) {
    p_feat <- tryCatch(
      FeaturePlot3.grid(obj,
        features = feat_genes, pt.size = FP_PT_SIZE,
        color = FP_COLOR, ncol = FP_NCOL),
      error = function(e) {
        cat("  FeaturePlot3.grid 不支持 color/ncol 参数，回退最小参数集:", conditionMessage(e), "\n")
        FeaturePlot3.grid(obj, features = feat_genes, pt.size = FP_PT_SIZE)
      })
    ggsave(file.path(OUT_DIR, "03_marker_UMAP.pdf"),
           p_feat, width = 5.5 * FP_NCOL,
           height = 5.5 * ceiling(length(feat_genes) / FP_NCOL), limitsize = FALSE)
    cat("  已保存图表: 03_marker_UMAP.pdf\n")
  }

  apply_celltype_map <- function(obj, MAP_PATH, CLUSTER_COL) {
    cat("检测到注释映射表 celltype_map.csv，正在应用...\n")
    map <- read.csv(MAP_PATH, stringsAsFactors = FALSE)
    if (!all(c("cluster", "celltype") %in% colnames(map))) {
      stop("celltype_map.csv 需要包含 cluster 和 celltype 两列")
    }
    map$celltype <- trimws(as.character(map$celltype))
    if ("top_hit" %in% colnames(map)) {
      empty_idx <- is.na(map$celltype) | map$celltype == ""
      if (sum(empty_idx) > 0) {
        map$celltype[empty_idx] <- map$top_hit[empty_idx]
        cat(sprintf("提示: %d 个 cluster 的 celltype 为空，已自动采用 top_hit 建议作为初始注释。\n", sum(empty_idx)))
        cat("      如需精调，请编辑 celltype_map.csv 中的 celltype 列后重新运行本脚本。\n")
      }
    }
    obj$celltype <- "Unknown"
    for (i in seq_len(nrow(map))) {
      if (is.na(map$celltype[i]) || map$celltype[i] == "") next
      obj$celltype[obj[[CLUSTER_COL]] == as.character(map$cluster[i])] <- map$celltype[i]
    }
    cat("注释结果（细胞类型分布）：\n")
    print(table(obj$celltype, useNA = "ifany"))
    unknown_n <- sum(obj$celltype == "Unknown")
    if (unknown_n > 0) cat(sprintf("注意: %d 个细胞未注释（Unknown），请在 celltype_map.csv 中补全后重跑\n", unknown_n))
    obj
  }

  write_map_template <- function(map_template, OUT_DIR, script_dir) {
    cat("尚未找到 celltype_map.csv（首次运行）。\n")
    cat("请先查看 output/02_score_heatmap.pdf 与 output/03_marker_UMAP.pdf，\n")
    cat("确认每个 cluster 的身份后，编辑 output/celltype_map_template.csv 中的 celltype 列，\n")
    cat("保存为 celltype_map.csv（与本脚本同目录），再重新运行本脚本即可应用注释。\n")
    write.csv(map_template, file.path(OUT_DIR, "celltype_map_template.csv"), row.names = FALSE, na = "")
    write.csv(map_template, file.path(script_dir, "celltype_map.csv"), row.names = FALSE, na = "")
    cat("  已生成注释模板: celltype_map.csv（脚本目录）与 celltype_map_template.csv（output/）\n")
  }

  plot_annotation_umap <- function(obj, OUT_DIR) {
    p_anno <- DimPlot2(obj, features = "celltype", label = TRUE, box = TRUE,
                       repel = TRUE, theme = NoLegend())
    ggsave(file.path(OUT_DIR, "04_annotation_UMAP.pdf"), p_anno, width = 9, height = 7)
    cat("  已保存图表: 04_annotation_UMAP.pdf\n")
    p_anno_grp <- DimPlot2(obj, features = "celltype", split.by = "group",
                           label = TRUE, box = TRUE, repel = TRUE)
    ggsave(file.path(OUT_DIR, "04_annotation_UMAP_by_group.pdf"), p_anno_grp, width = 16, height = 7)
    cat("  已保存图表: 04_annotation_UMAP_by_group.pdf\n")
  }

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

  plot_composition <- function(obj, OUT_DIR) {
    p_comp <- tryCatch({
      ClusterDistrBar(origin = obj$group, cluster = obj$celltype,
                      flip = TRUE, stack = COMP_STACK)
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

  derive_celltype_detailed <- function(obj, CLUSTER_COL) {
    ct <- as.character(obj$celltype)
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
    if (all(c("Score_CD8T1", "Score_CD4T1") %in% colnames(obj@meta.data))) {
      s8 <- obj$Score_CD8T1; s4 <- obj$Score_CD4T1
      isT <- ct %in% c("T cell", "Tcell", "T")
      ct[isT] <- ifelse(s8[isT] >= s4[isT] & pmax(s8[isT], s4[isT]) > 0,
                        "CD8T", ifelse(pmax(s8[isT], s4[isT]) > 0, "CD4T", "T cell"))
    }
    present <- intersect(TYPE_ORDER, unique(ct))
    obj$celltype_detailed <- factor(ct, levels = present)
    obj
  }

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

  get_manual_markers <- function(obj) {
    if (!exists("MANUAL_MARKERS") || length(MANUAL_MARKERS) == 0) return(list())
    lst <- MANUAL_MARKERS[intersect(TYPE_ORDER, names(MANUAL_MARKERS))]
    lapply(lst, function(g) g[g %in% rownames(obj)])
  }

  combine_auto_manual <- function(auto_lst, manual_lst) {
    common <- intersect(names(auto_lst), names(manual_lst))
    out <- list()
    for (ct in common) out[[ct]] <- c(auto_lst[[ct]], manual_lst[[ct]])
    out[intersect(TYPE_ORDER, names(out))]
  }

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

  run_three_dotplots <- function(obj, OUT_DIR) {
    auto_lst   <- get_auto_markers(obj, "celltype_detailed")
    manual_lst <- get_manual_markers(obj)
    plot_dot_generic(obj, auto_lst,                       "06_marker_dotplot_auto.pdf",     OUT_DIR)
    plot_dot_generic(obj, manual_lst,                     "06_marker_dotplot_manual.pdf",   OUT_DIR)
    plot_dot_generic(obj, combine_auto_manual(auto_lst, manual_lst),
                                                     "06_marker_dotplot_combined.pdf", OUT_DIR)
  }

  ## ---- 2. 读入对象与基础检查 ----
  if (!file.exists(INPUT_RDS)) stop("找不到输入 rds: ", INPUT_RDS)
  cat("正在读取 04_CD45_integrated.rds（约需 1~3 分钟）...\n")
  obj <- readRDS(INPUT_RDS)
  cat(sprintf("读入完成：%d 个基因 x %d 个细胞\n", nrow(obj), ncol(obj)))

  if (!CLUSTER_COL %in% colnames(obj@meta.data)) {
    stop("对象中找不到主 cluster 列 '", CLUSTER_COL, "'，请先运行 04_integration_clustering.R")
  }
  cat(sprintf("当前聚类（%s）簇数: %d\n", CLUSTER_COL, length(unique(obj[[CLUSTER_COL]]))))

  if (!"data" %in% Layers(obj, assay = DefaultAssay(obj))) {
    cat("对象缺少 data 层，执行 NormalizeData（LogNormalize）...\n")
    obj <- NormalizeData(obj, verbose = FALSE)
  }
  if (!"group" %in% colnames(obj@meta.data)) {
    cat("提示: 对象中无 'group' 列，组成图将退化为单一分组 'All'。\n")
    obj$group <- "All"
  }

  ## ---- 3. 细胞注释信息整理 ----
  CELLTYPE_MARKERS <- list(
    "T cell"      = c("CD3D", "CD3E", "CD4", "CD8A", "GZMK", "IL7R"),
    "B cell"      = c("CD19", "CD79A", "CD79B", "MS4A1"),
    "Plasma cell" = c("IGHG1", "MZB1", "SDC1", "CD79A", "PRDM1", "JCHAIN"),
    "NK cell"     = c("KLRD1", "GNLY", "NKG7", "CD247", "FCGR3A", "PRF1"),
    "Macrophage"  = c("LYZ", "CST3", "CD14", "CD68", "CD163", "C1QA"),
    "Neutrophil"  = c("CST3", "LYZ", "FCGR3B", "CSF3R", "MPO", "S100A8", "S100A9"),
    "Mast cell"   = c("TPSAB1", "TPSB2", "CPA3", "KIT", "CST3")
  )
  cat(sprintf("已整理 %d 种细胞类型的注释 marker 信息。\n", length(CELLTYPE_MARKERS)))

  MARKER_SETS <- list(
    Tcell    = c("CD3D", "CD3E", "CD3G"),
    CD8T     = c("CD8A", "CD8B"),
    CD4T     = c("CD4", "IL7R"),
    NK       = c("NKG7", "GNLY", "KLRD1", "NCAM1"),
    B        = c("MS4A1", "CD79A", "CD79B"),
    Plasma   = c("MZB1", "SDC1", "JCHAIN"),
    Mono     = c("LYZ", "FCGR3A", "CSF1R"),
    Macro    = c("CD68", "C1QA", "C1QB", "CSF1R"),
    DC       = c("CLEC9A", "ITGAX", "FCER1A", "BATF3"),
    Mast     = c("TPSAB1", "TPSB2", "CPA3"),
    Neutro   = c("FCGR3B", "S100A8", "S100A9")
  )

  FEATURE_GENES <- c("CD3D", "CD8A", "CD4", "NKG7", "MS4A1", "MZB1",
                     "LYZ", "CD68", "CLEC9A", "TPSAB1", "FCGR3B")

  ## ---- 4. 免疫 marker 打分 ----
  cat("\n==== 4. 免疫 marker 打分 ====\n")
  used_sets <- list()
  for (ct in names(MARKER_SETS)) {
    genes <- MARKER_SETS[[ct]]
    present <- genes[genes %in% rownames(obj)]
    if (length(present) == 0) {
      cat(sprintf("  警告: %s 的 marker 基因全部缺失，跳过该类型打分\n", ct))
      next
    }
    if (length(present) < length(genes)) {
      cat(sprintf("  %s: 缺失 %d 个基因（%s），使用剩余 %d 个\n",
                  ct, length(genes) - length(present),
                  paste(setdiff(genes, present), collapse = "/"), length(present)))
    }
    obj <- AddModuleScore(obj,
      features = list(present),
      name = paste0("Score_", ct),
      seed = SEED,
      assay = if (is.null(SCORE_ASSAY)) DefaultAssay(obj) else SCORE_ASSAY)
    used_sets[[ct]] <- present
    cat(sprintf("  %s 打分完成（%d 个基因）\n", ct, length(present)))
  }

  score_cols <- paste0("Score_", names(used_sets), "1")
  cat("打分列:", paste(score_cols, collapse = ", "), "\n")

  score_mat <- obj@meta.data %>%
    dplyr::group_by(.data[[CLUSTER_COL]]) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(score_cols), mean),
      n_cells = dplyr::n(),
      .groups = "drop"
    )
  colnames(score_mat)[2:(length(score_cols) + 1)] <- names(used_sets)

  write.csv(score_mat, file.path(OUT_DIR, "01_celltype_score_by_cluster.csv"), row.names = FALSE)
  cat("  已保存: 01_celltype_score_by_cluster.csv（cluster × 细胞类型平均得分）\n")

  mat_hm <- as.matrix(score_mat[, names(used_sets), drop = FALSE])
  rownames(mat_hm) <- score_mat[[CLUSTER_COL]]
  plot_score_heatmap(mat_hm, OUT_DIR)

  top_hit   <- apply(mat_hm, 1, function(x) names(used_sets)[which.max(x)])
  top_score <- apply(mat_hm, 1, max)

  ## ---- 5. marker UMAP 表达图 ----
  cat(sprintf("\n==== 5. marker UMAP 表达图 ====\n"))
  feat_genes <- FEATURE_GENES[FEATURE_GENES %in% rownames(obj)]
  cat(sprintf("  使用 %d 个代表性 marker 基因作图\n", length(feat_genes)))
  if (length(feat_genes) > 0) {
    plot_marker_umap(obj, feat_genes, OUT_DIR)
  }

  ## ---- 6. 应用注释 ----
  map_template <- data.frame(
    cluster  = score_mat[[CLUSTER_COL]],
    n_cells  = score_mat$n_cells,
    top_hit  = top_hit,
    celltype = NA_character_,
    stringsAsFactors = FALSE
  )
  cat("\n==== 6. 注释流程 ====\n")
  if (file.exists(MAP_PATH)) {
    obj <- apply_celltype_map(obj, MAP_PATH, CLUSTER_COL)
  } else {
    write_map_template(map_template, OUT_DIR, script_dir)
  }

  if (!"celltype" %in% colnames(obj@meta.data)) {
    obj$celltype <- top_hit[as.character(obj[[CLUSTER_COL]])]
    obj$celltype[is.na(obj$celltype)] <- "Unknown"
  }
  obj <- derive_celltype_detailed(obj, CLUSTER_COL)
  cat("celltype_detailed 分布:\n")
  print(table(obj$celltype_detailed, useNA = "ifany"))

  ## ---- 7. 注释后可视化与统计 ----
  cat("\n==== 7. 注释后可视化与统计 ====\n")
  plot_annotation_umap(obj, OUT_DIR)
  write_composition_stats(obj, OUT_DIR)
  plot_composition(obj, OUT_DIR)

  cat("\n==== 7d. 拆分三份 marker 点图（自动 / 手动 / 自动+经典）====\n")
  run_three_dotplots(obj, OUT_DIR)

  ## ---- 8. 保存注释对象与完成 ----
  saveRDS(obj, file.path(OUT_DIR, "05_annotated.rds"))
  cat(sprintf("  已保存: %s\n", file.path(OUT_DIR, "05_annotated.rds")))

  annotated <- "celltype" %in% colnames(obj@meta.data) && sum(obj$celltype != "Unknown", na.rm = TRUE) > 0
  cat("\n==== 05 细胞注释完成 ====\n")
  if (annotated) {
    cat(sprintf("已注释 %d 个细胞、%d 种细胞类型\n",
                ncol(obj), length(unique(obj$celltype[obj$celltype != "Unknown"]))))
  } else {
    cat("当前为首次运行：请查看打分热图与 marker UMAP 后填写 celltype_map.csv 并重跑。\n")
  }
  cat("\n输出文件清单（位于 output/05_celltype_annotation/ 目录）：\n")
  cat("  05_annotated.rds                   - 注释后的 Seurat 对象（含打分与 celltype 列）\n")
  cat("  01_celltype_score_by_cluster.csv   - cluster × 细胞类型平均得分矩阵\n")
  cat("  02_score_heatmap.pdf               - 打分热图（注释主要依据）\n")
  cat("  03_marker_UMAP.pdf                 - 关键 marker 的 UMAP 表达图\n")
  cat("  04_annotation_UMAP.pdf             - 注释后的 UMAP\n")
  cat("  04_annotation_UMAP_by_group.pdf    - 按分组的注释 UMAP\n")
  cat("  05_celltype_composition.pdf        - 细胞类型组成图\n")
  cat("  06_marker_dotplot_auto.pdf        - 自动聚类 top5 marker 点图\n")
  cat("  06_marker_dotplot_manual.pdf      - 手动指定经典 marker 点图\n")
  cat("  06_marker_dotplot_combined.pdf    - 自动+经典 marker 组合点图\n")
  cat("  02_celltype_summary.csv / 03_celltype_by_group.csv - 组成统计表\n")
  cat("  celltype_map_template.csv          - 注释映射模板（cluster → celltype 待填）\n")
  cat("\n05 celltype_annotation 完成。\n")
}


## ====================================================================
## SECTION 06 — PARAMETER SWEEP（DIMS × RESOLUTIONS × UMAP 参数扫描）
##   读取 03_CD45_positive.rds（GSE203115 CD45+ 子集）->
##   标准化/变量基因/缩放/PCA/Harmony(theta=5) ->
##   对每个 DIMS 构建邻域图 -> 多分辨率聚类 -> 多 UMAP 嵌入 ->
##   输出：每个 (DIMS, n.neighbors) 一个 PDF（6 页，每页一个 min.dist；
##         每页 8 分辨率网格：上行 cluster 着色 / 下行 group 着色）
##         + 1 份合成 sweep_combined.pdf（封面 + 全部网格页）
##         + cluster_count_by_dims_resolution.csv
##   运行：Rscript merged_pipeline_01_05.R RESUME_FROM=06
## ====================================================================
if (run_from("06")) {
  cat(sprintf("\n########## SECTION 06 — PARAMETER SWEEP ##########\n"))

  suppressPackageStartupMessages({
    library(Seurat); library(dplyr); library(ggplot2)
    library(patchwork); library(cowplot); library(harmony)
  })

  ## ---- 0. 扫描参数（用户指定）----
  SWEEP_DIMS           <- c(10, 20, 30, 40, 50)
  SWEEP_RESOLUTIONS    <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
  SWEEP_UMAP_NEIGHBORS <- c(10, 25, 50, 100, 150)
  SWEEP_UMAP_MIN_DIST  <- c(0.05, 0.1, 0.2, 0.4, 0.9, 1.5)

  ## ---- 固定参数（沿用 SECTION 04 硬编码，保证可比）----
  SWEEP_NFEATURES      <- 3000
  SWEEP_REGRESS_MT     <- TRUE
  SWEEP_NPCS           <- 50
  SWEEP_HARMONY_THETA  <- 5
  SWEEP_HARMONY_MAXITER<- 50
  SWEEP_HARMONY_SIGMA  <- 0.1
  SWEEP_UMAP_METHOD    <- "uwot"
  SWEEP_UMAP_METRIC    <- "cosine"
  SWEEP_UMAP_SPREAD    <- 1
  SWEEP_SEED           <- 123

  INPUT_RDS <- RDS_03
  OUT_DIR   <- file.path(script_dir, "sweep_output")
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

  if (!file.exists(INPUT_RDS)) {
    stop("找不到 ", INPUT_RDS, "，请先运行 STOP_AFTER=03 生成（或准备 GSE203115 子集）。")
  }
  obj <- readRDS(INPUT_RDS)
  cat(sprintf("读入: %d 基因 x %d 细胞\n", nrow(obj), ncol(obj)))
  cat("分组分布:\n"); print(table(obj$group, useNA = "ifany"))

  if (DefaultAssay(obj) != "RNA") DefaultAssay(obj) <- "RNA"
  if (length(Layers(obj, assay = "RNA")) > 1) obj <- JoinLayers(obj, assay = "RNA")
  if (!"sample" %in% colnames(obj@meta.data)) obj$sample <- as.character(obj$orig.ident)
  grp_col <- if ("group" %in% colnames(obj@meta.data)) "group" else "orig.ident"

  cat("标准化 / 变量基因 / 缩放 ...\n")
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = SWEEP_NFEATURES, verbose = FALSE)
  if (SWEEP_REGRESS_MT) {
    obj <- ScaleData(obj, vars.to.regress = "percent.mt", verbose = FALSE)
  } else {
    obj <- ScaleData(obj, verbose = FALSE)
  }
  cat("PCA ...\n")
  obj <- RunPCA(obj, npcs = SWEEP_NPCS, verbose = FALSE)

  if (requireNamespace("future", quietly = TRUE)) future::plan(future::sequential)
  options(future.globals.maxSize = 8 * 1024^3)
  set.seed(SWEEP_SEED)

  total_combos <- length(SWEEP_DIMS) * length(SWEEP_RESOLUTIONS) *
                  length(SWEEP_UMAP_NEIGHBORS) * length(SWEEP_UMAP_MIN_DIST)
  cat(sprintf("总组合数 = %d DIMS x %d res x %d neighbors x %d min.dist = %d\n",
              length(SWEEP_DIMS), length(SWEEP_RESOLUTIONS),
              length(SWEEP_UMAP_NEIGHBORS), length(SWEEP_UMAP_MIN_DIST), total_combos))

  # 流式生成（关键修复）：边计算边出图，绝不一次性把所有网格驻留内存。
  # 原先把 150 页 wrap_plots（每页 16 个 DimPlot = 共 2400 个 ggplot 对象）全存进 all_pages，
  # 再分别写入独立 PDF 与合成 PDF，峰值内存过高导致 R 段错误崩溃（日志块缓冲只刷到首个 UMAP）。
  # 现改为：合成 PDF 设备常开，对每个 (D,n,d) 实时 UMAP+构图，立即双写（独立 PDF + 合成 PDF），
  # 随后释放该 UMAP 嵌入，内存保持有界。
  summ_rows <- list()
  n_ind_pdf <- 0

  comb_path <- file.path(OUT_DIR, "sweep_combined.pdf")
  cat(sprintf("合成 %s ...\n", comb_path))
  pdf(comb_path, width = 34, height = 13)
  comb_dev  <- dev.cur()
  p_title <- ggdraw()
  p_title <- p_title + draw_text("DIMS x Resolution x UMAP 参数扫描",
                  x = 0.5, y = 0.86, size = 22, hjust = 0.5, fontface = "bold")
  p_title <- p_title + draw_text("数据: GSE203115 CD45+ 子集 (Harmony theta=5, cosine)",
                  x = 0.5, y = 0.76, size = 12, hjust = 0.5)
  p_title <- p_title + draw_text(paste0("DIMS            = ", paste(SWEEP_DIMS, collapse = ", ")),
                  x = 0.5, y = 0.66, size = 12, hjust = 0.5)
  p_title <- p_title + draw_text(paste0("RESOLUTIONS     = ", paste(SWEEP_RESOLUTIONS, collapse = ", ")),
                  x = 0.5, y = 0.58, size = 12, hjust = 0.5)
  p_title <- p_title + draw_text(paste0("UMAP_NEIGHBORS  = ", paste(SWEEP_UMAP_NEIGHBORS, collapse = ", ")),
                  x = 0.5, y = 0.50, size = 12, hjust = 0.5)
  p_title <- p_title + draw_text(paste0("UMAP_MIN_DIST   = ", paste(SWEEP_UMAP_MIN_DIST, collapse = ", ")),
                  x = 0.5, y = 0.42, size = 12, hjust = 0.5)
  p_title <- p_title + draw_text(paste0("共 ", length(SWEEP_DIMS) * length(SWEEP_UMAP_NEIGHBORS) * length(SWEEP_UMAP_MIN_DIST),
                  " 个网格页 (DIMS x n.neighbors x min.dist)；每页上行=cluster 着色，下行=group 着色。"),
                  x = 0.5, y = 0.32, size = 12, hjust = 0.5)
  print(p_title)

  for (D in SWEEP_DIMS) {
    dims <- 1:D
    cat(sprintf("\n--- DIMS = 1:%d ---\n", D))
    cat("  Harmony 整合 (theta=5) ...\n")
    obj <- RunHarmony(obj, group.by.vars = "sample", reduction.use = "pca",
                      reduction.save = "harmony", theta = SWEEP_HARMONY_THETA,
                      max_iter = SWEEP_HARMONY_MAXITER, sigma = SWEEP_HARMONY_SIGMA,
                      verbose = FALSE)
    cat("  FindNeighbors ...\n")
    obj <- FindNeighbors(obj, reduction = "harmony", dims = dims, verbose = FALSE)

    # 多分辨率聚类（同一邻域图，开销小）
    for (r in SWEEP_RESOLUTIONS) {
      obj <- FindClusters(obj, resolution = r, algorithm = 1, method = "igraph",
                          group.singletons = TRUE, verbose = FALSE)
      obj[[paste0("cluster_", r)]] <- as.character(obj$seurat_clusters)
      summ_rows[[length(summ_rows) + 1]] <-
        data.frame(dims = D, resolution = r,
                   n_clusters = length(unique(obj@meta.data[[paste0("cluster_", r)]])))
    }

    # 流式：(n, d) 实时 UMAP + 构图 + 双写（独立 PDF 与合成 PDF）+ 释放嵌入
    for (n in SWEEP_UMAP_NEIGHBORS) {
      ind_name <- sprintf("sweep_D%d_n%d.pdf", D, n)
      pdf(file.path(OUT_DIR, ind_name), width = 34, height = 13)
      ind_dev  <- dev.cur()
      for (d in SWEEP_UMAP_MIN_DIST) {
        key <- paste0("umap_D", D, "_n", n, "_d", d)
        # UMAP（cosine 失败则回退 euclidean，避免单点失败中断整轮扫描）
        obj <- tryCatch({
          RunUMAP(obj, reduction = "harmony", dims = dims,
                  umap.method = SWEEP_UMAP_METHOD, metric = SWEEP_UMAP_METRIC,
                  n.neighbors = n, min.dist = d, spread = SWEEP_UMAP_SPREAD,
                  seed.use = SWEEP_SEED, reduction.name = key, verbose = FALSE)
        }, error = function(e) {
          cat(sprintf("  [warn] cosine UMAP 失败 (n=%d d=%s)：%s；回退 euclidean\n",
                      n, d, conditionMessage(e)))
          RunUMAP(obj, reduction = "harmony", dims = dims,
                  n.neighbors = n, min.dist = d, spread = SWEEP_UMAP_SPREAD,
                  seed.use = SWEEP_SEED, reduction.name = key, verbose = FALSE)
        })
        # 构建 8 分辨率网格（上行 cluster / 下行 group）
        panels <- list()
        for (r in SWEEP_RESOLUTIONS) {
          res_col <- paste0("cluster_", r)
          p_cl <- DimPlot(obj, reduction = key, group.by = res_col,
                          label = TRUE, repel = TRUE, pt.size = 0.4) +
                  NoLegend() + ggtitle(paste0("cluster | res=", r))
          p_gr <- DimPlot(obj, reduction = key, group.by = grp_col,
                          pt.size = 0.4) + NoLegend() + ggtitle("group")
          panels[[paste0("c_", r)]] <- p_cl
          panels[[paste0("g_", r)]] <- p_gr
        }
        ord  <- c(paste0("c_", SWEEP_RESOLUTIONS), paste0("g_", SWEEP_RESOLUTIONS))
        grid <- wrap_plots(panels[ord], nrow = 2, ncol = 8) +
                plot_annotation(title = paste0("DIMS=1:", D,
                  "  |  UMAP n.neighbors=", n, "  min.dist=", d))
        dev.set(ind_dev); print(grid)   # 写入独立 PDF
        dev.set(comb_dev); print(grid)  # 写入合成 PDF
        obj[[key]] <- NULL               # 释放嵌入，控制内存
      }
      dev.set(ind_dev); dev.off()       # 关闭独立 PDF，回到合成设备
      n_ind_pdf <- n_ind_pdf + 1
      cat(sprintf("  已保存独立 PDF: %s\n", ind_name))
    }
  }
  dev.set(comb_dev); dev.off()
  cat(sprintf("  已保存合成 PDF: %s\n", comb_path))

  ## ---- 聚类数汇总 CSV（按 DIMS x resolution）----
  summ_df <- do.call(rbind, summ_rows)
  write.csv(summ_df, file.path(OUT_DIR, "cluster_count_by_dims_resolution.csv"), row.names = FALSE)

  cat(sprintf("\n完成 SECTION 06：%d 个独立 PDF + 1 份合成 PDF (%s)\n",
              n_ind_pdf, comb_path))
  cat(sprintf("总组合数 = %d\n", total_combos))
}


## ====================================================================
## 全流程结束
## ====================================================================
cat(sprintf("\n########## 合并流程 01->06 全部完成 ##########\n"))
cat(sprintf("各段关键产物：\n"))
cat(sprintf("  %s\n", RDS_01))
cat(sprintf("  %s\n", RDS_02))
cat(sprintf("  %s  (★ 关键断点)\n", RDS_03))
cat(sprintf("  %s\n", RDS_04))
cat(sprintf("  %s\n", RDS_05))
