# =====================================================================
# 01_initiation.R
# 单细胞数据分析流程 - 第一步（数据整理 + 读入 Seurat）
#
# 功能：
#   1) 将 HRS 开头的"未压缩"10x 文件夹（barcodes.tsv / genes.tsv / matrix.mtx）
#      gzip 成标准 10x 格式，并把 genes.tsv 改名为 features.tsv，
#      使全部样本与 HRR/GSE 文件夹（已是 .gz）格式统一。
#   2) 读取 data/样本分组信息.xlsx（Sample -> 分组/编号），作为 metadata。
#   3) 用 Seurat::Read10X 逐个读取全部样本，合并为一个 Seurat 对象并保存 .rds。
#
# 用法：本脚本位于 Rversion/01_initiation/ 下，数据 data/ 在其上层项目根目录。
#   无论从哪个工作目录调用均可（路径以脚本位置自动推算）：
#     Rscript Rversion/01_initiation/01_initiation.R
#   或在 RStudio 中 Source 本文件。
#
# 说明：脚本对 HRS 的压缩是「原地」且幂等的——已存在 .gz 的文件不会重复压缩，
#       可放心重复运行。运行前请关闭 Excel（避免 xlsx 被锁）。
# =====================================================================

## ---- 0. 配置（按需修改） ----
# 本脚本位于 Rversion/ 下；以下路径以「脚本自身所在目录」为基准计算，
# 因此无论从哪个工作目录调用，都能正确找到上层的 data/ 与同层的 output/。
script_dir <- (function() {
  f <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (is.null(f) || !nzchar(f)) {
    arg <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", grep("^--file=", arg, value = TRUE)[1])
  }
  if (is.null(f) || !nzchar(f)) return(getwd())
  dirname(normalizePath(f))
})()

DATA_DIR     <- file.path(script_dir, "..", "data")       # 样本文件夹（位于上层项目根目录）
XLSX_PATH    <- file.path(DATA_DIR, "样本分组信息.xlsx")  # 分组表
OUT_DIR      <- file.path(script_dir, "output")           # 输出目录（与脚本同层，自动创建）
OUT_RDS      <- file.path(OUT_DIR, "01_seurat_combined.rds")
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
cat(sprintf("发现 %d 个样本文件夹。\n", length(sample_dirs)))

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
  cat("没有需要新压缩的样本（HRS 此前已压缩或不存在）。\n")
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

# 分组表中有但无文件夹的样本（当前应无）
missing <- setdiff(meta$Sample, sample_dirs)
if (length(missing) > 0) {
  cat("警告: 分组表中有但无对应文件夹的样本（未读入）:",
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
cat("01initiation 完成。\n")
