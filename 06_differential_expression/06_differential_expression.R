# =====================================================================
# 06_differential_expression.R
# 单细胞数据分析流程 - 第六步（差异表达 + 富集分析）
#
# 功能：
#   1) 读取第五步注释后的对象（05_annotated.rds）
#   2) 差异表达分析：
#      a) FindAllMarkers：每个细胞类型 vs 其余细胞（该类型的 marker 基因）
#      b) FindMarkers：ypN0 vs ypN+ 总体比较（两组间的差异基因）
#      c) 可选：按细胞类型分别做 ypN0 vs ypN+ 比较（类型特异的分组差异）
#   3) 可视化：marker DotPlot（各类型 top marker）、火山图（分组比较）、
#      热图（top marker 表达）
#   4) 富集分析（AUCell 单细胞评分，参考 Lesson 4 教程，使用 SeuratExtend 包）：
#      GO（免疫过程）/ Hallmark 基因集富集，输出富集表与图
#
# 用法：
#   在 RStudio 中打开本文件并点击 Source，或在命令行执行：
#     "C:/Program Files/R/R-4.4.3/bin/Rscript.exe" 06_differential_expression.R
#   说明：所有路径按脚本自身所在目录自动推算，可从任意工作目录运行。
#
# 常用修改点：见下方 ## ---- 0. 配置 ---- 中的 CONFIG 参数。
# 参考：Lesson-3（DE/火山图）、Lesson-4（AUCell 富集）
# =====================================================================

## ---- 0. 配置（按需修改） ----

# 自动定位脚本所在目录（兼容 Rscript 与 RStudio 两种运行方式）：
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
INPUT_RDS <- file.path(script_dir, "..", "05_celltype_annotation", "output", "05_annotated.rds") # 第五步注释对象
OUT_DIR   <- file.path(script_dir, "output")                                                      # 本步输出目录（自动创建）

# ---- 分组与注释列 ----
GROUP_COL    <- "group"      # 分组列（ypN0 / ypN+）
CELLTYPE_COL <- "celltype"   # 注释列（05 步生成；未注释时回退到 cluster）
CLUSTER_COL  <- "cluster"    # 聚类列（回退用）

# ---- 差异表达参数 ----
DE_LOGFC_THRESHOLD  <- 0.25   # FindMarkers/FindAllMarkers 的 log2FC 阈值（过滤微小差异）
DE_MIN_PCT          <- 0.1    # 基因至少在多少比例的细胞中表达才纳入测试（min.pct）
DE_ONLY_POS         <- TRUE   # FindAllMarkers 是否只保留上调 marker（TRUE 常见）
DE_TOPMARKERS       <- 5      # 每个细胞类型取 top N marker 用于热图/点图
RUN_GROUP_DE        <- TRUE   # 是否按细胞类型分别做 ypN0 vs ypN+ 比较（较耗时）
DE_WORKERS          <- 4      # FindAllMarkers/FindMarkers 的并行 worker 数（每个 worker 会复制一份表达矩阵，别开太大）

# ---- 结果复用 ----
# 置 TRUE 时：若 output/ 中已存在对应结果 CSV，则直接读取、不再重算。
# 用途：只改了富集/绘图逻辑时，避免把 50 多分钟的 DE 重跑一遍；改了 DE 参数请置 FALSE。
REUSE_EXISTING      <- TRUE

# ---- 富集分析参数 ----
RUN_ENRICHMENT      <- TRUE   # 是否运行 AUCell 富集（需安装 SeuratExtend）
SPECIES             <- "human" # 物种（"human" / "mouse"）；SeuratExtend 通过 getOption("spe") 读取，不设置会直接报错
ENRICHMENT_CORES    <- 4      # 富集并行核心数
ENRICH_GO_PARENT    <- "immune_system_process" # GO 父项（免疫相关过程，Lesson 4 用法）
RUN_HALLMARK        <- TRUE   # 是否额外跑 MSigDB Hallmark（50 个基因集，肿瘤免疫解读性好）
ENRICH_CELLS_PER_CT <- 500    # AUCell 抽样上限：每个「细胞类型 × 分组」层最多取多少细胞（控制排名矩阵内存）
ENRICH_TOP_TERMS    <- 30     # GO 热图展示的通路数（按细胞类型间标准差取变异最大的前 N 个）
ENRICH_TOP_SIG      <- 50     # 显著通路清单最多写多少条（写论文时直接贴表）
ENRICH_MIN_CT_CONSENSUS <- 3  # 跨细胞类型一致变化的阈值：至少在多少个类型中同向显著才算"共识通路"

# ---- 图片参数 ----
IMG_DPI      <- 300           # PNG 分辨率
X_TEXT_ANGLE <- 45            # x 轴刻度文本倾角（度），与 05 步保持一致
X_TEXT_HJUST <- 1             # 倾角对应的水平对齐方式

# ---- 其他 ----
SEED <- 123                  # 随机种子（保持可复现）

## ---- 1. 依赖检查与加载 ----

# 通用依赖安装函数：ensure_pkg("包名", bioc=是否 Bioconductor 包)
ensure_pkg <- function(pkg, bioc = FALSE) {
  if (requireNamespace(pkg, quietly = TRUE)) return(TRUE)             # 已装则跳过
  cat(sprintf("  正在安装依赖包: %s ...\n", pkg))                    # 提示开始安装
  if (bioc) {                                                         # Bioc 包分支
    if (!requireNamespace("BiocManager", quietly = TRUE)) {           # 先确保 BiocManager
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, update = FALSE, ask = FALSE)            # 安装 Bioc 包
  } else {                                                            # CRAN 包分支
    install.packages(pkg, repos = "https://cloud.r-project.org")      # 安装 CRAN 包
  }
  requireNamespace(pkg, quietly = TRUE)                               # 返回是否安装成功
}

# 依次确保需要的包可用
for (p in c("Seurat", "dplyr", "tidyr", "ggplot2", "patchwork", "ggrepel")) ensure_pkg(p)

# 加载核心包
library(Seurat)    # 单细胞分析主包
library(dplyr)     # 数据整理
library(tidyr)     # 长宽格式转换
library(ggplot2)   # 绘图
library(patchwork) # 图组合

# 创建输出目录（已存在则不报错）
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 物种设置：SeuratExtend 的富集函数（GeneSetAnalysisGO / RenameGO 等）内部通过
# getOption("spe") 取物种，未设置时报 "Species not defined. Please set with
# options(spe = 'your_species')"。这里统一按配置写入，避免只在富集处临时设置而遗漏其它调用。
options(spe = SPECIES)

# 打印运行环境信息
cat(sprintf("工作目录: %s\n", getwd()))                       # 当前工作目录
cat(sprintf("输出目录: %s\n", normalizePath(OUT_DIR)))        # 输出目录
cat(sprintf("物种(spe): %s\n", SPECIES))                      # 富集用物种

## ---- 1.1 通用小工具 ----

# x 轴刻度文本倾斜 45°（与 05 步的 theme_x45 一致，避免样本/通路名重叠）
theme_x45 <- function(angle = X_TEXT_ANGLE, hjust = X_TEXT_HJUST) {
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = angle, hjust = hjust))
}

# 同时保存 PDF 与 PNG（同尺寸、同分辨率，便于插入文档与快速预览）
save_plot <- function(p, filename, width, height, dpi = IMG_DPI) {
  ggplot2::ggsave(file.path(OUT_DIR, paste0(filename, ".pdf")), p,
                  width = width, height = height, limitsize = FALSE)
  ggplot2::ggsave(file.path(OUT_DIR, paste0(filename, ".png")), p,
                  width = width, height = height, dpi = dpi, limitsize = FALSE)
  cat(sprintf("  已保存图表: %s.pdf / %s.png\n", filename, filename))
}

# 写出带 UTF-8 BOM 的 CSV（Excel 直接双击打开中文/特殊字符不乱码）。
# 注意：R 4.4.3 的 write.csv 不支持 fileEncoding="UTF-8-BOM"，故先写临时文件再补 BOM。
write_utf8_csv <- function(df, path, row.names = FALSE) {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  write.csv(df, tmp, row.names = row.names, fileEncoding = "UTF-8")
  body <- readBin(tmp, "raw", n = file.info(tmp)$size)
  con <- file(path, open = "wb")
  writeBin(as.raw(c(0xef, 0xbb, 0xbf)), con)
  writeBin(body, con)
  close(con)
}

# 读取 CSV（自动去掉可能存在的 UTF-8 BOM，否则首列名会变成 "\ufeffgene"）
read_csv_safe <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(x)[1] <- sub("^\ufeff", "", names(x)[1])
  x
}

# 物理内存检查（原 wmic 调用已被删除）。
# 原因：wmic.exe 在当前安全策略下会被拦截，导致进程异常；DE/富集计算密集是既定事实，
#   不再实时查询内存，由用户在运行前自行保证可用内存 >=8GB。

# Seurat v5 并行/全局对象大小限制
options(future.globals.maxSize = 8 * 1024^3)
set.seed(SEED)

# 并行 plan：让 FindAllMarkers / FindMarkers 在 8 万细胞上自动用多核计算。
# 原因：单机顺序 Wilcoxon 在 8 万细胞上需要数十分钟甚至数小时；开启 multisession 并行
#   通常能加速 3~6 倍（与 CPU 核数相关）。如果 presto 也可用，Seurat 会自动调用更快的
#   presto::wilcoauc 实现，叠加效果显著。
if (requireNamespace("future", quietly = TRUE)) {
  # worker 数取配置值与实际核数的较小者：每个 worker 都会复制一份表达矩阵，
  # 开太多会直接把内存吃满（8 万细胞 × 3.6 万基因的矩阵单个就有数 GB）。
  nw <- max(1L, min(DE_WORKERS, parallel::detectCores(logical = TRUE)))
  future::plan("multisession", workers = nw)
  cat(sprintf("已开启 future 并行（multisession, workers=%d）\n", nw))
}

## ---- 2. 读入对象与确定分组依据 ----

# 2.1 检查输入文件是否存在
if (!file.exists(INPUT_RDS)) stop("找不到输入 rds: ", INPUT_RDS)

# 2.2 读取第五步的注释对象（读取需 1~3 分钟）
cat("正在读取 05_annotated.rds（约需 1~3 分钟）...\n")
obj <- readRDS(INPUT_RDS)                          # 读入 Seurat 对象
cat(sprintf("读入完成：%d 个基因 x %d 个细胞\n", nrow(obj), ncol(obj))) # 打印维度

# 2.3 确定 DE 的分组列：优先用细胞类型注释（celltype），未注释时回退到 cluster
if (CELLTYPE_COL %in% colnames(obj@meta.data) &&
    length(unique(obj@meta.data[[CELLTYPE_COL]])) > 1 &&
    !all(obj@meta.data[[CELLTYPE_COL]] == "Unknown")) { # 已有效注释
  de_group <- CELLTYPE_COL                          # 用注释列
  cat(sprintf("检测到细胞类型注释（%s），DE 将以细胞类型为分组。\n", CELLTYPE_COL))
} else if (CLUSTER_COL %in% colnames(obj@meta.data)) { # 未注释则回退
  de_group <- CLUSTER_COL                           # 用聚类列
  cat("警告: 未检测到有效注释（可能 05 步尚未填写 celltype_map.csv），本次 DE 以 cluster 为分组。\n")
} else {
  stop("对象中既无注释列也无聚类列，无法进行 DE。") # 均缺失则终止
}
cat(sprintf("DE 分组列: %s（%d 个组）\n", de_group, length(unique(obj@meta.data[[de_group]])))) # 打印组数

# 2.4 检查分组列（ypN0/ypN+）是否存在（分组比较用）
if (!GROUP_COL %in% colnames(obj@meta.data)) {     # 分组列缺失
  stop("对象中找不到分组列 '", GROUP_COL, "'，无法进行 ypN0 vs ypN+ 比较") # 终止
}
cat("分组分布：\n")                                 # 打印分组分布
print(table(obj@meta.data[[GROUP_COL]], useNA = "ifany")) # ypN0/ypN+ 细胞数

# 2.5 将 DE 分组列设为 Seurat 当前 idents（FindMarkers/FindAllMarkers 基于 idents 分组）
Idents(obj) <- obj@meta.data[[de_group]]

## ---- 3. 差异表达分析 ----

cat("\n==== 3. 差异表达分析 ====\n")

# 3.1 FindAllMarkers：每个细胞类型 vs 其余细胞（该类型的特异性 marker 基因）
#     only.pos=TRUE 只保留上调基因（细胞类型 marker 通常关注上调）
f_markers <- file.path(OUT_DIR, "01_markers_by_celltype.csv")
if (REUSE_EXISTING && file.exists(f_markers)) {    # 复用已有结果（省去数十分钟计算）
  all_markers <- read_csv_safe(f_markers)
  cat(sprintf("  复用已有结果: 01_markers_by_celltype.csv（共 %d 行 marker）\n", nrow(all_markers)))
} else {
  cat("运行 FindAllMarkers（每类型 vs 其余，可能较慢）...\n")
  all_markers <- FindAllMarkers(obj,
    only.pos = DE_ONLY_POS,                        # 只保留上调 marker
    logfc.threshold = DE_LOGFC_THRESHOLD,          # log2FC 阈值
    min.pct = DE_MIN_PCT,                          # 最小表达细胞比例
    verbose = FALSE)                               # 静默模式
  write_utf8_csv(all_markers, f_markers)           # 写 marker 表
  cat(sprintf("  已保存: 01_markers_by_celltype.csv（共 %d 行 marker）\n", nrow(all_markers))) # 提示
}

# 3.2 每个细胞类型取 top N marker（按平均 log2FC 排序），用于后续热图/点图
top_markers <- all_markers %>%
  dplyr::group_by(cluster) %>%                     # 按细胞类型分组
  dplyr::slice_max(order_by = avg_log2FC, n = DE_TOPMARKERS) %>% # 取 logFC 最高的 N 个
  dplyr::ungroup()                                 # 取消分组
cat("各细胞类型 top marker：\n")                    # 打印提示
print(as.data.frame(top_markers[, c("cluster", "gene", "avg_log2FC", "p_val_adj")])) # 打印 top marker

# 3.3 分组比较（ypN0 vs ypN+ 总体）：FindMarkers 需要 idents 含分组标签，
#     先临时切换 idents 到分组列，比较后恢复
f_overall <- file.path(OUT_DIR, "02_DE_ypN0_vs_ypNplus.csv")
if (REUSE_EXISTING && file.exists(f_overall)) {    # 复用已有结果
  de_group_overall <- read_csv_safe(f_overall)
  cat(sprintf("  复用已有结果: 02_DE_ypN0_vs_ypNplus.csv（%d 个差异基因）\n", nrow(de_group_overall)))
} else {
  cat("\n运行 ypN0 vs ypN+ 总体差异比较...\n")
  Idents(obj) <- obj@meta.data[[GROUP_COL]]        # idents 切换为分组
  de_group_overall <- FindMarkers(obj,
    ident.1 = setdiff(unique(obj@meta.data[[GROUP_COL]]), "ypN+")[1], # 第一组（ypN0）
    ident.2 = "ypN+",                              # 第二组（ypN+）
    logfc.threshold = DE_LOGFC_THRESHOLD,          # log2FC 阈值
    min.pct = DE_MIN_PCT,                          # 最小表达比例
    verbose = FALSE)                               # 静默模式
  de_group_overall$gene <- rownames(de_group_overall) # 基因名列
  de_group_overall <- de_group_overall %>%         # 调整列顺序（基因放前面）
    dplyr::select(gene, dplyr::everything())        # 基因列移到最前
  write_utf8_csv(de_group_overall, f_overall)      # 写结果
  cat(sprintf("  已保存: 02_DE_ypN0_vs_ypNplus.csv（%d 个差异基因）\n", nrow(de_group_overall))) # 提示
}

# 3.4 按细胞类型分别做 ypN0 vs ypN+ 比较（类型特异的分组差异，可选）
if (RUN_GROUP_DE) {                                # 开关打开时
  f_by_ct <- file.path(OUT_DIR, "03_DE_by_celltype_group.csv")
  if (REUSE_EXISTING && file.exists(f_by_ct)) {    # 复用已有结果
    de_by_ct_all <- read_csv_safe(f_by_ct)
    cat("\n复用已有结果: 03_DE_by_celltype_group.csv\n")
    for (ct in names(table(de_by_ct_all$celltype))) { # 打印各类型条目数
      cat(sprintf("  %s: %d 个差异基因\n", ct, sum(de_by_ct_all$celltype == ct)))
    }
  } else {
  cat("\n按细胞类型分别做 ypN0 vs ypN+ 比较（RUN_GROUP_DE=TRUE）...\n")
  Idents(obj) <- obj@meta.data[[de_group]]         # idents 切回细胞类型
  de_by_ct <- list()                               # 存储各类型结果
  for (ct in unique(obj@meta.data[[de_group]])) {  # 遍历细胞类型
    cells_ct <- WhichCells(obj, idents = ct)       # 该类型的细胞
    grp_vec <- obj@meta.data[cells_ct, GROUP_COL]  # 这些细胞的分组
    if (length(unique(grp_vec)) < 2) {             # 该类型只有单组时跳过
      cat(sprintf("  跳过 %s（无两个分组）\n", ct))
      next
    }
    sub <- subset(obj, idents = ct)                # 取该类型子集
    Idents(sub) <- sub@meta.data[[GROUP_COL]]      # 子集内按分组设置 idents
    res <- FindMarkers(sub,                        # 组内差异比较
      ident.1 = "ypN0", ident.2 = "ypN+",          # 两组
      logfc.threshold = DE_LOGFC_THRESHOLD,        # 阈值
      min.pct = DE_MIN_PCT,                        # 最小表达比例
      verbose = FALSE)                             # 静默模式
    res$gene <- rownames(res)                      # 基因列
    res$celltype <- ct                             # 标注细胞类型
    de_by_ct[[ct]] <- res                          # 存入列表
    cat(sprintf("  %s: %d 个差异基因\n", ct, nrow(res))) # 打印行数
  }
  de_by_ct_all <- dplyr::bind_rows(de_by_ct)       # 合并所有类型结果
  write_utf8_csv(de_by_ct_all, f_by_ct)            # 写长表
  cat("  已保存: 03_DE_by_celltype_group.csv（按细胞类型 × 分组的差异表）\n") # 提示
  }
}

# 3.5 恢复 idents 为细胞类型（供后续可视化使用）
Idents(obj) <- obj@meta.data[[de_group]]

## ---- 4. 可视化 ----

cat("\n==== 4. 可视化 ====\n")

# 4.1 marker DotPlot：展示各细胞类型 top marker 的表达（点大小=表达比例，颜色=平均表达）
if (nrow(top_markers) > 0 && "gene" %in% colnames(top_markers)) { # 有 marker 时
  feat_order <- unique(top_markers$gene)           # 基因顺序（按类型排列）
  p_dot <- DotPlot(obj,                            # 点图
    features = feat_order,                         # 展示的基因
    group.by = de_group,                           # 按细胞类型分组
    cols = c("lightgrey", "#E64B35")) +            # 低-高表达配色
    theme_minimal() +                              # 简洁主题
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # 旋转基因标签
    labs(title = paste0("Top ", DE_TOPMARKERS, " markers by ", de_group)) # 图标题
  save_plot(p_dot, "04_marker_dotplot",             # 保存点图（PDF + PNG 同尺寸）
            width  = max(8, length(feat_order) * 0.5), # 宽度随基因数自适应
            height = max(6, length(unique(obj@meta.data[[de_group]])) * 0.5)) # 高度随类型数自适应
}

# 4.2 火山图：ypN0 vs ypN+ 总体比较（x=log2FC，y=-log10(调整p值)），
#     标记显著差异基因（|logFC|>阈值 且 p.adj<0.05）
if (nrow(de_group_overall) > 0) {                  # 有结果时
  vol <- de_group_overall                          # 引用结果
  vol$sig <- ifelse(vol$p_val_adj < 0.05 & vol$avg_log2FC > DE_LOGFC_THRESHOLD, "up",   # 显著上调
             ifelse(vol$p_val_adj < 0.05 & vol$avg_log2FC < -DE_LOGFC_THRESHOLD, "down", "ns")) # 显著下调/不显著
  # 标注 top 基因（显著中 |logFC| 最大的 15 个）
  top_genes <- vol %>%
    dplyr::filter(sig != "ns") %>%                 # 只取显著基因
    dplyr::slice_max(order_by = abs(avg_log2FC), n = 15) # |logFC| 最大的 15 个
  p_vol <- ggplot(vol, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) + # 火山图
    geom_point(size = 0.5, alpha = 0.6) +          # 散点
    scale_color_manual(values = c(up = "#E64B35", down = "#3C5488", ns = "grey70")) + # 红=上调,蓝=下调
    ggrepel::geom_text_repel(                      # 基因名标注（防重叠）
      data = top_genes, aes(label = gene),         # 标注 top 基因
      size = 2.5, max.overlaps = 20) +             # 标注大小与重叠控制
    geom_vline(xintercept = c(-DE_LOGFC_THRESHOLD, DE_LOGFC_THRESHOLD), linetype = "dashed", color = "grey40") + # FC 阈值线
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") + # p 阈值线
    theme_minimal() +                              # 简洁主题
    labs(title = "Volcano plot: ypN0 vs ypN+",     # 图标题
         x = "avg log2FC", y = "-log10(adjusted p-value)", color = "Regulation") # 坐标轴
  save_plot(p_vol, "05_volcano_ypN0_vs_ypNplus",   # 保存火山图
            width = 9, height = 7)                 # 画布尺寸
}

# 4.3 top marker 表达热图（抽样版）：每个类型抽最多 100 细胞、画 top marker 的 Z-score 表达
if (nrow(top_markers) > 0 && "gene" %in% colnames(top_markers)) { # 有 marker 时
  set.seed(SEED)                                   # 固定随机种子
  cells_use <- unlist(lapply(                      # 每类型抽样细胞
    unique(obj@meta.data[[de_group]]),             # 遍历类型
    function(ct) {                                 # 每类型
      cells_ct <- WhichCells(obj, idents = ct)     # 该类型细胞
      if (length(cells_ct) > 100) sample(cells_ct, 100) else cells_ct # 最多 100 个
    }))
  obj_hm <- subset(obj, cells = cells_use)         # 抽样子集（控制热图规模）
  p_hm <- DoHeatmap(obj_hm,                        # 热图
    features = unique(top_markers$gene),           # top marker 基因
    group.by = de_group,                           # 按类型分组
    angle = 45, size = 3) +                        # 标签角度与字号
    scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35") + # 蓝-白-红
    theme(legend.position = "right")               # 图例位置
  save_plot(p_hm, "06_marker_heatmap",             # 保存热图（每类型抽样 100 细胞）
            width = 12, height = 8)
}

## ---- 5. 富集分析（AUCell 通路活性，参考 Lesson 4 教程）----

# 与"先找差异基因再做 ORA 富集"不同，AUCell 在**单细胞层面**给每个细胞算出
# 每条通路的连续活性得分（AUC），于是可以直接回答：
#   "干扰素通路在 ypN+ 的 T 细胞里是否整体高于 ypN0？"（Wilcoxon 检验）
# 输出两组结果：通路 × 细胞类型的平均活性（热图），以及通路 × 细胞类型的两组差异（热图 + 检验表）。
if (RUN_ENRICHMENT) {                              # 开关打开时
  cat("\n==== 5. 富集分析（AUCell）====\n")

  # 5.1 确保 SeuratExtend 可用（GitHub 包，需 remotes 安装；失败则降级跳过）
  if (!requireNamespace("SeuratExtend", quietly = TRUE)) { # 未安装时
    cat("正在安装 SeuratExtend（GitHub 包，需数分钟）...\n") # 提示
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org") # 装 remotes
    tryCatch(                                     # 安装可能失败（网络/依赖）
      remotes::install_github("huayc09/SeuratExtend", upgrade = "never", quiet = TRUE), # 从 GitHub 安装（仓库见 Lesson-4 教程）
      error = function(e) cat("SeuratExtend 安装失败:", conditionMessage(e), "\n"))
  }

  if (!requireNamespace("SeuratExtend", quietly = TRUE)) { # 装不上时
    warning("SeuratExtend 不可用，跳过富集分析（DE 结果已正常输出）")
  } else {
    # 5.2 载入注释数据集。
    #     这是之前 GO 富集失败的第二个根因：GO_Data / Genesets_data 是
    #     SeuratExtendData 里的**懒加载数据集**，而 SeuratExtend 的内部函数
    #     （GetChilrenGO 等）是在全局环境里按名字找 GO_Data 的；不显式 data() 载入，
    #     即使把 options(spe=...) 设对了，也仍会报 "object 'GO_Data' not found"。
    #     第一个根因（Species not defined）已在配置区用 options(spe=SPECIES) 解决。
    cat("载入 GO 注释数据集 GO_Data（约 450MB）...\n")
    utils::data("GO_Data", package = "SeuratExtendData", envir = globalenv())

    # 5.3 分层抽样：AUCell 要为**每个细胞**建一张全基因排名矩阵，8 万细胞的排名矩阵
    #     常驻内存可达数 GB。这里按「细胞类型 × 分组」分层抽样（默认每层 500 个），
    #     内存可控，同时保证每个"类型 × 分组"格子都有足够细胞做统计检验。
    set.seed(SEED)
    ct_vec  <- as.character(obj@meta.data[[de_group]])
    grp_vec <- as.character(obj@meta.data[[GROUP_COL]])
    strat   <- interaction(ct_vec, grp_vec, sep = "|", drop = TRUE)
    # 先按位置抽样，再映射成细胞名（subset() 对细胞名最稳妥）
    cells_use <- colnames(obj)[unlist(lapply(levels(strat), function(k) { # 逐层抽样
      idx <- which(strat == k)
      if (length(idx) > ENRICH_CELLS_PER_CT) sample(idx, ENRICH_CELLS_PER_CT) else idx
    }))]
    cat(sprintf("富集子集细胞数: %d（按 %s × %s 分层，每层上限 %d）\n",
                length(cells_use), de_group, GROUP_COL, ENRICH_CELLS_PER_CT))
    obj_enr <- subset(obj, cells = cells_use)       # 抽样子集（只用于富集评分）
    ct_enr  <- factor(as.character(obj_enr@meta.data[[de_group]]),
                      levels = levels(obj@meta.data[[de_group]])) # 保持与 05 一致的类型顺序
    grp_enr <- as.character(obj_enr@meta.data[[GROUP_COL]])
    print(table(ct_enr, grp_enr, useNA = "ifany"))  # 打印各"类型×分组"格子的细胞数

    # 5.4 通用汇总函数：把「通路 × 细胞」AUC 矩阵整理成表格与热图
    #     auc_mat   : 行=通路，列=细胞的 AUC 矩阵
    #     term_names: 通路的可读名称（与 auc_mat 行名一一对应）
    #     file_prefix: 输出文件名前缀，如 "07_GO" / "11_Hallmark"
    #     top_n     : 热图最多展示多少条通路（按细胞类型间变异或 |Δ| 排序取前 N）
    summarize_auc <- function(auc_mat, term_names, file_prefix, top_n) {
      ct_levels <- levels(ct_enr)

      # (a) 通路 × 细胞类型 的平均活性
      mean_ct <- do.call(cbind, lapply(ct_levels, function(l) { # 逐类型求均值
        rowMeans(auc_mat[, ct_enr == l, drop = FALSE])
      }))
      colnames(mean_ct) <- ct_levels
      rownames(mean_ct) <- rownames(auc_mat)
      tbl_ct <- data.frame(term_id = rownames(auc_mat), term = term_names,
                           round(mean_ct, 5), check.names = FALSE)
      write_utf8_csv(tbl_ct, file.path(OUT_DIR, paste0(file_prefix, "_AUC_by_celltype.csv")))
      cat(sprintf("  已保存: %s_AUC_by_celltype.csv（%d 条通路 × %d 个细胞类型）\n",
                  file_prefix, nrow(auc_mat), length(ct_levels)))

      # 热图 1：取细胞类型间标准差最大的前 top_n 条通路，按行做 z-score
      n_show <- min(top_n, nrow(mean_ct))
      keep   <- head(order(apply(mean_ct, 1, stats::sd), decreasing = TRUE), n_show)
      mat_z  <- t(scale(t(mean_ct[keep, , drop = FALSE]))) # 行方向 z-score
      mat_z[is.na(mat_z)] <- 0                            # 标准差为 0 的行置 0
      p_hm <- SeuratExtend::Heatmap(mat_z, color_scheme = "BuRd",
                                    angle = X_TEXT_ANGLE, hjust = X_TEXT_HJUST) +
        labs(title = paste0("AUCell pathway activity (z-score): ", file_prefix),
             x = "cell type", y = "pathway")
      save_plot(p_hm, paste0(file_prefix, "_AUC_heatmap_by_celltype"),
                width = 10, height = max(6, n_show * 0.35))

      # (b) 通路 × 细胞类型 的 ypN0 vs ypN+ 差异（Wilcoxon + BH 校正）
      cmp_list <- list()
      for (l in ct_levels) {                              # 逐类型检验
        i0 <- which(ct_enr == l & grp_enr == "ypN0")      # 该类型的 ypN0 细胞
        i1 <- which(ct_enr == l & grp_enr == "ypN+")      # 该类型的 ypN+ 细胞
        if (length(i0) < 10 || length(i1) < 10) {         # 细胞太少不做检验
          cat(sprintf("    跳过 %s：某组细胞数不足 10（ypN0=%d, ypN+=%d）\n", l, length(i0), length(i1)))
          next
        }
        sub <- auc_mat[, c(i0, i1), drop = FALSE]         # 只取这两组细胞
        # 很多小通路的 AUCell 得分在该类型里恒为 0（该通路基因在这些细胞中都不表达），
        # 秩和检验对它无意义且会返回 NA，这里先剔除，保证 p_adj 列干净可用。
        sd_row <- base::apply(sub, 1, stats::sd)
        keep_r <- is.finite(sd_row) & sd_row > 0
        if (!any(keep_r)) {                               # 该类型所有通路都是常数
          cat(sprintf("    跳过 %s：所有通路得分均为常数，无法检验\n", l))
          next
        }
        sub <- sub[keep_r, , drop = FALSE]
        m0  <- rowMeans(sub[, seq_along(i0), drop = FALSE])
        m1  <- rowMeans(sub[, length(i0) + seq_along(i1), drop = FALSE])
        pv  <- base::apply(sub, 1, function(v)            # 逐通路做秩和检验
          suppressWarnings(stats::wilcox.test(v[seq_along(i0)],
                                              v[length(i0) + seq_along(i1)],
                                              exact = FALSE)$p.value))
        cmp_list[[l]] <- data.frame(celltype = l,
                                    term_id  = rownames(sub),
                                    term     = term_names[keep_r],
                                    AUC_ypN0   = round(m0, 5),
                                    AUC_ypNplus = round(m1, 5),
                                    diff_ypNplus_minus_ypN0 = round(m1 - m0, 5),
                                    p_value  = pv, stringsAsFactors = FALSE)
      }
      cmp <- dplyr::bind_rows(cmp_list)
      if (nrow(cmp) == 0) return(invisible(NULL))          # 全部类型被跳过
      cmp$p_adj <- stats::p.adjust(cmp$p_value, method = "BH") # BH 多重校正
      # 方向标签：ypN+ 高 = 红色（上调，= ypN+ 减去 ypN0 的差为正）；ypN0 高 = 蓝色
      cmp$direction <- ifelse(cmp$p_adj < 0.05 & cmp$diff_ypNplus_minus_ypN0 > 0, "ypN+ higher",
                       ifelse(cmp$p_adj < 0.05 & cmp$diff_ypNplus_minus_ypN0 < 0, "ypN0 higher", "ns"))
      cmp <- cmp[order(-abs(cmp$diff_ypNplus_minus_ypN0)), ]
      write_utf8_csv(cmp, file.path(OUT_DIR, paste0(file_prefix, "_AUC_ypN0_vs_ypNplus.csv")))
      cat(sprintf("  已保存: %s_AUC_ypN0_vs_ypNplus.csv（%d 条 通路×细胞类型 检验）\n",
                  file_prefix, nrow(cmp)))

      # (b1) 显著通路清单：按 |Δ| 排序取前 ENRICH_TOP_SIG 条，写论文时可直接贴表
      top_sig <- head(cmp[cmp$p_adj < 0.05, ], ENRICH_TOP_SIG)
      write_utf8_csv(top_sig, file.path(OUT_DIR, paste0(file_prefix, "_AUC_top_significant.csv")))
      cat(sprintf("  已保存: %s_AUC_top_significant.csv（最多 %d 条）\n",
                  file_prefix, ENRICH_TOP_SIG))

      # (b2) 火山图：x=ΔAUC(ypN+ - ypN0)，y=-log10(p_adj)，色=方向；
      #      标注 top 15 |Δ| 显著通路（"细胞类型: 通路名"）。
      #      颜色约定：ypN+ 高 = 红（与 05/06 已有火山图一致），ypN0 高 = 蓝。
      if (nrow(cmp) >= 5) {
        cmp$neglogp <- -log10(pmax(cmp$p_adj, .Machine$double.xmin))
        p_vol <- ggplot2::ggplot(cmp, ggplot2::aes(x = diff_ypNplus_minus_ypN0,
                                                   y = neglogp,
                                                   color = direction)) +
          ggplot2::geom_point(size = 1.2, alpha = 0.7) +
          ggplot2::scale_color_manual(values = c("ypN+ higher" = "#E64B35",
                                                "ypN0 higher" = "#3C5488",
                                                "ns"          = "grey70")) +
          ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
          ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
          ggplot2::theme_minimal() +
          ggplot2::labs(title = paste0("Volcano of pathway activity: ypN+ vs ypN0 — ", file_prefix),
                        x = "ΔAUC (ypN+ minus ypN0)",
                        y = "-log10(adjusted p)", color = "Direction")
        if (nrow(top_sig) > 0) {
          # 标签数据单独构造。
          # 注意：必须带上 direction 列——全局 aes 里有 color = direction，
          # 而 geom_text_repel 默认 inherit.aes=TRUE，会把 color=direction 也继承过来；
          # 标签数据里没有这一列就会报 "Aesthetics are not valid data columns: colour = direction"。
          ts <- head(top_sig, 15)
          label_df <- data.frame(
            diff_ypNplus_minus_ypN0 = ts$diff_ypNplus_minus_ypN0,
            neglogp                 = -log10(pmax(ts$p_adj, .Machine$double.xmin)),
            direction               = ts$direction,   # 供继承的 color 映射使用
            label                   = paste0(ts$celltype, ": ",
                                             ifelse(nchar(ts$term) > 30,
                                                    paste0(substr(ts$term, 1, 28), ".."),
                                                    ts$term)),
            stringsAsFactors = FALSE)
          p_vol <- p_vol + ggrepel::geom_text_repel(
            data = label_df,
            ggplot2::aes(label = label),
            size = 2.2, max.overlaps = 25, min.segment.length = 0, show.legend = FALSE)
        }
        save_plot(p_vol, paste0(file_prefix, "_AUC_volcano_by_group"),
                  width = 11, height = 7)
      }

      # 热图 2：ΔAUC（ypN+ 减 ypN0），取跨类型 |Δ| 最大的前 top_n 条通路
      d_mat <- do.call(cbind, lapply(ct_levels, function(l) { # 组装 term × celltype 的 Δ 矩阵
        s <- cmp[cmp$celltype == l, ]
        if (nrow(s) == 0) return(rep(NA_real_, nrow(auc_mat)))
        s$diff_ypNplus_minus_ypN0[match(rownames(auc_mat), s$term_id)]
      }))
      colnames(d_mat) <- ct_levels
      rownames(d_mat) <- rownames(auc_mat)
      d_mat <- d_mat[, colSums(is.na(d_mat)) < nrow(d_mat), drop = FALSE] # 去掉整体缺失的类型
      if (ncol(d_mat) == 0) return(invisible(NULL))
      # 只保留至少在一种细胞类型中被检验过的通路（全 NA 行没有信息，直接排在最后并剔除）
      d_max <- base::apply(abs(d_mat), 1, function(v)
        if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
      keep2 <- head(order(d_max, decreasing = TRUE, na.last = NA),
                    min(top_n, sum(!is.na(d_max))))
      p_diff <- SeuratExtend::Heatmap(d_mat[keep2, , drop = FALSE], color_scheme = "BuRd",
                                      angle = X_TEXT_ANGLE, hjust = X_TEXT_HJUST) +
        labs(title = paste0("AUCell activity difference (ypN+ minus ypN0): ", file_prefix),
             x = "cell type", y = "pathway")
      save_plot(p_diff, paste0(file_prefix, "_AUC_diff_heatmap_by_group"),
                width = 10, height = max(6, length(keep2) * 0.35))

      # (b3) 跨细胞类型一致的通路（"共识通路"）：
      #      在 ≥ ENRICH_MIN_CT_CONSENSUS 个细胞类型中同向显著的通路，
      #      比起单类型显著更能代表全局组学差异；写论文时这块最稳。
      sig <- cmp[cmp$p_adj < 0.05 & cmp$direction != "ns", , drop = FALSE]
      if (nrow(sig) > 0) {
        cons <- sig %>%
          dplyr::group_by(term_id, term) %>%
          dplyr::summarise(
            n_ypNplus_higher = sum(direction == "ypN+ higher"),
            n_ypN0_higher    = sum(direction == "ypN0 higher"),
            n_consistent     = dplyr::first(n_ypNplus_higher) + dplyr::first(n_ypN0_higher), # 故意 N+ + N0，等于跨类型显著总次数
            mean_diff        = mean(diff_ypNplus_minus_ypN0),
            median_padj      = median(p_adj),
            .groups = "drop")
        # 把 n_consistent 替换为同向最大数（一个通路可能同时在 ypN+ 高 和 ypN0 高 类型里都显著）
        cons$n_consistent <- pmax(cons$n_ypNplus_higher, cons$n_ypN0_higher)
        cons$consistent_direction <- ifelse(cons$n_ypNplus_higher >= cons$n_ypN0_higher,
                                            "ypN+ higher", "ypN0 higher")
        cons <- cons[cons$n_consistent >= ENRICH_MIN_CT_CONSENSUS, , drop = FALSE]
        cons <- cons[order(-cons$n_consistent, -abs(cons$mean_diff)), , drop = FALSE]
        write_utf8_csv(cons, file.path(OUT_DIR, paste0(file_prefix, "_AUC_consensus_by_celltype.csv")))
        cat(sprintf("  已保存: %s_AUC_consensus_by_celltype.csv（%d 条共识通路，阈值≥%d 类型）\n",
                    file_prefix, nrow(cons), ENRICH_MIN_CT_CONSENSUS))
        # 共识热图：取 |mean_diff| 最大的前 min(top_n, nrow(cons)) 条
        if (nrow(cons) > 0) {
          n_cs <- min(nrow(cons), top_n)
          cons_top <- head(cons, n_cs)
          d_cons <- d_mat[match(cons_top$term_id, rownames(d_mat)), , drop = FALSE]
          rownames(d_cons) <- cons_top$term  # 用可读名（GO 已 RenameGO；Hallmark 已是名字）
          p_cs <- SeuratExtend::Heatmap(d_cons, color_scheme = "BuRd",
                                        angle = X_TEXT_ANGLE, hjust = X_TEXT_HJUST) +
            labs(title = paste0("Consensus pathways (>= ", ENRICH_MIN_CT_CONSENSUS,
                                " cell types same direction): ", file_prefix),
                 x = "cell type", y = "pathway")
          save_plot(p_cs, paste0(file_prefix, "_AUC_consensus_heatmap"),
                    width = 10, height = max(6, n_cs * 0.35))
        }
      }
      invisible(list(mean_ct = mean_ct, cmp = cmp))
    }

    # 5.5 带降级的 AUC 计算：先用 slot="counts"（AUCell 惯例），若对象里
    #     counts 层不可用（如 Seurat v5 分层/被拆分），自动改用 slot="data" 重试。
    run_auc <- function(call_fun, label) {
      tryCatch(call_fun("counts"),
        error = function(e) {
          cat(sprintf("  %s 用 slot=counts 失败（%s），改用 slot=data 重试\n",
                      label, conditionMessage(e)))
          call_fun("data")
        })
    }

    # 5.6 GO 富集（免疫系统过程父项）
    cat(sprintf("运行 GO 富集（parent=%s, nCores=%d）...\n", ENRICH_GO_PARENT, ENRICHMENT_CORES))
    go_mat <- tryCatch(
      run_auc(function(sl) SeuratExtend::GeneSetAnalysisGO(obj_enr,
                parent = ENRICH_GO_PARENT, slot = sl,
                nCores = ENRICHMENT_CORES, export_to_matrix = TRUE), "GO"),
      error = function(e) { warning("GO 富集失败: ", conditionMessage(e)); NULL })
    if (!is.null(go_mat)) {                        # 成功时
      saveRDS(go_mat, file.path(OUT_DIR, "07_GO_AUC_matrix.rds")) # 原始「通路 × 细胞」矩阵
      cat(sprintf("  GO 通路数（过滤后）: %d\n", nrow(go_mat)))
      go_names <- tryCatch(                        # GO ID → 可读名称
        SeuratExtend::RenameGO(rownames(go_mat), add_id = FALSE, add_n_gene = FALSE, spe = SPECIES),
        error = function(e) rownames(go_mat))      # 改名失败则保留 ID
      summarize_auc(go_mat, go_names, "07_GO", ENRICH_TOP_TERMS)
    }

    # 5.7 Hallmark 富集（MSigDB 50 个标志性基因集，肿瘤免疫解读性最好）
    if (RUN_HALLMARK) {
      cat("运行 Hallmark 富集（MSigDB 50 基因集）...\n")
      utils::data("Genesets_data", package = "SeuratExtendData", envir = globalenv()) # 载入基因集数据
      hallmark <- Genesets_data[[SPECIES]][["GSEA"]][["hallmark gene sets"]]          # 取 hallmark 集合
      hallmark <- lapply(hallmark[grepl("^HALLMARK_", names(hallmark))], unique)      # 去重
      names(hallmark) <- sub("^HALLMARK_", "", names(hallmark))                       # 去掉前缀，标签更短
      hm_mat <- tryCatch(
        run_auc(function(sl) SeuratExtend::GeneSetAnalysis(obj_enr,
                  genesets = hallmark, title = "hallmark", slot = sl,
                  nCores = ENRICHMENT_CORES, export_to_matrix = TRUE), "Hallmark"),
        error = function(e) { warning("Hallmark 富集失败: ", conditionMessage(e)); NULL })
      if (!is.null(hm_mat)) {
        saveRDS(hm_mat, file.path(OUT_DIR, "11_Hallmark_AUC_matrix.rds"))
        cat(sprintf("  Hallmark 通路数（过滤后）: %d\n", nrow(hm_mat)))
        summarize_auc(hm_mat, rownames(hm_mat), "11_Hallmark", nrow(hm_mat)) # 50 条全部展示
      }
    }

    rm(obj_enr)                                    # 释放富集子集内存
    if (exists("Genesets_data")) rm(Genesets_data) # 释放基因集数据
    invisible(gc())                                # 主动回收
  }
}

## ---- 6. 完成 ----

# 打印最终结果摘要。
# 注意：这里对 obj 做了判空保护——如果富集步骤出现异常把对象置空，
# 摘要本身不应该再崩一次（此前正是因此出现过
# "no applicable method for @ applied to an object of class NULL"）。
cat("\n==== 06 差异表达与富集完成 ====\n")
if (!is.null(obj)) {
  cat(sprintf("DE 分组: %s（%d 个组）\n", de_group,
              length(unique(obj@meta.data[[de_group]])))) # 分组信息
  cat(sprintf("细胞总数: %d\n", ncol(obj)))
} else {
  cat("警告: obj 已为空，跳过统计摘要（请检查上方富集步骤的报错）\n")
}
cat("\n输出文件清单（位于 output/ 目录）：\n")
cat("  01_markers_by_celltype.csv              - 各细胞类型 marker 基因表\n")
cat("  02_DE_ypN0_vs_ypNplus.csv               - ypN0 vs ypN+ 总体差异基因表\n")
cat("  03_DE_by_celltype_group.csv             - 按细胞类型的组间差异表（可选）\n")
cat("  04_marker_dotplot.pdf/.png              - top marker 点图\n")
cat("  05_volcano_ypN0_vs_ypNplus.pdf/.png     - 分组比较火山图\n")
cat("  06_marker_heatmap.pdf/.png              - top marker 表达热图\n")
cat("  07_GO_AUC_matrix.rds                    - GO：通路 × 细胞 的 AUCell 得分矩阵\n")
cat("  07_GO_AUC_by_celltype.csv               - GO：通路 × 细胞类型 平均活性\n")
cat("  07_GO_AUC_heatmap_by_celltype.pdf/.png  - GO 活性热图（z-score）\n")
cat("  07_GO_AUC_ypN0_vs_ypNplus.csv           - GO：通路 × 细胞类型 两组差异（Wilcoxon + BH）\n")
cat("  07_GO_AUC_diff_heatmap_by_group.pdf/.png- GO ΔAUC 热图（ypN+ 减 ypN0）\n")
cat("  11_Hallmark_AUC_matrix.rds              - Hallmark：通路 × 细胞 的 AUCell 得分矩阵\n")
cat("  11_Hallmark_AUC_by_celltype.csv         - Hallmark：通路 × 细胞类型 平均活性\n")
cat("  11_Hallmark_AUC_heatmap_by_celltype.*   - Hallmark 活性热图（z-score）\n")
cat("  11_Hallmark_AUC_ypN0_vs_ypNplus.csv     - Hallmark：通路 × 细胞类型 两组差异\n")
cat("  11_Hallmark_AUC_diff_heatmap_by_group.* - Hallmark ΔAUC 热图（ypN+ 减 ypN0）\n")
cat("\n06_differential_expression 完成。\n")
