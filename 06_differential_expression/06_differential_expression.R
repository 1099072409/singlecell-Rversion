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

# ---- 富集分析参数 ----
RUN_ENRICHMENT      <- TRUE   # 是否运行 AUCell 富集（需安装 SeuratExtend，较耗时）
ENRICHMENT_CORES    <- 2      # 富集并行核心数（Windows 建议小值，如 2）
ENRICH_GO_PARENT    <- "immune_system_process" # GO 父项（免疫相关过程，Lesson 4 用法）

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

# 打印运行环境信息
cat(sprintf("工作目录: %s\n", getwd()))                       # 当前工作目录
cat(sprintf("输出目录: %s\n", normalizePath(OUT_DIR)))        # 输出目录

# 检查可用物理内存并给出提示（Windows 下通过 wmic 查询；查询失败时自动跳过）
tryCatch({
  w  <- system("wmic OS get FreePhysicalMemory /value", intern = TRUE) # 查询可用内存(KB)
  fr <- as.numeric(sub(".*=", "", w[grepl("FreePhysicalMemory", w)]))  # 提取数值
  if (length(fr) == 1 && !is.na(fr)) {                                 # 查询成功时
    free_gb <- fr / 1024^2                                             # 转 GB
    cat(sprintf("当前可用物理内存: %.1f GB\n", free_gb))                # 打印可用内存
    if (free_gb < 8) {                                                 # 低于 8GB 时提示
      cat("提示: 本步 DE/富集计算密集，建议可用内存 >=8GB。\n")
    }
  }
}, error = function(e) invisible(NULL))                                 # 查询失败静默跳过

# Seurat v5 并行/全局对象大小限制
options(future.globals.maxSize = 8 * 1024^3)
set.seed(SEED)

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
cat("运行 FindAllMarkers（每类型 vs 其余，可能较慢）...\n")
all_markers <- FindAllMarkers(obj,
  only.pos = DE_ONLY_POS,                          # 只保留上调 marker
  logfc.threshold = DE_LOGFC_THRESHOLD,            # log2FC 阈值
  min.pct = DE_MIN_PCT,                            # 最小表达细胞比例
  verbose = FALSE)                                 # 静默模式
write.csv(all_markers, file.path(OUT_DIR, "01_markers_by_celltype.csv"), row.names = FALSE) # 写 marker 表
cat(sprintf("  已保存: 01_markers_by_celltype.csv（共 %d 行 marker）\n", nrow(all_markers))) # 提示

# 3.2 每个细胞类型取 top N marker（按平均 log2FC 排序），用于后续热图/点图
top_markers <- all_markers %>%
  dplyr::group_by(cluster) %>%                     # 按细胞类型分组
  dplyr::slice_max(order_by = avg_log2FC, n = DE_TOPMARKERS) %>% # 取 logFC 最高的 N 个
  dplyr::ungroup()                                 # 取消分组
cat("各细胞类型 top marker：\n")                    # 打印提示
print(as.data.frame(top_markers[, c("cluster", "gene", "avg_log2FC", "p_val_adj")])) # 打印 top marker

# 3.3 分组比较（ypN0 vs ypN+ 总体）：FindMarkers 需要 idents 含分组标签，
#     先临时切换 idents 到分组列，比较后恢复
cat("\n运行 ypN0 vs ypN+ 总体差异比较...\n")
Idents(obj) <- obj@meta.data[[GROUP_COL]]          # idents 切换为分组
de_group_overall <- FindMarkers(obj,
  ident.1 = setdiff(unique(obj@meta.data[[GROUP_COL]]), "ypN+")[1], # 第一组（ypN0）
  ident.2 = "ypN+",                                # 第二组（ypN+）
  logfc.threshold = DE_LOGFC_THRESHOLD,            # log2FC 阈值
  min.pct = DE_MIN_PCT,                            # 最小表达比例
  verbose = FALSE)                                 # 静默模式
de_group_overall$gene <- rownames(de_group_overall) # 基因名列
de_group_overall <- de_group_overall %>%           # 调整列顺序（基因放前面）
  dplyr::select(gene, dplyr::everything())          # 基因列移到最前
write.csv(de_group_overall, file.path(OUT_DIR, "02_DE_ypN0_vs_ypNplus.csv"), row.names = FALSE) # 写结果
cat(sprintf("  已保存: 02_DE_ypN0_vs_ypNplus.csv（%d 个差异基因）\n", nrow(de_group_overall))) # 提示

# 3.4 按细胞类型分别做 ypN0 vs ypN+ 比较（类型特异的分组差异，可选）
if (RUN_GROUP_DE) {                                # 开关打开时
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
  write.csv(de_by_ct_all, file.path(OUT_DIR, "03_DE_by_celltype_group.csv"), row.names = FALSE) # 写长表
  cat("  已保存: 03_DE_by_celltype_group.csv（按细胞类型 × 分组的差异表）\n") # 提示
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
  ggsave(file.path(OUT_DIR, "04_marker_dotplot.pdf"), # 保存点图
         p_dot, width = max(8, length(feat_order) * 0.5), height = max(6, length(unique(obj@meta.data[[de_group]])) * 0.5), # 尺寸自适应
         limitsize = FALSE)                        # 允许超出默认尺寸
  cat("  已保存图表: 04_marker_dotplot.pdf\n")     # 提示
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
  ggsave(file.path(OUT_DIR, "05_volcano_ypN0_vs_ypNplus.pdf"), # 保存火山图
         p_vol, width = 9, height = 7)             # 画布尺寸
  cat("  已保存图表: 05_volcano_ypN0_vs_ypNplus.pdf\n") # 提示
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
  ggsave(file.path(OUT_DIR, "06_marker_heatmap.pdf"), # 保存热图
         p_hm, width = 12, height = 8, limitsize = FALSE) # 尺寸
  cat("  已保存图表: 06_marker_heatmap.pdf（每类型抽样 100 细胞）\n") # 提示
}

## ---- 5. 富集分析（AUCell，参考 Lesson 4 教程）----

# 富集基于 SeuratExtend 包（封装了 AUCell 单细胞富集评分），
# 对每个细胞类型计算 GO（免疫过程）/Hallmark 基因集的富集程度。
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
  if (requireNamespace("SeuratExtend", quietly = TRUE)) { # 安装成功时
    # 5.2 GO 富集（免疫系统过程父项）：对每个细胞计算各 GO term 的 AUCell 得分
    cat(sprintf("运行 GO 富集（parent=%s）...\n", ENRICH_GO_PARENT))
    obj <- tryCatch(
      SeuratExtend::GeneSetAnalysisGO(obj,         # GO 富集分析
        parent = ENRICH_GO_PARENT,                 # GO 父项（免疫过程）
        nCores = ENRICHMENT_CORES),                # 并行核心
      error = function(e) {                        # 失败时
        warning("GO 富集失败: ", conditionMessage(e)) # 警告
        NULL                                       # 返回 NULL
      })
    if (!is.null(obj)) {                           # 成功时
      go_res <- obj@misc$AUCell$GO[[ENRICH_GO_PARENT]] # 提取富集结果
      # 结果可能是 data.frame/列表：尝试保存为 CSV，失败则打印结构
      tryCatch({
        if (is.data.frame(go_res)) {               # data.frame 直接写
          write.csv(go_res, file.path(OUT_DIR, "07_enrichment_GO_immune.csv"), row.names = FALSE)
          cat("  已保存: 07_enrichment_GO_immune.csv\n")
        } else {                                   # 列表/其他结构时打印概览
          cat("GO 富集结果结构：\n")
          print(utils::str(go_res, max.level = 1))
          saveRDS(go_res, file.path(OUT_DIR, "07_enrichment_GO_immune.rds")) # 存为 rds 备用
          cat("  已保存: 07_enrichment_GO_immune.rds（原始结构）\n")
        }
      }, error = function(e) warning("GO 富集结果保存失败: ", conditionMessage(e)))
    }
  } else {                                         # 安装失败时
    warning("SeuratExtend 不可用，跳过富集分析（DE 结果已正常输出）")
  }
}

## ---- 6. 完成 ----

# 打印最终结果摘要
cat("\n==== 06 差异表达与富集完成 ====\n")
cat(sprintf("DE 分组: %s（%d 个组）\n", de_group, length(unique(obj@meta.data[[de_group]])))) # 分组信息
cat("\n输出文件清单（位于 output/ 目录）：\n")
cat("  01_markers_by_celltype.csv            - 各细胞类型 marker 基因表\n")
cat("  02_DE_ypN0_vs_ypNplus.csv             - ypN0 vs ypN+ 总体差异基因表\n")
cat("  03_DE_by_celltype_group.csv           - 按细胞类型的组间差异表（可选）\n")
cat("  04_marker_dotplot.pdf                 - top marker 点图\n")
cat("  05_volcano_ypN0_vs_ypNplus.pdf        - 分组比较火山图\n")
cat("  06_marker_heatmap.pdf                 - top marker 表达热图\n")
cat("  07_enrichment_GO_immune.csv/.rds      - GO 免疫过程富集结果（可选）\n")
cat("\n06_differential_expression 完成。\n")
