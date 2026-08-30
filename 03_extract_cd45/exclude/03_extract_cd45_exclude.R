# =====================================================================
# 03_extract_cd45_exclude.R
# 单细胞数据分析流程 - 第三步（提取 CD45+ 免疫细胞）【exclude 变体】
#   与 03_extract_cd45.R 逻辑、参数、输出完全一致，仅在读入对象后剔除
#   EXCLUDE_SAMPLES 指定的样本，结果输出到本脚本所在目录（Rversion/03_extract_cd45/exclude/）。
#
# 功能：
#   1) 读取第二步质控后的 Seurat 对象（02_seurat_qc.rds）
#   2) 基于 PTPRC（CD45）基因的 counts 表达量判定 CD45+ 细胞：
#      PTPRC 的 counts > 0（检测到转录本）即判定为 CD45+（阈值可调）
#   3) 对"提取前（全部细胞）"与"提取后（CD45+ 细胞）"分别输出 QC 图，
#      沿用 02 步的 4 张图（小提琴图、散点图、提取后小提琴图、
#      按样本细胞数柱状图），并补充两张图：
#         - PTPRC 表达分布直方图（展示判定阈值依据）
#         - 各样本 CD45+ 比例柱状图（按分组着色）
#   4) 输出 CD45+ 细胞的 Seurat 对象（03_CD45_positive.rds）及统计表
#      （CD45_summary_by_sample.csv），并打印分组层面的汇总
#
# 用法：
#   在 RStudio 中打开本文件并点击 Source，或在命令行执行：
#     "C:/Program Files/R/R-4.4.3/bin/Rscript.exe" 03_extract_cd45.R
#   说明：所有路径按脚本自身所在目录自动推算，可从任意工作目录运行。
#
# 常用修改点：见下方 ## ---- 0. 配置 ---- 中的 CONFIG 参数。
# 参考：Rversion/02_quality_control/02_quality_control.R（脚本风格与 QC 图范式）
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
# 本脚本位于 Rversion/03_extract_cd45/exclude/ 下：
#   - 输入：向上两级回到 Rversion，再进入 02_quality_control/output
#   - 输出：脚本自身所在目录（即 exclude 文件夹，数据与本脚本同处一室）
INPUT_RDS <- file.path(script_dir, "..", "..", "02_quality_control", "output", "02_seurat_qc.rds") # 第二步质控后的 rds
OUT_DIR   <- script_dir                                                                       # 输出到 exclude 文件夹（脚本所在目录）

# ---- CD45+ 判定参数 ----
CD45_GENE      <- "PTPRC"   # CD45 标志基因名（人类；小鼠为 Ptprc，如切换物种请修改）
CD45_MIN_COUNTS <- 0        # PTPRC counts 大于该值即判定为 CD45+（0 = 检测到转录本即阳性）
# 提示：阈值设定遵循"参数选择先有结果"的规范——首次运行后查看
#       output/05_PTPRC_expression_distribution.pdf（PTPRC 表达分布直方图）与
#       output/CD45_summary_by_sample.csv（各样本 CD45+ 比例），如发现背景噪音
#       （低表达细胞过多）再回来调高 CD45_MIN_COUNTS（如 1 或 2）重新提取。

# ---- 其他 ----
SEED <- 123                 # 随机种子（保持可复现）

# ---- 剔除样本（exclude 分析专属配置）----
# 在本分析中需从样本集合中剔除以下样本（按 meta.data 的 sample_id 匹配）。
# 注：用户给定标签 "ypN07、011、012、013、014"，实际 sample_id 为
#     ypN07 / ypN011 / ypN012 / ypN013 / ypN014（对象中无字面 "011" 等标签）。
EXCLUDE_SAMPLES <- c("ypN07", "ypN011", "ypN012", "ypN013", "ypN014")

## ---- 1. 依赖检查与加载 ----

# 通用依赖安装函数：ensure_pkg("包名")，未安装时自动从 CRAN 安装
ensure_pkg <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) return(TRUE)             # 已装则跳过
  cat(sprintf("  正在安装依赖包: %s ...\n", pkg))                    # 提示开始安装
  install.packages(pkg, repos = "https://cloud.r-project.org")        # 从官方镜像安装
  requireNamespace(pkg, quietly = TRUE)                               # 返回是否安装成功
}

# 依次确保需要的包可用（Seurat/dplyr/ggplot2/patchwork 均已安装，patchwork 为图组合依赖）
for (p in c("Seurat", "dplyr", "ggplot2", "patchwork")) ensure_pkg(p)

# 加载核心包
library(Seurat)    # 单细胞分析主包
library(dplyr)     # 数据整理（group_by/summarise 等）
library(ggplot2)   # 绘图
library(patchwork) # 图组合

# 创建输出目录（已存在则不报错）
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 打印运行环境信息，便于排查问题
cat(sprintf("脚本目录: %s\n", script_dir))                    # 脚本所在目录
cat(sprintf("工作目录: %s\n", getwd()))                       # 当前工作目录
cat(sprintf("输出目录: %s\n", normalizePath(OUT_DIR)))        # 输出目录

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

# Seurat v5 并行/全局对象大小限制（为后续大对象操作预留空间）
options(future.globals.maxSize = 8 * 1024^3)

## ---- 2. 读入对象并检查 CD45 基因 ----

# 2.1 检查输入文件是否存在
if (!file.exists(INPUT_RDS)) stop("找不到输入 rds: ", INPUT_RDS)

# 2.2 读取第二步质控后的 Seurat 对象（451MB 左右，读取需 1~3 分钟）
cat("正在读取 02_seurat_qc.rds（约需 1~3 分钟）...\n")
obj <- readRDS(INPUT_RDS)                          # 读入 Seurat 对象
cat(sprintf("读入完成：%d 个基因 x %d 个细胞\n", nrow(obj), ncol(obj))) # 打印维度

# 2.3 确保默认 assay 为 RNA（后续操作都基于它）
if (DefaultAssay(obj) != "RNA") {                  # 若默认 assay 不是 RNA
  warning("默认 assay 不是 RNA，已切换为 RNA")      # 提示后切换
  DefaultAssay(obj) <- "RNA"
}

# 2.4 统一样本标识列：obj$sample = 样本名（02 产物通常已含此列，此处兜底，缺失时用 orig.ident 生成）
if (!"sample" %in% colnames(obj@meta.data)) {      # 输入对象没有 sample 列时
  obj$sample <- as.character(obj$orig.ident)       # 用 orig.ident（样本文件夹名）生成
}

# 2.5 统一 Layer 结构：若存在多个 layer（merge 产物），JoinLayers 合并为单一 counts
#     （02 步已处理过，此处作为兜底，保证 LayerData 取数正确）
cat(sprintf("当前 RNA layer: %s\n", paste(Layers(obj, assay = "RNA"), collapse = ", ")))
if (length(Layers(obj, assay = "RNA")) > 1) {      # 若存在多个 layer
  cat("检测到多个 layer，执行 JoinLayers 合并...\n")
  obj <- JoinLayers(obj, assay = "RNA")            # 合并所有 layer
}

# 2.6 检查 CD45 标志基因 PTPRC 是否存在于基因列表中（不存在则终止，避免后续下标越界）
if (!CD45_GENE %in% rownames(obj)) {               # 基因名不在行名中
  hit <- grep("PTPRC", rownames(obj), value = TRUE, ignore.case = TRUE) # 尝试大小写不敏感查找
  if (length(hit) > 0) {                           # 找到近似名（如大小写不同）
    CD45_GENE <- hit[1]                            # 采用找到的基因名
    cat(sprintf("提示：基因名大小写与配置不同，已改用 '%s'\n", CD45_GENE))
  } else {
    stop("数据中找不到 CD45 标志基因 PTPRC，无法提取 CD45+ 细胞") # 完全找不到则终止
  }
}
cat(sprintf("CD45 标志基因: %s\n", CD45_GENE))     # 打印实际使用的基因名

## ---- 2.7 剔除指定样本（exclude 分析专属步骤）----
# 依据 EXCLUDE_SAMPLES 配置，从 meta.data 的 sample_id 列中剔除对应样本的全部细胞。
# 仅当指定的样本 ID 全部存在于对象中时才执行，否则中止并报错，避免静默误删。
if (length(EXCLUDE_SAMPLES) > 0) {
  if (!"sample_id" %in% colnames(obj@meta.data)) {  # 缺少 sample_id 列时无法匹配
    stop("对象 meta.data 中找不到 sample_id 列，无法按样本剔除。请检查 02 产物。")
  }
  all_sid <- sort(unique(obj$sample_id))            # 对象中实际存在的全部 sample_id
  missing  <- setdiff(EXCLUDE_SAMPLES, all_sid)     # 配置中存在但对象中不存在的 ID
  if (length(missing) > 0) {                        # 存在不存在的 ID → 中止
    stop(sprintf("剔除列表中存在对象中不存在的样本 ID: %s\n  对象中已有的 sample_id: %s",
                 paste(missing, collapse = ", "), paste(all_sid, collapse = ", ")))
  }
  cat(sprintf("\n==== 2.7 剔除样本（exclude 分析）====\n"))
  cat(sprintf("剔除前细胞数: %d（共 %d 个样本）\n", ncol(obj), length(all_sid)))
  cat(sprintf("待剔除样本（%d 个）: %s\n", length(EXCLUDE_SAMPLES), paste(EXCLUDE_SAMPLES, collapse = ", ")))
  keep_cells <- colnames(obj)[!(obj$sample_id %in% EXCLUDE_SAMPLES)]  # 保留细胞（不在剔除列表）
  obj <- subset(obj, cells = keep_cells)            # 仅删除细胞，保留基因/层结构/meta 列
  cat(sprintf("剔除后细胞数: %d（剩余 %d 个样本）\n", ncol(obj), length(setdiff(all_sid, EXCLUDE_SAMPLES))))
  cat(sprintf("保留样本: %s\n", paste(setdiff(all_sid, EXCLUDE_SAMPLES), collapse = ", ")))
}

## ---- 3. 判定 CD45+ 细胞 ----

# 3.1 提取 PTPRC 基因的 counts 表达量（稀疏矩阵按行取子集后转稠密向量）
counts_mat <- LayerData(obj, assay = "RNA", layer = "counts") # 取整个 counts 矩阵（dgCMatrix）
ptprc_counts <- as.numeric(counts_mat[CD45_GENE, ])          # 该基因在所有细胞中的 counts 向量
rm(counts_mat); gc()                                          # 及时释放大矩阵，节省内存

# 3.2 写回 meta.data：保留原始 counts 值与 CD45 阳性标记
obj$CD45_counts   <- ptprc_counts                             # PTPRC 的原始 counts（便于追溯/后续调阈值）
obj$CD45_positive <- ptprc_counts > CD45_MIN_COUNTS           # 判定：counts > 阈值 即 CD45+（逻辑向量）

# 3.3 打印 CD45+ 细胞总体统计
cat("\n==== CD45+ 判定结果 ====\n")
cat(sprintf("全部细胞: %d\n", ncol(obj)))                     # 提取前细胞数
cat(sprintf("CD45+ 细胞: %d（%.2f%%）\n",                     # 阳性细胞数与比例
            sum(obj$CD45_positive), 100 * mean(obj$CD45_positive)))
cat(sprintf("CD45- 细胞: %d\n", sum(!obj$CD45_positive)))     # 阴性细胞数

# 3.3b 打印 PTPRC 表达量的分布统计（供判断阈值合理性）：
#      若 CD45+ 细胞中大量细胞仅 1~2 个 PTPRC 转录本（背景噪音），可考虑调高 CD45_MIN_COUNTS
cat("PTPRC 表达量分布（仅对 CD45+ 细胞）：\n")                # 打印提示
cat(sprintf("  counts 分位数 [5%%, 25%%, 50%%, 75%%, 95%%]: %s\n",  # 分位数
            paste(quantile(ptprc_counts[obj$CD45_positive], c(0.05, 0.25, 0.5, 0.75, 0.95)), collapse = ", ")))
cat(sprintf("  counts==1 的 CD45+ 细胞占比: %.2f%%\n",        # 单转录本细胞占比（背景指标）
            100 * mean(ptprc_counts[obj$CD45_positive] == 1)))

# 3.4 分组（ypN0/ypN+）层面的 CD45+ 汇总，供快速查看
grp_summary <- obj@meta.data %>%
  dplyr::group_by(group) %>%                                     # 按分组统计
  dplyr::summarise(
    cells_total  = dplyr::n(),                                   # 总细胞数
    cells_CD45pos = sum(CD45_positive),                          # CD45+ 细胞数
    pct_CD45pos  = 100 * mean(CD45_positive),                    # CD45+ 比例(%)
    .groups      = "drop"
  )
cat("按分组汇总：\n")
print(as.data.frame(grp_summary))                                # 打印分组层面结果

## ---- 4. 按样本统计 CD45+ 细胞数并写 CSV ----

# 按样本（sample）+ 分组（group）+ 样本编号（sample_id）汇总 CD45+ 提取前后细胞数
cd45_summary <- obj@meta.data %>%
  dplyr::group_by(sample, group, sample_id) %>%                  # 按三列分组
  dplyr::summarise(
    cells_total   = dplyr::n(),                                  # 提取前细胞数（=CD45- + CD45+）
    cells_CD45pos = sum(CD45_positive),                          # CD45+ 细胞数
    pct_CD45pos   = 100 * mean(CD45_positive),                   # CD45+ 比例(%)
    .groups       = "drop"
  )
write.csv(cd45_summary, file.path(OUT_DIR, "CD45_summary_by_sample.csv"), row.names = FALSE) # 写统计表
cat("  已保存: CD45_summary_by_sample.csv\n")                    # 提示保存成功

## ---- 5. QC 图（提取前：全部细胞）----

# 通用保存 PDF 函数：把 ggplot 对象 p 输出到 OUT_DIR 下的 file，尺寸 w×h 英寸
save_pdf <- function(p, file, w, h) {
  pdf(file.path(OUT_DIR, file), width = w, height = h)          # 打开 PDF 设备
  print(p)                                                       # 打印图形
  dev.off()                                                      # 关闭设备（必须，否则文件不完整）
  cat(sprintf("  已保存图表: %s\n", file))                       # 提示保存成功
}

cat("\n==== 5. 绘制 QC 图（提取前：全部细胞）====\n")

# 5.1 提取前小提琴图：5 个 QC 指标在各样本中的分布（与 02 步同款）
p_vln_before <- VlnPlot(obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), # 5 个指标
  group.by = "sample",                          # 按样本分组
  pt.size = 0,                                  # 不画点，只看分布
  ncol = 5) +                                   # 5 列排布
  theme(axis.text.x = element_text(angle = 45, hjust = 1))       # x 轴标签旋转 45° 防重叠
save_pdf(p_vln_before, "01_QC_violin_before_filter.pdf", 18, 5)

# 5.2 提取前散点图：总计数 vs 基因数（与 02 步同款，直观显示细胞质量分布）
p_scatter <- FeatureScatter(obj,
  feature1 = "nCount_RNA",                      # x 轴：总 UMI 计数
  feature2 = "nFeature_RNA",                    # y 轴：检测基因数
  group.by = "sample")                          # 按样本着色
save_pdf(p_scatter, "02_QC_scatter_before_filter.pdf", 8, 6)

## ---- 6. 提取 CD45+ 细胞并保存 ----

# 6.1 按 CD45_positive 标记提取 CD45+ 细胞（subset 只删细胞，不改基因、不改 layer 结构）
obj_cd45 <- subset(obj, subset = CD45_positive == TRUE)
cat(sprintf("\n提取完成：%d 个细胞 -> %d 个 CD45+ 细胞（占比 %.2f%%）\n",
            ncol(obj), ncol(obj_cd45), 100 * ncol(obj_cd45) / ncol(obj))) # 打印提取结果

# 6.2 保存 CD45+ 细胞对象（meta.data 完整保留 sample/group/sample_id/CD45_counts/CD45_positive 等信息）
saveRDS(obj_cd45, file.path(OUT_DIR, "03_CD45_positive.rds"))
cat(sprintf("  已保存: %s\n", file.path(OUT_DIR, "03_CD45_positive.rds")))

## ---- 7. QC 图（提取后：CD45+ 细胞）----

cat("\n==== 7. 绘制 QC 图（提取后：CD45+ 细胞）====\n")

# 7.1 提取后小提琴图：3 个核心 QC 指标在 CD45+ 细胞中的分布（与 02 步同款）
p_vln_after <- VlnPlot(obj_cd45,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),     # 3 个核心指标
  group.by = "sample",                          # 按样本分组
  pt.size = 0,                                  # 不画点
  ncol = 3) +                                   # 3 列排布
  theme(axis.text.x = element_text(angle = 45, hjust = 1))      # 旋转标签
save_pdf(p_vln_after, "03_QC_violin_after_filter.pdf", 14, 5)

# 7.2 提取前后细胞数柱状图：按样本编号展示 before(全部)/after(CD45+) 对比（与 02 步同款思路）
sample_id_levels <- unique(cd45_summary$sample_id)               # 用汇总表中的样本编号顺序
cell_count_long <- rbind(
  data.frame(sample_id = cd45_summary$sample_id, group = cd45_summary$group,
             status = "before", count = cd45_summary$cells_total, stringsAsFactors = FALSE),  # 提取前
  data.frame(sample_id = cd45_summary$sample_id, group = cd45_summary$group,
             status = "after",  count = cd45_summary$cells_CD45pos, stringsAsFactors = FALSE)  # 提取后(CD45+)
)
cell_count_long$sample_id <- factor(cell_count_long$sample_id, levels = sample_id_levels)      # 固定 x 轴顺序
p_bar <- ggplot(cell_count_long, aes(x = sample_id, y = count, fill = status)) + # 柱状图
  geom_col(position = "dodge", width = 0.7) +    # 分组柱状图
  scale_fill_manual(values = c(before = "grey60", after = "#E64B35")) + # 灰色=全部细胞，红色=CD45+（提取）
  theme_minimal() +                              # 简洁主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +        # 旋转 x 轴标签
  labs(title = "Cell count before/after CD45+ extraction by sample", # 图标题（英文避免字体问题）
       x = "Sample ID", y = "Cell count", fill = "Status")          # 坐标轴与图例标签
save_pdf(p_bar, "04_QC_cell_count_by_sample.pdf", 14, 6)

## ---- 8. 补充图 ----

# 8.1 PTPRC 表达分布直方图：展示 CD45 判定的依据
#     横轴为 log10(counts+1)，红色虚线为判定阈值（counts>0 即 log10(1)=0 右侧为阳性）
cat("\n==== 8. 绘制补充图 ====\n")
p_ptprc <- ggplot(obj@meta.data, aes(x = log10(CD45_counts + 1))) +  # 全部细胞的 PTPRC 表达分布
  geom_histogram(bins = 60, fill = "steelblue", color = "white") +   # 直方图
  geom_vline(xintercept = log10(CD45_MIN_COUNTS + 1),                # 阈值竖线位置
             linetype = "dashed", color = "red", linewidth = 0.8) +  # 红色虚线
  annotate("text", x = log10(CD45_MIN_COUNTS + 1) + 0.15,            # 阈值标注文字
           y = Inf, label = paste0("CD45+ threshold (counts>", CD45_MIN_COUNTS, ")"),
           vjust = 1.5, hjust = 0, color = "red", size = 3.5) +      # 文字样式
  theme_minimal() +                                                  # 简洁主题
  labs(title = "PTPRC (CD45) expression distribution",               # 图标题
       x = "log10(PTPRC counts + 1)", y = "Cell count")              # 坐标轴标签
save_pdf(p_ptprc, "05_PTPRC_expression_distribution.pdf", 9, 5)

# 8.2 各样本 CD45+ 比例柱状图：按分组着色，直观比较不同样本的免疫细胞占比
p_prop <- ggplot(cd45_summary, aes(x = sample_id, y = pct_CD45pos, fill = group)) + # 柱状图
  geom_col(width = 0.7) +                                            # 柱子
  scale_fill_manual(values = c(ypN0 = "#4DBBD5", `ypN+` = "#E64B35")) + # ypN0=蓝, ypN+=红
  theme_minimal() +                                                  # 简洁主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +         # 旋转 x 轴标签
  labs(title = "CD45+ cell proportion by sample",                    # 图标题
       x = "Sample ID", y = "CD45+ proportion (%)", fill = "Group")  # 坐标轴与图例标签
save_pdf(p_prop, "06_CD45_proportion_by_sample.pdf", 14, 6)

# 8.3 每个样本 CD45+ / CD45- 堆叠柱状图：直接展示各样本中免疫细胞的阴阳分层
#     柱高 = 该样本全部细胞数，红段 = CD45+（免疫细胞），灰段 = CD45-（非免疫细胞），
#     一眼看出每个样本里免疫浸润的绝对规模及其阴阳构成（与 06 图互补：06 看比例，本图看绝对数）
cd45_stack_long <- rbind(
  data.frame(sample_id = cd45_summary$sample_id, group = cd45_summary$group,
             status = "CD45-", count = cd45_summary$cells_total - cd45_summary$cells_CD45pos,
             stringsAsFactors = FALSE),   # CD45- 细胞数（全部 - 阳性）
  data.frame(sample_id = cd45_summary$sample_id, group = cd45_summary$group,
             status = "CD45+", count = cd45_summary$cells_CD45pos,
             stringsAsFactors = FALSE)    # CD45+ 细胞数（免疫细胞）
)
cd45_stack_long$sample_id <- factor(cd45_stack_long$sample_id, levels = sample_id_levels) # 固定 x 轴顺序
p_stack <- ggplot(cd45_stack_long, aes(x = sample_id, y = count, fill = status)) +  # 堆叠柱状图
  geom_col(width = 0.7) +                                                 # 默认 position="stack" 堆叠
  scale_fill_manual(values = c("CD45-" = "#90A4AE", "CD45+" = "#E64B35")) + # 灰蓝=阴性，红=阳性
  theme_minimal() +                                                       # 简洁主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +              # 旋转 x 轴标签防重叠
  labs(title = "CD45+ / CD45- cell count by sample (stacked)",            # 图标题
       x = "Sample ID", y = "Cell count", fill = "CD45 status")           # 坐标轴与图例标签
save_pdf(p_stack, "07_CD45_pos_neg_stack_by_sample.pdf", 14, 6)

## ---- 9. 完成 ----

# 打印最终结果摘要
cat("\n==== CD45+ 提取完成 ====\n")
cat(sprintf("全部细胞: %d\n", ncol(obj)))                            # 提取前
cat(sprintf("CD45+ 细胞: %d（%.2f%%）\n", ncol(obj_cd45), 100 * ncol(obj_cd45) / ncol(obj))) # 提取后
cat("\n输出文件清单（位于 exclude/ 目录）：\n")
cat("  03_CD45_positive.rds                  - CD45+ 细胞 Seurat 对象\n")
cat("  CD45_summary_by_sample.csv            - 按样本/分组/编号的 CD45+ 统计表\n")
cat("  01_QC_violin_before_filter.pdf        - 提取前（全部细胞）QC 小提琴图\n")
cat("  02_QC_scatter_before_filter.pdf       - 提取前计数-基因数散点图\n")
cat("  03_QC_violin_after_filter.pdf         - 提取后（CD45+ 细胞）QC 小提琴图\n")
cat("  04_QC_cell_count_by_sample.pdf        - 提取前后细胞数柱状图\n")
cat("  05_PTPRC_expression_distribution.pdf  - PTPRC 表达分布直方图（判定阈值依据）\n")
cat("  06_CD45_proportion_by_sample.pdf      - 各样本 CD45+ 比例柱状图\n")
cat("  07_CD45_pos_neg_stack_by_sample.pdf    - 各样本 CD45+ / CD45- 堆叠柱状图\n")
cat("\n03_extract_cd45 完成。\n")

