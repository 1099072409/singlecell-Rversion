# 校验 04_CD45_integrated.rds 完整性与内部一致性
# 用途：确认最终落盘对象可正常读取、聚类数与 CSV 一致、降维完整
cat("正在加载 04_CD45_integrated.rds（约需 2~5 分钟）...\n")
t0 <- Sys.time()
obj <- readRDS("A:/Workbuddy/singlecell/Rversion/04_integration_clustering/output/04_CD45_integrated.rds")
cat(sprintf("加载耗时: %.1f 分钟\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat(sprintf("维度: %d 基因 x %d 细胞\n", nrow(obj), ncol(obj)))
cat(sprintf("样本数: %d\n", length(unique(obj$sample))))
cat("分组分布:\n"); print(table(obj$group, useNA = "ifany"))
cat("降维: ", paste(names(obj@reductions), collapse = ", "), "\n")
cat(sprintf("res=0.4 簇数: %d\n", length(unique(obj$cluster_0.4))))
cat(sprintf("主 cluster 列簇数: %d\n", length(unique(obj$cluster))))
cat("主 cluster 分布:\n"); print(sort(table(obj$cluster), decreasing = TRUE))
# 与 01_clustering_summary.csv 比对
csv <- read.csv("A:/Workbuddy/singlecell/Rversion/04_integration_clustering/output/01_clustering_summary.csv", stringsAsFactors = FALSE)
cat(sprintf("CSV 行数: %d（不含表头 %d 个簇）\n", nrow(csv), length(unique(csv$cluster))))
tab <- sort(table(obj$cluster), decreasing = TRUE)
csv_tab <- setNames(csv$n_cells, as.character(csv$cluster))
idx <- intersect(names(tab), names(csv_tab))
if (length(idx) > 0 && all(tab[idx] == csv_tab[idx])) {
  cat("一致性检查: rds 主 cluster 细胞数与 CSV 完全一致 ✔\n")
} else {
  cat("一致性检查: rds 与 CSV 不一致 ✘（见上）\n")
}
cat("VERIFY_DONE\n")
