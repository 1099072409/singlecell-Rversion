# manual_markers.R
# 手写经典 marker 字典，仅供 05_celltype_annotation.R 的 08_marker_dotplot_manual.pdf
# 点图使用（用于"验证手动注释对不对"），不参与任何自动注释判定。
#
# 2026-09-02 按最新手写注释体系重写：键名必须与 05 脚本里的 MANUAL_CELLTYPE_MAP
# 完全一致（T / NK / B / Plasma / Myeloid / Mast / T_proliferating /
# B_proliferating / Epithelial / Fibroblast / Endothelial），
# 否则该类型不会出现在点图里。
#
# 说明：
#   * Myeloid 是合并大类（单核 / 巨噬 / DC / 中性粒），marker 取各亚群代表性基因；
#   * T_proliferating / B_proliferating 是增殖状态而非独立谱系，
#     marker = 细胞周期基因 + 各自谱系基因（用于区分增殖细胞的来源）；
#   * 数据中不存在的基因会在运行时被自动过滤，不影响出图。

MANUAL_MARKERS <- list(
  "T"               = c("CD3D", "CD3E", "CD3G", "TRAC", "CD2", "IL7R"),
  "NK"              = c("NKG7", "GNLY", "KLRD1", "PRF1", "FCGR3A", "NCAM1"),
  "B"               = c("CD19", "MS4A1", "CD79A", "CD79B", "IGHD"),
  "Plasma"          = c("MZB1", "SDC1", "JCHAIN", "PRDM1", "IGHG1", "XBP1"),
  # 注意：FCGR3A(CD16) 已归到 NK，Myeloid 不再重复列出，否则 DotPlot2 会去重并告警
  "Myeloid"         = c("LYZ", "CST3", "CD68", "CD14", "ITGAM", "CSF1R",
                        "S100A8", "S100A9"),
  "Mast"            = c("TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"),
  "T_proliferating" = c("MKI67", "TOP2A", "PCNA", "CDK1", "TYMS", "STMN1",
                        "CD3D"),
  "B_proliferating" = c("MKI67", "TOP2A", "PCNA", "CDK1", "TYMS", "STMN1",
                        "MS4A1"),
  "Epithelial"      = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"),
  "Fibroblast"      = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "ACTA2"),
  "Endothelial"     = c("PECAM1", "VWF", "CDH5", "PLVAP", "CLDN5")
)
