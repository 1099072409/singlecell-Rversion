# manual_markers.R
# 从 细胞注释.md 整理出的经典 marker 字典，供 05_celltype_annotation.R 点图使用。
# 已校正拼写（FCG3RA->FCGR3A, TYOBP->TYROBP, RF4->XBP1），并统一为人类 HGNC 符号。
# 缺失基因会在运行时被过滤，不影响出图。

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
