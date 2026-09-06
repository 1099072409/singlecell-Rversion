library(Seurat)
library(SeuratExtend)
library(dplyr)

options(max.print = 50, spe = "human")

# figure 2 ----------------------------------------------------------------

?GeneSetAnalysisGO
pbmc <- GeneSetAnalysisGO(pbmc, parent = "immune_system_process")
matr <- pbmc@misc$AUCell$GO$immune_system_process
matr <- FilterGOTerms(matr, n.min = 5, change.name = T)
rownames(matr) <- gsub("antigen processing and presentation","Ag proc. and pres.", rownames(matr))

Heatmap(CalcStats(matr, f = pbmc$cluster, order = "p", n = 3), lab_fill = "zscore", plot.margin = margin(l = 20)) +
  theme(legend.position = "bottom")
# install.packages("svglite")
ggsave("figure v1/fig2 heatmap.svg", width = 8.5, height = 5.5)

# WaterfallPlot(matr, f = pbmc$cluster, ident.1 = "B cell", ident.2 = "CD8 T cell", top.n = 10, color = "p")

# pbmc <- GeneSetAnalysis(pbmc, genesets = hall50$human)
# matr <- pbmc@misc$AUCell$genesets
# Heatmap(CalcStats(matr, f = pbmc$cluster), lab_fill = "zscore")

paths <- c("GO:0030183")
cells <- colnames(pbmc)[pbmc$cluster %in% c("B cell","CD4 T Memory","CD8 T cell","DC","Mono CD14")]
matr <- pbmc@misc$AUCell$GO$immune_system_process[paths,cells,drop=F]
matr <- RenameGO(matr)
VlnPlot2(matr, f = pbmc$cluster[cells], ncol = 1, stat.method = "wilcox.test")
ggsave("figure v1/fig2 violin.svg", width = 4, height = 4)

GSEAplot(
  pbmc, ident.1 = "CD4 T Naive", title = "Hallmark_IFN_gamma_Response",
  geneset = hall50$human$HALLMARK_INTERFERON_GAMMA_RESPONSE)
ggsave("figure v1/fig2 gsea.svg", width = 5, height = 4)

# fig 3 -------------------------------------------------------------------

mye_small <- readRDS("~/Documents/SeuratExtend-pbmc/upload/pbmc10k_mye_small_velocyto.rds")
loom_path <- "~/Documents/SeuratExtend-pbmc/upload/pbmc10k_mye_small.loom"
scenic_loom_path <- "~/Documents/SeuratExtend-pbmc/upload/pbmc3k_small_pyscenic_integrated-output.loom"

# Set up the path for saving the AnnData object in the HDF5 (h5ad) format
adata_path <- file.path(tempdir(), "mye_small.h5ad")

# Integrate Seurat Object and velocyto loom information into one AnnData object.
# This object will be stored at the specified path.
scVelo.SeuratToAnndata(
  mye_small, # The downloaded example Seurat object
  filename = adata_path, # Path where the AnnData object will be saved
  velocyto.loompath = loom_path, # Path to the loom file
  prefix = "sample1_", # Prefix for cell IDs in the Seurat object
  postfix = "-1" # Postfix for cell IDs in the Seurat object
)
scVelo.Plot(color = "cluster", save = "umap1.png", figsize = c(5,4), palette = color_pro(3, 2))

mye_small <- Palantir.RunDM(mye_small)
DimPlot2(mye_small, reduction = "ms", group.by = "cluster", label = FALSE, theme = NoAxes(), cols = "light")
ggsave("figure v1/fig3 dm.svg", height = 4, width = 5)

mye_small <- Palantir.Pseudotime(mye_small, start_cell = "sample1_GAGAGGTAGCAGTACG-1")
ps <- mye_small@misc$Palantir$Pseudotime
head(ps)
colnames(ps)[3:4] <- c("fate1", "fate2")
mye_small@meta.data[,colnames(ps)] <- ps
DimPlot2(mye_small, features = colnames(ps), reduction = "ms",
         cols = list(Entropy = "D"), theme = NoAxes())
ggsave("figure v1/fig3 dm ps.svg", height = 5, width = 6)

GeneTrendCurve.Palantir(mye_small, features = c("CD14", "FCGR3A"), pseudotime.data = ps)
ggsave("figure v1/fig3 ps curve.svg", height = 3, width = 6)

GeneTrendHeatmap.Palantir(
  mye_small,
  features = c("CD14", VariableFeatures(mye_small)[1:10]),
  pseudotime.data = ps,
  magic = FALSE,
  lineage = "fate1"
)
ggsave("figure v1/fig3 ps heatmap.png", height = 3, width = 6)

mye_small <- Palantir.Magic(mye_small)

# Creates a new assay "magic" in the Seurat object
mye_small <- NormalizeData(mye_small)
DimPlot2(mye_small, features = c("CD14", "magic_CD14", "FLT3", "magic_FLT3"), theme = NoAxes())
ggsave("figure v1/fig3 magic.svg", height = 5, width = 6)

pbmc <- ImportPyscenicLoom(scenic_loom_path, seu = pbmc)
tf_auc <- pbmc@misc$SCENIC$RegulonsAUC

DimPlot2(
  pbmc,
  features = c("ETS1", "ATF3", "tf_ETS1", "tf_ATF3"),
  cols = list("tf_ETS1" = "D", "tf_ATF3" = "D"),
  theme = NoAxes()
)
ggsave("figure v1/fig3 scenic umap.svg", height = 5, width = 6)

DefaultAssay(pbmc) <- "TF"

# Creating a waterfall plot to compare regulon activity between monocytes and CD8 T cells
WaterfallPlot(
  pbmc,
  features = rownames(pbmc),  # Using all available TFs in the "TF" assay
  ident.1 = "Mono CD14",      # First group of cells
  ident.2 = "CD8 T cell",     # Second group of cells
  exp.transform = FALSE,      # Disable transformation of expression data
  top.n = 20                  # Display the top 20 most differentially active TFs
)
ggsave("figure v1/fig3 scenic waterfall.svg", height = 4, width = 6.5)

# SCENIC2

SeuratExtend::Seu2Loom(mye_small, "mye_small.loom")
RunScenic("mye_small.loom", spe = "human")

setwd("~/Documents/SeuratExtend-manuscript")
mye_small <- ImportPyscenicLoom("Scenic_project/output/pyscenic_integrated-output.loom", seu = mye_small)

DefaultAssay(mye_small) <- "TF"
WaterfallPlot(
  mye_small,
  features = rownames(mye_small),  # Using all available TFs in the "TF" assay
  ident.1 = "Mono CD14",      # First group of cells
  ident.2 = "Mono FCGR3A",     # Second group of cells
  exp.transform = FALSE,      # Disable transformation of expression data
  top.n = 20                  # Display the top 20 most differentially active TFs
)
ggsave("figure v1/fig3 scenic waterfall2.svg", height = 4, width = 6.5)

DimPlot2(
  mye_small,
  features = c("CREB5", "magic_CREB5", "tf_CREB5", "POU2F2", "magic_POU2F2", "tf_POU2F2"),
  cols = list("tf_CREB5" = "D", "tf_POU2F2" = "D"),
  theme = NoAxes(), pt.size = 0.8
)
ggsave("figure v1/fig3 magic+scenic.svg", height = 5, width = 9)

# Cytotrace

DefaultAssay(mye_small) <- "RNA"
gs <- mye_small@misc$SCENIC$Regulons
vg <- VariableFeatures(mye_small)[1:1000]
gs <- lapply(gs, function(x) intersect(x, vg))
gs <- gs[lengths(gs) > 0]
toptf <- CalcStats(GetAssayData(mye_small, assay = "TF"), f = mye_small$cluster, order = "p", n = 5)
toptf <- toptf[6:15,]
gs <- gs[names(gs) %in% rownames(toptf)]
net <- melt(gs)
write.table(net, "network.csv", sep = ",", quote = F, row.names = F)

tf <- unique(net$L1)
node <- data.frame(
  nodename = unlist(net) %>% unique(),
  type = "gene"
)
node$show.name <- node$nodename
node[node$nodename %in% tf, "type"] <- "TF"
node[!node$nodename %in% tf, "show.name"] <- ""

ge <- CalcStats(mye_small, features = unique(net$value), group.by = "cluster")
ge[,"nodename"] <- rownames(ge)
node <- left_join(node, ge, by = "nodename")
toptf[,"nodename"] <- rownames(toptf)
node[node$type == "TF",] <- left_join(node[node$type == "TF",1:3], toptf, by = "nodename")
write.table(node, "node.csv", sep = ",", quote = F, row.names = F)

# fig 4 -------------------------------------------------------------------

DimPlot2(pbmc, theme = NoAxes() + NoLegend(), cols = NULL)
library(hyc)
show_col2(list("color_pro" = color_pro(9)), ncol = 9)
p1 <- DimPlot2(pbmc, theme = NoAxes() + NoLegend())
p2 <- ClusterDistrBar(rep("C",ncol(pbmc)), pbmc$cluster) + NoAxes() + NoLegend()
ggarrange(p1, p2, ncol = 1, heights = c(3,1))

par_list <- list(
  "Indistinct" = pal_brewer(palette = "Greens")(9),
  "Overly Saturated" = c("#ccffaa","#c00bff","#cfdb00","#0147ee","#f67900","#1b002c","#00e748","#e30146","#ffb1e8"),
  "Lively" = c("#ff2026","#cf5d00","#ffd03f","#649f00","#a3f83d","#82cc58","#6645fe","#d8009c","#ff43a2"),
  "Default" = color_pro(9),
  "Light" = color_pro(9, 2),
  "Red" = color_pro(9, 3),
  "Yellow" = color_pro(9, 4),
  "Green" = color_pro(9, 5),
  "Blue" = color_pro(9, 6),
  "Purple" = color_pro(9, 7)
)

p_list <- list()
for (i in 1:10) {
  title <- names(par_list)[i]
  cols <- par_list[[i]]
  p1 <- DimPlot2(pbmc, cols = cols, theme = NoAxes() + NoLegend(), combine = F, pt.size = 0.8)
  p1 <- p1[[1]] + labs(title = title) + theme(title = element_text(size = 13))
  p2 <- show_col2(par_list[i], ncol = 9) + theme(strip.text.y.left = element_blank())
  pc <- plot_grid(p1, p2, ncol = 1, rel_heights = c(5,1))
  p_list[[i]] <- pc
}
p <- plot_grid(plotlist = p_list, ncol = 5)
ggsave("figure v1/fig4 dimplots.svg", plot = p, width = 12, height = 6)

FeaturePlot3.grid(pbmc, features = c("CD3D", "CD14", "CD79A"), pt.size = 0.5)
ggsave("figure v1/fig4 umap ryb.svg", width = 2.8, height = 3)
FeaturePlot3(pbmc, "CD3D", "CD14", "CD79A")
ggsave("figure v1/fig4 umap ryb legend.svg", width = 4, height = 5)
FeaturePlot3.grid(pbmc, features = c("CD3D", "CD14", "CD79A"), color = "rgb", pt.size = 1)
ggsave("figure v1/fig4 umap rgb.svg", width = 2.8, height = 3)
FeaturePlot3(pbmc, "CD3D", "CD14", "CD79A", color = "rgb")
ggsave("figure v1/fig4 umap rgb legend.svg", width = 4, height = 5)

# Fig 5 ---------

gc.mel <- readRDS("~/Documents/Essen/2023-8-8 melanoma integration/GC_BT_harmony_singlet_all_cells_annotated.rds")
gc.mel.go <- readRDS("~/Documents/Essen/2023-8-8 melanoma integration/2023-8-9 gc.mel.go.rds")

rn <- c(
  "Activated_CD8_Tcells" = "T/NK Cells",
  "CD4_Tcells" = "T/NK Cells",
  "Cycling_CD8_Tcells" = "T/NK Cells",
  "Cytotoxic_CD8_Tcells" = "T/NK Cells",
  "Dysfuntional_CD8_Tcells" = "T/NK Cells",
  "Memory_CD8_Tcells" = "T/NK Cells",
  "Tregs" = "T/NK Cells",
  "NK" = "T/NK Cells",

  "M2 Macrophages" = "Macrophages",
  "Macrophages_CXCL10" = "Macrophages",
  "Macrophages_FOS" = "Macrophages",
  "Macrophages_necrosis" = "Macrophages",
  "Macrophages_PLA2G2D" = "Macrophages",

  "Monocytes_CD14" = "Monocytes",
  "Monocytes_CD16" = "Monocytes",
  "Monocytes_TNFSF13" = "Monocytes",

  "B_cells" = "B Cells",
  "Plasma_cells" = "B Cells",

  "DC1" = "DC",
  "DC2" = "DC",
  "pDC" = "DC",

  "CAF_APCDD1" = "CAF",
  "CAF_CHI3L1" = "CAF",
  "CAF_PLA2G2A" = "CAF",
  "CAF_POSTN" = "CAF",

  "BEC" = "EC/Peri",
  "LEC" = "EC/Peri",
  "Pericyte" = "EC/Peri",

  "Melanocytic_OXPHOS" = "Malignant",
  "Melanoma_immune_like" = "Malignant",
  "Patient_specific_A" = "Malignant",
  "Patient_specific_B" = "Malignant",
  "Doublets" = "Malignant",
  "Mitochondrial" = "Malignant",
  "Mitotic" = "Malignant",
  "Neural_like" = "Malignant",
  "Stress(hypoxia)" = "Malignant",
  "Trans_reg" = "Malignant",
  "Interferon_alpha_beta" = "Malignant",
  "Antigen_presentation" = "Malignant",
  "Mesenchymal" = "Malignant",

  "Melanophages" = "Macrophages",

  "Mast_cells" = "Others"
)
gc.mel$cluster_main <- plyr::revalue(gc.mel$Detailed_cluster_all_5, rn)
gc.mel$cluster_main[is.na(gc.mel$cluster_main)] <- "Others"
# sort(unique(gc.mel$Detailed_cluster_all_5))
DimPlot2(gc.mel, features = "cluster_main", label = T, theme = NoAxes())
gc.mel$orig.ident <- stringr::str_extract(colnames(gc.mel), "(?<=_)[0-9]+")
gc.mel$patient <- paste0("P",sprintf("%02d", as.numeric(gc.mel$orig.ident)))
gc.mel <- RunUMAP(gc.mel, dims = 1:16, reduction.key = "UMAP2_", reduction.name = "umap2")
DimPlot2(gc.mel, label = T, reduction = "umap2", group.by = c("patient","cluster_main"), theme = NoAxes())

gc.mel.go$cluster_main <- plyr::revalue(gc.mel.go$Detailed_cluster_all_5, rn)
gc.mel.go$cluster_main[is.na(gc.mel.go$cluster_main)] <- "Others"
gc.mel.go$patient <- paste0("P",sprintf("%02d", as.numeric(gc.mel.go$orig.ident)))
DimPlot2(gc.mel.go, label = T, group.by = c("patient","cluster_main"), theme = NoAxes())

# tmp <- FindMarkers(gc.mel, ident.1 = "P40", group.by = "patient", logfc.threshold = 2)
# VlnPlot2(gc.mel, group.by = "patient", features = c("IGLV7-46"), pt = F)

# Subset
# gc.mel2 <- subset(gc.mel, cells = colnames(gc.mel)[
#   gc.mel$cluster_main != "Others" &
#     !gc.mel$patient %in% c("P43","P45")
# ])
gc.mel2 <- subset(gc.mel, cells = colnames(gc.mel)[
  gc.mel$cluster_main != "Others"
])

gc.mel2 <- RunUMAP(gc.mel2, dims = 1:16, reduction.key = "UMAP2_", reduction.name = "umap2")
DimPlot2(gc.mel2, reduction = "umap2", group.by = c("patient","cluster_main"), theme = NoAxes() + NoLegend(), cols = "light", label = T)
ggsave("figure v1/fig5 umap no harmony filtered.png", width = 9.5, height = 5)
# gc.mel2.old <- gc.mel2

# gc.mel.go2 <- subset(gc.mel.go, cells = colnames(gc.mel)[
#   gc.mel$cluster_main != "Others" &
#     !gc.mel$patient %in% c("P43","P45")
# ])
gc.mel.go2 <- subset(gc.mel.go, cells = colnames(gc.mel)[
  gc.mel$cluster_main != "Others"
])
gc.mel.go2 <- FindVariableFeatures(gc.mel.go2)
gc.mel.go2 <- ScaleData(gc.mel.go2)
gc.mel.go2 <- RunPCA(gc.mel.go2)
gc.mel.go2 <- RunUMAP(gc.mel.go2, dims = 1:16)
DimPlot2(gc.mel.go2, label = T, group.by = c("patient","cluster_main"), theme = NoAxes() + NoLegend())
ggsave("figure v1/fig5 umap GO no harmony filtered.png", width = 9.5, height = 5)
saveRDS(gc.mel.go2, file = "rds/2024-7-2 gc.mel.go2.rds")
# gc.mel.go2.old <- gc.mel.go2

mk.b <- FindMarkers(gc.mel.go2, ident.1 = "B Cells", group.by = "cluster_main", only.pos = T, logfc.threshold = 0.1)
library(dplyr)
mk.b <- arrange(mk.b, desc(avg_log2FC))
RenameGO(mk.b)
gs.mk.b <- head(rownames(mk.b))
RenameGO(gs.mk.b)
RenameGO(c("GO:0019815","GO:0042611"))
DimPlot2(gc.mel.go2, features = c("GO:0019815","GO:0042611"), ncol = 1, theme = NoAxes() + NoLegend())
ggsave("figure v1/fig5 umap GO feature.png", width = 4.5, height = 9)

# marker genes
pbmc <- GeneSetAnalysis(pbmc, genesets = PanglaoDB_data$marker_list_human, title = "panglao")
auc <- pbmc@misc$AUCell$panglao

df2 <- CalcStats(auc, f = pbmc$cluster, order = "p", n = 3)
pbmc$unnamed_clusters <- paste0("C",1:9)[pbmc$cluster]
df1 <- CalcStats(auc, f = pbmc$unnamed_clusters, order = "p", n = 5)
Heatmap(df1, lab_fill = "zscore", angle = 0, hjust = 0.5)
ggsave("figure v1/fig5 heatmap markers.svg", width = 5.5, height = 6)
DimPlot2(pbmc, features = c("unnamed_clusters","cluster"), cols = list(unnamed_clusters = "light"), ncol = 1,theme = NoAxes() + NoLegend(), label = T, box = T, label.color = "black", repel = T)
ggsave("figure v1/fig5 umap markers.svg", width = 3.5, height = 7.5)
Heatmap(df2, lab_fill = "zscore")

# fig5 added --------------------------------------------------------------

table(gc.mel.go@meta.data$Detailed_cluster_all_5)
DimPlot2(gc.mel.go2, group.by = "cluster_main")
DimPlot2(gc.mel.go2, features = "Detailed_cluster_all_5", cells = gc.mel.go2$cluster_main == "Malignant")
gc.mel.go2$malig_patient_specific <- gc.mel.go2$Detailed_cluster_all_5
gc.mel.go2$malig_patient_specific[!gc.mel.go2$malig_patient_specific %in% c("Patient_specific_A", "Patient_specific_B")] <- "Others"
cells_malig <- CellSelector(DimPlot(gc.mel.go2))
cells_malig <- intersect(cells_malig, colnames(gc.mel.go2)[gc.mel.go2$cluster_main == "Malignant"])
# DimPlot2(gc.mel.go2, features = "Detailed_cluster_all_5", cells = cells_malig)
DimPlot2(gc.mel.go2, features = "malig_patient_specific", cells = cells_malig, cols = c("lightgrey", color_pro(2)), order = T)

gc.mel.go3 <- subset(gc.mel.go2, cells = cells_malig)
mark1 <- FindMarkers(gc.mel.go3, ident.1 = "Patient_specific_A", group.by = "malig_patient_specific", logfc.threshold = 4, only.pos = T)
rownames(mark1) %>% RenameGO()
DimPlot2(gc.mel.go3, head(rownames(mark1)), theme = NoAxes())
mark2 <- FindMarkers(gc.mel.go3, ident.1 = "Patient_specific_B", group.by = "malig_patient_specific", logfc.threshold = 3, only.pos = T, min.pct = 0.2)
rownames(mark2) %>% RenameGO()
gc.mel.go3$malig2 <- gc.mel.go3$malig_patient_specific
cells_selected <- CellSelector(DimPlot(gc.mel.go3))
gc.mel.go3$malig2[cells_selected] <- "selected"
mark3 <- FindMarkers(gc.mel.go3, ident.1 = "Patient_specific_B", ident.2 = "selected", group.by = "malig2", logfc.threshold = 4, only.pos = T)
mark4 <- FindMarkers(gc.mel.go3, ident.1 = c("Patient_specific_B", "selected"), group.by = "malig2", logfc.threshold = 4, only.pos = T)
DimPlot2(gc.mel.go3, c("malig_patient_specific",head(rownames(mark3), 15)), theme = NoAxes() )
RenameGO(rownames(mark3))
RenameGO(rownames(mark4))

# mark2 <- arrange(mark2, desc(mark2))
# DimPlot2(gc.mel.go3, "malig2", theme = NoAxes() )
DimPlot2(gc.mel.go3, head(rownames(mark1)), theme = NoAxes(), cols = c("lightblue","red"))

# 1k, 3k and 5k variable pathways
# gc.mel.go2 <- FindVariableFeatures(gc.mel.go2, nfeatures = 1000)
# gc.mel.go2 <- FindVariableFeatures(gc.mel.go2, nfeatures = 3000)
gc.mel.go2 <- FindVariableFeatures(gc.mel.go2, nfeatures = 5000)
hvg <- VariableFeatures(gc.mel.go2)
hvg <- FilterGOTerms(hvg, n.min = 10)
hvg <- FilterGOTerms(hvg, n.min = 15)
VariableFeatures(gc.mel.go2) <- hvg[1:1000]
gc.mel.go2 <- ScaleData(gc.mel.go2)
gc.mel.go2 <- RunPCA(gc.mel.go2)
gc.mel.go2 <- RunUMAP(gc.mel.go2, dims = 1:10)
gc.mel.go2 <- RunUMAP(gc.mel.go2, dims = 1:16)
DimPlot2(gc.mel.go2, label = T, theme = NoAxes(), features = c("patient","cluster_main"))
# ggsave("fig5 discussion/5. hvg-5000-dims-16.png", width = 12, height = 5)
ggsave("fig5 discussion/5. hvg-1000-dims-16-nmin-15.png", width = 12, height = 5)

# fig5 redo ---------------------------------------------------------------

gc.mel <- readRDS("~/Documents/Essen/2023-8-8 melanoma integration/GC_BT_harmony_singlet_all_cells_annotated.rds")
# subset features of functional related genes

go.base <- c("GO:0008150","GO:0005575","GO:0003674")
# SearchDatabase("molecular_fun") %>% names
all.genes1 <- Reduce(union, GO_Data$human$GO2Gene[go.base])
all.genes2 <- unique(unlist(Reactome_Data$human$Path2Gene))
all.genes3 <- unique(unlist(Genesets_data$human$GSEA$`KEGG gene sets`))
all.genes4 <- unique(unlist(Genesets_data$human$GSEA$`BioCarta gene sets`))
all.genes <- Reduce(union, list(all.genes1, all.genes2, all.genes3, all.genes4))

all.genes <- intersect(all.genes, rownames(gc.mel))
setdiff(rownames(gc.mel), all.genes)

gc.mel2 <- CreateSeuratObject(
  GetAssayData(gc.mel, slot = "counts")[all.genes,],
  meta.data = gc.mel@meta.data)
genesets <- c(GO_Data$human$GO2Gene, Reactome_Data$human$Path2Gene, Genesets_data$human$GSEA$`KEGG gene sets`, Genesets_data$human$GSEA$`BioCarta gene sets`)
genesets <- genesets[lengths(genesets) >= 15]
gc.mel2 <- GeneSetAnalysis(gc.mel2, genesets = genesets, nCores = 4)

rn <- c(
  "Activated_CD8_Tcells" = "T/NK Cells",
  "CD4_Tcells" = "T/NK Cells",
  "Cycling_CD8_Tcells" = "T/NK Cells",
  "Cytotoxic_CD8_Tcells" = "T/NK Cells",
  "Dysfuntional_CD8_Tcells" = "T/NK Cells",
  "Memory_CD8_Tcells" = "T/NK Cells",
  "Tregs" = "T/NK Cells",
  "NK" = "T/NK Cells",

  "M2 Macrophages" = "Macrophages",
  "Macrophages_CXCL10" = "Macrophages",
  "Macrophages_FOS" = "Macrophages",
  "Macrophages_necrosis" = "Macrophages",
  "Macrophages_PLA2G2D" = "Macrophages",

  "Monocytes_CD14" = "Monocytes",
  "Monocytes_CD16" = "Monocytes",
  "Monocytes_TNFSF13" = "Monocytes",

  "B_cells" = "B Cells",
  "Plasma_cells" = "B Cells",

  "DC1" = "DC",
  "DC2" = "DC",
  "pDC" = "DC",

  "CAF_APCDD1" = "CAF",
  "CAF_CHI3L1" = "CAF",
  "CAF_PLA2G2A" = "CAF",
  "CAF_POSTN" = "CAF",

  "BEC" = "EC/Peri",
  "LEC" = "EC/Peri",
  "Pericyte" = "EC/Peri",

  "Melanocytic_OXPHOS" = "Malignant",
  "Melanoma_immune_like" = "Malignant",
  "Patient_specific_A" = "Malignant",
  "Patient_specific_B" = "Malignant",
  "Doublets" = "Malignant",
  "Mitochondrial" = "Malignant",
  "Mitotic" = "Malignant",
  "Neural_like" = "Malignant",
  "Stress(hypoxia)" = "Malignant",
  "Trans_reg" = "Malignant",
  "Interferon_alpha_beta" = "Malignant",
  "Antigen_presentation" = "Malignant",
  "Mesenchymal" = "Malignant",

  "Melanophages" = "Macrophages",

  "Mast_cells" = "Others"
)
gc.mel.go$cluster_main <- plyr::revalue(gc.mel.go$Detailed_cluster_all_5, rn)
gc.mel.go$cluster_main[is.na(gc.mel.go$cluster_main)] <- "Others"
# sort(unique(gc.mel$Detailed_cluster_all_5))
cells <- colnames(gc.mel)[gc.mel$cluster_main != "Others"]
gc.mel.go <- CreateSeuratObject(
  counts = gc.mel2@misc$AUCell$genesets[,cells], meta.data = gc.mel2@meta.data[cells,]
)
gc.mel.go$cluster_main <- plyr::revalue(gc.mel.go$Detailed_cluster_all_5, rn)
gc.mel.go$cluster_main[is.na(gc.mel.go$cluster_main)] <- "Others"
gc.mel.go$orig.ident <- stringr::str_extract(colnames(gc.mel.go), "(?<=_)[0-9]+")
gc.mel.go$patient <- paste0("P",sprintf("%02d", as.numeric(gc.mel.go$orig.ident)))
gc.mel.go@assays$RNA$data <- gc.mel.go@assays$RNA$counts
gc.mel.go <- FindVariableFeatures(gc.mel.go, nfeatures = 3000)
gc.mel.go <- ScaleData(gc.mel.go)
gc.mel.go <- RunPCA(gc.mel.go)
gc.mel.go <- RunUMAP(gc.mel.go, dims = 1:10)

DimPlot2(gc.mel.go, features = c("patient","cluster_main"), label = T, theme = NoAxes())

gc.mel.go$malig_patient_specific <- gc.mel.go$Detailed_cluster_all_5
gc.mel.go$malig_patient_specific[!gc.mel.go$malig_patient_specific %in% c("Patient_specific_A", "Patient_specific_B")] <- "Others"
cells_malig <- CellSelector(DimPlot(gc.mel.go))
cells_malig <- intersect(cells_malig, colnames(gc.mel.go)[gc.mel.go$cluster_main == "Malignant"])
#
# DimPlot2(gc.mel.go, features = "malig_patient_specific", cells = cells_malig, cols = c("lightgrey", color_pro(2)), order = T)

DimPlot2(gc.mel.go3, "malig_patient_specific", theme = NoAxes(), cols = c("lightgrey", color_pro(2)), order = T, pt.size = 1)
ggsave("fig5 discussion/6. fig5c-1.svg", width = 5, height = 3.5)
gc.mel.go3 <- subset(gc.mel.go, cells = cells_malig)
mark1 <- FindMarkers(gc.mel.go3, ident.1 = "Patient_specific_A", group.by = "malig_patient_specific", logfc.threshold = 4, only.pos = T)
rownames(mark1) %>% RenameGO() %>% RenameReactome()
DimPlot2(gc.mel.go3, head(rownames(mark1),9), theme = NoAxes())
DimPlot2(gc.mel.go3, c("R-HSA-5576886","GO:0014061"), theme = NoAxes() + NoLegend(), pt.size = 1)
ggsave("fig5 discussion/6. fig5c-2.png", width = 5.5, height = 3)
mark2 <- FindMarkers(gc.mel.go3, ident.1 = "Patient_specific_B", group.by = "malig_patient_specific", logfc.threshold = 3, only.pos = T, min.pct = 0.2)
rownames(mark2) %>% RenameGO() %>% RenameReactome()
DimPlot2(gc.mel.go3, head(rownames(mark2)), theme = NoAxes())
DimPlot2(gc.mel.go3, c("GO:0050711","GO:0046068"), theme = NoAxes() + NoLegend(), pt.size = 1)
ggsave("fig5 discussion/6. fig5c-3.png", width = 5.5, height = 3)

DimPlot2(gc.mel.go, label = T, group.by = c("patient","cluster_main"), theme = NoAxes() + NoLegend())
ggsave("fig5 discussion/6. fig5b umap GO no harmony filtered.png", width = 9.5, height = 5)

# gc.mel.go3 <- FindVariableFeatures(gc.mel.go3)
# gc.mel.go3 <- ScaleData(gc.mel.go3)
# gc.mel.go3 <- RunPCA(gc.mel.go3)
# gc.mel.go3 <- RunUMAP(gc.mel.go3, dims = 1:10)
# DimPlot2(gc.mel.go3)

saveRDS(gc.mel.go, "rds/2024-7-8 gc.mel.go.rds")
