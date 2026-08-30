
###############################################################################
# GEO scRNA-seq advanced one-click pipeline
# Author: ITWANGYANG
# Version: 2026-06-18-v6-hdWGCNA-publication-plots
#
# Goal:
#   GEO download -> robust format detection -> Seurat QC/integration -> annotation ->
#   markers/DE/pseudobulk -> enrichment -> multi-method pseudotime -> CellChat ->
#   hdWGCNA -> publication-style vector PDFs -> resumable checkpoints.
#
# Design principles:
#   1. Every heavy module has an RDS checkpoint.
#   2. Every optional module is wrapped by safe_run(). If it fails, the pipeline continues.
#   3. Missing optional packages never stop the whole analysis.
#   4. Results, logs, module status, and parameters are exported.
#
# Usage:
#   Rscript geo_scRNA_advanced_pipeline.R
#
# Main fields to edit are in CONFIG below.
###############################################################################

# =============================== CONFIG ======================================
CFG <- list(
  GSE = "GSE155468",
  species = "human",                    # human, mouse, rat, zebrafish
  work_root = "~/GEO_scRNA_pipeline",
  raw_dir = NULL,                        # NULL means work_dir/rawdata
  input_dir = NULL,                      # optional local directory, overrides GEO download if non-empty
  output_prefix = NULL,
  
  # Automatic grouping. Rules are applied to sample name + GEO metadata text.
  group_rules = list(
    # Regex rules. Sample/file-name hits are weighted much higher than GEO metadata hits.
    # This prevents datasets whose series title contains a disease word from turning all samples into Disease.
    Control = c("^con[0-9._-]*$", "^ctrl", "\\bcontrol\\b", "\\bnormal\\b", "\\bhealthy\\b", "\\bsham\\b", "\\bvehicle\\b", "\\bwildtype\\b", "\\bwt\\b", "\\bbaseline\\b"),
    Disease = c("^taa[0-9._-]*$", "^ataa", "\\bataa\\b", "\\btaa\\b", "\\bdisease\\b", "\\bcase\\b", "\\btumou?r\\b", "\\bcancer\\b", "\\blesion\\b", "\\baneurysm\\b", "\\bko\\b", "\\bmutant\\b"),
    Treatment = c("\\btreated\\b", "\\btreatment\\b", "\\bdrug\\b", "\\bstim", "\\btherapy\\b", "\\bpost\\b", "\\bresponder\\b"),
    NonResponder = c("nonresponder", "non_responder", "non-responder", "\\bnr\\b")
  ),
  min_group_confidence = 1,
  unknown_group_name = "Unknown",
  
  # Basic Seurat processing
  assay = "RNA",
  min_cells = 3,
  min_features = 200,
  qc_nmads = 3,
  max_percent_mt = NULL,                 # NULL means adaptive MAD only; set e.g. 20 if needed
  nfeatures = 3000,
  dims = 1:30,
  resolution = 0.5,
  use_sctransform = FALSE,
  use_harmony = TRUE,
  use_fastmnn = FALSE,
  seed = 123,
  
  # Downsampling guardrails for heavy modules
  max_cells_for_umap_plot = 150000,
  max_cells_for_pseudotime = 12000,
  max_cells_per_celltype_cellchat = 700,
  max_cells_for_cellchat_total = 25000,
  max_celltypes_cellchat = 18,
  max_celltypes_hdwgcna = 4,
  min_cells_per_celltype_hdwgcna = 500,
  min_samples_for_pseudobulk = 2,
  
  # Analysis toggles
  run_annotation = TRUE,
  run_markers = TRUE,
  run_de = TRUE,
  run_pseudobulk = TRUE,
  run_enrichment = TRUE,
  run_pseudotime = TRUE,
  run_cellchat = TRUE,
  run_hdwgcna = TRUE,
  run_pathway_scores = TRUE,
  
  # Pseudotime settings
  trajectory_celltypes = NULL,           # NULL = auto choose major non-immune or all; or c("Fibroblast")
  trajectory_methods = c("monocle3", "slingshot", "monocle2"),
  root_group_priority = c("Control", "Normal", "Healthy"),
  root_celltype_keywords = c("stem", "progenitor", "basal", "naive", "control", "normal"),
  n_pseudotime_genes = 300,
  n_pseudotime_heatmap_genes = 120,
  n_pseudotime_bins = 100,
  monocle2_official_gene_test = TRUE,   # run differentialGeneTest(~sm.ns(Pseudotime)) when possible
  monocle2_run_beam = TRUE,             # run BEAM(branch_point=...) when a branched trajectory exists
  monocle2_branch_point = 1,            # default branch point for BEAM
  monocle2_max_test_genes = 3000,       # cap genes for official Monocle2 tests to avoid huge runs
  monocle2_num_cores = 1,               # safer on macOS/R.app; increase when using terminal/Rscript
  go_show_category = 12,                # fewer categories reduces label overlap
  go_label_wrap = 42,
  
  # CellChat settings
  cellchat_min_cells = 30,
  cellchat_database_category = NULL,     # NULL = all; examples: "Secreted Signaling", "ECM-Receptor"
  cellchat_top_interactions = 20,       # bubble plot only shows top ligand-receptor pairs by probability
  cellchat_top_pathways = 5,            # number of pathway circle plots
  
  # hdWGCNA settings
  hdwgcna_fraction = 0.05,
  hdwgcna_metacell_k = 25,
  hdwgcna_metacell_max_shared = 10,
  hdwgcna_soft_power = NULL,             # NULL = auto
  hdwgcna_target_celltypes = NULL,       # NULL = auto; or c("Fibroblasts", "Smooth muscle cells")
  hdwgcna_priority_celltype_keywords = c("fibro", "smooth", "muscle", "endo", "mono", "macro", "epithel", "t cell", "b cell"),
  hdwgcna_gene_select = "fraction",      # "fraction" is robust for GEO datasets
  hdwgcna_network_type = "signed",
  hdwgcna_soft_power_r2 = 0.8,
  hdwgcna_default_soft_power = 6,
  hdwgcna_min_module_size = 30,
  hdwgcna_merge_cut_height = 0.25,
  hdwgcna_min_cells_per_sample_ct = 25,  # sample coverage guardrail for target cell-type choice
  hdwgcna_min_samples_per_ct = 2,
  hdwgcna_try_full_object_metacells = TRUE,  # official-style: metacells on full object, SetDatExpr on target CT
  hdwgcna_try_subset_fallback = TRUE,        # fallback: subset target CT first, then build metacells
  hdwgcna_force_legacy_assay = TRUE,         # CRITICAL for Seurat v5 + older hdWGCNA: cast Assay5 -> Assay
  hdwgcna_patch_getassay_slot = TRUE,        # extra compatibility shim: slot -> layer for legacy calls
  hdwgcna_rescue_plain_wgcna = TRUE,         # final fallback: plain WGCNA metacell network if hdWGCNA internals fail
  hdwgcna_rescue_plot_plain_dendrogram = FALSE, # plain WGCNA dendrogram is usually crowded; hide by default
  hdwgcna_publication_plots = TRUE,        # add official-style hdWGCNA downstream visualizations
  hdwgcna_max_modules_to_plot = 8,         # cap module network plots to keep PDFs readable
  hdwgcna_hub_genes_per_module = 8,        # hub genes shown in dot/network plots
  hdwgcna_module_feature_ncol = 3,         # grid columns for module feature plots
  
  # Output
  pdf_family = "Times",
  save_png_copy = FALSE,
  force_recompute = FALSE,
  stop_on_critical_failure = FALSE,
  auto_install_missing = FALSE,          # FALSE by default for reproducibility
  checkpoint_version = "v7_hdwgcna_rescue_publication_plots",  # new checkpoint folder; avoids old bad/corrupt CellChat/pseudotime caches
  validate_pdf_output = TRUE,          # remove broken/tiny PDF files and log warnings
  strict_plot_device = TRUE            # always close exact graphics device used by each plot
)
# =============================================================================

# ============================== ENV SETUP ====================================
set.seed(CFG$seed)
options(stringsAsFactors = FALSE)
options(timeout = 7200)
options(future.globals.maxSize = 16 * 1024^3)
Sys.setenv(LANGUAGE = "en")

`%||%` <- function(x, y) if (!is.null(x)) x else y

work_dir <- file.path(path.expand(CFG$work_root), paste0(CFG$GSE, "_advanced_scRNA"))
if (!is.null(CFG$output_prefix)) work_dir <- file.path(path.expand(CFG$work_root), CFG$output_prefix)
raw_dir <- CFG$raw_dir %||% file.path(work_dir, "rawdata")
if (!is.null(CFG$input_dir)) raw_dir <- path.expand(CFG$input_dir)
out_dir <- file.path(work_dir, "results")
fig_dir <- file.path(out_dir, "figures_pdf")
tab_dir <- file.path(out_dir, "tables")
rdata_dir <- file.path(out_dir, paste0("RData_", CFG$checkpoint_version))
log_dir <- file.path(out_dir, "logs")
report_dir <- file.path(out_dir, "reports")
for (d in c(work_dir, raw_dir, out_dir, fig_dir, tab_dir, rdata_dir, log_dir, report_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(log_dir, paste0("pipeline_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
status_file <- file.path(report_dir, "module_status.csv")
param_file <- file.path(report_dir, "parameters.csv")

log_msg <- function(..., level = "INFO") {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] [", level, "] ", paste0(..., collapse = ""))
  message(msg)
  cat(msg, "\n", file = log_file, append = TRUE)
}

module_status <- data.frame(module = character(), status = character(), message = character(), time = character())
add_status <- function(module, status, message = "") {
  row <- data.frame(module = module, status = status, message = as.character(message), time = as.character(Sys.time()))
  assign("module_status", rbind(get("module_status", envir = .GlobalEnv), row), envir = .GlobalEnv)
  try(write.csv(get("module_status", envir = .GlobalEnv), status_file, row.names = FALSE), silent = TRUE)
}

write.csv(data.frame(parameter = names(CFG), value = vapply(CFG, function(x) paste(capture.output(str(x)), collapse = " "), character(1))), param_file, row.names = FALSE)

grDevices::pdf.options(family = CFG$pdf_family)

# ============================= PACKAGE UTILS =================================
cran_pkgs <- c("Seurat", "Matrix", "data.table", "ggplot2", "patchwork", "dplyr", "tidyr", "stringr", "tibble", "ggrepel", "pheatmap", "RColorBrewer", "viridisLite")
bioc_pkgs <- c("SingleCellExperiment", "SummarizedExperiment", "SingleR", "celldex", "GEOquery", "scDblFinder", "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db", "monocle", "monocle3", "slingshot", "tradeSeq", "edgeR", "limma", "zellkonverter", "ComplexHeatmap")
extra_pkgs <- c("harmony", "CellChat", "hdWGCNA", "WGCNA", "SeuratDisk")

pkg_ok <- function(pkg) requireNamespace(pkg, quietly = TRUE)
load_pkg <- function(pkg) {
  ok <- pkg_ok(pkg)
  if (ok) suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  ok
}

install_if_requested <- function(pkg) {
  if (pkg_ok(pkg)) return(TRUE)
  if (!isTRUE(CFG$auto_install_missing)) return(FALSE)
  log_msg("Trying to install missing package: ", pkg)
  ok <- tryCatch({
    if (pkg %in% bioc_pkgs) {
      if (!pkg_ok("BiocManager")) install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
    pkg_ok(pkg)
  }, error = function(e) {
    log_msg("Install failed for ", pkg, ": ", conditionMessage(e), level = "WARN")
    FALSE
  })
  ok
}

for (p in c("Seurat", "Matrix", "data.table", "ggplot2", "patchwork", "dplyr", "tidyr", "stringr", "tibble")) {
  if (!install_if_requested(p) && p %in% c("Seurat", "Matrix", "data.table", "ggplot2", "dplyr")) {
    stop("Required package missing: ", p, ". Please install it first.")
  }
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

HAS <- function(pkg) pkg_ok(pkg)
seurat_v5 <- utils::packageVersion("Seurat") >= "5.0.0"

# ============================== GENERAL UTILS ================================
macaron <- c("#F7B2BD", "#A0E7E5", "#FBE7C6", "#B4F8C8", "#FFAEBC", "#C3AED6", "#FFD6A5", "#A8D8EA", "#FCBAD3", "#BFE9C9", "#F6DFEB", "#D5A6BD", "#AEC6CF", "#FFB7B2", "#C7CEEA", "#E2F0CB", "#FFDAC1", "#B5EAD7", "#D7BDE2", "#AED6F1", "#A9DFBF", "#F9E79F", "#F5CBA7", "#FADBD8")
macaron_pal <- function(n) if (n <= length(macaron)) macaron[seq_len(n)] else grDevices::colorRampPalette(macaron)(n)
macaron_grad <- grDevices::colorRampPalette(c("#5DADE2", "#F7DC6F", "#EC7063"))(100)

species_norm <- tolower(CFG$species)
mt_pattern <- switch(species_norm, mouse = "^mt-", rat = "^Mt-", zebrafish = "^mt-", human = "^MT-", "^MT-")
orgdb_pkg <- switch(species_norm, mouse = "org.Mm.eg.db", rat = NULL, zebrafish = NULL, human = "org.Hs.eg.db", "org.Hs.eg.db")

clean_name <- function(x) {
  x <- basename(x)
  x <- gsub("\\.(txt|csv|tsv|mtx|h5|h5ad|loom|rds|RDS|rdata|RData)(\\.gz)?$", "", x, ignore.case = TRUE)
  x <- gsub("^GSM[0-9]+[_-]?", "", x, ignore.case = TRUE)
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+$", "", x)
  if (!nzchar(x)) x <- paste0("sample_", sample.int(999999, 1))
  x
}

safe_file_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

get_assay_layer <- function(obj, assay = NULL, layer = "data") {
  assay <- assay %||% DefaultAssay(obj)
  if (seurat_v5) GetAssayData(obj, assay = assay, layer = layer) else GetAssayData(obj, assay = assay, slot = layer)
}

join_layers_safe <- function(obj) {
  if (seurat_v5) {
    obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
  }
  obj
}

safe_assay_names <- function(obj) {
  tryCatch(names(obj@assays), error = function(e) character())
}

safe_set_default_assay <- function(obj, assay) {
  if (is.null(obj) || !inherits(obj, "Seurat")) return(obj)
  assay <- assay %||% CFG$assay
  available <- safe_assay_names(obj)
  if (assay %in% available) {
    obj <- tryCatch({ DefaultAssay(obj) <- assay; obj }, error = function(e) {
      log_msg("DefaultAssay set failed for ", assay, ": ", conditionMessage(e), level = "WARN")
      obj
    })
  } else {
    log_msg("DefaultAssay skipped; assay not found: ", assay, level = "WARN")
  }
  obj
}

# hdWGCNA currently has releases/functions that still call Seurat::GetAssayData(slot=...).
# In SeuratObject >=5, slot= is defunct for Assay5.  The safest fix is to give hdWGCNA
# a legacy Assay object with counts/data/scale.data, while keeping the main Seurat object untouched.
force_legacy_assay_for_hdwgcna <- function(obj, assay = CFG$assay) {
  if (is.null(obj) || !inherits(obj, "Seurat")) return(obj)
  if (!isTRUE(CFG$hdwgcna_force_legacy_assay)) return(obj)
  assay <- assay %||% DefaultAssay(obj)
  obj <- join_layers_safe(obj)
  obj <- safe_set_default_assay(obj, assay)
  if (!assay %in% safe_assay_names(obj)) return(obj)
  obj <- tryCatch({
    # Ensure data layer/slot exists before casting; hdWGCNA SetDatExpr usually needs normalized data.
    dat <- tryCatch(get_assay_layer(obj, assay = assay, layer = "data"), error = function(e) NULL)
    if (is.null(dat) || nrow(dat) == 0 || ncol(dat) == 0) obj <- NormalizeData(obj, assay = assay, verbose = FALSE)
    if (inherits(obj[[assay]], "Assay5") || inherits(obj[[assay]], "StdAssay")) {
      obj[[assay]] <- as(obj[[assay]], Class = "Assay")
      log_msg("hdWGCNA compatibility: cast ", assay, " Assay5/StdAssay to legacy Assay for slot= support")
    }
    obj
  }, error = function(e) {
    log_msg("hdWGCNA legacy assay cast failed: ", conditionMessage(e), level = "WARN")
    obj
  })
  obj
}

# Extra safety net: when packages call GetAssayData(slot=) on Seurat v5 objects,
# map slot -> layer. This is deliberately conservative and only used around hdWGCNA.
patch_getassay_slot_compat <- function() {
  if (!isTRUE(CFG$hdwgcna_patch_getassay_slot)) return(invisible(FALSE))
  ns <- asNamespace("SeuratObject")
  if (isTRUE(get0(".chatgpt_getassay_patch_applied", envir = .GlobalEnv, ifnotfound = FALSE))) return(invisible(TRUE))
  patch_one <- function(method_name) {
    if (!exists(method_name, envir = ns, inherits = FALSE)) return(FALSE)
    original <- get(method_name, envir = ns)
    patched <- switch(method_name,
                      "GetAssayData.Seurat" = function(object, assay = NULL, layer = NULL, slot = NULL, ...) {
                        if (is.null(layer) && !is.null(slot)) layer <- slot
                        assay <- assay %||% SeuratObject::DefaultAssay(object)
                        SeuratObject::GetAssayData(object = object[[assay]], layer = layer %||% "data", ...)
                      },
                      "GetAssayData.StdAssay" = function(object, layer = NULL, slot = NULL, ...) {
                        if (is.null(layer) && !is.null(slot)) layer <- slot
                        SeuratObject::LayerData(object = object, layer = layer %||% "data", ...)
                      },
                      NULL
    )
    if (is.null(patched)) return(FALSE)
    tryCatch({
      unlockBinding(method_name, ns); assign(method_name, patched, envir = ns); lockBinding(method_name, ns); TRUE
    }, error = function(e) FALSE)
  }
  ok <- patch_one("GetAssayData.Seurat") | patch_one("GetAssayData.StdAssay")
  assign(".chatgpt_getassay_patch_applied", ok, envir = .GlobalEnv)
  if (ok) log_msg("Applied SeuratObject GetAssayData slot->layer compatibility patch for hdWGCNA")
  invisible(ok)
}

is_valid_pdf <- function(f) {
  if (!file.exists(f)) return(FALSE)
  if (isTRUE(file.info(f)$size < 800)) return(FALSE)
  sig <- tryCatch(readBin(f, what = "raw", n = 4), error = function(e) raw())
  length(sig) == 4 && identical(rawToChar(sig), "%PDF")
}

draw_plot_object <- function(x) {
  if (is.null(x)) return(invisible(NULL))
  if (inherits(x, c("gg", "ggplot", "patchwork"))) {
    print(x)
  } else if (inherits(x, c("Heatmap", "HeatmapList")) && HAS("ComplexHeatmap")) {
    ComplexHeatmap::draw(x, heatmap_legend_side = "right", annotation_legend_side = "right")
  } else if (inherits(x, c("grob", "gTree", "gtable"))) {
    grid::grid.draw(x)
  } else if (inherits(x, "recordedplot")) {
    grDevices::replayPlot(x)
  } else {
    try(print(x), silent = TRUE)
  }
  invisible(x)
}

close_exact_device <- function(dev_id) {
  open_devs <- grDevices::dev.list()
  if (!is.null(open_devs) && dev_id %in% as.integer(open_devs)) {
    try(grDevices::dev.off(which = dev_id), silent = TRUE)
  }
}

save_pdf <- function(plot_obj, file, w = 8, h = 6) {
  f <- file.path(fig_dir, file)
  ok <- FALSE
  grDevices::pdf(f, width = w, height = h, family = CFG$pdf_family, useDingbats = FALSE, onefile = TRUE)
  dev_id <- grDevices::dev.cur()
  on.exit(close_exact_device(dev_id), add = TRUE)
  ok <- tryCatch({ draw_plot_object(plot_obj); TRUE }, error = function(e) {
    log_msg("Figure failed: ", file, " | ", conditionMessage(e), level = "WARN")
    FALSE
  })
  close_exact_device(dev_id)
  if ((!ok || (isTRUE(CFG$validate_pdf_output) && !is_valid_pdf(f))) && file.exists(f)) {
    unlink(f)
    ok <- FALSE
    log_msg("Removed invalid PDF: ", file, level = "WARN")
  }
  if (ok) log_msg("Saved figure: ", file)
  invisible(ok)
}

save_pdf_expr <- function(file, expr, w = 8, h = 6) {
  f <- file.path(fig_dir, file)
  ok <- FALSE
  grDevices::pdf(f, width = w, height = h, family = CFG$pdf_family, useDingbats = FALSE, onefile = TRUE)
  dev_id <- grDevices::dev.cur()
  on.exit(close_exact_device(dev_id), add = TRUE)
  ok <- tryCatch({
    res <- force(expr)
    if (!is.null(res)) draw_plot_object(res)
    TRUE
  }, error = function(e) {
    log_msg("Figure failed: ", file, " | ", conditionMessage(e), level = "WARN")
    FALSE
  })
  close_exact_device(dev_id)
  if ((!ok || (isTRUE(CFG$validate_pdf_output) && !is_valid_pdf(f))) && file.exists(f)) {
    unlink(f)
    ok <- FALSE
    log_msg("Removed invalid PDF: ", file, level = "WARN")
  }
  if (ok) log_msg("Saved figure: ", file)
  invisible(ok)
}

save_table <- function(x, file) {
  f <- file.path(tab_dir, file)
  tryCatch({ write.csv(x, f, row.names = FALSE); log_msg("Saved table: ", file); TRUE }, error = function(e) {
    log_msg("Table save failed: ", file, " | ", conditionMessage(e), level = "WARN")
    FALSE
  })
}

safe_run <- function(module, expr, default = NULL, critical = FALSE) {
  log_msg("Start module: ", module)
  ans <- tryCatch({
    val <- force(expr)
    add_status(module, "OK", "completed")
    log_msg("Finish module: ", module)
    val
  }, error = function(e) {
    msg <- conditionMessage(e)
    add_status(module, if (critical) "FAILED_CRITICAL" else "SKIPPED_OR_FAILED", msg)
    log_msg("Module failed: ", module, " | ", msg, level = if (critical) "ERROR" else "WARN")
    if (critical && isTRUE(CFG$stop_on_critical_failure)) stop(e)
    default
  })
  ans
}

cache_run <- function(name, expr, force = CFG$force_recompute, default = NULL, critical = FALSE) {
  f <- file.path(rdata_dir, paste0(name, ".rds"))
  if (file.exists(f) && !isTRUE(force)) {
    cached <- tryCatch(readRDS(f), error = function(e) NULL)
    if (!is.null(cached)) {
      if (grepl("combined|pathway_scores", name) && !inherits(cached, "Seurat")) {
        log_msg("Checkpoint rejected because it is not a Seurat object: ", name, level = "WARN")
      } else {
        log_msg("Checkpoint load: ", name)
        return(cached)
      }
    } else {
      log_msg("Checkpoint read failed or NULL, recomputing: ", name, level = "WARN")
    }
  }
  log_msg("Start module: ", name)
  ok <- TRUE
  out <- tryCatch({
    val <- force(expr)
    add_status(name, "OK", "completed")
    log_msg("Finish module: ", name)
    val
  }, error = function(e) {
    ok <<- FALSE
    msg <- conditionMessage(e)
    add_status(name, if (critical) "FAILED_CRITICAL" else "SKIPPED_OR_FAILED", msg)
    log_msg("Module failed: ", name, " | ", msg, level = if (critical) "ERROR" else "WARN")
    if (critical && isTRUE(CFG$stop_on_critical_failure)) stop(e)
    default
  })
  if (ok && !is.null(out)) {
    tryCatch({ saveRDS(out, f); log_msg("Checkpoint saved: ", name) }, error = function(e) log_msg("Checkpoint save failed: ", name, " | ", conditionMessage(e), level = "WARN"))
  }
  out
}

is_outlier_mad <- function(x, nmads = 3, type = c("both", "lower", "higher"), log = FALSE, batch = NULL) {
  type <- match.arg(type)
  xx <- if (log) log10(x + 1) else x
  flag <- rep(FALSE, length(x))
  grp <- if (is.null(batch)) rep("all", length(x)) else as.character(batch)
  for (g in unique(grp)) {
    idx <- grp == g
    med <- stats::median(xx[idx], na.rm = TRUE)
    md <- stats::mad(xx[idx], na.rm = TRUE)
    if (!is.finite(md) || md == 0) md <- 1e-8
    if (type %in% c("both", "lower")) flag[idx] <- flag[idx] | (xx[idx] < med - nmads * md)
    if (type %in% c("both", "higher")) flag[idx] <- flag[idx] | (xx[idx] > med + nmads * md)
  }
  flag[is.na(flag)] <- TRUE
  flag
}

# ============================== GEO METADATA =================================
fetch_geo_metadata <- function(gse) {
  if (!HAS("GEOquery")) return(data.frame())
  safe_run("fetch_geo_metadata", {
    suppressPackageStartupMessages(library(GEOquery))
    gset <- GEOquery::getGEO(gse, GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = FALSE)
    if (length(gset) > 1) gset <- gset[[1]] else gset <- gset[[1]]
    pd <- Biobase::pData(gset)
    pd$geo_accession <- rownames(pd)
    pd
  }, default = data.frame())
}

assign_group_from_text <- function(sample, meta_text = "") {
  # Scientific auto-grouping strategy:
  # 1) file/sample name is the most trustworthy for GEO supplemental files;
  # 2) GEO metadata is useful but can contain series-level disease words affecting all samples;
  # 3) regex patterns are supported and sample-name hits are weighted higher.
  sample_text <- tolower(paste(sample, collapse = " "))
  meta_text <- tolower(paste(meta_text, collapse = " "))
  score_one <- function(txt, keys) {
    sum(vapply(keys, function(k) isTRUE(grepl(k, txt, ignore.case = TRUE, perl = TRUE)), logical(1)))
  }
  sample_scores <- vapply(CFG$group_rules, function(keys) score_one(sample_text, keys), numeric(1))
  meta_scores <- vapply(CFG$group_rules, function(keys) score_one(meta_text, keys), numeric(1))
  scores <- sample_scores * 10 + meta_scores
  if (length(scores) == 0 || max(scores, na.rm = TRUE) < CFG$min_group_confidence) {
    return(list(group = CFG$unknown_group_name, confidence = 0, matched = ""))
  }
  top <- names(scores)[which.max(scores)]
  keys <- CFG$group_rules[[top]]
  matched_sample <- keys[vapply(keys, function(k) grepl(k, sample_text, ignore.case = TRUE, perl = TRUE), logical(1))]
  matched_meta <- keys[vapply(keys, function(k) grepl(k, meta_text, ignore.case = TRUE, perl = TRUE), logical(1))]
  matched <- paste(unique(c(paste0("sample:", matched_sample), paste0("meta:", matched_meta))), collapse = ";")
  list(group = top, confidence = unname(scores[[top]]), matched = matched)
}

# ============================== DATA READING =================================
decompress_archives <- function(dir_path) {
  files <- list.files(dir_path, recursive = TRUE, full.names = TRUE)
  for (f in files) {
    bn <- basename(f)
    if (grepl("\\.tar$|\\.tar\\.gz$|\\.tgz$", bn, ignore.case = TRUE)) {
      safe_run(paste0("untar_", bn), untar(f, exdir = dirname(f)), default = NULL)
    } else if (grepl("\\.zip$", bn, ignore.case = TRUE)) {
      safe_run(paste0("unzip_", bn), unzip(f, exdir = dirname(f)), default = NULL)
    } else if (grepl("\\.gz$", bn, ignore.case = TRUE) && !grepl("\\.(mtx|tsv|txt|csv)\\.gz$", bn, ignore.case = TRUE)) {
      if (HAS("R.utils")) {
        safe_run(paste0("gunzip_", bn), R.utils::gunzip(f, remove = FALSE, overwrite = FALSE), default = NULL)
      }
    }
  }
  invisible(TRUE)
}

download_geo_supp <- function(gse, dir_path) {
  if (length(list.files(dir_path, recursive = TRUE, full.names = TRUE)) > 0) {
    log_msg("Raw directory is not empty; skip GEO download.")
    return(TRUE)
  }
  ok <- FALSE
  if (HAS("GEOquery")) {
    ok <- safe_run("GEOquery_getGEOSuppFiles", {
      GEOquery::getGEOSuppFiles(gse, baseDir = dir_path, makeDirectory = TRUE)
      TRUE
    }, default = FALSE)
  }
  if (!isTRUE(ok)) {
    dest <- file.path(dir_path, paste0(gse, "_RAW.tar"))
    url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=", gse, "&format=file")
    ok <- safe_run("download_RAW_tar", {
      utils::download.file(url, dest, mode = "wb", method = "libcurl")
      file.exists(dest) && file.info(dest)$size > 1024
    }, default = FALSE)
  }
  decompress_archives(dir_path)
  ok
}

looks_like_10x_dir <- function(d) {
  fs <- basename(list.files(d, full.names = FALSE))
  any(grepl("matrix\\.mtx(\\.gz)?$", fs, ignore.case = TRUE)) &&
    any(grepl("(features|genes)\\.tsv(\\.gz)?$", fs, ignore.case = TRUE)) &&
    any(grepl("barcodes\\.tsv(\\.gz)?$", fs, ignore.case = TRUE))
}

find_10x_dirs <- function(root) {
  dirs <- unique(dirname(list.files(root, pattern = "matrix\\.mtx(\\.gz)?$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)))
  dirs[vapply(dirs, looks_like_10x_dir, logical(1))]
}

read_dense_matrix <- function(f) {
  log_msg("Reading dense matrix: ", f)
  sep <- ifelse(grepl("\\.csv(\\.gz)?$", f, ignore.case = TRUE), ",", "\t")
  dt <- data.table::fread(f, sep = sep, data.table = FALSE, check.names = FALSE, showProgress = FALSE)
  if (ncol(dt) < 2) dt <- data.table::fread(f, data.table = FALSE, check.names = FALSE, showProgress = FALSE)
  gene_col <- 1
  genes <- as.character(dt[[gene_col]])
  mat <- as.matrix(dt[, -gene_col, drop = FALSE])
  mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  rownames(mat) <- make.unique(genes)
  colnames(mat) <- make.unique(colnames(dt)[-gene_col])
  Matrix::Matrix(mat, sparse = TRUE)
}

read_r_object_as_seurat <- function(f) {
  obj <- NULL
  if (grepl("\\.rds$", f, ignore.case = TRUE)) {
    obj <- readRDS(f)
  } else {
    e <- new.env(parent = emptyenv())
    load(f, envir = e)
    nms <- ls(e)
    if (length(nms) == 0) stop("RData file contains no object")
    # Prefer Seurat object, then SingleCellExperiment, then matrix-like object.
    for (nm in nms) if (inherits(get(nm, envir = e), "Seurat")) obj <- get(nm, envir = e)
    if (is.null(obj)) for (nm in nms) if (inherits(get(nm, envir = e), "SingleCellExperiment")) obj <- get(nm, envir = e)
    if (is.null(obj)) obj <- get(nms[1], envir = e)
  }
  if (inherits(obj, "Seurat")) return(obj)
  if (inherits(obj, "SingleCellExperiment")) return(Seurat::as.Seurat(obj, counts = "counts", data = NULL))
  if (inherits(obj, "matrix") || inherits(obj, "dgCMatrix") || inherits(obj, "data.frame")) {
    mat <- if (inherits(obj, "data.frame")) as.matrix(obj) else obj
    return(CreateSeuratObject(Matrix::Matrix(mat, sparse = TRUE), min.cells = CFG$min_cells, min.features = CFG$min_features))
  }
  stop("Unsupported R object class: ", paste(class(obj), collapse = ","))
}

read_h5ad_as_seurat <- function(f) {
  if (!HAS("zellkonverter")) stop("Package zellkonverter is required for h5ad files")
  sce <- zellkonverter::readH5AD(f)
  Seurat::as.Seurat(sce, counts = "X", data = NULL)
}

read_loom_as_seurat <- function(f) {
  if (!HAS("SeuratDisk")) stop("Package SeuratDisk is required for loom files")
  h5s <- sub("\\.loom$", ".h5seurat", f, ignore.case = TRUE)
  SeuratDisk::Convert(f, dest = "h5seurat", overwrite = TRUE)
  SeuratDisk::LoadH5Seurat(h5s)
}

discover_input_specs <- function(dir_path) {
  decompress_archives(dir_path)
  all_files <- list.files(dir_path, recursive = TRUE, full.names = TRUE)
  all_files <- all_files[!grepl("\\.(tar|tgz|zip)$", all_files, ignore.case = TRUE)]
  specs <- list()
  
  # 10x directories
  for (d in find_10x_dirs(dir_path)) {
    s <- clean_name(d)
    specs[[paste0("10xdir__", s, "__", length(specs) + 1)]] <- list(type = "10xdir", sample = s, path = d)
  }
  
  # 10x h5 / h5ad / loom / R objects
  for (f in all_files[grepl("\\.h5$", all_files, ignore.case = TRUE)]) {
    if (grepl("filtered_feature_bc_matrix|raw_feature_bc_matrix|matrix|10x|counts|feature", basename(f), ignore.case = TRUE)) {
      s <- clean_name(f)
      specs[[paste0("h5__", s, "__", length(specs) + 1)]] <- list(type = "h5", sample = s, path = f)
    }
  }
  for (f in all_files[grepl("\\.h5ad$", all_files, ignore.case = TRUE)]) {
    s <- clean_name(f)
    specs[[paste0("h5ad__", s, "__", length(specs) + 1)]] <- list(type = "h5ad", sample = s, path = f)
  }
  for (f in all_files[grepl("\\.loom$", all_files, ignore.case = TRUE)]) {
    s <- clean_name(f)
    specs[[paste0("loom__", s, "__", length(specs) + 1)]] <- list(type = "loom", sample = s, path = f)
  }
  for (f in all_files[grepl("\\.(rds|RDS|rdata|RData)$", all_files, ignore.case = TRUE)]) {
    s <- clean_name(f)
    specs[[paste0("robj__", s, "__", length(specs) + 1)]] <- list(type = "robj", sample = s, path = f)
  }
  
  # Dense matrices; exclude 10x sidecar files
  dense <- all_files[grepl("\\.(txt|csv|tsv)(\\.gz)?$", all_files, ignore.case = TRUE) &
                       !grepl("barcodes|features|genes|matrix\\.mtx|annotation|metadata|meta|sample", basename(all_files), ignore.case = TRUE)]
  for (f in dense) {
    s <- clean_name(f)
    specs[[paste0("dense__", s, "__", length(specs) + 1)]] <- list(type = "dense", sample = s, path = f)
  }
  
  # Remove child dense files already covered by 10x dirs by keeping all, but reader will be safe.
  specs
}

read_spec_to_seurat <- function(sp) {
  obj <- switch(sp$type,
                "10xdir" = {
                  x <- Read10X(data.dir = sp$path)
                  if (is.list(x)) x <- x[[which.max(vapply(x, nrow, integer(1)))]]
                  CreateSeuratObject(x, project = sp$sample, min.cells = CFG$min_cells, min.features = CFG$min_features)
                },
                "h5" = {
                  if (!HAS("hdf5r")) log_msg("hdf5r not installed; Read10X_h5 may still work if dependencies are available", level = "WARN")
                  x <- Read10X_h5(sp$path)
                  if (is.list(x)) x <- x[[which.max(vapply(x, nrow, integer(1)))]]
                  CreateSeuratObject(x, project = sp$sample, min.cells = CFG$min_cells, min.features = CFG$min_features)
                },
                "h5ad" = read_h5ad_as_seurat(sp$path),
                "loom" = read_loom_as_seurat(sp$path),
                "robj" = read_r_object_as_seurat(sp$path),
                "dense" = CreateSeuratObject(read_dense_matrix(sp$path), project = sp$sample, min.cells = CFG$min_cells, min.features = CFG$min_features),
                stop("Unsupported spec type: ", sp$type)
  )
  obj$sample <- if (!"sample" %in% colnames(obj@meta.data)) sp$sample else as.character(obj$sample)
  obj$source_file <- sp$path
  obj$source_type <- sp$type
  obj
}

# ============================== MAIN BUILD ===================================
build_combined <- function() {
  if (is.null(CFG$input_dir)) download_geo_supp(CFG$GSE, raw_dir)
  specs <- discover_input_specs(raw_dir)
  if (length(specs) == 0) stop("No supported scRNA input files were detected in: ", raw_dir)
  spec_table <- do.call(rbind, lapply(specs, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
  save_table(spec_table, "00_detected_input_specs.csv")
  
  geo_meta <- fetch_geo_metadata(CFG$GSE)
  if (nrow(geo_meta) > 0) save_table(geo_meta, "00_geo_metadata_raw.csv")
  
  seu_list <- list()
  read_status <- list()
  for (nm in names(specs)) {
    sp <- specs[[nm]]
    obj <- safe_run(paste0("read_", nm), read_spec_to_seurat(sp), default = NULL)
    if (is.null(obj) || !inherits(obj, "Seurat") || ncol(obj) == 0) {
      read_status[[nm]] <- data.frame(spec = nm, sample = sp$sample, type = sp$type, status = "failed")
      next
    }
    obj <- join_layers_safe(obj)
    sample <- unique(as.character(obj$sample))[1]
    meta_text <- ""
    if (nrow(geo_meta) > 0) {
      gidx <- grepl(sample, apply(geo_meta, 1, paste, collapse = " "), ignore.case = TRUE)
      if (any(gidx)) meta_text <- paste(apply(geo_meta[gidx, , drop = FALSE], 1, paste, collapse = " "), collapse = " ")
    }
    gp <- assign_group_from_text(sample, meta_text)
    obj$group <- gp$group
    obj$group_confidence <- gp$confidence
    obj$group_matched_keywords <- gp$matched
    obj$orig_sample <- sample
    colnames(obj) <- paste0(make.unique(rep(sample, ncol(obj))), "__", colnames(obj))
    seu_list[[sample]] <- obj
    read_status[[nm]] <- data.frame(spec = nm, sample = sample, type = sp$type, status = "ok", n_genes = nrow(obj), n_cells = ncol(obj), group = gp$group, confidence = gp$confidence, matched = gp$matched)
    log_msg("Read sample: ", sample, " | cells=", ncol(obj), " | genes=", nrow(obj), " | group=", gp$group)
  }
  if (length(read_status) > 0) save_table(rbind.fill2(read_status), "00_read_status.csv")
  if (length(seu_list) == 0) stop("All input readers failed; check raw files and logs.")
  combined <- if (length(seu_list) == 1) seu_list[[1]] else merge(seu_list[[1]], y = seu_list[-1], project = CFG$GSE)
  combined <- join_layers_safe(combined)
  combined
}

rbind.fill2 <- function(x) {
  if (length(x) == 0) return(data.frame())
  cols <- unique(unlist(lapply(x, names)))
  do.call(rbind, lapply(x, function(d) {
    miss <- setdiff(cols, names(d))
    for (m in miss) d[[m]] <- NA
    d[, cols, drop = FALSE]
  }))
}

combined <- cache_run("01_combined_raw", build_combined(), critical = TRUE)
combined <- safe_set_default_assay(combined, CFG$assay)
combined <- join_layers_safe(combined)

# =============================== QC MODULE ===================================
run_qc <- function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = ifelse(species_norm == "human", "^RP[SL]", "^Rp[sl]"))
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = ifelse(species_norm == "human", "^HB[ABDEGQZ]", "^Hb[ab]"))
  
  obj$discard_low_features <- is_outlier_mad(obj$nFeature_RNA, CFG$qc_nmads, "lower", log = TRUE, batch = obj$sample) | obj$nFeature_RNA <= CFG$min_features
  obj$discard_low_counts <- is_outlier_mad(obj$nCount_RNA, CFG$qc_nmads, "lower", log = TRUE, batch = obj$sample)
  obj$discard_high_mt <- is_outlier_mad(obj$percent.mt, CFG$qc_nmads, "higher", batch = obj$sample)
  if (!is.null(CFG$max_percent_mt)) obj$discard_high_mt <- obj$discard_high_mt | obj$percent.mt > CFG$max_percent_mt
  obj$doublet_call <- "not_tested"
  
  save_pdf(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.hb"), group.by = "sample", pt.size = 0, ncol = 5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)), "01_QC_violin_before_filter.pdf", 18, 5)
  save_pdf(FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "sample"), "02_QC_scatter_before_filter.pdf", 8, 6)
  
  if (HAS("scDblFinder") && HAS("SingleCellExperiment") && HAS("SummarizedExperiment")) {
    obj <- safe_run("doublet_detection_scDblFinder", {
      sce <- SingleCellExperiment::SingleCellExperiment(list(counts = get_assay_layer(obj, layer = "counts")))
      sce$sample <- obj$sample
      sce <- scDblFinder::scDblFinder(sce, samples = sce$sample)
      cls <- as.character(SummarizedExperiment::colData(sce)$scDblFinder.class)
      obj$doublet_call <- cls
      obj
    }, default = obj)
  }
  
  obj$discard_doublet <- obj$doublet_call %in% c("doublet", "Doublet")
  obj$discard <- obj$discard_low_features | obj$discard_low_counts | obj$discard_high_mt | obj$discard_doublet
  qc_summary <- obj@meta.data %>%
    group_by(sample, group) %>%
    summarise(cells_before = n(), discarded = sum(discard), kept = sum(!discard), median_features = median(nFeature_RNA), median_counts = median(nCount_RNA), median_mt = median(percent.mt), .groups = "drop")
  save_table(qc_summary, "01_QC_summary_by_sample.csv")
  obj <- subset(obj, subset = discard == FALSE)
  save_pdf(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "sample", pt.size = 0, ncol = 3) + theme(axis.text.x = element_text(angle = 45, hjust = 1)), "03_QC_violin_after_filter.pdf", 14, 5)
  obj
}

combined <- cache_run("02_combined_qc", run_qc(combined), default = combined, critical = FALSE)
combined <- safe_set_default_assay(combined, CFG$assay)

# =========================== NORMALIZE / INTEGRATE ===========================
run_integration <- function(obj) {
  obj <- join_layers_safe(obj)
  obj <- safe_set_default_assay(obj, CFG$assay)
  if (isTRUE(CFG$use_sctransform) && HAS("sctransform")) {
    obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = FALSE)
    obj <- safe_set_default_assay(obj, "SCT")
  } else {
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = CFG$nfeatures, verbose = FALSE)
    obj <- ScaleData(obj, vars.to.regress = "percent.mt", verbose = FALSE)
  }
  obj <- RunPCA(obj, npcs = max(CFG$dims), verbose = FALSE)
  reduction_use <- "pca"
  if (isTRUE(CFG$use_harmony) && length(unique(obj$sample)) > 1 && HAS("harmony")) {
    obj <- safe_run("RunHarmony", {
      harmony::RunHarmony(obj, group.by.vars = "sample", reduction.use = "pca", reduction.save = "harmony", verbose = FALSE)
    }, default = obj)
    if ("harmony" %in% names(obj@reductions)) reduction_use <- "harmony"
  }
  obj <- RunUMAP(obj, reduction = reduction_use, dims = CFG$dims, verbose = FALSE)
  obj <- RunTSNE(obj, reduction = reduction_use, dims = CFG$dims, check_duplicates = FALSE)
  obj <- FindNeighbors(obj, reduction = reduction_use, dims = CFG$dims, verbose = FALSE)
  obj <- FindClusters(obj, resolution = CFG$resolution, verbose = FALSE)
  obj$cluster <- as.character(obj$seurat_clusters)
  obj@misc$pipeline$main_reduction <- reduction_use
  obj
}

combined <- cache_run("03_combined_integrated", run_integration(combined), default = combined, critical = FALSE)
combined <- safe_set_default_assay(combined, if (isTRUE(CFG$use_sctransform) && "SCT" %in% safe_assay_names(combined)) "SCT" else CFG$assay)

# ============================= BASIC FIGURES =================================
plot_umaps <- function(obj) {
  n_sample <- length(unique(obj$sample)); n_group <- length(unique(obj$group)); n_cluster <- length(unique(obj$seurat_clusters))
  save_pdf(DimPlot(obj, group.by = "sample", cols = macaron_pal(n_sample)) + ggtitle("UMAP by sample"), "04_UMAP_sample.pdf", 8, 6)
  save_pdf(DimPlot(obj, group.by = "group", cols = macaron_pal(n_group)) + ggtitle("UMAP by inferred group"), "05_UMAP_group.pdf", 8, 6)
  save_pdf(DimPlot(obj, group.by = "seurat_clusters", label = TRUE, repel = TRUE, cols = macaron_pal(n_cluster)) + ggtitle("UMAP by cluster") + NoLegend(), "06_UMAP_cluster.pdf", 8, 6)
  save_pdf((DimPlot(obj, group.by = "sample", cols = macaron_pal(n_sample)) | DimPlot(obj, group.by = "group", cols = macaron_pal(n_group)) | DimPlot(obj, label = TRUE, cols = macaron_pal(n_cluster)) + NoLegend()), "07_UMAP_overview.pdf", 18, 6)
}
safe_run("basic_umap_figures", plot_umaps(combined), default = NULL)

# =============================== ANNOTATION ==================================
marker_db <- list(
  human = list(
    T_cells = c("CD3D", "CD3E", "TRAC", "IL7R"),
    CD8_T_cells = c("CD8A", "CD8B", "NKG7", "GZMB"),
    NK_cells = c("NKG7", "GNLY", "KLRD1", "PRF1"),
    B_cells = c("MS4A1", "CD79A", "CD79B", "BANK1"),
    Plasma_cells = c("MZB1", "JCHAIN", "SDC1", "XBP1"),
    Monocytes_Macrophages = c("LYZ", "CD68", "AIF1", "LST1"),
    Dendritic_cells = c("FCER1A", "CLEC10A", "ITGAX", "LILRA4"),
    Endothelial_cells = c("PECAM1", "VWF", "KDR", "CLDN5"),
    Fibroblasts = c("DCN", "LUM", "COL1A1", "COL1A2"),
    Smooth_muscle_cells = c("ACTA2", "MYH11", "TAGLN", "CNN1"),
    Epithelial_cells = c("EPCAM", "KRT8", "KRT18", "KRT19"),
    Mast_cells = c("TPSAB1", "TPSB2", "CPA3", "KIT"),
    Neutrophils = c("S100A8", "S100A9", "FCGR3B", "CSF3R"),
    Pericytes = c("RGS5", "PDGFRB", "CSPG4", "MCAM")
  ),
  mouse = list(
    T_cells = c("Cd3d", "Cd3e", "Trac", "Il7r"),
    CD8_T_cells = c("Cd8a", "Cd8b1", "Nkg7", "Gzmb"),
    NK_cells = c("Nkg7", "Klrd1", "Prf1", "Gzma"),
    B_cells = c("Ms4a1", "Cd79a", "Cd79b", "Bank1"),
    Plasma_cells = c("Mzb1", "Jchain", "Sdc1", "Xbp1"),
    Monocytes_Macrophages = c("Lyz2", "Adgre1", "Aif1", "Lpl"),
    Dendritic_cells = c("Fcer1a", "Itgax", "Clec10a", "Siglech"),
    Endothelial_cells = c("Pecam1", "Vwf", "Kdr", "Cldn5"),
    Fibroblasts = c("Dcn", "Lum", "Col1a1", "Col1a2"),
    Smooth_muscle_cells = c("Acta2", "Myh11", "Tagln", "Cnn1"),
    Epithelial_cells = c("Epcam", "Krt8", "Krt18", "Krt19"),
    Mast_cells = c("Tpsab1", "Tpsb2", "Cpa3", "Kit"),
    Neutrophils = c("S100a8", "S100a9", "Csf3r", "Mpo"),
    Pericytes = c("Rgs5", "Pdgfrb", "Cspg4", "Mcam")
  )
)

score_marker_annotation <- function(obj) {
  db <- marker_db[[ifelse(species_norm == "mouse", "mouse", "human")]]
  db <- lapply(db, function(gs) intersect(gs, rownames(obj)))
  db <- db[vapply(db, length, integer(1)) >= 2]
  if (length(db) == 0) {
    obj$celltype_marker <- "Unknown"
    return(obj)
  }
  obj <- AddModuleScore(obj, features = db, name = "marker_score_", assay = DefaultAssay(obj), search = FALSE)
  score_cols <- grep("^marker_score_", colnames(obj@meta.data), value = TRUE)
  if (length(score_cols) == 0) { obj$celltype_marker <- "Unknown"; return(obj) }
  score_mat <- obj@meta.data[, score_cols, drop = FALSE]
  colnames(score_mat) <- names(db)
  obj$celltype_marker <- colnames(score_mat)[max.col(score_mat, ties.method = "first")]
  obj
}

run_annotation <- function(obj) {
  obj$celltype_singleR <- NA_character_
  if (HAS("SingleR") && HAS("celldex") && HAS("SummarizedExperiment")) {
    obj <- safe_run("SingleR_annotation", {
      suppressPackageStartupMessages({ library(SingleR); library(celldex) })
      ref_list <- list()
      lab_list <- list()
      if (species_norm == "mouse") {
        ref <- celldex::MouseRNAseqData()
        ref_list <- list(MouseRNAseq = ref); lab_list <- list(ref$label.main)
      } else {
        hpca <- celldex::HumanPrimaryCellAtlasData()
        bp <- celldex::BlueprintEncodeData()
        ref_list <- list(HPCA = hpca, BlueprintEncode = bp); lab_list <- list(hpca$label.main, bp$label.main)
      }
      dat <- get_assay_layer(obj, assay = DefaultAssay(obj), layer = "data")
      pred <- SingleR(test = dat, ref = ref_list, labels = lab_list, clusters = obj$seurat_clusters)
      c2t <- setNames(pred$labels, rownames(pred))
      obj$celltype_singleR <- unname(c2t[as.character(obj$seurat_clusters)])
      ann <- data.frame(cluster = rownames(pred), SingleR_label = pred$labels, stringsAsFactors = FALSE)
      save_table(ann, "02_cluster_annotation_singleR.csv")
      obj
    }, default = obj)
  }
  obj <- score_marker_annotation(obj)
  obj$celltype <- ifelse(!is.na(obj$celltype_singleR) & nzchar(obj$celltype_singleR), obj$celltype_singleR, obj$celltype_marker)
  obj$celltype[is.na(obj$celltype) | !nzchar(obj$celltype)] <- "Unknown"
  obj$celltype <- as.character(obj$celltype)
  # Preserve one final label within each cluster by cluster majority.
  # Base R avoids dplyr::slice / Bioconductor S4 method conflicts.
  tmp <- as.data.frame(table(seurat_clusters = obj$seurat_clusters, celltype = obj$celltype), stringsAsFactors = FALSE)
  tmp <- tmp[tmp$Freq > 0, , drop = FALSE]
  tmp <- tmp[order(tmp$seurat_clusters, -tmp$Freq), , drop = FALSE]
  tmp <- tmp[!duplicated(tmp$seurat_clusters), , drop = FALSE]
  c2t <- setNames(as.character(tmp$celltype), as.character(tmp$seurat_clusters))
  obj$celltype <- unname(c2t[as.character(obj$seurat_clusters)])
  save_table(data.frame(cell = colnames(obj), obj@meta.data[, c("sample", "group", "seurat_clusters", "celltype", "celltype_singleR", "celltype_marker"), drop = FALSE]), "02_cell_metadata_annotation.csv")
  n_ct <- length(unique(obj$celltype))
  save_pdf(DimPlot(obj, group.by = "celltype", label = TRUE, repel = TRUE, cols = macaron_pal(n_ct)) + ggtitle("Cell types"), "08_UMAP_celltype.pdf", 9, 7)
  save_pdf(DimPlot(obj, group.by = "celltype", split.by = "group", label = TRUE, repel = TRUE, cols = macaron_pal(n_ct)) + ggtitle("Cell types by group"), "09_UMAP_celltype_by_group.pdf", 16, 7)
  db <- marker_db[[ifelse(species_norm == "mouse", "mouse", "human")]]
  canonical <- unique(intersect(unlist(db), rownames(obj)))
  canonical <- canonical[seq_len(min(length(canonical), 60))]
  if (length(canonical) >= 3) {
    save_pdf(DotPlot(obj, features = canonical, group.by = "celltype") + RotatedAxis() + ggtitle("Canonical marker dotplot"), "10_marker_dotplot.pdf", 14, 8)
    save_pdf(FeaturePlot(obj, features = canonical[seq_len(min(12, length(canonical)))], ncol = 4, order = TRUE), "11_marker_featureplot.pdf", 14, 10)
  }
  obj
}

if (isTRUE(CFG$run_annotation)) combined <- cache_run("04_combined_annotated", run_annotation(combined), default = combined, critical = FALSE) else combined$celltype <- as.character(combined$seurat_clusters)
if (!"celltype" %in% colnames(combined@meta.data)) combined$celltype <- as.character(combined$seurat_clusters)
combined <- safe_set_default_assay(combined, if (isTRUE(CFG$use_sctransform) && "SCT" %in% safe_assay_names(combined)) "SCT" else CFG$assay)

# ============================== MARKERS / DE =================================
run_markers <- function(obj) {
  if (is.null(obj) || !inherits(obj, "Seurat")) stop("Input object is not a Seurat object")
  Idents(obj) <- "celltype"
  mk <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  save_table(mk, "03_celltype_markers.csv")
  if (nrow(mk) > 0) {
    top <- mk %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>% pull(gene) %>% unique()
    top <- intersect(top, rownames(obj))
    if (length(top) >= 3) save_pdf(DoHeatmap(subset(obj, downsample = 60), features = top, group.by = "celltype", raster = FALSE) + ggtitle("Top markers per cell type"), "12_marker_heatmap.pdf", 12, 14)
  }
  mk
}
all_markers <- if (isTRUE(CFG$run_markers)) cache_run("05_all_markers", run_markers(combined), default = data.frame()) else data.frame()

volcano_plot <- function(df, title, fc = 0.25, padj = 0.05) {
  df <- df %>% mutate(sig = case_when(p_val_adj < padj & avg_log2FC > fc ~ "Up", p_val_adj < padj & avg_log2FC < -fc ~ "Down", TRUE ~ "NS"), neglogp = -log10(p_val_adj + 1e-300))
  p <- ggplot(df, aes(avg_log2FC, neglogp, color = sig)) + geom_point(size = 0.8, alpha = 0.7) +
    scale_color_manual(values = c(Up = "#EC7063", Down = "#5DADE2", NS = "grey82")) +
    geom_vline(xintercept = c(-fc, fc), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(padj), linetype = "dashed", color = "grey50") +
    labs(title = title, x = "Average log2 fold change", y = "-log10 adjusted p value") + theme_bw(base_size = 12) + theme(legend.title = element_blank())
  if (HAS("ggrepel")) {
    lab <- df %>% filter(sig != "NS") %>% arrange(desc(abs(avg_log2FC))) %>% head(15)
    p <- p + ggrepel::geom_text_repel(data = lab, aes(label = gene), color = "black", size = 3, max.overlaps = 50)
  }
  p
}

run_de <- function(obj) {
  if (is.null(obj) || !inherits(obj, "Seurat")) stop("Input object is not a Seurat object")
  groups <- setdiff(unique(as.character(obj$group)), CFG$unknown_group_name)
  if (length(groups) < 2) {
    log_msg("DE skipped because fewer than 2 non-Unknown groups were detected", level = "WARN")
    return(list(overall = data.frame(), by_celltype = data.frame()))
  }
  g1 <- groups[1]; g2 <- groups[2]
  Idents(obj) <- "group"
  overall <- safe_run("DE_overall_FindMarkers", {
    r <- FindMarkers(obj, ident.1 = g2, ident.2 = g1, logfc.threshold = 0, min.pct = 0.1)
    r$gene <- rownames(r)
    r$contrast <- paste(g2, "vs", g1)
    r
  }, default = data.frame())
  if (nrow(overall) > 0) {
    save_table(overall, "04_DE_overall_FindMarkers.csv")
    save_pdf(volcano_plot(overall, paste0("Overall DE: ", g2, " vs ", g1)), "13_volcano_overall.pdf", 8, 7)
  }
  obj$ct_grp <- paste(obj$celltype, obj$group, sep = "__")
  Idents(obj) <- "ct_grp"
  L <- list()
  for (ct in unique(obj$celltype)) {
    i1 <- paste0(ct, "__", g1); i2 <- paste0(ct, "__", g2)
    if (sum(Idents(obj) == i1) >= 20 && sum(Idents(obj) == i2) >= 20) {
      r <- safe_run(paste0("DE_celltype_", safe_file_name(ct)), {
        rr <- FindMarkers(obj, ident.1 = i2, ident.2 = i1, logfc.threshold = 0, min.pct = 0.1)
        rr$gene <- rownames(rr); rr$celltype <- ct; rr$contrast <- paste(g2, "vs", g1); rr
      }, default = data.frame())
      if (nrow(r) > 0) L[[ct]] <- r
    }
  }
  byct <- if (length(L) > 0) do.call(rbind, L) else data.frame()
  if (nrow(byct) > 0) {
    save_table(byct, "05_DE_by_celltype_FindMarkers.csv")
    dd <- byct %>% mutate(sig = case_when(p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "Up", p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "Down", TRUE ~ "NS"), neglogp = -log10(p_val_adj + 1e-300))
    save_pdf(ggplot(dd, aes(avg_log2FC, neglogp, color = sig)) + geom_point(size = 0.45, alpha = 0.6) +
               scale_color_manual(values = c(Up = "#EC7063", Down = "#5DADE2", NS = "grey82")) +
               facet_wrap(~celltype, scales = "free") + theme_bw(base_size = 10) +
               labs(title = paste0("Cell-type DE: ", g2, " vs ", g1), x = "Average log2FC", y = "-log10 adjusted p"), "14_volcano_by_celltype.pdf", 16, 12)
  }
  list(overall = overall, by_celltype = byct)
}

de_results <- if (isTRUE(CFG$run_de)) cache_run("06_DE_results", run_de(combined), default = list(overall = data.frame(), by_celltype = data.frame())) else list(overall = data.frame(), by_celltype = data.frame())

# ============================ CELL COMPOSITION ===============================
run_composition <- function(obj) {
  if (is.null(obj) || !inherits(obj, "Seurat")) stop("Input object is not a Seurat object")
  meta <- obj@meta.data
  sample_group <- meta %>% distinct(sample, group)
  tab <- as.data.frame(table(sample = meta$sample, celltype = meta$celltype), stringsAsFactors = FALSE)
  tab <- left_join(tab, sample_group, by = "sample")
  total <- tab %>% group_by(sample) %>% summarise(total = sum(Freq), .groups = "drop")
  tab <- left_join(tab, total, by = "sample") %>% mutate(prop = Freq / total)
  save_table(tab, "06_celltype_proportion_by_sample.csv")
  ptest <- tab %>% group_by(celltype) %>% summarise(p_value = tryCatch(wilcox.test(prop ~ group, data = cur_data())$p.value, error = function(e) NA_real_), .groups = "drop") %>% mutate(p_adj = p.adjust(p_value, method = "BH"))
  save_table(ptest, "07_celltype_proportion_wilcoxon.csv")
  grp_tab <- as.data.frame(prop.table(table(group = meta$group, celltype = meta$celltype), 1), stringsAsFactors = FALSE)
  save_pdf(ggplot(grp_tab, aes(group, Freq, fill = celltype)) + geom_col(color = "white", linewidth = 0.1) + scale_fill_manual(values = macaron_pal(length(unique(meta$celltype)))) + theme_classic(base_size = 13) + labs(title = "Cell-type composition", y = "Proportion", x = NULL), "15_celltype_composition_stacked.pdf", 8, 7)
  save_pdf(ggplot(tab, aes(group, prop, fill = group)) + geom_boxplot(outlier.shape = NA, alpha = 0.65) + geom_jitter(width = 0.15, size = 1.2) + facet_wrap(~celltype, scales = "free_y") + scale_fill_manual(values = macaron_pal(length(unique(tab$group)))) + theme_bw(base_size = 11) + theme(legend.position = "none") + labs(title = "Cell-type proportion per sample", y = "Proportion", x = NULL), "16_celltype_composition_boxplot.pdf", 14, 11)
  list(proportion = tab, test = ptest)
}
composition_results <- cache_run("07_cell_composition", run_composition(combined), default = NULL)

# ============================== PSEUDOBULK ===================================
run_pseudobulk_de <- function(obj) {
  if (is.null(obj) || !inherits(obj, "Seurat")) stop("Input object is not a Seurat object")
  if (!HAS("edgeR") || !HAS("limma")) {
    log_msg("edgeR or limma missing; pseudobulk DE skipped", level = "WARN")
    return(data.frame())
  }
  counts <- get_assay_layer(obj, assay = CFG$assay, layer = "counts")
  meta <- obj@meta.data
  out <- list()
  for (ct in unique(meta$celltype)) {
    cells <- rownames(meta)[meta$celltype == ct]
    if (length(cells) < 50) next
    md <- meta[cells, , drop = FALSE]
    groups <- table(md$sample, md$group)
    sample_group <- data.frame(sample = rownames(groups), group = colnames(groups)[max.col(groups, ties.method = "first")], stringsAsFactors = FALSE)
    if (length(unique(sample_group$group)) < 2) next
    if (min(table(sample_group$group)) < CFG$min_samples_for_pseudobulk) next
    samples <- unique(as.character(md$sample))
    mm <- do.call(cbind, lapply(samples, function(sid) Matrix::rowSums(counts[, rownames(md)[md$sample == sid], drop = FALSE])))
    colnames(mm) <- samples
    mm <- as.matrix(mm)
    sample_group <- sample_group[match(colnames(mm), sample_group$sample), ]
    dge <- edgeR::DGEList(counts = mm, group = sample_group$group)
    keep <- edgeR::filterByExpr(dge, group = sample_group$group)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    if (nrow(dge) < 50) next
    dge <- edgeR::calcNormFactors(dge)
    design <- model.matrix(~0 + group, data = sample_group)
    colnames(design) <- gsub("^group", "", colnames(design))
    dge <- edgeR::estimateDisp(dge, design)
    fit <- edgeR::glmQLFit(dge, design)
    glev <- colnames(design)
    contr <- paste0(glev[2], "-", glev[1])
    qlf <- edgeR::glmQLFTest(fit, contrast = limma::makeContrasts(contrasts = contr, levels = design))
    tt <- edgeR::topTags(qlf, n = Inf)$table
    tt$gene <- rownames(tt); tt$celltype <- ct; tt$contrast <- paste(glev[2], "vs", glev[1])
    out[[ct]] <- tt
  }
  res <- if (length(out) > 0) do.call(rbind, out) else data.frame()
  if (nrow(res) > 0) save_table(res, "08_pseudobulk_edgeR_by_celltype.csv")
  res
}
pseudobulk_results <- if (isTRUE(CFG$run_pseudobulk)) cache_run("08_pseudobulk_DE", run_pseudobulk_de(combined), default = data.frame()) else data.frame()

# =============================== ENRICHMENT ==================================
run_enrichment <- function(de_list_or_df, prefix = "DE") {
  if (!HAS("clusterProfiler") || is.null(orgdb_pkg) || !HAS(orgdb_pkg)) {
    log_msg("clusterProfiler or OrgDb missing; enrichment skipped", level = "WARN")
    return(data.frame())
  }
  suppressPackageStartupMessages(library(clusterProfiler))
  OrgDb <- get(orgdb_pkg, envir = asNamespace(orgdb_pkg))
  df <- if (is.list(de_list_or_df) && "by_celltype" %in% names(de_list_or_df)) de_list_or_df$by_celltype else de_list_or_df
  if (is.null(df) || nrow(df) == 0 || !"gene" %in% colnames(df)) return(data.frame())
  ct_col <- if ("celltype" %in% colnames(df)) "celltype" else NULL
  if (is.null(ct_col)) df$celltype <- "overall"
  out <- list()
  for (ct in unique(df$celltype)) {
    up <- df %>% filter(celltype == ct)
    padj_col <- intersect(c("p_val_adj", "FDR", "adj.P.Val"), colnames(up))[1]
    lfc_col <- intersect(c("avg_log2FC", "logFC", "log2FoldChange"), colnames(up))[1]
    if (is.na(padj_col) || is.na(lfc_col)) next
    genes <- up %>% filter(.data[[padj_col]] < 0.05, .data[[lfc_col]] > 0.25) %>% pull(gene) %>% unique()
    if (length(genes) < 10) next
    eg <- safe_run(paste0("bitr_", safe_file_name(ct)), clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb), default = data.frame())
    if (nrow(eg) < 5) next
    ego <- safe_run(paste0("enrichGO_", safe_file_name(ct)), clusterProfiler::enrichGO(eg$ENTREZID, OrgDb = OrgDb, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE), default = NULL)
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      tt <- as.data.frame(ego); tt$celltype <- ct; out[[ct]] <- tt
      p_go <- tryCatch(
        clusterProfiler::dotplot(ego, showCategory = CFG$go_show_category, label_format = CFG$go_label_wrap),
        error = function(e) clusterProfiler::dotplot(ego, showCategory = CFG$go_show_category)
      )
      p_go <- p_go + ggtitle(paste0(prefix, " GO BP: ", ct)) +
        theme_bw(base_size = 11) +
        theme(
          axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 9),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.margin = margin(8, 28, 8, 8)
        )
      save_pdf(p_go, paste0("17_GO_dotplot_", safe_file_name(ct), ".pdf"), 10, max(6, CFG$go_show_category * 0.45 + 3))
    }
  }
  res <- if (length(out) > 0) do.call(rbind, out) else data.frame()
  if (nrow(res) > 0) save_table(res, paste0("09_", prefix, "_GO_enrichment.csv"))
  res
}
enrichment_results <- if (isTRUE(CFG$run_enrichment)) cache_run("09_enrichment", run_enrichment(de_results, "DE"), default = data.frame()) else data.frame()

# ============================== PATHWAY SCORES ===============================
run_pathway_scores <- function(obj) {
  if (is.null(obj) || !inherits(obj, "Seurat")) stop("Input object is not a Seurat object")
  sigs <- list()
  if (species_norm == "human") {
    sigs <- list(
      Inflammation = c("IL1B", "TNF", "CXCL8", "CCL2", "NFKBIA", "PTGS2"),
      Interferon = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1"),
      Fibrosis = c("COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN", "ACTA2"),
      Hypoxia = c("HIF1A", "VEGFA", "LDHA", "ENO1", "SLC2A1", "CA9"),
      CellCycle = c("MKI67", "TOP2A", "PCNA", "STMN1", "UBE2C", "MCM5")
    )
  } else {
    sigs <- list(
      Inflammation = c("Il1b", "Tnf", "Cxcl2", "Ccl2", "Nfkbia", "Ptgs2"),
      Interferon = c("Isg15", "Ifit1", "Ifit2", "Ifit3", "Mx1", "Oas1a"),
      Fibrosis = c("Col1a1", "Col1a2", "Col3a1", "Fn1", "Postn", "Acta2"),
      Hypoxia = c("Hif1a", "Vegfa", "Ldha", "Eno1", "Slc2a1", "Car9"),
      CellCycle = c("Mki67", "Top2a", "Pcna", "Stmn1", "Ube2c", "Mcm5")
    )
  }
  sigs <- lapply(sigs, intersect, y = rownames(obj)); sigs <- sigs[vapply(sigs, length, integer(1)) >= 2]
  if (length(sigs) == 0) return(obj)
  obj <- AddModuleScore(obj, features = sigs, name = "PathwayScore_", search = FALSE)
  score_cols <- grep("^PathwayScore_", colnames(obj@meta.data), value = TRUE)
  if (length(score_cols) == 0) return(obj)
  names(score_cols) <- names(sigs)[seq_along(score_cols)]
  for (i in seq_along(score_cols)) obj[[paste0(names(score_cols)[i], "_score")]] <- obj@meta.data[[score_cols[i]]]
  long <- obj@meta.data %>% select(sample, group, celltype, ends_with("_score")) %>% pivot_longer(cols = ends_with("_score"), names_to = "signature", values_to = "score")
  save_table(long, "10_pathway_scores_long.csv")
  save_pdf(ggplot(long, aes(celltype, score, fill = group)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~signature, scales = "free_y") + theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = "Curated pathway/module scores", x = NULL, y = "Module score"), "18_pathway_scores_by_celltype.pdf", 16, 10)
  for (sig in names(sigs)) {
    col <- paste0(sig, "_score")
    save_pdf(FeaturePlot(obj, features = col, order = TRUE) + scale_color_gradientn(colors = macaron_grad) + ggtitle(sig), paste0("19_score_UMAP_", safe_file_name(sig), ".pdf"), 7, 6)
  }
  obj
}
if (isTRUE(CFG$run_pathway_scores)) combined <- cache_run("10_pathway_scores", run_pathway_scores(combined), default = combined)

# ========================== ADVANCED PSEUDOTIME ==============================
choose_trajectory_cells <- function(obj) {
  meta <- obj@meta.data
  if (!is.null(CFG$trajectory_celltypes)) {
    cells <- rownames(meta)[meta$celltype %in% CFG$trajectory_celltypes]
  } else {
    # Prefer large stromal/epithelial/endothelial groups for biologically meaningful state transitions.
    priority <- c("fibro", "smooth", "muscle", "endo", "epithel", "mono", "macro", "t cell", "b cell")
    ct_counts <- sort(table(meta$celltype), decreasing = TRUE)
    chosen <- names(ct_counts)[vapply(tolower(names(ct_counts)), function(ct) any(grepl(paste(priority, collapse = "|"), ct)), logical(1))]
    if (length(chosen) == 0) chosen <- names(ct_counts)[seq_len(min(3, length(ct_counts)))]
    cells <- rownames(meta)[meta$celltype %in% chosen]
  }
  if (length(cells) > CFG$max_cells_for_pseudotime) cells <- sample(cells, CFG$max_cells_for_pseudotime)
  cells
}

select_root_cells <- function(obj, cells) {
  meta <- obj@meta.data[cells, , drop = FALSE]
  root <- character()
  for (g in CFG$root_group_priority) {
    root <- rownames(meta)[tolower(meta$group) == tolower(g)]
    if (length(root) > 0) break
  }
  if (length(root) == 0) {
    key <- paste(CFG$root_celltype_keywords, collapse = "|")
    root <- rownames(meta)[grepl(key, tolower(meta$celltype))]
  }
  if (length(root) == 0) {
    # fallback: cells from smallest cluster id
    first_cluster <- sort(unique(as.character(meta$seurat_clusters)))[1]
    root <- rownames(meta)[as.character(meta$seurat_clusters) == first_cluster]
  }
  root <- intersect(root, cells)
  if (length(root) > 300) root <- sample(root, 300)
  root
}

plot_pseudotime_common <- function(df, prefix) {
  if (!all(c("Pseudotime", "celltype", "group") %in% colnames(df))) return(NULL)
  save_pdf(ggplot(df, aes(Pseudotime, fill = group)) + geom_density(alpha = 0.55, color = NA) + scale_fill_manual(values = macaron_pal(length(unique(df$group)))) + theme_bw(base_size = 13) + labs(title = paste0(prefix, " pseudotime density by group"), y = "Density"), paste0("20_", prefix, "_pseudotime_density_group.pdf"), 8, 5)
  save_pdf(ggplot(df, aes(celltype, Pseudotime, fill = celltype)) + geom_violin(scale = "width", color = "grey40") + geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.75) + scale_fill_manual(values = macaron_pal(length(unique(df$celltype)))) + theme_bw(base_size = 11) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = paste0(prefix, " pseudotime by cell type"), x = NULL), paste0("21_", prefix, "_pseudotime_violin_celltype.pdf"), 10, 6)
}

pseudotime_gene_dynamics <- function(obj, cells, pt, prefix) {
  pt <- pt[is.finite(pt)]
  common <- intersect(names(pt), cells)
  if (length(common) < 100) return(data.frame())
  vars <- intersect(VariableFeatures(obj), rownames(obj))
  if (length(vars) < 100) vars <- rownames(obj)[seq_len(min(3000, nrow(obj)))]
  expr <- get_assay_layer(obj, assay = DefaultAssay(obj), layer = "data")[vars, common, drop = FALSE]
  rho <- apply(as.matrix(expr), 1, function(x) suppressWarnings(cor(x, pt[common], method = "spearman", use = "pairwise.complete.obs")))
  pval <- apply(as.matrix(expr), 1, function(x) suppressWarnings(tryCatch(cor.test(x, pt[common], method = "spearman")$p.value, error = function(e) NA_real_)))
  res <- data.frame(gene = names(rho), spearman_rho = as.numeric(rho), p_value = as.numeric(pval), p_adj = p.adjust(pval, "BH")) %>% arrange(p_adj, desc(abs(spearman_rho)))
  save_table(res, paste0("11_", prefix, "_pseudotime_gene_correlation.csv"))
  top <- res %>% filter(is.finite(spearman_rho)) %>% arrange(desc(abs(spearman_rho))) %>% head(CFG$n_pseudotime_heatmap_genes) %>% pull(gene)
  top <- intersect(top, rownames(obj))
  if (length(top) >= 5 && HAS("pheatmap")) {
    ord_cells <- common[order(pt[common])]
    e <- as.matrix(get_assay_layer(obj, assay = DefaultAssay(obj), layer = "data")[top, ord_cells, drop = FALSE])
    bins <- cut(seq_along(ord_cells), breaks = CFG$n_pseudotime_bins, labels = FALSE)
    sm <- t(apply(e, 1, function(x) tapply(x, bins, mean, na.rm = TRUE)))
    sm <- t(scale(t(sm))); sm[!is.finite(sm)] <- 0
    save_pdf_expr(paste0("22_", prefix, "_pseudotime_gene_heatmap.pdf"), {
      pheatmap::pheatmap(sm, cluster_cols = FALSE, cluster_rows = TRUE, show_colnames = FALSE, color = macaron_grad, fontsize_row = ifelse(nrow(sm) <= 60, 6, 4), main = paste0(prefix, " dynamic genes"))
    }, 9, 11)
  }
  # Representative gene smooth curves
  curve_genes <- head(top, 12)
  if (length(curve_genes) >= 3) {
    edf <- as.data.frame(t(as.matrix(get_assay_layer(obj, assay = DefaultAssay(obj), layer = "data")[curve_genes, common, drop = FALSE])))
    edf$Pseudotime <- pt[common]
    edf$cell <- common
    long <- edf %>% pivot_longer(cols = all_of(curve_genes), names_to = "gene", values_to = "expression")
    save_pdf(ggplot(long, aes(Pseudotime, expression)) + geom_point(alpha = 0.15, size = 0.25) + geom_smooth(method = "loess", se = FALSE, linewidth = 0.7) + facet_wrap(~gene, scales = "free_y", ncol = 4) + theme_bw(base_size = 10) + labs(title = paste0(prefix, " genes along pseudotime")), paste0("23_", prefix, "_gene_curves.pdf"), 12, 8)
  }
  res
}

run_monocle3_pt <- function(obj, cells, root_cells) {
  if (!HAS("monocle3") || !HAS("SingleCellExperiment")) stop("monocle3 or SingleCellExperiment missing")
  suppressPackageStartupMessages(library(monocle3))
  sub <- subset(obj, cells = cells)
  counts <- get_assay_layer(sub, assay = CFG$assay, layer = "counts")
  gene_meta <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
  cds <- monocle3::new_cell_data_set(counts, cell_metadata = sub@meta.data, gene_metadata = gene_meta)
  cds <- monocle3::preprocess_cds(cds, num_dim = max(CFG$dims))
  cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP")
  cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP")
  cds <- monocle3::learn_graph(cds, use_partition = TRUE)
  root_cells <- intersect(root_cells, colnames(cds))
  if (length(root_cells) > 0) cds <- monocle3::order_cells(cds, root_cells = root_cells) else cds <- monocle3::order_cells(cds)
  pt <- monocle3::pseudotime(cds)
  md <- as.data.frame(colData(cds)); md$Pseudotime <- pt[rownames(md)]; md$cell <- rownames(md)
  save_table(md, "12_monocle3_pseudotime_metadata.csv")
  save_pdf(monocle3::plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = FALSE, label_leaves = TRUE, label_branch_points = TRUE) + ggtitle("Monocle3 pseudotime"), "24_monocle3_trajectory_pseudotime.pdf", 8, 7)
  save_pdf(monocle3::plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster = FALSE, label_leaves = TRUE, label_branch_points = TRUE) + ggtitle("Monocle3 trajectory by cell type"), "25_monocle3_trajectory_celltype.pdf", 9, 7)
  save_pdf(monocle3::plot_cells(cds, color_cells_by = "group", label_groups_by_cluster = FALSE, label_leaves = TRUE, label_branch_points = TRUE) + ggtitle("Monocle3 trajectory by group"), "26_monocle3_trajectory_group.pdf", 8, 7)
  plot_pseudotime_common(md, "monocle3")
  dyn <- pseudotime_gene_dynamics(obj, cells, pt, "monocle3")
  list(cds = cds, metadata = md, dynamic_genes = dyn)
}

run_slingshot_pt <- function(obj, cells, root_cells) {
  if (!HAS("slingshot") || !HAS("SingleCellExperiment")) stop("slingshot or SingleCellExperiment missing")
  suppressPackageStartupMessages({ library(SingleCellExperiment); library(slingshot) })
  sub <- subset(obj, cells = cells)
  counts <- get_assay_layer(sub, assay = CFG$assay, layer = "counts")
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts, logcounts = get_assay_layer(sub, assay = DefaultAssay(sub), layer = "data")), colData = sub@meta.data)
  if ("umap" %in% names(sub@reductions)) reducedDims(sce)$UMAP <- Embeddings(sub, "umap")
  if (!"UMAP" %in% names(reducedDims(sce))) stop("UMAP reduction missing for Slingshot")
  start.clus <- NULL
  if (length(root_cells) > 0) {
    root_clusters <- sub$seurat_clusters[match(root_cells, colnames(sub))]
    start.clus <- names(sort(table(root_clusters), decreasing = TRUE))[1]
  }
  sce <- slingshot::slingshot(sce, clusterLabels = "seurat_clusters", reducedDim = "UMAP", start.clus = start.clus)
  ptm <- slingshot::slingPseudotime(sce)
  pt <- apply(ptm, 1, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))
  md <- as.data.frame(colData(sce)); md$Pseudotime <- pt[rownames(md)]; md$cell <- rownames(md)
  save_table(md, "13_slingshot_pseudotime_metadata.csv")
  plot_df <- as.data.frame(reducedDims(sce)$UMAP); colnames(plot_df) <- c("UMAP_1", "UMAP_2"); plot_df$Pseudotime <- md$Pseudotime; plot_df$celltype <- md$celltype; plot_df$group <- md$group
  save_pdf(ggplot(plot_df, aes(UMAP_1, UMAP_2, color = Pseudotime)) + geom_point(size = 0.45) + scale_color_gradientn(colors = macaron_grad, na.value = "grey80") + theme_classic() + ggtitle("Slingshot pseudotime"), "27_slingshot_trajectory_pseudotime.pdf", 8, 7)
  save_pdf(ggplot(plot_df, aes(UMAP_1, UMAP_2, color = celltype)) + geom_point(size = 0.45) + scale_color_manual(values = macaron_pal(length(unique(plot_df$celltype)))) + theme_classic() + ggtitle("Slingshot by cell type"), "28_slingshot_trajectory_celltype.pdf", 9, 7)
  plot_pseudotime_common(md, "slingshot")
  dyn <- pseudotime_gene_dynamics(obj, cells, pt, "slingshot")
  list(sce = sce, metadata = md, dynamic_genes = dyn)
}

run_monocle2_pt <- function(obj, cells, root_cells) {
  if (!HAS("monocle")) stop("monocle2 package missing")
  suppressPackageStartupMessages(library(monocle))
  sub <- subset(obj, cells = cells)
  ex <- get_assay_layer(sub, assay = CFG$assay, layer = "counts")
  pd <- new("AnnotatedDataFrame", data = sub@meta.data)
  fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(ex), row.names = rownames(ex)))
  cds <- monocle::newCellDataSet(as(ex, "sparseMatrix"), phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, expressionFamily = monocle::negbinomial.size())
  cds <- monocle::estimateSizeFactors(cds)
  cds <- tryCatch(monocle::detectGenes(cds, min_expr = 0.1), error = function(e) cds)
  
  # Try official dispersion estimation first. New dplyr versions may break old Monocle2 internals;
  # if that happens, keep the pipeline running and still do DDRTree + correlation fallback.
  cds <- safe_run("monocle2_estimateDispersions", {
    monocle::estimateDispersions(cds)
  }, default = cds)
  
  og <- intersect(VariableFeatures(obj), rownames(cds))
  if (length(og) < 100) og <- rownames(cds)[seq_len(min(2000, nrow(cds)))]
  cds <- monocle::setOrderingFilter(cds, og)
  cds <- monocle::reduceDimension(cds, max_components = 2, reduction_method = "DDRTree", norm_method = "log")
  cds <- tryCatch(monocle::orderCells(cds), error = function(e) {
    log_msg("Monocle2 orderCells failed; using DDRTree-coordinate fallback pseudotime: ", conditionMessage(e), level = "WARN")
    rd <- t(monocle::reducedDimS(cds)); pseudo <- rd[, 1]
    pseudo <- pseudo - min(pseudo, na.rm = TRUE)
    if (max(pseudo, na.rm = TRUE) > 0) pseudo <- pseudo / max(pseudo, na.rm = TRUE)
    Biobase::pData(cds)$Pseudotime <- pseudo
    Biobase::pData(cds)$State <- factor(kmeans(rd, centers = min(3, nrow(rd)), nstart = 10)$cluster)
    cds
  })
  
  # Orient trajectory toward root cells when possible.
  root_cells <- intersect(root_cells, colnames(cds))
  if (length(root_cells) > 0 && "State" %in% colnames(Biobase::pData(cds))) {
    root_state <- names(sort(table(Biobase::pData(cds)[root_cells, "State"]), decreasing = TRUE))[1]
    cds <- tryCatch(monocle::orderCells(cds, root_state = root_state), error = function(e) cds)
  }
  
  pdat <- Biobase::pData(cds); pdat$cell <- rownames(pdat)
  save_table(pdat, "14_monocle2_pseudotime_metadata.csv")
  save_pdf(monocle::plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.6) + scale_color_gradientn(colors = macaron_grad) + ggtitle("Monocle2 DDRTree pseudotime"), "29_monocle2_trajectory_pseudotime.pdf", 8, 7)
  save_pdf(monocle::plot_cell_trajectory(cds, color_by = "State", cell_size = 0.6) + scale_color_manual(values = macaron_pal(length(unique(pdat$State)))) + ggtitle("Monocle2 trajectory states"), "29A_monocle2_trajectory_state.pdf", 8, 7)
  save_pdf(monocle::plot_cell_trajectory(cds, color_by = "celltype", cell_size = 0.6) + scale_color_manual(values = macaron_pal(length(unique(pdat$celltype)))) + ggtitle("Monocle2 by cell type"), "30_monocle2_trajectory_celltype.pdf", 9, 7)
  save_pdf(monocle::plot_cell_trajectory(cds, color_by = "group", cell_size = 0.6) + scale_color_manual(values = macaron_pal(length(unique(pdat$group)))) + ggtitle("Monocle2 by group"), "30A_monocle2_trajectory_group.pdf", 8, 7)
  plot_pseudotime_common(pdat, "monocle2")
  
  pt <- setNames(pdat$Pseudotime, rownames(pdat))
  dyn <- pseudotime_gene_dynamics(obj, cells, pt, "monocle2")
  
  official_pt <- data.frame()
  beam_res <- data.frame()
  
  if (isTRUE(CFG$monocle2_official_gene_test)) {
    test_genes <- intersect(og, rownames(cds))
    if (length(test_genes) > CFG$monocle2_max_test_genes) {
      # Prefer variable genes with non-zero expression and keep runtime bounded.
      test_genes <- test_genes[seq_len(CFG$monocle2_max_test_genes)]
    }
    official_pt <- safe_run("monocle2_differentialGeneTest_smns_Pseudotime", {
      dg <- monocle::differentialGeneTest(cds[test_genes, ], fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = CFG$monocle2_num_cores)
      dg <- dg[order(dg$qval, dg$pval), , drop = FALSE]
      dg$gene <- rownames(dg)
      save_table(dg, "14A_monocle2_differentialGeneTest_smns_Pseudotime.csv")
      dg
    }, default = data.frame())
    
    if (nrow(official_pt) > 0) {
      top_pt <- official_pt %>% filter(is.finite(qval)) %>% arrange(qval) %>% head(CFG$n_pseudotime_heatmap_genes) %>% pull(gene)
      top_pt <- intersect(top_pt, rownames(cds))
      if (length(top_pt) >= 5) {
        # Official Monocle2 heatmap first; if it breaks, use the self-drawn fallback already made above.
        save_pdf_expr("29B_monocle2_official_plot_pseudotime_heatmap.pdf", {
          monocle::plot_pseudotime_heatmap(cds[top_pt, ], num_clusters = 4, cores = CFG$monocle2_num_cores, show_rownames = FALSE, return_heatmap = FALSE)
          invisible(NULL)
        }, 9, 11)
        curve_genes <- head(top_pt, 9)
        save_pdf(monocle::plot_genes_in_pseudotime(cds[curve_genes, ], color_by = "celltype", ncol = 3) + ggtitle("Top Monocle2 pseudotime genes"), "29C_monocle2_official_genes_in_pseudotime.pdf", 12, 9)
      }
    }
  }
  
  if (isTRUE(CFG$monocle2_run_beam)) {
    branch_point <- CFG$monocle2_branch_point
    beam_genes <- intersect(og, rownames(cds))
    if (length(beam_genes) > CFG$monocle2_max_test_genes) beam_genes <- beam_genes[seq_len(CFG$monocle2_max_test_genes)]
    beam_res <- safe_run(paste0("monocle2_BEAM_branch_point_", branch_point), {
      br <- monocle::BEAM(cds[beam_genes, ], branch_point = branch_point, cores = CFG$monocle2_num_cores)
      br <- br[order(br$qval, br$pval), , drop = FALSE]
      br$gene <- rownames(br)
      save_table(br, paste0("14B_monocle2_BEAM_branch_point_", branch_point, ".csv"))
      br
    }, default = data.frame())
    
    if (nrow(beam_res) > 0) {
      top_beam <- beam_res %>% filter(is.finite(qval)) %>% arrange(qval) %>% head(CFG$n_pseudotime_heatmap_genes) %>% pull(gene)
      top_beam <- intersect(top_beam, rownames(cds))
      if (length(top_beam) >= 5) {
        save_pdf_expr(paste0("29D_monocle2_BEAM_branched_heatmap_bp", branch_point, ".pdf"), {
          monocle::plot_genes_branched_heatmap(cds[top_beam, ], branch_point = branch_point, num_clusters = 4, cores = CFG$monocle2_num_cores, show_rownames = FALSE, return_heatmap = FALSE)
          invisible(NULL)
        }, 10, 12)
        save_pdf(monocle::plot_genes_branched_pseudotime(cds[head(top_beam, 6), ], branch_point = branch_point, color_by = "celltype", ncol = 3) + ggtitle(paste0("BEAM branch genes, branch point ", branch_point)), paste0("29E_monocle2_BEAM_branched_gene_curves_bp", branch_point, ".pdf"), 12, 8)
      }
    }
  }
  
  list(cds = cds, metadata = pdat, dynamic_genes = dyn, official_differentialGeneTest = official_pt, BEAM = beam_res)
}

run_pseudotime_all <- function(obj) {
  if (is.null(obj) || !inherits(obj, "Seurat")) stop("Input object is not a Seurat object")
  cells <- choose_trajectory_cells(obj)
  if (length(cells) < 100) stop("Too few cells for pseudotime analysis")
  root_cells <- select_root_cells(obj, cells)
  save_table(data.frame(cell = cells, root = cells %in% root_cells, celltype = obj$celltype[cells], group = obj$group[cells]), "11_pseudotime_selected_cells.csv")
  save_pdf(DimPlot(obj, cells.highlight = cells, cols.highlight = "#EC7063", sizes.highlight = 0.35) + NoLegend() + ggtitle(paste0("Trajectory selected cells: n=", length(cells))), "20_pseudotime_selected_cells.pdf", 8, 7)
  res <- list(selected_cells = cells, root_cells = root_cells)
  if ("monocle3" %in% CFG$trajectory_methods) res$monocle3 <- cache_run("11_pseudotime_monocle3", run_monocle3_pt(obj, cells, root_cells), default = NULL)
  if ("slingshot" %in% CFG$trajectory_methods) res$slingshot <- cache_run("12_pseudotime_slingshot", run_slingshot_pt(obj, cells, root_cells), default = NULL)
  if ("monocle2" %in% CFG$trajectory_methods) res$monocle2 <- cache_run("13_pseudotime_monocle2", run_monocle2_pt(obj, cells, root_cells), default = NULL)
  res
}

pseudotime_results <- if (isTRUE(CFG$run_pseudotime)) cache_run("11_pseudotime_all", run_pseudotime_all(combined), default = list()) else list()

# =============================== CELLCHAT ====================================
prepare_cellchat_object <- function(obj) {
  meta <- obj@meta.data
  # Select top cell types to keep network readable and computationally stable.
  ct_keep <- names(sort(table(meta$celltype), decreasing = TRUE))[seq_len(min(CFG$max_celltypes_cellchat, length(unique(meta$celltype))))]
  cells <- rownames(meta)[meta$celltype %in% ct_keep]
  # Stratified downsampling per cell type.
  cells2 <- unlist(lapply(split(cells, meta[cells, "celltype"]), function(v) if (length(v) > CFG$max_cells_per_celltype_cellchat) sample(v, CFG$max_cells_per_celltype_cellchat) else v), use.names = FALSE)
  if (length(cells2) > CFG$max_cells_for_cellchat_total) cells2 <- sample(cells2, CFG$max_cells_for_cellchat_total)
  subset(obj, cells = cells2)
}

run_cellchat_one <- function(obj, label = "all") {
  if (!HAS("CellChat")) stop("CellChat missing")
  suppressPackageStartupMessages(library(CellChat))
  sub <- prepare_cellchat_object(obj)
  if (length(unique(sub$celltype)) < 2) stop("Need at least two cell types for CellChat")
  if (min(table(sub$celltype)) < CFG$cellchat_min_cells) log_msg("Some cell types have fewer than cellchat_min_cells; CellChat will filter weak groups", level = "WARN")
  data.input <- get_assay_layer(sub, assay = CFG$assay, layer = "data")
  meta <- sub@meta.data
  cellchat <- CellChat::createCellChat(object = data.input, meta = meta, group.by = "celltype")
  db.use <- if (species_norm == "mouse") CellChat::CellChatDB.mouse else CellChat::CellChatDB.human
  if (!is.null(CFG$cellchat_database_category)) db.use <- CellChat::subsetDB(db.use, search = CFG$cellchat_database_category)
  cellchat@DB <- db.use
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = CFG$cellchat_min_cells)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)
  cellchat <- CellChat::netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  comm <- CellChat::subsetCommunication(cellchat)
  save_table(comm, paste0("15_cellchat_communication_", safe_file_name(label), ".csv"))
  group_size <- as.numeric(table(cellchat@idents))
  save_pdf_expr(paste0("31_cellchat_circle_count_", safe_file_name(label), ".pdf"), {
    CellChat::netVisual_circle(cellchat@net$count, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, title.name = paste0("CellChat interaction count: ", label))
  }, 8, 8)
  save_pdf_expr(paste0("32_cellchat_circle_weight_", safe_file_name(label), ".pdf"), {
    CellChat::netVisual_circle(cellchat@net$weight, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, title.name = paste0("CellChat interaction strength: ", label))
  }, 8, 8)
  # CellChat heatmaps are often ComplexHeatmap/grid objects, not ggplot objects.
  # Do not add ggtitle() directly; draw with draw_plot_object() inside a tightly controlled PDF device.
  save_pdf_expr(paste0("33_cellchat_heatmap_", safe_file_name(label), ".pdf"), {
    ht <- CellChat::netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Reds")
    draw_plot_object(ht)
    invisible(NULL)
  }, 9, 8)
  
  lr_col <- intersect(c("interaction_name", "interaction_name_2"), colnames(comm))[1]
  pairLR.use <- NULL
  if (!is.na(lr_col) && "prob" %in% colnames(comm)) {
    comm_top <- comm[order(comm$prob, decreasing = TRUE), , drop = FALSE]
    top_lr <- unique(as.character(comm_top[[lr_col]]))
    top_lr <- top_lr[nzchar(top_lr) & !is.na(top_lr)]
    top_lr <- head(top_lr, CFG$cellchat_top_interactions)
    pairLR.use <- data.frame(interaction_name = top_lr, stringsAsFactors = FALSE)
    save_table(head(comm_top, CFG$cellchat_top_interactions), paste0("15A_cellchat_top", CFG$cellchat_top_interactions, "_communication_", safe_file_name(label), ".csv"))
  }
  bubble_plot <- tryCatch({
    if (!is.null(pairLR.use) && nrow(pairLR.use) > 0) {
      CellChat::netVisual_bubble(cellchat, pairLR.use = pairLR.use, remove.isolate = TRUE)
    } else {
      CellChat::netVisual_bubble(cellchat, remove.isolate = TRUE)
    }
  }, error = function(e) {
    log_msg("CellChat top-pair bubble failed; fallback to default bubble: ", conditionMessage(e), level = "WARN")
    CellChat::netVisual_bubble(cellchat, remove.isolate = TRUE)
  })
  if (inherits(bubble_plot, c("gg", "ggplot"))) bubble_plot <- bubble_plot + ggtitle(paste0("CellChat top ", CFG$cellchat_top_interactions, " ligand-receptor pairs: ", label)) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  save_pdf(bubble_plot, paste0("34_cellchat_bubble_", safe_file_name(label), ".pdf"), 12, 8)
  
  save_pdf_expr(paste0("35_cellchat_signaling_role_heatmap_", safe_file_name(label), ".pdf"), {
    ht2 <- CellChat::netAnalysis_signalingRole_heatmap(cellchat, pattern = "all")
    draw_plot_object(ht2)
    invisible(NULL)
  }, 10, 8)
  
  pathways <- tryCatch(cellchat@netP$pathways, error = function(e) character())
  if (length(pathways) > 0) {
    for (pw in head(pathways, CFG$cellchat_top_pathways)) {
      save_pdf_expr(paste0("36_cellchat_pathway_circle_", safe_file_name(label), "_", safe_file_name(pw), ".pdf"), {
        CellChat::netVisual_aggregate(cellchat, signaling = pw, layout = "circle")
      }, 8, 8)
    }
  }
  cellchat
}

run_cellchat_all <- function(obj) {
  if (is.null(obj) || !inherits(obj, "Seurat")) stop("Input object is not a Seurat object")
  res <- list()
  res$all <- cache_run("14_cellchat_all", run_cellchat_one(obj, "all"), default = NULL)
  groups <- setdiff(unique(as.character(obj$group)), CFG$unknown_group_name)
  if (length(groups) >= 2) {
    group_objs <- list()
    for (g in groups) {
      if (sum(obj$group == g) >= 200) {
        group_objs[[g]] <- cache_run(paste0("15_cellchat_group_", safe_file_name(g)), run_cellchat_one(subset(obj, subset = group == g), g), default = NULL)
      }
    }
    group_objs <- group_objs[!vapply(group_objs, is.null, logical(1))]
    res$by_group <- group_objs
    if (length(group_objs) >= 2 && HAS("CellChat")) {
      merged <- safe_run("CellChat_compare_groups", {
        cellchat_merged <- CellChat::mergeCellChat(group_objs, add.names = names(group_objs))
        save_pdf_expr("37_cellchat_compare_interactions_count.pdf", {
          print(CellChat::compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2), measure = "count"))
        }, 7, 6)
        save_pdf_expr("38_cellchat_compare_interactions_weight.pdf", {
          print(CellChat::compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2), measure = "weight"))
        }, 7, 6)
        save_pdf_expr("39_cellchat_diff_interaction_count.pdf", {
          CellChat::netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "count")
        }, 8, 8)
        save_pdf_expr("40_cellchat_diff_interaction_weight.pdf", {
          CellChat::netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "weight")
        }, 8, 8)
        cellchat_merged
      }, default = NULL)
      res$merged_compare <- merged
    }
  }
  res
}
cellchat_results <- if (isTRUE(CFG$run_cellchat)) cache_run("14_cellchat_results", run_cellchat_all(combined), default = list()) else list()

# ================================ hdWGCNA ====================================
# v4 design notes:
# 1) Primary strategy follows the hdWGCNA tutorial pattern: run metacell construction
#    on the full Seurat object grouped by celltype + sample, then use SetDatExpr()
#    for the target cell type. This preserves sample structure and reduces sparsity.
# 2) If the primary strategy fails because a GEO dataset has unusual metadata or too
#    few target cells in some samples, fallback to a target-celltype subset strategy.
# 3) We use argument-introspection so different hdWGCNA versions using soft_power / power
#    or slightly different formal arguments do not break the whole pipeline.

call_pkg_fun <- function(pkg, fun, args = list()) {
  f <- get(fun, envir = asNamespace(pkg))
  fn <- names(formals(f))
  args <- args[names(args) %in% fn]
  do.call(f, args)
}

plot_hdwgcna_object <- function(p, file, w = 9, h = 6) {
  save_pdf_expr(file, {
    draw_plot_object(p)
    invisible(NULL)
  }, w, h)
}

safe_get_modules <- function(obj) {
  tryCatch(hdWGCNA::GetModules(obj), error = function(e) data.frame())
}

safe_get_hubs <- function(obj, n_hubs = 30) {
  tryCatch(hdWGCNA::GetHubGenes(obj, n_hubs = n_hubs), error = function(e) data.frame())
}

safe_plot_hdwgcna <- function(expr, file, w = 9, h = 6) {
  save_pdf_expr(file, {
    p <- tryCatch(force(expr), error = function(e) {
      log_msg("hdWGCNA figure failed: ", file, " | ", conditionMessage(e), level = "WARN")
      NULL
    })
    if (!is.null(p)) draw_plot_object(p)
    invisible(NULL)
  }, w, h)
}


# Publication-style hdWGCNA plotting helpers. These do not replace the official
# hdWGCNA plots; they add readable summaries and hide overly crowded fallback
# dendrograms. All functions are fail-soft.
hdwgcna_available <- function(fun) {
  HAS("hdWGCNA") && exists(fun, envir = asNamespace("hdWGCNA"), inherits = FALSE)
}

get_top_modules <- function(mods, n = CFG$hdwgcna_max_modules_to_plot) {
  if (is.null(mods) || nrow(mods) == 0 || !"module" %in% colnames(mods)) return(character())
  tab <- sort(table(mods$module), decreasing = TRUE)
  names(tab) <- as.character(names(tab))
  names(tab)[!tolower(names(tab)) %in% c("grey", "gray")][seq_len(min(n, sum(!tolower(names(tab)) %in% c("grey", "gray"))))]
}

plot_module_size_bar <- function(mods, label, prefix = "48_hdwgcna_module_size") {
  if (is.null(mods) || nrow(mods) == 0 || !"module" %in% colnames(mods)) return(invisible(FALSE))
  df <- as.data.frame(table(module = mods$module), stringsAsFactors = FALSE)
  df <- df[!tolower(df$module) %in% c("grey", "gray"), , drop = FALSE]
  if (nrow(df) == 0) return(invisible(FALSE))
  df <- df[order(df$Freq, decreasing = TRUE), , drop = FALSE]
  df$module <- factor(df$module, levels = rev(df$module))
  p <- ggplot(df, aes(module, Freq, fill = module)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.2) +
    coord_flip() +
    scale_fill_manual(values = macaron_pal(nrow(df))) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none", axis.text.y = element_text(size = 9), plot.title = element_text(face = "bold", hjust = 0.5)) +
    labs(title = paste0("hdWGCNA module sizes: ", label), x = NULL, y = "Genes per module")
  save_pdf(p, paste0(prefix, "_", safe_file_name(label), ".pdf"), 8, max(5, 0.35 * nrow(df) + 2))
}

plot_hub_gene_dotplot <- function(hubs, label, prefix = "49_hdwgcna_top_hub_genes") {
  if (is.null(hubs) || nrow(hubs) == 0) return(invisible(FALSE))
  gene_col <- intersect(c("gene_name", "gene", "Gene", "symbol"), colnames(hubs))[1]
  module_col <- intersect(c("module", "Module", "color"), colnames(hubs))[1]
  kme_col <- intersect(c("kME", "kME_table", "kME_value", "eigengene_connectivity", "connectivity"), colnames(hubs))[1]
  if (is.na(gene_col) || is.na(module_col)) return(invisible(FALSE))
  if (is.na(kme_col)) {
    numeric_cols <- names(hubs)[vapply(hubs, is.numeric, logical(1))]
    kme_col <- numeric_cols[1]
  }
  if (is.na(kme_col)) return(invisible(FALSE))
  df <- hubs
  df$gene_plot <- as.character(df[[gene_col]])
  df$module_plot <- as.character(df[[module_col]])
  df$kME_plot <- as.numeric(df[[kme_col]])
  df <- df[is.finite(df$kME_plot) & !tolower(df$module_plot) %in% c("grey", "gray"), , drop = FALSE]
  if (nrow(df) == 0) return(invisible(FALSE))
  df <- do.call(rbind, lapply(split(df, df$module_plot), function(z) {
    z <- z[order(abs(z$kME_plot), decreasing = TRUE), , drop = FALSE]
    head(z, CFG$hdwgcna_hub_genes_per_module)
  }))
  df$gene_plot <- factor(df$gene_plot, levels = rev(unique(df$gene_plot)))
  p <- ggplot(df, aes(x = abs(kME_plot), y = gene_plot, color = module_plot, size = abs(kME_plot))) +
    geom_point(alpha = 0.95) +
    facet_wrap(~module_plot, scales = "free_y", ncol = 2) +
    scale_color_manual(values = macaron_pal(length(unique(df$module_plot)))) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", strip.text = element_text(face = "bold"), axis.text.y = element_text(size = 8), plot.title = element_text(face = "bold", hjust = 0.5)) +
    labs(title = paste0("Top hdWGCNA hub genes: ", label), x = "|kME|", y = NULL)
  save_pdf(p, paste0(prefix, "_", safe_file_name(label), ".pdf"), 10, max(6, 0.22 * nrow(df) + 3))
}

plot_me_by_group <- function(sw, label) {
  mes <- tryCatch(hdWGCNA::GetMEs(sw), error = function(e) NULL)
  if (is.null(mes) || nrow(mes) < 2) return(invisible(FALSE))
  md <- sw@meta.data[intersect(rownames(mes), rownames(sw@meta.data)), , drop = FALSE]
  mes <- mes[rownames(md), , drop = FALSE]
  if (!"group" %in% colnames(md)) return(invisible(FALSE))
  df <- cbind(md[, intersect(c("sample", "group", "celltype"), colnames(md)), drop = FALSE], as.data.frame(mes))
  me_cols <- colnames(mes)[seq_len(min(ncol(mes), CFG$hdwgcna_max_modules_to_plot))]
  if (length(me_cols) == 0) return(invisible(FALSE))
  long <- df %>% tidyr::pivot_longer(cols = all_of(me_cols), names_to = "module", values_to = "eigengene")
  p <- ggplot(long, aes(group, eigengene, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, linewidth = 0.3) +
    geom_jitter(width = 0.15, size = 0.6, alpha = 0.45) +
    facet_wrap(~module, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = macaron_pal(length(unique(long$group)))) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1), strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", hjust = 0.5)) +
    labs(title = paste0("Module eigengenes by group: ", label), x = NULL, y = "Module eigengene")
  save_pdf(p, paste0("50_hdwgcna_module_eigengene_by_group_", safe_file_name(label), ".pdf"), 11, 8)
}

plot_hdwgcna_official_networks <- function(sw, label, mods) {
  top_mods <- get_top_modules(mods, CFG$hdwgcna_max_modules_to_plot)
  if (length(top_mods) == 0) return(invisible(FALSE))
  # Official UMAP of module genes, when available.
  if (hdwgcna_available("ModuleUMAPPlot")) {
    safe_plot_hdwgcna(call_pkg_fun("hdWGCNA", "ModuleUMAPPlot", list(seurat_obj = sw, wgcna_name = attr(sw, "hdwgcna_wgcna_name"))), paste0("51_hdwgcna_module_umap_", safe_file_name(label), ".pdf"), 8, 7)
  }
  # Official hub-gene network; cap modules to avoid spaghetti plots.
  if (hdwgcna_available("HubGeneNetworkPlot")) {
    safe_plot_hdwgcna(call_pkg_fun("hdWGCNA", "HubGeneNetworkPlot", list(seurat_obj = sw, mods = top_mods, n_hubs = min(6, CFG$hdwgcna_hub_genes_per_module), n_other = 2, edge_prop = 0.35, vertex.label.cex = 0.55, wgcna_name = attr(sw, "hdwgcna_wgcna_name"))), paste0("52_hdwgcna_hub_gene_network_", safe_file_name(label), ".pdf"), 10, 8)
  }
  # Per-module network plots, only top modules to keep readable.
  if (hdwgcna_available("ModuleNetworkPlot")) {
    for (m in head(top_mods, 4)) {
      safe_plot_hdwgcna(call_pkg_fun("hdWGCNA", "ModuleNetworkPlot", list(seurat_obj = sw, mods = m, n_inner = 10, n_outer = 15, wgcna_name = attr(sw, "hdwgcna_wgcna_name"))), paste0("53_hdwgcna_module_network_", safe_file_name(label), "_", safe_file_name(m), ".pdf"), 8, 8)
    }
  }
  if (hdwgcna_available("ModuleTopologyBarplot")) {
    safe_plot_hdwgcna(call_pkg_fun("hdWGCNA", "ModuleTopologyBarplot", list(seurat_obj = sw, mods = top_mods, wgcna_name = attr(sw, "hdwgcna_wgcna_name"))), paste0("54_hdwgcna_module_topology_barplot_", safe_file_name(label), ".pdf"), 10, 7)
  }
  if (hdwgcna_available("ModuleTopologyHeatmap")) {
    safe_plot_hdwgcna(call_pkg_fun("hdWGCNA", "ModuleTopologyHeatmap", list(seurat_obj = sw, mods = top_mods, wgcna_name = attr(sw, "hdwgcna_wgcna_name"))), paste0("55_hdwgcna_module_topology_heatmap_", safe_file_name(label), ".pdf"), 9, 8)
  }
}

choose_hdwgcna_celltypes <- function(obj) {
  meta <- obj@meta.data
  if (!"celltype" %in% colnames(meta)) stop("celltype column missing; run annotation before hdWGCNA")
  if (!is.null(CFG$hdwgcna_target_celltypes)) {
    targets <- intersect(CFG$hdwgcna_target_celltypes, unique(as.character(meta$celltype)))
    if (length(targets) == 0) stop("None of hdwgcna_target_celltypes exist in object$celltype")
    return(head(targets, CFG$max_celltypes_hdwgcna))
  }
  ct_counts <- sort(table(meta$celltype), decreasing = TRUE)
  ct_counts <- ct_counts[!grepl("unknown|doublet|low|ambient|red blood|eryth|platelet", tolower(names(ct_counts)))]
  if (length(ct_counts) == 0) ct_counts <- sort(table(meta$celltype), decreasing = TRUE)
  
  sample_cov <- sapply(names(ct_counts), function(ct) {
    tab <- table(meta$sample[meta$celltype == ct])
    sum(tab >= CFG$hdwgcna_min_cells_per_sample_ct)
  })
  eligible <- names(ct_counts)[as.numeric(ct_counts) >= CFG$min_cells_per_celltype_hdwgcna & sample_cov >= CFG$hdwgcna_min_samples_per_ct]
  if (length(eligible) == 0) {
    eligible <- names(ct_counts)[as.numeric(ct_counts) >= CFG$min_cells_per_celltype_hdwgcna]
  }
  if (length(eligible) == 0) stop("No cell type has enough cells/sample coverage for hdWGCNA")
  
  key <- paste(CFG$hdwgcna_priority_celltype_keywords, collapse = "|")
  priority_bonus <- ifelse(grepl(key, tolower(eligible)), 10, 1)
  score <- log1p(as.numeric(ct_counts[eligible])) * log1p(sample_cov[eligible]) * priority_bonus
  targets <- eligible[order(score, decreasing = TRUE)]
  
  sel_table <- data.frame(
    celltype = eligible,
    cells = as.numeric(ct_counts[eligible]),
    samples_with_enough_cells = as.numeric(sample_cov[eligible]),
    priority_bonus = priority_bonus,
    selection_score = as.numeric(score),
    stringsAsFactors = FALSE
  )
  sel_table <- sel_table[order(sel_table$selection_score, decreasing = TRUE), ]
  save_table(sel_table, "16_hdwgcna_target_celltype_selection.csv")
  head(targets, CFG$max_celltypes_hdwgcna)
}

choose_hdwgcna_reduction <- function(obj) {
  red <- if ("harmony" %in% names(obj@reductions)) "harmony" else if ("pca" %in% names(obj@reductions)) "pca" else NULL
  if (is.null(red)) {
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = max(CFG$dims), verbose = FALSE)
    red <- "pca"
  }
  obj <- force_legacy_assay_for_hdwgcna(obj, CFG$assay)
  list(obj = obj, reduction = red)
}

select_hdwgcna_power <- function(obj, label) {
  if (!is.null(CFG$hdwgcna_soft_power)) return(CFG$hdwgcna_soft_power)
  ptab <- tryCatch(hdWGCNA::GetPowerTable(obj), error = function(e) data.frame())
  if (nrow(ptab) > 0) {
    save_table(ptab, paste0("16_hdwgcna_power_table_", safe_file_name(label), ".csv"))
    power_col <- intersect(c("Power", "power"), colnames(ptab))[1]
    r2_col <- intersect(c("SFT.R.sq", "SFT.R.sq.", "Rsquared", "R.sq", "sft_r2"), colnames(ptab))[1]
    if (!is.na(power_col) && !is.na(r2_col)) {
      ok <- which(is.finite(ptab[[r2_col]]) & ptab[[r2_col]] >= CFG$hdwgcna_soft_power_r2)
      if (length(ok) > 0) return(ptab[[power_col]][ok[1]])
      best <- which.max(ptab[[r2_col]])
      if (length(best) == 1 && is.finite(ptab[[power_col]][best])) return(ptab[[power_col]][best])
    }
  }
  CFG$hdwgcna_default_soft_power
}

prepare_hdwgcna_object <- function(obj, ct, mode = c("full", "subset")) {
  mode <- match.arg(mode)
  sw <- obj
  if (mode == "subset") {
    sw <- subset(sw, subset = celltype == ct)
    if (ncol(sw) < CFG$min_cells_per_celltype_hdwgcna) stop("Too few cells after subset for hdWGCNA: ", ct)
  }
  patch_getassay_slot_compat()
  sw <- join_layers_safe(sw)
  sw <- force_legacy_assay_for_hdwgcna(sw, CFG$assay)
  sw <- safe_set_default_assay(sw, CFG$assay)
  rr <- choose_hdwgcna_reduction(sw)
  sw <- rr$obj; red <- rr$reduction
  
  wname <- paste0("wgcna_", safe_file_name(ct), "_", mode)
  sw <- call_pkg_fun("hdWGCNA", "SetupForWGCNA", list(
    seurat_obj = sw,
    gene_select = CFG$hdwgcna_gene_select,
    fraction = CFG$hdwgcna_fraction,
    wgcna_name = wname
  ))
  
  metacell_args <- list(
    seurat_obj = sw,
    group.by = c("celltype", "sample"),
    reduction = red,
    k = CFG$hdwgcna_metacell_k,
    max_shared = CFG$hdwgcna_metacell_max_shared,
    ident.group = "celltype",
    assay = CFG$assay,
    min_cells = min(CFG$hdwgcna_metacell_k, CFG$hdwgcna_min_cells_per_sample_ct)
  )
  sw <- call_pkg_fun("hdWGCNA", "MetacellsByGroups", metacell_args)
  sw <- call_pkg_fun("hdWGCNA", "NormalizeMetacells", list(seurat_obj = sw, assay = CFG$assay, slot = "data", layer = "data"))
  
  # Official target-celltype extraction. In subset mode this still works because all
  # metacells belong to the selected cell type.
  sw <- call_pkg_fun("hdWGCNA", "SetDatExpr", list(
    seurat_obj = sw,
    group_name = ct,
    group.by = "celltype",
    assay = CFG$assay,
    slot = "data",
    layer = "data",
    wgcna_name = wname
  ))
  attr(sw, "hdwgcna_target_celltype") <- ct
  attr(sw, "hdwgcna_mode") <- mode
  attr(sw, "hdwgcna_reduction") <- red
  attr(sw, "hdwgcna_wgcna_name") <- wname
  sw
}

run_hdwgcna_build_network <- function(sw, label) {
  sw <- call_pkg_fun("hdWGCNA", "TestSoftPowers", list(
    seurat_obj = sw,
    networkType = CFG$hdwgcna_network_type,
    wgcna_name = attr(sw, "hdwgcna_wgcna_name")
  ))
  
  sp_plot <- tryCatch(hdWGCNA::PlotSoftPowers(sw), error = function(e) NULL)
  if (!is.null(sp_plot)) {
    if (is.list(sp_plot) && !inherits(sp_plot, c("gg", "ggplot", "patchwork"))) sp_plot <- patchwork::wrap_plots(sp_plot, ncol = 2)
    save_pdf(sp_plot, paste0("41_hdwgcna_soft_power_", safe_file_name(label), ".pdf"), 10, 7)
  }
  
  power <- select_hdwgcna_power(sw, label)
  log_msg("hdWGCNA selected soft power for ", label, ": ", power)
  
  # Current hdWGCNA uses soft_power; some older snippets used power. Keep both in
  # the candidate list and only pass supported formal arguments.
  sw <- call_pkg_fun("hdWGCNA", "ConstructNetwork", list(
    seurat_obj = sw,
    tom_name = safe_file_name(label),
    soft_power = power,
    power = power,
    overwrite_tom = TRUE,
    networkType = CFG$hdwgcna_network_type,
    minModuleSize = CFG$hdwgcna_min_module_size,
    mergeCutHeight = CFG$hdwgcna_merge_cut_height,
    setDatExpr = FALSE,
    wgcna_name = attr(sw, "hdwgcna_wgcna_name")
  ))
  
  sw <- call_pkg_fun("hdWGCNA", "ModuleEigengenes", list(
    seurat_obj = sw,
    group.by.vars = intersect(c("sample", "group"), colnames(sw@meta.data)),
    wgcna_name = attr(sw, "hdwgcna_wgcna_name")
  ))
  sw <- call_pkg_fun("hdWGCNA", "ModuleConnectivity", list(
    seurat_obj = sw,
    wgcna_name = attr(sw, "hdwgcna_wgcna_name")
  ))
  sw <- tryCatch(call_pkg_fun("hdWGCNA", "ResetModuleNames", list(
    seurat_obj = sw,
    new_name = paste0(safe_file_name(label), "_M"),
    wgcna_name = attr(sw, "hdwgcna_wgcna_name")
  )), error = function(e) sw)
  sw
}

plot_export_hdwgcna <- function(sw, label) {
  mods <- safe_get_modules(sw)
  if (nrow(mods) > 0) save_table(mods, paste0("16_hdwgcna_modules_", safe_file_name(label), ".csv"))
  hubs <- safe_get_hubs(sw, n_hubs = 30)
  if (nrow(hubs) > 0) save_table(hubs, paste0("16_hdwgcna_hubgenes_", safe_file_name(label), ".csv"))
  
  # Official hdWGCNA plots plus readable summaries. Dendrogram is kept but de-emphasized;
  # publication figures should usually prioritize module feature plots, module eigengenes,
  # hub-gene networks, and module-trait correlations.
  safe_plot_hdwgcna(hdWGCNA::PlotDendrogram(sw, main = paste0("hdWGCNA module dendrogram: ", label)), paste0("42_hdwgcna_official_dendrogram_", safe_file_name(label), ".pdf"), 10, 6)
  safe_plot_hdwgcna(hdWGCNA::PlotKMEs(sw, ncol = 4, text_size = 2.5), paste0("43_hdwgcna_kME_hubgenes_", safe_file_name(label), ".pdf"), 12, 9)
  safe_plot_hdwgcna(call_pkg_fun("hdWGCNA", "ModuleFeaturePlot", list(seurat_obj = sw, features = "hMEs", order = TRUE, restrict_range = FALSE, ncol = CFG$hdwgcna_module_feature_ncol, wgcna_name = attr(sw, "hdwgcna_wgcna_name"))), paste0("44_hdwgcna_module_featureplot_", safe_file_name(label), ".pdf"), 12, 10)
  safe_plot_hdwgcna(hdWGCNA::ModuleCorrelogram(sw), paste0("45_hdwgcna_module_correlogram_", safe_file_name(label), ".pdf"), 8, 8)
  plot_module_size_bar(mods, label)
  plot_hub_gene_dotplot(hubs, label)
  plot_me_by_group(sw, label)
  if (isTRUE(CFG$hdwgcna_publication_plots)) plot_hdwgcna_official_networks(sw, label, mods)
  
  # Module-trait correlation using official hdWGCNA functions first.
  sw <- safe_run(paste0("hdWGCNA_ModuleTraitCorrelation_", label), {
    if ("group" %in% colnames(sw@meta.data)) sw$group_fac <- factor(sw$group)
    traits <- intersect(c("group_fac", "nFeature_RNA", "nCount_RNA", "percent.mt"), colnames(sw@meta.data))
    if (length(traits) == 0) stop("No usable traits for ModuleTraitCorrelation")
    sw2 <- call_pkg_fun("hdWGCNA", "ModuleTraitCorrelation", list(
      seurat_obj = sw,
      traits = traits,
      group.by = "celltype",
      wgcna_name = attr(sw, "hdwgcna_wgcna_name")
    ))
    safe_plot_hdwgcna(hdWGCNA::PlotModuleTraitCorrelation(sw2, label = "fdr", text_size = 2, high_color = "#EC7063", mid_color = "white", low_color = "#5DADE2"), paste0("46_hdwgcna_module_trait_corr_", safe_file_name(label), ".pdf"), 9, 10)
    sw2
  }, default = sw)
  
  # Manual fallback correlation if official plotting/table is not available.
  if (HAS("pheatmap")) {
    mes <- tryCatch(hdWGCNA::GetMEs(sw), error = function(e) NULL)
    if (!is.null(mes) && nrow(mes) > 2) {
      md <- sw@meta.data[rownames(mes), , drop = FALSE]
      traits <- data.frame(row.names = rownames(md))
      if ("group" %in% colnames(md)) traits$group_numeric <- as.numeric(factor(md$group))
      if ("nFeature_RNA" %in% colnames(md)) traits$nFeature_RNA <- md$nFeature_RNA
      if ("nCount_RNA" %in% colnames(md)) traits$nCount_RNA <- md$nCount_RNA
      if ("percent.mt" %in% colnames(md)) traits$percent_mt <- md$percent.mt
      if (ncol(traits) > 0) {
        cors <- suppressWarnings(cor(mes, traits, use = "pairwise.complete.obs", method = "spearman"))
        save_table(data.frame(module = rownames(cors), cors, check.names = FALSE), paste0("16_hdwgcna_manual_module_trait_corr_", safe_file_name(label), ".csv"))
        save_pdf_expr(paste0("47_hdwgcna_manual_module_trait_heatmap_", safe_file_name(label), ".pdf"), {
          pheatmap::pheatmap(cors, color = macaron_grad, main = paste0("Manual module-trait Spearman: ", label))
        }, 8, 8)
      }
    }
  }
  sw
}

plain_wgcna_metacell_rescue <- function(obj, ct) {
  # Robust rescue layer for SeuratObject v5 + hdWGCNA/WGCNA version conflicts.
  # It does NOT call WGCNA::blockwiseModules, because some R/WGCNA combinations fail
  # with unused arguments such as weights.x / cosine. Instead it builds a compact
  # metacell correlation network with base R and exports publication-style summaries.
  if (!isTRUE(CFG$hdwgcna_rescue_plain_wgcna)) return(NULL)
  label <- paste0(safe_file_name(ct), "_rescue")
  sub <- subset(obj, subset = celltype == ct)
  if (ncol(sub) < CFG$min_cells_per_celltype_hdwgcna) stop("Too few cells for robust WGCNA rescue: ", ct)
  sub <- join_layers_safe(sub)
  sub <- safe_set_default_assay(sub, CFG$assay)
  
  data_mat <- get_assay_layer(sub, assay = CFG$assay, layer = "data")
  count_mat <- get_assay_layer(sub, assay = CFG$assay, layer = "counts")
  vars <- intersect(VariableFeatures(obj), rownames(data_mat))
  if (length(vars) < 500) {
    det <- Matrix::rowSums(count_mat > 0)
    vars <- names(sort(det, decreasing = TRUE))[seq_len(min(3000, length(det)))]
  }
  vars <- head(intersect(vars, rownames(data_mat)), 2500)
  if (length(vars) < 200) stop("Too few genes for robust WGCNA rescue: ", ct)
  data_mat <- data_mat[vars, , drop = FALSE]
  meta <- sub@meta.data
  
  # Build sample/group metacells; if coverage is too small, make sequential bins.
  groups <- interaction(meta$sample, meta$group, drop = TRUE, sep = "__")
  names(groups) <- colnames(data_mat)
  group_counts <- table(groups)
  keep_groups <- names(group_counts)[group_counts >= 20]
  if (length(keep_groups) >= 4) {
    keep_cells <- names(groups)[groups %in% keep_groups]
    data_mat <- data_mat[, keep_cells, drop = FALSE]
    groups <- droplevels(groups[keep_cells])
  } else {
    set.seed(CFG$seed)
    n_bins <- min(80, max(6, floor(ncol(data_mat) / 25)))
    groups <- factor(paste0("metacell_", cut(seq_len(ncol(data_mat)), breaks = n_bins, labels = FALSE)))
    names(groups) <- colnames(data_mat)
  }
  
  metacell_expr <- t(sapply(levels(groups), function(g) Matrix::rowMeans(data_mat[, groups == g, drop = FALSE])))
  metacell_expr <- as.data.frame(metacell_expr, check.names = FALSE)
  # Drop unstable genes and metacells.
  metacell_expr <- metacell_expr[, colSums(is.finite(as.matrix(metacell_expr))) == nrow(metacell_expr), drop = FALSE]
  sds <- apply(metacell_expr, 2, stats::sd, na.rm = TRUE)
  metacell_expr <- metacell_expr[, is.finite(sds) & sds > 0, drop = FALSE]
  if (ncol(metacell_expr) < 100 || nrow(metacell_expr) < 4) stop("Metacell expression matrix too small for robust WGCNA rescue: ", ct)
  
  # Use base correlation only to avoid WGCNA::blockwiseModules compatibility errors.
  zdat <- scale(as.matrix(metacell_expr))
  zdat[!is.finite(zdat)] <- 0
  gene_cor <- suppressWarnings(stats::cor(zdat, use = "pairwise.complete.obs", method = "pearson"))
  gene_cor[!is.finite(gene_cor)] <- 0
  diag(gene_cor) <- 1
  
  # Pick a conservative power from WGCNA::pickSoftThreshold when possible, otherwise default.
  power <- CFG$hdwgcna_soft_power %||% CFG$hdwgcna_default_soft_power
  fit <- data.frame()
  if (HAS("WGCNA")) {
    fit <- tryCatch({
      suppressPackageStartupMessages(library(WGCNA))
      sft <- WGCNA::pickSoftThreshold(metacell_expr, powerVector = c(1:10, seq(12, 20, 2)), networkType = CFG$hdwgcna_network_type, verbose = 0)
      sft$fitIndices
    }, error = function(e) data.frame())
    if (nrow(fit) > 0) {
      save_table(fit, paste0("16_hdwgcna_rescue_power_table_", label, ".csv"))
      ok <- which(is.finite(fit$SFT.R.sq) & fit$SFT.R.sq >= CFG$hdwgcna_soft_power_r2)
      power <- if (length(ok) > 0) fit$Power[ok[1]] else fit$Power[which.max(fit$SFT.R.sq)]
      if (!is.finite(power) || length(power) == 0) power <- CFG$hdwgcna_default_soft_power
    }
  }
  
  adjacency <- abs(gene_cor)^power
  diag(adjacency) <- 1
  # Cluster genes by correlation distance. This is a stable rescue substitute for TOM clustering.
  dist_gene <- as.dist(1 - adjacency)
  hc <- stats::hclust(dist_gene, method = "average")
  k <- min(12, max(4, round(ncol(metacell_expr) / 350)))
  mod_id <- stats::cutree(hc, k = k)
  module_names <- paste0("R", sprintf("%02d", mod_id))
  small <- names(which(table(module_names) < CFG$hdwgcna_min_module_size))
  module_names[module_names %in% small] <- "grey"
  modules <- data.frame(gene = colnames(metacell_expr), module = module_names, stringsAsFactors = FALSE)
  save_table(modules, paste0("16_hdwgcna_rescue_modules_", label, ".csv"))
  
  # Module eigengenes as first PC per module, with average fallback.
  module_levels <- setdiff(sort(unique(modules$module)), c("grey", "gray"))
  ME <- sapply(module_levels, function(m) {
    genes <- modules$gene[modules$module == m]
    if (length(genes) == 1) return(as.numeric(scale(metacell_expr[[genes]])))
    x <- scale(as.matrix(metacell_expr[, genes, drop = FALSE])); x[!is.finite(x)] <- 0
    pc <- tryCatch(stats::prcomp(x, center = FALSE, scale. = FALSE)$x[, 1], error = function(e) rowMeans(x))
    as.numeric(scale(pc))
  })
  if (is.null(dim(ME))) ME <- matrix(ME, ncol = 1, dimnames = list(rownames(metacell_expr), module_levels[1]))
  ME <- as.data.frame(ME, check.names = FALSE)
  rownames(ME) <- rownames(metacell_expr)
  save_table(data.frame(metacell = rownames(ME), ME, check.names = FALSE), paste0("16_hdwgcna_rescue_module_eigengenes_", label, ".csv"))
  
  hubs <- do.call(rbind, lapply(module_levels, function(m) {
    genes <- modules$gene[modules$module == m]
    if (length(genes) == 0) return(NULL)
    kme <- suppressWarnings(stats::cor(as.matrix(metacell_expr[, genes, drop = FALSE]), ME[[m]], use = "pairwise.complete.obs"))
    data.frame(module = m, gene = genes, kME = as.numeric(kme), stringsAsFactors = FALSE)
  }))
  if (!is.null(hubs) && nrow(hubs) > 0) {
    hubs <- hubs[order(hubs$module, -abs(hubs$kME)), , drop = FALSE]
    save_table(hubs, paste0("16_hdwgcna_rescue_hubgenes_", label, ".csv"))
  }
  
  # Clean, publication-style rescue figures.
  plot_module_size_bar(modules, label, prefix = "48_hdwgcna_rescue_module_size")
  plot_hub_gene_dotplot(hubs, label, prefix = "49_hdwgcna_rescue_top_hub_genes")
  
  # Module eigengene correlation heatmap.
  if (HAS("pheatmap") && ncol(ME) >= 2) {
    save_pdf_expr(paste0("50_hdwgcna_rescue_ME_correlation_heatmap_", label, ".pdf"), {
      pheatmap::pheatmap(stats::cor(ME, use = "pairwise.complete.obs"), color = macaron_grad,
                         main = paste0("Module eigengene correlation: ", ct), fontsize_row = 9, fontsize_col = 9)
    }, 8, 8)
  }
  
  # Module eigengenes by group/sample if sample/group metacells exist.
  mc_meta <- data.frame(metacell = rownames(metacell_expr), stringsAsFactors = FALSE)
  parts <- strsplit(mc_meta$metacell, "__", fixed = TRUE)
  mc_meta$sample <- vapply(parts, function(z) z[1] %||% NA_character_, character(1))
  mc_meta$group <- vapply(parts, function(z) z[2] %||% NA_character_, character(1))
  if (any(!is.na(mc_meta$group)) && ncol(ME) > 0) {
    df <- cbind(mc_meta, ME)
    top_me <- colnames(ME)[seq_len(min(ncol(ME), CFG$hdwgcna_max_modules_to_plot))]
    long <- tidyr::pivot_longer(df, cols = all_of(top_me), names_to = "module", values_to = "eigengene")
    p <- ggplot(long, aes(group, eigengene, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.75, linewidth = 0.3) +
      geom_jitter(width = 0.13, size = 1.1, alpha = 0.65) +
      facet_wrap(~module, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = macaron_pal(length(unique(long$group)))) +
      theme_bw(base_size = 11) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1), strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", hjust = 0.5)) +
      labs(title = paste0("Rescue module eigengenes by group: ", ct), x = NULL, y = "Module eigengene")
    save_pdf(p, paste0("51_hdwgcna_rescue_module_eigengene_by_group_", label, ".pdf"), 11, 8)
  }
  
  # Gene-network map using MDS of correlation distance for top hub genes only.
  top_hubs <- if (!is.null(hubs) && nrow(hubs) > 0) {
    do.call(rbind, lapply(split(hubs, hubs$module), function(z) head(z[order(abs(z$kME), decreasing = TRUE), , drop = FALSE], 25)))
  } else data.frame()
  if (nrow(top_hubs) >= 20) {
    genes_plot <- unique(top_hubs$gene)
    d2 <- as.dist(1 - abs(gene_cor[genes_plot, genes_plot, drop = FALSE]))
    xy <- tryCatch(as.data.frame(stats::cmdscale(d2, k = 2)), error = function(e) data.frame())
    if (nrow(xy) == length(genes_plot)) {
      colnames(xy) <- c("Dim1", "Dim2")
      xy$gene <- genes_plot
      xy <- merge(xy, top_hubs[, c("gene", "module", "kME")], by = "gene", all.x = TRUE)
      p <- ggplot(xy, aes(Dim1, Dim2, color = module, size = abs(kME))) +
        geom_point(alpha = 0.9) +
        ggrepel::geom_text_repel(aes(label = gene), size = 2.4, max.overlaps = 60, show.legend = FALSE) +
        scale_color_manual(values = macaron_pal(length(unique(xy$module)))) +
        theme_void(base_size = 12) +
        theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5)) +
        labs(title = paste0("Rescue hub-gene network map: ", ct), color = "Module", size = "|kME|")
      save_pdf(p, paste0("52_hdwgcna_rescue_hub_gene_network_map_", label, ".pdf"), 10, 8)
    }
  }
  
  attr(modules, "hdwgcna_mode") <- "robust_rescue_base_correlation"
  attr(modules, "hdwgcna_reduction") <- "metacell_base_correlation"
  attr(modules, "hdwgcna_power") <- power
  modules
}

run_hdwgcna_one <- function(obj, ct) {
  if (!HAS("hdWGCNA") || !HAS("WGCNA")) stop("hdWGCNA or WGCNA missing")
  suppressPackageStartupMessages({ library(hdWGCNA); library(WGCNA) })
  label <- safe_file_name(ct)
  
  # Primary: official-style full object metacells + SetDatExpr target CT.
  primary <- NULL
  if (isTRUE(CFG$hdwgcna_try_full_object_metacells)) {
    primary <- safe_run(paste0("hdWGCNA_primary_full_metacells_", label), {
      sw <- prepare_hdwgcna_object(obj, ct, mode = "full")
      sw <- run_hdwgcna_build_network(sw, paste0(label, "_full"))
      plot_export_hdwgcna(sw, paste0(label, "_full"))
    }, default = NULL)
  }
  if (!is.null(primary) && inherits(primary, "Seurat")) return(primary)
  
  # Fallback: subset target cell type before metacells. This is less global, but often
  # rescues GEO datasets where some non-target cell types or samples break metacell construction.
  if (isTRUE(CFG$hdwgcna_try_subset_fallback)) {
    fallback <- safe_run(paste0("hdWGCNA_fallback_subset_", label), {
      sw <- prepare_hdwgcna_object(obj, ct, mode = "subset")
      sw <- run_hdwgcna_build_network(sw, paste0(label, "_subset"))
      plot_export_hdwgcna(sw, paste0(label, "_subset"))
    }, default = NULL)
    if (!is.null(fallback) && inherits(fallback, "Seurat")) return(fallback)
  }
  
  rescue <- safe_run(paste0("hdWGCNA_plain_WGCNA_rescue_", label), {
    plain_wgcna_metacell_rescue(obj, ct)
  }, default = NULL)
  if (!is.null(rescue)) return(rescue)
  
  stop("Both hdWGCNA strategies failed for cell type: ", ct)
}

run_hdwgcna_all <- function(obj) {
  if (is.null(obj) || !inherits(obj, "Seurat")) stop("Input object is not a Seurat object")
  targets <- choose_hdwgcna_celltypes(obj)
  log_msg("hdWGCNA target cell types: ", paste(targets, collapse = ", "))
  res <- list()
  for (ct in targets) {
    res[[ct]] <- cache_run(paste0("16_hdwgcna_", safe_file_name(ct), "_v7"), run_hdwgcna_one(obj, ct), default = NULL)
  }
  res <- res[!vapply(res, is.null, logical(1))]
  
  # Summary table across successful networks.
  if (length(res) > 0) {
    summary <- do.call(rbind, lapply(names(res), function(ct) {
      mods <- if (inherits(res[[ct]], "Seurat")) safe_get_modules(res[[ct]]) else if (is.data.frame(res[[ct]]) && all(c("gene", "module") %in% colnames(res[[ct]]))) res[[ct]] else data.frame()
      data.frame(
        celltype = ct,
        mode = attr(res[[ct]], "hdwgcna_mode") %||% attr(res[[ct]], "hdwgcna_mode") %||% NA,
        reduction = attr(res[[ct]], "hdwgcna_reduction") %||% NA,
        modules = if (nrow(mods) > 0) length(setdiff(unique(mods$module), "grey")) else NA_integer_,
        genes_in_modules = if (nrow(mods) > 0) sum(mods$module != "grey", na.rm = TRUE) else NA_integer_,
        stringsAsFactors = FALSE
      )
    }))
    save_table(summary, "16_hdwgcna_network_summary.csv")
  }
  res
}
hdwgcna_results <- if (isTRUE(CFG$run_hdwgcna)) cache_run("16_hdwgcna_results_v7", run_hdwgcna_all(combined), default = list()) else list()

# ============================= FINAL EXPORT ==================================
export_final <- function(obj) {
  if (is.null(obj) || !inherits(obj, "Seurat")) stop("Input object is not a Seurat object")
  obj <- join_layers_safe(obj)
  pbmc <- obj
  cellAnn <- obj$celltype
  clusterAnn <- as.data.frame(table(seurat_clusters = obj$seurat_clusters, celltype = obj$celltype), stringsAsFactors = FALSE)
  clusterAnn <- clusterAnn[clusterAnn$Freq > 0, , drop = FALSE]
  clusterAnn <- clusterAnn[order(clusterAnn$seurat_clusters, -clusterAnn$Freq), , drop = FALSE]
  clusterAnn <- clusterAnn[!duplicated(clusterAnn$seurat_clusters), , drop = FALSE]
  final_rds <- file.path(out_dir, "Seurat_final.rds")
  final_rdata <- file.path(out_dir, "Seurat_final.RData")
  saveRDS(pbmc, final_rds)
  save(pbmc, cellAnn, clusterAnn, all_markers, de_results, pseudobulk_results, enrichment_results, composition_results, pseudotime_results, cellchat_results, hdwgcna_results, file = final_rdata)
  log_msg("Final Seurat RDS saved: ", final_rds)
  log_msg("Final RData saved: ", final_rdata)
  
  # Generate simple HTML-like text report.
  figs <- list.files(fig_dir, pattern = "\\.pdf$", full.names = FALSE)
  tabs <- list.files(tab_dir, pattern = "\\.csv$", full.names = FALSE)
  report <- file.path(report_dir, "analysis_summary.txt")
  cat("GEO scRNA-seq advanced pipeline summary\n", file = report)
  cat("======================================\n\n", file = report, append = TRUE)
  cat("GSE: ", CFG$GSE, "\n", sep = "", file = report, append = TRUE)
  cat("Species: ", CFG$species, "\n", sep = "", file = report, append = TRUE)
  cat("Cells after QC: ", ncol(obj), "\n", sep = "", file = report, append = TRUE)
  cat("Genes: ", nrow(obj), "\n", sep = "", file = report, append = TRUE)
  cat("Samples: ", paste(unique(obj$sample), collapse = ", "), "\n", sep = "", file = report, append = TRUE)
  cat("Groups: ", paste(unique(obj$group), collapse = ", "), "\n", sep = "", file = report, append = TRUE)
  cat("Cell types: ", paste(unique(obj$celltype), collapse = ", "), "\n\n", sep = "", file = report, append = TRUE)
  cat("Module status file: ", status_file, "\n", sep = "", file = report, append = TRUE)
  cat("Figures generated: ", length(figs), "\n", sep = "", file = report, append = TRUE)
  cat(paste0(" - ", figs, collapse = "\n"), "\n\n", file = report, append = TRUE)
  cat("Tables generated: ", length(tabs), "\n", sep = "", file = report, append = TRUE)
  cat(paste0(" - ", tabs, collapse = "\n"), "\n", file = report, append = TRUE)
  log_msg("Summary report saved: ", report)
  TRUE
}

safe_run("final_export", export_final(combined), default = NULL, critical = FALSE)

log_msg("Pipeline completed.")
log_msg("Results directory: ", out_dir)
log_msg("PDF figures: ", fig_dir)
log_msg("Tables: ", tab_dir)
log_msg("Checkpoints: ", rdata_dir)
log_msg("Logs: ", log_dir)
###############################################################################
