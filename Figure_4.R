library(Seurat)
library(harmony)
library(slingshot)
library(ggplot2)
library(ggrastr)
library(ggpubr)
library(viridis)
library(MetBrewer)      
library(WGCNA)
library(hdWGCNA)
library(doParallel)
library(SingleCellExperiment)
library(scuttle)
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggrepel)
library(ggtext)
library(ComplexHeatmap)
library(doParallel)
library(enrichR)
library(GeneOverlap)
library(tidyverse)
library(cowplot)
library(patchwork)

big_text_theme <- theme(
  axis.text = element_text(size = 14),
  axis.title = element_text(size = 16, face = "bold"),
  strip.text = element_text(size = 16, face = "bold"),
  legend.text = element_text(size = 13),
  legend.title = element_text(size = 14, face = "bold"),
  plot.title = element_text(size = 18, face = "bold", margin = margin(b = 12))
)

# -----------------------------------------------------------------------
# Load PER-CELL-LINE objects produced in Figure_4.R. IMPORTANT: these are
# the FRESHLY-RECOMPUTED per-cell-line UMAP embeddings (built and saved at
# the end of the Figure 4a loop in Figure_4.R), NOT the earlier pooled-
# embedding-then-subset objects. Slingshot trajectory inference below
# depends directly on UMAP coordinates, so using the same per-line
# embeddings that Figure 4a itself visualizes keeps Figure 4 and Figure 5
# consistent - each cell line gets its own independent embedding and
# downstream trajectory/WGCNA analysis, matching Figure 4's approach.
# -----------------------------------------------------------------------
Cortical_lineage_list <- readRDS("~/project/IPSC_2025_Data/merged_IPSC_derived_pallial_lineages_by_line_umap")
Hem_lineage_list <- readRDS("~/project/IPSC_2025_Data/merged_IPSC_derived_hem_lineages_by_line_umap")
cell_lines <- names(Cortical_lineage_list)

# Result containers, keyed by cell line
curves_cort_list <- list()
curves_hem_list  <- list()
p_pseudotime_cort_list <- list()
p_pseudotime_hem_list  <- list()

for (cl in cell_lines) {
  
  Cortical_lineage <- Cortical_lineage_list[[cl]]
  Hem_lineage <- Hem_lineage_list[[cl]]
  
  ENDS <- c("DL_ExN", "UL_ExN", "A1 Astrocyte", "A2 Astrocyte")
  set.seed(1)
  Cortical_lineages <- as.SlingshotDataSet(getLineages(
    data           = Cortical_lineage@reductions$umap@cell.embeddings,
    clusterLabels  = Cortical_lineage$Celltype2,
    dist.method    = "slingshot",
    end.clus       = ENDS,
    start.clus     = "RG"))
  
  ENDS <- c("CRN", "Epithelial")
  Hem_lineages <- as.SlingshotDataSet(getLineages(
    data           = Hem_lineage@reductions$umap@cell.embeddings,
    clusterLabels  = Hem_lineage$Celltype2,
    dist.method    = "slingshot",
    end.clus       = ENDS,
    start.clus     = "Hem_RG"))
  
  curves_cort <- as.SlingshotDataSet(getCurves(
    data          = Cortical_lineages,
    thresh        = 1e-1,
    stretch       = 1e-1,
    allow.breaks  = F,
    approx_points = 150
  ))
  
  curves_hem <- as.SlingshotDataSet(getCurves(
    data          = Hem_lineages,
    thresh        = 1e-1,
    stretch       = 1e-1,
    allow.breaks  = F,
    approx_points = 150
  ))
  
  curves_cort_list[[cl]] <- curves_cort
  curves_hem_list[[cl]]  <- curves_hem
  
  pseudotime_cort <- as.data.frame(slingPseudotime(curves_cort, na = FALSE))
  pseudotime_cort$Pseudotime <- apply(pseudotime_cort, 1, max)
  
  pseudotime_hem <- as.data.frame(slingPseudotime(curves_hem, na = FALSE))
  pseudotime_hem$Pseudotime <- apply(pseudotime_hem, 1, max)
  
  Cortical_lineage$pseudotime <- pseudotime_cort$Pseudotime
  Cortical_lineage$dp_pseudotime <- ifelse(Cortical_lineage$Celltype2 %in% c("RG", "IPC_ExN", "DL_ExN"), pseudotime_cort$Lineage2, NA)
  Cortical_lineage$up_pseudotime <- ifelse(Cortical_lineage$Celltype2 %in% c("RG", "IPC_ExN", 'UL_ExN'), pseudotime_cort$Lineage4 , NA)
  Cortical_lineage$A1_pseudotime <- ifelse(Cortical_lineage$Celltype2 %in% c("RG", "IPC_ExN", 'A1 Astrocyte'), pseudotime_cort$Lineage3, NA)
  Cortical_lineage$A2_pseudotime <- ifelse(Cortical_lineage$Celltype2 %in% c("RG", 'A2 Astrocyte'), pseudotime_cort$Lineage1, NA)
  
  Hem_lineage$pseudotime <- pseudotime_hem$Pseudotime
  Hem_lineage$crn_pseudotime <- ifelse(Hem_lineage$Celltype %in% c("Hem_RG", "CRN"), pseudotime_hem$Lineage1 , NA)
  Hem_lineage$epi_pseudotime <- ifelse(Hem_lineage$Celltype %in% c("Hem_RG", "Epithelial"), pseudotime_hem$Lineage2 , NA)
  
  Cortical_lineage$UMAP1 <- Cortical_lineage@reductions$umap@cell.embeddings[,1]
  Cortical_lineage$UMAP2 <- Cortical_lineage@reductions$umap@cell.embeddings[,2]
  Hem_lineage$UMAP1 <- Hem_lineage@reductions$umap@cell.embeddings[,1]
  Hem_lineage$UMAP2 <- Hem_lineage@reductions$umap@cell.embeddings[,2]
  
  # Figure 5a panels for this cell line
  p1 <- Cortical_lineage@meta.data %>%
    ggplot(aes(x=UMAP1, y=UMAP2, color=dp_pseudotime)) +
    ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
    coord_equal() + scale_color_gradientn(colors=plasma(256), na.value='grey') +
    umap_theme() + big_text_theme + ggtitle(paste0(cl, ": DP"))
  p2 <- Cortical_lineage@meta.data %>%
    ggplot(aes(x=UMAP1, y=UMAP2, color=up_pseudotime)) +
    ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
    coord_equal() + scale_color_gradientn(colors=viridis(256), na.value='grey') +
    umap_theme() + big_text_theme + ggtitle(paste0(cl, ": UP"))
  p3 <- Cortical_lineage@meta.data %>%
    ggplot(aes(x=UMAP1, y=UMAP2, color=A1_pseudotime)) +
    ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
    coord_equal() + scale_color_gradientn(colors=inferno(256), na.value='grey') +
    umap_theme() + big_text_theme + ggtitle(paste0(cl, ": A1"))
  p4 <- Cortical_lineage@meta.data %>%
    ggplot(aes(x=UMAP1, y=UMAP2, color=A2_pseudotime)) +
    ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
    coord_equal() + scale_color_gradientn(colors=mako(256), na.value='grey') +
    umap_theme() + big_text_theme + ggtitle(paste0(cl, ": A2"))
  p5 <- Hem_lineage@meta.data %>%
    ggplot(aes(x=UMAP1, y=UMAP2, color=crn_pseudotime)) +
    ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
    coord_equal() + scale_color_gradientn(colors=inferno(256), na.value='grey') +
    umap_theme() + big_text_theme + ggtitle(paste0(cl, ": CRN"))
  p6 <- Hem_lineage@meta.data %>%
    ggplot(aes(x=UMAP1, y=UMAP2, color=epi_pseudotime)) +
    ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
    coord_equal() + scale_color_gradientn(colors=mako(256), na.value='grey') +
    umap_theme() + big_text_theme + ggtitle(paste0(cl, ": Epi"))
  
  p_pseudotime_cort_list[[cl]] <- (p1 + p2) / (p3 + p4)
  p_pseudotime_hem_list[[cl]]  <- (p5 + p6)
  
  Cortical_lineage_list[[cl]] <- Cortical_lineage
  Hem_lineage_list[[cl]] <- Hem_lineage
}

# -----------------------------------------------------------------------
# Figure 5a: ONE combined figure, one row per cell line
# -----------------------------------------------------------------------
fig4a_cort_combined <- wrap_plots(p_pseudotime_cort_list, ncol = 1) +
  plot_annotation(title = "Cortical Lineage By Slingshot Pseudotime, per Cell Line")
ggsave("~/project/IPSC_2025_Data/Figure5a_cortical_pseudotime.tiff",
       plot = fig4a_cort_combined, device = "tiff",
       width = 12, height = 10 * length(cell_lines), dpi = 300, limitsize = FALSE)

fig4a_hem_combined <- wrap_plots(p_pseudotime_hem_list, ncol = 1) +
  plot_annotation(title = "Hem Lineage By Slingshot Pseudotime, per Cell Line")
ggsave("~/project/IPSC_2025_Data/Figure5a_hem_pseudotime.tiff",
       plot = fig4a_hem_combined, device = "tiff",
       width = 10, height = 5 * length(cell_lines), dpi = 300, limitsize = FALSE)


# -----------------------------------------------------------------------
# Custom kME plot: hub genes on the X-AXIS (not in a side list box like
# hdWGCNA's default PlotKMEs()), faceted by module. Manually implements a
# "reorder within facet" bar layout (so genes are sorted by kME within
# each module's panel) without requiring the tidytext package.
# -----------------------------------------------------------------------
plot_kme_custom <- function(obj, wgcna_name, n_hubs = 10, title = "") {
  modules_df <- GetModules(obj, wgcna_name = wgcna_name)
  mod_levels <- setdiff(levels(modules_df$module), "grey")
  
  # kME_<x> columns are keyed by each module's original WGCNA COLOR (e.g.
  # kME_turquoise), NOT by the module's current display name - renaming
  # (ResetModuleNames) only changes the `module` label column itself.
  # Look up each module's color from modules_df$color (present regardless
  # of renaming) to find its kME column; fall back to the display name in
  # case a given hdWGCNA version DOES key kME columns by name.
  has_color_col <- "color" %in% colnames(modules_df)
  
  hub_list <- lapply(mod_levels, function(m) {
    sub_df <- modules_df %>% filter(module == m)
    if (nrow(sub_df) == 0) return(NULL)
    
    kme_col <- NULL
    if (has_color_col) {
      this_color <- unique(sub_df$color)[1]
      candidate <- paste0("kME_", this_color)
      if (candidate %in% colnames(modules_df)) kme_col <- candidate
    }
    if (is.null(kme_col)) {
      candidate <- paste0("kME_", m)
      if (candidate %in% colnames(modules_df)) kme_col <- candidate
    }
    if (is.null(kme_col)) return(NULL)
    
    sub_df %>%
      arrange(desc(.data[[kme_col]])) %>%
      slice_head(n = n_hubs) %>%
      transmute(gene = gene_name, kME = .data[[kme_col]], module = m)
  })
  hub_df <- bind_rows(hub_list)
  if (nrow(hub_df) == 0) {
    warning("plot_kme_custom(): no hub genes found for wgcna_name='", wgcna_name,
            "'. Available GetModules() columns: ", paste(colnames(modules_df), collapse = ", "),
            " | module levels: ", paste(mod_levels, collapse = ", "),
            " -- Returning a blank placeholder plot.")
    return(ggplot() + theme_void() + ggtitle(paste0(title, " (no hub genes found)")))
  }
  
  # manual reorder-within-facet: combine gene+module into a unique factor
  # level, ordered by kME within each module, then strip the module suffix
  # back off for the displayed axis labels.
  hub_df <- hub_df %>%
    group_by(module) %>%
    mutate(gene_key = factor(paste(gene, module, sep = "___"),
                             levels = paste(gene[order(kME)], module, sep = "___"))) %>%
    ungroup()
  
  ggplot(hub_df, aes(x = gene_key, y = kME, fill = module)) +
    geom_col() +
    facet_wrap(~module, scales = "free_x", nrow = 2) +
    scale_x_discrete(labels = function(x) sub("___.*", "", x)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold")
    ) +
    labs(x = "Hub gene", y = "kME", title = title)
}

# -----------------------------------------------------------------------
# Figure 5b: WGCNA modules, run SEPARATELY within each cell line
#
# STAGE 1 (SLOW): module-finding. This is the expensive part - once it
# completes, a checkpoint is saved below and everything after it can be
# re-run independently without repeating this loop.
# -----------------------------------------------------------------------

for (cl in cell_lines) {
  
  Cortical_lineage <- Cortical_lineage_list[[cl]]
  Hem_lineage <- Hem_lineage_list[[cl]]
  
  cl_wgcna <- paste0("trajectory_", cl)
  
  # --- Cortical lineage WGCNA ---
  Cortical_lineage <- SetupForWGCNA(Cortical_lineage, gene_select = "fraction", fraction = 0.05, wgcna_name = cl_wgcna)
  Cortical_lineage <- MetacellsByGroups(
    Cortical_lineage, group.by = c("Celltype2", "SampleID"),
    reduction = 'harmony', k = 10, max_shared = 10, ident.group = 'Celltype2'
  )
  Cortical_lineage <- NormalizeMetacells(Cortical_lineage)
  cl_cores <- makeCluster(max(1, parallel::detectCores() - 1))
  registerDoParallel(cl_cores)
  Cortical_lineage <- SetDatExpr(
    Cortical_lineage,
    group_name = unique(Cortical_lineage@misc[[cl_wgcna]]$wgcna_metacell_obj$Celltype2),
    group.by='Celltype2', assay = 'RNA', layer = 'data'
  )
  Cortical_lineage <- TestSoftPowers(Cortical_lineage)
  Cortical_lineage <- ConstructNetwork(Cortical_lineage, tom_name = cl_wgcna, overwrite_tom = TRUE)
  Cortical_lineage <- ModuleEigengenes(Cortical_lineage, group.by.vars = c("SampleID"))
  Cortical_lineage <- ModuleConnectivity(Cortical_lineage)
  Cortical_lineage <- ResetModuleNames(Cortical_lineage, new_name = paste0(cl, "-Pallial-M"))
  stopCluster(cl_cores)
  
  # --- Hem lineage WGCNA ---
  cl_cores <- makeCluster(max(1, parallel::detectCores() - 1))
  registerDoParallel(cl_cores)
  Hem_lineage <- SetupForWGCNA(Hem_lineage, gene_select = "fraction", fraction = 0.05, wgcna_name = cl_wgcna)
  Hem_lineage <- MetacellsByGroups(
    Hem_lineage, group.by = c("Celltype2", "SampleID"),
    reduction = 'harmony', k = 10, max_shared = 10, ident.group = 'Celltype2'
  )
  Hem_lineage <- NormalizeMetacells(Hem_lineage)
  Hem_lineage <- SetDatExpr(
    Hem_lineage,
    group_name = unique(Hem_lineage@misc[[cl_wgcna]]$wgcna_metacell_obj$Celltype2),
    group.by='Celltype2', assay = 'RNA', layer = 'data'
  )
  Hem_lineage <- TestSoftPowers(Hem_lineage)
  Hem_lineage <- ConstructNetwork(Hem_lineage, tom_name = cl_wgcna, overwrite_tom = TRUE)
  Hem_lineage <- ModuleEigengenes(Hem_lineage, group.by.vars = c("SampleID"))
  Hem_lineage <- ModuleConnectivity(Hem_lineage)
  Hem_lineage <- ResetModuleNames(Hem_lineage, new_name = paste0(cl, "-Hem-M"))
  stopCluster(cl_cores)
  
  Cortical_lineage_list[[cl]] <- Cortical_lineage
  Hem_lineage_list[[cl]] <- Hem_lineage
}

# -----------------------------------------------------------------------
# CHECKPOINT: module-finding (the slow part - TestSoftPowers /
# ConstructNetwork / ModuleEigengenes / ModuleConnectivity) is complete
# at this point, with plain numeric module names (M1, M2, ...). Save here
# so any future re-run of renaming/relabeling/replotting below can start
# from this checkpoint instead of re-running the expensive construction.
# -----------------------------------------------------------------------
saveRDS(Cortical_lineage_list, "~/project/IPSC_2025_Data/checkpoint_pallial_modules_found_by_line")
saveRDS(Hem_lineage_list, "~/project/IPSC_2025_Data/checkpoint_hem_modules_found_by_line")

# =========================================================================
# EVERYTHING BELOW THIS POINT IS FAST (renaming + plotting only) and can be
# RE-RUN INDEPENDENTLY without repeating the module-finding above.
#
# If Cortical_lineage_list / Hem_lineage_list are not already in your R
# session (e.g. you're starting a fresh session today), uncomment these two
# lines to reload the checkpoint saved above instead of re-running the
# WGCNA construction loop:
#

print_top_go_terms <- function(obj, wgcna_name, database = "GO_Biological_Process_2023", n_terms = 8) {
  et <- GetEnrichrTable(obj, wgcna_name = wgcna_name)
  et %>%
    filter(db == database, module != "grey") %>%
    group_by(module) %>%
    arrange(P.value, .by_group = TRUE) %>%
    slice_head(n = n_terms) %>%
    select(module, Term, P.value, Adjusted.P.value, Genes) %>%
    ungroup()
}

for (cl in cell_lines) {
  cl_wgcna <- paste0("trajectory_", cl)
  
  cat("\n\n==========", cl, ": Pallial ==========\n")
  cort_top <- print_top_go_terms(Cortical_lineage_list[[cl]], cl_wgcna)
  print(cort_top, n = Inf)
  
  cat("\n\n==========", cl, ": Hem ==========\n")
  hem_top <- print_top_go_terms(Hem_lineage_list[[cl]], cl_wgcna)
  print(hem_top, n = Inf)
}

Cortical_lineage_list <- readRDS("~/project/IPSC_2025_Data/checkpoint_pallial_modules_found_by_line")
Hem_lineage_list <- readRDS("~/project/IPSC_2025_Data/checkpoint_hem_modules_found_by_line")
# cell_lines <- names(Cortical_lineage_list)
# =========================================================================
Cortical_lineage_list_org <- Cortical_lineage_list
Hem_lineage_list_org <- Hem_lineage_list


Cortical_lineage_list <- Cortical_lineage_list_org 
Hem_lineage_list <- Hem_lineage_list_org  
p_kme_cort_list <- list()
p_kme_hem_list  <- list()

# -----------------------------------------------------------------------
# Explicit, hand-specified module renaming: reorders/relabels modules so
# their NUMBER reflects a common functional category across cell lines
# (1 = translation, 2 = cell cycle, 3 = polarization, 4 = axon development,
# 5 = synapse activity for Pallial; 1 = translation, 2 = epitheliogenesis,
# 3 = neuron differentiation for Hem), matched to JHC1's ordering. "UM" =
# unique module (a module with no clear cross-line functional counterpart).
# Where one original module's function split across two new categories in
# a way that didn't cleanly fit a single number, a "-1"/"-2" suffix is used
# (e.g. M4-1, M4-2) rather than inventing a new top-level number.
# -----------------------------------------------------------------------
pallial_rename_maps <- list(
  JHC1 = c(
    "JHC1-Pallial-M2" = "JHC1-Pallial-M1",
    "JHC1-Pallial-M6" = "JHC1-Pallial-M2",
    "JHC1-Pallial-M3" = "JHC1-Pallial-M4-1",
    "JHC1-Pallial-M5" = "JHC1-Pallial-M4-2",
    "JHC1-Pallial-M1" = "JHC1-Pallial-M5",
    "JHC1-Pallial-M4" = "JHC1-Pallial-UM1"
  ),
  `KOLF2.1` = c(
    "KOLF2.1-Pallial-M1" = "KOLF2.1-Pallial-M1",
    "KOLF2.1-Pallial-M3" = "KOLF2.1-Pallial-M3",
    "KOLF2.1-Pallial-M4" = "KOLF2.1-Pallial-M5-1",
    "KOLF2.1-Pallial-M6" = "KOLF2.1-Pallial-M5-2",
    "KOLF2.1-Pallial-M2" = "KOLF2.1-Pallial-UM1",
    "KOLF2.1-Pallial-M5" = "KOLF2.1-Pallial-UM2"
  ),
  O2C3 = c(
    "O2C3-Pallial-M1" = "O2C3-Pallial-M1",
    "O2C3-Pallial-M4" = "O2C3-Pallial-M2",
    "O2C3-Pallial-M2" = "O2C3-Pallial-M3",
    "O2C3-Pallial-M5" = "O2C3-Pallial-M4",
    "O2C3-Pallial-M3" = "O2C3-Pallial-M5-1",
    "O2C3-Pallial-M7" = "O2C3-Pallial-M5-2",
    "O2C3-Pallial-M6" = "O2C3-Pallial-UM1"
  )
)

hem_rename_maps <- list(
  JHC1 = c(
    "JHC1-Hem-M1" = "JHC1-Hem-M1",
    "JHC1-Hem-M2" = "JHC1-Hem-M2",
    "JHC1-Hem-M3" = "JHC1-Hem-M3",
    "JHC1-Hem-M4" = "JHC1-Hem-UM1"
  ),
  `KOLF2.1` = c(
    "KOLF2.1-Hem-M1" = "KOLF2.1-Hem-M1",
    "KOLF2.1-Hem-M2" = "KOLF2.1-Hem-UM1",
    "KOLF2.1-Hem-M3" = "KOLF2.1-Hem-M3-1",
    "KOLF2.1-Hem-M4" = "KOLF2.1-Hem-M2-1",
    "KOLF2.1-Hem-M5" = "KOLF2.1-Hem-M2-2",
    "KOLF2.1-Hem-M6" = "KOLF2.1-Hem-M3-2"
  ),
  O2C3 = c(
    "O2C3-Hem-M1" = "O2C3-Hem-M1",
    "O2C3-Hem-M2" = "O2C3-Hem-M2",
    "O2C3-Hem-M3" = "O2C3-Hem-M3-1",
    "O2C3-Hem-M4" = "O2C3-Hem-M3-2",
    "O2C3-Hem-M5" = "O2C3-Hem-UM1",
    "O2C3-Hem-M6" = "O2C3-Hem-M3-3",
    "O2C3-Hem-M7" = "O2C3-Hem-M3-4"
  )
)

# -----------------------------------------------------------------------
# Module renaming: pure text relabeling via SetModules() - NOT
# ResetModuleNames()'s named-list mode.
#
# IMPORTANT (found via diagnostic): ResetModuleNames(new_name = <named
# list>) does NOT do a pure text relabel - it silently RECOMPUTES kME
# values (max observed difference 0.595 for the same gene in the same,
# UNRENAMED module M1->M1), corrupting downstream hub-gene selection and
# GO-term enrichment. A pure reorder via SetModules() (same names, just
# different factor level order) was verified to preserve kME EXACTLY
# (diff = 0), so the approach below only ever does simple, provably
# value-safe text substitution: relabel the `module` factor AND rename
# the matching kME_<old-name> column headers (pure colname substitution,
# touches no numeric values) - never a recomputation.
# -----------------------------------------------------------------------
rename_modules_native <- function(obj, wgcna_name, rename_map) {
  ordered_old_names <- names(rename_map)
  
  modules_df <- GetModules(obj, wgcna_name = wgcna_name)
  current_levels <- setdiff(levels(modules_df$module), "grey")
  if (!setequal(current_levels, ordered_old_names)) {
    stop(
      "rename_modules_native(): current module names don't match the ",
      "rename map's expected names - this object may already have been ",
      "renamed in a previous run. Reload it fresh from the Stage 1 ",
      "checkpoint before re-running Stage 2.\n",
      "  Current module levels: ", paste(current_levels, collapse = ", "), "\n",
      "  Rename map expects:    ", paste(ordered_old_names, collapse = ", ")
    )
  }
  
  new_names_ordered <- unname(rename_map[ordered_old_names])
  
  # Relabel the module factor values (pure text substitution) AND reorder
  # to the desired functional sequence in the same step.
  old_char <- as.character(modules_df$module)
  new_char <- ifelse(old_char %in% names(rename_map), unname(rename_map[old_char]), old_char)
  modules_df$module <- factor(new_char, levels = c(new_names_ordered, "grey"))
  
  # Rename the matching kME_<old-name> column HEADERS to kME_<new-name> -
  # a pure colname substitution that does not touch any numeric values.
  old_kme_cols <- paste0("kME_", ordered_old_names)
  new_kme_cols <- paste0("kME_", new_names_ordered)
  has_col <- old_kme_cols %in% colnames(modules_df)
  colnames(modules_df)[match(old_kme_cols[has_col], colnames(modules_df))] <- new_kme_cols[has_col]
  
  obj <- SetModules(obj, modules_df, wgcna_name = wgcna_name)
  obj
}

# GetMEs()'s eigengene matrix is a SEPARATE cached object that SetModules()
# does not touch (confirmed via diagnostic - it is unaffected regardless of
# any module-table renaming), so its columns must be renamed manually. This
# is a pure colname substitution (no values touched), safe for the same
# reason the module-table rename above is safe. Keyed on the RAW/Stage-1
# names, which is exactly what rename_map's keys already are.
rename_ME_columns <- function(MEs, rename_map) {
  colnames(MEs) <- ifelse(colnames(MEs) %in% names(rename_map),
                          unname(rename_map[colnames(MEs)]),
                          colnames(MEs))
  MEs
}

for (cl in cell_lines) {
  
  Cortical_lineage <- Cortical_lineage_list[[cl]]
  Hem_lineage <- Hem_lineage_list[[cl]]
  cl_wgcna <- paste0("trajectory_", cl)
  
  Cortical_lineage <- rename_modules_native(Cortical_lineage, cl_wgcna, pallial_rename_maps[[cl]])
  Hem_lineage <- rename_modules_native(Hem_lineage, cl_wgcna, hem_rename_maps[[cl]])
  
  # GetMEs() does NOT read from the modules table SetModules() just updated
  # - it reads from a SEPARATE internal cache (@misc[[wgcna_name]]$MEs /
  # $hMEs), confirmed via hdWGCNA's own source. PlotModuleTrajectory() and
  # FindDMEs() call GetMEs() internally, so renaming a local copy and
  # cbinding into meta.data (as before) is not enough - those functions
  # would still see the stale internal cache. Fix: rename BOTH the
  # non-harmonized (MEs) and harmonized (hMEs) versions and write them
  # back into the object's real internal cache using hdWGCNA's own
  # (unexported but fully functional) SetMEs().
  MEs_cort_raw <- rename_ME_columns(GetMEs(Cortical_lineage, harmonized = FALSE), pallial_rename_maps[[cl]])
  Cortical_lineage <- hdWGCNA:::SetMEs(Cortical_lineage, MEs_cort_raw, harmonized = FALSE, wgcna_name = cl_wgcna)
  MEs_cort <- rename_ME_columns(GetMEs(Cortical_lineage, harmonized = TRUE), pallial_rename_maps[[cl]])
  Cortical_lineage <- hdWGCNA:::SetMEs(Cortical_lineage, MEs_cort, harmonized = TRUE, wgcna_name = cl_wgcna)
  
  MEs_hem_raw <- rename_ME_columns(GetMEs(Hem_lineage, harmonized = FALSE), hem_rename_maps[[cl]])
  Hem_lineage <- hdWGCNA:::SetMEs(Hem_lineage, MEs_hem_raw, harmonized = FALSE, wgcna_name = cl_wgcna)
  MEs_hem <- rename_ME_columns(GetMEs(Hem_lineage, harmonized = TRUE), hem_rename_maps[[cl]])
  Hem_lineage <- hdWGCNA:::SetMEs(Hem_lineage, MEs_hem, harmonized = TRUE, wgcna_name = cl_wgcna)
  
  meta_cort <- Cortical_lineage@meta.data
  Cortical_lineage@meta.data <- cbind(meta_cort, MEs_cort)
  meta_hem <- Hem_lineage@meta.data
  Hem_lineage@meta.data <- cbind(meta_hem, MEs_hem)
  
  split_pseudotime <- function(df, colname) {
    plus_col  <- paste0(colname, "_plus_pseudotime")
    minus_col <- paste0(colname, "_minus_pseudotime")
    df[[plus_col]]  <- ifelse(df$Protocol == "plus", df[[colname]], NA)
    df[[minus_col]] <- ifelse(df$Protocol == "minus", df[[colname]], NA)
    df
  }
  md_cort <- Cortical_lineage@meta.data
  for (pt in c("dp_pseudotime", "up_pseudotime", "A1_pseudotime", "A2_pseudotime")) md_cort <- split_pseudotime(md_cort, pt)
  md_hem <- Hem_lineage@meta.data
  for (pt in c("crn_pseudotime", "epi_pseudotime")) md_hem <- split_pseudotime(md_hem, pt)
  Cortical_lineage@meta.data <- md_cort
  Hem_lineage@meta.data <- md_hem
  
  # Custom hub-genes-on-x-axis kME plots (Item 2), replacing PlotKMEs().
  # No display_map needed anymore - GetModules() already returns the final
  # fancy names directly.
  p_kme_cort_list[[cl]] <- plot_kme_custom(Cortical_lineage, cl_wgcna, n_hubs = 10,
                                           title = paste0(cl, ": Cortical kME"))
  p_kme_hem_list[[cl]] <- plot_kme_custom(Hem_lineage, cl_wgcna, n_hubs = 10,
                                          title = paste0(cl, ": Hem kME"))
  
  Cortical_lineage_list[[cl]] <- Cortical_lineage
  Hem_lineage_list[[cl]] <- Hem_lineage
}

fig4b_cort_combined <- wrap_plots(p_kme_cort_list, ncol = 1)
ggsave("~/project/IPSC_2025_Data/Figure_5b_cort.png",
       plot = fig4b_cort_combined, device = "png",
       width = 24, height = 8 * length(cell_lines), dpi = 300, limitsize = FALSE)
fig4b_hem_combined <- wrap_plots(p_kme_hem_list, ncol = 1)
ggsave("~/project/IPSC_2025_Data/Figure_5b_hem.png",
       plot = fig4b_hem_combined, device = "png",
       width = 20, height = 8 * length(cell_lines), dpi = 300, limitsize = FALSE)


# -----------------------------------------------------------------------
# Figure 5c: module eigengene trajectories, per cell line, combined
# -----------------------------------------------------------------------
protocol_colors <- met.brewer("Lakota", n = 2, type = 'discrete')
protocol_labels <- c("minus", "plus")

fig4c_cort_list <- list()
fig4c_hem_list  <- list()

for (cl in cell_lines) {
  Cortical_lineage <- Cortical_lineage_list[[cl]]
  Hem_lineage <- Hem_lineage_list[[cl]]
  
  # extra_margin: uniform margins around each panel. strip_fix blanks the
  # module facet labels (leaving the panels themselves) and enlarges the
  # axis tick text. Both are applied AFTER big_text_theme in each panel's
  # chain, and strip.text=element_blank() is repeated as a trailing theme()
  # on each panel (below) so the S7/ggplot2-4.0 theme-merge can't un-blank
  # it via big_text_theme's own strip.text setting.
  extra_margin <- theme(plot.margin = margin(t = 15, r = 15, b = 10, l = 15))
  #strip_fix <- theme(
  #  strip.text = element_blank(),
  #  strip.background = element_blank(),
  #  axis.text = element_text(size = 25)   # enlarge x/y tick text on all Fig4c panels
  #)
  
  strip_fix <- theme(
    strip.text = element_text(size = 11),      # small, plain labels (no bold, no grey box)
    strip.background = element_blank(),
    axis.text = element_text(size = 18)
  )
  
  # Fixed facet columns for BOTH lineages so every module sub-panel is the
  # same physical size regardless of Hem vs Pallial (Hem previously used
  # ncol=5, Pallial ncol=7, which made their facets different sizes).
  facet_ncol <- 7
  blank_strip <- theme(strip.text = element_blank())   # trailing override (S7 merge-proof)
  
  p_dp <- PlotModuleTrajectory(Cortical_lineage, pseudotime_col = c("dp_pseudotime_plus_pseudotime", "dp_pseudotime_minus_pseudotime"), group_colors = protocol_colors, ncol = facet_ncol) +
    scale_color_manual(values = protocol_colors, labels = protocol_labels) +
    ggtitle(paste0(cl, ": DP lineage eigengene")) + big_text_theme + extra_margin + strip_fix + NoLegend() + blank_strip
  p_up <- PlotModuleTrajectory(Cortical_lineage, pseudotime_col = c("up_pseudotime_plus_pseudotime", "up_pseudotime_minus_pseudotime"), group_colors = protocol_colors, ncol = facet_ncol) +
    scale_color_manual(values = protocol_colors, labels = protocol_labels) +
    ggtitle(paste0(cl, ": UP lineage eigengene")) + big_text_theme + extra_margin + strip_fix + blank_strip
  p_A1 <- PlotModuleTrajectory(Cortical_lineage, pseudotime_col = c("A1_pseudotime_plus_pseudotime", "A1_pseudotime_minus_pseudotime"), group_colors = protocol_colors, ncol = facet_ncol) +
    scale_color_manual(values = protocol_colors, labels = protocol_labels) +
    ggtitle(paste0(cl, ": A1 lineage eigengene")) + big_text_theme + extra_margin + strip_fix + NoLegend() + blank_strip
  p_A2 <- PlotModuleTrajectory(Cortical_lineage, pseudotime_col = c("A2_pseudotime_plus_pseudotime", "A2_pseudotime_minus_pseudotime"), group_colors = protocol_colors, ncol = facet_ncol) +
    scale_color_manual(values = protocol_colors, labels = protocol_labels) +
    ggtitle(paste0(cl, ": A2 lineage eigengene")) + big_text_theme + extra_margin + strip_fix + blank_strip
  
  p_crn <- PlotModuleTrajectory(Hem_lineage, pseudotime_col = c("crn_pseudotime_plus_pseudotime", "crn_pseudotime_minus_pseudotime"), group_colors = protocol_colors, ncol = facet_ncol) +
    scale_color_manual(values = protocol_colors, labels = protocol_labels) +
    ggtitle(paste0(cl, ": CRN lineage eigengene")) + big_text_theme + extra_margin + strip_fix + blank_strip
  p_epi <- PlotModuleTrajectory(Hem_lineage, pseudotime_col = c("epi_pseudotime_plus_pseudotime", "epi_pseudotime_minus_pseudotime"), group_colors = protocol_colors, ncol = facet_ncol) +
    scale_color_manual(values = protocol_colors, labels = protocol_labels) +
    ggtitle(paste0(cl, ": Epithelial lineage eigengene")) + big_text_theme + extra_margin + strip_fix + NoLegend() + blank_strip
  
  fig4c_cort_list[[cl]] <- (p_dp + p_up) / (p_A1 + p_A2)
  fig4c_hem_list[[cl]]  <- (p_crn + p_epi)
}

fig4c_cort_combined_heights <- sapply(cell_lines, function(cl) {
  n_mod <- length(pallial_rename_maps[[cl]])
  2 * ceiling(n_mod / 7)   # 2 stacked row-groups (dp+up, A1+A2), each needing ceiling(n_mod/7) internal rows
})
fig4c_cort_combined <- wrap_plots(fig4c_cort_list, ncol = 1, heights = fig4c_cort_combined_heights) +
  plot_annotation(
    title = "Effect of Protocol on Cortical Lineage Modules, per Cell Line",
    theme = theme(plot.title = element_text(size = 24, face = "bold", margin = margin(b = 20, t = 10)))
  )
ggsave("~/project/IPSC_2025_Data/Figure_5c_cort.png",
       plot = fig4c_cort_combined, device = "png",
       width = 30, height = 11 * sum(fig4c_cort_combined_heights) / length(cell_lines), dpi = 300, limitsize = FALSE)

fig4c_hem_combined_heights <- sapply(cell_lines, function(cl) {
  n_mod <- length(hem_rename_maps[[cl]])
  ceiling(n_mod / 7)   # 1 row-group (crn+epi), same facet_ncol=7 as Pallial so facet size matches
})
fig4c_hem_combined <- wrap_plots(fig4c_hem_list, ncol = 1, heights = fig4c_hem_combined_heights) +
  plot_annotation(
    title = "Effect of Protocol on Hem Lineage Modules, per Cell Line",
    theme = theme(plot.title = element_text(size = 24, face = "bold", margin = margin(b = 20, t = 10)))
  )
ggsave("~/project/IPSC_2025_Data/Figure_5c_hem.png",
       plot = fig4c_hem_combined, device = "png",
       width = 30, height = 11 * sum(fig4c_hem_combined_heights) / length(cell_lines), dpi = 300, limitsize = FALSE)


# -----------------------------------------------------------------------
# Figure 5d: DME lollipop plots, per cell line, combined
# -----------------------------------------------------------------------
run_dme <- function(obj, lineage_col, lineage_value, label, wgcna_name) {
  group1 <- obj@meta.data %>% subset(.[[lineage_col]] == lineage_value & Protocol == "plus") %>% rownames()
  group2 <- obj@meta.data %>% subset(.[[lineage_col]] == lineage_value & Protocol == "minus") %>% rownames()
  res <- FindDMEs(obj, barcodes1 = group1, barcodes2 = group2, test.use = "wilcox",
                  pseudocount.use = 0.01, wgcna_name = wgcna_name)
  res$lineage <- label
  res
}

fig4d_cort_list <- list()
fig4d_hem_list  <- list()

for (cl in cell_lines) {
  
  Cortical_lineage <- Cortical_lineage_list[[cl]]
  Hem_lineage <- Hem_lineage_list[[cl]]
  cl_wgcna <- paste0("trajectory_", cl)
  
  Cortical_lineage$dp_lineage <- ifelse(!is.na(Cortical_lineage$dp_pseudotime), "dp_ExN_lineage", NA)
  Cortical_lineage$up_lineage <- ifelse(!is.na(Cortical_lineage$up_pseudotime), "up_ExN_lineage", NA)
  Cortical_lineage$A1_lineage <- ifelse(!is.na(Cortical_lineage$A1_pseudotime), "A1_Astrocyte_lineage", NA)
  Cortical_lineage$A2_lineage <- ifelse(!is.na(Cortical_lineage$A2_pseudotime), "A2_Astrocyte_lineage", NA)
  
  DMEs_Cortical <- bind_rows(
    run_dme(Cortical_lineage, "dp_lineage", "dp_ExN_lineage", "DP_ExN_lineage", cl_wgcna),
    run_dme(Cortical_lineage, "up_lineage", "up_ExN_lineage", "UP_ExN_lineage", cl_wgcna),
    run_dme(Cortical_lineage, "A1_lineage", "A1_Astrocyte_lineage", "A1_Astrocyte_lineage", cl_wgcna),
    run_dme(Cortical_lineage, "A2_lineage", "A2_Astrocyte_lineage", "A2_Astrocyte_lineage", cl_wgcna)
  )
  
  facet_order <- c("DP_ExN_lineage", "UP_ExN_lineage", "A1_Astrocyte_lineage", "A2_Astrocyte_lineage")
  # DMEs_Cortical$module comes from FindDMEs(), which now reports the
  # object's ACTUAL current module names directly - already the final
  # fancy display names, since ResetModuleNames() renamed them at the
  # source (no separate relabeling step needed here anymore).
  # module_order is sourced from the object's actual module levels (not a
  # hardcoded count), preserving the functional order already established
  # by rename_modules_native()'s reordering step.
  module_order <- setdiff(levels(GetModules(Cortical_lineage, wgcna_name = cl_wgcna)$module), "grey")
  plot_df <- DMEs_Cortical %>%
    filter(module %in% module_order) %>%   # drop any leftover "grey"/unassigned rows
    mutate(sig_flag = p_val_adj < 0.05,
           module = factor(module, levels = module_order),
           lineage = factor(lineage, levels = facet_order))
  
  fig4d_cort_list[[cl]] <- ggplot(plot_df, aes(x = avg_log2FC, y = module)) +
    geom_segment(aes(x = 0, xend = avg_log2FC, y = module, yend = module), color = "grey60", linewidth = 0.6) +
    geom_point(aes(color = module), size = 3, data = subset(plot_df, sig_flag)) +
    geom_text(data = subset(plot_df, !sig_flag), aes(label = "×"), color = "black", size = 5, fontface = "bold") +
    facet_wrap(~ lineage, scales = "free_y") +
    coord_cartesian(xlim = c(-4, 4)) +   # fixed x-scale, same across ALL cell lines (Pallial)
    theme_minimal(base_size = 14) +
    theme(strip.text = element_text(size = 16, face = "bold"), axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12),
          plot.title = element_text(size = 18, face = "bold", margin = margin(b = 12)),
          plot.margin = margin(t = 20, r = 10, b = 10, l = 10)) +
    labs(title = paste0(cl, ": DME Lollipop Plot by lineage"), x = "avg_log2FC", y = "module")
  
  Hem_lineage$crn_lineage <- ifelse(!is.na(Hem_lineage$crn_pseudotime), "crn_lineage", NA)
  Hem_lineage$epi_lineage <- ifelse(!is.na(Hem_lineage$epi_pseudotime), "epi_lineage", NA)
  
  DMEs_Hem <- bind_rows(
    run_dme(Hem_lineage, "crn_lineage", "crn_lineage", "CRN_lineage", cl_wgcna),
    run_dme(Hem_lineage, "epi_lineage", "epi_lineage", "Epithelial_lineage", cl_wgcna)
  )
  
  facet_order_hem <- c("Epithelial_lineage", "CRN_lineage")
  # Same as above: source module order directly from the object's actual
  # (already fancy-renamed) module levels.
  module_order_hem <- setdiff(levels(GetModules(Hem_lineage, wgcna_name = cl_wgcna)$module), "grey")
  plot_df_hem <- DMEs_Hem %>%
    filter(module %in% module_order_hem) %>%
    mutate(sig_flag = p_val_adj < 0.05,
           module = factor(module, levels = module_order_hem),
           lineage = factor(lineage, levels = facet_order_hem))
  
  fig4d_hem_list[[cl]] <- ggplot(plot_df_hem, aes(x = avg_log2FC, y = module)) +
    geom_segment(aes(x = 0, xend = avg_log2FC, y = module, yend = module), color = "grey60", linewidth = 0.6) +
    geom_point(aes(color = module), size = 3, data = subset(plot_df_hem, sig_flag)) +
    geom_text(data = subset(plot_df_hem, !sig_flag), aes(label = "×"), color = "black", size = 5, fontface = "bold") +
    facet_wrap(~ lineage, scales = "free_y") +
    coord_cartesian(xlim = c(-1.25, 1.25)) +   # fixed x-scale, same across ALL cell lines (Hem)
    theme_minimal(base_size = 14) +
    theme(strip.text = element_text(size = 16, face = "bold"), axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12),
          plot.title = element_text(size = 18, face = "bold", margin = margin(b = 12)),
          plot.margin = margin(t = 20, r = 10, b = 10, l = 10)) +
    labs(title = paste0(cl, ": DME Lollipop Plot by lineage"), x = "avg_log2FC", y = "module")
}

fig4d_cort_combined <- wrap_plots(fig4d_cort_list, ncol = 1)
ggsave("~/project/IPSC_2025_Data/Figure_5d_cort_lineage.tiff",
       plot = fig4d_cort_combined, device = "tiff",
       width = 8, height = 5 * length(cell_lines), dpi = 300, limitsize = FALSE)

fig4d_hem_combined <- wrap_plots(fig4d_hem_list, ncol = 1)
ggsave("~/project/IPSC_2025_Data/Figure_5d_hem_lineage.tiff",
       plot = fig4d_hem_combined, device = "tiff",
       width = 8, height = 5 * length(cell_lines), dpi = 300, limitsize = FALSE)



Cortical_lineage_list_org
Hem_lineage_list_org

# -----------------------------------------------------------------------
# Supplementary Figure 6: Enrichr, per cell line, combined
# -----------------------------------------------------------------------


# NEW helper - insert just above build_supp_fig7_panel()
build_original_order_table <- function(rename_map_cl) {
  orig_names   <- names(rename_map_cl)
  fancy_names  <- unname(rename_map_cl)
  numeric_part <- as.numeric(sub(".*-M(\\d+)$", "\\1", orig_names))
  ord <- order(numeric_part)
  tibble(orig = orig_names[ord], fancy = fancy_names[ord])
}

# REPLACES the old build_supp_fig7_panel() definition
build_supp_fig7_panel <- function(obj, cl_wgcna, order_table, category_lookup, title) {
  module_order_fancy <- order_table$fancy
  module_order_label <- order_table$orig
  
  dotplot <- EnrichrDotPlot(
    obj, mods = "all", database = "GO_Biological_Process_2023",
    n_terms = 5, term_size = 8, p_adj = FALSE, wgcna_name = cl_wgcna
  ) +
    scale_color_stepsn(colors = rev(viridis::magma(256))) +
    scale_x_discrete(limits = module_order_fancy, labels = module_order_label) +
    ggtitle(title) +
    big_text_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  strip <- plot_category_strip(module_order_fancy, category_lookup)
  
  dotplot / strip + plot_layout(heights = c(10, 1))
}

# REPLACES the old for-loop body (note the two new order_table lines)
fig_s7_cort_list <- list()
fig_s7_hem_list  <- list()
for (cl in cell_lines) {
  cl_wgcna <- paste0("trajectory_", cl)
  
  cort_order <- build_original_order_table(pallial_rename_maps[[cl]])
  hem_order  <- build_original_order_table(hem_rename_maps[[cl]])
  
  fig_s7_cort_list[[cl]] <- build_supp_fig7_panel(
    Cortical_lineage_list[[cl]], cl_wgcna, cort_order, module_categories_pallial, paste0(cl, ": Cortical")
  )
  fig_s7_hem_list[[cl]] <- build_supp_fig7_panel(
    Hem_lineage_list[[cl]], cl_wgcna, hem_order, module_categories_hem, paste0(cl, ": Hem")
  )
}

# ggsave() calls stay exactly as before - unchanged
fig_s7_cort_combined <- wrap_plots(fig_s7_cort_list, ncol = 1)
ggsave("~/project/IPSC_2025_Data/Supplmental_Figure6_pallial_lineages.png",
       plot = fig_s7_cort_combined, device = "png",
       width = 10, height = 11 * length(cell_lines), dpi = 300, limitsize = FALSE)

fig_s7_hem_combined <- wrap_plots(fig_s7_hem_list, ncol = 1)
ggsave("~/project/IPSC_2025_Data/Supplmental_Figure6_hem_lineages.png",
       plot = fig_s7_hem_combined, device = "png",
       width = 10, height = 11 * length(cell_lines), dpi = 300, limitsize = FALSE)


install.packages("ggalluvial")

library(ggalluvial)
library(ggfittext)
library(patchwork)
library(scales)
library(stringr)
library(purrr)
library(tibble)
library(dplyr)

# =========================================================================
# 1. Category lookup tables, keyed by the EXACT final module name.
#    NOTE: several tags (M4-1/M4-2, M5-1/M5-2, M3-1..M3-4) are kept as
#    DISTINCT categories rather than collapsed back to one base label -
#    the enrichment showed these are genuinely different biology, not the
#    same theme split at two intensities. Modules with no term below
#    Adjusted P < 0.05 are labeled "weak enrichment" rather than given a
#    confident biological name.
# =========================================================================
module_categories_pallial <- c(
  "JHC1-Pallial-M1"     = "Translation",
  "JHC1-Pallial-M2"     = "Cell cycle",
  "JHC1-Pallial-M4-1"   = "Ciliogenesis",
  "JHC1-Pallial-M4-2"   = "Axon/neurite outgrowth",
  "JHC1-Pallial-M5"     = "Axon guidance & synapse activity",
  "JHC1-Pallial-UM1"    = "Unique module",
  
  "KOLF2.1-Pallial-M1"    = "Translation",
  "KOLF2.1-Pallial-M3"    = "Apical-basal polarity",
  "KOLF2.1-Pallial-M5-1"  = "Axon guidance & development",
  "KOLF2.1-Pallial-M5-2"  = "Synapse assembly & activity",
  "KOLF2.1-Pallial-UM1"   = "Unique module (weak enrichment)",
  "KOLF2.1-Pallial-UM2"   = "Unique module (immune-like)",
  
  "O2C3-Pallial-M1"    = "Translation",
  "O2C3-Pallial-M2"    = "Cell cycle",
  "O2C3-Pallial-M3"    = "Apical-basal polarity",
  "O2C3-Pallial-M4"    = "Adhesion / planar polarity",
  "O2C3-Pallial-M5-1"  = "Axon guidance & development",
  "O2C3-Pallial-M5-2"  = "Synapse assembly & activity",
  "O2C3-Pallial-UM1"   = "Unique module (lipid/immune)"
)

module_categories_hem <- c(
  "JHC1-Hem-M1"   = "Translation",
  "JHC1-Hem-M2"   = "Ciliogenesis",
  "JHC1-Hem-M3"   = "Neuron differentiation",
  "JHC1-Hem-UM1"  = "Unique module (weak enrichment)",
  
  "KOLF2.1-Hem-M1"    = "Translation",
  "KOLF2.1-Hem-M2-1"  = "Ciliogenesis",
  "KOLF2.1-Hem-M2-2"  = "Unique module (weak enrichment)",
  "KOLF2.1-Hem-M3-1"  = "Neuron differentiation",
  "KOLF2.1-Hem-M3-2"  = "Unique module (weak enrichment)",
  "KOLF2.1-Hem-UM1"   = "Neuron differentiation (unmerged)",
  
  "O2C3-Hem-M1"    = "Translation",
  "O2C3-Hem-M2"    = "Ciliogenesis",
  "O2C3-Hem-M3-1"  = "Neuron differentiation",
  "O2C3-Hem-M3-2"  = "Neuron differentiation",
  "O2C3-Hem-M3-3"  = "Neuron differentiation (Wnt+)",
  "O2C3-Hem-M3-4"  = "Neuron differentiation (Wnt-)",
  "O2C3-Hem-UM1"   = "Unique module (interferon-like)"
)

all_categories <- unique(c(module_categories_pallial, module_categories_hem))
category_colors <- setNames(scales::hue_pal()(length(all_categories)), all_categories)

# =========================================================================
# 2. Category strip - a thin colored bar under the dotplot's x-axis,
#    mapping each module (in the SAME order as the dotplot) to its category
# =========================================================================
plot_category_strip <- function(module_order, category_lookup) {
  strip_df <- tibble(
    module = factor(module_order, levels = module_order),
    Category = factor(category_lookup[module_order], levels = all_categories)
  )
  ggplot(strip_df, aes(x = module, y = 1, fill = Category)) +
    geom_tile(height = 1, color = "white", linewidth = 0.3) +
    scale_fill_manual(values = category_colors, drop = FALSE, name = "Category") +
    scale_x_discrete(limits = module_order) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 9),
          plot.margin = margin(t = 0, r = 5, b = 5, l = 5))
}


# =========================================================================
# 3. Grouping-logic alluvial plot, using the refined per-module categories
#    (joined against the rename maps you already have in your script)
# =========================================================================
build_module_grouping_df <- function(rename_maps, category_lookup, lineage_label) {
  map_dfr(names(rename_maps), function(cl) {
    m <- rename_maps[[cl]]
    tibble(CellLine = cl, Lineage = lineage_label,
           OldModule = names(m), NewModule = unname(m))
  }) %>%
    mutate(Category = category_lookup[NewModule])
}

pallial_grouping_df <- build_module_grouping_df(pallial_rename_maps, module_categories_pallial, "Pallial")
hem_grouping_df     <- build_module_grouping_df(hem_rename_maps,     module_categories_hem,     "Hem")

plot_module_grouping_alluvial <- function(df, title) {
  df$Category <- factor(df$Category, levels = all_categories)
  ggplot(df, aes(axis1 = CellLine, axis2 = OldModule, axis3 = Category)) +
    geom_alluvium(aes(fill = Category), width = 1/4, alpha = 0.85) +
    geom_stratum(width = 1/4, fill = "grey95", color = "grey30") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    scale_x_discrete(limits = c("Cell line", "Original module", "Functional category"),
                     expand = c(0.08, 0.08)) +
    scale_fill_manual(values = category_colors, drop = FALSE) +
    labs(title = title, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid = element_blank(), legend.position = "right",
          legend.text = element_text(size = 8))
}

p_grouping_pallial <- plot_module_grouping_alluvial(pallial_grouping_df, "Pallial module grouping logic")
p_grouping_hem     <- plot_module_grouping_alluvial(hem_grouping_df,     "Hem module grouping logic")

fig_s7_grouping_combined <- p_grouping_pallial / p_grouping_hem +
  plot_annotation(title = "Module renaming: original modules grouped by refined GO term theme")

ggsave("~/project/IPSC_2025_Data/Supplmental_Figure6_grouping_logic.png",
       plot = fig_s7_grouping_combined, device = "png",
       width = 12, height = 12, dpi = 300, limitsize = FALSE)

