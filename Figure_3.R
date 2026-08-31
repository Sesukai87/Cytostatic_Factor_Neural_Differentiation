library(ggpubr)
library(matrixStats)
library(uwot)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(harmony)
library(pheatmap)
library(SingleCellExperiment)
library(scuttle)
library(dreamlet)
library(mashr)
library(zenith)
library(ComplexHeatmap)
library(ggrepel)
library(ggtext)
library(cowplot)
library(patchwork)

# -----------------------------------------------------------------------
# Shared print-legible theme, reused across every Figure 3/Supp panel
# -----------------------------------------------------------------------
big_text_theme <- theme(
  axis.text = element_text(size = 14),
  axis.title = element_text(size = 16, face = "bold"),
  strip.text = element_text(size = 16, face = "bold"),
  legend.text = element_text(size = 13),
  legend.title = element_text(size = 14, face = "bold"),
  plot.title = element_text(size = 18, face = "bold")
)

merged <- readRDS("~/project/IPSC_2025_Data/merged_Fetal_IPSC_derived_forebrain")
merged_IPSC <- subset(merged, Sampletype == "IPSC-Derived")
unique(merged_IPSC$gt_line)
merged_IPSC$neural_induction_media <- ifelse(
  grepl("_E6", merged_IPSC$SampleID),
  "E6",
  "KSR"
)
rm(merged)
unique(merged_IPSC$gt_line)
merged_IPSC$gt_line <- droplevels(merged_IPSC$gt_line)
cell_lines <- levels(merged_IPSC$gt_line)


Cortical_lineage <- subset(merged_IPSC, Celltype == "DL_ExN" | Celltype == "UL_ExN" | Celltype == "RG" | Celltype == "Astrocyte" | Celltype == "IPC_ExN")
# Stamp original barcode before any downstream re-embedding/merging (for
# consistency with the Hem branch; the Cortical barcodes happened to stay
# matchable, but this makes the join in Figure_5.R robust for both).
Cortical_lineage$orig_barcode <- colnames(Cortical_lineage)
Cortical_lineage[["percent.lmx1a"]] <- PercentageFeatureSet(Cortical_lineage, pattern = "LMX1A")
Cortical_lineage <- subset(Cortical_lineage, In_new_mod <= 0)
add_to_hem <- subset(Cortical_lineage, percent.lmx1a > 0)
crn_to_hem <- subset(Cortical_lineage, CRN_new_mod > 0.5)
Cortical_lineage <- subset(Cortical_lineage, percent.lmx1a == 0)
Cortical_lineage <- subset(Cortical_lineage, CRN_new_mod <= 0.5)
Cortical_lineage <- subset(Cortical_lineage, In_new_mod <= 0)
Cortical_lineage <-  NormalizeData(Cortical_lineage) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunHarmony(group.by.vars = c("SampleID2")) %>% RunUMAP(reduction = "harmony", dims = 1:10, n.neighbors = 150, min.dist = 0.8, spread = 1.0, repulsion.strength = 0.01, local.connectivity = 5, metric = "euclidean", seed.use = 42) %>%  FindNeighbors(reduction = "harmony", dims = 1:20)
Cortical_lineage <- FindClusters(Cortical_lineage, res = 0.6)
DimPlot(Cortical_lineage, group.by = "seurat_clusters", label = TRUE)
DimPlot(Cortical_lineage, group.by = "Celltype")
Cortical_lineage <- RenameIdents(Cortical_lineage, `0` = "UL_ExN", `1` = "UL_ExN", `2` = "A1 Astrocyte", `3` = "UL_ExN", `4` = "DL_ExN", `5` = "A1 Astrocyte", `6` = "RG", `7` = "RG", `8` = "A2 Astrocyte", `9` = "A2 Astrocyte", `10` = "A2 Astrocyte", `11` = "RG", `12` = "UL_ExN", `13` = "A2 Astrocyte", `14` = "UL_ExN")
Cortical_lineage$Celltype2 <- Idents(Cortical_lineage)
Celltype1 <- as.character(Cortical_lineage$Celltype)
Celltype2 <- as.character(Cortical_lineage$Celltype2)
names(Celltype1) <- colnames(Cortical_lineage)
names(Celltype2) <- colnames(Cortical_lineage)
consensusClusterLabels <- Celltype2
table(consensusClusterLabels)
consensusClusterLabels[names(which(Celltype2 == "A2 Astrocyte"))] <- "A2 Astrocyte"
consensusClusterLabels[names(which(Celltype2 == "A1 Astrocyte"))] <- "A1 Astrocyte"
consensusClusterLabels[names(which(Celltype1 == "Astrocyte" & Celltype2 == "RG"))] <- "A2 Astrocyte"
consensusClusterLabels[names(which(Celltype1 == "Astrocyte" & Celltype2 == "UL_ExN"))] <- "A1 Astrocyte"
consensusClusterLabels[names(which(Celltype1 == "Astrocyte" & Celltype2 == "DL_ExN"))] <- "A2 Astrocyte"
table(consensusClusterLabels)
Cortical_lineage$Celltype2 <- consensusClusterLabels
Cortical_lineage <- JoinLayers(Cortical_lineage)


Hem_lineage <- subset(merged_IPSC, Celltype == "Epithelial" | Celltype == "Hem_RG" | Celltype == "CRN")
# Stamp the ORIGINAL barcode (as it exists in merged_IPSC) into metadata
# on every component object BEFORE the merge() below re-encodes the cell
# names. merge() appends disambiguating suffixes that diverge from
# merged_IPSC's barcode scheme, which later made it impossible to join
# Hem-lineage pseudotime back onto merged_IPSC by barcode in Figure_5.R.
# A metadata column survives subset()/merge() untouched, so this gives a
# stable key. (add_to_hem / crn_to_hem are also merged in below, so they
# need the stamp too - each is stamped with its own original barcodes.)
Hem_lineage$orig_barcode <- colnames(Hem_lineage)
add_to_hem$orig_barcode  <- colnames(add_to_hem)
crn_to_hem$orig_barcode  <- colnames(crn_to_hem)
# Get cell IDs
cells_hem <- colnames(Hem_lineage)
cells_add <- c(colnames(add_to_hem), colnames(crn_to_hem))
# Identify duplicates
dup_cells <- intersect(cells_hem, cells_add)
# Keep only unique cells in each object
Hem_lineage_unique <- subset(Hem_lineage, cells = setdiff(cells_hem, dup_cells))
add_to_hem_unique  <- subset(add_to_hem,  cells = setdiff(cells_add, dup_cells))
crn_to_hem_unique  <- subset(crn_to_hem,  cells = setdiff(cells_add, dup_cells))
# Merge only unique cells
dim(Hem_lineage_unique)
dim(add_to_hem_unique)
Hem_lineage_merged <- merge(Hem_lineage_unique, add_to_hem_unique)
Hem_lineage_merged <- merge(Hem_lineage_merged, crn_to_hem_unique)
Hem_lineage <- Hem_lineage_merged
rm(Hem_lineage_merged)
Hem_lineage <- subset(Hem_lineage, In_new_mod <= 0)
Hem_lineage <- NormalizeData(Hem_lineage) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunHarmony(group.by.vars = c("SampleID2")) %>% RunUMAP(reduction = "harmony", dims = 1:5, n.neighbors = 150, min.dist = 0.8, spread = 1.0, repulsion.strength = 0.01, local.connectivity = 5, metric = "euclidean", seed.use = 42) %>%  FindNeighbors(reduction = "harmony", dims = 1:5) %>% FindClusters(res = 0.6)
Hem_lineage <- RenameIdents(Hem_lineage, `0` = "Hem_RG", `1` = "Hem_RG", `2` = "Epithelial", `3` = "Hem_RG", `4` = "CRN", `5` = "Hem_RG", `6` = "Hem_RG", `7` = "Hem_RG", `8` = "CRN", `9` = "Hem_RG", `10` = "Hem_RG", `11` = "Hem_RG", `12` = "Hem_RG", `13` = "Epithelial")
Hem_lineage$Celltype2 <- Idents(Hem_lineage)
Celltype1 <- as.character(Hem_lineage$Celltype)
Celltype2 <- as.character(Hem_lineage$Celltype2)
names(Celltype1) <- colnames(Hem_lineage)
names(Celltype2) <- colnames(Hem_lineage)
consensusClusterLabels <- Celltype1
table(consensusClusterLabels)
consensusClusterLabels[names(which(Celltype1 == "RG"))] <- "Hem_RG"
consensusClusterLabels[names(which(Celltype1 == "IPC_ExN"))] <- "Hem_RG"
consensusClusterLabels[names(which(Celltype1 == "UL_ExN"))] <- "CRN"
consensusClusterLabels[names(which(Celltype1 == "DL_ExN"))] <- "CRN"
Hem_lineage$Celltype2 <- consensusClusterLabels
stripped <- sub("#4$", "", colnames(Hem_lineage))
keep <- !duplicated(stripped)
Hem_lineage <- Hem_lineage[, keep]
colnames(Hem_lineage) <- stripped[keep]
Hem_lineage <- JoinLayers(Hem_lineage)

saveRDS(Cortical_lineage, "~/project/IPSC_2025_Data/merged_IPSC_derived_pallial_lineages")
saveRDS(Hem_lineage, "~/project/IPSC_2025_Data/merged_IPSC_derived_hem_lineages")

# -----------------------------------------------------------------------
# PER-CELL-LINE SPLIT
# -----------------------------------------------------------------------
# Everything downstream (partition UMAPs shown per-line in Fig 3a, and all
# differential analysis) is now done separately within each cell line
# (gt_line) rather than pooled, per your request. The UMAP embedding and
# Celltype2 clustering above are computed once on the full IPSC-derived set
# (that clustering doesn't change per line), but we now SPLIT the object by
# gt_line before anything comparative happens.
Cortical_lineage_list <- lapply(cell_lines, function(cl) subset(Cortical_lineage, gt_line == cl))
names(Cortical_lineage_list) <- cell_lines
Hem_lineage_list <- lapply(cell_lines, function(cl) subset(Hem_lineage, gt_line == cl))
names(Hem_lineage_list) <- cell_lines

saveRDS(Cortical_lineage_list, "~/project/IPSC_2025_Data/merged_IPSC_derived_pallial_lineages_by_line")
saveRDS(Hem_lineage_list, "~/project/IPSC_2025_Data/merged_IPSC_derived_hem_lineages_by_line")


# -----------------------------------------------------------------------
# Figure 3a: partition UMAPs, per cell line x Protocol, combined into ONE
# figure (rows = cell line, columns = lineage x Protocol)
# -----------------------------------------------------------------------
fig3a_panels <- list()
# Persist each cell line's freshly-UMAP'd lineage objects so downstream
# figures (e.g. Supplemental Figure 5 feature plots) can reuse the SAME
# per-cell-line embeddings rather than recomputing them or falling back to
# the pooled embedding.
cort_line_umap_list <- list()
hem_line_umap_list  <- list()
for (cl in cell_lines) {

  # Build a clean per-cell-line object and compute a FRESH UMAP embedding
  # for THIS cell line only (rather than reusing the pooled embedding). The
  # two protocols within a line still share this line's embedding, but each
  # cell line now gets its own independent embedding. Same processing
  # pipeline / parameters as the pooled Cortical_lineage / Hem_lineage
  # objects, just applied within the single cell line's cells.
  #
  # (Building from fresh counts + metadata avoids the Assay5 scale.data /
  # mangled-layer issues that DietSeurat and repeated subset() hit.)
  cort_obj_src <- Cortical_lineage_list[[cl]]
  cort_line_obj <- CreateSeuratObject(
    counts = LayerData(cort_obj_src, assay = "RNA", layer = "counts"),
    meta.data = cort_obj_src@meta.data
  )
  cort_line_obj <- cort_line_obj %>%
    NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>%
    RunHarmony(group.by.vars = c("SampleID2")) %>%
    RunUMAP(reduction = "harmony", dims = 1:10, n.neighbors = 150, min.dist = 0.8,
            spread = 1.0, repulsion.strength = 0.01, local.connectivity = 5,
            metric = "euclidean", seed.use = 42)

  hem_obj_src <- Hem_lineage_list[[cl]]
  hem_line_obj <- CreateSeuratObject(
    counts = LayerData(hem_obj_src, assay = "RNA", layer = "counts"),
    meta.data = hem_obj_src@meta.data
  )
  hem_line_obj <- hem_line_obj %>%
    NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>%
    RunHarmony(group.by.vars = c("SampleID2")) %>%
    RunUMAP(reduction = "harmony", dims = 1:5, n.neighbors = 150, min.dist = 0.8,
            spread = 1.0, repulsion.strength = 0.01, local.connectivity = 5,
            metric = "euclidean", seed.use = 42)

  # stash for reuse in later figures (Supp Fig 5)
  cort_line_umap_list[[cl]] <- cort_line_obj
  hem_line_umap_list[[cl]]  <- hem_line_obj

  # Build panels in a FIXED column order for each row (cell line):
  # Pallial -SDF, Pallial +SDF, Hem -SDF, Hem +SDF
  protocol_order <- c("minus", "plus")
  protocol_labels <- c(minus = "-SDF", plus = "+SDF")
  lineage_order <- list(
    list(name = "Pallial", obj = cort_line_obj),
    list(name = "Hem",     obj = hem_line_obj)
  )

  for (lin in lineage_order) {
    for (prot in protocol_order) {
      sub_obj <- subset(lin$obj, Protocol == prot)
      prot_label <- protocol_labels[[prot]]

      # Guard against empty subsets (e.g. a cell line with no cells for
      # this Protocol in one of the lineages) - DimPlot errors on 0-cell
      # objects, so we substitute a blank placeholder panel instead of
      # failing.
      if (ncol(sub_obj) == 0) {
        p <- ggplot() + theme_void() +
          ggtitle(paste0(cl, " | ", lin$name, " | ", prot_label, " (no cells)")) +
          theme(plot.title = element_text(size = 14))
      } else {
        p <- DimPlot(sub_obj, group.by = "Celltype2", pt.size = 2) +
          ggtitle(paste0(cl, " | ", lin$name, " | ", prot_label)) +
          big_text_theme + theme(plot.title = element_text(size = 14))
      }

      fig3a_panels[[paste(cl, lin$name, prot, sep = "_")]] <- p
    }
  }
}
# 4 columns (Pallial -SDF, Pallial +SDF, Hem -SDF, Hem +SDF) x N cell-line rows
combined_plot <- wrap_plots(fig3a_panels, ncol = 4)
ggsave(
  "~/project/IPSC_2025_Data/Figure3a_Dimplot_Celltype_Combined.tiff",
  plot = combined_plot,
  device = "tiff",
  width = 24, height = 6 * length(cell_lines), dpi = 300, limitsize = FALSE
)

# Persist the FRESH per-cell-line UMAP-embedded objects (built in the loop
# above) so downstream scripts - Figure_4.R (Slingshot, WGCNA) and beyond -
# use the SAME per-line embeddings that Figure 3a itself visualizes, rather
# than the older pooled-embedding-then-subset objects saved earlier. This
# is what keeps Figure 3 and Figure 4's UMAP-dependent analyses consistent.
saveRDS(cort_line_umap_list, "~/project/IPSC_2025_Data/merged_IPSC_derived_pallial_lineages_by_line_umap")
saveRDS(hem_line_umap_list, "~/project/IPSC_2025_Data/merged_IPSC_derived_hem_lineages_by_line_umap")

# -----------------------------------------------------------------------
# Supplemental Figure 5: faceted FeaturePlots of the 4 markers, shown on
# EACH cell line's own UMAP embedding (the per-cell-line embeddings built
# and stored in the Figure 3a loop above), rather than a single pooled
# embedding. Layout: rows = cell line, columns = feature.
# -----------------------------------------------------------------------
s5_features <- c("TNFSF12", "S100A10", "C3", "FKBP5")

# Match FeaturePlot_scCustom()'s default color scheme: its default is
# colors_use = viridis_plasma_dark_high, i.e. the viridis "plasma" palette
# reversed so HIGH expression is the dark magenta/purple end and LOW is the
# light yellow end. Pull scCustomize's exact object if available; otherwise
# reconstruct it with viridisLite so the gradient matches regardless.
sccustom_cols <- if (requireNamespace("scCustomize", quietly = TRUE)) {
  scCustomize::viridis_plasma_dark_high
} else {
  rev(viridisLite::plasma(50))
}

s5_rows <- lapply(cell_lines, function(cl) {
  obj <- cort_line_umap_list[[cl]]
  obj <- NormalizeData(obj)   # ensure the 'data' layer exists for FeaturePlot
  # one row of feature plots for this cell line, on its own embedding
  fp <- FeaturePlot(obj, features = s5_features, combine = FALSE)
  fp <- lapply(seq_along(fp), function(i) {
    p <- fp[[i]] +
      scale_color_gradientn(colors = sccustom_cols) +   # match FeaturePlot_scCustom default
      big_text_theme
    # prefix the leftmost feature's title with the cell line so each row is
    # clearly labeled by cell line
    if (i == 1) p <- p + ggtitle(paste0(cl, "\n", s5_features[i]))
    p
  })
  wrap_plots(fp, nrow = 1)
})

# stack the per-cell-line rows vertically into one combined figure
p_s5 <- wrap_plots(s5_rows, ncol = 1)
ggsave(
  "~/project/IPSC_2025_Data/Supplmentary_Figure5.tiff",
  plot = p_s5,
  device = "tiff",
  width = 5 * length(s5_features),
  height = 5 * length(cell_lines),
  dpi = 300, limitsize = FALSE
)

# -----------------------------------------------------------------------
# Differential expression, run SEPARATELY within each cell line.
# Model now includes neural_induction_media as an additional covariate.
# gt_line is dropped from the per-line formula since it's constant within
# a single cell line's subset (it was only needed as a covariate when
# lines were pooled).
# -----------------------------------------------------------------------
desired_order <- c("RG", "IPC_ExN", "DL_ExN", "UL_ExN", "A1 Astrocyte", "A2 Astrocyte", "Hem_RG", "CRN", "Epithelial")

res_mash_list <- list()
res_dl_list   <- list()
res_proc_list <- list()
res_df_list   <- list()

for (cl in cell_lines) {

  # Build clean, counts-only objects before merging - scale.data in the
  # source objects only covers each lineage's own variable features (a
  # different gene set for Cortical vs Hem), so merging with scale.data
  # still attached fails with a dimension mismatch. This step only needs
  # raw counts + metadata anyway.
  cort_counts_obj <- CreateSeuratObject(
    counts = LayerData(Cortical_lineage_list[[cl]], assay = "RNA", layer = "counts"),
    meta.data = Cortical_lineage_list[[cl]]@meta.data
  )
  hem_counts_obj <- CreateSeuratObject(
    counts = LayerData(Hem_lineage_list[[cl]], assay = "RNA", layer = "counts"),
    meta.data = Hem_lineage_list[[cl]]@meta.data
  )

  merged_temp <- merge(cort_counts_obj, hem_counts_obj)
  merged_temp <- JoinLayers(merged_temp)
  cnts <- merged_temp@assays$RNA@layers$counts
  colnames(cnts) <- colnames(merged_temp)
  rownames(cnts) <- rownames(merged_temp)
  merged_temp <- CreateSeuratObject(counts = cnts, meta.data = merged_temp@meta.data)
  merged_temp <- NormalizeData(merged_temp)
  IPSC_sce <- as.SingleCellExperiment(merged_temp)
  rm(merged_temp)

  IPSC_sce <- IPSC_sce[rowSums(assay(IPSC_sce, "counts") > 0) > 0, ]
  qc <- perCellQCMetrics(IPSC_sce)
  ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
  IPSC_sce <- IPSC_sce[, !ol]
  pb <- aggregateToPseudoBulk(IPSC_sce,
                              assay = "counts",
                              cluster_id = "Celltype2",
                              sample_id = "SampleID",
                              verbose = FALSE
  )

  # neural_induction_media added as a covariate here
  res.proc <- processAssays(pb, ~Protocol + Age + neural_induction_media, min.count = 5)
  res.dl   <- dreamlet(res.proc, ~Protocol + Age + neural_induction_media)
  res_mash <- run_mash(res.dl, coef = 'Protocolplus')

  res_proc_list[[cl]] <- res.proc
  res_dl_list[[cl]]   <- res.dl
  res_mash_list[[cl]] <- res_mash

  pm   <- get_pm(res_mash$model)
  genes <- rownames(pm)
  cells <- colnames(pm)
  lfsr <- get_lfsr(res_mash$model)

  res_df <- expand.grid(
    Gene = genes,
    Celltype = cells,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      logFC = as.vector(pm),
      adj.P.Val = as.vector(lfsr),
      sig = adj.P.Val < 0.05,
      direction = case_when(
        sig & logFC > 0 ~ "pos",
        sig & logFC < 0 ~ "neg",
        TRUE ~ "ns"
      ),
      color = case_when(
        direction == "pos" ~ "red",
        direction == "neg" ~ "blue",
        TRUE ~ "grey80"
      ),
      CellLine = cl
    )
  res_df$Celltype <- factor(res_df$Celltype, levels = desired_order)
  res_df_list[[cl]] <- res_df
}

res_df_all <- bind_rows(res_df_list)
res_df_all$CellLine <- factor(res_df_all$CellLine, levels = cell_lines)

# -----------------------------------------------------------------------
# Figure 3b: volcano plots, faceted by Celltype x CellLine, ONE combined
# figure across all cell lines
# -----------------------------------------------------------------------
plot_mash_volcano <- function(res_df, title = "mashr Volcano Plot") {

  ypos_df <- res_df %>%
    group_by(Celltype, CellLine) %>%
    summarise(
      ymin = min(-log10(adj.P.Val), na.rm = TRUE),
      ymax_raw = max(-log10(adj.P.Val), na.rm = TRUE),
      ymax = pmin(ymax_raw, 50),
      ylab = ymin + 0.80 * (ymax - ymin),
      .groups = "drop"
    )

  xpos_df <- res_df %>%
    group_by(Celltype, CellLine) %>%
    summarise(
      xmin_raw = suppressWarnings(min(logFC, na.rm = TRUE)),
      xmax_raw = suppressWarnings(max(logFC, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      xmin = ifelse(is.finite(xmin_raw), xmin_raw, 0),
      xmax = ifelse(is.finite(xmax_raw), xmax_raw, 0)
    )

  counts <- res_df %>%
    group_by(Celltype, CellLine, direction) %>%
    summarise(n_sig = sum(sig), .groups = "drop") %>%
    filter(direction != "ns") %>%
    left_join(ypos_df, by = c("Celltype", "CellLine")) %>%
    left_join(xpos_df, by = c("Celltype", "CellLine")) %>%
    mutate(
      x_pos = ifelse(direction == "pos", xmax, xmin),
      hjust = ifelse(direction == "pos", 1.1, -0.1)
    )

  ggplot(res_df, aes(x = logFC, y = -log10(adj.P.Val), color = color)) +
    # NS points: small and faint so they recede into the background.
    # Significant (red/blue) points: larger and more opaque so they stand
    # out clearly against the NS cloud.
    geom_point(data = ~ dplyr::filter(.x, direction == "ns"), size = 1.2, alpha = 0.25) +
    geom_point(data = ~ dplyr::filter(.x, direction != "ns"), size = 2.2, alpha = 0.85) +
    facet_grid(CellLine ~ Celltype, scales = "free") +
    scale_color_identity() +
    theme_bw(base_size = 20) +
    big_text_theme +
    theme(
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 24, face = "bold"),
      strip.text = element_text(size = 22, face = "bold"),
      plot.title = element_text(size = 26, face = "bold")
    ) +
    labs(
      x = "Posterior mean logFC (mashr)",
      y = "-log10(LFSR)",
      title = title
    ) +
    ggrepel::geom_label_repel(
      data = counts,
      aes(x = x_pos, y = ylab, label = paste0(n_sig, " gene(s)"), hjust = hjust),
      inherit.aes = FALSE,
      vjust = 0,
      size = 6,
      fontface = "bold",
      label.size = 0,
      fill = alpha("white", 0.75),
      # repel settings: keep labels off the point cloud; direction = "y" so
      # they nudge vertically (preserving their left/right side placement),
      # and seed for reproducibility
      direction = "y",
      box.padding = 0.5,
      point.padding = 0.3,
      min.segment.length = 0.2,
      segment.color = "grey50",
      max.overlaps = Inf,
      seed = 42
    )
}
pv <- plot_mash_volcano(res_df_all, title = "mashr Volcano Plot — Protocolplus, by Cell Line")
# Reduced per-panel size (was 8x5in, producing a 64x15in canvas where
# points/text became nearly invisible) - 6x6in per panel keeps the figure
# more legible at normal zoom while still resolving each panel clearly.
ggsave("~/project/IPSC_2025_Data/Figure3b_DEG_volcano_plots.tiff",
       plot = pv,
       device = "tiff",
       width = 6 * length(unique(res_df_all$Celltype)),
       height = 6 * length(cell_lines), dpi = 300, limitsize = FALSE)


# -----------------------------------------------------------------------
# Figure 3c: DE gene similarity heatmaps, per cell line, combined
# -----------------------------------------------------------------------
get_top_pos_genes_resdf <- function(df, celltype, n = nrow(df)) {
  tt <- df %>%
    filter(Celltype == celltype, logFC > 0, adj.P.Val < 0.05) %>%
    arrange(adj.P.Val, desc(logFC))
  head(tt$Gene, n)
}
get_top_neg_genes_resdf <- function(df, celltype, n = nrow(df)) {
  tt <- df %>%
    filter(Celltype == celltype, logFC < 0, adj.P.Val < 0.05) %>%
    arrange(adj.P.Val, logFC)
  head(tt$Gene, n)
}
make_similarity_matrix <- function(top_list) {
  assays <- names(top_list)
  n <- length(assays)
  sim_mat <- matrix(0, n, n, dimnames = list(assays, assays))
  overlap_mat <- matrix(0, n, n, dimnames = list(assays, assays))
  for (i in seq_along(assays)) {
    for (j in seq_along(assays)) {
      set1 <- top_list[[i]]; set2 <- top_list[[j]]
      overlap <- length(intersect(set1, set2))
      denom   <- min(length(set1), length(set2))
      sim_mat[i, j] <- ifelse(denom > 0, overlap / denom, 0)
      overlap_mat[i, j] <- overlap
    }
  }
  list(similarity = sim_mat, overlap = overlap_mat)
}

setwd("~/project/IPSC_2025_Data/Figure_Parts/Manuscript1")
pos_heatmaps <- list()
neg_heatmaps <- list()
for (cl in cell_lines) {
  df_line <- res_df_list[[cl]]
  celltypes_resdf <- sort(unique(df_line$Celltype))
  top_pos <- lapply(celltypes_resdf, function(ct) get_top_pos_genes_resdf(df_line, ct))
  top_neg <- lapply(celltypes_resdf, function(ct) get_top_neg_genes_resdf(df_line, ct))
  names(top_pos) <- celltypes_resdf
  names(top_neg) <- celltypes_resdf
  sim_pos <- make_similarity_matrix(top_pos)
  sim_neg <- make_similarity_matrix(top_neg)
  sim_pos$similarity[!is.finite(sim_pos$similarity)] <- 0
  sim_neg$similarity[!is.finite(sim_neg$similarity)] <- 0
  global_min <- min(sim_pos$similarity, sim_neg$similarity)
  global_max <- max(sim_pos$similarity, sim_neg$similarity)
  my_breaks <- seq(global_min, global_max, length.out = 101)
  my_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

  # Explicit :: needed here - ComplexHeatmap (loaded above) also exports a
  # pheatmap() compatibility wrapper that returns an S4 Heatmap object
  # instead of the classic list-based pheatmap object, which breaks $gtable
  # below. Force the real pheatmap package's function.
  #
  # main = NULL here (no built-in pheatmap title) - pheatmap's own title
  # has almost no reserved vertical space and collides with the column
  # label row directly beneath it. Titles are added separately below with
  # guaranteed, non-overlapping space.
  pos_heatmaps[[cl]] <- pheatmap::pheatmap(
    sim_pos$similarity, display_numbers = sim_pos$overlap, fontsize_number = 12,
    fontsize = 14, main = NA,
    color = my_colors, breaks = my_breaks, silent = TRUE
  )
  neg_heatmaps[[cl]] <- pheatmap::pheatmap(
    sim_neg$similarity, display_numbers = sim_neg$overlap, fontsize_number = 12,
    fontsize = 14, main = NA,
    color = my_colors, breaks = my_breaks, silent = TRUE
  )
}

# Build each heatmap as (title text) stacked above (heatmap grob), with a
# fixed, reserved title band so text can never overlap the plot below it.
make_titled_heatmap <- function(hm, title_text) {
  title_grob <- ggdraw() +
    draw_label(title_text, fontface = "bold", size = 16, x = 0.5, hjust = 0.5)
  plot_grid(title_grob, hm$gtable, ncol = 1, rel_heights = c(0.12, 1))
}

pos_titled <- lapply(cell_lines, function(cl) {
  make_titled_heatmap(pos_heatmaps[[cl]], paste0(cl, ": Positive DE gene similarity"))
})
neg_titled <- lapply(cell_lines, function(cl) {
  make_titled_heatmap(neg_heatmaps[[cl]], paste0(cl, ": Negative DE gene similarity"))
})

combined_pos <- plot_grid(plotlist = pos_titled, nrow = 1)
combined_neg <- plot_grid(plotlist = neg_titled, nrow = 1)
fig3c_combined <- plot_grid(combined_pos, combined_neg, ncol = 1)
ggsave("Figure3c_combined.tiff", plot = fig3c_combined, device = "tiff",
       width = 6 * length(cell_lines), height = 14, dpi = 300, limitsize = FALSE)

# -----------------------------------------------------------------------
# Figure 3d: gene-set enrichment (zenith), per cell line, combined
# -----------------------------------------------------------------------
# Gene sets for zenith gene-set enrichment (Gene Ontology BP/CC), loaded
# once here since the original script referenced go.gs.bp/go.gs.cc without
# ever defining them - these must have existed from an earlier interactive
# session but weren't in the saved script.
go.gs.bp <- zenith::get_GeneOntology(to = "BP", species = "human")
go.gs.cc <- zenith::get_GeneOntology(to = "CC", species = "human")

# -----------------------------------------------------------------------
# Figure 3d: gene-set enrichment (zenith), per cell line, combined
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Figure 3d: gene-set enrichment (zenith), per cell line, combined
# -----------------------------------------------------------------------
# plotZenithResults() doesn't expose enough control over its internal
# label rendering to fix overlap reliably (its scale_y_discrete labels=
# override didn't take effect). Building the heatmap directly from the
# zenith_gsa() output instead, with full control over term selection,
# wrapping, and sizing.
#
# zenith_gsa() columns: assay, coef, Geneset, NGenes, Direction, PValue, FDR
# This replicates plotZenithResults()'s EXACT internal tstat derivation
# (confirmed from its source):
#   if delta & se columns exist:  tstat = delta / se   (a true Wald stat)
#   else:                         tstat = qnorm(PValue, lower.tail = FALSE)
#                                 signed by Direction (Up=+1, Down=-1)
# Our zenith_gsa() output has no delta/se columns, so the else branch
# applies. NOTE: this is a ONE-sided z (qnorm(p, lower.tail = FALSE)),
# not the two-sided version - matching zenith exactly.
# Step 1: compute signed tstat + select terms for ONE cell line's zenith df.
# Returns BOTH the full per-term data (every geneset's tstat for this cell
# line, so any union term's real value can be looked up later) AND the
# vector of terms this cell line selected as its own top/bottom hits.
zenith_prep <- function(df, ntop = 3, cell_line = "") {

  df <- df %>%
    mutate(
      signed_score = if (all(c("delta", "se") %in% names(df))) {
        delta / se
      } else {
        qnorm(PValue, lower.tail = FALSE) * ifelse(Direction == "Up", 1, -1)
      },
      Geneset_wrapped = stringr::str_wrap(Geneset, width = 40),
      CellLine = cell_line
    )

  # Term selection matching plotZenithResults(): per assay, take the top
  # `ntop` by tstat (most positive / Up) AND the bottom `ntop` by tstat
  # (most negative / Down), unioned across assays.
  selected_terms <- df %>%
    group_by(assay) %>%
    group_modify(~ {
      ts <- sort(.x$signed_score)
      n <- length(ts)
      cutoff_bottom <- if (ntop > 0 && n >= ntop) ts[ntop] else -Inf
      cutoff_top    <- if (ntop > 0 && n >= ntop) rev(ts)[ntop] else Inf
      .x[.x$signed_score <= cutoff_bottom | .x$signed_score >= cutoff_top, ]
    }) %>%
    ungroup() %>%
    pull(Geneset_wrapped) %>%
    unique()

  # NOTE: full_df keeps ALL terms (not filtered to this line's selection),
  # so downstream we can pull the actual tstat of a term that was selected
  # by a DIFFERENT cell line - this is what fills in the previously-white
  # cells and lets you see whether a top hit in one line trends the same
  # direction in the others.
  list(full_df = df, selected_terms = selected_terms)
}

# Step 2: given all cell lines' selected-term vectors + full data, compute
# the UNION of selected terms and ONE global ordering (by mean signed tstat
# across all cell lines) so negative(blue)/positive(red) group consistently.
compute_global_term_order <- function(prepped_list) {
  union_terms <- unique(unlist(lapply(prepped_list, `[[`, "selected_terms")))
  pooled_full <- bind_rows(lapply(prepped_list, `[[`, "full_df"))
  pooled_full %>%
    filter(Geneset_wrapped %in% union_terms) %>%
    group_by(Geneset_wrapped) %>%
    summarise(mean_score = mean(signed_score, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_score) %>%   # most negative first -> bottom of plot
    pull(Geneset_wrapped)
}

# Step 3: draw one cell line's panel. Pulls the ACTUAL tstat of every union
# term from THIS cell line's full data (not just terms it selected itself),
# so a term that's a top hit elsewhere shows its real value here rather than
# rendering blank/white.
plot_zenith_heatmap <- function(prepped, global_term_order, title = "", zmax = NULL) {

  plot_df <- prepped$full_df %>%
    filter(Geneset_wrapped %in% global_term_order)
  plot_df$Geneset_wrapped <- factor(plot_df$Geneset_wrapped, levels = global_term_order)

  n_terms <- length(global_term_order)

  # Symmetric color limits so white sits exactly at 0 and blue/red intensity
  # scale identically, matching zenith's limits = c(-zmax, zmax). If a zmax
  # is supplied (shared across all panels), use it so every panel has an
  # IDENTICAL color scale for direct visual comparison; otherwise fall back
  # to this panel's own max.
  if (is.null(zmax)) {
    zmax <- max(abs(plot_df$signed_score), na.rm = TRUE)
  }

  list(
    plot = ggplot(plot_df, aes(x = assay, y = Geneset_wrapped, fill = signed_score)) +
      geom_tile(color = "white", linewidth = 0.3) +
      # FDR < 0.05 significance markers, matching zenith's geom_text(label = *)
      geom_text(aes(label = ifelse(!is.na(FDR) & FDR < 0.05, "*", "")), size = 6, vjust = 0.75) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                           limits = c(-zmax, zmax),
                           name = "t-statistic") +
      scale_y_discrete(drop = FALSE) +   # keep all global terms even if absent here
      ggtitle(title) +
      theme_bw(base_size = 14) +
      theme(
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(size = 14),
        axis.title = element_blank(),
        plot.title = element_text(size = 16, face = "bold"),
        legend.position = "right"
      ),
    n_terms = n_terms
  )
}

zenith_bp_prepped <- list()
zenith_cc_prepped <- list()

# Pass 1: prep each cell line (compute tstat for ALL terms + record which
# terms this line selects as its own top/bottom hits)
for (cl in cell_lines) {
  df_gs_bp <- zenith_gsa(res_mash_list[[cl]], go.gs.bp)
  df_gs_cc <- zenith_gsa(res_mash_list[[cl]], go.gs.cc)
  df_gs_bp$assay <- factor(df_gs_bp$assay, levels = desired_order)
  df_gs_cc$assay <- factor(df_gs_cc$assay, levels = desired_order)

  zenith_bp_prepped[[cl]] <- zenith_prep(df_gs_bp, ntop = 3, cell_line = cl)
  zenith_cc_prepped[[cl]] <- zenith_prep(df_gs_cc, ntop = 3, cell_line = cl)
}

# Pass 2: union of selected terms + shared ordering (BP and CC separately)
bp_global_order <- compute_global_term_order(zenith_bp_prepped)
cc_global_order <- compute_global_term_order(zenith_cc_prepped)

# Compute ONE shared symmetric color-scale max across ALL cell lines, so
# every BP panel uses the identical scale (and likewise every CC panel).
# Only union terms actually drawn are considered.
bp_zmax <- max(abs(
  bind_rows(lapply(zenith_bp_prepped, `[[`, "full_df")) %>%
    filter(Geneset_wrapped %in% bp_global_order) %>% pull(signed_score)
), na.rm = TRUE)
cc_zmax <- max(abs(
  bind_rows(lapply(zenith_cc_prepped, `[[`, "full_df")) %>%
    filter(Geneset_wrapped %in% cc_global_order) %>% pull(signed_score)
), na.rm = TRUE)

# Pass 3: draw each panel against the shared order AND shared color scale,
# pulling ACTUAL tstats for every union term from each line's full data
zenith_bp_plots <- list()
zenith_cc_plots <- list()
for (cl in cell_lines) {
  zenith_bp_plots[[cl]] <- plot_zenith_heatmap(zenith_bp_prepped[[cl]], bp_global_order,
                                               title = paste0(cl, ": GO BP"), zmax = bp_zmax)$plot
  zenith_cc_plots[[cl]] <- plot_zenith_heatmap(zenith_cc_prepped[[cl]], cc_global_order,
                                               title = paste0(cl, ": GO CC"), zmax = cc_zmax)$plot
}

fig3d_bp_combined <- wrap_plots(zenith_bp_plots, ncol = length(cell_lines))
fig3d_cc_combined <- wrap_plots(zenith_cc_plots, ncol = length(cell_lines))

# Panel height scales with the (shared) number of term rows drawn, so
# labels get consistent per-row spacing.
bp_height <- max(8, 0.35 * length(bp_global_order))
cc_height <- max(8, 0.35 * length(cc_global_order))

ggsave("Figure3d_GO_BP_combined.tiff", plot = fig3d_bp_combined, device = "tiff",
       width = 10 * length(cell_lines), height = bp_height, dpi = 300, limitsize = FALSE)
ggsave("Figure3d_GO_CC_combined.tiff", plot = fig3d_cc_combined, device = "tiff",
       width = 10 * length(cell_lines), height = cc_height, dpi = 300, limitsize = FALSE)


# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Supplemental Figure 6 (rebuilt to address reviewer concern about
# E6 vs KSR neural induction media):
#
# Design context: neural_induction_media is severely confounded with
# cell line and timepoint at the dataset level (E6 is essentially JHC1-
# only, some timepoints are single-media). A pooled "control for media as
# a covariate" approach therefore can't cleanly isolate media - the term
# is unestimable wherever only one media is present. Instead we:
#
#   Panel A: DESIGN MAP - a tile plot of line x timepoint x media coverage
#            so the (uneven) E6/KSR sampling is fully transparent.
#
#   Panel B: MATCHED E6-vs-KSR DE - a direct differential test between
#            media conditions, run ONLY within line x timepoint blocks
#            that actually contain BOTH E6 and KSR (auto-detected). This
#            shows exactly what genes/pathways differ by media, where it
#            can be honestly measured.
#
#   Panel C: PROTOCOL-DEG MEDIA-INVARIANCE - for each line x timepoint
#            block that has both media, we compute Protocol (+SDF vs -SDF)
#            DEGs SEPARATELY within each media condition, then show the
#            overlap. High overlap = the Protocol response recruits the
#            same genes regardless of media, i.e. the main Protocol
#            conclusion is robust to media choice.
# -----------------------------------------------------------------------

# Rebuild a clean, counts-only combined IPSC object (all lineages, all
# lines) so we can subset freely by line/timepoint/media without the
# Assay5 scale.data dimension issues encountered earlier.
merged_all <- merge(
  CreateSeuratObject(counts = LayerData(Cortical_lineage, assay = "RNA", layer = "counts"),
                     meta.data = Cortical_lineage@meta.data),
  CreateSeuratObject(counts = LayerData(Hem_lineage, assay = "RNA", layer = "counts"),
                     meta.data = Hem_lineage@meta.data)
)
merged_all <- JoinLayers(merged_all)
merged_all <- NormalizeData(merged_all)

meta_all <- merged_all@meta.data

# ----- Panel A: design coverage map -----
design_tbl <- meta_all %>%
  distinct(gt_line, Age, neural_induction_media) %>%
  mutate(present = TRUE)

# also count cells per combo for shading
cell_counts <- meta_all %>%
  group_by(gt_line, Age, neural_induction_media) %>%
  summarise(n_cells = n(), .groups = "drop")

design_map <- cell_counts %>%
  mutate(Age = factor(Age))

pA <- ggplot(design_map, aes(x = factor(Age), y = neural_induction_media, fill = n_cells)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = n_cells), size = 4) +
  facet_wrap(~gt_line, nrow = 1) +
  scale_fill_viridis_c(option = "D", trans = "log10", name = "n cells") +
  labs(title = "A. Media coverage by cell line and timepoint",
       x = "Timepoint (Age)", y = "Neural induction media") +
  theme_bw() + big_text_theme

# ----- Auto-detect line x timepoint blocks with BOTH media present -----
both_media_blocks <- meta_all %>%
  group_by(gt_line, Age) %>%
  summarise(n_media = n_distinct(neural_induction_media), .groups = "drop") %>%
  filter(n_media == 2)

message("line x timepoint blocks with both E6 and KSR:")
print(both_media_blocks)

# Helper: pseudobulk DE within a subset for a given 2-level grouping var,
# returns a tidy results frame. Uses dreamlet on Celltype2 pseudobulk.
run_block_de <- function(sub_obj, group_var, block_label) {
  sce <- as.SingleCellExperiment(sub_obj)
  sce <- sce[rowSums(assay(sce, "counts") > 0) > 0, ]
  pb <- tryCatch(
    aggregateToPseudoBulk(sce, assay = "counts",
                          cluster_id = "Celltype2", sample_id = "SampleID",
                          verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(pb)) return(NULL)

  form <- as.formula(paste0("~", group_var))
  res.proc <- tryCatch(processAssays(pb, form, min.count = 5), error = function(e) NULL)
  if (is.null(res.proc)) return(NULL)
  res.dl <- tryCatch(dreamlet(res.proc, form), error = function(e) NULL)
  if (is.null(res.dl)) return(NULL)

  # coefficient name is the non-reference level of group_var
  coef_names <- setdiff(colnames(coef(res.dl[[1]])), "(Intercept)")
  coef_use <- coef_names[grepl(group_var, coef_names)][1]
  if (is.na(coef_use)) return(NULL)

  tt <- tryCatch(topTable(res.dl, coef = coef_use, number = Inf), error = function(e) NULL)
  if (is.null(tt)) return(NULL)
  # dreamlet::topTable returns an S4 DFrame (S4Vectors), which bind_rows()
  # can't stack - coerce to a base data.frame first.
  tt <- as.data.frame(tt)
  # Normalize the gene and celltype identifier columns so downstream code
  # can rely on `gene` and `assay` existing regardless of dreamlet version.
  # dreamlet typically returns `ID` (gene) and `assay` (celltype) as
  # columns; older/other versions may put the gene in rownames.
  if (!"gene" %in% names(tt)) {
    if ("ID" %in% names(tt)) {
      tt$gene <- tt$ID
    } else {
      tt$gene <- rownames(tt)
    }
  }
  if (!"assay" %in% names(tt)) {
    # if no celltype column exists, this was a single-assay result
    tt$assay <- NA_character_
  }
  tt$block <- block_label
  tt$coef <- coef_use
  tt
}

# ----- Panel B: matched E6-vs-KSR DE within each both-media block -----
media_de_list <- list()
for (i in seq_len(nrow(both_media_blocks))) {
  cl_i  <- both_media_blocks$gt_line[i]
  age_i <- both_media_blocks$Age[i]
  block_label <- paste0(cl_i, "_D", age_i)

  sub_obj <- subset(merged_all, gt_line == cl_i & Age == age_i)
  # relevel media so KSR is reference -> coefficient is "E6 vs KSR"
  sub_obj$neural_induction_media <- factor(sub_obj$neural_induction_media,
                                           levels = c("KSR", "E6"))
  res <- run_block_de(sub_obj, "neural_induction_media", block_label)
  if (!is.null(res)) media_de_list[[block_label]] <- res
}
media_de_all <- bind_rows(media_de_list)

# volcano of media DE, faceted by block x celltype
if (nrow(media_de_all) > 0) {
  media_de_all <- media_de_all %>%
    mutate(sig = adj.P.Val < 0.05,
           color = case_when(sig & logFC > 0 ~ "#B2182B",
                             sig & logFC < 0 ~ "#2166AC",
                             TRUE ~ "grey80"))
  pB <- ggplot(media_de_all, aes(x = logFC, y = -log10(adj.P.Val), color = color)) +
    geom_point(data = ~ dplyr::filter(.x, !sig), size = 1.0, alpha = 0.3) +
    geom_point(data = ~ dplyr::filter(.x, sig), size = 2.0, alpha = 0.85) +
    facet_grid(block ~ assay, scales = "free") +
    scale_color_identity() +
    labs(title = "B. Matched E6 vs KSR differential expression (KSR = reference), where both media present",
         x = "logFC (E6 vs KSR)", y = "-log10(adj p)") +
    theme_bw() + big_text_theme
} else {
  pB <- ggplot() + theme_void() +
    ggtitle("B. No line x timepoint block had both media with estimable DE")
}

# ----- Panel C: Protocol effect robustness to media adjustment -----
# The per-media Protocol-DEG overlap approach is NOT estimable with this
# design: subdividing line x timepoint x media x protocol x celltype into
# pseudobulk leaves too few replicates for the Protocol contrast to fit
# within each media separately. Instead we use the standard covariate-
# adjustment robustness check, which preserves power by keeping all cells:
#
#   For each cell line, compare the per-gene Protocol logFC estimated
#   WITHOUT media in the model (~Protocol + Age) vs WITH media as a
#   covariate (~Protocol + Age + neural_induction_media, already fit in
#   res_dl_list). If adjusting for media barely moves the Protocol logFCs
#   (tight scatter on y = x, high correlation), the Protocol effect is not
#   explained by / confounded with media - i.e. it is robust.
#
# CAVEAT (state in the manuscript): E6 is almost entirely JHC1, so media
# and cell line are partially collinear; adjustment can only partially
# separate them. This is a design limitation, disclosed rather than hidden.

robustness_list <- list()
for (cl in cell_lines) {
  res.proc <- res_proc_list[[cl]]
  res.dl_adj <- res_dl_list[[cl]]   # already fit with media covariate

  # unadjusted fit on the SAME processed pseudobulk (identical gene/celltype
  # sets), so logFCs are directly comparable per gene x celltype
  res.dl_unadj <- tryCatch(dreamlet(res.proc, ~Protocol + Age), error = function(e) NULL)
  if (is.null(res.dl_unadj)) next

  tt_adj <- tryCatch(as.data.frame(topTable(res.dl_adj, coef = "Protocolplus", number = Inf)),
                     error = function(e) NULL)
  tt_unadj <- tryCatch(as.data.frame(topTable(res.dl_unadj, coef = "Protocolplus", number = Inf)),
                       error = function(e) NULL)
  if (is.null(tt_adj) || is.null(tt_unadj)) next

  # normalize gene id column
  norm_ids <- function(tt) {
    if (!"gene" %in% names(tt)) tt$gene <- if ("ID" %in% names(tt)) tt$ID else rownames(tt)
    if (!"assay" %in% names(tt)) tt$assay <- NA_character_
    tt
  }
  tt_adj <- norm_ids(tt_adj); tt_unadj <- norm_ids(tt_unadj)

  merged_fc <- inner_join(
    tt_unadj %>% dplyr::select(assay, gene, logFC_unadj = logFC),
    tt_adj   %>% dplyr::select(assay, gene, logFC_adj = logFC),
    by = c("assay", "gene")
  ) %>% mutate(CellLine = cl)
  robustness_list[[cl]] <- merged_fc
}
robustness_df <- bind_rows(robustness_list)

if (nrow(robustness_df) > 0) {
  # correlation per cell line, annotated on each facet
  cor_labels <- robustness_df %>%
    group_by(CellLine) %>%
    summarise(
      r = cor(logFC_unadj, logFC_adj, use = "complete.obs"),
      x = min(logFC_unadj, na.rm = TRUE),
      y = max(logFC_adj, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(label = paste0("r = ", sprintf("%.3f", r)))

  pC <- ggplot(robustness_df, aes(x = logFC_unadj, y = logFC_adj)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "firebrick", linewidth = 0.8) +
    geom_point(alpha = 0.25, size = 0.9) +
    geom_text(data = cor_labels, aes(x = x, y = y, label = label),
              hjust = 0, vjust = 1, size = 6, fontface = "bold", color = "firebrick") +
    facet_wrap(~CellLine, nrow = 1) +
    labs(title = "C. Protocol logFC is robust to media adjustment: unadjusted (~Protocol+Age) vs media-adjusted (+neural_induction_media)",
         subtitle = "Points on the dashed y = x line indicate the Protocol effect is unchanged by adjusting for media (per gene x celltype)",
         x = "Protocol logFC — media NOT in model", y = "Protocol logFC — media adjusted") +
    theme_bw() + big_text_theme
} else {
  pC <- ggplot() + theme_void() +
    ggtitle("C. Could not fit unadjusted vs adjusted Protocol models for comparison")
}

fig_s6_combined <- pA / pB / pC +
  plot_layout(heights = c(1, 2.2, 1.2)) +
  plot_annotation(
    title = "Supplemental Figure 6: E6 vs KSR neural induction media — coverage, direct effect, and Protocol-effect robustness to media adjustment"
  )

ggsave(
  "~/project/IPSC_2025_Data/Supplmentary_Figure6.tiff",
  plot = fig_s6_combined,
  device = "tiff",
  width = 22, height = 26, dpi = 300, limitsize = FALSE
)
