library(SeuratExtend)
library(Seurat)
library(SeuratDisk)    
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggtext)
library(pheatmap)
library(purrr)
library(stringr)
library(ComplexHeatmap)   
library(Matrix)
library(patchwork)
library(cowplot)
library(magick)

big_text_theme <- theme(
  axis.text = element_text(size = 14),
  axis.title = element_text(size = 16, face = "bold"),
  strip.text = element_text(size = 16, face = "bold"),
  legend.text = element_text(size = 13),
  legend.title = element_text(size = 14, face = "bold"),
  plot.title = element_text(size = 18, face = "bold")
)

# Canonical lineage labels, used consistently across Fig 5A/B/C (matches
# the "<CellLine>: A1 lineage eigengene" style already used for Panel B)
lineage_labels <- c(
  dp  = "DP ExN lineage",
  up  = "UP ExN lineage",
  A1  = "A1 Astrocyte lineage",
  A2  = "A2 Astrocyte lineage",
  crn = "CRN lineage",
  epi = "Epithelial lineage"
)
lineage_order <- c("dp", "up", "A1", "A2", "epi", "crn")
lineage_order_labeled <- unname(lineage_labels[lineage_order])

# -----------------------------------------------------------------------
# ONE-TIME SETUP STEP - already completed. This block exports data for
# pySCENIC to run on in Python and produces merged_IPSC_aucell.loom (used
# by the block below). Since that loom file already exists and pySCENIC
# is NOT being re-run, this block is disabled to avoid regenerating an
# export nobody needs (this was the source of the SaveH5Seurat stall -
# SeuratDisk's SaveH5Seurat has known compatibility issues/hangs with
# Seurat v5 Assay5 objects, which is an additional reason not to re-run
# this unless you specifically need to redo the pySCENIC step itself).
#
# Only re-enable (set to TRUE) if you need to regenerate the h5ad export
# to re-run pySCENIC from scratch.
# -----------------------------------------------------------------------
RUN_PYSCENIC_EXPORT <- FALSE

if (RUN_PYSCENIC_EXPORT) {
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
  # Seurat v5 Assay5 objects store data in layers, not the old @counts slot -
  # use LayerData() instead of @assays$RNA@counts.
  cnts <- LayerData(merged_IPSC, assay = "RNA", layer = "counts")
  colnames(cnts) <- colnames(merged_IPSC)
  rownames(cnts) <- rownames(merged_IPSC)
  merged_IPSC <- CreateSeuratObject(counts = cnts, meta.data = merged_IPSC@meta.data, min.cells = 5)
  merged_IPSC[["RNA"]] <- as(object = merged_IPSC[["RNA"]], Class = "Assay")
  SaveH5Seurat(merged_IPSC, filename = "~/project/IPSC_2025_Data/merged_IPSC.h5Seurat")
  Convert("~/project/IPSC_2025_Data/merged_IPSC.h5Seurat", dest = "h5ad")
}


#After Runing pySCENIC (unchanged - SCENIC is not being re-run, per your
#note; we only re-partition the results per cell line downstream)
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
scenic_loom_path <- "~/project/IPSC_2025_Data/merged_IPSC_aucell.loom"
merged_IPSC <- ImportPyscenicLoom(scenic_loom_path, seu = merged_IPSC)
dim(merged_IPSC)
dim(merged_IPSC@misc$SCENIC$RegulonsAUC)



Cortical_lineage_list <- readRDS("~/project/IPSC_2025_Data/checkpoint_pallial_modules_found_by_line")
Hem_lineage_list <- readRDS("~/project/IPSC_2025_Data/checkpoint_hem_modules_found_by_line")



# Transfer metadata separately PER CELL LINE (barcodes are unique dataset-
# wide, so matching by rowname still correctly routes each cell's own
# line-specific pseudotime/lineage values into merged_IPSC).
for (cl in cell_lines) {
  Cortical_lineage <- Cortical_lineage_list[[cl]]
  Hem_lineage <- Hem_lineage_list[[cl]]
  
  # Match on the ORIGINAL barcode stamped into metadata (orig_barcode) in
  # Figure_3.R BEFORE the subset+re-embed branches re-encoded cell names.
  # The Cortical and Hem lineage objects are built via separate
  # subset+merge branches that append DIFFERENT disambiguating suffixes to
  # the cell names, so their barcodes no longer match merged_IPSC's (the
  # Hem branch in particular diverged completely - 0 barcode overlap). A
  # metadata column survives subset()/merge() untouched, so orig_barcode
  # holds each cell's name AS IT EXISTS in merged_IPSC, giving a stable,
  # exact key. rownames(merged_IPSC) ARE those original barcodes, so we
  # match orig_barcode directly against colnames(merged_IPSC).
  stopifnot("orig_barcode" %in% colnames(Hem_lineage@meta.data),
            "orig_barcode" %in% colnames(Cortical_lineage@meta.data))
  
  meta_hem <- Hem_lineage@meta.data[, c("orig_barcode", "Celltype2", grep("pseudotime", colnames(Hem_lineage@meta.data), value = TRUE), grep("lineage", colnames(Hem_lineage@meta.data), value = TRUE))]
  meta_hem$pseudotime <- NULL
  
  meta_cor <- Cortical_lineage@meta.data[, c("orig_barcode", "Celltype2", grep("pseudotime", colnames(Cortical_lineage@meta.data), value = TRUE), grep("lineage", colnames(Cortical_lineage@meta.data), value = TRUE))]
  meta_cor$pseudotime <- NULL
  
  meta_hem[] <- lapply(meta_hem, function(x) if (is.factor(x)) as.character(x) else x)
  meta_cor[] <- lapply(meta_cor, function(x) if (is.factor(x)) as.character(x) else x)
  
  # A cell can appear in BOTH the Cortical and Hem branches (or be
  # duplicated within one), so an orig_barcode could recur - drop any
  # colliding orig_barcode rather than arbitrarily picking one.
  ob_hem <- meta_hem$orig_barcode
  ob_cor <- meta_cor$orig_barcode
  dup_hem <- ob_hem %in% ob_hem[duplicated(ob_hem)]
  dup_cor <- ob_cor %in% ob_cor[duplicated(ob_cor)]
  if (sum(dup_hem) > 0) cat(cl, ": dropping", sum(dup_hem), "Hem cells with duplicate orig_barcode\n")
  if (sum(dup_cor) > 0) cat(cl, ": dropping", sum(dup_cor), "Cortical cells with duplicate orig_barcode\n")
  meta_hem <- meta_hem[!dup_hem, , drop = FALSE]
  meta_cor <- meta_cor[!dup_cor, , drop = FALSE]
  
  rownames(meta_hem) <- meta_hem$orig_barcode; meta_hem$orig_barcode <- NULL
  rownames(meta_cor) <- meta_cor$orig_barcode; meta_cor$orig_barcode <- NULL
  
  common_hem <- intersect(rownames(meta_hem), colnames(merged_IPSC))
  merged_IPSC@meta.data[common_hem, colnames(meta_hem)] <- meta_hem[common_hem, ]
  
  common_cor <- intersect(rownames(meta_cor), colnames(merged_IPSC))
  merged_IPSC@meta.data[common_cor, colnames(meta_cor)] <- meta_cor[common_cor, ]
  
  cat(cl, "- Hem cells matched:", length(common_hem), "/", nrow(meta_hem),
      "| Cortical cells matched:", length(common_cor), "/", nrow(meta_cor), "\n")
}
table(merged_IPSC$Celltype2)

merged_IPSC <- subset(merged_IPSC, Celltype2 == "A1 Astrocyte" | Celltype2 == "A2 Astrocyte" | Celltype2 == "CRN" | Celltype2 == "DL_ExN" | Celltype2 == "UL_ExN" | Celltype2 == "IPC_ExN" | Celltype2 == "RG" | Celltype2 == "Epithelial" | Celltype2 == "Hem_RG")
merged_IPSC@misc$SCENIC$RegulonsAUC <- merged_IPSC@misc$SCENIC$RegulonsAUC[colnames(merged_IPSC),]
tf_auc <- merged_IPSC@misc$SCENIC$RegulonsAUC
dim(tf_auc)
dim(merged_IPSC)
identical(rownames(tf_auc), colnames(merged_IPSC))
tf_auc <- tf_auc %>% mutate(across(everything(), ~ replace_na(., 0)))
identical(rownames(tf_auc), colnames(merged_IPSC))
tf_counts <- t(tf_auc)
# Properly CREATE the TF assay (as Assay5, matching Seurat v5's default
# object class used throughout this pipeline) rather than assigning
# directly into a nonexistent @assays$TF@counts slot, which was never
# valid regardless of Seurat version - there was no TF assay object here
# to assign a slot into at all.
tf_assay <- CreateAssay5Object(counts = tf_counts, data = tf_counts)
merged_IPSC[["TF"]] <- tf_assay
DefaultAssay(merged_IPSC) <- 'RNA'
dim(merged_IPSC)
DefaultAssay(merged_IPSC) <- 'TF'
dim(merged_IPSC)
merged_IPSC <- ScaleData(merged_IPSC)
merged_IPSC_tf <- merged_IPSC
DefaultAssay(merged_IPSC_tf) <- "TF"
merged_IPSC_tf@assays$RNA <- NULL

# -----------------------------------------------------------------------
# Figure 5a: Regulon differential "expression" (AUC), kept SEPARATE per
# cell line rather than collapsed across gt_line (previous version grouped
# by interaction(Age, gt_line) internally but then summarised away gt_line
# entirely when collapsing to one row per Regulon x Lineage - i.e. results
# were pooled across cell lines). Now gt_line is retained as its own column
# all the way through, and Age is the only dimension collapsed within a
# cell line.
# -----------------------------------------------------------------------
lineage_df <- merged_IPSC_tf@meta.data %>%
  tibble::rownames_to_column("cell") %>%        
  dplyr::select(
    cell, dp_pseudotime, up_pseudotime, A1_pseudotime, A2_pseudotime,
    crn_pseudotime, epi_pseudotime, Age, gt_line, Protocol
  ) %>%
  tidyr::pivot_longer(
    cols = ends_with("pseudotime"), names_to = "Lineage", values_to = "pseudotime"
  ) %>%
  dplyr::filter(!is.na(pseudotime)) %>%
  dplyr::mutate(Lineage = sub("_pseudotime", "", Lineage))

lineages <- unique(lineage_df$Lineage)

merged_IPSC_tf$dummy_ident <- "all_cells"
Idents(merged_IPSC_tf) <- "dummy_ident"

der_lineage_list <- lapply(lineages, function(lin) {
  
  cells_lin <- lineage_df %>% dplyr::filter(Lineage == lin)
  
  # split by Age x gt_line, but keep gt_line as a column we group by later
  # (previously: interaction(Age, gt_line) was used only to define FindMarkers
  # groups, then collapsed away entirely across BOTH Age and gt_line)
  groups <- cells_lin %>%
    mutate(group = interaction(Age, gt_line, drop = TRUE)) %>%
    split(.$group)
  
  group_results <- lapply(groups, function(g) {
    cells_plus  <- g %>% dplyr::filter(Protocol == "plus")  %>% pull(cell)
    cells_minus <- g %>% dplyr::filter(Protocol == "minus") %>% pull(cell)
    if (length(cells_plus) < 3 || length(cells_minus) < 3) return(NULL)
    
    res <- FindMarkers(
      object = merged_IPSC_tf[["TF"]],
      cells.1 = cells_plus, cells.2 = cells_minus,
      logfc.threshold = 0, min.pct = 0.25, layer = "data"   # Seurat v5: layer, not slot
    )
    res$gt_line <- unique(g$gt_line)[1]
    res
  })
  
  group_results <- group_results[!sapply(group_results, is.null)]
  if (length(group_results) == 0) return(NULL)
  
  group_results <- lapply(group_results, function(df) {
    df <- as.data.frame(df)
    df$Regulon <- rownames(df)
    df
  })
  
  combined <- bind_rows(group_results, .id = "Group")
  
  # Collapse across Age ONLY - gt_line stays as its own grouping variable
  collapsed <- combined %>%
    group_by(Regulon, gt_line) %>%
    summarise(
      median_p_val_adj = median(p_val_adj, na.rm = TRUE),
      max_p_val_adj    = max(p_val_adj, na.rm = TRUE),
      mean_log2FC      = mean(avg_log2FC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Lineage = lin)
  
  collapsed
})

names(der_lineage_list) <- lineages
der_lineage_results <- der_lineage_list %>%
  discard(is.null) %>%
  bind_rows()

res_df <- der_lineage_results %>%
  mutate(
    negLog10P = -log10(median_p_val_adj),
    sig = median_p_val_adj <= 0.05,
    direction = case_when(
      sig & mean_log2FC > 0 ~ "Positive",
      sig & mean_log2FC < 0 ~ "Negative",
      TRUE ~ "NotSig"
    ),
    LineageLabel = factor(unname(lineage_labels[Lineage]), levels = lineage_order_labeled)
  )

counts_df <- res_df %>%
  dplyr::filter(sig) %>%
  group_by(LineageLabel, gt_line, direction) %>%
  summarise(n = n(),
            y_pos = max(negLog10P, na.rm = TRUE) * 0.8,
            .groups = "drop") %>%
  mutate(x_pos = ifelse(direction == "Positive",
                        max(res_df$mean_log2FC, na.rm=TRUE) * 0.8,
                        min(res_df$mean_log2FC, na.rm=TRUE) * 0.8))

# Figure 5a - volcano, one panel per cell line (facet by Lineage only
# within each panel), combined via wrap_plots - matching the loop+combine
# style used in Fig 4c/4d and Fig 5b/5c.
fig5a_list <- list()
for (cl in cell_lines) {
  res_df_cl <- res_df %>% dplyr::filter(gt_line == cl)
  counts_df_cl <- counts_df %>% dplyr::filter(gt_line == cl)
  
  is_first  <- cl == cell_lines[1]                    # JHC1 - keeps a (bigger) title
  is_middle <- cl == cell_lines[2]                     # KOLF2.1 - shared y-axis title
  is_last   <- cl == cell_lines[length(cell_lines)]    # O2C3 - shared x-axis title
  
  p <- ggplot(res_df_cl, aes(x = mean_log2FC, y = negLog10P)) +
    geom_point(aes(color = direction), alpha = 0.6) +
    scale_color_manual(values = c("Positive" = "red", "Negative" = "blue", "NotSig" = "grey70")) +
    facet_wrap(~ LineageLabel, nrow = 1, scales = "free_y") +
    theme_bw(base_size = 13) +
    big_text_theme +
    labs(x = NULL, y = NULL, color = "Regulon class") +
    theme(
      axis.text = element_text(size = 18),
      strip.background = element_rect(fill = "grey90"),
      legend.position = "bottom"
    ) +
    geom_text(
      data = counts_df_cl,
      aes(x = x_pos, y = y_pos, label = n, color = direction),
      inherit.aes = FALSE, size = 10, fontface = "bold",
      show.legend = FALSE   # fixes legend showing a colored "a" glyph instead
      # of a dot - geom_text's color aes was being
      # merged into the same legend as geom_point's
    )
  
  if (is_first) {
    p <- p +
      labs(title = paste0(cl, ": Regulon DE by Lineage (significant counts labeled)")) +
      theme(plot.title = element_text(size = 24, face = "bold", margin = margin(b = 10)),
            plot.margin = margin(t = 25, r = 10, b = 10, l = 10))   # room so the bigger title isn't clipped
  }
  if (is_middle) {
    p <- p +
      labs(y = "-log10(median adj p-value)") +
      theme(axis.title.y = element_text(size = 22, face = "bold"))
  }
  if (is_last) {
    p <- p +
      labs(x = "Average log2 Fold Change (+SDF vs -SDF)") +
      theme(axis.title.x = element_text(size = 22, face = "bold"))
  }
  
  fig5a_list[[cl]] <- p
}

fig5a_combined <- wrap_plots(fig5a_list, ncol = 1)
ggsave("~/project/IPSC_2025_Data/Figure5a.png",
       plot = fig5a_combined, device = "png", bg = "white",
       width = 6 * length(lineage_order), height = 5 * length(cell_lines), dpi = 300, limitsize = FALSE)


#Figure 5b
merged_IPSC_tf$dp_lineage <- ifelse(!is.na(merged_IPSC_tf$dp_pseudotime), "dp_ExN_lineage", NA)
merged_IPSC_tf$up_lineage <- ifelse(!is.na(merged_IPSC_tf$up_pseudotime), "up_ExN_lineage", NA)
merged_IPSC_tf$A1_lineage <- ifelse(!is.na(merged_IPSC_tf$A1_pseudotime), "A1_Astrocyte_lineage", NA)
merged_IPSC_tf$A2_lineage <- ifelse(!is.na(merged_IPSC_tf$A2_pseudotime), "A2_Astrocyte_lineage", NA)
merged_IPSC_tf$crn_lineage <- ifelse(!is.na(merged_IPSC_tf$crn_pseudotime), "crn_lineage", NA)
merged_IPSC_tf$epi_lineage <- ifelse(!is.na(merged_IPSC_tf$epi_pseudotime), "epi_lineage", NA)

get_top_regulons <- function(df, alpha = 0.05, n = 5) {
  df %>%
    as.data.frame() %>%
    dplyr::filter(median_p_val_adj <= alpha) %>%
    mutate(direction = ifelse(mean_log2FC > 0, "Up", "Down")) %>%
    group_by(Lineage, gt_line, direction) %>%
    slice_max(order_by = abs(mean_log2FC), n = n, with_ties = FALSE) %>%
    ungroup()
}

top_regulons <- get_top_regulons(der_lineage_results, n = 10)
dim(top_regulons)

# -----------------------------------------------------------------------
# SHARED top-10-positive / top-10-negative regulons: a SINGLE fixed list
# (not computed separately per cell line), so every cell line's heatmap
# shows the exact same regulons in the exact same row order and rows line
# up for direct visual comparison.
#
# NOTE: requiring formal significance (median_p_val_adj <= 0.05) in EVERY
# SINGLE LINEAGE simultaneously was far too strict - zero regulons passed
# out of 262 candidates present in all 3 cell lines. Per request, this is
# now purely a RANKING (no significance gate, no fixed "top 10 of each
# line" ceiling) - a regulon qualifies as "shared" if it's present in all
# cell lines with a CONSISTENT effect direction across them (its overall
# Protocol effect, collapsed across lineage), and the top n by average
# log2FC magnitude are taken regardless of how far down the full ranked
# list that falls.
# -----------------------------------------------------------------------
get_shared_top_regulons <- function(df, cell_lines, n = 10) {
  per_line <- df %>%
    as.data.frame() %>%
    dplyr::group_by(Regulon, gt_line) %>%
    dplyr::summarise(
      mean_log2FC = mean(mean_log2FC, na.rm = TRUE),
      .groups = "drop"
    )
  
  shared <- per_line %>%
    dplyr::group_by(Regulon) %>%
    dplyr::filter(
      dplyr::n_distinct(gt_line) == length(cell_lines),    # present in every cell line
      dplyr::n_distinct(sign(mean_log2FC)) == 1             # same direction in every cell line
    ) %>%
    dplyr::summarise(overall_log2FC = mean(mean_log2FC), .groups = "drop")
  
  pos <- shared %>% dplyr::filter(overall_log2FC > 0) %>% dplyr::slice_max(overall_log2FC, n = n, with_ties = FALSE)
  neg <- shared %>% dplyr::filter(overall_log2FC < 0) %>% dplyr::slice_min(overall_log2FC, n = n, with_ties = FALSE)
  
  list(regulons = c(pos$Regulon, neg$Regulon), pos = pos$Regulon, neg = neg$Regulon)
}

shared_top <- get_shared_top_regulons(der_lineage_results, cell_lines, n = 10)
shared_regulon_order <- shared_top$regulons   # fixed row order used for EVERY cell line's heatmap
cat("Shared top-10 positive regulons:\n"); print(shared_top$pos)
cat("Shared top-10 negative regulons:\n"); print(shared_top$neg)

lineage_df <- merged_IPSC_tf@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::select(
    cell, dp_pseudotime, up_pseudotime, A1_pseudotime, A2_pseudotime,
    crn_pseudotime, epi_pseudotime, Protocol, gt_line
  ) %>%
  tidyr::pivot_longer(
    cols = ends_with("pseudotime"), names_to = "Lineage", values_to = "pseudotime"
  ) %>%
  dplyr::filter(!is.na(pseudotime)) %>%
  dplyr::mutate(Lineage = sub("_pseudotime", "", Lineage))

lineage_df <- lineage_df %>%
  mutate(Protocol_Lineage = paste(Protocol, Lineage, sep = "_"))

auc_mat <- GetAssayData(merged_IPSC_tf, assay = "TF", layer = "data")   # Seurat v5: layer, not slot

# AUC heatmap builder - uses the FIXED shared_regulon_order for every cell
# line (rather than deriving its own regulon list per call), and no longer
# sets its own title/main (a single shared title is drawn once above the
# combined row of heatmaps instead).
plot_auc_heatmap_lineage <- function(regulon_order, auc_df) {
  desired_order <- c(
    "minus_A1",  "plus_A1", "minus_A2",  "plus_A2",
    "minus_dp",  "plus_dp", "minus_up",  "plus_up",
    "minus_crn", "plus_crn", "minus_epi", "plus_epi"
  )
  regulons <- regulon_order[regulon_order %in% colnames(auc_df)]
  if (length(regulons) == 0) return(NULL)
  
  auc_df <- auc_df %>%
    dplyr::filter(Protocol_Lineage %in% desired_order) %>%
    mutate(Protocol_Lineage = factor(Protocol_Lineage, levels = desired_order))
  
  mat <- auc_df %>%
    arrange(Protocol_Lineage) %>%
    as.data.frame() %>%
    column_to_rownames("Protocol_Lineage") %>%
    dplyr::select(all_of(regulons)) %>%
    t()
  # Force BOTH the fixed row order (regulons) AND the fixed column order
  # (Protocol_Lineage combos, from desired_order) - some regulons OR some
  # Protocol/Lineage combos may be entirely absent for a given cell line
  # (e.g. zero cells of a given lineage under one protocol), which
  # previously shrank that cell line's column count silently and shifted
  # every subsequent column/spacer out of alignment with the other cell
  # lines' heatmaps. Reindexing to the full desired_order guarantees every
  # heatmap has the identical column count/order, with any truly-missing
  # combo rendered as a genuine NA (blank) cell instead of disappearing.
  mat_full <- matrix(NA_real_, nrow = length(regulon_order), ncol = length(desired_order),
                     dimnames = list(regulon_order, desired_order))
  common_regs <- intersect(rownames(mat), regulon_order)
  common_cols <- intersect(colnames(mat), desired_order)
  mat_full[common_regs, common_cols] <- mat[common_regs, common_cols]
  mat <- mat_full
  
  new_mat <- mat
  spacer_cols <- integer(0)
  for (i in seq_len(ncol(mat))) {
    if (i == 1) {
      new_mat <- mat[, i, drop = FALSE]
    } else {
      new_mat <- cbind(new_mat, mat[, i, drop = FALSE])
    }
    if (i %% 2 == 0 && i < ncol(mat)) {
      spacer_name <- paste0("spacer_", i)
      new_mat <- cbind(new_mat, rep(NA, nrow(mat)))
      colnames(new_mat)[ncol(new_mat)] <- spacer_name
      spacer_cols <- c(spacer_cols, ncol(new_mat))
    }
  }
  
  # Column labels: "-SDF"/"+SDF" instead of "minus"/"plus", and the
  # lineage codes translated to Figure 4's naming convention for
  # consistency, blank for spacer columns (rather than "spacer_1", etc).
  # NOTE: this renaming is applied ONLY to the display labels here - the
  # `desired_order` used above for filtering/leveling must stay as the
  # actual values in Protocol_Lineage ("minus_A1" etc.) or the filter step
  # silently matches zero rows (which is what caused the earlier error).
  lineage_display_map <- c(
    A1  = "A1_Astrocyte_lineage",
    A2  = "A2_Astrocyte_lineage",
    dp  = "DP_ExN_lineage",
    up  = "UP_ExN_lineage",
    crn = "CRN_lineage",
    epi = "Epithelial_lineage"
  )
  col_labels <- colnames(new_mat)
  is_spacer <- seq_along(col_labels) %in% spacer_cols
  for (i in seq_along(col_labels)) {
    if (is_spacer[i]) { col_labels[i] <- ""; next }
    parts <- strsplit(col_labels[i], "_", fixed = TRUE)[[1]]
    protocol <- parts[1]
    lineage_code <- paste(parts[-1], collapse = "_")
    protocol_label <- dplyr::case_when(
      protocol == "minus" ~ "-SDF",
      protocol == "plus"  ~ "+SDF",
      TRUE ~ protocol
    )
    lineage_label <- if (lineage_code %in% names(lineage_display_map)) {
      lineage_display_map[[lineage_code]]
    } else {
      lineage_code
    }
    col_labels[i] <- paste0(protocol_label, " ", lineage_label)
  }
  
  # Explicit :: needed - ComplexHeatmap (loaded elsewhere in this pipeline)
  # also exports a pheatmap() compatibility wrapper that returns an S4
  # Heatmap object instead of the classic list-based pheatmap object (with
  # a $gtable slot), which breaks any downstream $gtable access. Force the
  # real pheatmap package's function. main = NA (not NULL - NA is what
  # pheatmap's internal is.na() check requires to skip the title cleanly).
  pheatmap::pheatmap(
    new_mat, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
    na_col = "white", border_color = NA, fontsize = 18,
    labels_col = col_labels, fontsize_col = 20, angle_col = 315,
    main = NA, silent = TRUE
  )
}

heatmap_list <- list()
for (cl in cell_lines) {
  auc_df_cl <- lineage_df %>%
    filter(gt_line == cl) %>%
    dplyr::select(cell, Protocol_Lineage) %>%
    left_join(as.data.frame(t(auc_mat)) %>% tibble::rownames_to_column("cell"), by = "cell") %>%
    dplyr::group_by(Protocol_Lineage) %>%
    dplyr::summarise(across(-cell, mean), .groups = "drop")
  
  hm <- plot_auc_heatmap_lineage(shared_regulon_order, auc_df_cl)
  if (!is.null(hm)) {
    # Bigger cell-line label above each heatmap (own text only - the
    # shared "Average AUC..." title is added once, separately, below).
    hm$gtable <- gridExtra::arrangeGrob(hm$gtable, top = grid::textGrob(cl, gp = grid::gpar(fontsize = 22, fontface = "bold")))
    heatmap_list[[cl]] <- hm$gtable
  }
}

heatmap_row <- plot_grid(plotlist = heatmap_list, nrow = 1)

# Shared title drawn ONCE, in its own reserved band above the row of
# heatmaps (rather than repeated per-panel or crammed onto the middle
# panel's own label, which would risk colliding with that panel's cell-
# line name) - guaranteed non-overlapping space, same pattern used
# elsewhere in this pipeline for shared titles/legends.
shared_title <- ggdraw() +
  draw_label("Average AUC per Protocol × Lineage for Top Regulons",
             fontface = "bold", size = 22, x = 0.5, hjust = 0.5)

p_2 <- plot_grid(shared_title, heatmap_row, ncol = 1, rel_heights = c(0.08, 1))

# NOTE: R's png() device has known quirks compositing transparency in
# nested grid/cowplot viewports - fully-NA heatmap regions rendered BLACK
# instead of the specified na_col="white" when saved directly as png,
# even though an isolated pheatmap test confirmed na_col works correctly
# on its own. Saving as PNG first (which handles this correctly, verified
# via that same isolated test) and converting to png afterward via
# magick avoids the png() device's compositing issue entirely.
png_path_5b <- "~/project/IPSC_2025_Data/Figure5b.png"
png_path_5b <- "~/project/IPSC_2025_Data/Figure5b.png"
ggsave(png_path_5b,
       plot = p_2, device = "png", bg = "white",
       width = 10 * length(cell_lines), height = 20, dpi = 300, limitsize = FALSE)
magick::image_write(magick::image_read(png_path_5b), path = png_path_5b, format = "png")



# =========================================================================
# STANDALONE DIAGNOSTIC (does not modify Figure_5.R): find regulons that
# are "lineage_opposed" consistently across ALL 3 cell lines, split into
# Pallial lineages (dp/up = neuron, A1/A2 = astrocyte) and Hem lineages
# (crn = neuron, epi = epithelial). Recomputes cor_results2 per cell line
# from scratch (same logic as the existing Fig5c loop) rather than
# depending on anything the loop saved.
# =========================================================================

find_consistent_opposed <- function(df, lineage_group, group_label) {
  n_lines_total <- n_distinct(df$gt_line)
  
  # 1. Check for consistent lineage_opposed
  per_line_opposed <- df %>%
    dplyr::filter(Pseudotime %in% lineage_group) %>%
    dplyr::group_by(Regulon, gt_line) %>%
    dplyr::summarise(
      n_opposed = sum(color_group == "lineage_opposed", na.rm = TRUE),
      n_tested  = dplyr::n(),
      opposed_in_line = if (REQUIRE_ALL_LINEAGES) n_opposed == n_tested else n_opposed > 0,
      .groups = "drop"
    )
  
  consistent_opposed <- per_line_opposed %>%
    dplyr::group_by(Regulon) %>%
    dplyr::summarise(n_lines_opposed = sum(opposed_in_line), .groups = "drop") %>%
    dplyr::filter(n_lines_opposed == n_lines_total) %>%
    dplyr::pull(Regulon)
  
  # 2. Check for consistent lineage_aligned
  per_line_aligned <- df %>%
    dplyr::filter(Pseudotime %in% lineage_group) %>%
    dplyr::group_by(Regulon, gt_line) %>%
    dplyr::summarise(
      n_aligned = sum(color_group == "lineage_aligned", na.rm = TRUE),
      n_tested  = dplyr::n(),
      aligned_in_line = if (REQUIRE_ALL_LINEAGES) n_aligned == n_tested else n_aligned > 0,
      .groups = "drop"
    )
  
  consistent_aligned <- per_line_aligned %>%
    dplyr::group_by(Regulon) %>%
    dplyr::summarise(n_lines_aligned = sum(aligned_in_line), .groups = "drop") %>%
    dplyr::filter(n_lines_aligned == n_lines_total) %>%
    dplyr::pull(Regulon)
  
  # Print results for both directions
  cat("\n", group_label, "(", paste(lineage_group, collapse = ", "), ") - OPPOSED in all", n_lines_total, "cell lines:\n")
  if (length(consistent_opposed) == 0) cat("  (none found)\n") else print(consistent_opposed)
  
  cat("\n", group_label, "(", paste(lineage_group, collapse = ", "), ") - ALIGNED in all", n_lines_total, "cell lines:\n")
  if (length(consistent_aligned) == 0) cat("  (none found)\n") else print(consistent_aligned)
  
  # Return as a named list so downstream variables hold both directions
  invisible(list(opposed = consistent_opposed, aligned = consistent_aligned))
}

cat("\n=== PALLIAL LINEAGES ===\n")
# Capture the full list (both opposed and aligned) for neurons and astrocytes
pallial_neuron_res    <- find_consistent_opposed(cor_results_all, pallial_neuron_lineages, "Pallial neuron")
pallial_astrocyte_res <- find_consistent_opposed(cor_results_all, pallial_astrocyte_lineages, "Pallial astrocyte")

# Find the specific intersections you are looking for
neuron_aligned_astro_opposed <- intersect(pallial_neuron_res$aligned, pallial_astrocyte_res$opposed)
neuron_opposed_astro_aligned <- intersect(pallial_neuron_res$opposed, pallial_astrocyte_res$aligned)

# Print the final sets
cat("\n=== DIVERGENT REGULONS (PALLIAL) ===\n")
cat("Neuron-Aligned AND Astrocyte-Opposed:\n")
if (length(neuron_aligned_astro_opposed) == 0) {
  cat("  (none found)\n") 
} else {
  print(neuron_aligned_astro_opposed)
}

cat("\nNeuron-Opposed AND Astrocyte-Aligned:\n")
if (length(neuron_opposed_astro_aligned) == 0) {
  cat("  (none found)\n") 
} else {
  print(neuron_opposed_astro_aligned)
}




# -----------------------------------------------------------------------
# Figure 5c: regulon-pseudotime correlation, kept SEPARATE per cell line
# (previously run only on the "minus" subset pooled across all lines).
# Consistent lineage labels applied, matching Panels A/B.
# -----------------------------------------------------------------------
fig5c_list <- list()
legend_source_plot <- NULL   # captured from one panel, used to build ONE shared (larger) legend

for (cl in cell_lines) {
  
  merged_minus <- subset(merged_IPSC_tf, Protocol == "minus" & gt_line == cl)
  merged_minus@misc$SCENIC$RegulonsAUC <- merged_minus@misc$SCENIC$RegulonsAUC[colnames(merged_minus),]
  tf_auc <- merged_minus@misc$SCENIC$RegulonsAUC
  tf_auc <- tf_auc %>% mutate(across(everything(), ~ replace_na(., 0)))
  tf_counts <- t(tf_auc)
  # Seurat v5 Assay5: update layers via LayerData<-(), not @counts<-/@data<-
  LayerData(merged_minus, assay = "TF", layer = "counts") <- tf_counts
  LayerData(merged_minus, assay = "TF", layer = "data") <- tf_counts
  regulon_mat <- t(LayerData(merged_minus, assay = "TF", layer = "counts"))
  
  pseudotimes <- list(
    dp  = merged_minus$dp_pseudotime, up  = merged_minus$up_pseudotime,
    A1  = merged_minus$A1_pseudotime, A2  = merged_minus$A2_pseudotime,
    epi = merged_minus$epi_pseudotime, crn = merged_minus$crn_pseudotime
  )
  correlate_regulons <- function(mat, pseudotime) {
    apply(mat, 2, function(reg) {
      # Guard against too few finite paired observations (e.g. a lineage
      # with very few or zero cells within this cell line x protocol
      # subset, since pseudotime is NA for any cell not on that lineage) -
      # cor.test() errors with "not enough finite observations" below 2
      # pairs, and needs >=3 for a meaningful t-based p-value anyway.
      complete_idx <- is.finite(reg) & is.finite(pseudotime)
      if (sum(complete_idx) < 3) {
        # NOTE: names must exactly match the success branch below -
        # ct$estimate is already named "cor" (pearson), so c(cor=...)
        # produces the compound name "cor.cor", not "cor". Downstream code
        # references cor_results2$cor.cor directly, so this NA fallback
        # must match that naming exactly or apply()'s matrix assembly
        # could silently produce inconsistent column names.
        return(c(cor.cor = NA_real_, pval = NA_real_))
      }
      ct <- cor.test(reg[complete_idx], pseudotime[complete_idx], method = "pearson")
      c(cor = ct$estimate, pval = ct$p.value)
    }) %>% t() %>% as.data.frame()
  }
  cor_results <- lapply(names(pseudotimes), function(pt) {
    df <- correlate_regulons(regulon_mat, pseudotimes[[pt]])
    df$Regulon <- rownames(df)
    df$Pseudotime <- pt
    df
  }) %>% bind_rows()
  
  cor_results2 <- cor_results %>%
    left_join(
      der_lineage_results %>% filter(gt_line == cl) %>%
        dplyr::select(Regulon, Lineage, median_p_val_adj, max_p_val_adj, mean_log2FC),
      by = c("Regulon", "Pseudotime" = "Lineage")
    )
  
  threshold_cor <- 0
  cor_results2 <- cor_results2 %>%
    mutate(
      de_dir = case_when(
        !is.na(median_p_val_adj) & median_p_val_adj <= 0.05 & mean_log2FC > 0 ~ "up",
        !is.na(median_p_val_adj) & median_p_val_adj <= 0.05 & mean_log2FC < 0 ~ "down",
        TRUE ~ "ns"
      ),
      is_green = case_when(de_dir == "down" & cor.cor < -threshold_cor ~ TRUE, de_dir == "up" & cor.cor > threshold_cor ~ TRUE, TRUE ~ FALSE),
      is_red   = case_when(de_dir == "down" & cor.cor >  threshold_cor ~ TRUE, de_dir == "up" & cor.cor < -threshold_cor ~ TRUE, TRUE ~ FALSE),
      color_group = case_when(is_green ~ "lineage_aligned", is_red ~ "lineage_opposed", TRUE ~ "neutral"),
      PseudotimeLabel = factor(unname(lineage_labels[Pseudotime]), levels = lineage_order_labeled)
    )
  
  cor_results_plot <- cor_results2 %>% dplyr::filter(!is.na(mean_log2FC))
  green_counts <- cor_results2 %>% dplyr::filter(is_green) %>% count(PseudotimeLabel, name = "n_green")
  red_counts   <- cor_results2 %>% dplyr::filter(is_red) %>% count(PseudotimeLabel, name = "n_red")
  
  label_genes_epi_crn <- c("TCF7L1(+)", "SOX13(+)")
  #label_genes_dp_up_A1_A2 <- c("POU3F1(+)", "STAT3(+)", "NFIA(+)", "NFIX(+)", "NFIC(+)")
  label_genes_dp_up_A1_A2 <- c("POU3F1(+)", "E2F2(+)", "HMGA2(+)", "NFIX(+)", "NFIC(+)")
  
  is_first  <- cl == cell_lines[1]                    # JHC1
  is_middle <- cl == cell_lines[2]                     # KOLF2.1 - shared y-axis title
  is_last   <- cl == cell_lines[length(cell_lines)]    # O2C3 - shared x-axis title
  
  p <- ggplot(cor_results_plot, aes(x = cor.cor, y = mean_log2FC)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(aes(color = color_group), alpha = 0.7, size = 1.8) +
    scale_color_manual(values = c("lineage_aligned" = "forestgreen", "lineage_opposed" = "firebrick", "neutral" = "grey70")) +
    facet_wrap(~ PseudotimeLabel, scales = "free", nrow = 1) +
    geom_point(data = subset(cor_results_plot, Regulon %in% label_genes_dp_up_A1_A2 &
                               Pseudotime %in% c("dp", "up", "A1", "A2")), color = "black", size = 2.5) +
    geom_text_repel(data = subset(cor_results_plot, Regulon %in% label_genes_dp_up_A1_A2 &
                                    Pseudotime %in% c("dp", "up", "A1", "A2")),
                    aes(label = Regulon), size = 5.5, color = "black", max.overlaps = Inf,
                    box.padding = 0.4, point.padding = 0.3, segment.color = "black",
                    segment.size = 0.4, min.segment.length = 0) +
    geom_point(data = subset(cor_results_plot, Regulon %in% label_genes_epi_crn &
                               Pseudotime %in% c("crn", "epi")), color = "black", size = 2.5) +
    geom_text_repel(data = subset(cor_results_plot, Regulon %in% label_genes_epi_crn &
                                    Pseudotime %in% c("crn", "epi")),
                    aes(label = Regulon), size = 5.5, color = "black", max.overlaps = Inf,
                    box.padding = 0.4, point.padding = 0.3, segment.color = "black",
                    segment.size = 0.4, min.segment.length = 0) +
    geom_text(data = green_counts, aes(x = -Inf, y = Inf, label = paste0("aligned = ", n_green)),
              hjust = -0.1, vjust = 2.2, size = 6, color = "forestgreen", inherit.aes = FALSE) +
    geom_text(data = red_counts, aes(x = -Inf, y = Inf, label = paste0("opposed = ", n_red)),
              hjust = -0.1, vjust = 3.8, size = 6, color = "firebrick", inherit.aes = FALSE) +
    labs(
      x = NULL, y = NULL, color = "Regulon class",
      title = paste0(cl, ": Regulon alignment of synchronized response with lineage dynamics")
    ) +
    theme_bw(base_size = 13) + big_text_theme +
    theme(
      strip.background = element_rect(fill = "grey90"),
      axis.text = element_text(size = 15),
      plot.title = element_text(size = 20, face = "bold", margin = margin(b = 10)),
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10),   # room so bigger titles aren't clipped
      legend.position = "none"    # shared legend added once at combine time instead
    )
  
  if (is_middle) {
    p <- p + labs(y = "Lineage-level Average log2FC (+SDF vs -SDF protocol)") +
      theme(axis.title.y = element_text(size = 20, face = "bold"))
  }
  if (is_last) {
    p <- p + labs(x = "Pearson Correlation with pseudotime") +
      theme(axis.title.x = element_text(size = 20, face = "bold"))
  }
  
  # Capture a legend from ANY one panel (with a bigger legend theme applied)
  # to build the single shared legend added once at combine time.
  if (is.null(legend_source_plot)) {
    legend_source_plot <- p +
      theme(legend.position = "bottom",
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 18, face = "bold"),
            legend.key.size = unit(1.2, "cm"))
  }
  
  fig5c_list[[cl]] <- p
}

fig5c_legend <- cowplot::get_legend(legend_source_plot)
fig5c_combined <- plot_grid(
  wrap_plots(fig5c_list, ncol = 1),
  fig5c_legend,
  ncol = 1, rel_heights = c(1, 0.06)
)
ggsave("~/project/IPSC_2025_Data/Figure5c.png",
       plot = fig5c_combined, device = "png", bg = "white",
       width = 24, height = 8 * length(cell_lines), dpi = 300, limitsize = FALSE)

# -----------------------------------------------------------------------
# Final combined Figure 5: layout (A / C) | B - A stacked above C in a
# left column, B occupying a right column spanning the same total height.
# Since both columns share the SAME overall height by construction, A and
# C each naturally get half of B's height (i.e. B ends up ~2x the height
# of A and ~2x the height of C), without needing a special ratio.
# -----------------------------------------------------------------------
AC_stack <- cowplot::plot_grid(fig5a_combined, fig5c_combined, ncol = 1, rel_heights = c(1, 1),
                               labels = c("A", "C"), label_size = 24)
fig5_final <- cowplot::plot_grid(AC_stack, p_2, nrow = 1, rel_widths = c(1, 1),
                                 labels = c("", "B"), label_size = 24)
png_path_final <- "~/project/IPSC_2025_Data/Figure5_combined.png"
png_path_final <- "~/project/IPSC_2025_Data/Figure5_combined.png"
ggsave(png_path_final,
       plot = fig5_final, device = "png", bg = "white",
       width = 40, height = 24 * length(cell_lines), dpi = 300, limitsize = FALSE)
magick::image_write(magick::image_read(png_path_final), path = png_path_final, format = "png")
