
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


Cortical_lineage <- subset(merged_IPSC, Celltype == "DL_ExN" | Celltype == "UL_ExN" | Celltype == "RG" | Celltype == "Astrocyte" | Celltype == "IPC_ExN")
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



# Figure3a
p1_1 <- DimPlot(subset(Cortical_lineage, Protocol == "minus"), group.by = "Celltype2", pt.size = 2)
p1_2 <- DimPlot(subset(Cortical_lineage, Protocol == "plus"), group.by = "Celltype2", pt.size = 2)
p1_3 <- DimPlot(subset(Hem_lineage, Protocol == "minus"), group.by = "Celltype2", pt.size = 2)
p1_4 <- DimPlot(subset(Hem_lineage, Protocol == "plus"), group.by = "Celltype2", pt.size = 2)
combined_plot <- ggarrange(
  p1_1, p1_2,
  p1_3, p1_4,
  ncol = 2, nrow = 2,
  font.label = list(size = 20)
)
ggsave(
  "~/project/IPSC_2025_Data/Figure3a_Dimplot_Celltype_Combined.tiff",
  plot = combined_plot,
  device = "tiff",
  width = 20, height = 20, dpi = 300
)

# Supplemental Figure 4
p <- FeaturePlot(subset(Cortical_lineage, features = c("TNFSF12", "S100A10", "C3", "FKBP5"))
ggsave(
  "~/project/IPSC_2025_Data/Supplmentary_Figure4.tiff",
  plot = p,
  device = "tiff",
  width = 20, height = 20, dpi = 300
)

#Merge for differential expression analysis
merged_temp <- merge(Cortical_lineage, Hem_lineage)
merged_temp <- JoinLayers(merged_temp)
cnts <- merged_temp@assays$RNA@layers$counts
colnames(cnts) <- colnames(merged_temp)
rownames(cnts) <- rownames(merged_temp)
merged_temp <- CreateSeuratObject(counts = cnts, meta.data = merged_temp@meta.data)
merged_temp <- NormalizeData(merged_temp)
IPSC_sce <- as.SingleCellExperiment(merged_temp)
rm(merged_temp)

# only keep singlet cells with sufficient reads
IPSC_sce <- IPSC_sce[rowSums(assay(IPSC_sce, "counts") > 0) > 0, ]
# compute QC metrics
qc <- perCellQCMetrics(IPSC_sce)
# remove cells with few or many detected genes
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
IPSC_sce <- IPSC_sce[, !ol]
pb <- aggregateToPseudoBulk(IPSC_sce,
                            assay = "counts",
                            cluster_id = "Celltype2",
                            sample_id = "SampleID",
                            verbose = FALSE
)
assayNames(pb)
res.proc <- processAssays(pb, ~Protocol + Age + gt_line, min.count = 5)
res.dl <- dreamlet(res.proc, ~Protocol + Age + gt_line)
res_mash = run_mash(res.dl, coef='Protocolplus')
plotVolcano(res_mash)


#Figure3b
desired_order <- c("RG", "IPC_ExN", "DL_ExN", "UL_ExN","A1 Astrocyte", "A2 Astrocyte","Hem_RG", "CRN", "Epithelial")
pm   <- get_pm(res_mash$model)
genes <- rownames(pm)
cells <- colnames(pm)
psd <- get_psd(res_mash$model)
lfsr <- get_lfsr(res_mash$model)
prot_coef <- res_mash$coefList

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
    )
  )
res_df$Celltype <- factor(res_df$Celltype, levels = desired_order)

plot_mash_volcano <- function(res_df, title = "mashr Volcano Plot") {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # 1. Per-celltype y-axis label positions (avoid Inf)
  ypos_df <- res_df %>%
    group_by(Celltype) %>%
    summarise(
      ymin = min(-log10(adj.P.Val), na.rm = TRUE),
      ymax_raw = max(-log10(adj.P.Val), na.rm = TRUE),
      ymax = pmin(ymax_raw, 50),              # cap Inf or huge values
      ylab = ymin + 0.80 * (ymax - ymin),
      .groups = "drop"
    )
  
  # 2. Per-celltype x-axis label positions (avoid NA)
  xpos_df <- res_df %>%
    group_by(Celltype) %>%
    summarise(
      xmin_raw = suppressWarnings(min(logFC, na.rm = TRUE)),
      xmax_raw = suppressWarnings(max(logFC, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      xmin = ifelse(is.finite(xmin_raw), xmin_raw, 0),
      xmax = ifelse(is.finite(xmax_raw), xmax_raw, 0),
      fallback = !is.finite(xmin_raw) | !is.finite(xmax_raw)
    )
  
  # 3. Count significant genes per direction
  counts <- res_df %>%
    group_by(Celltype, direction) %>%
    summarise(n_sig = sum(sig), .groups = "drop") %>%
    filter(direction != "ns") %>%        # <--- remove NA rows
    left_join(ypos_df, by = "Celltype") %>%
    left_join(xpos_df, by = "Celltype") %>%
    mutate(
      x_pos = ifelse(direction == "pos", xmax, xmin),
      hjust = ifelse(direction == "pos", 1.1, -0.1)
    )
  
  # 4. Volcano plot
  ggplot(res_df, aes(x = logFC, y = -log10(adj.P.Val), color = color)) +
    geom_point(alpha = 0.6, size = 0.8) +
    facet_wrap(~ Celltype, scales = "free") +
    scale_color_identity() +
    theme_bw() +
    labs(
      x = "Posterior mean logFC (mashr)",
      y = "-log10(LFSR)",
      title = title
    ) +
    geom_text(
      data = counts,
      aes(
        x = x_pos,
        y = ylab,
        label = paste0(n_sig, " gene(s)"),
        hjust = hjust
      ),
      inherit.aes = FALSE,
      vjust = 0,
      size = 3,
      fontface = "bold"
    )
}
pv <- plot_mash_volcano(res_df, title = "mashr Volcano Plot — Protocolplus")
ggsave("~/project/IPSC_2025_Data/Figure3b_DEG_volcano_plots.tiff",
       plot = pv,
       device = "tiff",
       width = 10, height = 10, dpi = 300)


get_top_pos_genes_resdf <- function(df, celltype, n = nrow(df)) {
  tt <- df %>%
    filter(Celltype == celltype,
           logFC > 0,
           adj.P.Val < 0.05) %>%
    arrange(adj.P.Val, desc(logFC))
  head(tt$Gene, n)
}

# Figure2c
get_top_neg_genes_resdf <- function(df, celltype, n = nrow(df)) {
  tt <- df %>%
    filter(Celltype == celltype,
           logFC < 0,
           adj.P.Val < 0.05) %>%
    arrange(adj.P.Val, logFC)
  head(tt$Gene, n)
}
celltypes_resdf <- sort(unique(res_df$Celltype))
top_pos_genes_resdf <- lapply(celltypes_resdf, function(ct)
  get_top_pos_genes_resdf(res_df, ct))
top_neg_genes_resdf <- lapply(celltypes_resdf, function(ct)
  get_top_neg_genes_resdf(res_df, ct))
names(top_pos_genes_resdf) <- celltypes_resdf
names(top_neg_genes_resdf) <- celltypes_resdf
make_similarity_matrix <- function(top_list) {
  assays <- names(top_list)
  n <- length(assays)

  sim_mat <- matrix(0, n, n, dimnames = list(assays, assays))
  overlap_mat <- matrix(0, n, n, dimnames = list(assays, assays))
  
  for (i in seq_along(assays)) {
    for (j in seq_along(assays)) {
      set1 <- top_list[[i]]
      set2 <- top_list[[j]]
      
      overlap <- length(intersect(set1, set2))
      denom   <- min(length(set1), length(set2))   # Overlap coefficient denominator
      
      sim_mat[i, j] <- ifelse(denom > 0, overlap / denom, 0)
      overlap_mat[i, j] <- overlap
    }
  }
  
  list(similarity = sim_mat, overlap = overlap_mat)
}
pos_similarity_matrix_resdf <- make_similarity_matrix(top_pos_genes_resdf)
neg_similarity_matrix_resdf <- make_similarity_matrix(top_neg_genes_resdf)
sim_pos <- pos_similarity_matrix_resdf$similarity
sim_neg <- neg_similarity_matrix_resdf$similarity
sim_pos[!is.finite(sim_pos)] <- 0
sim_neg[!is.finite(sim_neg)] <- 0
global_min <- min(sim_pos, sim_neg)
global_max <- max(sim_pos, sim_neg)
global_min
global_max
my_breaks <- seq(global_min, global_max, length.out = 101)
my_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
p_sim_pos_1 <- pheatmap(
  sim_pos,
  display_numbers = pos_similarity_matrix_resdf$overlap,
  fontsize_number = 12,
  main = "Positive DE gene similarity (mashr, adj.P.Val < 0.05)",
  color = my_colors,
  breaks = my_breaks
)
p_sim_neg_2 <- pheatmap(
  sim_neg,
  display_numbers = neg_similarity_matrix_resdf$overlap,
  fontsize_number = 12,
  main = "Negative DE gene similarity (mashr, adj.P.Val < 0.05)",
  color = my_colors,
  breaks = my_breaks
)
setwd("~/project/IPSC_2025_Data/Figure_Parts/Manuscript1")
tiff("Figure_3c_1.tiff", width = 750, height = 750)
p_sim_pos_1
dev.off()
setwd("~/project/IPSC_2025_Data/Figure_Parts/Manuscript1")
tiff("Figure_3c_2.tiff", width = 750, height = 750)
p_sim_neg_2
dev.off()

#Figure3d
df_gs_bp = zenith_gsa(res_mash, go.gs.bp)
df_gs_cc = zenith_gsa(res_mash, go.gs.cc)
df_gs_cc$assay <- factor(df_gs_cc$assay, levels = desired_order)
df_gs_bp$assay <- factor(df_gs_bp$assay, levels = desired_order)
# Heatmap of results
setwd("~/project/IPSC_2025_Data/Figure_Parts/Manuscript1")
tiff("Figure_3d_1.tiff", width = 750, height = 1200)
plotZenithResults(df_gs_bp, 5, 5)
dev.off()
setwd("~/project/IPSC_2025_Data/Figure_Parts/Manuscript1")
tiff("Figure_3d_2.tiff", width = 750, height = 1200)
plotZenithResults(df_gs_cc, 5, 5)
dev.off()

#Supplementary Figure 5
vp.lst <- fitVarPart(res.proc, ~Protocol + Age + gt_line)
vp.global <- vp.lst %>%
  as.data.frame() %>%
  group_by(assay) %>%
  summarise(across(c(Protocol, Age, gt_line, Residuals),
                   ~mean(as.numeric(.), na.rm = TRUE)))
# Convert to plain data.frame
vp.global <- as.data.frame(vp.global)
# Set cell types as row names
rownames(vp.global) <- vp.global$assay
vp.global$assay <- NULL
# Ensure numeric
vp.global[] <- lapply(vp.global, as.numeric)
# Now plot
p <- plotPercentBars(vp.global)
ggsave(
  "~/project/IPSC_2025_Data/Supplmentary_Figure5.tiff",
  plot = p,
  device = "tiff",
  width = 10, height = 10, dpi = 300
)

