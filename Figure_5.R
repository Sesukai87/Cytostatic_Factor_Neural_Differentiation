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


#Saving File To Run pySCENIC
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
cnts <- merged_IPSC@assays$RNA@counts
colnames(cnts) <- colnames(merged_IPSC)
rownames(cnts) <- rownames(merged_IPSC)
merged_IPSC <- CreateSeuratObject(counts = cnts, meta.data = merged_IPSC@meta.data, min.cells = 5)
merged_IPSC[["RNA"]] <- as(object = merged_IPSC[["RNA"]], Class = "Assay")
SaveH5Seurat(merged_IPSC, filename = "~/project/IPSC_2025_Data/merged_IPSC.h5Seurat")
Convert("~/project/IPSC_2025_Data/merged_IPSC.h5Seurat", dest = "h5ad")


#After Runing pySCENIC
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
scenic_loom_path <- "~/project/IPSC_2025_Data/merged_IPSC_aucell.loom"
merged_IPSC <- ImportPyscenicLoom(scenic_loom_path, seu = merged_IPSC)
dim(merged_IPSC)
dim(merged_IPSC@misc$SCENIC$RegulonsAUC)


Cortical_lineage <- readRDS("~/project/IPSC_2025_Data/merged_IPSC_derived_pallial_lineages")
Hem_lineage <- readRDS("~/project/IPSC_2025_Data/merged_IPSC_derived_hem_lineages")

meta_hem <- Hem_lineage@meta.data[, c("Celltype2", grep("pseudotime", colnames(Hem_lineage@meta.data), value = TRUE), grep("lineage", colnames(Hem_lineage@meta.data), value = TRUE))]
meta_hem$pseudotime <- NULL

meta_cor <- Cortical_lineage@meta.data[, c("Celltype2", grep("pseudotime", colnames(Cortical_lineage@meta.data), value = TRUE), grep("lineage", colnames(Cortical_lineage@meta.data), value = TRUE))]
meta_cor$pseudotime <- NULL

meta_hem[] <- lapply(meta_hem, function(x) if (is.factor(x)) as.character(x) else x)
meta_cor[] <- lapply(meta_cor, function(x) if (is.factor(x)) as.character(x) else x)


# Hem lineage metadata transfer
common_hem <- intersect(rownames(meta_hem), colnames(merged_IPSC))
merged_IPSC@meta.data[common_hem, colnames(meta_hem)] <- meta_hem[common_hem, ]
table(merged_IPSC$Celltype2)
# Cortical lineage metadata transfer
common_cor <- intersect(rownames(meta_cor), colnames(merged_IPSC))
merged_IPSC@meta.data[common_cor, colnames(meta_cor)] <- meta_cor[common_cor, ]
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
merged_IPSC@assays$TF@counts <- tf_counts
merged_IPSC@assays$TF@data <- tf_counts
DefaultAssay(merged_IPSC) <- 'RNA'
dim(merged_IPSC)
DefaultAssay(merged_IPSC) <- 'TF'
dim(merged_IPSC)
merged_IPSC <- ScaleData(merged_IPSC)
merged_IPSC_tf <- merged_IPSC
DefaultAssay(merged_IPSC_tf) <- "TF"
merged_IPSC_tf@assays$RNA <- NULL

#Figure 5a
lineage_df <- merged_IPSC_tf@meta.data %>%
  tibble::rownames_to_column("cell") %>%        
  dplyr::select(
    cell,
    dp_pseudotime,
    up_pseudotime,
    A1_pseudotime,
    A2_pseudotime,
    crn_pseudotime,
    epi_pseudotime,
    Age,
    gt_line,
    Protocol
  ) %>%
  tidyr::pivot_longer(
    cols = ends_with("pseudotime"),
    names_to = "Lineage",
    values_to = "pseudotime"
  ) %>%
  dplyr::filter(!is.na(pseudotime)) %>%
  dplyr::mutate(Lineage = sub("_pseudotime", "", Lineage))

lineages <- unique(lineage_df$Lineage)

Idents(merged_IPSC_tf) <- "Protocol"

# Reset identities to avoid collisions
merged_IPSC_tf$dummy_ident <- "all_cells"
Idents(merged_IPSC_tf) <- "dummy_ident"

der_lineage_list <- lapply(lineages, function(lin) {
  
  cells_lin <- lineage_df %>% dplyr::filter(Lineage == lin)
  
  groups <- cells_lin %>%
    mutate(group = interaction(Age, gt_line, drop = TRUE)) %>%
    split(.$group)
  
  group_results <- lapply(groups, function(g) {
    
    cells_plus  <- g %>% dplyr::filter(Protocol == "plus")  %>% pull(cell)
    cells_minus <- g %>% dplyr::filter(Protocol == "minus") %>% pull(cell)
    
    if (length(cells_plus) < 3 || length(cells_minus) < 3) {
      return(NULL)
    }
    
    FindMarkers(
      object = merged_IPSC_tf[["TF"]],
      cells.1 = cells_plus,
      cells.2 = cells_minus,
      logfc.threshold = 0,
      min.pct = 0.25,
      slot = "data"
    )
  })
  
  group_results <- group_results[!sapply(group_results, is.null)]
  if (length(group_results) == 0) return(NULL)
  
  group_results <- lapply(group_results, function(df) {
    df <- as.data.frame(df)
    df$Regulon <- rownames(df)
    df
  })
  
  combined <- bind_rows(group_results, .id = "Group")
  
  collapsed <- combined %>%
    group_by(Regulon) %>%
    summarise(
      median_p_val_adj = median(p_val_adj, na.rm = TRUE),
      max_p_val_adj    = max(p_val_adj, na.rm = TRUE),
      mean_log2FC      = mean(avg_log2FC, na.rm = TRUE)
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
    )
  )

# Compute counts and per-facet y position
counts_df <- res_df %>%
  dplyr::filter(sig) %>%
  group_by(Lineage, direction) %>%
  summarise(n = n(),
            y_pos = max(negLog10P, na.rm = TRUE) * 0.8,   # per Celltype max
            .groups = "drop") %>%
  mutate(x_pos = ifelse(direction == "Positive", 
                        max(res_df$mean_log2FC, na.rm=TRUE) * 0.8,
                        min(res_df$mean_log2FC, na.rm=TRUE) * 0.8))


desired_order <- c("dp", "up", "A1", "A2", "epi", "crn")

# Make sure Celltype is a factor with the specified levels
res_df$Lineage    <- factor(res_df$Lineage, levels = desired_order)
counts_df$Lineage <- factor(counts_df$Lineage, levels = desired_order)

#Figure 5a
# Volcano plot
p_1 <- ggplot(res_df, aes(x = mean_log2FC, y = negLog10P)) +
  geom_point(aes(color = direction), alpha = 0.6) +
  scale_color_manual(values = c("Positive" = "red",
                                "Negative" = "blue",
                                "NotSig"   = "grey70")) +
  facet_wrap(~ Lineage, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Average log2 Fold Change",
    y = "-log10(median adj p-value)",
    title = "Faceted Volcano Plots by Lineage (significant counts labeled)"
  ) + theme_bw(base_size = 13) +
  theme(
    strip.text.x = element_text(size = 24, face = "bold"),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom"
  )+ geom_text(
    data = counts_df,
    aes(x = x_pos, y = y_pos, label = n, color = direction),
    inherit.aes = FALSE,
    size = 4,
    fontface = "bold"
  ) 
ggsave("~/project/IPSC_2025_Data/Figure5a.tiff",
       plot = p_1,
       device = "tiff",
       width = 20, height = 20, dpi = 300)


#Figure 5b
merged_IPSC_tf$dp_lineage <- ifelse(
  !is.na(merged_IPSC_tf$dp_pseudotime),
  "dp_ExN_lineage",
  NA
)
merged_IPSC_tf$up_lineage <- ifelse(
  !is.na(merged_IPSC_tf$up_pseudotime),
  "up_ExN_lineage",
  NA
)
merged_IPSC_tf$A1_lineage <- ifelse(
  !is.na(merged_IPSC_tf$A1_pseudotime),
  "A1_Astrocyte_lineage",
  NA
)
merged_IPSC_tf$A2_lineage <- ifelse(
  !is.na(merged_IPSC_tf$A2_pseudotime),
  "A2_Astrocyte_lineage",
  NA
)
merged_IPSC_tf$crn_lineage <- ifelse(
  !is.na(merged_IPSC_tf$crn_pseudotime),
  "crn_lineage",
  NA
)
merged_IPSC_tf$epi_lineage <- ifelse(
  !is.na(merged_IPSC_tf$epi_pseudotime),
  "epi_lineage",
  NA
)

get_top_regulons <- function(df, alpha = 0.05, n = 5) {
  df %>%
    as.data.frame() %>%
    dplyr::filter(median_p_val_adj <= alpha) %>%
    mutate(direction = ifelse(mean_log2FC > 0, "Up", "Down")) %>%
    group_by(Lineage, direction) %>%
    slice_max(order_by = abs(mean_log2FC), n = n, with_ties = FALSE) %>%
    ungroup()
}

top_regulons <- get_top_regulons(der_lineage_results, n = 10)
dim(top_regulons)

lineage_df <- merged_IPSC_tf@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::select(
    cell,
    dp_pseudotime, up_pseudotime,
    A1_pseudotime, A2_pseudotime,
    crn_pseudotime, epi_pseudotime,
    Protocol
  ) %>%
  tidyr::pivot_longer(
    cols = ends_with("pseudotime"),
    names_to = "Lineage",
    values_to = "pseudotime"
  ) %>%
  dplyr::filter(!is.na(pseudotime)) %>%
  dplyr::mutate(Lineage = sub("_pseudotime", "", Lineage))

lineage_df <- lineage_df %>%
  mutate(Protocol_Lineage = paste(Protocol, Lineage, sep = "_"))

auc_mat <- GetAssayData(merged_IPSC_tf, assay = "TF", slot = "data")

auc_df <- lineage_df %>%
  dplyr::select(cell, Protocol_Lineage) %>%
  left_join(
    as.data.frame(t(auc_mat)) %>% tibble::rownames_to_column("cell"),
    by = "cell"
  ) %>%
  dplyr::group_by(Protocol_Lineage) %>%
  dplyr::summarise(across(-cell, mean), .groups = "drop")


plot_auc_heatmap_lineage <- function(top_df, auc_df) {
  library(pheatmap)
  library(dplyr)
  library(tibble)
  
  # Desired Protocol_Lineage order
  desired_order <- c(
    "minus_A1",  "plus_A1",
    "minus_A2",  "plus_A2",
    "minus_dp",  "plus_dp",
    "minus_up",  "plus_up",
    "minus_crn", "plus_crn",
    "minus_epi", "plus_epi"
  )
  
  regulons <- unique(top_df$Regulon)
  
  # Keep only regulons present in AUC matrix
  regulons <- regulons[regulons %in% colnames(auc_df)]
  if (length(regulons) == 0) stop("No regulons matched between DE results and AUC matrix.")
  
  # Filter AUC df to only desired order (and drop missing)
  auc_df <- auc_df %>%
    dplyr::filter(Protocol_Lineage %in% desired_order) %>%
    mutate(Protocol_Lineage = factor(Protocol_Lineage, levels = desired_order))
  
  # Build matrix for heatmap
  mat <- auc_df %>%
    arrange(Protocol_Lineage) %>%
    as.data.frame() %>%
    column_to_rownames("Protocol_Lineage") %>%
    dplyr::select(all_of(regulons)) %>%
    t()
  
  # ---- INSERT WHITE SPACERS BETWEEN EACH PAIR ----
  spacer_cols <- list()
  
  for (i in seq(2, ncol(mat), by = 2)) {
    spacer_name <- paste0("spacer_", i)
    spacer_cols[[spacer_name]] <- rep(NA, nrow(mat))
  }
  
  new_mat <- mat
  for (i in seq_len(ncol(mat))) {
    if (i == 1) {
      new_mat <- mat[, i, drop = FALSE]
    } else {
      new_mat <- cbind(new_mat, mat[, i, drop = FALSE])
    }
    if (i %% 2 == 0) {
      spacer_name <- paste0("spacer_", i)
      new_mat <- cbind(new_mat, spacer_cols[[spacer_name]])
      colnames(new_mat)[ncol(new_mat)] <- spacer_name
    }
  }
  
  # ---- DEFAULT COLOR SCALE RESTORED ----
  pheatmap(
    new_mat,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    na_col = "white",      # spacers remain white
    border_color = NA,
    main = "Average AUC per Protocol × Lineage for Top Regulons"
  )
}

p_2 <- plot_auc_heatmap_lineage(top_regulons, auc_df)

ggsave("~/project/IPSC_2025_Data/Figure5b.tiff",
       plot = p_2,
       device = "tiff",
       width = 10, height = 20, dpi = 300)


#Figure 5c
merged_minus <- subset(merged_IPSC_tf, Protocol == "minus")
merged_minus@misc$SCENIC$RegulonsAUC <- merged_minus@misc$SCENIC$RegulonsAUC[colnames(merged_minus),]
dim(merged_minus)
tf_auc <- merged_minus@misc$SCENIC$RegulonsAUC
dim(tf_auc)
dim(merged_minus)
identical(rownames(tf_auc), colnames(merged_minus))
tf_auc <- tf_auc %>% mutate(across(everything(), ~ replace_na(., 0)))
identical(rownames(tf_auc), colnames(merged_minus))
tf_counts <- t(tf_auc)
merged_minus@assays$TF@counts <- tf_counts
merged_minus@assays$TF@data <- tf_counts
regulon_mat <- t(merged_minus@assays$TF@counts)
pseudotimes <- list(
  dp  = merged_minus$dp_pseudotime,
  up  = merged_minus$up_pseudotime,
  A1  = merged_minus$A1_pseudotime,
  A2  = merged_minus$A2_pseudotime,
  epi = merged_minus$epi_pseudotime,
  crn = merged_minus$crn_pseudotime
)
correlate_regulons <- function(mat, pseudotime) {
  apply(mat, 2, function(reg) {
    ct <- cor.test(reg, pseudotime, method = "pearson")
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
    der_lineage_results %>%
      dplyr::select(Regulon, Lineage, median_p_val_adj, max_p_val_adj, mean_log2FC),
    by = c("Regulon", "Pseudotime" = "Lineage")
  )

cor_results2 <- cor_results2 %>%
  mutate(cor.cor = ifelse(is.na(cor.cor) & !is.na(cor.cor), cor.cor, cor.cor))

threshold_cor <- 0   # or 0.2, 0.3, etc.

cor_results2 <- cor_results2 %>%
  mutate(
    de_dir = case_when(
      !is.na(median_p_val_adj) & median_p_val_adj <= 0.05 & mean_log2FC > 0 ~ "up",
      !is.na(median_p_val_adj) & median_p_val_adj <= 0.05 & mean_log2FC < 0 ~ "down",
      TRUE ~ "ns"
    ),
    is_green = case_when(
      de_dir == "down" & cor.cor < -threshold_cor ~ TRUE,
      de_dir == "up"   & cor.cor > threshold_cor ~ TRUE,
      TRUE ~ FALSE
    ),
    is_red = case_when(
      de_dir == "down" & cor.cor >  threshold_cor ~ TRUE,
      de_dir == "up"   & cor.cor < -threshold_cor ~ TRUE,
      TRUE ~ FALSE
    ),
    color_group = case_when(
      is_green ~ "lineage_aligned",
      is_red   ~ "lineage_opposed",
      TRUE     ~ "neutral"
    )
  )
cor_results_plot <- cor_results2 %>%
  dplyr::filter(!is.na(mean_log2FC))
green_counts <- cor_results2 %>%
  dplyr::filter(is_green) %>%
  count(Pseudotime, name = "n_green")
red_counts <- cor_results2 %>%
  dplyr::filter(is_red) %>%
  count(Pseudotime, name = "n_red")
cor_results_plot <- cor_results2 %>%
  dplyr::filter(!is.na(mean_log2FC))
# Define facet order
facet_order <- c("dp", "up", "A1", "A2", "crn", "epi")

# Apply ordering to Pseudotime
cor_results_plot$Pseudotime <- factor(cor_results_plot$Pseudotime, levels = facet_order)
green_counts$Pseudotime     <- factor(green_counts$Pseudotime,     levels = facet_order)
red_counts$Pseudotime       <- factor(red_counts$Pseudotime,       levels = facet_order)

# Label sets
label_genes_epi_crn       <- c("TCF7L1(+)", "OTX1(+)")
label_genes_dp_up_A1_A2   <- c("POU3F1(+)", "STAT3(+)", "NFIA(+)", "NFIX(+)", "NFIC(+)")

p_3 <- ggplot(cor_results_plot, aes(x = cor.cor, y = mean_log2FC)) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  
  geom_point(aes(color = color_group), alpha = 0.7, size = 1.8) +
  
  scale_color_manual(
    values = c(
      "lineage_aligned" = "forestgreen",
      "lineage_opposed" = "firebrick",
      "neutral"         = "grey70"
    )
  ) +
  
  # ⭐ Force 4×4 facet layout
  facet_wrap(~ Pseudotime, scales = "free", nrow = 4, ncol = 4) +
  
  # dp/up/A1/A2 labels
  geom_point(
    data = subset(cor_results_plot,
                  Regulon %in% label_genes_dp_up_A1_A2 &
                    Pseudotime %in% c("dp", "up", "A1", "A2")),
    color = "black", size = 2.5
  ) +
  geom_text_repel(
    data = subset(cor_results_plot,
                  Regulon %in% label_genes_dp_up_A1_A2 &
                    Pseudotime %in% c("dp", "up", "A1", "A2")),
    aes(label = Regulon),
    size = 4,
    color = "black",
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.4,
    min.segment.length = 0
  ) +
  
  # epi/crn labels
  geom_point(
    data = subset(cor_results_plot,
                  Regulon %in% label_genes_epi_crn &
                    Pseudotime %in% c("crn", "epi")),
    color = "black", size = 2.5
  ) +
  geom_text_repel(
    data = subset(cor_results_plot,
                  Regulon %in% label_genes_epi_crn &
                    Pseudotime %in% c("crn", "epi")),
    aes(label = Regulon),
    size = 4,
    color = "black",
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.4,
    min.segment.length = 0
  ) +
  
  geom_text(
    data = green_counts,
    aes(x = -Inf, y = Inf, label = paste0("aligned = ", n_green)),
    hjust = -0.1,
    vjust = 2.2,
    size = 3.8,
    color = "forestgreen",
    inherit.aes = FALSE
  ) +
  geom_text(
    data = red_counts,
    aes(x = -Inf, y = Inf, label = paste0("opposed = ", n_red)),
    hjust = -0.1,
    vjust = 3.8,
    size = 3.8,
    color = "firebrick",
    inherit.aes = FALSE
  ) +
  
  labs(
    x = "Pearson Correlation with pseudotime",
    y = "Lineage-level Average log2FC (plus vs minus protocol)",
    color = "Regulon class",
    title = "Regulon alignment of cytostatic response with lineage dynamics"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.text.x = element_text(size = 24, face = "bold"),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom"
  )
ggsave("~/project/IPSC_2025_Data/Figure5c.tiff",
       plot = p_3,
       device = "tiff",
       width = 10, height = 10, dpi = 300)