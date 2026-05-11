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



Cortical_lineage <- readRDS("~/project/IPSC_2025_Data/merged_IPSC_derived_pallial_lineages")
Hem_lineage <- readRDS("~/project/IPSC_2025_Data/merged_IPSC_derived_hem_lineages")

ENDS <- c("DL_ExN", "UL_ExN", "A1 Astrocyte", "A2 Astrocyte")
set.seed(1)
Cortical_lineages <- as.SlingshotDataSet(getLineages(
  data           = Cortical_lineage@reductions$umap@cell.embeddings,
  clusterLabels  = Cortical_lineage$Celltype2,
  dist.method    = "slingshot", # It can be: "simple", "scaled.full", "scaled.diag", "slingshot" or "mnn"
  end.clus       = ENDS, # You can also define the ENDS!
  start.clus     = "RG")) # define where to START the trajectories

ENDS <- c("CRN", "Epithelial")
Hem_lineages <- as.SlingshotDataSet(getLineages(
  data           = Hem_lineage@reductions$umap@cell.embeddings,
  clusterLabels  = Hem_lineage$Celltype2,
  dist.method    = "slingshot", # It can be: "simple", "scaled.full", "scaled.diag", "slingshot" or "mnn"
  end.clus       = ENDS, # You can also define the ENDS!
  start.clus     = "Hem_RG")) # define where to START the trajectories


Cortical_lineages
Hem_lineages


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


plot(Cortical_lineage@reductions$umap@cell.embeddings, col = Cortical_lineage$Celltype, pch = 16)
lines(curves_cort, lwd = 2, col = "red")

plot(Hem_lineage@reductions$umap@cell.embeddings, col = Hem_lineage$Celltype, pch = 16)
lines(curves_hem, lwd = 2, col = "red")




pseudotime_cort <- as.data.frame(slingPseudotime(curves_cort, na = FALSE))
pseudotime_cort$Pseudotime <- apply(pseudotime_cort, 1, max)
cellWeights_cort <- slingCurveWeights(curves_cort)
x <- rowMeans(pseudotime_cort)
x <- x / max(x)
o <- order(x)  

pseudotime_hem <- as.data.frame(slingPseudotime(curves_hem, na = FALSE))
pseudotime_hem$Pseudotime <- apply(pseudotime_hem, 1, max)
cellWeights_hem <- slingCurveWeights(curves_hem)
x <- rowMeans(pseudotime_hem)
x <- x / max(x)
o <- order(x) 

# separate pseudotime trajectories by the different cell subtypes in the lineage
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

#Figure 4a
p1 <- Cortical_lineage@meta.data %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=dp_pseudotime)) +
  ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
  coord_equal() +
  scale_color_gradientn(colors=plasma(256), na.value='grey') +
  umap_theme()
p1
p2 <- Cortical_lineage@meta.data %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=up_pseudotime)) +
  ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
  coord_equal() +
  scale_color_gradientn(colors=viridis(256), na.value='grey') +
  umap_theme()
p2
p3 <- Cortical_lineage@meta.data %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=A1_pseudotime)) +
  ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
  coord_equal() +
  scale_color_gradientn(colors=inferno(256), na.value='grey') +
  umap_theme()

p4 <- Cortical_lineage@meta.data %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=A2_pseudotime)) +
  ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
  coord_equal() +
  scale_color_gradientn(colors=mako(256), na.value='grey') +
  umap_theme()

p5 <- Hem_lineage@meta.data %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=crn_pseudotime)) +
  ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
  coord_equal() +
  scale_color_gradientn(colors=inferno(256), na.value='grey') +
  umap_theme()

p6 <- Hem_lineage@meta.data %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=epi_pseudotime)) +
  ggrastr::rasterise(geom_point(size=1), dpi=500, scale=0.75) +
  coord_equal() +
  scale_color_gradientn(colors=mako(256), na.value='grey') +
  umap_theme()

p_1_1 <- (p1 + p2) / (p3 + p4) + plot_annotation(title = "Cortical Lineage By Slingshot Pseudotime")
ggsave("~/project/IPSC_2025_Data/Figure4a_cortical_pseudotime.tiff",
       plot = p_1_1,
       device = "tiff",
       width = 10, height = 10, dpi = 300)

p_1_2 <- (p5 + p6) + plot_annotation(title = "Hem Lineage By Slingshot Pseudotime")
ggsave("~/project/IPSC_2025_Data/Figure4a_hem_pseudotime.tiff",
       plot = p_1_2,
       device = "tiff",
       width = 10, height = 5, dpi = 300)


#Figure4b
# set up the WGCNA experiment in the Seurat object
Cortical_lineage <- SetupForWGCNA(
  Cortical_lineage,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = 'trajectory'
)
# construct metacells 
Cortical_lineage <- MetacellsByGroups(
  Cortical_lineage,
  group.by = c("Celltype2", "SampleID"), # specify the columns in merged@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 10, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'Celltype2' # set the Idents of the metacell seurat object
)
Cortical_lineage <- NormalizeMetacells(Cortical_lineage)



# Use all but one core
cl <- makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)
# setup expression matrix
Cortical_lineage <- SetDatExpr(
  Cortical_lineage, # the name of the group of interest in the group.by column
  group_name = unique(Cortical_lineage@misc$trajectory$wgcna_metacell_obj$Celltype2),
  group.by='Celltype2', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)
# test soft power parameter
Cortical_lineage <- TestSoftPowers(Cortical_lineage)
# construct the co-expression network
Cortical_lineage <- ConstructNetwork(
  Cortical_lineage, 
  tom_name='trajectory', 
  overwrite_tom=TRUE
)
# compute module eigengenes & connectivity
Cortical_lineage <- ModuleEigengenes(Cortical_lineage, group.by.vars=c("SampleID"))
Cortical_lineage <- ModuleConnectivity(Cortical_lineage)
Cortical_lineage <- ResetModuleNames(
  Cortical_lineage,
  new_name = "Pallial-M"
)
stopCluster(cl)





cl <- makeCluster(parallel::detectCores() - 1)
registerDoParallel(cl)
# set up the WGCNA experiment in the Seurat object
Hem_lineage <- SetupForWGCNA(
  Hem_lineage,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = 'trajectory'
)
# construct metacells 
Hem_lineage <- MetacellsByGroups(
  Hem_lineage,
  group.by = c("Celltype2", "SampleID"), # specify the columns in merged@meta.data to group by
  reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 10, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'Celltype2' # set the Idents of the metacell seurat object
)
Hem_lineage <- NormalizeMetacells(Hem_lineage)
# setup expression matrix
Hem_lineage <- SetDatExpr(
  Hem_lineage, # the name of the group of interest in the group.by column
  group_name = unique(Hem_lineage@misc$trajectory$wgcna_metacell_obj$Celltype2),
  group.by='Celltype2', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)
# test soft power parameter
Hem_lineage <- TestSoftPowers(Hem_lineage)
# construct the co-expression network
Hem_lineage <- ConstructNetwork(
  Hem_lineage, 
  tom_name='trajectory', 
  overwrite_tom=TRUE
)
# compute module eigengenes & connectivity
Hem_lineage <- ModuleEigengenes(Hem_lineage, group.by.vars=c("SampleID"))
Hem_lineage <- ModuleConnectivity(Hem_lineage)
Hem_lineage <- ResetModuleNames(
  Hem_lineage,
  new_name = "Hem-M"
)
stopCluster(cl)

MEs_cort <- GetMEs(Cortical_lineage)
modules_cort <- GetModules(Cortical_lineage)
mods_cort <- levels(modules_cort$module)

MEs_hem <- GetMEs(Hem_lineage)
modules_hem <- GetModules(Hem_lineage)
mods_hem <- levels(modules_hem$module)


meta_cort <- Cortical_lineage@meta.data
Cortical_lineage@meta.data <- cbind(meta_cort, MEs_cort)
md_cort <- Cortical_lineage@meta.data

meta_hem <- Hem_lineage@meta.data
Hem_lineage@meta.data <- cbind(meta_hem, MEs_hem)
md_hem <- Hem_lineage@meta.data

split_pseudotime <- function(df, colname) {
  plus_col  <- paste0(colname, "_plus_pseudotime")
  minus_col <- paste0(colname, "_minus_pseudotime")
  
  df[[plus_col]]  <- ifelse(df$Protocol == "plus", df[[colname]], NA)
  df[[minus_col]] <- ifelse(df$Protocol == "minus", df[[colname]], NA)
  
  return(df)
}

# Apply to all four pseudotime columns
for (pt in c("dp_pseudotime", "up_pseudotime", "A1_pseudotime", "A2_pseudotime")) {
  md_cort <- split_pseudotime(md_cort, pt)
}

for (pt in c("crn_pseudotime", "epi_pseudotime")) {
  md_hem <- split_pseudotime(md_hem, pt)
}
#Should be 0
sum(is.na(md_cort$dp_pseudotime_plus_pseudotime)) - sum(is.na(md_cort$dp_pseudotime)) + sum(is.na(md_cort$dp_pseudotime_minus_pseudotime)) - sum(is.na(md_cort$dp_pseudotime)) == nrow(md_cort) - sum(is.na(md_cort$dp_pseudotime))
sum(is.na(md_hem$crn_pseudotime_plus_pseudotime)) - sum(is.na(md_hem$crn_pseudotime)) + sum(is.na(md_hem$crn_pseudotime_minus_pseudotime)) - sum(is.na(md_hem$crn_pseudotime)) == nrow(md_hem) - sum(is.na(md_hem$crn_pseudotime))

Cortical_lineage@meta.data <- md_cort
Hem_lineage@meta.data <- md_hem


#Figure5b
# Define colors and labels once
protocol_colors <- met.brewer("Lakota", n = 2, type = 'discrete')
protocol_labels <- c("minus", "plus")

# dp lineage
p_dp <- PlotModuleTrajectory(
  Cortical_lineage,
  pseudotime_col = c("dp_pseudotime_plus_pseudotime", "dp_pseudotime_minus_pseudotime"),
  group_colors = protocol_colors, ncol = 7
) +
  scale_color_manual(values = protocol_colors, labels = protocol_labels) +
  ggtitle("dp_lineage_eigengene_over_pseudotime") +
  theme(strip.text = element_text(size = 20))

# up lineage
p_up <- PlotModuleTrajectory(
  Cortical_lineage,
  pseudotime_col = c("up_pseudotime_plus_pseudotime", "up_pseudotime_minus_pseudotime"),
  group_colors = protocol_colors, ncol = 7
) +
  scale_color_manual(values = protocol_colors, labels = protocol_labels) +
  ggtitle("up_lineage_eigengene_over_pseudotime") +
  theme(strip.text = element_text(size = 20))
p_up
# A1 lineage
p_A1 <- PlotModuleTrajectory(
  Cortical_lineage,
  pseudotime_col = c("A1_pseudotime_plus_pseudotime", "A1_pseudotime_minus_pseudotime"),
  group_colors = protocol_colors, ncol = 7
) +
  scale_color_manual(values = protocol_colors, labels = protocol_labels) +
  ggtitle("A1_lineage_eigengene_over_pseudotime") +
  theme(strip.text = element_text(size = 20))

# A2 lineage
p_A2 <- PlotModuleTrajectory(
  Cortical_lineage,
  pseudotime_col = c("A2_pseudotime_plus_pseudotime", "A2_pseudotime_minus_pseudotime"),
  group_colors = protocol_colors, ncol = 7
) +
  scale_color_manual(values = protocol_colors, labels = protocol_labels) +
  ggtitle("A2_lineage_eigengene_over_pseudotime") +
  theme(strip.text = element_text(size = 20))


# CRN lineage
p_crn <- PlotModuleTrajectory(
  Hem_lineage,
  pseudotime_col = c("crn_pseudotime_plus_pseudotime", "crn_pseudotime_minus_pseudotime"),
  group_colors = protocol_colors, ncol = 5
) +
  scale_color_manual(values = protocol_colors, labels = protocol_labels) +
  ggtitle("crn_lineage_eigengene_over_pseudotime") +
  theme(strip.text = element_text(size = 20))

# Epi lineage
p_epi <- PlotModuleTrajectory(
  Hem_lineage,
  pseudotime_col = c("epi_pseudotime_plus_pseudotime", "epi_pseudotime_minus_pseudotime"),
  group_colors = protocol_colors, ncol = 5
) +
  scale_color_manual(values = protocol_colors, labels = protocol_labels) +
  ggtitle("epi_lineage_eigengene_over_pseudotime") +
  theme(strip.text = element_text(size = 20))

p_dp <- p_dp + NoLegend()
p_A1 <- p_A1 + NoLegend() 
p_A2 <- p_A2 + NoLegend()  
p_epi <- p_epi + NoLegend()  






p_2_1 <- PlotKMEs(Cortical_lineage, ncol = 10, text_size = 10, plot_widths = c(3, 10)) 
ggsave("~/project/IPSC_2025_Data/Figure_4b_cort.tiff",
       plot = p_2_1,
       device = "tiff",
       width = 50, height = 10, dpi = 300, limitsize = FALSE)
p_2_2 <- PlotKMEs(Hem_lineage, ncol = 10, text_size = 10, plot_widths = c(3, 10)) 
ggsave("~/project/IPSC_2025_Data/Figure_4b_hem.tiff", plot = p_2_2, device = "tiff", width = 40, height = 10, dpi = 300)


p_3_1 <- (p_dp + p_up) / (p_A1 + p_A2) + plot_annotation(title = "Figure5: Effect of 3i On Cortical Lineage Modules)")
ggsave("~/project/IPSC_2025_Data/Figure_4c_cort.tiff",
       plot = p_3_1,
       device = "tiff",
       width = 30, height = 10, dpi = 300)

p_3_2 <- (p_crn + p_epi) + plot_annotation(title = "Figure5: Effect of 3i On Hem Lineage Modules)")
ggsave("~/project/IPSC_2025_Data/Figure_4c_hem.tiff",
       plot = p_3_2,
       device = "tiff",
       width = 20, height = 5, dpi = 300)



Cortical_lineage$dp_lineage <- ifelse(
  !is.na(Cortical_lineage$dp_pseudotime),
  "dp_ExN_lineage",
  NA
)
Cortical_lineage$up_lineage <- ifelse(
  !is.na(Cortical_lineage$up_pseudotime),
  "up_ExN_lineage",
  NA
)
Cortical_lineage$A1_lineage <- ifelse(
  !is.na(Cortical_lineage$A1_pseudotime),
  "A1_Astrocyte_lineage",
  NA
)
Cortical_lineage$A2_lineage <- ifelse(
  !is.na(Cortical_lineage$A2_pseudotime),
  "A2_Astrocyte_lineage",
  NA
)

DMEs_Cortical <- list()
# Identify barcodes for each group
group1 <- Cortical_lineage@meta.data %>%
  subset(dp_lineage == "dp_ExN_lineage" & Protocol == "plus") %>%
  rownames()
group2 <- Cortical_lineage@meta.data %>%
  subset(dp_lineage == "dp_ExN_lineage" & Protocol == "minus") %>%
  rownames()

# Run DME test
dp_DMEs <- FindDMEs(
  Cortical_lineage,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use = "wilcox",
  pseudocount.use = 0.01,
  wgcna_name = "trajectory"
)

# Add cluster info
dp_DMEs$lineage <- rep("DP_ExN_lineage", times = nrow(dp_DMEs))
# Store result
DMEs_Cortical[["DP_ExN_lineage"]] <- dp_DMEs
# Identify barcodes for each group
group1 <- Cortical_lineage@meta.data %>%
  subset(up_lineage == "up_ExN_lineage" & Protocol == "plus") %>%
  rownames()
group2 <- Cortical_lineage@meta.data %>%
  subset(up_lineage == "up_ExN_lineage" & Protocol == "minus") %>%
  rownames()
# Run DME test
up_DMEs <- FindDMEs(
  Cortical_lineage,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use = "wilcox",
  pseudocount.use = 0.01,
  wgcna_name = "trajectory"
)
# Add cluster info
up_DMEs$lineage <- rep("UP_ExN_lineage", times = nrow(up_DMEs))
# Store result
DMEs_Cortical[["UP_ExN_lineage"]] <- up_DMEs
# Identify barcodes for each group
group1 <- Cortical_lineage@meta.data %>%
  subset(A1_lineage == "A1_Astrocyte_lineage" & Protocol == "plus") %>%
  rownames()
group2 <- Cortical_lineage@meta.data %>%
  subset(A1_lineage == "A1_Astrocyte_lineage" & Protocol == "minus") %>%
  rownames()
# Run DME test
A1_DMEs <- FindDMEs(
  Cortical_lineage,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use = "wilcox",
  pseudocount.use = 0.01,
  wgcna_name = "trajectory"
)
# Add cluster info
A1_DMEs$lineage <- rep("A1_Astrocyte_lineage", times = nrow(A1_DMEs))
# Store result
DMEs_Cortical[["A1_Astrocyte_lineage"]] <- A1_DMEs
# Identify barcodes for each group
group1 <- Cortical_lineage@meta.data %>%
  subset(A2_lineage == "A2_Astrocyte_lineage" & Protocol == "plus") %>%
  rownames()
group2 <- Cortical_lineage@meta.data %>%
  subset(A2_lineage == "A2_Astrocyte_lineage" & Protocol == "minus") %>%
  rownames()
# Run DME test
A2_DMEs <- FindDMEs(
  Cortical_lineage,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use = "wilcox",
  pseudocount.use = 0.01,
  wgcna_name = "trajectory"
)
# Add cluster info
A2_DMEs$lineage <- rep("A2_Astrocyte_lineage", times = nrow(A2_DMEs))
# Store result
DMEs_Cortical[["A2_Astrocyte_lineage"]] <- A2_DMEs
DMEs_Cortical <- dplyr::bind_rows(DMEs_Cortical)
# Desired facet order
facet_order <- c(
  "DP_ExN_lineage",
  "UP_ExN_lineage", "A1_Astrocyte_lineage", "A2_Astrocyte_lineage")
# Desired module order
module_order <- c(
  "Pallial-M1","Pallial-M2","Pallial-M3","Pallial-M4","Pallial-M5", "Pallial-M6")
# Prepare data
plot_df <- DMEs_Cortical %>%
  mutate(
    sig_flag = p_val_adj < 0.05,
    # enforce module order
    module = factor(module, levels = module_order),
    # enforce facet order
    lineage = factor(lineage, levels = facet_order)
  )
# Build faceted lollipop plot
p_4_1 <- ggplot(plot_df, aes(x = avg_log2FC, y = module)) +
  # Lollipop stems
  geom_segment(aes(x = 0, xend = avg_log2FC,
                   y = module, yend = module),
               color = "grey60", linewidth = 0.6) +
  # Points for significant modules
  geom_point(aes(color = module),
             size = 3,
             data = subset(plot_df, sig_flag)) +
  # X marks for non-significant modules
  geom_text(data = subset(plot_df, !sig_flag),
            aes(label = "×"),
            color = "black",
            size = 5,
            fontface = "bold") +
  # Facet by Celltype with your custom order
  facet_wrap(~ lineage, scales = "free_y") +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10)
  ) +
  labs(
    title = "DME Lollipop Plot by lineage",
    x = "avg_log2FC",
    y = "module"
  )
p_4_1



Hem_lineage$crn_lineage <- ifelse(
  !is.na(Hem_lineage$crn_pseudotime),
  "crn_lineage",
  NA
)
Hem_lineage$epi_lineage <- ifelse(
  !is.na(Hem_lineage$epi_pseudotime),
  "epi_lineage",
  NA
)
DMEs_Hem <- list()
# Identify barcodes for each group
group1 <- Hem_lineage@meta.data %>%
  subset(crn_lineage == "crn_lineage" & Protocol == "plus") %>%
  rownames()
group2 <- Hem_lineage@meta.data %>%
  subset(crn_lineage == "crn_lineage" & Protocol == "minus") %>%
  rownames()
# Run DME test
crn_DMEs <- FindDMEs(
  Hem_lineage,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use = "wilcox",
  pseudocount.use = 0.01,
  wgcna_name = "trajectory"
)
# Add cluster info
crn_DMEs$lineage <- rep("CRN_lineage", times = nrow(crn_DMEs))
# Store result
DMEs_Hem[["CRN_lineage"]] <- crn_DMEs
group1 <- Hem_lineage@meta.data %>%
  subset(epi_lineage == "epi_lineage" & Protocol == "plus") %>%
  rownames()
group2 <- Hem_lineage@meta.data %>%
  subset(epi_lineage == "epi_lineage" & Protocol == "minus") %>%
  rownames()
# Run DME test
epi_DMEs <- FindDMEs(
  Hem_lineage,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use = "wilcox",
  pseudocount.use = 0.01,
  wgcna_name = "trajectory"
)
# Add cluster info
epi_DMEs$lineage <- rep("Epithelial_lineage", times = nrow(epi_DMEs))
# Store result
DMEs_Hem[["Epithelial_lineage"]] <- epi_DMEs
DMEs_Hem <- dplyr::bind_rows(DMEs_Hem)
# Desired facet order
facet_order <- c(
  "Epithelial_lineage",
  "CRN_lineage")
# Desired module order
module_order <- c(
  "Hem-M1","Hem-M2","Hem-M3","Hem-M4",
  "Hem-M5")
# Prepare data
plot_df <- DMEs_Hem %>%
  mutate(
    sig_flag = p_val_adj < 0.05,
    # enforce module order
    module = factor(module, levels = module_order),
    # enforce facet order
    lineage = factor(lineage, levels = facet_order)
  )
# Build faceted lollipop plot
p_4_2 <- ggplot(plot_df, aes(x = avg_log2FC, y = module)) +
  # Lollipop stems
  geom_segment(aes(x = 0, xend = avg_log2FC,
                   y = module, yend = module),
               color = "grey60", linewidth = 0.6) +
  # Points for significant modules
  geom_point(aes(color = module),
             size = 3,
             data = subset(plot_df, sig_flag)) +
  # X marks for non-significant modules
  geom_text(data = subset(plot_df, !sig_flag),
            aes(label = "×"),
            color = "black",
            size = 5,
            fontface = "bold") +
  # Facet by Celltype with your custom order
  facet_wrap(~ lineage, scales = "free_y") +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10)
  ) +
  labs(
    title = "DME Lollipop Plot by lineage",
    x = "avg_log2FC",
    y = "module"
  )
p_4_2


ggsave("~/project/IPSC_2025_Data/Figure_4d_hem_lineage.tiff",
       plot = p_4_2,
       device = "tiff",
       width = 8, height = 5, dpi = 300)

ggsave("~/project/IPSC_2025_Data/Figure_4d_cort_lineage.tiff",
       plot = p_4_1,
       device = "tiff",
       width = 8, height = 5, dpi = 300)


#Supplementary Figure 6
theme_set(theme_cowplot())
set.seed(12345)
# define the enrichr databases to test
dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023')

# perform enrichment tests
Cortical_lineage <- RunEnrichr(
  Cortical_lineage,
  dbs=dbs,
  max_genes = 100 # use max_genes = Inf to choose all genes
)
Hem_lineage <- RunEnrichr(
  Hem_lineage,
  dbs=dbs,
  max_genes = 100 # use max_genes = Inf to choose all genes
)

# retrieve the output table
#enrich_cort_df <- GetEnrichrTable(Cortical_lineage)
#enrich_hem_df <- GetEnrichrTable(Hem_lineage)

# enrichr dotplot
p_s6_1 <- EnrichrDotPlot(
  Cortical_lineage,
  mods = "all", # use all modules (default)
  database = "GO_Biological_Process_2023", # this must match one of the dbs used previously
  n_terms=5, # number of terms per module
  term_size=8, # font size for the terms
  p_adj = FALSE # show the p-val or adjusted p-val?
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))

ggsave("~/project/IPSC_2025_Data/Supplmental_Figure6_pallial_lineages.tiff",
       plot = p_s6_1,
       device = "tiff",
       width = 10, height = 10, dpi = 300)

# enrichr dotplot
P_s6_2 <-EnrichrDotPlot(
  Hem_lineage,
  mods = "all", # use all modules (default)
  database = "GO_Biological_Process_2023", # this must match one of the dbs used previously
  n_terms=5, # number of terms per module
  term_size=8, # font size for the terms
  p_adj = FALSE # show the p-val or adjusted p-val?
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))

ggsave("~/project/IPSC_2025_Data/Supplmental_Figure6_hem_lineages.tiff",
       plot = p_s6_2,
       device = "tiff",
       width = 10, height = 10, dpi = 300)
