library(Seurat)
library(SeuratExtend)
library(harmony)
library(magrittr)
library(ggplot2)
library(dplyr)

merged <- readRDS("~/project/IPSC_2025_Data/merged_Fetal_IPSC_derived_forebrain")

merged <- NormalizeData(merged) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunHarmony(group.by.vars = "SampleID2") %>% RunUMAP(reduction = "harmony", dims = 1:30) 

#Figure 1b
p2 < - DimPlot2(merged, group.by = "Celltype")
p2

#Figure 1c
features <- c("RG_new_mod", "IPC_EN_new_mod", "ExN_new_mod", "DL_new_mod", "UL_new_mod", "In_new_mod", "CGE_new_mod", "MGE_new_mod", "LGE_new_mod", "Hem_new_mod", "CRN_new_mod", "Epithelial_new_mod", "Astrocyte_new_mod")

#remove fetal only clusters
merged <- subset(merged, Celltype != "OPC" & Celltype != "Vascular/Fibroblast" & Celltype != "Unknown" & Celltype != "Microglia")

desired_order <- c("RG", "IPC_ExN", "DL_ExN", "UL_ExN", "IPC_In", "CGE_In", "MGE_In", "LGE_In", "Hem_RG", "CRN", "Epithelial", "Astrocyte", "OPC")
# Make Consensus2 a factor with the desired order
merged$Celltype <- factor(merged2$Celltype, levels = desired_order)

p3 <- DotPlot(
  merged,
  features = features,
  group.by = "Celltype", split.by = "Sampletype", cols = c("green", "blue"))+ RotatedAxis()
p3