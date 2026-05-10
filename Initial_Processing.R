## Libraries
library(Seurat)
library(data.table)
library(ggplot2)
library(SingleCellExperiment)
library(DropletUtils)
library(scran)
library(SoupX)
library(DoubletFinder)
library(scCustomize)
library(dplyr)
library(ggpmisc)

## ------------------------------------------------------------------
## Load 10x data and souporcell annotations
## ------------------------------------------------------------------

## Read count tables (all sample sets)
minus_3W  <- Read10X("~/project/IPSC_2025_Data/3W_minus/filtered_feature_bc_matrix")
plus_3W   <- Read10X("~/project/IPSC_2025_Data/3W_plus/filtered_feature_bc_matrix")
minus_7W  <- Read10X("~/project/IPSC_2025_Data/7W_minus/filtered_feature_bc_matrix")
plus_7W   <- Read10X("~/project/IPSC_2025_Data/7W_plus/filtered_feature_bc_matrix")
minus_12W <- Read10X("~/project/IPSC_2025_Data/12W_minus/filtered_feature_bc_matrix")
plus_12W  <- Read10X("~/project/IPSC_2025_Data/12W_plus/filtered_feature_bc_matrix")
minus_17W <- Read10X("~/project/IPSC_2025_Data/17W_minus/filtered_feature_bc_matrix")
plus_17W  <- Read10X("~/project/IPSC_2025_Data/17W_plus/filtered_feature_bc_matrix")

minus_7W_2  <- Read10X("~/project/IPSC_2025_Data/7W_minus_2/filtered_feature_bc_matrix")
plus_7W_2   <- Read10X("~/project/IPSC_2025_Data/7W_plus_2/filtered_feature_bc_matrix")
minus_7W_E6   <- Read10X("~/project/IPSC_2025_Data/7W_minus_E6/filtered_feature_bc_matrix")
plus_7W_E6    <- Read10X("~/project/IPSC_2025_Data/7W_plus_E6/filtered_feature_bc_matrix")
minus_12W_E6  <- Read10X("~/project/IPSC_2025_Data/12W_minus_E6/filtered_feature_bc_matrix")
plus_12W_E6   <- Read10X("~/project/IPSC_2025_Data/12W_plus_E6/filtered_feature_bc_matrix")

## Read souporcell labels (ground truth)
Donor_minus_3W_gt  <- fread("~/project/IPSC_2025_Data/3W_minus/clusters.tsv")
Donor_plus_3W_gt   <- fread("~/project/IPSC_2025_Data/3W_plus/clusters.tsv")
Donor_minus_7W_gt  <- fread("~/project/IPSC_2025_Data/7W_minus/souporcell_gt/clusters.tsv")
Donor_plus_7W_gt   <- fread("~/project/IPSC_2025_Data/7W_plus/souporcell_gt/clusters.tsv")
Donor_minus_12W_gt <- fread("~/project/IPSC_2025_Data/12W_minus/souporcell_gt/clusters.tsv")
Donor_plus_12W_gt  <- fread("~/project/IPSC_2025_Data/12W_plus/souporcell_gt/clusters.tsv")
Donor_minus_17W_gt <- fread("~/project/IPSC_2025_Data/17W_minus/clusters.tsv")
Donor_plus_17W_gt  <- fread("~/project/IPSC_2025_Data/17W_plus/clusters.tsv")

Donor_minus_7W_2_gt  <- fread("~/project/IPSC_2025_Data/7W_minus_2/clusters.tsv")
Donor_plus_7W_2_gt   <- fread("~/project/IPSC_2025_Data/7W_plus_2/clusters.tsv")
Donor_minus_7W_E6_gt   <- fread("~/project/IPSC_2025_Data/7W_minus_E6/clusters.tsv")
Donor_plus_7W_E6_gt    <- fread("~/project/IPSC_2025_Data/7W_plus_E6/clusters.tsv")
Donor_minus_12W_E6_gt  <- fread("~/project/IPSC_2025_Data/12W_minus_E6/clusters.tsv")
Donor_plus_12W_E6_gt   <- fread("~/project/IPSC_2025_Data/12W_plus_E6/clusters.tsv")

## Check that barcodes match between matrices and souporcell labels
result <- identical(colnames(minus_3W),  Donor_minus_3W_gt$barcode);  print(result)
result <- identical(colnames(plus_3W),   Donor_plus_3W_gt$barcode);   print(result)
result <- identical(colnames(minus_7W),  Donor_minus_7W_gt$barcode);  print(result)
result <- identical(colnames(plus_7W),   Donor_plus_7W_gt$barcode);   print(result)
result <- identical(colnames(minus_12W), Donor_minus_12W_gt$barcode); print(result)
result <- identical(colnames(plus_12W),  Donor_plus_12W_gt$barcode);  print(result)
result <- identical(colnames(minus_17W), Donor_minus_17W_gt$barcode); print(result)
result <- identical(colnames(plus_17W),  Donor_plus_17W_gt$barcode);  print(result)

result <- identical(colnames(minus_7W_2),  Donor_minus_7W_2_gt$barcode);  print(result)
result <- identical(colnames(plus_7W_2),   Donor_plus_7W_2_gt$barcode);   print(result)
result <- identical(colnames(minus_7W_E6),   Donor_minus_7W_E6_gt$barcode);   print(result)
result <- identical(colnames(plus_7W_E6),    Donor_plus_7W_E6_gt$barcode);    print(result)
result <- identical(colnames(minus_12W_E6),  Donor_minus_12W_E6_gt$barcode);  print(result)
result <- identical(colnames(plus_12W_E6),   Donor_plus_12W_E6_gt$barcode);   print(result)

## ------------------------------------------------------------------
## Create Seurat objects and attach souporcell metadata
## ------------------------------------------------------------------

minus_3W <- CreateSeuratObject(minus_3W)
minus_3W$gt_doublet_status <- Donor_minus_3W_gt$status
minus_3W$gt_line           <- Donor_minus_3W_gt$assignment
table(minus_3W$gt_doublet_status, minus_3W$gt_line)
table(Donor_minus_3W_gt$status, Donor_minus_3W_gt$assignment)
minus_3W <- SetIdent(minus_3W, value = "gt_line")
minus_3W <- RenameIdents(minus_3W, `0` = "KOLF2.1", `1` = "O2C3", `2` = "JHC1")
minus_3W$gt_line <- Idents(minus_3W)
table(minus_3W$gt_doublet_status, minus_3W$gt_line)
table(minus_3W$gt_doublet_status, Idents(minus_3W))

plus_3W <- CreateSeuratObject(plus_3W)
plus_3W$gt_doublet_status <- Donor_plus_3W_gt$status
plus_3W$gt_line           <- Donor_plus_3W_gt$assignment
table(plus_3W$gt_doublet_status, plus_3W$gt_line)
table(Donor_plus_3W_gt$status, Donor_plus_3W_gt$assignment)
plus_3W <- SetIdent(plus_3W, value = "gt_line")
plus_3W <- RenameIdents(plus_3W, `0` = "KOLF2.1", `1` = "O2C3", `2` = "JHC1")
plus_3W$gt_line <- Idents(plus_3W)
table(plus_3W$gt_doublet_status, plus_3W$gt_line)
table(plus_3W$gt_doublet_status, Idents(plus_3W))

minus_7W <- CreateSeuratObject(minus_7W)
minus_7W$gt_doublet_status <- Donor_minus_7W_gt$status
minus_7W$gt_line           <- Donor_minus_7W_gt$assignment
table(minus_7W$gt_doublet_status, minus_7W$gt_line)
table(Donor_minus_7W_gt$status, Donor_minus_7W_gt$assignment)
minus_7W <- SetIdent(minus_7W, value = "gt_line")
minus_7W <- RenameIdents(minus_7W, `0` = "KOLF2.1", `1` = "O2C3")
minus_7W$gt_line <- Idents(minus_7W)
table(minus_7W$gt_doublet_status, minus_7W$gt_line)
table(minus_7W$gt_doublet_status, Idents(minus_7W))

plus_7W <- CreateSeuratObject(plus_7W)
plus_7W$gt_doublet_status <- Donor_plus_7W_gt$status
plus_7W$gt_line           <- Donor_plus_7W_gt$assignment
table(plus_7W$gt_doublet_status, plus_7W$gt_line)
table(Donor_plus_7W_gt$status, Donor_plus_7W_gt$assignment)
plus_7W <- SetIdent(plus_7W, value = "gt_line")
plus_7W <- RenameIdents(plus_7W, `0` = "KOLF2.1", `1` = "O2C3")
plus_7W$gt_line <- Idents(plus_7W)
table(plus_7W$gt_doublet_status, plus_7W$gt_line)
table(plus_7W$gt_doublet_status, Idents(plus_7W))

minus_12W <- CreateSeuratObject(minus_12W)
minus_12W$gt_doublet_status <- Donor_minus_12W_gt$status
minus_12W$gt_line           <- Donor_minus_12W_gt$assignment
table(minus_12W$gt_doublet_status, minus_12W$gt_line)
table(Donor_minus_12W_gt$status, Donor_minus_12W_gt$assignment)
minus_12W <- SetIdent(minus_12W, value = "gt_line")
minus_12W <- RenameIdents(minus_12W, `0` = "KOLF2.1", `1` = "O2C3")
minus_12W$gt_line <- Idents(minus_12W)
table(minus_12W$gt_doublet_status, minus_12W$gt_line)
table(minus_12W$gt_doublet_status, Idents(minus_12W))

plus_12W <- CreateSeuratObject(plus_12W)
plus_12W$gt_doublet_status <- Donor_plus_12W_gt$status
plus_12W$gt_line           <- Donor_plus_12W_gt$assignment
table(plus_12W$gt_doublet_status, plus_12W$gt_line)
table(Donor_plus_12W_gt$status, Donor_plus_12W_gt$assignment)
plus_12W <- SetIdent(plus_12W, value = "gt_line")
plus_12W <- RenameIdents(plus_12W, `0` = "KOLF2.1", `1` = "O2C3")
plus_12W$gt_line <- Idents(plus_12W)
table(plus_12W$gt_doublet_status, plus_12W$gt_line)
table(plus_12W$gt_doublet_status, Idents(plus_12W))

minus_17W <- CreateSeuratObject(minus_17W)
minus_17W$gt_doublet_status <- Donor_minus_17W_gt$status
minus_17W$gt_line           <- Donor_minus_17W_gt$assignment
table(minus_17W$gt_doublet_status, minus_17W$gt_line)
table(Donor_minus_17W_gt$status, Donor_minus_17W_gt$assignment)
minus_17W <- SetIdent(minus_17W, value = "gt_line")
minus_17W <- RenameIdents(minus_17W, `0` = "KOLF2.1", `1` = "O2C3")
minus_17W$gt_line <- Idents(minus_17W)
table(minus_17W$gt_doublet_status, minus_17W$gt_line)
table(minus_17W$gt_doublet_status, Idents(minus_17W))

plus_17W <- CreateSeuratObject(plus_17W)
plus_17W$gt_doublet_status <- Donor_plus_17W_gt$status
plus_17W$gt_line           <- Donor_plus_17W_gt$assignment
table(plus_17W$gt_doublet_status, plus_17W$gt_line)
table(Donor_plus_17W_gt$status, Donor_plus_17W_gt$assignment)
plus_17W <- SetIdent(plus_17W, value = "gt_line")
plus_17W <- RenameIdents(plus_17W, `0` = "KOLF2.1", `1` = "O2C3")
plus_17W$gt_line <- Idents(plus_17W)
table(plus_17W$gt_doublet_status, plus_17W$gt_line)
table(plus_17W$gt_doublet_status, Idents(plus_17W))

minus_7W_2 <- CreateSeuratObject(minus_7W_2)
minus_7W_2$gt_doublet_status <- Donor_minus_7W_2_gt$status
minus_7W_2$gt_line           <- Donor_minus_7W_2_gt$assignment
table(minus_7W_2$gt_doublet_status, minus_7W_2$gt_line)
table(Donor_minus_7W_2_gt$status, Donor_minus_7W_2_gt$assignment)
minus_7W_2 <- SetIdent(minus_7W_2, value = "gt_line")
minus_7W_2 <- RenameIdents(minus_7W_2, `0` = "KOLF2.1", `1` = "O2C3", `2` = "JHC1")
minus_7W_2$gt_line <- Idents(minus_7W_2)
table(minus_7W_2$gt_doublet_status, minus_7W_2$gt_line)
table(minus_7W_2$gt_doublet_status, Idents(minus_7W_2))

plus_7W_2 <- CreateSeuratObject(plus_7W_2)
plus_7W_2$gt_doublet_status <- Donor_plus_7W_2_gt$status
plus_7W_2$gt_line           <- Donor_plus_7W_2_gt$assignment
table(plus_7W_2$gt_doublet_status, plus_7W_2$gt_line)
table(Donor_plus_7W_2_gt$status, Donor_plus_7W_2_gt$assignment)
plus_7W_2 <- SetIdent(plus_7W_2, value = "gt_line")
plus_7W_2 <- RenameIdents(plus_7W_2, `0` = "KOLF2.1", `1` = "O2C3", `2` = "JHC1")
plus_7W_2$gt_line <- Idents(plus_7W_2)
table(plus_7W_2$gt_doublet_status, plus_7W_2$gt_line)
table(plus_7W_2$gt_doublet_status, Idents(plus_7W_2))

minus_7W_E6 <- CreateSeuratObject(minus_7W_E6)
minus_7W_E6$gt_doublet_status <- Donor_minus_7W_E6_gt$status
minus_7W_E6$gt_line           <- Donor_minus_7W_E6_gt$assignment
table(minus_7W_E6$gt_doublet_status, minus_7W_E6$gt_line)
table(Donor_minus_7W_E6_gt$status, Donor_minus_7W_E6_gt$assignment)
minus_7W_E6 <- SetIdent(minus_7W_E6, value = "gt_line")
minus_7W_E6 <- RenameIdents(minus_7W_E6, `0` = "KOLF2.1", `1` = "O2C3", `2` = "JHC1")
minus_7W_E6$gt_line <- Idents(minus_7W_E6)
table(minus_7W_E6$gt_doublet_status, minus_7W_E6$gt_line)
table(minus_7W_E6$gt_doublet_status, Idents(minus_7W_E6))

plus_7W_E6 <- CreateSeuratObject(plus_7W_E6)
plus_7W_E6$gt_doublet_status <- Donor_plus_7W_E6_gt$status
plus_7W_E6$gt_line           <- Donor_plus_7W_E6_gt$assignment
table(plus_7W_E6$gt_doublet_status, plus_7W_E6$gt_line)
table(Donor_plus_7W_E6_gt$status, Donor_plus_7W_E6_gt$assignment)
plus_7W_E6 <- SetIdent(plus_7W_E6, value = "gt_line")
plus_7W_E6 <- RenameIdents(plus_7W_E6, `0` = "KOLF2.1", `1` = "O2C3", `2` = "JHC1")
plus_7W_E6$gt_line <- Idents(plus_7W_E6)
table(plus_7W_E6$gt_doublet_status, plus_7W_E6$gt_line)
table(plus_7W_E6$gt_doublet_status, Idents(plus_7W_E6))

minus_12W_E6 <- CreateSeuratObject(minus_12W_E6)
minus_12W_E6$gt_doublet_status <- Donor_minus_12W_E6_gt$status
minus_12W_E6$gt_line           <- Donor_minus_12W_E6_gt$assignment
table(minus_12W_E6$gt_doublet_status, minus_12W_E6$gt_line)
table(Donor_minus_12W_E6_gt$status, Donor_minus_12W_E6_gt$assignment)
minus_12W_E6 <- SetIdent(minus_12W_E6, value = "gt_line")
minus_12W_E6 <- RenameIdents(minus_12W_E6, `0` = "KOLF2.1", `1` = "O2C3", `2` = "JHC1")
minus_12W_E6$gt_line <- Idents(minus_12W_E6)
table(minus_12W_E6$gt_doublet_status, minus_12W_E6$gt_line)
table(minus_12W_E6$gt_doublet_status, Idents(minus_12W_E6))

plus_12W_E6 <- CreateSeuratObject(plus_12W_E6)
plus_12W_E6$gt_doublet_status <- Donor_plus_12W_E6_gt$status
plus_12W_E6$gt_line           <- Donor_plus_12W_E6_gt$assignment
table(plus_12W_E6$gt_doublet_status, plus_12W_E6$gt_line)
table(Donor_plus_12W_E6_gt$status, Donor_plus_12W_E6_gt$assignment)
plus_12W_E6 <- SetIdent(plus_12W_E6, value = "gt_line")
plus_12W_E6 <- RenameIdents(plus_12W_E6, `0` = "KOLF2.1", `1` = "O2C3", `2` = "JHC1")
plus_12W_E6$gt_line <- Idents(plus_12W_E6)
table(plus_12W_E6$gt_doublet_status, plus_12W_E6$gt_line)
table(plus_12W_E6$gt_doublet_status, Idents(plus_12W_E6))

## ------------------------------------------------------------------
## Nuclear fraction and empty droplet removal
## ------------------------------------------------------------------
nf <- readRDS("~/project/IPSC_2025_Data/3W_minus/rds_nf")
minus_3W$nuclear_fraction <- nf$nuclear_fraction
meta <- as.data.frame(minus_3W@meta.data)
p <- ggplot(meta, aes(x = nuclear_fraction))
setwd("~/project/IPSC_2025_Data/3W_minus/")
tiff("minus_3W_nf_density_curve.tiff", width = 500, height = 500)
p + geom_density()
dev.off()

nf <- readRDS("~/project/IPSC_2025_Data/3W_plus/rds_nf")
plus_3W$nuclear_fraction <- nf$nuclear_fraction
meta <- as.data.frame(plus_3W@meta.data)
p <- ggplot(meta, aes(x = nuclear_fraction))
setwd("~/project/IPSC_2025_Data/3W_plus/")
tiff("plus_3W_nf_density_curve.tiff", width = 500, height = 500)
p + geom_density()
dev.off()

nf <- readRDS("~/project/IPSC_2025_Data/7W_minus/rds_nf")
minus_7W$nuclear_fraction <- nf$nuclear_fraction
meta <- as.data.frame(minus_7W@meta.data)
p <- ggplot(meta, aes(x = nuclear_fraction))
setwd("~/project/IPSC_2025_Data/7W_minus/")
tiff("minus_7W_nf_density_curve.tiff", width = 500, height = 500)
p + geom_density()
dev.off()

nf <- readRDS("~/project/IPSC_2025_Data/7W_plus/rds_nf")
plus_7W$nuclear_fraction <- nf$nuclear_fraction
meta <- as.data.frame(plus_7W@meta.data)
p <- ggplot(meta, aes(x = nuclear_fraction))
setwd("~/project/IPSC_2025_Data/7W_plus/")
tiff("plus_7W_nf_density_curve.tiff", width = 500, height = 500)
p + geom_density()
dev.off()

nf <- readRDS("~/project/IPSC_2025_Data/12W_minus/rds_nf")
minus_12W$nuclear_fraction <- nf$nuclear_fraction
meta <- as.data.frame(minus_12W@meta.data)
p <- ggplot(meta, aes(x = nuclear_fraction))
setwd("~/project/IPSC_2025_Data/12W_minus/")
tiff("minus_12W_nf_density_curve.tiff", width = 500, height = 500)
p + geom_density()
dev.off()

nf <- readRDS("~/project/IPSC_2025_Data/12W_plus/rds_nf")
plus_12W$nuclear_fraction <- nf$nuclear_fraction
meta <- as.data.frame(plus_12W@meta.data)
p <- ggplot(meta, aes(x = nuclear_fraction))
setwd("~/project/IPSC_2025_Data/12W_plus/")
tiff("plus_12W_nf_density_curve.tiff", width = 500, height = 500)
p + geom_density()
dev.off()

nf <- readRDS("~/project/IPSC_2025_Data/17W_minus/rds_nf")
minus_17W$nuclear_fraction <- nf$nuclear_fraction
meta <- as.data.frame(minus_17W@meta.data)
p <- ggplot(meta, aes(x = nuclear_fraction))
setwd("~/project/IPSC_2025_Data/17W_minus/")
tiff("minus_17W_nf_density_curve.tiff", width = 500, height = 500)
p + geom_density()
dev.off()

nf <- readRDS("~/project/IPSC_2025_Data/17W_plus/rds_nf")
plus_17W$nuclear_fraction <- nf$nuclear_fraction
meta <- as.data.frame(plus_17W@meta.data)
p <- ggplot(meta, aes(x = nuclear_fraction))
setwd("~/project/IPSC_2025_Data/17W_plus/")
tiff("plus_17W_nf_density_curve.tiff", width = 500, height = 500)
p + geom_density()
dev.off()

nf <- readRDS("~/project/IPSC_2025_Data/7W_minus_2/rds_nf")
minus_7W_2$nuclear_fraction <- nf$nuclear_fraction
meta <- as.data.frame(minus_7W_2@meta.data)
p <- ggplot(meta, aes(x = nuclear_fraction))
setwd("~/project/IPSC_2025_Data/7W_minus_2/")
tiff("minus_7W_2_nf_density_curve.tiff", width = 500, height = 500)
p + geom_density()
dev.off()

nf <- readRDS("~/project/IPSC_2025_Data/7W_plus_2/rds_nf")
plus_7W_2$nuclear_fraction <- nf$nuclear_fraction
meta <- as.data.frame(plus_7W_2@meta.data)
p <- ggplot(meta, aes(x = nuclear_fraction))
setwd("~/project/IPSC_2025_Data/7W_plus_2/")
tiff("plus_7W_2_nf_density_curve.tiff", width = 500, height = 500)
p + geom_density()
dev.off()


#Add Metadata

###############################################################################
## Helper: assign neural induction media automatically
###############################################################################
assign_media <- function(obj_name) {
  if (grepl("E6", obj_name)) "E6" else "KSR"
}

###############################################################################
## 3W (D21)
###############################################################################
minus_3W$Age      <- rep("D21",  ncol(minus_3W))
minus_3W$Protocol <- rep("minus", ncol(minus_3W))
minus_3W$neural_induction_media <- rep(assign_media("minus_3W"), ncol(minus_3W))
minus_3W$SampleID <- paste(minus_3W$Age, minus_3W$Protocol, minus_3W$gt_line,
                           minus_3W$neural_induction_media, sep = "_")

plus_3W$Age       <- rep("D21",  ncol(plus_3W))
plus_3W$Protocol  <- rep("plus", ncol(plus_3W))
plus_3W$neural_induction_media <- rep(assign_media("plus_3W"), ncol(plus_3W))
plus_3W$SampleID  <- paste(plus_3W$Age, plus_3W$Protocol, plus_3W$gt_line,
                           plus_3W$neural_induction_media, sep = "_")

###############################################################################
## 7W (D49)
###############################################################################
minus_7W$Age      <- rep("D49",  ncol(minus_7W))
minus_7W$Protocol <- rep("minus", ncol(minus_7W))
minus_7W$neural_induction_media <- rep(assign_media("minus_7W"), ncol(minus_7W))
minus_7W$SampleID <- paste(minus_7W$Age, minus_7W$Protocol, minus_7W$gt_line,
                           minus_7W$neural_induction_media, sep = "_")

plus_7W$Age       <- rep("D49",  ncol(plus_7W))
plus_7W$Protocol  <- rep("plus", ncol(plus_7W))
plus_7W$neural_induction_media <- rep(assign_media("plus_7W"), ncol(plus_7W))
plus_7W$SampleID  <- paste(plus_7W$Age, plus_7W$Protocol, plus_7W$gt_line,
                           plus_7W$neural_induction_media, sep = "_")

###############################################################################
## 7W_2 (D49) — duplicated samples get "_2"
###############################################################################
minus_7W_2$Age      <- rep("D49",  ncol(minus_7W_2))
minus_7W_2$Protocol <- rep("minus", ncol(minus_7W_2))
minus_7W_2$neural_induction_media <- rep(assign_media("minus_7W_2"), ncol(minus_7W_2))
minus_7W_2$SampleID <- paste(minus_7W_2$Age, minus_7W_2$Protocol, minus_7W_2$gt_line,
                             minus_7W_2$neural_induction_media, sep = "_")
minus_7W_2$SampleID <- paste0(minus_7W_2$SampleID, "_2")

plus_7W_2$Age       <- rep("D49",  ncol(plus_7W_2))
plus_7W_2$Protocol  <- rep("plus", ncol(plus_7W_2))
plus_7W_2$neural_induction_media <- rep(assign_media("plus_7W_2"), ncol(plus_7W_2))
plus_7W_2$SampleID  <- paste(plus_7W_2$Age, plus_7W_2$Protocol, plus_7W_2$gt_line,
                             plus_7W_2$neural_induction_media, sep = "_")
plus_7W_2$SampleID  <- paste0(plus_7W_2$SampleID, "_2")

###############################################################################
## 7W_E6 (D49)
###############################################################################
minus_7W_E6$Age      <- rep("D49",  ncol(minus_7W_E6))
minus_7W_E6$Protocol <- rep("minus", ncol(minus_7W_E6))
minus_7W_E6$neural_induction_media <- rep(assign_media("minus_7W_E6"), ncol(minus_7W_E6))
minus_7W_E6$SampleID <- paste(minus_7W_E6$Age, minus_7W_E6$Protocol, minus_7W_E6$gt_line,
                              minus_7W_E6$neural_induction_media, sep = "_")

plus_7W_E6$Age       <- rep("D49",  ncol(plus_7W_E6))
plus_7W_E6$Protocol  <- rep("plus", ncol(plus_7W_E6))
plus_7W_E6$neural_induction_media <- rep(assign_media("plus_7W_E6"), ncol(plus_7W_E6))
plus_7W_E6$SampleID  <- paste(plus_7W_E6$Age, plus_7W_E6$Protocol, plus_7W_E6$gt_line,
                              plus_7W_E6$neural_induction_media, sep = "_")

###############################################################################
## 12W (D84)
###############################################################################
minus_12W$Age      <- rep("D84",  ncol(minus_12W))
minus_12W$Protocol <- rep("minus", ncol(minus_12W))
minus_12W$neural_induction_media <- rep(assign_media("minus_12W"), ncol(minus_12W))
minus_12W$SampleID <- paste(minus_12W$Age, minus_12W$Protocol, minus_12W$gt_line,
                            minus_12W$neural_induction_media, sep = "_")

plus_12W$Age       <- rep("D84",  ncol(plus_12W))
plus_12W$Protocol  <- rep("plus", ncol(plus_12W))
plus_12W$neural_induction_media <- rep(assign_media("plus_12W"), ncol(plus_12W))
plus_12W$SampleID  <- paste(plus_12W$Age, plus_12W$Protocol, plus_12W$gt_line,
                            plus_12W$neural_induction_media, sep = "_")

###############################################################################
## 12W_E6 (D84)
###############################################################################
minus_12W_E6$Age      <- rep("D84",  ncol(minus_12W_E6))
minus_12W_E6$Protocol <- rep("minus", ncol(minus_12W_E6))
minus_12W_E6$neural_induction_media <- rep(assign_media("minus_12W_E6"), ncol(minus_12W_E6))
minus_12W_E6$SampleID <- paste(minus_12W_E6$Age, minus_12W_E6$Protocol, minus_12W_E6$gt_line,
                               minus_12W_E6$neural_induction_media, sep = "_")

plus_12W_E6$Age       <- rep("D84",  ncol(plus_12W_E6))
plus_12W_E6$Protocol  <- rep("plus", ncol(plus_12W_E6))
plus_12W_E6$neural_induction_media <- rep(assign_media("plus_12W_E6"), ncol(plus_12W_E6))
plus_12W_E6$SampleID  <- paste(plus_12W_E6$Age, plus_12W_E6$Protocol, plus_12W_E6$gt_line,
                               plus_12W_E6$neural_induction_media, sep = "_")

###############################################################################
## 17W (D119)
###############################################################################
minus_17W$Age      <- rep("D119", ncol(minus_17W))
minus_17W$Protocol <- rep("minus", ncol(minus_17W))
minus_17W$neural_induction_media <- rep(assign_media("minus_17W"), ncol(minus_17W))
minus_17W$SampleID <- paste(minus_17W$Age, minus_17W$Protocol, minus_17W$gt_line,
                            minus_17W$neural_induction_media, sep = "_")

plus_17W$Age       <- rep("D119", ncol(plus_17W))
plus_17W$Protocol  <- rep("plus", ncol(plus_17W))
plus_17W$neural_induction_media <- rep(assign_media("plus_17W"), ncol(plus_17W))
plus_17W$SampleID  <- paste(plus_17W$Age, plus_17W$Protocol, plus_17W$gt_line,
                            plus_17W$neural_induction_media, sep = "_")


## ------------------------------------------------------------------
## Ambient RNA estimation and filtering (all samples)
## ------------------------------------------------------------------

## minus_3W
a <- colnames(subset(minus_3W, nuclear_fraction < 0.3))
minus_3W_sce <- as.SingleCellExperiment(minus_3W)
x <- colnames(minus_3W)
b <- match(a, x)
ambient_profile_genes_minus_3W <- ambientProfileEmpty(counts(minus_3W_sce), known.empty = b, good.turing = FALSE, round = FALSE)
minus_3W_sce_agg <- aggregateAcrossCells(minus_3W_sce, ids = DataFrame(sample = minus_3W_sce$SampleID))
max.ambient_profile_genes_minus_3W <- ambientContribMaximum(counts(minus_3W_sce_agg), ambient_profile_genes_minus_3W, mode = "proportion")
contamination <- rowMeans(max.ambient_profile_genes_minus_3W, na.rm = TRUE)
top_genes <- names(sort(contamination, decreasing = TRUE)[1:50])
minus_3W <- NormalizeData(minus_3W)
minus_3W <- AddModuleScore(minus_3W, features = list(top_genes), name = "Amb_Genes")
meta <- minus_3W@meta.data
p <- ggplot(meta, aes(x = Amb_Genes1)) + geom_density() + scale_x_continuous(limits = c(-0.005, 0.005))
pb <- ggplot_build(p)
p + stat_valleys(
  data = pb[["data"]][[1]],
  aes(x = x, y = density),
  colour = "red",
  size = 3,
  geom = "text"
)
setwd("~/project/IPSC_2025_Data/3W_minus/")
tiff("GW3_minus_amb_gene_module_density_curve.tiff", width = 500, height = 500)
p + geom_vline(xintercept = -0.0001076, size = 0.5, color = "blue")
dev.off()
minus_3W_filt <- subset(minus_3W, Amb_Genes1 < -0.0001076 & nuclear_fraction > 0.4)
minus_3W_filt[["percent.mt"]] <- PercentageFeatureSet(minus_3W_filt, pattern = "^MT-")
minus_3W_filt$high_subsets_Mt_percent <- isOutlier(minus_3W_filt$percent.mt, type = "higher", min.diff = 0.5)
minus_3W_filt <- subset(minus_3W_filt, high_subsets_Mt_percent == "FALSE")

## plus_3W
a <- colnames(subset(plus_3W, nuclear_fraction < 0.3))
plus_3W_sce <- as.SingleCellExperiment(plus_3W)
x <- colnames(plus_3W)
b <- match(a, x)
ambient_profile_genes_plus_3W <- ambientProfileEmpty(counts(plus_3W_sce), known.empty = b, good.turing = FALSE, round = FALSE)
plus_3W_sce_agg <- aggregateAcrossCells(plus_3W_sce, ids = DataFrame(sample = plus_3W_sce$SampleID))
max.ambient_profile_genes_plus_3W <- ambientContribMaximum(counts(plus_3W_sce_agg), ambient_profile_genes_plus_3W, mode = "proportion")
contamination <- rowMeans(max.ambient_profile_genes_plus_3W, na.rm = TRUE)
top_genes <- names(sort(contamination, decreasing = TRUE)[1:25])
plus_3W <- NormalizeData(plus_3W)
plus_3W <- AddModuleScore(plus_3W, features = list(top_genes), name = "Amb_Genes")
meta <- plus_3W@meta.data
p <- ggplot(meta, aes(x = Amb_Genes1)) + geom_density() + scale_x_continuous(limits = c(-0.005, 0.005))
pb <- ggplot_build(p)
p + stat_valleys(
  data = pb[["data"]][[1]],
  aes(x = x, y = density),
  colour = "red",
  size = 3,
  geom = "text"
)
setwd("~/project/IPSC_2025_Data/3W_plus/")
tiff("GW3_plus_amb_gene_module_density_curve.tiff", width = 500, height = 500)
p + geom_vline(xintercept = -0.0001272, size = 0.5, color = "blue")
dev.off()
plus_3W_filt <- subset(plus_3W, Amb_Genes1 < -0.0001272 & nuclear_fraction > 0.4)
plus_3W_filt[["percent.mt"]] <- PercentageFeatureSet(plus_3W_filt, pattern = "^MT-")
plus_3W_filt$high_subsets_Mt_percent <- isOutlier(plus_3W_filt$percent.mt, type = "higher", min.diff = 0.5)
plus_3W_filt <- subset(plus_3W_filt, high_subsets_Mt_percent == "FALSE")

## minus_7W_2
a <- colnames(subset(minus_7W_2, nuclear_fraction < 0.2))
minus_7W_2_sce <- as.SingleCellExperiment(minus_7W_2)
x <- colnames(minus_7W_2)
b <- match(a, x)
ambient_profile_genes_minus_7W_2 <- ambientProfileEmpty(counts(minus_7W_2_sce), known.empty = b, good.turing = FALSE, round = FALSE)
minus_7W_2_sce_agg <- aggregateAcrossCells(minus_7W_2_sce, ids = DataFrame(sample = minus_7W_2_sce$SampleID))
max.ambient_profile_genes_minus_7W_2 <- ambientContribMaximum(counts(minus_7W_2_sce_agg), ambient_profile_genes_minus_7W_2, mode = "proportion")
contamination <- rowMeans(max.ambient_profile_genes_minus_7W_2, na.rm = TRUE)
top_genes <- names(sort(contamination, decreasing = TRUE)[1:18])
minus_7W_2 <- NormalizeData(minus_7W_2)
minus_7W_2 <- AddModuleScore(minus_7W_2, features = list(top_genes), name = "Amb_Genes")
meta <- minus_7W_2@meta.data
p <- ggplot(meta, aes(x = Amb_Genes1)) + geom_density() + scale_x_continuous(limits = c(-0.005, 0.005))
pb <- ggplot_build(p)
p + stat_valleys(
  data = pb[["data"]][[1]],
  aes(x = x, y = density),
  colour = "red",
  size = 3,
  geom = "text"
)
setwd("~/project/IPSC_2025_Data/7W_minus_KSR/")
tiff("7W_minus_KSR_amb_gene_module_density_curve.tiff", width = 500, height = 500)
p + geom_vline(xintercept = -0.0002838, size = 0.5, color = "blue")
dev.off()
minus_7W_2_filt <- subset(minus_7W_2, Amb_Genes1 < -0.0002838 & nuclear_fraction > 0.3)
minus_7W_2_filt[["percent.mt"]] <- PercentageFeatureSet(minus_7W_2_filt, pattern = "^MT-")
minus_7W_2_filt$high_subsets_Mt_percent <- isOutlier(minus_7W_2_filt$percent.mt, type = "higher", min.diff = 0.5)
minus_7W_2_filt <- subset(minus_7W_2_filt, high_subsets_Mt_percent == "FALSE")

## plus_7W_2
a <- colnames(subset(plus_7W_2, nuclear_fraction < 0.2))
plus_7W_2_sce <- as.SingleCellExperiment(plus_7W_2)
x <- colnames(plus_7W_2)
b <- match(a, x)
ambient_profile_genes_plus_7W_2 <- ambientProfileEmpty(counts(plus_7W_2_sce), known.empty = b, good.turing = FALSE, round = FALSE)
plus_7W_2_sce_agg <- aggregateAcrossCells(plus_7W_2_sce, ids = DataFrame(sample = plus_7W_2_sce$SampleID))
max.ambient_profile_genes_plus_7W_2 <- ambientContribMaximum(counts(plus_7W_2_sce_agg), ambient_profile_genes_plus_7W_2, mode = "proportion")
contamination <- rowMeans(max.ambient_profile_genes_plus_7W_2, na.rm = TRUE)
top_genes <- names(sort(contamination, decreasing = TRUE)[1:44])
plus_7W_2 <- NormalizeData(plus_7W_2)
plus_7W_2 <- AddModuleScore(plus_7W_2, features = list(top_genes), name = "Amb_Genes")
meta <- plus_7W_2@meta.data
p <- ggplot(meta, aes(x = Amb_Genes1)) + geom_density() + scale_x_continuous(limits = c(-0.005, 0.005))
pb <- ggplot_build(p)
p + stat_valleys(
  data = pb[["data"]][[1]],
  aes(x = x, y = density),
  colour = "red",
  size = 3,
  geom = "text"
)
setwd("~/project/IPSC_2025_Data/7W_plus_KSR/")
tiff("7W_plus_KSR_amb_gene_module_density_curve.tiff", width = 500, height = 500)
p + geom_vline(xintercept = -0.0001076, size = 0.5, color = "blue")
dev.off()
plus_7W_2_filt <- subset(plus_7W_2, Amb_Genes1 < -0.0001272 & nuclear_fraction > 0.2)
plus_7W_2_filt[["percent.mt"]] <- PercentageFeatureSet(plus_7W_2_filt, pattern = "^MT-")
plus_7W_2_filt$high_subsets_Mt_percent <- isOutlier(plus_7W_2_filt$percent.mt, type = "higher", min.diff = 0.5)
plus_7W_2_filt <- subset(plus_7W_2_filt, high_subsets_Mt_percent == "FALSE")

## minus_17W
a <- colnames(subset(minus_17W, nuclear_fraction < 0.3))
minus_17W_sce <- as.SingleCellExperiment(minus_17W)
x <- colnames(minus_17W)
b <- match(a, x)
ambient_profile_genes_minus_17W <- ambientProfileEmpty(counts(minus_17W_sce), known.empty = b, good.turing = FALSE, round = FALSE)
minus_17W_sce_agg <- aggregateAcrossCells(minus_17W_sce, ids = DataFrame(sample = minus_17W_sce$SampleID))
max.ambient_profile_genes_minus_17W <- ambientContribMaximum(counts(minus_17W_sce_agg), ambient_profile_genes_minus_17W, mode = "proportion")
contamination <- rowMeans(max.ambient_profile_genes_minus_17W, na.rm = TRUE)
top_genes <- names(sort(contamination, decreasing = TRUE)[1:36])
minus_17W <- AddModuleScore(minus_17W, features = list(top_genes), name = "Amb_Genes")
meta <- minus_17W@meta.data
p <- ggplot(meta, aes(x = Amb_Genes1)) + geom_density() + scale_x_continuous(limits = c(-0.005, 0.005))
pb <- ggplot_build(p)
p + stat_valleys(
  data = pb[["data"]][[1]],
  aes(x = x, y = density),
  colour = "red",
  size = 3,
  geom = "text"
)
setwd("~/project/IPSC_2025_Data/17W_minus/")
tiff("GW17_minus_amb_gene_module_density_curve.tiff", width = 500, height = 500)
p + geom_vline(xintercept = -0.0001663, size = 0.5, color = "blue")
dev.off()
minus_17W_filt <- subset(minus_17W, Amb_Genes1 < -0.0001468 & nuclear_fraction > 0.4)
minus_17W_filt[["percent.mt"]] <- PercentageFeatureSet(minus_17W_filt, pattern = "^MT-")
minus_17W_filt$high_subsets_Mt_percent <- isOutlier(minus_17W_filt$percent.mt, type = "higher", min.diff = 0.5)
minus_17W_filt <- subset(minus_17W_filt, high_subsets_Mt_percent == "FALSE")

## plus_17W
a <- colnames(subset(plus_17W, nuclear_fraction < 0.3))
plus_17W_sce <- as.SingleCellExperiment(plus_17W)
x <- colnames(plus_17W)
b <- match(a, x)
ambient_profile_genes_plus_17W <- ambientProfileEmpty(counts(plus_17W_sce), known.empty = b, good.turing = FALSE, round = FALSE)
plus_17W_sce_agg <- aggregateAcrossCells(plus_17W_sce, ids = DataFrame(sample = plus_17W_sce$SampleID))
max.ambient_profile_genes_plus_17W <- ambientContribMaximum(counts(plus_17W_sce_agg), ambient_profile_genes_plus_17W, mode = "proportion")
contamination <- rowMeans(max.ambient_profile_genes_plus_17W, na.rm = TRUE)
top_genes <- names(sort(contamination, decreasing = TRUE)[1:95])
plus_17W <- AddModuleScore(plus_17W, features = list(top_genes), name = "Amb_Genes")
meta <- plus_17W@meta.data
p <- ggplot(meta, aes(x = Amb_Genes1)) + geom_density() + scale_x_continuous(limits = c(-0.005, 0.005))
pb <- ggplot_build(p)
p + stat_valleys(
  data = pb[["data"]][[1]],
  aes(x = x, y = density),
  colour = "red",
  size = 3,
  geom = "text"
)
setwd("~/project/IPSC_2025_Data/17W_plus/")
tiff("GW17_plus_amb_gene_module_density_curve.tiff", width = 1000, height = 500)
p + geom_vline(xintercept = -8.806e-05, size = 0.5, color = "blue")
dev.off()
plus_17W_filt <- subset(plus_17W, Amb_Genes1 < -8.806e-05 & nuclear_fraction > 0.4)
plus_17W_filt[["percent.mt"]] <- PercentageFeatureSet(plus_17W_filt, pattern = "^MT-")
plus_17W_filt$high_subsets_Mt_percent <- isOutlier(plus_17W_filt$percent.mt, type = "higher", min.diff = 0.5)
plus_17W_filt <- subset(plus_17W_filt, high_subsets_Mt_percent == "FALSE")



## Simple mitochondrial filtering for Other sets
minus_7W_E6[["percent.mt"]] <- PercentageFeatureSet(minus_7W_E6, pattern = "^MT-")
minus_7W_E6$high_subsets_Mt_percent <- isOutlier(minus_7W_E6$percent.mt, type = "higher", min.diff = 0.5)
minus_7W_E6_filt <- subset(minus_7W_E6, high_subsets_Mt_percent == "FALSE")

plus_7W_E6[["percent.mt"]] <- PercentageFeatureSet(plus_7W_E6, pattern = "^MT-")
plus_7W_E6$high_subsets_Mt_percent <- isOutlier(plus_7W_E6$percent.mt, type = "higher", min.diff = 0.5)
plus_7W_E6_filt <- subset(plus_7W_E6, high_subsets_Mt_percent == "FALSE")

minus_12W_E6[["percent.mt"]] <- PercentageFeatureSet(minus_12W_E6, pattern = "^MT-")
minus_12W_E6$high_subsets_Mt_percent <- isOutlier(minus_12W_E6$percent.mt, type = "higher", min.diff = 0.5)
minus_12W_E6_filt <- subset(minus_12W_E6, high_subsets_Mt_percent == "FALSE")

plus_12W_E6[["percent.mt"]] <- PercentageFeatureSet(plus_12W_E6, pattern = "^MT-")
plus_12W_E6$high_subsets_Mt_percent <- isOutlier(plus_12W_E6$percent.mt, type = "higher", min.diff = 0.5)
plus_12W_E6_filt <- subset(plus_12W_E6, high_subsets_Mt_percent == "FALSE")

minus_7W_2[["percent.mt"]] <- PercentageFeatureSet(minus_7W_2, pattern = "^MT-")
minus_7W_2$high_subsets_Mt_percent <- isOutlier(minus_7W_2$percent.mt, type = "higher", min.diff = 0.5)
minus_7W_2_filt <- subset(minus_7W_2, high_subsets_Mt_percent == "FALSE")

plus_7W_2[["percent.mt"]] <- PercentageFeatureSet(plus_7W_2, pattern = "^MT-")
plus_7W_2$high_subsets_Mt_percent <- isOutlier(plus_7W_2$percent.mt, type = "higher", min.diff = 0.5)
plus_7W_2_filt <- subset(plus_7W_2, high_subsets_Mt_percent == "FALSE")

## ------------------------------------------------------------------
## DoubletFinder (all filtered objects)
## ------------------------------------------------------------------

## minus_3W_filt
minus_3W_filt <- NormalizeData(minus_3W_filt)
minus_3W_filt <- FindVariableFeatures(minus_3W_filt, selection.method = "vst", nfeatures = 2000)
minus_3W_filt <- ScaleData(minus_3W_filt)
minus_3W_filt <- RunPCA(minus_3W_filt)
minus_3W_filt <- RunUMAP(minus_3W_filt, dims = 1:30)
minus_3W_filt <- FindNeighbors(minus_3W_filt, dims = 1:30)
minus_3W_filt <- FindClusters(minus_3W_filt)
sweep.res.list_minus_3W_filt <- paramSweep(minus_3W_filt, PCs = 1:30, sct = FALSE)
sweep.stats_minus_3W_filt <- summarizeSweep(sweep.res.list_minus_3W_filt, GT = FALSE)
bcmvn_minus_3W_filt <- find.pK(sweep.stats_minus_3W_filt)
a <- (nrow(minus_3W_filt@meta.data) / 1000) * 0.004
homotypic.prop <- modelHomotypic(minus_3W_filt@meta.data$seurat_clusters)
nExp_poi <- round(a * nrow(minus_3W_filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
optimal.pk <- bcmvn_minus_3W_filt %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
minus_3W_filt <- doubletFinder(minus_3W_filt, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

## plus_3W_filt
plus_3W_filt <- NormalizeData(plus_3W_filt)
plus_3W_filt <- FindVariableFeatures(plus_3W_filt, selection.method = "vst", nfeatures = 2000)
plus_3W_filt <- ScaleData(plus_3W_filt)
plus_3W_filt <- RunPCA(plus_3W_filt)
plus_3W_filt <- RunUMAP(plus_3W_filt, dims = 1:30)
plus_3W_filt <- FindNeighbors(plus_3W_filt, dims = 1:30)
plus_3W_filt <- FindClusters(plus_3W_filt)
sweep.res.list_plus_3W_filt <- paramSweep(plus_3W_filt, PCs = 1:30, sct = FALSE)
sweep.stats_plus_3W_filt <- summarizeSweep(sweep.res.list_plus_3W_filt, GT = FALSE)
bcmvn_plus_3W_filt <- find.pK(sweep.stats_plus_3W_filt)
a <- (nrow(plus_3W_filt@meta.data) / 1000) * 0.004
homotypic.prop <- modelHomotypic(plus_3W_filt@meta.data$seurat_clusters)
nExp_poi <- round(a * nrow(plus_3W_filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
optimal.pk <- bcmvn_plus_3W_filt %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
plus_3W_filt <- doubletFinder(plus_3W_filt, PCs = 1:30, pN = 0.25, pK = 0.3, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

## minus_7W_filt
minus_7W_filt <- NormalizeData(minus_7W_filt)
minus_7W_filt <- FindVariableFeatures(minus_7W_filt, selection.method = "vst", nfeatures = 2000)
minus_7W_filt <- ScaleData(minus_7W_filt)
minus_7W_filt <- RunPCA(minus_7W_filt)
minus_7W_filt <- RunUMAP(minus_7W_filt, dims = 1:30)
minus_7W_filt <- FindNeighbors(minus_7W_filt, dims = 1:30)
minus_7W_filt <- FindClusters(minus_7W_filt)
sweep.res.list_minus_7W_filt <- paramSweep(minus_7W_filt, PCs = 1:30, sct = FALSE)
sweep.stats_minus_7W_filt <- summarizeSweep(sweep.res.list_minus_7W_filt, GT = FALSE)
bcmvn_minus_7W_filt <- find.pK(sweep.stats_minus_7W_filt)
homotypic.prop <- modelHomotypic(minus_7W_filt@meta.data$seurat_clusters)
nExp_poi <- round(0.08 * nrow(minus_7W_filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
optimal.pk <- bcmvn_minus_7W_filt %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
minus_7W_filt <- doubletFinder(minus_7W_filt, PCs = 1:30, pN = 0.25, pK = 0.26, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

## plus_7W_filt
plus_7W_filt <- NormalizeData(plus_7W_filt)
plus_7W_filt <- FindVariableFeatures(plus_7W_filt, selection.method = "vst", nfeatures = 2000)
plus_7W_filt <- ScaleData(plus_7W_filt)
plus_7W_filt <- RunPCA(plus_7W_filt)
plus_7W_filt <- RunUMAP(plus_7W_filt, dims = 1:30)
plus_7W_filt <- FindNeighbors(plus_7W_filt, dims = 1:30)
plus_7W_filt <- FindClusters(plus_7W_filt)
sweep.res.list_plus_7W_filt <- paramSweep(plus_7W_filt, PCs = 1:30, sct = FALSE)
sweep.stats_plus_7W_filt <- summarizeSweep(sweep.res.list_plus_7W_filt, GT = FALSE)
bcmvn_plus_7W_filt <- find.pK(sweep.stats_plus_7W_filt)
homotypic.prop <- modelHomotypic(plus_7W_filt@meta.data$seurat_clusters)
nExp_poi <- round(0.08 * nrow(plus_7W_filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
optimal.pk <- bcmvn_plus_7W_filt %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
plus_7W_filt <- doubletFinder(plus_7W_filt, PCs = 1:30, pN = 0.25, pK = 0.24, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

## minus_12W_filt
minus_12W_filt <- NormalizeData(minus_12W_filt)
minus_12W_filt <- FindVariableFeatures(minus_12W_filt, selection.method = "vst", nfeatures = 2000)
minus_12W_filt <- ScaleData(minus_12W_filt)
minus_12W_filt <- RunPCA(minus_12W_filt)
minus_12W_filt <- RunUMAP(minus_12W_filt, dims = 1:30)
minus_12W_filt <- FindNeighbors(minus_12W_filt, dims = 1:30)
minus_12W_filt <- FindClusters(minus_12W_filt)
sweep.res.list_minus_12W_filt <- paramSweep(minus_12W_filt, PCs = 1:30, sct = FALSE)
sweep.stats_minus_12W_filt <- summarizeSweep(sweep.res.list_minus_12W_filt, GT = FALSE)
bcmvn_minus_12W_filt <- find.pK(sweep.stats_minus_12W_filt)
homotypic.prop <- modelHomotypic(minus_12W_filt@meta.data$seurat_clusters)
nExp_poi <- round(0.072692 * nrow(minus_12W_filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
optimal.pk <- bcmvn_minus_12W_filt %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
minus_12W_filt <- doubletFinder(minus_12W_filt, PCs = 1:30, pN = 0.25, pK = 0.28, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

## plus_12W_filt
plus_12W_filt <- NormalizeData(plus_12W_filt)
plus_12W_filt <- FindVariableFeatures(plus_12W_filt, selection.method = "vst", nfeatures = 2000)
plus_12W_filt <- ScaleData(plus_12W_filt)
plus_12W_filt <- RunPCA(plus_12W_filt)
plus_12W_filt <- RunUMAP(plus_12W_filt, dims = 1:30)
plus_12W_filt <- FindNeighbors(plus_12W_filt, dims = 1:30)
plus_12W_filt <- FindClusters(plus_12W_filt)
sweep.res.list_plus_12W_filt <- paramSweep(plus_12W_filt, PCs = 1:30, sct = FALSE)
sweep.stats_plus_12W_filt <- summarizeSweep(sweep.res.list_plus_12W_filt, GT = FALSE)
bcmvn_plus_12W_filt <- find.pK(sweep.stats_plus_12W_filt)
homotypic.prop <- modelHomotypic(plus_12W_filt@meta.data$seurat_clusters)
nExp_poi <- round(0.08 * nrow(plus_12W_filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
optimal.pk <- bcmvn_plus_12W_filt %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
plus_12W_filt <- doubletFinder(plus_12W_filt, PCs = 1:30, pN = 0.25, pK = 0.08, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

## minus_17W_filt
minus_17W_filt <- NormalizeData(minus_17W_filt)
minus_17W_filt <- FindVariableFeatures(minus_17W_filt, selection.method = "vst", nfeatures = 2000)
minus_17W_filt <- ScaleData(minus_17W_filt)
minus_17W_filt <- RunPCA(minus_17W_filt)
minus_17W_filt <- RunUMAP(minus_17W_filt, dims = 1:30)
minus_17W_filt <- FindNeighbors(minus_17W_filt, dims = 1:30)
minus_17W_filt <- FindClusters(minus_17W_filt)
sweep.res.list_minus_17W_filt <- paramSweep(minus_17W_filt, PCs = 1:30, sct = FALSE)
sweep.stats_minus_17W_filt <- summarizeSweep(sweep.res.list_minus_17W_filt, GT = FALSE)
bcmvn_minus_17W_filt <- find.pK(sweep.stats_minus_17W_filt)
a <- (nrow(minus_17W_filt@meta.data) / 1000) * 0.004
homotypic.prop <- modelHomotypic(minus_17W_filt@meta.data$seurat_clusters)
nExp_poi <- round(a * nrow(minus_17W_filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
optimal.pk <- bcmvn_minus_17W_filt %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
minus_17W_filt <- doubletFinder(minus_17W_filt, PCs = 1:30, pN = 0.25, pK = 0.3, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

## plus_17W_filt
plus_17W_filt <- NormalizeData(plus_17W_filt)
plus_17W_filt <- FindVariableFeatures(plus_17W_filt, selection.method = "vst", nfeatures = 2000)
plus_17W_filt <- ScaleData(plus_17W_filt)
plus_17W_filt <- RunPCA(plus_17W_filt)
plus_17W_filt <- RunUMAP(plus_17W_filt, dims = 1:30)
plus_17W_filt <- FindNeighbors(plus_17W_filt, dims = 1:30)
plus_17W_filt <- FindClusters(plus_17W_filt)
sweep.res.list_plus_17W_filt <- paramSweep(plus_17W_filt, PCs = 1:30, sct = FALSE)
sweep.stats_plus_17W_filt <- summarizeSweep(sweep.res.list_plus_17W_filt, GT = FALSE)
bcmvn_plus_17W_filt <- find.pK(sweep.stats_plus_17W_filt)
a <- (nrow(plus_17W_filt@meta.data) / 1000) * 0.004
homotypic.prop <- modelHomotypic(plus_17W_filt@meta.data$seurat_clusters)
nExp_poi <- round(a * nrow(plus_17W_filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
optimal.pk <- bcmvn_plus_17W_filt %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
plus_17W_filt <- doubletFinder(plus_17W_filt, PCs = 1:30, pN = 0.25, pK = 0.27, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

## minus_7W_2_filt
minus_7W_2_filt <- NormalizeData(minus_7W_2_filt)
minus_7W_2_filt <- FindVariableFeatures(minus_7W_2_filt, selection.method = "vst", nfeatures = 2000)
minus_7W_2_filt <- ScaleData(minus_7W_2_filt)
minus_7W_2_filt <- RunPCA(minus_7W_2_filt)
minus_7W_2_filt <- RunUMAP(minus_7W_2_filt, dims = 1:30)
minus_7W_2_filt <- FindNeighbors(minus_7W_2_filt, dims = 1:30)
minus_7W_2_filt <- FindClusters(minus_7W_2_filt)
sweep.res.list_minus_7W_2_filt <- paramSweep(minus_7W_2_filt, PCs = 1:30, sct = FALSE)
sweep.stats_minus_7W_2_filt <- summarizeSweep(sweep.res.list_minus_7W_2_filt, GT = FALSE)
bcmvn_minus_7W_2_filt <- find.pK(sweep.stats_minus_7W_2_filt)
a <- (nrow(minus_7W_2_filt@meta.data) / 1000) * 0.004
homotypic.prop <- modelHomotypic(minus_7W_2_filt@meta.data$seurat_clusters)
nExp_poi <- round(a * nrow(minus_7W_2_filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
optimal.pk <- bcmvn_minus_7W_2_filt %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
minus_7W_2_filt <- doubletFinder(minus_7W_2_filt, PCs = 1:30, pN = 0.25, pK = 0.13, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

## plus_7W_2_filt
plus_7W_2_filt <- NormalizeData(plus_7W_2_filt)
plus_7W_2_filt <- FindVariableFeatures(plus_7W_2_filt, selection.method = "vst", nfeatures = 2000)
plus_7W_2_filt <- ScaleData(plus_7W_2_filt)
plus_7W_2_filt <- RunPCA(plus_7W_2_filt)
plus_7W_2_filt <- RunUMAP(plus_7W_2_filt, dims = 1:30)
plus_7W_2_filt <- FindNeighbors(plus_7W_2_filt, dims = 1:30)
plus_7W_2_filt <- FindClusters(plus_7W_2_filt)
sweep.res.list_plus_7W_2_filt <- paramSweep(plus_7W_2_filt, PCs = 1:30, sct = FALSE)
sweep.stats_plus_7W_2_filt <- summarizeSweep(sweep.res.list_plus_7W_2_filt, GT = FALSE)
bcmvn_plus_7W_2_filt <- find.pK(sweep.stats_plus_7W_2_filt)
a <- (nrow(plus_7W_2_filt@meta.data) / 1000) * 0.004
homotypic.prop <- modelHomotypic(plus_7W_2_filt@meta.data$seurat_clusters)
nExp_poi <- round(a * nrow(plus_7W_2_filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
optimal.pk <- bcmvn_plus_7W_2_filt %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
plus_7W_2_filt <- doubletFinder(plus_7W_2_filt, PCs = 1:30, pN = 0.25, pK = 0.3, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

## minus_7W_E6_filt
minus_7W_E6_filt <- NormalizeData(minus_7W_E6_filt)
minus_7W_E6_filt <- FindVariableFeatures(minus_7W_E6_filt, selection.method = "vst", nfeatures = 2000)
minus_7W_E6_filt <- ScaleData(minus_7W_E6_filt)
minus_7W_E6_filt <- RunPCA(minus_7W_E6_filt)
minus_7W_E6_filt <- RunUMAP(minus_7W_E6_filt, dims = 1:30)
minus_7W_E6_filt <- FindNeighbors(minus_7W_E6_filt, dims = 1:30)
minus_7W_E6_filt <- FindClusters(minus_7W_E6_filt)
sweep.res.list_minus_7W_E6_filt <- paramSweep(minus_7W_E6_filt, PCs = 1:30, sct = FALSE)
sweep.stats_minus_7W_E6_filt <- summarizeSweep(sweep.res.list_minus_7W_E6_filt, GT = FALSE)
bcmvn_minus_7W_E6_filt <- find.pK(sweep.stats_minus_7W_E6_filt)
a <- (nrow(minus_7W_E6_filt@meta.data) / 1000) * 0.004
homotypic.prop <- modelHomotypic(minus_7W_E6_filt@meta.data$seurat_clusters)
nExp_poi <- round(a * nrow(minus_7W_E6_filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
optimal.pk <- bcmvn_minus_7W_E6_filt %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
minus_7W_E6_filt <- doubletFinder(minus_7W_E6_filt, PCs = 1:30, pN = 0.25, pK = 0.22, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

## plus_7W_E6_filt
plus_7W_E6_filt <- NormalizeData(plus_7W_E6_filt)
plus_7W_E6_filt <- FindVariableFeatures(plus_7W_E6_filt, selection.method = "vst", nfeatures = 2000)
plus_7W_E6_filt <- ScaleData(plus_7W_E6_filt)
plus_7W_E6_filt <- RunPCA(plus_7W_E6_filt)
plus_7W_E6_filt <- RunUMAP(plus_7W_E6_filt, dims = 1:30)
plus_7W_E6_filt <- FindNeighbors(plus_7W_E6_filt, dims = 1:30)
plus_7W_E6_filt <- FindClusters(plus_7W_E6_filt)
sweep.res.list_plus_7W_E6_filt <- paramSweep(plus_7W_E6_filt, PCs = 1:30, sct = FALSE)
sweep.stats_plus_7W_E6_filt <- summarizeSweep(sweep.res.list_plus_7W_E6_filt, GT = FALSE)
bcmvn_plus_7W_E6_filt <- find.pK(sweep.stats_plus_7W_E6_filt)
a <- (nrow(plus_7W_E6_filt@meta.data) / 1000) * 0.004
homotypic.prop <- modelHomotypic(plus_7W_E6_filt@meta.data$seurat_clusters)
nExp_poi <- round(a * nrow(plus_7W_E6_filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
optimal.pk <- bcmvn_plus_7W_E6_filt %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
plus_7W_E6_filt <- doubletFinder(plus_7W_E6_filt, PCs = 1:30, pN = 0.25, pK = 0.03, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

## minus_12W_E6_filt
minus_12W_E6_filt <- NormalizeData(minus_12W_E6_filt)
minus_12W_E6_filt <- FindVariableFeatures(minus_12W_E6_filt, selection.method = "vst", nfeatures = 2000)
minus_12W_E6_filt <- ScaleData(minus_12W_E6_filt)
minus_12W_E6_filt <- RunPCA(minus_12W_E6_filt)
minus_12W_E6_filt <- RunUMAP(minus_12W_E6_filt, dims = 1:30)
minus_12W_E6_filt <- FindNeighbors(minus_12W_E6_filt, dims = 1:30)
minus_12W_E6_filt <- FindClusters(minus_12W_E6_filt)
sweep.res.list_minus_12W_E6_filt <- paramSweep(minus_12W_E6_filt, PCs = 1:30, sct = FALSE)
sweep.stats_minus_12W_E6_filt <- summarizeSweep(sweep.res.list_minus_12W_E6_filt, GT = FALSE)
bcmvn_minus_12W_E6_filt <- find.pK(sweep.stats_minus_12W_E6_filt)
a <- (nrow(minus_12W_E6_filt@meta.data) / 1000) * 0.004
homotypic.prop <- modelHomotypic(minus_12W_E6_filt@meta.data$seurat_clusters)
nExp_poi <- round(a * nrow(minus_12W_E6_filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
optimal.pk <- bcmvn_minus_12W_E6_filt %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
minus_12W_E6_filt <- doubletFinder(minus_12W_E6_filt, PCs = 1:30, pN = 0.25, pK = 0.16, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

## plus_12W_E6_filt
plus_12W_E6_filt <- NormalizeData(plus_12W_E6_filt)
plus_12W_E6_filt <- FindVariableFeatures(plus_12W_E6_filt, selection.method = "vst", nfeatures = 2000)
plus_12W_E6_filt <- ScaleData(plus_12W_E6_filt)
plus_12W_E6_filt <- RunPCA(plus_12W_E6_filt)
plus_12W_E6_filt <- RunUMAP(plus_12W_E6_filt, dims = 1:30)
plus_12W_E6_filt <- FindNeighbors(plus_12W_E6_filt, dims = 1:30)
plus_12W_E6_filt <- FindClusters(plus_12W_E6_filt)
sweep.res.list_plus_12W_E6_filt <- paramSweep(plus_12W_E6_filt, PCs = 1:30, sct = FALSE)
sweep.stats_plus_12W_E6_filt <- summarizeSweep(sweep.res.list_plus_12W_E6_filt, GT = FALSE)
bcmvn_plus_12W_E6_filt <- find.pK(sweep.stats_plus_12W_E6_filt)
a <- (nrow(plus_12W_E6_filt@meta.data) / 1000) * 0.004
homotypic.prop <- modelHomotypic(plus_12W_E6_filt@meta.data$seurat_clusters)
nExp_poi <- round(a * nrow(plus_12W_E6_filt@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
optimal.pk <- bcmvn_plus_12W_E6_filt %>% dplyr::filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
plus_12W_E6_filt <- doubletFinder(plus_12W_E6_filt, PCs = 1:30, pN = 0.25, pK = 0.3, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)

## ------------------------------------------------------------------
## SoupX
## ------------------------------------------------------------------

## minus_3W_filt SoupX
minus_3W_filt <- readRDS("~/project/IPSC_2025_Data/3W_minus/3W_minus_filt_after_df_seu")
counts <- minus_3W_filt@assays$RNA@layers$counts
colnames(counts) <- colnames(minus_3W_filt)
rownames(counts) <- rownames(minus_3W_filt)
write10xCounts(x = counts, path = "~/project/IPSC_2025_Data/3W_minus/filtered_feature_bc_matrix2")
toc <- Seurat::Read10X("~/project/IPSC_2025_Data/3W_minus/filtered_feature_bc_matrix2")
tod <- Seurat::Read10X("~/project/IPSC_2025_Data/3W_minus/raw_feature_bc_matrix")
common_genes <- intersect(rownames(toc), rownames(tod))
toc <- toc[common_genes, ]
tod <- tod[common_genes, ]
sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
sc <- estimateSoup(sc)
minus_3W_filt <- NormalizeData(minus_3W_filt)
minus_3W_filt <- FindVariableFeatures(minus_3W_filt)
minus_3W_filt <- ScaleData(minus_3W_filt)
minus_3W_filt <- RunPCA(minus_3W_filt)
minus_3W_filt <- RunUMAP(minus_3W_filt, dims = 1:30)
minus_3W_filt <- FindNeighbors(minus_3W_filt, dims = 1:30)
minus_3W_filt <- FindClusters(minus_3W_filt, resolution = 0.4)
minus_3W_filt$Clusters <- Idents(minus_3W_filt)
sc <- setClusters(sc, setNames(minus_3W_filt$Clusters, rownames(minus_3W_filt@meta.data)))
sc <- autoEstCont(sc)
setwd("~/project/IPSC_2025_Data/3W_minus/")
tiff("minus_3W_autoEstcon.tiff", width = 500, height = 500)
autoEstCont(sc)
dev.off()
out <- adjustCounts(sc, roundToInt = TRUE)
seu2 <- CreateSeuratObject(out)
seu2 <- NormalizeData(seu2)
seu2 <- FindVariableFeatures(seu2)
seu2 <- ScaleData(seu2)
seu2 <- RunPCA(seu2)
seu2 <- RunUMAP(seu2, dims = 1:30)
seu2 <- FindNeighbors(seu2, dims = 1:30)
seu2 <- FindClusters(seu2)
seu2@meta.data <- minus_3W_filt@meta.data
setwd("~/project/IPSC_2025_Data/3W_minus/")
tiff("3W_minus_before_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(minus_3W_filt, features = c("GLI3", "LMX1A", "DCX", "CFAP299"))
dev.off()
tiff("3W_minus_after_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(seu2, features = c("GLI3", "LMX1A", "DCX", "CFAP299"))
dev.off()
saveRDS(seu2, "~/project/IPSC_2025_Data/3W_minus/3W_minus_filt_after_df_and_soupX_auto_seu")

## plus_3W_filt SoupX
plus_3W_filt <- readRDS("~/project/IPSC_2025_Data/3W_plus/3W_plus_filt_after_df_seu")
counts <- plus_3W_filt@assays$RNA@layers$counts
colnames(counts) <- colnames(plus_3W_filt)
rownames(counts) <- rownames(plus_3W_filt)
write10xCounts(x = counts, path = "~/project/IPSC_2025_Data/3W_plus/filtered_feature_bc_matrix2")
toc <- Seurat::Read10X("~/project/IPSC_2025_Data/3W_plus/filtered_feature_bc_matrix2")
tod <- Seurat::Read10X("~/project/IPSC_2025_Data/3W_plus/raw_feature_bc_matrix")
common_genes <- intersect(rownames(toc), rownames(tod))
toc <- toc[common_genes, ]
tod <- tod[common_genes, ]
sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
sc <- estimateSoup(sc)
plus_3W_filt <- NormalizeData(plus_3W_filt)
plus_3W_filt <- FindVariableFeatures(plus_3W_filt)
plus_3W_filt <- ScaleData(plus_3W_filt)
plus_3W_filt <- RunPCA(plus_3W_filt)
plus_3W_filt <- RunUMAP(plus_3W_filt, dims = 1:30)
plus_3W_filt <- FindNeighbors(plus_3W_filt, dims = 1:30)
plus_3W_filt <- FindClusters(plus_3W_filt, resolution = 0.4)
plus_3W_filt$Clusters <- Idents(plus_3W_filt)
sc <- setClusters(sc, setNames(plus_3W_filt$Clusters, rownames(plus_3W_filt@meta.data)))
sc <- autoEstCont(sc)
setwd("~/project/IPSC_2025_Data/3W_plus/")
tiff("plus_3W_autoEstcon.tiff", width = 500, height = 500)
autoEstCont(sc)
dev.off()
out <- adjustCounts(sc, roundToInt = TRUE)
seu2_plus <- CreateSeuratObject(out)
seu2_plus <- NormalizeData(seu2_plus)
seu2_plus <- FindVariableFeatures(seu2_plus)
seu2_plus <- ScaleData(seu2_plus)
seu2_plus <- RunPCA(seu2_plus)
seu2_plus <- RunUMAP(seu2_plus, dims = 1:30)
seu2_plus <- FindNeighbors(seu2_plus, dims = 1:30)
seu2_plus <- FindClusters(seu2_plus)
seu2_plus@meta.data <- plus_3W_filt@meta.data
tiff("3W_plus_before_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(plus_3W_filt, features = c("GLI3", "LMX1A", "DCX", "CFAP299"))
dev.off()
tiff("3W_plus_after_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(seu2_plus, features = c("GLI3", "LMX1A", "DCX", "CFAP299"))
dev.off()
saveRDS(seu2_plus, "~/project/IPSC_2025_Data/3W_plus/3W_plus_filt_after_df_and_soupX_auto_seu")

minus_7W_filt <- readRDS("~/project/IPSC_2025_Data/7W_minus/7W_minus_filt_after_df_seu")
counts <- minus_7W_filt@assays$RNA@layers$counts
colnames(counts) <- colnames(minus_7W_filt)
rownames(counts) <- rownames(minus_7W_filt)
write10xCounts(x = counts, path = "~/project/IPSC_2025_Data/7W_minus/filtered_feature_bc_matrix2")
toc <- Seurat::Read10X("~/project/IPSC_2025_Data/7W_minus/filtered_feature_bc_matrix2")
tod <- Seurat::Read10X("~/project/IPSC_2025_Data/7W_minus/raw_feature_bc_matrix")
common_genes <- intersect(rownames(toc), rownames(tod))
toc <- toc[common_genes, ]
tod <- tod[common_genes, ]
sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
sc <- estimateSoup(sc)
minus_7W_filt <- NormalizeData(minus_7W_filt)
minus_7W_filt <- FindVariableFeatures(minus_7W_filt)
minus_7W_filt <- ScaleData(minus_7W_filt)
minus_7W_filt <- RunPCA(minus_7W_filt)
minus_7W_filt <- RunUMAP(minus_7W_filt, dims = 1:30)
minus_7W_filt <- FindNeighbors(minus_7W_filt, dims = 1:30)
minus_7W_filt <- FindClusters(minus_7W_filt, resolution = 0.8)
minus_7W_filt$Clusters <- Idents(minus_7W_filt)
sc <- setClusters(sc, setNames(minus_7W_filt$Clusters, rownames(minus_7W_filt@meta.data)))
sc <- autoEstCont(sc)
setwd("~/project/IPSC_2025_Data/7W_minus/")
tiff("7W_minus_autoEstcon.tiff", width = 500, height = 500)
autoEstCont(sc)
dev.off()
out <- adjustCounts(sc, roundToInt = TRUE)
seu2 <- CreateSeuratObject(out)
seu2 <- NormalizeData(seu2)
seu2 <- FindVariableFeatures(seu2)
seu2 <- ScaleData(seu2)
seu2 <- RunPCA(seu2)
seu2 <- RunUMAP(seu2, dims = 1:30)
seu2 <- FindNeighbors(seu2, dims = 1:30)
seu2 <- FindClusters(seu2)
seu2@meta.data <- minus_7W_filt@meta.data
tiff("7W_minus_before_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(minus_7W_filt, features = c("GLI3", "LMX1A", "TBR1", "GAD2"))
dev.off()
tiff("7W_minus_after_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(seu2, features = c("GLI3", "LMX1A", "TBR1", "GAD2"))
dev.off()
saveRDS(seu2, "~/project/IPSC_2025_Data/7W_minus/7W_minus_filt_after_df_and_soupX_auto_seu")

plus_7W_filt <- readRDS("~/project/IPSC_2025_Data/7W_plus/7W_plus_filt_after_df_seu")
counts <- plus_7W_filt@assays$RNA@layers$counts
colnames(counts) <- colnames(plus_7W_filt)
rownames(counts) <- rownames(plus_7W_filt)
write10xCounts(x = counts, path = "~/project/IPSC_2025_Data/7W_plus/filtered_feature_bc_matrix2")
toc <- Seurat::Read10X("~/project/IPSC_2025_Data/7W_plus/filtered_feature_bc_matrix2")
tod <- Seurat::Read10X("~/project/IPSC_2025_Data/7W_plus/raw_feature_bc_matrix")
common_genes <- intersect(rownames(toc), rownames(tod))
toc <- toc[common_genes, ]
tod <- tod[common_genes, ]
sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
sc <- estimateSoup(sc)
plus_7W_filt <- NormalizeData(plus_7W_filt)
plus_7W_filt <- FindVariableFeatures(plus_7W_filt)
plus_7W_filt <- ScaleData(plus_7W_filt)
plus_7W_filt <- RunPCA(plus_7W_filt)
plus_7W_filt <- RunUMAP(plus_7W_filt, dims = 1:30)
plus_7W_filt <- FindNeighbors(plus_7W_filt, dims = 1:30)
plus_7W_filt <- FindClusters(plus_7W_filt, resolution = 0.8)
plus_7W_filt$Clusters <- Idents(plus_7W_filt)
sc <- setClusters(sc, setNames(plus_7W_filt$Clusters, rownames(plus_7W_filt@meta.data)))
sc <- autoEstCont(sc)
setwd("~/project/IPSC_2025_Data/7W_plus/")
tiff("7W_plus_autoEstcon.tiff", width = 500, height = 500)
autoEstCont(sc)
dev.off()
out <- adjustCounts(sc, roundToInt = TRUE)
seu2 <- CreateSeuratObject(out)
seu2 <- NormalizeData(seu2)
seu2 <- FindVariableFeatures(seu2)
seu2 <- ScaleData(seu2)
seu2 <- RunPCA(seu2)
seu2 <- RunUMAP(seu2, dims = 1:30)
seu2 <- FindNeighbors(seu2, dims = 1:30)
seu2 <- FindClusters(seu2)
seu2@meta.data <- plus_7W_filt@meta.data
tiff("7W_plus_before_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(plus_7W_filt, features = c("GLI3", "LMX1A", "TBR1", "GAD2"))
dev.off()
tiff("7W_plus_after_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(seu2, features = c("GLI3", "LMX1A", "TBR1", "GAD2"))
dev.off()
saveRDS(seu2, "~/project/IPSC_2025_Data/7W_plus/7W_plus_filt_after_df_and_soupX_auto_seu")

minus_7W_E6_filt <- readRDS("~/project/IPSC_2025_Data/7W_minus_E6/7W_minus_E6_filt_after_df_seu")
counts <- minus_7W_E6_filt@assays$RNA@layers$counts
colnames(counts) <- colnames(minus_7W_E6_filt)
rownames(counts) <- rownames(minus_7W_E6_filt)
write10xCounts(x = counts, path = "~/project/IPSC_2025_Data/7W_minus_E6/filtered_feature_bc_matrix2")
toc <- Seurat::Read10X("~/project/IPSC_2025_Data/7W_minus_E6/filtered_feature_bc_matrix2")
tod <- Seurat::Read10X("~/project/IPSC_2025_Data/7W_minus_E6/raw_feature_bc_matrix")
common_genes <- intersect(rownames(toc), rownames(tod))
toc <- toc[common_genes, ]
tod <- tod[common_genes, ]
sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
sc <- estimateSoup(sc)
minus_7W_E6_filt <- NormalizeData(minus_7W_E6_filt)
minus_7W_E6_filt <- FindVariableFeatures(minus_7W_E6_filt)
minus_7W_E6_filt <- ScaleData(minus_7W_E6_filt)
minus_7W_E6_filt <- RunPCA(minus_7W_E6_filt)
minus_7W_E6_filt <- RunUMAP(minus_7W_E6_filt, dims = 1:30)
minus_7W_E6_filt <- FindNeighbors(minus_7W_E6_filt, dims = 1:30)
minus_7W_E6_filt <- FindClusters(minus_7W_E6_filt, resolution = 0.4)
minus_7W_E6_filt$Clusters <- Idents(minus_7W_E6_filt)
sc <- setClusters(sc, setNames(minus_7W_E6_filt$Clusters, rownames(minus_7W_E6_filt@meta.data)))
sc <- autoEstCont(sc)
setwd("~/project/IPSC_2025_Data/7W_minus_E6/")
tiff("7W_minus_E6_autoEstcon.tiff", width = 500, height = 500)
autoEstCont(sc)
dev.off()
out <- adjustCounts(sc, roundToInt = TRUE)
seu2 <- CreateSeuratObject(out)
seu2 <- NormalizeData(seu2)
seu2 <- FindVariableFeatures(seu2)
seu2 <- ScaleData(seu2)
seu2 <- RunPCA(seu2)
seu2 <- RunUMAP(seu2, dims = 1:30)
seu2 <- FindNeighbors(seu2, dims = 1:30)
seu2 <- FindClusters(seu2)
seu2@meta.data <- minus_7W_E6_filt@meta.data
tiff("7W_minus_E6_before_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(minus_7W_E6_filt, features = c("GLI3", "GAD2", "SATB2", "BCL11B"))
dev.off()
tiff("7W_minus_E6_after_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(seu2, features = c("GLI3", "GAD2", "SATB2", "BCL11B"))
dev.off()
saveRDS(seu2, "~/project/IPSC_2025_Data/7W_minus_E6/7W_minus_E6_filt_after_df_and_soupX_auto_seu")

plus_7W_E6_filt <- readRDS("~/project/IPSC_2025_Data/7W_plus_E6/7W_plus_E6_filt_after_df_seu")
counts <- plus_7W_E6_filt@assays$RNA@layers$counts
colnames(counts) <- colnames(plus_7W_E6_filt)
rownames(counts) <- rownames(plus_7W_E6_filt)
write10xCounts(x = counts, path = "~/project/IPSC_2025_Data/7W_plus_E6/filtered_feature_bc_matrix2")
toc <- Seurat::Read10X("~/project/IPSC_2025_Data/7W_plus_E6/filtered_feature_bc_matrix2")
tod <- Seurat::Read10X("~/project/IPSC_2025_Data/7W_plus_E6/raw_feature_bc_matrix")
common_genes <- intersect(rownames(toc), rownames(tod))
toc <- toc[common_genes, ]
tod <- tod[common_genes, ]
sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
sc <- estimateSoup(sc)
plus_7W_E6_filt <- NormalizeData(plus_7W_E6_filt)
plus_7W_E6_filt <- FindVariableFeatures(plus_7W_E6_filt)
plus_7W_E6_filt <- ScaleData(plus_7W_E6_filt)
plus_7W_E6_filt <- RunPCA(plus_7W_E6_filt)
plus_7W_E6_filt <- RunUMAP(plus_7W_E6_filt, dims = 1:30)
plus_7W_E6_filt <- FindNeighbors(plus_7W_E6_filt, dims = 1:30)
plus_7W_E6_filt <- FindClusters(plus_7W_E6_filt, resolution = 0.4)
plus_7W_E6_filt$Clusters <- Idents(plus_7W_E6_filt)
sc <- setClusters(sc, setNames(plus_7W_E6_filt$Clusters, rownames(plus_7W_E6_filt@meta.data)))
sc <- autoEstCont(sc)
setwd("~/project/IPSC_2025_Data/7W_plus_E6/")
tiff("7W_plus_E6_autoEstcon.tiff", width = 500, height = 500)
autoEstCont(sc)
dev.off()
out <- adjustCounts(sc, roundToInt = TRUE)
seu2 <- CreateSeuratObject(out)
seu2 <- NormalizeData(seu2)
seu2 <- FindVariableFeatures(seu2)
seu2 <- ScaleData(seu2)
seu2 <- RunPCA(seu2)
seu2 <- RunUMAP(seu2, dims = 1:30)
seu2 <- FindNeighbors(seu2, dims = 1:30)
seu2 <- FindClusters(seu2)
seu2@meta.data <- plus_7W_E6_filt@meta.data
tiff("7W_plus_E6_before_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(plus_7W_E6_filt, features = c("GLI3", "GAD2", "SATB2", "BCL11B"))
dev.off()
tiff("7W_plus_E6_after_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(seu2, features = c("GLI3", "GAD2", "SATB2", "BCL11B"))
dev.off()
saveRDS(seu2, "~/project/IPSC_2025_Data/7W_plus_E6/7W_plus_E6_filt_after_df_and_soupX_auto_seu")

minus_7W_2_filt <- readRDS("~/project/IPSC_2025_Data/7W_minus_2/7W_minus_2_filt_after_df_seu")
counts <- minus_7W_2_filt@assays$RNA@layers$counts
colnames(counts) <- colnames(minus_7W_2_filt)
rownames(counts) <- rownames(minus_7W_2_filt)
write10xCounts(x = counts, path = "~/project/IPSC_2025_Data/7W_minus_2/filtered_feature_bc_matrix2")
toc <- Seurat::Read10X("~/project/IPSC_2025_Data/7W_minus_2/filtered_feature_bc_matrix2")
tod <- Seurat::Read10X("~/project/IPSC_2025_Data/7W_minus_2/raw_feature_bc_matrix")
common_genes <- intersect(rownames(toc), rownames(tod))
toc <- toc[common_genes, ]
tod <- tod[common_genes, ]
sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
sc <- estimateSoup(sc)
minus_7W_2_filt <- NormalizeData(minus_7W_2_filt)
minus_7W_2_filt <- FindVariableFeatures(minus_7W_2_filt)
minus_7W_2_filt <- ScaleData(minus_7W_2_filt)
minus_7W_2_filt <- RunPCA(minus_7W_2_filt)
minus_7W_2_filt <- RunUMAP(minus_7W_2_filt, dims = 1:30)
minus_7W_2_filt <- FindNeighbors(minus_7W_2_filt, dims = 1:30)
minus_7W_2_filt <- FindClusters(minus_7W_2_filt, resolution = 0.4)
minus_7W_2_filt$Clusters <- Idents(minus_7W_2_filt)
sc <- setClusters(sc, setNames(minus_7W_2_filt$Clusters, rownames(minus_7W_2_filt@meta.data)))
sc <- autoEstCont(sc)
setwd("~/project/IPSC_2025_Data/7W_minus_2/")
tiff("7W_minus_2_autoEstcon.tiff", width = 500, height = 500)
autoEstCont(sc)
dev.off()
out <- adjustCounts(sc, roundToInt = TRUE)
seu2 <- CreateSeuratObject(out)
seu2 <- NormalizeData(seu2)
seu2 <- FindVariableFeatures(seu2)
seu2 <- ScaleData(seu2)
seu2 <- RunPCA(seu2)
seu2 <- RunUMAP(seu2, dims = 1:30)
seu2 <- FindNeighbors(seu2, dims = 1:30)
seu2 <- FindClusters(seu2)
seu2@meta.data <- minus_7W_2_filt@meta.data
tiff("7W_minus_2_before_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(minus_7W_2_filt, features = c("GLI3", "LMX1A", "DCX", "CFAP299"))
dev.off()
tiff("7W_minus_2_after_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(seu2, features = c("GLI3", "LMX1A", "DCX", "CFAP299"))
dev.off()
saveRDS(seu2, "~/project/IPSC_2025_Data/7W_minus_2/7W_minus_2_filt_after_df_and_soupX_auto_seu")

plus_7W_2_filt <- readRDS("~/project/IPSC_2025_Data/7W_plus_2/7W_plus_2_filt_after_df_seu")
counts <- plus_7W_2_filt@assays$RNA@layers$counts
colnames(counts) <- colnames(plus_7W_2_filt)
rownames(counts) <- rownames(plus_7W_2_filt)
write10xCounts(x = counts, path = "~/project/IPSC_2025_Data/7W_plus_2/filtered_feature_bc_matrix2")
toc <- Seurat::Read10X("~/project/IPSC_2025_Data/7W_plus_2/filtered_feature_bc_matrix2")
tod <- Seurat::Read10X("~/project/IPSC_2025_Data/7W_plus_2/raw_feature_bc_matrix")
common_genes <- intersect(rownames(toc), rownames(tod))
toc <- toc[common_genes, ]
tod <- tod[common_genes, ]
sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
sc <- estimateSoup(sc)
plus_7W_2_filt <- NormalizeData(plus_7W_2_filt)
plus_7W_2_filt <- FindVariableFeatures(plus_7W_2_filt)
plus_7W_2_filt <- ScaleData(plus_7W_2_filt)
plus_7W_2_filt <- RunPCA(plus_7W_2_filt)
plus_7W_2_filt <- RunUMAP(plus_7W_2_filt, dims = 1:30)
plus_7W_2_filt <- FindNeighbors(plus_7W_2_filt, dims = 1:30)
plus_7W_2_filt <- FindClusters(plus_7W_2_filt, resolution = 0.4)
plus_7W_2_filt$Clusters <- Idents(plus_7W_2_filt)
sc <- setClusters(sc, setNames(plus_7W_2_filt$Clusters, rownames(plus_7W_2_filt@meta.data)))
sc <- autoEstCont(sc)
setwd("~/project/IPSC_2025_Data/7W_plus_2/")
tiff("7W_plus_2_autoEstcon.tiff", width = 500, height = 500)
autoEstCont(sc)
dev.off()
out <- adjustCounts(sc, roundToInt = TRUE)
seu2 <- CreateSeuratObject(out)
seu2 <- NormalizeData(seu2)
seu2 <- FindVariableFeatures(seu2)
seu2 <- ScaleData(seu2)
seu2 <- RunPCA(seu2)
seu2 <- RunUMAP(seu2, dims = 1:30)
seu2 <- FindNeighbors(seu2, dims = 1:30)
seu2 <- FindClusters(seu2)
seu2@meta.data <- plus_7W_2_filt@meta.data
tiff("7W_plus_2_before_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(plus_7W_2_filt, features = c("GLI3", "LMX1A", "DCX", "CFAP299"))
dev.off()
tiff("7W_plus_2_after_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(seu2, features = c("GLI3", "LMX1A", "DCX", "CFAP299"))
dev.off()
saveRDS(seu2, "~/project/IPSC_2025_Data/7W_plus_2/7W_plus_2_filt_after_df_and_soupX_auto_seu")

minus_12W_filt <- readRDS("~/project/IPSC_2025_Data/12W_minus/12W_minus_filt_after_df_seu")
counts <- minus_12W_filt@assays$RNA@layers$counts
colnames(counts) <- colnames(minus_12W_filt)
rownames(counts) <- rownames(minus_12W_filt)
write10xCounts(x = counts, path = "~/project/IPSC_2025_Data/12W_minus/filtered_feature_bc_matrix2")
toc <- Seurat::Read10X("~/project/IPSC_2025_Data/12W_minus/filtered_feature_bc_matrix2")
tod <- Seurat::Read10X("~/project/IPSC_2025_Data/12W_minus/raw_feature_bc_matrix")
common_genes <- intersect(rownames(toc), rownames(tod))
toc <- toc[common_genes, ]
tod <- tod[common_genes, ]
sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
sc <- estimateSoup(sc)
minus_12W_filt <- NormalizeData(minus_12W_filt)
minus_12W_filt <- FindVariableFeatures(minus_12W_filt)
minus_12W_filt <- ScaleData(minus_12W_filt)
minus_12W_filt <- RunPCA(minus_12W_filt)
minus_12W_filt <- RunUMAP(minus_12W_filt, dims = 1:30)
minus_12W_filt <- FindNeighbors(minus_12W_filt, dims = 1:30)
minus_12W_filt <- FindClusters(minus_12W_filt, resolution = 0.8)
minus_12W_filt$Clusters <- Idents(minus_12W_filt)
sc <- setClusters(sc, setNames(minus_12W_filt$Clusters, rownames(minus_12W_filt@meta.data)))
sc <- autoEstCont(sc)
setwd("~/project/IPSC_2025_Data/12W_minus/")
tiff("12W_minus_autoEstcon.tiff", width = 500, height = 500)
autoEstCont(sc)
dev.off()
out <- adjustCounts(sc, roundToInt = TRUE)
seu2 <- CreateSeuratObject(out)
seu2 <- NormalizeData(seu2)
seu2 <- FindVariableFeatures(seu2)
seu2 <- ScaleData(seu2)
seu2 <- RunPCA(seu2)
seu2 <- RunUMAP(seu2, dims = 1:30)
seu2 <- FindNeighbors(seu2, dims = 1:30)
seu2 <- FindClusters(seu2)
seu2@meta.data <- minus_12W_filt@meta.data
tiff("12W_minus_before_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(minus_12W_filt, features = c("GLI3", "LMX1A", "TBR1", "GAD2"))
dev.off()
tiff("12W_minus_after_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(seu2, features = c("GLI3", "LMX1A", "TBR1", "GAD2"))
dev.off()
saveRDS(seu2, "~/project/IPSC_2025_Data/12W_minus/12W_minus_filt_after_df_and_soupX_auto_seu")

plus_12W_filt <- readRDS("~/project/IPSC_2025_Data/12W_plus/12W_plus_filt_after_df_seu")
counts <- plus_12W_filt@assays$RNA@layers$counts
colnames(counts) <- colnames(plus_12W_filt)
rownames(counts) <- rownames(plus_12W_filt)
write10xCounts(x = counts, path = "~/project/IPSC_2025_Data/12W_plus/filtered_feature_bc_matrix2")
toc <- Seurat::Read10X("~/project/IPSC_2025_Data/12W_plus/filtered_feature_bc_matrix2")
tod <- Seurat::Read10X("~/project/IPSC_2025_Data/12W_plus/raw_feature_bc_matrix")
common_genes <- intersect(rownames(toc), rownames(tod))
toc <- toc[common_genes, ]
tod <- tod[common_genes, ]
sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
sc <- estimateSoup(sc)
plus_12W_filt <- NormalizeData(plus_12W_filt)
plus_12W_filt <- FindVariableFeatures(plus_12W_filt)
plus_12W_filt <- ScaleData(plus_12W_filt)
plus_12W_filt <- RunPCA(plus_12W_filt)
plus_12W_filt <- RunUMAP(plus_12W_filt, dims = 1:30)
plus_12W_filt <- FindNeighbors(plus_12W_filt, dims = 1:30)
plus_12W_filt <- FindClusters(plus_12W_filt, resolution = 0.8)
plus_12W_filt$Clusters <- Idents(plus_12W_filt)
sc <- setClusters(sc, setNames(plus_12W_filt$Clusters, rownames(plus_12W_filt@meta.data)))
sc <- autoEstCont(sc)
setwd("~/project/IPSC_2025_Data/12W_plus/")
tiff("12W_plus_autoEstcon.tiff", width = 500, height = 500)
autoEstCont(sc)
dev.off()
out <- adjustCounts(sc, roundToInt = TRUE)
seu2 <- CreateSeuratObject(out)
seu2 <- NormalizeData(seu2)
seu2 <- FindVariableFeatures(seu2)
seu2 <- ScaleData(seu2)
seu2 <- RunPCA(seu2)
seu2 <- RunUMAP(seu2, dims = 1:30)
seu2 <- FindNeighbors(seu2, dims = 1:30)
seu2 <- FindClusters(seu2)
seu2@meta.data <- plus_12W_filt@meta.data
tiff("12W_plus_before_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(plus_12W_filt, features = c("GLI3", "LMX1A", "TBR1", "GAD2"))
dev.off()
tiff("12W_plus_after_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(seu2, features = c("GLI3", "LMX1A", "TBR1", "GAD2"))
dev.off()
saveRDS(seu2, "~/project/IPSC_2025_Data/12W_plus/12W_plus_filt_after_df_and_soupX_auto_seu")

minus_12W_E6_filt <- readRDS("~/project/IPSC_2025_Data/12W_minus_E6/12W_minus_E6_filt_after_df_seu")
counts <- minus_12W_E6_filt@assays$RNA@layers$counts
colnames(counts) <- colnames(minus_12W_E6_filt)
rownames(counts) <- rownames(minus_12W_E6_filt)
write10xCounts(x = counts, path = "~/project/IPSC_2025_Data/12W_minus_E6/filtered_feature_bc_matrix2")
toc <- Seurat::Read10X("~/project/IPSC_2025_Data/12W_minus_E6/filtered_feature_bc_matrix2")
tod <- Seurat::Read10X("~/project/IPSC_2025_Data/12W_minus_E6/raw_feature_bc_matrix")
common_genes <- intersect(rownames(toc), rownames(tod))
toc <- toc[common_genes, ]
tod <- tod[common_genes, ]
sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
sc <- estimateSoup(sc)
minus_12W_E6_filt <- NormalizeData(minus_12W_E6_filt)
minus_12W_E6_filt <- FindVariableFeatures(minus_12W_E6_filt)
minus_12W_E6_filt <- ScaleData(minus_12W_E6_filt)
minus_12W_E6_filt <- RunPCA(minus_12W_E6_filt)
minus_12W_E6_filt <- RunUMAP(minus_12W_E6_filt, dims = 1:30)
minus_12W_E6_filt <- FindNeighbors(minus_12W_E6_filt, dims = 1:30)
minus_12W_E6_filt <- FindClusters(minus_12W_E6_filt, resolution = 0.4)
minus_12W_E6_filt$Clusters <- Idents(minus_12W_E6_filt)
sc <- setClusters(sc, setNames(minus_12W_E6_filt$Clusters, rownames(minus_12W_E6_filt@meta.data)))
sc <- autoEstCont(sc)
setwd("~/project/IPSC_2025_Data/12W_minus_E6/")
tiff("12W_minus_E6_autoEstcon.tiff", width = 500, height = 500)
autoEstCont(sc)
dev.off()
out <- adjustCounts(sc, roundToInt = TRUE)
seu2 <- CreateSeuratObject(out)
seu2 <- NormalizeData(seu2)
seu2 <- FindVariableFeatures(seu2)
seu2 <- ScaleData(seu2)
seu2 <- RunPCA(seu2)
seu2 <- RunUMAP(seu2, dims = 1:30)
seu2 <- FindNeighbors(seu2, dims = 1:30)
seu2 <- FindClusters(seu2)
seu2@meta.data <- minus_12W_E6_filt@meta.data
tiff("12W_minus_E6_before_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(minus_12W_E6_filt, features = c("GLI3", "GAD2", "SATB2", "BCL11B"))
dev.off()
tiff("12W_minus_E6_after_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(seu2, features = c("GLI3", "GAD2", "SATB2", "BCL11B"))
dev.off()
saveRDS(seu2, "~/project/IPSC_2025_Data/12W_minus_E6/12W_minus_E6_filt_after_df_and_soupX_auto_seu")

plus_12W_E6_filt <- readRDS("~/project/IPSC_2025_Data/12W_plus_E6/12W_plus_E6_filt_after_df_seu")
counts <- plus_12W_E6_filt@assays$RNA@layers$counts
colnames(counts) <- colnames(plus_12W_E6_filt)
rownames(counts) <- rownames(plus_12W_E6_filt)
write10xCounts(x = counts, path = "~/project/IPSC_2025_Data/12W_plus_E6/filtered_feature_bc_matrix2")
toc <- Seurat::Read10X("~/project/IPSC_2025_Data/12W_plus_E6/filtered_feature_bc_matrix2")
tod <- Seurat::Read10X("~/project/IPSC_2025_Data/12W_plus_E6/raw_feature_bc_matrix")
common_genes <- intersect(rownames(toc), rownames(tod))
toc <- toc[common_genes, ]
tod <- tod[common_genes, ]
sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
sc <- estimateSoup(sc)
plus_12W_E6_filt <- NormalizeData(plus_12W_E6_filt)
plus_12W_E6_filt <- FindVariableFeatures(plus_12W_E6_filt)
plus_12W_E6_filt <- ScaleData(plus_12W_E6_filt)
plus_12W_E6_filt <- RunPCA(plus_12W_E6_filt)
plus_12W_E6_filt <- RunUMAP(plus_12W_E6_filt, dims = 1:30)
plus_12W_E6_filt <- FindNeighbors(plus_12W_E6_filt, dims = 1:30)
plus_12W_E6_filt <- FindClusters(plus_12W_E6_filt, resolution = 0.4)
plus_12W_E6_filt$Clusters <- Idents(plus_12W_E6_filt)
sc <- setClusters(sc, setNames(plus_12W_E6_filt$Clusters, rownames(plus_12W_E6_filt@meta.data)))
sc <- autoEstCont(sc)
setwd("~/project/IPSC_2025_Data/12W_plus_E6/")
tiff("12W_plus_E6_autoEstcon.tiff", width = 500, height = 500)
autoEstCont(sc)
dev.off()
out <- adjustCounts(sc, roundToInt = TRUE)
seu2 <- CreateSeuratObject(out)
seu2 <- NormalizeData(seu2)
seu2 <- FindVariableFeatures(seu2)
seu2 <- ScaleData(seu2)
seu2 <- RunPCA(seu2)
seu2 <- RunUMAP(seu2, dims = 1:30)
seu2 <- FindNeighbors(seu2, dims = 1:30)
seu2 <- FindClusters(seu2)
seu2@meta.data <- plus_12W_E6_filt@meta.data
tiff("12W_plus_E6_before_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(plus_12W_E6_filt, features = c("GLI3", "GAD2", "SATB2", "BCL11B"))
dev.off()
tiff("12W_plus_E6_after_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(seu2, features = c("GLI3", "GAD2", "SATB2", "BCL11B"))
dev.off()
saveRDS(seu2, "~/project/IPSC_2025_Data/12W_plus_E6/12W_plus_E6_filt_after_df_and_soupX_auto_seu")

minus_17W_filt <- readRDS("~/project/IPSC_2025_Data/17W_minus/17W_minus_filt_after_df_seu")
counts <- minus_17W_filt@assays$RNA@layers$counts
colnames(counts) <- colnames(minus_17W_filt)
rownames(counts) <- rownames(minus_17W_filt)
write10xCounts(x = counts, path = "~/project/IPSC_2025_Data/17W_minus/filtered_feature_bc_matrix2")
toc <- Seurat::Read10X("~/project/IPSC_2025_Data/17W_minus/filtered_feature_bc_matrix2")
tod <- Seurat::Read10X("~/project/IPSC_2025_Data/17W_minus/raw_feature_bc_matrix")
common_genes <- intersect(rownames(toc), rownames(tod))
toc <- toc[common_genes, ]
tod <- tod[common_genes, ]
sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
sc <- estimateSoup(sc)
minus_17W_filt <- NormalizeData(minus_17W_filt)
minus_17W_filt <- FindVariableFeatures(minus_17W_filt)
minus_17W_filt <- ScaleData(minus_17W_filt)
minus_17W_filt <- RunPCA(minus_17W_filt)
minus_17W_filt <- RunUMAP(minus_17W_filt, dims = 1:30)
minus_17W_filt <- FindNeighbors(minus_17W_filt, dims = 1:30)
minus_17W_filt <- FindClusters(minus_17W_filt, resolution = 0.8)
minus_17W_filt$Clusters <- Idents(minus_17W_filt)
sc <- setClusters(sc, setNames(minus_17W_filt$Clusters, rownames(minus_17W_filt@meta.data)))
sc <- autoEstCont(sc)
setwd("~/project/IPSC_2025_Data/17W_minus/")
tiff("17W_minus_autoEstcon.tiff", width = 500, height = 500)
autoEstCont(sc)
dev.off()
out <- adjustCounts(sc, roundToInt = TRUE)
seu2 <- CreateSeuratObject(out)
seu2 <- NormalizeData(seu2)
seu2 <- FindVariableFeatures(seu2)
seu2 <- ScaleData(seu2)
seu2 <- RunPCA(seu2)
seu2 <- RunUMAP(seu2, dims = 1:30)
seu2 <- FindNeighbors(seu2, dims = 1:30)
seu2 <- FindClusters(seu2)
seu2@meta.data <- minus_17W_filt@meta.data
tiff("17W_minus_before_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(minus_17W_filt, features = c("GLI3", "AQP4", "DCX", "CFAP299"))
dev.off()
tiff("17W_minus_after_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(seu2, features = c("GLI3", "AQP4", "DCX", "CFAP299"))
dev.off()
saveRDS(seu2, "~/project/IPSC_2025_Data/17W_minus/17W_minus_filt_after_df_and_soupX_auto_seu")

plus_17W_filt <- readRDS("~/project/IPSC_2025_Data/17W_plus/17W_plus_filt_after_df_seu")
counts <- plus_17W_filt@assays$RNA@layers$counts
colnames(counts) <- colnames(plus_17W_filt)
rownames(counts) <- rownames(plus_17W_filt)
write10xCounts(x = counts, path = "~/project/IPSC_2025_Data/17W_plus/filtered_feature_bc_matrix2")
toc <- Seurat::Read10X("~/project/IPSC_2025_Data/17W_plus/filtered_feature_bc_matrix2")
tod <- Seurat::Read10X("~/project/IPSC_2025_Data/17W_plus/raw_feature_bc_matrix")
common_genes <- intersect(rownames(toc), rownames(tod))
toc <- toc[common_genes, ]
tod <- tod[common_genes, ]
sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
sc <- estimateSoup(sc)
plus_17W_filt <- NormalizeData(plus_17W_filt)
plus_17W_filt <- FindVariableFeatures(plus_17W_filt)
plus_17W_filt <- ScaleData(plus_17W_filt)
plus_17W_filt <- RunPCA(plus_17W_filt)
plus_17W_filt <- RunUMAP(plus_17W_filt, dims = 1:30)
plus_17W_filt <- FindNeighbors(plus_17W_filt, dims = 1:30)
plus_17W_filt <- FindClusters(plus_17W_filt, resolution = 0.8)
plus_17W_filt$Clusters <- Idents(plus_17W_filt)
sc <- setClusters(sc, setNames(plus_17W_filt$Clusters, rownames(plus_17W_filt@meta.data)))
sc <- autoEstCont(sc)
setwd("~/project/IPSC_2025_Data/17W_plus/")
tiff("17W_plus_autoEstcon.tiff", width = 500, height = 500)
autoEstCont(sc)
dev.off()
out <- adjustCounts(sc, roundToInt = TRUE)
seu2 <- CreateSeuratObject(out)
seu2 <- NormalizeData(seu2)
seu2 <- FindVariableFeatures(seu2)
seu2 <- ScaleData(seu2)
seu2 <- RunPCA(seu2)
seu2 <- RunUMAP(seu2, dims = 1:30)
seu2 <- FindNeighbors(seu2, dims = 1:30)
seu2 <- FindClusters(seu2)
seu2@meta.data <- plus_17W_filt@meta.data
tiff("17W_plus_before_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(plus_17W_filt, features = c("GLI3", "LMX1A", "DCX", "CFAP299"))
dev.off()
tiff("17W_plus_after_soupX_markers.tiff", width = 1000, height = 1000)
FeaturePlot_scCustom(seu2, features = c("GLI3", "LMX1A", "DCX", "CFAP299"))
dev.off()
saveRDS(seu2, "~/project/IPSC_2025_Data/17W_plus/17W_plus_filt_after_df_and_soupX_auto_seu")

library(Seurat)

minus_3W_filt <- readRDS("~/project/IPSC_2025_Data/3W_minus/3W_minus_filt_after_df_and_soupX_auto_seu")
plus_3W_filt <- readRDS("~/project/IPSC_2025_Data/3W_plus/3W_plus_filt_after_df_and_soupX_auto_seu")
minus_7W_filt <- readRDS("~/project/IPSC_2025_Data/7W_minus/7W_minus_filt_after_df_and_soupX_auto_seu")
plus_7W_filt <- readRDS("~/project/IPSC_2025_Data/7W_plus/7W_plus_filt_after_df_and_soupX_auto_seu")
minus_7W_2_filt <- readRDS("~/project/IPSC_2025_Data/7W_minus_2/7W_minus_2_filt_after_df_and_soupX_auto_seu")
plus_7W_2_filt <- readRDS("~/project/IPSC_2025_Data/7W_plus_2/7W_plus_2_filt_after_df_and_soupX_auto_seu")
minus_7W_E6_filt <- readRDS("~/project/IPSC_2025_Data/7W_minus_E6/7W_minus_E6_filt_after_df_and_soupX_auto_seu")
plus_7W_E6_filt <- readRDS("~/project/IPSC_2025_Data/7W_plus_E6/7W_plus_E6_filt_after_df_and_soupX_auto_seu")
minus_12W_filt <- readRDS("~/project/IPSC_2025_Data/12W_minus/12W_minus_filt_after_df_and_soupX_auto_seu")
plus_12W_filt <- readRDS("~/project/IPSC_2025_Data/12W_plus/12W_plus_filt_after_df_and_soupX_auto_seu")
minus_12W_E6_filt <- readRDS("~/project/IPSC_2025_Data/12W_minus_E6/12W_minus_E6_filt_after_df_and_soupX_auto_seu")
plus_12W_E6_filt <- readRDS("~/project/IPSC_2025_Data/12W_plus_E6/12W_plus_E6_filt_after_df_and_soupX_auto_seu")
minus_17W_filt <- readRDS("~/project/IPSC_2025_Data/17W_minus/17W_minus_filt_after_df_and_soupX_auto_seu")
plus_17W_filt <- readRDS("~/project/IPSC_2025_Data/17W_plus/17W_plus_filt_after_df_and_soupX_auto_seu")
merged <- merge(
  minus_3W_filt,
  y = list(
    plus_3W_filt,
    minus_7W_filt,
    plus_7W_filt,
    minus_7W_2_filt,
    plus_7W_2_filt,
    minus_7W_E6_filt,
    plus_7W_E6_filt,
    minus_12W_filt,
    plus_12W_filt,
    minus_12W_E6_filt,
    plus_12W_E6_filt,
    minus_17W_filt,
    plus_17W_filt
  ),
  project = "IPSC_2025_merged"
)

merged <- NormalizeData(merged) %>% FindVariableFeatures() %>% ScaleData(merged_2, vars.to.regress = "neural_induction_media") %>% RunPCA() %>% RunHarmony(group.by.vars = c("gt_line", "Protocol", "neural_induction_media")) %>% RunUMAP(reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony") %>% FindNeighbors(reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 1)

merged <- RenameIdents(merged,
  `0`  = "Hem_RG",
  `1`  = "DL_ExN",
  `2`  = "Hem_RG",
  `3`  = "CGE_In",
  `4`  = "Unknown",
  `5`  = "MGE_In",
  `6`  = "Unknown",
  `7`  = "RG-mixed",
  `8`  = "UL_ExN",
  `9`  = "CGE_In",
  `10` = "CRN",
  `11` = "UL_ExN",
  `12` = "Hem_RG",
  `13` = "Epithelial",
  `14` = "Hem_RG",
  `15` = "RG",
  `16` = "MGE_In",
  `17` = "DL_ExN",
  `18` = "LGE_In",
  `19` = "DL_ExN",
  `20` = "MGE_In",
  `21` = "UL_ExN",
  `22` = "Unknown",
  `23` = "RG_mixed",
  `24` = "LGE_In",
  `25` = "RG-mixed",
  `26` = "Unknown",
  `27` = "Epithelial",
  `28` = "RG",
  `29` = "CRN",
  `30` = "Hem_RG",
  `31` = "Astrocyte",
  `32` = "oRG",
  `33` = "UL_ExN",
  `34` = "CRN",
  `35` = "MGE_In")

#manual section of small population of Astrocytes
iPSC_merged$ipsc_only_cluster <- Idents(iPSC_merged)
plot <- FeaturePlot_scCustom(iPSC_merged, features = "GFAP")
cells.located <- CellSelector(plot = plot)
iPSC_merged$ipsc_only_cluster[cells.located] <- "Astrocyte"

#manual section of small population of LGE_In
plot <- FeaturePlot_scCustom(iPSC_merged, features = "ISL1")
cells.located <- CellSelector(plot = plot)
iPSC_merged$ipsc_only_cluster[cells.located] <- "LGE_In"

marker_list <- list(
ExN = c("SATB2", "TAFA1", "FEZF2", "DOK5", "SLC17A7", "SLC17A6", "NEUROD1", "NEUROD4", "NEUROD2", "NEUROD6", "NEUROG1", "EOMES", "NEUROG2", "TBR1"),
In = c("GAD2", "GAD1", "SLC32A1", "DLX6-AS1", "DLX1", "DLX2", "DLX5", "DLX6"),
CRN = c("RELN", "EBF3", "LHX5", "LHX1", "TP73", "MAB21L1"),
UL = c("SATB2", "CUX2", "RORB", "POU3F2", "DOK5", "NRGN"),
DL = c("FEZF2", "BCL11B", "CRYM", "SEMA3E", "SORCS2", "HS3ST4", "ETV1"),  
CGE = c("ADARB2", "CALB2", "VIP", "CCK", "HTR3A", "NR2F1", "NR2F2"),
LGE = c("SP8", "SIX3", "ISL1", "ZNF503"),
MGE = c("LHX6", "NKX2-1", "SP9"),
IPC_EN = c("EOMES", "NEUROD4", "NEUROG2", "PPP1R17", "NEUROD1"),
Astrocyte = c("CD44", "GFAP", "S100B", "AQP4", "ALDOC"),
Hem = c("LMX1A", "WNT3A", "RSPO2"),
RG = c("PAX6", "EMX2", "SALL1", "GLI3", "HES1", "PTN", "SLC1A3", "SOX2", "HMGA2", "PCNA"),
Epithelial = c("DNAAF1", "VWA3A", "CFAP47","DEUP1", "SHISA8"),
#L6b = c("NXPH4", "CCN2", "NR4A2", "ST18"),
#Vascular = c("ESAM", "PECAM1", "PDGFRB")
)

merged <- JoinLayers(merged)

for (nm in names(marker_list)) {
  merged <- AddModuleScore(
    merged,
    features = list(marker_list[[nm]]),
    name = paste0(nm, "_new_mod")   # <-- use "new_mod" suffix
  )
  
  # Seurat will create columns like "RG_new_mod1"
  # Rename them to just "RG_new_mod"
  merged[[paste0(nm, "_new_mod")]] <- merged[[paste0(nm, "_new_mod1")]]
  merged[[paste0(nm, "_new_mod1")]] <- NULL
}

Celltype0 <- merged$ipsc_only_cluster
Celltype1 <- as.numeric(merged$ExN_new_mod)
Celltype2 <- as.numeric(merged$In_new_mod)
Celltype3 <- as.numeric(merged$IPC_EN_new_mod)
Celltype4 <- as.numeric(merged$RG_new_mod)
Celltype5 <- as.numeric(merged$DL_new_mod)
Celltype6 <- as.numeric(merged$UL_new_mod)
Celltype7 <- as.numeric(merged$CRN_new_mod)


names(Celltype0) <- colnames(merged)
names(Celltype1) <- colnames(merged)
names(Celltype2) <- colnames(merged)
names(Celltype3) <- colnames(merged)
names(Celltype4) <- colnames(merged)
names(Celltype5) <- colnames(merged)
names(Celltype6) <- colnames(merged)
names(Celltype7) <- colnames(merged)

consensusClusterLabels <- Celltype0
# All new labels you assign
new_labels <- c("IPC_In", "Unknown", "IPC_ExN")

# Expand levels on the object you will modify
consensusClusterLabels <- factor(
  consensusClusterLabels,
  levels = union(levels(consensusClusterLabels), new_labels)
)

consensusClusterLabels[names(which(Celltype0 == "RG-mixed" & Celltype2 > 0))] <- "IPC_In"
consensusClusterLabels[names(which(Celltype0 == "CRN" & Celltype7 <= 0))] <- "Unknown"
consensusClusterLabels[names(which(Celltype3 > 0))] <- "IPC_ExN"
consensusClusterLabels[names(which(Celltype0 == "UL_ExN" & Celltype6 <= 0 & Celltype2 > 0))] <- "CGE_In"
consensusClusterLabels[names(which(Celltype0 == "DL_ExN" & Celltype2 > 0))] <- "CGE_In"
consensusClusterLabels[consensusClusterLabels == "RG-mixed"] <- "RG"
consensusClusterLabels[consensusClusterLabels == "oRG"] <- "RG"

merged$ipsc_only_cluster_consensus <- consensusClusterLabels
saveRDS(merged, "~/project/IPSC_2025_Data/merged_IPSC_derived_forebrain")
