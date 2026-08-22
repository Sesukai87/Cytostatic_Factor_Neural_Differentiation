library(Seurat)
library(SeuratExtend)
library(harmony)
library(magrittr)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)

merged_fetal_ipsc <- readRDS("~/project/IPSC_2025_Data/merged_Fetal_IPSC_derived_forebrain")

merged_fetal_ipsc <- NormalizeData(merged_fetal_ipsc) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony(group.by.vars = "SampleID") %>%
  RunUMAP(reduction = "harmony", dims = 1:30)

# ---------------------------------------------------------------------------
# 1) Re-annotate fetal cells in merged_fetal_ipsc$Celltype using fetal$type
# ---------------------------------------------------------------------------
# `fetal` is the fetal-only subset object that carries the fine-grained
# `type` annotation. We copy those fine labels over into merged_fetal_ipsc
# for the matching cells, leaving the IPSC-derived cells' existing Celltype
# untouched.
#
# NOTE: merge() appended a numeric suffix (e.g. "_2") to fetal barcodes in
# merged_fetal_ipsc to disambiguate them from IPSC barcodes, so they no
# longer match colnames(fetal) exactly. We match by stripping that trailing
# "_<number>" suffix from the merged object's fetal-derived barcodes.

merged_fetal_ipsc$Celltype <- as.character(merged_fetal_ipsc$Celltype)

merged_cellnames <- colnames(merged_fetal_ipsc)
is_fetal_cell <- merged_fetal_ipsc$Sampletype == "Fetal"
stripped_barcode <- sub("_[0-9]+$", "", merged_cellnames)
fetal <- readRDS("~/project/11_29_25_merged_IPSC_Fetal")
fetal <- subset(fetal, Sampletype == "Fetal")
fetal_types <- as.character(fetal$type)
names(fetal_types) <- colnames(fetal)

# only attempt to match cells flagged as Fetal in Sampletype
match_idx <- is_fetal_cell & (stripped_barcode %in% names(fetal_types))
stopifnot(sum(is_fetal_cell) == sum(match_idx))  # sanity check: every fetal cell found a match

merged_fetal_ipsc$Celltype[match_idx] <- fetal_types[stripped_barcode[match_idx]]

# ---------------------------------------------------------------------------
# 2) Collapse the new fine fetal labels back into the original broad
#    Celltype groups, EXCEPT Microglia and Oligodendrocyte(-Immature),
#    which are kept as new, distinct categories.
# ---------------------------------------------------------------------------
# EDIT THIS LOOKUP TABLE if any of the assumed groupings below don't match
# your intended biology (a couple are flagged as best-guess assumptions).

celltype_map <- c(
  # Excitatory neurons - IT/upper-layer-like -> UL_ExN
  "EN-IT-Immature"           = "UL_ExN",
  "EN-L2_3-IT"                = "UL_ExN",
  "EN-L4-IT"                  = "UL_ExN",
  "EN-L5-IT"                  = "UL_ExN",
  "EN-L6-IT"                  = "UL_ExN",
  
  # Excitatory neurons - non-IT/deep-layer-like -> DL_ExN
  "EN-Non-IT-Immature"        = "DL_ExN",
  "EN-L5-ET"                  = "DL_ExN",
  "EN-L6b"                    = "DL_ExN",
  "EN-L5_6-NP"                = "DL_ExN",
  "EN-L6-CT"                  = "DL_ExN",
  # IPC (excitatory) -> IPC_ExN
  "IPC-EN"                    = "IPC_ExN",
  "EN-Newborn"                = "IPC_ExN",
  # Radial glia -> RG
  "RG-vRG"                    = "RG",
  "RG-oRG"                    = "RG",
  "RG-tRG"                    = "RG",
  
  # Interneurons - MGE-derived -> MGE_In
  "IN-MGE-Immature"           = "MGE_In",
  "IN-MGE-PV"                 = "MGE_In",
  "IN-MGE-SST"                = "MGE_In",
  
  # Interneurons - CGE-derived -> CGE_In
  "IN-CGE-Immature"           = "CGE_In",
  "IN-CGE-VIP"                = "CGE_In",
  "IN-CGE-SNCG"                = "CGE_In",
  "IN-CGE-LAMP5"               = "CGE_In",
  
  # Interneurons - LGE-derived -> LGE_In
  "IN-dLGE-Immature"          = "LGE_In",
  
  # Astrocytes -> Astrocyte
  "Astrocyte-Protoplasmic"    = "Astrocyte",
  "Astrocyte-Immature"        = "Astrocyte",
  "Astrocyte-Fibrous"         = "Astrocyte",
  "IPC-Glia"                  = "Astrocyte",
  
  # Cajal-Retzius -> CRN
  "Cajal-Retzius cell"        = "CRN",
  
  # Vascular -> Vascular/Fibroblast
  "Vascular"                  = "Vascular/Fibroblast",
  
  # Unchanged / already-final labels
  "OPC"                       = "OPC",
  "Unknown"                   = "Unknown",
  
  # NEW categories that should be KEPT as-is (not collapsed)
  "Microglia"                 = "Microglia",
  "Oligodendrocyte"           = "Oligodendrocyte",
  "Oligodendrocyte-Immature"  = "Oligodendrocyte"
)

# Sanity check: flag any fetal$type values not covered by the map
uncovered <- setdiff(unique(fetal_types), names(celltype_map))
if (length(uncovered) > 0) {
  warning("The following fetal$type values are not in celltype_map and will be left as-is: ",
          paste(uncovered, collapse = ", "))
}

remapped <- celltype_map[merged_fetal_ipsc$Celltype]
# keep original value where no mapping entry exists (e.g. original IPSC labels
# like IPC_In, Hem_RG, Epithelial that were never touched in step 1)
remapped[is.na(remapped)] <- merged_fetal_ipsc$Celltype[is.na(remapped)]
remapped <- unname(remapped)  # strip names (currently the lookup keys, not cell barcodes)
# so Seurat's `$<-` assigns by position, not by name-matching
merged_fetal_ipsc$Celltype <- remapped

desired_order <- c("RG", "IPC_ExN", "DL_ExN", "UL_ExN", "IPC_In", "CGE_In",
                   "MGE_In", "LGE_In", "Hem_RG", "CRN", "Epithelial",
                   "Astrocyte", "OPC", "Vascular/Fibroblast",
                   "Microglia", "Oligodendrocyte", "Unknown")
desired_order <- desired_order[desired_order %in% unique(merged_fetal_ipsc$Celltype)]
merged_fetal_ipsc$Celltype <- factor(merged_fetal_ipsc$Celltype, levels = desired_order)

# ---------------------------------------------------------------------------
# 3) Rename Protocol groups and build side-by-side DimPlots with a shared
#    Celltype color mapping across all three panels
# ---------------------------------------------------------------------------
protocol_char <- as.character(merged_fetal_ipsc$Protocol)
protocol_char[is.na(protocol_char)] <- "Fetal"
protocol_char[protocol_char == "minus"] <- "-SDF"
protocol_char[protocol_char == "plus"]  <- "+SDF"
merged_fetal_ipsc$Protocol <- factor(protocol_char, levels = c("Fetal", "-SDF", "+SDF"))
# NOTE: run this block in full (not just a patch on top of a prior partial
# run) any time Protocol needs correcting - if raw values other than
# "minus"/"plus"/NA slip through unmapped, they'll silently become NA here.
stopifnot(!anyNA(merged_fetal_ipsc$Protocol))

# Build one consistent color palette keyed by Celltype level so every panel
# uses the same color for the same cell type, even if a given group is
# missing some cell types entirely.
celltype_levels <- levels(merged_fetal_ipsc$Celltype)

# A maximally-distinct 17-color palette (default hue_pal's evenly-spaced
# hues are hard to tell apart at this many categories). Tries several
# sources in order of availability so it works without installing packages
# incompatible with your R version:
#   1. Seurat's built-in DiscretePalette("polychrome") - zero new deps,
#      Seurat is already loaded (36 maximally-distinct colors, we take 17).
#   2. pals::glasbey() - same glasbey palette Polychrome would give.
#   3. a hand-picked ColorBrewer-based fallback.
n_ct <- length(celltype_levels)
celltype_palette <- tryCatch(
  Seurat::DiscretePalette(n_ct, palette = "polychrome"),
  error = function(e) {
    if (requireNamespace("pals", quietly = TRUE)) {
      unname(pals::glasbey(n_ct))
    } else {
      c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628",
        "#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#E78AC3",
        "#A6D854","#FFD92F","#1B9E77","#D95F02","#7570B3")[seq_len(n_ct)]
    }
  }
)
celltype_colors <- setNames(celltype_palette, celltype_levels)

make_dimplot <- function(obj, protocol_label) {
  sub_obj <- subset(obj, Protocol == protocol_label)
  # Drop celltype levels with zero cells in THIS subset, so celltypes that
  # exist only in the fetal data (e.g. Microglia, OPC, Oligodendrocyte) do
  # not appear in the coloring of the iPSC-derived (-SDF / +SDF) panels
  # where they are absent. The shared celltype_colors mapping still keeps
  # each celltype's color consistent across panels.
  sub_obj$Celltype <- droplevels(sub_obj$Celltype)
  present_levels <- levels(sub_obj$Celltype)
  # order.by puts the RAREST celltypes last so they're drawn ON TOP of the
  # dense common populations - otherwise low-abundance fetal-only types
  # (Microglia, OPC, Oligodendrocyte) get buried under the big clusters and
  # look "missing" even when their points are present.
  rare_first <- names(sort(table(sub_obj$Celltype)))  # rarest ... commonest
  DimPlot(sub_obj, group.by = "Celltype", cols = celltype_colors[present_levels],
          order = rev(rare_first)) +   # DimPlot draws `order` last = on top
    ggtitle(protocol_label) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

# Build a MASTER legend from the FULL data (all celltypes present, including
# the fetal-only Microglia/OPC/Oligodendrocyte/Vascular types) so the shared
# legend lists every celltype - not just those in whichever single panel we
# happened to keep a legend on. Previously the legend was harvested from the
# +SDF panel, which legitimately lacks the fetal-only types, so they were
# missing from the legend entirely.
legend_src_plot <- DimPlot(merged_fetal_ipsc, group.by = "Celltype",
                           cols = celltype_colors) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_blank())
master_legend <- cowplot::get_legend(legend_src_plot)

# All three data panels drop their own legends; the master legend is added
# as a 4th column so it reflects every celltype.
p2_fetal <- make_dimplot(merged_fetal_ipsc, "Fetal") + NoLegend()
p2_minus <- make_dimplot(merged_fetal_ipsc, "-SDF") + NoLegend()
p2_plus  <- make_dimplot(merged_fetal_ipsc, "+SDF") + NoLegend()

# Combine the three panels + shared master legend. Using cowplot::plot_grid
# (rather than patchwork's `+`/`&`, which conflicts with Seurat's S4 `&`).
p2 <- cowplot::plot_grid(
  p2_fetal, p2_minus, p2_plus, master_legend,
  nrow = 1,
  rel_widths = c(1, 1, 1, 0.35)
)
p2

# ---------------------------------------------------------------------------
# 4) Dotplot of 3 canonical markers per Celltype, split by Protocol
#    (Seurat's DotPlot only supports a single split.by variable at a time,
#    so we build a custom dotplot: for each Celltype x Protocol x gene combo
#    we compute % expressing and average scaled expression ourselves.)
# ---------------------------------------------------------------------------

# Pick 3 canonical markers per Celltype from the provided module_list.
# Edit these choices if you'd like different markers per group.
marker_list <- list(
  RG                   = c("PAX6", "GLI3", "HES1"),
  IPC_ExN              = c("EOMES", "NEUROG2", "NEUROD4"),
  DL_ExN               = c("FEZF2", "BCL11B", "CRYM"),
  UL_ExN               = c("SATB2", "CUX2", "RORB"),
  IPC_In               = c("DLX1", "DLX2", "GAD1"),
  CGE_In               = c("ADARB2", "VIP", "NR2F2"),
  MGE_In               = c("LHX6", "NKX2-1", "SP9"),
  LGE_In               = c("SP8", "SIX3", "ISL1"),
  Hem_RG               = c("LMX1A", "WNT3A", "RSPO2"),
  CRN                  = c("RELN", "LHX5", "TP73"),
  Epithelial           = c("DNAAF1", "VWA3A", "DEUP1"),
  Astrocyte            = c("GFAP", "AQP4", "S100B"),
  OPC                  = c("PDGFRA", "OLIG1", "OLIG2"),
  "Vascular/Fibroblast" = c("PECAM1", "PDGFRB", "COL1A1"),
  Microglia            = c("AIF1", "CX3CR1", "PTPRC"),
  Oligodendrocyte      = c("MBP", "PLP1", "MOG")
  # NOTE: "Unknown" intentionally omitted here - it has no markers, but is
  # still shown as an empty labeled row in the dotplot below (see drop = FALSE)
)
marker_list <- marker_list[names(marker_list) %in% levels(merged_fetal_ipsc$Celltype)]

genes_present <- unlist(marker_list)
genes_present <- genes_present[genes_present %in% rownames(merged_fetal_ipsc)]
missing_genes <- setdiff(unlist(marker_list), genes_present)
if (length(missing_genes) > 0) {
  warning("These marker genes are missing from the object and will be dropped: ",
          paste(missing_genes, collapse = ", "))
}

# gene -> Celltype it represents, so we can later restrict each gene's row
# to only the Celltype/Protocol combos of interest if desired
gene_celltype <- unlist(lapply(names(marker_list), function(ct) {
  setNames(rep(ct, length(marker_list[[ct]])), marker_list[[ct]])
}))
gene_celltype <- gene_celltype[genes_present]

DefaultAssay(merged_fetal_ipsc) <- "RNA"
expr_mat <- GetAssayData(merged_fetal_ipsc, layer = "data")[genes_present, , drop = FALSE]

meta <- merged_fetal_ipsc@meta.data
meta$cell <- rownames(meta)

dot_df <- lapply(genes_present, function(g) {
  x <- expr_mat[g, ]
  df <- data.frame(
    cell = names(x),
    expr = as.numeric(x),
    Celltype = meta$Celltype[match(names(x), meta$cell)],
    Protocol = meta$Protocol[match(names(x), meta$cell)]
  )
  df %>%
    group_by(Celltype, Protocol) %>%
    summarise(
      pct.exp = 100 * mean(expr > 0),
      avg.exp = mean(expm1(expr)),
      .groups = "drop"
    ) %>%
    mutate(gene = g)
}) %>% bind_rows()

# Ensure every (Celltype, gene) pair we intend to show has a row for ALL
# three Protocol levels, even if that Celltype x Protocol combo has zero
# cells - those become explicit NA rows (drawn as small grey "x" marks
# below) rather than silently disappearing from the plot.
dot_df$gene <- factor(dot_df$gene, levels = genes_present)
# Include ALL desired Celltype levels (even ones with no markers, like
# "Unknown") so they still render as an empty labeled row via drop = FALSE
# in facet_grid below, rather than being silently dropped from the plot.
dot_df$Celltype <- factor(dot_df$Celltype, levels = desired_order)
dot_df$Protocol <- factor(dot_df$Protocol, levels = c("Fetal", "-SDF", "+SDF"))

dot_df <- dot_df %>%
  tidyr::complete(tidyr::nesting(Celltype, gene), Protocol)

# scale average expression per gene (z-score across Celltype x Protocol) so
# color is comparable across genes, matching Seurat's DotPlot behavior.
# NA avg.exp values (missing combos) propagate through scale() as NA.
dot_df <- dot_df %>%
  group_by(gene) %>%
  mutate(avg.exp.scaled = as.numeric(scale(avg.exp))) %>%
  ungroup()

p3 <- ggplot(dot_df, aes(x = gene, y = Protocol)) +
  geom_point(
    data = ~ dplyr::filter(.x, !is.na(pct.exp)),
    aes(size = pct.exp, color = avg.exp.scaled)
  ) +
  geom_point(
    data = ~ dplyr::filter(.x, is.na(pct.exp)),
    shape = 4, size = 2, color = "grey70", stroke = 0.6
  ) +
  facet_grid(Celltype ~ ., switch = "y", drop = FALSE) +
  scale_color_gradient(low = "lightgrey", high = "blue", name = "Avg.\nExpression", na.value = "grey70") +
  scale_size(range = c(0, 6), name = "Pct.\nExpressed") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 10),
    strip.text.y.left = element_text(angle = 0, size = 9, face = "bold"),
    strip.placement = "outside",
    panel.spacing = unit(0.1, "lines")
  ) +
  labs(x = NULL, y = NULL, caption = "x = Celltype x Protocol combination not present in data")
p3

# ---------------------------------------------------------------------------
# 5) Combine Figure 1b + 1c into a single high-resolution figure with
#    legible text (Figure 1a will be appended manually above afterwards)
# ---------------------------------------------------------------------------
fig1_combined <- plot_grid(
  p2, p3,
  ncol = 1,
  rel_heights = c(1, 1.6),
  labels = c("B", "C"),
  label_size = 16,
  label_y = c(1, 1.02)  # nudge "C" up slightly relative to panel C; tweak this value as needed
)

ggsave(
  filename = "Figure1_BC_combined.png",
  plot = fig1_combined,
  width = 20, height = 15, units = "in",
  dpi = 600, bg = "white"
)

ggsave(
  filename = "Figure1_BC_combined.pdf",
  plot = fig1_combined,
  width = 14, height = 16, units = "in",
  dpi = 600, bg = "white"
)
