library(rstatix)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(ggrepel)
library(ggtext)
library(ggsignif)
library(data.table)
library(hrbrthemes)   # theme_ipsum()
library(viridis)
library(cowplot)      # get_legend(), plot_grid() with panel alignment
library(ggh4x)        # guide_axis_nested() for the two-tier x-axis in Fig 2e

#Figure 2a
# -----------------------------
# 1. Prepare data
# -----------------------------
merged <- readRDS("~/project/IPSC_2025_Data/merged_Fetal_IPSC_derived_forebrain")
merged_IPSC <- subset(merged, Sampletype == "IPSC-Derived")
merged_IPSC$Age_line <- paste(merged_IPSC$Age, merged_IPSC$gt_line)
data <- as.data.table(merged_IPSC@meta.data)
# Count per SampleID2 / protocol / class
prop_data <- data[, .(n = .N), 
                  by = .(SampleID2, gt_line, Age, Protocol, Celltype, Age_line, neural_induction_media)]
# Compute frequencies within SampleID2
prop_data[, freq := n / sum(n), by = .(SampleID2)]

# -----------------------------
# 2. Filter out unwanted classes
# -----------------------------
remove_classes <- c("Vascular/Fibroblast", "Unknown", "OPC", "CGE_In", "MGE_In", "LGE_In", "IPC_In")
#(for supplmental figure 2)
#  remove_classes <- c("Vascular/Fibroblast", "Unknown", "OPC", "RG", "IPC_ExN", "DL_ExN", "UL_ExN", "Astrocyte", "Hem_RG", "Epithelial", "CRN") 
prop_data <- prop_data[!Celltype %in% remove_classes]
unique(prop_data$Celltype)
# -----------------------------
# Set facet order for Age_line
# -----------------------------

# -----------------------------
# Set class order
# -----------------------------
desired_order <- c(
  "RG", 
  "IPC_ExN",
  "DL_ExN",  
  "UL_ExN",
  #"IPC_In", "CGE_In", "LGE_In", "MGE_In"  #(for supplmental figure 2)
  "Astrocyte", 
  "Hem_RG", 
  "Epithelial",  
  "CRN"
)

desired_order <- c("RG", "IPC_ExN", "DL_ExN", "UL_ExN", "Astrocyte", "IPC_In", "CGE_In", "LGE_In", "MGE_In", "Hem_RG", "CRN", "Epithelial", "Unknown")


facet_order <-  c("D21 KOLF2.1", "D49 KOLF2.1",  "D84 KOLF2.1",  "D119 KOLF2.1",  "D21 O2C3",  "D49 O2C3",  "D84 O2C3",  "D119 O2C3", "D21 JHC1", "D49 JHC1", "D84 JHC1") 


prop_data[, Celltype := factor(Celltype, levels = desired_order)]
prop_data[, Age_line := factor(Age_line, levels = facet_order)]

# Protocol labels, consistent with Figure 1 (minus -> -SDF, plus -> +SDF)
prop_data[, Protocol := as.character(Protocol)]
prop_data[, Protocol := ifelse(Protocol == "minus", "-SDF",
                               ifelse(Protocol == "plus",  "+SDF", Protocol))]
prop_data[, Protocol := factor(Protocol, levels = c("-SDF", "+SDF"))]
stopifnot(!anyNA(prop_data$Protocol))

# -----------------------------
# Shared print-legible theme used for Fig 2a/2c
# -----------------------------
fig2ac_theme <- theme_ipsum() +
  theme(
    axis.text.x = element_text(size = 25, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 25, margin = margin(r = 5)),
    axis.title.x = element_text(size = 18, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 18, face = "bold", margin = margin(r = 18)),
    strip.text.x = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 22, face = "bold"),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 15)
  )

# -----------------------------
# Plot (Figure 2a)
# -----------------------------
# Y-scale: the lowest-max facet peaks around 0.05 and the highest around 0.8.
# A single fixed y-axis (no free_y) with a sqrt transform lets both the small
# and large clusters be read clearly on the same scale, instead of a linear
# axis that would flatten the low-frequency facets or a free_y axis that
# makes facets impossible to compare directly.
p1 <- ggplot(prop_data, aes(x = Celltype, y = freq, fill = Protocol)) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_fill_viridis_d(option = "E", end = 0.9) +
  scale_y_sqrt(breaks = c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8),
               labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~Age_line, nrow = 3, ncol = 4, scales = "fixed") +
  fig2ac_theme +
  xlab("Cell Class") +
  ylab("Frequency (sqrt scale)") +
  ggtitle("Class × Protocol Composition by Age")
p1

#Supplementary Figure 1
exp.mat <- read.table(file="~/project/IPSC_2025_Data/nestorawa_forcellcycle_expressionMatrix.txt",    header = TRUE, as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
merged_IPSC <- CellCycleScoring(merged_IPSC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
merged_IPSC$Phase_Celltype <- paste(merged_IPSC$Phase, merged_IPSC$Celltype)
merged_IPSC <- SetIdent(merged_IPSC, value = "dataset")
data <- as.data.table(merged_IPSC@meta.data)
data$Age_line <- forcats::fct_relevel(data$Age_line, "D21 KOLF2.1", "D49 KOLF2.1",  "D84 KOLF2.1",  "D119 KOLF2.1",  "D21 O2C3",  "D49 O2C3",  "D84 O2C3",  "D119 O2C3", "D21 JHC1", "D49 JHC1", "D84 JHC1")
# Calculate proportions by dataset_Celltype (not just dataset)
data2 <- data[, .(n = .N), keyby = .(Age_line, Celltype, Protocol, Age, Phase_Celltype, Phase)][, freq := prop.table(n), by = .(Age_line)]
data2 <- as.data.frame(data2)
data2 <- subset(data2, Celltype == "Hem_RG" | Celltype == "RG")
p_s1 <- ggplot(data2, aes(fill = Phase, y = freq, x = Phase_Celltype)) + 
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_viridis(discrete = TRUE, option = "E") +
  scale_y_sqrt(breaks = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
               labels = scales::percent_format(accuracy = 1)) +
  ggtitle("Progenitor Phase Representation Per Age") +
  facet_wrap(~Age_line, nrow = 4, ncol = 4, scales = "fixed") +
  theme_ipsum() +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    strip.text.x = element_text(size = 25)
  ) +
  xlab("") +
  ylab("Frequency (sqrt scale)") +
  RotatedAxis() +
  scale_x_discrete(expand = expansion(mult = c(0.2, 0.2)))
ggsave("~/project/IPSC_2025_Data/Supplementary_Figure_1.png",
       plot = p_s1,
       device = "png",
       width = 20, height = 20, dpi = 300)

#Figure 2c
ICC_data <- read.csv("~/project/Quantification_Table (1).csv")
ICC_data$CTIP2.PAX6.. <- NULL
ICC_data$SATB2.PAX6.. <- NULL

# Drop blank padding rows (Excel/CSV export artifact - these have an empty
# Protocol and no real data, not missing measurements).
ICC_data <- ICC_data %>% filter(Protocol %in% c("minus", "plus"))
# Clean columns
ICC_data <- ICC_data %>%
  mutate(
    Age = Age..Days.,
    Age_line = paste(Age, Line)
  ) %>%
  select(Sample, Line, Protocol, Age, Age_line,
         PAX6 = PAX6..,
         CTIP2 = CTIP2..,
         SATB2 = SATB2..)
# Long format
ICC_long <- ICC_data %>%
  pivot_longer(
    cols = c(PAX6, CTIP2, SATB2),
    names_to = "Marker",
    values_to = "Percent"
  )
facet_order <- c(
  "21 KOLF2.1", "49 KOLF2.1", "84 KOLF2.1", "119 KOLF2.1",
  "21 O2C3",    "49 O2C3",    "84 O2C3",    "119 O2C3",
  "21 JHC1", "49 JHC1", "84 JHC1"
)
ICC_long$Age_line <- factor(ICC_long$Age_line, levels = facet_order)

# Protocol labels, consistent with Figure 1 (minus -> -SDF, plus -> +SDF)
ICC_long <- ICC_long %>%
  mutate(Protocol = case_when(
    Protocol == "minus" ~ "-SDF",
    Protocol == "plus"  ~ "+SDF",
    TRUE ~ Protocol
  )) %>%
  mutate(Protocol = factor(Protocol, levels = c("-SDF", "+SDF")))
stopifnot(!anyNA(ICC_long$Protocol))

# Wilcoxon rank-sum test (standardized across all Fig 2 panels, matching
# Figure 2e), with BH multiple-testing correction for consistency.
valid_panels <- ICC_long %>%
  group_by(Age_line, Marker) %>%
  summarise(n_groups = n_distinct(Protocol), .groups = "drop") %>%
  filter(n_groups == 2)
# Wilcoxon rank-sum test (standardized across all Fig 2 panels, matching
# Figure 2e), with BH multiple-testing correction for consistency.
#
# NOTE ON INTERPRETATION: several Age_line x Marker panels have only n=3 per
# Protocol group. A two-sided Wilcoxon rank-sum test with n=3 vs n=3 has only
# C(6,3)=20 possible rank orderings, so the smallest achievable p-value is
# 1/10 = 0.10 - it is mathematically impossible to reach p<0.05 at this
# sample size regardless of effect size. Significance stars would therefore
# misleadingly read as "no difference" for every comparison. Instead we show
# the exact p-value on each bracket so the n=3 floor is visible to the
# reader, rather than collapsing everything to "ns".
valid_panels <- ICC_long %>%
  group_by(Age_line, Marker) %>%
  summarise(n_groups = n_distinct(Protocol), .groups = "drop") %>%
  filter(n_groups == 2)

stat_test <- ICC_long %>%
  inner_join(valid_panels, by = c("Age_line", "Marker")) %>%
  group_by(Age_line, Marker) %>%
  rstatix::wilcox_test(Percent ~ Protocol) %>%   # group1/group2 = the two Protocol levels
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  mutate(p.label = paste0("p = ", formatC(p.adj, format = "f", digits = 2)))
# Panels where a Protocol group is entirely absent have no valid comparison
# to bracket, so they're simply left out of stat_test (no bracket drawn)
# rather than shown with a fabricated "ns".

y_positions <- ICC_long %>%
  group_by(Age_line, Marker) %>%
  summarise(y.position = max(Percent) * 1.10, .groups = "drop")

stat_test_full <- stat_test %>%
  left_join(y_positions, by = c("Age_line", "Marker"))

p3 <- ggplot(
  ICC_long,
  aes(x = Marker, y = Percent, fill = Protocol)
) +
  geom_boxplot(
    width = 0.7,
    outlier.size = 1.5,
    color = "black",
    position = position_dodge(width = 0.8)
  ) +
  stat_pvalue_manual(
    stat_test_full,
    label = "p.label",
    x = "Marker",
    position = position_dodge(width = 0.8),  # matches geom_boxplot dodge, so bracket
    # spans exactly the two Protocol boxes
    y.position = "y.position",
    tip.length = 0.01,
    size = 4,
    bracket.size = 0.5,
    inherit.aes = FALSE
  ) +
  scale_fill_viridis_d(option = "E", end = 0.9) +
  facet_wrap(~Age_line, scales = "fixed") +
  fig2ac_theme +
  labs(
    x = "Marker",
    y = "Percent Positive Cells",
    title = "ICC: -SDF vs +SDF (Wilcoxon, exact p) Within Each Marker and Age × Line"
  )
p3

# ===========================================================================
# ALTERNATIVE STAT-TEST VERSIONS OF FIGURE 2c (t-test and Wilcoxon), each
# showing */ns significance labels for ALL panels.
#
# NOTE ON TEST CHOICE: "Mann-Whitney U" and "Wilcoxon rank-sum" are the SAME
# test (wilcox.test on two independent groups) - there is no separate,
# less-sample-size-sensitive option between them. The genuinely different
# alternative is the parametric t-test, which (unlike the rank-sum test)
# CAN reach p < 0.05 at n=3 vs n=3. Both alternatives below use BH
# correction, matching the rest of Figure 2.
# ===========================================================================

make_fig2c <- function(test = c("wilcox", "t"), label_type = c("signif", "exact")) {
  test <- match.arg(test)
  label_type <- match.arg(label_type)
  
  valid_panels <- ICC_long %>%
    group_by(Age_line, Marker) %>%
    summarise(n_groups = n_distinct(Protocol), .groups = "drop") %>%
    filter(n_groups == 2)
  
  base <- ICC_long %>%
    inner_join(valid_panels, by = c("Age_line", "Marker")) %>%
    group_by(Age_line, Marker)
  
  st <- if (test == "wilcox") {
    base %>% rstatix::wilcox_test(Percent ~ Protocol)
  } else {
    base %>% rstatix::t_test(Percent ~ Protocol)
  }
  st <- st %>%
    rstatix::adjust_pvalue(method = "BH") %>%
    rstatix::add_significance("p.adj") %>%
    mutate(p.label = paste0("p = ", formatC(p.adj, format = "f", digits = 2)))
  
  y_positions <- ICC_long %>%
    group_by(Age_line, Marker) %>%
    summarise(y.position = max(Percent) * 1.10, .groups = "drop")
  st_full <- st %>% left_join(y_positions, by = c("Age_line", "Marker"))
  
  lab <- if (label_type == "signif") "p.adj.signif" else "p.label"
  test_name <- if (test == "wilcox") "Wilcoxon rank-sum (= Mann-Whitney U)" else "t-test"
  
  ggplot(ICC_long, aes(x = Marker, y = Percent, fill = Protocol)) +
    geom_boxplot(width = 1, outlier.size = 1.5, color = "black",
                 position = position_dodge(width = 0.8)) +
    stat_pvalue_manual(st_full, label = lab, x = "Marker",
                       position = position_dodge(width = 0.8),
                       y.position = "y.position", tip.length = 0.01,
                       size = 8, bracket.size = 0.5, inherit.aes = FALSE) +
    scale_fill_viridis_d(option = "E", end = 0.9) +
    facet_wrap(~Age_line, scales = "fixed") +
    fig2ac_theme +
    labs(x = "Marker", y = "Percent Positive Cells",
         title = paste0("ICC: -SDF vs +SDF (", test_name, ", ", label_type, ") by Marker and Age × Line"))
}

# t-test version with */ns stars (the alternative you'd revert to for 2c)
p3_ttest_signif <- make_fig2c(test = "t", label_type = "signif")
ggsave("~/project/Figure2c_ttest_signif.png",
       plot = p3_ttest_signif, device = "png", width = 20, height = 20, dpi = 300)

# Wilcoxon version with */ns stars (for completeness / comparison)
p3_wilcox_signif <- make_fig2c(test = "wilcox", label_type = "signif")
ggsave("~/project/Figure2c_wilcox_signif.png",
       plot = p3_wilcox_signif, device = "png", width = 20, height = 20, dpi = 300)

# ===========================================================================
# POWER ANALYSIS FOR FIGURE 2c NON-PARAMETRIC TEST
#
# Purpose: justify to reviewers what sample size WOULD be needed to detect
# significance with a non-parametric (Wilcoxon rank-sum / Mann-Whitney U)
# test in these comparisons - supporting the decision to use a t-test for
# 2c while keeping Wilcoxon for the better-powered Figure 2e.
#
# Approach: the Wilcoxon rank-sum test has no simple closed-form power
# formula, so we (a) estimate the standardized effect size (Cohen's d) from
# the observed -SDF vs +SDF data per panel, then (b) use the asymptotic
# relative efficiency (ARE) adjustment - the Wilcoxon test needs ~1/0.955 =
# 1.047x the t-test's sample size under normality (Pitman ARE = 3/pi). We
# compute the t-test n per group via pwr, then inflate by the ARE factor to
# approximate the Wilcoxon requirement. This is the standard practical way
# to size a Wilcoxon test.
# ===========================================================================
suppressPackageStartupMessages(library(pwr))

ARE_WILCOX_NORMAL <- pi / 3   # ~1.047; Wilcoxon n = t-test n * (pi/3) under normality

effect_sizes_2c <- ICC_long %>%
  inner_join(valid_panels, by = c("Age_line", "Marker")) %>%
  group_by(Age_line, Marker) %>%
  summarise(
    mean_minus = mean(Percent[Protocol == "-SDF"], na.rm = TRUE),
    mean_plus  = mean(Percent[Protocol == "+SDF"], na.rm = TRUE),
    sd_pooled  = sqrt(mean(c(
      var(Percent[Protocol == "-SDF"], na.rm = TRUE),
      var(Percent[Protocol == "+SDF"], na.rm = TRUE)
    ), na.rm = TRUE)),
    n_minus = sum(Protocol == "-SDF"),
    n_plus  = sum(Protocol == "+SDF"),
    .groups = "drop"
  ) %>%
  mutate(
    cohens_d = abs(mean_minus - mean_plus) / sd_pooled,
    # t-test n per group for 80% power at alpha=0.05 (two-sided)
    n_per_group_ttest = mapply(function(d) {
      if (!is.finite(d) || d == 0) return(NA_real_)
      tryCatch(ceiling(pwr.t.test(d = d, sig.level = 0.05, power = 0.80,
                                  type = "two.sample")$n),
               error = function(e) NA_real_)
    }, cohens_d),
    # approximate Wilcoxon n per group via ARE inflation
    n_per_group_wilcox = ceiling(n_per_group_ttest * ARE_WILCOX_NORMAL)
  )

# Save the power table for the reviewer response
write.csv(effect_sizes_2c,
          "~/project/IPSC_2025_Data/Figure2c_power_analysis.csv",
          row.names = FALSE)

# Visual summary: required Wilcoxon n per group vs observed effect size
p_power_2c <- ggplot(effect_sizes_2c %>% filter(is.finite(cohens_d), !is.na(n_per_group_wilcox)),
                     aes(x = cohens_d, y = n_per_group_wilcox)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_hline(yintercept = 3, linetype = "dashed", color = "firebrick") +
  annotate("text", x = Inf, y = 3, label = "current n = 3", hjust = 1.1, vjust = -0.5,
           color = "firebrick", size = 4) +
  scale_y_log10() +
  labs(
    title = "Figure 2c power analysis: Wilcoxon/Mann-Whitney n per group needed for 80% power",
    subtitle = "Each point = one Marker × Age×Line comparison; y on log scale. Points far above the n=3 line are underpowered for a rank-sum test.",
    x = "Observed standardized effect size (Cohen's d, -SDF vs +SDF)",
    y = "Required n per group (Wilcoxon, 80% power, α=0.05)"
  ) +
  theme_bw() + fig2ac_theme
ggsave("~/project/IPSC_2025_Data/Figure2c_power_analysis.png",
       plot = p_power_2c, device = "png", width = 14, height = 10, dpi = 300)

# -----------------------------
# Combine Figure 2a + 2c horizontally with ONE shared legend in the middle
# -----------------------------
shared_legend <- cowplot::get_legend(
  p1 + theme(legend.position = "right", legend.box.margin = margin(0, 0, 0, 0))
)

p1_noleg <- p1 + theme(legend.position = "none", plot.title = element_blank())
p3_noleg <- p3 + theme(legend.position = "none", plot.title = element_blank())

fig2_ac_combined <- cowplot::plot_grid(
  p1_noleg, shared_legend, p3_noleg,
  nrow = 1,
  rel_widths = c(1, 0.18, 1),
  labels = c("A", "", "C"),
  label_size = 20
)

ggsave("~/project/IPSC_2025_Data/Figure2_AC_combined.png",
       plot = fig2_ac_combined,
       device = "png",
       width = 26, height = 12, dpi = 300)


#Figure2e

############################################
### UNIVERSAL THEME (NO CLIPPING, BIG TEXT, HORIZONTAL MARKER STRIPS)
############################################

big_theme <- theme_bw() +
  theme(
    # bottom tier of the nested x-axis (Protocol: -SDF / +SDF)
    axis.text.x = element_text(size = 16, angle = 40, hjust = 1, vjust = 1),
    # upper tier of the nested x-axis (Cell Line), added by ggh4x
    ggh4x.axis.nesttext.x = element_text(size = 18, face = "bold", angle = 0),
    ggh4x.axis.nestline.x = element_line(linewidth = 0.6, color = "black"),
    axis.text.y = element_text(size = 28, margin = margin(r = 12)),
    axis.title.y = element_text(size = 30, margin = margin(r = 25)),
    strip.text.x = element_text(size = 30, margin = margin(t = 10, b = 10)),
    # Marker strips shown horizontally (not rotated) on the left of each row
    strip.text.y.left = element_text(size = 26, angle = 0, hjust = 1,
                                     margin = margin(r = 10, l = 10)),
    strip.placement = "outside",
    legend.position = "top",
    legend.text = element_text(size = 26),
    panel.spacing.x = unit(1.2, "lines"),
    plot.margin = margin(t = 30, r = 40, b = 30, l = 80),
    strip.background = element_rect(fill = "grey90", color = "black")
  )

############################################
### FUNCTION TO SAFELY COMPUTE Y-POSITIONS
############################################
# Cap well below the axis ceiling (110) so significance brackets/labels have
# clear headroom above the tallest box/point and never get clipped.

safe_ypos <- function(df) {
  yp <- max(df$Pct, na.rm = TRUE)
  if (!is.finite(yp)) yp <- 50
  pmin(yp * 1.15, 92)
}

############################################
### SHARED PANEL BUILDER
############################################
# Renames Protocol to match Figure 1 (-SDF/+SDF), builds a Cell Line / Protocol
# combined variable for a two-tier x-axis (bottom tier = Protocol, upper tier
# = Cell Line, via ggh4x::guide_axis_nested splitting on " / "), runs a
# Wilcoxon rank-sum test (BH-adjusted, consistent with Fig 2a/2c), and returns
# a plot with FIXED facet scales / y-limits so every panel in the combined
# figure shares identical scaling (no per-panel distortion).

cell_line_levels <- c("JHC1", "KOLF2.1J", "O2C3")

build_fig2e_panel <- function(df_long, test = c("wilcox", "t")) {
  test <- match.arg(test)
  
  df_long <- df_long %>%
    mutate(
      Protocol = case_when(
        Protocol == "minus" ~ "-SDF",
        Protocol == "plus"  ~ "+SDF",
        TRUE ~ as.character(Protocol)
      ),
      Protocol = factor(Protocol, levels = c("-SDF", "+SDF")),
      Cell_Protocol = paste(Cell.Line, Protocol, sep = " / "),
      Cell_Protocol = factor(
        Cell_Protocol,
        levels = paste(rep(cell_line_levels, each = 2), c("-SDF", "+SDF"), sep = " / ")
      )
    )
  stopifnot(!anyNA(df_long$Protocol))
  
  stats_base <- df_long %>%
    filter(!is.na(Pct)) %>%
    group_by(Cell.Line, Marker, Age)
  stats <- if (test == "wilcox") {
    stats_base %>% wilcox_test(Pct ~ Protocol)
  } else {
    stats_base %>% t_test(Pct ~ Protocol)
  }
  stats <- stats %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    left_join(
      df_long %>%
        group_by(Cell.Line, Marker, Age) %>%
        summarise(y.position = safe_ypos(cur_data()), .groups = "drop"),
      by = c("Cell.Line", "Marker", "Age")
    ) %>%
    mutate(
      xmin = paste(Cell.Line, "-SDF", sep = " / "),
      xmax = paste(Cell.Line, "+SDF", sep = " / ")
    )
  
  ggplot(df_long, aes(x = Cell_Protocol, y = Pct)) +
    geom_boxplot(aes(fill = Protocol), width = 0.7, outlier.shape = NA) +
    geom_point(aes(color = Protocol),
               position = position_jitter(width = 0.08),
               alpha = 0.6, size = 4) +
    facet_grid(Marker ~ Age, scales = "fixed", switch = "y") +
    stat_pvalue_manual(
      stats,
      label = "p.adj.signif",
      xmin = "xmin", xmax = "xmax",
      y.position = "y.position",
      size = 8,
      bracket.size = 0.7,
      tip.length = 0.01
    ) +
    scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25)) +
    scale_x_discrete(guide = ggh4x::guide_axis_nested(delim = " / ")) +
    big_theme +
    coord_cartesian(clip = "off") +
    labs(x = "", y = "Percent (%)")
}

############################################
### 7W: SATB2/CTIP2/SOX9/KI67 Percentages
############################################

df1 <- read.csv("~/project/7W_cell_counts.csv")

df_long <- df1 %>%
  pivot_longer(
    cols = c(SATB2.Pct, KI67.Pct, SOX9.Pct, CTIP2.Pct),
    names_to = "Marker",
    values_to = "Pct"
  ) %>%
  mutate(
    Marker = recode(Marker,
                    "SATB2.Pct" = "SATB2+",
                    "CTIP2.Pct" = "CTIP2+",
                    "SOX9.Pct"  = "SOX9+",
                    "KI67.Pct"  = "KI67+"
    ),
    Marker = factor(Marker, levels = c("SATB2+","CTIP2+","SOX9+","KI67+"))
  )

df_long_p1 <- df_long
p5_1 <- build_fig2e_panel(df_long)

############################################
### 7W: KI67 Double‑Positive Percentages
############################################

df_long <- df1 %>%
  mutate(across(
    c(`CTIP2.KI67..CTIP2.`,`SATB2.KI67..SATB2.`,`SOX9.KI67..SOX9.`),
    ~ suppressWarnings(as.numeric(trimws(.)))
  )) %>%
  pivot_longer(
    cols = c(`CTIP2.KI67..CTIP2.`,`SATB2.KI67..SATB2.`,`SOX9.KI67..SOX9.`),
    names_to = "Marker",
    values_to = "Pct"
  ) %>%
  mutate(
    Marker = recode(
      Marker,
      "CTIP2.KI67..CTIP2." = "CTIP2+KI67+/CTIP2+",
      "SATB2.KI67..SATB2." = "SATB2+KI67+/SATB2+",
      "SOX9.KI67..SOX9."   = "SOX9+KI67+/SOX9+"
    ),
    Marker = factor(Marker, levels = c("SATB2+KI67+/SATB2+","CTIP2+KI67+/CTIP2+","SOX9+KI67+/SOX9+"))
  )

df_long_p2 <- df_long
p5_2 <- build_fig2e_panel(df_long)

############################################
### 12W: SATB2/CTIP2/SOX9/KI67 Percentages
############################################

df2 <- read.csv("~/project/12W_cell_counts.csv")

df_long <- df2 %>%
  pivot_longer(
    cols = c(SATB2.Pct, KI67.Pct, SOX9.Pct, CTIP2.Pct),
    names_to = "Marker",
    values_to = "Pct"
  ) %>%
  mutate(
    Marker = recode(Marker,
                    "SATB2.Pct" = "SATB2+",
                    "CTIP2.Pct" = "CTIP2+",
                    "SOX9.Pct"  = "SOX9+",
                    "KI67.Pct"  = "KI67+"
    ),
    Marker = factor(Marker, levels = c("SATB2+","CTIP2+","SOX9+","KI67+"))
  )

df_long_p3 <- df_long
p5_3 <- build_fig2e_panel(df_long)

############################################
### 12W: KI67 Double‑Positive Percentages
############################################

df_long <- df2 %>%
  mutate(across(
    c(`CTIP2.KI67..CTIP2.`,`SATB2.KI67..SATB2.`,`SOX9.KI67..SOX9.`),
    ~ suppressWarnings(as.numeric(trimws(.)))
  )) %>%
  pivot_longer(
    cols = c(`CTIP2.KI67..CTIP2.`,`SATB2.KI67..SATB2.`,`SOX9.KI67..SOX9.`),
    names_to = "Marker",
    values_to = "Pct"
  ) %>%
  mutate(
    Marker = recode(
      Marker,
      "CTIP2.KI67..CTIP2." = "CTIP2+KI67+/CTIP2+",
      "SATB2.KI67..SATB2." = "SATB2+KI67+/SATB2+",
      "SOX9.KI67..SOX9."   = "SOX9+KI67+/SOX9+"
    ),
    Marker = factor(Marker, levels = c("SATB2+KI67+/SATB2+","CTIP2+KI67+/CTIP2+","SOX9+KI67+/SOX9+"))
  )

df_long_p4 <- df_long
p5_4 <- build_fig2e_panel(df_long)

############################################
### 2×2 Combined Figure
############################################
# Using cowplot::plot_grid with align = "hv" / axis = "tblr" instead of
# ggarrange: this aligns each panel's actual plotting area (not just its
# outer bounding box), which is what was causing the stretching/distortion
# when sub-panels had different facet-grid dimensions and axis label widths.
# All four panels also now share identical facet scales/limits (set inside
# build_fig2e_panel), so proportions are visually comparable across panels.

p5 <- cowplot::plot_grid(
  p5_1, p5_3,
  p5_2, p5_4,
  ncol = 2, nrow = 2,
  align = "hv",
  axis = "tblr",
  labels = c("7W", "12W", "", ""),
  label_size = 30
)

ggsave("~/project/IPSC_2025_Data/Figure2e_wilcox.png",
       plot = p5,
       device = "png",
       width = 28, height = 26, dpi = 300)

# ---- Alternative Figure 2e using t-tests (parametric) instead of Wilcoxon ----
# Provided so you can compare; 2e's larger n per group means Wilcoxon is
# already reasonably powered here (unlike 2c), so Wilcoxon remains the
# recommended default for 2e. Both show */ns significance labels.
p5_1_t <- build_fig2e_panel(df_long_p1, test = "t")
p5_2_t <- build_fig2e_panel(df_long_p2, test = "t")
p5_3_t <- build_fig2e_panel(df_long_p3, test = "t")
p5_4_t <- build_fig2e_panel(df_long_p4, test = "t")

p5_t <- cowplot::plot_grid(
  p5_1_t, p5_3_t,
  p5_2_t, p5_4_t,
  ncol = 2, nrow = 2,
  align = "hv",
  axis = "tblr",
  labels = c("7W", "12W", "", ""),
  label_size = 30
)
ggsave("~/project/IPSC_2025_Data/Figure2e_ttest.png",
       plot = p5_t,
       device = "png",
       width = 28, height = 26, dpi = 300)






#Supplementary Figure 4

library(tidyverse)
library(rstatix)
library(ggpubr)

############################################
### PROTOCOL COLORS
############################################

protocol_colors <- c(
  "no3i_noTF" = "#999999",   # grey  — baseline/reference
  "no3i_TF"   = "#56B4E9",   # blue  — TF only
  "3i_noTF"   = "#E69F00",   # amber — 3i only
  "3i_TF"     = "#CC79A7"    # pink  — 3i + TF
)

############################################
### PRISM-STYLE THEME
############################################

prism_theme <- theme(
  # White background, no gridlines
  panel.background = element_rect(fill = "white", color = NA),
  plot.background  = element_rect(fill = "white", color = NA),
  panel.grid       = element_blank(),
  # Axis lines: left and bottom only, bold
  axis.line.x      = element_line(color = "black", linewidth = 1),
  axis.line.y      = element_line(color = "black", linewidth = 1),
  axis.ticks       = element_line(color = "black", linewidth = 0.8),
  axis.ticks.length = unit(0.2, "cm"),
  # Text sizing
  axis.text.x  = element_text(angle = 45, hjust = 1, size = 24, color = "black"),
  axis.text.y  = element_text(size = 24, color = "black"),
  axis.title.y = element_text(size = 28, face = "bold", margin = margin(r = 15)),
  # Facet strips
  strip.background = element_blank(),
  strip.text.x = element_text(size = 24, face = "bold", margin = margin(b = 8)),
  # Legend
  legend.position  = "top",
  legend.title     = element_blank(),
  legend.text      = element_text(size = 22),
  legend.key       = element_rect(fill = "white", color = NA),
  legend.key.size  = unit(1.2, "cm"),
  # Title
  plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
  plot.margin = margin(t = 20, r = 40, b = 20, l = 40)
)

############################################
### LOAD DATA
############################################

df <- read.csv("~/project/07_25_26_different_TF_3i_condition_comparisonss.csv",
               check.names = FALSE)

# Rename columns explicitly (make.names can't distinguish +/- in these headers)
colnames(df) <- c("Image", "Cell.Line", "Protocol",
                  "KI67_Pct", "TUBB3_Pct",
                  "KI67pos_TUBB3pos_Pct",
                  "KI67pos_TUBB3neg_Pct",
                  "KI67neg_TUBB3pos_Pct")

# Pivot to long format
df_long <- df %>%
  pivot_longer(
    cols = c(KI67_Pct, TUBB3_Pct, KI67pos_TUBB3pos_Pct,
             KI67pos_TUBB3neg_Pct, KI67neg_TUBB3pos_Pct),
    names_to = "Marker",
    values_to = "Pct"
  ) %>%
  mutate(
    Marker = recode(Marker,
                    "KI67_Pct"              = "KI67+",
                    "TUBB3_Pct"             = "TUBB3+",
                    "KI67pos_TUBB3pos_Pct"  = "KI67+TUBB3+",
                    "KI67pos_TUBB3neg_Pct"  = "KI67+TUBB3-",
                    "KI67neg_TUBB3pos_Pct"  = "KI67-TUBB3+"
    ),
    Marker = factor(Marker, levels = c("KI67+", "TUBB3+",
                                       "KI67+TUBB3+", "KI67+TUBB3-", "KI67-TUBB3+")),
    Protocol = factor(Protocol, levels = c("no3i_noTF", "no3i_TF", "3i_noTF", "3i_TF"))
  )

############################################
### FUNCTION: Plot one cell line
############################################

plot_cell_line <- function(data, cell_line_name, test = c("wilcox", "t")) {
  test <- match.arg(test)
  
  dl <- data %>% filter(Cell.Line == cell_line_name)
  
  # Pairwise test: each protocol vs no3i_noTF (reference).
  # NOTE: "Mann-Whitney U" and "Wilcoxon rank-sum" are the same test; the
  # genuinely different alternative is the parametric t-test (test = "t").
  stats_base <- dl %>% group_by(Marker)
  stats <- if (test == "wilcox") {
    stats_base %>% wilcox_test(Pct ~ Protocol, ref.group = "no3i_noTF")
  } else {
    stats_base %>% t_test(Pct ~ Protocol, ref.group = "no3i_noTF")
  }
  stats <- stats %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    # Compute y positions above the tallest bar+SEM or data point
    left_join(
      dl %>%
        group_by(Marker) %>%
        summarise(
          bar_top = max({
            tapply(Pct, Protocol, function(x)
              mean(x, na.rm = TRUE) + sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
          }, na.rm = TRUE),
          pt_top = max(Pct, na.rm = TRUE),
          ymax   = max(bar_top, pt_top),
          .groups = "drop"
        ) %>% select(Marker, ymax),
      by = "Marker"
    ) %>%
    group_by(Marker) %>%
    mutate(
      rank = row_number(),
      y.position = ymax * (1.08 + 0.10 * rank)
    ) %>%
    ungroup()
  
  # Compute mean ± SEM for bar heights and error bars
  dl_summary <- dl %>%
    group_by(Protocol, Marker) %>%
    summarise(
      mean_pct = mean(Pct, na.rm = TRUE),
      sem      = sd(Pct, na.rm = TRUE) / sqrt(sum(!is.na(Pct))),
      .groups  = "drop"
    )
  
  ggplot(dl_summary, aes(x = Protocol, y = mean_pct, fill = Protocol)) +
    # Mean bars (no outline for clean Prism look)
    geom_col(width = 0.7, color = "black", linewidth = 0.5) +
    # SEM error bars
    geom_errorbar(aes(ymin = mean_pct, ymax = mean_pct + sem),
                  width = 0.25, linewidth = 0.8) +
    # Individual data points on top
    geom_point(data = dl, aes(x = Protocol, y = Pct, fill = Protocol),
               shape = 21, color = "black", size = 3.5, stroke = 0.6,
               position = position_jitter(width = 0.12, seed = 42),
               alpha = 0.7) +
    facet_wrap(~ Marker, nrow = 1, scales = "free_y") +
    stat_pvalue_manual(
      stats,
      label = "p.adj.signif",
      xmin = "group1", xmax = "group2", inherit.aes = FALSE,
      y.position = "y.position",
      size = 10,
      bracket.size = 1.2,
      tip.length = 0.01
    ) +
    scale_fill_manual(values = protocol_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    prism_theme +
    coord_cartesian(clip = "off") +
    labs(x = "", y = "Percent (%)", title = cell_line_name)
}

############################################
### GENERATE PLOTS
############################################

p_JHC1   <- plot_cell_line(df_long, "JHC1")
p_O2C3   <- plot_cell_line(df_long, "O2C3")
p_KOLF21 <- plot_cell_line(df_long, "KOLF2.1")


p_combined <- ggarrange(
  p_JHC1, p_O2C3, p_KOLF21,
  ncol = 1, nrow = 3,
  common.legend = TRUE, legend = "top"
)

ggsave("~/project/Supplemental_Figure4b_wilcox.png",
       plot = p_combined, device = "png",
       width = 24, height = 28, dpi = 300)

# ---- Alternative Supplemental Figure 4b using t-tests (parametric) ----
# Same pairwise-vs-reference design, BH-adjusted, */ns labels. Provided so
# reviewers can compare parametric vs rank-based results for these
# 4-protocol comparisons.
p_JHC1_t   <- plot_cell_line(df_long, "JHC1",    test = "t")
p_O2C3_t   <- plot_cell_line(df_long, "O2C3",    test = "t")
p_KOLF21_t <- plot_cell_line(df_long, "KOLF2.1", test = "t")

p_combined_t <- ggarrange(
  p_JHC1_t, p_O2C3_t, p_KOLF21_t,
  ncol = 1, nrow = 3,
  common.legend = TRUE, legend = "top"
)

ggsave("~/project/Supplemental_Figure4b_ttest.png",
       plot = p_combined_t, device = "png",
       width = 24, height = 28, dpi = 300)
