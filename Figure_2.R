library(rstatix)
library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(ggrepel)
library(ggtext)
library(ggsignif)

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

# -----------------------------
# Plot
# -----------------------------
p1 <- ggplot(prop_data, aes(x = Celltype, y = freq, fill = Protocol)) +
  geom_col(position = position_dodge(width = 0.8)) +
  scale_fill_viridis_d(option = "E", end = 0.9) +
  facet_wrap(~Age_line, nrow = 3, ncol = 4, scales = "free_y") +
  theme_ipsum() +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    strip.text.x = element_text(size = 22),
    legend.position = "right"
  ) +
  xlab("Cell Class") +
  ylab("Frequency") +
  ggtitle("Class × Protocol Composition by Age")

ggsave("~/project/IPSC_2025_Data/Figure2a.tiff",
       plot = p1,
       device = "tiff",
       width = 20, height = 20, dpi = 300)

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
  ggtitle("Progenitor Phase Representation Per Age") +
  facet_wrap(~Age_line, nrow = 4, ncol = 4) +
  theme_ipsum() +
  theme(
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    strip.text.x = element_text(size = 25)
  ) +
  xlab("") +
  RotatedAxis() +
  scale_x_discrete(expand = expansion(mult = c(0.2, 0.2)))
ggsave("~/project/IPSC_2025_Data/Supplementary_Figure_1.tiff",
       plot = p_s1,
       device = "tiff",
       width = 20, height = 20, dpi = 300)

#Figure 2c
ICC_data <- read.csv("~/project/IPSC_2025_Data/Quantification_Table (1).csv")
ICC_data$CTIP2.PAX6.. <- NULL
ICC_data$SATB2.PAX6.. <- NULL
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


#Wilcox
valid_panels <- ICC_long %>%
  group_by(Age_line, Marker) %>%
  summarise(n_groups = n_distinct(Protocol), .groups = "drop") %>%
  filter(n_groups == 2)
stat_test <- ICC_long %>%
  inner_join(valid_panels, by = c("Age_line", "Marker")) %>%
  group_by(Age_line, Marker) %>%
  rstatix::wilcox_test(Percent ~ Protocol) %>%
  rstatix::add_significance("p")
stat_test <- stat_test %>%
  mutate(
    group1 = Marker,
    group2 = Marker
  )
y_positions <- ICC_long %>%
  group_by(Age_line, Marker) %>%
  summarise(y.position = max(Percent) * 1.10, .groups = "drop")
stat_test <- stat_test %>%
  left_join(y_positions, by = c("Age_line", "Marker"))
p_icc_stats <- ggplot(
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
    stat_test,
    label = "p.signif",
    y.position = "y.position",
    tip.length = 0.01,
    inherit.aes = FALSE
  ) +
  scale_fill_viridis_d(option = "E", end = 0.9) +
  facet_wrap(~Age_line, scales = "free_y") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    strip.text.x = element_text(size = 18),
    legend.position = "right"
  ) +
  labs(
    x = "Marker",
    y = "Percent Positive Cells",
    title = "ICC: Plus vs Minus Comparison Within Each Marker and Age × Line"
  )
p_icc_stats


#t-test (used in main figure)
valid_panels <- ICC_long %>%
  group_by(Age_line, Marker) %>%
  summarise(n_groups = n_distinct(Protocol), .groups = "drop") %>%
  filter(n_groups == 2)
stat_test <- ICC_long %>%
  inner_join(valid_panels, by = c("Age_line", "Marker")) %>%
  group_by(Age_line, Marker) %>%
  t_test(Percent ~ Protocol) %>%
  add_significance("p")
missing_panels <- ICC_long %>%
  group_by(Age_line, Marker) %>%
  summarise(n_groups = n_distinct(Protocol), .groups = "drop") %>%
  filter(n_groups < 2) %>%
  mutate(
    p = NA,
    p.signif = "ns"
  )
stat_test_full <- bind_rows(stat_test, missing_panels) %>%
  mutate(
    group1 = Marker,
    group2 = Marker
  )
y_positions <- ICC_long %>%
  group_by(Age_line, Marker) %>%
  summarise(y.position = max(Percent) * 1.10, .groups = "drop")
stat_test_full <- stat_test_full %>%
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
    label = "p.signif",
    y.position = "y.position",
    tip.length = 0.01,
    inherit.aes = FALSE
  ) +
  scale_fill_viridis_d(option = "E", end = 0.9) +
  facet_wrap(~Age_line, scales = "free_y") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    strip.text.x = element_text(size = 18),
    legend.position = "right"
  ) +
  labs(
    x = "Marker",
    y = "Percent Positive Cells",
    title = "ICC: Plus vs Minus (t‑test) Within Each Marker and Age × Line"
  )
p3
ggsave("~/project/IPSC_2025_Data/Figure2c.tiff",
       plot = p3,
       device = "tiff",
       width = 20, height = 20, dpi = 300)


#Figure2e

############################################
### UNIVERSAL THEME (NO CLIPPING, BIG TEXT)
############################################

big_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 28),
    axis.text.y = element_text(size = 28, margin = margin(r = 12)),
    axis.title.y = element_text(size = 30, margin = margin(r = 25)),
    strip.text.x = element_text(size = 30, margin = margin(t = 10, b = 10)),
    strip.text.y = element_text(size = 30, margin = margin(r = 10, l = 10)),
    legend.position = "top",
    legend.text = element_text(size = 26),
    plot.margin = margin(t = 20, r = 40, b = 20, l = 80),
    strip.background = element_rect(fill = "grey90", color = "black")
  )

############################################
### FUNCTION TO SAFELY COMPUTE Y-POSITIONS
############################################

safe_ypos <- function(df) {
  yp <- max(df$Pct, na.rm = TRUE)
  if (!is.finite(yp)) yp <- 50
  pmin(yp * 1.12, 98)
}

############################################
### 7W: SATB2/CTIP2/SOX9/KI67 Percentages
############################################

df1 <- read.csv("~/project/IPSC_2025_Data/7W_cell_counts.csv")

df_long <- df1 %>%
  pivot_longer(
    cols = c(SATB2.Pct, KI67.Pct, SOX9.Pct, CTIP2.Pct),
    names_to = "Marker",
    values_to = "Pct"
  ) %>%
  mutate(
    Marker = recode(Marker,
                    "SATB2.Pct" = "SATB2+Pct",
                    "CTIP2.Pct" = "CTIP2+Pct",
                    "SOX9.Pct"  = "SOX9+Pct",
                    "KI67.Pct"  = "KI67+Pct"
    ),
    Cell_Protocol = paste(Cell.Line, Protocol),
    Cell_Protocol = factor(Cell_Protocol,
                           levels = c("JHC1 minus","JHC1 plus",
                                      "KOLF2.1J minus","KOLF2.1J plus",
                                      "O2C3 minus","O2C3 plus")
    ),
    Protocol = factor(Protocol, levels = c("minus","plus")),
    Marker = factor(Marker, levels = c("SATB2+Pct","CTIP2+Pct","SOX9+Pct","KI67+Pct"))
  )

stats <- df_long %>%
  group_by(Cell.Line, Marker, Age) %>%
  wilcox_test(Pct ~ Protocol) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  left_join(
    df_long %>%
      group_by(Cell.Line, Marker, Age) %>%
      summarise(y.position = safe_ypos(cur_data()), .groups = "drop"),
    by = c("Cell.Line","Marker","Age")
  ) %>%
  mutate(
    xmin = paste(Cell.Line,"minus"),
    xmax = paste(Cell.Line,"plus")
  )

p5_1 <- ggplot(df_long, aes(x = Cell_Protocol, y = Pct)) +
  geom_boxplot(aes(fill = Protocol), width = 0.7, outlier.shape = NA) +
  geom_point(aes(color = Protocol),
             position = position_jitter(width = 0.08),
             alpha = 0.6, size = 4) +
  facet_grid(Marker ~ Age, scales = "fixed") +
  stat_pvalue_manual(
    stats,
    label = "p.adj.signif",
    xmin = "xmin", xmax = "xmax",
    y.position = "y.position",
    size = 14,
    bracket.size = 1.4,
    tip.length = 0.01
  ) +
  scale_y_continuous(limits = c(0, 100)) +
  big_theme +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "Percent (%)")

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
      "CTIP2.KI67..CTIP2." = "CTIP2+KI67+/CTIP2+Pct",
      "SATB2.KI67..SATB2." = "SATB2+KI67+/SATB2+Pct",
      "SOX9.KI67..SOX9."   = "SOX9+KI67+/SOX9+Pct"
    ),
    Cell_Protocol = paste(Cell.Line, Protocol),
    Cell_Protocol = factor(Cell_Protocol,
                           levels = c("JHC1 minus","JHC1 plus",
                                      "KOLF2.1J minus","KOLF2.1J plus",
                                      "O2C3 minus","O2C3 plus")
    ),
    Protocol = factor(Protocol, levels = c("minus","plus"))
  )

stats <- df_long %>%
  filter(!is.na(Pct)) %>%
  group_by(Cell.Line, Marker, Age) %>%
  wilcox_test(Pct ~ Protocol) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  left_join(
    df_long %>%
      group_by(Cell.Line, Marker, Age) %>%
      summarise(y.position = safe_ypos(cur_data()), .groups = "drop"),
    by = c("Cell.Line","Marker","Age")
  ) %>%
  mutate(
    xmin = paste(Cell.Line,"minus"),
    xmax = paste(Cell.Line,"plus")
  )

p5_2 <- ggplot(df_long, aes(x = Cell_Protocol, y = Pct)) +
  geom_boxplot(aes(fill = Protocol), width = 0.7, outlier.shape = NA) +
  geom_point(aes(color = Protocol),
             position = position_jitter(width = 0.08),
             alpha = 0.6, size = 4) +
  facet_grid(Marker ~ Age, scales = "free_y") +
  stat_pvalue_manual(
    stats,
    label = "p.adj.signif",
    xmin = "xmin", xmax = "xmax",
    y.position = "y.position",
    size = 14,
    bracket.size = 1.4,
    tip.length = 0.01
  ) +
  scale_y_continuous(limits = c(0, 100)) +
  big_theme +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "Percent (%)")

############################################
### 12W: SATB2/CTIP2/SOX9/KI67 Percentages
############################################

df2 <- read.csv("~/project/IPSC_2025_Data/12W_cell_counts.csv")

df_long <- df2 %>%
  pivot_longer(
    cols = c(SATB2.Pct, KI67.Pct, SOX9.Pct, CTIP2.Pct),
    names_to = "Marker",
    values_to = "Pct"
  ) %>%
  mutate(
    Marker = recode(Marker,
                    "SATB2.Pct" = "SATB2+Pct",
                    "CTIP2.Pct" = "CTIP2+Pct",
                    "SOX9.Pct"  = "SOX9+Pct",
                    "KI67.Pct"  = "KI67+Pct"
    ),
    Cell_Protocol = paste(Cell.Line, Protocol),
    Cell_Protocol = factor(Cell_Protocol,
                           levels = c("JHC1 minus","JHC1 plus",
                                      "KOLF2.1J minus","KOLF2.1J plus",
                                      "O2C3 minus","O2C3 plus")
    ),
    Protocol = factor(Protocol, levels = c("minus","plus")),
    Marker = factor(Marker, levels = c("SATB2+Pct","CTIP2+Pct","SOX9+Pct","KI67+Pct"))
  )

stats <- df_long %>%
  group_by(Cell.Line, Marker, Age) %>%
  wilcox_test(Pct ~ Protocol) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  left_join(
    df_long %>%
      group_by(Cell.Line, Marker, Age) %>%
      summarise(y.position = safe_ypos(cur_data()), .groups = "drop"),
    by = c("Cell.Line","Marker","Age")
  ) %>%
  mutate(
    xmin = paste(Cell.Line,"minus"),
    xmax = paste(Cell.Line,"plus")
  )

p5_3 <- ggplot(df_long, aes(x = Cell_Protocol, y = Pct)) +
  geom_boxplot(aes(fill = Protocol), width = 0.7, outlier.shape = NA) +
  geom_point(aes(color = Protocol),
             position = position_jitter(width = 0.08),
             alpha = 0.6, size = 4) +
  facet_grid(Marker ~ Age, scales = "fixed") +
  stat_pvalue_manual(
    stats,
    label = "p.adj.signif",
    xmin = "xmin", xmax = "xmax",
    y.position = "y.position",
    size = 14,
    bracket.size = 1.4,
    tip.length = 0.01
  ) +
  scale_y_continuous(limits = c(0, 100)) +
  big_theme +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "Percent (%)")

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
      "CTIP2.KI67..CTIP2." = "CTIP2+KI67+/CTIP2+Pct",
      "SATB2.KI67..SATB2." = "SATB2+KI67+/SATB2+Pct",
      "SOX9.KI67..SOX9."   = "SOX9+KI67+/SOX9+Pct"
    ),
    Cell_Protocol = paste(Cell.Line, Protocol),
    Cell_Protocol = factor(Cell_Protocol,
                           levels = c("JHC1 minus","JHC1 plus",
                                      "KOLF2.1J minus","KOLF2.1J plus",
                                      "O2C3 minus","O2C3 plus")
    ),
    Protocol = factor(Protocol, levels = c("minus","plus"))
  )

stats <- df_long %>%
  filter(!is.na(Pct)) %>%
  group_by(Cell.Line, Marker, Age) %>%
  wilcox_test(Pct ~ Protocol) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  left_join(
    df_long %>%
      group_by(Cell.Line, Marker, Age) %>%
      summarise(y.position = safe_ypos(cur_data()), .groups = "drop"),
    by = c("Cell.Line","Marker","Age")
  ) %>%
  mutate(
    xmin = paste(Cell.Line,"minus"),
    xmax = paste(Cell.Line,"plus")
  )

p5_4 <- ggplot(df_long, aes(x = Cell_Protocol, y = Pct)) +
  geom_boxplot(aes(fill = Protocol), width = 0.7, outlier.shape = NA) +
  geom_point(aes(color = Protocol),
             position = position_jitter(width = 0.08),
             alpha = 0.6, size = 4) +
  facet_grid(Marker ~ Age, scales = "free_y") +
  stat_pvalue_manual(
    stats,
    label = "p.adj.signif",
    xmin = "xmin", xmax = "xmax",
    y.position = "y.position",
    size = 14,
    bracket.size = 1.4,
    tip.length = 0.01
  ) +
  scale_y_continuous(limits = c(0, 100)) +
  big_theme +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "Percent (%)")

############################################
### 2×2 Combined Figure
############################################

p5 <- ggarrange(
  p5_1, p5_3,
  p5_2, p5_4,
  ncol = 2, nrow = 2,
  font.label = list(size = 30)
)

ggsave("~/project/IPSC_2025_Data/Figure2e.tiff",
       plot = p5,
       device = "tiff",
       width = 20, height = 20, dpi = 300)



