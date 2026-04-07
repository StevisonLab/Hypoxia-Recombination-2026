# ============================================================================
# SET WORKING DIRECTORY
# ============================================================================
setwd("~/GitHub/Hypoxia-Recombination-2026/")

# ============================================================================
# LOAD LIBRARIES
# ============================================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(grid)
library(gtable)

# ============================================================================
# LOAD ALL COMPARISON FILES
# ============================================================================
Hypoxia_vs_Control_all <- read.csv("inputs/HypoxiavsControl.all_KEGGenrich_significant.csv")
Normoxia_vs_Control_all <- read.csv("inputs/NormoxiavsControl.all_KEGGenrich_significant.csv")
Normoxia_vs_Hypoxia_all <- read.csv("inputs/NormoxiavsHypoxia.all_KEGGenrich_significant.csv")

# ============================================================================
# FUNCTION TO CREATE UP/DOWN DATAFRAMES WITH LABELS
# ============================================================================
prepare_comparison <- function(df, comparison_label) {
  # Up-regulated
  up_df <- df %>%
    filter(Up > 0) %>%
    mutate(
      Direction = "Up-Regulated",
      Comparison = comparison_label
    )
  
  # Down-regulated
  down_df <- df %>%
    filter(Down > 0) %>%
    mutate(
      Direction = "Down-Regulated",
      Comparison = comparison_label
    )
  
  # Combine
  return(bind_rows(up_df, down_df))
}

# ============================================================================
# APPLY FUNCTION TO EACH COMPARISON
# ============================================================================
up_down_combined <- bind_rows(
  prepare_comparison(Hypoxia_vs_Control_all, "Chamber 8% vs Control 21%"),
  prepare_comparison(Normoxia_vs_Control_all, "Chamber 21% vs Control 21%"),
  # prepare_comparison(Normoxia_vs_Hypoxia_all, "Chamber 21% vs Chamber 8%")
)

# ============================================================================
# FUNCTION TO PREPARE DATA FOR PLOTTING
# ============================================================================
prepare_and_order <- function(df) {
  df <- df %>%
    mutate(
      neg_log10_padj = -log10(padj),
      GeneRatio_numeric = sapply(
        strsplit(as.character(GeneRatio), "/"),
        function(x) as.numeric(x[1]) / as.numeric(x[2])
      ),
      Description_short = ifelse(
        nchar(Description) > 50,
        paste0(substr(Description, 1, 47), "..."),
        Description
      )
    )
  
  term_order <- df %>%
    group_by(Description_short) %>%
    summarise(max_sig = max(neg_log10_padj, na.rm = TRUE), .groups = "drop") %>%
    arrange(max_sig) %>%
    pull(Description_short)
  
  df <- df %>%
    mutate(Description_short = factor(Description_short, levels = term_order))
  
  return(df)
}

up_down_combined <- prepare_and_order(up_down_combined)

# ============================================================================
# FACTORIZE FOR FACETING
# ============================================================================
up_down_combined <- up_down_combined %>%
  mutate(
    Direction = factor(Direction, levels = c("Up-Regulated", "Down-Regulated")),
    Comparison = factor(Comparison, levels = c(
      "Chamber 21% vs Control 21%",
      "Chamber 21% vs Chamber 8%",
      "Chamber 8% vs Control 21%"
    ))
  )

prop_up   = up_down_combined$Up / up_down_combined$Count  
prop_down = up_down_combined$Down / up_down_combined$Count

# Filter to assign pathway to Up or Down if >=50% of genes in that direction
up_down_filtered <- up_down_combined %>%
  rowwise() %>%
  mutate(
    prop_up   = Up / Count,
    prop_down = Down / Count,
    # Assign dominant direction only if >=50%
    Direction_filtered = case_when(
      prop_up >= 0.5   ~ "Up-Regulated",
      prop_down >= 0.5 ~ "Down-Regulated",
      TRUE             ~ NA_character_  # pathways with no clear majority are dropped
    )
  ) %>%
  filter(!is.na(Direction_filtered)) %>%
  ungroup() %>%
  mutate(
    Direction = factor(Direction_filtered, levels = c("Up-Regulated", "Down-Regulated"))
  )
# ============================================================================
# CREATE COMBINED DOTPLOT — DOT SIZE = Up / Down GENES
# ============================================================================
p_combined <- ggplot() +
  # UP-REGULATED DOTS (RED)
  geom_point(
    data = up_down_filtered %>% filter(Direction == "Up-Regulated"),
    aes(
      x = GeneRatio_numeric,
      y = Description_short,
      size = Count,       # use total count
      color = neg_log10_padj
    ),
    alpha = 1
  ) +
  scale_color_gradientn(
    colors = c("#fcae91", "#fb6a4a", "#de2d26", "#a50f15"),
    name = "-log10(padj) (Up)",
    guide = guide_colorbar(order = 1)
  ) +
  ggnewscale::new_scale_color() +
  # DOWN-REGULATED DOTS (BLUE)
  geom_point(
    data = up_down_filtered %>% filter(Direction == "Down-Regulated"),
    aes(
      x = GeneRatio_numeric,
      y = Description_short,
      size = Count,       # use total count
      color = neg_log10_padj
    ),
    alpha = 1
  ) +
  scale_color_gradientn(
    colors = c("#B0E2FF", "#87CEFA", "#3182BD", "#27408B"),
    name = "-log10(padj) (Down)",
    guide = guide_colorbar(order = 2)
  ) +
  facet_grid(Direction ~ Comparison, scales = "free_y", space = "free_y") +
  scale_size_continuous(name = "Gene Count", range = c(3, 12)) +
  scale_x_continuous(name = "Gene Ratio", labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "KEGG Pathway Enrichment — Up & Down Regulated",
    subtitle = "Dot size = Gene Count  •  Color = −log10(padj)",
    y = "KEGG Pathway"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(color = "black", fill = "grey95", linewidth = 0.8),
    strip.text.x = element_text(size = 11, face = "bold"),
    strip.text.y = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 9, face = "bold", color = "black"),
    axis.text.x = element_text(size = 9, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    panel.spacing.y = unit(0.8, "cm"),
    panel.spacing.x = unit(0.5, "cm")
  )

# ============================================================================
# COLOR THE STRIP LABELS (UP = RED, DOWN = BLUE)
# ============================================================================
g <- ggplotGrob(p_combined)
strip_r_idx <- which(grepl("strip-r", g$layout$name))

for (i in seq_along(strip_r_idx)) {
  strip_color <- if (i == 1) "#C0392B" else "#2980B9"
  strip_grob <- g$grobs[[strip_r_idx[i]]]
  for (j in seq_along(strip_grob$grobs)) {
    inner <- strip_grob$grobs[[j]]
    for (k in seq_along(inner$children)) {
      child <- inner$children[[k]]
      if (inherits(child, "rect")) inner$children[[k]]$gp <- gpar(fill = strip_color, col = strip_color)
      if (inherits(child, "text") || inherits(child, "titleGrob")) {
        inner$children[[k]]$gp <- modifyList(as.list(inner$children[[k]]$gp), list(col = "white")) %>% do.call(gpar, .)
      }
    }
    strip_grob$grobs[[j]] <- inner
  }
  g$grobs[[strip_r_idx[i]]] <- strip_grob
}

# ============================================================================
# AUTO HEIGHT & SAVE FIGURE
# ============================================================================
n_up <- length(levels(up_down_combined$Description_short[up_down_combined$Direction == "Up-Regulated"]))
n_down <- length(levels(up_down_combined$Description_short[up_down_combined$Direction == "Down-Regulated"]))
fig_h <- max(12, (n_up + n_down) * 0.28 + 5)

png(
  "outputs/dotplot_COMBINED_KEGG_allcomparisons.png",
  width = 12,
  height = 6,
  units = "in",
  res = 300,
  bg = "white"
)
grid.draw(g)
dev.off()

cat("✓ Combined KEGG dot plot saved: dotplot_COMBINED_KEGG_allcomparisons.png\n")

cat("  Figure size: 20 x", round(fig_h, 1), "inches | 300 DPI\n")
