# ============================================================================
# DOT PLOT: GO Pathway Enrichment — Up vs Down Classification
# Uses .all_GOenrich files; classifies terms by majority direction
# Combined Up + Down panel plot | High DPI | No transparency
# ============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)   # for combining plots side-by-side

# ============================================================================
# PART 1: LOAD FILES
# ============================================================================

setwd("~/GitHub/Hypoxia-Recombination-2026/")

Hypoxia_vs_Control  <- read.csv("inputs/HypoxiavsControl.all_GOenrich_significant.csv")
Normoxia_vs_Control <- read.csv("inputs/NormoxiavsControl.all_GOenrich_significant.csv")

# ============================================================================
# PART 2: ADD COMPARISON LABELS
# ============================================================================

Hypoxia_vs_Control$Comparison  <- "Chamber 8% vs Control 21%"
Normoxia_vs_Control$Comparison <- "Chamber 21% vs Control 21%"

# ============================================================================
# PART 3: COMBINE, FILTER (padj < 0.05), CLASSIFY DIRECTION
# ============================================================================
# Classification rule:
#   - If Up   > 50% of Count  → "Up"
#   - If Down > 50% of Count  → "Down"
#   - Otherwise               → "Mixed" (excluded from both plots)
# ============================================================================

all_data <- bind_rows(
  Hypoxia_vs_Control,
  Normoxia_vs_Control,
  #  Normoxia_vs_Hypoxia
) %>%
  filter(padj < 0.05) %>%
  mutate(
    # Ensure numeric
    Up   = as.numeric(Up),
    Down = as.numeric(Down),
    Count = as.numeric(Count),
    
    # Classify direction
    Direction = case_when(
      Up   / Count > 0.5 ~ "Up",
      Down / Count > 0.5 ~ "Down",
      TRUE               ~ "Mixed"
    )
  )

cat("Total significant entries:   ", nrow(all_data), "\n")
cat("Classified as Up:            ", sum(all_data$Direction == "Up"), "\n")
cat("Classified as Down:          ", sum(all_data$Direction == "Down"), "\n")
cat("Classified as Mixed (excluded):", sum(all_data$Direction == "Mixed"), "\n\n")

up_data   <- filter(all_data, Direction == "Up")
down_data <- filter(all_data, Direction == "Down")

# ============================================================================
# PART 4: SHARED PREP FUNCTION
# ============================================================================

prepare_dotplot_data <- function(df) {
  df %>%
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
      ),
      Comparison = factor(Comparison, levels = c(
        "Chamber 21% vs Control 21%",
        "Chamber 21% vs Chamber 8%",
        "Chamber 8% vs Control 21%"
      ))
    )
}

# ============================================================================
# PART 5: ORDER GO TERMS BY MAX SIGNIFICANCE (per direction separately)
# ============================================================================

get_term_order <- function(df) {
  df %>%
    group_by(Description_short) %>%
    summarise(max_sig = max(neg_log10_padj, na.rm = TRUE), .groups = "drop") %>%
    arrange(max_sig) %>%
    pull(Description_short)
}

# ============================================================================
# PART 6: PLOT FUNCTION
# ============================================================================

make_dotplot <- function(df, direction, title_label, show_legend = TRUE) {
  
  df <- prepare_dotplot_data(df)
  term_order <- get_term_order(df)
  df <- df %>%
    mutate(Description_short = factor(Description_short, levels = term_order))
  
  color_scale <- if (direction == "Up") {
    scale_color_gradientn(
      colors = c("#FCBBA1", "#FC4E2A", "#800026"),
      name = "-log10(padj)"
    )
  } else {
    scale_color_gradientn(
      colors = c("#BDD7E7", "#2171B5", "#08306B"),
      name = "-log10(padj)"
    )
  }
  
  p <- ggplot(
    df,
    aes(
      x     = GeneRatio_numeric,
      y     = Description_short,
      size  = Count,
      color = neg_log10_padj
    )
  ) +
    geom_point(alpha = 1) +
    
    facet_wrap(~Comparison, ncol = 3) +
    
    color_scale +
    
    scale_size_continuous(
      name  = "Gene Count",
      range = c(3, 12)
    ) +
    
    scale_x_continuous(
      name   = "Gene Ratio",
      labels = scales::percent_format(accuracy = 1)
    ) +
    
    labs(
      title = title_label,
      y     = "GO Pathway"
    ) +
    
    theme_minimal(base_size = 20) +
    theme(
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
      strip.background = element_rect(color = "black", fill = "grey95", linewidth = 0.8),
      strip.text       = element_text(size = 14, face = "bold"),
      axis.text.y      = element_text(size = 16, face = "bold"),
      axis.text.x      = element_text(size = 16, face = "bold"),
      axis.title       = element_text(size = 14, face = "bold"),
      plot.title       = element_text(size = 16, face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      legend.title     = element_text(face = "bold"),
      legend.position  = if (show_legend) "right" else "none"
    )
  
  n_terms <- length(term_order)
  fig_h   <- max(10, n_terms * 0.30 + 4)
  
  return(list(plot = p, height = fig_h, n_terms = n_terms))
}

# ============================================================================
# PART 7: GENERATE INDIVIDUAL PLOTS
# ============================================================================

up_result   <- make_dotplot(up_data,   "Up",   "GO Enrichment — Up-Regulated Genes (>50% Up)")
down_result <- make_dotplot(down_data, "Down", "GO Enrichment — Down-Regulated Genes (>50% Down)")

# Save individual plots
ggsave(
  "outputs/dotplot_UP_regulated_GO.png",
  up_result$plot,
  width  = 20,
  height = up_result$height,
  dpi    = 300,
  bg     = "white"
)
cat("✓ Saved: dotplot_UP_regulated_GO.png  [20 x", round(up_result$height, 1), "in]\n")

ggsave(
  "outputs/dotplot_DOWN_regulated_GO.png",
  down_result$plot,
  width  = 20,
  height = down_result$height,
  dpi    = 300,
  bg     = "white"
)
cat("✓ Saved: dotplot_DOWN_regulated_GO.png  [20 x", round(down_result$height, 1), "in]\n")

# ============================================================================
# PART 8: COMBINED STACKED PLOT (Up on top, Down below)
# ============================================================================
# Use patchwork to stack vertically; heights are proportional to term count

combined <- up_result$plot / down_result$plot +
  plot_layout(heights = c(up_result$n_terms, down_result$n_terms)) +
  plot_annotation(
    title   = "GO Pathway Enrichment — Up & Down Regulated Genes",
    theme   = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )
  )

combined_height <- up_result$height + down_result$height

ggsave(
  "outputs/dotplot_COMBINED_UP_DOWN_GO.png",
  combined,
  width  = 20,
  height = combined_height,
  dpi    = 300,
  bg     = "white"
)
cat("✓ Saved: dotplot_COMBINED_UP_DOWN_GO.png  [20 x", round(combined_height, 1), "in]\n")

cat("\nDone! 3 files saved.\n")
