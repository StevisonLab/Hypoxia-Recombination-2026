setwd("~/GitHub/Hypoxia-Recombination-2026/")

# ============================================================================
# DOT PLOT: Separate Up & Down Regulated GO Pathway Plots
# High DPI | No transparency | Separate ordering per direction
# ============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================================
# PART 1: LOAD FILES
# ============================================================================

Hypoxia_vs_Control_up       <- read.csv("inputs/HypoxiavsControl.up_GOenrich_significant.csv")
Normoxia_vs_Control_up      <- read.csv("inputs/NormoxiavsControl.up_GOenrich_significant.csv")

Hypoxia_vs_Control_down     <- read.csv("inputs/HypoxiavsControl.down_GOenrich_significant.csv")
Normoxia_vs_Control_down    <- read.csv("inputs/NormoxiavsControl.down_GOenrich_significant.csv")

# ============================================================================
# PART 2: ADD LABELS
# ============================================================================

label_and_tag <- function(df, comparison, direction) {
  df$Comparison <- comparison
  df$Direction  <- direction
  return(df)
}

Hypoxia_vs_Control_up    <- label_and_tag(Hypoxia_vs_Control_up,    "Chamber 8% vs Control 21%",  "Up")
Normoxia_vs_Control_up   <- label_and_tag(Normoxia_vs_Control_up,   "Chamber 21% vs Control 21%", "Up")

Hypoxia_vs_Control_down  <- label_and_tag(Hypoxia_vs_Control_down,  "Chamber 8% vs Control 21%",  "Down")
Normoxia_vs_Control_down <- label_and_tag(Normoxia_vs_Control_down, "Chamber 21% vs Control 21%", "Down")

# ============================================================================
# PART 3: COMBINE AND FILTER (padj < 0.05)
# ============================================================================

up_data <- bind_rows(
  Hypoxia_vs_Control_up,
  Normoxia_vs_Control_up,
) %>% filter(padj < 0.05)

down_data <- bind_rows(
  Hypoxia_vs_Control_down,
  Normoxia_vs_Control_down,
) %>% filter(padj < 0.05)

cat("Up-regulated significant entries:  ", nrow(up_data), "\n")
cat("Down-regulated significant entries:", nrow(down_data), "\n")

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
    arrange(max_sig) %>%   # ascending so most significant appears at top of plot
    pull(Description_short)
}

# ============================================================================
# PART 6: SHARED PLOT FUNCTION
# ============================================================================

make_dotplot <- function(df, direction, dot_color, title_label) {
  
  # Prepare data
  df <- prepare_dotplot_data(df)
  
  # Get term ordering
  term_order <- get_term_order(df)
  df <- df %>%
    mutate(Description_short = factor(Description_short, levels = term_order))
  
  # Auto-scale figure height
  n_terms <- length(term_order)
  fig_h   <- max(10, n_terms * 0.30 + 4)
  
  p <- ggplot(
    df,
    aes(
      x    = GeneRatio_numeric,
      y    = Description_short,
      size = Count,
      color = neg_log10_padj   # color intensity = significance
    )
  ) +
    geom_point(alpha = 1) +    # NO transparency
    
    facet_wrap(~Comparison, ncol = 3) +
    
    # Color gradient for significance
    scale_color_gradientn(
      colors = if (direction == "Up") {
        c("#FCBBA1", "#FC4E2A", "#800026")   # light red → deep red
      } else {
        c("#BDD7E7", "#2171B5", "#08306B")   # light blue → deep blue
      },
      name = "-log10(padj)"
    ) +
    
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
    
    theme_minimal(base_size = 20) +        # increased base font size
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
      legend.title     = element_text(face = "bold", size = 14),
      legend.text      = element_text(size = 13),
      legend.position  = "right"
    )
  
  return(list(plot = p, height = fig_h))
}

# ============================================================================
# PART 7: GENERATE AND SAVE PLOTS
# ============================================================================

# --- UP REGULATED PLOT ---
up_result <- make_dotplot(
  df          = up_data,
  direction   = "Up",
  dot_color   = "red",
  title_label = "GO Pathway Enrichment — Up-Regulated Genes"
)

ggsave(
  "outputs/dotplot_UP_regulated_GO.png",
  up_result$plot,
  width  = 20,
  height = up_result$height,
  dpi    = 300,
  bg     = "white"
)

cat("✓ Up-regulated dot plot saved: dotplot_UP_regulated_GO.png\n")
cat("  Figure size: 20 x", round(up_result$height, 1), "inches | 300 DPI\n")

# --- DOWN REGULATED PLOT ---
down_result <- make_dotplot(
  df          = down_data,
  direction   = "Down",
  dot_color   = "blue",
  title_label = "GO Pathway Enrichment — Down-Regulated Genes"
)

ggsave(
  "outputs/dotplot_DOWN_regulated_GO.png",
  down_result$plot,
  width  = 20,
  height = down_result$height,
  dpi    = 300,
  bg     = "white"
)

cat("✓ Down-regulated dot plot saved: dotplot_DOWN_regulated_GO.png\n")
cat("  Figure size: 20 x", round(down_result$height, 1), "inches | 300 DPI\n")
