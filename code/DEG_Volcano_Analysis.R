setwd("~/GitHub/Hypoxia-Recombination-2026/")

# Load required libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)
library(patchwork)

cond_hypoxia  <- "Chamber 8%"
cond_normoxia <- "Chamber 21%"
cond_control  <- "Control 21%"

# Import your CSV files
Hypoxia_vs_Control_deg <- read.csv("inputs/HypoxiavsControl_deg.csv", row.names = 1)
Normoxia_vs_Control_deg <- read.csv("inputs/NormoxiavsControl_deg.csv", row.names = 1)
Normoxia_vs_Hypoxia_deg <- read.csv("inputs/NormoxiavsHypoxia_deg.csv", row.names = 1)

# --------------------------
# Legend setup
# --------------------------
legend_colors <- c(
  "#c51b7d",
  "#01665e",
  "#bf812d",
  "gray70"
)

names(legend_colors) <- c(
  paste0("Upregulated in ", cond_hypoxia),
  paste0("Upregulated in ", cond_normoxia),
  paste0("Upregulated in ", cond_control),
  "Not significant"
)

# --------------------------
# Global x- and y-limits
# --------------------------
all_log2FC <- c(Hypoxia_vs_Control_deg$log2FoldChange,
                Normoxia_vs_Control_deg$log2FoldChange,
                Normoxia_vs_Hypoxia_deg$log2FoldChange)

all_neglog10p <- c(-log10(Hypoxia_vs_Control_deg$padj),
                   -log10(Normoxia_vs_Control_deg$padj),
                   -log10(Normoxia_vs_Hypoxia_deg$padj))

x_limits <- c(floor(min(all_log2FC, na.rm = TRUE)), ceiling(max(all_log2FC, na.rm = TRUE)))
y_limits <- c(0, max(all_neglog10p, na.rm = TRUE) * 1.1)

# --------------------------
# Volcano wrapper function
# --------------------------
create_volcano_plot <- function(deg_data, title, condition_right, condition_left,
                                pval_threshold = 0.05,log2fc_threshold = 1, top_genes = 20,
                                ylim = NULL, xlim = NULL) {
  
  deg_data <- deg_data %>%
    mutate(
      log2FC = log2FoldChange,
      neglog10p = -log10(padj),
      significance = case_when(
        log2FC > log2fc_threshold  & padj < pval_threshold ~ paste0("Upregulated in ", condition_right),
        log2FC < -log2fc_threshold & padj < pval_threshold ~ paste0("Upregulated in ", condition_left),
        TRUE ~ "Not significant"
      )
    )
  
  deg_data$significance <- factor(deg_data$significance, levels = names(legend_colors))
  
  deg_data <- deg_data %>% filter(!is.na(log2FC) & !is.na(neglog10p) & !is.infinite(neglog10p))
  
  genes_to_label <- deg_data %>%
    filter(significance != "Not significant") %>%
    arrange(desc(abs(log2FC) * neglog10p)) %>%
    head(top_genes)
  
  p <- ggplot(deg_data, aes(x = log2FC, y = neglog10p)) +
    geom_point(aes(color = significance), alpha = 0.6, size = 4.5) +
    scale_color_manual(values = legend_colors, limits = names(legend_colors), drop = FALSE) +
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), 
               linetype = "dashed", color = "gray40")+
    geom_text_repel(
      data = genes_to_label,
      aes(label = gene_name),
      size = 5,
      box.padding = 0.5,
      max.overlaps = 20,
      segment.color = "gray50"
    ) +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value",
      color = "Legend"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 13, face = "bold"),
      legend.key.size = unit(1.5, "cm"),
      axis.text = element_text(size = 13),
      axis.title = element_text(size = 15, face = "bold")
    )
  
  if (!is.null(ylim) | !is.null(xlim)) {
    p <- p + coord_cartesian(ylim = ylim, xlim = xlim)
  }
  
  return(p)
}

# --------------------------
# Create individual volcano plots
# --------------------------
vol_HC <- create_volcano_plot(
  Hypoxia_vs_Control_deg,
  paste(cond_hypoxia, "vs", cond_control),
  cond_hypoxia, cond_control,
  ylim = y_limits, xlim = x_limits
)

vol_NC <- create_volcano_plot(
  Normoxia_vs_Control_deg,
  paste(cond_normoxia, "vs", cond_control),
  cond_normoxia, cond_control,
  ylim = y_limits, xlim = x_limits
)

vol_NH <- create_volcano_plot(
  Normoxia_vs_Hypoxia_deg,
  paste(cond_normoxia, "vs", cond_hypoxia),
  cond_normoxia, cond_hypoxia,
  ylim = y_limits, xlim = x_limits
)

# Save individual plots
ggsave("outputs/volcano_Hypoxia_vs_Control_padj.png", vol_HC, width = 8, height = 6, dpi = 300)
ggsave("outputs/volcano_Normoxia_vs_Control_padj.png", vol_NC, width = 8, height = 6, dpi = 300)
ggsave("outputs/volcano_Normoxia_vs_Hypoxia_padj.png", vol_NH, width = 8, height = 6, dpi = 300)

# --------------------------
# Combined plots
# --------------------------
library(cowplot)

# Strip legends for combined plots
vol_HC_noleg <- vol_HC + theme(legend.position = "none")
vol_NC_noleg <- vol_NC + theme(legend.position = "none")
vol_NH_noleg <- vol_NH + theme(legend.position = "none")

# Dummy plot for shared legend
dummy <- ggplot(
  data.frame(x = 1:4, y = 1:4, label = factor(names(legend_colors), levels = names(legend_colors))),
  aes(x = x, y = y, color = label)
) +
  geom_point(size = 6) +
  scale_color_manual(values = legend_colors, drop = FALSE) +
  guides(color = guide_legend(title = "Legend", override.aes = list(size = 6))) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 13, face = "bold"),
    legend.key.size = unit(1.5, "cm")
  )

leg <- get_legend(dummy)
stopifnot(!inherits(leg, "zeroGrob"))

# Combined horizontal
combined_horizontal <- plot_grid(
  vol_HC_noleg, vol_NC_noleg, vol_NH_noleg,
  ncol = 3,
  align = "hv"
)

final_horizontal <- plot_grid(
  combined_horizontal, leg,
  ncol = 2,
  rel_widths = c(3, 0.5)
)

ggsave("outputs/volcano_combined_horizontal.png", final_horizontal, width = 22, height = 6, dpi = 300)

# Combined vertical
combined_vertical <- plot_grid(
  vol_HC_noleg, vol_NC_noleg, vol_NH_noleg,
  ncol = 1,
  align = "hv"
)

final_vertical <- plot_grid(
  combined_vertical, leg,
  ncol = 2,
  rel_widths = c(3, 0.5)
)

ggsave("outputs/volcano_combined_vertical.png", final_vertical, width = 12, height = 18, dpi = 300)
