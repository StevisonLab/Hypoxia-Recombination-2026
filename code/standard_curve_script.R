setwd("~/GitHub/Hypoxia-Recombination-2026/")

# ============================================================
# Figure S2
#   (a) digital egg-count standard curves
#   (b) egg-to-adult viability vs. per-vial egg density by cross type
#
# (a) Type II (major axis, MA) regression of each observer's digital egg
# count on the microscope count for the same vials (calibration
# dataset of Schneider et al. 2025). The MA slope and intercept are
# used in 2_hypoxia_dataset_compilation_script.R to correct the
# digital counts: corrected = (digital - intercept) / slope.
# ============================================================

library(lmodel2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

pac_comparison <- read.csv("inputs/S1_hypoxia_standard_curves.csv", fileEncoding = "UTF-8-BOM")

EH_regression <- lmodel2(EH_digital_egg_count ~ SK_microscope_egg_count, data = pac_comparison)
EH_regression
CH_regression <- lmodel2(CH_digital_egg_count ~ SK_microscope_egg_count, data = pac_comparison)
CH_regression

get_ma <- function(model) {
  res <- model$regression.results
  res[res$Method == "MA", c("Intercept", "Slope")]
}
ma_coefs <- bind_rows(CH = get_ma(CH_regression), EH = get_ma(EH_regression), .id = "observer")
ma_coefs

# Check against the coefficients hard-coded in 2_hypoxia_dataset_compilation_script.R
hard_coded <- data.frame(observer  = c("CH", "EH"),
                         Intercept = c(2.738337, 2.383755),
                         Slope     = c(0.8962081, 0.9459906))
stopifnot(isTRUE(all.equal(ma_coefs$Intercept, hard_coded$Intercept, tolerance = 1e-5)),
          isTRUE(all.equal(ma_coefs$Slope,     hard_coded$Slope,     tolerance = 1e-5)))

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
observer_colors <- c(CH = "#c0392b", EH = "#3266ad")

counts_long <- pac_comparison %>%
  pivot_longer(c(CH_digital_egg_count, EH_digital_egg_count),
               names_to = "observer", values_to = "digital_count") %>%
  mutate(observer = sub("_digital_egg_count", "", observer))

x_range <- range(pac_comparison$SK_microscope_egg_count, na.rm = TRUE)
ma_lines <- ma_coefs %>%
  crossing(x = x_range) %>%
  mutate(y = Intercept + Slope * x)

eq_labels <- ma_coefs %>%
  mutate(label = sprintf("%s: y = %.3fx + %.2f", observer, Slope, Intercept),
         x = x_range[1],
         y = max(counts_long$digital_count, na.rm = TRUE) * c(1, 0.93))

panel_a <- ggplot(counts_long, aes(x = SK_microscope_egg_count, y = digital_count,
                                     color = observer)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(size = 2, alpha = 0.8) +
  geom_line(data = ma_lines, aes(x = x, y = y), linewidth = 0.9) +
  geom_text(data = eq_labels, aes(x = x, y = y, label = label),
            hjust = 0, size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = observer_colors) +
  coord_equal() +
  labs(x = "Microscope egg count",
       y = "Digital egg count",
       color = "Observer") +
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "inside", legend.position.inside = c(0.88, 0.15))

# (b) Viability vs. egg density, on the angular-transformed scale analyzed in
# the ANCOVA (one point per vial per F2 sex; datasets written by script 9).
# NRR_2 = marker-free crosses; NRR_1 = multiply-marked crosses.
cross_colors <- c("Marker-free" = "#1b7837", "Multiply-marked" = "#762a83")

viability <- bind_rows(
  read.csv("outputs/figure_three_b_marker_free_survival_dataset.csv"),
  read.csv("outputs/figure_three_b_multiply_marked_survival_dataset.csv")
) %>%
  mutate(Cross_type = recode(Markers, NRR_2 = "Marker-free", NRR_1 = "Multiply-marked"))

panel_b <- ggplot(viability, aes(x = Density, y = Arcsine, color = Cross_type)) +
  geom_point(size = 1, alpha = 0.35) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.9) +
  scale_color_manual(values = cross_colors) +
  labs(x = "Eggs per vial",
       y = "Egg-to-adult viability (angular transformed)",
       color = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(color = "black"),
        legend.position = "inside", legend.position.inside = c(0.8, 0.1))

figure_S2 <- panel_a + panel_b +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold"))

ggsave("outputs/figure_S2_standard_curves.pdf", figure_S2, width = 10, height = 5,
       device = cairo_pdf)
ggsave("outputs/figure_S2_standard_curves.png", figure_S2, width = 10, height = 5, dpi = 300)
