setwd("~/Desktop/2026_hypoxia")
pac_comparison <- read.csv("S1_hypoxia_standard_curves.csv")
library(lmodel2)
EH_regression <- lmodel2(EH_digital_egg_count ~ SK_microscope_egg_count, data=pac_comparison)
EH_regression
CH_regression <- lmodel2(CH_digital_egg_count ~ SK_microscope_egg_count, data=pac_comparison)
CH_regression