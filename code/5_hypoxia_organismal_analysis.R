#setwd("Raw Data/") # Please use Raw Data as working directory

# Perform Regression and write ANOVA table for supplemental material
marker_free_viability_anova.lm <- lm(Viability ~ Brood*Sex*Treatment, data=marker_free_viability_anova_input)
marker_free_viability_anova.table <- anova(marker_free_viability_anova.lm)
marker_free_viability_regression_coefficients <- summary(marker_free_viability_anova.lm)
marker_free_viability_emmeans_table <- emmeans(marker_free_viability_anova.lm, ~Treatment*Sex)
marker_free_contrasts_table<-pairs(marker_free_viability_emmeans_table)
write.csv(marker_free_viability_anova.table,"../outputs/table_S1_marker_free_viability_regression.csv",row.names = TRUE)


#INPUT FILE NOT AVAILABLE
# Perform Regression and write ANOVA table for supplemental material
#multiply_marked_viability_anova.lm <- lm(Viability ~ Brood + Sex + Cross, data=cleaned_multiply_marked_viability_anova_input)
#multiply_marked_viability_anova.table <- anova(multiply_marked_viability_anova.lm)
#multiply_marked_viability_regression_coefficients <- summary(multiply_marked_viability_anova.lm)
#rite.csv(multiply_marked_viability_anova.table,"../outputs/table_S2_multiply_marked_viability_regression.csv",row.names = TRUE)
