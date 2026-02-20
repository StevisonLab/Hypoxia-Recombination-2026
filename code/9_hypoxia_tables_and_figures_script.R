library(car)
library(emmeans)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpattern)
options(contrasts = c("contr.sum", "contr.poly")) 

figure_two_marker_free_egg_count <- marker_free_viability_dataset[,-c(7,8,9,10,11,13,14,15)]
figure_two_marker_free_fecundity <- figure_two_marker_free_egg_count$transformed_mean_marker_free_total_egg_count/5
figure_two_marker_free_egg_dataset <- cbind(figure_two_marker_free_egg_count,figure_two_marker_free_fecundity)
names(figure_two_marker_free_egg_dataset) <- c("Markers", "Cross", "Brood", "Vial", "Treatment", "Observer", "Egg_Count", "Fecundity")

figure_two_multiply_marked_egg_count <- cleaned_multiply_marked_viability_dataset[,-c(7:139,141:143)]
figure_two_multiply_marked_fecundity <- figure_two_multiply_marked_egg_count$transformed_mean_multiply_marked_total_egg_count/5
figure_two_multiply_marked_egg_dataset <- cbind(figure_two_multiply_marked_egg_count,figure_two_multiply_marked_fecundity)
names(figure_two_multiply_marked_egg_dataset) <- c("Markers", "Cross", "Brood", "Vial", "Treatment", "Observer", "Egg_Count", "Fecundity")

figure_two_fecundity_dataset <- rbind(figure_two_marker_free_egg_dataset,figure_two_multiply_marked_egg_dataset)
figure_two_fecundity_input <- as.data.frame(figure_two_fecundity_dataset)
figure_two_fecundity_input$Brood <- factor(figure_two_fecundity_input$Brood)
#figure_two_fecundity_input$Brood <- as.numeric(figure_two_fecundity_input$Brood)
write.csv(figure_two_fecundity_dataset,"../outputs/figure_two_fecundity_dataset.csv",row.names = TRUE)

#Type I sum of squares is sequential testing of model effects use only for balanced designs
#Type II sum of squares is conditional testing of model effects for unbalanced designs assuming no interactions
#Type III sum of squares is conditional testing of model effects for unbalanced designs when interactions are tested
#options(contrasts = c("contr.sum", "contr.poly")) #Contrasts settings necessary for type III sum of squares
options(contrasts = c("contr.treatment", "contr.poly")) #Default contrasts for R (type I sum of squares)

figure_two_fecundity_anova.lm <- lm(figure_two_fecundity_input$Fecundity ~ Brood*Markers*Treatment, data=figure_two_fecundity_input)
figure_two_fecundity_anova.table <- anova(figure_two_fecundity_anova.lm)
figure_two_fecundity_emmeans_table <- emmeans(figure_two_fecundity_anova.lm, ~Markers+Treatment)
figure_two_fecundity_contrasts_table <- pairs(figure_two_fecundity_emmeans_table)
write.csv(figure_two_fecundity_anova.table,"../outputs/figure_two_fecundity_anova_table.csv",row.names = TRUE)
write.csv(figure_two_fecundity_emmeans_table,"../outputs/figure_two_fecundity_emmeans_table.csv",row.names = TRUE)
write.csv(figure_two_fecundity_contrasts_table,"../outputs/figure_two_fecundity_contrasts_table.csv",row.names = TRUE)


figure_two_panel_a_estimates <- as.data.frame(figure_two_fecundity_emmeans_table)
figure_two_panel_a_order<-c(4,1,6,3,5,2)
figure_two_panel_a_markers <- c("Multiply-marked cross","Marker-free cross","Multiply-marked cross","Marker-free cross","Multiply-marked cross","Marker-free cross")
figure_two_panel_a_pattern <- c("stripe","none","stripe","none","stripe","none")
figure_two_panel_a_estimates <- as.data.frame(figure_two_fecundity_emmeans_table)
figure_two_panel_a_estimates <-figure_two_panel_a_estimates[,-c(1,4,5)]
figure_two_panel_a_values <- cbind(figure_two_panel_a_markers,figure_two_panel_a_order,figure_two_panel_a_estimates,figure_two_panel_a_pattern)
names(figure_two_panel_a_values) <- c("Markers", "Order", "Treatment", "emmean", "lower.CL", "upper.CL","Hatching")
figure_two_panel_a_multiply_marked_only <- data.frame(figure_two_panel_a_values[-c(2,4,6),])

ggplot(figure_two_panel_a_values, aes(x = reorder(Treatment,Order), y = emmean, fill=Treatment, pattern = Hatching)) +
  geom_bar(stat = "identity",color = "black", width = 0.9) +
  geom_col_pattern(
    color = "black",        # Outline color
    pattern_fill = "white", # Stripe color (fill color of the pattern)
    pattern_colour = "white", # Stripe outline color (if needed)
    pattern_density = 0.01 # Lower value for narrower stripes
  ) +
  scale_pattern_identity() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  facet_wrap(~ factor(Markers)) +
  coord_cartesian(ylim = c(0, 50)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +  
  scale_fill_manual(values = c("Control" = "#bf812d", "Normoxia" = "#01665e", "Hypoxia" = "#c51b7d")) +  theme_minimal() +
  theme(legend.position="none",
        panel.border = element_blank(), # Remove the default panel border (rectangle)
        axis.line.x = element_line(color = "black", size = 0.5), # Add bottom x-axis line
        axis.line.y = element_line(color = "black", size = 0.5), # Add left y-axis line
        axis.ticks = element_line(color = "black"),
        panel.spacing.x = unit(0, "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "F1 Female Fecundity", x = "Experimental Treatment")

# Code for building a data frame with sex-specific egg-to-adult survival from marker-free crosses from output file "F_marker_free_viability_dataset.csv"
figure_two_marker_free_survival <- marker_free_viability_dataset[,-c(7:11,15)]
figure_two_marker_free_survival_female_fecundity <- figure_two_marker_free_survival$transformed_mean_marker_free_total_egg_count/5
figure_two_marker_free_survival_female_proportions <- figure_two_marker_free_survival$derived_marker_free_female_viability

# Angular transformation of proportions, adjusting all values in excess of 1 to be 1.570796 on transformed scale, and back transformations
figure_two_marker_free_survival_female_arcsine <- asin(sqrt(figure_two_marker_free_survival_female_proportions))
figure_two_marker_free_survival_female_arcsine[is.nan(figure_two_marker_free_survival_female_arcsine)] <- 1.570796
figure_two_marker_free_survival_female_inverse <- (sin(figure_two_marker_free_survival_female_arcsine))^2

# Compiling the data, transformations, and covariates in a labeled data.frame
figure_two_marker_free_survival_female_vector <- rep("Female", 252)
figure_two_marker_free_survival_female_dataset <- cbind(figure_two_marker_free_survival,figure_two_marker_free_survival_female_fecundity,figure_two_marker_free_survival_female_arcsine,figure_two_marker_free_survival_female_inverse,figure_two_marker_free_survival_female_vector)
names(figure_two_marker_free_survival_female_dataset) <- c("Markers", "Cross", "Brood", "Vial", "Treatment", "Observer", "Density", "Female_Proportion", "Male_Proportion", "Fecundity", "Arcsine", "Inverse", "Sex")

# Code for building a data frame with sex-specific egg-to-adult survival from marker-free crosses from output file "F_marker_free_viability_dataset.csv"
figure_two_marker_free_survival_male_fecundity <- figure_two_marker_free_survival$transformed_mean_marker_free_total_egg_count/5
figure_two_marker_free_survival_male_proportions <- figure_two_marker_free_survival$derived_marker_free_male_viability

# Angular transformation of proportions, adjusting all values in excess of 1 to be 1.570796 on transformed scale, and back transformations
figure_two_marker_free_survival_male_arcsine <- asin(sqrt(figure_two_marker_free_survival_male_proportions))
figure_two_marker_free_survival_male_arcsine[is.nan(figure_two_marker_free_survival_male_arcsine)] <- 1.570796
figure_two_marker_free_survival_male_inverse <- (sin(figure_two_marker_free_survival_male_arcsine))^2

# Compiling the data, transformations, and covariates in a labeled data.frame
figure_two_marker_free_survival_male_vector <- rep("Male", 252)
figure_two_marker_free_survival_male_dataset <- cbind(figure_two_marker_free_survival,figure_two_marker_free_survival_male_fecundity,figure_two_marker_free_survival_male_arcsine,figure_two_marker_free_survival_male_inverse,figure_two_marker_free_survival_male_vector)
names(figure_two_marker_free_survival_male_dataset) <- c("Markers", "Cross", "Brood", "Vial", "Treatment", "Observer", "Density", "Female_Proportion", "Male_Proportion", "Fecundity", "Arcsine", "Inverse", "Sex")

# Combining the female and male data.frames with the option of coding brood as continuous versus categorical
figure_two_marker_free_survival_dataset <- rbind(figure_two_marker_free_survival_female_dataset,figure_two_marker_free_survival_male_dataset)
figure_two_marker_free_survival_input <- as.data.frame(figure_two_marker_free_survival_dataset)
figure_two_marker_free_survival_input$Brood <- factor(figure_two_marker_free_survival_input$Brood)
#figure_two_marker_free_survival_input$Brood <- as.numeric(figure_two_marker_free_survival_input$Brood)
write.csv(figure_two_marker_free_survival_dataset,"../outputs/figure_two_marker_free_survival_dataset.csv",row.names = TRUE)

#Type I sum of squares is sequential testing of model effects use only for balanced designs
#Type II sum of squares is conditional testing of model effects for unbalanced designs assuming no interactions
#Type III sum of squares is conditional testing of model effects for unbalanced designs when interactions are tested
options(contrasts = c("contr.sum", "contr.poly")) #Contrasts settings necessary for type III sum of squares
#options(contrasts = c("contr.treatment", "contr.poly")) #Default contrasts for R (type I sum of squares

# Performing ANCOVA on Egg-to-Adult survival with F2 density per vial as a covariate
figure_two_marker_free_survival_ancova.lm <- aov(Arcsine ~ Density+Sex*Brood*Treatment, data=figure_two_marker_free_survival_input)
figure_two_marker_free_survival_ancova.table <- Anova(figure_two_marker_free_survival_ancova.lm, type="III")
figure_two_marker_free_survival_emmeans_table <- emmeans(figure_two_marker_free_survival_ancova.lm, ~Treatment)
figure_two_marker_free_survival_contrasts_table <- pairs(figure_two_marker_free_survival_emmeans_table)
write.csv(figure_two_marker_free_survival_ancova.table,"../outputs/figure_two_marker_free_survival_ancova_table.csv",row.names = TRUE)
write.csv(figure_two_marker_free_survival_emmeans_table,"../outputs/figure_two_marker_free_survival_emmeans_table.csv",row.names = TRUE)
write.csv(figure_two_marker_free_survival_contrasts_table,"../outputs/figure_two_marker_free_survival_contrasts_table.csv",row.names = TRUE)

# Converting estimated marginal means from ANCOVA into a data.frame for figure two
figure_two_panel_c_transformed <- as.data.frame(figure_two_marker_free_survival_emmeans_table)
figure_two_panel_c_treatments <- as.factor(figure_two_panel_c_transformed$Treatment)
figure_two_panel_c_transformed_estimates <- as.data.frame(figure_two_panel_c_transformed[,-c(1,3,4)])
figure_two_panel_c_back_transformed_estimates <- (sin(figure_two_panel_c_transformed_estimates))^2
figure_two_panel_c_values <- cbind(c(1,3,2),figure_two_panel_c_treatments,figure_two_panel_b_back_transformed_estimates)
names(figure_two_panel_c_values) <- c("Order", "Treatment", "emmean", "lower.CL", "upper.CL")

# Code for a four panel figure two, commented out
# ggplot(figure_two_panel_c_values, aes(x = reorder(Treatment,Order), y = emmean, fill=Treatment)) + geom_bar(stat = "identity",color = "black") + coord_cartesian(ylim = c(0, 1)) + geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) + scale_fill_manual(values = c("Control" = "#bf812d", "Normoxia" = "#01665e", "Hypoxia" = "#c51b7d")) + theme_minimal() + theme(legend.position="none") + labs(title = "Marker-free crosses", y = "F2 Egg-to-Adult Survival", x = "Experimental Treatment")

# Code for building a data frame with sex-specific egg-to-adult survival from multiply marked crosses from output file "6_multiply_marked_viability_dataset.csv"
figure_two_multiply_marked_survival <- cleaned_multiply_marked_viability_dataset[,-c(7:139,143)]
figure_two_multiply_marked_survival_female_fecundity <- figure_two_multiply_marked_survival$transformed_mean_multiply_marked_total_egg_count/5
figure_two_multiply_marked_survival_female_proportions <- figure_two_multiply_marked_survival$derived_multiply_marked_female_viability

# Angular transformation of proportions, adjusting all values in excess of 1 to be 1.570796 on transformed scale, and back transformations
figure_two_multiply_marked_survival_female_arcsine <- asin(sqrt(figure_two_multiply_marked_survival_female_proportions))
figure_two_multiply_marked_survival_female_arcsine[is.nan(figure_two_multiply_marked_survival_female_arcsine)] <- 1.570796
figure_two_multiply_marked_survival_female_inverse <- (sin(figure_two_multiply_marked_survival_female_arcsine))^2

# Compiling the data, transformations, and covariates in a labeled data.frame
figure_two_multiply_marked_survival_female_vector <- rep("Female", 247)
figure_two_multiply_marked_survival_female_dataset <- cbind(figure_two_multiply_marked_survival,figure_two_multiply_marked_survival_female_fecundity,figure_two_multiply_marked_survival_female_arcsine,figure_two_multiply_marked_survival_female_inverse,figure_two_multiply_marked_survival_female_vector)
names(figure_two_multiply_marked_survival_female_dataset) <- c("Markers", "Cross", "Brood", "Vial", "Treatment", "Observer", "Density", "Female_Proportion", "Male_Proportion", "Fecundity", "Arcsine", "Inverse", "Sex")

# Code for building a data frame with sex-specific egg-to-adult survival from multiply marked crosses from output file "6_multiply_marked_viability_dataset.csv"
figure_two_multiply_marked_survival_male_fecundity <- figure_two_multiply_marked_survival$transformed_mean_multiply_marked_total_egg_count/5
figure_two_multiply_marked_survival_male_proportions <- figure_two_multiply_marked_survival$derived_multiply_marked_male_viability

# Angular transformation of proportions, adjusting all values in excess of 1 to be 1.570796 on transformed scale, and back transformations
figure_two_multiply_marked_survival_male_arcsine <- asin(sqrt(figure_two_multiply_marked_survival_male_proportions))
figure_two_multiply_marked_survival_male_arcsine[is.nan(figure_two_multiply_marked_survival_male_arcsine)] <- 1.570796
figure_two_multiply_marked_survival_male_inverse <- (sin(figure_two_multiply_marked_survival_male_arcsine))^2

# Compiling the data, transformations, and covariates in a labeled data.frame
figure_two_multiply_marked_survival_male_vector <- rep("Male", 247)
figure_two_multiply_marked_survival_male_dataset <- cbind(figure_two_multiply_marked_survival,figure_two_multiply_marked_survival_male_fecundity,figure_two_multiply_marked_survival_male_arcsine,figure_two_multiply_marked_survival_male_inverse,figure_two_multiply_marked_survival_male_vector)
names(figure_two_multiply_marked_survival_male_dataset) <- c("Markers", "Cross", "Brood", "Vial", "Treatment", "Observer", "Density", "Female_Proportion", "Male_Proportion", "Fecundity", "Arcsine", "Inverse", "Sex")

# Combining the female and male data.frames with the option of coding brood as continuous versus categorical
figure_two_multiply_marked_survival_dataset <- rbind(figure_two_multiply_marked_survival_female_dataset,figure_two_multiply_marked_survival_male_dataset)
figure_two_multiply_marked_survival_input <- as.data.frame(figure_two_multiply_marked_survival_dataset)
figure_two_multiply_marked_survival_input$Brood <- factor(figure_two_multiply_marked_survival_input$Brood)
#figure_two_multiply_marked_survival_input$Brood <- as.numeric(figure_two_multiply_marked_survival_input$Brood)
write.csv(figure_two_multiply_marked_survival_dataset,"../outputs/figure_two_multiply_marked_survival_dataset.csv",row.names = TRUE)

#Type I sum of squares is sequential testing of model effects use only for balanced designs
#Type II sum of squares is conditional testing of model effects for unbalanced designs assuming no interactions
#Type III sum of squares is conditional testing of model effects for unbalanced designs when interactions are tested
options(contrasts = c("contr.sum", "contr.poly")) #Contrasts settings necessary for type III sum of squares
#options(contrasts = c("contr.treatment", "contr.poly")) #Default contrasts for R (type I sum of squares

# Performing ANCOVA on Egg-to-Adult survival with F2 density per vial as a covariate
figure_two_multiply_marked_survival_ancova.lm <- aov(Arcsine ~ Density+Sex*Brood*Treatment, data=figure_two_multiply_marked_survival_input)
figure_two_multiply_marked_survival_ancova.table <- Anova(figure_two_multiply_marked_survival_ancova.lm, type="III")
#ANOVA equivalent if the nonsignificant density covariate is removed, commented out, no difference observed in inferential statistics
#figure_two_multiply_marked_survival_anova.lm <- aov(Arcsine ~ Sex*Brood*Treatment, data=figure_two_multiply_marked_survival_input)
#figure_two_multiply_marked_survival_anova.table <- Anova(figure_two_multiply_marked_survival_anova.lm, type="III")
#figure_two_multiply_marked_survival_emmeans_table <- emmeans(figure_two_multiply_marked_survival_anova.lm, ~Treatment)
#figure_two_multiply_marked_survival_contrasts_table <- pairs(figure_two_multiply_marked_survival_emmeans_table)
figure_two_multiply_marked_survival_emmeans_table <- emmeans(figure_two_multiply_marked_survival_ancova.lm, ~Treatment)
figure_two_multiply_marked_survival_contrasts_table <- pairs(figure_two_multiply_marked_survival_emmeans_table)
write.csv(figure_two_multiply_marked_survival_ancova.table,"../outputs/figure_two_multiply_marked_survival_ancova_table.csv",row.names = TRUE)
write.csv(figure_two_multiply_marked_survival_emmeans_table,"../outputs/figure_two_multiply_marked_survival_emmeans_table.csv",row.names = TRUE)
write.csv(figure_two_multiply_marked_survival_contrasts_table,"../outputs/figure_two_multiply_marked_survival_contrasts_table.csv",row.names = TRUE)

# Converting estimated marginal means from ANCOVA into a data.frame for figure two
figure_two_panel_d_transformed <- as.data.frame(figure_two_multiply_marked_survival_emmeans_table)
figure_two_panel_d_treatments <- as.factor(figure_two_panel_d_transformed$Treatment)
figure_two_panel_d_transformed_estimates <- as.data.frame(figure_two_panel_d_transformed[,-c(1,3,4)])
figure_two_panel_d_back_transformed_estimates <- (sin(figure_two_panel_d_transformed_estimates))^2
figure_two_panel_d_values <- cbind(c(4,6,5),figure_two_panel_b_treatments,figure_two_panel_d_back_transformed_estimates)
names(figure_two_panel_d_values) <- c("Order", "Treatment", "emmean", "lower.CL", "upper.CL")

# Code for a four panel figure two, commented out
# ggplot(figure_two_panel_d_values, aes(x = reorder(Treatment,Order), y = emmean, fill=Treatment)) + geom_bar(stat = "identity",color = "black") + coord_cartesian(ylim = c(0, 1)) + geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +  scale_fill_manual(values = c("Control" = "#bf812d", "Normoxia" = "#01665e", "Hypoxia" = "#c51b7d")) +  theme_minimal() + theme(legend.position="none") + labs(title = "Multiply-marked crosses", y = "F2 Egg-to-Adult Survival", x = "Experimental Treatment")

figure_two_panels_cd_values_combined <- rbind(figure_two_panel_c_values,figure_two_panel_d_values)
figure_two_panel_b_markers <- c("Marker-free cross","Marker-free cross","Marker-free cross","Multiply-marked cross","Multiply-marked cross","Multiply-marked cross")
figure_two_panel_b_pattern <- c("none","none","none","stripe","stripe","stripe")
figure_two_panel_b_values <- cbind(figure_two_panel_b_markers,figure_two_panels_cd_values_combined,figure_two_panel_b_pattern)
names(figure_two_panel_b_values) <- c("Markers", "Order", "Treatment", "emmean", "lower.CL", "upper.CL", "Hatching")

ggplot(figure_two_panel_b_values, aes(x = reorder(Treatment,Order), y = emmean, fill=Treatment, pattern = Hatching)) +
  geom_bar(stat = "identity",color = "black", width = 0.9) +
  geom_col_pattern(
    color = "black",        # Outline color
    pattern_fill = "white", # Stripe color (fill color of the pattern)
    pattern_colour = "white", # Stripe outline color (if needed)
    pattern_density = 0.01 # Lower value for narrower stripes
  ) +
  scale_pattern_identity() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  facet_wrap(~ factor(Markers)) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +  
  scale_fill_manual(values = c("Control" = "#bf812d", "Normoxia" = "#01665e", "Hypoxia" = "#c51b7d")) +  theme_minimal() +
  theme(legend.position="none",
        panel.border = element_blank(), # Remove the default panel border (rectangle)
        axis.line.x = element_line(color = "black", size = 0.5), # Add bottom x-axis line
        axis.line.y = element_line(color = "black", size = 0.5), # Add left y-axis line
        axis.ticks = element_line(color = "black"),
        panel.spacing.x = unit(0, "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "F2 Egg-to-Adult Survival", x = "Experimental Treatment")

#Type I sum of squares is sequential testing of model effects use only for balanced designs
#Type II sum of squares is conditional testing of model effects for unbalanced designs assuming no interactions
#Type III sum of squares is conditional testing of model effects for unbalanced designs when interactions are tested
options(contrasts = c("contr.sum", "contr.poly")) #Contrasts settings necessary for type III sum of squares
#options(contrasts = c("contr.treatment", "contr.poly")) #Default contrasts for R (type I sum of squares)

MLE_H3_looped_vials_outputs <- read.csv("21_output_for_anova.csv", header=TRUE)
figure_three_length_dataset <- as.data.frame(MLE_H3_looped_vials_outputs)
figure_three_length_dataset$brood_id <- as.factor(figure_three_length_dataset$brood_id)
figure_three_length.lm <- lm(MLE_X_length ~ brood_id*treatment_id*data_collector, data=figure_three_length_dataset)
figure_three_length_anova_table <- anova(figure_three_length.lm)
figure_three_length_emmeans_table <- emmeans(figure_three_length.lm, ~treatment_id)
figure_three_length_contrasts_table <- pairs(figure_three_length_emmeans_table)
write.csv(figure_three_length_anova_table,"../outputs/figure_three_length_anova_table.csv",row.names = TRUE)
write.csv(figure_three_length_emmeans_table,"../outputs/figure_three_length_emmeans_table.csv",row.names = TRUE)
write.csv(figure_three_length_contrasts_table,"../outputs/figure_three_length_contrasts_table.csv",row.names = TRUE)


#figure_three_panel_a_estimates <- as.data.frame(figure_three_length_emmeans_table)
#figure_three_panel_a_estimates <- figure_three_panel_a_estimates[,-c(3,4)]
#figure_three_panel_a_order<-c(1,3,2)
#figure_three_panel_a_pattern <- c("stripe","stripe","stripe")
#figure_three_panel_a_values <- cbind(figure_three_panel_a_order,figure_three_panel_a_estimates,figure_three_panel_a_pattern)
#names(figure_three_panel_a_values) <- c("Order", "Treatment", "emmean", "lower.CL", "upper.CL","Hatching")

figure_three_panel_a_order <- factor(figure_three_length_dataset$treatment_id, levels = c("Control", "Normoxia", "Hypoxia"))
figure_three_panel_a_order <- as.numeric(figure_three_panel_a_order)
figure_three_panel_a_pattern <- rep("stripe", 247)
figure_three_panel_a_values <-cbind(figure_three_length_dataset,figure_three_panel_a_order,figure_three_panel_a_pattern)

ggplot(figure_three_panel_a_values, aes(x = reorder(treatment_id,figure_three_panel_a_order), y = MLE_X_length, fill = treatment_id, pattern = figure_three_panel_a_pattern)) + # Map variables to aesthetics, use 'fill' for color by group
  geom_boxplot() +
  coord_cartesian(ylim = c(25, 125)) +
  scale_fill_manual(values = c("Control" = "#bf812d", "Normoxia" = "#01665e", "Hypoxia" = "#c51b7d")) + 
  geom_boxplot_pattern(
    pattern = 'stripe',
    color = "black",        # Outline color
    pattern_size = 0.1,
    pattern_fill = "white", # Stripe color (fill color of the pattern)
    pattern_colour = "white", # Stripe outline color (if needed)
    pattern_spacing = 0.02,
    pattern_density = 0.1# Lower value for narrower stripes
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_minimal() +
  theme(legend.position="none",
        panel.border = element_blank(), # Remove the default panel border (rectangle)
        axis.line.x = element_line(color = "black", size = 0.5), # Add bottom x-axis line
        axis.line.y = element_line(color = "black", size = 0.5), # Add left y-axis line
        axis.ticks = element_line(color = "black"),
        panel.spacing.x = unit(0, "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "X Chromosome Genetic Length (cM)", x = "Experimental Treatment")


figure_two_combined_survival_dataset <- rbind(figure_two_marker_free_survival_input,figure_two_multiply_marked_survival_input)
#Type I sum of squares is sequential testing of model effects use only for balanced designs
#Type II sum of squares is conditional testing of model effects for unbalanced designs assuming no interactions
#Type III sum of squares is conditional testing of model effects for unbalanced designs when interactions are tested
options(contrasts = c("contr.sum", "contr.poly")) #Contrasts settings necessary for type III sum of squares
#options(contrasts = c("contr.treatment", "contr.poly")) #Default contrasts for R (type I sum of squares)


# Performing ANCOVA on Egg-to-Adult survival with F2 density per vial as a covariate
figure_two_combined_survival_ancova.lm <- aov(Arcsine ~ Density*Sex*Brood*Markers*Treatment, data=figure_two_combined_survival_dataset)
figure_two_combined_survival_ancova.table <- Anova(figure_two_combined_survival_ancova.lm, type="III")
figure_two_combined_survival_emmeans_table <- emmeans(figure_two_combined_survival_ancova.lm, ~Density*Markers*Treatment)
figure_two_combined_survival_contrasts_table <- pairs(figure_two_combined_survival_emmeans_table)
write.csv(figure_two_combined_survival_ancova.table,"../outputs/figure_two_combined_survival_ancova_table.csv",row.names = TRUE)
write.csv(figure_two_combined_survival_emmeans_table,"../outputs/figure_two_combined_survival_emmeans_table.csv",row.names = TRUE)
write.csv(figure_two_combined_survival_contrasts_table,"../outputs/figure_two_combined_survival_contrasts_table.csv",row.names = TRUE)
