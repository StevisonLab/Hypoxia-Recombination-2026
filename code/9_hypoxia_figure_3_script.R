setwd("~/Desktop/2026_hypoxia/inputs") # Please use "inputs" as working directory

library(car)
library(emmeans)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpattern)
options(contrasts = c("contr.sum", "contr.poly")) 

#Figure 3 Panel A F1 Fecundity

figure_three_a_marker_free_egg_count <- marker_free_viability_dataset[,-c(7,8,9,10,11,13,14,15)]
figure_three_a_marker_free_fecundity <- figure_three_a_marker_free_egg_count$transformed_mean_marker_free_total_egg_count/5
figure_three_a_marker_free_egg_dataset <- cbind(figure_three_a_marker_free_egg_count,figure_three_a_marker_free_fecundity)
names(figure_three_a_marker_free_egg_dataset) <- c("Markers", "Cross", "Brood", "Vial", "Treatment", "Observer", "Egg_Count", "Fecundity")

figure_three_a_multiply_marked_egg_count <- cleaned_multiply_marked_viability_dataset[,-c(7:139,141:147)]
figure_three_a_multiply_marked_fecundity <- figure_three_a_multiply_marked_egg_count$transformed_mean_multiply_marked_total_egg_count/5
figure_three_a_multiply_marked_egg_dataset <- cbind(figure_three_a_multiply_marked_egg_count,figure_three_a_multiply_marked_fecundity)
names(figure_three_a_multiply_marked_egg_dataset) <- c("Markers", "Cross", "Brood", "Vial", "Treatment", "Observer", "Egg_Count", "Fecundity")

figure_three_a_fecundity_dataset <- rbind(figure_three_a_marker_free_egg_dataset,figure_three_a_multiply_marked_egg_dataset)
figure_three_a_fecundity_input <- as.data.frame(figure_three_a_fecundity_dataset)
figure_three_a_fecundity_input$Brood <- factor(figure_three_a_fecundity_input$Brood)
#figure_three_a_fecundity_input$Brood <- as.numeric(figure_three_a_fecundity_input$Brood)
write.csv(figure_three_a_fecundity_dataset,"../outputs/figure_three_a_fecundity_dataset.csv",row.names = TRUE)

#Type I sum of squares is sequential testing of model effects use only for balanced designs
#Type II sum of squares is conditional testing of model effects for unbalanced designs assuming no interactions
#Type III sum of squares is conditional testing of model effects for unbalanced designs when interactions are tested
options(contrasts = c("contr.sum", "contr.poly")) #Contrasts settings necessary for type III sum of squares
#options(contrasts = c("contr.treatment", "contr.poly")) #Default contrasts for R (type I sum of squares)

figure_three_a_fecundity_anova.lm <- lm(figure_three_a_fecundity_input$Fecundity ~ Brood*Markers*Treatment, data=figure_three_a_fecundity_input)
figure_three_a_fecundity_anova.table <- Anova(figure_three_a_fecundity_anova.lm, type="III")
figure_three_a_fecundity_emmeans_table <- emmeans(figure_three_a_fecundity_anova.lm, ~Markers+Treatment)
figure_three_a_fecundity_contrasts_table <- pairs(figure_three_a_fecundity_emmeans_table)
write.csv(figure_three_a_fecundity_anova.table,"../outputs/figure_three_a_fecundity_anova_table.csv",row.names = TRUE)
write.csv(figure_three_a_fecundity_emmeans_table,"../outputs/figure_three_a_fecundity_emmeans_table.csv",row.names = TRUE)
write.csv(figure_three_a_fecundity_contrasts_table,"../outputs/figure_three_a_fecundity_contrasts_table.csv",row.names = TRUE)


figure_three_a_estimates <- as.data.frame(figure_three_a_fecundity_emmeans_table)
figure_three_a_order<-c(4,1,6,3,5,2)
figure_three_a_markers <- c("Multiply-marked cross","Marker-free cross","Multiply-marked cross","Marker-free cross","Multiply-marked cross","Marker-free cross")
figure_three_a_pattern <- c("stripe","none","stripe","none","stripe","none")
figure_three_a_estimates <- as.data.frame(figure_three_a_fecundity_emmeans_table)
figure_three_a_estimates <-figure_three_a_estimates[,-c(1,4,5)]
figure_three_a_values <- cbind(figure_three_a_markers,figure_three_a_order,figure_three_a_estimates,figure_three_a_pattern)
names(figure_three_a_values) <- c("Markers", "Order", "Treatment", "emmean", "lower.CL", "upper.CL","Hatching")
figure_three_a_multiply_marked_only <- data.frame(figure_three_a_values[-c(2,4,6),])

ggplot(figure_three_a_values, aes(x = reorder(Treatment,Order), y = emmean, fill=Treatment, pattern = Hatching)) +
  geom_bar(stat = "identity",color = "black", width = 0.9) +
  geom_col_pattern(
    color = "black",        # Outline color
    pattern_fill = "white", # Stripe color (fill color of the pattern)
    pattern_colour = "white", # Stripe outline color (if needed)
    pattern_density = 0.01 # Lower value for narrower stripes
  ) +
  scale_pattern_identity() +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  facet_wrap(~ factor(Markers)) +
  coord_cartesian(ylim = c(0, 60)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +  
  scale_fill_manual(values = c("Control" = "#bf812d", "Normoxia" = "#01665e", "Hypoxia" = "#c51b7d")) +  theme_minimal() +
  theme(legend.position="none",
        panel.border = element_blank(), # Remove the default panel border (rectangle)
        axis.line.x = element_line(color = "black", size = 0.5), # Add bottom x-axis line
        axis.line.y = element_line(color = "black", size = 0.5), # Add left y-axis line
        axis.ticks = element_line(color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        panel.spacing.x = unit(0, "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Daily Individual Female Fecundity", x = "Experimental Treatment")

# Figure 3 Panel B Survival
# Code for building a data frame with sex-specific egg-to-adult survival from marker-free crosses from output file "F_marker_free_viability_dataset.csv"
figure_three_b_marker_free_survival <- marker_free_viability_dataset[,-c(7:11,15)]
figure_three_b_marker_free_survival_female_fecundity <- figure_three_b_marker_free_survival$transformed_mean_marker_free_total_egg_count/5
figure_three_b_marker_free_survival_female_proportions <- figure_three_b_marker_free_survival$derived_marker_free_female_viability

# Angular transformation of proportions, adjusting all values in excess of 1 to be 1.570796 on transformed scale, and back transformations
figure_three_b_marker_free_survival_female_arcsine <- asin(sqrt(figure_three_b_marker_free_survival_female_proportions))
figure_three_b_marker_free_survival_female_arcsine[is.nan(figure_three_b_marker_free_survival_female_arcsine)] <- 1.570796
figure_three_b_marker_free_survival_female_inverse <- (sin(figure_three_b_marker_free_survival_female_arcsine))^2

# Compiling the data, transformations, and covariates in a labeled data.frame
figure_three_b_marker_free_survival_female_vector <- rep("Female", 252)
figure_three_b_marker_free_survival_female_dataset <- cbind(figure_three_b_marker_free_survival,figure_three_b_marker_free_survival_female_fecundity,figure_three_b_marker_free_survival_female_arcsine,figure_three_b_marker_free_survival_female_inverse,figure_three_b_marker_free_survival_female_vector)
names(figure_three_b_marker_free_survival_female_dataset) <- c("Markers", "Cross", "Brood", "Vial", "Treatment", "Observer", "Density", "Female_Proportion", "Male_Proportion", "Fecundity", "Arcsine", "Inverse", "Sex")

# Code for building a data frame with sex-specific egg-to-adult survival from marker-free crosses from output file "F_marker_free_viability_dataset.csv"
figure_three_b_marker_free_survival_male_fecundity <- figure_three_b_marker_free_survival$transformed_mean_marker_free_total_egg_count/5
figure_three_b_marker_free_survival_male_proportions <- figure_three_b_marker_free_survival$derived_marker_free_male_viability

# Angular transformation of proportions, adjusting all values in excess of 1 to be 1.570796 on transformed scale, and back transformations
figure_three_b_marker_free_survival_male_arcsine <- asin(sqrt(figure_three_b_marker_free_survival_male_proportions))
figure_three_b_marker_free_survival_male_arcsine[is.nan(figure_three_b_marker_free_survival_male_arcsine)] <- 1.570796
figure_three_b_marker_free_survival_male_inverse <- (sin(figure_three_b_marker_free_survival_male_arcsine))^2

# Compiling the data, transformations, and covariates in a labeled data.frame
figure_three_b_marker_free_survival_male_vector <- rep("Male", 252)
figure_three_b_marker_free_survival_male_dataset <- cbind(figure_three_b_marker_free_survival,figure_three_b_marker_free_survival_male_fecundity,figure_three_b_marker_free_survival_male_arcsine,figure_three_b_marker_free_survival_male_inverse,figure_three_b_marker_free_survival_male_vector)
names(figure_three_b_marker_free_survival_male_dataset) <- c("Markers", "Cross", "Brood", "Vial", "Treatment", "Observer", "Density", "Female_Proportion", "Male_Proportion", "Fecundity", "Arcsine", "Inverse", "Sex")

# Combining the female and male data.frames with the option of coding brood as continuous versus categorical
figure_three_b_marker_free_survival_dataset <- rbind(figure_three_b_marker_free_survival_female_dataset,figure_three_b_marker_free_survival_male_dataset)
figure_three_b_marker_free_survival_input <- as.data.frame(figure_three_b_marker_free_survival_dataset)
figure_three_b_marker_free_survival_input$Brood <- factor(figure_three_b_marker_free_survival_input$Brood)
#figure_three_b_marker_free_survival_input$Brood <- as.numeric(figure_three_b_marker_free_survival_input$Brood)
write.csv(figure_three_b_marker_free_survival_dataset,"../outputs/figure_three_b_marker_free_survival_dataset.csv",row.names = TRUE)

# Code for building a data frame with sex-specific egg-to-adult survival from multiply marked crosses from output file "6_multiply_marked_viability_dataset.csv"
figure_three_b_multiply_marked_survival <- cleaned_multiply_marked_viability_dataset[,-c(7:139,143:147)]
figure_three_b_multiply_marked_survival_female_fecundity <- figure_three_b_multiply_marked_survival$transformed_mean_multiply_marked_total_egg_count/5
figure_three_b_multiply_marked_survival_female_proportions <- figure_three_b_multiply_marked_survival$derived_multiply_marked_female_viability

# Angular transformation of proportions, adjusting all values in excess of 1 to be 1.570796 on transformed scale, and back transformations
figure_three_b_multiply_marked_survival_female_arcsine <- asin(sqrt(figure_three_b_multiply_marked_survival_female_proportions))
figure_three_b_multiply_marked_survival_female_arcsine[is.nan(figure_three_b_multiply_marked_survival_female_arcsine)] <- 1.570796
figure_three_b_multiply_marked_survival_female_inverse <- (sin(figure_three_b_multiply_marked_survival_female_arcsine))^2

# Compiling the data, transformations, and covariates in a labeled data.frame
figure_three_b_multiply_marked_survival_female_vector <- rep("Female", 247)
figure_three_b_multiply_marked_survival_female_dataset <- cbind(figure_three_b_multiply_marked_survival,figure_three_b_multiply_marked_survival_female_fecundity,figure_three_b_multiply_marked_survival_female_arcsine,figure_three_b_multiply_marked_survival_female_inverse,figure_three_b_multiply_marked_survival_female_vector)
names(figure_three_b_multiply_marked_survival_female_dataset) <- c("Markers", "Cross", "Brood", "Vial", "Treatment", "Observer", "Density", "Female_Proportion", "Male_Proportion", "Fecundity", "Arcsine", "Inverse", "Sex")

# Code for building a data frame with sex-specific egg-to-adult survival from multiply marked crosses from output file "6_multiply_marked_viability_dataset.csv"
figure_three_b_multiply_marked_survival_male_fecundity <- figure_three_b_multiply_marked_survival$transformed_mean_multiply_marked_total_egg_count/5
figure_three_b_multiply_marked_survival_male_proportions <- figure_three_b_multiply_marked_survival$derived_multiply_marked_male_viability

# Angular transformation of proportions, adjusting all values in excess of 1 to be 1.570796 on transformed scale, and back transformations
figure_three_b_multiply_marked_survival_male_arcsine <- asin(sqrt(figure_three_b_multiply_marked_survival_male_proportions))
figure_three_b_multiply_marked_survival_male_arcsine[is.nan(figure_three_b_multiply_marked_survival_male_arcsine)] <- 1.570796
figure_three_b_multiply_marked_survival_male_inverse <- (sin(figure_three_b_multiply_marked_survival_male_arcsine))^2

# Compiling the data, transformations, and covariates in a labeled data.frame
figure_three_b_multiply_marked_survival_male_vector <- rep("Male", 247)
figure_three_b_multiply_marked_survival_male_dataset <- cbind(figure_three_b_multiply_marked_survival,figure_three_b_multiply_marked_survival_male_fecundity,figure_three_b_multiply_marked_survival_male_arcsine,figure_three_b_multiply_marked_survival_male_inverse,figure_three_b_multiply_marked_survival_male_vector)
names(figure_three_b_multiply_marked_survival_male_dataset) <- c("Markers", "Cross", "Brood", "Vial", "Treatment", "Observer", "Density", "Female_Proportion", "Male_Proportion", "Fecundity", "Arcsine", "Inverse", "Sex")

# Combining the female and male data.frames with the option of coding brood as continuous versus categorical
figure_three_b_multiply_marked_survival_dataset <- rbind(figure_three_b_multiply_marked_survival_female_dataset,figure_three_b_multiply_marked_survival_male_dataset)
figure_three_b_multiply_marked_survival_input <- as.data.frame(figure_three_b_multiply_marked_survival_dataset)
figure_three_b_multiply_marked_survival_input$Brood <- factor(figure_three_b_multiply_marked_survival_input$Brood)
#figure_three_b_multiply_marked_survival_input$Brood <- as.numeric(figure_three_b_multiply_marked_survival_input$Brood)
write.csv(figure_three_b_multiply_marked_survival_dataset,"../outputs/figure_three_b_multiply_marked_survival_dataset.csv",row.names = TRUE)


#Here is the problem
figure_three_b_combined_survival_dataset <- rbind(figure_three_b_marker_free_survival_input,figure_three_b_multiply_marked_survival_input)
#Type I sum of squares is sequential testing of model effects use only for balanced designs
#Type II sum of squares is conditional testing of model effects for unbalanced designs assuming no interactions
#Type III sum of squares is conditional testing of model effects for unbalanced designs when interactions are tested
options(contrasts = c("contr.sum", "contr.poly")) #Contrasts settings necessary for type III sum of squares
#options(contrasts = c("contr.treatment", "contr.poly")) #Default contrasts for R (type I sum of squares)

# Performing ANCOVA on Egg-to-Adult survival with F2 density per vial as a covariate
figure_three_b_combined_survival_ancova.lm <- aov(Arcsine ~ Density+Density*Markers+Sex*Brood*Markers*Treatment, data=figure_three_b_combined_survival_dataset)
figure_three_b_combined_survival_ancova.table <- Anova(figure_three_b_combined_survival_ancova.lm, type="III")
figure_three_b_combined_survival_emmeans_table <- emmeans(figure_three_b_combined_survival_ancova.lm, ~Markers+Treatment)
figure_three_b_combined_survival_contrasts_table <- pairs(figure_three_b_combined_survival_emmeans_table)
write.csv(figure_three_b_combined_survival_ancova.table,"../outputs/figure_three_b_combined_survival_ancova_table.csv",row.names = TRUE)
write.csv(figure_three_b_combined_survival_emmeans_table,"../outputs/figure_three_b_combined_survival_emmeans_table.csv",row.names = TRUE)
write.csv(figure_three_b_combined_survival_contrasts_table,"../outputs/figure_three_b_combined_survival_contrasts_table.csv",row.names = TRUE)

figure_three_b_combined_estimates <- as.data.frame(figure_three_b_combined_survival_emmeans_table)
figure_three_b_combined_order <- c(4,1,6,3,5,2)
figure_three_b_combined_treatments <- c("Control","Control","Hypoxia","Hypoxia","Normoxia","Normoxia")
figure_three_b_combined_markers <- c("Multiply-marked cross","Marker-free cross","Multiply-marked cross","Marker-free cross","Multiply-marked cross","Marker-free cross")
figure_three_b_combined_pattern <- c("stripe","none","stripe","none","stripe","none")
figure_three_b_combined_estimates <- as.data.frame(figure_three_b_combined_survival_emmeans_table)
figure_three_b_combined_estimates <-figure_three_b_combined_estimates[,-c(1,2,4,5)]
figure_three_b_combined_back_transformed_estimates <- (sin(figure_three_b_combined_estimates))^2
figure_three_b_combined_values <- cbind(figure_three_b_combined_markers, figure_three_b_combined_order, figure_three_b_combined_treatments, figure_three_b_combined_back_transformed_estimates, figure_three_b_combined_pattern)
names(figure_three_b_combined_values) <- c("Markers", "Order", "Treatment", "emmean", "lower.CL", "upper.CL","Hatching")

ggplot(figure_three_b_combined_values, aes(x = reorder(Treatment,Order), y = emmean, fill=Treatment, pattern = Hatching)) +
  geom_bar(stat = "identity",color = "black", width = 0.9) +
  geom_col_pattern(
    color = "black",        # Outline color
    pattern_fill = "white", # Stripe color (fill color of the pattern)
    pattern_colour = "white", # Stripe outline color (if needed)
    pattern_density = 0.01 # Lower value for narrower stripes
  ) +
  scale_pattern_identity() +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  facet_wrap(~ factor(Markers)) +
  coord_cartesian(ylim = c(0, 1.25)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +  
  scale_fill_manual(values = c("Control" = "#bf812d", "Normoxia" = "#01665e", "Hypoxia" = "#c51b7d")) +  theme_minimal() +
  theme(legend.position="none",
        panel.border = element_blank(), # Remove the default panel border (rectangle)
        axis.line.x = element_line(color = "black", size = 0.5), # Add bottom x-axis line
        axis.line.y = element_line(color = "black", size = 0.5), # Add left y-axis line
        axis.ticks = element_line(color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        panel.spacing.x = unit(0, "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "Proportion of Eggs Surviving to Eclosion", x = "Experimental Treatment")

#Figure 3 Panel C Recombination

#Type I sum of squares is sequential testing of model effects use only for balanced designs
#Type II sum of squares is conditional testing of model effects for unbalanced designs assuming no interactions
#Type III sum of squares is conditional testing of model effects for unbalanced designs when interactions are tested
options(contrasts = c("contr.sum", "contr.poly")) #Contrasts settings necessary for type III sum of squares
#options(contrasts = c("contr.treatment", "contr.poly")) #Default contrasts for R (type I sum of squares)

MLE_H3_looped_vials_outputs <- read.csv("21_output_for_anova.csv", header=TRUE)
figure_three_c_length_dataset <- as.data.frame(MLE_H3_looped_vials_outputs)
figure_three_c_length_dataset$brood_id <- as.factor(figure_three_c_length_dataset$brood_id)
figure_three_c_length.lm <- lm(MLE_X_length ~ brood_id*treatment_id*data_collector, data=figure_three_c_length_dataset)
figure_three_c_length_anova_table <- Anova(figure_three_c_length.lm, type="III")
figure_three_c_length_emmeans_table <- emmeans(figure_three_c_length.lm, ~treatment_id)
figure_three_c_length_contrasts_table <- pairs(figure_three_c_length_emmeans_table)
write.csv(figure_three_c_length_anova_table,"../outputs/figure_three_c_length_anova_table.csv",row.names = TRUE)
write.csv(figure_three_c_length_emmeans_table,"../outputs/figure_three_c_length_emmeans_table.csv",row.names = TRUE)
write.csv(figure_three_c_length_contrasts_table,"../outputs/figure_three_c_length_contrasts_table.csv",row.names = TRUE)

#Commented out dataframes for using the estimated marginal means in figure three panel c instead of per vial estimates
#figure_three_panel_a_estimates <- as.data.frame(figure_three_length_emmeans_table)
#figure_three_panel_a_estimates <- figure_three_panel_a_estimates[,-c(3,4)]
#figure_three_panel_a_order<-c(1,3,2)
#figure_three_panel_a_pattern <- c("stripe","stripe","stripe")
#figure_three_panel_a_values <- cbind(figure_three_panel_a_order,figure_three_panel_a_estimates,figure_three_panel_a_pattern)
#names(figure_three_panel_a_values) <- c("Order", "Treatment", "emmean", "lower.CL", "upper.CL","Hatching")

figure_three_c_order <- factor(figure_three_c_length_dataset$treatment_id, levels = c("Control", "Normoxia", "Hypoxia"))
figure_three_c_order <- as.numeric(figure_three_c_order)
figure_three_c_pattern <- rep("stripe", 247)
figure_three_c_values <-cbind(figure_three_c_length_dataset,figure_three_c_order,figure_three_c_pattern)

figure_three_c_plot <- ggplot(figure_three_c_values, aes(x = reorder(treatment_id,figure_three_c_order), y = MLE_X_length, fill = treatment_id, pattern = figure_three_c_pattern)) + # Map variables to aesthetics, use 'fill' for color by group
  geom_boxplot() +
  coord_cartesian(ylim = c(25, 150)) +
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
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme_minimal() +
  theme(legend.position="none",
        panel.border = element_blank(), # Remove the default panel border (rectangle)
        axis.line.x = element_line(color = "black", size = 0.5), # Add bottom x-axis line
        axis.line.y = element_line(color = "black", size = 0.5), # Add left y-axis line
        axis.ticks = element_line(color = "black"),
        axis.title.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        panel.spacing.x = unit(0, "line"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(y = "X Chromosome Genetic Length (cM)", x = "Experimental Treatment")

#Exported at the size of the panel placed in the Figure 3 PowerPoint (467 x 660 px at 96 ppi); treatment labels and brackets are added there
ggsave("../outputs/figure_three_c_boxplot.png", figure_three_c_plot, width = 467/96, height = 660/96, dpi = 400, bg = "white")
ggsave("../outputs/figure_three_c_boxplot.pdf", figure_three_c_plot, width = 467/96, height = 660/96)




MLE_H3_single_interval_outputs <- read.csv("21_output_for_anova_single_interval.csv", header=TRUE)
figure_three_c_single_interval_dataset <- as.data.frame(MLE_H3_single_interval_outputs)
figure_three_c_single_interval_dataset$brood_id <- as.factor(figure_three_c_single_interval_dataset$brood_id)

figure_three_c_x1.lm <- lm(MLE_x1 ~ brood_id*treatment_id*data_collector, data=figure_three_c_single_interval_dataset)
figure_three_c_x1_anova_table <- Anova(figure_three_c_x1.lm, type="III")
figure_three_c_x1_emmeans_table <- emmeans(figure_three_c_x1.lm, ~treatment_id)
figure_three_c_x1_contrasts_table <- pairs(figure_three_c_x1_emmeans_table)
write.csv(figure_three_c_x1_anova_table,"../outputs/figure_three_c_x1_anova_table.csv",row.names = TRUE)
write.csv(figure_three_c_x1_emmeans_table,"../outputs/figure_three_c_x1_emmeans_table.csv",row.names = TRUE)
write.csv(figure_three_c_x1_contrasts_table,"../outputs/figure_three_c_x1_contrasts_table.csv",row.names = TRUE)

figure_three_c_x2.lm <- lm(MLE_x2 ~ brood_id*treatment_id*data_collector, data=figure_three_c_single_interval_dataset)
figure_three_c_x2_anova_table <- Anova(figure_three_c_x2.lm, type="III")
figure_three_c_x2_emmeans_table <- emmeans(figure_three_c_x2.lm, ~treatment_id)
figure_three_c_x2_contrasts_table <- pairs(figure_three_c_x2_emmeans_table)
write.csv(figure_three_c_x2_anova_table,"../outputs/figure_three_c_x2_anova_table.csv",row.names = TRUE)
write.csv(figure_three_c_x2_emmeans_table,"../outputs/figure_three_c_x2_emmeans_table.csv",row.names = TRUE)
write.csv(figure_three_c_x2_contrasts_table,"../outputs/figure_three_c_x2_contrasts_table.csv",row.names = TRUE)

figure_three_c_x3.lm <- lm(MLE_x3 ~ brood_id*treatment_id*data_collector, data=figure_three_c_single_interval_dataset)
figure_three_c_x3_anova_table <- Anova(figure_three_c_x3.lm, type="III")
figure_three_c_x3_emmeans_table <- emmeans(figure_three_c_x3.lm, ~treatment_id)
figure_three_c_x3_contrasts_table <- pairs(figure_three_c_x3_emmeans_table)
write.csv(figure_three_c_x3_anova_table,"../outputs/figure_three_c_x3_anova_table.csv",row.names = TRUE)
write.csv(figure_three_c_x3_emmeans_table,"../outputs/figure_three_c_x3_emmeans_table.csv",row.names = TRUE)
write.csv(figure_three_c_x3_contrasts_table,"../outputs/figure_three_c_x3_contrasts_table.csv",row.names = TRUE)

figure_three_c_x4.lm <- lm(MLE_x4 ~ brood_id*treatment_id*data_collector, data=figure_three_c_single_interval_dataset)
figure_three_c_x4_anova_table <- Anova(figure_three_c_x4.lm, type="III")
figure_three_c_x4_emmeans_table <- emmeans(figure_three_c_x4.lm, ~treatment_id)
figure_three_c_x4_contrasts_table <- pairs(figure_three_c_x4_emmeans_table)
write.csv(figure_three_c_x4_anova_table,"../outputs/figure_three_c_x4_anova_table.csv",row.names = TRUE)
write.csv(figure_three_c_x4_emmeans_table,"../outputs/figure_three_c_x4_emmeans_table.csv",row.names = TRUE)
write.csv(figure_three_c_x4_contrasts_table,"../outputs/figure_three_c_x4_contrasts_table.csv",row.names = TRUE)

figure_three_c_x5.lm <- lm(MLE_x5 ~ brood_id*treatment_id*data_collector, data=figure_three_c_single_interval_dataset)
figure_three_c_x5_anova_table <- Anova(figure_three_c_x5.lm, type="III")
figure_three_c_x5_emmeans_table <- emmeans(figure_three_c_x5.lm, ~treatment_id)
figure_three_c_x5_contrasts_table <- pairs(figure_three_c_x5_emmeans_table)
write.csv(figure_three_c_x5_anova_table,"../outputs/figure_three_c_x5_anova_table.csv",row.names = TRUE)
write.csv(figure_three_c_x5_emmeans_table,"../outputs/figure_three_c_x5_emmeans_table.csv",row.names = TRUE)
write.csv(figure_three_c_x5_contrasts_table,"../outputs/figure_three_c_x5_contrasts_table.csv",row.names = TRUE)

#Crossover interference per vial (MLE_p), same model as total genetic length
figure_three_c_interference.lm <- lm(MLE_p ~ brood_id*treatment_id*data_collector, data=figure_three_c_length_dataset)
figure_three_c_interference_anova_table <- Anova(figure_three_c_interference.lm, type="III")
figure_three_c_interference_emmeans_table <- emmeans(figure_three_c_interference.lm, ~treatment_id)
figure_three_c_interference_contrasts_table <- pairs(figure_three_c_interference_emmeans_table)
write.csv(figure_three_c_interference_anova_table,"../outputs/figure_three_c_interference_anova_table.csv",row.names = TRUE)
write.csv(figure_three_c_interference_emmeans_table,"../outputs/figure_three_c_interference_emmeans_table.csv",row.names = TRUE)
write.csv(figure_three_c_interference_contrasts_table,"../outputs/figure_three_c_interference_contrasts_table.csv",row.names = TRUE)



