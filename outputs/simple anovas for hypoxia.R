setwd("~/Desktop/2026_hypoxia/outputs") # Please use "inputs" as working directory

# Read in csv for some simple anovas
vial_output_for_anova <- read.csv("21_output_for_anova.csv", header=TRUE)

# Convert variable "Cross" to a categorical factor
vial_output_for_anova$cross_id <- factor(vial_output_for_anova$cross_id)
vial_output_for_anova$data_collector <- factor(vial_output_for_anova$data_collector)
vial_output_for_anova$treatment_id <- factor(vial_output_for_anova$treatment_id)
vial_output_for_anova$brood_id <- factor(vial_output_for_anova$brood_id)

# Perform ANOVA on individual vial estimates fitting an intercept only model
vials_pooled_anova_input <- as.data.frame(vial_output_for_anova)
vials_pooled_anova.lm <- lm(vials_pooled_anova_input$MLE_X_length ~ brood_id*data_collector*treatment_id, data=vials_pooled_anova_input)
vials_pooled_anova.table <- anova(vials_pooled_anova.lm)
vials_pooled_emmeans_table <- emmeans(vials_pooled_anova.lm, ~treatment_id)
vials_pooled_contrast_table <- pairs(vials_pooled_emmeans_table)

write.csv(vials_pooled_anova.table,"../outputs/table_S10_vials_pooled_anova.csv",row.names = TRUE)
write.csv(vials_pooled_contrast_table,"../outputs/vials_pooled_contrast_table.csv",row.names = TRUE)


pooled_by_cross_anova_input <- read.csv("pooled_by_cross_anova_input.csv", header=TRUE)

# Perform ANOVA on cross pooled estimates fitting an intercept only model
cross_pooled_anova_input <- as.data.frame(pooled_by_cross_anova_input)
cross_pooled_anova.lm <- lm(cross_pooled_anova_input$MLE_X_length ~ treatment_id*data_collector, data=cross_pooled_anova_input)
cross_pooled_anova_anova.table <- anova(cross_pooled_anova.lm)
cross_pooled_emmeans_table <- emmeans(cross_pooled_anova.lm, ~treatment_id)
cross_pooled_contrast_table <- pairs(cross_pooled_emmeans_table)
write.csv(cross_pooled_anova_anova.table,"../outputs/table_S11_cross_pooled_anova.csv",row.names = TRUE)
write.csv(cross_pooled_contrast_table,"../outputs/cross_pooled_contrast_table.csv",row.names = TRUE)



pooled_by_brood_anova_input <- read.csv("pooled_by_brood_anova_input.csv", header=TRUE)

# Perform ANOVA on brood pooled estimates fitting an intercept only model
pooled_by_brood_anova_input$cross_id <- factor(pooled_by_brood_anova_input$cross_id)
pooled_by_brood_anova_input$data_collector <- factor(pooled_by_brood_anova_input$data_collector)
pooled_by_brood_anova_input$treatment_id <- factor(pooled_by_brood_anova_input$treatment_id)
pooled_by_brood_anova_input$brood_id <- factor(pooled_by_brood_anova_input$brood_id)

brood_pooled_anova_input <- as.data.frame(pooled_by_brood_anova_input)
brood_pooled_anova.lm <- lm(brood_pooled_anova_input$MLE_X_length ~ treatment_id, data=brood_pooled_anova_input)
brood_pooled_anova_anova.table <- anova(brood_pooled_anova.lm)
brood_pooled_emmeans_table <- emmeans(brood_pooled_anova.lm, ~treatment_id)
brood_pooled_contrast_table <- pairs(brood_pooled_emmeans_table)
write.csv(brood_pooled_anova_anova.table,"../outputs/table_S12_brood_pooled_anova.csv",row.names = TRUE)
write.csv(brood_pooled_contrast_table,"../outputs/brood_pooled_contrast_table.csv",row.names = TRUE)





fecundity_anova <- read.csv("fecundity_anova_input.csv", header=TRUE)

fecundity_anova <- as.data.frame(fecundity_anova)
fecundity_anova$cross_id <- factor(fecundity_anova$cross_id)
fecundity_anova$data_collector <- factor(fecundity_anova$data_collector)
fecundity_anova$treatment_id <- factor(fecundity_anova$treatment_id)
fecundity_anova$brood_id <- factor(fecundity_anova$brood_id)

fecundity_anova.lm <- lm(fecundity_anova$egg_count ~ markers*treatment_id*brood_id, data=fecundity_anova)
fecundity_anova.table <- anova(fecundity_anova.lm)
fecundity_anova_emmeans_table <- emmeans(fecundity_anova.lm, ~treatment_id)
fecundity_anova_contrasts_table <-pairs(fecundity_anova_emmeans_table)
write.csv(fecundity_anova.table,"../outputs/fecundity_anova.csv",row.names = TRUE)
write.csv(fecundity_anova_contrasts_table,"../outputs/fecundity_contrasts.csv",row.names = TRUE)



viability_anova <- read.csv("viability_anova_input.csv", header=TRUE)

viability_anova <- as.data.frame(viability_anova)
viability_anova$cross_id <- factor(viability_anova$cross_id)
viability_anova$data_collector <- factor(viability_anova$data_collector)
viability_anova$treatment_id <- factor(viability_anova$treatment_id)
viability_anova$brood_id <- factor(viability_anova$brood_id)

viability_anova.lm <- lm(viability_anova$viability ~ markers*treatment_id*brood_id*sex, data=viability_anova)
viability_anova.table <- anova(viability_anova.lm)
viability_anova_emmeans_table <- emmeans(viability_anova.lm, ~treatment_id)
viability_anova_contrasts_table <-pairs(viability_anova_emmeans_table)
write.csv(viability_anova.table,"../outputs/viability_anova.csv",row.names = TRUE)
write.csv(viability_anova_contrasts_table,"../outputs/viability_contrasts.csv",row.names = TRUE)
