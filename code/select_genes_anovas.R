setwd("~/Desktop/2026_hypoxia/inputs") # Please use "inputs" as working directory

library(car)
library(emmeans)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpattern)
options(contrasts = c("contr.sum", "contr.poly")) 

#Select gene ANOVAs

select_genes_input <- read.csv("hypoxia_select_genes.csv", header=TRUE)

Hph_anova.lm <- lm(select_genes_input$Hph ~ Treatment*Day, data=select_genes_input)
Hph_anova.table <- anova(Hph_anova.lm)
Hph_emmeans_table <- emmeans(Hph_anova.lm, ~Treatment)
Hph_contrasts_table <- pairs(Hph_emmeans_table)
write.csv(Hph_anova.table,"../outputs/Hph_anova_table.csv",row.names = TRUE)
write.csv(Hph_emmeans_table,"../outputs/Hph_emmeans_table.csv",row.names = TRUE)
write.csv(Hph_contrasts_table,"../outputs/Hph_contrasts_table.csv",row.names = TRUE)

figure_two_select_gene_order <- factor(select_genes_input$Treatment, levels = c("Control", "Normoxia", "Hypoxia"))
figure_two_select_gene_order <- as.numeric(figure_two_select_gene_order)
figure_two_select_gene_values <-cbind(select_genes_input,figure_two_select_gene_order)

ggplot(figure_two_select_gene_values, aes(x = reorder(Treatment,figure_two_select_gene_order), y = Hph, fill = Treatment)) + # Map variables to aesthetics, use 'fill' for color by group
  geom_boxplot() +
  coord_cartesian(ylim = c(30, 60)) +
  scale_fill_manual(values = c("Control" = "#bf812d", "Normoxia" = "#01665e", "Hypoxia" = "#c51b7d")) + 
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
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "fatiga", y = "Expression (fpkm)", x = "Experimental Treatment")



sima_anova.lm <- lm(select_genes_input$sima ~ Treatment*Day, data=select_genes_input)
sima_anova.table <- anova(sima_anova.lm)
sima_emmeans_table <- emmeans(sima_anova.lm, ~Treatment)
sima_contrasts_table <- pairs(sima_emmeans_table)
write.csv(sima_anova.table,"../outputs/sima_anova_table.csv",row.names = TRUE)
write.csv(sima_emmeans_table,"../outputs/sima_emmeans_table.csv",row.names = TRUE)
write.csv(sima_contrasts_table,"../outputs/sima_contrasts_table.csv",row.names = TRUE)

ggplot(figure_two_select_gene_values, aes(x = reorder(Treatment,figure_two_select_gene_order), y = sima, fill = Treatment)) + # Map variables to aesthetics, use 'fill' for color by group
  geom_boxplot() +
  coord_cartesian(ylim = c(20, 50)) +
  scale_fill_manual(values = c("Control" = "#bf812d", "Normoxia" = "#01665e", "Hypoxia" = "#c51b7d")) + 
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
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "sima", y = "Expression (fpkm)", x = "Experimental Treatment")

tgo_anova.lm <- lm(select_genes_input$tgo ~ Treatment*Day, data=select_genes_input)
tgo_anova.table <- anova(tgo_anova.lm)
tgo_emmeans_table <- emmeans(tgo_anova.lm, ~Treatment)
tgo_contrasts_table <- pairs(tgo_emmeans_table)
write.csv(tgo_anova.table,"../outputs/tgo_anova_table.csv",row.names = TRUE)
write.csv(tgo_emmeans_table,"../outputs/tgo_emmeans_table.csv",row.names = TRUE)
write.csv(tgo_contrasts_table,"../outputs/tgo_contrasts_table.csv",row.names = TRUE)

ggplot(figure_two_select_gene_values, aes(x = reorder(Treatment,figure_two_select_gene_order), y = tgo, fill = Treatment)) + # Map variables to aesthetics, use 'fill' for color by group
  geom_boxplot() +
  coord_cartesian(ylim = c(20,50)) +
  scale_fill_manual(values = c("Control" = "#bf812d", "Normoxia" = "#01665e", "Hypoxia" = "#c51b7d")) + 
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
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "tango", y = "Expression (fpkm)", x = "Experimental Treatment")

Vhl_anova.lm <- lm(select_genes_input$Vhl ~ Treatment*Day, data=select_genes_input)
Vhl_anova.table <- anova(Vhl_anova.lm)
Vhl_emmeans_table <- emmeans(Vhl_anova.lm, ~Treatment)
Vhl_contrasts_table <- pairs(Vhl_emmeans_table)
write.csv(Vhl_anova.table,"../outputs/Vhl_anova_table.csv",row.names = TRUE)
write.csv(Vhl_emmeans_table,"../outputs/Vhl_emmeans_table.csv",row.names = TRUE)
write.csv(Vhl_contrasts_table,"../outputs/Vhl_contrasts_table.csv",row.names = TRUE)

ggplot(figure_two_select_gene_values, aes(x = reorder(Treatment,figure_two_select_gene_order), y = Vhl, fill = Treatment)) + # Map variables to aesthetics, use 'fill' for color by group
  geom_boxplot() +
  coord_cartesian(ylim = c(0,15)) +
  scale_fill_manual(values = c("Control" = "#bf812d", "Normoxia" = "#01665e", "Hypoxia" = "#c51b7d")) + 
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
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Vhl", y = "Expression (fpkm)", x = "Experimental Treatment")


ERR_anova.lm <- lm(select_genes_input$ERR ~ Treatment*Day, data=select_genes_input)
ERR_anova.table <- anova(ERR_anova.lm)
ERR_emmeans_table <- emmeans(ERR_anova.lm, ~Treatment)
ERR_contrasts_table <- pairs(ERR_emmeans_table)
write.csv(ERR_anova.table,"../outputs/ERR_anova_table.csv",row.names = TRUE)
write.csv(ERR_emmeans_table,"../outputs/ERR_emmeans_table.csv",row.names = TRUE)
write.csv(ERR_contrasts_table,"../outputs/ERR_contrasts_table.csv",row.names = TRUE)


ggplot(figure_two_select_gene_values, aes(x = reorder(Treatment,figure_two_select_gene_order), y = ERR, fill = Treatment)) + # Map variables to aesthetics, use 'fill' for color by group
  geom_boxplot() +
  coord_cartesian(ylim = c(0,15)) +
  scale_fill_manual(values = c("Control" = "#bf812d", "Normoxia" = "#01665e", "Hypoxia" = "#c51b7d")) + 
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
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "ERR", y = "Expression (fpkm)", x = "Experimental Treatment")



