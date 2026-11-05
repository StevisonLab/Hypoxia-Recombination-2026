# Per-vial H3 counting-model fits: same model as 8_hypoxia_pooling_analysis.R, applied to the 247 vials
# in outputs/13_multi_locus_individual_vials_dataset.csv. Each fit is written to its own file named by
# vial ID, so estimates cannot be separated from their vial. Run from "inputs"; vials already fit are skipped.
# Usage: Rscript ../code/8b_hypoxia_per_vial_fits.R <worker> <n_workers>   (about 3 min per vial per worker)
# Then run 8c_hypoxia_build_anova_inputs.R.
library(dfoptim)
out_dir <- "../outputs/21_per_vial_fits"
dir.create(out_dir, showWarnings = FALSE)
source("../code/4_hypoxia_multi_locus_custom_functions.R")
MLE_H3_looped_vials_outputs <- matrix(nrow=0, ncol=39)
colnames(MLE_H3_looped_vials_outputs)<-c("ID",
                                         "MLE_X_length",
                                         "MLE_x1",
                                         "MLE_x2",
                                         "MLE_x3",
                                         "MLE_x4",
                                         "MLE_x5",
                                         "MLE_p",
                                         "MLE_y1",
                                         "MLE_y2",
                                         "MLE_y3",
                                         "MLE_y4",
                                         "MLE_y5",
                                         "MLE_v1f+",
                                         "MLE_v2f+",
                                         "MLE_v3f+",
                                         "MLE_v4f+",
                                         "MLE_v5f+",
                                         "MLE_v6f+",
                                         "MLE_v1f-",
                                         "MLE_v2f-",
                                         "MLE_v3f-",
                                         "MLE_v4f-",
                                         "MLE_v5f-",
                                         "MLE_v6f-",
                                         "MLE_v1m+",
                                         "MLE_v2m+",
                                         "MLE_v3m+",
                                         "MLE_v4m+",
                                         "MLE_v5m+",
                                         "MLE_v6m+",
                                         "MLE_v1m-",
                                         "MLE_v2m-",
                                         "MLE_v3m-",
                                         "MLE_v4m-",
                                         "MLE_v5m-",
                                         "MLE_v6m-",
                                         "lnL",
                                         "Convergence")

vials <- read.csv("../outputs/13_multi_locus_individual_vials_dataset.csv", check.names = FALSE)[, -1]
multi_locus_full_dataset <- vials
args <- as.integer(commandArgs(TRUE))
todo <- which(!file.exists(file.path(out_dir, paste0(vials$ID, ".csv"))))
todo <- todo[seq(args[1], length(todo), by = args[2])]
for (g in todo) {
  single_experimental_unit <- as.numeric(multi_locus_full_dataset[g,])
  observed_count <- c(single_experimental_unit[130],single_experimental_unit[2:129])
  
  
  
  # Create matrix to store all optimized outputs of Uniform.CxCoM.Likelihood
  # Log likelihoods calculated for each combination of y & p model parameters
  MLE_H3_looped_vials_internal <- matrix(nrow=0, ncol=39)
  colnames(MLE_H3_looped_vials_internal)<-c("ID",
                                            "MLE_X_length",
                                            "MLE_x1",
                                            "MLE_x2",
                                            "MLE_x3",
                                            "MLE_x4",
                                            "MLE_x5",
                                            "MLE_p",
                                            "MLE_y1",
                                            "MLE_y2",
                                            "MLE_y3",
                                            "MLE_y4",
                                            "MLE_y5",
                                            "MLE_v1f+",
                                            "MLE_v2f+",
                                            "MLE_v3f+",
                                            "MLE_v4f+",
                                            "MLE_v5f+",
                                            "MLE_v6f+",
                                            "MLE_v1f-",
                                            "MLE_v2f-",
                                            "MLE_v3f-",
                                            "MLE_v4f-",
                                            "MLE_v5f-",
                                            "MLE_v6f-",
                                            "MLE_v1m+",
                                            "MLE_v2m+",
                                            "MLE_v3m+",
                                            "MLE_v4m+",
                                            "MLE_v5m+",
                                            "MLE_v6m+",
                                            "MLE_v1m-",
                                            "MLE_v2m-",
                                            "MLE_v3m-",
                                            "MLE_v4m-",
                                            "MLE_v5m-",
                                            "MLE_v6m-",
                                            "lnL",
                                            "Convergence")
  
  # For loop Uniform.CxCoM.Likelihood for all possible combinations of p
  # Output model parameters and likelihoods in matrix MLE_H3_full_outputs
  for (v in 1:max_p_value) {
    
    # Define parameters starting points for y, v, and p in single run of Nelder-Mead algorithm
    # When bounding Nelder-Mead using package dfoptim, starting points cannot be at the boundaries
    parameters <- as.numeric(c(y_uniform[v,1],
                               y_uniform[v,2],
                               y_uniform[v,3],
                               y_uniform[v,4],
                               y_uniform[v,5],
                               0.738,
                               0.741,
                               0.739,
                               0.772,
                               0.817,
                               0.766,
                               0.619,
                               0.615,
                               0.617,
                               0.584,
                               0.540,
                               0.591,
                               0.606,
                               0.606,
                               0.619,
                               0.673,
                               0.748,
                               0.657,
                               0.479,
                               0.479,
                               0.466,
                               0.412,
                               0.337,
                               0.428))
    
    p <- p_uniform[v,]
    
    # Optimize CxCoM model parameters y and v with bounded Nelder-Mead method
    # Upper limit on y's is arbitrary, whereas upper and lower limit on v's is
    # given by definition of egg-to-adult viability (cannot exceed 0-1 range)
    # nmkb stops with an error if the likelihood returns NaN; log that p and go on to the next one
    MLE_H3_Vial_CxCoM_Table <- tryCatch(nmkb(par = parameters,
                                    fn = CxCoM.Likelihood.H3,
                                    control = list(maxfeval = 50000),
                                    lower = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                    upper = c(9,9,9,9,9,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)),
                                    error = function(e) { cat("FAILED", vials$ID[g], "p =", p[1], conditionMessage(e), "\n"); NULL })
    if (is.null(MLE_H3_Vial_CxCoM_Table)) next
    
    # Create output vector of y, p, likelihood, convergence code
    MLE_Vial_outputs <- c(g,
                          (((MLE_H3_Vial_CxCoM_Table$par[1]+MLE_H3_Vial_CxCoM_Table$par[2]+MLE_H3_Vial_CxCoM_Table$par[3]+MLE_H3_Vial_CxCoM_Table$par[4]+MLE_H3_Vial_CxCoM_Table$par[5])/(2*p[1]))*100),
                          MLE_H3_Vial_CxCoM_Table$par[1]/(2*p[1]),
                          MLE_H3_Vial_CxCoM_Table$par[2]/(2*p[1]),
                          MLE_H3_Vial_CxCoM_Table$par[3]/(2*p[1]),
                          MLE_H3_Vial_CxCoM_Table$par[4]/(2*p[1]),
                          MLE_H3_Vial_CxCoM_Table$par[5]/(2*p[1]),
                          p[1],
                          MLE_H3_Vial_CxCoM_Table$par[1],
                          MLE_H3_Vial_CxCoM_Table$par[2],
                          MLE_H3_Vial_CxCoM_Table$par[3],
                          MLE_H3_Vial_CxCoM_Table$par[4],
                          MLE_H3_Vial_CxCoM_Table$par[5],
                          MLE_H3_Vial_CxCoM_Table$par[6],
                          MLE_H3_Vial_CxCoM_Table$par[7],
                          MLE_H3_Vial_CxCoM_Table$par[8],
                          MLE_H3_Vial_CxCoM_Table$par[9],
                          MLE_H3_Vial_CxCoM_Table$par[10],
                          MLE_H3_Vial_CxCoM_Table$par[11],
                          MLE_H3_Vial_CxCoM_Table$par[12],
                          MLE_H3_Vial_CxCoM_Table$par[13],
                          MLE_H3_Vial_CxCoM_Table$par[14],
                          MLE_H3_Vial_CxCoM_Table$par[15],
                          MLE_H3_Vial_CxCoM_Table$par[16],
                          MLE_H3_Vial_CxCoM_Table$par[17],
                          MLE_H3_Vial_CxCoM_Table$par[18],
                          MLE_H3_Vial_CxCoM_Table$par[19],
                          MLE_H3_Vial_CxCoM_Table$par[20],
                          MLE_H3_Vial_CxCoM_Table$par[21],
                          MLE_H3_Vial_CxCoM_Table$par[22],
                          MLE_H3_Vial_CxCoM_Table$par[23],
                          MLE_H3_Vial_CxCoM_Table$par[24],
                          MLE_H3_Vial_CxCoM_Table$par[25],
                          MLE_H3_Vial_CxCoM_Table$par[26],
                          MLE_H3_Vial_CxCoM_Table$par[27],
                          MLE_H3_Vial_CxCoM_Table$par[28],
                          MLE_H3_Vial_CxCoM_Table$par[29],
                          MLE_H3_Vial_CxCoM_Table$value,
                          MLE_H3_Vial_CxCoM_Table$convergence)
    
    # Add vector of y, p, likelihood, code to MLE_Uniform_CxCoM
    MLE_H3_looped_vials_internal <- rbind(MLE_H3_looped_vials_internal, MLE_Vial_outputs)
  }
  
  # Take output matrix MLE_Uniform_CxCoM and sort by likelihood
  # Print the sorted table for records, make sure to RENAME FILE
  MLE_H3_vials <- as.data.frame(MLE_H3_looped_vials_internal)
  df4 <- MLE_H3_vials[order(MLE_H3_vials$lnL),]
  MLE_H3_looped_vials_outputs<- rbind(MLE_H3_looped_vials_outputs, df4[1,])
  result <- cbind(vial_unique_id = vials$ID[g], df4[1, ])
  write.csv(result, file.path(out_dir, paste0(vials$ID[g], ".csv")), row.names = FALSE)
}
