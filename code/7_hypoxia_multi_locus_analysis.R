#setwd("Raw Data/") # Please use Raw Data as working directory

# Begin algorithm for fitting CxCoM H0 to full experiment dataset
# Create matrix to store all optimized outputs of CxCoM.Likelihood.H0
# Log likelihoods calculated for each combination of y & p model parameters
MLE_H0_full_outputs <- matrix(nrow=0, ncol=39)
colnames(MLE_H0_full_outputs)<-c("ID",
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
# Output model parameters and likelihoods in matrix MLE_H0_full_outputs
for (v in 1:max_p_value) {
  
  # Define parameters starting points for y, v, and p in single run of Nelder-Mead algorithm
  # When bounding Nelder-Mead using package dfoptim, starting points cannot be at the boundaries
  parameters <- as.numeric(c(y_uniform[v,1],
                             y_uniform[v,2],
                             y_uniform[v,3],
                             y_uniform[v,4],
                             y_uniform[v,5]))
  p <- p_uniform[v,]
  
  # Optimize CxCoM model parameters y and v with bounded Nelder-Mead method
  # Upper limit on y's is arbitrary, whereas upper and lower limit on v's is
  # given by definition of egg-to-adult viability (cannot exceed 0-1 range)
  MLE_H0_Full_CxCoM_Table <- nmkb(par = parameters,
                                  fn = CxCoM.Likelihood.H0,
                                  control = list(maxfeval = 50000),
                                  lower = c(0,0,0,0,0),
                                  upper = c(9,9,9,9,9))
  
  # Create output vector of y, p, likelihood, convergence code
  MLE_H0_outputs <- c(0,
                      (((MLE_H0_Full_CxCoM_Table$par[1]+MLE_H0_Full_CxCoM_Table$par[2]+MLE_H0_Full_CxCoM_Table$par[3]+MLE_H0_Full_CxCoM_Table$par[4]+MLE_H0_Full_CxCoM_Table$par[5])/(2*p[1]))*100),
                      MLE_H0_Full_CxCoM_Table$par[1]/(2*p[1]),
                      MLE_H0_Full_CxCoM_Table$par[2]/(2*p[1]),
                      MLE_H0_Full_CxCoM_Table$par[3]/(2*p[1]),
                      MLE_H0_Full_CxCoM_Table$par[4]/(2*p[1]),
                      MLE_H0_Full_CxCoM_Table$par[5]/(2*p[1]),
                      p[1],
                      MLE_H0_Full_CxCoM_Table$par[1],
                      MLE_H0_Full_CxCoM_Table$par[2],
                      MLE_H0_Full_CxCoM_Table$par[3],
                      MLE_H0_Full_CxCoM_Table$par[4],
                      MLE_H0_Full_CxCoM_Table$par[5],
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      MLE_H0_Full_CxCoM_Table$value,
                      MLE_H0_Full_CxCoM_Table$convergence)
  
  # Add vector of y, p, likelihood, code to MLE_Uniform_CxCoM
  MLE_H0_full_outputs <- rbind(MLE_H0_full_outputs, MLE_H0_outputs)
}

# Take output matrix MLE_Uniform_CxCoM and sort by likelihood
# Print the sorted table for records, make sure to RENAME FILE
MLE_H0_full<-as.data.frame(MLE_H0_full_outputs)
df0 <- MLE_H0_full[order(MLE_H0_full$lnL),]
write.csv(df0,"../Summaries and Outputs/17_multi_locus_full_experiment_H0_mle_output.csv",row.names = TRUE)


# Begin algorithm for fitting CxCoM H1 to full experiment dataset
# Create matrix to store all optimized outputs of CxCoM.Likelihood.H1
# Log likelihoods calculated for each combination of y & p model parameters
MLE_H1_full_outputs <- matrix(nrow=0, ncol=39)
colnames(MLE_H1_full_outputs)<-c("ID",
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
# Output model parameters and likelihoods in matrix MLE_H1_full_outputs
for (v in 1:max_p_value) {
  
  # Define parameters starting points for y, v, and p in single run of Nelder-Mead algorithm
  # When bounding Nelder-Mead using package dfoptim, starting points cannot be at the boundaries
  parameters <- as.numeric(c(y_uniform[v,1],
                             y_uniform[v,2],
                             y_uniform[v,3],
                             y_uniform[v,4],
                             y_uniform[v,5],
                             0.5,
                             0.5))
  
  p <- p_uniform[v,]
  
  # Optimize CxCoM model parameters y and v with bounded Nelder-Mead method
  # Upper limit on y's is arbitrary, whereas upper and lower limit on v's is
  # given by definition of egg-to-adult viability (cannot exceed 0-1 range)
  MLE_H1_Full_CxCoM_Table <- nmkb(par = parameters,
                                  fn = CxCoM.Likelihood.H1,
                                  control = list(maxfeval = 50000),
                                  lower = c(0,0,0,0,0,0,0),
                                  upper = c(9,9,9,9,9,1,1))
  
  # Create output vector of y, p, likelihood, convergence code
  MLE_H1_outputs <- c(1,
                      (((MLE_H1_Full_CxCoM_Table$par[1]+MLE_H1_Full_CxCoM_Table$par[2]+MLE_H1_Full_CxCoM_Table$par[3]+MLE_H1_Full_CxCoM_Table$par[4]+MLE_H1_Full_CxCoM_Table$par[5])/(2*p[1]))*100),
                      MLE_H1_Full_CxCoM_Table$par[1]/(2*p[1]),
                      MLE_H1_Full_CxCoM_Table$par[2]/(2*p[1]),
                      MLE_H1_Full_CxCoM_Table$par[3]/(2*p[1]),
                      MLE_H1_Full_CxCoM_Table$par[4]/(2*p[1]),
                      MLE_H1_Full_CxCoM_Table$par[5]/(2*p[1]),
                      p[1],
                      MLE_H1_Full_CxCoM_Table$par[1],
                      MLE_H1_Full_CxCoM_Table$par[2],
                      MLE_H1_Full_CxCoM_Table$par[3],
                      MLE_H1_Full_CxCoM_Table$par[4],
                      MLE_H1_Full_CxCoM_Table$par[5],
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      MLE_H1_Full_CxCoM_Table$value,
                      MLE_H1_Full_CxCoM_Table$convergence)
  
  # Add vector of y, p, likelihood, code to MLE_Uniform_CxCoM
  MLE_H1_full_outputs <- rbind(MLE_H1_full_outputs, MLE_H1_outputs)
}

# Take output matrix MLE_Uniform_CxCoM and sort by likelihood
# Print the sorted table for records, make sure to RENAME FILE
MLE_H1_full<-as.data.frame(MLE_H1_full_outputs)
df1 <- MLE_H1_full[order(MLE_H1_full$lnL),]
write.csv(df1,"../Summaries and Outputs/18_multi_locus_full_experiment_H1_mle_output.csv",row.names = TRUE)

# Begin algorithm for fitting CxCoM H2 to full experiment dataset
# Create matrix to store all optimized outputs of CxCoM.Likelihood.H2
# Log likelihoods calculated for each combination of y & p model parameters
MLE_H2_full_outputs <- matrix(nrow=0, ncol=39)
colnames(MLE_H2_full_outputs)<-c("ID",
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
# Output model parameters and likelihoods in matrix MLE_H2_full_outputs
for (v in 1:max_p_value) {
  
  # Define parameters starting points for y, v, and p in single run of Nelder-Mead algorithm
  # When bounding Nelder-Mead using package dfoptim, starting points cannot be at the boundaries
  parameters <- as.numeric(c(y_uniform[v,1],
                             y_uniform[v,2],
                             y_uniform[v,3],
                             y_uniform[v,4],
                             y_uniform[v,5],
                             0.619,
                             0.615,
                             0.617,
                             0.584,
                             0.540,
                             0.591,
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
  MLE_H2_Full_CxCoM_Table <- nmkb(par = parameters,
                                  fn = CxCoM.Likelihood.H2,
                                  control = list(maxfeval = 50000),
                                  lower = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                  upper = c(9,9,9,9,9,1,1,1,1,1,1,1,1,1,1,1,1))
  
  # Create output vector of y, p, likelihood, convergence code
  MLE_H2_outputs <- c(2,
                      (((MLE_H2_Full_CxCoM_Table$par[1]+MLE_H2_Full_CxCoM_Table$par[2]+MLE_H2_Full_CxCoM_Table$par[3]+MLE_H2_Full_CxCoM_Table$par[4]+MLE_H2_Full_CxCoM_Table$par[5])/(2*p[1]))*100),
                      MLE_H2_Full_CxCoM_Table$par[1]/(2*p[1]),
                      MLE_H2_Full_CxCoM_Table$par[2]/(2*p[1]),
                      MLE_H2_Full_CxCoM_Table$par[3]/(2*p[1]),
                      MLE_H2_Full_CxCoM_Table$par[4]/(2*p[1]),
                      MLE_H2_Full_CxCoM_Table$par[5]/(2*p[1]),
                      p[1],
                      MLE_H2_Full_CxCoM_Table$par[1],
                      MLE_H2_Full_CxCoM_Table$par[2],
                      MLE_H2_Full_CxCoM_Table$par[3],
                      MLE_H2_Full_CxCoM_Table$par[4],
                      MLE_H2_Full_CxCoM_Table$par[5],
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      MLE_H2_Full_CxCoM_Table$par[6],
                      MLE_H2_Full_CxCoM_Table$par[7],
                      MLE_H2_Full_CxCoM_Table$par[8],
                      MLE_H2_Full_CxCoM_Table$par[9],
                      MLE_H2_Full_CxCoM_Table$par[10],
                      MLE_H2_Full_CxCoM_Table$par[11],
                      1,
                      1,
                      1,
                      1,
                      1,
                      1,
                      MLE_H2_Full_CxCoM_Table$par[12],
                      MLE_H2_Full_CxCoM_Table$par[13],
                      MLE_H2_Full_CxCoM_Table$par[14],
                      MLE_H2_Full_CxCoM_Table$par[15],
                      MLE_H2_Full_CxCoM_Table$par[16],
                      MLE_H2_Full_CxCoM_Table$par[17],
                      MLE_H2_Full_CxCoM_Table$value,
                      MLE_H2_Full_CxCoM_Table$convergence)
  
  # Add vector of y, p, likelihood, code to MLE_Uniform_CxCoM
  MLE_H2_full_outputs <- rbind(MLE_H2_full_outputs, MLE_H2_outputs)
}

# Take output matrix MLE_Uniform_CxCoM and sort by likelihood
# Print the sorted table for records, make sure to RENAME FILE
MLE_H2_full<-as.data.frame(MLE_H2_full_outputs)
df2 <- MLE_H2_full[order(MLE_H2_full$lnL),]
write.csv(df2,"../Summaries and Outputs/19_multi_locus_full_experiment_H2_mle_output.csv",row.names = TRUE)

# Begin algorithm for fitting CxCoM H3 to full experiment dataset
# Create matrix to store all optimized outputs of CxCoM.Likelihood.H3
# Log likelihoods calculated for each combination of y & p model parameters
MLE_H3_full_outputs <- matrix(nrow=0, ncol=39)
colnames(MLE_H3_full_outputs)<-c("ID",
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
  MLE_H3_Full_CxCoM_Table <- nmkb(par = parameters,
                                  fn = CxCoM.Likelihood.H3,
                                  control = list(maxfeval = 50000),
                                  lower = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                                  upper = c(9,9,9,9,9,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
  
  # Create output vector of y, p, likelihood, convergence code
  MLE_H3_outputs <- c(3,
                      (((MLE_H3_Full_CxCoM_Table$par[1]+MLE_H3_Full_CxCoM_Table$par[2]+MLE_H3_Full_CxCoM_Table$par[3]+MLE_H3_Full_CxCoM_Table$par[4]+MLE_H3_Full_CxCoM_Table$par[5])/(2*p[1]))*100),
                      MLE_H3_Full_CxCoM_Table$par[1]/(2*p[1]),
                      MLE_H3_Full_CxCoM_Table$par[2]/(2*p[1]),
                      MLE_H3_Full_CxCoM_Table$par[3]/(2*p[1]),
                      MLE_H3_Full_CxCoM_Table$par[4]/(2*p[1]),
                      MLE_H3_Full_CxCoM_Table$par[5]/(2*p[1]),
                      p[1],
                      MLE_H3_Full_CxCoM_Table$par[1],
                      MLE_H3_Full_CxCoM_Table$par[2],
                      MLE_H3_Full_CxCoM_Table$par[3],
                      MLE_H3_Full_CxCoM_Table$par[4],
                      MLE_H3_Full_CxCoM_Table$par[5],
                      MLE_H3_Full_CxCoM_Table$par[6],
                      MLE_H3_Full_CxCoM_Table$par[7],
                      MLE_H3_Full_CxCoM_Table$par[8],
                      MLE_H3_Full_CxCoM_Table$par[9],
                      MLE_H3_Full_CxCoM_Table$par[10],
                      MLE_H3_Full_CxCoM_Table$par[11],
                      MLE_H3_Full_CxCoM_Table$par[12],
                      MLE_H3_Full_CxCoM_Table$par[13],
                      MLE_H3_Full_CxCoM_Table$par[14],
                      MLE_H3_Full_CxCoM_Table$par[15],
                      MLE_H3_Full_CxCoM_Table$par[16],
                      MLE_H3_Full_CxCoM_Table$par[17],
                      MLE_H3_Full_CxCoM_Table$par[18],
                      MLE_H3_Full_CxCoM_Table$par[19],
                      MLE_H3_Full_CxCoM_Table$par[20],
                      MLE_H3_Full_CxCoM_Table$par[21],
                      MLE_H3_Full_CxCoM_Table$par[22],
                      MLE_H3_Full_CxCoM_Table$par[23],
                      MLE_H3_Full_CxCoM_Table$par[24],
                      MLE_H3_Full_CxCoM_Table$par[25],
                      MLE_H3_Full_CxCoM_Table$par[26],
                      MLE_H3_Full_CxCoM_Table$par[27],
                      MLE_H3_Full_CxCoM_Table$par[28],
                      MLE_H3_Full_CxCoM_Table$par[29],
                      MLE_H3_Full_CxCoM_Table$value,
                      MLE_H3_Full_CxCoM_Table$convergence)
  
  # Add vector of y, p, likelihood, code to MLE_Uniform_CxCoM
  MLE_H3_full_outputs <- rbind(MLE_H3_full_outputs, MLE_H3_outputs)
}
