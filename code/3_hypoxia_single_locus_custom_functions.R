# custom script for goodness-of-fit G-test for marker viability effects under six models of viability
# H0: no excess experimental mortality in multiply-marked crosses above and beyond the marker-free cross
# H1: excess experimental mortality in multiply-marked crosses is random with respect to marker alleles
# H2: excess experimental mortality in multiply-marked crosses is solely due to mutant alleles at marker loci
# H3: excess experimental mortality in multiply-marked crosses is due to viability effects associated with both mutant and wildtype alleles at marker loci
# unattributable experimental mortality (v_x) estimated as 0.921 from external marker free experimental evidence

single_locus_analysis <- function(){# enter egg count sample size (n)
  n <- 10952
  
  # enter observed and unobserved counts
  observed_female_wildtype_adult_count <- single_locus_input[1]
  observed_male_wildtype_adult_count   <- single_locus_input[2]
  observed_female_mutant_adult_count   <- single_locus_input[3]
  observed_male_mutant_adult_count     <- single_locus_input[4]
  observed_failed_development_count    <- single_locus_input[5]
  
  # solve linear system of equations to estimate viability effects
  
  # H0: no excess experimental mortality (no coefficients to be estimated)
  
  # H1: excess experimental mortality in multiply-marked crosses is random with respect to marker alleles (1 parameter to be estimated)
  # note that 2521.698 is (1/4)*n*v_x = (1/4)*10952*0.921 and 865.208 is n*(1-v_x) = 10952*.079
  H1_a <- c(-2521.698, -2521.698, -2521.698, -2521.698, 10086.79)
  H1_y <- c(observed_female_wildtype_adult_count-2521.698,
            observed_male_wildtype_adult_count-2521.698,
            observed_female_mutant_adult_count-2521.698,
            observed_male_mutant_adult_count-2521.698,
            observed_failed_development_count-865.208)
  H1_matrix <- as.data.frame(cbind(H1_y, H1_a))
  H1_model <- lm(H1_y ~ H1_a + 0, data=H1_matrix)
  H1_vector <- as.vector(coef(H1_model))
  
  # H2: excess experimental mortality in multiply-marked crosses is solely due to mutant alleles at marker loci (2 parameters to be estimated)
  H2_a <- c(0, 0, -2521.698, 0, 2521.698)
  H2_b <- c(0, 0, 0, -2521.698, 2521.698)
  H2_y <- c(observed_female_wildtype_adult_count-2521.698,
            observed_male_wildtype_adult_count-2521.698,
            observed_female_mutant_adult_count-2521.698,
            observed_male_mutant_adult_count-2521.698,
            observed_failed_development_count-865.208)
  H2_matrix <- as.data.frame(cbind(H2_y, H2_a, H2_b))
  H2_model <- lm(H2_y ~ H2_a + H2_b + 0, data=H2_matrix)
  H2_vector <- as.vector(coef(H2_model))
  
  # H3: excess experimental mortality in multiply-marked crosses is due to viability effects associated with both mutant and wildtype alleles at marker loci
  H3_a <- c(-2521.698, 0, 0, 0, 2521.698)
  H3_b <- c(0, -2521.698, 0, 0, 2521.698)
  H3_c <- c(0, 0, -2521.698, 0, 2521.698)
  H3_d <- c(0, 0, 0, -2521.698, 2521.698)
  H3_y <- c(observed_female_wildtype_adult_count-2521.698,
            observed_male_wildtype_adult_count-2521.698,
            observed_female_mutant_adult_count-2521.698,
            observed_male_mutant_adult_count-2521.698,
            observed_failed_development_count-865.208)
  H3_matrix <- as.data.frame(cbind(H3_y, H3_a, H3_b, H3_c, H3_d))
  H3_model <- lm(H3_y ~ H3_a + H3_b + H3_c + H3_d + 0, data=H3_matrix)
  H3_vector <- as.vector(coef(H3_model))
  
  # save outputs as viability effects for goodness-of-fit tests
  H0_v_x <- 0.921
  H1_v_x <- 0.921
  H1_v_1 <- 1-H1_vector[1]
  H2_v_x <- 0.921
  H2_v_2 <- 1-H2_vector[1]
  H2_v_3 <- 1-H2_vector[2]
  H3_v_x <- 0.921
  H3_v_4 <- 1-H3_vector[1]
  H3_v_5 <- 1-H3_vector[2]
  H3_v_6 <- 1-H3_vector[3]
  H3_v_7 <- 1-H3_vector[4]
  
  # expected counts under specified hypotheses/viability models
  # H0: no experimental mortality
  H0_expected_female_wildtype_count <- (n*H0_v_x)/4
  H0_expected_male_wildtype_count   <- (n*H0_v_x)/4
  H0_expected_female_mutant_count   <- (n*H0_v_x)/4
  H0_expected_male_mutant_count     <- (n*H0_v_x)/4
  H0_expected_failed_development    <- n*(1-H0_v_x)
  
  # H1: experimental mortality is random with respect to marker alleles
  H1_expected_female_wildtype_count <- (n*H1_v_x*H1_v_1)/4
  H1_expected_male_wildtype_count   <- (n*H1_v_x*H1_v_1)/4
  H1_expected_female_mutant_count   <- (n*H1_v_x*H1_v_1)/4
  H1_expected_male_mutant_count     <- (n*H1_v_x*H1_v_1)/4
  H1_expected_failed_development    <- n*(1-(H1_v_x*H1_v_1))
  
  # H2: experimental mortality has both a random component and marker associated effects
  H2_expected_female_wildtype_count <- (n*H2_v_x)/4
  H2_expected_male_wildtype_count   <- (n*H2_v_x)/4
  H2_expected_female_mutant_count   <- (n*H2_v_x*H2_v_2)/4
  H2_expected_male_mutant_count     <- (n*H2_v_x*H2_v_3)/4
  H2_expected_failed_development    <- (n*(4-(H2_v_x*(2+H2_v_2+H2_v_3))))/4
  
  # H3: experimental mortality has both a random component and sex-specific marker associated effects
  H3_expected_female_wildtype_count <- (n*H3_v_x*H3_v_4)/4
  H3_expected_male_wildtype_count   <- (n*H3_v_x*H3_v_5)/4
  H3_expected_female_mutant_count   <- (n*H3_v_x*H3_v_6)/4
  H3_expected_male_mutant_count     <- (n*H3_v_x*H3_v_7)/4
  H3_expected_failed_development    <- (n*(4-(H3_v_x*(H3_v_4+H3_v_5+H3_v_6+H3_v_7))))/4
  
  # Goodness-of-fit G-test for each hypothesis/viability model using parameters fit by method of least squares with 0 intercept
  G_statistic_H0 <- 2*sum(observed_female_wildtype_adult_count*(log((observed_female_wildtype_adult_count/H0_expected_female_wildtype_count))),
                          observed_male_wildtype_adult_count*(log((observed_male_wildtype_adult_count/H0_expected_male_wildtype_count))),
                          observed_female_mutant_adult_count*(log((observed_female_mutant_adult_count/H0_expected_female_mutant_count))),
                          observed_male_mutant_adult_count*(log((observed_male_mutant_adult_count/H0_expected_male_mutant_count))),
                          observed_failed_development_count*(log((observed_failed_development_count/H0_expected_failed_development))))
  
  G_statistic_H1 <- 2*sum(observed_female_wildtype_adult_count*(log((observed_female_wildtype_adult_count/H1_expected_female_wildtype_count))),
                          observed_male_wildtype_adult_count*(log((observed_male_wildtype_adult_count/H1_expected_male_wildtype_count))),
                          observed_female_mutant_adult_count*(log((observed_female_mutant_adult_count/H1_expected_female_mutant_count))),
                          observed_male_mutant_adult_count*(log((observed_male_mutant_adult_count/H1_expected_male_mutant_count))),
                          observed_failed_development_count*(log((observed_failed_development_count/H1_expected_failed_development))))
  
  G_statistic_H2 <- 2*sum(observed_female_wildtype_adult_count*(log((observed_female_wildtype_adult_count/H2_expected_female_wildtype_count))),
                          observed_male_wildtype_adult_count*(log((observed_male_wildtype_adult_count/H2_expected_male_wildtype_count))),
                          observed_female_mutant_adult_count*(log((observed_female_mutant_adult_count/H2_expected_female_mutant_count))),
                          observed_male_mutant_adult_count*(log((observed_male_mutant_adult_count/H2_expected_male_mutant_count))),
                          observed_failed_development_count*(log((observed_failed_development_count/H2_expected_failed_development))))
  
  G_statistic_H3 <- 2*sum(observed_female_wildtype_adult_count*(log((observed_female_wildtype_adult_count/H3_expected_female_wildtype_count))),
                          observed_male_wildtype_adult_count*(log((observed_male_wildtype_adult_count/H3_expected_male_wildtype_count))),
                          observed_female_mutant_adult_count*(log((observed_female_mutant_adult_count/H3_expected_female_mutant_count))),
                          observed_male_mutant_adult_count*(log((observed_male_mutant_adult_count/H3_expected_male_mutant_count))),
                          observed_failed_development_count*(log((observed_failed_development_count/H3_expected_failed_development))))
  
  # Calculate p-value for G-tests
  # Note H3 has 0 degrees of freedom to conduct this test
  H0_p_value <- pchisq(G_statistic_H0, 4, lower.tail=FALSE)
  H1_p_value <- pchisq(G_statistic_H1, 3, lower.tail=FALSE)
  H2_p_value <- pchisq(G_statistic_H2, 2, lower.tail=FALSE)
  
  # Build a dataframe to export as supplementary tables
  female_wildtype_viability <- c(1, H1_v_1, 1, H3_v_4)
  male_wildtype_viability <- c(1, H1_v_1, 1, H3_v_5)
  female_mutant_viability <- c(1, H1_v_1, H2_v_2, H3_v_6)
  male_mutant_viability <- c(1, H1_v_1, H2_v_3, H3_v_7)
  G_statistic <- c(G_statistic_H0,G_statistic_H1,G_statistic_H2,NA)
  degrees_of_freedom <- c(4,3,2,NA)
  p_value <- c(H0_p_value,H1_p_value,H2_p_value,NA)
  
  single_locus_output <- data.frame(female_wildtype_viability=female_wildtype_viability,
                                    male_wildtype_viability=male_wildtype_viability,
                                    female_mutant_viability=female_mutant_viability,
                                    male_mutant_viability=male_mutant_viability,
                                    G_statistic=G_statistic,
                                    degrees_of_freedom=degrees_of_freedom,
                                    p_value=p_value,
                                    row.names =c("H0:", "H1:", "H2:", "H3:"))
  return(single_locus_output)}