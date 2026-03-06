
library(MASS)
library(matrixcalc)
library(mbend)
library(Matrix)
library(Rmosek)
library(ggplot2)
library(lme4)
library(bartCause)
library(grf)
library(parallel)
library(dplyr)
library(future)
library(future.apply)

source('functions_sim1.r')

###########################################################
#                                                         
#       Case 1
#
###########################################################

run_iteration_1 <- function(iter){
  
  # 2026.03.04 Updated
  MSE_CATE <- function(tau_hat, tau, cluster_id){
    # pooled: MSE within each cluster, then average across clusters
    clusters <- split(seq_along(tau_hat), cluster_id)
    cluster_mses <- sapply(clusters, function(idx) mean((tau_hat[idx] - tau[idx])^2))
    return(mean(cluster_mses))
  }
  
  pooled_value <- function(policy, Y1, Y0, cluster_id){
    # V(d) using P_n: J^{-1} sum_j n_j^{-1} sum_i [d_ij * Y1_ij + (1-d_ij) * Y0_ij]
    clusters <- split(seq_along(policy), cluster_id)
    cluster_values <- sapply(clusters, function(idx) {
      mean(policy[idx] * Y1[idx] + (1 - policy[idx]) * Y0[idx])
    })
    return(mean(cluster_values))
  }
  
  # Generate data
  df <- create_multileveldata_D1(cluster_num = 300, # number of clusters
                                    cluster_size_min = 10,
                                    cluster_size_max = 50,
                                    R_var = 0.6653, # 0.6653 variance of residual
                                    V_var = 0.8,  # 0.8 variance of cluster effect in selection
                                    U_var = 1.2, # 1.2 variance of cluster effect in outcome
                                    tau_var = 0, # 0 variance of cluster effect in treatment
                                    U_var_coef = 1,
                                    clustereffect=TRUE)
  
  # True beta
  # to be consistent with the design 1, the coefficient of S2 changed from -0.2 to 0.8
  beta_true <- c(0.8, 0.3, 0.2, 0.2, -0.5, 
                 0.5, 0.4, 0.2, 0.1, 0.8, -0.2)
  
  
  # Prepare data
  dat0 <- df[, c("id", "X11", "X12", "X13", "S1", "X14", "X21", "X22", "X23", "S2", "S2_2", "A", "Y")]
  dat0$id <- as.character(dat0$id)
  tau <- df$tau
  ate <- mean(tau)
  
  # 2026.03.04 Updated:
  C_id <- dat0$id
  
  # true policy
  policy_true <- df$policy
  
  # true counterfactuals
  Y1 <- df$Y1
  Y0 <- df$Y0
  
  # true value
  #value_true <- mean(ifelse(policy_true > 0, Y1, Y0))
  value_true <- pooled_value(policy_true, Y1, Y0, C_id)
  
  # generate random policy
  A_random <- rbinom(nrow(dat0), 1, 0.42)
  # value_random <- mean(ifelse(A_random >0, Y1, Y0))
  value_random <- pooled_value(A_random, Y1, Y0, C_id)
  
  # get the relative utility(RU) for true
  RU_true <- (value_true - value_random)/value_random
  
  # Set constraints
  delta_set <- c(20, 20)
  
  # GLMM models fitting
  glmer_Control <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
  glmm_out <- GLMM_model(data = dat0,
                         outcome = "Y",
                         treatment = "A",
                         cluster = "id",
                         n_AGQ = 2,
                         fixed_intercept = TRUE,
                         random_slope = FALSE,
                         glmer_Control = glmer_Control)
  
  
  # Setting 1: Original cluster, old_version = FALSE
  ml_fr_out_setting1 <- fairCATE_multilevel(data = dat0,
                                            sensitive = c("S1", "S2", "S2_2"),
                                            legitimate = NULL,
                                            fairness = c("tau~S1", "tau~S2"),
                                            treatment = "A",
                                            outcome = "Y",
                                            cluster = "id",
                                            multicategorical = NULL,
                                            outcome.LMM = glmm_out$outcome.LMM,
                                            ps.GLMM = glmm_out$ps.GLMM,
                                            delta = delta_set,
                                            ps.trim = "Sturmer.1")
  
  beta_hat_setting1 <- ml_fr_out_setting1$beta.hat
  tau_hat_setting1 <- ml_fr_out_setting1$tau_hat
  
  MSE_beta_setting1 <- mean((beta_hat_setting1 - beta_true)^2)
  #MSE_tau_setting1 <- mean((tau_hat_setting1 - tau)^2)
  MSE_tau_setting1 <- MSE_CATE(tau_hat_setting1, tau, C_id)
  
  # get the estimated policy
  policy_est_setting1 <- ml_fr_out_setting1$policy_hat
  
  # get the value of setting 1
  # value_setting1 <- mean(ifelse(policy_est_setting1 > 0, Y1, Y0))
  value_setting1 <- pooled_value(policy_est_setting1, Y1, Y0, C_id)
  
  # get the relative utility(RU) for setting 1
  RU_setting1 <- (value_setting1 - value_random)/value_random
  
  # get the accuracy of estimated policy comparing to true policy
  policy_accuracy_setting1 <- mean(policy_est_setting1 == policy_true)
  
  # Setting 2: No cluster, old_version = FALSE
  dat0_no_cluster <- dat0 %>% mutate(id = "1")
  ml_fr_out_setting2 <- fairCATE_multilevel(data = dat0_no_cluster,
                                            sensitive = c("S1", "S2", "S2_2"),
                                            legitimate = NULL,
                                            fairness = c("tau~S1", "tau~S2"),
                                            treatment = "A",
                                            outcome = "Y",
                                            cluster = "id",
                                            multicategorical = NULL,
                                            outcome.LMM = glmm_out$outcome.LMM,
                                            ps.GLMM = glmm_out$ps.GLMM,
                                            delta = delta_set,
                                            ps.trim = "Sturmer.1")
  
  beta_hat_setting2 <- ml_fr_out_setting2$beta.hat
  tau_hat_setting2 <- ml_fr_out_setting2$tau_hat
  
  MSE_beta_setting2 <- mean((beta_hat_setting2 - beta_true)^2)
  # MSE_tau_setting2 <- mean((tau_hat_setting2 - tau)^2)
  MSE_tau_setting2 <- MSE_CATE(tau_hat_setting2, tau, C_id)
  
  # get the estimated policy
  policy_est_setting2 <- ml_fr_out_setting2$policy_hat
  
  # get the value of setting 2
  #value_setting2 <- mean(ifelse(policy_est_setting2 > 0, Y1, Y0))
  value_setting2 <- pooled_value(policy_est_setting2, Y1, Y0, C_id)
  
  # get the relative utility(RU) for setting 2
  RU_setting2 <- (value_setting2 - value_random)/value_random
  
  # get the accuracy of estimated policy comparing to true policy
  policy_accuracy_setting2 <- mean(policy_est_setting2 == policy_true)
  
  # Setting 3: Totally cluster-ignored
  glm_out <- GLM_model(data = dat0_no_cluster,
                       outcome = "Y",
                       treatment = "A",
                       cluster = "id")
  
  outcome.LM <- glm_out$outcome.LM
  ps.GLM <- glm_out$ps.GLM
  
  ml_fr_out <- fairCATE_multilevel_cluster_ignored(data = dat0_no_cluster,
                                                   sensitive = c("S1", "S2"),
                                                   legitimate = NULL,
                                                   fairness = c("tau~S1", "tau~S2"),
                                                   treatment = "A",
                                                   outcome = "Y",
                                                   cluster = "id",
                                                   multicategorical = NULL,
                                                   outcome.LMM = outcome.LM,
                                                   ps.GLMM = ps.GLM,
                                                   delta = c(20, 20),
                                                   ps.trim = "Sturmer.1")
  
  beta_hat_setting3 <- ml_fr_out$beta.hat
  tau_hat_setting3 <- ml_fr_out$tau_hat
  MSE_beta_setting3 <- mean((beta_hat_setting3 - beta_true)^2)
  # MSE_tau_setting3 <- mean((tau_hat_setting3 - tau)^2)
  MSE_tau_setting3 <- MSE_CATE(tau_hat_setting3, tau, C_id)
  
  # get the estimated policy
  policy_est_setting3 <- ml_fr_out$policy_hat
  
  # get the value of setting 3
  #value_setting3 <- mean(ifelse(policy_est_setting3 > 0, Y1, Y0))
  value_setting3 <- pooled_value(policy_est_setting3, Y1, Y0, C_id)
  
  # get the relative utility(RU) for setting 3
  RU_setting3 <- (value_setting3 - value_random)/value_random
  
  # get the accuracy of estimated policy comparing to true policy
  policy_accuracy_setting3 <- mean(policy_est_setting3 == policy_true)
  
  return(list(MSE_beta_setting1 = MSE_beta_setting1,
              MSE_tau_setting1 = MSE_tau_setting1,
              policy_accuracy_setting1 = policy_accuracy_setting1,
              RU_setting1 = RU_setting1,
              MSE_beta_setting2 = MSE_beta_setting2,
              MSE_tau_setting2 = MSE_tau_setting2,
              policy_accuracy_setting2 = policy_accuracy_setting2,
              RU_setting2 = RU_setting2,
              MSE_beta_setting3 = MSE_beta_setting3,
              MSE_tau_setting3 = MSE_tau_setting3,
              policy_accuracy_setting3 = policy_accuracy_setting3,
              RU_setting3 = RU_setting3,
              RU_true = RU_true))
}


# Set up parallel processing
plan(multisession, workers = 8)

total_time0 <- Sys.time()
# Run the simulation for each combination of U_var and tau_var
results <- future_lapply(1:500, run_iteration_1)


# Convert list of lists into a data frame
results_df1 <- do.call(rbind, lapply(results, function(res) {
  data.frame(
    MSE_beta_setting1 = res$MSE_beta_setting1,
    MSE_beta_setting2 = res$MSE_beta_setting2,
    MSE_beta_setting3 = res$MSE_beta_setting3,
    
    MSE_tau_setting1 = res$MSE_tau_setting1,
    MSE_tau_setting2 = res$MSE_tau_setting2,
    MSE_tau_setting3 = res$MSE_tau_setting3,
    
    policy_accuracy_setting1 = res$policy_accuracy_setting1,
    policy_accuracy_setting2 = res$policy_accuracy_setting2,
    policy_accuracy_setting3 = res$policy_accuracy_setting3,
    
    RU_true = res$RU_true,
    RU_setting1 = res$RU_setting1,
    RU_setting2 = res$RU_setting2,
    RU_setting3 = res$RU_setting3
  )
}))

# save the results_df3 as a rdata file
save(results_df1, file = "results/results_df_case_1.rdata")

# get the averaged results
output_mean <- colMeans(results_df1, na.rm = TRUE)

