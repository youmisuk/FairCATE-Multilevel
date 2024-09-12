# date: 2024.08.12
# purpose:
#     - this file is based on the Design_2_parallel_V2.r version
#     - this file add the propensity score in to the causal forest model.

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

source("data_generating_process_multilevel.r")
source("fairCATE_multilevel_functions.r")


run_iteration <- function(iter) {
  
  ###################################################################
  #                                                                 #
  #                     1.0 In each Repetition                      #
  #                                                                 #
  ###################################################################
  
  library(dplyr)
  MSE_CATE <- function(tau_hat, tau){
    return(mean((tau_hat - tau)^2))
  }

  # 
  # ------------------ #
  # 1.1 Generate Data  #
  # ------------------ #
  
  # 1.1.1, generate simulated data
  
  df <- create_multileveldata_D2(cluster_num = 300, 
                                    cluster_size = 25, 
                                    E_var = 0.6653, # var of residual
                                    R_var = 1.95, # var of cluster effect in selection
                                    U_var = 0.0776, # var of cluster effect in outcome
                                    clustereffect=TRUE)
  
  
  ATE_true <- mean(df$tau)
  value_true <- mean(ifelse(df$policy > 0, df$Y1, df$Y0))
  dat0 <- df[,c("id","X11","X12","X13","S1","L","X21","X22","X23","S2","A","Y","S_is_1","S_is_2","S_is_3")]
  
  # 1.1.2, generate a random treatment A in the proportion of 0.42, similar to df$A
  dat0$id <- as.character(dat0$id)
  
  # PAN: 2024.07.20
  # mutate the dat0 with intersectional sensitive variables
  dat0_is <- dat0 %>% select(-S1, -S2)
  dat0 <- dat0 %>% select(-S_is_1, -S_is_2, -S_is_3)
  
  # create the reference indicator but not to be included into the dat0_is
  S_is_0 <- ifelse(dat0_is$S_is_1 == 0 & dat0_is$S_is_2 == 0 & dat0_is$S_is_3 == 0, 1, 0)
  
  # extract all the S related columns to save ram and runnning time
  S1 <- dat0$S1
  S2 <- dat0$S2
  S_is_1 <- dat0_is$S_is_1
  S_is_2 <- dat0_is$S_is_2
  S_is_3 <- dat0_is$S_is_3
  policy <- df$policy
  Y1 <- df$Y1
  Y0 <- df$Y0
  tau <- df$tau
  
  A_random <- rbinom(nrow(dat0), 1, 0.42)
  value_random <- mean(ifelse(A_random > 0, Y1, Y0))
  
  # get the relative utility
  RU_true <- (value_true - value_random ) / value_random
  
  
  unfairness_1 <- abs(mean(policy[S_is_1==1]) - mean(policy[S_is_1==0]))
  
  unfairness_2 <- abs(mean(policy[S_is_2==1]) - mean(policy[S_is_2==0]))
  
  unfairness_3 <- abs(mean(policy[S_is_3==1]) - mean(policy[S_is_3==0]))
  
  # get the average unfairness
  average_unfairness_true <- (unfairness_1 + unfairness_2 + unfairness_3) / 3
  
  # ------------------ #
  #   1.2 Run BART     #
  # ------------------ #
  
  # 1.2.1 fit BART causal model for CATE
  # out.BART <- bartc(response = dat0$Y, treatment = dat0$A, 
  #                   confounders = dat0[, c("X11","X12","X13","S1","L","X21","X22","X23","S2")], 
  #                   method.rsp = "bart", method.trt = "bart", 
  #                   p.scoreAsCovariate = TRUE, 
  #                   use.rbrt = FALSE, 
  #                   keepTrees = FALSE)
  
  out.BART <- bartc(response = dat0_is$Y, treatment = dat0_is$A, 
                    confounders = dat0_is[, c("X11","X12","X13","S_is_1","L","X21","X22","X23","S_is_2","S_is_3")], 
                    method.rsp = "bart", method.trt = "bart", 
                    p.scoreAsCovariate = TRUE, 
                    use.rbrt = FALSE, 
                    keepTrees = FALSE)
  
  cate.BART <- fitted(out.BART, type = "icate", sample = "all") 
  
  # 1.2.2 retrieve the OTR, value, and ATE
  OTR_BART <- as.numeric(cate.BART > 0)
  value_BART <- mean(ifelse(OTR_BART > 0, Y1, Y0))
  ATE_BART <- mean(cate.BART)
  
  # 1.2.3 get the Relative Utility
  RU_BART <- (value_BART - value_random ) / value_random
  
  # 1.2.4 get the Average Unfairness
  
  unfairness_1 <- abs(mean(OTR_BART[S_is_1==1]) - mean(OTR_BART[S_is_1==0]))
  
  unfairness_2 <- abs(mean(OTR_BART[S_is_2==1]) - mean(OTR_BART[S_is_2==0]))
  
  unfairness_3 <- abs(mean(OTR_BART[S_is_3==1]) - mean(OTR_BART[S_is_3==0]))
  
  # get the average unfairness
  average_unfairness_bart <- (unfairness_1 + unfairness_2 + unfairness_3) / 3
  
  # PAN: 2024.07.29
  # get the CATE_MSE
  MSE_CATE_BART <- MSE_CATE(cate.BART, tau)
  # get the CATE unfairness
  unfairness_1 <- abs(mean(cate.BART[S_is_1==1]) - mean(cate.BART[S_is_1==0]))
  unfairness_2 <- abs(mean(cate.BART[S_is_2==1]) - mean(cate.BART[S_is_2==0]))
  unfairness_3 <- abs(mean(cate.BART[S_is_3==1]) - mean(cate.BART[S_is_3==0]))
  
  average_unfairness_bart_CATE <- (unfairness_1 + unfairness_2 + unfairness_3) / 3
  
  # Remove the variables
  rm(out.BART, cate.BART)
  
  # Run garbage collection
  gc()
  
  # ------------------ #
  #   1.3 Run CausalF  #
  # ------------------ #
  
  # 1.3.1 implement Causal Forests with 
  # X = confounders; Y = outcome; W = treatment
  # out.cf <- causal_forest(X = dat0[, c("X11","X12","X13","S1","L","X21","X22","X23","S2")], 
  #                         Y = dat0$Y, 
  #                         W = dat0$A) 
  
  out.cf <- causal_forest(X = dat0_is[, c("X11","X12","X13","S_is_1","L","X21","X22","X23","S_is_2","S_is_3")], 
                          Y = dat0_is$Y, 
                          W = dat0_is$A) 
  
  # get individual treatment effect estimates, \tau_ij
  cate.cf <- predict(out.cf, type="vector", estimate.variance = FALSE) 
  
  # 1.3.2 retrieve the estimated OTR, value, and ATE
  OTR_cf <- as.numeric(cate.cf$predictions > 0)
  value_cf <- mean(ifelse(OTR_cf > 0, Y1, Y0))
  cate_cf_pre <- cate.cf$predictions
  
  # 1.3.3 get the Relative Utility
  RU_cf <- (value_cf - value_random ) / value_random
  
  unfairness_1 <- abs(mean(OTR_cf[S_is_1==1]) - mean(OTR_cf[S_is_1==0]))
  
  unfairness_2 <- abs(mean(OTR_cf[S_is_2==1]) - mean(OTR_cf[S_is_2==0]))
  
  unfairness_3 <- abs(mean(OTR_cf[S_is_3==1]) - mean(OTR_cf[S_is_3==0]))
  
  # get the average unfairness
  average_unfairness_cf <- (unfairness_1 + unfairness_2 + unfairness_3) / 3
  
  # PAN: 2024.07.29
  # get the CATE_MSE
  MSE_CATE_cf <- MSE_CATE(cate_cf_pre, tau)
  # get the CATE unfairness
  unfairness_1 <- abs(mean(cate_cf_pre[S_is_1==1]) - mean(cate_cf_pre[S_is_1==0]))
  unfairness_2 <- abs(mean(cate_cf_pre[S_is_2==1]) - mean(cate_cf_pre[S_is_2==0]))
  unfairness_3 <- abs(mean(cate_cf_pre[S_is_3==1]) - mean(cate_cf_pre[S_is_3==0]))
  
  average_unfairness_cf_CATE <- (unfairness_1 + unfairness_2 + unfairness_3) / 3
  
  # Remove the variables
  rm(out.cf, cate.cf)
  
  # 1.4.0.0 model fitting
  glmer_Control <- glmerControl(optimizer = "bobyqa",
                                optCtrl = list(maxfun=100000))
  
  # model fitting
  glmm_out <- GLMM_model(data = dat0,
                         outcome = "Y",
                         treatment = "A",
                         cluster = "id",
                         n_AGQ = 2,
                         fixed_intercept = FALSE,
                         glmer_Control = glmer_Control)
  
  glmm_out_is <- GLMM_model(data = dat0_is,
                            outcome = "Y",
                            treatment = "A",
                            cluster = "id",
                            n_AGQ = 2,
                            fixed_intercept = FALSE,
                            glmer_Control = glmer_Control)
  
  # 1.4.0.1 get the unconstrained results to calculate the UG and UD
  ml_fr_out <- fairCATE_multilevel(data = dat0,
                                   sensitive = c("S1","S2"),
                                   legitimate = NULL,
                                   fairness = c("tau~S1","tau~S2"),
                                   treatment = "A",
                                   outcome = "Y",
                                   cluster = "id",
                                   multicategorical = NULL,
                                   outcome.LMM = glmm_out$outcome.LMM,
                                   ps.GLMM = glmm_out$ps.GLMM,
                                   fixed_intercept = FALSE,
                                   delta = c(20,20),
                                   ps.trim="Sturmer.1")
  
  # 1.4.0.4 get the unconstrained results to calculate the UG and UD
  ml_fr_out_is <- fairCATE_multilevel(data = dat0_is,
                                      sensitive = c("S_is_1","S_is_2","S_is_3"),
                                      legitimate = NULL,
                                      fairness = c("tau~S_is_1","tau~S_is_2","tau~S_is_3"),
                                      treatment = "A",
                                      outcome = "Y",
                                      cluster = "id",
                                      multicategorical = NULL,
                                      outcome.LMM = glmm_out_is$outcome.LMM,
                                      ps.GLMM = glmm_out_is$ps.GLMM,
                                      fixed_intercept = FALSE,
                                      delta = c(20, 20, 20),
                                      ps.trim="Sturmer.1")
  
  # PAN: 2024.08.12
  out.cf <- causal_forest(X = dat0_is[, c("X11","X12","X13","S_is_1","L",
                                          "X21","X22","X23","S_is_2","S_is_3")], 
                          Y = dat0_is$Y, 
                          W = dat0_is$A,
                          W.hat = ml_fr_out_is$ps_scores) 
  
  # get individual treatment effect estimates, \tau_ij
  cate.cf <- predict(out.cf, type="vector", estimate.variance = FALSE) 
  
  # 1.3.2 retrieve the estimated OTR, value, and ATE
  OTR_cf <- as.numeric(cate.cf$predictions > 0)
  value_cf <- mean(ifelse(OTR_cf > 0, Y1, Y0))
  cate_cf_pre <- cate.cf$predictions
  
  # 1.3.3 get the Relative Utility
  RU_cf_ps <- (value_cf - value_random ) / value_random
  
  unfairness_1 <- abs(mean(OTR_cf[S_is_1==1]) - mean(OTR_cf[S_is_1==0]))
  
  unfairness_2 <- abs(mean(OTR_cf[S_is_2==1]) - mean(OTR_cf[S_is_2==0]))
  
  unfairness_3 <- abs(mean(OTR_cf[S_is_3==1]) - mean(OTR_cf[S_is_3==0]))
  
  # get the average unfairness
  average_unfairness_cf_ps <- (unfairness_1 + unfairness_2 + unfairness_3) / 3
  
  # PAN: 2024.07.29
  # get the CATE_MSE
  MSE_CATE_cf_ps <- MSE_CATE(cate_cf_pre, tau)
  # get the CATE unfairness
  unfairness_1 <- abs(mean(cate_cf_pre[S_is_1==1]) - mean(cate_cf_pre[S_is_1==0]))
  unfairness_2 <- abs(mean(cate_cf_pre[S_is_2==1]) - mean(cate_cf_pre[S_is_2==0]))
  unfairness_3 <- abs(mean(cate_cf_pre[S_is_3==1]) - mean(cate_cf_pre[S_is_3==0]))
  
  average_unfairness_cf_CATE_ps <- (unfairness_1 + unfairness_2 + unfairness_3) / 3
  
  # Remove the variables
  rm(out.cf, cate.cf)
  
  
  # --------------------- #
  #   1.4 Run ML_FRCATE   #
  # --------------------- #
  
  
  tau_hat_init <- ml_fr_out$tau_hat
  # get the OTR, value
  OTR_ml_fr_init <- as.numeric(tau_hat_init > 0)
  value_ml_fr_init <- mean(ifelse(OTR_ml_fr_init > 0, Y1, Y0))
  
  # condition1: E[d=1|S_is_1=1] - E[d=1|S_is_1=0]
  unfairness_1_init <- abs(mean(OTR_ml_fr_init[S_is_1==1]) - mean(OTR_ml_fr_init[S_is_1==0]))
  
  # condition2: E[d=1|S_is_2=1] - E[d=1|S_is_0=1]
  unfairness_2_init <- abs(mean(OTR_ml_fr_init[S_is_2==1]) - mean(OTR_ml_fr_init[S_is_2==0]))
  
  # condition3: E[d=1|S_is_3=1] - E[d=1|S_is_0=1]
  unfairness_3_init <- abs(mean(OTR_ml_fr_init[S_is_3==1]) - mean(OTR_ml_fr_init[S_is_3==0]))
  
  # get the average unfairness
  average_unfairness_ml_fr_init_1 <- (unfairness_1_init + unfairness_2_init + unfairness_3_init) / 3
  
  # 1.4.0.3 model fitting for the intersectional data
  

  tau_hat_init_is <- ml_fr_out_is$tau_hat
  # get the OTR, value
  OTR_ml_fr_init_is <- as.numeric(tau_hat_init_is > 0)
  value_ml_fr_init_is <- mean(ifelse(OTR_ml_fr_init_is > 0, Y1, Y0))
  
  # PAN: 2024.07.21
  # Revised the unfairness calculation for intersectional fairness
  # condition1: E[d=1|S_is_1=1] - E[d=1|S_is_1=0]
  unfairness_1_init_is <- abs(mean(OTR_ml_fr_init_is[S_is_1==1]) - mean(OTR_ml_fr_init_is[S_is_1==0]))
  
  # condition2: E[d=1|S_is_2=1] - E[d=1|S_is_0=1]
  unfairness_2_init_is <- abs(mean(OTR_ml_fr_init_is[S_is_2==1]) - mean(OTR_ml_fr_init_is[S_is_2==0]))
  
  # condition3: E[d=1|S_is_3=1] - E[d=1|S_is_0=1]
  unfairness_3_init_is <- abs(mean(OTR_ml_fr_init_is[S_is_3==1]) - mean(OTR_ml_fr_init_is[S_is_3==0]))
  
  # get the average unfairness
  average_unfairness_ml_fr_init_is <- (unfairness_1_init_is + unfairness_2_init_is + unfairness_3_init_is) / 3
  
  # set the range of the deltas
  delta_set <- c(0.0001,0.05, 
                 seq(0.1, 2, by = 0.1), 
                 seq(2,4, by = 0.5))
  
  # create an empty set to store the results
  value_ml_fr_indv_set <- c()
  RU_ml_fr_indv_set <- c()
  AU_ml_fr_indv_set <- c()
  
  value_ml_fr_indv_cluster_set <- c()
  RU_ml_fr_indv_cluster_set <- c()
  AU_ml_fr_indv_cluster_set <- c()
  
  FURG_indv_set <- c()
  FUTR_indv_set <- c()
  FURG_indv_cluster_set <- c()
  FUTR_indv_cluster_set <- c()
  
  # create an empty set to store the results from intersectional sensitive variables
  value_ml_fr_is_set <- c()
  RU_ml_fr_is_set <- c()
  AU_ml_fr_is_set <- c()
  
  FURG_is_set <- c()
  FUTR_is_set <- c()
  
  # PAN: 2024.07.29
  RU_ml_fr_indv_CATE_set <- c()
  AU_ml_fr_indv_CATE_set <- c()
  
  RU_ml_fr_indv_cluster_CATE_set <- c()
  AU_ml_fr_indv_cluster_CATE_set <- c()
  
  RU_ml_fr_is_CATE_set <- c()
  AU_ml_fr_is_CATE_set <- c()
  
  # using a for loop to get all the results
  for (delta in delta_set) {
    
    # 1.4.1 
    # fair constraints on the individual level only
    ml_fr_out <- fairCATE_multilevel(data = dat0,
                                     sensitive = c("S1"),
                                     legitimate = NULL,
                                     fairness = c("tau~S1"),
                                     treatment = "A",
                                     outcome = "Y",
                                     cluster = "id",
                                     multicategorical = NULL,
                                     outcome.LMM = glmm_out$outcome.LMM,
                                     ps.GLMM = glmm_out$ps.GLMM,
                                     fixed_intercept = FALSE,
                                     delta = c(delta),
                                     ps.trim="Sturmer.1")
    
    # ml_fr_out <- ML_FRCATE_alter_v10(data = dat0_is,
    #                                  data_alter = dat0,
    #                                  sensitive = c("S1"),
    #                                  legitimate = NULL,
    #                                  fairness = c("tau~S1"),
    #                                  treatment = "A",
    #                                  outcome = "Y",
    #                                  cluster = "id",
    #                                  multicategorical = NULL,
    #                                  outcome.LMM = glmm_out_is$outcome.LMM,
    #                                  ps.GLMM = glmm_out_is$ps.GLMM,
    #                                  delta = c(delta),
    #                                  ps.trim="Sturmer.1")
    
    tau_hat <- ml_fr_out$tau_hat
    # get the OTR, value
    OTR_ml_fr <- as.numeric(tau_hat > 0)
    value_ml_fr_indv <- mean(ifelse(OTR_ml_fr > 0, Y1, Y0))
    
    # get the Relative Utility
    RU_ml_fr_indv <- (value_ml_fr_indv - value_random ) / value_random
    
    # condition1: E[d=1|S_is_1=1] - E[d=1|S_is_0=1]
    unfairness_1 <- abs(mean(OTR_ml_fr[S_is_1==1]) - mean(OTR_ml_fr[S_is_1==0]))
    
    # condition2: E[d=1|S_is_2=1] - E[d=1|S_is_0=1]
    unfairness_2 <- abs(mean(OTR_ml_fr[S_is_2==1]) - mean(OTR_ml_fr[S_is_2==0]))
    
    # condition3: E[d=1|S_is_3=1] - E[d=1|S_is_0=1]
    unfairness_3 <- abs(mean(OTR_ml_fr[S_is_3==1]) - mean(OTR_ml_fr[S_is_3==0]))
    
    # get the average unfairness
    average_unfairness_ml_fr_indv <- (unfairness_1 + unfairness_2 + unfairness_3) / 3
    
    # calculate the UG
    UG_indv <- (value_ml_fr_indv - value_ml_fr_init) / (value_ml_fr_init - value_random)
    if(UG_indv == 0){
      UG_indv <- -0.01
    }
    
    
    # calculate the UD
    # PAN: 2024.07.25
    # revised the UD
    UD_indv_1 <- (unfairness_1_init - unfairness_1) / unfairness_1_init
    UD_indv_2 <- (unfairness_2_init - unfairness_2) / unfairness_2_init
    UD_indv_3 <- (unfairness_3_init - unfairness_3) / unfairness_3_init
    UD_indv <- (UD_indv_1 + UD_indv_2 + UD_indv_3) / 3
    
    FURG_indv <- UG_indv + UD_indv
    
    FUTR_indv <- - UD_indv / UG_indv
    
    # PAN: 2024.07.29
    # get the CATE_MSE
    MSE_CATE_ml_fr_indv <- MSE_CATE(tau_hat, tau)
    # get the CATE unfairness
    unfairness_1 <- abs(mean(tau_hat[S_is_1==1]) - mean(tau_hat[S_is_1==0]))
    unfairness_2 <- abs(mean(tau_hat[S_is_2==1]) - mean(tau_hat[S_is_2==0]))
    unfairness_3 <- abs(mean(tau_hat[S_is_3==1]) - mean(tau_hat[S_is_3==0]))
    
    average_unfairness_ml_fr_indv_CATE <- (unfairness_1 + unfairness_2 + unfairness_3) / 3
    
    # 1.4.2 
    # fair constraints on both individual and cluster levels
    ml_fr_out <- fairCATE_multilevel(data = dat0,
                                     sensitive = c("S1","S2"),
                                     legitimate = NULL,
                                     fairness = c("tau~S1", "tau~S2"),
                                     treatment = "A",
                                     outcome = "Y",
                                     cluster = "id",
                                     multicategorical = NULL,
                                     outcome.LMM = glmm_out$outcome.LMM,
                                     ps.GLMM = glmm_out$ps.GLMM,
                                     fixed_intercept = FALSE,
                                     delta = c(delta, delta),
                                     ps.trim="Sturmer.1")
    
    # ml_fr_out <- ML_FRCATE_alter_v10(data = dat0_is,
    #                                  data_alter = dat0,
    #                                  sensitive = c("S1","S2"),
    #                                  legitimate = NULL,
    #                                  fairness = c("tau~S1", "tau~S2"),
    #                                  treatment = "A",
    #                                  outcome = "Y",
    #                                  cluster = "id",
    #                                  multicategorical = NULL,
    #                                  outcome.LMM = glmm_out_is$outcome.LMM,
    #                                  ps.GLMM = glmm_out_is$ps.GLMM,
    #                                  delta = c(delta, delta),
    #                                  ps.trim="Sturmer.1")
    # 
    tau_hat <- ml_fr_out$tau_hat
    # get the OTR, value
    OTR_ml_fr <- as.numeric(tau_hat > 0)
    value_ml_fr_indv_cluster <- mean(ifelse(OTR_ml_fr > 0, df$Y1, df$Y0))
    
    # get the Relative Utility
    RU_ml_fr_indv_cluster <- (value_ml_fr_indv_cluster - value_random ) / value_random
    
    # ---------------------------------------------------------------
    # condition1: E[d=1|S_is_1=1] - E[d=1|S_is_0=1]
    unfairness_1 <- abs(mean(OTR_ml_fr[S_is_1==1]) - mean(OTR_ml_fr[S_is_1==0]))
    
    # condition2: E[d=1|S_is_2=1] - E[d=1|S_is_0=1]
    unfairness_2 <- abs(mean(OTR_ml_fr[S_is_2==1]) - mean(OTR_ml_fr[S_is_2==0]))
    
    # condition3: E[d=1|S_is_3=1] - E[d=1|S_is_0=1]
    unfairness_3 <- abs(mean(OTR_ml_fr[S_is_3==1]) - mean(OTR_ml_fr[S_is_3==0]))
    
    # get the average unfairness
    average_unfairness_ml_fr_indv_cluster <- (unfairness_1 + unfairness_2 +unfairness_3) / 3
    
    # calculate the UG
    UG_indv_cluster <- (value_ml_fr_indv_cluster - value_ml_fr_init) / (value_ml_fr_init - value_random)
    
    if(UG_indv_cluster == 0){
      UG_indv_cluster <- -0.01
    }
    
    # calculate the UD
    UD_1 <- (unfairness_1_init - unfairness_1) / unfairness_1_init
    UD_2 <- (unfairness_2_init - unfairness_2) / unfairness_2_init
    UD_3 <- (unfairness_3_init - unfairness_3) / unfairness_3_init
    
    UD_indv_cluster <- (UD_1 + UD_2 + UD_3) / 3
    
    FURG_indv_cluster <- UG_indv_cluster + UD_indv_cluster
    
    FUTR_indv_cluster <- - UD_indv_cluster / UG_indv_cluster
    
    # PAN: 2024.07.29
    # get the CATE_MSE
    MSE_CATE_ml_fr_indv_cluster <- MSE_CATE(tau_hat, tau)
    # get the CATE unfairness
    unfairness_1 <- abs(mean(tau_hat[S_is_1==1]) - mean(tau_hat[S_is_1==0]))
    unfairness_2 <- abs(mean(tau_hat[S_is_2==1]) - mean(tau_hat[S_is_2==0]))
    unfairness_3 <- abs(mean(tau_hat[S_is_3==1]) - mean(tau_hat[S_is_3==0]))
    
    average_unfairness_ml_fr_indv_cluster_CATE <- (unfairness_1 + unfairness_2 + unfairness_3) / 3
    
    
    # 1.4.3 
    # fairness constrains on the inter-sectional sensitive variables
    
    ml_fr_out_is <- fairCATE_multilevel(data = dat0_is,
                                        sensitive = c("S_is_1","S_is_2","S_is_3"),
                                        legitimate = NULL,
                                        fairness = c("tau~S_is_1","tau~S_is_2","tau~S_is_3"),
                                        treatment = "A",
                                        outcome = "Y",
                                        cluster = "id",
                                        multicategorical = NULL,
                                        outcome.LMM = glmm_out_is$outcome.LMM,
                                        ps.GLMM = glmm_out_is$ps.GLMM,
                                        fixed_intercept = FALSE,
                                        delta = c(delta, delta, delta),
                                        ps.trim="Sturmer.1")
    
    tau_hat_is <- ml_fr_out_is$tau_hat
    # get the OTR, value
    OTR_ml_fr_is <- as.numeric(tau_hat_is > 0)
    value_ml_fr_is <- mean(ifelse(OTR_ml_fr_is > 0, Y1, Y0))
    
    # get the Relative Utility
    RU_ml_fr_is <- (value_ml_fr_is - value_random ) / value_random
    
    # condition1: E[d=1|S_is_1=1] - E[d=1|S_is_0=1]
    unfairness_1_is <- abs(mean(OTR_ml_fr_is[S_is_1==1]) - mean(OTR_ml_fr_is[S_is_1==0]))
    
    # condition2: E[d=1|S_is_2=1] - E[d=1|S_is_0=1]
    unfairness_2_is <- abs(mean(OTR_ml_fr_is[S_is_2==1]) - mean(OTR_ml_fr_is[S_is_2==0]))
    
    # condition3: E[d=1|S_is_3=1] - E[d=1|S_is_0=1]
    unfairness_3_is <- abs(mean(OTR_ml_fr_is[S_is_3==1]) - mean(OTR_ml_fr_is[S_is_3==0]))
    
    # get the average unfairness
    average_unfairness_ml_fr_is <- (unfairness_1_is + unfairness_2_is + unfairness_3_is) / 3
    
    # calculate the UG
    UG_is <- (value_ml_fr_is - value_ml_fr_init_is) / (value_ml_fr_init_is - value_random)
    
    if(UG_is == 0){
      UG_is <- -0.01
    }
    
    # calculate the UD
    UD_1_is <- (unfairness_1_init_is - unfairness_1_is) / unfairness_1_init_is
    UD_2_is <- (unfairness_2_init_is - unfairness_2_is) / unfairness_2_init_is
    UD_3_is <- (unfairness_3_init_is - unfairness_3_is) / unfairness_3_init_is
    
    UD_is <- (UD_1_is + UD_2_is + UD_3_is) / 3
    
    FURG_is <- UG_is + UD_is
    
    FUTR_is <- - UD_is / UG_is
    
    # PAN: 2024.07.29
    # get the CATE_MSE
    MSE_CATE_ml_fr_is <- MSE_CATE(tau_hat_is, tau)
    # get the CATE unfairness
    unfairness_1 <- abs(mean(tau_hat_is[S_is_1==1]) - mean(tau_hat_is[S_is_1==0]))
    unfairness_2 <- abs(mean(tau_hat_is[S_is_2==1]) - mean(tau_hat_is[S_is_2==0]))
    unfairness_3 <- abs(mean(tau_hat_is[S_is_3==1]) - mean(tau_hat_is[S_is_3==0]))
    
    average_unfairness_ml_fr_is_CATE <- (unfairness_1 + unfairness_2 + unfairness_3) / 3
    
    
    # store the results
    value_ml_fr_indv_set <- c(value_ml_fr_indv_set, value_ml_fr_indv)
    RU_ml_fr_indv_set <- c(RU_ml_fr_indv_set, RU_ml_fr_indv)
    AU_ml_fr_indv_set <- c(AU_ml_fr_indv_set, average_unfairness_ml_fr_indv)
    
    value_ml_fr_indv_cluster_set <- c(value_ml_fr_indv_cluster_set, value_ml_fr_indv_cluster)
    RU_ml_fr_indv_cluster_set <- c(RU_ml_fr_indv_cluster_set, RU_ml_fr_indv_cluster)
    AU_ml_fr_indv_cluster_set <- c(AU_ml_fr_indv_cluster_set, average_unfairness_ml_fr_indv_cluster)
    
    FURG_indv_set <- c(FURG_indv_set, FURG_indv)
    FUTR_indv_set <- c(FUTR_indv_set, FUTR_indv)
    FURG_indv_cluster_set <- c(FURG_indv_cluster_set, FURG_indv_cluster)
    FUTR_indv_cluster_set <- c(FUTR_indv_cluster_set, FUTR_indv_cluster)
    
    value_ml_fr_is_set <- c(value_ml_fr_is_set, value_ml_fr_is)
    RU_ml_fr_is_set <- c(RU_ml_fr_is_set, RU_ml_fr_is)
    AU_ml_fr_is_set <- c(AU_ml_fr_is_set, average_unfairness_ml_fr_is)
    
    FURG_is_set <- c(FURG_is_set, FURG_is)
    FUTR_is_set <- c(FUTR_is_set, FUTR_is)
    
    # PAN: 2024.07.29
    RU_ml_fr_indv_CATE_set <- c(RU_ml_fr_indv_CATE_set, MSE_CATE_ml_fr_indv)
    AU_ml_fr_indv_CATE_set <- c(AU_ml_fr_indv_CATE_set, average_unfairness_ml_fr_indv_CATE)
    
    RU_ml_fr_indv_cluster_CATE_set <- c(RU_ml_fr_indv_cluster_CATE_set, MSE_CATE_ml_fr_indv_cluster)
    AU_ml_fr_indv_cluster_CATE_set <- c(AU_ml_fr_indv_cluster_CATE_set, average_unfairness_ml_fr_indv_cluster_CATE)
    
    RU_ml_fr_is_CATE_set <- c(RU_ml_fr_is_CATE_set, MSE_CATE_ml_fr_is)
    AU_ml_fr_is_CATE_set <- c(AU_ml_fr_is_CATE_set, average_unfairness_ml_fr_is_CATE)
    
  }
  # store the results into the list
  list(value_true = value_true,
       value_random = value_random,
       value_BART = value_BART,
       value_cf = value_cf,
       RU_true = RU_true,
       RU_BART = RU_BART,
       RU_cf = RU_cf,
       RU_cf_ps = RU_cf_ps,
       RU_BART_CATE = MSE_CATE_BART,
       RU_cf_CATE = MSE_CATE_cf,
       RU_cf_CATE_ps = MSE_CATE_cf_ps,
       average_unfairness_true = average_unfairness_true,
       average_unfairness_bart = average_unfairness_bart,
       average_unfairness_cf = average_unfairness_cf,
       average_unfairness_cf_ps = average_unfairness_cf_ps,
       average_unfairness_bart_CATE = average_unfairness_bart_CATE,
       average_unfairness_cf_CATE = average_unfairness_cf_CATE,
       average_unfairness_cf_CATE_ps = average_unfairness_cf_CATE_ps,
       value_ml_fr_indv_set = value_ml_fr_indv_set,
       RU_ml_fr_indv_set = RU_ml_fr_indv_set,
       AU_ml_fr_indv_set = AU_ml_fr_indv_set,
       value_ml_fr_indv_cluster_set = value_ml_fr_indv_cluster_set,
       RU_ml_fr_indv_cluster_set = RU_ml_fr_indv_cluster_set,
       AU_ml_fr_indv_cluster_set = AU_ml_fr_indv_cluster_set,
       FURG_indv_set = FURG_indv_set,
       FUTR_indv_set = FUTR_indv_set,
       FURG_indv_cluster_set = FURG_indv_cluster_set,
       FUTR_indv_cluster_set = FUTR_indv_cluster_set,
       value_ml_fr_is_set = value_ml_fr_is_set,
       RU_ml_fr_is_set = RU_ml_fr_is_set,
       AU_ml_fr_is_set = AU_ml_fr_is_set,
       FURG_is_set = FURG_is_set,
       FUTR_is_set = FUTR_is_set,
       RU_ml_fr_indv_CATE_set = RU_ml_fr_indv_CATE_set,
       AU_ml_fr_indv_CATE_set = AU_ml_fr_indv_CATE_set,
       RU_ml_fr_indv_cluster_CATE_set = RU_ml_fr_indv_cluster_CATE_set,
       AU_ml_fr_indv_cluster_CATE_set = AU_ml_fr_indv_cluster_CATE_set,
       RU_ml_fr_is_CATE_set = RU_ml_fr_is_CATE_set,
       AU_ml_fr_is_CATE_set = AU_ml_fr_is_CATE_set)
}


time_0 <- Sys.time()

cl <- makeCluster(detectCores() - 1)
clusterEvalQ(cl, {
  library(MASS)
  library(matrixcalc)
  library(mbend)
  library(Matrix)
  library(Rmosek)
  library(lme4)
  library(bartCause)
  library(grf)
  library(dplyr)
  source("data_generating_process_multilevel.r")
  source("fairCATE_multilevel_functions.r")
})

clusterExport(cl, c("run_iteration", "create_multileveldata_D2", "GLMM_model", 
                    "fairCATE_multilevel"))
results <- parLapply(cl, 1:20, run_iteration)
stopCluster(cl)
# save the results
save(results, file = "results/Simulation_Design_2_results.rda")

time_1 <- Sys.time()

time_1-time_0


