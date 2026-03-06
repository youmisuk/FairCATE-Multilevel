library(MASS)
library(matrixcalc)
library(mbend)
library(Matrix)
library(Rmosek)
library(lme4)
library(fastDummies)



#####################################################
#                                                   #
#           DGP for design 1                        #
#           version: v7, 2025.04.27                 #
#                                                   #
#####################################################

create_multileveldata_D1 <- function(cluster_num = 300, # number of clusters
                                    cluster_size_min = 10,
                                    cluster_size_max = 50,
                                    R_var = .6653, # variance of residual
                                    V_var = 0.8,  # variance of cluster effect in selection
                                    U_var = 0.0776, # variance of cluster effect in outcome
                                    U_V_same = FALSE, # whether the cluster effect in selection and outcome are the same
                                    tau_var = 0, # variance of the cluster level tau
                                    U_var_coef = 1, # coefficient before the U_var
                                    clustereffect=TRUE){
  
  # DGP v7 version of the multilevel data generation function    
  #   - adding the cluster level random effect in tau            
  #   - adding repeat process to ensure the S1 in each cluster   
  #     has at least one 0 and one 1  
  
  # ::::: 1) generate the number of observations in each cluster :::::
  
  balanced <- FALSE
  
  while (balanced == FALSE) {
    
    J <- cluster_num # level-2 unit, the number of cluster
    
    # v4 updated part
    n.clus <- round(runif(J, min = cluster_size_min , max = cluster_size_max ), 0)
    
    id <- as.factor(rep(1:J, times=n.clus))              # cluster id
    
    # ::::: 2) generate level-2 covariates, W: X21,X22,X23 :::::
    # cluster level sensitive variable, should be historically disadvantaged minority group
    S2 <- rbinom(J, size = 1, prob = 0.3)
    
    S2_2 <- rbinom(J, size = 1, prob = 0.6)
    
    # cluster level covariate 1
    X21 <- rnorm(J, 0, 1)  
    
    # cluster level covariate 2 influenced by S2
    X22 <- rnorm(J, -0.1*S2 + 0.5, 0.5)
    
    # cluster level covariate 3 influenced by S2, like the Government's funding for the 
    # protected groups
    
    X23 <- rnorm(J, - 0.3 * S2, 0.5)
    
    names(X21) <- names(X22) <- names(X23) <- names(S2) <- names(S2_2) <- levels(id) 
    
    # ::::: 3) generate level-1 covariates, Xs with cluster-specific means :::::
    
    totalN <- length(id)
    # the protected group might be the minority group
    # S1 <- rbinom(totalN, size = 1, prob = 0.4)  
    
    # 2025.03.28
    # Assign individual-level S1 probability based on their cluster's S2 status
    S2_for_each_individual <- S2[id]  # get the cluster-level S2 for each individual
    
    # Repeat S1 generation until all clusters have at least one 0 and one 1
    repeat {
      
      # Define conditional probabilities
      prob_S1_given_S2 <- ifelse(S2_for_each_individual == 1, 0.6, 0.4)
      
      # Generate S1
      S1 <- rbinom(totalN, size = 1, prob = prob_S1_given_S2)
      
      # Check if each cluster contains both 0 and 1 in S1
      library(dplyr)
      s1_check <- data.frame(cluster = id, S1 = S1) %>%
        group_by(cluster) %>%
        summarise(has_0 = any(S1 == 0), has_1 = any(S1 == 1)) %>%
        summarise(all_ok = all(has_0 & has_1)) %>%
        pull(all_ok)
      
      # Exit the loop if the condition is satisfied
      if (s1_check) break
      
    }
    
    # individual level covariate 1
    X11 <- rnorm(totalN, 0, 1)
    
    # individual level covariate 2 influence by S1
    X12 <- rnorm(totalN, -2 * S1 + 1, 1)
    
    # individual level covariate 3 influnced by S1
    X13 <- rnorm(totalN, -0.5 * S1, 1)
    
    # legitimate variable, the variable you wanna control for conditional statistical disparity
    L_continuous <- rnorm(totalN, -S1 + 0.25, 1)
    
    # discretize the L_continuous
    X14 <- ifelse(L_continuous > 0.1, 1, 0)
    
    pop <- data.frame(id, X11, X12, X13, S1, X14, 
                      X21=X21[id], X22=X22[id], X23=X23[id], S2=S2[id], S2_2=S2_2[id])
    
    # ::::: 4) generate selection probabilities and potential outcome ::::::::
    
    E <- rnorm(totalN, 0, sqrt(R_var))   # error terms for pot.   
    
    if (clustereffect == FALSE) {
      
      # Include all the covariates as moderators except the S2
      pop$lps <- -0.7 + 0.3*(pop$X11 + pop$X12 + pop$X13 + pop$S1 + pop$X14) + 
        0.3*(pop$X21 + pop$X22 + pop$X23 - pop$S2 - pop$S2_2)
      
      pop$Y0  <- 4 + 0.4*(pop$X11  + pop$X12 + pop$X13 + pop$S1 + pop$X21 - 
                            pop$X22 + pop$X23 - pop$S2 + pop$X14) + E
      
      pop$Y1  <- pop$Y0 + 0.8 + 0.3 * pop$X11 + 0.2 * pop$X12 + 0.2 * pop$X13 - 
        0.5*pop$S1 + 0.5*pop$X21 + 0.4*pop$X22 + 0.2*pop$X23 +  0.1*pop$X14 - 0.2*pop$S2 - 0.2*pop$S2_2
      
    } else {
      # adding the cluster-specific random effect and increase the variance
      R_j <- rnorm(J, 0, sqrt(V_var)) # level-2 cluster effect in selection
      
      if (U_V_same){
        U_j <- R_j # level-2 cluster effect in outcome
      } else {
        U_j <- rnorm(J, 0, sqrt(U_var)) # level-2 cluster effect in outcome
      }
      
      tau_var_j <- rnorm(J, 0, sqrt(tau_var)) # level-2 cluster effect in tau
      
      pop$R_j <- R_j[id]
      pop$U_j <- U_j[id]
      pop$tau_var_j <- tau_var_j[id]
      
      pop$lps <- -0.1 + 0.2*(pop$X11 + pop$X12 + pop$X13 + pop$S1 + pop$X14) +
        -0.2*(pop$X21 + pop$X22 + pop$X23 - pop$S2 + pop$S2_2)  + pop$R_j
      
      # pop$lps <- -1.7 + 0.3*(pop$X11 + pop$X12 + pop$X13 + pop$S1 + pop$X14) + 
      #   0.3*(pop$X21 + pop$X22 + pop$X23 - pop$S2)  + pop$R_j
      # 
      
      pop$Y0  <- 4 + 0.4*(pop$X11  + pop$X12 + pop$X13 + pop$S1 + pop$X21 -
                            pop$X22 + pop$X23 - pop$S2 + pop$X14) + E + U_var_coef*pop$U_j
      
      pop$Y1  <- pop$Y0 + 0.8 + 0.3 * pop$X11 + 0.2 * pop$X12 + 0.2 * pop$X13 - 
        0.5*pop$S1 + 0.5*pop$X21 + 0.4*pop$X22 + 0.2*pop$X23 +  0.1*pop$X14 + 0.8*pop$S2 - 0.2*pop$S2_2 + pop$tau_var_j
    }
    # changed  S2's beta from -0.2 to +0.8
    
    pop$tau <- pop$Y1 - pop$Y0
    pop$policy <- as.numeric(pop$tau > 0)
    
    pop$ps <- 1 / (1 + exp(-pop$lps))   # propensity scores
    pop$A <- rbinom(totalN, 1, pop$ps) # treatment indicator  
    
    # 2025.03.28
    # adding criterion to make sure that in each cluster, there are at least one treated and one control
    
    cluster_ids <- unique(pop$id)
    pop$Y <- with(pop, ifelse(A == 1, Y1, Y0))
    
    # detect whether there is any cluster with single treatment assignment
    detecter <- lapply(cluster_ids, function(clust) {
      ifelse (length(unique(pop$A[pop$id == clust])) == 2, TRUE, FALSE)
    })
    
    balanced <- all(unlist(detecter)) # if all TRUE,balanced = TRUE, loop stops
    message("Balanced check not pass, re-generate the data...")
  }
  
  message("Balanced check pass")
  
  return(pop)
}



#####################################################
#                                                   #
#           GLMM model fitting function             #
#                                                   #
#####################################################

GLMM_model <- function(data,
                       outcome,
                       treatment,
                       cluster,
                       multicategorical=NULL,
                       n_AGQ = 1,
                       fixed_intercept = TRUE,
                       random_slope = FALSE, # for the outcome model
                       glmer_Control=NULL){
  # note: 2025.03.26
  # the fixed intercept only influence the outcome model
  # the propensity model always has a fixed intercept.
  
  
  # change the variable name
  if(outcome != "Y"){names(data)[which(names(data) == outcome)] <- "Y"}
  if(treatment != "A"){names(data)[which(names(data) == treatment)] <- "A"}
  
  # if(legitimate != "L"){names(data)[which(names(data) == legitimate)] <- "L"}
  
  if(!is.null(multicategorical)){
    if(!(all(sapply(data[, multicategorical], is.factor)))){
      stop("The variable(s) in multicategorical should be in factor format. \n
         You don't need to include the binary coded factor variable. \n
         Only whose category level greater than 2 should be included.\n")
    }
  }
  
  if(class(data[,cluster]) != "character"){
    stop("In case of misusing of nominal scale of Cluster variable as numeric values, please convert it as character first!")}
  
  if (fixed_intercept){
    # Handling the multi-categorical variables, i.e., varaibles specified in reference
    df_mc_dummied <- as.data.frame(matrix(1, nrow = nrow(data), ncol = 1))
    
    # Convert them into a dummy coded dataframe 
    for (mc_var in multicategorical) {
      df_mc <- as.data.frame(data[,mc_var])
      names(df_mc) <-  mc_var
      df_mc <- dummy_cols(df_mc, select_columns = c(mc_var),
                          ignore_na = T, remove_selected_columns = T)
      # drop the last level within each variable to avoid Singularity matrix
      df_mc <- df_mc[,-ncol(df_mc)]
      # attach the dummy coded variable into the dataframe
      df_mc_dummied <- cbind(df_mc_dummied, df_mc)
    }
    
    df_mc_dummied <- df_mc_dummied[,-1]
    
    # exclude the non-covariates
    exclude_vars <- c(cluster,multicategorical,"A")
    include_vars <- setdiff(names(data),c("Y", exclude_vars))
    
    # merge the dummied multicategorical covariataes matrix with the other covariates
    b.mat <- as.matrix(cbind(1,cbind(data[include_vars],df_mc_dummied)))
    
    A <- data$A
    C <- data[, cluster]
    Y <- data$Y
    
    # convert b.mat to a data.frame
    b.mat.df <- as.data.frame(b.mat)
    
    # fit the LMM model with interaction terms
    message("Fitting the LMM model for the outcome with interaction terms....")
    
    # include the interaction between the A and covariates
    # using [-1] to drop the intercept term
    fixed_effects <- paste("A *", names(b.mat.df)[-1], collapse = " + ")
    
    # Combine the fixed effects with the random effects
    if (random_slope){
      formula <- as.formula(paste("Y ~ ", fixed_effects, "+ (1 + A| C)"))
    } else{
      formula <- as.formula(paste("Y ~ ", fixed_effects, "+ (1 | C)"))
    }
    
    # Fit the model using the lmer function
    outcome.LMM <- lmer(formula, data = cbind(Y, A, C, b.mat.df))
    
    message("...Finished!")
    
    # fit the GLMM model
    message("Fitting the GLMM model for propensity score (this might run over 4 mins)....")
    ps.GLMM <- glmer(A ~ 0 + b.mat + (1|C), family="binomial",
                     nAGQ = n_AGQ,control = glmer_Control) # propensity score
    message("...Finished!")
    
  } else {
    
    # Handling the multi-categorical variables, i.e., varaibles specified in reference
    df_mc_dummied <- as.data.frame(matrix(1, nrow = nrow(data), ncol = 1))
    
    # Convert them into a dummy coded dataframe 
    for (mc_var in multicategorical) {
      df_mc <- as.data.frame(data[,mc_var])
      names(df_mc) <-  mc_var
      df_mc <- dummy_cols(df_mc, select_columns = c(mc_var),
                          ignore_na = T, remove_selected_columns = T)
      # drop the last level within each variable to avoid Singularity matrix
      df_mc <- df_mc[,-ncol(df_mc)]
      # attach the dummy coded variable into the dataframe
      df_mc_dummied <- cbind(df_mc_dummied, df_mc)
    }
    df_mc_dummied <- df_mc_dummied[,-1]
    
    # exclude the non-covariates
    exclude_vars <- c(cluster,multicategorical,"A")
    include_vars <- setdiff(names(data),c("Y", exclude_vars))
    
    # merge the dummied multicategorical covairtaes matrix with the other covariates
    b.mat <- as.matrix(cbind(1,cbind(data[include_vars],df_mc_dummied)))
    
    A <- data$A
    C <- data[, cluster]
    Y <- data$Y
    
    # convert b.mat to a data.frame
    b.mat.df <- as.data.frame(b.mat)
    
    # fit the LMM model with interaction terms
    message("Fitting the LMM model for the outcome with interaction terms....")
    
    # Include the interaction between the A and covariates
    # using [-1] to drop the intercept term
    fixed_effects <- paste("A *", names(b.mat.df)[-1], collapse = " + ")
    
    if (random_slope){
      formula <- as.formula(paste("Y ~ ", fixed_effects, "+ (1 + A| C)"))
    } else{
      formula <- as.formula(paste("Y ~ ", fixed_effects, "+ (1 | C)"))
    }
    
    # Fit the model using the lmer function
    outcome.LMM <- lmer(formula, data = cbind(Y, A, C, b.mat.df))
    
    message("...Finished!")
    # fit the GLMM model
    message("Fitting the GLMM model for propensity score (this might run over 4 mins)....")
    ps.GLMM <- glmer(A ~ 0 + b.mat + (1|C), family="binomial",
                     nAGQ = n_AGQ,control = glmer_Control) # propensity score
    message("...Finished!")
  }
  
  
  return(list(
    outcome.LMM = outcome.LMM,
    ps.GLMM = ps.GLMM))
}



#####################################################
#                                                   #
#           GLM model fitting function              #
#                                                   #
#####################################################

GLM_model <- function(data,
                      outcome,
                      treatment,
                      cluster,
                      multicategorical=NULL,
                      fixed_intercept = TRUE){
  # this function is used to fit linear regression for outcome model
  # and logistic regression for propensity score model without considering the cluster effect
  # Note: the fixed_intercept argument influnece the outcome model only.
  # version 1: 2025.03.26
  
  # change the variable name
  if(outcome != "Y"){names(data)[which(names(data) == outcome)] <- "Y"}
  if(treatment != "A"){names(data)[which(names(data) == treatment)] <- "A"}
  
  # if(legitimate != "L"){names(data)[which(names(data) == legitimate)] <- "L"}
  
  if(!is.null(multicategorical)){
    if(!(all(sapply(data[, multicategorical], is.factor)))){
      stop("The variable(s) in multicategorical should be in factor format. \n
     You don't need to include the binary coded factor variable. \n
     Only whose category level greater than 2 should be included.\n")
    }
  }
  
  if(class(data[,cluster]) != "character"){
    stop("In case of misusing of nominal scale of Cluster variable as numeric values, please convert it as character first!")}
  
  
  
  # Handling the multi-categorical variables, i.e., varaibles specified in reference
  df_mc_dummied <- as.data.frame(matrix(1, nrow = nrow(data), ncol = 1))
  
  # Convert them into a dummy coded dataframe 
  for (mc_var in multicategorical) {
    df_mc <- as.data.frame(data[,mc_var])
    names(df_mc) <-  mc_var
    df_mc <- dummy_cols(df_mc, select_columns = c(mc_var),
                        ignore_na = T, remove_selected_columns = T)
    # drop the last level within each variable to avoid Singularity matrix
    df_mc <- df_mc[,-ncol(df_mc)]
    # attach the dummy coded variable into the dataframe
    df_mc_dummied <- cbind(df_mc_dummied, df_mc)
  }
  df_mc_dummied <- df_mc_dummied[,-1]
  
  # exclude the non-covariates
  exclude_vars <- c(cluster,multicategorical,"A")
  include_vars <- setdiff(names(data),c("Y", exclude_vars))
  
  # merge the dummied multicategorical covairtaes matrix with the other covariates
  b.mat <- as.matrix(cbind(cbind(data[include_vars],df_mc_dummied)))
  
  A <- data$A
  C <- data[, cluster]
  Y <- data$Y
  
  # convert b.mat to a data.frame
  b.mat.df <- as.data.frame(b.mat)
  head(b.mat.df)
  
  if(fixed_intercept){
    message("Fixed intercept is included in the model.")
    # fit the LM model with interaction terms
    message("Fitting the LM model for the outcome with interaction terms....")
    
    # Include the interaction between the A and covariates
    # using [-1] to drop the intercept term
    fixed_effects <- paste("A *", names(b.mat.df), collapse = " + ")
    
    # Combine the fixed effects with the random effects
    formula <- as.formula(paste("Y ~ ", fixed_effects))
    
    # Fit the model using the lmer function
    outcome.LM <- lm(formula, data = cbind(Y, A, C, b.mat.df))
    
    
    message("...Finished!")
    
    # fit the logistic model
    message("Fitting the GLM model for propensity score ....")
    ps.GLM <- glm(A ~ b.mat, family="binomial") # propensity score
    
    summary(ps.GLM)
    message("...Finished!")
    
  } else {
    message("Fixed intercept is excluded in the model.")
    # Include the interaction between the A and covariates
    # using [-1] to drop the intercept term
    fixed_effects <- paste("A *", names(b.mat.df), collapse = " + ")
    
    # Combine the fixed effects with the random effects
    formula <- as.formula(paste("Y ~ 0 + ", fixed_effects))
    
    # Fit the model using the lmer function
    outcome.LM <- lm(formula, data = cbind(Y, A, C, b.mat.df))
    
    message("...Finished!")
    
    # fit the logistic model
    message("Fitting the GLM model for propensity score ....")
    ps.GLM <- glm(A ~ b.mat, family="binomial") # propensity score
    
    summary(ps.GLM)
    message("...Finished!")
  }
  return(list(
    outcome.LM = outcome.LM,
    ps.GLM = ps.GLM))
}


###############################################################################
#                                                                             #
#           multilevel fair CATE function ignoring the cluster effect         #
#                                                                             #
###############################################################################



fairCATE_multilevel_cluster_ignored<- function(data,
                                               sensitive,
                                               legitimate=NULL,
                                               fairness,
                                               treatment,
                                               outcome,
                                               cluster,
                                               multicategorical=NULL,
                                               outcome.LMM,
                                               ps.GLMM,
                                               fixed_intercept = TRUE,
                                               random_slope = FALSE,
                                               delta,
                                               ps.trim="Sturmer.1"){
  
  # this function is used after GLM_model function, the argument `fixed_intercept`
  # must match with each other.
  # version 1: 2025.03.26
  
  # change the variable name
  if(outcome != "Y"){names(data)[which(names(data) == outcome)] <- "Y"}
  if(treatment != "A"){names(data)[which(names(data) == treatment)] <- "A"}
  
  if(length(fairness) != length(delta)){
    stop("fariness and delta must have the same length!")}
  
  if(!is.null(multicategorical)){
    if(!(all(sapply(data[, multicategorical], is.factor)))){
      stop("The variable(s) in multicategorical should be in factor format. \n
         You don't need to include the binary coded factor variable. \n
         Only whose category level greater than 2 should be included.")
    }
  }
  
  if(class(data[,cluster]) != "character"){
    stop("In case of misusing of nominal scale of Cluster variable as numeric values, please convert it as character first!")}
  
  merged_fair_specification <- paste0(fairness, collapse = " ")
  
  if (grepl("\\|", merged_fair_specification) && is.null(legitimate)){
    stop("It seems you plan to run on the conditional statistical disparity. \n
         But you forget to specify the legitimate variable into the argument. \n
         Please specify the legitimate variable(s) for the conditional statistical parity!")
  }
  
  # define the expit function
  expit <- function(x){ exp(x)/(1+exp(x))}
  
  # define the number of clusters
  m <- length(unique(data[,cluster]))
  
  # Handling the multicategorical variables, i.e., variables specified in reference
  df_mc_dummied <- as.data.frame(matrix(1, nrow = nrow(data), ncol = 1))
  
  # Convert them into a dummy coded dataframe 
  for (mc_var in multicategorical) {
    df_mc <- as.data.frame(data[,mc_var])
    names(df_mc) <-  mc_var
    df_mc <- dummy_cols(df_mc, select_columns = c(mc_var),
                        ignore_na = T, remove_selected_columns = T)
    # drop the last level within each variable to avoid Singularity matrix
    df_mc <- df_mc[,-ncol(df_mc)]
    # attach the dummy coded variable into the dataframe
    df_mc_dummied <- cbind(df_mc_dummied, df_mc)
  }
  
  df_mc_dummied <- df_mc_dummied[,-1]
  
  # exclude the non-covariates
  exclude_vars <- c(cluster,multicategorical,"A")
  include_vars <- setdiff(names(data),c("Y", exclude_vars))
  
  if (fixed_intercept){
    # merge the dummied multicategorical covariataes matrix with the other covariates
    b.mat <- as.matrix(cbind(1,cbind(data[include_vars], df_mc_dummied)))
  } else {
    b.mat <- as.matrix(cbind(data[include_vars], df_mc_dummied))
  }
  
  A <- data$A
  C <- data[, cluster]
  
  # cluster dummy matrix
  Cl.mat <- t(do.call(rbind, lapply(1:m, function(i) ifelse(C == i, 1, 0)))) 
  colnames(Cl.mat) <- paste0("C", 1:m)
  
  Y <- data$Y
  
  # create an empty L matrix
  L.mat <- matrix(nrow = dim(data)[1],ncol = 0)
  if (is.null(legitimate)) {
    L.mat <- NULL
    dat <- cbind(Cl.mat, data)
  } else {
    # combine all legitimate matrix into one unified L matrix
    for (i in 1:length(legitimate)) {
      # PAN: support adding the cluster ID as legitimate variable
      if (legitimate[i] == cluster){
        Cl.mat_ <- Cl.mat
        colnames(Cl.mat_) <- paste0(legitimate[i],"_", 1:(ncol(Cl.mat_)))
        L.mat <- cbind(L.mat, Cl.mat_)
      } else{
        L_ <- data[, legitimate[i]]
        L.mat_ <- model.matrix(~as.factor(L_)+0)
        # adding support for multicategorical legitimate 
        colnames(L.mat_) <- paste0(legitimate[i],"_", 0:(nlevels(as.factor(L_))-1)) #### 
        L.mat <- cbind(L.mat,L.mat_)
      }
    }
    dat <- cbind(Cl.mat, data, L.mat)
  }
  
  
  cluster_vec <- unique(data[, cluster])
  
  if (fixed_intercept){
    # create the interaction term in b.matrix
    interaction_A1 <- b.mat[, -1] * 1
    interaction_A0 <- b.mat[, -1] * 0
    
    # note: b.mat has the first column of 1
    outcome.A1.estimate <- as.numeric( cbind(1, 1, b.mat[,-1], interaction_A1)%*%coef(outcome.LMM)) 
    outcome.A0.estimate <- as.numeric( cbind(1, 0, b.mat[,-1], interaction_A0)%*%coef(outcome.LMM)) 
    log_odds <- b.mat %*% coef(ps.GLMM)
  } else {
    # create the interaction term in b.matrix
    interaction_A1 <- b.mat * 1
    interaction_A0 <- b.mat * 0
    
    # note: b.mat has the first column of 1
    outcome.A1.estimate <- as.numeric( cbind(1, b.mat, interaction_A1)%*%coef(outcome.LMM))
    outcome.A0.estimate <- as.numeric( cbind(0, b.mat, interaction_A0)%*%coef(outcome.LMM)) 
    log_odds <- cbind(1, b.mat) %*% coef(ps.GLMM)
  }
  
  # estimate of Pr(A=1|X_i,S_i) for all i
  message("Getting the estimated propensity scores...")
  ps.estimate <- 1 / (1 + exp(-log_odds))
  
  # PS trimming - removing the individuals with extreme PS values
  idx.exclude <- NULL
  if (ps.trim == "Sturmer.1") {
    # the common range method ver.1 by St??rmer et al. Am J Epidemiol 2021;190:1659???1670.
    idx.exclude <- c(which((ps.estimate < min(ps.estimate[A==1])) == TRUE),
                     which((max(ps.estimate[A==0]) < ps.estimate) == TRUE)
    )
  } 
  
  if (ps.trim == "Sturmer.2") {
    # the common range method ver.2 by St??rmer et al. Am J Epidemiol 2010;172:843???854.
    idx.exclude <- c(which((ps.estimate < quantile(ps.estimate[dat$A==1], probs = 0.05)) == TRUE),
                     which((quantile(ps.estimate[dat$A==0], probs = 0.95) < ps.estimate) == TRUE)
    )
  }
  
  message("Getting the cluster's weights matrix...")
  # in case we want to reweight on a cluster level later on
  phi.b.hat <- matrix(NA,nrow=m, ncol=ncol(b.mat)) 
  for(jj in 1:m){
    # extract the corresponding school id
    cluster_id <- cluster_vec[jj]
    if (is.null(idx.exclude)) {
      phi.b.hat[jj,] <-
        colMeans( (( A*(Y-outcome.A1.estimate) / ps.estimate - 
                       (1-A)*(Y-outcome.A0.estimate) / (1-ps.estimate) + 
                       (outcome.A1.estimate-outcome.A0.estimate) ) * b.mat)[C==cluster_id,] )
    } else {
      mat_sub <- (as.numeric(( A*(Y-outcome.A1.estimate) / ps.estimate - 
                                 (1-A)*(Y-outcome.A0.estimate) / (1-ps.estimate) + 
                                 (outcome.A1.estimate-outcome.A0.estimate) )) * b.mat)[C==cluster_id & !(1:length(C) %in% idx.exclude),]
      if ( "numeric" %in% class(mat_sub)) {
        # this condition handles the cluster with only one observation
        phi.b.hat[jj,] <- mat_sub}
      else if(all(is.na(mat_sub))==T){
        # this condition handles the cluster whose all observations are with trimmed propensity
        # will give a weight as 0
        phi.b.hat[jj,] <- rep(0,ncol(b.mat))
      } else{
        colmean <- colMeans(mat_sub)
        phi.b.hat[jj,] <- colmean
      }
    }
  }
  phi.b.hat <- colMeans(phi.b.hat)
  
  #########################################################################
  # Optimization                                                          #
  #########################################################################
  
  message("Optimizing using the given fariness condition...")
  
  Q.mat <- (t(b.mat) %*% b.mat)/nrow(b.mat)
  C.mat <- NULL
  
  k <- ncol(b.mat)
  N.tot <- nrow(b.mat)
  
  
  delta.v <- c()
  # based on the input fariness argument to build the delta vector
  if (is.list(delta)){
    delta.v <- unlist(delta)
  } else {
    for (i in 1:length(fairness)) {
      if (grepl("\\|",fairness[i])){
        
        # to make it works on the multicategorical legitimate variables
        L_temp <- sub(".*\\|", "", fairness[i])
        n_level <- nlevels(as.factor(data[,L_temp]))
        delta.v <- c(delta.v, rep(delta[i],n_level))
      } else{
        delta.v <- c(delta.v, delta[i])
      }
    }
  }
  
  prob <- list(sense="min")
  prob$c <- c(-phi.b.hat, 1)
  prob$bx <- rbind(blx=c(rep(-Inf,k),0), bux=c(rep(Inf,k+1)))
  
  # using whether has "|" to choose the conditional SP or Global SP
  for (fair in fairness) {
    
    if (grepl("\\|", fair)){ 
      # this is the conditional SP
      # extract the name of legitimate variable
      L_temp <- sub(".*\\|", "", fair)
      # extract the dummies legitimate variable names
      nm.l.factors <- colnames(L.mat)[grep(paste0("^",L_temp), colnames(L.mat))]
      
      # extract the sensitive variable
      s.nm <- sub(".*~(.*?)\\|.*", "\\1", fair)
      for (l.nm in nm.l.factors) {
        
        temp.C <- cbind((dat[,l.nm]*(1-dat[,s.nm])) %*% b.mat/(N.tot*mean(dat[,l.nm]*(1-dat[,s.nm])))
                        - (dat[,l.nm]*dat[,s.nm]) %*% b.mat/(N.tot*mean(dat[,l.nm]*dat[,s.nm])), 0)
        
        C.mat <- rbind(C.mat, temp.C) # C.mat's row counts are n.constraints
      }
    } else { 
      # this is the global statistical parity
      s.nm <- sub(".*~", "", fair)
      # global statistical parity
      temp.C <- cbind((1-dat[,s.nm]) %*% b.mat/(N.tot*mean(1 - dat[,s.nm])) - dat[,s.nm] %*% b.mat/(N.tot*mean(dat[,s.nm])),0)
      C.mat <- rbind(C.mat, temp.C)
    }
  }
  
  prob$A <- Matrix(C.mat)
  
  prob$bc <- rbind(blc=-delta.v,
                   buc=delta.v)
  
  # Cholesky factorization
  Q <- suppressWarnings(chol(Q.mat, TRUE))
  r <- attr(Q, 'rank')
  
  Q[-(1:r), -(1:r)] <- 0
  oo <- order(attr(Q, 'pivot'))
  unpivQ.mat <- Q[, oo]
  
  prob$F <- Matrix(rbind(c(rep(0,k), 1),
                         rep(0,k+1),
                         cbind(unpivQ.mat, as(matrix(0,k,1), "dgCMatrix"))
  ), sparse = TRUE)
  prob$g <- c(0, 1, rep(0,k))
  prob$cones <- matrix(list("RQUAD", k+2, NULL), nrow=3, ncol=1)
  rownames(prob$cones) <- c("type","dim","conepar")
  r <- mosek(prob, list(verbose=0))
  if (!r$sol$itr$solsta=="OPTIMAL") {
    warning("optimal solution not found")
  }
  beta.hat <- r$sol$itr$xx[1:k]
  tau.hat <- b.mat %*% beta.hat
  
  # unfairness 
  message("Fairness restrained CATE calculation finished...")
  message("Extract the unfairness...")
  
  # adding a vector to load all the unfairness
  unfairness <- c()
  
  for (fair in fairness) {
    
    if (grepl("\\|", fair)){ 
      # this is the conditional SP
      # extract the name of legitimate variable
      L_temp <- sub(".*\\|", "", fair)
      # extract the dummies legitimate variable names
      nm.l.factors <- colnames(L.mat)[grep(paste0("^",L_temp), colnames(L.mat))]
      # extract the sensitive variable
      s.nm <- sub(".*~(.*?)\\|.*", "\\1", fair)
      for (l.nm in nm.l.factors) {
        ufair <- abs(mean(tau.hat[which(dat[,s.nm]==1 & dat[,l.nm]==1),]) - mean(tau.hat[which(dat[,s.nm]==0 & dat[,l.nm]==1),]))
        print(paste0("tau~",s.nm,"|",l.nm,"=1", ": ", ufair))
        unfairness <- c(unfairness,ufair)
      }
    } else { 
      # this is the global statistical parity
      s.nm <- sub(".*~", "", fair)
      # global statistical parity
      ufair <- abs(mean(tau.hat[which(dat[,s.nm]==1),]) - mean(tau.hat[which(dat[,s.nm]==0),]))
      print(paste0("tau~",s.nm,": ",ufair))
      unfairness <- c(unfairness,ufair)
    }
  }
  message("...All set!")
  
  cf_outcome <- data.frame(
    A1 = outcome.A1.estimate,
    A0 = outcome.A0.estimate
  )
  
  # PAN: 2025.04.15
  # add the estimated policy
  policy_hat <- ifelse(tau.hat > 0, 1, 0)
  
  return(list(tau_hat = tau.hat,
              unfair = unfairness,
              cf_outcomes = cf_outcome,
              ps_scores = ps.estimate,
              phi.b.hat = phi.b.hat,
              beta.hat = beta.hat,
              policy_hat = policy_hat))
}

###############################################################################
#                                                                             #
#           multilevel fair CATE function                                     #
#           version: V8, 2025.04.27                                           #
#                                                                             #
###############################################################################

fairCATE_multilevel<-function(data,
                              sensitive,
                              legitimate=NULL,
                              fairness,
                              treatment,
                              outcome,
                              cluster,
                              multicategorical=NULL,
                              outcome.LMM,
                              ps.GLMM,
                              cf_outcomes=NULL,
                              ps_scores=NULL,
                              delta,
                              ps.trim="Sturmer.1",
                              tau_var_j = NULL){
  # PAN: 2025.04.20 v8 is based on v7
  # Based on Prof.Kim's latest suggestions on constraints, I revised the C.mat part
  # It should works well for both individual or cluster levels fairness constraints
  # Additionally, the conditional SP should work well too.
  
  
  # PAN: 2025.04.16  v7 is based on V6
  # adding the tau_var_j to the function for a correctly specified projection if
  # we add the random slope into the DGP
  # after predicting the counterfactuals and propensity scores, I augmented the b.mat 
  # to b.mat.aug with tau_var_j column. There will be no influence if it is NULL.
  # You need to specify the tau_var_j like tau_var_j = data$tau_var_j
  
  # PAN: 2025.04.14
  # predict the ps using random effect, edited by Prof.suk
  
  # PAN: 2025.04.06
  # 1) v4 version is based on v3;
  # 2) add two more argument to indicate the sensitives' levels for constructing the uf;
  # 3) modify the uf function to align with proposed cluster adjusted method;
  # 4) revised the counterfactuls prediction using both fixed and random effects. Previous
  #    version uses the fixed effects only.
  
  
  # 2025.04.05.
  # 1) v3 version is based on v2;
  # 2) keep using cluster adjusted Q.mat all the time;
  # 3) predict the counterfatucls using both the fixed effect and random effect;
  # 4) modify the uf function to align with proposed cluster adjusted method;
  
  
  # define the expit function
  expit <- function(x){ exp(x)/(1+exp(x))}
  
  # define the number of clusters
  m <- length(unique(data[,cluster]))
  
  # Handling the multicategorical variables, i.e., variables specified in reference
  df_mc_dummied <- as.data.frame(matrix(1, nrow = nrow(data), ncol = 1))
  
  # Convert them into a dummy coded dataframe 
  for (mc_var in multicategorical) {
    df_mc <- as.data.frame(data[,mc_var])
    names(df_mc) <-  mc_var
    df_mc <- dummy_cols(df_mc, select_columns = c(mc_var),
                        ignore_na = T, remove_selected_columns = T)
    # drop the last level within each variable to avoid Singularity matrix
    df_mc <- df_mc[,-ncol(df_mc)]
    # attach the dummy coded variable into the dataframe
    df_mc_dummied <- cbind(df_mc_dummied, df_mc)
  }
  
  df_mc_dummied <- df_mc_dummied[,-1]
  
  # exclude the non-covariates
  exclude_vars <- c(cluster,multicategorical,"A")
  include_vars <- setdiff(names(data),c("Y", exclude_vars))
  
  # merge the dummied multicategorical covariataes matrix with the other covariates
  b.mat <- as.matrix(cbind(1,cbind(data[include_vars],df_mc_dummied)))
  
  A <- data$A
  C <- data[, cluster]
  
  # cluster dummy matrix
  Cl.mat <- t(do.call(rbind, lapply(1:m, function(i) ifelse(C == i, 1, 0)))) 
  colnames(Cl.mat) <- paste0("C", 1:m)
  
  Y <- data$Y
  
  # create an empty L matrix
  L.mat <- matrix(nrow = dim(data)[1],ncol = 0)
  if (is.null(legitimate)) {
    L.mat <- NULL
    dat <- cbind(Cl.mat, data)
  } else {
    # combine all legitimate matrix into one unified L matrix
    for (i in 1:length(legitimate)) {
      # PAN: support adding the cluster ID as legitimate variable
      if (legitimate[i] == cluster){
        Cl.mat_ <- Cl.mat
        colnames(Cl.mat_) <- paste0(legitimate[i],"_", 1:(ncol(Cl.mat_)))
        L.mat <- cbind(L.mat, Cl.mat_)
      } else{
        L_ <- data[, legitimate[i]]
        L.mat_ <- model.matrix(~as.factor(L_)+0)
        # adding support for multicategorical legitimate 
        colnames(L.mat_) <- paste0(legitimate[i],"_", 0:(nlevels(as.factor(L_))-1)) #### 
        L.mat <- cbind(L.mat,L.mat_)
      }
    }
    dat <- cbind(Cl.mat, data, L.mat)
  }
  
  
  cluster_vec <- unique(data[, cluster])
  ps.var <- as.numeric(VarCorr(ps.GLMM)) # variance of random effect in the propensity score
  
  if (ps.var==0){warning("There is no random effect in the propensity score model!")}
  
  # create the interaction term in b.matrix
  interaction_A1 <- b.mat[, -1] * 1
  interaction_A0 <- b.mat[, -1] * 0
  
  # the old method, consider the fixed effect only
  # outcome.A1.estimate <- as.numeric( cbind(1, 1, b.mat[,-1], interaction_A1)%*%fixef(outcome.LMM) ) # estimate of E(Y|A=1,X,S)
  # outcome.A0.estimate <- as.numeric( cbind(1, 0, b.mat[,-1], interaction_A0)%*%fixef(outcome.LMM) ) # estimate of E(Y|A=0,X,S)
  
  # 2025.04.05: the new version, consider the random effect (customizable)
  # check the fitted object
  newdata_A1 <- as.data.frame(cbind(1, b.mat[,-1], interaction_A1))
  newdata_A1$C <- C
  names(newdata_A1)[1] <- "A" 
  
  newdata_A0 <- as.data.frame(cbind(0, b.mat[,-1], interaction_A0))
  newdata_A0$C <- C
  names(newdata_A0)[1] <- "A"
  
  if (is.null(cf_outcomes)){
    # get individualized predictions that take cluster-specific deviations into account
    outcome.A1.estimate <- predict(outcome.LMM, newdata = newdata_A1, re.form = NULL)
    outcome.A0.estimate <- predict(outcome.LMM, newdata = newdata_A0, re.form = NULL)
  } else {
    # use the user specified counterfactuals
    outcome.A1.estimate <- cf_outcomes$A1
    outcome.A0.estimate <- cf_outcomes$A0
  }
  
  
  # estimate of Pr(A=1|X,S,random effect=u)
  # ps.given.re <- function(x,u){ 
  #   expit( sum(x*as.numeric(fixef(ps.GLMM))) + u )
  # }
  # 
  
  
  # # estimate of Pr(A=1|X,S) = \int Pr(A=1|X,S,random effect=u)*f(u) du
  # ps.marginal <- function(x) {
  #   # adding a condition to skip the integration if the ps.var = 0
  #   if (ps.var == 0) {
  #     return(ps.given.re(x, 0))
  #   } else {
  #     integrate(function(u) {
  #       ps.given.re(x, u) * dnorm(u, mean=0, sd=sqrt(ps.var))
  #     },
  #     lower = -5 * sqrt(ps.var),
  #     upper = 5 * sqrt(ps.var))$value
  #   }
  # }
  # 
  # # estimate of Pr(A=1|X_i,S_i) for all i
  # message("Getting the estimated propensity scores...")
  # ps.estimate <- 
  #   sapply(1:(nrow(data)),function(ii){ps.marginal(b.mat[ii,])}) 
  
  
  ## YS updated ##
  if (is.null(ps_scores)) {
    ps.estimate <-  predict(ps.GLMM, type="response", re.form=NULL)
  } else {
    # use the user specified propensity scores
    ps.estimate <- ps_scores
  }
  
  
  
  # PS trimming - removing the individuals with extreme PS values
  idx.exclude <- NULL
  if (ps.trim == "Sturmer.1") {
    # the common range method ver.1 by St??rmer et al. Am J Epidemiol 2021;190:1659???1670.
    idx.exclude <- c(which((ps.estimate < min(ps.estimate[A==1])) == TRUE),
                     which((max(ps.estimate[A==0]) < ps.estimate) == TRUE)
    )
  } 
  
  if (ps.trim == "Sturmer.2") {
    # the common range method ver.2 by St??rmer et al. Am J Epidemiol 2010;172:843???854.
    idx.exclude <- c(which((ps.estimate < quantile(ps.estimate[dat$A==1], probs = 0.05)) == TRUE),
                     which((quantile(ps.estimate[dat$A==0], probs = 0.95) < ps.estimate) == TRUE)
    )
  }
  
  # PAN: 2025.04.16 
  # augment the b.mat with tau_var_j
  
  
  b.mat.aug <- cbind(b.mat, tau_var_j)
  
  message("Getting the cluster's weights matrix...")
  # in case we want to reweight on a cluster level later on
  phi.b.hat <- matrix(NA,nrow=m, ncol=ncol(b.mat.aug)) 
  for(jj in 1:m){
    # extract the corresponding school id
    cluster_id <- cluster_vec[jj]
    if (is.null(idx.exclude)) {
      phi.b.hat[jj,] <-
        colMeans( (( A*(Y-outcome.A1.estimate) / ps.estimate - 
                       (1-A)*(Y-outcome.A0.estimate) / (1-ps.estimate) + 
                       (outcome.A1.estimate-outcome.A0.estimate) ) * b.mat.aug)[C==cluster_id,] )
    } else {
      mat_sub <- (( A*(Y-outcome.A1.estimate) / ps.estimate - 
                      (1-A)*(Y-outcome.A0.estimate) / (1-ps.estimate) + 
                      (outcome.A1.estimate-outcome.A0.estimate) ) * b.mat.aug)[C==cluster_id & !(1:length(C) %in% idx.exclude),]
      if ( "numeric" %in% class(mat_sub)) {
        # this condition handles the cluster with only one observation
        phi.b.hat[jj,] <- mat_sub}
      else if(all(is.na(mat_sub))==T){
        # this condition handles the cluster whose all observations are with trimmed propensity
        # will give a weight as 0
        phi.b.hat[jj,] <- rep(0,ncol(b.mat))
      } else{
        colmean <- colMeans(mat_sub)
        phi.b.hat[jj,] <- colmean
      }
    }
  }
  phi.b.hat <- colMeans(phi.b.hat)
  
  #########################################################################
  # Optimization                                                          #
  #########################################################################
  
  message("Optimizing using the given fariness condition...")
  
  
  # 2025.03.29: the new version of Q.mat using cluster adjusted
  # 2025.04.05: using the cluster adjusted Q.mat all the time
  Q.mat <- matrix(0, ncol(b.mat.aug), ncol(b.mat.aug))
  
  for (jj in 1:m){ # across all clusters
    # extract the corresponding school id
    cluster_id <- cluster_vec[jj]
    
    # extract the cluster_j's b.mat
    # 2025.05.03: adding drop=FALSE to avoid the matrix drop
    b.mat_j <- b.mat.aug[C==cluster_id, ,drop= FALSE]
    
    # get the cluster_j's Q_mat and add it to the global Q_mat
    Q.mat <- Q.mat + t(b.mat_j) %*% b.mat_j / nrow(b.mat_j)
    
  }
  
  # get the final cluster adjusted Q.mat
  Q.mat <- Q.mat / m
  
  C.mat <- NULL
  
  k <- ncol(b.mat.aug)
  N.tot <- nrow(b.mat.aug)
  
  
  delta.v <- c()
  # based on the input fariness argument to build the delta vector
  if (is.list(delta)){
    delta.v <- unlist(delta)
  } else {
    for (i in 1:length(fairness)) {
      if (grepl("\\|",fairness[i])){
        
        # to make it works on the multicategorical legitimate variables
        L_temp <- sub(".*\\|", "", fairness[i])
        n_level <- nlevels(as.factor(data[,L_temp]))
        delta.v <- c(delta.v, rep(delta[i],n_level))
      } else{
        delta.v <- c(delta.v, delta[i])
      }
    }
  }
  
  prob <- list(sense="min")
  prob$c <- c(-phi.b.hat, 1)
  prob$bx <- rbind(blx=c(rep(-Inf,k),0), bux=c(rep(Inf,k+1)))
  
  
  # using whether has "|" to choose the conditional SP or Global SP
  for (fair in fairness) {
    
    if (grepl("\\|", fair)){ 
      
      # this is the conditional SP
      # extract the name of legitimate variable
      L_temp <- sub(".*\\|", "", fair)
      # extract the dummies legitimate variable names
      nm.l.factors <- colnames(L.mat)[grep(paste0("^",L_temp), colnames(L.mat))]
      
      # extract the sensitive variable
      s.nm <- sub(".*~(.*?)\\|.*", "\\1", fair)
      
      s <- dat[[s.nm]]
      
      for (l.nm in nm.l.factors) {
        
        L  <- dat[[l.nm]]
        clusters <- split(seq_len(nrow(dat)), C)
        
        p0 <- mean(sapply(clusters, function(idx) mean((1-s[idx])*(L[idx]==1))))
        p1 <- mean(sapply(clusters, function(idx) mean(    s[idx] *(L[idx]==1))))
        
        Uf_by_cluster <- t(sapply(clusters, function(idx){
          uf <- ((1-s[idx])*(L[idx]==1))/p0 - (s[idx]*(L[idx]==1))/p1
          colMeans(uf * b.mat.aug[idx, , drop = FALSE])
        }))
        
        C.mat <- rbind(C.mat, c(colMeans(Uf_by_cluster), 0))
        print(C.mat)
      }
    } else { 
      
      # this is the global statistical parity
      s.nm <- sub(".*~", "", fair)
      s <- dat[[s.nm]]
      clusters <- split(seq_len(nrow(dat)), C)   # list of row indices by cluster
      
      ## denominators  P_n(1-S)  and  P_n(S)
      p0 <- mean(sapply(clusters, function(idx) mean(1 - s[idx])))
      p1 <- mean(sapply(clusters, function(idx) mean(s[idx])))
      
      ## per‑cluster contribution to P_n{ u_f W^T }
      Uf_by_cluster <- t(sapply(clusters, function(idx) {
        uf_vec <- (1 - s[idx]) / p0 - s[idx] / p1                # u_f(O_ij)
        colMeans(uf_vec * b.mat.aug[idx, , drop = FALSE])            # (1/n_j)∑ u_f Wᵀ
      }))
      
      C_row <- colMeans(Uf_by_cluster)          # (1/J)∑ cluster means
      C.mat <- rbind(C.mat, c(C_row, 0))        # append, then continue building
    }
  }
  
  prob$A <- Matrix(C.mat)
  
  prob$bc <- rbind(blc=-delta.v,
                   buc=delta.v)
  
  # Cholesky factorization
  Q <- suppressWarnings(chol(Q.mat, TRUE))
  r <- attr(Q, 'rank')
  
  Q[-(1:r), -(1:r)] <- 0
  oo <- order(attr(Q, 'pivot'))
  unpivQ.mat <- Q[, oo]
  
  prob$F <- Matrix(rbind(c(rep(0,k), 1),
                         rep(0,k+1),
                         cbind(unpivQ.mat, as(matrix(0,k,1), "dgCMatrix"))
  ), sparse = TRUE)
  prob$g <- c(0, 1, rep(0,k))
  prob$cones <- matrix(list("RQUAD", k+2, NULL), nrow=3, ncol=1)
  rownames(prob$cones) <- c("type","dim","conepar")
  r <- mosek(prob, list(verbose=0))
  if (!r$sol$itr$solsta=="OPTIMAL") {
    warning("optimal solution not found")
  }
  beta.hat <- r$sol$itr$xx[1:k]
  tau.hat <- b.mat.aug %*% beta.hat
  
  # unfairness 
  message("Fairness restrained CATE calculation finished...")
  message("Extract the unfairness...")
  
  # adding a vector to load all the unfairness
  unfairness <- c()
  
  for (fair in fairness) {
    
    if (grepl("\\|", fair)){ 
      # this is the conditional SP
      # extract the name of legitimate variable
      L_temp <- sub(".*\\|", "", fair)
      # extract the dummies legitimate variable names
      nm.l.factors <- colnames(L.mat)[grep(paste0("^",L_temp), colnames(L.mat))]
      # extract the sensitive variable
      s.nm <- sub(".*~(.*?)\\|.*", "\\1", fair)
      for (l.nm in nm.l.factors) {
        ufair <- abs(mean(tau.hat[which(dat[,s.nm]==1 & dat[,l.nm]==1),]) - mean(tau.hat[which(dat[,s.nm]==0 & dat[,l.nm]==1),]))
        print(paste0("tau~",s.nm,"|",l.nm,"=1", ": ", ufair))
        unfairness <- c(unfairness,ufair)
      }
    } else { 
      # this is the global statistical parity
      s.nm <- sub(".*~", "", fair)
      # global statistical parity
      ufair <- abs(mean(tau.hat[which(dat[,s.nm]==1),]) - mean(tau.hat[which(dat[,s.nm]==0),]))
      print(paste0("tau~",s.nm,": ",ufair))
      unfairness <- c(unfairness,ufair)
    }
  }
  message("...All set!")
  cf_outcome <- data.frame(
    A1 = outcome.A1.estimate,
    A0 = outcome.A0.estimate
  )
  
  # PAN: 2025.04.15
  # add the estimated policy
  policy_hat <- ifelse(tau.hat > 0, 1, 0)
  
  return(list(tau_hat = tau.hat,
              unfair = unfairness,
              cf_outcomes = cf_outcome,
              ps_scores = ps.estimate,
              phi.b.hat = phi.b.hat,
              beta.hat = beta.hat,
              C.mat = C.mat,
              policy_hat = policy_hat))
}