# This file contains four functions:
#     1) GLMM_model(): GLMM model fitting function
#     2) fairCATE_multilevel(): multilevel fair CATE function
#     3) create_multileveldata_D1(): Data Generating function for design 1
#     4) create_multileveldata_D2(): Data Generating function for design 2



library(MASS)
library(matrixcalc)
library(mbend)
library(Matrix)
library(Rmosek)
library(lme4)
library(fastDummies)



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
                       glmer_Control=NULL){
  
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
    formula <- as.formula(paste("Y ~ ", fixed_effects, "+ (1 | C)"))
    
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
    
    # Combine the fixed effects with the random effects
    formula <- as.formula(paste("Y ~ 0 +", fixed_effects, "+ (1 | C)"))
    
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
#           multilevel fair CATE function           #
#                                                   #
#####################################################

fairCATE_multilevel <- function(data,
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
                                delta,
                                ps.trim="Sturmer.1"){
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
  
  if (fixed_intercept){
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
    
    # get the counterfactual outcome    
    design_matrix_A1 <- cbind(b.mat, interaction_A1)
    design_matrix_A0 <- cbind(b.mat, interaction_A0)
    
    outcome.A1.estimate <- as.numeric( cbind(1, 1, b.mat[,-1], interaction_A1)%*%fixef(outcome.LMM) ) # estimate of E(Y|A=1,X,S)
    outcome.A0.estimate <- as.numeric( cbind(1, 0, b.mat[,-1], interaction_A0)%*%fixef(outcome.LMM) ) # estimate of E(Y|A=0,X,S)
    
    # estimate of Pr(A=1|X,S,random effect=u)
    ps.given.re <- function(x,u){ 
      expit( sum(x*as.numeric(fixef(ps.GLMM))) + u )
    }
    
    # estimate of Pr(A=1|X,S) = \int Pr(A=1|X,S,random effect=u)*f(u) du
    ps.marginal <- function(x) {
      # adding a condition to skip the integration if the ps.var = 0
      if (ps.var == 0) {
        return(ps.given.re(x, 0))
      } else {
        integrate(function(u) {
          ps.given.re(x, u) * dnorm(u, mean=0, sd=sqrt(ps.var))
        },
        lower = -5 * sqrt(ps.var),
        upper = 5 * sqrt(ps.var))$value
      }
    }
    
    # estimate of Pr(A=1|X_i,S_i) for all i
    message("Getting the estimated propensity scores...")
    ps.estimate <- 
      sapply(1:(nrow(data)),function(ii){ps.marginal(b.mat[ii,])}) 
    
    # PS trimming - removing the individuals with extreme PS values
    idx.exclude <- NULL
    if (ps.trim == "Sturmer.1") {
      # the common range method ver.1 by Stürmer et al. Am J Epidemiol 2021;190:1659–1670.
      idx.exclude <- c(which((ps.estimate < min(ps.estimate[A==1])) == TRUE),
                       which((max(ps.estimate[A==0]) < ps.estimate) == TRUE)
      )
    } 
    
    if (ps.trim == "Sturmer.2") {
      # the common range method ver.2 by Stürmer et al. Am J Epidemiol 2010;172:843–854.
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
        mat_sub <- (( A*(Y-outcome.A1.estimate) / ps.estimate - 
                        (1-A)*(Y-outcome.A0.estimate) / (1-ps.estimate) + 
                        (outcome.A1.estimate-outcome.A0.estimate) ) * b.mat)[C==cluster_id & !(1:length(C) %in% idx.exclude),]
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
        # support adding the cluster ID as legitimate variable
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
    
    # get the counterfactual outcome
    design_matrix_A1 <- cbind(b.mat, interaction_A1)
    design_matrix_A0 <- cbind(b.mat, interaction_A0)
    
    outcome.A1.estimate <- as.numeric( cbind(1,b.mat[,-1], interaction_A1)%*%fixef(outcome.LMM) ) # estimate of E(Y|A=1,X,S)
    outcome.A0.estimate <- as.numeric( cbind(0,b.mat[,-1], interaction_A0)%*%fixef(outcome.LMM) ) # estimate of E(Y|A=0,X,S)
    
    # estimate of Pr(A=1|X,S,random effect=u)
    ps.given.re <- function(x,u){ 
      expit( sum(x*as.numeric(fixef(ps.GLMM))) + u )
    }
    
    # estimate of Pr(A=1|X,S) = \int Pr(A=1|X,S,random effect=u)*f(u) du
    ps.marginal <- function(x) {
      # adding a condition to skip the integration if the ps.var = 0
      if (ps.var == 0) {
        return(ps.given.re(x, 0))
      } else {
        integrate(function(u) {
          ps.given.re(x, u) * dnorm(u, mean=0, sd=sqrt(ps.var))
        },
        lower = -5 * sqrt(ps.var),
        upper = 5 * sqrt(ps.var))$value
      }
    }
    
    # estimate of Pr(A=1|X_i,S_i) for all i
    message("Getting the estimated propensity scores...")
    ps.estimate <- 
      sapply(1:(nrow(data)),function(ii){ps.marginal(b.mat[ii,])}) 
    
    # PS trimming - removing the individuals with extreme PS values
    idx.exclude <- NULL
    if (ps.trim == "Sturmer.1") {
      # the common range method ver.1 by Stürmer et al. Am J Epidemiol 2021;190:1659–1670.
      idx.exclude <- c(which((ps.estimate < min(ps.estimate[A==1])) == TRUE),
                       which((max(ps.estimate[A==0]) < ps.estimate) == TRUE)
      )
    } 
    if (ps.trim == "Sturmer.2") {
      # the common range method ver.2 by Stürmer et al. Am J Epidemiol 2010;172:843–854.
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
        mat_sub <- (( A*(Y-outcome.A1.estimate) / ps.estimate - 
                        (1-A)*(Y-outcome.A0.estimate) / (1-ps.estimate) + 
                        (outcome.A1.estimate-outcome.A0.estimate) ) * b.mat)[C==cluster_id & !(1:length(C) %in% idx.exclude),]
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
        # this is the global SP
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
  }
  
  return(list(tau_hat = tau.hat,
              unfair = unfairness,
              cf_outcomes = cf_outcome,
              ps_scores = ps.estimate))
}


#####################################################
#                                                   #
#     Data Generating function for design 1         #
#                                                   #
#####################################################

# ::: Data Generating Models :::
create_multileveldata_D1 <- function(cluster_num = 150, # number of clusters
                                     cluster_size = 25,
                                     R_var = 18, # variance of residual
                                     V_var = 0.0001,  # variance of cluster effect in selection
                                     U_var = 0.0001, # variance of cluster effect in outcome
                                     clustereffect=FALSE) {
  
  # ::::: 1) generate the number of observations in each cluster :::::
  
  J <- cluster_num # level-2 unit, the number of cluster
  n.clus <- cluster_size # level-1 unit, cluster sizes
  
  id <- as.factor(rep(1:J, each=n.clus))              # cluster id
  
  # ::::: 2) generate level-2 covariates, W: X21,X22,X23 :::::
  # cluster level sensitive variable, should be historically disadvantaged minority group
  S2 <- rbinom(J, size = 1, prob = 0.3)
  
  # cluster level covariate 1
  X21 <- rnorm(J, 0, 1)  
  
  # cluster level covariate 2 influenced by S2
  X22 <- rnorm(J, -S2 + 0.5, 1)
  
  # cluster level covariate 3 influenced by S2, like the Government's funding for the 
  # protected groups
  
  X23 <- rnorm(J, - 0.3 * S2, 1)
  
  names(X21) <- names(X22) <- names(X23) <- names(S2) <- levels(id) 
  
  # ::::: 3) generate level-1 covariates, Xs with cluster-specific means :::::
  
  totalN <- length(id)
  # the protected group might be the minority group
  S1 <- rbinom(totalN, size = 1, prob = 0.4)  

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
  
  pop <- data.frame(id, X11, X12, X13, S1, X14, X21=X21[id], X22=X22[id], X23=X23[id], S2=S2[id])
  
  # ::::: 4) generate selection probabilities and potential outcome ::::::::
  
  E <- rnorm(totalN, 0, sqrt(R_var))   # error terms for pot.   
  
  if (clustereffect == FALSE) {
    
    # Include all the covariates as moderators except the S2
    pop$lps <- -0.7 + 0.3*(pop$X11 + pop$X12 + pop$X13 + pop$S1 + pop$X14) + 
      0.3*(pop$X21 + pop$X22 + pop$X23 - pop$S2)
    
    pop$Y0  <- 4 + 0.4*(pop$X11  + pop$X12 + pop$X13 + pop$S1 + pop$X21 - 
                          pop$X22 + pop$X23 - pop$S2 + pop$X14) + E
    
    pop$Y1  <- pop$Y0 + 0.8 + 0.3 * pop$X11 + 0.2 * pop$X12 + 0.2 * pop$X13 - 
      0.5*pop$S1 + 0.5*pop$X21 + 0.4*pop$X22 + 0.2*pop$X23 +  0.1*pop$X14 - 0.2*pop$S2
    
  } else {
    # adding the cluster-specific random effect and increase the variance
    R_j <- rnorm(J, 0, sqrt(V_var)) # level-2 cluster effect in selection
    U_j <- rnorm(J, 0, sqrt(U_var)) # level-2 cluster effect in outcome
    pop$R_j <- R_j[id]
    pop$U_j <- U_j[id]
    
    pop$lps <- -0.7 + 0.3*(pop$X11 + pop$X12 + pop$X13 + pop$S1 + pop$X14) + 
      0.3*(pop$X21 + pop$X22 + pop$X23 - pop$S2)  + pop$R_j
    
    pop$Y0  <- 4 + 0.4*(pop$X11  + pop$X12 + pop$X13 + pop$S1 + pop$X21 - 
                          pop$X22 + pop$X23 - pop$S2 + pop$X14) + E + pop$U_j
    
    pop$Y1  <- pop$Y0 + 0.8 + 0.3 * pop$X11 + 0.2 * pop$X12 + 0.2 * pop$X13 - 
      0.5*pop$S1 + 0.5*pop$X21 + 0.4*pop$X22 + 0.2*pop$X23 +  0.1*pop$X14 - 0.2*pop$S2
  }
  
  
  pop$tau <- pop$Y1 - pop$Y0
  pop$policy <- as.numeric(pop$tau > 0)
  
  pop$ps <- 1 / (1 + exp(-pop$lps))   # propensity scores
  pop$A <- rbinom(totalN, 1, pop$ps) # treatment indicator  
  pop$Y <- with(pop, ifelse(A == 1, Y1, Y0))
  
  return(pop)
}


#####################################################
#                                                   #
#     Data Generating function for design 2         #
#                                                   #
#####################################################

# ::: Data Generating Models: Design 2 :::
create_multileveldata_D2 <- function(cluster_num = 150, # number of clusters
                                     cluster_size = 25,
                                     R_var = 18, # variance of residual
                                     V_var = 0.0001,  # variance of cluster effect in selection
                                     U_var = 0.0001, # variance of cluster effect in outcome
                                     clustereffect=FALSE) {
  
  # ::::: 1) generate the number of observations in each cluster :::::
  
  J <- cluster_num # level-2 unit, the number of cluster
  n.clus <- cluster_size # level-1 unit, cluster sizes
  
  id <- as.factor(rep(1:J, each=n.clus))              # cluster id
  totalN <- length(id)
  
  # ::::: 2) generate level-2 covariates, W: X21,X22,X23 :::::
  # cluster level sensitive variable, should be historically disadvantaged minority group
  S2 <- rbinom(J, size = 1, prob = 0.3)
  # the protected group might be the minority group
  S1 <- rbinom(totalN, size = 1, prob = 0.4)  
  
  # cluster level covariate 1
  X21 <- rnorm(J, 0, 1)  
  
  # cluster level covariate 2
  X22 <- rnorm(J, 0.5, 1)
  
  # cluster level covariate 3
  X23 <- rnorm(J, -1, 1)
  
  names(X21) <- names(X22) <- names(X23) <- names(S2) <- levels(id) 
  
  S2 <- S2[id]
  
  # create the intersectional groups
  S_is_0 = ifelse(S1 == 0 & S2 == 0, 1, 0)
  S_is_1 = ifelse(S1 == 1 & S2 == 0, 1, 0)
  S_is_2 = ifelse(S1 == 0 & S2 == 1, 1, 0)
  S_is_3 = ifelse(S1 == 1 & S2 == 1, 1, 0)
  
  # ::::: 3) generate level-1 covariates, Xs with cluster-specific means :::::
  
  # individual level covariate 1
  # X11 <- rnorm(totalN, 0, 1)
  X11 <- rnorm(totalN, -S_is_2, 1)
  
  # individual level covariate 2 influence by S1_is_3
  X12 <- rnorm(totalN, -2 * S_is_3, 1)
  
  # individual level covariate 3 influnced by S1
  X13 <- rnorm(totalN, S_is_1, 1)
  
  # legitimate variable, the variable you wanna control for conditional statistical disparity
  L_continuous <- rnorm(totalN, -S_is_1 + 0.25, 1)
  # discretize the L_continuous
  X14 <- ifelse(L_continuous > 0.1, 1, 0)
  
  pop <- data.frame(id, X11, X12, X13, S1, X14, 
                    X21=X21[id], X22=X22[id], X23=X23[id], S2=S2,
                    S_is_0, S_is_1, S_is_2, S_is_3)
  
  
  # ::::: 4) generate selection probabilities and potential outcome ::::::::
  
  E <- rnorm(totalN, 0, sqrt(R_var))   # error terms for pot.   
  
  if (clustereffect == FALSE) {
    
    # Include all the covariates as moderators except the S2
    # Add more difference from the intersectional conditions
    pop$lps <- -0.7 + 0.3*(pop$X11 + pop$X12 + pop$X13 + pop$S1 + pop$X14) + 
      0.3*(pop$X21 + pop$X22 + pop$X23 - pop$S2)
    
    # model for Y0
    pop$Y0 <- 4 +
      0.4 * (pop$X11 + pop$X12 + pop$X13 + pop$X21 - pop$X22 + pop$X23 + pop$X14) +
      -1 * pop$S_is_1 +
      0.5 * pop$S_is_2 +
      - 0.5 * pop$S_is_3 + E
    
    # model for Y1
    pop$Y1 <- pop$Y0 + 0.8 +
      0.3 * pop$X11 +
      0.2 * pop$X12 +
      0.2 * pop$X13 +
      - 1 * pop$S_is_1 +
      0.5 * pop$S_is_2 +
      - 0.5 * pop$S_is_3 +
      0.5 * pop$X21 +
      0.4 * pop$X22 +
      0.2 * pop$X23 +
      0.1 * pop$X14 
    
  } else {
    # add the cluster-specific random effect and increase the variance
    R_j <- rnorm(J, 0, sqrt(V_var)) # level-2 cluster effect in selection
    U_j <- rnorm(J, 0, sqrt(U_var)) # level-2 cluster effect in outcome
    pop$R_j <- R_j[id]
    pop$U_j <- U_j[id]
    
    # add more difference from the intersectional conditions
    pop$lps <- -0.7 + 0.3*(pop$X11 + pop$X12 + pop$X13 + pop$S1 + pop$X14) + 
      0.3*(pop$X21 + pop$X22 + pop$X23 - pop$S2)  + pop$R_j
    
    # model for Y0
    pop$Y0 <- 4 +
      0.4 * (pop$X11 + pop$X12 + pop$X13 + pop$X21 - pop$X22 + pop$X23 + pop$X14) +
      -1 * pop$S_is_1 +
      0.5 * pop$S_is_2 +
      - 0.5 * pop$S_is_3 + pop$U_j + E
    
    # model for Y1
    pop$Y1 <- pop$Y0 + 0.8 +
      0.3 * pop$X11 +
      0.2 * pop$X12 +
      0.2 * pop$X13 +
      - 1 * pop$S_is_1 +
      0.5 * pop$S_is_2 +
      - 0.5 * pop$S_is_3 +
      0.5 * pop$X21 +
      0.4 * pop$X22 +
      0.2 * pop$X23 +
      0.1 * pop$X14 
  }
  
  
  pop$tau <- pop$Y1 - pop$Y0
  pop$policy <- as.numeric(pop$tau > 0)
  
  pop$ps <- 1 / (1 + exp(-pop$lps))   # propensity scores
  pop$A <- rbinom(totalN, 1, pop$ps) # treatment indicator  
  pop$Y <- with(pop, ifelse(A == 1, Y1, Y0))
  
  return(pop)
}
