library(MASS)
library(matrixcalc)
library(mbend)
library(Matrix)
library(Rmosek)
library(lme4)
library(fastDummies)

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
    
    # merge the dummied multicategorical covairtaes matrix with the other covariates
    b.mat <- as.matrix(cbind(1,cbind(data[include_vars],df_mc_dummied)))
    
    A <- data$A
    C <- data[, cluster]
    Y <- data$Y
    
    # PAN:2024.07.10
    # revise the fitting process by adding all interaction terms
    # convert b.mat to a data.frame
    b.mat.df <- as.data.frame(b.mat)
    
    # PAN:2024.07.10
    # fit the LMM model with interaction terms
    message("Fitting the LMM model for the outcome with interaction terms....")
    
    # PAN: 2024.07.11
    # to avoid the interaction between the A and C, we need manually specify the lmm formula
    # fixed_effects <- paste("A *",include_vars, collapse = " + ")
    # PAN: 2024.07.12
    # including the inetraction between the A and the dummied multicategorical variables
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
    
    # PAN:2024.07.10
    # revise the fitting process by adding all interaction terms
    # convert b.mat to a data.frame
    b.mat.df <- as.data.frame(b.mat)
    
    # PAN:2024.07.10
    # fit the LMM model with interaction terms
    message("Fitting the LMM model for the outcome with interaction terms....")
    
    # PAN: 2024.07.11
    # to avoid the interaction between the A and C, we need manually specify the lmm formula
    # fixed_effects <- paste("A *",include_vars, collapse = " + ")
    # PAN: 2024.07.12
    # including the inetraction between the A and the dummied multicategorical variables
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
    
    # get the counterfactual outcome
    # PAN: 2024.07.11
    # create the interaction term in b.matrix
    interaction_A1 <- b.mat[, -1] * 1
    interaction_A0 <- b.mat[, -1] * 0
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
      # PAN: 2024.07.03 
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
          # PAN: this condition handles the cluster with only one observation
          phi.b.hat[jj,] <- mat_sub}
        else if(all(is.na(mat_sub))==T){
          # PAN: this condition handles the cluster whose all observations are with trimmed propensity
          # PAN: will give a weight as 0
          phi.b.hat[jj,] <- rep(0,ncol(b.mat))
        } else{
          colmean <- colMeans(mat_sub)
          phi.b.hat[jj,] <- colmean
        }
      }
    }
    phi.b.hat <- colMeans(phi.b.hat)
    
    ####################################################################################
    # Optimization
    ####################################################################################
    
    ####
    # we want to achieve unfairness at S1 for conditional statistical parity and S1 for global statistical parity.
    
    message("Optimizing using the given fariness condition...")
    
    Q.mat <- (t(b.mat) %*% b.mat)/nrow(b.mat)
    C.mat <- NULL
    
    k <- ncol(b.mat)
    N.tot <- nrow(b.mat)
    
    
    delta.v <- c()
    # PAN: based on the input fariness argument to build the delta vector
    if (is.list(delta)){
      delta.v <- unlist(delta)
    } else {
      for (i in 1:length(fairness)) {
        if (grepl("\\|",fairness[i])){
          # 2024.05.29
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
    
    #### PAN revision ####
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
        # this is the global SP
        s.nm <- sub(".*~", "", fair)
        # global statistical parity
        temp.C <- cbind((1-dat[,s.nm]) %*% b.mat/(N.tot*mean(1 - dat[,s.nm])) - dat[,s.nm] %*% b.mat/(N.tot*mean(dat[,s.nm])),0)
        C.mat <- rbind(C.mat, temp.C)
      }
    }
    
    prob$A <- Matrix(C.mat)
    # PAN
    prob$bc <- rbind(blc=-delta.v,
                     buc=delta.v)
    # Cholesky factorization
    Q <- suppressWarnings(chol(Q.mat, TRUE))
    r <- attr(Q, 'rank')
    
    # if (r < nrow(x)) Q[(r+1):nrow(x), (r+1):nrow(x)] <- 0
    Q[-(1:r), -(1:r)] <- 0
    oo <- order(attr(Q, 'pivot'))
    unpivQ.mat <- Q[, oo]
    # all.equal(crossprod(unpivQ.mat), Q.mat)
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
    
    # unfairness #### PAN revision ####
    message("Fairness restrained CATE calculation finished...")
    message("Extract the unfairness...")
    
    # 2024.06.18
    # adding a vector to load all the unfariness
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
    
    # get the counterfactual outcome
    # PAN: 2024.07.11
    # create the interaction term in b.matrix
    interaction_A1 <- b.mat[, -1] * 1
    interaction_A0 <- b.mat[, -1] * 0
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
      # PAN: 2024.07.03 
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
          # PAN: this condition handles the cluster with only one observation
          phi.b.hat[jj,] <- mat_sub}
        else if(all(is.na(mat_sub))==T){
          # PAN: this condition handles the cluster whose all observations are with trimmed propensity
          # PAN: will give a weight as 0
          phi.b.hat[jj,] <- rep(0,ncol(b.mat))
        } else{
          colmean <- colMeans(mat_sub)
          phi.b.hat[jj,] <- colmean
        }
      }
    }
    phi.b.hat <- colMeans(phi.b.hat)
    
    ####################################################################################
    # Optimization
    ####################################################################################
    
    ####
    # we want to achieve unfairness at S1 for conditional statistical parity and S1 for global statistical parity.
    
    message("Optimizing using the given fariness condition...")
    
    Q.mat <- (t(b.mat) %*% b.mat)/nrow(b.mat)
    C.mat <- NULL
    
    k <- ncol(b.mat)
    N.tot <- nrow(b.mat)
    
    
    delta.v <- c()
    # PAN: based on the input fariness argument to build the delta vector
    if (is.list(delta)){
      delta.v <- unlist(delta)
    } else {
      for (i in 1:length(fairness)) {
        if (grepl("\\|",fairness[i])){
          # 2024.05.29
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
    
    #### PAN revision ####
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
        # this is the global SP
        s.nm <- sub(".*~", "", fair)
        # global statistical parity
        temp.C <- cbind((1-dat[,s.nm]) %*% b.mat/(N.tot*mean(1 - dat[,s.nm])) - dat[,s.nm] %*% b.mat/(N.tot*mean(dat[,s.nm])),0)
        C.mat <- rbind(C.mat, temp.C)
      }
    }
    
    prob$A <- Matrix(C.mat)
    # PAN
    prob$bc <- rbind(blc=-delta.v,
                     buc=delta.v)
    # Cholesky factorization
    Q <- suppressWarnings(chol(Q.mat, TRUE))
    r <- attr(Q, 'rank')
    
    # if (r < nrow(x)) Q[(r+1):nrow(x), (r+1):nrow(x)] <- 0
    Q[-(1:r), -(1:r)] <- 0
    oo <- order(attr(Q, 'pivot'))
    unpivQ.mat <- Q[, oo]
    # all.equal(crossprod(unpivQ.mat), Q.mat)
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
    
    # unfairness #### PAN revision ####
    message("Fairness restrained CATE calculation finished...")
    message("Extract the unfairness...")
    
    # 2024.06.18
    # adding a vector to load all the unfariness
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