# date:2024.07.27
# purpose:
#   - this file is to fine tune the v1 version to create large difference!
#   - the v1 is good, but the difference is not large enough




# ::: Data Generating Models :::
create_multileveldata_D1 <- function(cluster_num = 150, # number of clusters
                                        cluster_size = 25,
                                        E_var = 18, # variance of residual
                                        R_var = 0.0001,  # variance of cluster effect in selection
                                        U_var = 0.0001, # variance of cluster effect in outcome
                                        clustereffect=FALSE) {
  
  # ::::: 1) generate the number of observations in each cluster :::::
  
  J <- cluster_num # level-2 unit, the number of cluster
  n.clus <- cluster_size # level-1 unit, cluster sizes
  
  id <- as.factor(rep(1:J, each=n.clus))              # cluster id
  
  # ::::: 2) generate level-2 covariates, W: X21,X22,X23 :::::
  # PAN: cluster level sensitive variable, should be historically disadvantaged minority group
  S2 <- rbinom(J, size = 1, prob = 0.3)
  
  # cluster level covariate 1
  X21 <- rnorm(J, 0, 1)  
  
  # cluster level covariate 2 influnced by S2
  X22 <- rnorm(J, -S2 + 0.5, 1)
  
  # PAN: 2024.07.11
  # cluster level covariate 3 influnced by S2, like the Government's funding for the 
  # protected groups
  
  X23 <- rnorm(J, - 0.3 * S2, 1)
  
  names(X21) <- names(X22) <- names(X23) <- names(S2) <- levels(id) 
  
  # ::::: 3) generate level-1 covariates, Xs with cluster-specific means :::::
  
  totalN <- length(id)
  # PAN: the protected group might be the minority group
  S1 <- rbinom(totalN, size = 1, prob = 0.4)  
  # PAN: 
  # individual level covariate 1
  X11 <- rnorm(totalN, 0, 1)
  
  # individual level covariate 2 influence by S1
  X12 <- rnorm(totalN, -2 * S1 + 1, 1)
  
  # PAN: 2024.07.11
  # individual level covariate 3 influnced by S1
  X13 <- rnorm(totalN, -0.5 * S1, 1)
  
  # legitimate variable, the variable you wanna control for conditional statistical disparity
  L_continuous <- rnorm(totalN, -S1 + 0.25, 1)
  # discretize the L_continuous
  L <- ifelse(L_continuous > 0.1, 1, 0)
  
  pop <- data.frame(id, X11, X12, X13, S1, L, X21=X21[id], X22=X22[id], X23=X23[id], S2=S2[id])
  
  # ::::: 4) generate selection probabilities and potential outcome ::::::::
  
  E <- rnorm(totalN, 0, sqrt(E_var))   # error terms for pot.   
  
  if (clustereffect == FALSE) {
    
    # PAN:2024.07.03
    # I include all the covariates as moderators except the S2
    # This should be same to the Professor's version.
    pop$lps <- -0.7 + 0.3*(pop$X11 + pop$X12 + pop$X13 + pop$S1 + pop$L) + 
      0.3*(pop$X21 + pop$X22 + pop$X23 - pop$S2)
    
    pop$Y0  <- 4 + 0.4*(pop$X11  + pop$X12 + pop$X13 + pop$S1 + pop$X21 - 
                          pop$X22 + pop$X23 - pop$S2 + pop$L) + E
    
    pop$Y1  <- pop$Y0 + 0.8 + 0.3 * pop$X11 + 0.2 * pop$X12 + 0.2 * pop$X13 - 
      0.5*pop$S1 + 0.5*pop$X21 + 0.4*pop$X22 + 0.2*pop$X23 +  0.1*pop$L - 0.2*pop$S2
    
  } else {
    # R_j <- rnorm(totalN, 0, 0.25) # level-2 cluster effect in selection
    # U_j <- rnorm(totalN, 0, 0.5) # level-2 cluster effect in outcome
    # 2024.07.02 
    # adding the cluster-specific random effect and increase the variance
    R_j <- rnorm(J, 0, sqrt(R_var)) # level-2 cluster effect in selection
    U_j <- rnorm(J, 0, sqrt(U_var)) # level-2 cluster effect in outcome
    pop$R_j <- R_j[id]
    pop$U_j <- U_j[id]
    
    pop$lps <- -0.7 + 0.3*(pop$X11 + pop$X12 + pop$X13 + pop$S1 + pop$L) + 
      0.3*(pop$X21 + pop$X22 + pop$X23 - pop$S2)  + pop$R_j
    
    pop$Y0  <- 4 + 0.4*(pop$X11  + pop$X12 + pop$X13 + pop$S1 + pop$X21 - 
                          pop$X22 + pop$X23 - pop$S2 + pop$L) + E + pop$U_j
    
    pop$Y1  <- pop$Y0 + 0.8 + 0.3 * pop$X11 + 0.2 * pop$X12 + 0.2 * pop$X13 - 
      0.5*pop$S1 + 0.5*pop$X21 + 0.4*pop$X22 + 0.2*pop$X23 +  0.1*pop$L - 0.2*pop$S2
  }
  
  
  pop$tau <- pop$Y1 - pop$Y0
  pop$policy <- as.numeric(pop$tau > 0)
  
  pop$ps <- 1 / (1 + exp(-pop$lps))   # propensity scores
  pop$A <- rbinom(totalN, 1, pop$ps) # treatment indicator  
  pop$Y <- with(pop, ifelse(A == 1, Y1, Y0))
  
  return(pop)
}

# ::: Data Generating Models: Design 2 :::
create_multileveldata_D2 <- function(cluster_num = 150, # number of clusters
                                        cluster_size = 25,
                                        E_var = 18, # variance of residual
                                        R_var = 0.0001,  # variance of cluster effect in selection
                                        U_var = 0.0001, # variance of cluster effect in outcome
                                        clustereffect=FALSE) {
  
  # ::::: 1) generate the number of observations in each cluster :::::
  
  J <- cluster_num # level-2 unit, the number of cluster
  n.clus <- cluster_size # level-1 unit, cluster sizes
  
  id <- as.factor(rep(1:J, each=n.clus))              # cluster id
  totalN <- length(id)
  # ::::: 2) generate level-2 covariates, W: X21,X22,X23 :::::
  # PAN: cluster level sensitive variable, should be historically disadvantaged minority group
  S2 <- rbinom(J, size = 1, prob = 0.3)
  # PAN: the protected group might be the minority group
  S1 <- rbinom(totalN, size = 1, prob = 0.4)  
  
  # cluster level covariate 1
  X21 <- rnorm(J, 0, 1)  
  
  # cluster level covariate 2 influnced by S2
  # X22 <- rnorm(J, -S2 + 0.5, 1)
  X22 <- rnorm(J, 0.5, 1)
  
  # PAN: 2024.07.11
  # cluster level covariate 3 influnced by S2, like the Government's funding for the 
  # protected groups
  
  # X23 <- rnorm(J, 0.3 * S2, 1)
  X23 <- rnorm(J, -1, 1)
  
  names(X21) <- names(X22) <- names(X23) <- names(S2) <- levels(id) 
  
  S2 <- S2[id]
  
  # create the intersectional groups
  S_is_0 = ifelse(S1 == 0 & S2 == 0, 1, 0)
  S_is_1 = ifelse(S1 == 1 & S2 == 0, 1, 0)
  S_is_2 = ifelse(S1 == 0 & S2 == 1, 1, 0)
  S_is_3 = ifelse(S1 == 1 & S2 == 1, 1, 0)
  
  # ::::: 3) generate level-1 covariates, Xs with cluster-specific means :::::
  
  # PAN: 
  # individual level covariate 1
  # X11 <- rnorm(totalN, 0, 1)
  X11 <- rnorm(totalN, -S_is_2, 1)
  
  # individual level covariate 2 influence by S1_is_3
  #X12 <- rnorm(totalN, -2 * S_is_3 + 1, 1)
  X12 <- rnorm(totalN, -2 * S_is_3, 1)
  
  # PAN: 2024.07.11
  # individual level covariate 3 influnced by S1
  # X13 <- rnorm(totalN, -0.5 * S_is_1 + 1, 1)
  X13 <- rnorm(totalN, S_is_1, 1)
  
  # legitimate variable, the variable you wanna control for conditional statistical disparity
  L_continuous <- rnorm(totalN, -S_is_1 + 0.25, 1)
  # discretize the L_continuous
  L <- ifelse(L_continuous > 0.1, 1, 0)
  
  pop <- data.frame(id, X11, X12, X13, S1, L, 
                    X21=X21[id], X22=X22[id], X23=X23[id], S2=S2,
                    S_is_0, S_is_1, S_is_2, S_is_3)
  
  
  # ::::: 4) generate selection probabilities and potential outcome ::::::::
  
  E <- rnorm(totalN, 0, sqrt(E_var))   # error terms for pot.   
  
  if (clustereffect == FALSE) {
    
    # PAN:2024.07.03
    # I include all the covariates as moderators except the S2
    # This should be same to the Professor's version.
    # PAN:2024.07.26
    # based on design 1, I only add more difference from the intersectional conditions
    pop$lps <- -0.7 + 0.3*(pop$X11 + pop$X12 + pop$X13 + pop$S1 + pop$L) + 
      0.3*(pop$X21 + pop$X22 + pop$X23 - pop$S2)
    
    # Y0 construction
    pop$Y0 <- 4 +
      0.4 * (pop$X11 + pop$X12 + pop$X13 + pop$X21 - pop$X22 + pop$X23 + pop$L) +
      -1 * pop$S_is_1 +
      0.5 * pop$S_is_2 +
      - 0.5 * pop$S_is_3 + E
    
    # Y1 construction
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
      0.1 * pop$L 
    
  } else {
    # 2024.07.02 
    # adding the cluster-specific random effect and increase the variance
    R_j <- rnorm(J, 0, sqrt(R_var)) # level-2 cluster effect in selection
    U_j <- rnorm(J, 0, sqrt(U_var)) # level-2 cluster effect in outcome
    pop$R_j <- R_j[id]
    pop$U_j <- U_j[id]
    
    # PAN:2024.07.26
    # based on design 1, I only add more difference from the intersectional conditions
    pop$lps <- -0.7 + 0.3*(pop$X11 + pop$X12 + pop$X13 + pop$S1 + pop$L) + 
      0.3*(pop$X21 + pop$X22 + pop$X23 - pop$S2)  + pop$R_j
    
    # Y0 construction
    pop$Y0 <- 4 +
      0.4 * (pop$X11 + pop$X12 + pop$X13 + pop$X21 - pop$X22 + pop$X23 + pop$L) +
      -1 * pop$S_is_1 +
      0.5 * pop$S_is_2 +
      - 0.5 * pop$S_is_3 + pop$U_j + E
    
    # Y1 construction
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
      0.1 * pop$L 
  }
  
  
  pop$tau <- pop$Y1 - pop$Y0
  pop$policy <- as.numeric(pop$tau > 0)
  
  pop$ps <- 1 / (1 + exp(-pop$lps))   # propensity scores
  pop$A <- rbinom(totalN, 1, pop$ps) # treatment indicator  
  pop$Y <- with(pop, ifelse(A == 1, Y1, Y0))
  
  return(pop)
}
