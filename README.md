# FairCATE-Multilevel  

## 0 Table of Contents

1. [Overview](#1-overview)
2. [R files](#2-r-files)  
    2.1 [Preparation](#21-preparation)  
    2.2 [Running the R files](#22-running-the-r-files)  
        2.2.1 [functions.R](#221-functionsr)  
        2.2.2 [simulation_design_1.R and simulation_design_2.R](#222-simulation_design_1r-and-simulation_design_2r)  
        2.2.3 [simulation_results_analysis.R and simulation_results_visualizations.R](#223-simulation_results_analysisr-and-simulation_results_visualizationsr)

## 1 Overview  

**Title**: Fair and Robust Estimation of Heterogeneous Treatment Effects for Optimal Policies in Multilevel Studies

**Authors**: Youmi Suk, Chan Park, Chenguang Pan, Kwangho Kim

Recently, there have been growing efforts in developing fair algorithms for treatment effect estimation and optimal treatment recommendations to mitigate discriminatory biases against disadvantaged groups. While most of this work has primarily focused on addressing discrimination due to individual-level attributes (e.g., race/ethnicity), it overlooks the broader impact of societal structures and cultural norms (e.g., structural racism) beyond the individual level. In this paper, we formalize the concept of multilevel fairness for estimating heterogeneous treatment effects to improve fairness in optimal policies. Specifically, we propose a general framework for the estimation of conditional average treatment effects under multilevel fairness constraints that incorporate individual-level sensitive variables, cluster-level sensitive variables, and their combinations. Using this framework, we analyze the trade-off between fairness and the maximum achievable utility by the optimal policy. We evaluate the effectiveness of our framework through a simulation study and a real data study on math course-taking plans using data from the High School Longitudinal Study of 2009.

For more details of our proposed methods, see [our paper](https://osf.io/preprints/psyarxiv/xz3jw). 
  
## 2 R files  
  
### 2.1 Preparation  
Before running the R files, please be sure all required R packages are installed successfully.  

```r
install.packages("Rmosek")
install.packages("MASS")
install.packages("matrixcalc")
install.packages("mbend")
install.packages("Matrix")
install.packages("ggplot2")
install.packages("lme4")
install.packages("bartCause")
install.packages("grf")
install.packages("parallel")
install.packages("dplyr")
install.packages("fastDummies")
```

Please note that installing `Rmosek` package is a little bit different from the traditional way. Please check [this link](https://docs.mosek.com/latest/rmosek/install-interface.html) for instructions and [this link](https://www.mosek.com/products/academic-licenses/) for requesting a free license.  
  
### 2.2 Running the R files  

#### 2.2.1 `functions.R`  
This file contains four functions for GLMM model fitting, estimating fair CATE on multilevel data, and simulated data generating:  

**`GLMM_model()`**   

**Description**  
to fit the outcome model and treatment model using GLMM method.

  
**Usage**
```r
GLMM_model <- function(data,
                       outcome,
                       treatment,
                       cluster,
                       multicategorical=NULL,
                       n_AGQ = 1,
                       fixed_intercept = TRUE,
                       glmer_Control=NULL)
```

**Arguments**  

- `data`: a dataframe. You do not need to convert the multicategorical covariates into dummy (one-hot) codings before model fitting. We recommend users to normalize all the continuous covariates in advance for faster convergence in GLMM model fitting.  
- `outcome`: a string to indicate the name of the outcome variable.  
- `treatment`: a string to indicate the name of the treatment variable. Treatment must be a binary variable for the current version.   
- `cluster`: a string to indicate the name of the cluster id.  
- `multicategroical`: a string array to specify all the multicategorical variables. Use `NULL` by default to indicate that there is no multicategorical variable in your dataframe. If specified, this function will automatically convert these multicategorical variables into one-hot codings. The multicategorical variable(s) must be in factor format.  
- `nAGQ`: integer scalar, same to the argument in `glmer()` function from the `lme4` packages. Values greater than 1 produce greater accuracy in the evaluation of the log-likelihood at the expense of speed.  
- `fixed_intercept`: logical. Whether to include a fixed grand intercept in the outcome model. Using `TRUE` by default. This may depend on your research setting.  
- `gmler_Control`: a list of correct class, resulting from `glmerControl()` containing control parameters. Same as the argument `control` from the `lme4` package. It is used to solve the non-convergence issue in fitting the GLMM model.  

**Example**  
```r
# customize the optimizer and the maximal estimation iterations
glmer_Control <- glmerControl(optimizer = "bobyqa",
                              optCtrl = list(maxfun=100000))
  
# model fitting
glmm_out <- GLMM_model(data = dat0, outcome = "Y",
                         treatment = "A",
                         cluster = "id",
                         n_AGQ = 2,
                         fixed_intercept = FALSE,
                         glmer_Control = glmer_Control)

# check the fitted objects
summary(glmm_out$outcome.LMM)
summary(glmm_out$ps.GLMM)
```

**`fairCATE_multilevel()`**  

**Description:**  
to estimate the CATE with fairness constraints. It will return:  
- `tau.hat`: the estimated CATEs.
- `unfair`: the unfairness on each condition you specified in the argument `fairness`.
- `cf_outcomes`: the estimated counterfactual outcomes from the fitted outcome model.
- `ps_scores`: the estimated propensity scores from the fitted treatment model.

  
**Usage**
```r
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
                                ps.trim="Sturmer.1")
```

  
**Argument**  
- `data`: a dataframe. You do not need to convert the multicategorical covariates into dummy (one-hot) codings before model fitting. We recommend users to normalize all the continuous covariates in advance for faster convergence in GLMM model fitting.  
- `sensitive`: a string array to indicate the sensitive variable(s).
- `legitimate`: a string array to indicate the legitimate variable(s) for conditional statistical disparity. Using NULL by default for statistical disparity.
- `fairness`: a formula style string array to indicate the fairness constraints. It should start with `tau ~`, where we use `tau` to represent the CATE. For example, suppose your sensitive variable is `S1` and the legitimate variable is `L1`, your fairness constraint (conditional statistical disparity) should be written as `tau ~ S1|L1`, which means "get the CATE with fairness constraint on S1 after controlling for L1". The spaces between any letters or symbols in this `fairness` specification do not matter.
- `treatment`: a string to indicate the name of the treatment variable. Treatment must be a binary variable for the current version.
- `outcome`: a string to indicate the name of the outcome variable.
- `cluster`: a string to indicate the name of the cluster id.
- `multicategroical`: a string array to specify all the multicategorical variables. Use `NULL` by default to indicate that there is no multicategorical variable in your dataframe. If specified, this function will automatically convert these multicategorical variables into one-hot codings. The multicategorical variable(s) must be in factor format.
- `outcome.LMM`: a fitted outcome model returned from the function `GLMM_model`.
- `ps.GLMM`: a fitted treatment model returned from the function `GLMM_model`.
- `fixed_intercept`: logical. Whether to include a fixed grand intercept in the outcome model. Using `TRUE` by default. This may depend on your research setting.
- `delta`: a numeric array or a list to indicate the (un)fairness tolerance level $\delta$. Usually, you can use 20 to indicate no fairness constraint (i.e., $\delta = \infty$) and 0.0001 for the most strict fairness constraint (i.e., $\delta = 0$). For example, if you set the `fairness` argument to be `c("tau~S1", "tau~S2")` and you want to give fairness constraints on both two conditions, your `delta` should be `delta = c(0.0001, 0.0001)`.
- `ps.trim`: a string to choose the trimming method for propensity scores. It should be either "Sturmer.1" (the default) or "Sturmer.2". "Sturmer.1" is the common range method ver.1 by [Stürmer et al. Am J Epidemiol 2021;190:1659–1670](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8327194/). "Sturmer.2" is the common range method ver.2 by [Stürmer et al. Am J Epidemiol 2010;172:843–854](https://pubmed.ncbi.nlm.nih.gov/20716704/).

**Example**  
```r
# customize the optimizer and the maximal estimation iterations
glmer_Control <- glmerControl(optimizer = "bobyqa",
                              optCtrl = list(maxfun=100000))
  
# model fitting
glmm_out <- GLMM_model(data = dat0, outcome = "Y",
                         treatment = "A",
                         cluster = "id",
                         n_AGQ = 2,
                         fixed_intercept = FALSE,
                         glmer_Control = glmer_Control)

# estimate CATEs with fairness constraints
ml_fr_out <- fairCATE_multilevel(data = dat0,
                                 sensitive = c("S1","S2"),
                                 legitimate = c("L1"),
                                 fairness = c("tau~S1|L1","tau~S2"),
                                 treatment = "A",
                                 outcome = "Y",
                                 cluster = "id",
                                 multicategorical = NULL,
                                 outcome.LMM = glmm_out$outcome.LMM,
                                 ps.GLMM = glmm_out$ps.GLMM,
                                 fixed_intercept = FALSE,
                                 delta = c(20,20),
                                 ps.trim="Sturmer.1")

# retrieve the estimated CATEs
cates_est <- ml_fr_out$tau_hat

# retrieve the details of unfairness
unfair_detail <- ml_fr_out$unfair

# retrieve the estimated counterfactual outcomes and propensity scores
cf_out <- ml_fr_out$cf_outcomes
ps_out <- ml_fr_out$ps_scores

```


**`create_multileveldata_D1()`**  

**Description**  
  
To generate the multilevel (individual-level and cluster-level) data. If setting the argument `clustereffect=TRUE`, it will return a dataframe with:  
- indivdiaul-level covariates: `X11`, `X12`, `X13`, and `X14`.
- individual-level sensitive variable: `S1`.
- cluster-level covariates: `X21`, `X22`, and `X23`.
- cluster-level sensitive variable: `S2`.
- cluster ID:`id`.
- true CATEs: `tau`.
- true optimal treatment regime: `policy`.
- propensity scores: `ps`.
- observed treatment assignments: `A`.
- the observed outcome: `Y`.


**Usage**  
```r
create_multileveldata_D1 <- function(cluster_num = 150,
                                     cluster_size = 25,
                                     E_var = 18, 
                                     R_var = 0.0001,
                                     U_var = 0.0001, 
                                     clustereffect=FALSE)
```

**Argument**  
- `cluster_num`: an integer to define the number of clusters.
- `cluster_size`: an integer to define the cluster size.
- `E_var`: a float to define the variance of residual in the outcome model.
- `R_var`: a float to define the variance of cluster effect in treatment.
- `U_var`: a float to define the variance of cluster effect in outcome.
- `clustereffect`: logical. To determine whether to generate multilevel data.

**Example**  
```r
# Generate 7500 observations nested in 300 clusters, with the conditional ICC of .105
# for the outcome model, and .372 for the treatment model.
df <- create_multileveldata_D1(cluster_num = 300, 
                                cluster_size = 25, 
                                E_var = 0.6653, # var of residual
                                R_var = 1.95, # var of cluster effect in selection
                                U_var = 0.0776, # var of cluster effect in outcome
                                clustereffect=TRUE)

# retrieve the true heterogeneous treatment effect
tau <- df$tau

# retrieve the true ATE
ATE_true <- mean(df$tau)

# retrieve the true optimal regimes
OTRs_true <- df$policy

# retrieve the true Value under OTRs
value_true <- mean(ifelse(df$policy > 0, df$Y1, df$Y0))
```

**`create_multileveldata_D2()`**  

**Description**  
Based on the function `create_multileveldata_D1`, to generate multilevel data with intersectionality. If setting the argument `clustereffect=TRUE`, it will return a dataframe with:  
- indivdiaul-level covariates: `X11`, `X12`, `X13`, and `X14`.
- individual-level sensitive variable: `S1`.
- cluster-level covariates: `X21`, `X22`, and `X23`.
- cluster-level sensitive variable: `S2`.
- intersectional sensitive variable: `S_is_0`, `S_is_1`, `S_is_2`, and `S_is_3` for the intersectional groups $(S1=0, S2=0)$, $(S1=1, S2=0)$, $(S1=0, S2=1)$, and $(S1=1, S2=1)$, respectively.
- cluster ID:`id`.
- true CATEs: `tau`.
- true optimal treatment regime: `policy`.
- propensity scores: `ps`.
- observed treatment assignments: `A`.
- the observed outcome: `Y`.

**Usage**  
```r
create_multileveldata_D2 <- function(cluster_num = 150,
                                     cluster_size = 25,
                                     R_var = 18, 
                                     V_var = 0.0001,
                                     U_var = 0.0001, 
                                     clustereffect=FALSE)
```

**Argument**  
- `cluster_num`: an integer to define the number of clusters.
- `cluster_size`: an integer to define the cluster size.
- `E_var`: a float to define the variance of residual in the outcome model.
- `R_var`: a float to define the variance of cluster effect in treatment.
- `U_var`: a float to define the variance of cluster effect in outcome.
- `clustereffect`: logical. To determine whether to generate multilevel data.

**Example**  
```r
# Generate 7500 observations nested in 300 clusters, with the conditional ICC of .105
# for the outcome model, and .372 for the treatment model.
df <- create_multileveldata_D2(cluster_num = 300, 
                                    cluster_size = 25, 
                                    R_var = 0.6653, # var of residual
                                    V_var = 1.95, # var of cluster effect in selection
                                    U_var = 0.0776, # var of cluster effect in outcome
                                    clustereffect=TRUE)

# retrieve the true heterogeneous treatment effect
tau <- df$tau

# retrieve the true ATE
ATE_true <- mean(df$tau)

# retrieve the true optimal regimes
OTRs_true <- df$policy

# retrieve the true Value under OTRs
value_true <- mean(ifelse(df$policy > 0, df$Y1, df$Y0))

# subset the generated dataframe for model building
dat0 <- df[,c("id","X11","X12","X13","S1","X14","X21","X22",
              "X23","S2","A","Y","S_is_1","S_is_2","S_is_3")]

# mutate the dat0 with intersectional sensitive variables
library(dplyr)
dat0_is <- dat0 %>% select(-S1, -S2)
```
#### 2.2.2 `simulation_design_1.R` and `simulation_design_2.R`  

These two R files use parallel computing techniques to accelerate the simulation study process. Running either of them for 500 repetitions will be finished in about 2 hours on our platform, which is as follows:  
- OS: Windows 11 Pro
- CPU: i7 14700F with PL1 250W, PL2 250W
- RAM: 96 GB
- R version: 4.4.1 (2024-06-14 ucrt) -- "Race for Your Life"

A large RAM (at least 64 GB) is recommended. R will borrow space from your hard drive if RAM is full, which costs much more running time like 5 hours on 32GB RAM compared to 2 hours on 96GB RAM.  

These two scripts run well on Apple's Mac systems.

Please scroll down to the last few lines of code and replace `20` with the repetition number you'd like in this line of code:  
```r
results <- parLapply(cl, 1:20, run_iteration)
```

For details about our simulation design, please refer to section **5 Simulation Study** in [our paper](https://osf.io/preprints/psyarxiv/xz3jw).  

#### 2.2.3 `simulation_results_analysis.R` and `simulation_results_visualizations.R`  

After running the simulation scripts in section 2.2.2, please run these two files for results analysis and visualizations. `simulation_results_analysis.R` will return the results similar to `Table 1` and `Table 2` in [our paper](https://osf.io/preprints/psyarxiv/xz3jw). And, `simulation_results_visualizations.R` will return the plots similar to `Figure 1`, `Figure 2`, and `Figure S1` in  [our paper](https://osf.io/preprints/psyarxiv/xz3jw).


