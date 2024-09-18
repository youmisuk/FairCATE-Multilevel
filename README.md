# FairCATE-Multilevel  
  
## 1.0 Overview  
To be updated...  
  
## 2.0 R files  
  
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

Please note that installing `Rmosek` package is a little bit different from the traditional way. Please check [this link](https://docs.mosek.com/latest/rmosek/install-interface.html) for instruction and [this link](https://www.mosek.com/products/academic-licenses/) for requesting a free license.  
  
### 2.2 Running the R files  

#### 2.2.1 `functions.R`  
This file contains four functions for GLMM model fitting, estimating fair CATE on multilevel data, and simulated data generating:  

**`GLMM_model()`**   

**Description:** fitting the outcome model and treatment model using GLMM method.  

  
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

- `data`: a dataframe. You do not need to convert the multicategorical covariates into dummy (one-hot) codings before model fiiting. We recommend users to normalize all the continuous covariates in advance for a faster convergence in GLMM model fitting.  
- `outcome`: a string to indicate the name of outcome vairable.  
- `treatment`: a string to indicate the name of treatment variable. Treatment must be a binary variable for the current version.   
- `cluster`: a string to indicate the name of the cluster id.  
- `multicategroical`: a string array to specify all the multicategorical variables. Using `NULL` by default to indicate that there is no multicategorical variable in your dataframe. If specified, this function will automatically convert these multicategorical variables into one-hot codings. The multicategorical variable(s) must be in factor format.  
- `nAGQ`: integer scalar, same to the argument in `glmer()` function from the `lme4` packages. Values greater than 1 produce than 1 produce greater accuracy in the evaluation of the log-likelihood at the expense of speed.  
- `fixed_intercept`: logical. Whether to include a fixed grand intercept in the outcome model. Using `TRUE` by default. This may depend on your research setting.  
- `gmler_Control`: a list of correct class, resulting from `glmerControl()` containing control parameters. Same to the argument `control` from the `lme4` package. It is used to solve the non-convergence issue in fitting the GLMM model.  

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

**Description:** to estimate the CATE with fairness contraint.  

  
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
- `data`: a dataframe. You do not need to convert the multicategorical covariates into dummy (one-hot) codings before model fiiting. We recommend users to normalize all the continuous covariates in advance for a faster convergence in GLMM model fitting.  
- `sensitive`: a string array to indicate the sensitive variable(s).
- `legitimate`: a string array to indicate the legitimate variable(s) for conditional stastistical disparity. Using NULL by default for statistical disparity.
- `fairness`: a formula style string array to indicate the fairness constraints. It should start with `tau ~`, where we use `tau` to represent the CATE. For example, supporse your sensitive variable is `S1` and legitimate variable is `L1`, you fairness constraint (conditional statistical disparity) should be written as `tau ~ S1|L1`, which means "get the CATE with fariness constraint on S1 after controlling for L1". The spaces between any letters or symbols in this `fairness` specification does not matter.
- `treatment`: a string to indicate the name of treatment variable. Treatment must be a binary variable for the current version.
- `outcome`: a string to indicate the name of outcome vairable.
- `cluster`: a string to indicate the name of the cluster id.
- `multicategroical`: a string array to specify all the multicategorical variables. Using `NULL` by default to indicate that there is no multicategorical variable in your dataframe. If specified, this function will automatically convert these multicategorical variables into one-hot codings. The multicategorical variable(s) must be in factor format.
- `outcome.LMM`: a fitted outcome model returned from the function `GLMM_model`.
- `ps.GLMM`: a fitted treatment model returned from the function `GLMM_model`.
- `fixed_intercept`: logical. Whether to include a fixed grand intercept in the outcome model. Using `TRUE` by default. This may depend on your research setting.
- `delta`: a numeric array or a list to indicate the (un)fairness tolerance level $\delta$. Usually, you can use 20 to indicate no fairness constraint (i.e., $\delta = \infty$) and 0.0001 for most strict fariness constraint (i.e., $\delta = 0$). For example, if you set `fairness` argument to be `c("tau~S1", "tau~S2")` and you want to give fairness constraints on both two conditions, your `delta` should be `delta = c(0.0001, 0.0001)`.
- `ps.trim`: a string to choose the trimming method for propensity scores. It should be either "Sturmer.1" (the default) or "Sturmer.2". "Sturmer.1" is the common range method ver.1 by Stürmer et al. Am J Epidemiol 2021;190:1659–1670. "Sturmer.2" is the common range method ver.2 by Stürmer et al. Am J Epidemiol 2010;172:843–854.


