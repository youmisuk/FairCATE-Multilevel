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

**GLMM_model()** 
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

