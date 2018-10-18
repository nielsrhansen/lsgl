## High Dimensional Linear Multiple-Response Regression

Linear multiple-response regression with feature and parameter selection
using sparse group lasso. Suitable for high dimensional problems.

This is the **release** of R package **lsgl version 1.3.7**.

### R-package Overview

This package implements procedures for working with linear
multiple-response regression models using sparse group lasso. This
includes procedures for fitting and cross validating sparse models in a
high dimensional setup. See the [Quick Start (Predict airline ticket
prices for multiple airlines)](quick-start.md) for an example of a
traditional workflow consisting of 1) model selection and assessment
using cross validation, 2) estimation of a final model and 3) using the
selected model for carrying out predictions on new data.

![The multiple lasso estimator and the least squares
estimate](https://raw.github.com/nielsrhansen/lsgl/master/fig1.png)

> Comparison of the multiple lasso estimator and least squares estimate
> on simulated data with 50 samples, 50 features and 25 groups. See the
> lsgl example in the package, i.e. run example(lsgl).

**Package highlights:**

  - Feature and parameter selection
  - Fast coordinate gradient descent algorithm
  - Suitable for high dimensional multiclass classification
  - Support for lasso, group lasso and sparse group lasso
  - Supports custom grouping of features
  - Supports sample weighting
  - Supports individual weighting of the group and parameter penalties

The penalized maximum likelihood estimator for the linear
multiple-response regression model is computed using a coordinate
gradient descent algorithm via the
[sglOptim](https://github.com/nielsrhansen/sglOptim) optimizer. Use of
parallel computing for cross validation and subsampling is supported
through the [foreach](https://cran.r-project.org/package=foreach) and
[doParallel](https://cran.r-project.org/package=doParallel) packages.

### Installation

Get the released version from CRAN:

``` r
install.packages("lsgl")
```

Install the release candidate from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("vincent-dk/sglOptim")
devtools::install_github("vincent-dk/lsgl")
```

Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("vincent-dk/sglOptim", ref = "develop")
devtools::install_github("vincent-dk/lsgl", ref = "develop")
```

### Minimal Example

``` r
library(lsgl)

# Load sone data
data(AirlineTicketPrices)

# Setup 2 parallel units
cl <- makeCluster(2)
registerDoParallel(cl)

# Do 10-fold cross validation on 100 models with increasing complexity, using the 2 parallel units
fit.cv <- lsgl::cv(
  x = X,
  y = Y,
  alpha = 0.5,
  lambda = 0.01,
  use_parallel = TRUE
)
```

    ## 
    ## Running lsgl 10 fold cross validation 
    ## 
    ##  Samples:  Features:  Models:  Groups:  Parameters: 
    ##        337        412        6      412       2.472k

``` r
stopCluster(cl)

# Print information about models
# and cross validation errors
fit.cv
```

    ## 
    ## Call:
    ## lsgl::cv(x = X, y = Y, alpha = 0.5, lambda = 0.01, use_parallel = TRUE)
    ## 
    ## Models:
    ## 
    ##  Index:  Lambda:  Features:  Parameters:  Error: 
    ##        1    1.000        2.9         17.4     132
    ##       20    0.413          4           24     103
    ##       40    0.163       10.5         60.9      78
    ##       60    0.064       14.7         83.8      66
    ##       80    0.025       33.7        167.7      58
    ##      100    0.010       48.6          215      51
    ## 
    ## Best model:
    ## 
    ##  Index:  Lambda:  Features:  Parameters:  Error: 
    ##      100     0.01       48.6          215      51

### Documentation

  - R package documentation
  - [Quick Start (Predict airline ticket prices for multiple
    airlines)](quick-start.md)

### Author

Martin Vincent

### License

GPL (\>=2)
