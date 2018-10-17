# R/iprior: An R package for I-prior regression

[![Build Status](https://travis-ci.org/haziqj/iprior.svg?branch=master)](https://travis-ci.org/haziqj/iprior)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/haziqj/iprior?branch=master&svg=true)](https://ci.appveyor.com/project/haziqj/iprior)
[![Coverage Status](https://img.shields.io/codecov/c/github/haziqj/iprior/master.svg)](https://codecov.io/gh/haziqj/iprior)
[![CRAN_Status_Badge_version_ago](http://www.r-pkg.org/badges/version-ago/iprior)](https://cran.r-project.org/package=iprior)
[![CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/iprior)](https://cran.r-project.org/package=iprior)

Based on the manuscript entitled "Regression and Classification with I-priors" by Wicher Bergsma (2018, [arXiv:1707.00274](https://arxiv.org/abs/1707.00274)). 
In a general regression setting, priors can be assigned to the regression function in a vector space framework, and the posterior estimate of the regression function obtained. 
An I-prior is defined as Gaussian with some prior mean (usually zero) and covariance kernel proportional to the Fisher information for the regression function.

This package performs regression modelling using I-priors in R. 
It is intuitively designed to be similar to `lm`, making use of familiar syntactical conventions and S3 methods, with both formula and non-formula based input. 
The package estimates these parameters using direct log-likelihood maximisation, the expectation-maximisation (EM) algorithm, or a combination of both.
While the main interest of I-prior modelling is prediction, inference is also possible, e.g. via log-likelihood ratio tests.

For installation instructions and some examples of I-prior modelling, continue reading below. 
The package is documented with help files, and the [vignette](http://phd.haziqj.ml/iprior_paper.pdf) provides an introduction to the concept of I-priors and also to using the package.

## Installation

Install R/iprior either by downloading the latest CRAN release

```r
install.packages("iprior")
library(iprior)
```

or the developmental version from this GitHub repository. R/iprior makes use of several C++ code, so as a prerequisite, you must have a working C++ compiler. To get it:

-   On Windows, install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
-   On Mac, install Xcode from the app store.
-   On Linux, `sudo apt-get install r-base-dev` or similar.

The easiest way to then install from this repo is by using the [devtools](https://github.com/hadley/devtools) package. 
Install this first.

``` r
install.packages("devtools")
```

Then, run the following code to install and attach the `iprior` package.

``` r
devtools::install_github("haziqj/iprior", build_vignettes = TRUE)
library(iprior)
```
[//]: # (*Note: The option `build_vignettes = TRUE` builds the package vignettes for viewing, but takes slightly longer. Set `build_vignettes = FALSE`, or remove this option entirely, to skip building the vignettes.*)

## Syntax

To fit an I-prior model to `mod` regressing `y` against `x`, where these are contained in the data frame `dat`, the following syntax are equivalent.

``` r
mod <- iprior(y ~ x, data = dat)     # formula based input
mod <- iprior(y = dat$y, x = dat$x)  # non-formula based input
```

The call to `iprior()` can be accompanied by several other model specification arguments, including choosing the RKHS, hyperparameters to estimate, estimation method, and others. 
Control options for the estimation method of choice is done through the option `control = list()`. 
Find the full list of options by typing `?iprior` in R.

## Resources

View the package vignette by typing `browseVignettes("iprior")` in R or visiting this [link](http://phd.haziqj.ml/iprior_paper.pdf). 
This package is part of the PhD project entitled "Regression Modelling using priors depending on Fisher information covariance kernels (I-priors)" by Haziq Jamil [[link](http://phd.haziqj.ml)].
