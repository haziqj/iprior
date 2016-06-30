# R/iprior: An R package for I-prior regression

[![Build Status](https://travis-ci.org/haziqjamil/iprior.svg?branch=master)](https://travis-ci.org/haziqjamil/iprior)

>**[v0.5.1](https://github.com/haziqjamil/iprior/releases/tag/v0.5.1) UPDATED: 28/6/16. Two-stage fitting now possible with `kernL()`.**

Based on the manuscript entitled "Objective Bayes regression using I-priors" by Wicher Bergsma [2016, unpublished]. In a linear regression setting, priors can be assigned to the regression function using a vector space framework, and the posterior estimate of the regression function obtained. I-priors are a class of such priors based on the principle of maximum entropy. 

This package performs linear regression modelling using I-priors in R. It is intuitively designed to be similar to `lm`, with both formula and non-formula based input. The parameters of an I-prior model are the scale parameters of the reproducing kernel Hilbert space over the set of covariates, `lambda`, and the standard deviation of model errors, `sigma`. While the main interest of I-prior modelling is prediction, inference is also possible, e.g. via log-likelihood ratio tests.

For installation instructions and some examples of I-prior modelling, continue reading below. A full list of options, along with some other useful documentation, can be found in the [wiki](https://github.com/haziqjamil/iprior/wiki/). 

## Installation

R/iprior makes use of several C++ code, so as a prerequisite, you must have a working C++ compiler. To get it:

* On Windows, install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
* On Mac, install Xcode from the app store.
* On Linux, `sudo apt-get install r-base-dev` or similar.

Install R/iprior from this GitHub repository. To do this, you must first install the [devtools](https://github.com/hadley/devtools) package.

```{r devtool_install, eval=FALSE}
install.packages("devtools")
library(devtools)
```

Then, run the following code to install and attach the `iprior` package.

```{r iprior_install, eval=FALSE}
install_github("haziqjamil/iprior")
library(iprior)
```

## Syntax

To fit an I-prior model to `mod` regressing `y` against `x`, where these are contained in the data frame `dat`, the following syntax are equivalent.

```{r syntax, eval=FALSE}
mod <- iprior(y ~ x, data=dat)    #formula based input
mod <- iprior(y=dat$y, x=dat$x)   #non-formula based input
```

The call to `iprior()` can be accompanied by model options in the form of `model=list()`, such as choosing the RKHS, number of scale parameters, and others. Control options for the EM algorithm fit is done through the option `control=list()`. For details, see this [wiki page](https://github.com/haziqjamil/iprior/wiki/List-of-options).

## Example 1: Multiple linear regression

We will be analysing Brownlee's stack loss plant data, which is available in R built-in. For more information and a description of the dataset, consult the help section `?stackloss`.

```{r data1}
str(stackloss)
```
```r
## 'data.frame':    21 obs. of  4 variables:
##  $ Air.Flow  : num  80 80 75 62 62 62 62 62 58 58 ...
##  $ Water.Temp: num  27 27 25 24 22 23 24 24 23 18 ...
##  $ Acid.Conc.: num  89 88 90 87 87 87 93 93 87 80 ...
##  $ stack.loss: num  42 37 37 28 18 18 19 20 15 14 ...
```

We can fit a multiple regression model on the dataset, regressing `stack.loss` against the other three variables. The I-prior for our function lives in a "straight line" RKHS which we call the *Canonical* RKHS. We fit an I-prior model as follows:

```{r mod1}
mod.iprior <- iprior(stack.loss ~ ., data=stackloss)
```
```r
## Iteration 0:    Log-likelihood = -145.52941 .......
## Iteration 100:  Log-likelihood = -56.607013 ........
## Iteration 200:  Log-likelihood = -56.354227 ........
## Iteration 300:  Log-likelihood = -56.347994 .......
## Iteration 384:  Log-likelihood = -56.347910 
## EM complete.
```

The `iprior` package estimates the model by an EM algorithm, and by default prints reports for every 100 iterations completed. Several options are available to tweak this by supplying a list of control options (see the [wiki](https://github.com/haziqjamil/iprior/wiki/List-of-options)).

The summary output was designed to look similar to an `lm` output. The only differences are the inclusion of RKHS information, EM convergence report and the final log-likelihood value.

```{r summary_mod1}
summary(mod.iprior)
```
```r
## 
## Call:
## iprior(formula = stack.loss ~ ., data = stackloss)
## 
## RKHS used:
## Canonical (Air.Flow, Water.Temp, Acid.Conc.) 
## with multiple scale parameters.
## 
## Residuals:
##    Min. 1st Qu.  Median 3rd Qu.    Max. 
## -7.4230 -1.7040 -0.3862  1.8700  5.7680 
## 
##                 Estimate    S.E.      z P[|Z>z|]    
## (Intercept)      17.5238  0.6711 26.113   <2e-16 ***
## lam1.Air.Flow     0.0408  0.0250  1.634    0.102    
## lam2.Water.Temp   0.2227  0.1372  1.623    0.104    
## lam3.Acid.Conc.  -0.0123  0.0104 -1.182    0.237    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## EM converged to within 1e-07 tolerance. No. of iterations: 384
## Standard deviation of errors: 3.075 with S.E.: 0.5046
## T2 statistic: 1.835 on ??? degrees of freedom.
## Log-likelihood value: -56.34791
##
```

The object `mod.iprior` is of class `iprior` and contains a bunch of things of interest that we can extract. Of interest, among other things, might be fitted values, `fitted(mod.iprior)`, residuals, `residuals(mod.iprior)` and the model coefficients, `coef(mod.iprior)`.

To compare the I-prior model against a regular linear regression model, we could look at the fitted versus residual plot.

![eg1plot](/images/frontpage/eg1plot.png)

## Example 2: Multilevel modelling

The High School and Beyond is a national longitudinal survey of of students from public and private high schools in the United States, with information such as students' cognitive and non-cognitive skills, high school experiences, work experiences and future plans collected. Papers such as Raudenbush and Bryk (2002) and Raudenbush et. al. (2004) had analyzed this particular dataset, as mentioned in Rabe-Hesketh and Skrondal (2008).

```{r data2}
data(hsbsmall)
str(hsbsmall)
```
```r
## 'data.frame':    661 obs. of  3 variables:
##  $ mathach : num  16.663 -2.155 0.085 18.804 2.409 ...
##  $ ses     : num  0.322 0.212 0.682 -0.148 -0.468 0.842 0.072 0.332 -0.858 0.902 ...
##  $ schoolid: Factor w/ 16 levels "1374","1433",..: 1 1 1 1 1 1 1 1 1 1 ...
```

This dataset contains the variables mathach, a measure of mathematics achievement; ses, the socioeconomic status of the students based on parental education, occupation and income; and schoolid, the school identifier for the students. The original dataset contains 160 groups with varying number of observations per group (n=7185 in total). However, this smaller set contains only 16 randomly selected schools, such that the total sample size is n=661. This was mainly done for computational reasons to illustrate this example.

We fit an I-prior model, with the aim of predicting `mathach` from `ses`, with the assumption that students' `ses` varied by `schoolid`. 

```{r mod2}
(mod.iprior <- iprior(mathach ~ ses + schoolid + ses:schoolid, data=hsbsmall))
```
```r
## Iteration 0:    Log-likelihood = -3716.6891 .....
## Iteration 69:   Log-likelihood = -2137.7988 
## EM complete.
## 
## Call:
## iprior(formula = mathach ~ ses + schoolid + ses:schoolid, data = hsbsmall)
## 
## RKHS used: Canonical & Pearson, with multiple scale parameters.
## 
## 
## Parameter estimates:
## (Intercept)     lambda1     lambda2         psi 
## 13.68325416  0.41779890  0.13230405  0.02804779
##
```

On a technical note, the vector space for functions over the set of nominal-type variables (such as `schoolid`) is called the *Pearson* RKHS. 

A plot of fitted lines, one for each school, is produced using the `plot()` function. The option `plots="fitted"` produces the plot of interest, but there are other options for this as well.

```{r plot2}
plot(mod.iprior, plots="fitted")
```

![eg2plot](/images/frontpage/eg2plot.png)

## Example 3: One-dimensional smoothing

Instead of just "straight-line" regression functions, we could also smoothed curves. The vector space which contains such curves is called the *Fractional Brownian Motion* RKHS with a Hurst coefficient `Hurst`, which defaults to 0.5. The Hurst coefficient can be thought of as a smoothing parameter, but this is treated as a fixed parameter in the `iprior` package.

Consider a simulated set of points in `datfbm` which were generated from a mixed Gaussian distribution `fx <- function(x) 65*dnorm(x, mean=2) + 35*dnorm(x,mean=7,sd=1.5)`. 

```{r data3}
data(datfbm)
str(datfbm)
```
```r
## 'data.frame':    100 obs. of  2 variables:
##  $ y: num  3.85 8.78 6.75 6.99 4.59 ...
##  $ x: num  -0.168 0.133 0.255 0.412 0.424 ...
```

To illustrate one-dimensional smoothing, we fit an iprior model and produce the plot of fitted values.

```{r mod3}
mod.iprior <- iprior(y ~ x, data=datfbm, model=list(kernel="FBM"), control=list(silent=T))
summary(mod.iprior)
plot(mod.iprior, plots="fitted")
```
```r
## 
## Call:
## iprior(formula = y ~ x, data = datfbm)
## 
## RKHS used:
## Fractional Brownian Motion with Hurst coef. 0.5 (x) 
## with a single scale parameter.
## 
## Residuals:
##    Min. 1st Qu.  Median 3rd Qu.    Max. 
## -3.4880 -1.1460 -0.1998  1.0970  3.1390 
## 
##             Estimate   S.E.      z  P[|Z>z|]    
## (Intercept)   9.9961 0.1694 58.993 < 2.2e-16 ***
## lambda        5.8861 1.2468  4.721 < 2.2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## EM converged to within 1e-07 tolerance. No. of iterations: 500
## Standard deviation of errors: 1.694 with S.E.: 0.1354
## T2 statistic: 16.05 on ??? degrees of freedom.
## Log-likelihood value: -222.8153
##
```

![eg3plot](/images/frontpage/eg3plot.png)
