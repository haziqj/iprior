# R/iprior: An R package for I-prior regression

>**[v0.4.7](https://github.com/haziqjamil/iprior/releases/tag/v0.4.7) NEW: Fractional Brownian Motion kernel**

>**WARNING: The I-prior package is currently not optimised for large datasets yet. You might encounter debilitating slowness for `n>1000`. This is mainly due to the matrix multiplication and data storing process when the EM initialises. See issue [#20](https://github.com/haziqjamil/iprior/issues/20).**

> *For installation instructions and a basic functionality tutorial, read on below. Have a look at the [wiki](https://github.com/haziqjamil/iprior/wiki/) for further guidance on topics such as the full list of modelling options, modelling interaction effects, using the predicted log-likelihood feature, and the `update()` function.*

Based on manuscript entitled "Regression modelling with I-priors" by Wicher Bergsma [2014, unpublished]. This package performs linear regression modelling like `lm`. It is formula based, but also takes vectors and matrices. It enjoys all the methods of `lm` like `summary`, `coef`, `predict`, and so on.

Currently, either a single scale parameter or several individual scale parameters for each covariate is supported. This is done by calling the `iprior()` function with option `one.lam=T` and `one.lam=F` (default) respectively. Future updates will hopefully see finer control of the scale parameter for each predictor. A full list of options can be found in the [wiki](https://github.com/haziqjamil/iprior/wiki/List-of-options).

`iprior()` function now recognises categorical covariates in the data frame (as factors). This enables predictors which are categorical, ordinal, or even multi-level modelling (categorical group indicators). As of [v0.4](https://github.com/haziqjamil/iprior/releases/tag/v0.4), the programme supports modelling **interaction effects**. For a tutorial on this, have a look at [this page](https://github.com/haziqjamil/iprior/wiki/Modelling-interactions). 

By default, the Canonical kernel (straight lines) is used for continuous variables, but this can be changed to the [FBM kernel](https://github.com/haziqjamil/iprior/releases/tag/v0.4.6) for smoothing (see below, or the [wiki](https://github.com/haziqjamil/iprior/wiki/List-of-options)). Also, a [small writeup](https://github.com/haziqjamil/iprior/wiki/Using-FBM-kernel-for-regression-with-functional-covariates) on using the FBM kernel for the Tecator dataset.

## Installation
Install R/iprior from this GitHub repository. To do this, you must first install the [devtools](https://github.com/hadley/devtools) package.

```r
install.packages("devtools")
library(devtools)
```

Then install R/iprior

```r
install_github("haziqjamil/iprior")
library(iprior)
```

## Example use 1
We will be analysing Brownlee's stack loss plant data, which is available in R built-in. For more information and a description of the dataset, consult the help section `?stackloss`.

```r
summary(stackloss)
```

An I-prior model is fitted as follows
```r
mod.iprior <- iprior(stack.loss ~ ., data=stackloss)
mod.iprior
```

By default, `iprior` prints report iterations for every 100 iterations completed. You can change this to silent mode by adding `silent=T` to the call. Other default options are `delt=1e-03`, `report.int=100`, and `maxit=50000`.

The summary was designed to look like the summary for `lm`. Parameter estimates are given, along with their standard errors Wald's statistic and corresponding p-values. Note that this may be improved in the future (maybe a log-likelihood ratio test?), especially for the parameter `psi` which is obviously positive.
```r
summary(mod.iprior)
```

The list `mod.iprior` contains a bunch of useful things which we can extract.
```r
ls(mod.iprior)
```

Of interest might be the fitted values
```r
fit.iprior <- fitted(mod.iprior)
```
Or even predicted values for a new data point
```r
newdata <- data.frame(Air.Flow=72, Water.Temp=20, Acid.Conc.=85)
predict(mod.iprior, newdata)
```

We can also compare the I-prior fit to a normal OLS fit.
```r
mod.lm <- lm(stack.loss ~ ., data=stackloss)
fit.lm <- fitted(mod.lm)
plot(fit.lm, fit.iprior, type="n", xlab="OLS predicted values", ylab="I-prior fitted values", main="Comparison between I-prior and classical regression predicted values")
text(fit.lm, fit.iprior, pch=as.character(1:length(fit.lm)), col=1:length(fit.lm), cex=0.7)
abline(a=0, b=1)
```

A comparison of the value of `psi = 1/sigma^2`
```r
sigma.iprior <- 1/sqrt(mod.iprior$psi)
sigma.lm <- summary(mod.lm)$sigma
sigma.iprior; sigma.lm
```

## Example use 2
Here is an example of the Pearson kernel in action. The data set is the classic anatomical data of domestic cats. We would like to predict heart weight from body weight and sex. The dataframe `cats` contains the variable `Sex` which in this instance is categorical. In R, this is recognised as being a factor.

```r
data(cats, package="MASS")
summary(cats)
is.factor(cats$Sex)
```

The `iprior` function recognises categorical variables automatically and applies a Pearson kernel accordingly.
```r
mod.iprior <- iprior(Hwt ~ Bwt + Sex, data=cats)
mod.iprior
```

The summary now gives more information regarding the I-prior model, including which RKHS were used, the final log-likelihood value, whether or not the EM converged and how many iterations were performed.
```r
summary(mod.iprior)
```

## Example use 3
We can also use the Fractional Brownian Motion (FBM) kernel to perform 1-dimensional smoothing. The package contains a simulated dataset of points called `datfbm` which were generated from a mixed Gaussian distribution function `fx <- function(x) 65*dnorm(x, mean=2)+ 35*dnorm(x,mean=7,sd=1.5)`.

```r
data(datfbm)
attach(datfbm)
plot(x, y, cex=0.5)
lines(x, fx(x), type="l", col=3, lty=2)
```

To use the FBM kernel, we include it in the call options:

```r
mod.iprior <- iprior(y~x, datfbm, kernel="FBM")
summary(mod.iprior)
```

Currently, the Hurst coefficient is not estimated in the EM procedure (it is treated as fixed). We can change the value of the Hurst coefficient by including `gam=<Hurst.value>` in the call option, where `Hurst.value` is between 0 and 1, e.g. `gam=0.5`. By default, if no Hurst coefficient is specified, `gam=0.5` is used.

Here's a plot of the fitted values:

```r
yhat <- fitted(mod.iprior)
lines(x, yhat, type="l", col=2)
detach(datfbm)
```

![FBMplot](/images/Rplot7.jpg)

