# R/iprior: An R package for I-prior regression

>**[v1.3.2 (interactions)](https://github.com/haziqjamil/iprior/releases/tag/v1.3.2) WARNING: `iprior()` (mainly the Pearson kernel generation) may be slow for very large datasets (n > 10,000) - see this [comment](https://github.com/haziqjamil/iprior/commit/87f72554a0a35f6ad5a07ab108ff367f1599d251#commitcomment-15054556
). This version fixes an error in the `predict()` functionality.**

> *For installation instructions and basic functionality tutorial, read on below. For a specific tutorial on modelling interactions (v1.3.0 and later), click [here](README-INTERACTIONS.md).*

Based on manuscript entitled "Regression modelling with I-priors" by Wicher Bergsma [2014, unpublished]. This package performs linear regression modelling like `lm`. It is formula based, but also takes vectors and matrices. It enjoys all the methods of `lm` like `summary`, `coef`, and so on.

Currently, either one single scale parameter or individual scale parameters for each covariate is supported. This is done by calling the `iprior()` function with option `one.lam=T` and `one.lam=F` (default) respectively. Future updates will hopefully see finer control of the scale parameter for each predictor.

`iprior()` function now recognises categorical covariates in the data frame (as factors). This enables predictors which are categorical, ordinal, or even multi-level modelling (categorical group indicators). The current version also supports modelling **interaction effects**. For a tutorial on this, have a look at [this page](README-INTERACTIONS.md).

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
