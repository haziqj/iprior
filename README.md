# R/iprior: An R package for I-prior regression

>**[v1.1.0 (mult-lambda)](https://github.com/haziqjamil/iprior/releases/tag/v1.1.0) Only centred canonical kernel supported**

Based on manuscript entitled "Regression modelling with I-priors" by Wicher Bergsma [2014, unpublished]. This package performs linear regression modelling like `lm`. It is formula based, but also takes vectors and matrices. It enjoys all the methods of `lm` like `summary`, `coef`, and so on.

Currently, either one single scale parameter or individual scale parameters for each covariate is supported. This is done by calling the `iprior` function with option `one.lam=T` and `one.lam=F` (default) respectively. Future updates will hopefully see finer control of the scale parameter for each predictor.

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

## Example use
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
