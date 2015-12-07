# R/iprior: An R package for I-prior regression

**[v1.0.0 (one-lambda)](https://github.com/haziqjamil/iprior/releases/tag/v1.0.0) Only centred canonical kernel and a single scale parameter able to be used**

Based on manuscript entitled "Regression modelling with I-priors" by Wicher Bergsma [2014, unpublished]. This package performs linear regression modelling like `lm`. It is formula based, but also takes vectors and matrices. It enjoys all the methods of `lm` like `summary`, `coef`, and so on.

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
We will be analysing the classic Fisher cats data, to predict heart weight from body weight.

```r
data(cats, package="MASS")
summary(cats
```

An I-prior model is fitted as follows
```r
mod.iprior <- iprior(Hwt~Bwt, data=cats)
mod.iprior
```

By default, `iprior` prints report iterations for every 1000 iterations completed. You can change this to silent mode by adding `silent=T` to the call. Other default options are `delt=1e-05`, `report.int=1000`, and `maxit=50000`.

The summary was designed to look like the summary for `lm`. Parameter estimates are given, along with their standard errors Wald's statistic and corresponding p-values. Note that this may be improved in the future (maybe a log-likelihood ratio test?), especially for the parameter `psi` which is obviously positive.
```r
summary(mod.iprior)
```

The list `mod.iprior` contains a bunch of useful things which we can extract.
```r
ls(mod.iprior)
```

Of interest might be the fitted values (or even predicted values for new data sets)
```r
fit.iprior <- fitted(mod.iprior)
fitnew.iprior <- predict(mod.iprior, newdata)
```

We can also compare the I-prior fit to a normal OLS fit.
```r
mod.lm <- lm(Hwt~Bwt, data=cats)
fit.lm <- fitted(mod.lm)
plot(fit.lm, fit.iprior, type="n", xlab="Classical regression model estimates", ylab="I-prior estimates", main="Comparison between I-prior and classical regression predicted values")
text(fit.lm, fit.iprior, pch=as.character(1:length(fit.lm)), col=1:length(fit.lm), cex=0.7)
abline(a=0, b=1)
```

A comparison of the value of `psi = 1/sigma^2`
```r
sigma.iprior <-  1/sqrt(mod.iprior$psi)
sigma.lm <- summary(mod.lm)$sigma
sigma.iprior; sigma.lm
```
