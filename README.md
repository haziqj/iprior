# R/iprior: An R package for I-prior regression

>**[v0.5](https://github.com/haziqjamil/iprior/releases/tag/v0.5) UPDATED: 28/6/16. Speed bumps and plotting feature.**

Based on the manuscript entitled "Objective Bayes regression using I-priors" by Wicher Bergsma [2016, unpublished]. In a linear regression setting, priors can be assigned to the regression function using a vector space framework, and the posterior estimate of the regression function obtained. I-priors are a class of such priors based on the principle of maximum entropy. 

This package performs linear regression modelling using I-priors in R. It is intuitively designed to be similar to `lm`, with both formula and non-formula based input. The parameters of an I-prior model are the scale parameters of the reproducing kernel Hilbert space over the set of covariates, `lambda`, and the standard deviation of model errors,`sigma`. While the main interest of I-prior modelling is prediction, inference is also possible, e.g. via log-likelihood ratio tests.

For installation instructions and some examples of I-prior modelling, continue reading below. A full list of options, along with some other useful documentation, can be found in the [wiki](https://github.com/haziqjamil/iprior/wiki/). 

## Installation

Install R/iprior from this GitHub repository. To do this, you must first install the [devtools](https://github.com/hadley/devtools) package.

```{r devtool_install, eval=FALSE}
install.packages("devtools")
library(devtools)
```

Then, run the following code to install R/iprior and load the `iprior` library.

```{r iprior_install, eval=FALSE}
install_github("haziqjamil/iprior")
library(iprior)
```

## Example 1: Multiple linear regression

We will be analysing Brownlee's stack loss plant data, which is available in R built-in. For more information and a description of the dataset, consult the help section `?stackloss`.

```{r data1}
str(stackloss)
```

We can fit a multiple regression model on the dataset, regressing `stack.loss` against the other three variables. The I-prior for our function lives in a "straight line" RKHS which we call the *Canonical* RKHS. We fit an I-prior model as follows:

```{r mod1}
mod.iprior <- iprior(stack.loss ~ ., data=stackloss)
```

The `iprior` package estimates the model by an EM algorithm, and by default prints reports for every 100 iterations completed. Several options are available to tweak this by supplying a list of control options (see the [wiki](https://github.com/haziqjamil/iprior/wiki/List-of-options)).

The summary output was designed to look similar to an `lm` output. The only differences are the inclusion of RKHS information, EM convergence report and the final log-likelihood value.

```{r summary_mod1}
summary(mod.iprior)
```

The object `mod.iprior` is of class `iprior` and contains a bunch of things of interest that we can extract. Of interest, among other things, might be fitted values, `fitted(mod.iprior)`, residuals, `residuals(mod.iprior)` and the model coefficients, `coef(mod.iprior)`.

To compare the I-prior model against a regular linear regression model, we could look at the fitted versus residual plot.

```{r plot1, echo=FALSE}
plot(y=residuals(mod.iprior), x=fitted(mod.iprior), pch=19, ylab="Residuals", xlab="Fitted values", main="Fitted vs. Residuals")
mod.lm <- lm(stack.loss~., stackloss)
points(y=residuals(mod.lm), x=fitted(mod.lm),, col=2)
abline(0,0, col="gray", lty=2)
legend("bottomleft", legend=c("iprior", "lm"), pch=c(19,1), col=c(1,2))
```

## Example 2: Multilevel modelling

The High School and Beyond is a national longitudinal survey of of students from public and private high schools in the United States, with information such as students' cognitive and non-cognitive skills, high school experiences, work experiences and future plans collected. Papers such as Raudenbush and Bryk (2002) and Raudenbush et. al. (2004) had analyzed this particular dataset, as mentioned in Rabe-Hesketh and Skrondal (2008).

```{r data2}
data(hsbsmall)
str(hsbsmall)
```

This dataset contains the variables mathach, a measure of mathematics achievement; ses, the socioeconomic status of the students based on parental education, occupation and income; and schoolid, the school identifier for the students. The original dataset contains 160 groups with varying number of observations per group (n=7185 in total). However, this smaller set contains only 16 randomly selected schools, such that the total sample size is n=661. This was mainly done for computational reasons to illustrate this example.

We fit an I-prior model, with the aim of predicting `mathach` from `ses`, with the assumption that students' `ses` varied by `schoolid`. 

```{r mod2}
(mod.iprior <- iprior(mathach ~ ses + schoolid + ses:schoolid, data=hsbsmall))
```

A plot of fitted lines, one for each school, is produced using the `plot()` function. The option `plots="fitted"` produces the plot of interest, but there are other options for this as well.

```{r plot2}
plot(mod.iprior, plots="fitted")
```

## Example 3: One-dimensional smoothing

Instead of just "straight-line" regression functions, we could also smoothed curves. The vector space which contains such curves is called the *Fractional Brownian Motion* RKHS with a Hurst coefficient `Hurst`, which defaults to 0.5. The Hurst coefficient can be thought of as a smoothing parameter, but this is treated as a fixed parameter in the `iprior` package.

Consider a simulated set of points in `datfbm` which were generated from a mixed Gaussian distribution `fx <- function(x) 65*dnorm(x, mean=2)+ 35*dnorm(x,mean=7,sd=1.5)`. 

```{r data3}
data(datfbm)
str(datfbm)
```

To illustrate one-dimensional smoothing, we fit an iprior model and produce the plot of fitted values.

```{r mod3}
mod.iprior <- iprior(y ~ x, data=datfbm, model=list(kernel="FBM"), control=list(silent=T))
summary(mod.iprior)
plot(mod.iprior, plots="fitted")
```

![FBMplot](/images/frontpage/Rplot3.png)

