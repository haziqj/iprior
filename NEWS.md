# v0.7.2
* Updated `get_kern_matrix()` accessor function.
* Fixed bug in delta method conversion of standard errors in direct optimisation.
* Updated references and README.

# v0.7.1

* Modification to centering of SE and polynomial kernels.
* Added option `train.samp` and `test.samp` to `kernL()` and `iprior()` to easily split training and test samples for cross-validation.
* Added a function to perform k-fold cross validation experiments for I-prior models.
* Fixed minor bug in `iprior_em_closed()` which caused lambda to expand together with the number of iterations.
* Fixed incorrect calculation of polynomial kernel.
* Removed all legacy functions.
* Updated vignette.
* Added vignette for cross-validation function.

# v0.7.0

* **This udpate provides a complete redesign of the internals of the package. There are more kernels supported, new estimation methods, and plots are done using the `ggplot2` package.**
* Enhanced the methods and calculations for the linear (canonical) kernel, the fractional Brownian motion kernel, and the Pearson kernel.
* Added support for the squared exponential kernel and the `d`-degree polynomial kernel with offset `c`.
* Newly redesigned kernel loader function `kernL()`, while still keeping support for the legacy `.kernL()` function - although there are plans to phase out this in favour of the new one.
* There is now a `summary` method for `ipriorKernel2` objects. 
* The legacy kernels `Canonical`, `FBM` and `Pearson` are now referred to as `linear`, `fbm` and `pearson`, but there is backward compatability with the old references. 
* `parsm` option for interactions has been removed - it's hardly likely that this is ever useful.
* `rootkern` option for Gaussian process regression has been removed. Should use specialised GPR software for this and keep this package for I-priors only.
* `order` option to specify higher order terms has been removed in favour of polynomial kernels.
* The package now supports the following estimation methods: 
    1. Direct minimisation of the marginal deviance;
    2. EM algorithm (efficient closed-form version and the "regular" version);
    3. Combination of direct and EM methods;
    4. A fixed estimation method to obtain the posterior regression function without estimating any hyperparameters; and
    5. The Nystrom kernel approximation method.
* Parallel restarts is supported via `control = list(restarts = TRUE)`. By default it will use the maximum number of available cores to fit the model in parallel from different random initial values.
* New plot functions added: `plot_fitted()`, `plot_predict()`, and `plot_iter()`.
* Updated documentation throughout.
* New vignette added which gives an overview of regression modelling using I-priors.

# v0.6.5

* Updated documentation.
* Edit FBM kernel. Corrected a mistake. Initially for multivariate `x` then  `H(x) = H1(x[1]) + ... + H_p(x[p])`. This is only true for Canonical kernel. Now correctly applies the FBM kernel using the norm function on each multivariate `x_i`.
* Added support for Gaussian process regression with the currently available kernels.
* Fixed memory leak in FBM kernel function. Also made Canonical kernel function more efficient.
* While linear I-prior models can perform classification tasks, one cannot obtain estimation of probabilities for the classes. This is the motivation behind the [`iprobit`] (https://github.com/haziqjamil/iprobit) package. By using a probit link, the I-prior methodology is extended to categorical responses.
* Most functions written here can be used by I-prior probit models in the `iprobit` package. Added support for categorical response kernel loading.
* Exported some helper functions like `is.ipriorKernel()` and `is.ipriorMod()`.

# v0.6.4

* Fixed "override warning" bug in kernel loader when multiple Hurst coefficients used.
* Updated documentation for `iprior()` and `kernL()`.
* Trimmed down the size of `ipriorMod` objects by not saving `Psql`, `Sl`, `Hlam.mat`, and `VarY.inv`. Although these are no longer stored within an `ipriorMod` object, they can still be retrieved via the functions `Hlam()` and `vary()`.
* Fixed a bug with `ipriorOptim()` or `fbmOptim()` whereby standard errors could not be calculated.
* Added new features to `fbmOptim()`: Ability to specify an interval to search for, and also the maximum number of iterations for the initial EM step.

# v0.6.3

* Changed some code to match JSS paper.
* Commented on the line where Pearson kernels are always used for factor-type variables. Should this always be the case?
* Added control option to set intercept at a fixed value.
* Added (hidden) options for `str()` when printing `ipriorKernel` objects.
* Added  `fbmOptim()` function to find optimum Hurst coefficient for fitting FBM I-prior models.
* Added new way to specify Hurst coefficient using the syntax `kernel = "FBM,<value>"`.
* Wrote vignette manual guide which details how to calculate the matrices required for the closed form estimate of `lambda`.
* Removed the T2 statistic from the `summary()` output for now.

# v0.6.2

* Fix for the installation error (#26) on old R releases (prior to 3.3.0). This error was caused by the generic S3 method `sigma()` not being available from the `stats` package prior to R v3.3.0. 

# v0.6.1

* Several bug fixes and cleanups makes this a CRAN-ready release.

# v0.6

* Added documentation for the package.

# v0.5.1

* Added multi-stage model fitting via `kernL()`.

# v0.5

* Massive improvement to the EM engine which brings about speed improvements.
* Added a plotting feature.
 
# v0.4.7

* Bug fixes.
 
# v0.4.6
 
* Added support for Fractional Brownian Motion kernel (i.e. smoothing models).
 
# v0.4.5
 
* Added the 'predicted log-likelihood feature' in the EM reporting.
* WARNING: The I-prior package is currently not optimised for large datasets yet. You might encounter debilitating slowness for `n > 1000`. This is mainly due to the matrix multiplication and data storing process when the EM initialises. See issue #20.
 
# v0.4.4

* More bug fixes. 
 
# v0.4.3
 
* Fixed an error in the `predict()` functionality.
 
# v0.4.2
 
* Added progress feedback reporting feature for the EM algorithm.

# v0.4.1

* Improved Pearson kernel generation, but still requires tweaking.
 
# v0.4
 
* Added support for Pearson kernels (i.e. regression with categorical variables)

# v0.3

* Major bug fixes.
 
# v0.2

* Multiple scale parameters supported.

# v0.1

* First useful release.
* Only centred canonical kernel and a single scale parameter able to be used.
