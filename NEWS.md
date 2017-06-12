# v0.6.4.9005

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
