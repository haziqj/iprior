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
