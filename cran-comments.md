## Fixes
This version addresses `These (still?) have files with vignette-like filenames with no identifiable vignette engine (as spotted by tools::pkgVignettes(check = TRUE)), see below: can you pls fix as necessary?` as instructed by K Hornik.

## Test environments
* local macOS High Sierra 10.13, R 3.4.3
* macOS Sierra 10.12.6 on travis-ci, R 3.3.3 
* Ubuntu 14.04.5 LTS on travis-ci, R devel, release and oldrel (with valgrind)
* Windows Server 2012 R2 x64 build 9600 on appveyor-ci, R devel, release and oldrel 
* win-builder (devel)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTEs:

* checking installed package size ... NOTE
  installed size is  10.3Mb
  sub-directories of 1Mb or more:
    libs   8.6Mb 
  
  *This NOTE appears only with linux systems (i.e., it does not appear with macOS or win-builder). The size is due to a .so file resulting from the usage of RcppEigen.*

* Uses the superseded package: 'doSNOW'

  *Use of the 'doSNOW' package as opposed to the 'doParallel' package is required due to the support of the printed txtProgressBar in the 'doSNOW' package.*
  
## Downstream dependencies
There are currently no downstream dependencies for this package.
