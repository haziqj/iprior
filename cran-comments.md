## Fixes
This version addresses failed build under R-devel with error message: `In the numeric method for all.equal(), a numeric scale argument is now checked to be positive and allowed to be of length > 1`.

## Test environments
* local macOS Mojave 10.14, R 3.5.1
* macOS High Sierra 10.13.3 on travis-ci, R 3.4.4
* Ubuntu 14.04.5 LTS on travis-ci, R devel, release and oldrel (with valgrind)
* Windows Server 2012 R2 x64 build 9600 on appveyor-ci, R devel, release and oldrel 
* win-builder (devel)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTEs:

* checking installed package size ... NOTE
  installed size is  10.1Mb
  sub-directories of 1Mb or more:
    libs   8.7Mb
  
  *This NOTE appears only with linux systems (i.e., it does not appear with macOS or win-builder). The size is due to a .so file resulting from the usage of RcppEigen.*
  
* Uses the superseded package: 'doSNOW'

  *Use of the 'doSNOW' package as opposed to the 'doParallel' package is required due to the support of the printed txtProgressBar in the 'doSNOW' package.*

## Downstream dependencies
There are currently no downstream dependencies for this package.
