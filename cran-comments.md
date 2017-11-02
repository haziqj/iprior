## Test environments
* local macOS High Sierra 10.13, R 3.4.2
* Ubuntu 14.04.5 (devel, release and oldrel on travis-ci)
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking installed package size ... NOTE
  installed size is  9.1Mb
  sub-directories of 1Mb or more:
    libs   8.6Mb 
  
  *This NOTE appears only with linux systems (i.e., it does not appear with macOS or win-builder). The size is due to a .so file resulting from the usage of RcppEigen.*

## Downstream dependencies
There are currently no downstream dependencies for this package.
