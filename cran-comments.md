## Test environments
* local OS X 10.12.5, R 3.4.0
* Ubuntu 12.04.5 (on travis-ci), R 3.4.0
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking installed package size ... NOTE
  installed size is  10.1Mb
  sub-directories of 1Mb or more:
    libs   9.5Mb 
  
  *This NOTE appears only with linux systems (i.e., it does not appear with OS X or win-builder). The size is due to a .so file resulting from the usage of RcppEigen.*

## Downstream dependencies
There are currently no downstream dependencies for this package.
