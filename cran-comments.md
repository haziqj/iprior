## Test environments
* local OS X 10.11.6, R 3.3.1
* ubuntu 12.04.5 (on travis-ci), R 3.3.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* This NOTE appears only with linux systems (i.e., it does not appear with OS X or win-builder):
  
  checking installed package size ... NOTE
  installed size is  8.2Mb
  sub-directories of 1Mb or more:
    libs   7.7Mb 

  The size is due to a .so file resulting from the usage of RcppEigen.
