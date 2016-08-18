## Test environments
* local OS X 10.11.6, R 3.3.1
* ubuntu 12.04.5 (on travis-ci), R 3.3.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTEs:

* checking installed package size ... NOTE
  installed size is  8.2Mb
  sub-directories of 1Mb or more:
    libs   7.7Mb 
  
  *This NOTE appears only with linux systems (i.e., it does not appear with OS X or win-builder). The size is due to a .so file resulting from the usage of RcppEigen.*

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Haziq Jamil <haziq.jamil@gmail.com>'
  
  New submission
  
  Possibly mis-spelled words in DESCRIPTION:
    Bergsma (8:101)
    Wicher (8:94)
    iprior (8:834)
    lm (8:671)
    
  The Title field should be in title case, current version then in title case: 
  'Linear Regression using I-priors'
  'Linear Regression using I-Priors'

  *This is my first submission. Wicher Bergsma is the name of a colleague whose paper I have cited in the Description field. iprior and lm are R functions which I have mentioned in the Description field. 'I-prior' is a proper noun which has already been correcly capitalised.*
