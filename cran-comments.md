## Fixes
This version fixes the Rd metadata check NOTEs regarding PKGNAME alias.

## Test environments
r-lib GH action
- {os: macos-latest,   r: 'release'}
- {os: windows-latest, r: 'release'}
- {os: ubuntu-latest,   r: 'devel', http-user-agent: 'release'}
- {os: ubuntu-latest,   r: 'release'}
- {os: ubuntu-latest,   r: 'oldrel-1'}

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 1 NOTEs:

* checking installed package size ... NOTE
  installed size is 19.7Mb
  sub-directories of 1Mb or more:
    libs  18.2Mb
  
  *This NOTE appears only with linux systems (i.e., it does not appear with macOS or win-builder). The size is due to a .so file resulting from the usage of RcppEigen.*

## Downstream dependencies
There are currently no downstream dependencies for this package.
