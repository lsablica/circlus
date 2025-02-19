## Test environments
* local Ubuntu 24.04 LTS, R 4.3.3
* win-builder 
* GitHub Actions (ubuntu-20.04, macOS-latest, windows-latest) https://github.com/lsablica/circlus/actions/

## R CMD check results
There were no ERRORs or WARNINGs.

It seems that on Ubuntu architectures, the CHECK returns one NOTE because the libs and data sub-directories are above the 1MB threshold.
My understanding is that this inflation of the libs sub-directory is due to the use of Rcpp and RcppArmadillo. Without the speed up gained from those C++ functions, this package would become impractical.
The dataset in data folder has been compressed using the most efficient method available in R.

## Downstream dependencies
none
