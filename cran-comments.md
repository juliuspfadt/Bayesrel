

## R CMD CHECK results
1 WARNING:
The package was formerly archived due to a warning originating from Rcpp (https://cran-archive.r-project.org/web/checks/2022/2022-11-11_check_results_Bayesrel.html). This warning will likely be generated again with this submission. I ask you to ignore it, as it originates from Rcpp and is already fixed (https://github.com/RcppCore/Rcpp/commit/c05fd1bd6b3e3e4fee8ad60be0402ed5e8adbab6). As soon as Rcpp is udpated on CRAN the warning should disappear. On my machine, I have checked the package with the updated Rcpp package from github (https://github.com/RcppCore/Rcpp/commit/97ae61797528b16d6d672c1cf2c2b1f1eaeef8f0) and found no warnings. 

## Submission
Officially this is a new submission, as the package was archived on CRAN.
Practically, this is a resubmission. In this version I have:

* updated the version number

* simplified an example for a function

* added -DARMA_DONT_USE_FORTRAN_HIDDEN_ARGS to Makevars.in to avoid the mismatch warnings of runs with LTO


