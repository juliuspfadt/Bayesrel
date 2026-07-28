## Submission of Bayesrel 0.7.9

This is a patch release fixing two bugs in `bomegas()`:

* the `callback` argument was accepted but never invoked during Gibbs sampling,
* `model.type = "bi-factor"` together with `param.out = TRUE` errored because the
  array for the posterior factor variances was allocated without room for the
  general factor.

## Test environments

* local macOS 15 (arm64), R 4.5.2
* win-builder (devel and release)
* macOS builder (r-release, arm64)
* R-hub: linux, macos, macos-arm64, windows

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies on CRAN.
