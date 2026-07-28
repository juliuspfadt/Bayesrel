## Submission of Bayesrel 0.7.9

This is a patch release fixing two bugs in `bomegas()`:

* the `callback` argument was accepted but never invoked during Gibbs sampling,
  so external progress reporting never advanced,
* `model.type = "bi-factor"` together with `param.out = TRUE` errored because the
  array holding the posterior factor variances was allocated without room for the
  general factor.

It also drops `knitr` and `rmarkdown` from Suggests, which were declared but
unused (the package has no vignettes).

## Test environments

* local macOS 15 (arm64), R 4.5.2
* GitHub Actions: ubuntu-24.04 (R-devel, R-release, R-oldrel-1),
  macOS (R-release), Windows (R-release)
* R-hub v2: linux (R-devel), windows (R-devel), macos (R-devel)
* win-builder (R-devel and R-release)

## R CMD check results

0 errors | 0 warnings | 0 notes

All R-hub and GitHub Actions platforms returned Status: OK.

## Reverse dependencies

There are no reverse dependencies on CRAN.
