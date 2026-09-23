## Submission of Bayesrel 0.8.0

This release fixes `omegasCFA()` reporting a zero-width confidence interval
(point estimate as both bounds) when the underlying lavaan fit's information
matrix could not be inverted, which read as a perfectly precise estimate. The
interval is now `NA` in that case, and the function warns when the factor
model did not converge or the solution is inadmissible. `fit.measures = TRUE`
no longer aborts the whole call when lavaan refuses to compute fit measures
for a non-converged model.

The returned object gained `diagnostics` (`converged`, `admissible`,
`se.available`) and `lavaan.fit`.

## Test environments

* local macOS 15 (arm64), R 4.5.2
* GitHub Actions: ubuntu-24.04 (R-devel, R-release, R-oldrel-1),
  macOS (R-release), Windows (R-release): OK
* R-hub v2: linux (R-devel), macos-arm64 (R-devel), windows (R-devel): OK
* win-builder R-devel (2026-09-21 r90579 ucrt): OK
* win-builder R-release (R 4.6.1, 2026-06-24 ucrt): OK

## R CMD check results

0 errors | 0 warnings | 0 notes.

All R-hub and GitHub Actions platforms returned Status: OK.

## Reverse dependencies

There are no reverse dependencies on CRAN.
