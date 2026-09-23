# Bayesrel 0.8.0

## Bug fixes

* `omegasCFA()` no longer reports a zero-width confidence interval when the standard
  errors of the omega coefficients could not be computed. lavaan returns the point
  estimate as both interval bounds when the information matrix cannot be inverted,
  which reads as a perfectly precise estimate. The interval is now `NA` and a warning
  is issued. This affects in particular the second-order model with two group factors,
  where the general factor is not identified: the general factor has one loading per
  group factor, but the group factors supply only a single correlation. The point
  estimates are still returned.

* `omegasCFA()` no longer fails entirely when `fit.measures = TRUE` and the fit measures
  cannot be computed, which lavaan refuses for a model that did not converge. The
  coefficients are returned and the fit measures are omitted.

* `omegasCFA()` now warns when the factor model did not converge, and when the solution
  is inadmissible, for instance because of a negative variance estimate. Previously a
  non-converged model returned plausible-looking coefficients without any indication.

## New features

* The object returned by `omegasCFA()` gained `diagnostics`, a list with the entries
  `converged`, `admissible`, and `se.available`, and `lavaan.fit`, holding the fitted
  lavaan object. `print()` reports the diagnostics when any of them fail.

# Bayesrel 0.7.9

## Bug fixes

* `bomegas()` now calls the function supplied to `callback` during the Gibbs
  sampling loops. Previously the argument was accepted but never invoked, so
  external progress reporting never advanced.

* `bomegas()` no longer errors for `model.type = "bi-factor"` in combination
  with `param.out = TRUE`. The array holding the posterior samples of the
  factor variances was allocated without room for the general factor.
