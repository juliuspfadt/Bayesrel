# Bayesrel 0.7.9

## Bug fixes

* `bomegas()` now calls the function supplied to `callback` during the Gibbs
  sampling loops. Previously the argument was accepted but never invoked, so
  external progress reporting never advanced.

* `bomegas()` no longer errors for `model.type = "bi-factor"` in combination
  with `param.out = TRUE`. The array holding the posterior samples of the
  factor variances was allocated without room for the general factor.
