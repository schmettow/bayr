# bayr 0.9.8

First release prepared for CRAN.

## New features

* `clu()` now works for `lme4` models (`glmerMod`), tidying them with
  `broom.mixed` so a model structure can be checked quickly.

## Bug fixes

* `md_coef()` and `frm_coef()` returned `character(0)` instead of the
  formatted estimate.

## Changes

* Requires `brms (>= 2.16.0)`.
* Removed the unused `MCMCglmm` methods.
* Prediction for `brms` models now uses the public
  `brms::posterior_predict()`.
* Documented the intended overloading of the `brms`/`rstanarm` methods and
  the resulting load-order caveat (see `README.md`).

## Internal

* Rebuilt the test suite around pre-computed fixtures; no model is fitted
  during `R CMD check`.
