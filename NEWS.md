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
* Renamed `expand_grid()` to `expand_grid_df()` so that `library(bayr)` no
  longer masks `tidyr::expand_grid()`.
* Dropped the `plyr` dependency; `discard_all_na()` and `discard_redundant()`
  now use tidy implementations.
* Replaced superseded tidyverse calls (`gather()`/`spread()`,
  `mutate_all()`/`transmute_all()`, `sample_n()`) with their current
  equivalents, and narrowed the dplyr/tidyr imports to the functions used.
* Prediction for `brms` models now uses the public
  `brms::posterior_predict()`.
* Documented the intended overloading of the `brms`/`rstanarm` methods and
  the resulting load-order caveat (see `README.md`).

## Internal

* Rebuilt the test suite around pre-computed fixtures; no model is fitted
  during `R CMD check`.
