## Test environments

* local: Windows 11 x64, R 4.4.0 (2024-04-24), 64-bit

<!-- Add win-builder (R-release, R-devel) results here before submitting. -->

## R CMD check results

0 errors | 0 warnings | 0 notes

`R CMD check --as-cran` on the built tarball adds only the standard
"New submission" note. (Runs in the local sandbox also report two environment
artefacts — `unable to verify current time` and `README.md`/`NEWS.md` cannot be
checked without `pandoc` — neither of which occurs on CRAN's check machines.)

## Notes

* The package deliberately registers S3 methods for classes owned by other
  packages: `predict()`/`coef()` for `brmsfit` and `stanreg`, and `clu()` for
  `lme4`'s `glmerMod`. This is intended behaviour kept for backwards
  compatibility; which implementation wins depends on the order in which the
  packages are loaded (documented in the README).
* `Suggests`-only packages (`brms`, `rstanarm`, `lme4`, `broom.mixed`) are used
  conditionally, guarded by `requireNamespace()` or
  `testthat::skip_if_not_installed()`.

## Downstream dependencies

There are currently no downstream dependencies for this package.
