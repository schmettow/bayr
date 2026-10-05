# bayr

Tidy and unified reporting of Bayesian regression models.

`bayr` extracts posterior draws from Bayesian regression models into a common
long format, so that parameter estimates, predictions and information criteria
can be reported with tidy processing chains and `knitr`/Markdown output.
Wrapper functions cover `brms` and `rstanarm` models; `lme4` fits can be
summarised quickly via `broom.mixed`.

## Installation

```r
# from GitHub
remotes::install_github("schmettow/bayr")
```

## Usage

```r
library(rstanarm)
library(bayr)

fit <- stan_glm(mpg ~ wt + cyl, data = mtcars)

posterior(fit)   # long-format draws         -> tbl_post
clu(fit)         # center-lower-upper table  -> tbl_clu
fixef(fit)       # fixed effects             -> tbl_coef
grpef(fit)       # group-level SDs           -> tbl_coef
predict(fit)     # tidy predictions          -> tbl_predicted
```

## Overloaded `brms` and `rstanarm` methods — load order matters

`bayr` deliberately registers S3 methods for classes owned by other packages,
namely `predict()` and `coef()` for `brmsfit` and `stanreg` objects. That
overloading is what makes the tidy model-level shortcuts above work, and it is
**intended behaviour**, kept for backwards compatibility with existing code.

R has exactly one slot per S3 method, so when `brms` or `rstanarm` is loaded
*after* `bayr`, it takes the slot back and R reports, for example:

```
Registered S3 methods overwritten by 'brms':
  method          from
  coef.brmsfit    bayr
  predict.brmsfit bayr
```

Which implementation you get therefore depends on the order in which the
packages are loaded:

| Load order | `predict(fit)` / `coef(fit)` |
| --- | --- |
| `library(brms)` / `library(rstanarm)` **before** `library(bayr)` | `bayr`'s tidy output (`tbl_predicted`, `tbl_coef`) |
| `library(bayr)` **before** `brms` / `rstanarm` | the model package's own output — for a `stan_glmer` fit, `rstanarm::predict()` even errors and asks for `posterior_predict()` |

**If you want `bayr`'s methods, load `bayr` last.** The same caveat applies to
`fixef()` and `ranef()`: those names are also exported by `nlme`/`lme4`, so
attaching `lme4` after `bayr` masks `bayr`'s generics as well.

The tidy pipeline itself is independent of the load order: `bayr::post_pred()`
and bayr's own generics (`clu()`, `grpef()`, `posterior()`, `tbl_post()`, and
`predict()` on a `tbl_post_pred`) always behave the same, whichever package was
loaded last. In particular, the following always returns a `tbl_predicted`:

```r
predict(post_pred(fit))
```

## License

GPL-3
