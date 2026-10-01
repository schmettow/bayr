# bayr — CRAN readiness TODO

Prepared from a hands-on audit of the current `master` tree (`Version: 0.9.8`).
This is **not** a style opinion list: every item below was produced by actually
running the tools on this repo, and the quoted messages are copied from the
check output.

## How this was assessed

- `R CMD build .` — succeeded (`bayr_0.9.8.tar.gz`).
- `R CMD check` on the **unmodified** sources stops immediately:
  ```
  * checking package dependencies ... ERROR
  Package suggested but not available: 'broom.mixed'
  Namespace dependency missing from DESCRIPTION Imports/Depends entries: 'knitr'
  ```
- To see the remaining problems, check was re-run on a throwaway copy with
  `knitr` added to `Imports` and `broom.mixed` removed, forcing
  `_R_CHECK_FORCE_SUGGESTS_=false`:
  ```
  Status: 1 ERROR, 6 WARNINGs, 3 NOTEs
  ```
- `R CMD check --as-cran` adds the CRAN-incoming NOTEs listed in Phase 2.
- Toolchain used: R 4.4.0 (Windows), roxygen2 7.3.1, rcmdcheck available.
  `broom.mixed` is **not installed** here; `mascutils` (used only by tests) is.

## Current status at a glance

| Check | Result |
|---|---|
| Package installs | ✅ fixed in Phase 0 (removed `library(tidyverse)`; globals moved to `R/globals.R`) |
| Dependencies declared | ✅ fixed in Phase 0 (`knitr`+`plyr/purrr/rlang/tibble` added; phantom `broom.mixed` dropped; `broom` added to `Suggests`) |
| Examples run | ✅ fixed in Phase 1 |
| Tests run | ❌ still Phase 3 (fixture location + dev script); test deps now declared except `tidyverse` |
| Docs code/doc match | ✅ fixed in Phase 1 |
| S3 consistency | ✅ fixed in Phase 1 |
| `R code for possible problems` | ✅ fixed in Phase 1 |
| DESCRIPTION metadata | ⚠️ as-cran NOTEs (title case, "This package", stale `Date`, invalid URL) |

Target for submission: **0 ERROR, 0 WARNING, 0 NOTE**.

**After Phase 0:** `checking package dependencies ... OK` and the package
installs cleanly.

**After Phase 1 (this pass):** `Status: 1 WARNING, 1 NOTE` with tests skipped
(tests are Phase 3). `checking examples ... OK`, and the S3, codoc,
`Rd \usage` and "possible problems" checks are all OK. The remaining WARNING
(`tidyverse` in tests) and NOTE (`require(broom)` + `brms:::predict.brmsfit`)
are the two deliberately skipped items / Phase 3 work.

---

## Phase 0 — Blockers ✅ DONE

- [x] **Removed the top-level `library(tidyverse)` from `R/package.R`.** That
  file is deleted; the de-duplicated `utils::globalVariables(...)` call now
  lives in a new `R/globals.R`. The commented-out `library(tidyverse)` and dead
  `globalVariables` blocks were also removed from the other `R/*.R` files.
  (This install-time call was the reason installation failed under `--as-cran`:
  `Error in library(tidyverse) : there is no package called 'tidyverse'`.)
- [x] **Declared `knitr` in `DESCRIPTION` `Imports`.** `NAMESPACE` already had
  `importFrom(knitr,knit_print)`; the code also uses `knitr::kable()` /
  `knitr::asis_output()`.
- [x] **Dropped the phantom `broom.mixed` from `Suggests`** (unused anywhere in
  `R/` or `tests/`).
- [x] **Added the actually-used packages to `Imports`:** `knitr`, `tibble`,
  `rlang`, `purrr`, `plyr`. `broom` (needed only by `clu.glmerMod`) was added to
  `Suggests`. `Imports`/`Suggests` are now sorted alphabetically.

Verified on a fresh `R CMD build` + `R CMD check`:

```
* checking package dependencies ... OK
* checking whether package 'bayr' can be installed ... OK
```

> Note: the old `'library' or 'require' call to 'broom'` and `':::'` findings
> remain as a **NOTE** (not a WARNING) — they are tracked in Phase 1.5, not
> here. `checking dependencies in R code` now only reports the `:::`-related
> items plus the `require(broom)` call.

---

## Phase 1 — `R CMD check` ERRORs / WARNINGs — ✅ DONE (2 items skipped)

### 1.1 Examples must run ✅

- [x] **`go_first` / `go_arrange` example** (`R/dplyr_extensions.R`,
  `man/go_first.Rd`). Fixed by replacing `data_frame()` with `tibble::tibble()`
  and by using the tidyselect form the implementation actually supports:
  `go_first(D, y)`, `go_first(D, y:z)`, `go_arrange(D, y)`. `select()` rejects
  the `~` formula shorthand in every dplyr version, so the documented formula
  form was stale docs rather than a working API; `@param ...` was updated to
  match. Verified: `checking examples ... OK`.
  > ⚠️ **Decision to confirm:** I kept the current function behaviour and fixed
  > the docs. If the formula API (`~y`) was actually intended, the function
  > needs to be changed instead — see "Open questions".
- [x] **`update_by` example** (`R/dplyr_extensions.R`, `man/update_by.Rd`).
  Fixed `mutate_by()` → `update_by()`. The example is now self-contained
  (`tibble::tibble()`, `dplyr::mutate()`, `dplyr::if_else()`, no `%>%`),
  because examples run with only `bayr` attached and cannot see
  imported-but-not-exported helpers such as `tribble` or `%>%`.
- [x] No example needs a fitted model, so no `\dontrun{}` / `\donttest{}`
  wrapping was needed.

### 1.2 Codoc mismatches ✅

All four fixed in the roxygen source and regenerated with
`roxygen2::roxygenise()` (RoxygenNote bumped to 7.3.1): the `as_tbl_obs`
methods now take `...`; the hand-written `@usage` overrides in `post_pred`,
`posterior`, `as_tbl_obs` and `re_scores` were removed; `@param newdata` was
added to `post_pred`. `checking for code/documentation mismatches ... OK`.

### 1.2 Codoc mismatches — `\usage` must match the code (`WARNING`)

Regenerate docs (`devtools::document()`) after fixing. Concrete mismatches:

- [ ] `as_tbl_obs.Rd`: code `function(x, ...)`, docs `function(x)` — the
  methods `as_tbl_obs.data.frame` / `as_tbl_obs.tbl_df` omit `...`.
- [ ] `post_pred.Rd`: docs say `scale = "obs"`, `function(model, scale,
  model_name, thin = 1)`; code is
  `function(model, scale = "resp", model_name = deparse(substitute(model)), newdata = NULL, thin = 1, ...)`.
  Remove the hand-written `@usage` and let roxygen generate it; add `@param newdata`.
- [ ] `posterior.Rd`: hand-written `@usage posterior(model, shape, ...)`
  hides `thin`, `type`, `model_name`, and the `shape = "long"` default.
  Remove the `@usage` override and document all arguments.
- [ ] `re_scores.Rd`: docs omit the `type = "ranef"` argument.

### 1.3 Rd `\usage` sections ✅

All entries resolved by removing stale `@usage` tags and bringing `@param` in
line with the signatures (`@param ic` → `ic_list`, `@param rounding` → `round`,
`@param filter` → `by`, `@param model` → `object` in `fixef_ml`, plus
`@param x` / `@param ...` / `@param scale` / `@param model_name` where methods
and arguments were previously undocumented). The not-implemented placeholders
`join.tbl_coef*` / `seperate.tbl_coef` and their `man/*.Rd` were deleted.
`checking Rd \usage sections ... OK`.

### 1.4 S3 generic / method consistency ✅

All methods now match their generics: `as_tbl_obs.*` gained `...`;
`clu.tbl_post` and `coef.tbl_post` gained `...`; `clu.data.frame`,
`clu.glmerMod` and `coef.data.frame` use `object`; `discard_redundant.*` use
`D` + `...`; `predict.tbl_post_pred/brmsfit/stanreg` use `object`;
`tbl_post.data.frame` uses `model`.
`checking S3 generic/method consistency ... OK`.

### 1.5 Dependencies in R code — ✅ mostly DONE (2 items skipped)

- [x] Removed all `bayr:::` self-calls (`AllCols`, `prep_print_tbl_post`,
  `tbl_post.data.frame`) and replaced all three `base:::print.data.frame`
  calls with `print.data.frame()`.
- [x] Replaced `brms:::fixef.brmsfit` / `brms:::ranef.brmsfit` with the
  exported `brms::fixef()` / `brms::ranef()` (both are exported).
- [ ] ⏭️ **SKIPPED — `require(broom)` in `clu.glmerMod()`.** Whether to keep
  the `glmerMod` backend at all is still open; the `glmerMod` tidier lives in
  `broom.mixed` (not installed here), not in `broom`. Decision needed.
- [ ] ⏭️ **SKIPPED — `brms:::predict.brmsfit`.** The obvious replacement,
  `stats::predict()`, would recurse into `bayr`'s own `predict.brmsfit` method.
  Needs either `brms::posterior_predict()` (verify the returned shape first) or
  an accepted `:::` NOTE. Decision needed.

### 1.6 `R code for possible problems` ✅

`R/globals.R` now declares the full set of NSE globals; the remaining
unqualified helpers were namespace-qualified (`stringr::str_c`,
`purrr::map_dfr`, `stringr::str_extract`); `@importFrom stats formula` and
`@importFrom stats sd` were added.
`checking R code for possible problems ... OK`.

### 1.7 Unstated test dependencies — ⚠️ PARTIAL

- [x] Added `testthat (>= 3.0.0)` to `Suggests` and
  `Config/testthat/edition: 3` to `DESCRIPTION`.
- [ ] Remaining: `'library' or 'require' call not declared from: 'tidyverse'`
  (the tests also use `mascutils`). Deliberately left to Phase 3, which
  rewrites the tests and their fixtures.

---

## Phase 2 — DESCRIPTION metadata (as-cran NOTEs)

- [ ] **Title case.** Current: `tidy and unified reporting of Bayesian regression models`.
  Required: `Tidy and Unified Reporting of Bayesian Regression Models`.
- [ ] **Description must not start with "This package"/package name.**
  Also fix typos: `varios` → `various`, `an unified` → `a unified`. Consider
  adding a reference/`<doi:...>` and confirming brms/rstanarm support matches
  the code (the docs claim MCMCglmm/stanfit/glmerMod too — see below).
- [ ] **Use `Authors@R`** instead of `Author` + `Maintainer`, e.g.
  `Authors@R: person("Martin", "Schmettow", email = "schmettow@web.de", role = c("aut", "cre"))`,
  and drop the standalone `Author`/`Maintainer` fields.
- [ ] **Remove the `Date:` field** (as-cran: *"The Date field is over a month old"*;
  `R CMD build` adds the canonical build timestamp).
- [ ] **Fix the URL**: `http://github.com/schmettow/bayr` redirects;
  as-cran reports *"URL moved to https://github.com/schmettow/bayr"*. Use https
  for both URLs.
- [ ] **Remove `LazyData: TRUE`** (no `data/` directory) — the build already
  prints `Omitted 'LazyData' from DESCRIPTION` and check emits
  `NOTE: 'LazyData' is specified without a 'data' directory`.
- [ ] Keep `Imports` sorted and versioned consistently; bump `RoxygenNote`
  after regenerating with the installed roxygen2 (7.3.1). Confirm
  `Depends: R (>= 3.6.1)` is still the floor you want (the code needs
  rlang/tidyr/dplyr features but no newer base R).

---

## Phase 3 — Tests and fixtures

The tests cannot pass as-is during `R CMD check`:

- [ ] **`tests/prepare_test_models.R` runs during check** (R CMD check sources
  every `.R` in `tests/`). It builds `brms`/`rstanarm` models, reads
  `"tests/Pumps.csv"` (wrong path from the `tests/` working directory), and
  needs `mascutils`/`readr`. Move it to `data-raw/` (and add `^data-raw$` to
  `.Rbuildignore`) so it is never part of the shipped package.
- [ ] **`tests/testthat/*.R` load `"M_1.Rda"`** with a path that does not match
  the file location (`M_1.Rda` currently sits at package root and is shipped
  at the tarball top level). testthat runs with the working directory set to
  `tests/testthat/`, so the file is not found.
  - Either move the fixture to `tests/testthat/` and load it relative to the
    test, or ship it under `inst/extdata/` and load via `system.file()`.
- [ ] **Remove the non-CRAN test dependency `mascutils`.** It is only used to
  load fixture models. Replace with base `load()`/`readr` (add `readr` to
  `Suggests`) or drop it.
- [ ] **Avoid fitting Bayesian models inside `R CMD check`** — it is slow and
  fragile on CRAN machines. Prefer:
  - small pre-computed fixtures (`.rds`/`.Rda`) of a *tbl_post*-like object so
    the pure-data functions are tested without Stan, and
  - `testthat::skip_on_cran()` around anything that must run a sampler.
- [ ] Re-instate the currently commented-out tests once fixtures are stable, or
  delete them — half-dead test files give false confidence.
- [ ] Add `tests/testthat/test-*.R` coverage for the exported helpers that have
  no tests yet (`md_coef`/`frm_coef`, `z_trans`, `rescale_*`, `reorder_levels`,
  `discard_redundant`, `update_by`, `left_union`).

---

## Phase 4 — Housekeeping / artifacts

- [x] **Deleted the not-implemented placeholders** `join.tbl_coef`,
  `join.tbl_coefcomp`, `seperate.tbl_coef` (`R/coefficient_extraction.R`) and
  their `man/join.tbl_coef.Rd`, `man/seperate.tbl_coef.Rd` (done in Phase 1.3).
- [ ] **Silence roxygen2 7.3.1's unregistered-S3-method warnings:**
  `knit_print.tbl_obs_old`, `knit_print.tbl_post_old`
  (`R/markup_helpers.R`) and `mtx_post_pred.data.frame`
  (`R/postpred_extraction.R`) look like S3 methods but are not registered.
  Either delete the dead `*_old` functions or register them if they are meant
  to be used.
- [ ] **Decide the fate of the dead/partial backends.** `NAMESPACE` registers
  methods for `MCMCglmm`, `stanfit` and `glmerMod`, but:
  - there is **no `tbl_post.MCMCglmm`**, so `clu.MCMCglmm`/`coef.MCMCglmm`/
    `fixef.MCMCglmm` can never work; and
  - `MCMCglmm`, `rstan`/`stanfit` are not in `Suggests`, and the docs say
    support is "currently: brms and rstanarm".
  Either implement + declare these backends, or remove the methods to match the
  documented scope. Same question for `clu.glmerMod` (needs broom.mixed).
- [ ] **`cran-comments.Rmd` is stale** (dated 2016, "Windows 7", "R 3.2.3",
  GPL noise) and is excluded from the build. Replace with a current
  `cran-comments.md` that reflects the actual test environments and check
  results for the submission.
- [ ] **Tidy `.Rbuildignore`**: it lists `^.*\.Rproj$` twice and omits the new
  artifacts. Suggested additions:
  ```
  ^cran-comments\.(Rmd|md)$
  ^data-raw$
  ^M_1\.Rda$          # only if not moved into tests/inst
  ^\.github$
  ^.*\.tar\.gz$
  ^.*\.Rcheck$
  ^TODO\.md$
  ```
- [ ] **Add `NEWS.md`** (CRAN likes a changelog; makes release notes easy) and a
  `README.md` / `README.Rmd` for the GitHub landing page.
- [ ] **Add a `LICENSE` note?** Not required for `GPL-3`, but consider
  `License: GPL (>= 3)` for explicitness. Current `GPL-3` is accepted.
- [ ] **Optional but strong for a submission: a getting-started vignette**
  (`vignettes/bayr.Rmd`) showing `posterior()` → `fixef()` → `md_coef()` on a
  small model. This is often the first thing CRAN reviewers/community look for.

---

## Phase 5 — Optional quality improvements (not blocking)

- [ ] `expand_grid()` (`R/dplyr_extensions.R`) shadows `tidyr::expand_grid`,
  which is imported wholesale. `library(bayr)` will emit a masking message.
  Consider renaming (e.g. `expand_grid_df`) or dropping it.
- [ ] Replace superseded APIs so examples/tests don't emit deprecation noise:
  `tidyr::gather`/`spread` → `pivot_longer`/`pivot_wider`;
  `mutate_all`/`transmute_all` → `across()`;
  `arrange_` → `arrange()` + `across()`;
  `dplyr::sample_n` → `slice_sample()`;
  `one_of()` → `any_of()`;
  `data_frame()` → `tibble()`.
- [ ] `discard_all_na()` uses `plyr::aaply(as.matrix(D), 2, ...)`, which
  coerces mixed-type tibbles to character. Prefer a pure-dplyr
  `select(where(~ !all(is.na(.x))))` and drop the `plyr` dependency entirely.
- [ ] Replace bare `T`/`F` with `TRUE`/`FALSE` throughout (they are ordinary
  variables and a common source of bugs).
- [ ] Consider narrowing `@import dplyr` / `@import tidyr` to targeted
  `@importFrom` to reduce the chance of future name clashes.
- [ ] Point `URL:` at a pkgdown site and add `BugReports:` (already present) —
  keep the GitHub links consistent.

---

## Suggested verification workflow (repeat until clean)

```sh
# from the package root
Rscript -e 'if (!requireNamespace("roxygen2", quietly=TRUE)) install.packages("roxygen2")'
Rscript -e 'devtools::document()'        # regenerate NAMESPACE + man/
Rscript -e 'devtools::build()'           # produces bayr_<ver>.tar.gz
_R_CHECK_FORCE_SUGGESTS_=false R CMD check --as-cran bayr_0.9.8.tar.gz
```

Also install `broom.mixed` (or remove it — see Phase 0) before the final run so
the dependency check is exercised honestly. For a real submission, run
`--as-cran` on **both** Windows and Linux, and on the current R release and
R-devel.

## Open questions for you

1. **Scope of backends** — is MCMCglmm / rstan (`stanfit`) / `glmerMod` support
   intended, or should the package officially support only `brms` + `rstanarm`
   as the DESCRIPTION says? This determines whether we implement and declare
   those packages, or delete the dead methods.
2. **`go_first`/`go_arrange` API** — ✅ resolved provisionally in Phase 1: the
   `~` examples were stale (never worked with `select()`), so the docs now use
   the bare-name/tidyselect form (`go_first(D, y)`). Confirm, or the function
   will be changed to accept `~y` instead.
3. **Test strategy** — OK to ship small pre-computed `.rda`/`.rds` fixtures and
   skip sampler-fitting on CRAN?
4. **`broom.mixed`** — remove the phantom `Suggests` entry, or wire up the
   `glmerMod` tidier properly?
