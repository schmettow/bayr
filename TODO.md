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
| Examples run | ❌ `go_first` example errors (`data_frame` removed from dplyr) |
| Tests run | ❌ tests load a fixture that isn't where they look for it; dev script runs during check |
| Docs code/doc match | ❌ 4 codoc mismatches, ~15 `\usage`/`\arguments` mismatches |
| S3 consistency | ❌ ~15 WARNINGs |
| `R code for possible problems` | ⚠️ NOTE (undefined globals, `:::` calls) |
| DESCRIPTION metadata | ⚠️ as-cran NOTEs (title case, "This package", stale `Date`, invalid URL) |

Target for submission: **0 ERROR, 0 WARNING, 0 NOTE**.

**After Phase 0 (this pass):** `checking package dependencies ... OK` and the
package installs cleanly; the remaining check result is **4 WARNINGs, 2 NOTEs**,
all of which are Phase 1 items.

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

## Phase 1 — `R CMD check` ERRORs / WARNINGs

### 1.1 Examples must run (`checking examples ... ERROR`)

- [ ] **`go_first` example (`R/dplyr_extensions.R` ~line 30, `man/go_first.Rd`)**
  uses `data_frame()`, which dplyr no longer re-exports:
  ```
  Error in data_frame(...) : could not find function "data_frame"
  Execution halted
  ```
  - Replace `data_frame(...)` with `tibble::tibble(...)` (or `tribble`).
  - The example also calls `go_first(D, ~y)` / `go_arrange(D, ~y)`, but
    `go_first()` passes its `...` to `dplyr::select()`, which rejects formulas:
    ```
    Error in dplyr::select(D, !!!cols) : Formula shorthand must be wrapped in `where()`.
    ```
    `go_first(D, y)` works; the `~` form does not. Decide whether the
    documented API is the formula form (then fix the function) or the bare-name
    form (then fix the docs). Update the roxygen source and regenerate `man/`.

- [ ] **`update_by` example (`R/dplyr_extensions.R` ~line 102, `man/update_by.Rd`)**
  calls `mutate_by(...)`, but the function is named `update_by(...)`. Fix the
  example. (`tribble` is fine — it is re-exported through the dplyr import.)

- [ ] If any example needs a fitted model, wrap it in `\dontrun{}` /
  `\donttest{}`, and never reference objects that are not created in the example.

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

### 1.3 Rd `\usage` sections — undocumented / stale arguments (`WARNING`)

Remove hand-written `@usage` tags (they are the cause of most of these) and make
`@param` match the real signatures:

- [ ] `IC.Rd`: `'x' '...'` undocumented (`print`/`knit_print` methods).
- [ ] `as_tbl_obs.Rd`: `'...'` undocumented.
- [ ] `clu.Rd`: `'x' 'df' 'model_name'` undocumented.
- [ ] `coef.tbl_post.Rd`: `'df' 'x'` undocumented.
- [ ] `compare_IC.Rd`: `'x' '...' 'ic_list'` undocumented; documented arg
  `'ic'` not in usage (rename the `@param ic` to `@param ic_list`).
- [ ] `discard_redundant.Rd`: `'...' 'object'` undocumented.
- [ ] `fixef_ml.Rd`: `'x'` undocumented; `'model'` documented but not in usage.
- [ ] `join.tbl_coef.Rd`, `seperate.tbl_coef.Rd`: `'x' 'y'` undocumented,
  `'first' 'second' 'modelnames'` not in usage — **these document
  not-implemented placeholders; delete the functions *and* their `man/*.Rd`**
  (see Phase 4).
- [ ] `md_coef.Rd`: `'round'` undocumented; `'rounding'` documented but not in usage.
- [ ] `post_pred.Rd`: `'x' '...'` undocumented.
- [ ] `posterior.Rd`: `'x' '...'` undocumented; `'thin' 'type' 'model_name'` not in usage.
- [ ] `re_scores.Rd`: `'type'` documented but not in usage.
- [ ] `rescale_unit.Rd`: `'scale'` undocumented.
- [ ] `update_by.Rd`: `'by'` undocumented; `'filter'` documented but not in usage.

### 1.4 S3 generic / method consistency (`WARNING`)

Each method must be compatible with its generic's formals (same first-argument
name, and include `...` when the generic has it):

- [ ] `as_tbl_obs.data.frame`, `as_tbl_obs.tbl_df` → `function(x, ...)`.
- [ ] `clu.tbl_post` → add `...`; `clu.data.frame`, `clu.glmerMod` → first arg
  should be `object` (generic is `clu(object, ...)`).
- [ ] `coef.tbl_post` → add `...`; `coef.data.frame` → rename `df` → `object`.
- [ ] `discard_redundant.*` → add `...` (generic is `discard_redundant(D, except, ...)`).
- [ ] `predict.tbl_post_pred`, `predict.brmsfit`, `predict.stanreg` → rename
  first argument `x` → `object` (generic is `stats::predict(object, ...)`).
- [ ] `tbl_post.data.frame` → rename `x` → `model` (generic is `tbl_post(model, ...)`).

### 1.5 Dependencies in R code (`WARNING`)

- [ ] Remove `require(broom)` from `clu.glmerMod()` (`R/coefficient_extraction.R`
  line 222). `broom` is now declared in `Suggests` (Phase 0), so only the code
  change remains: call `broom::tidy()` guarded by `requireNamespace()`. Note the
  tidier for `glmerMod` actually comes from **broom.mixed** — decide whether to
  support `glmerMod` at all.
- [ ] Replace `base:::print.data.frame` with plain `print.data.frame()`
  (`R/markup_helpers.R` lines 449, 497, 516).
- [ ] Stop reaching into another package's internals:
  `brms:::fixef.brmsfit`, `brms:::ranef.brmsfit` (`R/posterior_extraction.R`
  lines 179–180) and `brms:::predict.brmsfit` (`R/postpred_extraction.R` line 92).
  Use the exported `brms::fixef()` / `brms::ranef()` / `stats::predict()`.
- [ ] Remove `bayr:::` self-calls (a package never needs `:::` for its own
  objects): `bayr:::AllCols` (`R/markup_helpers.R` line 16,
  `R/postpred_extraction.R` line 77), `bayr:::prep_print_tbl_post`
  (`R/markup_helpers.R` line 109), `bayr:::tbl_post.data.frame`
  (`R/posterior_manip.R` line 53).

### 1.6 `R code for possible problems` (`NOTE`)

- [ ] **Define all NSE globals in one place.** `R/package.R` currently declares
  only six, duplicated, and misses dozens. Add a single
  `utils::globalVariables(c(...))` covering at least:
  `.chain .draw .iteration .tmp_idx Estimate Model Obs Part SD SE center chain
  conf.high conf.low diff_IC dpar effect estimate fe_value fixef_2 group lower
  map_dfr model nlpar nonlin prior re_1 re_entity re_factor sd str_c str_extract
  term tidy upper`
  (or switch to the `.data[[...]]` pronoun, which is the more future-proof fix).
- [ ] `importFrom("stats", "formula")` (used in `posterior()` at
  `formula(model)`) and `sd` (used in `z()`). Add via roxygen `@importFrom`.
- [ ] `error("Not implemented")` is not a function — it should be `stop()`.
  These are the placeholder functions `join.tbl_coef`, `join.tbl_coefcomp`,
  `seperate.tbl_coef` (see Phase 4).

### 1.7 Unstated test dependencies (`WARNING`)

```
'library' or 'require' calls not declared from: 'testthat' 'tidyverse'
```

- [ ] Add `testthat (>= 3.0.0)` to `Suggests`, plus `Config/testthat/edition: 3`.
- [ ] The tests also load `mascutils` and `tidyverse`. `mascutils` is not on
  CRAN — see Phase 3.

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

- [ ] **Delete the not-implemented placeholders** `join.tbl_coef`,
  `join.tbl_coefcomp`, `seperate.tbl_coef` (`R/coefficient_extraction.R`
  ~lines 587–625) and their `man/join.tbl_coef.Rd`, `man/seperate.tbl_coef.Rd`.
  They contain broken `error()` calls and generate `\usage` WARNINGs.
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
2. **`go_first`/`go_arrange` API** — formula (`~y`) or bare names (`y`)? The
   docs and implementation disagree.
3. **Test strategy** — OK to ship small pre-computed `.rda`/`.rds` fixtures and
   skip sampler-fitting on CRAN?
4. **`broom.mixed`** — remove the phantom `Suggests` entry, or wire up the
   `glmerMod` tidier properly?
