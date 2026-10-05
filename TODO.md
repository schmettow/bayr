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
| Dependencies declared | ✅ fixed (Phases 0–1: `knitr`+`plyr/purrr/rlang/tibble` added to `Imports`; `broom.mixed` declared in `Suggests` for the `glmerMod` tidier) |
| Examples run | ✅ fixed in Phase 1 |
| Tests run | ✅ fixed in Phase 3 (`FAIL 0 | PASS 82`; no samplers during check) |
| Docs code/doc match | ✅ fixed in Phase 1 |
| S3 consistency | ✅ fixed in Phase 1 |
| `R code for possible problems` | ✅ fixed in Phase 1 |
| DESCRIPTION metadata | ✅ fixed in Phase 2 (as-cran incoming NOTE is only `Maintainer:` + `New submission`) |
| Quality issues (Phase 5) | ✅ fixed (masking, superseded verbs, many-to-many warnings, `plyr`, imports, pkgdown scaffold) |

Target for submission: **0 ERROR, 0 WARNING, 0 NOTE**.

**After Phase 0:** `checking package dependencies ... OK` and the package
installs cleanly.

**After Phase 1:** `Status: 1 WARNING, 1 NOTE` with tests skipped (tests are
Phase 3). `checking examples ... OK`, and the S3, codoc, `Rd \usage` and
"possible problems" checks are all OK. The remaining WARNING (`tidyverse` in
tests) is Phase 3; the NOTE is now only `brms:::predict.brmsfit`.

**After Phase 2:** plain check unchanged at `Status: 1 WARNING, 1 NOTE`
(tests skipped). With `--as-cran` the incoming-feasibility NOTE contains only
`Maintainer:` and `New submission`; the extra as-cran NOTEs are `M_1.Rda` at
top level (Phase 3), `brms:::predict.brmsfit` (Phase 1.5 skip) and `unable to
verify current time` (environmental, not a package issue). (The `M_1.Rda` and
`brms:::predict.brmsfit` NOTEs have since been resolved — see Phase 3 and
Phase 1.5.)

> **Environment note:** while wiring up `glmerMod`, the user library was
> upgraded (`dplyr` 1.2.1, `rlang` 1.3.0, `broom.mixed` 0.2.9.7 installed via
> Posit PPM binaries). In `dplyr` 1.2.1 `arrange_()` is now *defunct* (errors
> instead of warning), which broke the `go_arrange` example — fixed in Phase 5.
> Other superseded verbs (`sample_n`, `transmute_all`, `one_of`, `gather`,
> `spread`) still work with deprecation warnings; they were replaced in
> Phase 5 as well.

**After Phase 3:** full `R CMD check` (tests enabled) runs in ~50 s with
`checking tests ... OK` (`FAIL 0 | WARN 6 | PASS 96`).

**After removing the last `:::` (this pass):** plain `R CMD check` reports
**`Status: OK`** — no ERRORs, WARNINGs or NOTEs. With `--as-cran` the remaining
`3 NOTEs` are all non-package: `New submission` (boilerplate), `unable to verify
current time` (environmental) and `README.md ... cannot be checked without
'pandoc'` (pandoc is not installed in this sandbox; CRAN check machines have
it). The 6 test WARNs are dplyr's many-to-many join notices (Phase 5).

**After Phase 4:** the full check *with the vignette built* (pandoc, see below)
and the tests enabled is still **`Status: OK`**. `--as-cran` reports only
environment artefacts (`qpdf`, `README.md`/`NEWS.md` without pandoc on `PATH`,
`unable to verify current time`) plus the boilerplate `New submission`.

**After Phase 5:** `R CMD check --no-manual` on the rebuilt tarball (tests and
vignette enabled) is still **`Status: OK`**; the `WARN 6` many-to-many join
notices are gone (`FAIL 0 | WARN 0 | PASS 96`). `plyr` is no longer a
dependency, and dplyr/tidyr are imported function-by-function (see
`R/bayr-package.R`).

> ⚠️ The pkgdown site (`_pkgdown.yml`, `.github/workflows/pkgdown.yaml`) is
> scaffolded but not yet deployed. `DESCRIPTION` therefore keeps the live book
> URL; `--as-cran` reports `https://schmettow.github.io/bayr/` as a 404 as long
> as it is listed. After the first push has deployed GitHub Pages, switch the
> URL to the pkgdown site (see Phase 5).

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
  `rlang`, `purrr`, `plyr`. `Imports`/`Suggests` are sorted alphabetically.
  (The `glmerMod` tidier ended up needing `broom.mixed`, not `broom` — see
  Phase 1.5.)

Verified on a fresh `R CMD build` + `R CMD check`:

```
* checking package dependencies ... OK
* checking whether package 'bayr' can be installed ... OK
```

> Note: the old `'library' or 'require' call to 'broom'` and `':::'` findings
> were tracked to Phase 1.5 and are now resolved (the `glmerMod` tidier uses
> `broom.mixed`; only the `brms:::predict.brmsfit` NOTE remains).

---

## Phase 1 — `R CMD check` ERRORs / WARNINGs ✅ DONE

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

### 1.5 Dependencies in R code ✅

- [x] Removed all `bayr:::` self-calls (`AllCols`, `prep_print_tbl_post`,
  `tbl_post.data.frame`) and replaced all three `base:::print.data.frame`
  calls with `print.data.frame()`.
- [x] Replaced `brms:::fixef.brmsfit` / `brms:::ranef.brmsfit` with the
  exported `brms::fixef()` / `brms::ranef()` (both are exported).
- [x] **`clu.glmerMod()` — kept, now uses `broom.mixed`.** Replaced
  `require(broom)` with `requireNamespace("broom.mixed")` plus a clear error
  message, and the unqualified `tidy()` with `broom.mixed::tidy(conf.int = TRUE)`
  (the empty `conf.level = ` argument was a bug and was dropped). `broom` was
  swapped for `broom.mixed` in `Suggests`. Runtime-tested on a `glmer()`
  binomial model; returns fixed, SD and random-effect rows as expected.
- [x] **`brms:::predict.brmsfit` → `brms::posterior_predict()`** (Option B).
  `stats::predict()` would have recursed into `bayr`'s own `predict.brmsfit`,
  and the unexported brms method could change without notice. brms's own
  `predict.brmsfit` is a thin wrapper over the public `posterior_predict()`
  (`prepare_predictions()` + `posterior_predict(..., summary = summary)`) and
  both return the same matrix shape, so the call now passes
  `ndraws = n_draws, summary = FALSE`. This removes the last `:::` NOTE.
  `Suggests: brms` floor raised to `(>= 2.16.0)` (already required by the
  `as_draws_df()` call in `tbl_post.brmsfit()`).

### 1.6 `R code for possible problems` ✅

`R/globals.R` now declares the full set of NSE globals; the remaining
unqualified helpers were namespace-qualified (`stringr::str_c`,
`purrr::map_dfr`, `stringr::str_extract`); `@importFrom stats formula` and
`@importFrom stats sd` were added.
`checking R code for possible problems ... OK`.

### 1.7 Unstated test dependencies ✅

- [x] Added `testthat (>= 3.0.0)` to `Suggests` and
  `Config/testthat/edition: 3` to `DESCRIPTION`.
- [x] The `tidyverse` / `mascutils` calls were removed in Phase 3;
  `checking for unstated dependencies in 'tests' ... OK`.

---

## Phase 2 — DESCRIPTION metadata ✅ DONE

- [x] **Title case:** `Tidy and Unified Reporting of Bayesian Regression Models`.
- [x] **Description rewritten** so it no longer starts with "This package";
  typos fixed (`varios` → `various`, `an unified` → `a unified`) and software
  names quoted (`'brms'`, `'rstanarm'`, `'lme4'`, `'knitr'`, `'Markdown'`).
  A reference was added: Kruschke and Liddell (2018),
  <doi:10.1177/2515245918771329>. The wording now matches the confirmed
  backend scope (`glmerMod` included via `broom.mixed`).
- [x] **`Authors@R`** added (`aut` + `cre`); the standalone `Author` /
  `Maintainer` fields were removed.
- [x] **`Date:` removed** (`R CMD build` supplies the build timestamp).
- [x] **URLs fixed to https** (website and GitHub).
- [x] **`LazyData: TRUE` removed** (there is no `data/` directory).
- [x] `Imports`/`Suggests` sorted; `RoxygenNote` is 7.3.1.
  `Depends: R (>= 3.6.1)` left unchanged — no newer base-R features are used.
- [x] Verified with `R CMD check --as-cran`: the incoming-feasibility NOTE now
  contains only the expected `Maintainer:` and `New submission` lines.
  (The DOI was accepted and the URL check is clean.)

---

## Phase 3 — Tests and fixtures ✅ DONE

The test suite was rebuilt around a pre-computed fixture; no model is fitted
during `R CMD check`.

- [x] **`tests/prepare_test_models.R` moved to `data-raw/prepare_test_models.R`**
  (together with `Pumps.csv`) and rewritten to read from `data-raw/` and write
  the fixture to `tests/testthat/M_1.Rda`. `data-raw/` is in `.Rbuildignore`,
  so the script is never shipped and never runs during check.
- [x] **Fixture relocated** to `tests/testthat/M_1.Rda` and re-compressed with
  `compress = "xz"` (1.9 MB → 0.9 MB); tests load it via
  `testthat::test_path("M_1.Rda")`. This also removes the non-standard
  top-level file reported by `--as-cran`.
- [x] **`mascutils` (and `tidyverse`) removed from the tests.** They now need
  only `brms`, `rstanarm`, `lme4` and `broom.mixed` (all in `Suggests`), each
  guarded with `skip_if_not_installed()`.
- [x] **No samplers run during check.** The extraction tests load the
  pre-computed fits; everything else runs on a dependency-free synthetic
  `tbl_post` built in `helper-fixtures.R` (`tbl_post_fixture()`).
- [x] **Dead, commented-out tests deleted** and replaced with focused coverage:

| File | Covers |
|---|---|
| `test-posterior-extraction.R` | `tbl_post()` validation, `posterior()` and `fixef()`/`grpef()`/`ranef()`/`clu()` on real `brms`/`rstanarm` fits |
| `test-par-tables.R` | `clu()`/`coef()`/`fixef()`/`grpef()`/`ranef()`/`re_scores()`, print/`knit_print` methods, `lme4` + `broom.mixed` path |
| `test-helpers.R` | `z_trans()`, `rescale_*()`, `reorder_levels()`, `discard_redundant()`, `discard_all_na()`, `update_by()`, `left_union()`, `go_first()`/`go_arrange()`, `expand_grid()` |
| `test-md-coef.R` | `md_coef()`/`frm_coef()` formatting, selection and error paths |

Bugs the new tests found and that were fixed in this pass:

- [x] **`md_coef()`/`frm_coef()` always returned `character(0)`** — the
  accumulator was initialised with `as.character()` (length 0) and
  `stringr::str_c()` then recycled to zero length. `R/markup_helpers.R` now
  starts from `""`.
- [x] **`brms::fixef()`/`brms::ranef()` caused infinite recursion** in
  `extr_brms_par()`: `brms` and `bayr` both register `fixef.brmsfit`, so
  `brms::fixef()` dispatched back into `bayr`'s method → `tbl_post()` →
  `extr_brms_par()` → ... `extr_brms_par()` now derives the fixed-effect names
  from the draw names (`b_*`) and the unused `ranef()` call was dropped.
- [x] Deprecated tidyselect usage surfaced by the tests replaced:
  `select(cols)` → `select(all_of(cols))` in `pre_print_tbl_clu()` and
  `one_of()` → `any_of()` in `discard_all_na()`.

---

## Phase 4 — Housekeeping / artifacts ✅ DONE

- [x] **Deleted the not-implemented placeholders** `join.tbl_coef`,
  `join.tbl_coefcomp`, `seperate.tbl_coef` (`R/coefficient_extraction.R`) and
  their `man/join.tbl_coef.Rd`, `man/seperate.tbl_coef.Rd` (done in Phase 1.3).
- [x] **Silenced roxygen2's unregistered-S3-method warnings** by deleting the
  three dead functions: `knit_print.tbl_obs_old`, `knit_print.tbl_post_old`
  (`R/markup_helpers.R`) and `mtx_post_pred.data.frame`
  (`R/postpred_extraction.R`). The latter was also broken — it validated against
  the posterior `AllCols` schema instead of `Cols_pp` — and was unreachable
  because no `mtx_post_pred.data.frame` method was registered.
- [x] **MCMCglmm backend removed** (`clu`/`coef`/`fixef`/`ranef`/`grpef`
  methods, their NAMESPACE registrations and the commented-out
  `predicted.MCMCglmm` prototype). There was no `tbl_post.MCMCglmm`, so the
  methods could never work; the `@param object` docs now read
  "(brms, rstanarm)".
- [x] **`stanfit` methods removed** (`clu`/`coef`/`fixef`/`ranef`/`grpef`, the
  commented-out `predicted.stanfit` prototype and the now-unused
  `importFrom(stats, fitted)`). There was no `tbl_post.stanfit`, so they could
  never work. With that, the backend scope is settled: `brms`, `rstanarm` and
  `glmerMod` (Phase 6).
- [x] **`cran-comments.Rmd` replaced by a current `cran-comments.md`**
  (`git mv`, content rewritten): local test environment, the 0/0/0 check
  results, a note on the intended S3 overloading and on the conditional use of
  `Suggests` packages, and the downstream-dependencies statement. A marker is
  left for the win-builder / R-devel results to be added before submitting.
- [x] **`.Rbuildignore` tidied** (Phase 2 pass): the duplicated `^.*\.Rproj$`
  was removed and `^cran-comments\.(Rmd|md)$`, `^AGENTS\.md$`, `^TODO\.md$`,
  `^data-raw$` and `^\.github$` were added. The fixture now lives under
  `tests/` (Phase 3), so no `M_1.Rda` entry is needed.
- [x] **`README.md` added**, including the note that overloading the
  `brms`/`rstanarm` `predict()`/`coef()` methods is intended behaviour and that
  the load order decides which implementation wins.
- [x] **`NEWS.md` added** with the 0.9.8 release notes (new `clu()` support for
  lme4, the `md_coef()` fix, the `MCMCglmm` removal, the brms floor and the
  rebuilt tests).
- [x] **License changed to MIT**: `License: MIT + file LICENSE` in
  `DESCRIPTION` plus the `LICENSE` file
  (`YEAR: 2016-2026`, `COPYRIGHT HOLDER: Martin Schmettow`), and the README
  updated.
  > ⚠️ This relicenses the package from GPL-3 to MIT — please confirm you hold
  > the rights to all of the code (DESCRIPTION lists a single author; first
  > commit 2016-03-07).
- [x] **Getting-started vignette added** (`vignettes/bayr.Rmd`, "Getting
  started with bayr"), following the *fitting regression models* section of
  <https://schmettow.github.io/New_Stats/gsr.html#fitting>: simulated data →
  `as_tbl_obs()` → two `rstanarm::stan_glm()` fits → `clu()`, `fixef()`,
  `md_coef()`, `posterior()`, `predict()`. The chunks are skipped when
  `rstanarm` is unavailable; `VignetteBuilder: knitr` and `Suggests: rmarkdown`
  were added. Builds to `inst/doc/bayr.html` and is rebuilt cleanly during
  check.

---

## Phase 5 — Optional quality improvements ✅ DONE

- [x] **`expand_grid()` renamed to `expand_grid_df()`** (`R/dplyr_extensions.R`,
  `man/expand_grid_df.Rd`): `library(bayr)` after `tidyr` no longer masks
  `tidyr::expand_grid()`. The old export is gone; the rename is documented in
  `NEWS.md`. This is a deliberate API change.
- [x] `arrange_` → `arrange(!!!rlang::syms(...))` in `go_arrange()`
  (`arrange_()` became **defunct** in dplyr 1.2.1, so this was blocking).
- [x] `one_of()` → `any_of()` in `discard_all_na()` (Phase 3).
- [x] **Superseded verbs replaced:** `tidyr::gather()`/`spread()` →
  `pivot_longer()`/`pivot_wider()` (`tbl_post.brmsfit`, `tbl_post.stanreg`,
  `fixef_ml`, and the dead wide-format branch in `posterior()`);
  `transmute_all(z)` → `mutate(across(everything(), z))` (`z_trans`);
  `dplyr::sample_n()` → `slice_sample()` (all `print()` / `knit_print()`
  methods).
- [x] **Many-to-many join warning removed** (`tbl_post.brmsfit`). The
  `expand.grid(...) %>% full_join(type_patterns, by = "type")` scaffold was a
  parameter × pattern crossing in disguise (`"shape"` maps to three regex
  entries) and dplyr 1.2.1 flagged it on every `posterior()` call. It is now
  `tidyr::crossing(par_order["parameter"], type_patterns)` plus the same
  `str_detect()` filter — no `relationship =` override and no dplyr version
  bump needed.
- [x] **`discard_all_na()` / `discard_redundant.default()` rewritten with tidy
  tools** (`select(where(any_not_na))`, `purrr::map_lgl()` + `n_distinct()`)
  instead of `plyr::aaply(as.matrix(D), 2, ...)`, which coerced mixed-type
  tibbles to character. `plyr` was dropped from `Imports`.
- [x] **Bare `T`/`F` replaced with `TRUE`/`FALSE`** in all `R/*.R` files
  (function defaults, `na.rm =`, `remove =`, `row.names =`, plus the
  commented-out prototypes).
- [x] **Wholesale `@import dplyr` / `@import tidyr` / `@import assertthat`
  narrowed to function-level `@importFrom`** in the new `R/bayr-package.R`
  (dplyr incl. `%>%`, rlang `enquo`/`enquos`/`quos`, tibble, tidyr,
  assertthat); `NAMESPACE` contains no `import()` calls any more. This
  surfaced one new codetools NOTE (`count` used as a column in `join_by()`),
  fixed by adding `count` to `utils::globalVariables()`.
- [x] **pkgdown scaffold added:** `_pkgdown.yml` (bootstrap 5) and
  `.github/workflows/pkgdown.yaml` (build + deploy to `gh-pages` on push);
  `^_pkgdown\.yml$`/`^docs$` added to `.Rbuildignore`. `URL:`/`BugReports:`
  stay on the **live** URLs (book + GitHub) for now: `--as-cran` checks URLs
  and flags the not-yet-deployed pkgdown site as `Status: 404`. After the
  first push has deployed the site, replace
  `https://schmettow.github.io/New_Stats/` with
  `https://schmettow.github.io/bayr/` in `DESCRIPTION` and re-run
  `R CMD check --as-cran`.

---

## Phase 6 — Complete the `glmerMod` interface (backlog, post-CRAN)

`clu.glmerMod()` is kept as a working stub: it tidies `lme4` fits via
`broom.mixed`, so a model structure can be checked quickly without refitting
the same model in `brms`/`rstanarm`. The rest of the interface is not
implemented yet.

- [ ] `coef.glmerMod()` / `fixef.glmerMod()` — population-level effects with
  intervals (`broom.mixed::tidy(conf.int = TRUE)`), mapped to the `tbl_coef`
  column scheme (`model`, `type`, `nonlin`, `fixef`, `re_factor`, `re_entity`,
  `center`, `lower`, `upper`).
- [ ] `grpef.glmerMod()` — group-level SDs/correlations, currently emitted by
  `clu.glmerMod()` as `type == "sd"` / `"cor"` rows; standardise to `grpef`.
- [ ] `ranef.glmerMod()` — random-effect modes; `lme4::ranef()` returns a named
  list of data frames (one per grouping factor) and needs reshaping to the long
  `tbl_coef` scheme (see the commented-out `clu.glm` prototype).
- [ ] Make `assert_clu()` accept the lme4 output (`re_entity` is currently a
  character column there) and check whether `assert_coef()` is needed.
- [ ] Verify the `print()` / `knit_print()` paths for the produced tables.
- [ ] Tests: extend the single `clu()` smoke test on `lme4::cbpp`
  (`tests/testthat/test-par-tables.R`) to cover the new methods.
- [ ] Docs: mention `glmerMod` in `clu()`'s `@param object`; consider a "quick
  model checks with lme4" section in the vignette.
- [ ] Decide whether `posterior()` / `post_pred()` should support `glmerMod` at
  all (they are Bayesian-specific; currently not supported).

---

## Suggested verification workflow (repeat until clean)

```sh
# from the package root
Rscript -e 'if (!requireNamespace("roxygen2", quietly=TRUE)) install.packages("roxygen2")'
Rscript -e 'devtools::document()'        # regenerate NAMESPACE + man/
Rscript -e 'devtools::build()'           # produces bayr_<ver>.tar.gz
_R_CHECK_FORCE_SUGGESTS_=false R CMD check --as-cran bayr_0.9.8.tar.gz
```

`broom.mixed` is now installed and declared in `Suggests`. For a real
submission, run
`--as-cran` on **both** Windows and Linux, and on the current R release and
R-devel.

> `checking for future file timestamps ... NOTE: unable to verify current
> time` is an environment artefact (no reliable time source in this sandbox),
> not a package issue.

> **Building the vignette needs pandoc.** It is not on `PATH` in this sandbox,
> but RStudio bundles one — point `rmarkdown` at it instead of installing
> pandoc by setting
> `RSTUDIO_PANDOC="C:/Users/SchmettowM/AppData/Local/Programs/RStudio/resources/app/bin/quarto/bin/tools"`
> before `R CMD build` / `R CMD check`. Two local-only artefacts then remain:
> `README.md`/`NEWS.md` cannot be checked without pandoc *on `PATH`*, and
> `'qpdf' is needed for checks on size reduction of PDFs`. CRAN check machines
> have both installed.

## Open questions for you

1. ~~Scope of backends~~ — ✅ resolved: `MCMCglmm` and `stanfit` methods removed
   (neither could work without a `tbl_post` method), `glmerMod` kept with its
   completion tracked in Phase 6.
2. **`go_first`/`go_arrange` API** — ✅ resolved provisionally in Phase 1: the
   `~` examples were stale (never worked with `select()`), so the docs now use
   the bare-name/tidyselect form (`go_first(D, y)`). Confirm, or the function
   will be changed to accept `~y` instead.
3. **Test strategy** — ✅ resolved provisionally in Phase 3: the suite ships a
   pre-computed `brms`/`rstanarm` fixture (`tests/testthat/M_1.Rda`) plus a
   dependency-free synthetic `tbl_post`, and no sampler runs during
   `R CMD check`.
4. ~~`broom.mixed`~~ — ✅ resolved: declared in `Suggests` and used for the
   `glmerMod` tidier.
