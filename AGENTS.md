# AGENTS.md

Project-specific rules for the `bayr` R package.

## 1. Always return a commit message after changes

Every run that modifies files MUST end with a proposed Git commit message, in
addition to the normal summary of the changes. Write it to the repo's existing
conventions:

- Separate the subject from the body with a blank line.
- Keep the subject to ~50 characters, capitalize it, use the imperative mood,
  and do not end it with punctuation.
- Wrap the body at 72 characters, keep it short, and only include it when it
  adds information the subject does not.
- Do not include the raw diff, file listing, or meta-commentary in the message.
- Propose the message only; do not commit or create branches unless explicitly
  asked.

Present it in a fenced block so it can be copied directly, for example:

```
Fix dependency declarations for CRAN

Add knitr, plyr, purrr, rlang and tibble to Imports and drop the
unused broom.mixed. Remove the install-time library() call that
broke installation.
```

## 2. Use strictly tidy R

Write code in tidyverse style and follow tidy-data principles.

- Prefer `dplyr`, `tidyr`, `stringr`, `purrr`, `rlang` and `tibble` over base R
  equivalents (no `apply`/`sapply` loops, `reshape`, `grep`, etc.).
- Follow tidy-data rules: one variable per column, one observation per row, and
  long format via `pivot_longer()`/`pivot_wider()` instead of matrix-style data.
- Chain transformations with pipes (use `%>%` to match the surrounding code)
  rather than deeply nested calls.
- Use tidy evaluation (`{{ }}`) and the `.data`/`.env` pronouns inside functions
  instead of leaving bare column names that trigger "no visible binding" NOTEs.
- Avoid superseded/deprecated APIs: `across()` not `mutate_all()`/`transmute_all()`,
  `slice_sample()` not `sample_n()`, `pivot_*()` not `gather()`/`spread()`,
  `any_of()` not `one_of()`, `tibble()` not `data_frame()`, `arrange()` not
  `arrange_()`.
- Return tibbles (not data frames or matrices) from user-facing functions.
- Use `snake_case` and stay consistent with the existing `bayr` column scheme
  (`model`, `chain`, `iter`, `parameter`, `value`, `type`, `fixef`, `grpef`,
  `ranef`, `re_factor`, `re_entity`, ...).
- Use `TRUE`/`FALSE`, never `T`/`F`.
- Never call `library()`/`require()` in package code; use `pkg::fun` or
  `@importFrom` and regenerate `NAMESPACE` with roxygen2.

## 3. Ambiguity in multiple subtasks

When commanded to perform multiple subtasks, 
- always perform them sequentially, one after the other, rather than in parallel.
- when ambiguity arises for one or more subtasks, skip these tasks and add a note after the remaining tasks. Outline the ambiguity and make suggestions how to resolve it.
