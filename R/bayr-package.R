#' bayr: Tidy and Unified Reporting of Bayesian Regression Models
#'
#' Provides a unified, tidy interface for reporting Bayesian regression models
#' in the spirit of the Bayesian New Statistics. Wrapper functions extract
#' parameter estimates and posterior distributions from 'brms' and 'rstanarm'
#' models, and tidy summaries from 'lme4' fits, into a common long format.
#'
#' @keywords internal
#' @importFrom assertthat assert_that has_name
#' @importFrom dplyr "%>%" across all_of any_of arrange bind_cols bind_rows
#'   case_when distinct everything filter full_join group_by if_else join_by
#'   left_join matches mutate n n_distinct rename right_join row_number select
#'   slice slice_sample starts_with summarize ungroup where
#' @importFrom rlang enquo enquos quos
#' @importFrom tibble as_tibble tibble tribble
#' @importFrom tidyr crossing extract pivot_longer pivot_wider
"_PACKAGE"
