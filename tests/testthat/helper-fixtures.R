## A small, dependency-free tbl_post used by the pure-data tests. It mimics the
## output of posterior() for a model with two fixed effects, one group-level
## standard deviation and one random effect.
tbl_post_fixture <- function(model = "m1", n_iter = 10) {
	pars <- tibble::tribble(
		~parameter,            ~type,    ~fixef,      ~re_factor,    ~re_entity,    ~order,
		"b_Intercept",         "fixef",  "Intercept", NA_character_, NA_character_, 1L,
		"b_GroupB",            "fixef",  "GroupB",    NA_character_, NA_character_, 2L,
		"sd_Part__Intercept",  "grpef",  "Intercept", "Part",        NA_character_, 3L,
		"r_Part[1,Intercept]", "ranef",  "Intercept", "Part",        "1",           4L
	)

	out <- tidyr::expand_grid(pars, chain = 1:2, iter = seq_len(n_iter))
	out$model <- model
	out$nonlin <- NA_character_
	out$value <- out$iter / n_iter + out$order / 100
	out <- out[, c("model", "chain", "iter", "order", "parameter", "type",
								 "nonlin", "fixef", "re_factor", "re_entity", "value")]
	bayr::tbl_post(out)
}

## The pre-computed brms/rstanarm fits used by the extraction tests
## (see data-raw/prepare_test_models.R). Loading them needs brms and rstanarm,
## so callers must skip_if_not_installed() first.
load_model_fixture <- function() {
	env <- new.env(parent = emptyenv())
	load(testthat::test_path("M_1.Rda"), envir = env)
	env
}
