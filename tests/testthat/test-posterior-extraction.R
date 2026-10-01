test_that("tbl_post() accepts valid data and rejects invalid data", {
	expect_s3_class(tbl_post_fixture(), "tbl_post")
	expect_error(tbl_post(tibble::tibble(a = 1)))
})

test_that("posterior() returns a long-format tbl_post for brms and rstanarm", {
	skip_if_not_installed("brms")
	skip_if_not_installed("rstanarm")

	fx  <- load_model_fixture()
	P_b <- posterior(fx$M_1_b, model_name = "M_1_b")
	P_s <- posterior(fx$M_1_s, model_name = "M_1_s")

	for (P in list(P_b, P_s)) {
		expect_s3_class(P, "tbl_post")
		expect_s3_class(P, "tbl_df")
		expect_true(all(c("model", "chain", "iter", "order", "parameter", "type",
											"nonlin", "fixef", "re_factor", "re_entity", "value")
										%in% names(P)))
		expect_true("fixef" %in% P$type)
	}

	expect_identical(unique(P_b$model), "M_1_b")
	expect_identical(unique(P_s$model), "M_1_s")
})

test_that("coefficient extraction works on real fits", {
	skip_if_not_installed("brms")
	skip_if_not_installed("rstanarm")

	fx <- load_model_fixture()

	for (fit in list(fx$M_1_b, fx$M_1_s)) {
		expect_s3_class(fixef(fit), "tbl_coef")
		expect_s3_class(clu(fit), "tbl_clu")
		expect_gt(nrow(fixef(fit)), 0)
		expect_gt(nrow(grpef(fit)), 0)
		expect_gt(nrow(ranef(fit)), 0)
	}
})
