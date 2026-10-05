test_that("post_pred() and predict() return tidy predictions", {
	skip_if_not_installed("brms")
	skip_if_not_installed("rstanarm")

	fx <- load_model_fixture()

	for (fit in list(fx$M_1_b, fx$M_1_s)) {
		PP <- post_pred(fit, thin = 10)

		expect_s3_class(PP, "tbl_post_pred")
		expect_true(all(c("model", "Obs", "chain", "iter", "scale", "value")
										%in% names(PP)))
		expect_gt(nrow(PP), 0)

		PR <- predict(PP)

		expect_s3_class(PR, "tbl_predicted")
		expect_true(all(c("model", "Obs", "center", "lower", "upper")
										%in% names(PR)))
		expect_true(all(PR$lower <= PR$center & PR$center <= PR$upper))
	}
})

test_that("the model-level predict method is built on post_pred()", {
	skip_if_not_installed("brms")

	fx <- load_model_fixture()
	PR <- bayr:::predict.brmsfit(fx$M_1_b, thin = 10)

	expect_s3_class(PR, "tbl_predicted")
	expect_gt(nrow(PR), 0)
})
