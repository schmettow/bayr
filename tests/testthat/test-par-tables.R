test_that("clu() returns a tidy CLU table", {
	P <- tbl_post_fixture()
	C <- clu(P)

	expect_s3_class(C, "tbl_clu")
	expect_true(all(c("model", "parameter", "type", "center", "lower", "upper")
									%in% names(C)))
	expect_false("order" %in% names(C))
	expect_true(all(C$lower <= C$center & C$center <= C$upper))
	expect_equal(nrow(C), 4)
	expect_equal(attr(C, "interval"), .95)
})

test_that("coef() and friends extract the requested parameter types", {
	P <- tbl_post_fixture()

	expect_equal(nrow(coef(P, type = "fixef")), 2)
	expect_equal(nrow(coef(P, type = "ranef")), 1)
	expect_setequal(fixef(P)$fixef, c("Intercept", "GroupB"))
	expect_equal(nrow(grpef(P)), 1)
	expect_equal(nrow(ranef(P)), 1)
	expect_s3_class(fixef(P), "tbl_coef")
})

test_that("coefficient intervals widen with the requested level", {
	P <- tbl_post_fixture()
	narrow <- coef(P, type = "fixef", interval = .5)
	wide   <- coef(P, type = "fixef", interval = 1)

	expect_true(all(narrow$lower >= wide$lower))
	expect_true(all(narrow$upper <= wide$upper))
})

test_that("re_scores() adds the population-level effect to random effects", {
	P <- tbl_post_fixture()
	S <- re_scores(P)

	re <- P[P$type == "ranef", ]
	fe <- P[P$type == "fixef" & P$fixef == "Intercept", ]

	expect_s3_class(S, "tbl_post")
	expect_true(all(S$type == "ranef"))
	expect_equal(sort(S$value), sort(fe$value + re$value))
})

test_that("print and knit_print methods work", {
	P <- tbl_post_fixture()

	expect_output(print(P), "MCMC posterior")
	expect_output(print(clu(P)), "credibility limits")
	expect_output(print(fixef(P)), "Coefficient estimates")
	expect_match(as.character(knitr::knit_print(clu(P))), "credibility limits")
	expect_match(as.character(knitr::knit_print(fixef(P))), "Coefficient estimates")
})

test_that("clu() tidies lme4 models via broom.mixed", {
	skip_if_not_installed("lme4")
	skip_if_not_installed("broom.mixed")

	m <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
									 family = binomial, data = lme4::cbpp)
	C <- clu(m)

	expect_s3_class(C, "tbl_clu")
	expect_true(any(C$type == "fixef"))
	expect_true(any(C$type == "ranef"))
	expect_true(all(c("center", "lower", "upper") %in% names(C)))
})

test_that("clu() validates pre-computed tables", {
	expect_s3_class(clu(clu(tbl_post_fixture())), "tbl_clu")
	expect_error(clu(data.frame(a = 1)))
})
