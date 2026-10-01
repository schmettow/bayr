test_that("md_coef() and frm_coef() format coefficient estimates", {
	P <- tbl_post_fixture()
	C <- coef(P, type = "fixef")

	## Intercept draws are 0.11, 0.21, ..., 1.01 twice (both chains) -> median 0.56,
	## 2.5%/97.5% quantiles 0.11/1.01
	expect_equal(md_coef(C, fixef == "Intercept"), "0.56")
	expect_equal(md_coef(C, fixef == "Intercept", center = FALSE, interval = TRUE),
						 " [0.11, 1.01]_{CI95}")
	expect_match(md_coef(C, fixef == "Intercept", interval = TRUE),
						 "^0\\.56 \\[[0-9.-]+, [0-9.-]+\\]_\\{CI95\\}$")
	expect_equal(frm_coef(C, fixef == "Intercept"), "$0.56 [0.11, 1.01]_{CI95}$")
})

test_that("md_coef() can select a row and negate estimates", {
	P <- tbl_post_fixture()
	C <- coef(P, type = "fixef")

	expect_equal(md_coef(C, row = 1), "0.56")
	expect_equal(md_coef(C, row = 1, neg = TRUE), "-0.56")
})

test_that("md_coef() validates its input and selection", {
	C <- coef(tbl_post_fixture(), type = "fixef")

	expect_error(md_coef(C, fixef == "nope"), "does not exist")
	expect_error(md_coef(C, type == "fixef"), "not unique")
	expect_error(md_coef(C, row = 99), "between")
	expect_error(md_coef(data.frame(a = 1)), "coefficient table required")
})
