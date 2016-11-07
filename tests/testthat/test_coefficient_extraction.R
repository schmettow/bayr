library(testthat)
library(bayr)

context("coefficient extraction")

load("bayr_test.Rda")
C_ <- list() # coefficient tales

test_that("group coefficients can be extracted without issues",{
	expect_silent(C_$Stroop_1$brms <<-
									M$Stroop_1$brms %>% grpef())
	expect_silent(C_$Stroop_1$mcgl <<-
									M$Stroop_1$mcgl %>% grpef())
	expect_silent(C_$Stroop_2$brms <<-
									M$Stroop_2$brms %>% grpef())
	expect_silent(C_$Stroop_2$mcgl <<-
									M$Stroop_2$mcgl %>% grpef())
})


test_that("fixed effects coefficients can be extracted without issues",{
	expect_silent(C_$Stroop_1$brms <<-
									M$Stroop_1$brms %>% fixef())
	expect_silent(C_$Stroop_1$mcgl <<-
									M$Stroop_1$mcgl %>% fixef())
	expect_silent(C_$Stroop_2$brms <<-
									M$Stroop_2$brms %>% fixef())
	expect_silent(C_$Stroop_2$mcgl <<-
									M$Stroop_2$mcgl %>% fixef())
})


test_that("random effects coefficients can be extracted without issues",{
	expect_silent(C_$Stroop_1$brms <<-
									M$Stroop_1$brms %>% ranef())
	expect_silent(C_$Stroop_1$mcgl <<-
									M$Stroop_1$mcgl %>% ranef())
	expect_silent(C_$Stroop_2$brms <<-
									M$Stroop_2$brms %>% ranef())
	expect_silent(C_$Stroop_2$mcgl <<-
									M$Stroop_2$mcgl %>% ranef())
})



test_that("tbl_coef object inherits from tbl_df (dplyr)",{
	expect_is(C_$Stroop_1$mcgl,
						class = c("tbl_coef", "tbl_df"))
	expect_is(C_$Stroop_1$brms,
						class = c("tbl_coef", "tbl_df"))
	expect_is(C_$Stroop_2$mcgl,
						class = c("tbl_coef", "tbl_df"))
	expect_is(C_$Stroop_2$brms,
						class = c("tbl_coef", "tbl_df"))
})





# test_that("MCMCglmm and brms return same order of fixed effects",{
# 	expect_equal(fixef.posterior(P$Stroop_1$mcgl) %>% select(parameter),
# 							 fixef.posterior(P$Stroop_1$brms) %>% select(parameter))
# })
#
# test_that("MCMCglmm and brms return same number of group-level effects",{
# 	expect_equal(P$Stroop_1$mcgl %>%  grpef() %>% select(parameter) %>% nrow(),
# 							 P$Stroop_1$brms %>%  grpef() %>% select(parameter) %>% nrow())
# })
#


## TODO:
# * create test with three-way group effects
# * re-instantiate the last two tests without use of fixef and grpef
