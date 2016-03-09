library(testthat)
library(bayr)

context("posterior extraction")

load("bayr_test.Rda")
P <- list()

test_that("posterior objects can be extracted without issues",{
	expect_silent(P$Stroop_1$brms <<- M$Stroop_1$brms %>% posterior())
	expect_silent(P$Stroop_1$mcgl <<- M$Stroop_1$mcgl %>% posterior())
	expect_silent(P$Stroop_2$brms <<- M$Stroop_2$brms %>% posterior())
	expect_silent(P$Stroop_2$mcgl <<- M$Stroop_2$mcgl %>% posterior())
})



test_that("posterior object inherits from tbl_df (dplyr)",{
	expect_is(P$Stroop_1$mcgl, class = c("posterior", "tbl_df"))
	expect_is(P$Stroop_1$brms, class = c("posterior", "tbl_df"))
	expect_is(P$Stroop_2$mcgl, class = c("posterior", "tbl_df"))
	expect_is(P$Stroop_2$brms, class = c("posterior", "tbl_df"))
})

test_that("posterior returns chain iter parameter value type order",{
	expect_equal(P$Stroop_1$brms %>% names(),
							 c("chain","iter","parameter","value","type","order"))
	expect_equal(P$Stroop_1$mcgl %>% names(),
							 c("chain","iter","parameter","value","type","order"))
	expect_equal(P$Stroop_2$brms %>% names(),
							 c("chain","iter","parameter","value","type","order"))
	expect_equal(P$Stroop_2$mcgl %>% names(),
							 c("chain","iter","parameter","value","type","order"))
})


test_that("posterior returns a grpef parameter resid ",{
	expect_true("resid" %in% unique(P$Stroop_1$mcgl)$parameter)
	expect_true("resid" %in% unique(P$Stroop_1$brms)$parameter)
	expect_true("resid" %in% unique(P$Stroop_2$mcgl)$parameter)
	expect_true("resid" %in% unique(P$Stroop_2$brms)$parameter)
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
