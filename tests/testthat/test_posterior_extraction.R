library(testthat)
library(bayr)
library(tidyverse)
library(rstanarm)
library(brms)
library(mascutils)

context("posterior extraction")

load("M_1.Rda")

M_1 <- list(M_1_s, M_1_b)

test_that("low-level posterior extraction works without issues",{
	expect_silent(bayr:::tbl_post.brmsfit(M_1_b))
	expect_silent(bayr:::tbl_post.stanreg(M_1_s))
	expect_silent(bayr:::tbl_post(M_1_b))
})

test_that("user-level posterior extraction works without issues",{
	expect_silent(P_1_b <<- bayr::posterior(M_1_b))
	expect_silent(P_1_s <<- bayr::posterior(M_1_s))
})


test_that("posterior object inherits from tbl_df (dplyr)",{
	expect_is(P_1_s, class = c("tbl_post", "tbl_df"))
	expect_is(P_1_b, class = c("tbl_post", "tbl_df"))
})

test_that("posterior object print",{
	expect_silent(bayr:::prep_print_tbl_post(P_1_s))
	expect_output(bayr:::print.tbl_post(P_1_s), regexp = "tbl_post")
	expect_silent(bayr:::knit_print.tbl_post(P_1_s))
})

test_that("extracting random scores",{
	expect_silent(R_1 <<- re_scores(P_1_s))
	expect_s3_class(R_1, "tbl_post")
})



# test_that("posterior returns chain iter parameter value type order",{
# 	expect_equal(P_$Stroop_1$brms %>% names(),
# 							 c("chain","iter","parameter","value","type","order"))
# 	expect_equal(P_$Stroop_1$mcgl %>% names(),
# 							 c("chain","iter","parameter","value","type","order"))
# 	expect_equal(P_$Stroop_2$brms %>% names(),
# 							 c("chain","iter","parameter","value","type","order"))
# 	expect_equal(P_$Stroop_2$mcgl %>% names(),
# 							 c("chain","iter","parameter","value","type","order"))
# })
#
#
# test_that("posterior returns a grpef parameter resid ",{
# 	expect_true("resid" %in% unique(P_$Stroop_1$mcgl)$parameter)
# 	expect_true("resid" %in% unique(P_$Stroop_1$brms)$parameter)
# 	expect_true("resid" %in% unique(P_$Stroop_2$mcgl)$parameter)
# 	expect_true("resid" %in% unique(P_$Stroop_2$brms)$parameter)
# })


# test_that("extraction of wide posterior without issues",{
# 		expect_silent(P_$Stroop_1$brms_w <<- M$Stroop_1$brms %>% posterior(shape = "wide"))
# 		expect_silent(P_$Stroop_1$mcgl_w <<- M$Stroop_1$mcgl %>% posterior(shape = "wide"))
# 		expect_silent(P_$Stroop_2$brms_w <<- M$Stroop_2$brms %>% posterior(shape = "wide"))
# 		expect_silent(P_$Stroop_2$mcgl_w <<- M$Stroop_2$mcgl %>% posterior(shape = "wide"))
# })



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
