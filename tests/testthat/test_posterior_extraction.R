library(testthat)
library(bayr)
library(dplyr)
library(tidyr)
library(stringr)
library(MCMCglmm)
library(brms)

context("posterior extraction")

purp.mcmc = F

## archiving system
M <- list()
P <- list() #posterior distributions
if(file.exists("bayr_test.Rda")) load("bayr_test.Rda")

## redo MCMC runs
if(purp.mcmc){
	library(MCMCglmm)
	library(brms)
	load("Book.Rda")
	load("Classic_linear_models.Rda")
	load("Linear_mixed-effects_models.Rda")

	## Stroop intercept-only model
	M$Stroop_1$mcgl <-
		MCMCglmm(RT ~ 1,
						 data = Dlme$Stroop,
						 pr = T,
						 pl = T,
						 nitt = 150, thin = 5, burnin = 50)
	# summary(M$Stroop_1$mcgl)

	M$Stroop_1$brms <-
		brm(RT ~ 1,
				data = Dlme$Stroop,
				iter = 150, warmup = 50, thin = 5)

	## Stroop mixed-effects model
	M$Stroop_2$mcgl <-
		MCMCglmm(RT ~ Condition + age + trial,
						 random = ~ idh(Condition):Participant,
						 data = Dlme$Stroop,
						 pr = T,
						 pl = T,
						 nitt = 450, thin = 5, burnin = 150)
	# summary(M$Stroop_1$mcgl)

	M$Stroop_2$brms <-
		brm(RT ~ Condition + age + trial + (1 + Condition|Participant),
				data =
					Dlme$Stroop,
				iter = 150, warmup = 50, thin = 5)
	names(M)
	save(M, D, Dlme, file = "bayr_test.Rda")
}

load("bayr_test.Rda")


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
