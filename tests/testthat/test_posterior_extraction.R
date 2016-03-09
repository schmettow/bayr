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

	## Stroop
	M$Stroop_1$mcgl <-
		Dlme$Stroop %>%
		MCMCglmm(RT ~ Condition + age + trial,
						 random = ~ idh(Condition):Participant,
						 data = .,
						 pr = T,
						 pl = T,
						 nitt = 550, thin = 5, burnin = 50)
	summary(M$Stroop_1$mcgl)

	M$Stroop_1$brms <-
		Dlme$Stroop %>%
		brm(RT ~ Condition + age + trial + (1 + Condition|Participant),
				data = .,
				pr = T,
				iter = 550, warmup = 50, thin = 5)
	save(M, D, Dlme, file = "bayr_test.Rda")
}

load("bayr_test.Rda")

P$Stroop_1$brms <-
	M$Stroop_1$brms %>% posterior()

P$Stroop_1$mcgl <-
	M$Stroop_1$mcgl %>% posterior()


test_that("brms.posterior returns posterior object, inherited from data_frame",{
	expect_is(P$Stroop_1$brms, class = c("posterior", "tbl_df"))
})

test_that("MCMCglmm.posterior returns posterior object, inherited from data_frame",{
	expect_is(P$Stroop_1$mcgl, "posterior", "tbl_df")
})

test_that("brms.posterior returns chain iter parameter value type order",{
	expect_equal(P$Stroop_1$brms %>% names(),
							 c("chain","iter","parameter","value","type","order"))
})

test_that("MCMCglmm.posterior returns chain iter parameter value type order",{
	expect_equal(P$Stroop_1$mcgl %>% names(),
							 c("chain","iter","parameter","value","type","order"))
})

test_that("MCMCglmm.posterior returns a grpef parameter resid ",{
	expect_true("resid" %in% unique(P$Stroop_1$mcgl)$parameter)
})

test_that("brms.posterior returns a grpef parameter resid ",{
	expect_true("resid" %in% unique(P$Stroop_1$brms)$parameter)
})


test_that("MCMCglmm and brms return same order of fixed effects",{
	expect_equal(fixef.posterior(P$Stroop_1$mcgl) %>% select(parameter),
							 fixef.posterior(P$Stroop_1$brms) %>% select(parameter))
})

test_that("MCMCglmm and brms return same number of group-level effects",{
	expect_equal(P$Stroop_1$mcgl %>%  grpef() %>% select(parameter) %>% nrow(),
							 P$Stroop_1$brms %>%  grpef() %>% select(parameter) %>% nrow())
})



## TODO:
# * create test with three-way group effects
