library(dplyr)
library(MCMCglmm)
library(brms)

## archiving system
M <- list()

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
