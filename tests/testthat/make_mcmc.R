library(dplyr)
library(MCMCglmm)
library(brms)
library(rstanarm)

setwd("tests/testthat/")
## archiving system
M <- list()

load("Book.Rda")
load("Classic_linear_models.Rda")
load("Linear_mixed-effects_models.Rda")
load("bayr_test.Rda")

## Stroop intercept-only model
M$Stroop_1$mcgl <-
	MCMCglmm(RT ~ 1,
					 data = Dlme$Stroop,
					 pr = T,
					 pl = T,
					 nitt = 150, thin = 5, burnin = 50)

M$Stroop_1$brms <-
	brm(RT ~ 1,
			data = Dlme$Stroop,
			iter = 150, warmup = 50, thin = 5)

M$Stroop_1$arm <-
	stan_glm(RT ~ 1,
			data = Dlme$Stroop,
			algorithm = "fullrank",
			iter = 150, warmup = 50, thin = 5)



## Stroop mixed-effects
M$Stroop_2$mcgl <-
	MCMCglmm(RT ~ Condition + age + trial,
					 random = ~ idh(Condition):Participant,
					 data = Dlme$Stroop,
					 pr = T,
					 pl = T,
					 nitt = 450, thin = 5, burnin = 150)

M$Stroop_2$brms <-
	brm(RT ~ Condition + age + trial + (1 + Condition|Participant),
			data = Dlme$Stroop,
			iter = 150, warmup = 50, thin = 5)
names(M)

M$Stroop_2$arm <-
	stan_lmer(RT ~ Condition + age + trial + (1 + Condition|Participant),
			data = Dlme$Stroop,
			iter = 150)


## Stroop ANOVA
M$Stroop_3$mcgl <-
	MCMCglmm(RT ~ Condition + age,
					 data = Dlme$Stroop,
					 pr = T,
					 pl = T,
					 nitt = 150, thin = 5, burnin = 50)

M$Stroop_3$brms <-
	brm(RT ~ Condition + age,
			data = Dlme$Stroop,
			iter = 150, warmup = 50, thin = 5)

M$Stroop_3$arm <-
	stan_glm(RT ~ Condition + age,
					 data = Dlme$Stroop,
					 algorithm = "fullrank",
					 iter = 150, warmup = 50, thin = 5)




names(M$Stroop_1)
names(M$Stroop_2)
names(M$Stroop_3)

save(M, D, Dlme, file = "bayr_test.Rda")
