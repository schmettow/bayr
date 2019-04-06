## Prepare models for testing

library(tidyverse)
library(rstanarm)
library(brms)
library(mascutils)

options(mc.cores = 2)

## GLMM

Ipump <-
	read_csv("tests/Pumps.csv") %>%
	filter(Part <= 10, Task <= 5) %>%
	as_tbl_obs()


F_1 <- formula(ToT ~ Design * session + (1 + Design|Part) + (1|Task))

M_1_b <- brm(F_1, family = "Gaussian", data = Ipump, chains = 2, iter = 200)

M_1_s <- stan_glmer(F_1, data = Ipump, chains = 2, iter = 200)

load("M_1.Rda")

save(M_1_b, M_1_s, file = "M_1.Rda")

P_1_s <- bayr:::tbl_post.stanreg(model = M_1_s)
P_1_b <- bayr:::tbl_post.brmsfit(model = M_1_b)

bayr:::prep_print_tbl_post(P_1_s)
bayr:::prep_print_tbl_post(P_1_b)

bayr:::print.tbl_post(P_1_s)
bayr:::print.tbl_post(P_1_b)

bayr:::knit_print.tbl_post(P_1_s)
bayr:::knit_print.tbl_post(P_1_b)

P_1_s



	# NLM

