library(dplyr)
library(rstanarm)
library(brms)
library(bayr)

thisdir = getwd()
cases_dir = paste0(GDRIVE, "Aktenkoffer/Publications/New_Stats/Cases/")
setwd(cases_dir)
load("CUE8.Rda")
load("BrowsingAB.Rda")
setwd(thisdir)
load(paste0(GDRIVE, "PUBLICATIONS/Publicatie these Stefan Huijser/Data/Lap15.Rda"))

P_1_brms <- posterior(CUE8$M_1)
P_1_brms
fixef(P_1_brms)
grpef(P_1_brms)
ranef(P_1_brms)

P_1_rstn <-
	posterior(CUE8$M_1_rstanarm)

P_1_rstn
fixef(P_1_rstn)
grpef(P_1_rstn)
ranef(P_1_rstn)






## Lap15: non-linear

P_2_brms <-
	posterior(Lap15$M_2_Dur)

P_2_brms
fixef(P_2_brms)
grpef(P_2_brms)
ranef(P_2_brms)


### Case BrowsingAB

M_5 <- stan_glm(ToT ~ Design*age, data = BrowsingAB$BAB1, iter = 200, chains =2)
tbl_post.stanreg(M_5)

M_6 <-
	BrowsingAB$BAB5 %>%
	stan_glmer(ToT ~ Design*age + (1|Task), data = ., iter = 100, chains =2)
tbl_post.stanreg(M_6)

M_7 <- brm(ToT ~ Design*age, data = BrowsingAB$BAB1, iter = 200, chains =2)
tbl_post.brmsfit(M_7)

M_8 <-
	BrowsingAB$BAB5 %>%
	mascutils::z_score(age) %>%
	brm(ToT ~ Design*zage + (1|Task), data = ., iter = 100, chains =2)

P_8 <- tbl_post.brmsfit(M_8)
unique(P_8$parameter)


M_9 <-
	Lap15$D %>%
	brm(Lap15$F_nonlinear_Dur,
		nonlinear = Lap15$F_ME_z,
		prior = Lap15$F_priors,
		chains = 2,
		iter = 200,
		data = .)

P_9 <- tbl_post.brmsfit(M_9)
unique(P_9$parameter)


posterior(M_7)
posterior(M_8)
posterior(M_9)


P_78 <-	bind_rows(posterior(M_7), posterior(M_8))
fixef(P_78)
ranef(P_78)
grpef(P_78)
coef.tbl_post(P_78, type = c("fixef", "ranef", "grpef"))
