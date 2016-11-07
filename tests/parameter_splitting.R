library(tidyverse)
library(bayr)
# library(rstanarm)
# library(brms)
# library(bayr)
# library(stringr)

GDRIVE = "c:/Users/martin/Google Drive/"
thisdir = getwd()
cases_dir = paste0(GDRIVE, "Aktenkoffer/Publications/New_Stats/Cases/")
setwd(cases_dir)
load("CUE8.Rda")
load("BrowsingAB.Rda")
setwd(thisdir)
load(paste0(GDRIVE, "PUBLICATIONS/Publicatie these Stefan Huijser/Data/Lap15.Rda"))
load(paste0(GDRIVE, "THESES/MT Daan Keeris/Data/DK.Rda"))



## Testing parameter splitting on a bunch of models

posterior(Lap15$M_2_Dur)

posterior(Lap15$M_1_Dur)

posterior(M_$MathurRepl_3)
parnames(M_$MathurRepl_3)

posterior(M_$MathurRepl_4)
parnames(M_$MathurRepl_3)

bayr:::extr_brms_par(CUE8$M_1)
posterior(CUE8$M_1)

## New method for fes and res

parnames(CUE8$M_1)
brms:::fixef.brmsfit(CUE8$M_1)

parnames(Lap15$M_3)
brms:::fixef.brmsfit(Lap15$M_3)

## Testing tbl_post.data.frame
## (recasting tbl_post)

P_recast <-
	posterior(M_$MathurRepl_3) %>%
	mutate(anno_1 = "x") %>%
	posterior()

class(P_recast)[1] == "tbl_post"
P_recast ## anno_1 in user annos

## Testing coef siblings

coef(Lap15$P_2_Dur)
coef(Lap15$M_2_Dur)
fixef(Lap15$M_2_Dur)
grpef(Lap15$M_2_Dur)
ranef(Lap15$M_2_Dur)

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



## Uncanny: underscores in var names



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


