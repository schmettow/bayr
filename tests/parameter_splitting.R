
library(dplyr)
library(rstanarm)
library(brms)

cases_dir = "../../../../Publications/New_Stats/Cases/"
load(paste0(cases_dir,"CUE8.Rda"))
load("../../../../../PUBLICATIONS/Publicatie these Stefan Huijser/Data/M_2_Dur_corr.Rda")

summary(CUE8$M_1)
summary(CUE8$M_1_rstanarm)


P_1_brms <- posterior(CUE8$M_1)
P_1_brms %>%
	select_(.dots = AllCols)


P_1_rstn <- posterior(CUE8$M_1_rstanarm)
P_1_rstn %>%
	select_(.dots = AllCols)


grpef(P_1_brms)

P_1_rstn

P_1_brms %>%
	distinct(type, nonlin, fixef, re_factor, re_unit) %>%
	group_by(type, nonlin, fixef, re_factor) %>%
	summarize(units = n())

P_1_rstn %>%
	distinct(type, nonlin, fixef, re_factor, re_unit) %>%
	group_by(type, nonlin, fixef, re_factor) %>%
	summarize(units = n())


## Lap15: non-linear


P_2_brms <-
	posterior(Lap15$M_2_Dur)


## Testing posterior()
rm("posterior")

P_3_rstn_p <- bayr::posterior(CUE8$M_1_rstanarm)
P_3_rstn_t <- bayr:::tbl_post.stanreg(CUE8$M_1_rstanarm)

## The ultimate test case for brms
## is the Lap15 data with z-transformed preds, as it contains:
## nonlin, fixef, ranef, unit and correlations between ranef

P_4 <- posterior(M_2_Dur_corr)
unique(P_4$type)

