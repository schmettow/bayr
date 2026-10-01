## Dev script: refit the brms/rstanarm models used as the test fixture
## (tests/testthat/M_1.Rda).
##
## Run manually from the package root:
##   Rscript data-raw/prepare_test_models.R
##
## Not part of the shipped package: data-raw/ is listed in .Rbuildignore.
## Requires brms, rstanarm and a working C++ toolchain (Stan compilation).
## The tests themselves never fit models; they only load the fixture.

if (requireNamespace("pkgload", quietly = TRUE)) {
	pkgload::load_all(".", quiet = TRUE)
} else {
	library(bayr)
}

library(brms)
library(rstanarm)

options(mc.cores = 2)

## GLMM

Ipump <- read.csv("data-raw/Pumps.csv")
Ipump <- Ipump[Ipump$Part <= 10 & Ipump$Task <= 5, ]
Ipump <- as_tbl_obs(Ipump)

F_1 <- ToT ~ Design * session + (1 + Design | Part) + (1 | Task)

M_1_b <- brm(F_1, family = "Gaussian", data = Ipump, chains = 2, iter = 200)
M_1_s <- stan_glmer(F_1, data = Ipump, chains = 2, iter = 200)

## Smoke checks of the extraction pipeline
print(fixef(M_1_b))
print(fixef(M_1_s))
print(clu(M_1_b))
print(clu(M_1_s))

save(M_1_b, M_1_s, file = "tests/testthat/M_1.Rda", compress = "xz")
