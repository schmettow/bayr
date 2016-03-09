library(modeest)
library(dplyr)
library(tidyr)
library(stringr)

## dplyr is used with NSE, which gives "no visible binding for global variable errors"
utils::globalVariables(names = c("type", "parameter", "value", "new_name", "iter", "pattern"))


################ COEF ###############################


#' Coefficient extraction
#'
#' summary table of fixed, random or group-level coefficients from posterior
#'
#' @param object posterior distribution object (or MCMCglmm/brmsfit directly)
#' @param type type of effect default: fixef (grpef, ranef)
#' @param estimate function for computing the location estimate (posterior mode)
#' @param ... ignored
#' @return coefficient table with parameter name, location and CI
#'
#' The standard location function is the posterior mode computed
#' by modeest::shorth
#'
#' @author Martin Schmettow
#' @import dplyr
#' @importFrom modeest shorth
#' @importFrom nlme fixef
#' @importFrom nlme ranef
#' @importFrom stats coef
#' @export

coef.posterior <-
	function(object, type = "fixef", estimate = shorth, ...) {
		object %>%
		filter(type == type) %>%
		group_by(parameter, order) %>%
		summarize(location = estimate(value),
							"l-95% CI" = quantile(value, .025),
							"u-95% CI" = quantile(value, .975)) %>%
		ungroup() %>%
		arrange(order) %>%
		select(-order)
}

coef.MCMCglmm <-
	function(object, estimate = shorth, ...)
		posterior(object) %>% fixef(estimate = estimate, ...)

coef.brms <-
	function(object, estimate = shorth, ...)
		posterior(object) %>% fixef(estimate = estimate, ...)


################ COEF ###############################

#' @rdname coef.posterior
#' @export

fixef.posterior <-
	function(object, estimate = shorth, ...)
		coef(object, type = "fixef", estimate = estimate, ...)

#' @rdname coef.posterior
#' @export

fixef.MCMCglmm <-
	function(object, estimate = shorth, ...)
	posterior(object) %>% fixef(estimate = estimate, ...)


#' @rdname coef.posterior
#' @export

fixef.brmsfit <-
	function(object, estimate = shorth, ...)
	posterior(object) %>% fixef(estimate = estimate, ...)



############## RANEF ##############

#' @rdname coef.posterior
#' @export

ranef.posterior <-
	function(object, estimate = shorth, ...)
		coef(object, type = "ranef", estimate = estimate, ...)


#' @rdname coef.posterior
#' @export

ranef.MCMCglmm <-
	function(object, estimate = shorth, ...)
	posterior(object) %>% ranef(estimate = estimate)

#' @rdname coef.posterior
#' @export

ranef.brmsfit <-
	function(object, estimate = shorth, ...)
	posterior(object) %>% ranef(estimate = estimate)




##################### GRPEF ######################

grpef <- function(object, estimate = shorth, ...) UseMethod("grpef", object)

#' @rdname coef.posterior
#' @export

grpef.posterior <-
	function(object, estimate = shorth, ...)
		coef(object, type = "grpef", estimate = estimate, ...)


#' @rdname coef.posterior
#' @export

grpef.MCMCglmm <-
	function(object, estimate = shorth, ...)
	posterior(object) %>% ranef(estimate = estimate)

#' @rdname coef.posterior
#' @export

grpef.brmsfit <-
	function(object, estimate = shorth, ...)
	posterior(object) %>% ranef(estimate = estimate)


# TODO
# 1) support for the other stan
# 2) function resid()
# 3) function cor()
