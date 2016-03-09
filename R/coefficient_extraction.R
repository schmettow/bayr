library(modeest)
library(dplyr)
library(tidyr)
library(stringr)

## dplyr is used with NSE, which gives "no visible binding for global variable errors"
utils::globalVariables(names = c("type", "parameter", "value", "new_name", "iter", "pattern",
																 "tbl_coef"))


################ COEF ###############################


#' Coefficient extraction
#'
#' summary table of fixed, random or group-level coefficients from posterior
#'
#' @param object tbl_post (brms, MCMCglmm) object holding the posterior in long format
#' @param type type of coefficient: fixef (grpef, ranef)
#' @param estimate function for computing the center estimate (posterior mode)
#' @param interval credibility interval: .95
#' @param ... ignored
#' @return coefficient table with parameter name, estimate and interval
#'
#' The standard center function is the posterior mode computed
#' by modeest::shorth
#'
#' @author Martin Schmettow
#' @import dplyr
#' @importFrom modeest shorth
#' @importFrom nlme fixef
#' @importFrom nlme ranef
#' @importFrom stats coef
#' @export

coef.tbl_post <-
	function(object, type = "fixef", estimate = shorth, interval = .95, ...) {
		lower <- (1-interval)/2
		upper <- 1-((1-interval)/2)
		partype <- type

		tbl_coef <-
			object %>%
			filter(type == partype) %>%
			group_by(parameter, order) %>%
			summarize(center = estimate(value),
								lower = quantile(value, lower),
								upper = quantile(value, upper)) %>%
			ungroup() %>%
			arrange(order) %>%
			select(-order)

		class(tbl_coef) <- append("tbl_coef", class(tbl_coef))
		attr(tbl_coef, "estimate") <- bquote(estimate)
		attr(tbl_coef, "interval") <- interval
		attr(tbl_coef, "lower") <- lower
		attr(tbl_coef, "upper") <- upper
		attr(tbl_coef, "type") <- type
		return(tbl_coef)
	}

#' @export


print.tbl_coef <-
	function(x, digits = NULL, title = T, footnote = T, ...){
		types <- data_frame(type = c("fixef", "ranef", "grpef"),
											 name = c("fixed effects","random effects", "group effects"))
		name <-
			types %>%
			filter(type == attr(x, "type")) %>%
			select(name) %>%
			as.character()

		if(title) cat(name, " coefficients", "\n***\n")
		print.data.frame(x, digits = digits, row.names = F)
		if(footnote) cat("\n*\nestimate with ",
										 attr(x, "interval")*100,
										 "% credibility limits")

		cat("\n")
		invisible(x)
	}


coef.MCMCglmm <-
	function(object, estimate = shorth, ...)
		tbl_post(object) %>% fixef(estimate = estimate, ...)

coef.brms <-
	function(object, estimate = shorth, ...)
		tbl_post(object) %>% fixef(estimate = estimate, ...)


################ FIXEF ###############################

#' @rdname coef.tbl_post
#' @export

fixef.tbl_post <-
	function(object, estimate = shorth, ...)
		coef(object, type = "fixef", estimate = estimate, ...)

#' @rdname coef.tbl_post
#' @export

fixef.MCMCglmm <-
	function(object, estimate = shorth, ...)
		tbl_post(object) %>% fixef(estimate = estimate, ...)


#' @rdname coef.tbl_post
#' @export

fixef.brmsfit <-
	function(object, estimate = shorth, ...)
		tbl_post(object) %>% fixef(estimate = estimate, ...)


############## RANEF ##############

#' @rdname coef.tbl_post
#' @export

ranef.tbl_post <-
	function(object, estimate = shorth, ...)
		coef(object, type = "ranef", estimate = estimate, ...)


#' @rdname coef.tbl_post
#' @export

ranef.MCMCglmm <-
	function(object, estimate = shorth, ...)
		tbl_post(object) %>% ranef(estimate = estimate)

#' @rdname coef.tbl_post
#' @export

ranef.brmsfit <-
	function(object, estimate = shorth, ...)
		tbl_post(object) %>% ranef(estimate = estimate)




##################### GRPEF ######################

grpef <- function(object, estimate = shorth, ...) UseMethod("grpef", object)

#' @rdname coef.tbl_post
#' @export

grpef.tbl_post <-
	function(object, estimate = shorth, ...)
		coef(object, type = "grpef", estimate = estimate, ...)


#' @rdname coef.tbl_post
#' @export

grpef.MCMCglmm <-
	function(object, estimate = shorth, ...)
		tbl_post(object) %>% ranef(estimate = estimate)

#' @rdname coef.tbl_post
#' @export

grpef.brmsfit <-
	function(object, estimate = shorth, ...)
		tbl_post(object) %>% ranef(estimate = estimate)


# TODO
# 1) support for the other stan
# 2) function resid()
# 3) function cor()
