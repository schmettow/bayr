library(modeest)
library(dplyr)
library(tidyr)
library(stringr)

## dplyr is used with NSE, which gives "no visible binding for global variable errors"
utils::globalVariables(names = c("type", "parameter", "value", "new_name", "iter", "pattern","tbl_coef"))


################ COEF ###############################


#' Coefficient extraction
#'
#' summary table of fixed, random or group-level coefficients and fitted values (eta,
#' only stanfit models) from posterior
#'
#' @param object tbl_post (brms, MCMCglmm) object holding the posterior in long format
#' @param model model
#' @param type type of coefficient: fixef (grpef, ranef)
#' @param mean function (identity)
#' @param estimate function for computing the center estimate (posterior mode)
#' @param interval credibility interval: .95
#' @param ... ignored
#' @return coefficient table with parameter name, estimate and interval
#'
#' The standard center function is the posterior median
#'
#' @author Martin Schmettow
#' @import dplyr
#' @importFrom nlme fixef
#' @importFrom nlme ranef
#' @importFrom stats coef median fitted quantile
#' @importFrom knitr knit_print
#' @export

coef.tbl_post <-
	function(object,
					 model = unique(object$model),
					 type = c("fixef", "disp", "grpef", "shape"), ## maybe deprecate for ~filter
					 mean.func = identity,
					 estimate = median,
					 interval = .95) {
		lower <- (1-interval)/2
		upper <- 1-((1-interval)/2)
		partype <- type

		tbl_coef <-
			object %>%
			filter(type %in% partype) %>%
			mutate(value = mean.func(value)) %>%
			group_by(model, parameter, type, nonlin, fixef, re_factor, re_entity, order) %>%
			summarize(center = estimate(value),
								lower = quantile(value, lower),
								upper = quantile(value, upper)) %>%
			ungroup() %>%
			arrange(model, order) %>%
			select(-order)

		class(tbl_coef) <- append("tbl_coef", class(tbl_coef))
		attr(tbl_coef, "estimate") <- bquote(estimate)
		attr(tbl_coef, "interval") <- interval
		attr(tbl_coef, "lower") <- lower
		attr(tbl_coef, "upper") <- upper
		attr(tbl_coef, "type") <- type
		return(tbl_coef)
	}



# coef <-
#  	function(object, estimate = median, ...) UseMethod("coef", object)



#' @rdname coef.tbl_post
#' @export

coef.data.frame <-
	## dirty hack! partype columns not preserved
	function(df, ...) {
		if(! all(c("parameter", "center", "lower", "upper") %in% names(df))) stop("not a valid tbl_coef, some columns missing")
		out = df
		class(out) = append("tbl_coef", class(out))
		out
	}


#' @rdname coef.tbl_post
#' @export

coef.MCMCglmm <-
	function(object, estimate = median, ...)
		tbl_post(object) %>% coef(estimate = estimate, ...)

#' @rdname coef.tbl_post
#' @export

coef.brmsfit <-
	function(object, estimate = median, ...)
		tbl_post(object) %>% coef(estimate = estimate, ...)

#' @rdname coef.tbl_post
#' @export

coef.stanfit <-
	function(object, estimate = median, ...)
		tbl_post(object) %>% coef(estimate = estimate, ...)

#' @rdname coef.tbl_post
#' @export

coef.stanreg <-
	function(object, estimate = median, ...)
		tbl_post(object) %>% coef(estimate = estimate, ...)


################ FIXEF ###############################

#' @rdname coef.tbl_post
#' @export

fixef <-
	function(object, estimate = median, ...) UseMethod("fixef", object)

#' @rdname coef.tbl_post
#' @export

fixef.tbl_post <-
	function(object, estimate = median, ...)
		coef(object, type = "fixef", estimate = estimate, ...) %>%
	select(-parameter)

#' @rdname coef.tbl_post
#' @export

fixef.MCMCglmm <-
	function(object, estimate = median, ...)
		tbl_post(object) %>% fixef(estimate = estimate, ...)


#' @rdname coef.tbl_post
#' @export

fixef.brmsfit <-
	function(object, estimate = median, ...)
		tbl_post(object) %>% fixef(estimate = estimate, ...)


#' @rdname coef.tbl_post
#' @export

fixef.stanfit <-
	function(object, estimate = median, ...)
		tbl_post(object) %>% fixef(estimate = estimate, ...)

#' @rdname coef.tbl_post
#' @export

fixef.stanreg <-
	function(object, estimate = median, ...)
		posterior(object) %>% fixef(estimate = estimate, ...)

############## RANEF ##############

#' @rdname coef.tbl_post
#' @export


ranef <-
	function(object, estimate = median, ...) UseMethod("ranef", object)



#' @rdname coef.tbl_post
#' @export

ranef.tbl_post <-
	function(object, estimate = median, ...)
		coef(object, type = "ranef", estimate = estimate, ...) %>%
	select(-parameter)


#' @rdname coef.tbl_post
#' @export

ranef.MCMCglmm <-
	function(object, estimate = median, ...)
		tbl_post(object) %>% ranef(estimate = estimate)

#' @rdname coef.tbl_post
#' @export

ranef.brmsfit <-
	function(object, estimate = median, ...)
		tbl_post(object) %>% ranef(estimate = estimate)


#' @rdname coef.tbl_post
#' @export

ranef.stanfit <-
	function(object, estimate = median, ...)
		tbl_post(object) %>% ranef(estimate = estimate)

#' @rdname coef.tbl_post
#' @export

ranef.stanreg <-
	function(object, estimate = median, ...)
		tbl_post(object) %>% ranef(estimate = estimate)


##################### GRPEF ######################


#' @rdname coef.tbl_post
#' @export


grpef <-
	function(object, estimate = median, ...) UseMethod("grpef", object)

#' @rdname coef.tbl_post
#' @export

grpef.tbl_post <-
	function(object, estimate = median, ...)
		coef(object, type = "grpef", estimate = estimate, ...) %>%
	select(-parameter)


#' @rdname coef.tbl_post
#' @export

grpef.MCMCglmm <-
	function(object, estimate = median, ...)
		tbl_post(object) %>% grpef(estimate = estimate)

#' @rdname coef.tbl_post
#' @export

grpef.brmsfit <-
	function(object, estimate = median, ...)
		tbl_post(object) %>% grpef(estimate = estimate)


#' @rdname coef.tbl_post
#' @export

grpef.stanfit <-
	function(object, estimate = median, ...)
		tbl_post(object) %>% grpef(estimate = estimate)


#' @rdname coef.tbl_post
#' @export

grpef.stanreg <-
	function(object, estimate = median, ...)
		tbl_post(object) %>% grpef(estimate = estimate)

#################### PRINT #######################

#' @rdname coef.tbl_post
#' @export

print.tbl_coef <- function(x, ...) {
	tab <- x
	if(nrow(tab) > 1)	{
		tab <- mascutils::discard_redundant(tab)
	} else if(tab$fixef[1] == "Intercept"){
		#		warning("Intercept model")
		tab <- select(tab, fixef, center, lower, upper)
	} else {
		tab <- mascutils::discard_all_na(x)
	}
	cap <-
		paste0("Estimates with ",
					 attr(x, "interval")*100,
					 "% credibility limits")
	out <-
		knitr::kable(tab, caption = cap)
	print(out)
	invisible(tab)
}

#' @rdname coef.tbl_post
#' @export

knit_print.tbl_coef <- function(x, ...) {
	tab <- x
	if(nrow(tab) > 1)	{
		tab <- mascutils::discard_redundant(tab)
	} else if(tab$fixef[1] == "Intercept"){
		#		warning("Intercept model")
		tab <- select(tab, fixef, center, lower, upper)
	} else {
		tab <- mascutils::discard_all_na(x)
	}
	cap <-
		paste0("Estimates with ",
					 attr(x, "interval")*100,
					 "% credibility limits")
	out <-	knitr::kable(tab, caption = cap)
	print(out)
	invisible(tab)
}




print.tbl_coef_EATME <-
	function(x, digits = NULL, title = F, footnote = T,
					 kable = bayr:::by_knitr()){
		out <- mascutils::discard_all_na(x)
		if(nrow(out) > 1)	{
			out <- mascutils::discard_redundant(out)}
		else{
			out <- select(out, -model, -type)
		}
		types <-
			data_frame(type = c("fixef",
													"ranef",
													"grpef",
													"fitted"),
								 title_text = c("fixed effects",
								 							 "random effects",
								 							 "factor-level variation (sd)",
								 							 "fitted values (linear predictor)"))
		title_text <- ""
		if(title)
			title_text <-
			types %>%
			filter(type %in% attr(x, "type")) %>%
			select(title_text) %>%
			as.character()

		footnote_text <-
			paste0("\n*\nestimate with ",
						 attr(x, "interval")*100,
						 "% credibility limits")

		if(kable) { ## prepared for knitr table output, not yet working
			print(knitr::kable(out, caption = title_text))
		} else {
			if(title) cat(stringr::str_c(title_text, sep = "|", collapse = T), "\n***\n")
			print.data.frame(out, digits = 3, row.names = F)
			if(footnote) cat(footnote_text)
			cat("\n")
			invisible(out)
		}
		return(out)
	}





########################### FUTURE ##########################

#' Joining two coefficent tables (NI)
#'
#' NOT IMPLEMENTED
#' creates an overview table of point estimates for models with overlapping
#' parameter sets
#'
#' @param first coefficient table tbl_coef or tbl_compcoef
#' @param second coefficient table tbl_coef
#' @param modelnames list of model names
#' @return data_frame with parameter names and point estimates
#'
#' Returns a comparative coefficient table (tbl_coefcomp).
#' Models are shown in columns, Parameters ordered by the most complex model
#' Models ordered by number of parameters (nesting?)
#' Model names are taken via lazy eval, or can be given.
#'
#' @author Martin Schmettow

join.tbl_coef <-
	function(x, y) {
		error("Not implemented")
		out <- data_frame()
		class(out) <- c(class(out), "tbl_coefcomp")
	}



join.tbl_coefcomp <-
	function(x, y) {
		error("Not implemented")
		out <- data_frame()
		class(out) <- c(class(out), "tbl_coefcomp")
	}


#' separating factors an group means model
#'
#' NOT IMPLEMENTED
#' In a model with interaction effects only (group means contrasts)
#' coefficient names are split into factor variables (and cleaned)
#'
#' @param first coefficient table tbl_coef or tbl_compcoef
#' @param second coefficient table tbl_coef
#' @param modelnames list of model names
#' @return data_frame with two or more factor variables, center and CI.
#'
#' Returns a data_frame with parameter plit into factors and cleaned.
#' Prepares for intercation plots
#'
#' @author Martin Schmettow

seperate.tbl_coef <-
	function(x, y) {
		error("Not implemented")
		out <- data_frame()
		class(out) <- c("tbl_df")
	}

#' @rdname coef.tbl_post
#' @export



# resid_plot_1 <-
# 	function(object){
# 		data_frame(fitted = fitted(object)[,1],
# 							 residual = residuals(object)[,1]) %>%
# 			ggplot(aes(x = fitted, y = residual)) +
# 			geom_quantile()
#
# 	}



# TODO
# 1) support for the other stan
# 2) function resid()
# 3) function cor()
