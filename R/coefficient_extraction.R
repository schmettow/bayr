#library(tidyverse)
#
# ## dplyr is used with NSE, which gives "no visible binding for global variable errors"
# utils::globalVariables(names = c("type", "parameter", "value", "new_name", "iter", "pattern","tbl_coef"))





################ CLU ###############################

#' center-lower-upper table
#'
#' returns a summary table for estimates with median (center) and 95% limits'
#' @param object tbl_post (brms, MCMCglmm) object holding the posterior in long format
#' @param model name of model
#' @param mean.func function (identity)
#' @param estimate function for computing the center estimate (posterior mode)
#' @param interval credibility interval: .95
#' @param ... ignored
#' @return CLU table (tbl_clu) with parameter name, estimate and interval.
#' @author Martin Schmettow
#' @import dplyr
#' @import lme4
#' @import broom.mixed
#' @import assertthat
#' @importFrom nlme fixef
#' @importFrom nlme ranef
#' @importFrom stats coef median fitted quantile
#' @importFrom knitr knit_print
#' @export



clu <-
  	function(object, ...) UseMethod("clu", object)

#' @rdname clu
#' @export


assert_clu <- function(x){
	assert_names(x, model, parameter, type, center, lower, upper)
	assert_key(x, model, parameter)
}

#' @rdname clu
#' @export


clu.tbl_post <- function(object,
												 model = unique(object$model),
												 mean.func = identity,
												 estimate = median,
												 interval = .95){
	lower <- (1-interval)/2
	upper <- 1-((1-interval)/2)
	tbl_clu <-
		object %>%
		mutate(value = mean.func(value)) %>%
		group_by(model, parameter, type, nonlin, fixef, re_factor, re_entity, order) %>%
		summarize(center = estimate(value),
							lower = quantile(value, lower),
							upper = quantile(value, upper)) %>%
		ungroup() %>%
		arrange(model, order) %>%
		select(-order)
	class(tbl_clu) <- append("tbl_clu", class(tbl_clu))
	attr(tbl_clu, "estimate") <- bquote(estimate)
	attr(tbl_clu, "interval") <- interval
	attr(tbl_clu, "lower") <- lower
	attr(tbl_clu, "upper") <- upper
	return(tbl_clu)
}




# clu.tbl_post.old <-
# 	function(object,
# 					 model = unique(object$model),
# 					 type = NA,
# 					 mean.func = identity,
# 					 center = median,
# 					 interval = .95) {
# 		lower <- (1-interval)/2
# 		upper <- 1-((1-interval)/2)
#
# 		if(!is.na(type)) object <- filter(object, type %in% partype)
#
# 		tbl_clu <-
# 			object %>%
# 			mutate(value = mean.func(value)) %>%
# 			group_by(model, parameter, type, nonlin, fixef, re_factor, re_entity, order) %>%
# 			summarize(center = center(value),
# 								lower = quantile(value, lower),
# 								upper = quantile(value, upper)) %>%
# 			ungroup() %>%
# 			arrange(model, order) %>%
# 			select(-order)
#
# 		class(tbl_clu) <- append("tbl_clu", class(tbl_clu))
# 		attr(tbl_clu, "estimate") <- bquote(estimate)
# 		attr(tbl_clu, "interval") <- interval
# 		attr(tbl_clu, "lower") <- lower
# 		attr(tbl_clu, "upper") <- upper
# 		attr(tbl_clu, "type") <- type
# 		return(tbl_clu)
# 	}



# clu_1.tbl_post <- function(tbl_post){
# 	out <-
# 		tbl_post %>%
# 		group_by(parameter) %>%
# 		summarize(center = median(value),
# 							lower = quantile(value, .05, na.rm = T),
# 							upper = quantile(value, .95, na.rm = T))
# 	class(out) = append("tbl_clu", class(out))
# 	out
# }


#' @rdname clu
#' @export

clu.data.frame <-
	## dirty hack! partype columns not preserved
	function(df, ...) {
		assert_clu(df)
		out = df
		class(out) = append("tbl_clu", class(out))
		out
	}


#' @rdname clu
#' @export

clu.MCMCglmm <-
	function(object, ...)
		tbl_post(object) %>% clu()

#' @rdname clu
#' @export

clu.brmsfit <-
	function(object, ...)
		tbl_post(object) %>% clu()

#' @rdname clu
#' @export

clu.stanfit <-
	function(object, ...)
		tbl_post(object) %>% clu()

#' @rdname clu
#' @export

clu.stanreg <-
	function(object, ...)
		tbl_post(object) %>% clu()

#' @rdname clu
#' @export

clu.glmerMod <-
	function(model,
					 mean.func = identity,
					 model_name = as.character(deparse(substitute(model))),
					 interval = .95){
		lower <- (1-interval)/2
		upper <- 1-(1-interval)/2

		pop_level <-
			model %>%
			tidy(conf.int = T, conf.level = ) %>%
			rename(parameter = term,
						 re_factor = group,
						 type = effect,
						 center = estimate,
						 lower = conf.low,
						 upper = conf.high) %>%
			mutate(model = model_name,
						 center = mean.func(center),
						 lower = mean.func(lower),
						 upper = mean.func(upper),
						 type = if_else(type == "fixed",
						 							 "fixef",
						 							 str_extract(parameter,
						 							 						"^cor|^sd")))

		ranefs <-
			nlme::ranef(model) %>%
			map_dfr(as_tibble, rownames = "re_entity", .id = "re_factor") %>%
			pivot_longer(!starts_with("re_"),
									 names_to = "fixef",
									 values_to = "center") %>%
			mutate(model = model_name,
						 parameter = str_c(re_factor,"_",re_entity,".",fixef),
						 type = "ranef",
						 lower = NA,
						 upper = NA)

		out <-
			bind_rows(pop_level,
								ranefs) %>%
			filter(!is.na(center))

		class(out) <- c("tbl_clu", class(out))
		attr(out, "interval") <- interval
		out
	}


################ COEF ###############################


#' Coefficient extraction
#'
#' summary table of fixed, random or group-level coefficients and fitted values (eta,
#' only stanfit models) from posterior
#'
#' @param object tbl_post (brms, MCMCglmm) object holding the posterior in long format
#' @param model model
#' @param type type of coefficient: fixef (grpef, ranef)
#' @param mean.func function (identity)
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
					 type = c("fixef", "ranef"), ## maybe deprecate for ~filter
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



################ FIXEF MULTILEVEL ###############################

#' Multi-level coefficient table
#'
#' produces a CLU table of population-level effects together with
#' random effects standard deviation for all random factor levels.
#'
#' @param object tbl_post (brms, rstanarm) object holding the posterior in long format
#' @param model model
#' @param ... passed on to grpef and fixef
#' @return coefficient table with standard deviations
#'
#'
#' @author Martin Schmettow
#' @import dplyr
#' @export


fixef_ml <-
	function(object, ...){
		model <- posterior(object)
		T_grpef <-
			grpef(model, ...) %>%
			select(model, fixef, re_factor, SD = center) %>%
			mutate(re_factor = str_c("SD_", re_factor)) %>%
			spread(re_factor, SD)

		out <-
			bayr::fixef(model, ...)  %>%
			left_join(T_grpef, by = c("model", "fixef")) %>%
			mascutils::discard_redundant()
		class(out) <- append("tbl_fixef_ml", class(out))
		out
	}







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
