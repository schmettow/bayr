library(modeest)
library(dplyr)
library(tidyr)
library(stringr)

## dplyr is used with NSE, which gives "no visible binding for global variable errors"
utils::globalVariables(names = c("type", "parameter", "value", "new_name", "iter", "pattern"))

#' MCMC chain extraction
#'
#' MCMC chains are  extracted from a Bayesian (regression) object
#' and returned as a posterior object, which is in long format
#' (chain, iter, parameter, value, type, order). Parameters are classified
#' as fixef, ranef, grpef and named after a common scheme.
#'
#' @usage posterior(model, ...)
#' @param model Bayesian model object
#' @param ... ignored
#' @return posterior object with MCMC chain in long format
#'
#' The MCMC chains are extracted from the model and stored in a
#' common format. Different to most internal representations, and
#' mcmc objects from package coda, the chains are stored in a
#' long table with the following columns:
#'
#'
#' @author Martin Schmettow
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @export


posterior <-
	function (model, ...) {
	UseMethod("posterior", model)
}

#' @describeIn posterior extraction from MCMCglmm
#' @export

posterior.MCMCglmm <-
	function(model, ...) {
	parameters <-
		bind_rows(data_frame(parameter = model$X@Dimnames[[2]],
												 type = "fixef") %>%
								mutate(new_name = str_replace(parameter, "\\(Intercept\\)", "Intercept")),
							data_frame(parameter = model$Z@Dimnames[[2]],
												 type = "ranef") %>%
								mutate(new_name = parameter),
							data_frame(parameter = colnames(model$VCV),
												 type = "grpef") %>%
								mutate(new_name = str_replace(parameter, "(.*)\\.(.*)", "\\2_\\1"),
											 new_name = str_replace(parameter, "units", "resid"))
		) %>%
		mutate(order = row_number())

	fixed <-
		as.data.frame(model$Sol) %>%
		as_data_frame() %>%
		mutate(iter = row_number()) %>%
		gather(parameter, value, -iter) %>%
		mutate(parameter = as.character(parameter))

	random <-
		as.data.frame(model$VCV) %>%
		as_data_frame() %>%
		mutate(iter = row_number()) %>%
		gather(parameter, value, -iter) %>%
		mutate(value = sqrt(value),
					 parameter = as.character(parameter))

	out <-
		bind_rows(fixed, random) %>%
		mutate(chain = as.factor(1)) %>%
		full_join(parameters, by = "parameter") %>%
		select(chain, iter, parameter = new_name, value, type, order)

	class(out) <- append(class(out), "posterior")

	return(out)
}


#' @describeIn posterior extraction from brmsfit
#' @export

posterior.brmsfit <-
	function(model, ...){
	samples <-
		brms::posterior_samples(model, add_chain = T) %>%
		as_data_frame()

	par_order <-
		data_frame(parameter = colnames(samples)) %>%
		mutate(order = row_number())

	type_patterns <-
		data_frame(pattern = c("^b_", "^sd_", "^r_", "^sigma_", "^lp__"),
							 type = c("fixef", "grpef", "ranef", "grpef", "diag"))
	type_mapping <-
		expand.grid(type = type_patterns$type,
								parameter = par_order$parameter) %>%
		mutate(parameter = as.character(parameter),
					 type = as.character(type)) %>%
		full_join(type_patterns, by = "type") %>%
		mutate(match = str_detect(parameter, pattern)) %>%
		filter(match) %>%
		select(parameter, type, pattern)

	out <-
		samples %>%
		mutate(iter = row_number()) %>%
		gather(parameter, value, -iter, -chain) %>%
		mutate(parameter = as.character(parameter)) %>%
		full_join(type_mapping, by = "parameter") %>%
		full_join(par_order, by = "parameter") %>%
		arrange(order, chain, iter) %>%
		mutate(parameter = str_replace(parameter, "^sigma_(.*)", "sigma_resid")) %>%
		mutate(parameter = str_replace(parameter, pattern, "")) %>%
		select(-pattern)

	class(out) <-
		append(class(out), "posterior")

	return(out)
}

################ FIXEF ###############################


#' Fixed effects
#'
#' Tabular summary of fixed effects coefficients from posterior
#'
#' @usage fixef(posterior, loc.func = shorth, ...)
#' @param posterior posterior distribution object (or MCMCglmm/brmsfit directly)
#' @param loc.func function for computing the location
#' @param ... ignored
#' @return coefficient table with parameter name, location and CI
#'
#' The standard location function is the posterior mode computed
#' by modeest::shorth
#'
#' @author Martin Schmettow
#' @import dplyr
#' @import tidyr
#' @importFrom modeest shorth
#' @export

fixef <-
	function (posterior, loc.func = shorth, ...) {
	UseMethod("fixef", posterior)
}

#' @describeIn fixef fixed effects extraction
#' @export


fixef.posterior <-
	function(posterior, loc.func = shorth, ...) {
	posterior %>%
		filter(type == "fixef") %>%
		group_by(parameter, order) %>%
		summarize(location = loc.func(value),
							"l-95% CI" = quantile(value, .025),
							"u-95% CI" = quantile(value, .975)) %>%
		ungroup() %>%
		arrange(order) %>%
		select(-order)
}

#' @describeIn fixef fixed effects extraction
#' @export


fixef.MCMCglmm <-
	function(posterior, loc.func = shorth, ...)
	posterior(posterior) %>% fixef(loc.func = loc.func, ...)

#' @describeIn fixef fixed effects extraction
#' @export

fixef.brmsfit <-
	function(posterior, loc.func = shorth, ...)
	posterior(posterior) %>% fixef(loc.func = loc.func, ...)



############## RANEF ##############

#' Random effects
#'
#' Tabular summary of random effects levels from posterior
#'
#' @usage ranef(posterior, loc.func = shorth, ...)
#' @param posterior posterior distribution object
#' @param loc.func function for computing the location
#' @param ... ignored
#' @return coefficient table with parameter name, location and CI
#'
#' The standard location function is the posterior mode computed
#' by modeest::shorth
#'
#' @author Martin Schmettow
#' @importFrom modeest shorth
#' @export

ranef <-
	function (posterior, loc.func = shorth, ...) {
	UseMethod("ranef", posterior)
}


#' @describeIn ranef random effects effects extraction
#' @export

ranef.posterior <-
	function(posterior, loc.func = shorth, ...) {
	posterior %>%
		filter(type == "ranef") %>%
		group_by(parameter, order) %>%
		summarize(location = loc.func(value),
							"l-95% CI" = quantile(value, .025),
							"u-95% CI" = quantile(value, .975)) %>%
		arrange(order) %>%
		select(-order)
}

#' @describeIn ranef random effects effects extraction
#' @export

ranef.MCMCglmm <-
	function(posterior, loc.func = shorth, ...)
	posterior(posterior) %>% ranef(loc.func = loc.func)

#' @describeIn ranef random effects effects extraction
#' @export

ranef.brmsfit <-
	function(posterior, loc.func = shorth, ...)
	posterior(posterior) %>% ranef(loc.func = loc.func)




##################### GRPEF ######################

#' Group-level effects
#'
#' Tabular summary of group level coefficients from posterior
#'
#' @usage grpef(posterior, loc.func = shorth, ...)
#' @param posterior posterior distribution object
#' @param loc.func function for computing the location
#' @param ... ignored
#' @return coefficient table with parameter name, location and CI
#'
#' The standard location function is the posterior mode computed
#' by modeest::shorth
#'
#' @author Martin Schmettow
#' @importFrom modeest shorth
#' @export

grpef <-
	function (posterior, loc.func = shorth, ...) {
	UseMethod("grpef", posterior)
}

#' @describeIn grpef group effects extraction
#' @export

grpef.posterior <-
	function(posterior, loc.func = shorth, ...) {
	posterior %>%
		filter(type == "grpef") %>%
		group_by(parameter, order) %>%
		summarize(location = loc.func(value),
							"l-95% CI" = quantile(value, .025),
							"u-95% CI" = quantile(value, .975))
}

#' @describeIn grpef group effects extraction
#' @export

grpef.MCMCglmm <-
	function(posterior, loc.func = shorth, ...)
	posterior(posterior) %>% ranef(loc.func = loc.func)

#' @describeIn grpef group effects extraction
#' @export

grpef.brmsfit <-
	function(posterior, loc.func = shorth, ...)
	posterior(posterior) %>% ranef(loc.func = loc.func)


# TODO
# 1) support for the other stan
# 2) function resid()
# 3) function cor()
