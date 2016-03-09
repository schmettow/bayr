library(modeest)
library(dplyr)
#library(tidyr)
#library(stringr)

## dplyr is used with NSE, which gives "no visible binding for global variable errors"
utils::globalVariables(names = c("type", "parameter", "value", "new_name", "iter", "pattern"))

#' posterior extraction (c
#'
#' MCMC chains are  extracted from a Bayesian (regression) object
#' and returned as a posterior object, which is in long format
#' (chain, iter, parameter, value, type, order). Parameters are classified
#' as fixef, ranef, grpef and named after a common scheme.
#'
#' @usage tbl_post(model, ...)
#' @param model Bayesian model object
#' @param ... ignored
#' @return tbl_post object with MCMC chain in long format
#'
#' The MCMC chains are extracted from the model and stored in a
#' common format. Different to most internal representations, and
#' mcmc objects from package coda, the chains are stored in a
#' long table with the following columns:
#'
#'
#' @author Martin Schmettow
#' @import dplyr
#' @export


tbl_post <-
	function (model, ...) {
		UseMethod("tbl_post", model)
	}

#' @rdname tbl_post
#' @export

posterior <-
	function(model, ...){
		tbl_post(model, ...)
	}


#' @rdname tbl_post
#' @export

tbl_post.MCMCglmm <-
	function(model, ...) {
		# building a parameter catalogue
		parameters <-
			bind_rows(data_frame(parameter = model$X@Dimnames[[2]],
													 type = "fixef") %>%
									mutate(new_name = stringr::str_replace(parameter, "\\(Intercept\\)", "Intercept")),
								data_frame(parameter = colnames(model$VCV),
													 type = "grpef") %>%
									mutate(new_name = stringr::str_replace(parameter, "(.*)\\.(.*)", "\\2_\\1"),
												 new_name = stringr::str_replace(parameter, "units", "resid"))
			)
		## extract random effects if these are present
		if(!is.null(model$Z))	parameters <-
				parameters %>%
				bind_rows(data_frame(parameter = model$Z@Dimnames[[2]],
														 type = "ranef") %>%
										mutate(new_name = parameter))

		parameters <-
			parameters %>%
			mutate(order = row_number())

		fixed <-
			as.data.frame(model$Sol) %>%
			as_data_frame() %>%
			mutate(iter = row_number()) %>%
			tidyr::gather(parameter, value, -iter) %>%
			mutate(parameter = as.character(parameter))

		random <-
			as.data.frame(model$VCV) %>%
			as_data_frame() %>%
			mutate(iter = row_number()) %>%
			tidyr::gather(parameter, value, -iter) %>%
			mutate(value = sqrt(value),
						 parameter = as.character(parameter))

		out <-
			bind_rows(fixed, random) %>%
			mutate(chain = as.factor(1)) %>%
			full_join(parameters, by = "parameter") %>%
			select(chain, iter, parameter = new_name, value, type, order)

		class(out) <- append(class(out), "tbl_post")

		return(out)
	}


#' @rdname tbl_post
#' @export

tbl_post.brmsfit <-
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
			mutate(match = stringr::str_detect(parameter, pattern)) %>%
			filter(match) %>%
			select(parameter, type, pattern)

		out <-
			samples %>%
			mutate(iter = row_number()) %>%
			tidyr::gather(parameter, value, -iter, -chain) %>%
			mutate(parameter = as.character(parameter)) %>%
			full_join(type_mapping, by = "parameter") %>%
			full_join(par_order, by = "parameter") %>%
			arrange(order, chain, iter) %>%
			mutate(parameter = stringr::str_replace(parameter, "^sigma_(.*)", "sigma_resid")) %>%
			mutate(parameter = stringr::str_replace(parameter, pattern, "")) %>%
			select(-pattern)

		class(out) <-
			append(class(out), "tbl_post")

		return(out)
	}

