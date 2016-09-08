#library(dplyr)
#library(tidyr)
#library(stringr)

## dplyr is used with NSE, which gives "no visible binding for global variable errors"
utils::globalVariables(names = c("type", "parameter", "value",
																 "new_name", "iter", "pattern"))


ParameterIDCols = list("parameter", "type", "fixef", "nonlin", "re_factor", "re_unit")
AllCols = append(append(list("model", "chain", "iter", "parameter", "order", "type"), ParameterIDCols), "value")


#' posterior extraction
#'
#' MCMC chains are  extracted from a Bayesian (regression) object
#' and returned as a posterior object, which is in long format
#' (chain, iter, parameter, value, type, order). Parameters are classified
#' as fixef, ranef, grpef and named after a common scheme.
#'
#' @usage posterior(model, shape, ...)
#' @param model Bayesian model object
#' @param shape return tbl_post in long shape or tbl_df in wide shape
#' @param thin thinning factor
#' @param type select parameter types to store ("fixef", "grpef", "ranef")
#' @return tbl_post object with MCMC chain in long format or tbl.df in wide shape
#'
#' The MCMC chains are extracted from the model and stored in a
#' common format. Different to most internal representations, and
#' mcmc objects from package coda, the chains are stored in a
#' long table with the following columns:
#'
#' chain  iter   parameter   value  type  order
#' (fctr) (int)  (chr)       (dbl)  (chr) (int)
#'
#'
#' @author Martin Schmettow
#' @import dplyr
#' @export


posterior <-
	function(model = CUE8$M_1_rstanarm,
					 shape = "long",
					 thin = 1,
					 type = c("fixef", "grpef", "ranef"), ...){
		if (!(shape %in%  c("wide", "long"))) warning("shape must be either wide or long")

		model_name <- deparse(substitute(model))

		post <-
			tbl_post(model, ...)

		type_regex <- stringr::str_c(type, collapse = "|")

		out <-
			post %>%
			filter(!iter %% thin) %>%
			filter(stringr::str_detect(type, type_regex)) %>%
			arrange(chain, iter, type, parameter) %>%
			mutate(model = model_name) %>%
			select_(.dots = AllCols)

		class(out) <- append("tbl_post", class(out))
		attr(out, which = "formula") <- formula(model)

		if (shape == "wide"){
			out <-
			out %>%
			#mutate(parameter = str_c(type, parameter, sep = "_")) %>%
			select(-order, -type) %>%
			tidyr::spread(parameter, value)
		}

		## Attention: Attributes are not preserved!
		return(out)
	}


#' @rdname posterior
#' @export

tbl_post <-
	function (model, ...) {
		UseMethod("tbl_post", model)
	}

#' @rdname posterior
#' @export

print.tbl_post_ <-
	function(x){
		n_iter <- length(unique(x$iter))
		n_chain <- length(unique(x$chain))
		fixef <-
			filter(x, type == "fixef") %>%
			distinct(parameter)
		grpef <-
			filter(x, type == "grpef") %>%
			distinct(parameter)
		ranef <-
			filter(x, type == "ranef") %>%
			distinct(parameter)
		disp <-
			filter(x, type == "disp") %>%
			distinct(parameter)
		frm <-
			formula.tools:::as.character.formula(attr(x, "formula"))
		cat("posterior samples in ", n_chain, " chains with ",
				n_iter, " samples each\n")
		cat(frm, "\n\n")
		cat("fixed: ", fixef$parameter, "\n")
		cat("group: ", grpef$parameter, "\n")
#		cat("ranef: ", ranef$parameter, "\n")
		cat("disp: ", disp$parameter, "\n")
		cat("\n")
		invisible(x)
	}


#' @rdname posterior
#' @export



tbl_post.data.frame <-
	function(model, ...) {
		out <- select_(model, dots = AllCols)
		class(out) <- append("tbl_post", class(out))
		out
	}



#' @rdname posterior
#' @export

tbl_post.MCMCglmm <-
	function(model, ...) {
		warning("MCMCglmm is out-of-sync: parameter splitting not implemented. Class tbl_post_v1")
		# building a parameter catalogue
		parameters <-
			bind_rows(data_frame(parameter = model$X@Dimnames[[2]],
													 type = "fixef") %>%
									mutate(new_name = stringr::str_replace(parameter,
																												 "\\(Intercept\\)", "Intercept")),
								data_frame(parameter = colnames(model$VCV),
													 type = "grpef") %>%
									mutate(new_name = stringr::str_replace(parameter,
																												 "(.*)\\.(.*)", "\\2_\\1"),
												 new_name = stringr::str_replace(parameter,
												 																"units", "resid"))
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

		class(out) <-
			append("tbl_post_v1", class(out))


		return(out)
	}


#' @rdname posterior
#' @export

tbl_post.brmsfit <-
	function(model, ...){
		samples <-
			brms::posterior_samples(model, add_chain = T) %>%
			as_data_frame()

		par_order <-
			data_frame(parameter = colnames(samples)) %>%
			mutate(order = row_number()) %>%
			filter(!parameter %in% c("chain", "iter"))

		type_patterns <-
			data_frame(pattern = c("^b.","^b_", "^sd_", "^r_", "^sigma", "^lp__", "^cor_"),
								 type = c("fixef", "fixef", "grpef", "ranef", "disp", "diag", "cor"))

		type_mapping <-
			expand.grid(type = unique(type_patterns$type),
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
			mutate(parameter = stringr::str_replace(parameter, "^lp__(.*)", "lp__log_density")) %>%
			mutate(parameter = stringr::str_replace(parameter, pattern, "")) %>%
			select(-pattern) %>%
			distinct() ### <--- dirty hack. it seems that fixed-effects models iterations get doubled

		## separating parameters

		out_id <-
			out %>%
			distinct(parameter, type)

		out_re <-
			out_id %>%
			filter(type == "ranef") %>%
			tidyr::extract(parameter,
										 into = c("nonlin","re_factor", "re_unit", "fixef"),
										 "^(.*_){0,1}(.+)\\[(.+),(.+)\\]$",
										 remove = F) %>%
			select_(.dots = ParameterIDCols)

		out_fe <-
			out_id %>%
			filter(type == "fixef") %>%
			mutate(nonlin = NA,
						 re_factor = NA,
						 re_unit = NA,
						 fixef = parameter) %>%
			select_(.dots = ParameterIDCols)

		out_ge <-
			out_id %>%
			filter(type == "grpef") %>%
			tidyr::extract(parameter,
										 into = c("nonlin","re_factor", "fixef"),
										 "^(.*_){0,1}(.+)_(.+)$",
										 remove = F) %>%
			mutate(re_unit = NA) %>%
			select_(.dots = ParameterIDCols)



		## putting it back together

		out <-
			bind_rows(out_fe, out_re) %>%
			bind_rows(out_ge) %>%
			right_join(out, by = c("parameter", "type")) %>%
			mutate(model = NA) %>%
			select_(.dots = AllCols)



		class(out) <-
			append("tbl_post", class(out))


		return(out)
	}

#' @rdname posterior
#' @export


tbl_post.stanfit <-
	function(model, ...){
		samples <-
			rstan::extract(model) %>%
			as.data.frame() %>%
			as_data_frame() #%>%
			# select(dplyr:::starts_with("b_"),
			# 			 dplyr:::starts_with("sd_"),
			# 			 dplyr:::starts_with("r_"),
			# 			 dplyr:::starts_with("sigma_"),
			# 			 dplyr:::starts_with("lp__"),
			# 			 dplyr:::starts_with("cor_"),
			# 			 dplyr:::starts_with("eta."))

		par_order <-
			data_frame(parameter = colnames(samples)) %>%
			mutate(order = row_number()) %>%
			filter(!parameter %in% c("chain", "iter"))

		type_patterns <-
			data_frame(pattern = c("^b_", "^b.","^sd_", "^r_", "^sigma_",
														 "^lp__", "^cor_", "eta."),
								 type = c("fixef", "fixef", "grpef", "ranef", "disp",
								 				 "diag","cor", "fitted"))

		type_mapping <-
			expand.grid(type = unique(type_patterns$type),
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
			tidyr::gather(parameter, value, -iter) %>%
			mutate(parameter = as.character(parameter)) %>%
			full_join(type_mapping, by = "parameter") %>%
			full_join(par_order, by = "parameter") %>%
			arrange(order, iter) %>%
			mutate(parameter = stringr::str_replace(parameter,
																							"^sigma_(.*)", "sigma_resid")) %>%
			mutate(parameter = stringr::str_replace(parameter,
																							"^lp__(.*)", "lp__log_density")) %>%
			mutate(parameter = stringr::str_replace(parameter, pattern, "")) %>%
			select(-pattern)

		class(out) <-
			append("tbl_post", class(out))

		return(out)
	}



#' @rdname posterior
#' @export


tbl_post.stanreg <-
	function(model, ...){
		if("lm" %in% class(model)) {
			samples <-
				rstanarm:::as.data.frame.stanreg(model) %>%
				as_data_frame()

			out <-
				samples %>%
				mutate(chain = NA,
							 iter = row_number()) %>%
				tidyr::gather("parameter", "value", -chain, -iter) %>%
				left_join(data_frame(parameter = colnames(samples)) %>%
										mutate(order = row_number(),
													 type = ifelse(parameter %in% c("sigma","shape"),
													 							"disp", "fixef")),
									by = "parameter") %>%
				mutate(parameter = stringr::str_replace(parameter,
																								"\\(Intercept\\)", "Intercept"),
							 parameter = stringr::str_replace(parameter,
							 																 "sigma", "sigma_resid"))

		}

		if("lmerMod" %in% class(model)) {
			out <-
				out %>%
				mutate(type = ifelse(stringr::str_detect(parameter, "^b"),
														 "ranef", type))

		}

		## separating parameters

		out_id <-
			out %>%
			distinct(parameter, type)

		out_re <-
			out_id %>%
			filter(type == "ranef") %>%
			tidyr::extract(parameter,
										 into = c("fixef", "re_factor", "re_unit"),
										 "^b\\[(.+) (.+):(.+)\\]$",
										 remove = F) %>%
			mutate(nonlin = NA) %>%
			select_(.dots = ParameterIDCols)

		out_fe <-
			out_id %>%
			filter(type == "fixef") %>%
			mutate(nonlin = NA,
						 fixef = parameter,
						 re_factor = NA,
						 re_unit = NA) %>%
			select_(.dots = ParameterIDCols)

		## putting it back together

		out <-
			bind_rows(out_fe, out_re) %>%
			right_join(out, by = c("parameter", "type")) %>%
			mutate(model = NA) %>%
			select_(.dots = AllCols)


		class(out) <-
			append("tbl_post", class(out))

		return(out)
	}




