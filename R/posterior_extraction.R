#library(dplyr)
#library(tidyr)
#library(stringr)

## dplyr is used with NSE, which gives "no visible binding for global variable errors"
utils::globalVariables(names = c("type", "parameter", "value",
																 "new_name", "iter", "pattern"))


ParameterIDCols = list("parameter", "type", "fixef", "nonlin", "re_factor", "re_entity")
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
	function(model,
					 shape = "long",
					 thin = 1,
					 type = c("fixef", "grpef", "ranef", "disp", "cor"), ...){
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
			stop("Wide format currently not implemented")
			out <-
				out %>%
				mutate(parameter = str_c(model, nonlin, fixef, re_factor, re_entity, sep = "_")) %>%
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

print.tbl_post <-
	## TODO: refactor code to work with brms:parnames
	function(tbl_post, kable = by_knitr(), ...){
		n_iter <- length(unique(tbl_post$iter))
		n_chain <- length(unique(tbl_post$chain))
		effects <-
			tbl_post %>%
			filter(type %in% c("fixef", "ranef", "grpef")) %>%
			distinct(model, type, fixef, nonlin, re_factor, re_entity) %>%
			group_by(model, type, nonlin, fixef, re_factor) %>%
			summarize(entities = n()) %>%
			ungroup() %>%
			mascutils::discard_all_na()

		# corr <-
		# 	tbl_post %>%
		# 	filter(type == "corr") %>%
		# 	distinct(type, re_factor, fixef_1, fixef_2,ranef_1, ranef_2)

		disp <-
			filter(tbl_post, type == "disp") %>%
			distinct(parameter) %>%
			mascutils::discard_all_na()

		# frm <-
		# 	formula.tools:::as.character.formula(attr(tbl_post, "formula"))

		if(kable){ ## prepared for knitr table output, not yet working
			cat("tbl_post: ", n_iter, " samples in ", n_chain, " chains\n\n")
			cat("Effects: \n")
			print(knitr::kable(effects))
			cat("\nDispersion: \n")
			print(knitr::kable(disp))
		} else {
			cat("tbl_post: ", n_iter, " samples in ", n_chain, " chains\n\n")
			#		cat(frm, "\n\n")
			cat("Effects: \n")
			print.data.frame(effects, row.names = F)
			cat("\nDispersion: \n")
			print.data.frame(disp, row.names = F)
		}
		invisible(tbl_post)
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

tbl_post.brmsfit <-
	function(model, ...){
		samples <-
			brms::posterior_samples(model, add_chain = T) %>%
			as_data_frame()

		out <-
			samples %>%
			mutate(iter = row_number()) %>%
			tidyr::gather(parameter, value, -iter, -chain) %>%
			mutate(parameter = as.character(parameter))

		par_order <-
			data_frame(parameter = colnames(samples)) %>%
			mutate(order = row_number()) %>%
			filter(!parameter %in% c("chain", "iter"))

		type_patterns <-
			data_frame(pattern = c("^b.","^b_", "^sd_", "^r_",
														 "^sigma$", "^shape$", "^lp__", "^cor_"),
								 type = c("fixef", "fixef", "grpef",
								 				 "ranef", "disp", "disp", "diag", "cor"))

		type_mapping <-
			expand.grid(type = unique(type_patterns$type),
									parameter = par_order$parameter) %>%
			mutate(parameter = as.character(parameter),
						 type = as.character(type)) %>%
			full_join(type_patterns, by = "type") %>%
			mutate(match = stringr::str_detect(parameter, pattern)) %>%
			filter(match) %>%
			select(parameter, type, pattern)

		par_all <-
			out %>%
			distinct(parameter) %>%
			full_join(type_mapping, by = "parameter") %>%
			full_join(par_order, by = "parameter") %>%
			# mutate(parameter = stringr::str_replace(parameter,
			# 																				"^sigma_(.*)", "sigma_resid")) %>%
			# mutate(parameter = stringr::str_replace(parameter,
			# 																				"^lp__(.*)", "lp__log_density")) %>%
			# mutate(parameter = stringr::str_replace(parameter, pattern, "")) %>%
			select(-pattern)

		## we assume there is always at least one fixef,
		## otherwise would be creating an empty par_out to start with

		if(any(par_all$type == "fixef")){
			par_fe <-
				par_all %>%
				filter(type == "fixef") %>%
				tidyr::extract(parameter,
											 into = c("nonlin","fixef"),
											 "^b_(?:(.+)_)?+(.+)$",
											 fill = "left",
											 remove = F) %>%
				mutate(re_factor = NA,
							 re_entity = NA) %>%
				select_(.dots = ParameterIDCols)
			par_out <- par_fe
		}

		if(any(par_all$type == "disp")){
			par_disp <-
				par_all %>%
				filter(type == "disp") %>%
				mutate(nonlin = NA,
							 fixef = NA,
							 re_factor = NA,
							 re_entity = NA) %>%
				select_(.dots = ParameterIDCols)
			par_out <- bind_rows(par_out, par_disp)
		}

		if(any(par_all$type == "ranef")){
			par_re <-
				par_all %>%
				filter(type == "ranef") %>%
				tidyr::extract(parameter,
											 into = c("re_1", "re_entity","fixef"),
											 "^r_(.*)\\[(.+),(.+)\\]$",
											 remove = F) %>%
				tidyr::extract(re_1,
											 into = c("re_factor","nonlin"),
											 "^(.*)(?:__(.+))$",
											 fill = "left",
											 remove = F) %>%
				select_(.dots = ParameterIDCols)
			par_out <- bind_rows(par_out, par_re)
		}

		if(any(par_all$type == "grpef")){
			par_ge <-
				par_all %>%
				filter(type == "grpef") %>%
				tidyr::extract(parameter,
											 into = c("re_factor","coef"),
											 "^sd_(.+)__(.+)$",
											 remove = F) %>%
				tidyr::extract(coef,
											 into = c("nonlin","fixef"),
											 "^(?:(.*)_)?(.+)$",
											 fill = "right",
											 remove = F) %>%
				mutate(re_entity = NA) %>%
				select_(.dots = ParameterIDCols)
			par_out <- bind_rows(par_out, par_ge)
		}


		## TODO: correlations

		if(any(par_all$type == "cor")){
			par_cor <-
				par_all %>%
				filter(type == "cor") %>%
				mutate(nonlin = NA,
							 fixef = NA,
							 re_factor = NA,
							 re_entity = NA) %>%
				select_(.dots = ParameterIDCols)
			par_out <- bind_rows(par_out, par_cor)
		}

		if(any(par_all$type == "diag")){
			par_diag <-
				par_all %>%
				filter(type == "diag") %>%
				mutate(nonlin = NA,
							 fixef = NA,
							 re_factor = NA,
							 re_entity = NA) %>%
				select_(.dots = ParameterIDCols)
			par_out <- bind_rows(par_out, par_diag)
		}

		par_out <-
			left_join(par_all,
								par_out,
								by = c("parameter", "type")) %>%
			distinct()



		## putting it back together

		out <-
			full_join(out, par_out, by = "parameter") %>%
			mutate(model = NA) %>% #,
						 #re_factor = stringr::str_replace(re_factor, "_$", ""),
						 #nonlin = stringr::str_replace(nonlin, "_$", "")) %>%
			select_(.dots = AllCols)



		class(out) <-
			append("tbl_post", class(out))


		return(out)
	}



#' @rdname posterior
#' @export


tbl_post.stanreg <-
	function(model, ...){

		## classifying parameters lm

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

		## classifying parameters lmer

		if("lmerMod" %in% class(model)) {
			out <-
				out %>%
				mutate(type = ifelse(stringr::str_detect(parameter, "^b\\[.*\\]"),
														 "ranef", type))

		}

		## annotating parameters lm

		par_all <-
			out %>%
			distinct(parameter, type)


		par_fe <-
			par_all %>%
			filter(type == "fixef") %>%
			mutate(nonlin = NA,
						 fixef = parameter,
						 re_factor = NA,
						 re_entity = NA) %>%
			select_(.dots = ParameterIDCols)

		par_disp <-
			par_all %>%
			filter(type == "dispa") %>%
			mutate(nonlin = NA,
						 fixef = NA,
						 re_factor = NA,
						 re_entity = NA) %>%
			select_(.dots = ParameterIDCols)

		par_out <- bind_rows(par_fe, par_disp)

		## annotating parameters lmer

		if("lmerMod" %in% class(model)) {
			par_re <-
				par_all %>%
				filter(type == "ranef") %>%
				tidyr::extract(parameter,
											 into = c("fixef", "re_factor", "re_entity"),
											 "^b\\[(.+) (.+):(.+)\\]$",
											 remove = F) %>%
				mutate(nonlin = NA) %>%
				select_(.dots = ParameterIDCols)

			par_out <-
				bind_rows(par_out, par_re)
		}



		## putting it back together

		out <-
			par_out %>%
			right_join(out, by = c("parameter", "type")) %>%
			mutate(model = NA) %>%
			select_(.dots = AllCols)


		class(out) <-
			append("tbl_post", class(out))

		return(out)
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
												 																"entities", "resid"))
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
