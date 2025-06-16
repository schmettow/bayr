#library(tidyverse)

## dplyr is used with NSE, which gives "no visible binding for global variable errors"
#utils::globalVariables(names = c("type", "parameter", "value",
#																 "new_name", "iter", "pattern"))

ParameterIDCols = c("parameter", "type", "nonlin", "fixef", "re_factor", "re_entity")
AllCols = c("model", "chain", "iter", "order", ParameterIDCols, "value")



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
#' @param model_name provides a name for the model
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
#' @importFrom knitr knit_print
#' @export


posterior <-
	function(model,
					 shape = "long",
					 thin = 1,
					 type = c("fixef", "grpef", "ranef", "disp", "shape","cor"),
					 model_name = NA_character_, ...){


		if (!(shape %in%  c("wide", "long"))) warning("shape must be either wide or long")

		if(is.na(model_name)) model_name <- deparse(substitute(model))

		post <-
			tbl_post(model, ...) %>%
			as_tibble()

		if(! all(as.character(AllCols) %in% names(post)))
			stop("not a valid tbl_post")

		type_regex <- stringr::str_c(type, collapse = "|")

		out <-
			post %>%
			filter(!iter %% thin) %>%
			filter(stringr::str_detect(type, type_regex)) %>%
			arrange(chain, iter, type, parameter) %>%
			dplyr::mutate(model = model,
						 parameter = parameter,
						 type = type,
						 fixef = fixef,
						 re_factor = re_factor,
						 re_entity = re_entity)

		## this prohibits model overwriting (tbl_post.data.frame), needs testing
		if(any(is.na(out$model))) out$model = model_name

		class(out) <- append("tbl_post", class(out))
		attr(out, which = "formula") <- formula(model)

		if (shape == "wide"){
			stop("Wide format currently not implemented")
			out <-
				out %>%
				dplyr::mutate(parameter = stringr::str_c(model, nonlin, fixef, re_factor, re_entity, sep = "_")) %>%
				dplyr::select(-order, -type) %>%
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

assert_tbl_post <-
	function (x, ...) {
		assert_names(x, model, chain, iter, parameter, type, value, order)
		assert_key(x, model, chain, iter, parameter)
	}


# tidy_samples.brmsfit <-	function(model, ...){
# 	model = M_1_b
#
# 	samples <-
# 		tidybayes::tidy_draws(model) %>%
# 		gather(parameter, value, -.chain, -.iteration, -.draw)
# }



# annotate.brmsfit <- function(model){
#
# 	pars <- insight::find_parameters(model)
# 	annotations <- list()
#
# 	annotations <-
# 		tibble(parameter = pars$conditional,
# 									type = "fixef",
# 									nonlin = NA_character_,
# 									re_factor = NA_character_,
# 									re_entity = NA_character_) %>%
# 		mutate(fixef = stringr::str_remove(parameter, "^b_"),
# 					 fixef = stringr::str_replace(fixef, ".", ":")) %>%
# 		select(all_of(ParameterIDCols))
#
# 	ranef <-
# 		tibble(parameter = pars$random,
# 					 type = "ranef",
#
# 					 nonlin = NA_character_) %>%
# 		filter(stringr::str_detect(parameter, "^r_")) %>%
# 		tidyr::separate(col = parameter,
# 										into = c("r_", "re_factor", "re_entity", "fixef"),
# 										remove = F,
# 										extra = "merge") %>%
# 		mutate(fixef = stringr::str_remove(fixef, ".$")) %>%
# 		select(all_of(ParameterIDCols))
#
# 	disp <-
#
# 		bind_rows(fixef, ranef)
# 		out
# }


#' @rdname posterior
#' @export


tbl_post.data.frame <-
	## IDEA: write methods for
	## - identifying user_annos (all user annos)
	## - registering user annos (explicit user annos)
	## - keep attribute user_annos (keep user annos)
	function(x, ...) {
		assert_tbl_post(x)
		out <- x
		class(out) <- append("tbl_post", class(out))
		out
	}

## extracting some parameter annotations from brms model prior


extr_brms_par <-
	function(model){

		# model <- M_tot

		## use fixef and ranef parnames for check
		pn_fe <- rownames(brms:::fixef.brmsfit(model))
		try(pn_re <- names(brms:::ranef.brmsfit(model)), silent = T)


		pars <-
			model$prior %>%
			select(-prior, fixef = coef, re_factor = group, nonlin = nlpar, dist = dpar) %>%
			filter(class %in% c("b", "sd"))

		pars$type <-
			case_when(
				pars$class == "b" &! pars$dist == "" ~ "shape",
				pars$class == "b"                    ~ "fixef",
				pars$class == "sd"                   ~ "disp")

		pars$parameter <-
			case_when(
				pars$class == "b" & pars$nonlin == ""  ~
					stringr::str_c(pars$class, "_", pars$fixef),
				pars$class == "b" & pars$nonlin != ""  ~
					stringr::str_c(pars$class, "_",pars$nonlin,"_", pars$fixef),
				pars$class == "sd" & pars$nonlin == "" ~
					stringr::str_c(pars$class, "_", pars$re_factor, "__", pars$fixef),
				pars$class == "sd" & pars$nonlin != "" ~
					stringr::str_c(pars$class, "_", pars$re_factor, "__", pars$nonlin,"_", pars$fixef)
			)

		## hack for missing b_Intercept (CUE8)

		if("Intercept" %in% pn_fe) {
			pars <- pars %>% bind_rows(
				tibble(
					parameter = "b_Intercept",
					type = "fixef",
					nonlin =NA_character_,
					fixef = "Intercept",
					re_factor = NA_character_,
					re_entity = NA_character_
				)) %>%
				select(all_of(ParameterIDCols))
		}

		pars %>%
			filter(!stringr::str_detect(parameter, "_$")) %>%
			mutate(re_entity = NA_character_) %>%
			#mutate_all(funs(ifelse(. == "", NA_character_, .))) %>%
			mutate(across(where(is.character), ~ ifelse(.x == "", NA_character_, .x))) %>%
			distinct() %>%
			select(all_of(ParameterIDCols))


	}



#' @rdname posterior
#' @export

tbl_post.brmsfit <-
	function(model, ...){

		#model = M_1_b
		#model = M_tot

		samples <-
			brms::as_draws_df(model) %>%
			as_tibble() %>%
			rename(iter = .iteration,
						 chain = .chain) %>%
			select(-.draw)

		samples_long <-
			samples %>%
			tidyr::gather(parameter, value, -iter, -chain) %>%
			mutate(parameter = as.character(parameter))

		par_order <-
			tibble(parameter = colnames(samples)) %>%
			mutate(order = row_number()) %>%
			filter(!parameter %in% c("chain", "iter"))

		type_patterns <-
			tribble(~pattern,   ~type,
							#"^b.",      "fixef",
							"^b_",      "fixef",
							"^sd_",     "grpef",
							"^r_",      "ranef",
							"^sigma",   "disp",
							"^shape$",  "shape",
							"^beta$",   "shape",
							"^phi$",    "shape",
							"^lp__",    "mcmc",
							"^cor_",    "cor")

		type_mapping <-
			expand.grid(type = unique(type_patterns$type),
									parameter = par_order$parameter) %>%
			mutate(parameter = as.character(parameter),
						 type = as.character(type)) %>%
			full_join(type_patterns, by = "type") %>%
			mutate(match = stringr::str_detect(parameter, pattern)) %>%
			filter(match) %>%
			select(parameter, type, pattern)

		brms_pars <-
			extr_brms_par(model) %>%
			select(-type)

		par_all <-
			samples_long %>%
			distinct(parameter) %>%
			full_join(type_mapping, by = "parameter") %>%
			full_join(brms_pars, by = "parameter") %>%
			distinct() %>%
			select(all_of(ParameterIDCols))

		par_out <-
			filter(par_all, type %in% c("fixef", "disp", "grpef"))

		#if(any(par_all$type == "disp", na.rm = T)){
		par_disp <-
			par_all %>%
			filter(type == "disp") %>%
			mutate(nonlin = NA_character_,
						 fixef = NA_character_,
						 re_factor = NA_character_,
						 re_entity = NA_character_) %>%
			select(all_of(ParameterIDCols))
		par_out <- bind_rows(par_out, par_disp)
		#}

		if(any(par_all$type == "ranef", na.rm = T)){
			par_re <-
				par_all %>%
				select(parameter, type) %>%
				filter(type == "ranef") %>%
				tidyr::extract(parameter,
											 into = c("re_1", "re_entity","fixef"),
											 "^r_(.*)\\[(.+),(.+)\\]$",
											 remove = F) %>%
				tidyr::extract(re_1,
											 into = c("re_factor","nonlin"),
											 "^(.+)(?:__(.+))$",
											 #fill = "left",
											 remove = F) %>%
				## hotfix for when there is no nonlin
				mutate(re_factor = ifelse(is.na(re_factor), re_1, re_factor)) %>%
				select(all_of(ParameterIDCols))
			par_out <- bind_rows(par_out, par_re)
		}

			par_cor <-
				par_all %>%
				filter(type == "cor") %>%
				mutate(nonlin = NA_character_,
							 fixef = NA_character_,
							 re_factor = NA_character_,
							 re_entity = NA_character_) %>%
				select(all_of(ParameterIDCols))
			par_out <- bind_rows(par_out, par_cor)

			par_diag <-
				par_all %>%
				filter(type == "diag") %>%
				mutate(nonlin = NA_character_,
							 fixef = NA_character_,
							 re_factor = NA_character_,
							 re_entity = NA_character_) %>%
				select(all_of(ParameterIDCols))
			par_out <- bind_rows(par_out, par_diag)

			par_shape <-
				par_all %>%
				filter(type == "shape") %>%
				mutate(nonlin = NA_character_,
							 fixef = NA_character_,
							 re_factor = NA_character_,
							 re_entity = NA_character_) %>%
				select(all_of(ParameterIDCols))
			par_out <- bind_rows(par_out, par_shape)

		par_out <-
			par_out %>%
			distinct(parameter, type, nonlin, fixef, re_factor, re_entity) %>%
			full_join(par_order, by = "parameter")




		## joining with samples

		out <-
			full_join(samples_long, par_out, by = "parameter") %>%
			mutate(model = NA_character_) %>% #,
			select(all_of(AllCols))



		class(out) <-
			append("tbl_post", class(out))

		return(out)
	}






#' @rdname posterior
#' @export


tbl_post.stanreg <-
	function(model, ...){

#		model <- M_1_s

		samples <-
			as.data.frame(model) %>%
			as_tibble() %>%
			mutate(chain = NA_integer_,
						 iter = row_number()) %>%
			tidyr::gather("parameter", "value", -chain, -iter) %>%
			mutate(parameter = stringr::str_replace(parameter, "\\(Intercept\\)", "Intercept"),
						 parameter = stringr::str_replace(parameter, "sigma", "sigma_resid"))

		parnames <-
			samples %>%
			select(parameter) %>%
			distinct() %>%
			mutate(type = case_when(
				stringr::str_detect(.$parameter, "sigma_resid|reciprocal_dispersion") ~ "disp",
				stringr::str_detect(.$parameter, "^b\\[.*\\]") ~ "ranef",
				stringr::str_detect(.$parameter, "^Sigma\\[.*\\]") ~ "grpef",
				stringr::str_detect(.$parameter, "^shape") ~ "shape",
				TRUE ~ "fixef"),
				order = row_number())

		par_fe <-
			parnames %>%
			filter(type == "fixef") %>%
			mutate(nonlin = NA_character_,
						 fixef = parameter,
						 re_factor = NA_character_,
						 re_entity = NA_character_) %>%
			select(all_of(ParameterIDCols))

		par_disp <-
			parnames %>%
			filter(type == "disp") %>%
			mutate(nonlin = NA_character_,
						 fixef = NA_character_,
						 re_factor = NA_character_,
						 re_entity = NA_character_) %>%
			select(all_of(ParameterIDCols))

		par_shape <-
			parnames %>%
			filter(type == "shape") %>%
			mutate(nonlin = NA_character_,
						 fixef = NA_character_,
						 re_factor = NA_character_,
						 re_entity = NA_character_) %>%
			select(all_of(ParameterIDCols))

		par_out <- bind_rows(par_fe, par_disp, par_shape)

		if("lmerMod" %in% class(model)) {
			par_re <-
				parnames %>%
				filter(type == "ranef") %>%
				tidyr::extract(parameter,
											 into = c("fixef", "re_factor", "re_entity"),
											 "^b\\[(.+) (.+):(.+)\\]$",
											 remove = F) %>%
				mutate(nonlin = NA_character_) %>%
				select(all_of(ParameterIDCols))


			par_ge <-
				parnames %>%
				filter(type == "grpef") %>%
				tidyr::extract(parameter,
											 into = c("re_factor", "fixef", "fixef_2"),
											 "^Sigma\\[(.+?):(.+),(.+)\\]$",
											 remove = F) %>%
				## pulling variance and correlations apart
				mutate(fixef_2 = if_else(fixef_2 == "(Intercept)", "Intercept", fixef_2),
							 type = if_else(fixef == fixef_2, "grpef", "corr")) %>%
				mutate(nonlin = NA_character_,
							 re_entity = NA_character_) %>%
				select(all_of(c(ParameterIDCols, "fixef_2")))
			par_out <- bind_rows(par_out, par_re, par_ge)
		}

		par_out <-
			parnames %>%
			select(-type) %>%  ## necessary, because above, type is changed (grpef, corr)
			left_join(par_out, by = c("parameter"))


		## putting it all together
		out <-
			par_out %>%
			right_join(samples, by = "parameter") %>%
			mutate(model = NA_character_) %>%
			update_by(type == "grpef",
													 value = sqrt(value)) %>%
			select(all_of(AllCols))

		# ## splitting by type for post-processing
		#
		# splt_type <-
		# 	split(samples, samples$type)
		#
		# ## transforming var to sd
		#
		# splt_type$grpef <-
		# 	splt_type$grpef %>%
		# 	mutate(value = sqrt(value))
		#
		# ## putting types together
		#
		# samples <-
		# 	bind_rows(splt_type)



		class(out) <-
			append("tbl_post", class(out))

		return(out)
	}


