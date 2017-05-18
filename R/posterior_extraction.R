#library(dplyr)
#library(tidyr)
#library(stringr)

## dplyr is used with NSE, which gives "no visible binding for global variable errors"
utils::globalVariables(names = c("type", "parameter", "value",
																 "new_name", "iter", "pattern"))


ParameterIDCols = list("parameter", "type", "nonlin", "fixef", "re_factor", "re_entity")
AllCols = append(append(list("model", "chain", "iter", "order"), ParameterIDCols), "value")


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
#' @export


posterior <-
	function(model,
					 shape = "long",
					 thin = 1,
					 type = c("fixef", "grpef", "ranef", "disp", "cor"),
					 model_name = NA, ...){

		if (!(shape %in%  c("wide", "long"))) warning("shape must be either wide or long")

		if(is.na(model_name)) model_name <- deparse(substitute(model))

		post <-
			tbl_post(model, ...)

		if(! all(as.character(bayr:::AllCols) %in% names(post)))
			stop("not a valid tbl_post")

		type_regex <- stringr::str_c(type, collapse = "|")

		out <-
			post %>%
			filter(!iter %% thin) %>%
			filter(stringr::str_detect(type, type_regex)) %>%
			arrange(chain, iter, type, parameter)

		## this prohibits model overwriting (tbl_post.data.frame), needs testing
		if(any(is.na(out$model))) out$model = model_name

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
	## TODO: add formula and corr
	function(tbl_post, kable = by_knitr(), ...){
		n_iter <- length(unique(tbl_post$iter))
		n_chain <- length(unique(tbl_post$chain))

		user_annos <- setdiff(names(tbl_post),
													as.character(bayr:::AllCols))

		effects <-
			tbl_post %>%
			filter(type %in% c("fixef", "ranef", "grpef")) %>%
			distinct(model, parameter, type, fixef, nonlin, re_factor, re_entity) %>%
			mutate(parameter = ifelse(type == "ranef", "", parameter)) %>%
			group_by(model, parameter, type, nonlin, fixef, re_factor) %>%
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
			cat("tbl_post: ", n_iter, " samples in ", n_chain, " chains\n")
			cat("Effects: \n")
			print(knitr::kable(effects))
			cat("\nDispersion: \n")
			print(knitr::kable(disp))
		} else {

			cat("** tbl_post: ", n_iter, " samples in ", n_chain, " chains\n\n")
			#		cat(frm, "\n\n")

			cat("** Effects: \n")
			print.data.frame(effects, row.names = F)

			cat("\n** Dispersion: \n")
			print.data.frame(disp, row.names = F)

			cat("\n** User annotations: \n", user_annos)

		}
		invisible(tbl_post)
	}


#' @rdname posterior
#' @export



tbl_post.data.frame <-
	## IDEA: write methods for
	## - identifying user_annos (all user annos)
	## - registering user annos (explicit user annos)
	## - keep attribute user_annos (keep user annos)
	function(df, ...) {
		if(! all(as.character(bayr:::AllCols) %in% names(df))) stop("not a valid tbl_post, some columns missing")
		out = df
		class(out) = append("tbl_post", class(out))
		out
	}

## extracting some parameter annotations from brms model prior

extr_brms_par <-
	function(model){
		## use fixef and ranef parnames for check
		pn_fe <- rownames(brms:::fixef.brmsfit(model))
		try(pn_re <- names(brms:::ranef.brmsfit(model)), silent = T)


		pars <-
			model$prior %>%
			select(-prior, fixef = coef, re_factor = group, nonlin = nlpar) %>%
			filter(class %in% c("b", "sd"))

		pars$type <-
			case_when(
				pars$class == "b"  ~ "fixef",
				pars$class == "sd" ~ "grpef")

		pars$parameter <-
			case_when(
				pars$class == "b" & pars$nonlin == "" ~
					stringr::str_c(pars$class, "_", pars$fixef),
				pars$class == "b" & pars$nonlin != "" ~
					stringr::str_c(pars$class, "_",pars$nonlin,"_", pars$fixef),
				pars$class == "sd" & pars$nonlin == "" ~
					stringr::str_c(pars$class, "_", pars$re_factor, "__", pars$fixef),
				pars$class == "sd" & pars$nonlin != "" ~
					stringr::str_c(pars$class, "_", pars$re_factor, "__", pars$nonlin,"_", pars$fixef)
			)

		## hack for missing b_Intercept (CUE8)

		if("Intercept" %in% pn_fe) {
			pars <- pars %>% bind_rows(
				data_frame(
					parameter = "b_Intercept",
					type = "fixef",
					nonlin =NA,
					fixef = "Intercept",
					re_factor = NA,
					re_entity = NA
				)) %>%
			select_(.dots = bayr:::ParameterIDCols)
		}

		pars %>%
			filter(!stringr::str_detect(parameter, "_$")) %>%
			mutate(re_entity = NA) %>%
			mutate_all(funs(ifelse(. == "", NA, .))) %>%
			distinct() %>%
			select_(.dots = bayr:::ParameterIDCols)


	}



#' @rdname posterior
#' @export

tbl_post.brmsfit <-
	function(model, ...){

		samples <-
			brms::posterior_samples(model, add_chain = T) %>%
			as_data_frame()

		samples_long <-
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
														 "^sigma", "^shape$", "^phi$", "^lp__", "^cor_"),
								 type = c("fixef", "fixef", "grpef",
								 				 "ranef", "disp", "disp", "disp","diag", "cor"))

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
			# mutate(type = ifelse(!is.na(brms_type), brms.type, type)) %>%
			# select(-pattern) %>%
			distinct() %>%
			# type.y is from brms_extr
			#mutate(type = if_else(is.na(type.y), type.x, NA_character_)) #%>%
			select_(.dots = bayr:::ParameterIDCols)
			# mutate(parameter = stringr::str_replace(parameter,
			# 																				"^sigma_(.*)", "sigma_resid")) %>%
			# mutate(parameter = stringr::str_replace(parameter,
			# 																				"^lp__(.*)", "lp__log_density")) %>%
			# mutate(parameter = stringr::str_replace(parameter, pattern, "")) %>%

		## we assume there is always at least one fixef,
		## otherwise would be creating an empty par_samples_long to start with

		# if(any(par_all$type == "fixef")){
		# 	par_fe <-
		# 		par_all %>%
		# 		filter(type == "fixef") %>%
		# 		tidyr::extract(parameter,
		# 									 into = c("nonlin","fixef"),
		# 									 "^b_(?:(.+)_)?+(.+)$",
		# 									 fill = "left",
		# 									 remove = F) %>%
		# 		mutate(re_factor = NA,
		# 					 re_entity = NA) %>%
		# 		select_(.dots = ParameterIDCols)
		# 	par_out <- par_fe
		# }
		par_out <-
			filter(par_all, type %in% c("fixef", "grpef"))

		#if(any(par_all$type == "disp", na.rm = T)){
			par_disp <-
				par_all %>%
				filter(type == "disp") %>%
				mutate(nonlin = NA,
							 fixef = NA,
							 re_factor = NA,
							 re_entity = NA) %>%
				select_(.dots = ParameterIDCols)
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
											 fill = "left",
											 remove = F) %>%
				## hotfix for when there is no nonlin
				mutate(re_factor = ifelse(is.na(re_factor), re_1, re_factor)) %>%
				select_(.dots = ParameterIDCols)
			par_out <- bind_rows(par_out, par_re)
		}

		# if(any(par_all$type == "grpef")){
		# 	par_ge <-
		# 		par_all %>%
		# 		filter(type == "grpef") %>%
		# 		tidyr::extract(parameter,
		# 									 into = c("re_factor","coef"),
		# 									 "^sd_(.+)__(.+)$",
		# 									 remove = F) %>%
		# 		tidyr::extract(coef,
		# 									 into = c("nonlin","fixef"),
		# 									 "^(?:(.*)_)?(.+)$",
		# 									 fill = "right",
		# 									 remove = F) %>%
		# 		mutate(re_entity = NA) %>%
		# 		select_(.dots = ParameterIDCols)
		# 	par_out <- bind_rows(par_out, par_ge)
		# }


		## TODO: correlations,
		## gives strange error "missing values where T/F needed"
		## with Mathur_Repl_4

		# if(any(par_all$type == "cor", na.rm = T)){
			par_cor <-
				par_all %>%
				filter(type == "cor") %>%
				mutate(nonlin = NA,
							 fixef = NA,
							 re_factor = NA,
							 re_entity = NA) %>%
				select_(.dots = ParameterIDCols)
			par_out <- bind_rows(par_out, par_cor)
		#}

		# if(any(par_all$type == "diag", na.rm = T)){
			par_diag <-
				par_all %>%
				filter(type == "diag") %>%
				mutate(nonlin = NA,
							 fixef = NA,
							 re_factor = NA,
							 re_entity = NA) %>%
				select_(.dots = ParameterIDCols)
			par_out <- bind_rows(par_out, par_diag)
		#}

		# par_out <-
		# 	left_join(par_all,
		# 						par_out,
		# 						by = c("parameter", "type")) %>%
		par_out <-
			par_out %>%
			distinct_(.dots = ParameterIDCols) %>%
			full_join(par_order, by = "parameter")




		## joining with samples

		out <-
			full_join(samples_long, par_out, by = "parameter") %>%
			mutate(model = NA) %>% #,
			#re_factor = stringr::str_replace(re_factor, "_$", ""),
			#nonlin = stringr::str_replace(nonlin, "_$", "")) %>%
			select_(.dots = AllCols)



		class(out) <-
			append("tbl_post", class(out))

		return(out)
	}



tbl_post_old.brmsfit <-
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
														 "^sigma", "^shape$", "^phi$", "^lp__", "^cor_"),
								 type = c("fixef", "fixef", "grpef",
								 				 "ranef", "disp", "disp", "disp","diag", "cor"))

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
		# model <- M_cue8_rstn
		samples <-
			rstanarm:::as.data.frame.stanreg(model) %>%
			as_data_frame() %>%
			mutate(chain = NA,
						 iter = row_number()) %>%
			tidyr::gather("parameter", "value", -chain, -iter) %>%
			mutate(parameter = stringr::str_replace(parameter, "\\(Intercept\\)", "Intercept"),
						 parameter = stringr::str_replace(parameter, "sigma", "sigma_resid"))

		parnames <-
			samples %>%
			select(parameter) %>%
			distinct() %>%
			mutate(type = case_when(
				stringr::str_detect(.$parameter, "sigma_resid") ~ "disp",
				stringr::str_detect(.$parameter, "^b\\[.*\\]") ~ "ranef",
				stringr::str_detect(.$parameter, "^Sigma\\[.*\\]") ~ "grpef",
				TRUE ~ "fixef"),
				order = row_number())

		par_fe <-
			parnames %>%
			filter(type == "fixef") %>%
			mutate(nonlin = NA,
						 fixef = parameter,
						 re_factor = NA,
						 re_entity = NA) %>%
			select_(.dots = bayr:::ParameterIDCols)

		par_disp <-
			parnames %>%
			filter(type == "disp") %>%
			mutate(nonlin = NA,
						 fixef = NA,
						 re_factor = NA,
						 re_entity = NA) %>%
			select_(.dots = bayr:::ParameterIDCols)

		par_out <- bind_rows(par_fe, par_disp)
		if("lmerMod" %in% class(model)) {
			par_re <-
				parnames %>%
				filter(type == "ranef") %>%
				tidyr::extract(parameter,
											 into = c("fixef", "re_factor", "re_entity"),
											 "^b\\[(.+) (.+):(.+)\\]$",
											 remove = F) %>%
				mutate(nonlin = NA) %>%
				select_(.dots = bayr:::ParameterIDCols)


			par_ge <-
				parnames %>%
				filter(type == "grpef") %>%
				tidyr::extract(parameter,
											 into = c("re_factor", "fixef"),
											 "^Sigma\\[(.+):(.+),(.+)\\]$",
											 remove = F) %>%
				mutate(nonlin = NA,
							 re_entity = NA) %>%
				select_(.dots = bayr:::ParameterIDCols)
			par_out <- bind_rows(par_out, par_re, par_ge)
		}
			par_out <-
				par_out %>%
				right_join(parnames, by = c("parameter","type"))


		## putting it back together
		## and transforming var to sd
		out <-
			par_out %>%
			right_join(samples, by = "parameter") %>%
			mutate(model = NA) %>%
			mutate(value = ifelse(type == "grpef", sqrt(value), value)) %>%
			select_(.dots = bayr:::AllCols)


		class(out) <-
			append("tbl_post", class(out))

		return(out)
	}




tbl_post.stanreg.old <-
	function(model, ...){

		## classifying parameters lm
model <- M_cue8_rstn


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
				mutate(type = case_when(
					stringr::str_detect(parameter, "^b\\[.*\\]") ~ "ranef",
					stringr::str_detect(parameter, "^Sigma\\[.*\\]") ~ "grpef",
					TRUE ~ type))

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
			filter(type == "disp") %>%
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


		if("lmerMod" %in% class(model)) {
			par_ge <-
				par_all %>%
				filter(type == "grpef") %>%
				tidyr::extract(parameter,
											 into = c("re_entity", "fixef"),
											 "^Sigma\\[(.+):(.+),(.+)\\]$",
											 remove = F) %>%
				mutate(nonlin = NA) %>%
				select_(.dots = ParameterIDCols)

			par_out <-
				bind_rows(par_out, par_re, par_ge)
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
