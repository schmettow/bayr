#library(tidyverse)

## dplyr is used with NSE, which gives "no visible binding for global variable errors"
#utils::globalVariables(names = c("type", "parameter", "value", "new_name", "iter", "pattern"))

Cols_pp = list("model", "Obs", "chain", "iter", "scale", "value") ## columns for a postpred object to be valid



#' posterior predictive extraction
#'
#' MCMC predicted values are  extracted from a Bayesian (regression) object
#' and returned as a tbl_post_pred object
#'
#' @usage post_pred(model, scale = "obs", model_name, thin = 1)
#' @param model Bayesian model object
#' @param newdata new data to predict from
#' @param scale "response" or "lin_pred"
#' @param model_name provides a name for the model
#' @param thin thinning factor
#' @return tbl_postpred object with MCMC draws
#'
#' chains are stored in a
#' long table with the following columns:
#'
#' chain  iter   Obs   value  type  order
#' (fctr) (int)  (int)       (dbl)  (chr) (int)
#'
#'
#' @author Martin Schmettow
#' @import dplyr
#' @export


post_pred <-
	function(model,
					 scale = "obs",
					 model_name = deparse(substitute(model)),
					 newdata = NULL,
					 thin = 1, ...){

		post_matrix <- mtx_post_pred(model,
																 newdata = newdata,
																 thin = thin, ...)
		out <- post_matrix %>%
			as.data.frame() %>%
			as_tibble() %>%
			mutate(model = model_name,
						 chain = 1, ## fixme
						 scale = "resp", ## fixme
						 iter = 1:n()) %>%
			pivot_longer(matches("\\d+"),
									 names_to = "Obs",
									 values_to = "value") %>%
			mutate(Obs = as.integer(stringr::str_replace(Obs, "^V", ""))) %>%
			select(model, Obs, chain, iter, scale, value) %>%
			arrange(model, Obs, chain, iter, scale)

		#assert_names(out, model, Obs, chain, iter, scale, value)
		class(out) <-
			append("tbl_post_pred", class(out))
		out
	}

#' @rdname post_pred
#' @export

mtx_post_pred <-
	function (model, ...) {
		UseMethod("mtx_post_pred", model)
	}


mtx_post_pred.data.frame <-
	## IDEA: write methods for
	## - identifying user_annos (all user annos)
	## - registering user annos (explicit user annos)
	## - keep attribute user_annos (keep user annos)
	function(df, model_name, thin = 1, ...) {
		if(! all(as.character(bayr:::AllCols) %in% names(df)))
			stop("not a valid tbl_post_pred, some columns missing")
		out = df
		class(out) = append("tbl_post_pred", class(out))
		out
	}

#' @rdname post_pred
#' @export

mtx_post_pred.brmsfit <-
	function(model, model_name, newdata = NULL, thin = 1, ...){
		n_iter <- brms::ndraws(model)
		n_draws <- round(n_iter/thin, 0)
		#draws <- sort(sample.int(n_iter, n_draws, replace = F))
		brms:::predict.brmsfit(model, newdata = newdata, nsamples = n_draws, summary = F)

	}


#' @rdname post_pred
#' @export


mtx_post_pred.stanreg <-
	function(model, model_name, newdata = NULL, thin = 1, ...){
		n_iter <- sum(model$stanfit@sim$n_save)
		n_draws <- round(n_iter/thin, 0)
		rstanarm::posterior_predict(model, newdata = newdata, draws = n_draws)
	}


# tbl_post_pred_old.generic <-
# 	function(sample_matrix, model_name = model_name){
# 		out <-
# 			sample_matrix %>%
# 			as.data.frame() %>%
# 			as_tibble() %>%
# 			mutate(model = model_name,
# 						 chain = 1, ## fixme
# 						 scale = "resp", ## fixme
# 						 iter = 1:n()) %>%
# 			pivot_longer(names_to = "Obs",
# 										value = "value", -model, -iter, -chain, -scale) %>%
# 			mutate(Obs = as.integer(stringr::str_replace(Obs, "^V", "")))  ## making observations integer (Vx)
#
# 		class(out) <-
# 			append("tbl_post_pred", class(out))
# 		out
# 	}
#
# ## Totaler Bullshit! Merged into post_pred.
# tbl_post_pred.generic <-
# 	function(sample_matrix, model_name = model_name){
# 		out <-
# 			sample_matrix %>%
# 			as.data.frame() %>%
# 			as_tibble() %>%
# 			mutate(model = model_name,
# 						 chain = 1, ## fixme
# 						 scale = "resp", ## fixme
# 						 iter = 1:n()) %>%
# 			pivot_longer(matches("^V\\d+"),
# 									 names_to = "Obs",
# 									 values_to = "value") %>%
# 			mutate(Obs = as.integer(stringr::str_replace(Obs, "^V", ""))) %>%
# 			select(model, Obs, chain, iter, scale, value) %>%
# 			arrange(model, Obs, chain, iter, scale)
# 		class(out) <-
# 			append("tbl_post_pred", class(out))
# 		out
# 	}
