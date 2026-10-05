

Cols_pp = list("model", "Obs", "chain", "iter", "scale", "value") ## columns for a postpred object to be valid



#' posterior predictive extraction
#'
#' MCMC predicted values are extracted from a Bayesian (regression) object
#' and returned as a tidy tbl_post_pred object.
#'
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
#' @export


post_pred <-
	function(model,
					 scale = "resp",
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
						 scale = scale, ## fixme
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



#' @rdname post_pred
#' @export

mtx_post_pred.brmsfit <-
	function(model, model_name, newdata = NULL, thin = 1, ...){
		n_iter <- brms::ndraws(model)
		n_draws <- round(n_iter/thin, 0)
		#draws <- sort(sample.int(n_iter, n_draws, replace = FALSE))
		brms::posterior_predict(model, newdata = newdata, ndraws = n_draws, summary = FALSE)

	}


#' @rdname post_pred
#' @export


mtx_post_pred.stanreg <-
	function(model, model_name, newdata = NULL, thin = 1, ...){
		n_iter <- sum(model$stanfit@sim$n_save)
		n_draws <- round(n_iter/thin, 0)
		rstanarm::posterior_predict(model, newdata = newdata, draws = n_draws)
	}

