library(modeest)
library(dplyr)
library(tidyr)
library(stringr)

#library(dplyr)
#library(tidyr)
#library(stringr)

## dplyr is used with NSE, which gives "no visible binding for global variable errors"
utils::globalVariables(names = c("type", "parameter", "value", "new_name", "iter", "pattern"))

Cols_pp = list("model", "chain", "iter", "Obs", "scale", "value") ## columns for a postpred object to be valid



#' posterior predictive extraction
#'
#' MCMC predicted values are  extracted from a Bayesian (regression) object
#' and returned as a tbl_post_pred object
#'
#' @usage post_pred(model, scale = "obs", model_name, thin)
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
#' chain  iter   parameter   value  type  order
#' (fctr) (int)  (chr)       (dbl)  (chr) (int)
#'
#'
#' @author Martin Schmettow
#' @import dplyr
#' @export


post_pred <-
	function(model,
					 scale = "obs",
					 model_name = NA,
					 thin = 1, ...){
		if(is.na(model_name)) model_name <- deparse(substitute(model))

		post <- tbl_post_pred(model, ...)

		if(! all(as.character(bayr:::Cols_pp) %in% names(post)))
			stop("not a valid tbl_post_pred")

		out <- post %>%
			select_(.dots = bayr:::Cols_pp) %>%
			arrange_(.dots = bayr:::Cols_pp) %>%
			mutate(model = model_name)

		class(out) <- append("tbl_post_pred", class(out))
		return(out)
	}


#' @rdname post_pred
#' @export

tbl_post_pred <-
	function (model, ...) {
		UseMethod("tbl_post_pred", model)
	}

#' @rdname post_pred
#' @export

print.tbl_post_pred <-
	function(tbl_post_pred, kable = by_knitr(), ...){
		n_iter <- length(unique(tbl_post_pred$iter))
		n_chain <- length(unique(tbl_post_pred$chain))
		n_Obs <- length(unique(tbl_post_pred$Obs))
		scales <- unique(tbl_post_pred$scale)

		# if(kable){ ## prepared for knitr table output, not yet working
		# 	cat("tbl_post_pred: ", n_iter, " samples in ", n_chain, " chains\n")
		# 	cat("Effects: \n")
		# 	print(knitr::kable(effects))
		# 	cat("\nDispersion: \n")
		# 	print(knitr::kable(disp))
		# } else {
		#
		cat("** tbl_post_pred : ",
				n_iter, " samples in ", n_chain, " chains\n** Observations :  ", n_Obs, "\n\n")

		tbl_post_pred %>%
			sample_n(5) %>% print.data.frame()

		invisible(tbl_post_pred)
	}


#' @rdname post_pred
#' @export

tbl_post_pred.data.frame <-
	## IDEA: write methods for
	## - identifying user_annos (all user annos)
	## - registering user annos (explicit user annos)
	## - keep attribute user_annos (keep user annos)
	function(df, ...) {
		if(! all(as.character(bayr:::AllCols) %in% names(df)))
			stop("not a valid tbl_post_pred, some columns missing")
		out = df
		class(out) = append("tbl_post_pred", class(out))
		out
	}



tbl_post.generic <-
	function(sample_matrix){
		out <-
			sample_matrix %>%
			as.data.frame() %>%
			as_data_frame() %>%
			mutate(model = NA,
						 chain = 1, ## fixme
						 scale = "resp", ## fixme
						 iter = 1:n()) %>%
			tidyr::gather(key = "Obs",
										value = "value", -model, -iter, -chain, -scale) %>%
			mutate(Obs = as.integer(stringr::str_replace(Obs, "^V", ""))) %>%  ## making observations integer (Vx)
			select_(.dots = bayr:::Cols_pp) %>%
			arrange_(.dots = bayr:::Cols_pp)

		class(out) <-
			append("tbl_post_pred", class(out))
		out
	}

#' @rdname post_pred
#' @export

tbl_post_pred.brmsfit <-
	function(model, newdata = NULL, ...){
		brms:::predict.brmsfit(model, newdata = newdata, summary = F) %>%
			bayr:::tbl_post.generic()
	}


#' @rdname post_pred
#' @export


tbl_post_pred.stanreg <-
	function(model, newdata = NULL, ...){
		rstanarm::posterior_predict(model, newdata = newdata) %>%
			bayr:::tbl_post.generic()
	}


