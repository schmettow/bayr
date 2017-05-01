library(modeest)
library(dplyr)
library(tidyr)
library(stringr)

## dplyr is used with NSE, which gives "no visible binding for global variable errors"
utils::globalVariables(names = c("type", "parameter", "value", "new_name", "iter", "pattern","tbl_coef"))


################ PREDICTED ###############################


#' Posterior predictions
#'
#' summary table of predicted values from predictive posterior
#'
#' @param object tbl_post_pred object holding the predictive posterior in long format
#' @param model model
#' @param scale linpred or resp
#' @param center function for computing the center estimate (median)
#' @param interval credibility interval: .95
#' @param ... ignored
#' @return coefficient table with center and interval estimates per obs
#'
#' The standard center function is the posterior median
#'
#' @author Martin Schmettow
#' @import dplyr
#' @importFrom stats median quantile predict
#' @export

predict.tbl_post_pred <-
	function(object,
					 model = unique(object$model),
					 scale = c("resp"),
					 center =  median,
					 interval = .95, ...) {
		lower <- (1-interval)/2
		upper <- 1-((1-interval)/2)

		tbl_predicted <-
			object %>%
			group_by(Obs) %>%
			summarize(center = center(value),
								lower = quantile(value, lower),
								upper = quantile(value, upper)) %>%
			ungroup() %>%
			arrange(Obs)

		class(tbl_predicted) <- append("tbl_predicted",
																	 class(tbl_predicted))
		attr(tbl_predicted, "center") <- bquote(center)
		attr(tbl_predicted, "interval") <- interval
		attr(tbl_predicted, "lower") <- lower
		attr(tbl_predicted, "upper") <- upper
		return(tbl_predicted)
	}

#' @rdname predict.tbl_post_pred
#' @export
#
predict.brmsfit <-
	function(object, center =  median, ...)
		tbl_post_pred(object) %>% predict(center =  center, ...)



#' @rdname predict.tbl_post_pred
#' @export
#
predict.stanreg <-
	function(object, center =  median, ...)
		tbl_post_pred(object) %>% predict(center =  center, ...)


# predicted.MCMCglmm <-
# 	function(object, center =  median, ...)
# 		tbl_post_pred(object) %>% coef(center =  estimate, ...)


# predicted.stanfit <-
# 	function(object, center =  median, ...)
# 		tbl_post_pred(object) %>% coef(center =  estimate, ...)

