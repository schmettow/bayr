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
#' @param x regression model or tbl_postpred
#' @param model model
#' @param scale linpred or resp
#' @param center function for computing the center estimate (median)
#' @param interval credibility interval: .95
#' @param ... passing parameters to tbl_postpred (e.g. thin)
#' @return coefficient table with center and interval estimates per obs
#'
#' The standard center function is the posterior median
#'
#' @author Martin Schmettow
#' @import dplyr
#' @importFrom stats median quantile predict
#' @export

predict.tbl_post_pred <-
	function(x,
					 model = unique(x$model),
					 scale = c("resp"),
					 center =  median,
					 interval = .95, ...) {
		lower <- (1-interval)/2
		upper <- 1-((1-interval)/2)

		tbl_predicted <-
			x %>%
			group_by(model, Obs) %>%
			summarize(center = center(value),
								lower = quantile(value, lower),
								upper = quantile(value, upper)) %>%
			ungroup() %>%
			arrange(Obs, model)

		class(tbl_predicted) <- append("tbl_predicted",
																	 class(tbl_predicted))
		attr(tbl_predicted, "center") <- bquote(center)
		attr(tbl_predicted, "interval") <- interval
		attr(tbl_predicted, "lower") <- lower
		attr(tbl_predicted, "upper") <- upper
		attr(tbl_predicted, "scale") <- scale
		return(tbl_predicted)
	}



#' @rdname predict.tbl_post_pred
#' @export

predict <-
	function(x, ...) UseMethod("predict", x)


#' @rdname predict.tbl_post_pred
#' @export

predict.brmsfit <-
	function(x,
					 model = unique(x$model),
					 scale = c("resp"),
					 center =  median,
					 interval = .95, ...)
		tbl_post_pred(x, ...) %>% predict(model, scale, center, interval)



#' @rdname predict.tbl_post_pred
#' @export

predict.stanreg <-
	function(x,
					 model = unique(x$model),
					 scale = c("resp"),
					 center =  median,
					 interval = .95, ...)
		tbl_post_pred(x, ...) %>% predict(model, scale, center, interval)


#' @rdname predict.tbl_post_pred
#' @export

knit_print.tbl_predicted <- function(x, ...) {
	n_obs <- length(unique(x$Obs))
	cap <-
		stringr::str_c(n_obs, " predictions (scale: ", attr(x, "scale") ,") with ",
					 attr(x, "interval")*100, "% credibility limits (five shown below)")
	out <-	x %>%
		sample_n(5) %>%
		arrange(Obs, model) %>%
		mascutils::discard_redundant() %>%
		knitr::kable(caption = cap)

	print(out)
	invisible(x)
}


#' @rdname predict.tbl_post_pred
#' @export


print.tbl_predicted <- function(x, ...) {
	n_obs <- length(unique(x$Obs))
	cap <-
		stringr::str_c(n_obs, " predictions (scale: ", attr(x, "scale") ,") with ",
									 attr(x, "interval")*100, "% credibility limits (five shown below)")
	tab <-	x %>%
		arrange(Obs, model) %>%
		mascutils::discard_redundant() %>%
		sample_n(5)

	cat("** ", cap, "\n")
	print.data.frame(tab)
	invisible(x)
}




# predicted.MCMCglmm <-
# 	function(x, center =  median, ...)
# 		tbl_post_pred(x) %>% coef(center =  estimate, ...)


# predicted.stanfit <-
# 	function(x, center =  median, ...)
# 		tbl_post_pred(x) %>% coef(center =  estimate, ...)

