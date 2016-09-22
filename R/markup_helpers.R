library(modeest)
library(dplyr)
library(tidyr)
library(stringr)

## dplyr is used with NSE, which gives "no visible binding for global variable errors"
utils::globalVariables(names = c("type", "parameter", "value", "new_name", "iter", "pattern","tbl_coef"))


################ COEF ###############################


#' md_coef
#'
#' inline markup reporting of estimates with credibility limits
#'
#'
#' @param tbl_coef coefficient table
#' @param ... filter formulas, using annotations
#' @param row row number of parameter
#' @param center add center estimate (TRUE)
#' @param interval add interval (TRUE)
#' @param prefix add prefix term with coef name (not implemented)
#' @param rounding digits (2)
#' @return markdown string
#'
#' The parameter can alternatively selected by row number, or using a
#' dplyr-like selection formula. When the selection criteria are ambiguous
#' an error is thrown
#'
#' @author Martin Schmettow
#' @export

md_coef =
	function(tbl_coef, ...,
					 row = NULL,
					 center = T,
					 interval = T,
					 prefix = F,
					 round = 2) {

		filter_crit = list(...)
		if(!is.null(row)){
			if(!is.numeric(row)) stop("row must be integer")
			if(row < 1 || row > nrow(tbl_coef)) stop("row must be between ", 1, " and ", nrow(tbl_coef))
			if(length(filter_crit) > 0) warning("selection by row overides selection by filter")
			tbl_coef = slice(tbl_coef, row)
		} else {
			tbl_coef = filter_(tbl_coef, .dots = filter_crit)
		}

		if(nrow(tbl_coef) == 0) stop("mrk_coef: parameter does not exist")
		if(nrow(tbl_coef) > 1)  stop("mrk_coef: parameter is not unique, ", print(tbl_coef))

		out = as.character()

		if(prefix) warning("prefix not yet implemented")

		if(center) out =
			stringr::str_c(out, round(tbl_coef[[1, "center"]], round))

		if(interval) out =
			stringr::str_c(out, " [", round(tbl_coef[[1, "lower"]], round), ", ",
						round(tbl_coef[[1, "upper"]], round), "]_{CI",attr(tbl_coef, "interval") * 100 ,"}")
		out
	}

