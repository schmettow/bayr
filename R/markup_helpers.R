library(modeest)
library(dplyr)
library(tidyr)
library(stringr)

## dplyr is used with NSE, which gives "no visible binding for global variable errors"
utils::globalVariables(names = c("type", "parameter", "value", "new_name", "iter", "pattern","tbl_coef"))


################ COEF ###############################


#' mrk_coef
#'
#' inline markup reporting of estimates with credibility limits (formula)
#'
#'
#' @param tbl_coef coefficient table
#' @param coef coefficient (row number or name)
#' @param center add center estimate (TRUE)
#' @param interval add interval (TRUE)
#' @param prefix add prefix term with coef name (not implemented)
#' @param rounding digits (2)
#' @return markdown string
#'
#'
#' @author Martin Schmettow
#' @export

mrk_coef =
	function(tbl_coef,
					 coef,
					 center = T,
					 interval = T,
					 prefix = F,
					 round = 2) {
		if(is.character(coef)){
			coefpos = which(tbl_coef$parameter == coef)
			if(length(coefpos) == 0) stop("coefficient ", coef, " not found in table")
		} else {
			coefpos = coef
		}

		out = as.character()

		if(prefix) stop("prefix not yet implemented")

		if(center) out =
			str_c(out, round(tbl_coef[[coefpos, 2]], round))

		if(interval) out =
			str_c(out, " [", round(tbl_coef[[coefpos, 3]], round), ", ",
						round(tbl_coef[[coefpos, 4]], round), "]_{CI",attr(tbl_coef, "interval") * 100 ,"}")
		out
	}

