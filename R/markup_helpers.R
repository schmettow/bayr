library(modeest)
library(dplyr)
library(tidyr)
library(stringr)

## dplyr is used with NSE, which gives "no visible binding for global variable errors"
utils::globalVariables(names = c("type", "parameter", "value", "new_name", "iter", "pattern","tbl_coef"))


################ COEF ###############################


#' coefficient formula
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
#' @param neg output negative estimate
#' @param mean_fnc mean function (identity)
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
					 interval = F,
					 prefix = F,
					 round = 2,
					 mean_fnc = identity,
					 neg = F) {

		filter_crit = list(...)

		if(!c("tbl_coef") %in% class(tbl_coef)) stop("coefficient table required as input, e.g. coef(posterior)")
		if(!is.null(row)){
			if(!is.numeric(row)) stop("row must be integer")
			if(row < 1 || row > nrow(tbl_coef)) stop("row must be between ", 1, " and ", nrow(tbl_coef))
			if(length(filter_crit) > 0) warning("selection by row overides selection by filter")
			tbl_coef = slice(tbl_coef, row)
		} else {
			tbl_coef = filter_(tbl_coef, .dots = filter_crit)
		}

		if(nrow(tbl_coef) == 0) stop("md_coef: parameter does not exist")
		if(nrow(tbl_coef) > 1)  stop("md_coef: parameter is not unique, ", print(tbl_coef))

		out = as.character()

		if(prefix) warning("prefix not yet implemented")

		if(neg) {
			tbl_coef$center <- -tbl_coef$center
			tbl_coef$lower <- -tbl_coef$lower
			tbl_coef$upper <- -tbl_coef$upper
		}

		if(center) out =
			stringr::str_c(out, round(mean_fnc(tbl_coef[[1, "center"]]), round))

		if(interval) out =
			stringr::str_c(out,
										 " [", round(mean_fnc(tbl_coef[[1, "lower"]]), round), ", ",
										 round(mean_fnc(tbl_coef[[1, "upper"]]), round),
										 "]_{CI",attr(tbl_coef, "interval") * 100 ,"}")
		out
	}

#' @rdname md_coef
#' @export

frm_coef =
	function(tbl_coef, ...,
					 row = NULL,
					 center = T,
					 interval = T, # <--
					 prefix = F,
					 round = 2,
					 neg = F) {
		out = bayr::md_coef(tbl_coef, ...,
												row = row,
												center = center,
												interval = interval,
												prefix = prefix,
												round = round,
												neg = neg)
		stringr::str_c("$", out, "$")
	}

