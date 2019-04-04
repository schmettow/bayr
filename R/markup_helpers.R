library(tidyverse)

## dplyr is used with NSE, which gives "no visible binding for global variable errors"
utils::globalVariables(names = c("type", "parameter", "value", "new_name", "iter", "pattern","tbl_coef"))


#################### PRINT #######################

#' @rdname coef.tbl_post
#' @export

print.tbl_coef <- function(x, ...) {
	tab <- x
	if(nrow(tab) > 1)	{
		tab <- mascutils::discard_redundant(tab)
	} else if(tab$fixef[1] == "Intercept"){
		#		warning("Intercept model")
		tab <- select(tab, fixef, center, lower, upper)
	} else {
		tab <- mascutils::discard_all_na(x)
	}
	cap <-
		paste0("Estimates with ",
					 attr(x, "interval")*100,
					 "% credibility limits")
	out <-
		knitr::kable(tab, caption = cap)
	print(out)
	invisible(tab)
}


# print.tbl_coef_EATME <-
# 	function(x, digits = NULL, title = F, footnote = T,
# 					 kable = bayr:::by_knitr()){
# 		out <- mascutils::discard_all_na(x)
# 		if(nrow(out) > 1)	{
# 			out <- mascutils::discard_redundant(out)}
# 		else{
# 			out <- select(out, -model, -type)
# 		}
# 		types <-
# 			data_frame(type = c("fixef",
# 													"ranef",
# 													"grpef",
# 													"fitted"),
# 								 title_text = c("fixed effects",
# 								 							 "random effects",
# 								 							 "factor-level variation (sd)",
# 								 							 "fitted values (linear predictor)"))
# 		title_text <- ""
# 		if(title)
# 			title_text <-
# 			types %>%
# 			filter(type %in% attr(x, "type")) %>%
# 			select(title_text) %>%
# 			as.character()
#
# 		footnote_text <-
# 			paste0("\n*\nestimate with ",
# 						 attr(x, "interval")*100,
# 						 "% credibility limits")
#
# 		if(kable) { ## prepared for knitr table output, not yet working
# 			print(knitr::kable(out, caption = title_text))
# 		} else {
# 			if(title) cat(stringr::str_c(title_text, sep = "|", collapse = T), "\n***\n")
# 			print.data.frame(out, digits = 3, row.names = F)
# 			if(footnote) cat(footnote_text)
# 			cat("\n")
# 			invisible(out)
# 		}
# 		return(out)
# 	}


#' @rdname posterior
#' @export

print.tbl_post <-
	## TODO: add formula and corr
	function(tbl_post, ...){
		tbls <- prep_print_tbl_post(tbl_post)
		# frm <-
		# 	formula.tools:::as.character.formula(attr(tbl_post, "formula"))


		cat("** tbl_post: ", tbls$n_iter, " samples in ", tbls$n_chain, " chains\n\n")
		#		cat(frm, "\n\n")

		cat("** Effects: \n")
		print.data.frame(tbls$effects, row.names = F)

		cat("\n** Dispersion: \n")
		print.data.frame(tbls$disp, row.names = F)

		cat("\n** Shape: \n")
		print.data.frame(tbls$shape, row.names = F)

		cat("\n** Correlations: \n")
		print.data.frame(tbls$cor, row.names = F)

		cat("\n** User annotations: \n", tbls$user_annos)


		invisible(tbl_post)
	}

prep_print_tbl_post <-
	function(tbl_post){
		res <- list()

		res$n_iter <- length(unique(tbl_post$iter))
		res$n_chain <- length(unique(tbl_post$chain))

		res$user_annos <- setdiff(names(tbl_post),
													as.character(bayr:::AllCols))

		res$effects <-
			tbl_post %>%
			filter(type %in% c("fixef", "ranef", "grpef")) %>%
			distinct(model, parameter, type, fixef, nonlin, re_factor, re_entity) %>%
			mutate(parameter = ifelse(type == "ranef", "", parameter)) %>%
			group_by(model, parameter, type, nonlin, fixef, re_factor) %>%
			summarize(entities = n()) %>%
			ungroup() %>%
			mascutils::discard_all_na()

		res$cor <-
			tbl_post %>%
			filter(type == "cor") %>%
			distinct(parameter) %>%
			mascutils::discard_all_na()

		res$disp <-
			filter(tbl_post, type == "disp") %>%
			distinct(parameter) %>%
			mascutils::discard_all_na()

		res$shape <-
			filter(tbl_post, type == "shape") %>%
			distinct(parameter) %>%
			mascutils::discard_all_na()
		res
	}




#' @rdname post_pred
#' @export


print.tbl_post_pred <-
	function(x, ...){
		n_iter <- length(unique(x$iter))
		n_chain <- length(unique(x$chain))
		n_Obs <- length(unique(x$Obs))
		scales <- unique(x$scale)
		cap <- stringr::str_c("posterior predictions: ",
													n_iter, " samples in ", n_chain, " chains on ", n_Obs, " observations. (five shown below)")
		cat("**", cap, "\n\n")

		x %>%
			sample_n(5) %>%
			arrange(model, Obs, chain, iter) %>%
			print.data.frame()

		invisible(x)
	}


#' @rdname predict.tbl_post_pred
#' @export

print.tbl_predicted <-
	function(x, ...) {
		n_obs <- length(unique(x$Obs))
		cap <-
			stringr::str_c(n_obs, " predictions (scale: ", attr(x, "scale") ,") with ",
										 attr(x, "interval")*100, "% credibility limits (five shown below)")
		tab <-	x %>%
			sample_n(5) %>%
			arrange(Obs, model) %>%
			mascutils::discard_redundant() %>%
			mascutils::discard_all_na()

		cat("** ", cap, "\n")
		print.data.frame(tab)
		invisible(x)
	}




#################### KNIT_PRINT #######################


#' @rdname posterior
#' @export

knit_print.tbl_post <- function(x, ...) {
	tbls <- bayr:::prep_print_tbl_post(x)
	res <- paste0("\n\n** tbl_post: ", tbls$n_iter,
								" samples in ", tbls$n_chain, " chains\n\n",
								collapse ="\n")

	if(nrow(tbls$effects)) res <- c(res, knitr::kable(tbls$effects, format = "markdown",
																														cap = "Coefficients"), "\n")
	if(nrow(tbls$disp)) res <- c(res, knitr::kable(tbls$disp, format = "markdown",
																												 cap = "Dispersion"), "\n")
	if(nrow(tbls$shape)) res <- c(res, knitr::kable(tbls$shape, format = "markdown",
																													cap = "Shape"), "\n")
	if(nrow(tbls$cor)) res <- c(res, knitr::kable(tbls$cor, format = "markdown",
																												cap = "Correlations"), "\n")

	out <- paste0(res, collapse = "\n")

	knitr::asis_output(out)

	}



#' @rdname coef.tbl_post
#' @export



knit_print.tbl_coef <- function (x, ...)
{
	cap =
		paste0("Estimates with ",
					 attr(x, "interval")*100,
					 "% credibility limits")
	tab = mascutils::discard_redundant(x) %>%
		mascutils::discard_all_na()

	res = paste0(c("", "", knitr::kable(tab, caption = cap, format = "markdown", ...), "\n\n"),
							collapse = "\n")
	knitr::asis_output(res)
}


#' @rdname post_pred
#' @export

knit_print.tbl_post_pred <-
	function(x, ...){
		n_iter <- nrow(distinct(x, iter, chain))
		n_chain <- length(unique(x$chain))
		n_Obs <- length(unique(x$Obs))
		scales <- unique(x$scale)
		cap <- stringr::str_c("posterior predictions: ",
													n_iter, " samples in ", n_chain, " chains on ",
													n_Obs, " observations. (five shown below)")
		tab <- x %>%
			sample_n(5) %>%
			arrange(model, Obs, chain, iter) %>%
			mascutils::discard_redundant() %>%
			mascutils::discard_all_na()

		res = paste0(c("", "", knitr::kable(tab, caption = cap, format = "markdown", ...), "\n\n"),
								collapse = "\n")
		knitr::asis_output(res)
	}


#' @rdname predict.tbl_post_pred
#' @export

knit_print.tbl_predicted <-
	function(x, ...) {
	n_obs <- nrow(x)
	cap <- paste0(n_obs, " predictions (scale: ", attr(x, "scale") ,") with ",
									 attr(x, "interval")*100, "% credibility limits (five shown below)", collapse = "")
	tab <-	x %>%
		sample_n(5) %>%
		arrange(Obs, model) %>%
		mascutils::discard_redundant() %>%
		mascutils::discard_all_na()

	res = paste0(c("", "", knitr::kable(tab, caption = cap, format = "markdown", ...)), collapse = "\n")
	knitr::asis_output(res)
}

# knit_print.tbl_predicted






################ FORMULAS ###############################


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



