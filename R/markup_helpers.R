#library(tidyverse)

## dplyr is used with NSE, which gives "no visible binding for global variable errors"


#################### TBL_POST #######################


prep_print_tbl_post <-
	function(tbl_post){
		res <- list()
		res$n_model <- n_distinct(tbl_post$model)
		res$n_iter <- n_distinct(tbl_post$iter)
		res$n_chain <- n_distinct(tbl_post$chain)
		res$n_param <- nrow(distinct(tbl_post, model, parameter))
		res$user_annos <- setdiff(names(tbl_post),
															as.character(bayr:::AllCols))


		res$model <-
			tbl_post %>%
			group_by(model) %>%
			summarize(n_param = n_distinct(parameter),
								n_iter = n_distinct(iter))


		res$effects <-
			tbl_post %>%
			as_tibble() %>%
			filter(type %in% c("fixef", "ranef")) %>%
			distinct(model, parameter, type, fixef, nonlin, re_factor, re_entity) %>%
			mutate(parameter = ifelse(type == "ranef", "", parameter)) %>%
			group_by(model, parameter, type, nonlin, fixef, re_factor) %>%
			summarize(count = n()) %>%
			ungroup() %>%
			discard_all_na()

		res$disp <-
			tbl_post %>%
			as_tibble() %>%
			filter(type %in% c( "disp", "grpef")) %>%
			distinct(model, parameter, fixef, type)

		res$shape <-
			tbl_post %>%
			as_tibble() %>%
			filter(type == "shape") %>%
			distinct(model, parameter, type) ## also creates empty columns

		## Correlations are just counted
		res$cor <-
			tbl_post %>%
			as_tibble() %>%
			filter(type == "cor") %>%
			distinct(model, parameter, type) %>%
			group_by(model, type) %>%
			summarize(count = n()) %>%
			mutate(parameter = "")

		res$comb <-
			full_join(res$effects, res$disp) %>%
			full_join(res$shape) %>%
			full_join(res$cor)

		res
	}




#' @rdname posterior
#' @export

print.tbl_post <-
	function(x, ...){

		tbls <- prep_print_tbl_post(x)
		# frm <-
		# 	formula.tools:::as.character.formula(attr(tbl_post, "formula"))


		cat("MCMC posterior: ", tbls$n_iter, " samples of ",
				tbls$n_param ," parameters  in ", tbls$n_model, " model(s)\n\n")
		#		cat(frm, "\n\n")

		print.data.frame(tbls$comb, row.names = F)

		invisible(x)
	}


#' @rdname posterior
#' @export

knit_print.tbl_post <- function(x, ...) {
	tbls <- prep_print_tbl_post(x)
	tab <- tbls$comb
	cap <- paste0("MCMC posterior with ", tbls$n_iter, " samples of ",
								tbls$n_param ," parameters in ", tbls$n_model, " model(s)")
	kab <- knitr::kable(tab, caption = cap, format = "markdown", ...)
	out <- paste0(c("", "", kab, "\n\n"), collapse = "\n")

	knitr::asis_output(out)
}




knit_print.tbl_post_old <- function(x, ...) {
	tbls <- bayr:::prep_print_tbl_post(x)
	# res <- paste0("\n\n** tbl_post: ", tbls$n_iter,
	# 							" samples in ", tbls$n_chain, " chains\n\n",
	# 							collapse ="\n")

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






#################### POST_PRED #######################




# print.tbl_coef_EATME <-
# 	function(x, digits = NULL, title = F, footnote = T,
# 					 kable = bayr:::by_knitr()){
# 		out <- discard_all_na(x)
# 		if(nrow(out) > 1)	{
# 			out <- discard_redundant(out)}
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


#' @rdname post_pred
#' @export



print.tbl_post_pred <-
	function(x, ...){
		n_iter <- n_distinct(x$iter)
		n_chain <- n_distinct(x$chain)
		n_Obs <- n_distinct(x$Obs)
		scales <- unique(x$scale)
		cap <- stringr::str_c("posterior predictions: ",
													n_iter, " samples in ", n_chain, " chains on ", n_Obs, " observations. (five shown below)")
		cat("**", cap, "\n\n")

		x %>%
			sample_n(min(n_Obs, 5)) %>%
			arrange(model, Obs, chain, iter) %>%
			print.data.frame()

		invisible(x)
	}

#' @rdname post_pred
#' @export

print.tbl_predicted <-
	function(x, ...) {
		n_Obs <- n_distinct(x$Obs)
		cap <-
			stringr::str_c(n_Obs, " predictions (scale: ", attr(x, "scale") ,") with ",
										 attr(x, "interval")*100, "% credibility limits (five shown below)")
		tab <-	x %>%
			sample_n(min(n_Obs, 5)) %>%
			arrange(Obs, model) %>%
			discard_redundant() %>%
			discard_all_na()

		cat("** ", cap, "\n")
		print.data.frame(tab)
		invisible(x)
	}



#' @rdname post_pred
#' @export

knit_print.tbl_post_pred <-
	function(x, ...){
		n_iter <- nrow(distinct(x, iter, chain))
		n_chain <- n_distinct(x$chain)
		n_Obs <- n_distinct(x$Obs)
		scales <- unique(x$scale)
		cap <- stringr::str_c("posterior predictions: ",
													n_iter, " samples in ", n_chain, " chains on ",
													n_Obs, " observations. (five shown below)")
		tab <- x %>%
			sample_n(min(n_Obs, 5)) %>%
			arrange(model, Obs, chain, iter) %>%
			discard_redundant() %>%
			discard_all_na()

		res = paste0(c("", "", knitr::kable(tab, caption = cap, format = "markdown", ...), "\n\n"),
								 collapse = "\n")
		knitr::asis_output(res)
	}


#' @rdname post_pred
#' @export

knit_print.tbl_predicted <-
	function(x, ...) {
		n <- min(8, nrow(x))
		n_Obs <- nrow(x)
		cap <- paste0(n_Obs, " predictions (scale: ", attr(x, "scale") ,") with ",
									attr(x, "interval")*100, "% credibility limits (8 shown)", collapse = "")
		tab <-	x %>%
			sample_n(min(n_Obs, 8)) %>%
			arrange(Obs, model) %>%
			discard_redundant() %>%
			discard_all_na()

		res = paste0(c("", "", knitr::kable(tab, caption = cap, format = "markdown", ...)), collapse = "\n")
		knitr::asis_output(res)
	}

# knit_print.tbl_predicted




#################### COEF ######################

#' @rdname coef.tbl_post
#' @export

print.tbl_coef <- function(x, ...) {
	tab <- x
	if(nrow(tab) > 1)	{
		tab <- discard_redundant(tab)
	} else if(tab$fixef[1] == "Intercept"){
		#		warning("Intercept model")
		tab <- select(tab, fixef, center, lower, upper)
	} else {
		tab <- discard_all_na(x)
	}
	cap <-
		paste0("Coefficient estimates with ",
					 attr(x, "interval")*100,
					 "% credibility limits")
	out <-
		knitr::kable(tab, caption = cap)
	print(out)
	invisible(tab)
}




#' @rdname coef.tbl_post
#' @export



knit_print.tbl_coef <- function (x, ...)
{
	cap =
		paste0("Coefficient estimates with ",
					 attr(x, "interval")*100,
					 "% credibility limits")
	tab = discard_redundant(x) %>%
		discard_all_na()

	res = paste0(c("", "", knitr::kable(tab, caption = cap, format = "markdown", ...), "\n\n"),
							 collapse = "\n")
	knitr::asis_output(res)
}


#################### FIXEF_ML ######################

#' @rdname fixef_ml
#' @export

print.tbl_fixef_ml <- function(x, ...) {
	tab <- x
	cap <-	paste0("Population-level coefficients with random effects standard deviations")
	out <-
		knitr::kable(tab, caption = cap)
	print(out)
	invisible(tab)
}




#' @rdname fixef_ml
#' @export



knit_print.tbl_fixef_ml <- function (x, ...)
{
	cap <-	paste0("Population-level coefficients with random effects standard deviations")
	tab <- x %>%	discard_all_na()
	if(nrow(tab > 1)) tab <- tab %>% discard_redundant()
	res <- paste0(c("", "", knitr::kable(tab, caption = cap, format = "markdown", ...), "\n\n"),
							 collapse = "\n")
	knitr::asis_output(res)
}




#################### CLU #######################

pre_print_tbl_clu <- function(tbl_clu){
	cols <- c()
	types <- unique(tbl_clu$type)
	if(n_distinct(tbl_clu$model) > 1) cols <- union(cols, "model")
	cols <- union(cols, "parameter")
	if("fixef" %in% types)          cols <- union(cols, "fixef")
	if("grpef" %in% types)          cols <- union(cols, c("fixef", "re_factor"))
	if("ranef" %in% types)          cols <- union(cols, c("fixef", "re_factor", "re_entity"))
	if("disp"  %in% types)          cols <- union(cols, c("parameter"))
	cols <- c(cols, "center", "lower", "upper")
	out <- tbl_clu %>%
		select(cols)
	out <- out %>%	discard_all_na()
	if(nrow(out > 1)) out <- out %>% discard_redundant()
	out
}


#' @rdname clu
#' @export


# print.tbl_clu.old <- function(x, ...) {
# 	tab <- x
# 	if(nrow(tab) > 1)	{
# 		tab <- discard_redundant(tab)
# 	} else if(!is.null(tab$fixef) & tab$fixef[1] == "Intercept"){
# 		#		warning("Intercept model")
# 		tab <- select(tab, fixef, center, lower, upper)
# 	} else {
# 		tab <- discard_all_na(x)
# 	}
# 	cap <-
# 		paste0("Estimates with ",
# 					 attr(x, "interval")*100,
# 					 "% credibility limits")
# 	out <-
# 		knitr::kable(tab, caption = cap)
# 	print(out)
# 	invisible(tab)
# }


print.tbl_clu <- function(x, ...) {
	tab <- pre_print_tbl_clu(x)
	if(dim(tab)[2] == 0) warning("Empty CLU table. Do you need to adjust filters applied to posterior object?")
	cap <-
		paste0("Parameter estimates with ",
					 attr(x, "interval")*100,
					 "% credibility limits")
	out <-
		knitr::kable(tab, caption = cap)
	print(out)
	invisible(x)
}



#' @rdname clu
#' @export


knit_print.tbl_clu <- function (x, ...)
{
	tab <- pre_print_tbl_clu(x)
	cap =
		paste0("Parameter estimates with ",
					 attr(x, "interval")*100,
					 "% credibility limits")
	res = paste0(c("", "", knitr::kable(tab, caption = cap, format = "markdown", ...), "\n\n"),
							 collapse = "\n")
	knitr::asis_output(res)
}



################ OBS ############################

#' @rdname as_tbl_obs
#' @export

print.tbl_obs <- function(x, ...) {
	n <- min(8, nrow(x))
	tab <- dplyr::sample_n(x, n)
	if("Obs" %in% colnames(tab)) tab <- dplyr::arrange(tab, Obs)
	if("Part" %in% colnames(tab)) tab <- dplyr::arrange(tab, Part)
	cap <- stringr::str_c("Data set",": showing ", n, " of ", nrow(x), " observations")
	print(cap)
	base:::print.data.frame(tab)
	invisible(x)
}

#' @rdname as_tbl_obs
#' @export

knit_print.tbl_obs <- function(x, ...) {
	#data_set <- deparse(substitute(x))
	n <- min(8, nrow(x))
	cap <- paste0("Data set with ", ncol(x)," variables, showing ", n, " of ", nrow(x), " observations.")
	tab <- dplyr::sample_n(x, n)
	if("Obs" %in% colnames(tab)) tab <- dplyr::arrange(tab, Obs)
	if("Part" %in% colnames(tab)) tab <- dplyr::arrange(tab, Part, Obs)

	kab <- knitr::kable(tab, caption = cap, format =
												"markdown", ...)
	out <- paste0(c("", "", kab, "\n\n"), collapse = "\n")
	knitr::asis_output(out)
}



knit_print.tbl_obs_old <- function(x, ...) {
	#data_set <- deparse(substitute(x))
	n <- min(8, nrow(x))
	tab <- dplyr::sample_n(x, n)
	if("Obs" %in% colnames(tab)) tab <- dplyr::arrange(tab, Obs)
	if("Part" %in% colnames(tab)) tab <- dplyr::arrange(tab, Part)
	cap <- stringr::str_c("Data set with ", ncol(x)," variables. Showing ", n, " of ", nrow(x), " observations.")

	out <- knitr::kable(tab, caption = cap, format = "markdown")
	knitr::asis_output(out)
	#
	# res = paste(c("", "", knitr::kable(tab, format = "markdown", ...)),
	#             collapse = "\n")
	# knitr::asis_output(res)
	invisible(tab)
}


################# IC ########################

#' @rdname IC
#' @export

print.tbl_IC <- function(x, ...) {
	cap <- stringr::str_c("Estimated Information Criterion")
	base:::print.data.frame(x)
	invisible(x)
}

#' @rdname IC
#' @export

knit_print.tbl_IC <- function(x, ...) {
	cap <- paste0("Estimated Information Criterion")
	kab <- knitr::kable(x, caption = cap, format = "markdown", ...)
	out <- paste0(c("", "", kab, "\n\n"), collapse = "\n")
	knitr::asis_output(out)
}

#' @rdname compare_IC
#' @export

print.tbl_IC_comp <- function(x, ...) {
	cap <- stringr::str_c("Model ranking by predictive accuracy")
	base:::print.data.frame(x)
	invisible(x)
}

#' @rdname compare_IC
#' @export

knit_print.tbl_IC_comp <- function(x, ...) {
	cap <- paste0("Model ranking by predictive accuracy")
	kab <- knitr::kable(x, caption = cap, format =
												"markdown", ...)
	out <- paste0(c("", "", kab, "\n\n"), collapse = "\n")
	knitr::asis_output(out)
}



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

md_coef <- function(tbl_coef, ...,
										row = NULL,
										center = T,
										interval = F,
										prefix = F,
										round = 2,
										mean_fnc = identity,
										neg = F) {


	if(0 == sum(stringr::str_detect(class(tbl_coef),
																	"tbl_clu|tbl_coef")))
		stop("coefficient table required as input, e.g. coef(posterior)")
	if(!is.null(row)){
		if(!is.numeric(row)) stop("row must be integer")
		if(row < 1 || row > nrow(tbl_coef)) stop("row must be between ", 1, " and ", nrow(tbl_coef))
		#if(length(...) > 0) warning("selection by row overides selection by filter")
		tbl_coef = slice(tbl_coef, row)
	} else {
		tbl_coef = filter(tbl_coef, ...)
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
		out =  md_coef(tbl_coef, ...,
												row = row,
												center = center,
												interval = interval,
												prefix = prefix,
												round = round,
												neg = neg)
		stringr::str_c("$", out, "$")
	}



