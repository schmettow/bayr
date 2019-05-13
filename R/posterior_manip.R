library(tidyverse)

## dplyr is used with NSE, which gives "no visible binding for global variable errors"
utils::globalVariables(names = c("type", "parameter", "value",
																 "new_name", "iter", "pattern"))


DrawIDCols <- c("model", "chain", "iter")
ParameterIDCols = list("parameter", "type", "nonlin", "fixef", "re_factor", "re_entity")
AllCols = append(append(list("model", "chain", "iter", "order"), ParameterIDCols), "value")


#' total random effects scores
#'
#' sums over random effects and their respective population-level effects
#'
#' @usage re_scores(tbl_post, fixef, ranef)
#' @param tbl_post posterior (tbl_post)
#' @param re_factor random factor (char)
#' @return tbl_post
#'
#' Random coefficients are created as differences to the group mean.
#' In some situations one needs the total (or absolute) score.
#'
#' Example: It is typical for psychometric settings to think of item and
#' person scores. An average person has an IQ of 100,
#' but would have a zero in terms of random effects.
#'
#'
#' @author Martin Schmettow
#' @import dplyr
#' @export


re_scores <-
	function(tbl_post, re_factor){

		#tbl_post <- P_1 %>% filter(iter < 3)

		tbl_ranef <-
			tbl_post %>%
			as_tibble() %>%
			dplyr::filter(type == "ranef",
										re_factor == re_factor)
		#select(model, chain, iter, fixef, nonlin, re_factor, re_entity, value)

		tbl_fixef <-
			tbl_post %>%
			as_tibble() %>%
			dplyr::filter(type == "fixef") %>%
			rename(fe_value = value) %>%
			select(model, chain, iter, fixef, nonlin, fe_value)

		scores <-
			tbl_ranef %>%
			dplyr::left_join(tbl_fixef,
											 by = c("model", "chain", "iter", "fixef", "nonlin")) %>%
			dplyr::mutate(value = fe_value + value,
										type = "score") %>%
			dplyr::select(-fe_value) %>%
			bayr:::tbl_post.data.frame()

		return(scores)
	}



# re_scores <-
# 	function(tbl_post,
# 					 fixef,
# 					 re_factor,
# 					 type = "ranef"){
#
# 		# tbl_post <- P_1
# 		# fixef <- "^DesignNovel|Intercept"
# 		# re_factor <- "Part"
# 		# type = "score"
#
# 		fe <- fixef
# 		rf <- re_factor
# 		tp <- type
#
# 		tbl_ranef <-
# 			tbl_post %>%
# 			as_tibble() %>%
# 			dplyr::filter(stringr::str_detect(fixef, fe),
# 										re_factor == rf,
# 										type == "ranef")
#
# 		tbl_fixef <-
# 			tbl_post %>%
# 			as_tibble() %>%
# 			dplyr::filter(stringr::str_detect(fixef, fe),
# 						 type == "fixef") %>%
# 			mutate(fe_value = value) %>%
# 			select(model, chain, iter, fixef, fe_value)
#
# 		scores <-
# 			tbl_ranef %>%
# 			dplyr::left_join(tbl_fixef,
# 											 by = c("model", "chain", "iter", "fixef")) %>%
# 			dplyr::mutate(value = fe_value + value,
# 						 type = tp) %>%
# 			dplyr::select(-fe_value) %>%
# 			bayr:::tbl_post.data.frame()
# 		return(scores)
# 	}

# load("M_1.Rda")
# P_1 <- bayr::posterior(M_1_s)

# re_scores(P_1, "Intercept", "Part")

