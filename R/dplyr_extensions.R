library(tidyverse)
## dplyr is used with NSE, which gives "no visible binding for global variable errors"
utils::globalVariables(names = c("type", "parameter", "value", "new_name", "iter", "pattern","tbl_coef"))


################ DPLYR overloaded ###############################


#' computing new variables
#'
#' computing new variables, like dplyr::mutate does it.
#'
#' @param tbl_post posterior table
#'
#' NOT IMPLEMENTED
#' The code is copied from dyplr::mutate_ and adapted to work with
#' tbl_post. Internally, tbl_post is converted to wide format, then back to long
#'
#'
#' @author Martin Schmettow
#' @import dplyr
#' @import tidyr

# mutate.tbl_post <-
# 	function (.data, ..., .dots)
# 	{
# 		print(class(.data))
# 		dots <- lazyeval::all_dots(.dots, ..., all_named = TRUE)
# 		dplyr:::mutate_impl(.data, dots)
# 		out <- dplyr:::distinct_impl(.data, dots)
# 		class(out)[2:length(class(out))]
# 		out
# 	}
#
# D <- data_frame(x = 1:6, y = 3:8, z = 5:10)
# mymutate(D,s = x + y)



#' @rdname mutate.tbl_post

distinct.tbl_post <-
	function (.data, ..., .dots)
	{
		dots <- lazyeval::all_dots(.dots, ..., all_named = TRUE)
		out <- dplyr:::distinct_impl(.data, dots)
		class(out)[2:length(class(out))]
		out
	}


########################### FUTURE ##########################



# TODO
# overload mutate
