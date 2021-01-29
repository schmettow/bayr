#' tibble class for observations
#'
#' A class for tibbles with observations providing a compact output
#'
#' @usage as_tbl_obs(x)
#' @param x a data frame
#' @return tbl_obs
#'
#' The constructor makes any data.frame or tibble class tbl_obs.
#' The provided print method randomly selects a small number of
#' observations and puts it into a knitr::kable, which looks nice
#' on both, console and knitr output.
#'
#' @author Martin Schmettow
#' @importFrom knitr knit_print
#' @export



as_tbl_obs <- function(x, ...) UseMethod("as_tbl_obs", x)

#' @rdname as_tbl_obs
#' @export


as_tbl_obs.tbl_df <- function(x) {
	x <- mutate(x, Obs = row_number()) %>%
		mascutils::go_first(Obs)
	class(x) <- c("tbl_obs", class(x))
	x
}

#' @rdname as_tbl_obs
#' @export


as_tbl_obs.data.frame <- function(x) {
	x <-
		tibble::as_tibble(x) %>%
		as_tbl_obs.tbl_df()
	class(x) <- c("tbl_obs", class(x))
	x
}

