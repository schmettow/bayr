
#' Tibble-returning expand.grid
#'
#' works like expand.grid, but returns a tibble. Named `expand_grid_df()` to
#' avoid masking `tidyr::expand_grid()`.
#'
#' @param ... factors
#' @return tibble
#'
#' @author Martin Schmettow
#' @export


expand_grid_df <-
	function(...)
		expand.grid(stringsAsFactors = FALSE, ...) %>% tibble::as_tibble()


#' Moving column to first position
#'
#' The referred column becomes the first, keeping the order of the others
#'
#' @param D data frame
#' @param ... columns to move, using dplyr::select tidyselect semantics
#' @return data frame
#'
#' go_first and go_arrange expect a column specification, similar to
#' dplyr::select
#'
#' @examples
#' D <- tibble::tibble(x = 1:3, y = 6:4, z = c(8,9,7))
#' go_first(D, y)
#' go_first(D, y:z)
#' go_arrange(D, y)
#'
#'
#'
#' @author Martin Schmettow
#' @export


go_first <-
	function(D, ...){
		class <- class(D)
		cols <- quos(...)
		df1 <- dplyr::select(D, !!!cols)
		out <- left_union(df1, D)
		class(out) <- class
		out
	}


#' @rdname go_first
#' @export

go_arrange <-
	function(D, ...){
		class <- class(D)
		cols <- quos(...)
		df1 <- dplyr::select(D, !!!cols)
		out <- left_union(df1, D) %>%
			dplyr::arrange(!!!rlang::syms(names(df1)))
		class(out) <- class
		out
	}



#' adding columns that don't exist yet
#'
#'
#' @param df1 master data frame, where columns are kept
#' @param df2 data frame where additional columns are taken from
#' @return data frame
#'
#'
#'
#' @author Martin Schmettow
#' @export

left_union <-
	function(df1, df2){
		class <- class(df1)
		if(nrow(df1) != nrow(df2)) stop("dataframes have different number of rows")
		added_cols <- setdiff(names(df2), names(df1))
		out <- dplyr::bind_cols(df1, df2[,added_cols])
		class(out) <- class
		out
	}


#' update columns conditionally
#'
#'
#' @param D data frame
#' @param by predicate function (like dplyr::filter)
#' @param ... expressions
#' @return data frame
#'
#' Applies mutations to the filtered group, only.
#'
#' @examples
#' D <- tibble::tibble(group = c(1, 1, 2, 2), value = c(4, 9, -4, -9))
#'
#' dplyr::mutate(D, value = dplyr::if_else(group == 1, sqrt(value), value))
#' ## Produces NaNs, because sqrt() is evaluated before selection
#'
#' update_by(D, group == 1, value = sqrt(value))
#' ## sqrt() is only evaluated
#'
#' @author Martin Schmettow
#' @export



update_by <-
	function(D, by, ...){
		flt <- enquo(by)
		args <- enquos(...)

		missing_cols <- dplyr::setdiff(names(args), names(D))
		if (length(missing_cols > 0)) stop(paste0("\nColumn does not exist: ",
																							missing_cols))

		D <- dplyr::mutate(D, .tmp_idx = dplyr::row_number())

		mut <- D %>%
			dplyr::filter(!!flt) %>%
			dplyr::mutate(!!!args)

		rmn <- D %>%
			dplyr::filter(!(!!flt))

		bind_rows(mut, rmn) %>%
			dplyr::arrange(.tmp_idx) %>%
			dplyr::select(-.tmp_idx)
	}





#' Adding z-transformed scores
#'
#' The (simple) range of columns are z-transformed and added to the data frame (foo_z)
#'
#' @param D data frame
#' @param ... dplyr::select range of variables
#' @return data frame
#'
#'
#' @author Martin Schmettow
#' @importFrom stats sd
#' @export


z_trans <- function(D,  ...){
	col_spec <- quos(...)
	df_z <- dplyr::select(D, !!!col_spec) %>%
		mutate(across(everything(), z))
	names(df_z) <- stringr::str_c("z", names(df_z), sep = "")
	bind_cols(D, df_z)
}

z <- function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)

#' @rdname z_trans
#' @export

z_score <- z_trans


#' centered rescaling
#'
#' symmetrically scales a variable up or down, maintaining the center of the scale
#'
#' @param x numerical vector
#' @param scale rescale factor
#' @return numerical vector
#'
#'
#' @author Martin Schmettow
#' @export


rescale_centered <- function(x, scale = .999){
	mean_x <- mean(x, na.rm = TRUE)
	x_center <- x - mean_x
	x_shrink <- x_center * scale
	out <- x_shrink + mean_x
	out
}

#' rescaling to unit interval
#'
#' scales a variable to the interval [0;1]
#'
#' @param x numerical vector
#' @param lower lowest possible value
#' @param upper highest possible value
#' @param scale rescale factor
#' @return numerical vector
#'
#'
#' @author Martin Schmettow
#' @export


rescale_unit     <- function(x,
														 lower = min(x, na.rm = TRUE),
														 upper = max(x, na.rm = TRUE), scale = 1){
	x_to_zero <- x - lower
	x_to_one <- x_to_zero/(upper -lower)
	out <- rescale_centered(x_to_one, scale = scale)
	out <- x_to_one
	out
}

#' @rdname rescale_unit
#' @export

rescale_zero_one <- rescale_unit
## TODO:
# check out dplyr:::select_.data.frame for how to use dplyrs column expansion



#' Removing completely unmeasured variables
#'
#' all columns that are completely NA are discarded
#'
#' @param D data frame
#' @return data frame
#'
#'
#' @author Martin Schmettow
#' @export

discard_all_na <-
	function(D){
		select(D, where(any_not_na))
	}

all_na <-
	function(x) all(is.na(x))

any_not_na <-
	function(x) !all_na(x)




#' Re-order levels of a factor
#'
#' Re-orders the levels of a factor by a vector of new positions.
#'
#' @param x factor
#' @param positions vector of positions
#' @return factor
#'
#'
#' @author Martin Schmettow
#' @export


reorder_levels <- function(x, positions){
	levels <- unique(x)
	if(length(positions) != length(levels))
		stop(stringr::str_c("Number of positions (",length(positions),")
                           does not match number of levels (", length(levels) ,")"))
	missing <- setdiff(1:length(levels), positions)
	above <- setdiff(positions, 1:length(levels))
	if(length(missing) | length(above)) {
		stop(stringr::str_c("Invalid positions vector. Must contain exactly the values 1:", length(levels)))
	}
	out <- factor(x, levels = levels[positions])
	out
}

# x <- c("A", "B", "C")
# reorder_levels(x, c(3,2,4))

# discard_redundant <- function(D){
#   if(nrow(D) < 2) return(D)
#
#   keep <- purrr::map_lgl(D, ~ n_distinct(.x) > 1)
#   D[, keep]
# }


#' Removes variables that do not vary in value
#'
#' @param D data frame
#' @param except vector of column names to keep
#' @return data frame
#'
#' @export
#' @author Martin Schmettow

#' all variables that have a constant value are removed from a data frame
#'
#' @param D data frame
#' @param except vector of column names to keep
#' @param ... ignored
#' @return data frame
#'
#'
#' @author Martin Schmettow
#' @export

discard_redundant <-
	function(D, except, ...) UseMethod("discard_redundant", D)


#' @rdname discard_redundant
#' @export
#'
discard_redundant.default <- function(D, except = c(), ...){
	if(nrow(D) < 2) return(D)
	varies <- purrr::map_lgl(D, ~ n_distinct(.x) > 1)
	keep   <- names(D) %in% except | varies
	select(D, all_of(names(D)[keep]))
}


#' @rdname discard_redundant
#' @export
discard_redundant.tbl_clu <- function(D, except = c(), ...)
	as_tibble(D) %>% discard_redundant(except = c(except, "parameter", "center", "lower", "upper"))

#' @rdname discard_redundant
#' @export
discard_redundant.tbl_coef <- function(D, except = c(), ...)
	as_tibble(D) %>% discard_redundant(except = c(except, "parameter", "center", "lower", "upper"))

#' @rdname discard_redundant
#' @export
discard_redundant.tbl_post_pred <- function(D, except = c(), ...)
	as_tibble(D) %>% discard_redundant(except = c(except, "Obs","value"))

#' @rdname discard_redundant
#' @export
discard_redundant.tbl_predicted <- function(D, except = c(), ...)
	as_tibble(D) %>% discard_redundant(except = c(except, "Obs", "center", "lower", "upper"))

#' @rdname discard_redundant
#' @export
discard_redundant.tbl_post <- function(D, except = c(), ...)
	as_tibble(D) %>% discard_redundant(except = c(except, "parameter", "value"))

