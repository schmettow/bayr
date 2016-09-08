
## dplyr is used with NSE, which gives "no visible binding for global variable errors"
utils::globalVariables(names = c("type", "parameter", "value",
																 "new_name", "iter", "pattern"))

#' formula interface for MCMCglmm
#'
#' returns the formula of an MCMCglmm object
#'
#' @usage formula()
#' @param model MCMCglmm object
#'
#' only fixed part implemented by now
#'
#' @author Martin Schmettow
#' @export


formula.MCMCglmm <-
	function(model) model$Fixed$formula
