#### ASSERTIONS ####

assert_has_name <- function(var, x){
	x_name <- as.character(deparse(substitute(x)))
	assert_that(has_name(x, var),
							msg = str_c("Variable ", var," missing in ", x_name))
}

assert_names <- function(x, ...){
	var_names <- as.character(rlang::exprs(...))
	walk(var_names, assert_has_name, x)
}

assert_key <- function(x, ...){
	var_names <- as.character(rlang::exprs(...))
	assert_that(length(var_names) > 0,
							msg = "At least one key variable must be provided as argument")
	assert_names(x, ...)
	assert_that(nrow(x) == nrow(distinct(x, ...)),
							msg = str_c("<",str_c(var_names, collapse = ","),
													"> not unique. Duplicate combinations found."))
}
