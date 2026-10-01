## dplyr is used with NSE, which gives "no visible binding for global variable
## errors". All such global variables are declared here.
utils::globalVariables(names = c("type", "parameter", "value",
																 "new_name", "iter", "pattern", "tbl_coef"))
