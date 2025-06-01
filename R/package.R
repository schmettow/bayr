library(tidyverse)

## dplyr is used with NSE, which gives "no visible binding for global variable errors"
utils::globalVariables(names = c("type", "parameter", "value",
																 "new_name", "iter", "pattern","tbl_coef"))
utils::globalVariables(names = c("type", "parameter", "value",
																 "new_name", "iter", "pattern","tbl_coef"))

## for testing purposes using the NewStats infrastructure
## path_NewStats <- stringr::str_c("../../../../../Publications/New_Stats/")
## path_Cases <- stringr::str_c(path_NewStats, "Cases/")
