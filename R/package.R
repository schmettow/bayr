library(tidyverse)

## dplyr is used with NSE, which gives "no visible binding for global variable errors"
utils::globalVariables(names = c("type", "parameter", "value",
																 "new_name", "iter", "pattern"))

## for testing purposes using the NewStats infrastructure
path_NewStats <- str_c("../../../../../Publications/New_Stats/")
path_Cases <- str_c(path_NewStats, "Cases/")
