## dplyr is used with NSE, which gives "no visible binding for global variable
## errors". All such global variables are declared here.
utils::globalVariables(names = c(
	".chain", ".draw", ".iteration", ".tmp_idx",
	"Estimate", "Model", "Obs", "Part", "SD", "SE",
	"center", "chain", "conf.high", "conf.low",
	"diff_IC", "dpar", "effect", "estimate", "fe_value", "fixef_2",
	"group", "iter", "lower",
	"model", "new_name", "nlpar", "nonlin",
	"parameter", "pattern", "prior",
	"re_1", "re_entity", "re_factor",
	"term", "tbl_coef", "tidy", "type",
	"upper", "value"
))
