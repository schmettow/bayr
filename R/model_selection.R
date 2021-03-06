#' Information Criteria extraction
#'
#' Tidies the output of model evaluation commands from package Loo.
#'
#' @usage IC(ic)
#' @param ic loo.psis, kfold or waic object
#' @return tbl_df
#' @author Martin Schmettow
#' @import dplyr
#' @importFrom knitr knit_print
#' @export


IC <- function (ic){
	tbl_IC <-
		ic$estimates %>%
		as_tibble(rownames = "IC") %>%
		mutate(Model = attr(ic, "model_name")) %>%
		select(Model, IC, Estimate, SE)
	class(tbl_IC) <- append("tbl_IC", class(tbl_IC))
	tbl_IC
}

#' Information Criteria comparison
#'
#' Collects results from multiple Loo objects and
#' puts them into a tidy model comparison table.
#'
#' @usage compare_IC(ic_list)
#' @param ic loo.psis, kfold or waic object
#' @return tbl_df
#' @author Martin Schmettow
#' @import dplyr
#' @importFrom knitr knit_print
#' @export


compare_IC <- function(ic_list){
	tbl_IC_comp <-
		ic_list %>%
		purrr::map_df(IC) %>%
		filter(IC %in% c("looic", "waic", "kfoldic")) %>%
		group_by(IC) %>%
		mutate(diff_IC = Estimate - min(Estimate)) %>%
		ungroup() %>%
		arrange(IC, diff_IC)
	class(tbl_IC_comp) <- append("tbl_IC_comp", class(tbl_IC_comp))
	tbl_IC_comp
}
