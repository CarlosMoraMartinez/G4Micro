#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param var PARAM_DESCRIPTION
#' @param group PARAM_DESCRIPTION
#' @param sep PARAM_DESCRIPTION, Default: '-'
#' @param paired PARAM_DESCRIPTION, Default: F
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{arrange}}
#' @rdname getTestsForAllCombinations
#' @export 
#' @importFrom dplyr arrange
getTestsForAllCombinations <- function(var, group, sep="-", paired=F){
  comp <- combn(unique(group[!is.na(group)]), 2, simplify = F)
  res <- lapply(comp, function(combo){
    a <- var[group==combo[1]]
    b <- var[group==combo[2]]
    res <- data.frame(groups_compared = paste(combo, sep=sep, collapse=sep),
                      t_pval = t.test(a, b, paired=paired)$p.value,
                      wilcox_pval = wilcox.test(a, b, paired=paired)$p.value,
                      shapiro_test =  shapiro.test(c(a, b))$p.value,
                      bartlett_test = bartlett.test(var[group %in% combo],
                                                    group[group %in% combo])[[3]],
                      levene_test = leveneTest(var[group %in% combo],
                                               factor(group[group %in% combo]))[1, 3]
    )
  }) %>% bind_rows  %>% dplyr::arrange(groups_compared)
  return(res)
}
