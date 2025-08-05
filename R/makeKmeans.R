#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param datasc PARAM_DESCRIPTION
#' @param levs PARAM_DESCRIPTION
#' @param varnames PARAM_DESCRIPTION
#' @param SEED PARAM_DESCRIPTION, Default: 123
#' @param folds PARAM_DESCRIPTION, Default: c()
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{select}}
#' @rdname makeKmeans
#' @export 
#' @importFrom dplyr select
makeKmeans <- function(datasc, levs, varnames, SEED=123, folds=c()){
  library(stats)
  library(caret)
  set.seed(SEED)
  train_df <- datasc %>% dplyr::select(-class, -sample)  %>% dplyr::select(all_of(varnames))
  mod_kmeans <- kmeans(train_df, centers=length(levs), iter.max = 100, nstart=100)
  predict_kmeans <-levels(datasc$class)[mod_kmeans$cluster] %>% factor(levels=levs)
  confmat_kmeans <- confusionMatrix(predict_kmeans, datasc$class, positive = levs[2])
  return(list(confmat_no_l1o=confmat_kmeans, mod=mod_kmeans, predicted=predict_kmeans))
}
