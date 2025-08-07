#' @title K-means Clustering with Confusion Matrix Evaluation
#' @description
#' Performs K-means clustering on the given dataset using specified variables,
#' then evaluates cluster assignments against true class labels via a confusion matrix.
#'
#' @param datasc A data frame containing the dataset. Must include a factor column `class` with true class labels, and `sample`.
#' @param levs A character vector specifying the class levels (clusters) to identify; determines number of clusters and factor levels.
#' @param varnames A character vector of column names in `datasc` to use as features for clustering.
#' @param SEED An integer seed for reproducibility of clustering results. Default is 123.
#' @param folds Currently not used in the function; reserved for compatibility. Default is an empty vector.
#'
#' @return A list containing:
#' \describe{
#'   \item{confmat_no_l1o}{Confusion matrix comparing K-means cluster assignments to true class labels.}
#'   \item{mod}{The fitted `kmeans` clustering model object.}
#'   \item{predicted}{Factor vector of cluster assignments labeled by class levels.}
#' }
#'
#' @details
#' The function performs K-means clustering with the number of clusters equal to the number of class levels.
#' The cluster labels are mapped to the factor levels of the true classes to generate a confusion matrix.
#' Note that cluster labels are arbitrary and may not correspond to actual class labels in order.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   data(iris)
#'   iris$class <- factor(ifelse(iris$Species == "setosa", "A", "B"))
#'   iris$sample <- 1:nrow(iris)
#'   result <- makeKmeans(iris, levs=c("A", "B"), varnames=colnames(iris)[1:4])
#'   print(result$confmat_no_l1o)
#' }
#' }
#'
#' @seealso
#' \code{\link[stats]{kmeans}}, \code{\link[caret]{confusionMatrix}}, \code{\link[dplyr]{select}}
#'
#' @rdname makeKmeans
#' @export
#' @importFrom dplyr select
#' @importFrom caret confusionMatrix
#' @importFrom stats kmeans
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
