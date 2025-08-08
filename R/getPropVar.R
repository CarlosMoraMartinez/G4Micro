#' @title getPropVar
#' @description Returns a label with the percentage of variance explained by one or more principal components from a PCA object.
#' @param pca PCA object (result of prcomp or similar) from which the proportion of explained variance is extracted.
#' @param PC Numeric vector indicating the principal component(s) to query (e.g., 1 for PC1, 2 for PC2).
#' @return A character vector with labels for each specified principal component in the format "PCn (xx.xx%)".
#' @details The function extracts the proportion of variance explained from the importance matrix in the summary of the PCA object, rounds it to two decimals, and constructs a readable label for plotting or reporting.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   pca <- prcomp(USArrests, scale. = TRUE)
#'   getPropVar(pca, 1)        # "1 (62.00%)"
#'   getPropVar(pca, c(1,2))   # "1 (62.00%)" "2 (24.74%)"
#' }
#' }
#' @rdname getPropVar
#' @export
getPropVar <- function(pca, PC){
  xx <- summary(pca)
  percc<- round(100*xx$importance["Proportion of Variance" , PC],2) %>% as.character
  labb <- paste(PC, " (", percc, "%)", sep="", collapse="")
  return(labb)
}
