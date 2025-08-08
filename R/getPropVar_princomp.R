#' @title Generate Labels of Variance Explained from a princomp Object
#' @description Creates character labels for each principal component indicating the percentage of variance explained, based on a `princomp` PCA object.
#' @param pca A `princomp` class object resulting from PCA analysis.
#' @return A named character vector with labels like "PC1 (xx.xx%)" representing variance explained by each component.
#' @details The function calculates the proportion of variance explained by each principal component, rounds it to two decimals, and formats the result into readable labels for visualization or reporting.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   data(USArrests)
#'   pca <- princomp(USArrests, cor=TRUE)
#'   getPropVar_princomp(pca)
#' }
#' }
#' @rdname getPropVar_princomp
#' @export
getPropVar_princomp <- function(pca){
  pcts <- round(100*pca$sdev^2/sum(pca$sdev^2), 2) %>% as.character
  nnames <- gsub("Comp.", "PC", names(pca$sdev))
  labb <- paste(nnames, " (", pcts, "%)", sep="")
  names(labb)<- nnames
  return(labb)
}
