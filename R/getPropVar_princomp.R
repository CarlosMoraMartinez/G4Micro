#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param pca PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
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
