#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param pca PARAM_DESCRIPTION
#' @param PC PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname getPropVar
#' @export 
getPropVar <- function(pca, PC){
  xx <- summary(pca)
  percc<- round(100*xx$importance["Proportion of Variance" , PC],2) %>% as.character
  labb <- paste(PC, " (", percc, "%)", sep="", collapse="")
  return(labb)
}
