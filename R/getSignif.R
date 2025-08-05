#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param vector PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname getSignif
#' @export 
getSignif <- function(vector){
  x <- ifelse(vector < 0.001, "***",
              ifelse(vector < 0.01, "**",
                     ifelse(vector< 0.05, "*", "NS")
              )
  )
  return(x)
}
