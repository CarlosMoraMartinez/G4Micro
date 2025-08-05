#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param vector PARAM_DESCRIPTION
#' @param limits PARAM_DESCRIPTION, Default: c(0.05, 0.01, 0.001)
#' @param legends PARAM_DESCRIPTION, Default: c("", "*", "**", "***")
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname getSignif2
#' @export 
getSignif2 <- function(vector, limits=c(0.05, 0.01, 0.001),
                       legends = c("", "*", "**", "***")){
  x <- ifelse(vector < limits[3], legends[4],
              ifelse(vector < limits[2], legends[3],
                     ifelse(vector< limits[1], legends[2], legends[1])
              )
  )
  return(x)
}
