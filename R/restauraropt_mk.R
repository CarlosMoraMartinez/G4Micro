#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param opt PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname restauraropt_mk
#' @export 
restauraropt_mk <- function(opt){
  output <- opt$out
  restaurar <- function(optbad){
    optbad$out <- output
    return(opt)
  }
}
