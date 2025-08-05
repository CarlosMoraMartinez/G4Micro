#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param phobj PARAM_DESCRIPTION
#' @param vars PARAM_DESCRIPTION, Default: c("nreads")
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname updatePsWithLogs
#' @export 
updatePsWithLogs <- function(phobj, vars = c("nreads")){

  for(v in vars){
    newname <- paste0(v, "_log")
    sample_data(phobj)[[newname]] <- log(sample_data(phobj)[[v]] +1  )
  }
  return(phobj)
}
