#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param tab PARAM_DESCRIPTION
#' @param namedtab PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname nameProcesses
#' @export 
nameProcesses<- function(tab, namedtab){
  if(is.null(namedtab)){return(tab)}
  tab$process_name <- namedtab$long[match(rownames(tab), namedtab$short)]
  return(tab)
}
