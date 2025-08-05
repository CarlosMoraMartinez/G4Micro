#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param namelist PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname process_names
#' @export 
process_names <- function(namelist){
  gsub("-", ".", namelist) %>%
    gsub("\\[|\\]|\\(|\\)", "", ., perl=T)
}
