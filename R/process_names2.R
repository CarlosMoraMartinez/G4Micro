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
#' @rdname process_names2
#' @export 
process_names2 <- function(namelist){

  res <- gsub("_", " ", namelist) %>%
    sapply(\(x){
      if(x=="") return("")
      x <- strsplit(x, " ")[[1]]
      y <- toupper(strsplit(x[1], "")[[1]][1])
      y <- paste0(y, ". ", x[2])
      return(y)
    })
}
