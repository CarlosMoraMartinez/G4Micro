#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mat PARAM_DESCRIPTION
#' @param opt PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname defWriteMatAsDF
#' @export 
defWriteMatAsDF <- function(mat, opt, name){
  fname <-paste(opt$out, name, sep="/", collapse="/")
  df <- mat %>% as.data.frame(row.names = rownames(.)) %>%
    rownames_to_column("gene")
  write_tsv(df, fname)
  return(df)
}
