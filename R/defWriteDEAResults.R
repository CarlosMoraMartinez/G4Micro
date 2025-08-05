#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param res PARAM_DESCRIPTION
#' @param resLFC PARAM_DESCRIPTION
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
#' @seealso 
#'  \code{\link[dplyr]{mutate}}
#' @rdname defWriteDEAResults
#' @export 
#' @importFrom dplyr mutate
defWriteDEAResults <- function(res, resLFC, opt, name){
  fname <-paste(opt$out, name, sep="/", collapse="/")
  resdf <- res %>% as.data.frame(row.names = rownames(.)) %>%
    rownames_to_column("taxon") %>%
    dplyr::mutate(log2FoldChangeShrink = resLFC$log2FoldChange,
                  lfcSE_Shrink = resLFC$lfcSE, svalue = resLFC$svalue)
  write_tsv(resdf, fname)
  return(resdf)
}
