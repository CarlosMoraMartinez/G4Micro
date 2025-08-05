#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param resdf PARAM_DESCRIPTION
#' @param plim PARAM_DESCRIPTION, Default: 0.05
#' @param fc PARAM_DESCRIPTION, Default: 1
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{filter}}
#' @rdname filter_taxa_praw
#' @export 
#' @importFrom dplyr filter
filter_taxa_praw <- function(resdf, plim=0.05, fc=1){
  taxa <- resdf %>%
    dplyr::filter(pvalue <= plim & abs(log2FoldChangeShrink) >= log2(fc) ) %>%
    pull(taxon)
  return(taxa)
}
