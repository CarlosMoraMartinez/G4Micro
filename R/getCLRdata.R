#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param phobj PARAM_DESCRIPTION
#' @param prevalence_lim PARAM_DESCRIPTION, Default: 0.05
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname getCLRdata
#' @export 
getCLRdata <- function(phobj, prevalence_lim = 0.05){
  library(CoDaSeq)
  #1) Imputate zeros
  # cmultRepl expects samples in rows and taxa in columns
  tax_matrix <- getZCompositionImputedData(phobj, prevalence_lim)

  #2) CLR transform
  tax_matrix_clr <- codaSeq.clr(tax_matrix, samples.by.row=F, IQLR=TRUE)

  return(tax_matrix_clr)
}
