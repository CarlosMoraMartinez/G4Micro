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
#' @seealso 
#'  \code{\link[zCompositions]{cmultRepl}}
#' @rdname getZCompositionImputedData
#' @export 
#' @importFrom zCompositions cmultRepl
getZCompositionImputedData <- function(phobj, prevalence_lim = 0.05){
  library(zCompositions)
  otutab <- otu_table(phobj)
  #Umbral de prevalencia
  prevalence <- apply(otutab, MAR=1, FUN=function(x)sum(x>0)/length(x))
  otutab_filt <- otutab[prevalence > prevalence_lim ,]

  # CLR transform as in Xia, Sun & Xen book, pag 350

  #1) Imputate zeros
  # cmultRepl expects samples in rows and taxa in columns
  tax_matrix <- zCompositions::cmultRepl(X = t(otutab_filt), output="p-counts") %>% t

  return(tax_matrix)
}
