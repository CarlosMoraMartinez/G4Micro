#' @title Zero Imputation of Compositional Data Using zCompositions
#' @description
#' Filters taxa based on prevalence and imputes zeros using the multiplicative replacement method from the \pkg{zCompositions} package.
#'
#' @param phobj A phyloseq object containing the OTU table.
#' @param prevalence_lim Numeric threshold for prevalence filtering.
#' Taxa with prevalence (proportion of samples with count > 0)
#' below this value are removed. Default is 0.05.
#' @return A numeric matrix with zeros imputed using multiplicative replacement, with taxa in rows and samples in columns.
#'
#' @details
#' This function filters taxa based on prevalence, then imputes zeros using \code{cmultRepl} from \pkg{zCompositions}.
#' Note that \code{cmultRepl} expects samples as rows and taxa as columns, so the OTU table is transposed accordingly before and after imputation.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   imputed_matrix <- getZCompositionImputedData(phobj, prevalence_lim = 0.1)
#' }
#' }
#' @seealso
#'  \code{\link[zCompositions]{cmultRepl}}
#' @rdname getZCompositionImputedData
#' @export
#' @importFrom zCompositions cmultRepl
getZCompositionImputedData <- function(phobj, prevalence_lim = 0.05){

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
