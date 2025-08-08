#' @title CLR Transformation of Compositional Data with Zero Imputation
#' @description
#' Imputes zeros in compositional microbiome data and performs Centered Log-Ratio (CLR) transformation.
#'
#' @param phobj Phyloseq object containing abundance data.
#' @param prevalence_lim Numeric threshold for prevalence filtering to identify zeros to impute. Default is 0.05.
#' @return A numeric matrix with CLR-transformed data after zero imputation.
#'
#' @details
#' The function imputes zeros using \code{getZCompositionImputedData},
#' then applies the CLR transformation using the \code{codaSeq.clr} function from the \pkg{CoDaSeq} package.
#' Note that \code{codaSeq.clr} expects taxa in rows, so \code{samples.by.row} is set to FALSE.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Assume phobj is a phyloseq object already loaded
#'   clr_data <- getCLRdata(phobj, prevalence_lim = 0.1)
#' }
#' }
#' @seealso \code{\link[CoDaSeq]{codaSeq.clr}}, \code{\link{getZCompositionImputedData}}
#' @rdname getCLRdata
#' @export
#' @importFrom CoDaSeq codaSeq.clr
getCLRdata <- function(phobj, prevalence_lim = 0.05){
  library(CoDaSeq)
  #1) Imputate zeros
  # cmultRepl expects samples in rows and taxa in columns
  tax_matrix <- getZCompositionImputedData(phobj, prevalence_lim)

  #2) CLR transform
  tax_matrix_clr <- codaSeq.clr(tax_matrix, samples.by.row=F, IQLR=TRUE)

  return(tax_matrix_clr)
}
