#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param datamat PARAM_DESCRIPTION
#' @param tax_matrix_clr PARAM_DESCRIPTION
#' @param pcvar2retain PARAM_DESCRIPTION, Default: 0.9
#' @param outdir PARAM_DESCRIPTION, Default: ''
#' @param name PARAM_DESCRIPTION, Default: 'CCA_from_CLR'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname makeCCAfromMatrixVegan
#' @export 
makeCCAfromMatrixVegan <- function(datamat, tax_matrix_clr, pcvar2retain=0.9, outdir = "", name="CCA_from_CLR"){
  #Assumes matrices have the same number of rows (subjects) and that are ordered in the same way
  txpca <- prcomp(tax_matrix_clr %>% t)
  txsum <- txpca %>% summary()
  pc2use <- colnames(txsum$importance)[txsum$importance["Cumulative Proportion",] > pcvar2retain][1]
  txpcdata <- txpca$rotation[, 1:(which(colnames(txpca$rotation) == pc2use))]
  ccares <- CCorA(datamat, txpcdata)

  save(ccares, file = paste0(outdir, name, ".RData"))
  pdf(paste0(outdir, name, "_BiplotVegan.pdf"), width = 12, height = 12)
  biplot(ccares)
  dev.off()
  return(ccares)
}
