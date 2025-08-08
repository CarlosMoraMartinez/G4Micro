#' @title Perform Canonical Correlation Analysis (CCA) Using Vegan on CLR-transformed Data
#' @description
#' Conducts a canonical correlation analysis (CCA) between a data matrix and
#' principal components derived from a CLR-transformed taxonomy matrix.
#' This approach reduces dimensionality of the taxonomy matrix before
#' performing CCA, retaining components explaining a specified proportion
#' of variance.
#'
#' @param datamat Numeric matrix or data frame of samples (rows) by features (columns),
#' usually representing environmental or phenotypic variables.
#' @param tax_matrix_clr Numeric matrix of CLR-transformed taxonomic abundances,
#' with taxa/features as rows and samples as columns.
#' @param pcvar2retain Numeric value between 0 and 1 indicating the cumulative proportion
#' of variance in the taxonomy PCA to retain for downstream CCA. Default is 0.9.
#' @param outdir Character string specifying the output directory where results
#' (RData and PDF biplot) will be saved. Default is current working directory ("").
#' @param name Character string base name for output files. Default is "CCA_from_CLR".
#'
#' @return An object of class \code{CCorA} (from the \pkg{vegan} package)
#' containing results of the canonical correlation analysis.
#'
#' @details
#' The function assumes that \code{datamat} and \code{tax_matrix_clr} have the same
#' number of rows (samples) and are ordered identically. It performs PCA on the
#' transposed CLR taxonomic matrix to reduce dimensionality, retaining the number
#' of principal components needed to reach the specified variance proportion.
#' CCA is then performed between \code{datamat} and these principal components.
#' Results are saved as an RData file, and a biplot of the CCA is saved as a PDF.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Simulated example with dummy data:
#'   datamat <- matrix(rnorm(100), nrow=10, ncol=10)
#'   tax_matrix_clr <- matrix(rnorm(200), nrow=20, ncol=10)
#'   results <- makeCCAfromMatrixVegan(datamat, tax_matrix_clr, pcvar2retain=0.9, outdir=".", name="TestCCA")
#' }
#' }
#'
#' @seealso \code{\link[vegan]{CCorA}}, \code{\link[stats]{prcomp}}, \code{\link[graphics]{biplot}}
#' @rdname makeCCAfromMatrixVegan
#' @importFrom vegan CCorA
#' @importFrom stats prcomp biplot
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
