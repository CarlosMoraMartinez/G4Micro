#' @title Run Regularized Canonical Correlation (RCC) with mixOmics
#' @description Performs a regularized canonical correlation analysis (RCC) between microbiome abundance data
#'              and dietary data using the \code{\link[mixOmics]{rcc}} function.
#'
#' @param microbiome_mat A numeric matrix or data frame of microbiome features (e.g., taxa abundances),
#'                       with rows corresponding to samples and columns to taxa.
#' @param diet_mat A numeric matrix or data frame of dietary variables, with rows corresponding to samples
#'                 and columns to dietary features. Missing values are replaced with zeros before analysis.
#' @param tax2plot A character vector of taxa names (matching column names in \code{microbiome_mat})
#'                 to include in the analysis. Default: \code{c()} (all taxa are included).
#' @param meta2plot A character vector of dietary variables (matching column names in \code{diet_mat})
#'                  to include in the analysis. Default: \code{c()} (all variables are included).
#' @param outdir A character string specifying the output directory where results will be saved.
#'               Default: \code{'./'}.
#' @param name A character string specifying a label for the output file name.
#'             Default: \code{'all'}.
#'
#' @return An object of class \code{rcc} from \pkg{mixOmics}, containing the results of the canonical correlation analysis.
#'
#' @details The function filters out samples with only missing or zero dietary information,
#'          scales both microbiome and dietary matrices, removes microbiome features containing missing values
#'          after scaling, and then applies \code{\link[mixOmics]{rcc}} with the "shrinkage" method.
#'          The resulting RCC object is also saved as an \code{.RData} file in the specified \code{outdir}.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Example data
#'   set.seed(123)
#'   microbiome <- matrix(rpois(50, 10), nrow = 10, ncol = 5)
#'   colnames(microbiome) <- paste0("Taxon", 1:5)
#'   rownames(microbiome) <- paste0("Sample", 1:10)
#'
#'   diet <- matrix(runif(30, 0, 100), nrow = 10, ncol = 3)
#'   colnames(diet) <- paste0("DietVar", 1:3)
#'   rownames(diet) <- paste0("Sample", 1:10)
#'
#'   # Run RCC
#'   res <- makeRCC_MixOmics(microbiome, diet, tax2plot = c("Taxon1", "Taxon3"), name = "example")
#'   print(res)
#' }
#' }
#'
#' @seealso
#'  \code{\link[assertthat]{assert_that}},
#'  \code{\link[mixOmics]{rcc}}
#' @rdname makeRCC_MixOmics
#' @export
#' @importFrom assertthat assert_that
#' @importFrom mixOmics rcc
akeRCC_MixOmics <- function(microbiome_mat, diet_mat, tax2plot=c(), meta2plot = c(),
                             outdir="./", name="all"){

  assertthat::assert_that(all(tax2plot %in% colnames(microbiome_mat)))
  #Some NAs are already encoded as 0s
  diet_mat[is.na(diet_mat)] <- 0
  not_nas <- rowSums(diet_mat != 0) > 0
  not_nas <- names(not_nas)[not_nas > 0]

  microbiome_mat_filt <- microbiome_mat[not_nas, ]
  if(length(tax2plot)>0){
    microbiome_mat_filt <- microbiome_mat[, tax2plot]
  }

  diet_mat_filt <- diet_mat[not_nas, ]
  if(length(meta2plot)>0){
    diet_mat_filt <- diet_mat_filt[, meta2plot]
  }

  microbiome_mat_filt_scale <- scale(microbiome_mat_filt)
  with_nas <- apply(microbiome_mat_filt_scale, MAR=2, \(x)length(which(is.na(x)))) %>% sort
  with_nas <- with_nas[with_nas > 0]
  microbiome_mat_filt_scale <- microbiome_mat_filt_scale[, ! colnames(microbiome_mat_filt_scale) %in% names(with_nas)]

  diet_mat_filt_scale <- scale(diet_mat_filt)

  rcc_res <- mixOmics::rcc(microbiome_mat_filt_scale, diet_mat_filt_scale, method = "shrinkage")
  save(rcc_res, file = paste0(outdir, "mixomics_RCC_", name, ".RData"))

  return(rcc_res)
}
