#' @title Run scImpute on OTU Table
#' @description Performs dropout imputation on count data using the scImpute package.
#'
#' The function writes the OTU count table to a CSV file, runs scImpute with specified parameters,
#' and reads back the imputed count matrix. It returns the outlier samples identified and the imputed counts.
#'
#' @param otutab A count matrix or data frame of raw counts (e.g., OTU table) with samples as columns and features as rows.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{outliers}{Output from the scImpute function, typically sample or cell outlier information.}
#'   \item{imputed_counts}{A data frame of counts after imputation by scImpute.}
#' }
#'
#' @details
#' The function writes the input counts to a temporary CSV file `"test_input.csv"`, runs scImpute with default parameters
#' (including dropout threshold 0.5, Kcluster=2, and 4 cores), and reads the imputed counts from `"scimpute_count.csv"`.
#'
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Assuming 'counts' is your OTU count matrix:
#'   result <- calc_scImpute(counts)
#'   head(result$imputed_counts)
#' }
#' }
#'
#' @rdname calc_scImpute
#' @importFrom scImpute scimpute
#' @export
calc_scImpute <- function(otutab){
  xx <- otutab %>% as.data.frame()
  write.table(xx, file="test_input.csv", sep=",", quote=F, row.names = T)

  outlier_samples <- scImpute::scimpute(# full path to raw count matrix
    count_path ="test_input.csv",
    infile = "csv",           # format of input file
    outfile = "csv",          # format of output file
    out_dir = "./",           # full path to output directory
    labeled = FALSE,          # cell type labels not available
    drop_thre = 0.5,          # threshold set on dropout probability
    Kcluster = 2,             # 2 cell subpopulations
    ncores = 4)
  scimputed <- read.csv("scimpute_count.csv")
  return(list(outliers = outlier_samples, imputed_counts=scimputed))
}
