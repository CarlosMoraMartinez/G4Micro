#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param otutab PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname calc_scImpute
#' @export 
calc_scImpute <- function(otutab){
  xx <- otutab %>% as.data.frame()
  write.table(xx, file="test_input.csv", sep=",", quote=F, row.names = T)
  library(scImpute)
  outlier_samples <- scimpute(# full path to raw count matrix
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
