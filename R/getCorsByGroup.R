#' @title Compute Correlations for Two Groups Separately and Combined
#' @description
#' Calculates correlations between two matrices overall and within two specified groups,
#' then merges results and writes a combined table to disk.
#'
#' @param mat1 Numeric matrix or data.frame with samples as rows and variables as columns.
#' @param mat2 Numeric matrix or data.frame with samples as rows and variables as columns.
#' @param group1 Logical or integer vector indicating the subset of rows belonging to the first group.
#' @param group2 Logical or integer vector indicating the subset of rows belonging to the second group.
#' @param outdir Output directory path (should end with '/' if not empty).
#' @param name Output file name (e.g., "correlations.tsv").
#' @param levnames Character vector of length 2 specifying group names for labeling results. Default is c("PSO", "CTRL").
#' @return A list containing:
#'   \item{alldf}{Data frame combining overall and per-group correlation results.}
#'   \item{all}{Output of getCorrelations for all samples combined.}
#'   \item{gr1}{Output of getCorrelations for group1.}
#'   \item{gr2}{Output of getCorrelations for group2.}
#'
#' @details
#' The function calls \code{getCorrelations} for all samples together and separately for
#' each group defined by \code{group1} and \code{group2}. It then merges the results,
#' renames columns for clarity, writes a TSV file, and returns all data.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Assume mat1 and mat2 are matrices with matching rownames representing samples
#'   group1 <- sample(1:nrow(mat1), size = floor(nrow(mat1)/2))
#'   group2 <- setdiff(1:nrow(mat1), group1)
#'   results <- getCorsByGroup(mat1, mat2, group1, group2, outdir = "./results/", name = "cors.tsv")
#' }
#' }
#' @seealso \code{\link[dplyr]{select}}, \code{\link{getCorrelations}}
#' @rdname getCorsByGroup
#' @export
#' @importFrom dplyr select
getCorsByGroup <- function(mat1, mat2, group1, group2, outdir, name, levnames=c("PSO", "CTRL")){

  cors <- getCorrelations(mat1, mat2)
  cors_pso <- getCorrelations(mat1[group1,], mat2[group1,])
  cors_ct <- getCorrelations(mat1[group2,], mat2[group2,])

  cors_pso_df <- cors_pso[[1]] %>% dplyr::select(-var1, -var2)
  names(cors_pso_df) <- paste(levnames[1], names(cors_pso_df), sep="_")

  cors_ct_df <- cors_ct[[1]] %>% dplyr::select(-var1, -var2)
  names(cors_ct_df) <- paste(levnames[2], names(cors_ct_df), sep="_")


  cordf <- cbind(cors[[1]], cors_pso_df, cors_ct_df)

  oname <- paste0(outdir, name)
  write_tsv(cordf, file = oname)
  return(list(alldf=cordf, all=cors, gr1=cors_pso, gr2=cors_ct))
}
