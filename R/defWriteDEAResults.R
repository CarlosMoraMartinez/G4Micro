#' @title Write DESeq2 Differential Expression Results to File
#' @description Converts DESeq2 results and log fold change shrinkage results into a data frame,
#' adds shrinkage statistics, writes the table to a tab-separated file, and returns the data frame.
#' @param res A DESeq2 results object (typically from \code{DESeq2::results}).
#' @param resLFC A DESeq2 log fold change shrinkage results object (e.g., from \code{lfcShrink}).
#' @param opt A list of options, must include \code{out} specifying the output directory.
#' @param name A string used as the filename (without extension) for saving the results.
#' @param nested_dir Subdirectory within \code{opt$out} to write output.
#' @return A data frame containing the original DESeq2 results along with added columns for
#' log2 fold change shrinkage, standard error of shrinkage, and s-values. The same data frame
#' is also written to a tab-separated file at the specified output path.
#' @details The function converts the DESeq2 results object to a data frame, adds shrinkage
#' estimates from \code{resLFC}, then writes the combined table to a file named
#' \code{<opt$out>/<name>} (no file extension added). The returned data frame includes
#' columns: \code{taxon}, original DESeq2 statistics, \code{log2FoldChangeShrink},
#' \code{lfcSE_Shrink}, and \code{svalue}.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   res <- DESeq2::results(dds)
#'   resLFC <- DESeq2::lfcShrink(dds, coef=2, type="apeglm")
#'   opt <- list(out="results")
#'   defWriteDEAResults(res, resLFC, opt, "comparison1.tsv")
#' }
#' }
#' @seealso
#'  \code{\link[dplyr]{mutate}},
#'  \code{\link[DESeq2]{results}},
#'  \code{\link[DESeq2]{lfcShrink}}
#' @rdname defWriteDEAResults
#' @export
#' @importFrom dplyr mutate
defWriteDEAResults <- function(res, resLFC, opt, name, nested_dir = ""){
  if(nested_dir != ""){
    if(! dir.exists(paste0(opt$out, nested_dir))) dir.create(paste0(opt$out, nested_dir))
  }
  fname <-paste(opt$out, nested_dir, name, sep="/", collapse="/")
  resdf <- res %>% as.data.frame(row.names = rownames(.)) %>%
    rownames_to_column("taxon") %>%
    dplyr::mutate(log2FoldChangeShrink = resLFC$log2FoldChange,
                  lfcSE_Shrink = resLFC$lfcSE, svalue = resLFC$svalue)
  write_tsv(resdf, fname)
  return(resdf)
}
