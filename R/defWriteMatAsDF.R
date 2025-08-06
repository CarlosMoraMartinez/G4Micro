#' @title Write Matrix as Data Frame and Save as TSV File
#' @description Converts a matrix to a data frame with row names as a column, writes it to a TSV file, and returns the data frame.
#' @param mat A \code{DESeqResults} object, typically returned by \code{DESeq2::results()}, containing differential expression results.
#' @param opt A list of options containing an output directory path in \code{opt$out}.
#' @param name A character string specifying the output file name.
#' @return A data frame with the original matrix data and a "gene" column containing row names.
#' @details This function converts the input matrix to a data frame, adds a "gene" column with the original row names, writes the data frame to a tab-separated file in the specified output directory, and returns the data frame.
#'@examples
#' \dontrun{
#' if(interactive()){
#'  library(DESeq2)
#'   # assume dds is a DESeqDataSet already created and DESeq run
#'   res <- results(dds)
#'   opt <- list(out = tempdir())
#'   df <- defWriteMatAsDF(res, opt, "DE_results.tsv")
#'  }
#' }
#' @rdname defWriteMatAsDF
#' @export
defWriteMatAsDF <- function(mat, opt, name){
  fname <-paste(opt$out, name, sep="/", collapse="/")
  df <- mat %>% as.data.frame(row.names = rownames(.)) %>%
    rownames_to_column("gene")
  write_tsv(df, fname)
  return(df)
}
