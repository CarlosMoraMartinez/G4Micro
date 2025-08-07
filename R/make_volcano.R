#' @title Create and Save a Volcano Plot from Differential Expression Results
#' @description
#' Generates a volcano plot using the \code{EnhancedVolcano} package based on
#' differential expression results, and saves the plot as a PDF file.
#'
#' @param res A data frame or matrix containing differential expression results,
#'   including \code{log2FoldChange} and p-value columns.
#' @param opt A list or object containing analysis options, including:
#'   \itemize{
#'     \item \code{pval}: numeric, significance threshold for p-values.
#'     \item \code{fc}: numeric, fold-change cutoff.
#'     \item \code{out}: character, output directory to save the plot PDF.
#'   }
#' @param name Character. Filename for the output PDF volcano plot.
#' @param pcol Character. Name of the p-value column to use for significance (default \code{"pvalue"}).
#'
#' @return The \code{EnhancedVolcano} plot object (invisibly returned after saving the PDF).
#'
#' @details
#' The function attempts to create a volcano plot with labeled points from the row names
#' of \code{res}. It highlights significant genes or taxa based on p-value and fold-change
#' thresholds provided in \code{opt}. The plot is saved to a PDF file in the specified
#' output directory.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   opt <- list(pval = 0.05, fc = 2, out = "results")
#'   make_volcano(res = my_results, opt = opt, name = "volcano_plot.pdf")
#' }
#' }
#' @rdname make_volcano
#' @export
#' @importFrom EnhancedVolcano EnhancedVolcano
make_volcano <- function(res, opt, name, pcol="pvalue"){
  outname <- paste(opt$out, name, sep="/", collapse="/")
  pdf(outname, width = 12, height = 8)
  volc <- NULL
  try({
    volc <- EnhancedVolcano(res,
                            lab = rownames(res),
                            x = 'log2FoldChange', title = "", subtitle="",
                            y = pcol,
                            pCutoff = opt$pval,
                            pCutoffCol = pcol,
                            FCcutoff = log2(opt$fc),
                            legendPosition = 'right',
                            legendLabSize = 12,
                            legendIconSize = 3.0,
                            legendLabels=c('Not sig.',
                                           paste("FC > ", as.character(opt$fc), sep="", collapse=""),
                                           paste(pcol, " < ", as.character(opt$pval), sep="", collapse=""),
                                           paste("FC > ", as.character(opt$fc), " & ", pcol, " < ", as.character(opt$pval), sep="", collapse=""))
    )
  })
  print(volc)
  tmp <- dev.off()
  volc
}
