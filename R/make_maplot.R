#' @title Generate MA-Plot from DESeq2 Results and Save as PDF
#' @description
#' Creates an MA-plot from DESeq2 differential expression results using specified
#' p-value and fold-change cutoffs, saves the plot as a PDF file, and also displays it.
#' Horizontal lines indicate the fold-change thresholds.
#'
#' @param res A \code{DESeq2} results object containing differential expression statistics.
#' @param opt A list or object with plotting options, including:
#'   \itemize{
#'     \item \code{pval}: numeric, significance threshold for adjusted p-value (alpha).
#'     \item \code{fc}: numeric, fold-change cutoff (linear scale, e.g., 2 for 2-fold).
#'     \item \code{out}: character, output directory path for saving the PDF plot.
#'   }
#' @param name Character string specifying the filename of the saved plot PDF. Default is \code{"maplot.pdf"}.
#'
#' @return None. The function invisibly returns \code{NULL} after saving the plot and displaying it.
#'
#' @details
#' The MA-plot visualizes log2 fold changes versus mean normalized counts with
#' points colored based on adjusted p-values. Fold-change thresholds are marked
#' by horizontal dashed blue lines at \code{±log2(fc)}.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   opt <- list(pval = 0.05, fc = 2, out = "results")
#'   make_maplot(res, opt, name = "my_ma_plot.pdf")
#' }
#' }
#' @seealso
#'  \code{\link[DESeq2]{plotMA}}
#' @rdname make_maplot
#' @export
#' @importFrom DESeq2 plotMA
make_maplot<- function(res, opt, name="maplot.pdf"){
  outname <- paste(opt$out, name, sep="/", collapse="/")
  pdf(outname)
  ma <- DESeq2::plotMA(res, ylim=c(-3, 3), alpha=opt$pval)
  abline(h = c(-1*log2(opt$fc), log2(opt$fc)), col = "dodgerblue3", lwd =2, lty =2)
  tmp <- dev.off()
  DESeq2::plotMA(res, ylim=c(-3, 3), alpha=opt$pval)
  abline(h = c(-1*log2(opt$fc), log2(opt$fc)), col = "dodgerblue3", lwd =2, lty =2)

}
