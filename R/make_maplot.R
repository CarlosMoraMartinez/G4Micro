#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param res PARAM_DESCRIPTION
#' @param opt PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION, Default: 'maplot.pdf'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
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
