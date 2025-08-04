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
