
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param resdf PARAM_DESCRIPTION
#' @param met2use PARAM_DESCRIPTION
#' @param df2plot PARAM_DESCRIPTION
#' @param variable PARAM_DESCRIPTION, Default: 'condition'
#' @param opt PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION, Default: 'heatmap.pdf'
#' @param logscale PARAM_DESCRIPTION, Default: FALSE
#' @param ptype PARAM_DESCRIPTION, Default: 'padj'
#' @param w PARAM_DESCRIPTION, Default: 20
#' @param h PARAM_DESCRIPTION, Default: 14
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{filter}}
#' @rdname makeHeatmapFunctional
#' @export 
#' @importFrom dplyr filter
makeHeatmapFunctional <- function(resdf, met2use, df2plot,
                                  variable = "condition",
                                  opt,
                                  name = "heatmap.pdf",
                                  logscale=FALSE,
                                  ptype = "padj", w=20, h=14){
  default_annot_colors <- c("gray80", "dodgerblue3", "firebrick3", "green4", "pink4", "skyblue")
  outname <- paste(opt$out, name, sep="/", collapse="/")
  annot <- as.data.frame(met2use[, variable])
  names(annot) <- c(variable)
  rownames(annot) <- met2use$sampleID

  if(ptype == "padj"){
    taxa <- resdf %>% dplyr::filter(padj <= opt$pval &
                                      abs(log2FoldChange) >= log2(opt$fc) ) %>%
      rownames
  }else{
    taxa <- resdf %>% dplyr::filter(pvalue <= opt$pval &
                                      abs(log2FoldChange) >= log2(opt$fc) ) %>%
      rownames
  }

  if(length(taxa) < .GlobalEnv$opt$num_genes_default){
    #get only first n genes
    taxa <- resdf[order(resdf$pvalue), ]  %>%
      head(.GlobalEnv$opt$num_genes_default) %>%
      rownames
  }
  mat <- df2plot %>%
    filter(Pathway %in% taxa) %>%
    column_to_rownames("Pathway") %>%
    as.matrix
  if(logscale){
    mat <- log(mat + 1)
  }
  mat <- mat %>%
    t %>% scale %>% t

  cond_colors <- default_annot_colors[1:length(unique(annot$Condition))]
  names(cond_colors) <- unique(annot$Condition)
  annot_colors <- list(Condition=cond_colors)

  fontsize_row = 14 - nrow(mat) / 15
  fontsize_col = 14 - ncol(mat) / 15
  hm <- pheatmap(mat, cluster_rows=T,
                 show_rownames=nrow(mat) < 100,
                 annotation_col = annot,
                 fontsize_row = fontsize_row,
                 fontsize_col = fontsize_col,
                 cluster_cols=T, annotcluster_rowsation_col=annot,
                 annotation_colors = annot_colors
  )
  pdf(outname, width = w, height = h)
  print(hm)
  tmp <- dev.off()
  return(hm)
}
