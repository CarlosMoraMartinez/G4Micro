#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param resdf PARAM_DESCRIPTION
#' @param dds PARAM_DESCRIPTION
#' @param df2plot PARAM_DESCRIPTION
#' @param variable PARAM_DESCRIPTION, Default: 'condition'
#' @param opt PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION, Default: 'heatmap.pdf'
#' @param logscale PARAM_DESCRIPTION, Default: FALSE
#' @param ptype PARAM_DESCRIPTION, Default: 'padj'
#' @param w PARAM_DESCRIPTION, Default: 5
#' @param h PARAM_DESCRIPTION, Default: 4
#' @param trim_values PARAM_DESCRIPTION, Default: FALSE
#' @param italics_rownames PARAM_DESCRIPTION, Default: TRUE
#' @param taxalist PARAM_DESCRIPTION, Default: c()
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
#' @rdname makeHeatmap
#' @export 
#' @importFrom dplyr filter
makeHeatmap <- function(resdf, dds, df2plot,
                        variable = "condition",
                        opt,
                        name = "heatmap.pdf",
                        logscale=FALSE,
                        ptype = "padj", w=5, h=4,
                        trim_values=FALSE,
                        italics_rownames=TRUE, taxalist=c()){
  outname <- paste(opt$out, name, sep="/", collapse="/")
  annot <- as.data.frame(colData(dds)[variable])
  names(annot) <- c(variable)
  rownames(annot) <- colData(dds)$sampleID

  if(length(taxalist)==0){
    if(ptype == "padj"){
      taxa <- resdf %>% filter(padj <= opt$pval &
                                 !is.na(padj) &
                                 abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>%
        pull(taxon)
    }else{
      taxa <- resdf %>% filter(pvalue <= opt$pval &
                                 abs(log2FoldChangeShrink) >= log2(opt$fc) ) %>%
        pull(taxon)
    }
    if(length(taxa) < .GlobalEnv$opt$num_genes_default){
      #get only first n genes
      taxa <- resdf[order(resdf$pvalue), ]  %>%
        head(.GlobalEnv$opt$num_genes_default) %>%
        pull(taxon)
    }
  }else{
    taxa <- resdf$taxon[resdf$taxon %in% taxalist]
  }

  mat <- df2plot %>%
    dplyr::filter(gene %in% taxa) %>%
    column_to_rownames("gene") %>%
    as.matrix %>%
    t %>% scale %>% t
  if(logscale){
    mat <- log(mat + 1)
  }

  fontsize_row = 10 - nrow(mat) / 15

  mat_trim <- mat
  if(trim_values){
    mat_trim[mat> quantile(mat, 0.99)] <- quantile(mat, 0.99)
    mat_trim[mat < quantile(mat, 0.01)] <- quantile(mat, 0.01)
  }

  labels_row <- gsub("_", " ", rownames(mat_trim))
  if(italics_rownames){
    labels_row <- lapply(labels_row, \(x){
      bquote(italic(.(x)))
    }) %>% as.expression()
  }
  hm <- pheatmap(mat_trim, cluster_rows=T,
                 show_rownames=nrow(mat) < 120,
                 annotation_col = annot,
                 fontsize_row = fontsize_row,
                 border_color = NA,
                 labels_row = labels_row,
                 cluster_cols=T,
                 annotcluster_rowsation_col=annot
  )
  w <- w+0.05*ncol(mat)
  h <- h+0.05*nrow(mat)
  pdf(outname, width = w, height = h)
  print(hm)
  tmp <- dev.off()
}
