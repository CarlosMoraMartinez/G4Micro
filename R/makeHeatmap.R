#' @title Generate Heatmap of Selected Taxa or Genes from Differential Expression Results
#' @description
#' Creates and saves a heatmap of expression values for selected taxa or genes based on
#' significance criteria from differential expression analysis results. The heatmap
#' includes sample annotations and offers options for log-scaling, trimming extreme values,
#' and italicizing row names.
#'
#' @param resdf A data frame containing differential expression results, including
#'   p-values, adjusted p-values, log2 fold changes, and taxa or gene identifiers.
#' @param dds A \code{DESeqDataSet} object containing expression data and sample metadata.
#' @param df2plot A data frame or matrix of expression values with genes/taxa as rows
#'   and samples as columns. Must include a column named "gene" for filtering.
#' @param variable Character. The name of the sample metadata column to use for annotation
#'   (default is \code{"condition"}).
#' @param opt A list or object containing analysis options, including:
#'   \itemize{
#'     \item \code{pval}: numeric, significance threshold for p-value or adjusted p-value.
#'     \item \code{fc}: numeric, fold-change cutoff.
#'     \item \code{out}: character, output directory to save the heatmap PDF.
#'   }
#' @param name Character. Filename for the output PDF heatmap (default \code{"heatmap.pdf"}).
#' @param logscale Logical. Whether to log-transform the expression values before plotting (default \code{FALSE}).
#' @param ptype Character. Type of p-value to filter by: \code{"padj"} for adjusted p-values or \code{"pvalue"} (default \code{"padj"}).
#' @param w Numeric. Width of the output PDF in inches (default 5).
#' @param h Numeric. Height of the output PDF in inches (default 4).
#' @param trim_values Logical. Whether to trim expression values at the 1st and 99th percentiles to reduce outliers (default \code{FALSE}).
#' @param italics_rownames Logical. Whether to italicize the row labels (default \code{TRUE}).
#' @param taxalist Character vector. A custom list of taxa or genes to include; if empty, selection is based on significance criteria (default empty vector).
#'
#' @return None. The function saves a PDF heatmap file to the specified output directory.
#'
#' @details
#' The function filters taxa/genes based on p-value or adjusted p-value thresholds and
#' fold-change cutoffs. If the number of selected taxa is fewer than a default number
#' (from \code{opt$num_genes_default}), it selects the top taxa by p-value.
#' Expression values are optionally log-transformed and scaled per gene before plotting.
#' The heatmap is generated with sample annotations for the specified metadata variable.
#' Row labels can be italicized for better visualization of taxon names.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   opt <- list(pval = 0.05, fc = 2, out = "results", num_genes_default = 50)
#'   makeHeatmap(resdf, dds, expression_df,
#'               variable = "condition", opt = opt,
#'               name = "diff_heatmap.pdf", logscale = TRUE)
#' }
#' }
#'
#' @seealso
#' \code{\link[dplyr]{filter}}, \code{\link[pheatmap]{pheatmap}}
#'
#' @rdname makeHeatmap
#' @export
#' @importFrom pheatmap pheatmap
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
    if(length(taxa) < opt$num_genes_default){
      #get only first n genes
      taxa <- resdf[order(resdf$pvalue), ]  %>%
        head(opt$num_genes_default) %>%
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
