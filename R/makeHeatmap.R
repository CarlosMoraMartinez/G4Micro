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
#' @param check_taxa Check whether all taxa in \code{taxalist} exist or not.
#' @param max_hm_h Heatmap height. If 0, will be automatically determined depending on number of rows.
#' @param max_hm_w Max heatmap width. If 0, will be automatically determined depending on number of cols
#' @param annotcols_num Max heatmap height
#' @param annotpal_cat Name of Wes Anderson palette to annotate categorical variables. Default = \code{c("blue", "white", "red")}
#' @return None. The function saves a PDF heatmap file to the specified output directory. Default = "Darjeeling1"
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
#' @importFrom wesanderson wes_palette
#' @importFrom SummarizedExperiment colData
makeHeatmap <- function(resdf, dds, df2plot,
                        variable = "condition",
                        opt,
                        name = "heatmap.pdf",
                        logscale=FALSE,
                        ptype = "padj", w=5, h=4,
                        trim_values=FALSE,
                        italics_rownames=TRUE, taxalist=c(),
                        check_taxa = TRUE,
                        max_hm_h=16,
                        max_hm_w=7,
                        annotcols_num = c("blue", "white", "red"),
                        annotpal_cat = "Darjeeling1"){
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
    if(check_taxa){
      taxa <- resdf$taxon[resdf$taxon %in% taxalist]
    }else{
      taxa <- taxalist
    }
  }

  mat <- df2plot %>%
    dplyr::filter(gene %in% taxa) %>%
    column_to_rownames("gene") %>%
    as.matrix %>%
    t %>% scale %>% t
  if(logscale){
    mat <- log(mat + 1)
  }

  #fontsize_row = 10 - nrow(mat) / 15
  fontsize_row = 10 - nrow(mat) / 13
  fontsize_col = 10 - ncol(mat) / 13

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

  numcols <- max(sapply(annot %>% select_if(\(x) !is.numeric(x)), \(x)length(unique(x))))
  cat("NUM COLS: ", numcols, "\n")
  #cc <- ggsci::pal_npg()(numcols) #palette = "category10"
  cc <- grDevices::colorRampPalette( wesanderson::wes_palette(annotpal_cat))(numcols)

  user.colfn=colorRampPalette(cc)
  newcc <- user.colfn(numcols) # in case there are too many colors

  color_list_cols <- lapply(annot %>% select_if(\(x) !is.numeric(x)),
                            \(x) {y <-newcc[1:length(unique(x))]; names(y)<- unique(x); y})

  # Color list for continuous vars
  cont_colors <- lapply(annot %>% select_if(\(x)is.numeric(x)), function(x) {
    colorRampPalette(annotcols_num)(100)
  })
  names(cont_colors) <- names(annot %>% select_if(\(x) is.numeric(x)))
  color_list_cols <- append(color_list_cols, cont_colors)
  if(! all(rownames(mat_trim) %in% colnames(mat_trim))){
    hm <- pheatmap(mat_trim,
                   show_rownames=nrow(mat) < 120,
                   annotation_col = annot,
                   annotation_colors = color_list_cols,
                   fontsize_col = fontsize_col,
                   fontsize_row = fontsize_row,
                   #annotcluster_rowsation_col=annot,
                   border_color = NA,
                   labels_row = labels_row,
                   cluster_cols=T,
                   cluster_rows=T
    )
  }else{
    # Correlation heatmap
    hm <- pheatmap(mat_trim,
                   show_rownames=nrow(mat) < 120,
                   annotation_col = annot,
                   annotation_row = annot,
                   fontsize_row = fontsize_row,
                   fontsize_col = fontsize_row,
                   border_color = NA,
                   labels_row = labels_row,
                   cluster_cols=T,
                   cluster_rows=T
    )
  }
  #w <- w+0.05*ncol(mat)
  #h <- h+0.05*nrow(mat)

  if(max_hm_h == 0){
    h <- if(nrow(mat)>10) h+0.05*nrow(mat) else 7
  }else{
    h <- max_hm_h
  }
  if(max_hm_w == 0){
    w <- if(ncol(mat)>10) w+0.05*ncol(mat) else 7
  }else{
    w <- max_hm_w
  }
  pdf(outname, width = w, height = h)
  print(hm)
  tmp <- dev.off()
}
