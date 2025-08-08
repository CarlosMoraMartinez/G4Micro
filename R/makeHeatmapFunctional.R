
#' @title Heatmap Visualization of Differential Functional Pathways
#' @description Generates a clustered heatmap of pathway abundances for pathways identified as significantly differential according to provided statistical results. The heatmap is annotated by sample metadata and optionally log-transformed and scaled.
#'
#' @param resdf A data frame of differential analysis results, typically including columns such as `padj`, `pvalue`, and `log2FoldChange`. Row names should correspond to pathway identifiers.
#' @param met2use A metadata data frame containing sample annotations. Must include a column with sample IDs (named `sampleID`) and a grouping variable specified by `variable`.
#' @param df2plot A data frame of pathway abundance values with at least columns `Pathway` (pathway IDs) and sample columns matching those in `met2use$sampleID`.
#' @param variable Character string naming the metadata column to use for heatmap annotation (default: `"condition"`).
#' @param opt A list of options controlling filtering and plotting, expected to contain at least `pval` (p-value cutoff), `fc` (fold-change cutoff), and `out` (output directory path).
#' @param name Character string filename for the output PDF heatmap file (default: `"heatmap.pdf"`).
#' @param logscale Logical flag indicating whether to log-transform the abundance values prior to plotting (default: `FALSE`).
#' @param ptype Character string indicating which p-value column to use for filtering pathways: `"padj"` for adjusted p-values or `"pvalue"` for raw p-values (default: `"padj"`).
#' @param w Numeric width of the output PDF in inches (default: 20).
#' @param h Numeric height of the output PDF in inches (default: 14).
#'
#' @return A pheatmap object representing the heatmap plot. Additionally, saves the heatmap as a PDF file in the directory specified by `opt$out`.
#'
#' @details
#' This function filters the pathways based on statistical significance and fold-change thresholds defined in `opt`. It selects pathways with adjusted or raw p-values below the cutoff and absolute log2 fold-changes above the threshold.
#'
#' If the number of selected pathways is below a default minimum (`.GlobalEnv$opt$num_genes_default`), the function selects the top pathways by p-value to ensure a meaningful heatmap size.
#'
#' The abundance data for the selected pathways is extracted from `df2plot`, optionally log-transformed, and then row-scaled (z-score normalization) across samples.
#'
#' The heatmap is annotated with sample groupings from `met2use` and colored using a default palette.
#'
#' The plot is generated using the `pheatmap` package with clustering of both rows and columns. Font sizes for row and column labels are adjusted dynamically based on heatmap dimensions to improve readability.
#'
#' The heatmap is saved as a PDF file to the specified path.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Assuming 'resdf' is a differential results dataframe,
#'   # 'metadata' contains sample info,
#'   # and 'abundances' contains pathway abundances
#'   opt <- list(pval=0.05, fc=1.5, out="results")
#'   heatmap_plot <- makeHeatmapFunctional(resdf = resdf,
#'                                        met2use = metadata,
#'                                        df2plot = abundances,
#'                                        variable = "Condition",
#'                                        opt = opt,
#'                                        logscale = TRUE)
#' }
#' }
#'
#' @seealso
#' \code{\link[dplyr]{filter}}, \code{\link[pheatmap]{pheatmap}}
#'
#' @rdname makeHeatmapFunctional
#' @export
#' @importFrom dplyr filter
#' @importFrom pheatmap pheatmap
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
