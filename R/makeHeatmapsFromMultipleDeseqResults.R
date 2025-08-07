#' @title Generate Heatmaps from Multiple DESeq2 Results
#' @description Creates heatmaps visualizing significance levels (-10*log10 p-values and adjusted p-values)
#' across multiple differential abundance analysis (DAA) results, combining results from multiple contrasts.
#' Saves merged DAA tables and heatmaps to specified output directory.
#' @param deseq_results_list A named list of data frames, each containing DESeq2 results for a different variable/contrast.
#' @param daa_main A data frame containing the main DESeq2 results to highlight.
#' @param main_name Character. The name/key in the list to assign to the main DESeq2 results.
#' @param vars2plot Character vector of variable names indicating which DESeq2 results in the list to include in the heatmaps.
#' @param italics_rownames Logical, default TRUE. Whether to display row names (taxa) in italics on the heatmap.
#' @param pfilt Numeric, default 0.05. P-value cutoff for filtering taxa to display in the heatmaps.
#' @param pplot Numeric, default 0.05. P-value cutoff threshold above which p-values are replaced with NA for visualization.
#' @param name Character. Base name used for output files.
#' @param outdir Character. Output directory path where results and heatmaps will be saved.
#' @param w Numeric, default 5. Base width (in inches) of the output heatmap PDF.
#' @param h Numeric, default 4. Base height (in inches) of the output heatmap PDF.
#' @return None. The function saves heatmaps as PDF files and writes a merged DAA results TSV file.
#' @details
#' This function merges multiple DESeq2 results data frames by taxa, transforms p-values for visualization,
#' filters taxa by significance, and generates heatmaps showing p-values and adjusted p-values across contrasts.
#' Row names can be displayed in italics for better readability of taxon names.
#' Heatmaps are saved as PDF files with dynamically adjusted dimensions based on matrix size.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Assume deseq_list is a named list of DESeq2 results,
#'   # daa_main is a main DESeq2 result data.frame,
#'   # vars_to_plot is a vector of contrasts to include
#'   makeHeatmapsFromMultipleDeseqResults(
#'     deseq_results_list = deseq_list,
#'     daa_main = daa_main,
#'     main_name = "main_contrast",
#'     vars2plot = vars_to_plot,
#'     name = "output_heatmaps",
#'     outdir = "results/",
#'     pfilt = 0.05,
#'     pplot = 0.01
#'   )
#' }
#' }
#' @seealso
#'  \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{select}}
#' @rdname makeHeatmapsFromMultipleDeseqResults
#' @export
#' @importFrom pheatmap pheatmap
#' @importFrom dplyr mutate select
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tidyr pivot_wider
#' @importFrom readr write_tsv
#' @importFrom purrr map
#' @importFrom stringr str_replace
#' @importFrom stats setNames
makeHeatmapsFromMultipleDeseqResults <- function(deseq_results_list,
                                                 daa_main,
                                                 main_name,
                                                 vars2plot,
                                                 italics_rownames=T,
                                                 pfilt=0.05, pplot=0.05,
                                                 name, outdir, w=5, h=4){
  library(RColorBrewer)
  deseq_results_list <- deseq_results_list[vars2plot]
  deseq_results_list[[main_name]] <- daa_main
  all_daa_tab <- names(deseq_results_list) %>%
    map(\(vname) {
      deseq_results_list[[vname]] %>%
        dplyr::mutate(variable_controlled=vname)
    }) %>%
    bind_rows() %>%
    dplyr::mutate(variable_controlled = gsub("tratamiento_coagul_betabloq_etc", "tratamiento_coag", .$variable_controlled),
                  pvalue = ifelse(pvalue>pplot, NA, -10*log10(pvalue)), #non-significant p-values will be shown in gray
                  padj=ifelse(padj > pplot, NA, -10*log10(padj))
    ) %>%
    dplyr::select(taxon, variable_controlled, pvalue, padj, log2FoldChange, log2FoldChangeShrink) %>%
    pivot_wider( names_from = variable_controlled,
                 values_from = c(pvalue, padj, log2FoldChange, log2FoldChangeShrink))
  write_tsv(all_daa_tab, paste0(outdir, name, "_MergedCorrectedDAA.tsv"))
  lapply(c("padj", "pvalue"), \(vname){
    mat <- all_daa_tab %>%
      column_to_rownames("taxon") %>%
      dplyr::select(starts_with(vname)) %>%
      as.matrix
    colnames(mat) <- gsub(paste0(vname, "_"), "", colnames(mat))
    mat <- mat[mat[, main_name] >= -10*log10(pfilt) & !is.na(mat[, main_name]) , ]
    mat <- mat[order(mat[, main_name], decreasing = T), ]

    #mat <- mat[ !apply(mat, MAR=1, \(x) all(is.na(x))) ,]
    annot <- data.frame(main_name = mat[, main_name])
    names(annot) <- main_name
    mat <- mat[, colnames(mat) != main_name]
    colnames(mat) <- colnames(mat) %>% gsub("_", " ", .) %>% gsub("tratamiento", "treat.", .)

    labels_row <- gsub("_", " ", rownames(mat))
    if(italics_rownames){
      labels_row <- lapply(labels_row, \(x){
        bquote(italic(.(x)))
      }) %>% as.expression()
    }
    annoCol<-list(Group= colorRampPalette(rev(brewer.pal(n = 7, name =
                                                           "RdYlBu")))(100))#colorRamp(c("dodgerblue3", "firebrick4")))
    names(annoCol) <- main_name
    fontsize_row = 10 - nrow(mat) / 15
    fname <- paste0(outdir, name, '_',main_name, '_', vname, "_.pdf")
    pheatmap(mat, cluster_rows=F,
             show_rownames=nrow(mat) < 120,
             annotation_row = annot,
             fontsize_row = fontsize_row,
             filename = fname,
             width =  w+0.05*ncol(mat),
             height =  h+0.05*nrow(mat),
             color =  colorRampPalette(rev(brewer.pal(n = 7, name =
                                                        "RdYlBu")))(100),
             #legend_breaks = c(-10*log10(c(0.05, 0.001, 0.000001, 0.00000001)), max(mat)),
             #legend_names = as.character(c(0.05, 0.001, 0.000001, 0.00000001, max(mat[!is.na(mat)]))),
             #border_color = NA,
             labels_row = labels_row,
             na_col = "#DDDDDD",
             cluster_cols=T,
             annotation_colors = annoCol
             #annotcluster_rowsation_col=annot
    )

  })

}
