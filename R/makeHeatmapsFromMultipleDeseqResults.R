#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param deseq_results_list PARAM_DESCRIPTION
#' @param daa_main PARAM_DESCRIPTION
#' @param main_name PARAM_DESCRIPTION
#' @param vars2plot PARAM_DESCRIPTION
#' @param italics_rownames PARAM_DESCRIPTION, Default: T
#' @param pfilt PARAM_DESCRIPTION, Default: 0.05
#' @param pplot PARAM_DESCRIPTION, Default: 0.05
#' @param name PARAM_DESCRIPTION
#' @param outdir PARAM_DESCRIPTION
#' @param w PARAM_DESCRIPTION, Default: 5
#' @param h PARAM_DESCRIPTION, Default: 4
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{select}}
#' @rdname makeHeatmapsFromMultipleDeseqResults
#' @export 
#' @importFrom dplyr mutate select
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
