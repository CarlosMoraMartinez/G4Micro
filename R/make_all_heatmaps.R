#' @title Generate Heatmaps for All Differential Abundance Results
#' @description
#' Iterates over a list of differential abundance results and generates
#' heatmaps for features (taxa/genes) that are significant according to
#' raw p-value and adjusted p-value thresholds. Two heatmaps are generated
#' per result: one using raw p-values and another using adjusted p-values.
#'
#' @param dearesults A list of differential abundance result objects. Each object
#'   is expected to contain at least:
#'   \itemize{
#'     \item `resdf`: a data frame with columns `taxon`, `pvalue`, and `padj`.
#'     \item `contrast_vec`: a vector specifying the contrast used (first element is the variable name, next two are levels to compare).
#'     \item `nested_dir`: output directory path for saving heatmaps.
#'     \item `name`: identifier for the comparison.
#'   }
#' @param df2plot A data frame with abundance or expression data to plot,
#'   where the first column is `gene` and subsequent columns are sample IDs.
#' @param metadata A data frame with sample metadata, including the variable
#'   used for contrasts and a column `sampleID` to match `df2plot`.
#' @param vars2heatmap A vector of variable names (metadata columns) to annotate
#'   in the heatmap.
#' @param dds A `DESeqDataSet` or similar object containing the full dataset,
#'   used by `makeHeatmap` to retrieve normalized data.
#' @param opt A list of options, where `opt$pval` is the p-value threshold
#'   for selecting features to plot.
#' @param max_hm_h Maximum height of the heatmap (in inches). Default: `16`.
#' @param max_hm_w Max heatmap width. If 0, will be automatically determined depending on number of cols. Default: `7`.
#' @return Invisibly returns `NULL`. Heatmap PDF files are saved to the paths
#'   specified by `nested_dir` and `name` for each result.
#' @details
#' For each element in `dearesults`, this function:
#' \enumerate{
#'   \item Selects features (taxa/genes) with `pvalue < opt$pval` and `padj < opt$pval`.
#'   \item Filters `metadata` to include only the samples in the current contrast.
#'   \item Subsets `df2plot` to those samples.
#'   \item Calls `makeHeatmap()` twice — once for raw p-values, once for adjusted p-values.
#' }
#' Errors in heatmap generation are caught and reported without stopping the loop.
#' Open graphics devices are closed after each iteration.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   make_all_heatmaps(
#'     dearesults = my_dea_list,
#'     df2plot = abundance_table,
#'     metadata = sample_metadata,
#'     vars2heatmap = c("Condition", "Sex"),
#'     dds = my_dds,
#'     opt = list(pval = 0.05)
#'   )
#' }
#' }
#' @seealso
#'  \code{\link[dplyr]{filter}}
#' @rdname make_all_heatmaps
#' @export
#' @importFrom dplyr filter select pull
make_all_heatmaps<- function(dearesults, df2plot, metadata, vars2heatmap, dds, opt, max_hm_h=16, max_hm_w=7){
  for(singleres in dearesults){
    taxalist_praw <-singleres$resdf %>% dplyr::filter(pvalue < opt$pval) %>% pull(taxon) %>% unlist %>% unique
    taxalist_padj <-  singleres$resdf %>% dplyr::filter(padj < opt$pval) %>% pull(taxon) %>% unlist %>% unique
    if(singleres$is_numeric){
      samples <- metadata %>% dplyr::filter(!is.na(!!sym(singleres$contrast_name))) %>%
                                              pull(sampleID)
    }else{
      samples <- metadata %>% dplyr::filter(!!sym(singleres$contrast_vec[1]) %in% singleres$contrast_vec[2:3] ) %>%
        pull(sampleID)
    }
    df2plot2 <- df2plot %>% select(gene, all_of(samples))

    tryCatch(makeHeatmap(singleres$resdf, dds, df2plot2, vars2heatmap,
                         opt, name = paste0(singleres$nested_dir, singleres$name, "diff_ab_heatmap_rawpval.pdf"),
                         logscale = F, ptype="pvalue", trim_values = TRUE, taxalist=taxalist_praw, max_hm_h, max_hm_w),
             error=\(x) cat("Error makeHeatmap praw"))
    tryCatch(makeHeatmap(singleres$resdf, dds, df2plot2, vars2heatmap,
                         opt, name = paste0(singleres$nested_dir, singleres$name, "diff_ab_heatmap_adjpval.pdf"),
                         logscale = F, ptype="padj", trim_values = TRUE, taxalist=taxalist_padj, max_hm_h, max_hm_w),
             error=\(x) cat("Error makeHeatmap padj"))

    while(dev.cur() != 1) dev.off()
  }

}
