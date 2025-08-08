#' @title Filter Taxa by Raw p-value and Fold Change Threshold
#' @description
#' Filters taxa from a differential expression results data frame based on a significance threshold
#' (raw p-value) and a minimum fold change.
#'
#' @param resdf A data frame containing differential abundance or expression results. Must include
#'   columns `pvalue` (adjusted p-value), `log2FoldChangeShrink` (log2 fold change), and `taxon` (taxon names).
#' @param plim Numeric cutoff for p-value to select significant taxa. Default is 0.05.
#' @param fc Numeric minimum fold change threshold (non-log scale). Taxa with absolute log2 fold change
#'   greater or equal to `log2(fc)` are selected. Default is 1 (i.e., any change).
#'
#' @return A character vector of taxa names meeting the p-value and fold change criteria.
#'
#' @details
#' This function extracts taxa that pass both an unadjusted p-value cutoff and a fold change threshold.
#' Fold change is evaluated on a log2 scale, so the input `fc` is converted internally.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   sig_taxa <- filter_taxa_padj(resdf = my_results_df, plim = 0.01, fc = 2)
#' }
#' }
filter_taxa_praw <- function(resdf, plim=0.05, fc=1){
  taxa <- resdf %>%
    dplyr::filter(pvalue <= plim & abs(log2FoldChangeShrink) >= log2(fc) ) %>%
    pull(taxon)
  return(taxa)
}
