#' @title Generate Summary Tables from DESeq2 Results
#' @description
#' Summarizes differential expression results by categorizing taxa based on raw and adjusted p-values
#' and log2 fold change thresholds, and outputs contingency tables showing counts in each category.
#' The tables separate taxa into more frequent, less frequent, or equal abundance based on the
#' specified fold-change cutoff, for both raw p-values and adjusted p-values.
#'
#' @param res A \code{DESeq2} results object or data frame containing differential expression results,
#'   including columns \code{pvalue}, \code{padj}, and \code{log2FoldChange}.
#' @param opt A list or object containing analysis options:
#'   - \code{pval}: numeric, p-value cutoff threshold.
#'   - \code{fc}: numeric, fold-change cutoff (linear scale, e.g., 2 for 2-fold).
#'   - \code{out}: character, output directory where summary tables will be saved.
#'
#' @return A list with two contingency tables:
#'   \describe{
#'     \item{\code{restab}}{Table of counts by categories using raw p-values and fold changes.}
#'     \item{\code{restab_adj}}{Table of counts by categories using adjusted p-values and fold changes.}
#'   }
#'   Both tables show counts of taxa split by significance (p-value threshold) and direction/magnitude of fold change.
#'
#' @details
#' This function processes DESeq2 differential abundance results to create summary tables
#' that help interpret the number of taxa meeting significance and fold-change criteria.
#' It writes the contingency tables as TSV files in the specified output directory.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Assuming 'res' is a DESeq2 results object and 'opt' is a list of options:
#'   opt <- list(pval = 0.05, fc = 2, out = "./results/")
#'   summary_tables <- getSummaryTablesDeseq(res, opt)
#'   print(summary_tables$restab)
#'   print(summary_tables$restab_adj)
#' }
#' }
#' @seealso
#'  \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{select}}
#' @rdname getSummaryTablesDeseq
#' @export
#' @importFrom dplyr mutate select
getSummaryTablesDeseq <- function(res, opt){
  namestab <- c(paste("p < ", as.character(opt$pval), sep="", collapse="") ,
                paste("LFC > ", as.character(log2(opt$fc)), sep="", collapse="")
  )
  restab <- res %>% as.data.frame() %>%
    dplyr::mutate(a = pvalue <= opt$pval,
                  b = ifelse(log2FoldChange <= -log2(opt$fc),"less frequent",
                             ifelse( log2FoldChange >= log2(opt$fc), "more frequent", "equal"))
    ) %>%
    dplyr::select(a:b) %>%
    set_names(namestab) %>%
    table()
  rownames(restab) <- ifelse(rownames(restab) == "TRUE",
                             paste("p <= ", as.character(opt$pval), sep="", collapse=""),
                             paste("p > ", as.character(opt$pval), sep="", collapse=""))

  #restab %>% kable(caption="Number of taxons(species) per category using raw p-values")

  restab_adj <- res %>% as.data.frame() %>%
    dplyr::mutate(a = padj <= opt$pval,
                  b = ifelse(log2FoldChange <= -log2(opt$fc),"less frequent",
                             ifelse( log2FoldChange >= log2(opt$fc), "more frequent", "equal"))
    ) %>%
    dplyr::select(a:b) %>%
    set_names(namestab) %>%
    table()
  rownames(restab_adj) <- ifelse(rownames(restab_adj) == "TRUE",
                                 paste("p <= ", as.character(opt$pval), sep="", collapse=""),
                                 paste("p > ", as.character(opt$pval), sep="", collapse=""))

  #restab_adj %>% kable(caption="Number of taxons (species) per category using adjusted p-values")
  write_tsv(as.data.frame(restab), paste0(opt$out, "/num_diff_rawpval.tsv"))
  write_tsv(as.data.frame(restab_adj), paste0(opt$out, "/num_diff_adjpval.tsv"))
  return(list("restab"=restab, "restab_adj"=restab_adj))
}
