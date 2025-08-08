#' @title Generate a Formatted GT Table from Differential Expression Results
#' @description Creates a styled gt table from a results data frame filtered by specified genes.
#' Applies color gradients to p-value, adjusted p-value, and log2 fold change columns for easy visualization.
#' @param res A data.frame or tibble containing differential expression results. Must include row names representing taxa or genes and columns `pvalue`, `padj`, and `log2FoldChange`.
#' @param genes A character vector of gene or taxa names to filter the results and include in the table.
#' @param name A character string to be used as the table caption.
#' @return A `gt` table object with formatted colors and styling applied.
#' @details
#' The function filters the results to only include rows corresponding to the specified genes,
#' then creates a gt table with conditional coloring:
#' - pvalue and padj are colored from red ("firebrick3") for high values to blue ("dodgerblue2") for low values.
#' - log2FoldChange is colored from blue for negative values to red for positive values.
#' The table options adjust font sizes, background colors, padding, and width for clarity.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   library(DESeq2)
#'   res <- results(dds)  # DESeq2 result object
#'   genes_of_interest <- c("GeneA", "GeneB", "GeneC")
#'   gt_table <- getGTTableFromRes(res, genes_of_interest, "Differential Expression Table")
#'   print(gt_table)
#' }
#' }
#' @rdname getGTTableFromRes
#' @export
#' @import gt
getGTTableFromRes <- function(res, genes = c(), name = "DAA result"){

  gttable <- res %>% as.data.frame() %>%
    rownames_to_column("Taxa")
  if(length(genes) > 0){
    gttable <- gttable %>%
      dplyr::filter(Taxa %in% genes)
  }

  gttable <- gttable %>%
    dplyr::mutate(Taxa = gsub("_", " ", Taxa)) %>%
    dplyr::mutate(Taxa = paste0("*", Taxa, "*")) %>%
    gt() %>%
    tab_caption(name) %>%
    fmt_markdown(columns = "Taxa") %>%
    fmt_number(
      columns = everything(),
      decimals = 3
    ) %>%
    data_color(
      method = "numeric",
      palette = c("firebrick3", "dodgerblue2"),
      columns = `pvalue`
    ) %>%
    data_color(
      method = "numeric",
      palette = c("firebrick3", "dodgerblue2"),
      columns = `padj`
    ) %>%
    data_color(
      method = "numeric",
      palette = c("dodgerblue2", "firebrick3"),
      columns = `log2FoldChange`
    ) %>%
    tab_options(
      table.background.color = "white",
      column_labels.background.color = "lightgray",
      column_labels.font.size = px(16),
      table.font.size = px(12),
      data_row.padding = px(16),
      table.width = px(600)
    )
  return(gttable)
}
