#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param res PARAM_DESCRIPTION
#' @param genes PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname getGTTableFromRes
#' @export 
getGTTableFromRes <- function(res, genes, name){
  gttable <- res %>% as.data.frame() %>%
    rownames_to_column("Taxa") %>%
    filter(Taxa %in% genes) %>%
    gt() %>%
    tab_caption(name) %>%
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
