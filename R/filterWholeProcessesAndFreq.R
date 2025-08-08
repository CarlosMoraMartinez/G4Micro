#' @title Filter Functional Processes by Count and Frequency
#' @description
#' Filters a functional matrix data frame to retain only whole (non-combined) processes,
#' excludes unwanted categories, and retains rows with sufficient counts across samples.
#'
#' @param df A data frame or tibble containing functional pathways as rows
#'   and samples as columns. The first column should be named "Pathway".
#' @param opt A list containing filtering thresholds:
#'   \itemize{
#'     \item \code{mincount}: Minimum count threshold to consider a value as present.
#'     \item \code{minsampleswithcount}: Minimum number of samples that must have counts above \code{mincount}.
#'   }
#'
#' @return A filtered data frame retaining only rows meeting the criteria:
#'   - The \code{Pathway} does not contain the "|" character (i.e., whole processes).
#'   - The \code{Pathway} does not contain "UNINTEGRATED", "UNMAPPED", or "UNGROUPED".
#'   - Rows where the count is above \code{mincount} in at least \code{minsampleswithcount} samples.
#'
#' @details
#' This function removes aggregated or ambiguous pathways and filters rows based on count thresholds.
#' It assumes that \code{df} has pathways in the first column named "Pathway" and numeric counts in other columns.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   opt <- list(mincount = 5, minsampleswithcount = 3)
#'   df_filtered <- filterWholeProcessesAndFreq(df, opt)
#' }
#' }
#'
#' @seealso
#' \code{\link[dplyr]{filter}}
#'
#' @rdname filterWholeProcessesAndFreq
#' @export
#' @importFrom dplyr filter
filterWholeProcessesAndFreq <- function(df, opt){
  df2 <- df %>% dplyr::filter(!grepl("\\|", Pathway)) %>%
    dplyr::filter(!grepl("UNINTEGRATED|UNMAPPED|UNGROUPED", Pathway, perl=T))
  keep <- rowSums(df2 > opt$mincount) > opt$minsampleswithcount
  df2 <- df2[keep, ]
  return(df2)
}
