
#' @title Filter Species-Level Processes by Count and Frequency
#' @description
#' Filters a functional matrix data frame to retain only species-level processes
#' (identified by the presence of "|" in the pathway name), excludes unwanted categories,
#' and retains rows with sufficient counts across samples.
#'
#' @param df A data frame or tibble with functional pathways as rows (first column named "Pathway")
#'   and numeric counts or abundances in the subsequent columns.
#' @param opt A list containing filtering parameters:
#'   \itemize{
#'     \item \code{mincount}: Minimum count threshold to consider a value as present.
#'     \item \code{minsampleswithcount}: Minimum number of samples that must have counts above \code{mincount}.
#'   }
#'
#' @return A filtered data frame retaining only rows where:
#'   - \code{Pathway} contains "|" indicating species-level resolution.
#'   - \code{Pathway} does not contain "UNINTEGRATED", "UNMAPPED", or "UNGROUPED".
#'   - The count is above \code{mincount} in at least \code{minsampleswithcount} samples.
#'
#' @details
#' This function focuses on species-level pathway abundances, filtering out ambiguous or unassigned entries
#' and ensuring that pathways are sufficiently represented across samples.
#' It assumes the first column of \code{df} is named "Pathway" and contains pathway identifiers.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{filter}}
#' @rdname filterBySpeciesAndFreq
#' @export
#' @importFrom dplyr filter
filterBySpeciesAndFreq <- function(df, opt){
  df2 <- df %>% dplyr::filter(grepl("\\|", Pathway))%>%
    dplyr::filter(!grepl("UNINTEGRATED|UNMAPPED|UNGROUPED", Pathway, perl=T))
  keep <- rowSums(df2 > opt$mincount) > opt$minsampleswithcount
  return(df2[keep, ])
}
