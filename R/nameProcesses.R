#' @title Assign Descriptive Process Names to Table Rows
#' @description
#' This function adds a descriptive process name column to a data table based on a mapping table.
#' It matches row names of the input table (`tab`) to short process codes in `namedtab` and
#' appends the corresponding long descriptive names as a new column.
#' If the `namedtab` is `NULL`, the original table is returned unchanged.
#'
#' @param tab A data frame or matrix whose row names correspond to short process identifiers.
#' @param namedtab A data frame with at least two columns: `short` containing short process codes,
#' and `long` containing the corresponding descriptive process names.
#'
#' @return The original data frame `tab` with an additional column `process_name` containing
#' descriptive names matched by row names from `namedtab`. Rows without matches will have `NA`.
#'
#' @details
#' The function performs a vectorized match of `rownames(tab)` against the `short` column of `namedtab`.
#' It then adds the `long` names from `namedtab` to a new column `process_name` in `tab`.
#' This is useful to enhance tables with human-readable labels for processes or features.
#' If `namedtab` is `NULL`, the function simply returns the input table as is.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   tab <- data.frame(value = c(10, 15), row.names = c("p1", "p2"))
#'   namedtab <- data.frame(short = c("p1", "p2"), long = c("Process A", "Process B"))
#'   tab_named <- nameProcesses(tab, namedtab)
#'   print(tab_named)
#' }
#' }
#'
#' @rdname nameProcesses
#' @export
nameProcesses<- function(tab, namedtab){
  if(is.null(namedtab)){return(tab)}
  tab$process_name <- namedtab$long[match(rownames(tab), namedtab$short)]
  return(tab)
}
