#' @title Convert P-values to Significance Symbols
#' @description Translates numeric p-values into significance symbols for easier interpretation in plots or tables.
#' @param vector A numeric vector of p-values to be converted.
#' @param limits Numeric vector of cutoff thresholds for significance levels, sorted from largest to smallest. Default is c(0.05, 0.01, 0.001).
#' @param legends Character vector of symbols corresponding to significance levels, where the first element is for non-significant values and subsequent elements represent increasing significance. Default is c("", "*", "**", "***").
#' @return A character vector of the same length as \code{vector}, containing significance symbols.
#' @details
#' The function compares each p-value to the provided \code{limits} and assigns the corresponding symbol from \code{legends}. This is useful for annotating statistical test results visually.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   pvals <- c(0.2, 0.04, 0.005, 0.0005)
#'   getSignif2(pvals)
#'   # Output: ""  "*"  "**"  "***"
#' }
#' }
#' @rdname getSignif2
#' @export
getSignif2 <- function(vector, limits=c(0.05, 0.01, 0.001),
                       legends = c("", "*", "**", "***")){
  x <- ifelse(vector < limits[3], legends[4],
              ifelse(vector < limits[2], legends[3],
                     ifelse(vector< limits[1], legends[2], legends[1])
              )
  )
  return(x)
}
