#' @title Get Significance Stars from p-values
#' @description
#' Converts numeric p-values into significance levels represented as stars:
#' "***" for p < 0.001, "**" for p < 0.01, "*" for p < 0.05, and "NS" otherwise.
#'
#' @param vector A numeric vector of p-values.
#'
#' @return A character vector of the same length as \code{vector}, with significance stars or "NS" for each p-value.
#'
#' @details
#' This function is useful for annotating statistical test results with conventional significance symbols.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   pvals <- c(0.0005, 0.007, 0.03, 0.1)
#'   getSignif(pvals)
#' }
#' }
#'
#' @rdname getSignif
#' @export
getSignif <- function(vector){
  x <- ifelse(vector < 0.001, "***",
              ifelse(vector < 0.01, "**",
                     ifelse(vector< 0.05, "*", "NS")
              )
  )
  return(x)
}
