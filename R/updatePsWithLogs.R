#' @title Add Log-Transformed Metadata Variables to a Phyloseq Object
#'
#' @description
#' Adds log-transformed versions of specified numeric metadata variables
#' in a \code{\link[phyloseq]{phyloseq}} object. Each transformed variable is
#' named by appending \code{"_log"} to the original variable name and is
#' computed as \eqn{\log(x + 1)} to handle zero values safely.
#'
#' @param phobj A \code{\link[phyloseq]{phyloseq}} object containing
#'   sample metadata.
#' @param vars Character vector of metadata variable names to log-transform.
#'   Default is \code{c("nreads")}.
#'
#' @return A \code{\link[phyloseq]{phyloseq}} object with additional
#'   log-transformed variables added to its \code{\link[phyloseq]{sample_data}}.
#'
#' @details
#' This transformation is commonly used to stabilize variance and
#' normalize skewed distributions, especially for read counts and other
#' strictly non-negative metadata variables. The addition of 1 before
#' taking the logarithm prevents issues with zero values.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Add log-transformed read counts to a phyloseq object
#'   ps <- updatePsWithLogs(phobj, vars = c("nreads", "total_abundance"))
#' }
#' }
#' @rdname updatePsWithLogs
#' @export
updatePsWithLogs <- function(phobj, vars = c("nreads")){

  for(v in vars){
    newname <- paste0(v, "_log")
    sample_data(phobj)[[newname]] <- log(sample_data(phobj)[[v]] +1  )
  }
  return(phobj)
}
