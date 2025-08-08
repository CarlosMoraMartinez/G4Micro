#' @title Identify Valid Numeric Variables in a Data Frame
#' @description
#' Scans a data frame and determines which variables qualify as valid numeric variables
#' based on type, number of unique values, and missing value thresholds.
#'
#' @param df A `data.frame` containing the variables to evaluate.
#' @param max_nas Maximum allowed proportion of missing values (including "-" entries)
#'   for a variable to be considered valid. Default is `0` (no missing values allowed).
#' @param max_classes Minimum number of unique values a numeric variable must have
#'   to be considered valid. Default is `4`.
#'
#' @return
#' A `data.frame` with the following columns:
#' \itemize{
#'   \item `var` – variable name
#'   \item `is_number` – whether the variable is numeric and has at least `max_classes` unique values
#'   \item `meet_nas` – whether the variable meets the `max_nas` criterion
#'   \item `is_valid` – overall validity based on both criteria
#' }
#'
#' @details
#' A variable is considered a valid numeric if:
#' \enumerate{
#'   \item It is of numeric type.
#'   \item It contains at least `max_classes` unique values.
#'   \item The proportion of missing values (including "-") does not exceed `max_nas`.
#' }
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   age = c(25, 30, 35, 40, NA),
#'   score = c(1, 2, 2, 1, 1),
#'   text = c("A", "B", "C", "A", "B")
#' )
#' getValidNumeric(df, max_nas = 0.2, max_classes = 3)
#' }
#'
#' @seealso
#'   \code{\link[dplyr]{mutate}},
#'   \code{\link[base]{is.numeric}},
#'   \code{\link[base]{unique}}
#' @importFrom dplyr mutate
#' @export
getValidNumeric <- function(df, max_nas=0.0, max_classes=4){
  is_number <- function(x){
    (is.numeric(x) &
       (length(unique(x)) >= max_classes))
  }
  meet_nas <- function(x){
    sum(is.na(x) | (x== "-"))/length(x) <= max_nas
  }

  res <- data.frame(var = names(df),
                    is_number = sapply(df, is_number),
                    meet_nas = sapply(df, meet_nas)
  ) %>% dplyr::mutate(is_valid = is_number & meet_nas)
  return(res)
}
