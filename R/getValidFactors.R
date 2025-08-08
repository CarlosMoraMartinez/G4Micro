#' @title Identify Valid Factor Variables in a Data Frame
#' @description
#' Scans a data frame and determines which variables qualify as valid factors based on type,
#' missing values, class balance, and optional exclusion criteria.
#'
#' @param df A `data.frame` containing the variables to evaluate.
#' @param max_nas Maximum allowed proportion of missing values (including "-" entries)
#'   for a variable to be considered valid. Default is `0.2` (20\%).
#' @param min_pct Minimum proportion of the least frequent class for a factor to be considered balanced.
#'   Default is `0.1` (10\%).
#' @param max_classes Maximum number of unique values for a numeric variable to be considered categorical.
#'   Default is `5`.
#' @param exclude A character vector of variable names to exclude from consideration. Default is `c()`.
#'
#' @return
#' A `data.frame` with the following columns:
#' \itemize{
#'   \item `var` – variable name
#'   \item `is_factor` – whether the variable is a factor, character, or categorical numeric
#'   \item `meet_nas` – whether the variable meets the `max_nas` criterion
#'   \item `is_balanced` – whether the variable meets the `min_pct` balance criterion
#'   \item `not_excluded` – whether the variable is not in the `exclude` list
#'   \item `is_valid` – overall validity based on all above criteria
#' }
#'
#' @details
#' A variable is considered a valid factor if:
#' \enumerate{
#'   \item It is a factor, a character variable, or a numeric variable with
#'         no more than `max_classes` unique values.
#'   \item The proportion of missing values (including "-") does not exceed `max_nas`.
#'   \item It has more than one class and the least frequent class occurs in at least `min_pct` of cases.
#'   \item It is not in the exclusion list.
#' }
#'
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   group = rep(c("A", "B"), each = 5),
#'   score = c(1, 2, 2, 1, 1, 3, 3, 2, 1, 2),
#'   missing_var = c(NA, "-", NA, "-", "-", 1, 1, 1, 1, 1)
#' )
#' getValidFactors(df, max_nas = 0.3, min_pct = 0.2)
#' }
#'
#' @seealso
#'   \code{\link[dplyr]{mutate}},
#'   \code{\link[base]{is.factor}},
#'   \code{\link[base]{table}}
#' @importFrom dplyr mutate
#' @export
getValidFactors <- function(df, max_nas=0.2, min_pct=0.1, max_classes=5, exclude=c()){
  is_factor <- function(x){
    (is.factor(x) |
       is.character(x) |
       (is.numeric(x) & length(unique(x)) <= max_classes))
  }
  meet_nas <- function(x){
    sum(is.na(x) | (x== "-"))/length(x) <= max_nas
  }
  is_balanced <- function(x){
    tt <- table(x)
    tt <- tt/sum(tt)
    min(tt) >= min_pct & (length(tt) > 1)
  }


  res <- data.frame(var = names(df),
                    is_factor = sapply(df, is_factor),
                    meet_nas = sapply(df, meet_nas),
                    is_balanced = sapply(df, is_balanced),
                    not_excluded = (!names(df) %in% exclude)
  ) %>% dplyr::mutate(is_valid = is_factor & meet_nas & is_balanced & not_excluded)
  return(res)
}
