#' @title Apply String Replacements and Title-Case to Data Frame Columns
#' @description Iteratively applies string replacements (via \code{gsub}) to specified columns in a data frame,
#'              and optionally converts selected columns to title case.
#'
#' @param df A data frame to be modified.
#' @param replace_list A list of character vectors, each of length 3, where:
#'   \itemize{
#'     \item \code{[1]} = pattern (regular expression to match)
#'     \item \code{[2]} = replacement string
#'     \item \code{[3]} = column name in \code{df} where replacement should be applied
#'   }
#' @param totitlecase A character vector of column names in \code{df} that should be converted
#'                    to title case using \code{\link[tools]{toTitleCase}}. Default: \code{c()}.
#'
#' @return A modified data frame where the specified replacements have been applied, and selected
#'         columns converted to title case.
#'
#' @details The function uses \code{\link[purrr]{reduce}} to sequentially apply all replacements
#'          in \code{replace_list} to the corresponding columns in \code{df}.
#'          Regular expressions are supported for pattern matching in \code{gsub} with \code{perl = TRUE}.
#'          After replacements, the function optionally applies \code{tools::toTitleCase} to the
#'          columns listed in \code{totitlecase}.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   library(dplyr)
#'   df <- data.frame(Name = c("john doe", "jane-smith"),
#'                    City = c("new_york", "los-angeles"))
#'
#'   replacements <- list(
#'     c("_", " ", "City"),       # replace underscores with spaces in City
#'     c("-", " ", "City")        # replace dashes with spaces in City
#'   )
#'
#'   df_mod <- makeReplacementsDF(df, replacements, totitlecase = c("Name", "City"))
#'   print(df_mod)
#' }
#' }
#'
#' @seealso
#'  \code{\link[purrr]{reduce}},
#'  \code{\link[dplyr]{mutate}},
#'  \code{\link[tools]{toTitleCase}}
#' @rdname makeReplacementsDF
#' @export
#' @importFrom purrr reduce
#' @importFrom dplyr mutate
#' @importFrom tools toTitleCase
makeReplacementsDF <- function(df, replace_list, totitlecase=c()){
  dfmod <- purrr::reduce(replace_list, ~ .x %>%
                           dplyr::mutate(!!sym(.y[3]) := gsub(.y[1], .y[2], !!sym(.y[3]), perl=TRUE)),
                         .init = df) %>%
    dplyr::mutate(across(all_of(totitlecase), tools::toTitleCase))

  return(dfmod)
}
