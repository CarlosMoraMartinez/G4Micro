#' @title Map CAZy Families to CAZy Classes
#' @description
#' Given a character vector of CAZy family identifiers separated by pipes ("|"),
#' returns the corresponding CAZy class names based on the CAZy classification system.
#'
#' @param cazy_tt A character vector of CAZy family strings, where each string
#'   may contain one or multiple CAZy family codes separated by "|". Each family code
#'   starts with a two-letter prefix (e.g., "GH", "GT").
#'
#' @return A character vector with the CAZy class names corresponding to each family code
#'   extracted from the input vector. The output is of the same length as the total number
#'   of family codes found (splitting by "|").
#'
#' @details
#' The function extracts the first two letters of each CAZy family code to map it
#' to one of the main CAZy classes:
#' \itemize{
#'   \item GT: GlycosylTransferases
#'   \item GH: Glycoside Hydrolases
#'   \item PL: Polysaccharide Lyases
#'   \item CE: Carbohydrate Esterases
#'   \item AA: Auxiliary Activities
#'   \item CB: Carbohydrate-Binding Modules
#' }
#'
#' The input strings can contain multiple families separated by "|".
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   families <- c("GH13|GT2|PL9", "CE4", "AA3|CBM50")
#'   getCazyClass(families)
#' }
#' }
#'
#' @rdname getCazyClass
#' @export
getCazyClass <- function(cazy_tt){
  cazytypes <- c( "GT"="GlycosylTransferases",
                  "GH"="Glycoside Hydrolases",
                  "PL"="Polysaccharide Lyases",
                  "CE"="Carbohydrate Esterases",
                  "AA"="Auxiliary Activities",
                  "CB"="Carbohydrate-Binding Modules")
  split_list <- str_split(cazy_tt, pattern = "\\|") %>% unlist %>% substr(1,2)
  cazy_class <- cazytypes[split_list]
  return(cazy_class)
}
