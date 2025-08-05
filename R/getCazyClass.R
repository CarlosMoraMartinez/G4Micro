#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param cazy_tt PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
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
