
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param vars2venn PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION, Default: 'VennDiagram'
#' @param opt PARAM_DESCRIPTION
#' @param w PARAM_DESCRIPTION, Default: 5
#' @param h PARAM_DESCRIPTION, Default: 5
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname makeVenn
#' @export 
makeVenn <- function(vars2venn, name="VennDiagram", opt, w=5, h=5){
  gv <- ggvenn(
    vars2venn, columns = names(vars2venn),
    stroke_size = 0.5,
    stroke_color = C_NS,
    fill_color = c(C_CASE, C_CASE2, C_CTRL2, C_CTRL),show_elements = F
  )
  ggsave(filename = paste0(opt$out, name, ".pdf"), gv, width = w, height = h)
  return(gv)
}
