#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param bacorder PARAM_DESCRIPTION
#' @param plim_plot PARAM_DESCRIPTION
#' @param custom_colors PARAM_DESCRIPTION
#' @param col_ctrl PARAM_DESCRIPTION, Default: C_CTRL
#' @param col_case PARAM_DESCRIPTION, Default: C_CASE
#' @param col_ns PARAM_DESCRIPTION, Default: C_NS
#' @param col_other PARAM_DESCRIPTION, Default: C_OTHER
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname getTotalColors
#' @export 
getTotalColors <- function(bacorder, plim_plot, custom_colors,
                           col_ctrl = C_CTRL, col_case = C_CASE, col_ns = C_NS, col_other = C_OTHER){
  cols <- ifelse(bacorder$Estimate < 0, col_ctrl, col_case)
  if(is.null(custom_colors)){
    cols <- ifelse(bacorder$p.value <= plim_plot, cols, col_ns)
  }else{
    if(is.na(custom_colors$total)){
      cols <- ifelse(bacorder$p.value <= plim_plot, cols, col_ns)
    }else if(!custom_colors$total){
      cols <- rep(col_ns, nrow(bacorder))
    }
  }
  cols <- c(cols, col_case, col_other)
  return(cols)
}
