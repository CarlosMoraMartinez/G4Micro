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
