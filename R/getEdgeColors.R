getEdgeColors <- function(param, p.value, Estimate, plim_plot, custom_colors,
                          col_ctrl = C_CTRL, col_case = C_CASE, col_ns = C_NS){
  cols <- ifelse(Estimate < 0, col_ctrl, col_case)
  if(is.null(custom_colors)){
    cols <- ifelse(p.value > plim_plot, col_ns, cols)
  }else{
    cols <- ifelse(param=="b",
                   cols,
                   ifelse(grepl("^a", param),
                          ifelse(is.na(custom_colors$indirect),
                                 ifelse(p.value > plim_plot, col_ns, cols),
                                 ifelse(custom_colors$indirect, cols, col_ns)
                          ),
                          ifelse(is.na(custom_colors$direct),
                                 ifelse(p.value > plim_plot, col_ns, cols),
                                 ifelse(custom_colors$direct, cols, col_ns)
                          )
                   )
    )
  }
  return(cols)
}
