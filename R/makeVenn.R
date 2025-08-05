
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
