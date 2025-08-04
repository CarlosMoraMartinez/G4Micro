makePCoA <- function(phobj, pcoa.bray, evals, var2color="Condition", name = "Bray-Curtis", extradims= 2:5, labelsamples="sampleID"){
  gs <- lapply(extradims, FUN=function(axis){
    gg <- plot_ordination(phobj, pcoa.bray, color = var2color,
                          title = name, axes=c(1, axis)) +
      #coord_fixed(sqrt(evals[2] / evals[1])) +
      #scale_color_manual(values=palette2)+
      stat_ellipse(level=0.95, linetype=2, alpha = 0.8, na.rm = TRUE) +
      geom_point(size = 2) +
      geom_text_repel(aes_string(label = labelsamples)) +
      mytheme
    # if(!is.numeric(gg$data[, var2color])){
    #   gg <- gg + scale_color_lancet() + scale_fill_lancet()
    # }
    return(gg)
  })
  cw <- cowplot::plot_grid(plotlist=gs)
  return(cw)
}
