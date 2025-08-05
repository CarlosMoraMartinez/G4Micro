plotNetWithTopNames <- function(net_obj, nodeprops, net_class, ntop=10, outdir, name, add_labels=TRUE,
                                w=8, h=8){

  top_taxa <- nodeprops %>% dplyr::filter(class==net_class) %>% top_n(sum_measures, n=10) %>%
    pull(vertex)

  net_vnames <- names(V(net_obj$net))
  if(add_labels){
    ver_names <- ifelse(net_vnames %in% top_taxa, net_vnames, "") %>%
      process_names %>% process_names2
  }else{
    ver_names <- ""
  }
  fname <- paste0(outdir, name, ".pdf")
  V(net_obj$net)$label.cex <- 1.2
  V(net_obj$net)$label.face <- "italic"
  pdf(fname, width=w, height=h)
  p1 <- plot(net_obj$net, vertex.label.size=6, vertex.label=ver_names,
             vertex.label.color="black",
             vertex.label.dist=0.5,
             vertex.size=ifelse(net_vnames %in% top_taxa, 8, 4),
             vertex.color=net_obj$cols)
  dev.off()

}
