
makeTopTaxaPlot <- function(nodeprops, ntop, outdir, name,
                            w=8, h=8, manual_scale=c("#A1C6EA", "#FD8B2F", "#A5ABBD")){
  top_taxa <- nodeprops %>% dplyr::group_by(class) %>% top_n(sum_measures, n=ntop) %>%
    gather("measure", "value", closeness, betweenness, degree) %>%
    dplyr::mutate(vertex = gsub("_", " ", vertex))
  fname <- paste0(outdir, name, "_top", as.character(ntop), "Table.tsv")
  write_tsv(x = top_taxa, file = fname)
  gprops <- ggplot(top_taxa, aes(x=vertex, y=value, fill=color))+
    geom_col()+
    coord_flip()+
    facet_grid(class~measure, scales = "free")+
    theme_pubclean()+
    scale_fill_manual(values=manual_scale)+
    theme(axis.text.y=element_text(face="italic"))
  fname <- paste0(outdir, name, "_top", as.character(ntop), "Barplot.pdf")
  ggsave(filename = fname, gprops, width = w, height = h)
  return(gprops)
}
