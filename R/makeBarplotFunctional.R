makeBarplotFunctional <- function(tab, plim=0.001, name, outdir,
                                  include_longnames = FALSE,
                                  w=12, h=12){
  tab <- tab %>% rownames_to_column("Pathway")
  if(include_longnames & "process_name" %in% names(tab)){
    pnames <- gsub(" \\[[a-zA-Z0-9\\.\\]:\\- ]+", "", tab$process_name, perl=T) %>% gsub(" / ", "/", .)
    tab$Pathway <- paste(tab$Pathway, pnames, sep=": ")
  }
  tab2 <- tab %>% dplyr::filter(padj < plim) %>%
    dplyr::arrange(log2FoldChange) %>%
    dplyr::mutate(Pathway= factor(Pathway, levels = Pathway),
                  direction = ifelse(log2FoldChange<0, "less abundant", "more abundant"),
                  #direction = factor(direction, levels=legendnames),
                  `-10log(padj)` = -10*log10(padj)
    )


  g1 <- ggplot(tab2, aes(x=Pathway, y = log2FoldChange, fill=direction, col=direction ))+ #`-10log(padj)`
    geom_hline(yintercept = 0, linetype=3, col="gray80") +
    geom_bar(stat="identity") +
    coord_flip()+
    mytheme +
    scale_fill_manual(values = c("#A1C6EA", "#FD8B2F")) +
    scale_color_manual(values = c("#A1C6EA", "#FD8B2F")) +
    #scale_fill_gradient2(low="steelblue1", "#DDDDDD", high="tomato")+
    #scale_color_gradient2(low="steelblue1", "#DDDDDD", high="tomato")+
    theme(strip.text.y = element_text(size = 10,
                                      colour = "black", angle = 0, face = "bold"))+
    guides(fill=guide_legend(title="LFC"), colour=guide_legend(title="LFC"))+
    theme(panel.grid = element_blank())+
    ggtitle(paste0(name, " - proc. with padj<", as.character(plim)))+
    theme(axis.text.y = element_text(size = 12,
                                     colour = "black", angle = 0, face = "plain"))

  ggsave(filename = paste0(outdir, "/", name, ".pdf"), g1, width = w, height = h)
  return(g1)
}
