compareLFCContrats3 <- function(contrastlist, firstContrast,
                                contrastNamesOrdered, mainContrastName,
                                plim_select= 0.000001, plim_plot=0.05,
                                name2remove = "",
                                resdfname="resdf", outdir = "./", name="LFC_compare3", w=8, h=4,
                                scale_mode="fixed"){
  alldeatables <- map(names(contrastlist),
                      \(x) contrastlist[[x]][[resdfname]] %>% mutate(comparison=x)) %>%
    bind_rows()
  alldeatables <- rbind(alldeatables, firstContrast[[resdfname]] %>%
                          mutate(comparison=mainContrastName)) %>%
    mutate(taxon = gsub("_", " ", taxon))

  tax2plot <- alldeatables %>% dplyr::select(taxon, padj, comparison) %>%
    spread(key=comparison, value = padj) %>%
    filter_if(is.numeric, all_vars(. < plim_select)) %>%
    pull(taxon)

  tax2plot %>% length
  taxorder <- firstContrast[[resdfname]] %>% arrange(desc(log2FoldChangeShrink)) %>%
    mutate(taxon = gsub("_", " ", taxon)) %>%
    filter(taxon %in% tax2plot) %>% pull(taxon)

  alldeatables_filt <- alldeatables %>%
    dplyr::filter(taxon %in% taxorder) %>%
    dplyr::mutate(taxon=factor(taxon, levels=taxorder),
                  comparison=gsub(name2remove, "", comparison),
                  comparison=gsub("_", " ", comparison),
                  comparison=factor(comparison,
                                    levels=contrastNamesOrdered),
                  UpOrDown = ifelse(log2FoldChangeShrink<0, "Down", "Up"),
                  UpOrDownSig = ifelse(padj<=plim_plot & !is.na(padj), UpOrDown, "NS"),
                  UpOrDownSig = factor(UpOrDownSig, levels=c("Up", "NS", "Down"))
    )

  g1<-ggplot(alldeatables_filt, aes(x=taxon, y=log2FoldChangeShrink, fill=UpOrDownSig)) +
    facet_grid(. ~ comparison, scales=scale_mode) +
    geom_hline(yintercept = 0, linetype=3, col="gray80") +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c(C_CASE,C_CTRL))+
    mytheme +
    theme_classic()+
    theme(strip.text.x = element_text(size = 12,
                                      colour = "black", angle = 0, face = "bold"))+
    theme(axis.text.y = element_text(size = 12,
                                     colour = "black", angle = 0,
                                     face = "italic", hjust=1, vjust=0.3))+
    theme(axis.title.x = element_text(size = 12))+
    coord_flip() +
    thin_barplot_lines

  ggsave(filename = paste0(outdir, "/", name, "3.pdf"), g1, width = w, height = h)
  return(g1)
}
