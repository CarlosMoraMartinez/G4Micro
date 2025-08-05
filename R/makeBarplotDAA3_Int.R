makeBarplotDAA3_Int <- function(daalist, outdir, plim=0.05, name="",
                                cond_names=list()){

  if(length(cond_names) == 0) cond_names <- names(daalist)
  cond_names <- paste0(cond_names, "__padj")

  daatab <- lapply(names(daalist), \(x){daalist[[x]]$contrast <- x; return(daalist[[x]])}) %>% bind_rows() %>%
    dplyr::select(taxon, log2FoldChangeShrink, padj, contrast) %>%
    dplyr::mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    gather(key="param", value ="value", padj, log2FoldChangeShrink) %>%
    tidyr::unite(col="tmp", contrast, param, sep="__") %>%
    spread(tmp, value) %>%
    dplyr::mutate(
      AllSig = !!sym(cond_names[1]) <= plim & !!sym(cond_names[2]) <= plim & !!sym(cond_names[3]) <= plim & !!sym(cond_names[4]) <= plim,
      CondAdjusted = !!sym(cond_names[1]) <= plim & !!sym(cond_names[2]) <= plim,
      CondOnly = !!sym(cond_names[2])  > plim & !!sym(cond_names[1]) <= plim,
      CondAdjustedOnly = !!sym(cond_names[2]) <= plim & !!sym(cond_names[1]) > plim,
      BMIAdjusted = !!sym(cond_names[4]) <= plim & !!sym(cond_names[2]) > plim &  !!sym(cond_names[1]) > plim,

      #CondAdjusted = CondAdjusted & !AllSig,
      #CondOnly = CondAdjusted & !AllSig,
      #CondAdjustedOnly = CondAdjusted & !AllSig,
      #BMIAdjusted = CondAdjusted & !AllSig

    ) %>%
    gather(key="param", value="value", -taxon, -AllSig, -CondAdjusted, -CondOnly, -CondAdjustedOnly, -BMIAdjusted) %>%
    separate(param, into=c("Contrast", "param"), sep="__") %>%
    spread(key=param, value = value) %>%
    dplyr::mutate(Contrast= gsub("plus", "+", Contrast),
                  Contrast= gsub("Inter", "*", Contrast),
                  Contrast = gsub("_", " ", Contrast),
                  Contrast = gsub("adj", "adj.", Contrast),
                  Contrast = gsub("Depr", "Depr.", Contrast),
                  taxon = gsub("_", " ", taxon),
                  taxon = gsub("[\\[\\]]", "", taxon),
                  color = ifelse(padj < plim, ifelse(log2FoldChangeShrink < 0, "Down", "Up"), "NS"),
                  color = factor(color, levels=c("Down", "Up", "NS")))

  orderedlevels <- daatab %>% filter(Contrast == "D vs C") %>% arrange(log2FoldChangeShrink) %>% pull(taxon)
  contrastlevels <- c("D vs C", "D vs C adj. Sex*Age + BMI", "BMI", "BMI adj. Sex*Age + MDD")
  daatab <- daatab %>% mutate(taxon=factor(taxon, levels=orderedlevels),
                              Contrast = factor(Contrast, levels = contrastlevels))

  vars2separate <- c("AllSig", "CondAdjusted", "CondOnly", "CondAdjustedOnly", "BMIAdjusted")
  vars2separatenames <- c("All significant",
                          "Different in D vs C after adjusting",
                          "Not significant after correction",
                          "Significant only after correction",
                          "Significant BMI after MDD correction")
  plots <- map2(vars2separate, vars2separatenames,\(x, y){
    tmptab <- daatab %>% dplyr::filter(!!sym(x))
    g1 <- ggplot(tmptab, aes(x=taxon, y=log2FoldChangeShrink, fill=color))+
      facet_grid( . ~ Contrast) +
      geom_hline(yintercept = 0, linetype=2, col=C_NS) +
      scale_fill_manual(values=c(C_CTRL, C_CASE, C_NS)) +
      geom_col() +
      coord_flip() +
      theme_minimal() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(strip.text.y = element_text(size = 10,
                                        colour = "black", angle = 0, face = "italic")) +
      theme(axis.text.x = element_text(size = 8,
                                       colour = "black", angle = 0,
                                       face = "plain"))+
      theme(axis.text.y = element_text(size = 8,
                                       colour = "black", angle = 0,
                                       vjust = 0.5,
                                       face = "italic")) +
      theme(legend.position = 'none')+
      #theme(axis.text.y = element_blank())+
      xlab(y)+
      ylab("") +
      #ylim(c(-9, 9))+
      thin_barplot_lines
    return(list(plot=g1, size=length(unique(tmptab$taxon))))
  })
  names(plots) <- vars2separate
  cw <- cowplot::plot_grid(plotlist = list(plots$CondAdjusted$plot,
                                           plots$CondOnly$plot + theme(strip.text.x = element_blank()),
                                           plots$CondAdjustedOnly$plot+ theme(strip.text.x = element_blank()) + ylab("LFC")
  ),
  rel_heights = c(plots$CondAdjusted$size,
                  plots$CondOnly$size,
                  plots$CondAdjustedOnly$size),
  ncol=1)

  filename = paste0(outdir, "/", name, "_Barplot_LFC_all1.pdf")
  pdf(filename, w=8, h=18)
  print(cw)
  dev.off()

  ggsave(filename = paste0(outdir, "/", name, "_Barplot_LFC_all_BMI.pdf"), plots$BMIAdjusted$plot, w=8, h=8)
  write_tsv(daatab, file =  paste0(outdir, "/", name, "_Barplot_LFC_all1.tsv"))

  ## Cowplot including BMI
  cw <- cowplot::plot_grid(plotlist = list(plots$CondAdjusted$plot,
                                           plots$CondOnly$plot + theme(strip.text.x = element_blank()),
                                           plots$CondAdjustedOnly$plot+ theme(strip.text.x = element_blank()),
                                           plots$BMIAdjusted$plot + theme(strip.text.x = element_blank()) + ylab("LFC")
  ),
  rel_heights = c(plots$CondAdjusted$size, plots$CondOnly$size, plots$CondAdjustedOnly$size, plots$BMIAdjusted$size),
  ncol=1)

  filename = paste0(outdir, "/", name, "_Barplot_LFC_all2.pdf")
  pdf(filename, w=8, h=20)
  print(cw)
  dev.off()

  #cw <- cowplot::plot_grid(plotlist = list(plots$AllSig$plot,
  #                                         plots$CondAdjusted$plot + theme(strip.text.x = element_blank()),
  #                                         plots$CondOnly$plot + theme(strip.text.x = element_blank()),
  #                                         plots$CondAdjustedOnly$plot+ theme(strip.text.x = element_blank()),
  #                                         plots$BMIAdjusted$plot + theme(strip.text.x = element_blank()) + ylab("LFC")
  #),
  #rel_heights = c(plots$AllSig$size,
  #                plots$CondAdjusted$size,
  #                plots$CondOnly$size,
  #                plots$CondAdjustedOnly$size,
  #                plots$BMIAdjusted$size),
  #ncol=1)
  #
  #filename = paste0(outdir, "/", name, "_Barplot_LFC_all3.pdf")
  #pdf(filename, w=8, h=20)
  #print(cw)
  #dev.off()
  return(list(plots=plots, tab=daatab))
}
