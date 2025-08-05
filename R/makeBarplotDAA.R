#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param daalist PARAM_DESCRIPTION
#' @param outdir PARAM_DESCRIPTION
#' @param plim PARAM_DESCRIPTION, Default: 0.05
#' @param name PARAM_DESCRIPTION, Default: ''
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{filter}}
#'  \code{\link[tidyr]{unite}}
#'  \code{\link[cowplot]{plot_grid}}
#' @rdname makeBarplotDAA
#' @export 
#' @importFrom dplyr select mutate filter
#' @importFrom tidyr unite
#' @importFrom cowplot plot_grid
makeBarplotDAA <- function(daalist, outdir, plim=0.05, name=""){
  daatab <- lapply(names(daalist), \(x){daalist[[x]]$contrast <- x; return(daalist[[x]])}) %>% bind_rows() %>%
    dplyr::select(taxon, log2FoldChangeShrink, padj, contrast) %>%
    dplyr::mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    gather(key="param", value ="value", padj, log2FoldChangeShrink) %>%
    tidyr::unite(col="tmp", contrast, param, sep="__") %>%
    spread(tmp, value) %>%
    dplyr::mutate(
      CondAdjusted = D_vs_C_adj_BMIplusAge__padj <= plim & D_vs_C__padj <= plim,
      CondOnly = D_vs_C_adj_BMIplusAge__padj > plim & D_vs_C__padj <= plim,
      CondAdjustedOnly = D_vs_C_adj_BMIplusAge__padj <= plim & D_vs_C__padj > plim,
      BMIAdjusted = BMI_adj_DeprplusAge__padj <= plim,
      AgeAdjusted = Age_adj_DeprplusBMI__padj <= plim
    ) %>%
    gather(key="param", value="value", -taxon, -CondAdjusted, -CondOnly, -CondAdjustedOnly, -BMIAdjusted, -AgeAdjusted ) %>%
    separate(param, into=c("Contrast", "param"), sep="__") %>%
    spread(key=param, value = value) %>%
    dplyr::mutate(Contrast= gsub("plus", "+", Contrast),
                  Contrast = gsub("_", " ", Contrast),
                  Contrast = gsub("adj", "adj.", Contrast),
                  taxon = gsub("_", " ", taxon),
                  taxon = gsub("[\\[\\]]", "", taxon),
                  color = ifelse(padj < plim, ifelse(log2FoldChangeShrink < 0, "Down", "Up"), "NS"),
                  color = factor(color, levels=c("Down", "Up", "NS")))

  orderedlevels <- daatab %>% filter(Contrast == "D vs C") %>% arrange(log2FoldChangeShrink) %>% pull(taxon)
  contrastlevels <- c("D vs C", "D vs C adj. BMI+Age", "BMI adj. Depr+Age", "Age adj. Depr+BMI")
  daatab <- daatab %>% mutate(taxon=factor(taxon, levels=orderedlevels),
                              Contrast = factor(Contrast, levels = contrastlevels))

  vars2separate <- c("CondAdjusted", "CondOnly", "CondAdjustedOnly", "BMIAdjusted", "AgeAdjusted")
  vars2separatenames <- c("Different in D vs C after adjusting", "Not significant after BMI+Age correction",
                          "Significant only after correction", "Significant BMI after Depr+Age correction",
                          "Significant Age after Depr+BMI correction")
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
      ylim(c(-9, 9))+
      thin_barplot_lines
    return(list(plot=g1, size=length(unique(tmptab$taxon))))
  })
  names(plots) <- vars2separate
  cw <- cowplot::plot_grid(plotlist = list(plots$CondAdjusted$plot,
                                           plots$CondOnly$plot + theme(strip.text.x = element_blank()),
                                           plots$CondAdjustedOnly$plot+ theme(strip.text.x = element_blank()) + ylab("LFC")
  ),
  rel_heights = c(plots$CondAdjusted$size, plots$CondOnly$size, plots$CondAdjustedOnly$size),
  ncol=1)

  filename = paste0(outdir, "/", name, "_Barplot_LFC_all1.pdf")
  pdf(filename, w=8, h=16)
  print(cw)
  dev.off()
  ggsave(filename = paste0(outdir, "/", name, "_Barplot_LFC_all_BMI.pdf"), plots$BMIAdjusted$plot, w=8, h=12)
  ggsave(filename = paste0(outdir, "/", name, "_Barplot_LFC_all_Age.pdf"), plots$AgeAdjusted$plot, w=8, h=6)
  write_tsv(daatab, file =  paste0(outdir, "/", name, "_Barplot_LFC_all1.tsv"))
  return(list(plots=plots, tab=daatab))
}
