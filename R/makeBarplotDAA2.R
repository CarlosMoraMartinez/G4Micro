#' @title Create Barplots for Differential Abundance Analysis Results
#' @description
#' Generates barplots visualizing log2 fold changes of taxa across multiple contrasts
#' from a list of differential abundance analysis results. It highlights taxa significant
#' under different conditions including adjustments for BMI and depression status.
#'
#' @param daalist A named list of data frames, each representing differential abundance
#' results for a contrast. Each data frame must contain columns: `taxon`,
#' `log2FoldChangeShrink`, and `padj`.
#' @param outdir Character. Directory path where plots and tables will be saved.
#' @param plim Numeric. Adjusted p-value threshold for significance (default 0.05).
#' @param name Character. Prefix for output filenames (default "").
#'
#' @return A list with two elements:
#' \item{plots}{A list of ggplot2 objects for different significance categories.}
#' \item{tab}{A processed data frame with combined results and significance flags.}
#'
#' @details
#' This function merges differential abundance results from multiple contrasts,
#' computes significance categories based on p-value thresholds,
#' and generates multiple barplots showing log2 fold changes (LFC) colored by direction
#' (Up/Down/Not Significant). The plots are saved as PDF files in the specified output directory.
#'
#' The contrasts considered include:
#' - "D vs C" (Disease vs Control)
#' - "D vs C adj. BMI" (adjusted for BMI)
#' - "BMI" (BMI effect)
#' - "BMI adj. Depr." (BMI adjusted for Depression)
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming daalist is prepared as a named list of data.frames with required columns
#'   result <- makeBarplotDAA2(daalist = my_da_results, outdir = "plots", plim = 0.05, name = "my_analysis")
#'   # Access plots:
#'   print(result$plots$CondAdjusted$plot)
#' }
#' }
#'
#' @seealso
#' \code{\link[dplyr]{select}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{filter}},
#' \code{\link[tidyr]{unite}}, \code{\link[cowplot]{plot_grid}}, \code{\link[ggplot2]{ggplot}}
#'
#' @rdname makeBarplotDAA2
#' @export
#' @importFrom dplyr select mutate filter
#' @importFrom tidyr unite gather spread separate
#' @importFrom cowplot plot_grid
makeBarplotDAA2 <- function(daalist, outdir, plim=0.05, name=""){
  daatab <- lapply(names(daalist), \(x){daalist[[x]]$contrast <- x; return(daalist[[x]])}) %>% bind_rows() %>%
    dplyr::select(taxon, log2FoldChangeShrink, padj, contrast) %>%
    dplyr::mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    gather(key="param", value ="value", padj, log2FoldChangeShrink) %>%
    tidyr::unite(col="tmp", contrast, param, sep="__") %>%
    spread(tmp, value) %>%
    dplyr::mutate(
      CondAdjusted = D_vs_C_adj_BMI__padj <= plim & D_vs_C__padj <= plim,
      CondOnly = D_vs_C_adj_BMI__padj > plim & D_vs_C__padj <= plim,
      CondAdjustedOnly = D_vs_C_adj_BMI__padj <= plim & D_vs_C__padj > plim,
      BMIAdjusted = BMI_adj_Depr__padj <= plim & D_vs_C_adj_BMI__padj > plim &  D_vs_C__padj > plim
    ) %>%
    gather(key="param", value="value", -taxon, -CondAdjusted, -CondOnly, -CondAdjustedOnly, -BMIAdjusted) %>%
    separate(param, into=c("Contrast", "param"), sep="__") %>%
    spread(key=param, value = value) %>%
    dplyr::mutate(Contrast= gsub("plus", "+", Contrast),
                  Contrast = gsub("_", " ", Contrast),
                  Contrast = gsub("adj", "adj.", Contrast),
                  Contrast = gsub("Depr", "Depr.", Contrast),
                  taxon = gsub("_", " ", taxon),
                  taxon = gsub("[\\[\\]]", "", taxon),
                  color = ifelse(padj < plim, ifelse(log2FoldChangeShrink < 0, "Down", "Up"), "NS"),
                  color = factor(color, levels=c("Down", "Up", "NS")))

  orderedlevels <- daatab %>% filter(Contrast == "D vs C") %>% arrange(log2FoldChangeShrink) %>% pull(taxon)
  contrastlevels <- c("D vs C", "D vs C adj. BMI", "BMI", "BMI adj. Depr.")
  daatab <- daatab %>% mutate(taxon=factor(taxon, levels=orderedlevels),
                              Contrast = factor(Contrast, levels = contrastlevels))

  vars2separate <- c("CondAdjusted", "CondOnly", "CondAdjustedOnly", "BMIAdjusted")
  vars2separatenames <- c("Different in D vs C after adjusting", "Not significant after BMI correction",
                          "Significant only after correction", "Significant BMI after Depr correction")
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
  write_tsv(daatab, file =  paste0(outdir, "/", name, "_Barplot_LFC_all1.tsv"))

  ## Cowplot including BMI
  cw <- cowplot::plot_grid(plotlist = list(plots$CondOnly$plot + theme(strip.text.x = element_blank()),
                                           plots$CondAdjustedOnly$plot+ theme(strip.text.x = element_blank()),
                                           plots$BMIAdjusted$plot + theme(strip.text.x = element_blank()) + ylab("LFC")
  ),
  rel_heights = c(plots$CondOnly$size, plots$CondAdjustedOnly$size, plots$BMIAdjusted$size),
  ncol=1)

  filename = paste0(outdir, "/", name, "_Barplot_LFC_all2.pdf")
  pdf(filename, w=8, h=20)
  print(cw)
  dev.off()
  return(list(plots=plots, tab=daatab))
}
