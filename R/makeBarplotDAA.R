#' @title Generate and Save Barplots for Differential Abundance Analysis Results
#' @description Creates multiple barplots showing log2 fold changes of taxa across contrasts with significance filtering.
#'   The function combines results from a list of differential abundance analysis outputs, applies p-value filtering, annotates significance categories,
#'   and generates facetted barplots that are saved as PDF files along with a summary TSV table.
#' @param daalist A named list of data frames, each representing differential abundance analysis results for a specific contrast.
#'   Each data frame must contain columns: \code{taxon}, \code{log2FoldChangeShrink}, and \code{padj}.
#' @param outdir A character string specifying the directory where output files (PDFs and TSV) will be saved.
#' @param plim Numeric threshold for adjusted p-value significance filtering. Default is 0.05.
#' @param name Character string prefix for output file names. Default is an empty string.
#' @return A list with two elements:
#'   \item{plots}{A named list of \code{ggplot} objects representing different subsets of significant taxa and their fold changes.}
#'   \item{tab}{A processed data frame summarizing the differential abundance results, including significance categories and plotting colors.}
#' @details
#'   The function first combines and reformats the input list of DAA results, adds categories based on adjusted p-value thresholds,
#'   and creates multiple barplots to visualize taxa log2 fold changes across contrasts.
#'   Significant taxa are highlighted with colors indicating direction (Up/Down), and the plots are arranged vertically with cowplot.
#'   The results are saved as PDF files and a TSV summary in the specified output directory.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Example usage with a list of DAA results named by contrast
#'   daa_list <- list(
#'     "D_vs_C" = daa_result1,
#'     "D_vs_C_adj_BMIplusAge" = daa_result2,
#'     "BMI_adj_DeprplusAge" = daa_result3,
#'     "Age_adj_DeprplusBMI" = daa_result4
#'   )
#'   output <- makeBarplotDAA(daalist = daa_list, outdir = "results/plots", plim = 0.05, name = "MyStudy")
#'   print(output$plots$CondAdjusted)
#' }
#' }
#' @seealso
#' \code{\link[dplyr]{select}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{filter}}, \code{\link[tidyr]{unite}}, \code{\link[cowplot]{plot_grid}}, \code{\link[ggplot2]{ggplot}}
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
