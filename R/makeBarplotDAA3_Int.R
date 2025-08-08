#' @title Barplot Visualization of Differential Abundance Analysis Results
#' @description Generates barplots of log2 fold changes from differential abundance analyses (DAA) across multiple contrasts,
#' highlighting taxa that pass specific significance thresholds under different model adjustments.
#' @param daalist A named list of data.frames, each representing a DAA result for one contrast.
#' Each data.frame must contain columns: \code{taxon}, \code{log2FoldChangeShrink}, and \code{padj}.
#' @param outdir Character string specifying the directory where output plots and tables will be saved.
#' @param plim Numeric value for adjusted p-value cutoff used to determine significance. Default is 0.05.
#' @param name Character string prefix for output filenames. Default is an empty string.
#' @param cond_names Character vector or list of adjusted p-value column names used internally to evaluate significance across contrasts.
#' If empty (default), names are derived from \code{daalist}.
#' @return A list containing:
#' \item{plots}{A named list of ggplot objects for different significance subsets.}
#' \item{tab}{A data.frame with merged DAA results and added significance flags.}
#' @details
#' The function merges multiple DAA results, reshapes the data, and creates barplots of log2 fold changes for taxa.
#' It identifies taxa significant under different model conditions:
#' \itemize{
#'   \item \code{AllSig}: significant in all contrasts.
#'   \item \code{CondAdjusted}: significant in the first two contrasts.
#'   \item \code{CondOnly}: significant only in the first contrast.
#'   \item \code{CondAdjustedOnly}: significant only in the second contrast.
#'   \item \code{BMIAdjusted}: significant only after BMI adjustment.
#' }
#' The plots are saved as PDF files in \code{outdir}, and a TSV table with results is also written.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Example usage assuming daalist contains differential abundance result data.frames:
#'   # daalist <- list(
#'   #   "D_vs_C" = daa_res1,
#'   #   "D_vs_C_adj" = daa_res2,
#'   #   "BMI" = daa_res3,
#'   #   "BMI_adj" = daa_res4
#'   # )
#'   # makeBarplotDAA3_Int(daalist, outdir = "results/", plim = 0.05, name = "example")
#' }
#' }
#' @seealso
#' \code{\link[dplyr]{select}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{filter}},
#' \code{\link[tidyr]{unite}}, \code{\link[tidyr]{spread}}, \code{\link[cowplot]{plot_grid}}
#' @rdname makeBarplotDAA3_Int
#' @export
#' @importFrom dplyr select mutate filter pull
#' @importFrom tidyr unite gather separate spread
#' @importFrom cowplot plot_grid
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
