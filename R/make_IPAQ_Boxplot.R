#' @title Generate Boxplots of Alpha Diversity by IPAQ Activity Levels
#' @description This function generates boxplots of selected alpha diversity indices
#' across different physical activity levels based on the IPAQ
#' (International Physical Activity Questionnaire).
#' It returns two plots: one aggregated across all samples,
#' and another split by condition. It also performs pairwise significance
#' testing and annotates the plots accordingly.
#' @param divtab A data frame containing alpha diversity indices, sample metadata, and a column with IPAQ activity levels (e.g., "Low", "Mid", "High").
#' @param var A character string indicating the name of the column in `divtab` that contains IPAQ activity levels. Default is `"IPAQ_act_fisica"`.
#' @param test2show The statistical test to use for comparing groups (e.g., `"wilcox.test"` or `"t.test"`). Default is `"wilcox.test"`.
#' @param alpha_indices A character vector with the names of the alpha diversity indices to plot (e.g., `"Shannon"`, `"Chao1"`).
#' @param outdir A character string indicating the path to the output directory where PDF files will be saved.
#' @param name A character string used to personalize output file names (e.g., a cohort or experiment name).
#' @param signif_levels A named numeric vector indicating the significance thresholds for p-values to be displayed on the plots (e.g., `c("***"=0.001, "**"=0.01, "*"=0.05, "ns"=1.1)`). Default matches this example.
#' @param w Width (in inches) of the output PDF plots. Default is 10.
#' @param h Height (in inches) of the single-row plot. The condition-specific plot is twice this height. Default is 4.
#' @param correct_pvalues Logical; if `TRUE`, p-values are corrected for multiple testing using the Bonferroni method. Default is `TRUE`.
#'
#' @return A list containing two ggplot2 objects:
#' \describe{
#'   \item{g1}{Boxplot of alpha diversity indices by IPAQ activity level (all samples pooled).}
#'   \item{g2}{Boxplot of alpha diversity indices by IPAQ activity level split by condition.}
#' }
#'
#' @details
#' This function prepares the input table by filtering missing values in the
#' IPAQ activity column and reshaping it into long format.
#' Then, it creates two faceted boxplots for the provided alpha diversity
#' indices, using the `ggsignif` package to overlay statistical significance of
#' pairwise comparisons between IPAQ levels. P-values can be optionally corrected
#' using the Bonferroni method depending on the number of comparisons. The plots are
#' saved as PDF files in the specified output directory.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  data <- read.csv("diversity_table.csv")
#'   make_IPAQ_Boxplot(
#'     divtab = data,
#'     alpha_indices = c("Shannon", "Chao1"),
#'     outdir = "results/plots",
#'     name = "Cohort1"
#'   )
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{filter}}, \code{\link[dplyr]{select}}, \code{\link[dplyr]{mutate}}
#'  \code{\link[ggsignif]{stat_signif}}
#' @rdname make_IPAQ_Boxplot
#' @export
#' @importFrom dplyr filter select mutate
#' @importFrom ggsignif stat_signif
#' @importFrom ggpubr theme_pubclean
make_IPAQ_Boxplot <- function(divtab, var = "IPAQ_act_fisica",
                              test2show = "wilcox.test",
                              alpha_indices,
                              outdir,
                              name,
                              signif_levels=c("***"=0.001, "**"=0.01, "*"=0.05, "ns"=1.1),
                              w=10, h=4,
                              correct_pvalues=TRUE){
  divtab2 <- divtab %>% dplyr::filter(!is.na(IPAQ_act_fisica)) %>%
    dplyr::select(sampleID, Condition, all_of(c(var, alpha_indices))) %>%
    dplyr::mutate(IPAQ_act_fisica=factor(IPAQ_act_fisica, levels=c("Low", "Mid", "High"))) %>%
    gather("index", "value", all_of(alpha_indices)) %>%
    dplyr::mutate(index = factor(index, levels=alpha_indices))

  comp <- combn(unique(divtab2[, var]), 2, simplify = F) %>% lapply(as.character)
  num_comparisons <- length(alpha_indices)*length(comp)
  if(correct_pvalues & num_comparisons>1){
    signif_levels_bonferroni <- c(signif_levels[1:3]/num_comparisons, signif_levels[4])
  }else{
    signif_levels_bonferroni <- signif_levels
  }
  g1 <- ggplot(divtab2, aes(x=!!sym(var), y=value, fill=!!sym(var), col=!!sym(var))) +
    facet_wrap(. ~ index, scales = "free", nrow = 1) +
    geom_boxplot(alpha = 0.7) +
    labs(title = var, x = '') +
    ggpubr::theme_pubclean() +
    mytheme +
    ggsignif::stat_signif(test=test2show, na.rm=T, comparisons = comp,
                          step_increase=0.06,
                          tip_length = 0.01,
                          map_signif_level=signif_levels_bonferroni,
                          vjust=0.3,
                          color = "black"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_blank())
  ggsave(filename = paste0(outdir, "/alphadiv_IPAQ_", name, ".pdf"), g1, height = h, width = w )

  g2 <- ggplot(divtab2, aes(x=!!sym(var), y=value, fill=!!sym(var), col=!!sym(var))) +
    facet_wrap(Condition ~ index, scales = "free", nrow = 2) +
    geom_boxplot(alpha = 0.7) +
    labs(title = var, x = '') +
    ggpubr::theme_pubclean() +
    mytheme +
    ggsignif::stat_signif(test=test2show, na.rm=T, comparisons = comp,
                          step_increase=0.06,
                          tip_length = 0.01,
                          map_signif_level=signif_levels_bonferroni,
                          vjust=0.3,
                          color = "black"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_blank())
  ggsave(filename = paste0(outdir, "/alphadiv_IPAQ_Condition_", name, ".pdf"), g2, height = h*2, width = w )
  return(list(g1=g1, g2=g2))

}
