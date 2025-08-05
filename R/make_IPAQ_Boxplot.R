#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param divtab PARAM_DESCRIPTION
#' @param var PARAM_DESCRIPTION, Default: 'IPAQ_act_fisica'
#' @param test2show PARAM_DESCRIPTION, Default: 'wilcox.test'
#' @param alpha_indices PARAM_DESCRIPTION
#' @param outdir PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION
#' @param signif_levels PARAM_DESCRIPTION, Default: c(`***` = 0.001, `**` = 0.01, `*` = 0.05, ns = 1.1)
#' @param w PARAM_DESCRIPTION, Default: 10
#' @param h PARAM_DESCRIPTION, Default: 4
#' @param correct_pvalues PARAM_DESCRIPTION, Default: TRUE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{filter}}, \code{\link[dplyr]{select}}, \code{\link[dplyr]{mutate}}
#'  \code{\link[ggsignif]{stat_signif}}
#' @rdname make_IPAQ_Boxplot
#' @export 
#' @importFrom dplyr filter select mutate
#' @importFrom ggsignif stat_signif
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
    theme_pubclean() +
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
    theme_pubclean() +
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
