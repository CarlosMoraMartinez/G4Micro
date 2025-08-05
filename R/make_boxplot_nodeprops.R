#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param nodeprops PARAM_DESCRIPTION
#' @param outdir PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION
#' @param w PARAM_DESCRIPTION, Default: 12
#' @param h PARAM_DESCRIPTION, Default: 6
#' @param signif_levels PARAM_DESCRIPTION, Default: c(`***` = 0.001, `**` = 0.01, `*` = 0.05, ns = 1.1)
#' @param correct_pvals PARAM_DESCRIPTION, Default: TRUE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{select}}
#'  \code{\link[ggsignif]{stat_signif}}
#' @rdname make_boxplot_nodeprops
#' @export 
#' @importFrom dplyr mutate select
#' @importFrom ggsignif stat_signif
make_boxplot_nodeprops <- function(nodeprops, outdir, name, w=12, h=6,
                                   signif_levels=c("***"=0.001, "**"=0.01, "*"=0.05, "ns"=1.1),
                                   correct_pvals = TRUE){
  propslong <- nodeprops %>%
    dplyr::mutate(closeness = ifelse(closeness==1, NA, closeness)) %>%
    gather("property", "value", closeness, betweenness, degree)

  tests <- map(c("closeness", "betweenness", "degree"),
               \(x)getTestsForAllCombinations(nodeprops[, x], nodeprops[, "class"]) %>%
                 dplyr::mutate(measure=x) %>%
                 dplyr::select(measure, everything())
  ) %>% bind_rows
  fname <- paste0(outdir, "/", name, "_stats.tsv")
  write_tsv(tests, file = fname)

  num_comparisons <- length(unique(tests$groups_compared)>1)
  if(correct_pvals & num_comparisons>1){
    signif_levels <- c(signif_levels[1:3]/num_comparisons, signif_levels[4])
  }

  comp <- combn(unique(as.character(propslong$class)), 2, simplify = F)

  (g1 <- ggplot(propslong, aes(x=class, y=value, fill=class))+
      facet_wrap(~ property, scales="free")+
      geom_violin(alpha=0.5, outlier.shape = NA)+
      geom_boxplot(width=0.2, fill="darkgray",
                   outlier.shape = NA,
                   notchwidth = 0.5, notch=F)+
      ggsignif::stat_signif(test="wilcox.test", na.rm=T, comparisons = comp,
                            step_increase=0.05,
                            tip_length = 0.02,
                            map_signif_level=signif_levels,
                            vjust=0.3,
                            color = "black"
      ) +
      #coord_flip() +
      theme_pubclean()+
      theme(axis.text.x = element_text(size = 12,
                                       colour = "black", angle = 45,
                                       face = "plain", hjust=1, vjust=1))
  )

  fname <- paste0(outdir, name, "_boxplot.pdf")
  ggsave(filename = fname, g1, width = w, height = h)
  return(list(tests=tests, plot=g1))
}
