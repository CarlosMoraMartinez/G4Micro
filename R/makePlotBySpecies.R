#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param bacnames PARAM_DESCRIPTION
#' @param df_all PARAM_DESCRIPTION
#' @param outdir PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION
#' @param quantvar PARAM_DESCRIPTION, Default: 'IMC_log'
#' @param quantvar_name PARAM_DESCRIPTION, Default: 'log(IMC)'
#' @param corrmethod PARAM_DESCRIPTION, Default: 'pearson'
#' @param plim_plot PARAM_DESCRIPTION, Default: 0.05
#' @param w PARAM_DESCRIPTION, Default: 8
#' @param h PARAM_DESCRIPTION, Default: 10
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[assertthat]{assert_that}}
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{summarise}}
#'  \code{\link[broom]{reexports}}
#'  \code{\link[cowplot]{plot_grid}}
#' @rdname makePlotBySpecies
#' @export 
#' @importFrom assertthat assert_that
#' @importFrom dplyr select summarise
#' @importFrom broom glance
#' @importFrom cowplot plot_grid
makePlotBySpecies <- function(bacnames, df_all, outdir, name, quantvar="IMC_log",
                              quantvar_name = "log(IMC)",
                              corrmethod="pearson", plim_plot=0.05, w=8, h=10){
  assertthat::assert_that(all(bacnames %in% names(df_all)))
  vars2factor <- bacnames %>% gsub("_", " ", .) %>%
    gsub("sp ", "sp. ", .) #%>%
  # paste0("italic('", ., "')")

  vars2factor <- factor(vars2factor, levels = vars2factor)

  df2box <- df_all %>% dplyr::select(all_of(c(bacnames, "Condition", mediator_name))) %>%
    gather(key="taxon", value="abundance", bacnames) %>%
    mutate(taxon = gsub("_", " ", taxon) %>%
             gsub("sp ", "sp. ", .) %>%
             #paste0("italic('", ., "')") %>%
             factor(levels=vars2factor)) %>%
    mutate(Condition = ifelse(Condition == "Control", "C", "D"),
           Condition = factor(Condition, levels=c("C", "D")))
  dfmeans <- df2box %>% group_by(Condition, taxon) %>% dplyr::summarise(abundance=mean(abundance))
  g1 <- ggplot(df2box, aes(x=Condition, y=abundance, col=Condition))+
    facet_grid( taxon ~ .)+
    geom_boxplot(outlier.alpha = 0, notch = F, width=0.1, size=1.5, varwidth = T)+
    geom_point(data=dfmeans, aes(x=Condition, y=abundance, size=1.8)) +
    coord_flip() +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(strip.text.y = element_text(size = 10,
                                      colour = "black", angle = 0, face = "italic")) +
    theme(axis.text.x = element_text(size = 8,
                                     colour = "black", angle = 0,
                                     face = "plain"))+
    theme(axis.text.y = element_text(size = 8,
                                     colour = "black", angle = 0,
                                     face = "plain")) +
    #theme(axis.text.y = element_blank())+
    ylab("Taxon abundance")
  ggsave(filename = paste0(outdir, "/", name, ".pdf"), g1, width = w, height = h)

  mods <- do.call('rbind', map(bacnames, \(x){
    form <- as.formula(paste0(quantvar, " ~ ", x))
    mm <- lm(data=df_all, formula = form) %>% summary()
    mm2 <- cbind(broom::glance(mm),
                 data.frame(intercept=mm$coefficients[1],
                            slope=mm$coefficients[2],
                            taxon=x,
                            correlation = cor(df_all[, quantvar], df_all[, x], method = corrmethod, use="complete.obs"),
                            correlation_pvalue = cor.test(df_all[, quantvar], df_all[, x], method = corrmethod, use="complete.obs")$p.value
                 ))
    return(mm2)
  }
  )) %>% mutate(taxon = gsub("_", " ", taxon) %>%
                  gsub("sp ", "sp. ", .) %>%
                  #paste0("italic('", ., "')") %>%
                  factor(levels=vars2factor),
                color = ifelse(correlation_pvalue > plim_plot, C_NS,
                               ifelse(slope <0 , C_CTRL, C_CASE))
  )

  g2 <- ggplot(mods, aes(y=correlation, x=taxon))+
    facet_grid( taxon ~ ., scales="free")+
    geom_col(fill=mods$color)+
    coord_flip()+
    theme_minimal() +
    geom_hline(yintercept = 0, col=C_NS, linetype=2)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    #theme(strip.text.y = element_text(size = 10,
    #                                  colour = "black", angle = 0, face = "italic")) +
    theme(strip.text.y = element_blank())+
    theme(axis.text.x = element_text(size = 8,
                                     colour = "black", angle = 0,
                                     face = "plain")) +
    theme(axis.text.y = element_text(size = 12,
                                     colour = "black", angle = 0,
                                     vjust = 0.5,
                                     face = "italic"))+
    ylab(paste0(make_clean_names(corrmethod, case="big_camel"), " corr. ", quantvar_name))+
    thin_barplot_lines

  ggsave(filename = paste0(outdir, "/", name,'_',quantvar ,"_barplot.pdf"), g2, width = w, height = h)

  g3 <- g1 + theme(strip.text.y = element_blank()) + theme(legend.position = "none")
  cw <- cowplot::plot_grid(plotlist=list(g2, g3), nrow = 1, rel_widths = c(2,1))
  pdf(paste0(outdir, name, "_merged.pdf"), width = w, height = h)
  print(cw)
  dev.off()
  write_tsv(mods, file = paste0(outdir, name,'_',quantvar ,"_correlations.tsv"))
  return(list(box=g1, bars=g2, cw=cw, correlations=mods))
}
