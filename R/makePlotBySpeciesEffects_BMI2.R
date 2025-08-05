#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param bacnames PARAM_DESCRIPTION
#' @param df_all PARAM_DESCRIPTION
#' @param res_final PARAM_DESCRIPTION
#' @param outdir PARAM_DESCRIPTION
#' @param name PARAM_DESCRIPTION
#' @param w PARAM_DESCRIPTION, Default: 8
#' @param h PARAM_DESCRIPTION, Default: 10
#' @param plim_plot PARAM_DESCRIPTION, Default: 0.05
#' @param fix_limits PARAM_DESCRIPTION, Default: FALSE
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
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{summarise}}, \code{\link[dplyr]{filter}}
#'  \code{\link[cowplot]{plot_grid}}
#' @rdname makePlotBySpeciesEffects_BMI2
#' @export 
#' @importFrom assertthat assert_that
#' @importFrom dplyr select mutate summarise filter
#' @importFrom cowplot plot_grid
makePlotBySpeciesEffects_BMI2 <- function(bacnames, df_all, res_final, outdir, name,
                                          w=8, h=10, plim_plot=0.05, fix_limits=FALSE){

  titlesize = ifelse(h <=10, 12, h+2)
  assertthat::assert_that(all(bacnames %in% names(df_all)))
  vars2factor <- bacnames %>% gsub("_", " ", .) %>%
    gsub("sp ", "sp. ", .) #%>%
  # paste0("italic('", ., "')")

  vars2factor <- factor(vars2factor, levels = vars2factor)

  df2box <- df_all %>% dplyr::select(all_of(c(bacnames, "Condition", mediator_name))) %>%
    gather(key="taxon", value="abundance", bacnames) %>%
    dplyr::mutate(taxon = gsub("_", " ", taxon) %>%
                    gsub("[\\[\\]]", "", ., perl=T) %>%
                    gsub("sp ", "sp. ", .) %>%
                    #paste0("italic('", ., "')") %>%
                    factor(levels=vars2factor)) %>%
    dplyr::mutate(Condition = ifelse(Condition == "Control", "C", "D"),
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
    theme(axis.text.x = element_text(size = ifelse(h<=8, 8, 10),
                                     colour = "black", angle = 0,
                                     face = "plain"))+
    theme(axis.text.y = element_text(size = ifelse(h<=8, 8, 10),
                                     colour = "black", angle = 0,
                                     face = "plain")) +
    #theme(axis.text.y = element_blank())+
    theme(axis.title = element_text(size = titlesize)) +
    ylab("Taxon abundance")
  ggsave(filename = paste0(outdir, "/", name, ".pdf"), g1, width = w, height = h)


  res_mod <- res_final %>% dplyr::mutate(
    param = gsub("[0-9]+", "", param),
    Effect = ifelse(param=="a*b", "Indirect (a*b)",
                    ifelse(param=="cp", "Direct (c')",
                           ifelse(param=="cp+a*b", "Total (c' + a*b)",
                                  ifelse(param=="a", "Effects on BMI (a)", "Other")))),
  ) %>% dplyr::filter(Effect!="Other") %>%
    dplyr::mutate(taxon = gsub("_", " ", x_labels) %>%
                    gsub("[\\[\\]]", "", ., perl=T) %>%
                    gsub("sp ", "sp. ", .) %>%
                    #paste0("italic('", ., "')") %>%
                    factor(levels=vars2factor),
                  color = ifelse(p.value <= plim_plot, ifelse(Estimate< 0, "Neg", "Pos") , "NS"),
                  color = factor(color, levels = c("Neg", "Pos", "NS")),
                  color2 = ifelse(p.value <= plim_plot, ifelse(Estimate< 0, C_POS_EFF, C_NEG_EFF) , C_NS),
                  Estimate = 10*Estimate
    )
  val_lims = c(min(res_mod$Estimate), max(res_mod$Estimate))
  plots <- map(c("Total (c' + a*b)", "Direct (c')", "Effects on BMI (a)", "Indirect (a*b)"), \(xx){
    ptab <- res_mod %>% dplyr::filter(Effect == xx) %>%
      dplyr::mutate(taxon = factor(taxon, levels = rev(vars2factor)))
    g2 <- ggplot(ptab, aes(y=Estimate, x=taxon))+ #, fill=color
      geom_col(fill=ptab$color2)+
      coord_flip()+
      scale_fill_manual(values = c(C_CTRL, C_CASE, C_NS))+
      theme_minimal() +
      geom_hline(yintercept = 0, col=C_NS, linetype=2)+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      #theme(strip.text.y = element_text(size = 10,
      #                                  colour = "black", angle = 0, face = "italic")) +
      theme(strip.text.y = element_blank())+
      theme(axis.text.x = element_text(size =  ifelse(h<=10, 8, 10), # 8,
                                       colour = "black", angle = 0,
                                       face = "plain")) +
      theme(axis.text.y = element_text(size = titlesize, # 12,
                                       colour = "black", angle = 0,
                                       vjust = 0.5,
                                       face = "italic"))+
      theme(axis.title = element_text(size = titlesize))+
      ylab("")+
      thin_barplot_lines
    if(fix_limits){
      g2 <- g2 + ylim(val_lims)
    }
    return(g2)
  })

  ga = plots[[1]] + theme(legend.position = "none") + ylab("Total effects (c' + a*b)")
  gb = plots[[2]] + theme(legend.position = "none")  + theme(axis.text.y = element_blank()) + xlab("") + ylab("Direct Effects (c')")
  gc = plots[[3]] + theme(legend.position = "none")  + theme(axis.text.y = element_blank())+ xlab("") + ylab("Effects on BMI (a)")
  gd = plots[[4]] + theme(legend.position = "none")  + theme(axis.text.y = element_blank())+ xlab("") + ylab("Indirect Effects (a*b)")
  ge <- g1 + theme(legend.position = "none") + theme(strip.text.y = element_blank())
  cw <- cowplot::plot_grid(plotlist=list(ga, gb, gc, gd, ge), nrow = 1, rel_widths = c(2.5,1, 1, 1,1))
  pdf(paste0(outdir, name, "_BarplotEffectsmerged.pdf"), width = w, height = h)
  print(cw)
  dev.off()
  write_tsv(res_mod, file = paste0(outdir, name ,"_BarplotEffectsmerged.tsv"))
  return(list(box=g1, bars=plots, cw=cw, tab=res_mod))
}
