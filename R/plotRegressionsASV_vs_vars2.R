#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param df3 PARAM_DESCRIPTION
#' @param corrdf_annot PARAM_DESCRIPTION
#' @param vars2cor PARAM_DESCRIPTION
#' @param groupvar PARAM_DESCRIPTION, Default: 'Psoriasis'
#' @param group_levels PARAM_DESCRIPTION, Default: c("no", "yes")
#' @param top_n PARAM_DESCRIPTION, Default: 10
#' @param groupnames PARAM_DESCRIPTION, Default: c("Control", "Psoriasis")
#' @param ylabel PARAM_DESCRIPTION, Default: 'log(pg/mL)'
#' @param outdir PARAM_DESCRIPTION, Default: ''
#' @param name PARAM_DESCRIPTION, Default: 'regr.pdf'
#' @param w PARAM_DESCRIPTION, Default: 7
#' @param h PARAM_DESCRIPTION, Default: 14
#' @param opt PARAM_DESCRIPTION, Default: list()
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{desc}}
#'  \code{\link[cowplot]{plot_grid}}
#' @rdname plotRegressionsASV_vs_vars2
#' @export 
#' @importFrom dplyr mutate desc
#' @importFrom cowplot plot_grid
plotRegressionsASV_vs_vars2<- function(df3, corrdf_annot, vars2cor, groupvar="Psoriasis",
                                       group_levels = c("no", "yes"),
                                       top_n=10,
                                       groupnames = c("Control", "Psoriasis"),
                                       ylabel = "log(pg/mL)",
                                       outdir = "", name = "regr.pdf", w=7, h=14, opt=list()){
  library(ggpmisc)
  df4 <- df3 %>% gather("ASV", "Abundance", matches("^ASV"))
  df4[, groupvar] <- factor(df4[, groupvar], levels=group_levels)

  corrdf_annot <- corrdf_annot %>%
    dplyr::mutate(var1_g = paste0(sapply(Genus, FUN=function(x){
      gsub("\\[|\\]", "", strsplit(x, " ")[[1]][1], perl=T)
    }), ":", var1))
  df4$ASV_g <- corrdf_annot$var1_g[match(df4$ASV, corrdf_annot$var1)]
  plotlist <- list()
  for(var2cor in vars2cor){
    best_asvs <- corrdf_annot %>% filter(var2==var2cor) %>%
      arrange(dplyr::desc(PSO_spearman_cor + CTRL_spearman_cor)) %>%
      head(top_n) %>% pull(var1_g)
    df4b <- df4 %>%  filter(ASV_g %in% best_asvs) %>% dplyr::mutate(ASV_g = factor(ASV_g, levels=best_asvs))
    df5 <- df4b %>% filter(df4b[, metaname]==group_levels[1])
    df6 <-  df4b %>% filter(df4b[, metaname]==group_levels[2])

    gno<-ggplot(df5, aes_string(x="Abundance", y=paste0("`", var2cor, "`"))) +
      facet_wrap( ASV_g ~ . , scales="free", ncol=1) +
      geom_point(col="dodgerblue3", alpha=1, size=0.5) +
      geom_smooth(method="lm", col="dodgerblue3") +
      stat_poly_eq(use_label(c("R2", "p", "n")), #c("eq", "R2", "f", "p", "n")
                   method="lm", small.p=T, small.r=F, label.y=0.99,
                   col="dodgerblue3") +
      theme_pubr() +
      #scale_color_lancet() +
      ylab(ylabel) +
      xlab("CLR-transformed Abundance")+
      ggtitle(groupnames[1])


    gyes<-ggplot(df6, aes_string(x="Abundance", y=paste0("`", var2cor, "`")))+
      facet_wrap( ASV_g ~ . , scales="free", ncol=1) +
      geom_smooth(method="lm", col="firebrick3") +
      geom_point( col="firebrick3", alpha=1, size=0.5) +
      stat_poly_eq(use_label(c("R2", "p", "n")), #c("eq", "R2", "f", "p", "n")
                   method="lm", small.p=T, small.r=F, label.y=0.99,
                   col="firebrick3") +
      theme_pubr() +
      #scale_color_lancet() +
      ylab(ylabel)+
      xlab("CLR-transformed Abundance") +
      ggtitle(groupnames[2])

    plotlist[[var2cor]] <- cowplot::plot_grid(gno, gyes)

  }
  WriteManyPlots(plotlist, name, outdir, w=w, h=h, separate = F, opt=opt)
  return(plotlist)
}
