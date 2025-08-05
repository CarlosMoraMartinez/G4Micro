#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param df3 PARAM_DESCRIPTION
#' @param asvs PARAM_DESCRIPTION
#' @param vars2cor PARAM_DESCRIPTION
#' @param groupvar PARAM_DESCRIPTION, Default: 'Psoriasis'
#' @param group_levels PARAM_DESCRIPTION, Default: c("no", "yes")
#' @param groupnames PARAM_DESCRIPTION, Default: c("Control", "Psoriasis")
#' @param xlabel PARAM_DESCRIPTION, Default: 'log(pg/mL)'
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
#'  \code{\link[cowplot]{plot_grid}}
#' @rdname plotRegressionsASV_vs_vars
#' @export 
#' @importFrom cowplot plot_grid
plotRegressionsASV_vs_vars<- function(df3, asvs, vars2cor, groupvar="Psoriasis",
                                      group_levels = c("no", "yes"),
                                      groupnames = c("Control", "Psoriasis"),
                                      xlabel = "log(pg/mL)",
                                      outdir = "", name = "regr.pdf", w=7, h=14, opt=list()){
  library(ggpmisc)
  df4 <- df3 %>% gather("IL", "log_val", vars2cor)
  df4[, groupvar] <- factor(df4[, groupvar], levels=group_levels)
  df4[, "IL"] <- factor(df4$IL, levels=vars2cor)

  df5 <- df4 %>% filter(df4[, metaname]==group_levels[1])
  df6 <-  df4 %>% filter(df4[, metaname]==group_levels[2])

  plotlist <- list()
  for(ASV in asvs){
    gno<-ggplot(df5, aes_string(x="log_val", y=ASV)) +
      facet_wrap( IL ~ . , scales="free", ncol=1) +
      geom_point(col="dodgerblue3", alpha=1, size=0.5) +
      geom_smooth(method="lm", col="dodgerblue3") +
      stat_poly_eq(use_label(c("R2", "p", "n")), #c("eq", "R2", "f", "p", "n")
                   method="lm", small.p=T, small.r=F, label.y=0.99,
                   col="dodgerblue3") +
      theme_pubr() +
      #scale_color_lancet() +
      ylab(ASV) +
      xlab(xlabel)+
      ggtitle(groupnames[1])


    gyes<-ggplot(df6, aes_string(x="log_val", y=ASV))+
      facet_wrap( IL ~ . , scales="free", ncol=1) +
      geom_smooth(method="lm", col="firebrick3") +
      geom_point( col="firebrick3", alpha=1, size=0.5) +
      stat_poly_eq(use_label(c("R2", "p", "n")), #c("eq", "R2", "f", "p", "n")
                   method="lm", small.p=T, small.r=F, label.y=0.99,
                   col="firebrick3") +
      theme_pubr() +
      #scale_color_lancet() +
      ylab(ASV)+
      xlab(xlabel) +
      ggtitle(groupnames[2])

    plotlist[[ASV]] <- cowplot::plot_grid(gno, gyes)

  }
  WriteManyPlots(plotlist, name, outdir, w=6, h=12, separate = F, opt=opt)
  return(plotlist)
}
