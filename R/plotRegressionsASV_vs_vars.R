#' @title Plot Linear Regressions of ASVs vs Variables by Group
#' @description
#' Creates scatterplots with linear regression fits of ASVs against variables of interest,
#' faceted by variable, and separated by group.
#'
#' @param df3 Data frame containing ASVs and variables.
#' @param asvs Character vector of ASV column names to plot as response variables.
#' @param vars2cor Character vector of variable column names to use as predictors.
#' @param groupvar Name of the grouping variable in \code{df3}, default is "Psoriasis".
#' @param group_levels Levels of the grouping variable, default is c("no", "yes").
#' @param groupnames Names to use in plot titles for each group, default is c("Control", "Psoriasis").
#' @param xlabel Label for the x-axis, default is "log(pg/mL)".
#' @param outdir Output directory to save plots, default is "" (current directory).
#' @param name Filename for the saved plot PDF, default is "regr.pdf".
#' @param w Width of the output plot(s) in inches, default is 7.
#' @param h Height of the output plot(s) in inches, default is 14.
#' @param opt List of additional options to pass to the plotting/saving function.
#' @return A named list of ggplot objects, one for each ASV.
#' @details
#' The function creates separate regression plots for each ASV against all variables in \code{vars2cor},
#' faceted by variable, and plots them separately for the two groups defined in \code{groupvar}.
#' The regression line includes equation, R², p-value, and sample size using \pkg{ggpmisc}.
#' Plots for each group are combined side-by-side using \pkg{cowplot}.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   plotlist <- plotRegressionsASV_vs_vars(df, asvs = c("ASV1", "ASV2"), vars2cor = c("Var1", "Var2"))
#' }
#' }
#' @seealso
#' \code{\link[cowplot]{plot_grid}}, \code{\link[ggpmisc]{stat_poly_eq}}
#' @rdname plotRegressionsASV_vs_vars
#' @export
#' @importFrom cowplot plot_grid
#' @importFrom ggpmisc stat_poly_eq
#' @importFrom tidyr gather
#' @importFrom dplyr filter
#' @importFrom ggpubr theme_pubr
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
