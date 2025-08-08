#' @title Plot Regressions of ASV Abundances vs Variables
#' @description Generates regression plots between CLR-transformed ASV abundances and selected variables, stratified by groups.
#'
#' For each variable in `vars2cor`, the function selects the top `top_n` ASVs based on combined Spearman correlations
#' from `corrdf_annot`. It then plots regression lines and statistics separately for each group defined in `groupvar`.
#'
#' @param df3 A data.frame containing the ASV abundance data and metadata variables.
#' @param corrdf_annot A data.frame with correlation annotations including columns `Genus`, `var1`, `var2`, and group-specific Spearman correlations (`PSO_spearman_cor`, `CTRL_spearman_cor`).
#' @param vars2cor Character vector of variable names to correlate against ASVs.
#' @param groupvar Character string naming the grouping variable in `df3`. Default is `"Psoriasis"`.
#' @param group_levels Character vector defining the factor levels of `groupvar`. Default is `c("no", "yes")`.
#' @param top_n Integer specifying the number of top correlated ASVs to plot per variable. Default is `10`.
#' @param groupnames Character vector of length two for naming the groups in plots. Default is `c("Control", "Psoriasis")`.
#' @param ylabel Character string for the y-axis label in plots. Default is `"log(pg/mL)"`.
#' @param outdir Directory path where to save the output plots. Default is current directory `""`.
#' @param name File name for the combined plot output. Default is `"regr.pdf"`.
#' @param w Numeric width of the output plot in inches. Default is `7`.
#' @param h Numeric height of the output plot in inches. Default is `14`.
#' @param opt List of additional options passed to the plot saving function. Default is empty list `list()`.
#'
#' @return A named list of ggplot objects, one per variable in `vars2cor`.
#'
#' @details
#' The function transforms the data to long format for plotting, matches ASVs to genus-level annotations,
#' selects top correlated ASVs per variable, and generates separate regression plots per group with
#' linear model fits and polynomial equation statistics.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Assuming df3, corrdf_annot, and vars2cor are pre-loaded:
#'   plots <- plotRegressionsASV_vs_vars2(df3, corrdf_annot, vars2cor = c("IL6", "TNFa"))
#' }
#' }
#'
#' @seealso
#' \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{desc}}, \code{\link[cowplot]{plot_grid}}
#'
#' @rdname plotRegressionsASV_vs_vars2
#' @export
#' @importFrom dplyr mutate desc
#' @importFrom cowplot plot_grid
#' @importFrom ggpmisc stat_poly_eq
#' @importFrom ggpubr theme_pubr
plotRegressionsASV_vs_vars2<- function(df3, corrdf_annot, vars2cor, groupvar="Psoriasis",
                                       group_levels = c("no", "yes"),
                                       top_n=10,
                                       groupnames = c("Control", "Psoriasis"),
                                       ylabel = "log(pg/mL)",
                                       outdir = "", name = "regr.pdf", w=7, h=14, opt=list()){
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
