
plotBarplot_pred2 <- function(df, variable, grouping_var, wrap_var,
                        fname,
                        signif_levels=c("***"=0.001, "**"=0.01, "*"=0.05, "ns"=1.1),
                        num_comparisons = 1,
                        ylabel = "pg/mL",
                        correct_pvalues = TRUE, write=TRUE){
  if(correct_pvalues){
    signif_levels_bonferroni <- c(signif_levels[1:3]/num_comparisons, signif_levels[4])
  }else{
    signif_levels_bonferroni <- signif_levels
  }
  wrap_formula <- paste0(". ~ ",  wrap_var) %>% as.formula
  df[, grouping_var] <- as.factor(df[, grouping_var])
  comp <- combn(unique(as.character(df[, grouping_var])), 2, simplify = F)

  g1 <-  ggplot(df, aes_string(x=grouping_var, y=variable, fill=grouping_var, col=grouping_var)) +
    facet_wrap(wrap_formula, scales = "free", nrow=1) +
    #geom_bar(stat="identity") +
    stat_summary(fun.data=mean_sdl, geom="bar", width=0.8, alpha=0.8) +
    stat_summary(fun.data=mean_cl_boot, geom="errorbar", width=0.3, col="black") +
    geom_point()+
    #scale_color_manual(values = c("#ffafcc", "#90DBF4")) +
    #scale_fill_manual(values = c("#ffafcc", "#90DBF4")) +

    labs(title = variable, x = '') +
    #mytheme +
    theme_pubclean() +
    geom_signif(test="wilcox.test", na.rm=T, comparisons = comp,
                step_increase=0.03,
                tip_length = 0.01,
                map_signif_level=signif_levels_bonferroni,
                vjust=0.4,
                color = "black"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_blank()) +
    mytheme +
    ggtitle("Cytokine levels per group") +
    ylab(ylabel)
  if(length(unique(df[, grouping_var] )) == 2){
    g1 <- g1  +
      scale_color_manual(values =c("firebrick4", "dodgerblue4")) +
      scale_fill_manual(values =c("firebrick2", "dodgerblue2"))
  }else{
    g1 <- g1  +
      scale_color_lancet() +
      scale_fill_lancet()
  }
  if(write){
    ggsave(filename=fname, g1, width = 12, height = 6)
  }
  return(g1)
  #theme(axis.text.x = element_text(angle = 360, hjust = 0.5, size = 10))
} #Plots qualitative variables
