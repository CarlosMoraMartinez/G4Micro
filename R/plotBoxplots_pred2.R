
plotBoxplots_pred2 <- function(df, variable, grouping_var, wrap_var,
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
  g1 <-  ggplot(df, aes_string(x=grouping_var, y=variable, col=grouping_var)) +
    facet_wrap(wrap_formula, scales = "free") +
    geom_boxplot(alpha = 0.7, width=0.5, outlier.alpha=0) +
    geom_jitter(alpha=0.5)+
    #scale_color_manual(values = c("#ffafcc", "#90DBF4")) +
    #scale_fill_manual(values = c("#ffafcc", "#90DBF4")) +
    scale_color_lancet() +
    scale_fill_lancet() +
    labs(title = variable, x = '') +
    #mytheme +
    theme_pubclean() +
    ggtitle("") +
    geom_signif(test="wilcox.test", na.rm=T, comparisons = comp,
                step_increase=0.03,
                tip_length = 0.01,
                map_signif_level=signif_levels_bonferroni,
                vjust=0.4,
                color = "black"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_blank()) +
    ylab(ylabel)

  if(write){
    ggsave(filename=fname, g1, width = 6, height = 8)
  }
  return(g1)
  #theme(axis.text.x = element_text(angle = 360, hjust = 0.5, size = 10))
} #Plots qualitative variables

