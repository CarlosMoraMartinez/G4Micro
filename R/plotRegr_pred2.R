plotRegr_pred2 <- function(df, variable, x_var, wrap_var, col_var,
                     fname,
                     write=TRUE){

  wrap_formula <- paste0(". ~ ",  wrap_var) %>% as.formula
  df[, col_var] <- as.factor(df[, col_var])
  points_col <- ifelse(df[, col_var] == "no", "dodgerblue3", "firebrick3") #levels(df[, col_var])[1]
  g1 <-  ggplot(df, aes_string(x=x_var, y=variable)) +
    facet_wrap(wrap_formula, scales = "free") +
    geom_smooth(method="lm") +
    geom_point(alpha=0.6, col=points_col)+
    stat_poly_eq(use_label(c("R2", "p", "n")), #c("eq", "R2", "f", "p", "n")
                 method="lm", small.p=T, small.r=F, label.y=0.99)+
    #scale_color_manual(values = c("#ffafcc", "#90DBF4")) +
    #scale_fill_manual(values = c("#ffafcc", "#90DBF4")) +
    scale_color_lancet() +
    scale_fill_lancet() +
    labs(title = variable, x = '') +
    #mytheme +
    theme_pubclean() +
    ggtitle("") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text.x = element_blank()) +
    ggtitle(x_var)

  if(write){
    ggsave(filename=fname, g1, width = 6, height = 8)
  }
  return(cowplot::plot_grid(g1))
  #theme(axis.text.x = element_text(angle = 360, hjust = 0.5, size = 10))
} #Plots qualitative variables
