
plotAgainstAllFactors_Barplot_pred2<- function(alldf, groupvars, name, outdir, opt, correct_bonferroni=F, names_cyt=c()){

  plotlist <- list()
  for(grouping_var in groupvars){
    if(length(unique(alldf[, grouping_var])) < 2 | min(table(alldf[, grouping_var])) <2) next
    df <- alldf[, c(grouping_var, names_cyt)] %>% gather(IL, pg_mL, names_cyt)
    oname <- paste0(outdir, name)
    plotlist[[grouping_var]] <- plotBarplot(df, "pg_mL", grouping_var, "IL",
                                            fname=oname,
                                            signif_levels=c("***"=0.001, "**"=0.01, "*"=0.05, "ns"=1.1),
                                            num_comparisons = nrow(test_res),
                                            ylabel = "pg/mL",
                                            correct_pvalues = correct_bonferroni,
                                            write = F)
  }
  WriteManyPlots(plotlist, name = name, outdir =outdir, w = 12, h = 12, separate = F, opt)
  return(plotlist)
}
