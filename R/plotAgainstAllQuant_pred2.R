
plotAgainstAllQuant_pred2<- function(alldf, groupvars, name, outdir, col_var, opt, names_cyt=c()){

  plotlist <- list()
  for(grouping_var in groupvars){
    df <- alldf[, c(grouping_var, names_cyt, col_var)] %>% tidyr::gather(IL, pg_mL, names_cyt)
    oname <- paste0(outdir, name)
    plotlist[[grouping_var]] <- plotRegr(df, variable="pg_mL",
                                         x_var = grouping_var,
                                         wrap_var = "IL",
                                         col_var ="Psoriasis",
                                         fname=oname,
                                         write = F)
  }
  WriteManyPlots(plotlist, name = name, outdir =outdir, w = 12, h = 12, separate = F, opt)
  return(plotlist)
}
