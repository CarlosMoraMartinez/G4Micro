plotAllModelPredictions <- function(phname, all_model_results, opt,
                                    get_pcnames_from="padj_taxa_res",
                                    pca_name="padj_taxa_pcas",
                                    varname="Condition", w=16, h=16,
                                    pred_mode="l1o"){
  outdir <- paste0(opt$out, "/", phname)
  if(!file.exists(outdir)) dir.create(outdir)
  modelnames <- names(all_model_results[[phname]][[get_pcnames_from]]$models)
  modnames_order <- map_vec(modelnames, \(x)all_model_results[[phname]][[get_pcnames_from]]$models[[x]]$confmat$overall["Accuracy"]) %>% order(decreasing = T)
  modelnames <- modelnames[modnames_order]
  modplots <- map(modelnames,
                  \(modpred){
                    tryCatch(plotPrediction(phname, modpred, all_model_results, opt,
                                            get_pcnames_from, pca_name, varname, pred_mode), error=\(x)list())
                  })
  names(modplots) <- names(all_model_results[[phname]][[get_pcnames_from]]$models)
  modplots2 <- list()
  for(mn in names(modplots)){
    m <- modplots[[mn]]
    if(class(m)[1] == "list" & length(m)==0) next
    modplots2[[mn]] <- m
  }
  if(length(modplots2) > 1){
    cw <- cowplot::plot_grid(plotlist = modplots2)
    pdf(paste0(outdir, "/all_model_predictions.pdf"), width = w, height = h)
    print(cw)
    dev.off()
  }else{
    cat("FAILED prediction plots for", phname)
  }
  return(list(cw=cw, plots=modplots))
}
