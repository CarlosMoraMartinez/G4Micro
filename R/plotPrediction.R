
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param phname PARAM_DESCRIPTION
#' @param mod2plot PARAM_DESCRIPTION
#' @param all_model_results PARAM_DESCRIPTION
#' @param opt PARAM_DESCRIPTION
#' @param get_pcnames_from PARAM_DESCRIPTION, Default: 'padj_taxa_res'
#' @param pca_name PARAM_DESCRIPTION, Default: 'padj_taxa_pcas'
#' @param varname PARAM_DESCRIPTION, Default: 'Condition'
#' @param pred_mode PARAM_DESCRIPTION, Default: 'l1o'
#' @param w PARAM_DESCRIPTION, Default: 6
#' @param h PARAM_DESCRIPTION, Default: 4
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{mutate}}
#'  \code{\link[caret]{confusionMatrix}}
#' @rdname plotPrediction
#' @export 
#' @importFrom dplyr mutate
#' @importFrom caret confusionMatrix
plotPrediction<-function(phname, mod2plot, all_model_results, opt,
                         get_pcnames_from="padj_taxa_res",
                         pca_name="padj_taxa_pcas",
                         varname="Condition", pred_mode="l1o", w=6, h=4){
  outdir <- paste0(opt$out, "/", phname)
  if(!dir.exists(outdir)) dir.create(outdir)
  if(pred_mode=="l1o"){
    predictions <- all_model_results[[phname]][[get_pcnames_from]]$models[[mod2plot]]$preds
  }else{
    predictions <- all_model_results[[phname]][[get_pcnames_from]]$models[[mod2plot]]$preds_no_l1o
  }
  snames <- all_model_results[[phname]][[pca_name]][[varname]]$pca$x %>% rownames
  PCs <- all_model_results[[phname]][[get_pcnames_from]]$varnames
  pc_order <- gsub("PC", "", PCs) %>% as.numeric %>% order
  PCs <- PCs[pc_order]
  PCs_newnames <- getPCnamesFromAllresults(phname, all_model_results, get_pcnames_from, pca_name, varname) %>% names %>% sort
  PCs_newnames <- PCs_newnames[pc_order]

  df <- pcBoxplots[[phname]]$tab %>%
    dplyr::mutate(Predicted = predictions[match(sampleID, snames)]) %>%
    spread(key=PC, value=score)
  if(all(unique(df$Condition) %in% c("Control", "Depression", "Depr."))){
    df <- df %>% dplyr::mutate(Condition = fct_recode(Condition, "Control"="Control", "Depression"="Depr."),
                               Good = ifelse(Condition==Predicted, TRUE, FALSE),
                               size2 = ifelse(Good, 0, 1))
  }else{
    df <- df %>% dplyr::mutate(Predicted = gsub("Depression", "D", as.character(Predicted)),
                               Predicted = gsub("Control", "C", as.character(Predicted)),
                               Good = ifelse(Condition==Predicted, TRUE, FALSE),
                               size2 = ifelse(Good, 0, 1))
  }
  confmat <- caret::confusionMatrix(factor(df$Condition), factor(df$Predicted))

  gm <- ggplot(df, aes(x=!!sym(PCs_newnames[1]), y =!!sym(PCs_newnames[2]), col=Condition))+
    geom_point(size=2)+
    geom_point(col="black", alpha=df$size2, size=0.6)+
    mytheme +
    ggtitle(paste0(mod2plot, ", Acc=", as.character(round(confmat$overall["Accuracy"], 3))))
  ggsave(filename = paste0(outdir, "/", mod2plot, "_predictionPlot.pdf"), width = w, height = h)
  gm <- gm + theme(legend.position = 'none')
  return(gm)
}
