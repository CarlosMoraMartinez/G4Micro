
#' @title Plot Model Predictions on PCA Components
#' @description
#' Generates a scatter plot showing model predictions over the first two principal components (PCs).
#' Points are colored by their true condition, with misclassified points highlighted.
#' The plot includes the accuracy of the model in the title and is saved as a PDF.
#'
#' @param phname Character. Name of the phenotype or dataset subset to analyze.
#' @param mod2plot Character. Name of the model to plot predictions for (e.g., "KNN-K=5").
#' @param all_model_results List. Nested list containing PCA results, model predictions, and metadata.
#' @param opt List. Options list containing at least the output directory path as \code{opt$out}.
#' @param get_pcnames_from Character. Name of the sublist in \code{all_model_results} to get PCA variable names and models from. Default is "padj_taxa_res".
#' @param pca_name Character. Name of the PCA results list in \code{all_model_results}. Default is "padj_taxa_pcas".
#' @param varname Character. Name of the metadata variable to use for grouping/condition. Default is "Condition".
#' @param pred_mode Character. Prediction mode to use: "l1o" (leave-one-out) or other mode (e.g., "no_l1o"). Default is "l1o".
#' @param w Numeric. Width of the saved plot in inches. Default is 6.
#' @param h Numeric. Height of the saved plot in inches. Default is 4.
#'
#' @return A \code{ggplot} object showing predicted vs actual classes on the first two PCs.
#' The plot is saved as a PDF file in the output directory.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Retrieves model predictions and matches them to sample IDs.
#'   \item Extracts the first two principal components for samples.
#'   \item Highlights correct vs. incorrect predictions by point size and layering.
#'   \item Calculates a confusion matrix and shows the overall accuracy in the plot title.
#'   \item Saves the plot to the specified output folder.
#' }
#'
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
