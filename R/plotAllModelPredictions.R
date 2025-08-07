#' @title Plot All Model Predictions on PCA Components for a Phenotype
#' @description
#' Generates and saves combined plots of predictions from multiple models
#' for a given phenotype/dataset over PCA components.
#' The models are ordered by accuracy and plotted using \code{plotPrediction}.
#' The combined plot is saved as a PDF.
#'
#' @param phname Character. Name of the phenotype or dataset subset to analyze.
#' @param all_model_results List. Nested list containing PCA results, model predictions, and metadata.
#' @param opt List. Options list containing at least the output directory path as \code{opt$out}.
#' @param get_pcnames_from Character. Name of the sublist in \code{all_model_results} to get PCA variable names and models from. Default is "padj_taxa_res".
#' @param pca_name Character. Name of the PCA results list in \code{all_model_results}. Default is "padj_taxa_pcas".
#' @param varname Character. Name of the metadata variable to use for grouping/condition. Default is "Condition".
#' @param w Numeric. Width of the saved combined plot PDF in inches. Default is 16.
#' @param h Numeric. Height of the saved combined plot PDF in inches. Default is 16.
#' @param pred_mode Character. Prediction mode to use: "l1o" (leave-one-out) or other mode. Default is "l1o".
#'
#' @return A list with two elements:
#' \item{cw}{A \code{ggplot} object representing the combined grid of all model prediction plots.}
#' \item{plots}{A named list of individual \code{ggplot} objects or empty lists (if a plot failed).}
#'
#' @details
#' The function:
#' \itemize{
#'   \item Creates an output directory if it does not exist.
#'   \item Extracts all model names and orders them by their accuracy.
#'   \item Calls \code{plotPrediction} for each model, safely capturing errors.
#'   \item Combines successful plots into a grid using \code{cowplot::plot_grid}.
#'   \item Saves the combined plot as a PDF file.
#'   \item Returns the combined plot object and the list of individual plots.
#' }
#' If no plots are successfully generated, a message is printed instead.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[cowplot]{plot_grid}}
#'  \code{\link{plotPrediction}}
#' @rdname plotAllModelPredictions
#' @export
#' @importFrom cowplot plot_grid
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
