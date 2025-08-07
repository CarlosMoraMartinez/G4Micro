
#' @title Compare Model Accuracies Across Multiple Conditions for a Single Phyloseq Object
#' @description This function compares classification model accuracies across several conditions
#' (e.g., groups of selected features or experimental groups) for a given phyloseq object.
#' It aggregates model summaries from a list of results, formats them, and generates a line plot of accuracies.
#' A TSV file with all model summaries is also saved.
#'
#' @param phname Name of the phyloseq object (used to subset \code{all_model_results}).
#' @param condnames Character vector with the names of different feature selection conditions or experimental groups to compare.
#' @param all_model_results A nested list of model results. The top-level list is indexed by phyloseq object names,
#' and the second level is indexed by condition names (e.g., "padj_taxa_res", "raw_taxa_res"), each containing a \code{modummary} data frame.
#' @param name Prefix used to name the output plot file.
#' @param opt A list containing output options. Must include the \code{opt$out} directory path.
#' @param w Width of the output PDF file in inches. Default is 8.
#' @param h Height of the output PDF file in inches. Default is 12.
#'
#' @return A ggplot2 object representing the accuracy comparison plot. Additionally, a TSV file and a PDF plot are saved to disk.
#'
#' @details
#' The function loops through the specified conditions, extracts the model accuracy summaries for each,
#' appends metadata for plotting, merges all results into a single table, and generates a line plot
#' showing accuracy (Accuracy_l1out) across different models and conditions.
#'
#' The y-axis shows accuracy, the x-axis shows model type, and lines/points are grouped by condition.
#'
#' Output files are saved in a subdirectory named after the \code{phname}, within the path specified by \code{opt$out}.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{arrange}}, \code{\link[dplyr]{select}}, \code{\link[dplyr]{summarise}}
#'  \code{\link[ggpubr]{theme_pubr}}
#' @rdname makeLinePlotComparingSamePhobjModels_Cov
#' @export
#' @importFrom dplyr mutate arrange select summarise
#' @importFrom ggpubr theme_pubr
#' @importFrom ggsci scale_color_cosmic scale_fill_cosmic
makeLinePlotComparingSamePhobjModels_Cov<- function(phname, condnames,
                                                    all_model_results,
                                                    name, opt, w=8, h=12){
  outdir <- paste0(opt$out, phname)
  if(!file.exists(outdir)) dir.create(outdir)
  all_sig_tables <- map(condnames, \(name){
    all_model_results[[phname]][[name]]$modummary %>%
      dplyr::mutate(taxa_group=name, input_data=name)
  }) %>% bind_rows() %>% dplyr::arrange(desc(Accuracy_l1out)) %>%
    dplyr::select(input_data, model, taxa_group, everything())

  model_levels <- all_sig_tables %>% group_by(model) %>%
    dplyr::summarise(acc=max(Accuracy_l1out)) %>%
    dplyr::arrange(desc(acc)) %>% pull(model)
  all_sig_tables <- all_sig_tables %>% dplyr::mutate(model = factor(model, levels=model_levels))
  write_tsv(all_sig_tables, file = paste0(outdir, "/all_model_summaries.tsv"))

  (g1 <- ggplot(all_sig_tables, aes(x=model,
                                    y=Accuracy_l1out,
                                    col=input_data,
                                    fill=input_data,
                                    group=input_data))+
      #facet_grid(taxa_group~.)+
      geom_point()+
      geom_line()+
      scale_color_cosmic()+
      scale_fill_cosmic() +
      ggpubr::theme_pubr() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  )
  ggsave(paste0(outdir, "/all_model_accuracy_", name, ".pdf"), g1,
         width = w, height = h)
  return(g1)

}
