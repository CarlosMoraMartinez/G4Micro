
#' @title Compare Model Accuracies from Two Feature Selection Strategies Across Phyloseq Objects
#' @description This function compares model performance across two feature selection strategies, each performed on features selected from different phyloseq objects.
#' (e.g., padj-selected and raw p-value-selected taxa). It extracts accuracy results from model summaries,
#' merges them, and generates a faceted line plot comparing model performances for each group.
#' The resulting data and figure are saved to disk.
#'
#' @param all_model_results A named list of model results, each containing sublists for different feature selection strategies
#' (e.g., 'padj_taxa_res' and 'praw_taxa_res'), each with a `modummary` data frame.
#' @param opt A list containing output options, specifically `opt$out` for the output directory path.
#' @param models_name1 Character string specifying the name of the first model result type to compare (typically adjusted p-value taxa), Default: 'padj_taxa_res'
#' @param models_name2 Character string specifying the name of the second model result type to compare (typically raw p-value taxa), Default: 'praw_taxa_res'
#'
#' @return A ggplot object containing the faceted accuracy comparison line plot.
#'
#' @details This function assumes that each element in `all_model_results` contains two elements named by `models_name1` and `models_name2`, each of which includes a `modummary` data frame with at least the columns `model` and `Accuracy_l1out`. The plot compares accuracies for each selection method (padj vs raw) across input data sets and models.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{arrange}}, \code{\link[dplyr]{select}}
#'  \code{\link[ggpubr]{theme_pubr}}
#' @rdname makeLinePlotComparingPhobjs
#' @export
#' @importFrom dplyr mutate arrange select
#' @importFrom ggpubr theme_pubr
#' @importFrom ggsci scale_color_cosmic scale_fill_cosmic
makeLinePlotComparingPhobjs <- function(all_model_results, opt, models_name1="padj_taxa_res", models_name2="praw_taxa_res"){

  all_sig_tables <- names(all_model_results) %>% map(\(name){
    a <- all_model_results[[name]][[models_name1]]$modummary %>%
      dplyr::mutate(taxa_group="padj")
    b <- all_model_results[[name]][[models_name2]]$modummary %>%
      dplyr::mutate(taxa_group="praw")
    rbind(a, b) %>% dplyr::mutate(input_data=name)
  }) %>% bind_rows() %>% dplyr::arrange(desc(Accuracy_l1out)) %>%
    dplyr::select(input_data, model, taxa_group, everything())

  write_tsv(all_sig_tables, file = paste0(opt$out, "/all_model_summaries.tsv"))

  (g1 <- ggplot(all_sig_tables, aes(x=model,
                                    y=Accuracy_l1out,
                                    col=input_data,
                                    fill=input_data,
                                    group=input_data))+
      facet_grid(taxa_group~.)+
      geom_point()+
      geom_line()+
      scale_color_cosmic()+
      scale_fill_cosmic() +
      ggpubr::theme_pubr() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  )
  ggsave(paste0(opt$out, "/all_model_accuracy.pdf"), g1,
         width = 8, height = 8)
  return(g1)
}
