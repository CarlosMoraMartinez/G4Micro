
#' @title Generate Barplots of Principal Components and Log Fold Change
#' @description
#' Creates barplots comparing PCA rotation values of taxa with their corresponding log2 fold change and adjusted p-values.
#' The function adjusts PCA rotation scores by factors derived from boxplot results to align with log fold change directions,
#' merges these with differential abundance analysis (DAA) results, and plots them faceted by variable.
#' It saves the plots and processed data to disk.
#'
#' @param phname Character. Name of the phenotype or dataset subset to analyze.
#' @param all_model_results List. Nested list containing PCA results, model predictions, and metadata.
#' @param pcBoxplots List. List of outputs from \code{\link{makePCsBoxplot}} containing PCA factors for sign adjustment.
#' @param daa_all List. List containing differential abundance analysis results with log2 fold change and p-values.
#' @param opt List. Options list containing at least the output directory path as \code{opt$out}.
#' @param get_pcnames_from Character. Name of the sublist in \code{all_model_results} to get PCA variable names and models from. Default is "padj_taxa_res".
#' @param pca_name Character. Name of the PCA results list in \code{all_model_results}. Default is "padj_taxa_pcas".
#' @param varname Character. Name of the metadata variable to use for grouping/condition. Default is "Condition".
#' @param w Numeric. Width of the saved plot in inches. Default is 8.
#' @param h Numeric. Height of the saved plot in inches. Default is 14.
#'
#' @return A \code{ggplot} object with barplots of PCs, log2 fold changes, and adjusted p-values for taxa.
#' The plot is saved as a PDF and the data used to create the plot is saved as a TSV file in the output directory.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Retrieves PCA rotation scores and adjusts their sign based on factors from \code{pcBoxplots} to align with log fold change direction.
#'   \item Merges adjusted PCA rotations with differential abundance statistics from \code{daa_all}.
#'   \item Creates a long-format dataframe combining PCA components, log2 fold changes, and adjusted p-values (transformed as -10*log10(padj)).
#'   \item Plots faceted barplots with taxa on the y-axis and variable values on the x-axis, colored by significance or condition.
#'   \item Saves the plot and data tables in the specified output folder.
#' }
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{arrange}}
#'  \code{\link[tidyr]{gather}}
#' @rdname makePCBarplot
#' @export
#' @importFrom dplyr select mutate arrange
#' @importFrom tidyr gather
makePCBarplot <- function(phname, all_model_results, pcBoxplots, daa_all, opt,
                          get_pcnames_from="padj_taxa_res",
                          pca_name="padj_taxa_pcas",
                          varname="Condition", w=8, h=14){
  outdir <- paste0(opt$out, "/", phname)
  if(!dir.exists(outdir)) dir.create(outdir)
  predictions <- all_model_results[[phname]][[get_pcnames_from]]$models$`KNN-K=5`$preds_no_l1o
  PCs <- all_model_results[[phname]][[get_pcnames_from]]$varnames
  pc_order<- gsub("PC", "", PCs) %>% as.numeric %>% order
  PCs <- PCs[pc_order]
  PCs_newnames <- getPCnamesFromAllresults(phname, all_model_results, get_pcnames_from, pca_name, varname)
  PCs_newnames <- PCs_newnames[pc_order]
  daatab <- daa_all[[phname]]$resdf
  pcfactors <- pcBoxplots[[phname]]$pcfactors %>% column_to_rownames("PC") %>% dplyr::select(factor)

  df <- all_model_results$remove_tanda2[[pca_name]][[varname]]$pca$rotation %>%
    as.data.frame() %>%
    dplyr::select(all_of(PCs)) %>%
    rownames_to_column("taxon") %>%
    dplyr::mutate(across(all_of(PCs), \(x)x*pcfactors[cur_column(), 1])) ## Multiplicar por factor para que coincida con LFC

  dfmerged <- merge(df, daatab, by="taxon", all.x=T, all.y=F) %>%
    dplyr::arrange(desc(!!sym(PCs[1]))) %>%
    dplyr::mutate(taxon = gsub("_", " ", taxon),
                  taxon = gsub("[\\[\\]]", "", taxon),
                  taxon = factor(taxon, levels = taxon))

  newlevnames <- c( PCs_newnames, "LFC"="log2FoldChangeShrink","-10log(adj. p)"="padj")
  dflong <- dfmerged %>%
    dplyr::select(all_of(c("taxon", "padj", PCs, "log2FoldChangeShrink"))) %>%
    dplyr::mutate(padj = -10*log10(padj)) %>%
    tidyr::gather(key="variable", "value", -taxon) %>%
    dplyr::mutate(color = ifelse(variable == "padj",C_NS, ifelse(value < 0, C_CTRL, C_CASE)),
                  variable = fct_recode(variable, !!!newlevnames),
                  variable = factor(variable, levels = names(newlevnames))
    )

  gbars <- ggplot(dflong, aes(y=value, x=taxon, fill=color)) +
    facet_wrap(~ variable, nrow=1, scales = "free_x")+
    geom_col()+
    coord_flip() +
    theme_classic() +
    scale_fill_manual(values = c(C_CTRL, C_CTRL_LINK2, C_CASE))+
    theme(axis.text.y = element_text(size = 8,
                                     colour = "black", angle = 0,
                                     face = "italic"))+
    theme(axis.text.x = element_text(size = 10,
                                     colour = "black", angle = 0,
                                     face = "plain"))+
    theme(strip.text.x = element_text(size = 14,
                                      colour = "black", angle = 0, face = "plain")) +
    thin_barplot_lines +
    theme(legend.position="none")
  ggsave(filename = paste0(outdir, "/", get_pcnames_from, "_barplots_PCs_and_LFC.pdf"), gbars, width = w, height = h)
  write_tsv(dflong,  paste0(outdir, "/", get_pcnames_from, "_barplots_PCs_and_LFC_data.tsv"))
  return(gbars)
}
