#' @title Create Barplot of Significant Pathways with Log2 Fold Change
#' @description
#' Generates a horizontal barplot visualizing pathways or processes that are significantly
#' differentially abundant based on adjusted p-values and log2 fold changes.
#' The bars are colored by direction (more or less abundant), and optionally the pathway names
#' can be augmented with longer descriptive names if available.
#' The plot is saved as a PDF file to the specified output directory.
#'
#' @param tab A data frame containing differential abundance results with row names representing
#' pathways or features, and columns including `padj` (adjusted p-value), `log2FoldChange`, and optionally `process_name`.
#' @param plim Numeric threshold for adjusted p-value significance filtering. Only pathways with `padj` less than this
#' threshold will be plotted. Default is 0.001.
#' @param name Character string used as the plot title prefix and PDF filename.
#' @param outdir Character string specifying the directory path where the PDF plot will be saved.
#' @param include_longnames Logical indicating whether to append long descriptive names (from `process_name` column)
#' to pathway names in the plot. Default is FALSE.
#' @param w Numeric width of the output PDF in inches. Default is 12.
#' @param h Numeric height of the output PDF in inches. Default is 12.
#'
#' @return A ggplot2 object representing the generated barplot.
#'
#' @details
#' The function first converts the row names of the input data frame into a column named `Pathway`.
#' If `include_longnames` is TRUE and the `process_name` column exists, it cleans and appends descriptive
#' names to the pathway labels for better interpretability.
#'
#' It then filters the table to keep only pathways with an adjusted p-value below `plim` and orders them
#' by log2 fold change. The `Pathway` factor levels are set to maintain the sorted order for plotting.
#'
#' The barplot shows the log2 fold change on the x-axis, with bars colored to indicate whether the
#' pathway is more or less abundant (positive or negative fold change). The plot includes a horizontal
#' reference line at zero and is flipped to show pathway names on the y-axis for readability.
#'
#' The plot is saved as a PDF to the specified output directory using the provided name.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Assuming `df_results` contains differential abundance results with appropriate columns
#'   makeBarplotFunctional(df_results, plim = 0.01, name = "Differential_Pathways",
#'                         outdir = "./plots", include_longnames = TRUE)
#' }
#' }
#'
#' @seealso
#'  \code{\link[dplyr]{filter}}, \code{\link[dplyr]{arrange}}, \code{\link[dplyr]{mutate}}, \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{geom_bar}}
#'
#' @rdname makeBarplotFunctional
#' @export
#' @importFrom dplyr filter arrange mutate
#' @importFrom tibble rownames_to_column
#' @importFrom ggplot2 ggplot geom_hline geom_bar coord_flip scale_fill_manual scale_color_manual theme guides theme element_text ggtitle ggsave
makeBarplotFunctional <- function(tab, plim=0.001, name, outdir,
                                  include_longnames = FALSE,
                                  w=12, h=12){
  tab <- tab %>% rownames_to_column("Pathway")
  if(include_longnames & "process_name" %in% names(tab)){
    pnames <- gsub(" \\[[a-zA-Z0-9\\.\\]:\\- ]+", "", tab$process_name, perl=T) %>% gsub(" / ", "/", .)
    tab$Pathway <- paste(tab$Pathway, pnames, sep=": ")
  }
  tab2 <- tab %>% dplyr::filter(padj < plim) %>%
    dplyr::arrange(log2FoldChange) %>%
    dplyr::mutate(Pathway= factor(Pathway, levels = Pathway),
                  direction = ifelse(log2FoldChange<0, "less abundant", "more abundant"),
                  #direction = factor(direction, levels=legendnames),
                  `-10log(padj)` = -10*log10(padj)
    )


  g1 <- ggplot(tab2, aes(x=Pathway, y = log2FoldChange, fill=direction, col=direction ))+ #`-10log(padj)`
    geom_hline(yintercept = 0, linetype=3, col="gray80") +
    geom_bar(stat="identity") +
    coord_flip()+
    mytheme +
    scale_fill_manual(values = c("#A1C6EA", "#FD8B2F")) +
    scale_color_manual(values = c("#A1C6EA", "#FD8B2F")) +
    #scale_fill_gradient2(low="steelblue1", "#DDDDDD", high="tomato")+
    #scale_color_gradient2(low="steelblue1", "#DDDDDD", high="tomato")+
    theme(strip.text.y = element_text(size = 10,
                                      colour = "black", angle = 0, face = "bold"))+
    guides(fill=guide_legend(title="LFC"), colour=guide_legend(title="LFC"))+
    theme(panel.grid = element_blank())+
    ggtitle(paste0(name, " - proc. with padj<", as.character(plim)))+
    theme(axis.text.y = element_text(size = 12,
                                     colour = "black", angle = 0, face = "plain"))

  ggsave(filename = paste0(outdir, "/", name, ".pdf"), g1, width = w, height = h)
  return(g1)
}
