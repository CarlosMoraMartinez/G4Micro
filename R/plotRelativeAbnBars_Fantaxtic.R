#' @title Barplot of Top Taxa Using fantaxtic
#' @description Plots a barplot of the top `n` most abundant taxa at a given taxonomic level using the `fantaxtic` package, optionally grouped by a sample variable and colored with a custom Wes Anderson palette.
#' @param phobj A `phyloseq` object containing microbiome data.
#' @param variable Character string indicating the sample metadata variable used to facet the barplots (e.g., "Condition"). Default: 'Condition'
#' @param topn Integer indicating the number of top abundant taxa to include in the plot. Remaining taxa are grouped as "Other". Default: 15
#' @param tax_level Taxonomic level to aggregate by (e.g., "Genus", "Family", etc.). Default: 'Genus'
#' @param outname Output file name for saving the plot as a PDF. Default: 'GenusBarplotFx.pdf'
#' @param height Height of the output PDF plot in inches. Default: 7
#' @param width Width of the output PDF plot in inches. Default: 12
#' @param wespalette Name of the Wes Anderson color palette (from `wesanderson` package) to use for coloring the taxa. Default: 'AsteroidCity1'
#' @return A `ggplot2` object representing the stacked barplot of relative abundances.
#' @details This function uses `fantaxtic` to prepare and plot relative abundance data of the top
#' taxa across sample groups. The top `n` taxa are
#' extracted using `get_top_taxa`, and "Other" groups all remaining taxa.
#' The final barplot is faceted by the `variable`, and colored using a selected palette
#' from the `wesanderson` package. The plot is saved to PDF and also returned as a `ggplot`
#' object for further manipulation.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  data("GlobalPatterns") # example phyloseq dataset
#'
#'   # Plot top 10 genera faceted by SampleType, save as PDF
#'   p <- plotRelativeAbnBars_Fantaxtic(
#'     phobj = GlobalPatterns,
#'     variable = "SampleType",
#'     topn = 10,
#'     tax_level = "Genus",
#'     outname = "GlobalPatterns_GenusBarplot.pdf",
#'     height = 28,
#'     width = 10,
#'     wespalette = "Darjeeling1"
#'   )
#'   print(p)
#'  }
#' }
#' @rdname plotRelativeAbnBars_Fantaxtic
#' @export
#' @importFrom fantaxtic fantaxtic_bar get_top_taxa name_taxa
#' @importFrom wesanderson wes_palette
plotRelativeAbnBars_Fantaxtic <- function(phobj, variable="Condition", topn = 15,
                                          tax_level = "Genus",
                                          outname="GenusBarplotFx.pdf", height=7, width=12,
                                          wespalette="AsteroidCity1"){
  topntx <- fantaxtic::get_top_taxa(physeq_obj = phobj, n = topn, relative = T,
                         discard_other = T, other_label = "Other")

  mycolors <- colorRampPalette(wesanderson::wes_palette(wespalette))(topn)
  topntx <- fantaxtic::name_taxa(topntx, label = "", species = F, other_label = "Other")
  topntx <- fantaxtic::fantaxtic_bar(topntx, color_by = tax_level, label_by = tax_level,
                          facet_by = variable, grid_by = NULL,
                          other_color = "Grey") +
    mytheme +
    theme(axis.text.x = element_text(size = 10,
                                     colour = "black", angle = 90,
                                     face = "plain", hjust=1, vjust=1)) +
    scale_fill_manual(values = mycolors)
  #theme(strip =element_rect(fill="white"))+
  ggsave(filename = outname, topntx, height = height, width = width)
  return(topntx)
}
