#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param phobj PARAM_DESCRIPTION
#' @param variable PARAM_DESCRIPTION, Default: 'Condition'
#' @param topn PARAM_DESCRIPTION, Default: 15
#' @param tax_level PARAM_DESCRIPTION, Default: 'Genus'
#' @param outname PARAM_DESCRIPTION, Default: 'GenusBarplotFx.pdf'
#' @param height PARAM_DESCRIPTION, Default: 7
#' @param width PARAM_DESCRIPTION, Default: 12
#' @param wespalette PARAM_DESCRIPTION, Default: 'AsteroidCity1'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname plotRelativeAbnBars_Fantaxtic
#' @export 
plotRelativeAbnBars_Fantaxtic <- function(phobj, variable="Condition", topn = 15,
                                          tax_level = "Genus",
                                          outname="GenusBarplotFx.pdf", height=7, width=12,
                                          wespalette="AsteroidCity1"){
  library(fantaxtic)
  topntx <- get_top_taxa(physeq_obj = phobj, n = topn, relative = T,
                         discard_other = T, other_label = "Other")

  mycolors <- colorRampPalette(wes_palette(wespalette))(topn)
  topntx <- name_taxa(topntx, label = "", species = F, other_label = "Other")
  topntx <- fantaxtic_bar(topntx, color_by = tax_level, label_by = tax_level,
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
