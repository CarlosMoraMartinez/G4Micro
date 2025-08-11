
#' @title Create and Save a Sankey Plot Using Plotly
#' @description Generates an interactive Sankey diagram from given nodes and links data frames and saves it as an HTML file.
#' @param nodelist A data frame containing node information. Must include columns: \code{value} (node labels), \code{color} (node colors), \code{nodesize} (size values), and \code{xpos} (horizontal position of nodes).
#' @param linklist A data frame containing link information. Must include columns: \code{source} (index of source node), \code{target} (index of target node), \code{Size} (link values), \code{LFC} (link labels), and \code{color} (link colors).
#' @param name A character string specifying the plot title.
#' @param fname A character string specifying the filename (without extension) for the saved HTML widget.
#' @param outdir A character string specifying the directory path where the HTML file will be saved.
#' @param fontsize Numeric value specifying the font size for the plot title and labels. Default is 18.
#' @return Returns a \code{plotly} object representing the Sankey plot.
#' @details This function uses the \code{plotly} package to create an interactive Sankey diagram, with nodes and links colored according to the supplied data frames. It also truncates node labels longer than 60 characters by removing parentheses content to improve plot readability.
#' The resulting plot is saved as an HTML file in the specified output directory.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Example data frames nodelist and linklist must be created first
#'   plot <- makeSankeyPlot(nodelist, linklist, name = "My Sankey", fname = "mysankey", outdir = tempdir())
#'   print(plot)
#' }
#' }
#' @seealso
#' \code{\link[dplyr]{mutate}}, \code{\link[htmlwidgets]{saveWidget}}, \code{\link[plotly]{plot_ly}}
#' @rdname makeSankeyPlot
#' @export
#' @importFrom dplyr mutate
#' @importFrom htmlwidgets saveWidget
#'@importFrom plotly plot_ly layout
makeSankeyPlot <- function(nodelist, linklist, name, fname, outdir, fontsize=18){

  nodelist <- nodelist %>%
    dplyr::mutate(value  = ifelse(nchar(value) > 60, gsub("\\([^)]+\\)", "", value, perl=TRUE), value)
    )
  # From https://plotly.com/r/sankey-diagram/
  print(head(nodelist))
  fig <- plot_ly(
    type = "sankey",
    arrangement="fixed",
    domain = list(
      x =  c(0,1),
      y =  c(0,1)
    ),
    orientation = "h",
    valueformat = ".0f",
    valuesuffix = "TWh",
    node = list(
      label = nodelist$value,
      color = nodelist$color,
      value = nodelist$nodesize,
      x = nodelist$xpos,
      pad = 15,
      thickness = 15,
      line = list(
        color = "black",
        width = 0.5
      )
    ),
    link = list(
      source = linklist$source,
      target = linklist$target,
      value =  linklist$Size,
      label =  as.character(round(linklist$LFC, 2)),
      color =  linklist$color
    )
  )

  fig <- fig %>% layout(
    title = name,
    font = list(
      size = fontsize
    ),
    xaxis = list(showgrid = F, zeroline = F),
    yaxis = list(showgrid = F, zeroline = F)
  )

  #fig
  #tmp <- paste0(outdir, "/", fname, ".png")
  #save_image(fig, tmp)
  #use_virtualenv("r-reticulate")
  htmlwidgets::saveWidget(fig, file=paste0(outdir, fname, ".html"))
  return(fig)
}
