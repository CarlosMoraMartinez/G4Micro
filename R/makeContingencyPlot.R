#' @title Create and Save a Mosaic Plot for Two Categorical Variables
#' @description
#' Generates a mosaic plot showing the relationship between two categorical variables
#' and saves the resulting figure as a PDF.
#'
#' @param df A data frame containing the categorical variables to be plotted.
#' @param var1 Character string, name of the first categorical variable (displayed on the y-axis).
#' @param var2 Character string, name of the second categorical variable (displayed on the x-axis).
#' @param outdir Character string, directory path where the plot will be saved.
#' @param name Character string, base name for the output file (without extension).
#' @param w Numeric, plot width in inches. Default: 8.
#' @param h Numeric, plot height in inches. Default: 6.
#'
#' @return A \code{ggplot} object containing the mosaic plot.
#'
#' @details
#' The function uses \code{\link[ggmosaic]{geom_mosaic}} to generate a mosaic plot and
#' \code{\link[ggmosaic]{geom_mosaic_text}} to display cell proportions as labels.
#' The plot is styled with \code{theme_classic()} and a custom theme (\code{mytheme}).
#' The output PDF file is named as \code{paste0(outdir, "/", name, "_mosaic.pdf")}.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   df <- data.frame(
#'     Group = sample(c("A", "B"), 100, replace = TRUE),
#'     Outcome = sample(c("Yes", "No"), 100, replace = TRUE)
#'   )
#'   makeContingencyPlot(df, var1 = "Group", var2 = "Outcome",
#'                       outdir = tempdir(), name = "example_plot")
#' }
#' }
#'
#' @seealso
#'  \code{\link[ggmosaic]{geom_mosaic}},
#'  \code{\link[ggmosaic]{geom_mosaic_text}},
#'  \code{\link[ggplot2]{ggsave}}
#'
#' @rdname makeContingencyPlot
#' @export
#' @importFrom ggmosaic geom_mosaic geom_mosaic_text product
#' @importFrom ggplot2 ggplot aes theme_classic element_text xlab ylab ggsave
#' @importFrom rlang sym
makeContingencyPlot <- function(df, var1, var2, outdir, name, w=8, h=6){
  gmos <- ggplot(data = df) +
    geom_mosaic(aes(x = product(!!sym(var1), !!sym(var2)), fill=!!sym(var1)), na.rm=T) +
    geom_mosaic_text(aes(x = product(!!sym(var1), !!sym(var2)),
                         fill=!!sym(var1), label=after_stat(.wt)),
                     na.rm=T, as.label=T, size=6) +
    theme_classic() +
    mytheme +
    theme(axis.text.x=element_text(size=14),
          axis.title=element_text(size=14,face="bold"),
          axis.text.y = element_text(size = 14,
                                     colour = "black", angle = 90, face = "bold"),
          legend.position = 'none')+
    xlab(gsub("_", " ", var2))+
    ylab(gsub("_", " ", var1))
  fname <- paste0(outdir, "/", name, "_mosaic.pdf")
  ggsave(filename = fname, gmos, width = w, height = h)
  return(gmos)
}
