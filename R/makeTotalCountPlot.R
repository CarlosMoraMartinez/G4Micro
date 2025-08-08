#' @title Plot Total Read Counts by Group
#'
#' @description
#' Creates and saves a boxplot of total read counts (\code{"nreads"})
#' grouped by a specified metadata variable. The plot uses both
#' \code{\link[ggplot2]{geom_boxplot}} and \code{\link[ggplot2]{geom_jitter}}
#' to show distribution and individual sample points.
#'
#' @param metadata A \code{data.frame} or tibble containing sample metadata,
#'   including a numeric column named \code{"nreads"}.
#' @param variable Character string giving the name of a column in \code{metadata}
#'   to group samples by in the plot.
#' @param outdir Character string specifying the output directory where
#'   the PDF plot will be saved.
#' @param name Optional character string to prepend to the saved file name.
#'   Default is \code{""}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object representing the generated plot.
#'
#' @details
#' The resulting PDF will be saved in \code{outdir} with a file name in the
#' format:
#' \code{"<name>_<variable>.pdf"}
#' If \code{name} is empty, only the grouping variable name will be used.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Example metadata
#'   meta <- data.frame(
#'     group = rep(c("A", "B"), each = 5),
#'     nreads = c(1000, 1500, 1100, 1700, 1600, 2000, 1800, 1900, 2100, 2500)
#'   )
#'
#'   # Make plot and save to 'plots' folder
#'   makeTotalCountPlot(meta, variable = "group", outdir = "plots", name = "experiment1")
#' }
#' }
#'
#' @seealso
#' \code{\link[ggplot2]{ggplot}},
#' \code{\link[ggplot2]{geom_boxplot}},
#' \code{\link[ggplot2]{geom_jitter}},
#' \code{\link[ggplot2]{ggsave}}
#' @rdname makeTotalCountPlot
#' @export
makeTotalCountPlot <- function(metadata, variable, outdir, name = ""){
  g1 <- ggplot(metadata, aes_string(x=variable, y="nreads", col=variable, fill=variable)) +
    geom_boxplot(alpha=0.6, width=0.5) +
    #scale_color_lancet() +
    #scale_fill_lancet() +
    ylab("Total counts") +
    geom_jitter() + #size=3*samplesums$IMC/max(samplesums$IMC[!is.na(samplesums$IMC)])
    mytheme
  outname <- paste0(outdir,"/", name, '_', variable, ".pdf")
  ggsave(filename=outname, g1)
  return(g1)
}
