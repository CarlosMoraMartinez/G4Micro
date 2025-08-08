#' @title Generate and Save Rarefaction Curves from a Phyloseq Object
#'
#' @description
#' Calculates rarefaction curves for each sample in a
#' \code{\link[phyloseq]{phyloseq}} object using \code{\link[vegan]{rarecurve}},
#' saves the plots and summary tables to disk, and returns the data and plot.
#'
#' @param phyloseq_rawdata A \code{\link[phyloseq]{phyloseq}} object containing
#'   an OTU/ASV table. The table will be extracted and converted to a matrix
#'   for rarefaction analysis.
#' @param opt A list containing at least an element named \code{"out"} specifying
#'   an existing output directory where plots and tables will be saved.
#' @param add_hlines Logical indicating whether to add horizontal lines at the
#'   interpolated number of taxa for the minimum read depth across samples.
#'   Default is \code{TRUE}.
#' @param name Character string used as the base name for output files
#'   (tables, plots, and RData). Default is \code{"rarecurve"}.
#'
#' @return A named list with the following elements:
#' \itemize{
#'   \item \code{plot} — A \code{\link[ggplot2]{ggplot}} object of the rarefaction curves.
#'   \item \code{df} — A \code{data.frame} with the rarefaction curve data
#'     (taxa counts at each subsampling depth for each sample).
#'   \item \code{df_max} — A \code{data.frame} with the maximum read depth per sample.
#'   \item \code{min_reads_all} — The overall minimum read depth across samples.
#'   \item \code{intersection2} — A \code{data.frame} with the interpolated taxa
#'     counts at the overall minimum read depth.
#' }
#'
#' @details
#' Rarefaction curves visualize the accumulation of observed taxa as a function
#' of sequencing depth. This function computes them using
#' \code{\link[vegan]{rarecurve}} and also provides ggplot-based versions for
#' cleaner presentation. Summary tables include:
#' \enumerate{
#'   \item All rarefaction curve points (\code{df}).
#'   \item Final point for each sample (\code{df_max}).
#'   \item Interpolated taxa counts at the overall minimum read depth
#'         (\code{intersection2}).
#' }
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   library(phyloseq)
#'   # Example: run rarefaction on a phyloseq object 'ps'
#'   res <- makeRarefactionCurve(ps, opt = list(out = "results"))
#'   res$plot
#' }
#' }
#'
#' @seealso
#' \code{\link[vegan]{rarecurve}},
#' \code{\link[ggplot2]{ggplot}},
#' \code{\link[dplyr]{mutate}},
#' \code{\link[dplyr]{summarise}},
#' \code{\link[ggrepel]{geom_text_repel}},
#' \code{\link[readr]{write_tsv}}
#'
#' @rdname makeRarefactionCurve
#' @export
#'
#' @importFrom vegan rarecurve
#' @importFrom phyloseq otu_table
#' @importFrom dplyr mutate summarise group_by filter bind_rows slice_min
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_vline geom_hline xlab ylab ggtitle scale_x_continuous theme ggsave
#' @importFrom ggrepel geom_text_repel
#' @importFrom readr write_tsv
#' @importFrom stats setNames
makeRarefactionCurve <- function(phyloseq_rawdata, opt, add_hlines=TRUE, name="rarecurve"){
  otu2rare <- otu_table(phyloseq_rawdata)
  class(otu2rare) <- "matrix"
  png(paste0(opt$out, "/00_plot_rarecurve.png"), width = 15, height = 15, units = "cm", res = 300)

  rarplot <- rarecurve(t(otu2rare),
                       step = 10000,
                       #sample = 20000,
                       col = "blue",
                       cex = 0.75,
                       main = "Rarefaction curve")
  dev.off()

  #Table with all data
  rardf <- mapply(rarplot,  colnames(otu2rare), SIMPLIFY=FALSE,
                  FUN=function(x, s) data.frame(sample=s, nreads = attr(x, "Subsample"), txcount=x)) %>%
    bind_rows() %>%
    group_by(sample) %>%
    dplyr::mutate(last_point = ifelse(nreads == max(nreads), 2, 0))

  #Table with max reads in each sample, to plot points
  rarmax <- rardf %>% filter(last_point == 2)

  #Table with the intersection between min reads of all samples and each of the curves, to plot hlines
  MIN_READS_ALLSAMPLES <- min(rarmax$nreads)

  intersection <- rardf %>%
    dplyr::mutate(diff2min = abs( nreads - MIN_READS_ALLSAMPLES)) %>%
    slice_min(n=2, order_by=diff2min) %>%
    dplyr::mutate(prop_dist = 1- diff2min/sum(diff2min))
  intersection2 <- intersection %>% dplyr::summarise(nreads=MIN_READS_ALLSAMPLES,
                                                     txcount = sum(txcount*prop_dist))
  write_tsv(rardf, paste0(opt$out, "/00_plot_", name, "2.tsv"))
  write_tsv(rarmax, paste0(opt$out, "/00_plot_", name, "_max.tsv"))
  write_tsv(intersection2, paste0(opt$out, "/00_plot_", name, "_intersect.tsv"))

  ggrare <- ggplot(rardf, aes(col=sample, x=nreads, y=txcount)) +
    geom_line() +
    geom_point(data=rarmax)+
    geom_text_repel(data=rarmax, aes(label=sample))+
    geom_vline(xintercept=MIN_READS_ALLSAMPLES, col="gray40", linetype=2) +
    xlab("Sample Size") +
    ylab("Species") +
    ggtitle("Rarefaction curve") +
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
    mytheme +
    theme(legend.position="none")

  if(add_hlines) ggrare <- ggrare + geom_hline(data=intersection2,
                                               aes(yintercept=txcount, col=sample), linetype=2, size=0.2)

  ggsave(paste0(opt$out, "/00_plot_", name, "2.pdf"), ggrare)

  alldata <- list(
    plot = ggrare,
    df = rardf,
    df_max = rarmax,
    min_reads_all = MIN_READS_ALLSAMPLES,
    intersection2 = intersection2
  )
  save(alldata, file=paste0(opt$out, "/00_plot_", name, "_all.RData"))

  return(alldata)

}
