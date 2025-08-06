#' @title makeRarefactionCurve
#' @description Calculates a rarefaction curve for a phyloseq object
#' @param phyloseq_rawdata Input phyloseq object
#' @param opt Option list, including a field named 'out' with an existing output directory.
#' @param add_hlines If TRUE, plot horizontal lines with the depth of each sample, Default: TRUE
#' @param name name of the output files, including tables, R objects and plots saved to opt$out, Default: 'rarecurve'
#' @return A list with the following elements:
#' `plot`: ggplot object with the rarefaction curve
#' `df`: DataFrame with the rarefaction curve data (number of taxa detected for each number of reads for each sample)
#' `df_max`: DataFrame with the number of reads of each sample
#' `min_reads_all`: Overall minimum number of reads
#' `intersection2`: DataFrame with the number of taxa at the overall minimum number of reads, extrapolated from the curve
#' @details Calculates a rarefaction curve for a phyloseq object using the vegan \code{rarecurve()} function
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[vegan]{rarecurve}},
#'  \code{\link[dplyr]{mutate}},
#'  \code{\link[dplyr]{summarise}}
#' @rdname makeRarefactionCurve
#' @export
#' @importFrom vegan rarecurve
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
