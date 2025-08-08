#' @title Plot Raw and Normalized Counts from DESeq2 Data
#' @description Generates plots summarizing raw and normalized counts for DESeq2
#' experiment samples, including total counts per condition and the distribution
#' of normalized counts across samples.
#' @param dds A DESeqDataSet object containing the DESeq2 experiment data, including sample metadata.
#' @param raw_counts A matrix or data frame of raw counts with genes as rows and samples as columns.
#' @param norm_counts A matrix or data frame of normalized counts matching \code{raw_counts} dimensions.
#' @return NULL (plots are drawn to the current graphics device).
#' @details
#' This function creates three plots arranged in a single row:
#' \itemize{
#'   \item Total raw counts per sample by condition with sample labels.
#'   \item Total normalized counts per sample by condition with sample labels.
#'   \item Density distribution of normalized counts per sample on a log2 scale,
#'         truncated at the 95th percentile for visualization clarity.
#' }
#' The function uses the \code{Condition} column from \code{colData(dds)} to group samples.
#' It requires \code{ggplot2}, \code{ggrepel}, \code{dplyr}, and \code{gridExtra} packages.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   dds <- DESeq2::makeExampleDESeqDataSet()
#'   raw_counts <- counts(dds)
#'   norm_counts <- DESeq2::counts(dds, normalized=TRUE)
#'   plotCountsDeseq(dds, raw_counts, norm_counts)
#' }
#' }
#' @seealso
#' \code{\link[DESeq2]{plotCounts}},
#' \code{\link[ggplot2]{ggplot}},
#' \code{\link[dplyr]{mutate}},
#' \code{\link[ggrepel]{geom_text_repel}},
#' \code{\link[gridExtra]{grid.arrange}}
#' @rdname plotCountsDeseq
#' @export
#' @import DESeq2
#' @importFrom dplyr mutate filter
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot aes geom_point geom_text_repel geom_density scale_x_continuous guides theme_bw ggtitle
#' @importFrom gridExtra grid.arrange
plotCountsDeseq <- function(dds, raw_counts, norm_counts){
  d <- plotCounts(dds, gene=which.min(res$padj), intgroup="Condition", normalized = T,
                  returnData=TRUE)
  d <- raw_counts %>% colSums() %>% as.data.frame()
  d$sample <- rownames(d)
  names(d)[1] <- "count"
  d$Condition <- colData(dds)$Condition[match( d$sample, colData(dds)$sampleID)]

  g1 <- ggplot(d, aes(x=Condition, y=count, col=Condition, fill=Condition)) +
    geom_point(position=position_jitter(w=0.1,h=0), size=2) +
    #  scale_y_log10(breaks=c(25,100,400)) +
    geom_text_repel(aes(label=sample)) +
    ggtitle("Total raw counts by sample") +
    theme_bw()

  d2 <- norm_counts %>% colSums() %>% as.data.frame()
  d2$sample <- rownames(d)
  names(d2)[1] <- "count"
  d2$Condition <- colData(dds)$Condition[match( d2$sample, colData(dds)$sampleID)]

  g2 <- ggplot(d2, aes(x=Condition, y=count, col=Condition, fill=Condition)) +
    geom_point(position=position_jitter(w=0.1,h=0), size=2) +
    #  scale_y_log10(breaks=c(25,100,400)) +
    geom_text_repel(aes(label=sample)) +
    ggtitle("Total normalized counts by sample") +
    theme_bw()

  nc_exp <- raw_counts %>% as.data.frame() %>%
    dplyr::mutate(gene=rownames(raw_counts)) %>%
    gather(key=sample, value='norm_counts', colnames(norm_counts)) %>%
    filter(norm_counts > 0)

  g3 <- ggplot(nc_exp, aes(x=norm_counts, col=sample))+
    geom_density() +
    xlim(x=c(0,  quantile(nc_exp$norm_counts, c(0.95))[1] ))+
    ggtitle("Distribution of normalized counts by sample") +
    theme_bw() +
    scale_x_continuous(trans='log2') +
    guides(color=FALSE)

  grid.arrange(g1, g2, g3, nrow = 1)

}
