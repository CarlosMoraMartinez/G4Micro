#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dds PARAM_DESCRIPTION
#' @param raw_counts PARAM_DESCRIPTION
#' @param norm_counts PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{mutate}}
#' @rdname plotCountsDeseq
#' @export 
#' @importFrom dplyr mutate
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
