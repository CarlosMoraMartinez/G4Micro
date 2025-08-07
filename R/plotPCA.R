#' @title PCA Plot Generator from Normalized Counts
#' @description Performs PCA on a selected subset of genes and generates four scatter plots showing PC1 vs PC2–PC5, colored by a specified metadata variable.
#' @param counts A data frame with normalized counts where rows are taxa/genes and columns are samples. It must include a "gene" column.
#' @param design A data frame containing sample metadata. Must include a "sampleID" column and the column named in `plotvar`.
#' @param genes A character vector of gene/taxa names to retain for PCA.
#' @param plotvar The metadata variable used to color the PCA plots, Default: 'Condition'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{pca}: the result of \code{\link[stats]{prcomp}}.
#'   \item \code{plots}: a combined \code{\link[ggplot2]{ggplot}} object (via \code{\link[cowplot]{plot_grid}}) with four PCA plots (PC1 vs PC2, PC3, PC4, and PC5).
#' }
#' @details
#' The function filters the normalized count matrix to include only the provided genes, scales and transposes the matrix for PCA, and joins metadata for plotting.
#' If the resulting matrix has fewer than 3 rows (samples), a placeholder plot is returned via \code{\link{getNullPlot}}.
#' It uses `ggplot2` for plotting, `ggrepel` for labeling, and a custom theme (`mystyle`) for aesthetics.
#' The proportion of variance explained by each principal component is shown in axis labels (via \code{\link{getPropVar}}).
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   result <- plotPCA(counts, design, genes = c("geneA", "geneB", "geneC"), plotvar = "group")
#'   result$plots
#' }
#' }
#' @seealso
#'  \code{\link[dplyr]{filter}},
#'  \code{\link[dplyr]{select}},
#'  \code{\link[tibble]{column_to_rownames}},
#'  \code{\link[stats]{prcomp}},
#'  \code{\link[ggplot2]{ggplot}},
#'  \code{\link[ggrepel]{geom_text_repel}},
#'  \code{\link[cowplot]{plot_grid}},
#'  \code{\link{getPropVar}},
#'  \code{\link{getNullPlot}}
#' @rdname plotPCA
#' @export
#' @importFrom dplyr filter select
#' @importFrom cowplot plot_grid
#' @importFrom ggrepel geom_text_repel
plotPCA <- function(counts, design, genes, plotvar="Condition"){

  mat <- counts %>%
    dplyr::filter(gene %in% genes) %>%
    column_to_rownames("gene") %>%
    dplyr::select(design$sampleID) %>%
    as.matrix %>%
    t %>% scale

  if(nrow(mat) < 3){
    cw <- getNullPlot()
    return(cw)
  }
  pca <- prcomp(mat)
  df <- pca$x %>% as.data.frame %>%
    rownames_to_column("sample")
  df$group <- design[match(df$sample, design$sampleID), plotvar] %>% unlist

  legname <- gsub("_", " ", plotvar)
  legname <- gsub("\\.", "-", legname)

  g1 <- ggplot(df, aes(x=PC1, y=PC2, col=group)) +
    geom_point() +
    ggrepel::geom_text_repel(aes(label = sample)) +
    xlab(getPropVar(pca, "PC1")) +
    ylab(getPropVar(pca, "PC2")) +
    mystyle +
    guides(col=guide_legend(title=legname))
  g2 <- ggplot(df, aes(x=PC1, y=PC3, col=group)) +
    geom_point() +
    ggrepel::geom_text_repel(aes(label = sample)) +
    xlab(getPropVar(pca, "PC1")) +
    ylab(getPropVar(pca, "PC3")) +
    mystyle +
    guides(col=guide_legend(title=legname))
  g3 <- ggplot(df, aes(x=PC1, y=PC4, col=group)) +
    geom_point() +
    ggrepel::geom_text_repel(aes(label = sample)) +
    xlab(getPropVar(pca, "PC1")) +
    ylab(getPropVar(pca, "PC4")) +
    mystyle +
    guides(col=guide_legend(title=legname))
  g4 <- ggplot(df, aes(x=PC1, y=PC5, col=group)) +
    geom_point() +
    ggrepel::geom_text_repel(aes(label = sample)) +
    xlab(getPropVar(pca, "PC1")) +
    ylab(getPropVar(pca, "PC5")) +
    mystyle +
    guides(col=guide_legend(title=legname))
  cw <- cowplot::plot_grid(g1, g2, g3, g4)

  return(list(pca=pca, plots=cw))
}
