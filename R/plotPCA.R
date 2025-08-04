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
    geom_text_repel(aes(label = sample)) +
    xlab(getPropVar(pca, "PC1")) +
    ylab(getPropVar(pca, "PC2")) +
    mystyle +
    guides(col=guide_legend(title=legname))
  g2 <- ggplot(df, aes(x=PC1, y=PC3, col=group)) +
    geom_point() +
    geom_text_repel(aes(label = sample)) +
    xlab(getPropVar(pca, "PC1")) +
    ylab(getPropVar(pca, "PC3")) +
    mystyle +
    guides(col=guide_legend(title=legname))
  g3 <- ggplot(df, aes(x=PC1, y=PC4, col=group)) +
    geom_point() +
    geom_text_repel(aes(label = sample)) +
    xlab(getPropVar(pca, "PC1")) +
    ylab(getPropVar(pca, "PC4")) +
    mystyle +
    guides(col=guide_legend(title=legname))
  g4 <- ggplot(df, aes(x=PC1, y=PC5, col=group)) +
    geom_point() +
    geom_text_repel(aes(label = sample)) +
    xlab(getPropVar(pca, "PC1")) +
    ylab(getPropVar(pca, "PC5")) +
    mystyle +
    guides(col=guide_legend(title=legname))
  cw <- cowplot::plot_grid(g1, g2, g3, g4)

  return(list(pca=pca, plots=cw))
}
