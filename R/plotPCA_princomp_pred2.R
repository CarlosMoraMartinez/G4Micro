
plotPCA_princomp_pred2 <- function(dfraw, vars2pc, labvar = "PACIENTE", plotvars=c("Condition"),
                             transform = "scale", dims=2:5, name = "PCAs.pdf", outdir="", w=12, h=8){

  dfraw <- dfraw %>% column_to_rownames(labvar)
  mat <- dfraw[, vars2pc] %>%
    as.matrix

  if(transform=="scale"){
    mat <- mat %>% scale
  }else if(transform == "log"){
    mat <- log(mat + 1)
  }else if(transform == "sqrt"){
    mat <- mat %>% sqrt
  }else{
    mat <- mat
  }

  dfraw$sample <- rownames(dfraw)
  pca <- princomp(mat, cor=TRUE)
  df <- pca$scores %>% as.data.frame
  names(df) <- gsub("Comp.", "PC", names(df))
  df2 <- df %>%  rownames_to_column("sample") %>%
    merge( dfraw, by="sample")

  pcts <- getPropVar_princomp(pca)
  plotlist <- list()
  for(plotvar in plotvars){

    legname <- gsub("_", " ", plotvar)
    legname <- gsub("\\.", "-", legname)
    gs <- lapply(paste0("PC", dims), FUN=function(PC){
      g1 <- ggplot(df2, aes_string(x="PC1", y=PC, col=plotvar)) +
        geom_point() +
        geom_text_repel(aes(label = sample)) +
        xlab(pcts["PC1"]) +
        ylab(pcts[PC]) +
        mystyle +
        guides(col=guide_legend(title=legname))
      if(! is.numeric(df2[, plotvar])){
        g1 <- g1 + scale_color_lancet()
      }
      return(g1)
    })
    plotlist[[plotvar]] <- cowplot::plot_grid(plotlist=gs)
  }

  WriteManyPlots(plotlist,name = name, outdir = outdir, w = w, h=h)
  return(list(pca=pca, plots=plotlist))
}
