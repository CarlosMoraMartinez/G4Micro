#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dfraw PARAM_DESCRIPTION
#' @param vars2pc PARAM_DESCRIPTION
#' @param labvar PARAM_DESCRIPTION, Default: 'PACIENTE'
#' @param plotvars PARAM_DESCRIPTION, Default: c("Condition")
#' @param transform PARAM_DESCRIPTION, Default: 'scale'
#' @param dims PARAM_DESCRIPTION, Default: 2:5
#' @param name PARAM_DESCRIPTION, Default: 'PCAs.pdf'
#' @param outdir PARAM_DESCRIPTION, Default: ''
#' @param w PARAM_DESCRIPTION, Default: 12
#' @param h PARAM_DESCRIPTION, Default: 8
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[cowplot]{plot_grid}}
#' @rdname plotPCA_pred2
#' @export 
#' @importFrom cowplot plot_grid
plotPCA_pred2 <- function(dfraw, vars2pc, labvar = "PACIENTE", plotvars=c("Condition"),
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
  }#else{
  #   mat <- mat
  # }

  dfraw$sample <- rownames(dfraw)
  pca <- prcomp(mat, center = T, scale. = T)
  df <- pca$x %>% as.data.frame %>%
    rownames_to_column("sample")
  df2 <- merge(df, dfraw, by="sample")

  plotlist <- list()
  for(plotvar in plotvars){

    legname <- gsub("_", " ", plotvar)
    legname <- gsub("\\.", "-", legname)
    gs <- lapply(paste0("PC", dims), FUN=function(PC){
      g1 <- ggplot(df2, aes_string(x="PC1", y=PC, col=plotvar)) +
        geom_point() +
        geom_text_repel(aes(label = sample)) +
        xlab(getPropVar(pca, "PC1")) +
        ylab(getPropVar(pca, PC)) +
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
  dev.off()
  return(list(pca=pca, plots=plotlist))
}
