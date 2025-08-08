#' @title PCA Plotting Function with Multiple PCs and Grouping Variables
#' @description Performs PCA on selected variables of a data frame with optional transformations, then creates scatterplots of PC1 vs multiple other PCs colored by grouping variables. Saves plots to PDF.
#' @param dfraw A data.frame containing the raw data with samples as rows and variables as columns.
#' @param vars2pc Character vector of column names in `dfraw` to use for PCA.
#' @param labvar Name of the column in `dfraw` that contains sample labels (default: "PACIENTE").
#' @param plotvars Character vector of column names in `dfraw` used to color the PCA plots (default: "Condition").
#' @param transform Character specifying how to transform the data before PCA. Options: "scale" (standardize variables), "log" (log-transform +1), "sqrt" (square root), or other (no transformation). Default is "scale".
#' @param dims Integer vector specifying which principal components (PCs) to plot against PC1 (default: 2:5).
#' @param name Filename for the output PDF containing the plots (default: "PCAs.pdf").
#' @param outdir Directory where the PDF will be saved (default: current directory).
#' @param w Width of the saved plot in inches (default: 12).
#' @param h Height of the saved plot in inches (default: 8).
#' @return A list with two elements:
#' \itemize{
#'   \item \code{pca} - the result of the `prcomp` PCA object,
#'   \item \code{plots} - a named list of combined ggplot objects for each grouping variable.
#' }
#' @details
#' The function first subsets the data to the specified variables, applies the chosen transformation, and performs PCA with `prcomp`.
#' It then merges the PCA rotation results with metadata for plotting.
#' For each grouping variable, it creates a combined plot of PC1 against the PCs in `dims`.
#' The plots are saved into a multi-page PDF using a helper function `WriteManyPlots`.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   df <- data.frame(PACIENTE=letters[1:10], Var1=rnorm(10), Var2=rnorm(10), Condition=rep(c("A","B"),5))
#'   plotPCA2(df, vars2pc=c("Var1", "Var2"), labvar="PACIENTE", plotvars="Condition", transform="scale")
#' }
#' }
#' @seealso \code{\link[cowplot]{plot_grid}}
#' @rdname plotPCA2
#' @export
#' @importFrom cowplot plot_grid
plotPCA2 <- function(dfraw, vars2pc, labvar = "PACIENTE", plotvars=c("Condition"),
                     transform = "scale", dims=2:5, name = "PCAs.pdf", outdir="", w=12, h=8){

  dfraw <- dfraw %>% column_to_rownames(labvar)
  mat <- dfraw[, vars2pc] %>%
    as.matrix

  if(transform=="scale"){
    mat <- mat %>% scale %>% t
  }else if(transform == "log"){
    mat <- log(mat + 1) %>% t
  }else if(transform == "sqrt"){
    mat <- mat %>% sqrt %>% t
  }else{
    mat <- mat %>% t
  }

  dfraw$sample <- rownames(dfraw)
  pca <- prcomp(mat)
  df <- pca$rotation %>% as.data.frame %>%
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
      # if(! is.numeric(df2[, plotvar])){
      #   g1 <- g1 + scale_color_lancet()
      # }
      return(g1)
    })
    plotlist[[plotvar]] <- cowplot::plot_grid(plotlist=gs)
  }

  WriteManyPlots(plotlist,name = name, outdir = outdir, w = w, h=h)
  dev.off()
  return(list(pca=pca, plots=plotlist))
}
