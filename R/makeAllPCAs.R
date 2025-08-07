#' @title Generate and Save Multiple PCA Plots
#' @description This function computes and saves PCA plots for a set of specified metadata variables, using count data and gene selection. It also exports PCA components and rotations.
#' @param phobj A phyloseq object containing sample metadata and the OTU table.
#' @param counts_df A data frame or matrix with normalized counts (e.g., variance-stabilized or log-transformed counts).
#' @param genes A vector of gene/taxa names, corresponding to rownames of \code{counts_df} to include in the PCA.
#' @param vars2pca A character vector with the names of metadata variables to use for PCA visualization (e.g., treatment groups, batches).
#' @param opt A list containing an `out` element, which specifies the output directory path.
#' @param name A character string to use as the output file name prefix, Default: 'PCAs'
#' @return A named list of PCA plot objects (one per variable in `vars2pca`). Also writes `.RData`, `.tsv`, and `.pdf` files to disk.
#' @details
#' This function uses the sample metadata and the normalized count matrix to compute PCA and generate plots colored by metadata variables.
#' It also appends a log10 transformed total read count (`reads_log10_current`) to the metadata and includes it in the PCA plot list, to visually check for depth biases.
#'
#' PCA computation and plotting are delegated to a helper function `plotPCA`, which must return a list with components `$pca` (containing the PCA object) and `$plots` (a ggplot2 object).
#' If plotting fails for a specific variable, a fallback null plot is generated using `getNullPlot`.
#'
#' Outputs:
#' - A PDF with all PCA plots.
#' - A `.RData` file with all PCA objects.
#' - TSV tables with PCA coordinates and rotations for the first variable.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   pcas <- makeAllPCAs(phobj, counts_df, genes, vars2pca = c("group", "sex"), opt = list(out = "results/"))
#' }
#' }
#' @seealso
#'  \code{\link[dplyr]{mutate}},
#'  \code{\link[phyloseq]{sample_data}},
#'  \code{\link[phyloseq]{otu_table}},
#'  \code{\link[readr]{write_tsv}},
#'  \code{\link[stats]{prcomp}},
#'  \code{\link{plotPCA}},
#'  \code{\link{getNullPlot}}
#' @rdname makeAllPCAs
#' @export
#' @importFrom dplyr mutate
#' @importFrom phyloseq sample_data otu_table
#' @importFrom readr write_tsv
makeAllPCAs <- function(phobj, counts_df, genes, vars2pca, opt, name = "PCAs"){
  design <- sample_data(phobj)
  design<- data.frame(design)
  nreads <- otu_table(phobj) %>% colSums()

  design <- design %>%
    dplyr::mutate(reads_log10_current = log10(nreads[as.character(sampleID)]))

  vars2pca <- c(vars2pca,"reads_log10_current")

  pca_plots <- lapply(vars2pca,  FUN=function(vv, counts, design, genes){
    plotPCA(counts,design, genes, vv)
  },counts_df, design, genes)
  names(pca_plots) <- vars2pca

  save(pca_plots, file = paste0(opt$out, name, "PCA.RData"))
  write_tsv(as.data.frame(pca_plots[[1]]$pca$x), file = paste0(opt$out, name, "PCAnewvars.tsv"))
  write_tsv(as.data.frame(pca_plots[[1]]$pca$rotation), file = paste0(opt$out, name, "PCArotation.tsv"))

  pdf(paste0(opt$out, name, ".pdf"), width=16, height = 12)
  for(pl in vars2pca){
    tryCatch({print(pca_plots[[pl]][["plots"]])}, error = function(x){print(getNullPlot(opt, pl))})
  }
  tmp <- dev.off()

  return(pca_plots)
}
