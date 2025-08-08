#' @title Generate Correlation Heatmaps by Group
#' @description
#' Creates correlation heatmaps between two matrices (e.g., microbial abundances and metadata values) for all samples combined and for each group within a given metadata variable.
#' The function first generates a global correlation heatmap and then produces separate heatmaps for each level of a grouping variable, arranging them together for comparison.
#'
#' @param mat1 A numeric matrix or data frame (samples in rows, variables/features in columns). Typically microbial abundances or similar quantitative data.
#' @param mat2 A numeric matrix or data frame (samples in rows, variables/features in columns). Typically metadata variables, functional abundances, or other measurements to correlate with `mat1`.
#' @param metadata A data frame containing metadata for the samples. Must include both `var2plot` and `varwithnames` columns.
#' @param var2plot Character string specifying the name of the metadata column used to group samples (default: `"Psoriasis"`).
#' @param varwithnames Character string specifying the column in `metadata` that contains sample identifiers (default: `"sampleID"`).
#' @param cormethod Character string indicating the correlation method to use. Passed to \code{\link[stats]{cor}} (default: `"pearson"`).
#' @param pval Numeric value specifying the p-value threshold for including correlations in the heatmap (default: `0.05`).
#' @param select_asvs Character vector of ASVs or feature names to plot. Overrides the `pval` threshold if provided (default: empty vector).
#' @param outdir Character string specifying the directory where the output PDF file will be saved (default: `""` = current directory).
#' @param name Character string specifying the base name for output files (default: `"corrHeatmap"`).
#' @param clust Logical, whether to cluster rows and columns in the heatmap (default: `TRUE`).
#' @param order_rows Optional character vector specifying a custom row order for the heatmap (default: empty vector).
#' @param order_cols Optional character vector specifying a custom column order for the heatmap (default: empty vector).
#' @param annot_asvs Optional data frame with annotations for ASVs/features (default: `NULL`).
#'
#' @return A \code{\link[pheatmap]{pheatmap}} object representing the final combined correlation heatmap.
#'
#' @details
#' This function is designed to compare correlation patterns between two datasets across different groups defined in a metadata column.
#' The workflow:
#' 1. Generate a correlation heatmap for all samples.
#' 2. Generate separate heatmaps for each group (level of `var2plot`).
#' 3. Combine results into one figure, showing the global pattern followed by group-specific patterns.
#' Only significant correlations (based on `pval`) or features in `select_asvs` are plotted.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Example data
#'   set.seed(123)
#'   mat1 <- matrix(rnorm(50), nrow = 10)
#'   rownames(mat1) <- paste0("Sample", 1:10)
#'   colnames(mat1) <- paste0("ASV", 1:5)
#'
#'   mat2 <- matrix(rnorm(30), nrow = 10)
#'   rownames(mat2) <- paste0("Sample", 1:10)
#'   colnames(mat2) <- paste0("Var", 1:3)
#'
#'   metadata <- data.frame(
#'     sampleID = paste0("Sample", 1:10),
#'     Group = rep(c("A", "B"), each = 5)
#'   )
#'
#'   makeCorrelationHeatmapByGroup(
#'     mat1 = mat1,
#'     mat2 = mat2,
#'     metadata = metadata,
#'     var2plot = "Group",
#'     varwithnames = "sampleID",
#'     outdir = tempdir(),
#'     name = "example_corrHeatmap"
#'   )
#' }
#' }
#'
#' @seealso
#' \code{\link[dplyr]{select}}, \code{\link[pheatmap]{pheatmap}}
#' @rdname makeCorrelationHeatmapByGroup
#' @export
#' @importFrom dplyr select
#' @importFrom pheatmap pheatmap
makeCorrelationHeatmapByGroup <- function(mat1, mat2, metadata,var2plot="Psoriasis",
                                          varwithnames = "sampleID",
                                          cormethod="pearson",
                                          pval=0.05, select_asvs=c(),
                                          outdir="", name="corrHeatmap",
                                          clust=T, order_rows = c(),
                                          order_cols = c(),
                                          annot_asvs = NULL){

  hmall <- makeCorrelationHeatmap(mat1, mat2, cormethod=cormethod,
                                  pval=0.05, select_asvs=select_asvs, outdir=outdir,
                                  name=paste0(name, "_allsamples"), clust=T)
  default_annot_colors <- c("gray40", "dodgerblue3", "firebrick3", "green4", "pink4", "skyblue")
  levs <- levels(metadata[, var2plot])
  hmlevs <- lapply(levs, FUN=function(lev){
    samples2select <- metadata[ metadata[, var2plot]==lev, varwithnames]
    mat1_tmp <- mat1[rownames(mat1) %in% samples2select, ]
    mat2_tmp <- mat2[rownames(mat2) %in% samples2select, ]
    hmlist <- makeCorrelationHeatmap(mat1_tmp, mat2_tmp, cormethod=cormethod,
                                     pval=0.05, select_asvs=hmall$asvs, outdir=outdir,
                                     name=paste0(name,"_", var2plot,"_", lev ), clust=T)
    return(hmlist)
  })
  names(hmlevs)<-levs

  hm1 <- hmall$hm
  asvorder <- hm1$tree_col$labels[hm1$tree_col$order]
  varorder <- hm1$tree_row$labels[hm1$tree_row$order]

  newcormat <- hmall$cormat[varorder, asvorder]
  newptext <- hmall$cor_ptext[varorder, asvorder]
  annot_rows <- rep("All samples", nrow(newcormat))
  vspaces <- c()
  for(ll in levs){
    vspaces <- c(vspaces, nrow(newcormat))
    newcormat <- rbind(newcormat, hmlevs[[ll]]$cormat[varorder, asvorder])
    newptext <- rbind(newptext, hmlevs[[ll]]$cor_ptext[varorder, asvorder])
    annot_rows <- c(annot_rows,rep(paste0(var2plot, "-", ll), nrow(hmlevs[[ll]]$cor_ptext)))
  }

  # select_asvs: plot these ASVs. Overrides pval
  # pval: plot ASVs with any correlation >  than this threshold

  oname <- paste0(outdir, "/", name, ".pdf")
  fontsize_row = 16 - nrow(newcormat) / 15
  fontsize_col = 16 - ncol(newcormat) / 15

  annot_rows <- data.frame( Var=rownames(newcormat), Condition = annot_rows)
  rownames(annot_rows) <- paste(annot_rows$Condition, annot_rows$Var, sep="-")
  rownames(newcormat) <- rownames(annot_rows)
  rownames(newptext) <- rownames(annot_rows)

  cond_colors <- default_annot_colors[1:length(unique(annot_rows$Condition))]
  names(cond_colors) <- unique(annot_rows$Condition)
  annot_colors <- list(Condition=cond_colors)
  if(! is.null(annot_asvs)){annot_colors[["LFC"]] <- c("Up"="red", "-"="white", "Down"="blue")}

  #annot_asvs <- annot_asvs[colnames(newcormat), ]

  tmp <-1
  while(!is.null(tmp)){tmp <- tryCatch(dev.off(), error=function(x){})}
  hm <- pheatmap(newcormat, display_numbers = newptext,
                 #filename=oname,
                 #width = 20, height = 12,
                 fontsize_row = fontsize_row,
                 fontsize_number = 16,
                 fontsize_col = fontsize_col,
                 cluster_rows = F, cluster_cols = F,
                 gaps_row = vspaces,
                 annotation_row = annot_rows %>% dplyr::select(Condition),
                 labels_row = annot_rows$Var,
                 annotation_colors = annot_colors,
                 annotation_col = annot_asvs
  )
  #tmp <- tryCatch(dev.off(), error=function(x){})
  pdf(oname, width = 20, height = 10)
  print(hm)
  tmp <- dev.off()
  return(hm)
  #return(list(hm=hm, hm_allonly=hmall, hmbygroup=hmlevs))
}
