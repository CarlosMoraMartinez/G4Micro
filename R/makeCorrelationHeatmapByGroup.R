#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mat1 PARAM_DESCRIPTION
#' @param mat2 PARAM_DESCRIPTION
#' @param metadata PARAM_DESCRIPTION
#' @param var2plot PARAM_DESCRIPTION, Default: 'Psoriasis'
#' @param varwithnames PARAM_DESCRIPTION, Default: 'sampleID'
#' @param cormethod PARAM_DESCRIPTION, Default: 'pearson'
#' @param pval PARAM_DESCRIPTION, Default: 0.05
#' @param select_asvs PARAM_DESCRIPTION, Default: c()
#' @param outdir PARAM_DESCRIPTION, Default: ''
#' @param name PARAM_DESCRIPTION, Default: 'corrHeatmap'
#' @param clust PARAM_DESCRIPTION, Default: T
#' @param order_rows PARAM_DESCRIPTION, Default: c()
#' @param order_cols PARAM_DESCRIPTION, Default: c()
#' @param annot_asvs PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{select}}
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
