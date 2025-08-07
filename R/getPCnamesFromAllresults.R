
#' @title Generate Named Principal Component Labels with Explained Variance Percentages
#' @description
#' Extracts principal component (PC) names and attaches the percentage of variance explained
#' by each PC based on PCA results stored within a nested results list.
#'
#' This function is useful for labeling PCs in plots or tables with their variance contribution.
#'
#' @param phname Character. Name of the phyloseq object or dataset key to extract results from.
#' @param all_model_results Nested list containing all model and PCA results organized by phyloseq object names and condition names.
#' @param get_pcnames_from Character. Name of the list element containing the PC names to extract. Default is \code{"padj_taxa_res"}.
#' @param pca_name Character. Name of the list element containing PCA objects with variance summaries. Default is \code{"padj_taxa_pcas"}.
#' @param varname Character. Variable name inside the PCA list to use for accessing the PCA object. Default is \code{"Condition"}.
#'
#' @return Named character vector of PC names, where each name is appended with the percentage of variance explained in parentheses.
#'
#' @details
#' The function accesses the PC names stored in \code{all_model_results[[phname]][[get_pcnames_from]]$varnames} and
#' the PCA importance (variance explained) summary stored in
#' \code{all_model_results[[phname]][[pca_name]][[varname]]$pca}. It calculates the percentage of variance explained
#' by each PC and appends this to the PC name to generate labels like "PC1 (45.2%)".
#'
#' These labeled PCs can be used to annotate plots or tables for easier interpretation of PCA results.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname getPCnamesFromAllresults
#' @export
getPCnamesFromAllresults <- function(phname, all_model_results,
                                     get_pcnames_from="padj_taxa_res",
                                     pca_name="padj_taxa_pcas",
                                     varname="Condition"){
  PCs <- all_model_results[[phname]][[get_pcnames_from]]$varnames
  sumPCA <- summary(all_model_results[[phname]][[pca_name]][[varname]]$pca)$importance
  PC_names <- paste(PCs, ' (', round(100*sumPCA[2, PCs], 1),'%)', sep="")
  PCs_newnames <- PCs
  names(PCs_newnames) <- PC_names
  return(PCs_newnames)
}
