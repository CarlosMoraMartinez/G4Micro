
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param phname PARAM_DESCRIPTION
#' @param all_model_results PARAM_DESCRIPTION
#' @param get_pcnames_from PARAM_DESCRIPTION, Default: 'padj_taxa_res'
#' @param pca_name PARAM_DESCRIPTION, Default: 'padj_taxa_pcas'
#' @param varname PARAM_DESCRIPTION, Default: 'Condition'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
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
