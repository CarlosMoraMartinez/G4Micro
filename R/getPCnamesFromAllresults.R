
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
