#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param phobj PARAM_DESCRIPTION
#' @param seed PARAM_DESCRIPTION, Default: 123
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[phyloseq]{tax_glom}}, \code{\link[phyloseq]{subset_samples}}, \code{\link[phyloseq]{otu_table}}
#'  \code{\link[HMP]{Xdc.sevsample}}
#' @rdname testHMP
#' @export 
#' @importFrom phyloseq tax_glom subset_samples otu_table
#' @importFrom HMP Xdc.sevsample
testHMP <- function(phobj, seed=123){
  ps_phylum <- phyloseq::tax_glom(pre_phyloseq_filt, "Phylum")
  controls <- phyloseq::subset_samples(ps_phylum, Condition == "Control")
  cf <- phyloseq::subset_samples(ps_phylum, Condition == "Depression")
  #Output OTU tables
  control_otu <- data.frame(phyloseq::otu_table(controls))
  cf_otu <- data.frame(phyloseq::otu_table(cf))

  group_data <- list(control_otu, cf_otu)
  (xdc <- HMP::Xdc.sevsample(group_data))
}
