#' @title Test for Differences in Microbiome Composition Using HMP Xdc.sevsample at the Phylum level
#' @description Performs a multi-sample permutation test from the HMP package to assess differences in microbiome composition between two groups at the Phylum level.
#' @param phobj A `phyloseq` object containing microbiome count data with sample metadata including a "Condition" variable that has at least "Control" and "Depression" groups.
#' @param seed Integer seed for reproducibility of the permutation test. Default: 123
#' @return The test statistic and p-value from the \code{HMP::Xdc.sevsample} function assessing differences between groups.
#' @details
#' The function first agglomerates the phyloseq object at the Phylum level,
#' then subsets samples into two groups based on the "Condition" variable.
#' It extracts the OTU tables for each group and applies the \code{Xdc.sevsample} test
#' from the HMP package, which tests for equality of community composition distributions.
#'
#' The input assumes the `Condition` variable has values "Control" and "Depression".
#' @examples
#' \dontrun{
#' if(interactive()){
#'   library(phyloseq)
#'   library(HMP)
#'   data(GlobalPatterns) # Example dataset, replace with your own phyloseq object
#'
#'   # Add a dummy Condition variable for demonstration
#'   sample_data(GlobalPatterns)$Condition <- rep(c("Control", "Depression"), length.out = nsamples(GlobalPatterns))
#'
#'   test_result <- testHMP(GlobalPatterns, seed = 123)
#'   print(test_result)
#' }
#' }
#' @seealso
#'  \code{\link[phyloseq]{tax_glom}}, \code{\link[phyloseq]{subset_samples}}, \code{\link[phyloseq]{otu_table}}
#'  \code{\link[HMP]{Xdc.sevsample}}
#' @rdname testHMP
#' @export
#' @importFrom phyloseq tax_glom subset_samples otu_table
#' @importFrom HMP Xdc.sevsample
testHMP <- function(phobj, seed=123){
  ps_phylum <- phyloseq::tax_glom(phobj, "Phylum")
  controls <- phyloseq::subset_samples(ps_phylum, Condition == "Control")
  cf <- phyloseq::subset_samples(ps_phylum, Condition == "Depression")
  #Output OTU tables
  control_otu <- data.frame(phyloseq::otu_table(controls) %>% t)
  cf_otu <- data.frame(phyloseq::otu_table(cf) %>% t)

  group_data <- list(control_otu, cf_otu)
  xdc <- HMP::Xdc.sevsample(group_data)
  return(xdc)
}
