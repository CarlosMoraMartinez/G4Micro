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
