plotPhylumBoxplots <- function(phobj, var="Condition", outname="phylumBarplot.pdf", height=10, width=8, paired=F){
  ## Total abundance by phylum, apparently by summing over all ASVs in phylum
  ps_phylum <- phyloseq::tax_glom(phobj, "Phylum")
  phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]

  sample_data(ps_phylum)[, var] <- unlist(sample_data(ps_phylum)[, var]) %>% as.character
  comps <- combn(unique(unlist(sample_data(ps_phylum)[, var])), 2, simplify = F)
  signif_codes <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))

  g1 <- phyloseq::psmelt(ps_phylum) %>%
    ggplot(data = ., aes_string(x = var, y = "Abundance")) +
    facet_wrap(~ OTU, scales = "free")+
    geom_boxplot(outlier.shape  = NA) +
    geom_jitter(aes(color = OTU), height = 0, width = .2) +
    stat_compare_means(method="wilcox.test", comparisons = comps,
                       symnum.args = signif_codes, paired=paired) +
    labs(x = var, y = "Abundance\n") +
    theme_pubclean() +
    guides(color = FALSE)
  if(paired){
    g1 <- g1 + geom_line(aes(group = pacienteID), color = "gray20", linetype=1, size=0.1, alpha=0.25)
  }
  ggsave(outname, g1, width=width, height = height)
  return(g1)
}
