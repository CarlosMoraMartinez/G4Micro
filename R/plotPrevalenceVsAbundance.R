plotPrevalenceVsAbundance <- function(phobj, outname="phylumBarplot.pdf", height=10, width=12){
  pre_prevalence <- getRelAbundanceTab(phobj)
  g1 <- ggplot(pre_prevalence, aes(TotalAbundance, relative_prevalence, color = Phylum)) +
    # Agregamos una línea para nuestro umbral
    geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
    scale_x_log10() +  xlab("Total Abundance") +
    ylab("Prevalence/num. samples") +
    facet_wrap(~Phylum) +
    theme_pubclean() +
    guides(color = FALSE)
  ggsave(outname, g1, width=width, height = height)
  return(g1)
}
