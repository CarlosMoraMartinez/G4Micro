plotAbundanceFullPipeline <- function(phobj, interestvar, outdir, phname, levs, tops=c(5,10,15,20)){
  oname <- paste0(outdir, phname, "_phylumBarplot.pdf")
  g1 <- plotRelativeAbnBarsPhylum(phobj, interestvar, oname)

  oname <- paste0(outdir, phname, "_phylumBoxplots.pdf")
  g2 <-plotPhylumBoxplots(phobj, interestvar, oname, paired=F)

  oname <- paste0(outdir, phname, "_phylumBivariateTests.tsv")
  tests_tab <- getPhylumTests(phobj, interestvar, oname, paired=F)

  outname <- paste0(outdir, phname, "_PrevalenceOfPhyla.tsv")
  df_prevalence <- getRelAbundancesByPhylumAndVariable(phobj, interestvar, outname,
                                                       oldlevs=levs)
  togenplots <- list()
  for(topn in tops){
    oname <- paste0(outdir, phname, "_relAbund_byGenus_top", as.character(topn), ".pdf")
    togenplots[[as.character(topn)]] <- plotRelativeAbnBarsGenus(phobj, interestvar, topn=topn, outname=oname,
                                                                 oldlevs=levs,
                                                                 width = 20, height = 8)
  }

  oname <- paste0(outdir, phname, "_prevalence_byGenus.tsv")
  pre_prevalence <- getRelAbundancesByGenusAndVariable(phobj, interestvar, oname,
                                                       oldlevs=levs)

  oname <- paste0(outdir, phname, "_relAbund_byASV_ColByGenus_top", as.character(n_species), ".pdf")
  g3 <- plotRelativeAbnBars_Fantaxtic(phobj, interestvar, topn = n_species, tax_level="Genus", outname = oname)

  oname <- paste0(outdir, phname, "_PrevVsAbund.pdf")
  g4 <- plotPrevalenceVsAbundance(phobj, oname)

  resall <- list(
    prev_vs_abn=g4,
    sp_col_by_gen=g3,
    prevalence_tab_genus=pre_prevalence,
    prevalence_tab_phylum=df_prevalence,
    top_genus_plots=togenplots,
    phylum_abundance_tests = tests_tab,
    phylum_boxplots = g2,
    phylum_rel_bars = g1
  )
  return(resall)
}
