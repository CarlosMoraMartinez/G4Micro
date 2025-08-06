#' @title Run All Taxon Abundance Plotting Functions
#'
#' @description
#' Executes a complete pipeline to generate multiple visualizations and summary tables of microbial taxonomic abundances at different levels (Phylum, Genus, ASV), saving plots and tables to the specified output directory.
#'
#' @param phobj A \code{phyloseq} object containing microbial abundance and sample metadata.
#' @param interestvar A character string specifying the sample variable of interest (e.g., condition or group) for grouping and comparisons.
#' @param outdir A directory path where the output files (PDFs and TSVs) will be saved.
#' @param phname A short name or prefix to include in all output file names.
#' @param levs A named vector of factor levels to use for the grouping variable (e.g., to enforce group order).
#' @param tops A numeric vector specifying how many top genera to visualize. Default: \code{c(5, 10, 15, 20)}.
#'
#' @return
#' A named list containing:
#' \itemize{
#'   \item \code{prev_vs_abn}: A ggplot object showing prevalence vs abundance.
#'   \item \code{sp_col_by_gen}: A ggplot object of ASVs colored by Genus.
#'   \item \code{prevalence_tab_genus}: A data frame of genus-level prevalence.
#'   \item \code{prevalence_tab_phylum}: A data frame of phylum-level prevalence.
#'   \item \code{top_genus_plots}: A list of ggplot objects showing relative abundances by genus for different top-N thresholds.
#'   \item \code{phylum_abundance_tests}: A data frame with test results comparing phylum abundances between groups.
#'   \item \code{phylum_boxplots}: A ggplot object with phylum-level boxplots.
#'   \item \code{phylum_rel_bars}: A ggplot object with relative abundance barplots at the phylum level.
#' }
#'
#' @details
#' This function runs several plotting and analysis functions:
#' \itemize{
#'   \item Barplots and boxplots for phylum-level abundances.
#'   \item Prevalence and abundance summaries at phylum and genus levels.
#'   \item Genus-level barplots for multiple top-N values.
#'   \item ASV-level barplots colored by genus.
#'   \item Prevalence vs abundance scatterplot.
#' }
#' All outputs are written to the given directory using filenames that incorporate \code{phname}.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  plotAbundanceFullPipeline(
#'     phobj = ps_data,
#'     interestvar = "Group",
#'     outdir = "results/",
#'     phname = "gut_microbiome",
#'     levs = c("Control", "TreatmentA", "TreatmentB")
#'   )
#'  }
#' }
#' @seealso
#'  \code{\link{plotRelativeAbnBarsPhylum}},
#'  \code{\link{plotPhylumBoxplots}},
#'  \code{\link{getPhylumTests}},
#'  \code{\link{getRelAbundancesByPhylumAndVariable}},
#'  \code{\link{plotRelativeAbnBarsGenus}},
#'  \code{\link{getRelAbundancesByGenusAndVariable}},
#'  \code{\link{plotRelativeAbnBars_Fantaxtic}},
#'  \code{\link{plotPrevalenceVsAbundance}}
#' @rdname plotAbundanceFullPipeline
#' @export
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
