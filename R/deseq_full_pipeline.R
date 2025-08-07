#' @title Full DESeq2 Differential Expression Analysis Pipeline
#' @description
#' Runs a comprehensive differential expression analysis pipeline on a phyloseq object,
#' including filtering, DESeq2 result generation, annotation, visualization (MA plots, volcano plots, heatmaps),
#' and saving of annotated results.
#'
#' @param phobj A phyloseq object containing count data and taxonomy information.
#' @param name A character string used as a prefix for output directories and files.
#' @param vars2deseq Character vector of variable names to use in the DESeq2 design formula.
#' @param opt A list of options controlling output directories, filtering thresholds, p-value cutoffs, and other parameters.
#'
#' @return A list containing DESeq2 result data frames and associated objects produced by \code{getDeseqResults}.
#'
#' @details
#' This function creates output directories if needed, computes the minimum number of samples with counts based on frequency thresholds,
#' executes differential expression analysis through \code{getDeseqResults()}, annotates results with genus-level taxonomy,
#' saves annotated tables as TSV files, and generates diagnostic and visualization plots such as MA plots, volcano plots, and heatmaps.
#' It also performs intersection of taxa lists across contrasts when applicable.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   opt <- list(out = "results/", minfreq = 0.01, pval = 0.05)
#'   vars <- c("Condition", "Batch")
#'   res <- deseq_full_pipeline(physeq_object, "MyExperiment", vars, opt)
#' }
#' }
#' @seealso
#'  \code{\link{getDeseqResults}}, \code{\link{make_maplot}}, \code{\link{make_volcano}},
#'  \code{\link{makeHeatmap}}, \code{\link{getSummaryTablesDeseq}}
#'  \code{\link[DESeq2]{DESeq}}, \code{\link[DESeq2]{results}}, \code{\link[DESeq2]{plotDispEsts}},
#'  \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{select}}, \code{\link[dplyr]{arrange}}, \code{\link[dplyr]{filter}}
#' @rdname deseq_full_pipeline
#' @export
#' @importFrom dplyr mutate select arrange filter
#' @importFrom purrr map
#' @importFrom readr write_tsv
#' @importFrom phyloseq tax_table
#' @importFrom DESeq2 plotDispEsts
deseq_full_pipeline <- function(phobj, name, vars2deseq, opt){
  if(!dir.exists(paste0(opt$out, "DeSEQ2"))) dir.create(paste0(opt$out, "DeSEQ2"))
  outdir <- paste0(opt$out, "DeSEQ2/", name, "/")
  opt$reserva <- opt$out
  opt$out <- outdir
  if(!dir.exists(opt$out)) dir.create(opt$out)
  if(opt$minfreq > 0){
    opt$minsampleswithcount <- opt$minfreq*nsamples(phobj)
    cat("Minfreq: ", opt$minfreq, ", setting minsampleswithcount to ", opt$minsampleswithcount)
  }
  dearesults <- getDeseqResults(phobj, opt, name, variables = vars2deseq)

  list2env(dearesults, envir = environment())
  tax2annot <- tax_table(phobj)
  resdf_annot <- resdf %>%
    dplyr::mutate(Genus = data.frame(tax2annot[taxon, "Genus"])$Genus) %>%
    dplyr::select(Genus, everything()) %>%
    dplyr::arrange(pvalue)
  write_tsv(resdf_annot, file=paste0(opt$out, "/DEA_annot.tsv"))
  # resdf_annot %>% filter(pvalue < 0.05) %>% dplyr::select(Genus, taxon) %>%
  #   kable(caption="Differentially abundant ASVs at adjusted p-value < 0.05")
  tryCatch(make_maplot(res, opt, paste0(name, "_MAPlot-rawFC.pdf")),  error=\(x)cat("Error make_maplot"))
  tryCatch(make_maplot(resLFC, opt,  paste0(name, "_MAPlot-rawFC.pdf")),  error=\(x)cat("Error make_maplot"))
  tryCatch(make_maplot(resLFC_ape, opt,  paste0(name, "_MAPlot-rawFC-ape.pdf")),  error=\(x)cat("Error make_maplot"))
  tryCatch(make_maplot(resLFC_ashr, opt,  paste0(name, "_MAPlot-rawFC-ashr.pdf")),  error=\(x)cat("Error make_maplot"))
  pdf(paste0(opt$out, name, "_DispEsts.pdf"))
  plotDispEsts(dds, CV=T , ylim = c(1e-6, 1e1))
  dev.off()
  rtabs <- getSummaryTablesDeseq(res, opt)
  tryCatch(make_volcano(resLFC, opt, paste0(name, "volcano_rawfc_rawpval.pdf"), "pvalue"), error=\(x)cat("Error make_volcano"))
  tryCatch(make_volcano(res, opt, paste0(name, "volcano_rawfc_adjpval.pdf"), "padj"),  error=\(x)cat("Error make_volcano"))
  df2plot <- if(! nrow(vst_counts_df)){norm_counts_df}else{vst_counts_df}

  if(all_combos_done & length(all_contrasts) > 1){
    cat("All contrasts TRUE, intersecting Taxon list")
    taxalist_praw <- map(all_contrasts, \(x){
      x$resdf %>% dplyr::filter(pvalue < opt$pval) %>% pull(taxon)
    }) %>% unlist %>% unique
    taxalist_padj <- map(all_contrasts, \(x){
      x$resdf %>% dplyr::filter(padj < opt$pval) %>% pull(taxon)
    }) %>% unlist %>% unique
  }else{
    taxalist_praw = taxalist_padj = c()
  }
  tryCatch(makeHeatmap(resdf, dds, df2plot, vars2deseq,
                       opt, name = paste0(name, "diff_ab_heatmap_rawpval.pdf"),
                       logscale = F, ptype="pvalue", trim_values = TRUE, taxalist=taxalist_praw),
           error=\(x) cat("Error makeHeatmap praw"))
  tryCatch(makeHeatmap(resdf, dds, df2plot, vars2deseq,
                       opt, name = paste0(name, "diff_ab_heatmap_adjpval.pdf"),
                       logscale = F, ptype="padj", trim_values = TRUE, taxalist=taxalist_padj),
           error=\(x) cat("Error makeHeatmap padj"))
  opt$out <- opt$reserva
  return(dearesults)
}
