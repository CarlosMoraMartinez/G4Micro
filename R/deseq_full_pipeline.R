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
#' @param plot_all Whether to make plots for all contrasts or not. Default, TRUE.
#' @param doPoscounts Logical; whether to force the use of the `poscounts` normalization method. Default: FALSE.
#' @return A list containing DESeq2 result data frames and associated objects produced by \code{getDeseqResults}.
#' @param deseqname A character string specifying the subdirectory name for storing DESeq2 results.
#'   Default: \code{"DeSEQ2/"}.
#' @param vars2heatmap Character vector of variable names (from sample metadata) to annotate
#'   in generated heatmaps. Default: \code{c("Condition")}.
#' @param max_hm_h Numeric; maximum height of the heatmaps (in inches) when saving PDF output.
#'   Passed to \code{makeHeatmap}. Default: \code{16}.
#' @param max_hm_w Max heatmap width. If 0, will be automatically determined depending on number of cols. Default: `7
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
deseq_full_pipeline <- function(phobj, name, vars2deseq, opt, plot_all=TRUE, doPoscounts=FALSE,
                                deseqname = "DeSEQ2/", vars2heatmap=c("Condition"), max_hm_h=16, max_hm_w=7){
  restaurar <- restauraropt_mk(opt)
  opt <- restaurar(opt)
  if(!dir.exists(paste0(opt$out, deseqname))) dir.create(paste0(opt$out, deseqname))
  outdir <- paste0(opt$out, deseqname, name, "/")
  opt$out <- outdir
  if(!dir.exists(opt$out)) dir.create(opt$out)
  if(opt$minfreq > 0){
    opt$minsampleswithcount <- opt$minfreq*nsamples(phobj)
    cat("Minfreq: ", opt$minfreq, ", setting minsampleswithcount to ", opt$minsampleswithcount)
  }
  dearesults <- getDeseqResults(phobj, opt, name, variables = vars2deseq, doPoscounts=doPoscounts)

  #list2env(dearesults, envir = environment())
  if(plot_all){
    make_all_maplots(dearesults$all_contrasts, opt)
  }else{
    make_all_maplots(list(dearesults$all_contrasts[[1]]), opt)
  }

  df2plot <-  if(! nrow(dearesults$vst_counts_df)){dearesults$norm_counts_df}else{dearesults$vst_counts_df}
  if(plot_all){
    make_all_heatmaps(dearesults$all_contrasts, df2plot,
                      sample_data(phobj) %>% data.frame,
                      vars2heatmap, dearesults$dds, opt,
                      max_hm_h=max_hm_h, max_hm_w=max_hm_w)
  }

  ## Global heatmaps
  if(dearesults$all_combos_done & length(dearesults$all_contrasts) > 1){
    cat("All contrasts TRUE, intersecting Taxon list")
    taxalist_praw <- map(dearesults$all_contrasts, \(x){
      x$resdf %>% dplyr::filter(pvalue < opt$pval) %>% pull(taxon)
    }) %>% unlist %>% unique
    taxalist_padj <- map(dearesults$all_contrasts, \(x){
      x$resdf %>% dplyr::filter(padj < opt$pval) %>% pull(taxon)
    }) %>% unlist %>% unique
  }else{
    taxalist_praw = taxalist_padj = c()
  }

  if(length(vars2heatmap) == 0) vars2heatmap <- vars2deseq

  tryCatch(makeHeatmap(dearesults$resdf, dearesults$dds, df2plot, vars2heatmap,
                       opt, name = paste0(name, "diff_ab_heatmap_rawpval.pdf"),
                       logscale = F, ptype="pvalue", trim_values = TRUE, taxalist=taxalist_praw,
                       max_hm_h=max_hm_h),
           error=\(x) cat("Error makeHeatmap praw"))
  tryCatch(makeHeatmap(dearesults$resdf, dearesults$dds, df2plot, vars2heatmap,
                       opt, name = paste0(name, "diff_ab_heatmap_adjpval.pdf"),
                       logscale = F, ptype="padj", trim_values = TRUE, taxalist=taxalist_padj,
                       max_hm_h=max_hm_h),
           error=\(x) cat("Error makeHeatmap padj"))

  dfcorr <- df2plot %>% column_to_rownames("gene") %>%
    as.matrix %>% cor %>%
    as.data.frame %>% rownames_to_column("gene")
  tryCatch(makeHeatmap(dearesults$resdf, dearesults$dds, dfcorr, vars2heatmap,
                       opt, name = paste0(name, "corr_heatmap_rawpval.pdf"),
                       logscale = F, ptype="pvalue", trim_values = TRUE, taxalist=dfcorr$gene,
                       italics_rownames = FALSE, check_taxa = FALSE,
                       max_hm_h=max_hm_h),
           error=\(x) cat("Error makeHeatmap praw"))
  tryCatch(makeHeatmap(dearesults$resdf, dearesults$dds, dfcorr, vars2heatmap,
                       opt, name = paste0(name, "corr_heatmap_adjpval.pdf"),
                       logscale = F, ptype="padj", trim_values = TRUE, taxalist=dfcorr$gene,
                       italics_rownames = FALSE, check_taxa = FALSE,
                       max_hm_h=max_hm_h),
           error=\(x) cat("Error makeHeatmap padj"))

  tryCatch({
    pdf(paste0(opt$out, name, "_DispEsts.pdf"))
    plotDispEsts(dds, CV=T , ylim = c(1e-6, 1e1))
    dev.off()
  }, error=\(x)cat("Error plotDispEsts\n"))


  opt <- restaurar(opt)
  return(dearesults)
}
