#' @title Run DESeq2 Differential Abundance Analysis
#' @description This function prepares and runs DESeq2 differential expression analysis from a `phyloseq` object,
#' handling design formulas and multiple contrasts.
#' @param phobj A `phyloseq` object containing count data and sample metadata.
#' @param opt A list of options including output path (`out`),
#' filtering thresholds (`mincount`, `minsampleswithcount`, and optionally `minfreq` (which, if specified, overrides `minsampleswithcount`)).
#' @param name Optional string used as a prefix for output files. Default: ''.
#' @param variables A character vector with variables to use in the design formula (e.g., conditions to compare). Default: c("Condition").
#' @param formula A custom formula to override automatic construction from `variables`. Default: NULL.
#' @param doposcounts Logical; whether to force the use of the `poscounts` normalization method. Default: FALSE.
#' @return A named list with components:
#' \itemize{
#'   \item \code{dds}: DESeq2 object.
#'   \item \code{raw_counts}, \code{norm_counts}: Raw and normalized count matrices.
#'   \item \code{res}, \code{resLFC}, \code{resLFC_ape}, \code{resLFC_ashr}: Main DESeq2 results and shrinkage estimates.
#'   \item \code{resdf}, \code{resdf_ape}, \code{resdf_shr}: Results as data frames.
#'   \item \code{raw_df}, \code{norm_counts_df}, \code{vst_counts_df}: Count matrices written to files.
#'   \item \code{vstds}: Variance-stabilized transformed data.
#'   \item \code{all_contrasts}: List of all contrast results.
#'   \item \code{all_combos_done}: Logical flag. Indicates whether all the possible contrasts were performed or not.
#'   \item \code{options}: A list summarizing filtering parameters and whether pseudocounts were used.
#' }
#' @details The function supports both categorical and numeric variables.
#' It automatically generates pairwise contrasts between factor levels or handles numeric
#' predictors using likelihood ratio tests. If all count rows contain zeros for some samples, `poscounts` normalization
#' is applied even when the \code{doposcounts} argument is set to \code{FALSE}.
#' Intermediate and final results are written to files using a naming pattern defined by `name`.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  opt <- list(out="results/", mincount=10, minsampleswithcount=3, minfreq=0.001)
#'  getDeseqResults(phobj, opt, name="comparison1", variables=c("Treatment"))
#'  }
#' }
#' @seealso
#'   \code{\link{getDeseqContrastFromNumerical}},
#'   \code{\link{getDeseqContrastFromCategorical}},
#'   \code{\link{defWriteMatAsDF}},
#'   \code{\link[DESeq2]{DESeq}},
#'   \code{\link[DESeq2]{counts}},
#'   \code{\link[DESeq2]{results}},
#'   \code{\link[DESeq2]{resultsNames}},
#'   \code{\link[DESeq2]{varianceStabilizingTransformation}},
#'   \code{\link[phyloseq]{phyloseq_to_deseq2}}
#' @rdname getDeseqResults
#' @export
#' @importFrom phyloseq phyloseq_to_deseq2
#' @importFrom DESeq2 DESeq counts results resultsNames varianceStabilizingTransformation
getDeseqResults <- function(phobj, opt, name="", variables = c("Condition"), formula=NULL, doposcounts=FALSE){

  if(is.null(formula)){
    formula <- paste0("~ ", paste(variables, sep=" + ", collapse=" + ")) %>%
      as.formula
  }else{
    variables <- gsub("~", "", formula) %>%
      gsub(" ", "", .) %>%
      strsplit(split="\\+|\\*") %>%
      unlist
    formula <- as.formula(formula)
  }
  dds <- phyloseq_to_deseq2(phobj, design= formula)
  raw_counts <- counts(dds)
  if(opt$minsampleswithcount == 0){
    dds <- dds[rowSums(counts(dds)) >= opt$mincount,]
  }else{
    dds <- dds[rowSums(counts(dds) >= opt$mincount) >= opt$minsampleswithcount,]
  }
  filt_counts <- counts(dds)

  ##Add pseudocount if necessary
  anyNonZero <- raw_counts %>% apply(MAR=1, all) %>% any
  if(!anyNonZero | doposcounts){
    do_poscounts = TRUE
    dds <- DESeq(dds, betaPrior = F, sfType = "poscounts")
  }else{
    do_poscounts = FALSE
    dds <- DESeq(dds, betaPrior = F)
  }

  write_file(paste(resultsNames(dds), collapse="\t" ), file=paste0(opt$out, name, "_", "DEA_resultsNames.tsv"))

  design <- dds@colData %>% as.data.frame()
  all_combos_done <- TRUE

  all_combins <- map(variables, \(x){
    if(is.numeric(design[, x])){
      return(list(c(x, "NUMERIC")))
    }
    levs <- levels(design[, x] %>% unlist)
    combins <- lapply(combn(1:length(levs), 2, simplify = F), \(y)c(x, levs[y]))
  }) %>% flatten

  all_contrasts <- map(all_combins, \(lev_combin){
    if(lev_combin[2] == "NUMERIC"){
      getDeseqContrastFromNumerical(dds, lev_combin[1], opt, name)
    }else{
      getDeseqContrastFromCategorical(dds, lev_combin, opt, name)
    }
  })
  names(all_contrasts) <- lapply(all_combins, \(x)ifelse(x[2]=="NUMERIC", x[1], paste0(x[1],'_', x[3], '_vs_', x[2]) %>% gsub(" ", ".", .)))
  all_combos_done <- TRUE
  # Write raw counts
  rawc_df <- defWriteMatAsDF(raw_counts, opt, paste0(name, "_", "raw_counts.tsv") )
  #filtx_df <- defWriteMatAsDF(filt_counts, opt, "raw_counts_filtered.tsv")

  # Normalized counts
  norm_counts <- counts(dds, normalized = T)
  norm_counts_df <- defWriteMatAsDF(norm_counts, opt, paste0(name, "_", "norm_counts.tsv") )

  tryCatch({
    vstds <- varianceStabilizingTransformation(raw_counts, blind=F)
    vst_counts_df <- defWriteMatAsDF(vstds, opt, paste0(name, "_", "vst_counts.tsv") )
  }, error= function(x){
    vstds <<- NULL
    vst_counts_df <<- data.frame()
  })

  all_results_list <- list(
    "dds"=dds,
    "raw_counts"=raw_counts,
    "all_contrasts"=all_contrasts,
    "res"=all_contrasts[[1]]$res,
    "resLFC"=all_contrasts[[1]]$resLFC,
    "resLFC_ape"=all_contrasts[[1]]$resLFC_ape,
    "resLFC_ashr"=all_contrasts[[1]]$resLFC_ashr,
    "resdf"=all_contrasts[[1]]$resdf,
    "resdf_ape"=all_contrasts[[1]]$resdf_ape,
    "resdf_shr"=all_contrasts[[1]]$resdf_ashr,
    "raw_df" = rawc_df,
    "norm_counts"=norm_counts,
    "norm_counts_df"=norm_counts_df,
    "vstds"=vstds,
    "vst_counts_df"=vst_counts_df,
    "all_combos_done"=all_combos_done,
    "options"= list(mincount=opt$mincount, minsampleswithcount=opt$minsampleswithcount,
                    minfreq=opt$minfreq, poscount=do_poscounts)
  )
  save(file = paste0(opt$out, "DESEQ2_all_results_", name, ".R"), all_results_list)
  return(all_results_list)
}
