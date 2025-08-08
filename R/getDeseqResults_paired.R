#' @title Differential Expression Analysis for Paired Samples with DESeq2
#' @description
#' Performs differential expression analysis on paired samples using DESeq2,
#' applying a likelihood ratio test (LRT) to test the effect of specified
#' variables while controlling for individual effects. The function handles filtering,
#' pseudocount addition, multiple shrinkage methods for log fold changes,
#' and returns various DESeq2 results and normalized counts.
#'
#' @param phobj A phyloseq object containing OTU counts and sample metadata.
#' @param opt A list of options including filtering thresholds and output paths.
#'   Expected elements include:
#'   \itemize{
#'     \item \code{mincount}: minimum count threshold for filtering features.
#'     \item \code{minsampleswithcount}: minimum number of samples required to meet the count threshold.
#'     \item \code{minfreq}: minimum frequency threshold (not directly used here).
#'     \item \code{fc}: fold-change threshold used for log fold change shrinkage.
#'     \item \code{out}: output directory path for saving results.
#'   }
#' @param name Character string suffix to append to output file names. Default is an empty string.
#' @param variables Character vector of variable names to include as fixed effects in the model. Default is \code{c("Condition")}.
#' @param individual Character string naming the individual identifier variable used for pairing. Default is \code{"pacienteID"}.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{dds}: DESeqDataSet object from the full model.
#'   \item \code{dds_null}: DESeqDataSet object from the reduced model.
#'   \item \code{raw_counts}: matrix of raw counts.
#'   \item \code{res}: DESeq2 results from the likelihood ratio test.
#'   \item \code{res_null}: results from the null model.
#'   \item \code{resLFC}: log fold change shrinkage results (normal).
#'   \item \code{resLFC_ape}: LFC shrinkage with "apeglm" method.
#'   \item \code{resLFC_ashr}: LFC shrinkage with "ashr" method.
#'   \item \code{resdf}, \code{resdf_ape}, \code{resdf_shr}: data frames written by \code{defWriteDEAResults}.
#'   \item \code{raw_df}: data frame of raw counts written by \code{defWriteMatAsDF}.
#'   \item \code{norm_counts}: matrix of normalized counts.
#'   \item \code{norm_counts_df}: data frame of normalized counts.
#'   \item \code{vstds}: variance stabilized transformed counts matrix (or NULL if error).
#'   \item \code{vst_counts_df}: data frame of VST counts.
#'   \item \code{options}: list of options including filtering parameters and whether pseudocounts were added.
#' }
#'
#' @details
#' This function fits a paired design model for differential expression using DESeq2,
#' controlling for individual effects. It filters low-count features based on
#' specified thresholds, applies pseudocounts if necessary, and performs likelihood
#' ratio tests between full and reduced models. Multiple methods for log fold
#' change shrinkage are applied and results are saved to files.
#'
#' The function requires \code{defWriteDEAResults} and \code{defWriteMatAsDF} helper functions
#' for writing results to disk.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Assuming 'ps' is a phyloseq object and 'options' is a list with required parameters:
#'   results <- getDeseqResults_paired(
#'     phobj = ps,
#'     opt = list(mincount=10, minsampleswithcount=3, fc=1.5, out="./results/"),
#'     name = "experiment1",
#'     variables = c("Treatment"),
#'     individual = "PatientID"
#'   )
#'   print(results$resdf)
#' }
#' }
#'
#' @rdname getDeseqResults_paired
#' @export
#' @importFrom phyloseq otu_table sample_data
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq counts results lfcShrink resultsNames varianceStabilizingTransformation
#' @importFrom dplyr %>%
getDeseqResults_paired <- function(phobj, opt, name="", variables = c("Condition"), individual = "pacienteID"){
  formula <- paste0("~ ", paste(variables, sep=" + ", collapse=" + "), " + ",
                    paste(individual, sep=" + ", collapse=" + ")) %>%
    as.formula
  reduced_formula = paste0("~ ", paste(individual, sep=" + ", collapse=" + ")) %>%
    as.formula

  counts_table <- otu_table(phobj) %>% data.frame %>% as.matrix
  colnames(counts_table) <- gsub("^X", "", colnames(counts_table))
  design <- sample_data(phobj) %>% data.frame
  design <- design[colnames(counts_table) ,]

  dds <- DESeqDataSetFromMatrix(countData = counts_table, design = formula, colData = design)

  #dds <- phyloseq_to_deseq2(phobj, design= formula)
  raw_counts <- counts(dds)
  if(opt$minsampleswithcount == 0){
    dds <- dds[rowSums(counts(dds)) >= opt$mincount,]
  }else{
    dds <- dds[rowSums(counts(dds) >= opt$mincount) >= opt$minsampleswithcount,]
  }
  filt_counts <- counts(dds)

  ##Add pseudocount if necessary
  anyNonZero <- raw_counts %>% apply(MAR=1, all) %>% any
  if(!anyNonZero){
    do_poscounts = TRUE
    dds_null <- DESeq(dds, betaPrior = F, sfType = "poscounts")
    dds <- DESeq(dds, betaPrior = F, sfType = "poscounts", reduced=reduced_formula, test="LRT")
  }else{
    do_poscounts = FALSE
    dds_null <- DESeq(dds, betaPrior = F)
    dds <- DESeq(dds, betaPrior = F, reduced=reduced_formula, test="LRT")
  }

  levelscond <- levels(design[, variables[1]])
  res <- tryCatch(results(dds, contrast=c(variables[1], levelscond[2], levelscond[1])),
                  error=results(dds, contrast=list(paste0(variables[1], "_" ,levelscond[2], "_vs_",levelscond[1]) ))
  )
  res_null <- tryCatch(results(dds_null, contrast=c(variables[1], levelscond[2], levelscond[1])),
                       error=results(dds_null, contrast=list(paste0(variables[1], "_" ,levelscond[2], "_vs_",levelscond[1]) ))
  )

  resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="normal", lfcThreshold = log2(opt$fc)) #apeglm gives weird results
  resLFC_null <- lfcShrink(dds_null, coef=resultsNames(dds)[2], type="normal", lfcThreshold = log2(opt$fc)) #apeglm gives weird results
  resLFC_ape <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm", lfcThreshold = log2(opt$fc)) #apeglm gives weird results
  resLFC_ashr <- lfcShrink(dds, coef=resultsNames(dds)[2], type="ashr", lfcThreshold = log2(opt$fc)) #apeglm gives weird results
  resdf <- defWriteDEAResults(res, resLFC, opt, paste0(name, "_", "DEA_results_shrinkNormal.tsv"))
  resLFC_null <- defWriteDEAResults(res, resLFC, opt, paste0(name, "_", "DEA_results_shrinkNormal_againstNullModel.tsv"))
  resdf_ape <- defWriteDEAResults(res, resLFC_ape, opt, paste0(name, "_", "DEA_results_shrinkApe.tsv"))
  resdf_ashr <- defWriteDEAResults(res, resLFC_ashr, opt, paste0(name, "_", "DEA_results_shrinkAshr.tsv"))

  # Write raw counts
  rawc_df <- defWriteMatAsDF(raw_counts, opt, paste0(name, "_", "raw_counts.tsv") )
  #filtx_df <- defWriteMatAsDF(filt_counts, opt, "raw_counts_filtered.tsv")

  # Normalized counts
  #  dds <- estimateSizeFactors(dds)##  already done
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
    "dds_null"=dds_null,
    "raw_counts"=raw_counts,
    "res"=res,
    "res_null"=res_null,
    "resLFC"=resLFC,
    "resLFC_ape"=resLFC_ape,
    "resLFC_ashr"=resLFC_ashr,
    "resdf"=resdf,
    "resdf_ape"=resdf_ape,
    "resdf_shr"=resdf_ashr,
    "raw_df" = rawc_df,
    "norm_counts"=norm_counts,
    "norm_counts_df"=norm_counts_df,
    "vstds"=vstds,
    "vst_counts_df"=vst_counts_df,
    "options"= list(mincount=opt$mincount, minsampleswithcount=opt$minsampleswithcount,
                    minfreq=opt$minfreq, poscount=do_poscounts)
  )
  save(file = paste0(opt$out, "DESEQ2_all_results_", name, ".R"), all_results_list)
  return(all_results_list)
}
