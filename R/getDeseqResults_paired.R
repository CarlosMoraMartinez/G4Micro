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
