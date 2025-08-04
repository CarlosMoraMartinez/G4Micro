getDeseqResults_fromFormula <- function(phobj, opt, name="", formula= "~ Condition"){
  variables <- gsub("~", "", formula) %>% gsub(" ", "", .) %>% strsplit(split="\\+|\\*") %>% unlist
  dds <- phyloseq_to_deseq2(phobj, design= as.formula(formula))
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
