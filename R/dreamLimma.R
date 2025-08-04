dreamLimma <- function(phobj, opt, variables =  c("Condition"), individual="pacienteID", outdir = "", name = "dream"){
  formula <- paste0("~ 0 +", paste(variables, sep=" + ", collapse=" + ")) %>%
    as.formula
  formula_ind <- paste0("~ ", paste(variables, sep=" + ", collapse=" + "), " + (1|", individual, ")") %>%
    as.formula
  count_mat <- otu_table(phobj) %>% data.frame %>% as.matrix()
  count_mat <- count_mat[rowSums(count_mat  > opt$mincount) > opt$minsampleswithcount, ]
  colnames(count_mat)<- gsub("^X", "", colnames(count_mat))
  # Standard usage of limma/voom
  genes = DGEList( count_mat )
  genes = calcNormFactors( genes )

  metadata <- metadata[colnames(count_mat), ]
  #design = model.matrix(formula_ind, metadata)
  vobj_tmp = voomWithDreamWeights(genes, formula_ind, metadata)

  contrvec <- paste(variables[1] , levels(metadata[, variables[1]]), sep="")
  L = getContrast( vobj_tmp, formula_ind, metadata,
                   contrvec[2])

  #fitmm_1 = dream(vobj_tmp, formula_ind, metadata)
  fitmm_dream = dream(vobj_tmp, formula_ind, metadata, L)
  #result_dream1 <- variancePartition::topTable(fitmm_1, sort.by = "P", n = Inf, coef=contrvec[2])
  result_dream <- variancePartition::topTable(fitmm_dream, sort.by = "P", n = Inf, coef=contrvec[2])

  write.table(result_dream, file=paste0(outdir, "/", name, "_contrast.tsv"))
  #write.table(result_dream1, file=paste0(outdir, "/", name, ".tsv"))

  reslist <- list(result=result_dream, vobj_tmp=vobj_tmp, fitmm=fitmm_dream, L=L)
  oname <- paste0(outdir, "/" ,name, ".RData")
  save(reslist, file=oname)
  return(reslist)
}
