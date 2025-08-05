#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param phobj PARAM_DESCRIPTION
#' @param opt PARAM_DESCRIPTION
#' @param variables PARAM_DESCRIPTION, Default: c("Condition")
#' @param individual PARAM_DESCRIPTION, Default: 'pacienteID'
#' @param outdir PARAM_DESCRIPTION, Default: ''
#' @param name PARAM_DESCRIPTION, Default: 'basic_limma'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname basicLimma
#' @export 
basicLimma <- function(phobj, opt, variables =  c("Condition"), individual="pacienteID", outdir="", name="basic_limma"){
  formula <- paste0("~ 0 +", paste(variables, sep=" + ", collapse=" + ")) %>%
    as.formula
  #formula_ind <- paste0("~ ", paste(variables, sep=" + ", collapse=" + "), " + (1|", individual, ")") %>%
  as.formula
  count_mat <- otu_table(phobj) %>% data.frame %>% as.matrix()
  count_mat <- count_mat[rowSums(count_mat  > opt$mincount) > opt$minsampleswithcount, ]

  # Standard usage of limma/voom
  genes = DGEList( count_mat )
  genes = calcNormFactors( genes )

  snames <- gsub("^X", "", colnames(count_mat))
  metadata <- metadata[snames, ]
  design = model.matrix(formula, metadata)
  vobj_tmp = voom( genes, design, plot=FALSE)

  # apply duplicateCorrelation
  dupcor <- duplicateCorrelation(vobj_tmp,design,block=metadata[, individual])

  # run voom considering the duplicateCorrelation results
  # in order to compute more accurate precision weights
  # Otherwise, use the results from the first voom run
  vobj = voom( genes, design, plot=FALSE, block=metadata[, individual], correlation=dupcor$consensus)

  # But this step uses only the genome-wide average for the random effect
  fitDupCor <- lmFit(vobj, design, block=metadata$Individual, correlation=dupcor$consensus)

  # Fit Empirical Bayes for moderated t-statistics
  #fitDupCor_bay <- eBayes( fitDupCor )
  contr <- makeContrasts(ConditionAfter -ConditionBefore, levels =  colnames(coef(fitDupCor)))
  tmp <- contrasts.fit(fitDupCor, contr)
  tmp <- eBayes(tmp)
  result <- topTable(tmp, sort.by = "P", n = Inf)

  write.table(result, file=paste0(outdir, name, ".tsv"))
  reslist <- list(result=result, vobj=vobj, vobj_tmp = vobj_tmp, fitmm=tmp, contrast=contr)
  oname <- paste0(outdir, "/", name, ".RData")
  save(reslist, file=oname)
  return(reslist)
}
