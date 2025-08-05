#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param phobj PARAM_DESCRIPTION
#' @param opt PARAM_DESCRIPTION
#' @param variables PARAM_DESCRIPTION, Default: c("Condition")
#' @param individual PARAM_DESCRIPTION, Default: 'pacienteID'
#' @param outdir PARAM_DESCRIPTION, Default: ''
#' @param name PARAM_DESCRIPTION, Default: 'dream'
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[variancePartition]{topTable}}
#' @rdname dreamLimma_old
#' @export 
#' @importFrom variancePartition topTable
dreamLimma_old <- function(phobj, opt, variables =  c("Condition"), individual="pacienteID", outdir = "", name = "dream"){
  formula <- paste0("~ ", paste(variables, sep=" + ", collapse=" + ")) %>%
    as.formula
  formula_ind <- paste0("~ ", paste(variables, sep=" + ", collapse=" + "), " + (1|", individual, ")") %>%
    as.formula

  count_mat <- otu_table(phobj) %>% data.frame %>% as.matrix()
  count_mat <- count_mat[rowSums(count_mat  > opt$mincount) > opt$minsampleswithcount, ]
  colnames(count_mat)<- gsub("^X", "", colnames(count_mat))
  # Standard usage of limma/voom
  genes = DGEList( count_mat )
  genes = calcNormFactors( genes )

  #snames <- gsub("^X", "", colnames(count_mat))
  metadata <- metadata[colnames(count_mat), ]
  design = model.matrix(formula, metadata)
  vobj_tmp = voom( genes, design, plot=FALSE)
  #colnames(vobj_tmp) <- gsub("^X", "", colnames(vobj_tmp))
  dupcor <- duplicateCorrelation(vobj_tmp,design,block=metadata[, individual])
  vobj = voom( genes, design, plot=FALSE, block=metadata[, individual], correlation=dupcor$consensus)
  ##DREAM
  cl <- makeCluster(4)
  registerDoParallel(cl)

  # Get the contrast matrix for the hypothesis test
  L = getContrast( dupcor, formula_ind, metadata,
                   colnames(vobj_tmp$design)[2])
  fitmm = dream( vobj, formula_ind, metadata, L)
  fitmm_ebayes = eBayes( fitmm )


  # get results
  result <- variancePartition::topTable(fitmm_ebayes, sort.by = "P", n = Inf)

  write.table(result, file=paste0(outdir, "/", name, ".tsv"))
  reslist <- list(result=result, vobj_tmp=vobj_tmp, fitmm=fitmm, L=L)
  oname <- paste0(outdir, name, ".RData")
  save(reslist, file=oname)
  return(reslist)
}
