#' @title Differential expression analysis with variancePartition's dream and limma-voom
#' @description
#' Runs a repeated-measures differential expression analysis on count data
#' (e.g., microbiome or RNA-seq) using the `dream` method from the
#' \pkg{variancePartition} package, which extends `limma-voom` to handle
#' random effects such as subject IDs.
#'
#' @param phobj A \pkg{phyloseq} object containing an OTU/ASV count table
#' and associated sample metadata.
#' @param opt A list of options, must contain:
#'   \itemize{
#'     \item \code{mincount}: minimum count per feature to be kept.
#'     \item \code{minsampleswithcount}: minimum number of samples with counts above \code{mincount} to retain a feature.
#'   }
#' @param variables Character vector specifying the fixed-effect variables
#' to include in the model. Default: \code{c("Condition")}.
#' @param individual Name of the metadata column identifying individuals
#' for the random effect term. Default: \code{"pacienteID"}.
#' @param outdir Output directory for saving results. Default: \code{""}.
#' @param name Base name for output files. Default: \code{"dream"}.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{result}: topTable output from \code{dream}.
#'     \item \code{vobj_tmp}: voom-transformed object with precision weights.
#'     \item \code{fitmm}: fitted model object from \code{dream}.
#'     \item \code{L}: contrast matrix used in the analysis.
#'   }
#'
#' @details
#' The workflow follows these steps:
#' \enumerate{
#'   \item Extract counts from the \pkg{phyloseq} object and filter low-abundance features.
#'   \item Create an \pkg{edgeR} \code{DGEList} and normalize library sizes with \code{calcNormFactors}.
#'   \item Apply \code{voomWithDreamWeights} to transform counts to logCPM and estimate observation-level weights.
#'   \item Specify contrasts with \code{getContrast} and fit the model with \code{dream}, including random effects.
#'   \item Extract top results with \code{topTable} and save outputs as TSV and RData.
#' }
#' The \code{dream} method allows modeling of repeated measures (e.g., multiple samples from the same subject) by
#' adding a random effect term \code{(1|individual)} in the model formula.
#'
#' For details on the methodology, see:
#' \url{https://bioconductor.org/packages/release/bioc/vignettes/variancePartition/inst/doc/variancePartition.html}
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   library(phyloseq)
#'   library(variancePartition)
#'   ps <- phyloseq(otu_table(matrix(rpois(200, lambda=5), ncol=10), taxa_are_rows=TRUE),
#'                  sample_data(data.frame(Condition=rep(c("A","B"), each=5),
#'                                         pacienteID=rep(1:5, each=2))))
#'   opt <- list(mincount=5, minsampleswithcount=2)
#'   results <- dreamLimma(ps, opt, variables="Condition", individual="pacienteID",
#'                         outdir="results", name="dream_test")
#' }
#' }
#'
#' @seealso
#'  \code{\link[variancePartition]{dream}},
#'  \code{\link[variancePartition]{topTable}},
#'  \code{\link[variancePartition]{voomWithDreamWeights}},
#'  \code{\link[limma]{voom}},
#'  \code{\link[edgeR]{calcNormFactors}}
#' @rdname dreamLimma
#' @export
#' @importFrom phyloseq otu_table
#' @importFrom magrittr %>%
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom variancePartition voomWithDreamWeights getContrast dream topTable
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
