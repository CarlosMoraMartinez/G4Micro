#'  @title Differential Expression Analysis Using Limma-Voom with Duplicate Correlation
#' @description
#' Performs differential expression analysis on count data from a phyloseq object
#' using the limma-voom pipeline. This method models the mean-variance relationship of
#' RNA-seq or similar count data by applying voom transformation,
#' allowing the use of linear models for differential expression.
#'
#' To account for repeated measures or batch effects from individuals, this function
#' estimates within-individual correlations using the \code{duplicateCorrelation} method.
#' This approach fits a linear mixed model by incorporating a consensus correlation
#' structure as a random effect, improving the accuracy of variance estimation and
#' statistical testing.
#'
#' The model formula is constructed dynamically from specified covariates (e.g., experimental conditions).
#' Differential contrasts (e.g., "ConditionAfter - ConditionBefore") are specified and tested
#' with moderated t-statistics from empirical Bayes shrinkage.
#'
#' The function writes results to disk as a tab-separated values file and saves the entire analysis
#' environment for downstream inspection or visualization.
#'
#' @param phobj A phyloseq object containing OTU count data and sample metadata.
#' @param opt A list of options controlling filtering thresholds:
#'   \itemize{
#'     \item \code{mincount} Minimum count threshold for filtering OTUs.
#'     \item \code{minsampleswithcount} Minimum number of samples that must pass the count threshold.
#'   }
#' @param variables Character vector of covariate names to include in the design formula (default: \code{c("Condition")}).
#' @param individual Character string naming the individual/sample grouping variable to account for repeated measures (default: \code{"pacienteID"}).
#' @param outdir Character string specifying output directory path for saving results (default: \code{""}).
#' @param name Character string specifying the base name for output files (default: \code{"basic_limma"}).
#' @return A list containing:
#'   \itemize{
#'     \item \code{result}: A data frame with the differential expression test results sorted by p-value.
#'     \item \code{vobj}: The \code{voom}-transformed expression set after duplicateCorrelation adjustment.
#'     \item \code{vobj_tmp}: Initial \code{voom} object before duplicateCorrelation adjustment.
#'     \item \code{fitmm}: The linear model fit object after contrast and empirical Bayes steps.
#'     \item \code{contrast}: The contrast matrix used for hypothesis testing.
#'   }
#' @details
#' This function implements the limma-voom pipeline tailored for count data with repeated measures.
#'
#' Steps include:
#' \enumerate{
#'   \item Filtering low-count OTUs based on \code{opt$mincount} and \code{opt$minsampleswithcount}.
#'   \item Creating a design matrix from the specified variables.
#'   \item Normalizing counts with TMM (Trimmed Mean of M-values) via \code{calcNormFactors}.
#'   \item Applying the voom transformation to estimate the mean-variance trend and compute weights.
#'   \item Estimating the correlation between repeated samples using \code{duplicateCorrelation}.
#'   \item Re-running voom with the correlation and fitting a linear model with blocking for individuals.
#'   \item Defining contrasts for comparison (here, "ConditionAfter - ConditionBefore").
#'   \item Computing moderated t-statistics with empirical Bayes shrinkage.
#' }
#'
#' Output files include a TSV file with DE results and an RData file with all key objects.
#'
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'
#'   # Load or prepare phyloseq object with count data and sample metadata
#'   phobj <- phyloseq::your_phyloseq_object
#'
#'   # Define options for filtering
#'   opt <- list(mincount=5, minsampleswithcount=3)
#'
#'   # Run basic Limma differential expression analysis
#'   res <- basicLimma(phobj, opt, variables = c("Condition"), individual = "SampleID", outdir = "results/", name = "DE_Analysis")
#' }
#' }
#' @seealso
#' \code{\link[limma]{voom}},
#' \code{\link[limma]{duplicateCorrelation}},
#' \code{\link[limma]{lmFit}},
#' \code{\link[limma]{makeContrasts}},
#' \code{\link[limma]{contrasts.fit}},
#' \code{\link[limma]{eBayes}},
#' \code{\link[edgeR]{DGEList}},
#' \code{\link[edgeR]{calcNormFactors}},
#' \code{\link[phyloseq]{otu_table}}
#' @rdname basicLimma
#' @export
#' @importFrom magrittr %>%
#' @importFrom phyloseq otu_table
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom limma voom duplicateCorrelation lmFit makeContrasts contrasts.fit eBayes topTable

basicLimma <- function(phobj, opt, variables =  c("Condition"), individual="pacienteID", outdir="", name="basic_limma"){
  formula <- paste0("~ 0 +", paste(variables, sep=" + ", collapse=" + ")) %>%
    as.formula
  #formula_ind <- paste0("~ ", paste(variables, sep=" + ", collapse=" + "), " + (1|", individual, ")") %>%
  # as.formula
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
  # in order to compute more accurate weights
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
