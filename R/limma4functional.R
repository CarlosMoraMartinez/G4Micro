#' @title Differential Pathway Abundance Analysis Using limma
#' @description Performs differential abundance testing on pathway-level quantitative data using linear modeling and empirical Bayes moderation via the limma package.
#'
#' @param df2 A data frame of quantitative pathway abundances. The first column must be named "Pathway" and contain pathway identifiers. Subsequent columns represent samples, with column names matching sample IDs in `metad2`.
#' @param metad2 A data frame of sample metadata, with at least columns `sampleID` (matching `df2` columns) and the grouping variable specified in `interestvar`. May also include covariates.
#' @param interestvar Character string specifying the main factor variable in `metad2` to test differences across (default is `"Condition"`).
#' @param covars Character vector of covariate column names in `metad2` to adjust for in the model (default is empty vector, meaning no covariates).
#' @param form Optional custom model formula as a string or formula object for design matrix construction. If `NULL` (default), a formula is automatically created from `interestvar` and `covars`.
#' @param make_all_contrasts Logical, whether to compute all possible contrasts between factor levels rather than only the main contrast (default `FALSE`).
#'
#' @return A data frame of differential abundance statistics for the main contrast if `make_all_contrasts = FALSE`. If `TRUE`, returns a named list of such data frames, one per contrast.
#'
#' @details
#' This function analyzes functional pathway abundance data (or other similar expression-like matrices) to identify differentially abundant pathways between experimental groups.
#'
#' The procedure follows these steps:
#' \enumerate{
#'   \item Metadata rows with missing values (`NA`) in specified covariates are filtered out to ensure clean modeling.
#'   \item The abundance data frame (`df2`) is converted to a numeric matrix with pathways as rows and samples as columns. Samples are reordered to match the metadata order.
#'   \item A log-transformation (`log(x + 1)`) is applied to stabilize variance and reduce skewness.
#'   \item The main variable of interest (`interestvar`) in `metad2` is coerced to a factor.
#'   \item If no formula is provided, a design formula is built automatically, modeling the main variable without an intercept and including covariates if specified.
#'   \item A design matrix is generated from the formula and metadata, encoding the experimental design.
#'   \item A linear model is fit to the log-transformed data using limma's `lmFit`.
#'   \item A contrast matrix is constructed to compare two groups: specifically, `Depression` versus `Control` (hardcoded, modify as needed).
#'   \item The contrast is applied, and empirical Bayes moderation (`eBayes`) is performed to improve variance estimation.
#'   \item A top table of statistics is extracted, including log2 fold changes, average expression, moderated t-statistics, raw and adjusted p-values, and log-odds scores.
#'   \item If requested, the function computes differential statistics for all other contrasts between groups, returning a list of results.
#' }
#'
#' **Note:** The contrast definitions are currently hardcoded and may require modification for custom group names or additional contrasts.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Example input data preparation
#'   df2 <- data.frame(Pathway = c("PWY1", "PWY2"),
#'                     Sample1 = c(10, 5),
#'                     Sample2 = c(20, 15),
#'                     Sample3 = c(8, 3),
#'                     Sample4 = c(1, 13))
#'   metad2 <- data.frame(sampleID = c("Sample1", "Sample2", "Sample3", "Sample4"),
#'                       Condition = c("Control", "Depression", "Control", "Depression"))
#'
#'   # Run differential abundance analysis
#'   results <- limma4functional(df2, metad2, interestvar = "Condition")
#'   head(results)
#' }
#' }
#'
#' @seealso
#' \code{\link[dplyr]{filter}}, \code{\link[limma]{lmFit}}, \code{\link[limma]{makeContrasts}}, \code{\link[limma]{eBayes}}, \code{\link[limma]{topTable}}
#'
#' @rdname limma4functional
#' @export
#' @importFrom dplyr filter
#' @import limma
limma4functional <- function(df2, metad2, interestvar = "Condition", covars=c(), form=NULL,
                            make_all_contrasts = FALSE){
  for(covar in covars){
    metad2 <- metad2 %>% dplyr::filter(!is.na(metad2[, covar]))
  }
  rownames(df2) <- NULL
  expr <- df2 %>% column_to_rownames("Pathway") %>% as.matrix
  expr <- expr[, metad2$sampleID]
  expr <- log(expr+1)
  metad2[, interestvar] <- as.factor(metad2[, interestvar])

  if(is.null(form)){
    if(length(covars) > 0){
      form <- paste("~0 ", interestvar, paste(covars, collapse = ' + '),
                    sep = ' + ', collapse=" + ") %>%
        as.formula
    }else{
      form <- paste0("~0 + ", interestvar) %>% as.formula()
    }
  }else{
    form <- as.formula(form)
  }
  print(form)
  design <- model.matrix(form, metad2)
  colnames(design) <- gsub(interestvar, "", colnames(design), perl=F)
  if(grepl("\\*", as.character(form)[2])) colnames(design) <- make.names(colnames(design))

  #colnames(design) <- gsub("metad2\\$Condition", "", colnames(design), perl=F)
  fit <- lmFit(expr, design)
  cont.matrix <- makeContrasts(case_vs_control = Depression - Control,
                               levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  tt <- topTable(fit2, n=Inf, sort.by = "P", adjust.method = "BH")
  names(tt) <- c("log2FoldChange", "AveExpr", "t", "pvalue", "padj", "B")
  dim(tt)

  if(make_all_contrasts){
    othercontrasts <- colnames(design)[! colnames(design) %in% unique(metad2[, interestvar])]
    other_results <- list(main_contraast=tt)
    for(othc in othercontrasts){
      cont.matrix_oth <- makeContrasts(case_vs_control = IMC,
                                       levels = design)
      fit2_oth <- contrasts.fit(fit, cont.matrix_oth)
      fit2_oth <- eBayes(fit2_oth)
      tt_oth <- topTable(fit2_oth, n=Inf, sort.by = "P", adjust.method = "BH")
      names(tt_oth) <- c("log2FoldChange", "AveExpr", "t", "pvalue", "padj", "B")
      other_results[[othc]] <- tt_oth

    }
    return(other_results)
  }
  return(tt)
}
