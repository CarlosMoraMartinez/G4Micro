#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param df2 PARAM_DESCRIPTION
#' @param metad2 PARAM_DESCRIPTION
#' @param interestvar PARAM_DESCRIPTION, Default: 'Condition'
#' @param covars PARAM_DESCRIPTION, Default: c()
#' @param form PARAM_DESCRIPTION, Default: NULL
#' @param make_all_contrasts PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[dplyr]{filter}}
#' @rdname imma4functional
#' @export 
#' @importFrom dplyr filter
imma4functional <- function(df2, metad2, interestvar = "Condition", covars=c(), form=NULL,
                            make_all_contrasts = FALSE){
  library(limma)
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
