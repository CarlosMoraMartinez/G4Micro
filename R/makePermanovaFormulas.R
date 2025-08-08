#' @title Run Multiple PERMANOVAs Using Custom Formulas
#' @description
#' Applies multiple PERMANOVA tests using a list of formulas on the sample metadata of a phyloseq object.
#'
#' @param phobj A \code{phyloseq} object containing OTU/ASV table and sample metadata.
#' @param formulas A character vector of formulas (e.g., \code{"~ Condition"}, \code{"~ Condition + Time"}).
#' @param dist_method Dissimilarity method to use (e.g., \code{"bray"}, \code{"jaccard"}). Passed to \code{\link[phyloseq]{distance}}. Default: "bray".
#' @param seed Random seed for reproducibility of adonis2. Default: 123.
#' @param outname File name to save the summary table as a TSV file. Default: "permanovas_mult.tsv".
#' @param adonisby String passed to \code{\link[vegan]{adonis2}} in the \code{by} argument. Default: "terms"
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{res}: A data frame with PERMANOVA results (including p-values and adjusted p-values).
#'   \item \code{modelos}: A named list of \code{adonis2} model objects for each formula.
#' }
#'
#' @details
#' This function runs PERMANOVA tests using the \code{adonis2} function from \pkg{vegan}.
#' Sample metadata column names are normalized by removing accents and replacing spaces with underscores.
#' A summary table of results is saved as a TSV file. P-values are adjusted using the Benjamini-Hochberg method.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  library(phyloseq)
#'   data(GlobalPatterns)
#'   makePermanovaFormulas(GlobalPatterns,
#'                          formulas = c("~ SampleType", "~ SampleType + Country"))
#'  }
#' }
#' @seealso
#'  \code{\link[phyloseq]{distance}},
#'  \code{\link[vegan]{adonis2}},
#'  \code{\link[dplyr]{arrange}},
#'  \code{\link[stringi]{stri_trans_general}},
#'  \code{\link[readr]{write_tsv}},
#'  \code{\link[stats]{p.adjust}}
#'  \code{\link{adonis2table}}
#'
#' @rdname makePermanovaFormulas
#' @export
#' @importFrom phyloseq distance
#' @importFrom dplyr arrange
#' @importFrom stringi stri_trans_general
#' @importFrom vegan adonis2
makePermanovaFormulas <- function(phobj, formulas, dist_method = "bray", seed = 123,
                                  outname = "permanovas_mult.tsv", adonisby="terms"){
  ## From https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova

  meta <- sample_data(phobj)
  braydist <- phyloseq::distance(phobj, method = dist_method)
  sampledf <- data.frame(sample_data(phobj))
  names(sampledf) <- stri_trans_general(str = names(sampledf), id = "Latin-ASCII") %>%
    gsub(" ", "_", .)

  modelos <- list()
  res <- data.frame()
  for(form in formulas){
    set.seed(seed)
    mod1 <- adonis2(as.formula(form),
                    data = sampledf,
                    na.action=na.exclude,
                    by = adonisby)
    newtab <- adonis2table(mod1, adonisby = adonisby) %>%
      dplyr::mutate(model = form,
                    padj_bymodel = p.adjust(P, method = "BH"))
    res <-rbind(res, newtab)
    modelos[[form]] <- mod1
  }
  res <- res %>% dplyr::select(model, everything())
  res$padj <- p.adjust(res$P, method="BH")
  write_tsv(res, file=outname)
  return(list(res=res, modelos=modelos))

}
