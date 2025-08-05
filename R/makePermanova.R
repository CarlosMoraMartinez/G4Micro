#' @title FUNCTION_TITLE
#' @description
#' This function performs PERMANOVA, as implemented in the vegan::adonis2 function,
#' for all the variables included in the metadata of a phyloseq object,
#' except for those included in the 'exclude_vars' argument. Based on:
#' @param phobj Phyloseq Object
#' @param dist_method Distance method to pass to the phyloseq::distance function, Default: 'bray'
#' @param seed Random seed, Default: 123
#' @param exclude_vars PARAM_DESCRIPTION, Default: c("sampleID")
#' @param outname PARAM_DESCRIPTION, Default: 'permanovas.tsv'
#' @param disp_permutations PARAM_DESCRIPTION, Default: 1000
#' @return A DataFrame object with the PERMANOVA results for each variable in the phyloseq object metadata
#' @references
#' From https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html#permanova
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  makePermanova(phobj, dist_method = "bray",
#'                  seed = 123,
#'                  exclude_vars = c("sampleID", "patient_name"),
#'                  outname = "my_output_name")
#'  }
#' }
#' @seealso
#'  \code{\link[phyloseq]{distance}}, \code{\link[phyloseq]{prune_samples}}
#'  \code{\link[dplyr]{arrange}}
#' @rdname makePermanova
#' @export
#' @importFrom vegan permutest adonis2 betadisper capscale anova.cca
#' @importFrom phyloseq distance prune_samples
#' @importFrom dplyr arrange
#' @importFrom stringi stri_trans_general
makePermanova <- function(phobj, dist_method = "bray", seed = 123,
                          exclude_vars = c("sampleID"), outname = "permanovas.tsv", disp_permutations=1000){

  meta <- sample_data(phobj)
  braydist <- phyloseq::distance(phobj, method = dist_method)
  sampledf <- data.frame(sample_data(phobj))
  names(sampledf) <- stri_trans_general(str = names(sampledf), id = "Latin-ASCII") %>%
    gsub(" ", "_", .)

  vars2test <- names(sampledf)[! names(sampledf) %in% exclude_vars]

  res <- data.frame()
  for(var in vars2test){
    set.seed(seed)
    form <- paste0("braydist ~ ", var) %>% as.formula

    if(length(unique(sampledf[!is.na(sampledf[, var]) , var])) > 1){
      # Adonis test
      #cat(var, "\n")
      # get original data frame
      sampledf <- data.frame(sample_data(phobj))
      names(sampledf) <- stri_trans_general(str = names(sampledf), id = "Latin-ASCII") %>%
        gsub(" ", "_", .)

      s2use <- sampledf$sampleID[!is.na(sampledf[, var]) ]
      phobj_filt <- phyloseq::prune_samples(s2use, phobj)

      braydist <- phyloseq::distance(phobj_filt, method = dist_method)
      sampledf <- data.frame(sample_data(phobj_filt))
      names(sampledf) <- stri_trans_general(str = names(sampledf), id = "Latin-ASCII") %>%
        gsub(" ", "_", .)

      mod1 <- adonis2(form, data = sampledf, na.action=na.exclude, permutations = disp_permutations)

      perm_disp <- tryCatch({
        bd <- betadisper(braydist, sampledf[, var])
        permutest(bd, permutations = disp_permutations)
      },error=\(x) NULL )

      ancap <- tryCatch({
        cap <- capscale(form, data = sampledf)
        ancap <- anova.cca(cap, permutations = disp_permutations)
      },error=\(x) NULL )

      res <-rbind(res, adonis2table(mod1, perm_disp, ancap, var))
    }
  }
  res <- res %>%
    dplyr::arrange(P)
  res$padj <- p.adjust(res$P, method="BH")
  if("perm_disp_P" %in% names(res)){
    res$perm_disp_Padj <- p.adjust(res$perm_disp_P, method="BH")
  }
  if("capscaleanova_P" %in% names(res)){
    res$capscaleanova_Padj <- p.adjust(res$capscaleanova_P, method="BH")
  }
  write_tsv(res, file=outname)
  return(res)

}
