#' @title Plot Relative Abundance Barplot of Top Species Colored by Genus
#' @description Generates a stacked bar plot of the relative abundances of the top N most abundant ASVs (species-level),
#'              coloring each species by its genus. The samples are optionally clustered or ranked for better visualization.
#' @param phobj A phyloseq object containing OTU counts and taxonomy data.
#' @param variable Sample metadata variable used for grouping/faceting (not fully implemented in this snippet). Default: "Condition"
#' @param topn Number of top ASVs (species) to include in the plot, sorted by total abundance. Default: 15
#' @param outname File name to save the resulting plot as PDF. Default: "SpeciescolByGenusBarplot.pdf"
#' @param height Height of the saved plot in inches. Default: 8
#' @param width Width of the saved plot in inches. Default: 12
#' @param docluster Logical; whether to cluster samples and genera to reorder the plot axes. Default: TRUE
#' @return A ggplot2 object representing the stacked barplot.
#' @details
#' The function first normalizes counts to relative abundances per sample,
#' then selects the top N ASVs by total abundance.
#' It assigns unique genus names (with suffixes if needed) to distinguish genera with duplicates.
#' Samples and genera are optionally reordered by hierarchical clustering.
#' Finally, a stacked bar plot is created and saved to file.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   library(phyloseq)
#'   data(GlobalPatterns)
#'
#'   # Plot top 10 species colored by genus, save to a PDF
#'   plotRelativeAbnBarsSpecies_ColByGenus(GlobalPatterns, variable="SampleType", topn=10,
#'                                        outname="example_barplot.pdf", docluster=FALSE)
#' }
#' }
#' @seealso
#'  \code{\link[phyloseq]{transform_sample_counts}}
#'  \code{\link[dplyr]{filter}}, \code{\link[dplyr]{arrange}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{mutate_if}}
#'  \code{\link[tidyr]{gather}}
#' @rdname plotRelativeAbnBarsSpecies_ColByGenus
#' @export
#' @importFrom phyloseq transform_sample_counts
#' @importFrom dplyr filter arrange mutate mutate_if
#' @importFrom tidyr gather
plotRelativeAbnBarsSpecies_ColByGenus <- function(phobj, variable="Condition", topn = 15,
                                                  outname="SpeciescolByGenusBarplot.pdf", height=8,
                                                  width=12, docluster=T){
  library(RColorBrewer)
  ps_rel_abund = phyloseq::transform_sample_counts(phobj, function(x){x / sum(x)})
  aa <- otu_table(ps_rel_abund)
  bb <- rowSums(aa) %>% sort(decreasing = T)
  top_asv <- names(bb)[1:topn]
  taxdata <- tax_table(phobj) %>% data.frame() %>%
    rownames_to_column("ASV") %>%
    dplyr::filter(ASV %in% top_asv) %>%
    dplyr::arrange(Genus)
  taxdata$Genus2 <- sapply(1:nrow(taxdata), FUN=function(i, lista){
    aux <- table(lista[1:i])
    aux2 <- table(lista)[names(aux)]
    aux <- ifelse(aux2 > 1, as.character(aux), "")
    paste0(lista[i], as.character(aux[lista[i]] ))
  }, taxdata$Genus)


  otus_gen <- data.frame(otu_table(phobj) ) %>%
    rownames_to_column("ASV") %>%
    dplyr::mutate(is_top = ASV %in% taxdata$ASV,
                  Genus = taxdata$Genus2[match(ASV, taxdata$ASV)]
    ) %>%
    dplyr::mutate(Genus = ifelse(is.na(Genus), "Other", Genus)) %>%
    # group_by(Genus) %>%
    # summarise_if(is.numeric, sum) %>%
    filter(Genus != "Other") %>%
    dplyr::mutate(
      Phylum = taxdata$Phylum[match(ASV, taxdata$ASV)],
      Family = taxdata$Family[match(ASV, taxdata$ASV)],
    )

  otus_gen_gen_rel <- otus_gen %>%
    dplyr::mutate_if(is.numeric, function(x)x/sum(x))

  otus_df <- otus_gen_gen_rel %>%
    tidyr::gather(key="sampleID", value="Abundance", -ASV, -Genus, -Phylum, -Family, -is_top) %>%
    dplyr::mutate(sampleID = gsub("^X", "", sampleID))


  df_merged <- merge(otus_df, data.frame(sample_data(phobj)), by = "sampleID")

  if(! ocluster){
    ranked_genera <- otus_gen_gen_rel %>%
      dplyr::mutate(mean_rel = rowMeans(across(where(is.numeric)))) %>%
      top_n(topn, mean_rel) %>% dplyr::arrange(mean_rel) #desc(mean_rel)
    most_prev_genus <- ranked_genera$Genus[nrow(ranked_genera)]

    ranked_samples <- df_merged %>%
      filter(Genus == most_prev_genus) %>%
      dplyr::arrange(Abundance) %>%
      pull(sampleID)
  }else{
    # dmat <- otus_gen_gen_rel
    # rownames(dmat) <- dmat$Genus
    # dmat <- dmat %>%
    #   select_if(is.numeric) %>% dist
    # hcl <- hclust(dmat)
    # ranked_genera <-hcl$labels[hcl$order]
    ranked_genera <- taxdata$Genus2

    dmat <- otus_gen_gen_rel %>% select_if(is.numeric) %>% t %>% dist
    hcl <- hclust(dmat)
    ranked_samples <- gsub("^X", "", hcl$labels[hcl$order])
  }

  df_merged <- df_merged %>% dplyr::mutate(sampleID = factor(sampleID, levels=ranked_samples),
                                           Genus = factor(Genus, levels = ranked_genera)
  )

  #cols <- brewer.pal(n=nrow(ranked_genera), name="Set2") # No more colors than those in palette
  mycolors <- colorRampPalette(brewer.pal(8, "Spectral"))(length(ranked_genera))
  g1 <-ggplot(df_merged, aes(x=sampleID, y=Abundance, color = Genus,
                             fill = Genus)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "", y = "Relative Abundance\n") +
    facet_wrap(facet_form, scales = "free") +
    #scale_fill_uchicago() +
    #scale_color_uchicago() +
    #theme_pubclean() +
    scale_fill_manual(values = mycolors) +
    scale_color_manual(values = mycolors) +
    mytheme +
    theme(panel.background = element_blank(),
          # axis.text.x=element_text(angle=40, vjust = -0.5, hjust=0.5)
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.text=element_text(face="italic")
    ) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0))
  ggsave(outname, g1, width=width, height = height)
  return(g1)
}
