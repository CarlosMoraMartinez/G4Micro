#' @title Plot Relative Abundance Barplots by Phylum
#' @description
#' This function creates a stacked barplot of relative abundances at the phylum level across samples,
#' grouped by a categorical variable. Samples can be ordered by the abundance of the most prevalent phylum
#' or clustered hierarchically.
#'
#' @param phobj A \code{phyloseq} object containing taxonomic and sample data.
#' @param variable A string indicating the sample variable to facet the barplots by. Default: \code{"Condition"}
#' @param outname Output filename for the PDF plot. Default: \code{"phylumBarplot.pdf"}
#' @param height Height (in inches) of the output PDF plot. Default: \code{8}
#' @param width Width (in inches) of the output PDF plot. Default: \code{12}
#' @param ocluster Logical indicating whether to cluster the samples instead of ordering by abundance. Default: \code{FALSE}
#'
#' @return A \code{ggplot} object with the generated barplot, and a PDF file is saved to disk.
#'
#' @details
#' This function transforms the data to relative abundances using \code{transform_sample_counts}, then
#' creates a stacked barplot of phylum-level taxa using \code{plot_bar}. Samples can be sorted by
#' the abundance of the most prevalent phylum or by hierarchical clustering.
#'
#' Phyla are ordered from least to most prevalent globally to ensure consistent bar color mapping.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[phyloseq]{transform_sample_counts}},
#'  \code{\link[phyloseq]{plot_bar}}
#'  \code{\link[dplyr]{arrange}},
#'  \code{\link[dplyr]{summarise}}
#' @rdname plotRelativeAbnBarsPhylum
#' @export
#' @importFrom phyloseq transform_sample_counts plot_bar
#' @importFrom dplyr arrange summarise
#' @importFrom RColorBrewer brewer.pal
plotRelativeAbnBarsPhylum <- function(phobj, variable="Condition",
                                      outname="phylumBarplot.pdf",
                                      height=8, width=12, ocluster=F){
  ps_rel_abund = phyloseq::transform_sample_counts(phobj, function(x){x / sum(x)})
  facet_form <- as.formula(paste0(". ~ ", variable))


  g1 <- phyloseq::plot_bar(ps_rel_abund, fill = "Phylum") +
    geom_bar(aes(color = Phylum,
                 fill = Phylum),
             stat = "identity", position = "stack") +
    labs(x = "", y = "Relative Abundance\n") +
    facet_wrap(facet_form, scales = "free") +

    #scale_fill_uchicago() +
    #scale_color_uchicago() +
    #theme_pubclean() +
    mytheme +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme(panel.background = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

  df_prevalence <- getMeanRelAbundancesByPhylum(phobj)  %>%
    dplyr::arrange(prevalence)
  g1$data$Phylum <- factor(g1$data$Phylum, levels = df_prevalence$Phylum) #Sort phylums from most to least prevalent

  #Sort samples according to prevalence of the most prevalent taxon or by clustering
  if(! ocluster){
    most_prevalent <- df_prevalence$Phylum[nrow(df_prevalence)]
    sorted_samples <- g1$data %>% filter(Phylum == most_prevalent) %>%
      group_by(sampleID) %>% dplyr::summarise(Abundance = sum(Abundance)) %>%
      dplyr::arrange(Abundance) %>% pull(sampleID)
    g1$data$Sample <- factor(g1$data$Sample, levels = sorted_samples)
  }else{
    dmat <- ps_rel_abund %>% otu_table() %>% t %>% dist
    hcl <- hclust(dmat)
    sorted_samples <- gsub("^X", "", hcl$labels[hcl$order])
    g1$data$Sample <- factor(g1$data$Sample, levels = sorted_samples)
  }

  mycolors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Spectral"))(nrow(df_prevalence)) # %>% rev
  #mycolors[is.na(df_prevalence$Phylum) | df_prevalence$Phylum == "NA"] <- "#000000"

  g1 <- g1 +
    scale_fill_manual(values = mycolors) +
    scale_color_manual(values = mycolors)
  ggsave(outname, g1, width=width, height = height)
  return(g1)
}
