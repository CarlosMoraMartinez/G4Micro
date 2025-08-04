plotRelativeAbnBarsGenus <- function(phobj, variable="Condition", topn = 15,
                                     outname="phylumBarplot.pdf", height=8,
                                     width=12, ocluster=T,
                                     oldlevs=c("Control", "Depression"),
                                     wespalette="Darjeeling2"){
  library(RColorBrewer)

  df_prevalence <- getRelAbundancesByGenusAndVariable(phobj, variable, outname="", oldlevs=oldlevs)
  df_top <- df_prevalence %>% top_n(topn, total_rel_prevalence) %>%
    filter(!is.na(Genus)) %>%
    dplyr::arrange(desc(total_rel_prevalence))
  taxdata <- tax_table(phobj) %>% data.frame() %>% filter(Genus %in% df_top$Genus)

  gendata <- taxdata %>% group_by(Genus) %>% dplyr::summarise(Phylum=unique(Phylum),
                                                              Family=unique(Family))

  otus_gen <- data.frame(otu_table(phobj) ) %>%
    rownames_to_column("ASV") %>%
    dplyr::mutate(is_top = ASV %in% rownames(taxdata),
                  Genus = taxdata$Genus[match(ASV, rownames(taxdata))]
    ) %>%
    dplyr::mutate(Genus = ifelse(is.na(Genus), "Other", Genus)) %>%
    group_by(Genus) %>%
    dplyr::summarise_if(is.numeric, sum) %>%
    filter(Genus != "Other") %>%
    dplyr::mutate(
      Phylum = gendata$Phylum[match(Genus, gendata$Genus)],
      Family = gendata$Family[match(Genus, gendata$Genus)],
    )

  otus_gen_gen_rel <- otus_gen %>%
    dplyr::mutate_if(is.numeric, function(x)x/sum(x))

  otus_df <- otus_gen_gen_rel %>%
    tidyr::gather(key="sampleID", value="Abundance", -Genus, -Phylum, -Family) %>%
    dplyr::mutate(sampleID = gsub("^X", "", sampleID))

  ranked_genera <- otus_gen_gen_rel %>%
    dplyr::mutate(mean_rel = rowMeans(across(where(is.numeric)))) %>%
    #top_n(topn, mean_rel) %>%
    dplyr::arrange(mean_rel) #desc(mean_rel)
  most_prev_genus <- ranked_genera$Genus[nrow(ranked_genera)]

  df_merged <- merge(otus_df, data.frame(sample_data(phobj)), by = "sampleID") %>%
    dplyr::mutate(Genus = factor(Genus, levels = ranked_genera$Genus))

  if(! ocluster){
    ranked_samples <- df_merged %>%
      filter(Genus == most_prev_genus) %>%
      dplyr::arrange(Abundance) %>%
      pull(sampleID)
  }else{
    dmat <- otus_gen_gen_rel %>% select_if(is.numeric) %>% t %>% dist
    hcl <- hclust(dmat)
    ranked_samples <- gsub("^X", "", hcl$labels[hcl$order])
  }

  df_merged <- df_merged %>% dplyr::mutate(sampleID = factor(sampleID, levels=ranked_samples))

  #cols <- brewer.pal(n=nrow(ranked_genera), name="Set2") # No more colors than those in palette
  #mycolors <- colorRampPalette(brewer.pal(8, "Spectral"))(nrow(ranked_genera))
  mycolors <- colorRampPalette(wes_palette(wespalette))(nrow(ranked_genera))
  facet_form <- paste0( ". ~ ", variable) %>% as.formula

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
