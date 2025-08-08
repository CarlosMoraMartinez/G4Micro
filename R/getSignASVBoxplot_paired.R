#' @title Generate Paired Boxplots for Significant ASVs
#' @description Creates violin and boxplots for ASVs significantly associated with conditions in paired samples.
#' @param resdf_annot Data frame containing ASV differential abundance results with taxonomic annotation and p-values.
#' @param raw_counts Numeric matrix or data frame of raw counts for ASVs (rows = ASVs, columns = samples).
#' @param outdir Output directory path where the plot PDF will be saved. Default is current directory ("./").
#' @param name Output filename for the saved plot PDF. Default is "sign_asvs_boxplot.pdf".
#' @param mode Character specifying data transformation mode: "PROP" to convert counts to proportions by sample, otherwise raw counts used. Default is "PROP".
#' @param max_n Maximum number of significant ASVs (p-value < 0.05) to plot. Defaults to 10.
#' @return A ggplot2 object representing the generated plot.
#' @details The function filters ASVs by p-value significance, optionally converts counts to proportions, reshapes the data, and plots paired violin and boxplots for each ASV faceted by taxon. The paired nature is shown by lines connecting samples from the same patient across conditions.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   plot <- getSignASVBoxplot_paired(resdf_annot, raw_counts, outdir = "./plots", name = "asvs_plot.pdf")
#'   print(plot)
#' }
#' }
#' @seealso
#'  \code{\link[dplyr]{arrange}}, \code{\link[dplyr]{mutate}}
#'  \code{\link[tidyr]{gather}}
#' @rdname getSignASVBoxplot_paired
#' @export
#' @importFrom dplyr arrange mutate
#' @importFrom tidyr gather
getSignASVBoxplot_paired<-function(resdf_annot, raw_counts, outdir="./", name="sign_asvs_boxplot.pdf", mode="PROP", max_n=10){
  ASVs <- resdf_annot %>% dplyr::arrange(pvalue) %>% filter(pvalue < 0.05)  %>% unite("taxon2", taxon, Genus, sep=":", remove=F)
  if(nrow(ASVs) > max_n){ASVs <- ASVs[1:10,]}

  if(mode=="PROP"){
    prop_counts <- apply(raw_counts, MAR=2,function(x)x/sum(x))
  }else{
    prop_counts <- raw_counts
  }
  df2plot <- prop_counts[ASVs$taxon, ] %>%
    as.data.frame %>% rownames_to_column("ASV") %>%
    dplyr::mutate(Genus = ASVs$Genus,taxon=ASVs$taxon2) %>%
    tidyr::gather("sampleID", "raw_counts", -ASV, -Genus, -taxon ) %>%
    dplyr::mutate(patient = metadata$pacienteID[match(sampleID, metadata$sampleID)],
                  Condition = metadata$Condition[match(sampleID, metadata$sampleID)])
  df2plot$taxon=gsub(":$", "", df2plot$taxon)

  g1 <- ggplot(df2plot, aes(x = Condition, y=raw_counts, col=Condition, fill=Condition)) +
    facet_wrap(. ~ taxon , scales = "free")+
    geom_violin(alpha=0.1)+
    geom_boxplot(width=0.1, fill="darkgray", notchwidth = 0.5, notch=F)+
    geom_point(size=1)+
    geom_line(aes(group = patient), color = "gray20", linetype=1, size=0.1, alpha=0.25) +
    # scale_color_lancet()+
    # scale_fill_lancet() +
    theme_pubclean() +
    ylab("Proportion")

  oname <- paste0(outdir, "/", name)
  ggsave(filename=oname, g1, width=6, height = 8)
  return(g1)

}
