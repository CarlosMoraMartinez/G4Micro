
#' @title Prepare Data for Sankey Plot Visualization
#' @description
#' Prepares node and link tables for creating a Sankey plot that shows relationships
#' between biological processes and species, including log-fold changes and expression values.
#' Supports filtering by adjusted p-values, aggregating non-significant species, and customizing
#' plot aesthetics such as node colors and sizes.
#'
#' @param daa_proc_all Data frame with process-level differential abundance analysis results.
#'   Must include columns with adjusted p-values (`padj`) and log2 fold changes (`log2FoldChange`).
#' @param daa_by_sp_all Data frame with process-species level differential abundance results,
#'   including columns `Pathway`, `Species`, `padj`, `log2FoldChange`, `AveExpr`, and taxonomic info.
#' @param quant_by_proc Data frame with abundance quantifications by process across samples.
#'   Rows are processes, columns are sample IDs.
#' @param metadf Data frame with sample metadata, including a `Condition` column and `sampleID` matching `quant_by_proc` columns.
#' @param plim Numeric threshold for adjusted p-value cutoff to select significant processes (default 0.01).
#' @param plim_sp Numeric threshold for adjusted p-value cutoff to select significant species within processes (default 0.01).
#' @param include_others Logical, whether to aggregate non-significant species into a single 'others' category (default FALSE).
#' @param include_longnames Logical, whether to use long descriptive names for processes if available (default FALSE).
#' @param others_name Character string naming the aggregated non-significant taxa group (default "Aggr. taxa").
#' @param case_name Character string for the case group name in the metadata (default "Depression").
#' @param control_name Character string for the control group name in the metadata (default "Control").
#' @param outdir Character string path where output TSV files for nodes and links will be saved (default "~/").
#' @param name Character string prefix for output file names (default "test").
#' @param ntop Integer indicating the maximum number of top significant processes to select by adjusted p-value;
#'   if 0, all processes passing `plim` are included (default 0).
#'
#' @return A list with two elements:
#' \item{nodelist}{Data frame describing nodes with columns for node ID, label, position, color, and size.}
#' \item{linklist}{Data frame describing links between nodes with source, target, size, and color for Sankey plotting.}
#'
#' @details
#' This function filters and processes differential abundance results at the process and species levels
#' based on specified significance thresholds. It normalizes abundance data within conditions and processes,
#' creates color gradients based on log-fold changes, and constructs node/link tables ready for Sankey plot
#' visualization tools. It optionally aggregates low-significance species and can annotate processes
#' with longer descriptive names.
#'
#' Output node and link tables are saved as TSV files in the specified output directory with
#' filenames based on the provided name prefix.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   result <- prepareSankeyPlot(
#'     daa_proc_all = process_da_df,
#'     daa_by_sp_all = species_da_df,
#'     quant_by_proc = quant_df,
#'     metadf = metadata_df,
#'     plim = 0.05,
#'     include_others = TRUE,
#'     outdir = "results/",
#'     name = "mySankey"
#'   )
#'   # Use result$nodelist and result$linklist for plotting
#' }
#' }
#'
#' @seealso
#' \code{\link[dplyr]{filter}}, \code{\link[dplyr]{mutate}}, \code{\link[dplyr]{summarise_if}},
#' \code{\link[dplyr]{select}}, \code{\link[dplyr]{arrange}}, \code{\link[scales]{gradient_n_pal}}
#'
#' @rdname prepareSankeyPlot
#' @export
#' @importFrom dplyr filter mutate summarise_if select arrange group_by ungroup slice_min
#' @importFrom tidyr pivot_longer
#' @importFrom readr write_tsv
#' @importFrom scales gradient_n_pal
prepareSankeyPlot <- function(daa_proc_all,
                              daa_by_sp_all,
                              quant_by_proc,
                              metadf,
                              plim = 0.01,
                              plim_sp = 0.01,
                              include_others = FALSE,
                              include_longnames = FALSE,
                              others_name = "Aggr. taxa",
                              case_name = "Depression",
                              control_name = "Control",
                              outdir = "~/",
                              name = "test",
                              ntop = 0){
  ## Returns a list with
  ##  1) a data frame of nodes and
  ##  2) a data frame of links between nodes.
  # Get processes to plot

  daa_proc1 <- daa_proc_all %>%
    dplyr::filter(padj< plim) %>%
    rownames_to_column("Pathway")
  if(nrow(daa_proc1) > ntop & ntop > 0){
    daa_proc1 <- daa_proc_all %>%
      slice_min(order_by = padj, n = ntop) %>%
      rownames_to_column("Pathway")
  }
  # Get process-species pairs that are present in selected processes,
  # with or without aggregating non-significant species
  if(include_others){
    daa_sp <- daa_by_sp_all %>%
      dplyr::filter(Pathway %in% daa_proc1$Pathway) %>%
      #dplyr::filter(padj < plim_sp) #%>%
      dplyr::mutate(Species=ifelse(padj < plim_sp, Species, others_name)) %>%
      group_by(Pathway, Genus, Species) %>%
      dplyr::summarise_if(is.numeric, mean) %>%
      ungroup()
  }else{
    daa_sp <- daa_by_sp_all %>%
      dplyr::filter(Pathway %in% daa_proc1$Pathway) %>%
      dplyr::filter(padj < plim_sp) #%>%
  }
  #all(daa_proc1$Pathway %in% daa_sp$Pathway)

  # Get color scale for species nodes (limits according to process:species LFC)
  scale_species = scales::gradient_n_pal(colours=c(C_CTRL, C_WHITE, C_CASE),
                                         values=c(min(c(-1,min(daa_sp$log2FoldChange))),
                                                  0,
                                                  max(c(1,max(daa_sp$log2FoldChange)))))
  # Get color scale for process nodes (limits according to process LFC)
  scale_proc = scales::gradient_n_pal(colours=c(C_CTRL, C_WHITE, C_CASE),
                                      values=c(min(c(-1,min(daa_proc1$log2FoldChange))),
                                               0,
                                               max(c(1,max(daa_proc1$log2FoldChange)))))
  # Get color scale for process-species link (limits according to process:species LFC)
  scale_species_link = scales::gradient_n_pal(colours=c(C_CTRL_LINK, C_WHITE, C_CASE_LINK),
                                              values=c(min(c(-1,min(daa_sp$log2FoldChange))),
                                                       0,
                                                       max(c(1,max(daa_sp$log2FoldChange)))))

  part1 <- daa_sp %>% dplyr::select(Pathway, Species, AveExpr, log2FoldChange, padj) %>%
    #mutate(col = ifelse(padj > plim, C_NS, ifelse(log2FoldChange < 0, C_CTRL, C_CASE))) %>%
    dplyr::mutate(col=scale_species_link(log2FoldChange)) %>%
    dplyr::select(-padj) %>%
    group_by(Pathway) %>%
    dplyr::mutate(AveExpr = AveExpr/(sum(AveExpr)))
  names(part1) <- c("node1", "node2", "Size", "LFC", "color")

  tab_quant_all <- quant_by_proc %>% dplyr::filter(Pathway %in% daa_proc1$Pathway)

  part2 <- tab_quant_all %>%
    pivot_longer(2:ncol(tab_quant_all)) %>%
    dplyr::mutate(Condition = metadf$Condition[match(name, metadf$sampleID)]) %>%
    group_by(Condition, Pathway) %>%
    dplyr::summarise(value=mean(value)) %>%
    dplyr::mutate(LFC= 0,
                  color = ifelse(Condition == case_name, C_CASE_LINK, C_CTRL_LINK)) %>%
    ungroup() %>%
    group_by(Condition) %>%
    dplyr::mutate(value = value/sum(value)) %>%
    ungroup() %>% group_by(Pathway) %>% dplyr::mutate(value = value/sum(value))


  part2 <- part2 %>% select(Condition, Pathway, value, LFC, color) %>%
    dplyr::arrange(Condition)
  names(part2) <-  c("node1", "node2", "Size", "LFC", "color")

  nodes2fill <- daa_proc1$Pathway[!daa_proc1$Pathway %in% part1$node1]
  if(length(nodes2fill)>0){
    part3 <- data.frame(node1 = nodes2fill,
                        node2 = others_name,
                        Size=OTHERS_SIZE,
                        LFC=daa_proc1$log2FoldChange[match(nodes2fill, daa_proc1$Pathway)]) %>%
      dplyr::mutate(color = scale_species_link(LFC))

    full_v1 <- rbind(part2, part1, part3)
  }else{
    full_v1 <- rbind(part2, part1)
  }


  nodelist_1 <- data.frame(value = c(control_name, case_name),
                           value2match = c(control_name, case_name),
                           xpos = c(0,0),
                           color = c(C_CTRL, C_CASE),
                           nodesize = c(prop.table(table(metad$Condition))[control_name],
                                        prop.table(table(metad$Condition))[case_name]))

  tb_norm <- tab_quant_all %>% mutate_if(is.numeric, \(x)x/sum(x)) %>%
    dplyr::select(-Pathway) %>% rowSums()
  tb_norm <- tb_norm/sum(tb_norm)

  include_longnames = (include_longnames & "process_name" %in% names(daa_proc1) & "process_name" %in% names(daa_by_sp_all))
  nodelist_2 <- rbind(nodelist_1,
                      data.frame(
                        value = {if(include_longnames)paste(daa_proc1$Pathway, daa_proc1$process_name, sep=":") else daa_proc1$Pathway},
                        value2match = daa_proc1$Pathway,
                        xpos = 0.15,
                        color = scale_proc(daa_proc1$log2FoldChange),
                        #color = ifelse(daa_proc1$log2FoldChange < 0, C_CASE, C_CTRL),
                        nodesize = tb_norm[match(daa_proc1$Pathway, tab_quant_all$Pathway)]
                      ))
  species <- part1 %>% group_by(node2) %>%
    dplyr::summarise(LFC= mean(LFC),
                     SumExpr = sum(Size)) #%>%
  # dplyr::mutate(color = ifelse(LFC < 0, C_CASE, C_CTRL))
  scale_species_node = scales::gradient_n_pal(colours=c(C_CTRL, C_WHITE, C_CASE),
                                              values=c(min(c(-0.1, min(species$LFC))),
                                                       0,
                                                       max(c(0.1, max(species$LFC)))))
  nodelist_3 <- data.frame(
    value = species$node2,
    value2match = species$node2,
    xpos = 1,
    #color = species$color,
    color = scale_species_node(species$LFC),
    nodesize = species$SumExpr
  )

  nodelist <- rbind(nodelist_2, nodelist_3) %>% dplyr::mutate(nodenum=0:(n()-1))
  if(others_name %in% full_v1$node2 & ! others_name %in% nodelist$value){
    nodelist <- rbind(nodelist,
                      data.frame(value = others_name,
                                 value2match = others_name,
                                 xpos = 1,
                                 color=C_NS,
                                 nodesize=1, nodenum=nrow(nodelist))
    )
  }
  nodelist$color[nodelist$value==others_name] <- C_NS

  # nodelist <- full_v1 %>% select(node1, node2, color) %>% pivot_longer(node1:node2) %>%
  #   select(-name) %>% dplyr::distinct() %>%
  #   mutate(nodenum = 0:(n()-1))
  full_v1 <- full_v1 %>% dplyr::mutate(source=nodelist$nodenum[match(node1, nodelist$value2match)],
                                       target = nodelist$nodenum[match(node2, nodelist$value2match)])
  write_tsv(nodelist, file=paste0(outdir, "/", name, "_nodelist4sankey.tsv"))
  write_tsv(full_v1, file=paste0(outdir, "/", name, "_linklist4sankey.tsv"))
  return(list(nodelist=nodelist, linklist=full_v1))
}
