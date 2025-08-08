#' @title Generate MA and Volcano Plots for Multiple DESeq2 Contrasts
#' @description
#' Creates MA plots and volcano plots for a set of DESeq2 differential expression results
#' under different log-fold change shrinkage methods (raw, normal, apeglm, ashr).
#' For each contrast, plots are saved as PDF files in the corresponding results directory.
#'
#' @param all_contrasts
#' A list of contrast result objects.
#' Each element must contain:
#' - `res`: DESeq2 results table (no shrinkage)
#' - `resLFC`: results table with LFC shrinkage using `"normal"` (optional)
#' - `resLFC_ape`: results table with LFC shrinkage using `"apeglm"` (optional)
#' - `resLFC_ashr`: results table with LFC shrinkage using `"ashr"` (optional)
#' - `nested_dir`: output directory path for saving plots
#' - `name`: name of the contrast (used in file names)
#' @param opt
#' A list of plotting and analysis options to be passed to `getSummaryTablesDeseq()`
#' and plotting functions (`make_maplot()`, `make_volcano()`).
#'
#' @return
#' Invisibly returns `NULL`.
#' Side effects: generates and saves MA and volcano plots as PDF files in the
#' specified directories.
#'
#' @details
#' This function iterates over multiple DESeq2 contrast results and produces:
#' - **MA plots** for raw and shrunken LFC estimates
#' - **Volcano plots** for p-values and adjusted p-values
#' Each plot is generated for available shrinkage types (`normal`, `apeglm`, `ashr`).
#' Errors in individual plot generation are caught and reported without stopping the loop.
#' All open graphical devices are closed at the end of each iteration.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname make_all_maplots
#' @export
make_all_maplots <- function(all_contrasts, opt){
  for(singleres in all_contrasts){

    rtabs <- getSummaryTablesDeseq(singleres$res, opt, singleres$nested_dir)

    tryCatch(make_maplot(singleres$res, opt, paste0(singleres$nested_dir,singleres$name, "_MAPlot-rawFC.pdf")),
             error=\(x)cat("Error make_maplot rawFC\n"))

    tryCatch(make_volcano(singleres$res, opt, paste0(singleres$nested_dir,singleres$name, "volcano_rawfc_rawpval.pdf"), "pvalue"),
             error=\(x)cat("Error make_volcano raw FC, raw p-val\n"))

    tryCatch(make_volcano(singleres$res, opt, paste0(singleres$nested_dir,singleres$name, "volcano_rawfc_adjpval.pdf"), "padj"),
             error=\(x)cat("Error make_volcano  raw FC, adj p-val\n"))

    if(! is_empty(singleres$resLFC)){
      tryCatch(make_maplot(singleres$resLFC, opt,  paste0(singleres$nested_dir, singleres$name, "_MAPlot-rawFC-normal.pdf")),
               error=\(x)cat("Error make_maplot rawFC-normal\n"))

      tryCatch(make_volcano(singleres$resLFC, opt, paste0(singleres$nested_dir,singleres$name, "volcano_shnormfc_rawpval.pdf"), "pvalue"),
               error=\(x)cat("Error make_volcano Shrink normal, raw p-val\n"))

      tryCatch(make_volcano(singleres$resLFC, opt, paste0(singleres$nested_dir,singleres$name, "volcano_shnormfc_adjpval.pdf"), "padj"),
               error=\(x)cat("Error make_volcanoShrink normal, adj p-val \n"))
      }

    if(! is_empty(singleres$resLFC_ape)){

      tryCatch(make_maplot(singleres$resLFC_ape, opt,  paste0(singleres$nested_dir,singleres$name, "_MAPlot-rawFC-ape.pdf")),
               error=\(x)cat("Error make_maplot rawFC-apet\n"))

      tryCatch(make_volcano(singleres$resLFC_ape, opt, paste0(singleres$nested_dir,singleres$name, "volcano_shapefc_rawpval.pdf"), "pvalue"),
               error=\(x)cat("Error make_volcano Shrink ape, raw p-val\n"))

      tryCatch(make_volcano(singleres$resLFC_ape, opt, paste0(singleres$nested_dir,singleres$name, "volcano_shapefc_adjpval.pdf"), "padj"),
               error=\(x)cat("Error make_volcano Shrink ape, adj p-val\n"))
    }

    if(! is_empty(singleres$resLFC_ashr)){

      tryCatch(make_maplot(singleres$resLFC_ashr, opt,  paste0(singleres$nested_dir,singleres$name, "_MAPlot-rawFC-ashr.pdf")),
               error=\(x)cat("Error make_maplot rawFC-ashr\n"))

      tryCatch(make_volcano(singleres$resLFC_ashr, opt, paste0(singleres$nested_dir,singleres$name, "volcano_shashfc_rawpval.pdf"), "pvalue"),
               error=\(x)cat("Error make_volcano Shrink ashr, raw p-val\n"))

      tryCatch(make_volcano(singleres$resLFC_ashr, opt, paste0(singleres$nested_dir,singleres$name, "volcano_shashfc_adjpval.pdf"), "padj"),
               error=\(x)cat("Error make_volcano Shrink ashr, adj p-val\n"))
    }

    while(dev.cur() != 1) dev.off()
  }
}
