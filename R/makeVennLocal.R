
#' @title Create and Save a Venn Diagram Plot
#' @description Generates a Venn diagram plot from a named list of sets and saves it as a PDF file in the specified output directory.
#' @param vars2venn A named list of vectors or sets to visualize in the Venn diagram. Each element represents a set.
#' @param name Filename prefix for the saved plot PDF. Default: \code{"VennDiagram"}.
#' @param outdir Directory path where the PDF will be saved.
#' @param w Width of the saved PDF in inches. Default: \code{5}.
#' @param h Height of the saved PDF in inches. Default: \code{5}.
#' @return A \code{ggplot} object representing the Venn diagram.
#' @details
#' The function uses \code{ggvenn} to create a Venn diagram from the input sets, customizing stroke and fill colors. The plot is saved as a PDF file named with the specified prefix in the given directory.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   sets <- list(GroupA = c("a", "b", "c"), GroupB = c("b", "c", "d"))
#'   venn_plot <- makeVennLocal(sets, name="ExampleVenn", outdir="./plots/", w=6, h=6)
#'   print(venn_plot)
#' }
#' }
#' @rdname makeVennLocal
#' @export
#' @importFrom ggvenn ggvenn
makeVennLocal <- function(vars2venn, name="VennDiagram", outdir, w=5, h=5){

  gv <- ggvenn(
    vars2venn, columns = names(vars2venn),
    stroke_size = 0.5,
    stroke_color = C_NS,
    fill_color = c(C_CASE, C_CASE2, C_CTRL2, C_CTRL),show_elements = F
  )
  ggsave(filename = paste0(outdir, name, ".pdf"), gv, width = w, height = h)
  return(gv)
}
