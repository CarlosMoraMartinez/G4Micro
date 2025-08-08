
#' @title Generate a Venn Diagram Plot
#' @description Creates a Venn diagram visualization using the `ggvenn` package, saves it as a PDF file, and returns the plot object.
#' @param vars2venn A named list of vectors, each containing elements to be plotted as a set in the Venn diagram.
#' @param name A character string specifying the base name for the output PDF file. Default is `"VennDiagram"`.
#' @param opt A list containing output options. Must include an `out` element specifying the directory path where the PDF will be saved.
#' @param w Numeric value specifying the width of the output PDF in inches. Default is 5.
#' @param h Numeric value specifying the height of the output PDF in inches. Default is 5.
#' @return A ggplot object representing the generated Venn diagram.
#' @details The function uses fixed colors for the diagram sections and disables the display of set elements inside the plot. Adjust `fill_color` and other graphical parameters as needed.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   set1 <- c("A", "B", "C")
#'   set2 <- c("B", "C", "D")
#'   sets <- list(Group1 = set1, Group2 = set2)
#'   outdir <- list(out = "./plots/")
#'   dir.create(outdir$out, showWarnings = FALSE)
#'   venn_plot <- makeVenn(vars2venn = sets, name = "ExampleVenn", opt = outdir)
#'   print(venn_plot)
#' }
#' }
#' @rdname makeVenn
#' @export
#' @importFrom ggvenn ggvenn
makeVenn <- function(vars2venn, name="VennDiagram", opt, w=5, h=5){
  gv <- ggvenn(
    vars2venn, columns = names(vars2venn),
    stroke_size = 0.5,
    stroke_color = C_NS,
    fill_color = c(C_CASE, C_CASE2, C_CTRL2, C_CTRL),show_elements = F
  )
  ggsave(filename = paste0(opt$out, name, ".pdf"), gv, width = w, height = h)
  return(gv)
}
