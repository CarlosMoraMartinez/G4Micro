#' @title Generate a Placeholder Plot with an Error Message
#' @description Creates an empty ggplot object displaying an error message as text, useful when a plot cannot be generated.
#' @param opt Parameter reserved for future use or to match generic function signature (currently unused).
#' @param name Character string specifying the context or name related to the plot. Default is "var".
#' @param error Logical flag indicating the type of error message to display. If FALSE, the message includes the `name`; if TRUE, a generic error message is shown. Default is FALSE.
#' @return A ggplot2 plot object showing an error message as text.
#' @details This function helps to return a consistent placeholder plot when an error occurs during plotting, so that downstream code that expects a plot can still run without crashing.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Generate a placeholder plot for a failed plot named "variance"
#'   print(getNullPlot(NULL, name = "variance"))
#' }
#' }
#' @rdname getNullPlot
#' @export
getNullPlot <- function(opt, name="var", error=FALSE){
  if(!error){
    emptylab <- paste("Error plotting ", name, sep="", collapse="")
  }else{
    emptylab = "Error while plotting"
  }

  gempty <- ggplot() +
    geom_text(aes(x=1, y=1),label=emptylab, size=14)+
    theme_light() +
    theme(axis.text = element_blank(), axis.title = element_blank())
  return(gempty)
}
