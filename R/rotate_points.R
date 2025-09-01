#' @title Rotate a Set of 2D Points
#' @description Rotates a set of 2D points so that a specified element lies on the positive x-axis.
#'
#' @param mat A numeric matrix with two columns representing x and y coordinates of points.
#'            Each row corresponds to a point.
#' @param vnames A character vector of point names. Must be the same length as the number of rows in \code{mat}.
#' @param element_to_right A character string indicating which point (from \code{vnames}) should be rotated to lie on the positive x-axis.
#'
#' @return A numeric matrix of the same dimensions as \code{mat}, containing the rotated coordinates of the input points.
#'         Row names correspond to \code{vnames}.
#'
#' @details This function computes the angle required to rotate the specified \code{element_to_right} point onto the positive x-axis,
#' and applies this rotation to all points in \code{mat}. The rotation is performed counterclockwise using a standard 2D rotation matrix.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Example: rotate three points so that "B" lies on the positive x-axis
#'   mat <- matrix(c(1,1, -1,2, 0,-1), ncol = 2, byrow = TRUE)
#'   vnames <- c("A", "B", "C")
#'   rotated <- rotate_points(mat, vnames, "B")
#'   print(rotated)
#' }
#' }
#'
#' @seealso
#'  \code{\link[assertthat]{assert_that}}
#' @rdname rotate_points
#' @export
#' @importFrom assertthat assert_that
rotate_points <- function(mat, vnames, element_to_right){
  assertthat::assert_that(element_to_right %in% vnames, msg = "rotate_points: element_to_right to calculate angle not in vnames")
  rownames(mat) <- vnames

  rotation <- pi - atan2(mat[element_to_right, 2], mat[element_to_right, 1])

  R <- matrix(c(cos(rotation), -sin(rotation),
                sin(rotation),  cos(rotation)),
              nrow = 2)
  rotated <- mat %*% R
  return(rotated)
}
