#' @title ApplyBC
#'
#' @description Boundary constraint for element centroids based on coordinate points. For the x & y direction per centroid create matrix with boundary 1(unfixed) or 0(fixed).
#' @usage ApplyBC(meshP, BoundConx, BoundCony)
#'
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh
#' @param BoundConx Boundary constraint for nodes in the x-direction
#' @param BoundCony Boundary constraint for nodes in the y-direction
#'
#' @return A data frame with constraint parameters applied to each node in the x and y directions. Formatted for use in reduced element matrix.
#' \item{NodeKnownL}{Constraint parameters}
#'
#' @examples
#' # meshP = MeshPts$p #mesh points
#' # BoundConx = numeric(NROW(meshP))
#' # BoundCony = numeric(NROW(meshP))
#' # BoundConx[1:NROW(meshP)] = BoundCony[1:NROW(meshP)] = 1
#' # ApplyBC(meshP, BoundConx, BoundCony)
#'
#' @export

ApplyBC = function(meshP, BoundConx, BoundCony){
  n = 2*NROW(meshP) # m=col, n=row, z=element#
  BoundConxy = matrix(rbind(BoundConx, BoundCony), byrow = TRUE, ncol = 1)

  Apply_BC = function(BoundConxy, n){
    for (m in 1:n){if (BoundConxy[m] < 0){BoundConxy[m] = 0} else {next}}
    NodeKnownL = (1:n)*BoundConxy} #ONLY if loads are already listed as in the example

  NodeKnownL = Apply_BC(BoundConxy, n)}
