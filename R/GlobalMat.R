#' @title GlobalMat
#'
#' @description Generates global stiffness matrix - once established, the expanded element matrix must be combined to create the global structural stiffness matrix by adding the expanded matrices.
#'
#' @param meshT Matrix (3 x n) containing the number of the coordinate point that forms a given triangle within the mesh.
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param ExEM Expanded element matrix
#'
#' @return Produces large (n x n) global matrix
#' \item{GlobalMat}{Global matrix}
#'
#' @examples
#' # meshP = MeshPts$p
#' # meshT = MeshPts$T
#' # ExEM = ExEM #expanded element matrix
#' # GlobalMat(ExEM, meshT)
#'
#' @export

GlobalMat = function(meshP, meshT, ExEM){
  m= n= o= 2*NROW(meshP)
  z= NROW(meshT) # m=col, n=row, z=element#
  GlobalMat = matrix(0, nrow = o, ncol = o)
  for (m in 1:o){for (n in 1:o){
    cell = numeric(z)
    for (z in 1:z){
      cell[z] = ExEM[[z]][m, n]}
    GlobalMat[m,n] = sum(cell)}}
  return(GlobalMat)}
