#' @title ExpandEM
#'
#' @description Generates the expanded element matrix, which represents the contribution of individual finite elements towards the global structural matrix
#'
#' @usage ExpandEM(meshP, meshT, centroid, EMatrixlist)
#'
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param meshT Matrix (3 x n) containing the number of the coordinate point that forms a given triangle within the mesh.
#' @param centroid Matrix (2 x n) containing coordinate points of the centroid of each triangular element.
#' @param EMatrixlist EMPStress or EMPStrain generated from ElementMat function. List of element matrices.
#'
#' @return Produces large (n x n) matrix.
#' \item{ExpandedMat}{The expanded element matrix}
#'
#' @examples
#' # EMatrixlist = EMPStress
#' # meshP = MeshPts$p #nodal points
#' # meshT = MeshPts$T #triangulation list
#' # centroid = Centroids
#' # ExpandEM(meshP, meshT, centroid, EMatrixlist)
#'
#' @export


ExpandEM = function(meshP, meshT, centroid, EMatrixlist){
  m= n= z= o= NROW(meshT)
  TDOF = 2*NROW(meshP)

  Expanded_EM = function(meshT, TDOF, centroid, EMatrixlist){
    Expanded = matrix(0, nrow=TDOF, ncol = TDOF)
    r1=2*meshT[m,1]-1
    r2=2*meshT[m,1]
    r3=2*meshT[m,2]-1
    r4=2*meshT[m,2]
    r5=2*meshT[m,3]-1
    r6=2*meshT[m,3]
    Expanded[c(r1,r2,r3,r4,r5,r6),c(r1,r2,r3,r4,r5,r6)] = EMatrixlist[[m]][c(1,2,3,4,5,6), c(1,2,3,4,5,6)]
    return(Expanded)}

  #Run for acquiring individual expanded element matrix=Expanded
  ExpandedMat = list()
  for (m in 1:z){
    ExpandedMat[[m]] = Expanded_EM(meshT, TDOF, centroid, EMatrixlist)}
  return(ExpandedMat)}
