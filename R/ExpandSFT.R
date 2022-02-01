#' @title ExpandSFT
#'
#' @description Generates expanded surface force element matrix from SurfaceTraction function
#'
#' @param meshT Matrix (3 x n) containing the number of the coordinate point that forms a given triangle within the mesh.
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param SurfTrac List of surface forces.
#'
#' @return Produces a large (n x n) element matrix of surface forces.
#' \item{ExpandedSurf}{Expanded surface force element matrix.}
#'
#' @examples
#' # meshT = MeshPts$T
#' # meshP = MeshPts$p
#' # SurfTrac = SurfaceTraction #matrix
#' # ExpandSFT(meshP, meshT, SurfTrac)
#'
#' @export

#Expand Surface Force
ExpandSFT = function(meshP, meshT, SurfTrac){
  m= n= z= o= NROW(meshP)
  TDOF = 2*z

  ESFT = function(meshT, TDOF, SurfTrac){
    Expanded = matrix(0, nrow=TDOF)
    r1=2*meshT[m,1]-1
    r2=2*meshT[m,1]
    r3=2*meshT[m,2]-1
    r4=2*meshT[m,2]
    r5=2*meshT[m,3]-1
    r6=2*meshT[m,3]
    Expanded[c(r1,r2,r3,r4,r5,r6)] = SurfTrac[[m]]
    return(Expanded)}

  ExpandedSurf = list()
  for (m in 1:z){
    ExpandedSurf[[m]] = ESFT(meshT, TDOF, SurfTrac) }
  return(ExpandedSurf)}
