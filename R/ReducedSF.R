#' @title ReducedSF
#'
#' @description Reduced matrix of surface forces
#'
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param ExSurf Expanded surface matrix, output from ExpandSFT
#'
#' @return Produces a large matrix.
#' \item{RSF}{Produces a large, reduced surface force matrix}
#'
#' @examples
#' # meshP = MeshPts$p
#' # ExSurf = ExSurf #expanded surface matrix
#' # ReducedSF(meshP, ExSurf)
#'
#' @export

#Reduced Surface force
ReducedSF = function(meshP, ExSurf){
  TDOF = 2*NROW(meshP)
  z= NROW(meshP)

  RSF = matrix(0, nrow = TDOF)
  for (n in 1:TDOF){
    cell = numeric(TDOF)
    for (m in 1:z){
      cell[m] = ExSurf[[m]][n,]}
    RSF[n,] = sum(cell) }
  return(RSF)}
