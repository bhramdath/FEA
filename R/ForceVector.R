#' @title ForceVector
#'
#' @description Creates a matrix of loads in the x & y direction for each load unconstrained node.
#'
#' @param Fx Load vector for the x-direction
#' @param Fy Load vector for the y-direction
#' @param RSF If surface traction is present assign value as the ReducedSF matrix; if there is no surface traction set RSF = 0
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param NodeKnownL data frame with constraint parameters applied to each node in the x and y directions. Formatted for use in reduced element matrix. Generated from ApplyBC function.
#'
#' @return Produces a matrix with loading parameters for each node.
#' \item{ReducedFV}{Reduced force vector matrix containing the model load parameters.}
#'
#' @examples
#' # meshP = MeshPts$p #mesh points
#' # RSF = RSF #reduced surface traction. If none is present, RSF = 0
#' # Fx = 10
#' # Fy = 10
#' # NodeKnownL = Nodelist #boundary conditions at nodes
#' # ForceVector(Fx, Fy, RSF, meshP, NodeKnownL)
#'
#' @export

ForceVector = function(Fx, Fy, RSF, meshP, NodeKnownL){
  z= NROW(meshP) # m=col, n=row, z=element#
  Vec = matrix(rbind(Fx, Fy), nrow = z, ncol =2)
  F.vector = matrix((rbind(Vec[,1], Vec[,2])), ncol = 1) + RSF
  ReducedFV = matrix(F.vector[c(NodeKnownL)], ncol = 1)}
