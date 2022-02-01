#' @title FEMStress
#'
#' @description Creates a complete finite element model using stress for a given 2D mesh under specified boundary conditions (constrain and load).
#' @usage FEMStress(meshP, meshT, centroid, BoundConx, BoundCony, SFShear,
#' SFTensile, Length, area, Fx, Fy,  Y, Nu, Thick)
#'
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param meshT Matrix (3 x n) containing the number of the coordinate point that forms a given triangle within the mesh.
#' @param centroid Matrix (2 x n) containing coordinate points of the centroid of each triangular element.
#' @param BoundConx Boundary constraint for nodes in the x-direction
#' @param BoundCony Boundary constraint for nodes in the y-direction
#' @param SFTensile Magnitude of tensile surface traction; if there is no surface traction then SFTensile = 0
#' @param SFShear Magnitude of positive shear traction; if there is no surface traction then SFShear = 0
#' @param Length Truss length
#' @param Thick Triangle element thickness
#' @param area Triangle element area
#' @param Fx Load vector for the x-direction
#' @param Fy Load vector for the y-direction
#' @param Y Value of Young's (Elastic) modulus
#' @param Nu Value of Poisson's ratio
#' @param Thick Value of the thickness of the mesh, a value must be given.
#'
#' @return Completes the FEM to generate values of stress and strain; produces three (3) [3 x n] matrix.
#' \item{Strain}{Calculated strain. [x, y, tau]}
#' \item{Stress}{Calculated stress in pascals. [x, y, tau]}
#' \item{StressStrain}{Stress as calucated from strain. [x, y, tau]}
#'
#' @examples
#' Y = 20e5
#' Nu = 0.45
#' Thick = 0.001
#' DOF = 6
#'
#' # meshP = MeshPts$p
#' # meshT = MeshPts$T
#' # centroid = Centroids
#' # BoundConx = numeric(NROW(meshP))
#' # BoundCony = numeric(NROW(meshP))
#' # BoundConx[1:NROW(meshP)] = BoundCony[1:NROW(meshP)] = 1
#' # SFShear = 0 #can leave as 0
#' # SFTensile = 0 #can leave as 0
#' # Length = TrussLength
#' # area = Area
#' # Fx = 20
#' # Fy = 10
#' # FEAtest = FEAStress(meshT, meshP, centroid, BoundConx, BoundCony, SFShear, SFTensile, Length, area, Fx, Fy, Nu, Y, Thick)
#'
#' # Following the FEAtest, plot for value map:
#' # PlotVal = FEAtest$Stress
#' # Kol = 3
#' # a= 1; b= 2; c=3; d=4; e=5; f=6; g=7; h=8; i=9
#' # PlotSystem(meshP, meshT, PlotVal, Kol, a, b, c, d, e, f, g, h, i)
#'
#' @export

#Complete FEM
FEMStress = function(meshP, meshT, centroid, BoundConx, BoundCony, SFShear,
                             SFTensile, Length, area, Fx, Fy, Y, Nu, Thick){
  FEAt1 = ElementMat(meshP, meshT, Nu, Y, Thick)

  #Expanded element matrix represents the contribution of individual finite elements towards the global structural matrix
  ElementMatrixlist = FEAt1$EMPStress
  FEAt2 = ExpandEM(meshP, meshT, centroid, ElementMatrixlist)

  #Global stiffness matrix - once established expanded must be combined to create the global structural stiffness matrix (by adding the expanded matrices)
  ExEM = FEAt2 #expanded element matrix
  FEAt3 = GlobalMat(meshP, meshT, ExEM)

  #Boundary conditions for element centroids based on coordinate points. For the x & y direction per centroid create matrix with boundary 1(unfixed) or 0(fixed)
  FEAt4 = ApplyBC(BoundConx,BoundCony, meshP)

  #Reduced stiffness matrix - use boundary condition to reduce matrix to smaller form by removing systems that are bound.
  GMat = FEAt3
  NodeKnownL = FEAt4
  FEAt5 = ReducedEM(GMat, NodeKnownL)

  #Element Surface Traction - generates the column matric for the uniformly distributed, SFTensile and SFShear must be given values. load.
  FEAt6 = SurfaceTraction(meshP,SFTensile, SFShear, Length, Thick, area)

  #Expanded Surface Force
  SurfTrac = FEAt6 #surface traction
  FEAt7 = ExpandSFT(meshP, meshT, SurfTrac)

  #Reduced Surface force
  ExSurf = FEAt7
  FEAt8 = ReducedSF(meshP, ExSurf)

  #Reduced Load vector - values must be known or given in the problem statement and applied to the model.
  RSF = FEAt8 #reduced surface traction. If none is present RSF = 0
  NodeKnownL = FEAt4
  FEAt9 = ForceVector(Fx, Fy, RSF, meshP, NodeKnownL)

  #Global Nodal Displacement
  ForceV = FEAt9 #reduced force vector
  REM = FEAt5 #reduced element martix
  NodeKnownL = FEAt4
  FEAt10 = NodeDis(meshP, REM, ForceV, NodeKnownL)

  #Global and Local Forces
  GMat = FEAt3
  GlobalND = FEAt10$GlobalND
  ElementMatrixlist = FEAt1$EMPStress
  FEAt11 = GLForces(meshP, meshT, GMat, GlobalND, ElementMatrixlist)

  #Stresses
  GlobalND = FEAt10$GlobalND
  FEAt12 = LocalStress(meshP, meshT, Y, Nu, GlobalND)
}
