#' @title LocalStress
#'
#' @description Calculates local stress and strain for triangular elements of the mesh
#'
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param meshT Matrix (3 x n) containing the number of the coordinate point that forms a given triangle within the mesh.
#' @param Nu Value of Poisson's ratio
#' @param Y Value of Young's (Elastic) modulus
#' @param GlobalND Global nodal displacement, return from function NodeDis
#'
#' @return Completes FEM by calculating values of stress and strain, produces three (3) [3 x n] matrix.
#' \item{Strain}{Calculated strain. [x, y, tau]}
#' \item{Stress}{Calculated stress in pascals. [x, y, tau]}
#' \item{StressStrain}{Stress as calucated from strain. [x, y, tau]}
#'
#' @examples
#' # Y = 17e9
#' # Nu = 0.45
#' # meshP = MeshPts$p
#' # meshT = MeshPts$T
#' # GlobalND = GlobalND
#' # LocalStress(meshP, meshT, Y, Nu, GlobalND)
#'
#' @export

#Stresses
LocalStress = function(meshP, meshT, Y, Nu, GlobalND){
  m= n= o= NROW(meshP)*2
  z= NROW(meshT) # m=col, n=row, z=element#

  LStress = function(meshP, meshT, Y, Nu, GlobalND){
    x1 = meshP[meshT[m,1],1]
    y1 = meshP[meshT[m,1],2]
    x2 = meshP[meshT[m,2],1]
    y2 = meshP[meshT[m,2],2]
    x3 = meshP[meshT[m,3],1]
    y3 = meshP[meshT[m,3],2]

    A2 = x3*(y1-y2) + x1*(y2-y3) + x2*(-y1+y3)
    A = A2/2

    #Strain (E) = [B] *Ue
    B = matrix(c(((y2-y3)/A2), 0, ((y3-y1)/A2), 0, ((y1-y2)/A2), 0,
                 0, ((x3-x2)/A2), 0, ((x1-x3)/A2), 0, ((x2-x1)/A2),
                 ((x3-x2)/A2), ((y2-y3)/A2), ((x1-x3)/A2),((y3-y1)/A2), ((x2-x1)/A2),((y1-y2)/A2)), nrow = 3, byrow = TRUE)

    #Material properties of the element [D] differs based on the nature of the problem
    #for plane stress
    d1 = Y/(1-Nu^2)
    d2 = matrix(c(1, Nu, 0,
                  Nu, 1, 0,
                  0, 0, ((1-Nu)/2)), nrow = 3, byrow = TRUE)
    D1 = d1*d2

    #for plane strain
    d3 = (Y*(1-Nu))/((1+Nu)*(1-2*Nu))
    d4 = matrix(c(1, Nu/(1-Nu), 0,
                  Nu/(1-Nu), 1, 0,
                  0, 0, (1-2*Nu)/(2*(1-Nu))), nrow = 3, byrow = TRUE)
    D2 = d3 * d4

    r1=2*meshT[m,1]-1
    r2=2*meshT[m,1]
    r3=2*meshT[m,2]-1
    r4=2*meshT[m,2]
    r5=2*meshT[m,3]-1
    r6=2*meshT[m,3]

    #Local Strain
    Strain = B %*% GlobalND[c(r1, r2, r3, r4, r5, r6)]

    #Local stress with plane stress
    LocalStress1 = D1 %*% B %*% GlobalND[c(r1, r2, r3, r4, r5, r6)]

    #Local stress with plane strain
    LocalStress2 = D2 %*% B %*% GlobalND[c(r1, r2, r3, r4, r5, r6)]

    Out = matrix(c(Strain, LocalStress1, LocalStress2), ncol = 3, nrow = 3, byrow = TRUE)
    return(Out)}

  OutMat = list() #run for local stress calculation for each element. (Computed stresses are σ_xx, σ_yy, and Tau_xy)
  for (m in 1:z){
    OutMat[[m]] = LStress(meshP, meshT, Y, Nu, GlobalND)}

  Strain = matrix(0, ncol=3, nrow = z, byrow = TRUE)
  Stress1 = matrix(0, ncol=3, nrow = z, byrow = TRUE)
  Stress2 = matrix(0, ncol=3, nrow = z, byrow = TRUE)

  for (m in 1:z){
    Strain[m,] = OutMat[[m]][1,]
    Stress1[m,] = OutMat[[m]][2,]
    Stress2[m,] = OutMat[[m]][3,]}
  Rlist = list("Strain" = Strain, "Stress" = Stress1, "StressFromStrain" = Stress2)}
