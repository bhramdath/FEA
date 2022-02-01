#' @title ElementMat
#'
#' @description Generates an element stiffness matrix
#'
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param meshT Matrix (3 x n) containing the number of the coordinate point that forms a given triangle within the mesh.
#' @param Nu Value of Poisson's ratio
#' @param Y Value of Young's (Elastic) modulus
#' @param Thick Value of the thickness of the mesh, a positive value must be given.
#'
#' @return Generates initial element matrix needed for the finite element model.
#' \item{EMPStress}{An element matrix of the geometry under stress.}
#' \item{EMPStrain}{An element matrix of the geometry under strain.}
#'
#' @examples
#' # Y = 20 #Elastic modulus
#' # Nu = 0.45 #Poisson ratio
#' # Thick = 0.001
#' # DOF = 6
#' # meshP = MeshPts$p
#' # meshT = MeshPts$T
#' # test7 = ElementMat(meshP, meshT, Nu, Y, Thick)
#'
#' @export

ElementMat = function(meshP, meshT, Nu, Y, Thick){
  CST_EM = function(meshP, meshT, Nu, Y, Thick){
    x1 = meshP[meshT[m,1],1]
    y1 = meshP[meshT[m,1],2]
    x2 = meshP[meshT[m,2],1]
    y2 = meshP[meshT[m,2],2]
    x3 = meshP[meshT[m,3],1]
    y3 = meshP[meshT[m,3],2]

    A2 = x3*(y1-y2) + x2*(y3-y1) + x1*(y2-y3)
    A = A2/2

    B = matrix(c(((y2-y3)/A2), 0, ((y3-y1)/A2), 0, ((y1-y2)/A2), 0,
                 0, ((x3-x2)/A2), 0, ((x1-x3)/A2), 0, ((x2-x1)/A2),
                 ((x3-x2)/A2), ((y2-y3)/A2), ((x1-x3)/A2),((y3-y1)/A2), ((x2-x1)/A2),((y1-y2)/A2)), nrow = 3, byrow = TRUE)

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

    #Element matrix
    Emat1 = (Thick*A)* t(B) %*% D1 %*% B #Element matrix for plane stress
    Emat2 = (Thick*A)* t(B) %*% D2 %*% B #Element matrix for plane strain

    Rlist = list("Pstress" = Emat1, "Pstrain" = Emat2)
    return(Rlist) }

  #Run for acquiring individual element matrix
  m= n= z= o= NROW(meshT) # m=col, n=row, z=element#
  EMPStress = list()
  for (m in 1:z){
    test6 = CST_EM(meshP, meshT, Nu, Y, Thick)
    EMPStress[[m]] = test6$Pstress}

  #Run for acquiring individual element matrix
  m= n= z= o= NROW(meshT) # m=col, n=row, z=element#
  EMPStrain = list()
  for (m in 1:z){
    test6 = CST_EM(meshP, meshT, Nu, Y, Thick)
    EMPStrain[[m]] = test6$Pstrain}

  Rlist = list("EMPStress" = EMPStress, "EMPStrain" = EMPStrain)}
