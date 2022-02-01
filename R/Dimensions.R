#' @title Dimensions

#' @description Calculates dimensional values for each triangular element, including truss length & angles, distance from nodal point to centroid, aspect ratio of each triangle element, and area of the triangle.
#'
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param meshT Matrix (3 x n) containing the number of the coordinate point that forms a given triangle within the mesh.
#' @param centroid Matrix (2 x n) containing coordinate points of the centroid of each triangular element.
#'
#' @return Evaluation of triangle elements truss, angle, and area.
#' \item{Truss}{Nodal pairs that form each truss.}
#' \item{TrussLength}{Distance between each paired nodes forming a truss, its length.}
#' \item{Dist2Cent}{Shortest distance from truss to triangle centroid.}
#' \item{Truss angle}{Angles of the triangle created from truss meeting.}
#' \item{AspectRatio}{Aspect ratio of triangle elements.}
#' \item{Area}{Area within triangle elements.}
#'
#' @examples
#' # meshP = MeshPts$p #mesh points
#' # meshT = MeshPts$T #triangulation list
#' # centroid = Centroids
#' # Dimensions(meshP, meshT, centroid)
#'
#' @export

Dimensions = function(meshP, meshT, centroid){
  B <- function(meshP, meshT){
    B = cbind(rbind(meshP[meshT[m,1],1],meshP[meshT[m,1],2], meshP[meshT[m,2],1],meshP[meshT[m,2],2]),
              rbind(meshP[meshT[m,2],1],meshP[meshT[m,2],2], meshP[meshT[m,3],1],meshP[meshT[m,3],2]),
              rbind(meshP[meshT[m,3],1],meshP[meshT[m,3],2], meshP[meshT[m,1],1],meshP[meshT[m,1],2]))}
  #Truss lengths
  L <- function(meshP, meshT){
    L = cbind(sqrt(((meshP[meshT[m,1],2] - meshP[meshT[m,2],2])^2) + ((meshP[meshT[m,1],1] - meshP[meshT[m,2],1])^2)),
              sqrt(((meshP[meshT[m,2],2] - meshP[meshT[m,3],2])^2) + ((meshP[meshT[m,2],1] - meshP[meshT[m,3],1])^2)),
              sqrt(((meshP[meshT[m,3],2] - meshP[meshT[m,1],2])^2) + ((meshP[meshT[m,3],1] - meshP[meshT[m,1],1])^2)))}
  #Length from point to centroid
  L2 <- function(meshP, meshT, centroid){
    L2 = cbind(sqrt(((meshP[meshT[m,1],2] - centroid[m,2])^2) + ((meshP[meshT[m,1],1] - centroid[m,1])^2)),
               sqrt(((meshP[meshT[m,2],2] - centroid[m,2])^2) + ((meshP[meshT[m,2],1] - centroid[m,1])^2)),
               sqrt(((meshP[meshT[m,3],2] - centroid[m,2])^2) + ((meshP[meshT[m,3],1] - centroid[m,1])^2)))}

  #Truss angles
  Ang <- function(Length){
    Ang = cbind((180/pi) * acos((Length[m,2]^2 + Length[m,3]^2 - Length[m,1]^2) / (2*Length[m,2]*Length[m,3])),
                (180/pi) * acos((Length[m,1]^2 + Length[m,3]^2 - Length[m,2]^2) / (2*Length[m,1]*Length[m,3])),
                (180/pi) * acos((Length[m,1]^2 + Length[m,2]^2 - Length[m,3]^2) / (2*Length[m,1]*Length[m,2])))}

  #Aspect ratio
  AspectRatio <- function(Length){
    s = (1/2) * (Length[m,1]+Length[m,2]+Length[m,3]) #skewness
    AR = (Length[m,1]*Length[m,2]*Length[m,3])/(8*(s-Length[m,1])*(s-Length[m,2])*(s-Length[m,3]))
    return(AR) }

  v = NROW(centroid)
  Truss = list()   #Total truss coordinates
  Length = matrix(0, nrow = v, ncol = 3, byrow = TRUE)  #Truss lengths
  LengthCent = matrix(0, nrow = v, ncol = 3, byrow = TRUE)  #Length to triangle centroid
  Angle = matrix(0, nrow = v, ncol = 3, byrow = TRUE)  #Truss angles
  AspectR = numeric(v)
  for (m in 1:v){
    Truss[[m]] = B(meshP, meshT)
    Length[m,] = L(meshP, meshT)
    LengthCent[m,] = L2(meshP, meshT, centroid)
    Angle[m,] = Ang(Length)
    AspectR[m] = AspectRatio(Length) }

  #Area of triangle
  Area = function(meshP, meshT){
    v = NROW(meshT)
    S = numeric(v)
    for (m in 1:v){
      x1 = meshP[meshT[m,1],1]
      y1 = meshP[meshT[m,1],2]
      x2 = meshP[meshT[m,2],1]
      y2 = meshP[meshT[m,2],2]
      x3 = meshP[meshT[m,3],1]
      y3 = meshP[meshT[m,3],2]
      S[m] = abs((x1*y2 + x2*y3 + x3*y1 - y1*x2 - y2*x3 - y3*x1)/2)}
    return(S)}
  area = Area(meshP, meshT)


  Rlist = list("Truss" = Truss, "TrussLength" = Length, "Dist2Cent" = LengthCent, "Truss angle" = Angle, "AspectRatio" =AspectR, "Area" = area)}
