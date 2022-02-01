#' @title triangulate0
#'
#' @description Triangulation by Delaunayn algorithm. Automatically generates a triangular mesh for a geometry containing nodal points.
#'
#' @param u0 Matrix (2 x n) of node coordinates within the geometry.
#' @param edge Matrix (2 x n) of coordinate points on the perimeter of the geometry.
#'
#' @return Produces data for generated mesh.
#' \item{Meshpts}{Includes both new mesh coordinate points and triangulation of points.}
#' \item{Centroids}{Centroid positions for each triangle element.}
#'
#' @examples
#' # u0 = CleanedNodes
#' # edge = Line
#' # triangulate0(u0, edge)
#'
#' @export

triangulate0 = function(u0, edge){
  t = geometry::delaunayn(u0) #list of triangle nodes
  Tpmid=matrix(c(((u0[t[,1],1]+ u0[t[,2],1] + u0[t[,3],1])/3),
                 (u0[t[,1],2]+ u0[t[,2],2] + u0[t[,3],2])/3), ncol = 2) #triangle centroids

  #Graphical output of initial mesh
  Meshpt = geometry::trimesh(t, u0)

  Meshpt = geometry::trimesh(t, u0)
  data1 = c(1:NROW(Meshpt$p))
  data2 = c(1:NROW(Meshpt$T))
  text(Meshpt$p[,1],Meshpt$p[,2], data1, col = "red", cex= 0.6, font = 2)
  text(Tpmid[,1],Tpmid[,2], data2, col = "blue", cex= 0.6, font = 1)

  Rlist = list("MeshPts" = Meshpt, "Centroids" = Tpmid)}
