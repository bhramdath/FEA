#' @title ManualAdjust
#'
#' @description Allows for manual refinement of the triangular mesh generated based on given parameters. Will remove triangle elements that are identified in the input (loc).
#'
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param meshT Matrix (3 x n) containing the number of the coordinate point that forms a given triangle within the mesh.
#' @param edge Coordinate points of the initial geometry.
#' @param centroid Matrix (2 x n) of triangle elements.
#' @param loc String containing the number of the meshT matrix row of the triangle chosen to be removed.
#'
#' @importFrom graphics text
#' @importFrom stats na.omit
#'
#' @return Generates new mesh and centroid tables
#' \item{Meshpts}{Includes both new mesh coordinate points and triangulation of points.}
#' \item{Centroids}{Centroid positions for each triangle element.}
#'
#' @examples
#' # meshP = MeshPts$p
#' # meshT = MeshPts$T
#' # edge = Line
#' # centroid = tCentroids
#' # loc = c(11, 12, 13) #number of the centroid point of the triangle that is chosen for removal
#' # ManualAdjust(meshP, meshT, edge, centroid, loc)
#'
#' @export

ManualAdjust = function(meshP, meshT, edge, centroid, loc){
  total = numeric(NROW(meshT)); total[1:NROW(meshT)] = 1
  total[loc] = NA
  pd1 = meshT*total
  pd2 = na.omit(pd1)   #Remove triangles outside shape
  Meshpt = geometry::trimesh(pd2, meshP)     #Trimesh with inner regions removed
  centroid = matrix(c(((meshP[pd2[,1],1]+ meshP[pd2[,2],1] + meshP[pd2[,3],1])/3),
                      (meshP[pd2[,1],2]+ meshP[pd2[,2],2] + meshP[pd2[,3],2])/3), ncol = 2) #triangle centroids
  data1 = c(1:NROW(Meshpt$p))
  data2 = c(1:NROW(centroid))
  text(Meshpt$p[,1],Meshpt$p[,2], data1, col = "red", cex= 0.6, font = 2)
  text(centroid[,1],centroid[,2], data2, col = "blue", cex= 0.6, font = 1)

  Rlist = list("MeshPts" = Meshpt, "Centroids" = centroid)}
