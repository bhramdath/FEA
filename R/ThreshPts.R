#' @title ThreshPts
#'
#' @description Clean node distribution within or outside of geometry. Optional function for complex geometries.
#'
#' @param coords Nodal coordinates
#' @param thresh Threshold for point removal. Ranges include: 500000-50000000
#' @param edge Coordinate points of the initial geometry.
#'
#' @return Coordinate points of valid nodes.
#' \item{CleanedNodes}{Matrix of new nodes that abide by given threshold rules.}
#' \item{NodeReport}{Report identifying with nodes were kept and which were removed.}
#'
#' @examples
#' # coords = Within
#' # thresh = 500000
#' # edge = Line
#' # ThreshPts(coords, thresh, edge)
#'
#' @export

ThreshPts = function(coords, thresh, edge){
  u = coords
  v = NROW(u)
  dist1 <- as.data.frame(geosphere::dist2Line(u, edge))
  total = numeric(v)
  total1 = numeric(v)
  #rejection method for nodes
  for (m in 1:v){
    total[m] = mean(dist1$distance[m])
    if (total[m] > thresh){total1[m] = NA}
    else {total1[m] = 1}}

  fd1 = total1*u
  fd = na.omit(fd1)
  u0 = as.matrix(rbind(edge, fd), ncol = 2) #present nodes (must add original shape nodes to this)
  plot(u0, pch = 46, type = "p") #polygon coordinates

  Rlist = list("CleanedNodes" = u0, "NodeReport" = fd1)}
