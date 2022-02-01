#' @title SinglePoly
#'
#' @description Generates a mesh for polygon with a single continuous geometry
#'
#' @param x X-coordinates for geometry.
#' @param y Y-coordinates for geometry.
#' @param ptDS Density of points desired within the geometry.
#' @param ptDL Density of points desired at the perimeter of the geometry.
#'
#' @return Coordinate points of nodes distributed within and on the line of a given geometry.
#' \item{AllCoords}{all coordinate points distributed across the geometry.}
#' \item{Within}{all coordinate points within the geometry ONLY.}
#' \item{Line}{all coordinate points that lay on the perimeter of the geometry ONLY.}
#'
#' @examples
#' x= c(0.23, 0.61, 0.75, 0.61, 0.23, -0.23, -0.61, -0.75, -0.61, -0.23, -0.23)
#' y= c(0.62, 0.38, 0.51, -0.38, -0.62, -0.62, -0.38, -0.51, 0.38, 0.62, 0.65)
#' ptDS = 30
#' ptDL = 20
#' SinglePoly(x, y, ptDS, ptDL)
#'
#'@export

SinglePoly = function(x, y, ptDS, ptDL){
  poly <- sp::Polygon(matrix(rbind(x, y),  nrow = NROW(x), ncol =2, byrow = T))
  pts <- as.data.frame(sp::spsample(poly, n= ptDS, "regular"), pch = 3) #change n to reflect desired point density
  names(pts)[1] <- "x"
  names(pts)[2] <- "y"
  line1 = as.data.frame(sp::spsample(sp::SpatialLines(list(sp::Lines(sp::Line(sp::SpatialPoints(poly)), ID="a"))), n= ptDL, offset =c(0,1), "regular"), pch = 3)
  names(line1)[1] <- "x"
  names(line1)[2] <- "y"
  coords= rbind(pts,line1)
  plot(coords, pch = 46, type = "p") #polygon coordinates

  Rlist = list("AllCoords" = coords, "Within" = pts, "Line" = line1)}
