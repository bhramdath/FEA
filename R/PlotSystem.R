#' @title PlotSystem
#'
#' @description Generates heat map for given stress or strain on the geometry. Threshold values for the color must be assigned.
#'
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param meshT Matrix (3 x n) containing the number of the coordinate point that forms a given triangle within the mesh.
#' @param PlotVal Value to be plotted, either stress or strain, return from function LocalStress function.
#' @param Kol Which plot value column to plot (choose 1, 2, or 3)? 1 = X-stress/strain; 2 = Y-stress/strain; 3 = Shear stress/strain.
#' @param a Threshold for royal blue
#' @param b Threshold for sky blue
#' @param c Threshold for bluegreen
#' @param d Threshold for dark green
#' @param e Threshold for green
#' @param f Threshold for gold
#' @param g Threshold for orange
#' @param h Threshold for pink
#' @param i Threshold for red
#'
#' @return Plot of colored polygon with mesh colored based on the plot value
#'
#' @examples
#' # meshP = MeshPts$p #mesh points
#' # meshT = MeshPts$T #triangulation list
#' # PlotVal = FEAtest$Stress  #output value wanted in plot
#' # Kol =  3 #plotting shear stress/strain
#' # a= 1; b= 2; c=3; d=4; e=5; f=6; g=7; h=8; i=9 #color threshold values from min to max
#' # PlotSystem(meshP, meshT, PlotVal, Kol, a, b, c, d, e, f, g, h, i)
#'
#' @export

#Plot colored elements
PlotSystem = function(meshP, meshT, PlotVal, Kol, a, b, c, d, e, f, g, h, i){
  m= n= o= NROW(meshP)*2
  z= NROW(meshT) # m=col, n=row, z=element#

  ColorMap = matrix(nrow=z, ncol=3)
  for (m in 1:z){
    for (n in 1:3){
      if (PlotVal[m,n] == 0) {ColorMap[m,n] = "purple"}
      else if (PlotVal[m,n] < a) {ColorMap[m,n] = "royalblue"}
      else if (PlotVal[m,n] < b) {ColorMap[m,n] = "skyblue"}
      else if (PlotVal[m,n] < c) {ColorMap[m,n] = "darkturquoise"}
      else if (PlotVal[m,n] < d) {ColorMap[m,n] = "darkolivegreen3"}
      else if (PlotVal[m,n] < e) {ColorMap[m,n] = "forestgreen"}
      else if (PlotVal[m,n] < f) {ColorMap[m,n] = "gold"}
      else if (PlotVal[m,n] < g) {ColorMap[m,n] = "orange"}
      else if (PlotVal[m,n] < h) {ColorMap[m,n] = "deeppink"}
      else if (PlotVal[m,n] < i) {ColorMap[m,n] = "red"}
      else {ColorMap[m,n] = "black"}
    }}

  PolyCo = list()
  PolyColor= function(meshP, meshT){
    x1 = meshP[meshT[m,1],1]
    y1 = meshP[meshT[m,1],2]
    x2 = meshP[meshT[m,2],1]
    y2 = meshP[meshT[m,2],2]
    x3 = meshP[meshT[m,3],1]
    y3 = meshP[meshT[m,3],2]
    sp::Polygons(list(sp::Polygon(cbind(c(x1, x2, x3, x1), c(y1, y2, y3, y1)))), m)}
  for (m in 1:z){PolyCo[[m]]=PolyColor(meshP, meshT)}
  PolyCo2 = sp::SpatialPolygons(PolyCo)
  sp::plot(PolyCo2, col = ColorMap[,Kol]) }
