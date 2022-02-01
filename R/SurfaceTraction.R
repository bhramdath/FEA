#' @title SurfaceTraction.
#'
#' @description Element Surface Traction - generates the column matrix for uniformly distributed surface traction. If surface traction is not present, assign SFTensile and SFShear a value of 0.
#'
#' @param meshP Matrix (2 x n) containing coordinate points of the mesh nodes.
#' @param SFTensile Magnitude of tensile surface traction
#' @param SFShear Magnitude of positive shear traction
#' @param Length Truss length
#' @param Thick Triangle element thickness
#' @param area Triangle element area
#'
#' @return List of element matrices containing surface forces.
#' \item{SurfT}{List of surface forces for each element.}
#'
#' @examples
#' # SFShear = 0
#' # SFTensile = 0
#' # Length = TrussLength
#' # area = Area
#' # SurfaceTraction(meshP, SFTensile, SFShear, Length, Thick, area)
#'
#' @export

SurfaceTraction = function(meshP, SFTensile, SFShear, Length, Thick, area){
  m= n= z= o= NROW(meshP)
  TDOF = 2*z

  SFT = function(TDOF, SFTensile, SFShear, Length, Thick, area){
    L = min(Length[m,])
    b =  Thick
    a = area[m]
    px = SFTensile #magnitude of tensile surface traction
    py = SFShear #magnitude of positive shear traction
    EqLoad = matrix(c(0,0, px*L*b/2, py*L*b/2, px*L*b/2, py*L*b/2), byrow = TRUE)
    return(EqLoad)}

  SurfT = list()
  for (m in 1:z){
    SurfT[[m]] =  SFT(TDOF, SFTensile, SFShear, Length, Thick, area)}
  return(SurfT)}
