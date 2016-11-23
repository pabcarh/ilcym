#' p.area
#'
#' no details
#' 
#' @param ilat numeric vector, range of the latitude
#' @param R double, number of regions to divide
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
p.area<-function(ilat,R){
  R=R+1
  lats=seq(ilat[1],ilat[2],length.out=R)+0.0000000001
  mat.lat=matrix(NA,R-1,2)
  for(i in 2:R) 
  {
    mat.lat[i-1,]=c(lats[i-1],lats[i])
  }
  return(mat.lat)
}