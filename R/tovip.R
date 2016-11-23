#' tovip
#'
#' no details
#' 
#' @param v1 numeric vector
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
tovip<-function(v1)
{
  v1=na.omit(v1);v1=v1[v1!=0];tovip=length(v1)
  return(tovip)
}
