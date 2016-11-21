#' descum
#'
#' no details
#'
#' @param vec vector of counts
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
descum<-function(vec){dc<-c(1);dc[1]=vec[1];for(i in 2:length(vec)){dc[i]=vec[i]-vec[i-1]};return(dc)}
