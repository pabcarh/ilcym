#' media
#'
#' no details
#' 
#' @param vec numeric vector
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
media<-function(vec){if(length(vec[is.na(vec)])==length(vec)){med=0}else{med= mean(vec,na.rm=T)};return(med)}