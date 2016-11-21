#' summ
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

summ<-function(vec){if(length(vec[is.na(vec)])==length(vec)){suma=NA}else(suma=sum(vec,na.rm=T));return(suma)}
