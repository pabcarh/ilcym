#' recalc
#'
#' no details
#' 
#' @param ini numeric vector, initial values to model
#' @param niv double, iteration level
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
recalc<-function(ini,niv=1)
{
  vec<-list(1)
  for(i in 1:length(ini))
  {
    vec[[i]]<-ini[[i]]+ini[[i]]*runif(1,-(niv),niv)
  }
  names(vec)=names(ini)
  return(vec)
}
