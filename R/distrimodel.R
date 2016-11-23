#' distrimodel
#'
#' no details
#' 
#' @param vec numeric vector
#' @param sll integer value
#' @param hfeno list, objects about phenology
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
distrimodel<-function(vec,sll,hfeno)
{
  # probit
  if(hfeno$distri_dv[[sll]]=="probit")
  {
    WW<-pnorm(log(vec)*hfeno$slope_dv[[sll]])
  }
  # logit
  if(hfeno$distri_dv[[sll]]=="logit")
  {
    WW<-1/(1+exp(-(log(vec)*hfeno$slope_dv[[sll]])))
  }
  # cloglog
  if(hfeno$distri_dv[[sll]]=="cloglog")
  {
    WW<-1-exp(-exp(log(vec)*hfeno$slope_dv[[sll]]))
  }
  return(WW)
}