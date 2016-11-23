#' distrimodeldeim
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
distrimodeldeim<-function(vec,sll,hfeno) ## this is GML function
{
  # probit
  if(hfeno$distri_dv[[sll]]=="probit")
  {
    pdd<-pnorm(log(vec)*hfeno$slope_dv[[sll]])
  }
  # logit
  if(hfeno$distri_dv[[sll]]=="logit")
  {
    pdd<-1/(1+exp(-(log(vec)*hfeno$slope_dv[[sll]])))
  }
  # cloglog
  if(hfeno$distri_dv[[sll]]=="cloglog")
  {
    pdd<-1-exp(-exp(log(vec)*hfeno$slope_dv[[sll]]))
  }
  return(pdd)
}
