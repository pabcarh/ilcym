#' distrimodeldema
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
distrimodeldema<-function(vec,sll,hfeno)
{
  if(sll==1)
  {
    #hembras
    # probit
    if(hfeno$distri_snh=="probit")
    {
      pdd<-pnorm(log(vec)*hfeno$slope_snh)
    }
    # logit
    if(hfeno$distri_snh=="logit")
    {
      pdd<-1/(1+exp(-(log(vec)*hfeno$slope_snh)))
    }
    # cloglog
    if(hfeno$distri_snh=="cloglog")
    {
      pdd<-1-exp(-exp(log(vec)*hfeno$slope_snh))
    }
  }
  if(sll==2)
  {
    #machos
    # probit
    if(hfeno$distri_snm=="probit")
    {
      pdd<-pnorm(log(vec)*hfeno$slope_snm)
    }
    # logit
    if(hfeno$distri_snm=="logit")
    {
      pdd<-1/(1+exp(-(log(vec)*hfeno$slope_snm)))
    }
    # cloglog
    if(hfeno$distri_snm=="cloglog")
    {
      pdd<-1-exp(-exp(log(vec)*hfeno$slope_snm))
    }
  }
  return(pdd)
}
