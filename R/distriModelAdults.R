#' distriModelAdults
#'
#' no details
#' 
#' @param vec numeric vector
#' @param k integer value, non linear model (position)
#' @param sizeInmaduros integer value
#' @param hfeno list, objects about phenology
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
distriModelAdults<-function(vec,k, sizeInmaduros,hfeno){  # sl es 4 o 5
  if(k == (sizeInmaduros +1)){
    # probit
    if(hfeno$distri_snh=="probit"){
      WW<-pnorm(log(vec)*hfeno$slope_snh)
    }
    # logit
    if(hfeno$distri_snh=="logit"){
      WW<-1/(1+exp(-(log(vec)*hfeno$slope_snh)))
    }
    # cloglog
    if(hfeno$distri_snh=="cloglog"){
      WW<-1-exp(-exp(log(vec)*hfeno$slope_snh))
    }
  }
  
  if(k == (sizeInmaduros +2)){
    # probit
    if(hfeno$distri_snm=="probit"){
      WW<-pnorm(log(vec)*hfeno$slope_snm)
    }
    # logit
    if(hfeno$distri_snm=="logit"){
      WW<-1/(1+exp(-(log(vec)*hfeno$slope_snm)))
    }
    # cloglog
    if(hfeno$distri_snm=="cloglog"){
      WW<-1-exp(-exp(log(vec)*hfeno$slope_snm))
    }
  }
  return(WW)
}
