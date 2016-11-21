#' pr6mortal
#'
#' no details
#' 
#' @param modelm integer value, non linear model (position)
#' @param datm data frame, temperatures and mortality
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
pr6mortal<-function(modelm,datm)
{
  x<-datm[,1]
  y<-datm[,2]
  if(modelm==1)  modelo<-stats::lm(y~x+I(x^2))
  if(modelm==2)  modelo<-stats::lm(y~x+I(sqrt(x)))
  if(modelm==3)  modelo<-stats::lm(y~x+I(1/sqrt(x)))
  if(modelm==4)  modelo<-stats::lm(y~x+I(1/x^2))
  if(modelm==5)  modelo<-stats::lm(y~x+I(1/x))
  if(modelm==6)  modelo<-stats::lm(y~x+I(log(x)))
  coefi<-round(as.numeric(stats::coef(modelo)),3)
  for(i in 1:length(coefi))if(coefi[i]==0) coefi[i]<-0.001
  coefi<-c("a"=coefi[3],"b"=coefi[2],"c"=coefi[1])
  for (i in names(coefi))
  {
    temp <- coefi[[i]]
    storage.mode(temp) <- "double"
    assign(i, temp)
  }
  salidas<-list(ini=coefi)
  return(salidas)
}
