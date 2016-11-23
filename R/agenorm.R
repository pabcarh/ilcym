#' agenorm
#'
#' no details
#' 
#' @param estimor numeric vector, estimated coefficients of the nonlinear model
#' @param modelm integer value, non linear model (position)
#' @param xini double, initial value of iteration
#' @param ni integer value, number of iterations 
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
agenorm<-function(estimor,modelm,xini=0.1,ni=20)
{
  ind <- as.list(estimor)
  for (i in names(ind))
  {
    temp <- estimor[[i]]
    storage.mode(temp) <- "double"
    assign(i, temp)
  }
  
  if(modelm==1){  f<-function(x){1-exp(-(a*x+b*x^2+c*x^3))-0.5}; dx1x <- deriv(~ 1-exp(-(a*x+b*x^2+c*x^3))-0.5, "x") }
  if(modelm==2){  f<-function(p){qgamma(p,a,b)}}
  if(modelm==3){  f<-function(x){1/(1+exp(a+b*x))-0.5}; dx1x <- deriv(~ 1/(1+exp(a+b*x))-0.5, "x") }
  if(modelm==4){  f<-function(x){1-exp(-a*x^b)-0.5}; dx1x <- deriv(~ 1-exp(-a*x^b)-0.5, "x") }
  if(modelm==5){  f<-function(x){1-exp(-((x-a)/n)^b)-0.5}; dx1x <- deriv(~ 1-exp(-((x-a)/n)^b)-0.5, "x") }
  if(modelm==6){  f<-function(p){qweibull(p,a,b)}}
  
  ## Aplicamos Newton Raphson
  if(modelm==2 || modelm==6){x=f(0.5)}else
  {
    x=xini
    for(j in 1:ni)
    {
      de = eval(dx1x)
      xd = attr(de,"gradient")[,1]
      x  = x - f(x)/xd	
    }
  }
  
  return(list(xi=x))
}
