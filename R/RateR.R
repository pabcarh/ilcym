#' RateR
#'
#' no details
#' 
#' @param vec numeric vector
#' @param Table2 data frame, temperatures (tmin and tmax)
#' @param nmax double
#' @param steps integer value, between 1 and 48
#' @param ff formula, equation
#' @param paramR numeric vector, estimated coefficients of the nonlinear model

#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
RateR<-function(vec,Table2,nmax,steps,ff,paramR)
{
  for (i in names(paramR)){temp <- paramR[i];storage.mode(temp) <- "double";assign(i, temp)}
  i=0:(steps-1) 
  T1=((vec[3]-vec[2])/2)*cos(pi*(i+0.5)/steps) + (vec[3]+vec[2])/2
  if(vec[1]!=nmax){T2<-((vec[3]-Table2[vec[1]+1,2])/2)*cos(pi*(i+0.5)/steps) + (vec[3]+Table2[vec[1]+1,2])/2}else{
    T2<-((vec[3]-Table2[1,2])/2)*cos(pi*(i+0.5)/steps) + (vec[3]+Table2[1,2])/2}
  
  x=T1;r1=eval(ff[[3]]);r1[r1 < 0]=0;r1[r1 > 1]=1
  R1=(sum(r1))/steps
  x=T2;r1=eval(ff[[3]]);r1[r1 < 0]=0;r1[r1 > 1]=1
  R2=(sum(r1))/steps
  Rt=(R1+R2)/2
  return(Rt)
}
