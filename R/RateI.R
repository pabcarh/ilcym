#' RateI
#'
#' no details
#' 
#' @param vec numeric vector
#' @param Table2 data frame, temperatures (tmin and tmax)
#' @param Ki integer value, state position
#' @param parametrosc numeric vector, estimated coefficients of the nonlinear model about the development rate
#' @param parametrosm numeric vector, estimated coefficients of the nonlinear model about the mortality
#' @param funciont formula, development rate equation
#' @param funcionm formula, mortality equation
#' @param nmax double
#' @param steps integer value, between 1 and 48
#' @param J integer value, state evaluated (1: , 2:)
#' @param modelim integer value, position of the nonlinear model (development rate)

#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
RateI<-function(vec,Table2,Ki,parametrosc,parametrosm=NULL,funciont,funcionm=NULL,nmax,steps,J=NA,modelim)
{
  for (i in names(parametrosc)){temp <- parametrosc[i];storage.mode(temp) <- "double";assign(i, temp)}
  i=0:(steps-1)
  T1=((vec[3]-vec[2])/2)*cos(pi*(i+0.5)/steps) + (vec[3]+vec[2])/2
  if(vec[1]!=nmax){T2<-((vec[3]-Table2[vec[1]+1,2])/2)*cos(pi*(i+0.5)/steps) + (vec[3]+Table2[vec[1]+1,2])/2}else{
    T2<-((vec[3]-Table2[1,2])/2)*cos(pi*(i+0.5)/steps) + (vec[3]+Table2[1,2])/2}
  if((modelim[Ki]==1 || modelim[Ki]==2 || modelim[Ki]==3 ||
      modelim[Ki]==4 || modelim[Ki]==5 || modelim[Ki]==6 ||
      modelim[Ki]==7 || modelim[Ki]==8 || modelim[Ki]==9 ||
      modelim[Ki]==10 || modelim[Ki]==11 || modelim[Ki]==12 || modelim[Ki]==13 || modelim[Ki]==14) & is.na(J)){x = T1 + 273.15;x2 = T2 + 273.15}else{x = T1;x2 = T2}
  
  
  Rat=eval(funciont);if(!is.na(J)){Rat[Rat<=0]=0}
  Ratetot1<-(sum(Rat))/steps 
  x=x2;Rat=eval(funciont);if(!is.na(J)){Rat[Rat<=0]=0}
  Ratetot2<-(sum(Rat))/steps    
  Rate=(Ratetot1+Ratetot2)/2
  
  if(!is.null(funcionm))
  {
    for (i in names(parametrosm)){temp <- parametrosm[i];storage.mode(temp) <- "double";assign(i, temp)}
    x=T1;M1=eval(funcionm);M1[M1>1]=1; Mortality1 <- (sum(M1))/steps    
    x=T2;M2=eval(funcionm);M2[M2>1]=1; Mortality2 <- (sum(M2))/steps    
    Mortality=(Mortality1+Mortality2)/2
    return(c(Rate1=Ratetot1,Rate2=Rate,Mortality1=Mortality1,Mortality2=Mortality))
  }else{return(c(Rate1=Ratetot1,Rate2=Rate))}
}
