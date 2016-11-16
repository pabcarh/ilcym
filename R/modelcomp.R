#' modelcomp
#'
#' no details
#' 
#' @param test a character, dicotomic model selected('logit','probit', 'cloglog')
#' @param model a character, dicotomic model ('logit','probit', 'cloglog')
#' @param matri data frame, counts of insect by days and others variables
#' @param parametros numeric vector, the parameters values estimated
#' @param shap.estimados data frame
#' @param coefi numeric vector, linear regression parameters
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
modelcomp<-function(test,model,matri,parametros,shap.estimados,coefi){
  slope<-parametros[length(parametros)]
  est1<-cbind(shap.estimados)
  est1<-as.list(est1)
  for (i in names(est1)) {
    temp <- est1[[i]]
    storage.mode(temp) <- "double"
    assign(i, temp)
  }
  kelvin<-matri[,1]+273.15
  celsius<-matri[,1]
  time<-matri[,2]
  des1<-matri[,5]/matri[,6]
  inteest1<-parametros[1:(length(parametros)-1)]
  slop<-parametros[length(parametros)]
  if(test=="logit"){
    if(model==1)
    {
      meta <- function(kelvin,time,Ha, Hl,Tl,Hh,Th)
      {
        expr <- expression(1/(1+exp(-(log(time*((coefi[1]+coefi[2]*((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))) * (kelvin/(((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th))))) * exp((Ha/1.987) * ((1/((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))) - (1/kelvin))))/
                                            (1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((Hh/1.987) * ((1/Th) - (1/kelvin)))))*slope))))
        eval(expr)
      }
      jmeta<- function(kelvin,time, Ha, Hl,Tl,Hh,Th)
      {
        expr <- expression(1/(1+exp(-(log(time*((coefi[1]+coefi[2]*((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))) * (kelvin/(((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th))))) * exp((Ha/1.987) * ((1/((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))) - (1/kelvin))))/
                                            (1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((Hh/1.987) * ((1/Th) - (1/kelvin)))))*slope))))
        c(eval(D(expr, "Ha"  )),  eval(D(expr, "Hl" )),
          eval(D(expr, "Tl" )), eval(D(expr, "Hh" )), eval(D(expr, "Th" )))
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    if(model==2)
    {
      meta <- function(kelvin,time,To, Ha, Hl,Tl,Hh,Th)
      {
        expr <- expression(1/(1+exp(-(log(time*((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/
                                            (1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((Hh/1.987) * ((1/Th) - (1/kelvin)))))*slope))))
        eval(expr)
      }
      jmeta<- function(kelvin,time,To, Ha, Hl,Tl,Hh,Th)
      {
        expr <- expression(1/(1+exp(-(log(time*((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/
                                            (1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((Hh/1.987) * ((1/Th) - (1/kelvin)))))*slope))))
        c(eval(D(expr, "To")),  eval(D(expr, "Ha"  )),  eval(D(expr, "Hl" )),
          eval(D(expr, "Tl" )), eval(D(expr, "Hh" )), eval(D(expr, "Th" )))
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    if(model==3)
    {
      meta <- function(kelvin,time,p,To, Ha, Hl,Tl,Hh,Th)
      {
        expr <- expression(1/(1+exp(-(log(time*((p) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/
                                            (1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((Hh/1.987) * ((1/Th) - (1/kelvin)))))*slope))))
        eval(expr)
      }
      jmeta<- function(kelvin,time,p,To, Ha, Hl,Tl,Hh,Th)
      {
        expr <- expression(1/(1+exp(-(log(time*((p) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/
                                            (1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((Hh/1.987) * ((1/Th) - (1/kelvin)))))*slope))))
        c(eval(D(expr, "p")),eval(D(expr, "To")),  eval(D(expr, "Ha"  )),  eval(D(expr, "Hl" )),
          eval(D(expr, "Tl" )), eval(D(expr, "Hh" )), eval(D(expr, "Th" )))
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    if(model==4)
    {
      meta <- function(kelvin,time,p,To, Ha)  ## cambian sus parametros
      {
        expr <- expression(1/(1+exp(-(log(time*((p * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))))*slope))))  ## ingreso la funcion dentro de: ..time*(...))*slope..
        eval(expr)
      }
      jmeta<- function(kelvin,time,p,To, Ha)  ## cambian sus parametros
      {
        expr <- expression(1/(1+exp(-(log(time*((p * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))))*slope))))  ## lo mismo aqui
        c(eval(D(expr, "p"  )),  eval(D(expr, "To" )),eval(D(expr, "Ha"  )))
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    
    ######
    if(model==5)
    {
      meta <- function(kelvin,time,p,To, Ha, Tl, Hl) ## cambian sus parametros
      {
        expr <- expression(1/(1+exp(-(log(time*((p * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))))))*slope))))
        eval(expr)
      }
      jmeta<- function(kelvin,time,p,To, Ha, Tl, Hl) ## cambian sus parametros
      {
        expr <- expression(1/(1+exp(-(log(time*((p * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))))))*slope))))
        c(eval(D(expr, "p"  )),  eval(D(expr, "To" )),eval(D(expr, "Ha"  )),eval(D(expr, "Tl"  )),eval(D(expr, "Hl"  )))  ## segun los parametros
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    ######
    if(model==6)
    {
      meta <- function(kelvin,time,p,To, Ha, Hh, Th)  ## cambian sus parametros
      {
        expr <- expression(1/(1+exp(-(log(time*((p * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hh/1.987) * ((1/Th) - (1/kelvin))))))*slope))))  ## ingreso la funcion dentro de: ..time*(...))*slope..
        eval(expr)
      }
      jmeta<- function(kelvin,time,p,To, Ha, Hh, Th)  ## cambian sus parametros
      {
        expr <- expression(1/(1+exp(-(log(time*((p * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hh/1.987) * ((1/Th) - (1/kelvin))))))*slope))))  ## lo mismo aqui
        c(eval(D(expr, "p"  )),  eval(D(expr, "To" )),eval(D(expr, "Ha"  )),eval(D(expr, "Hh"  )),eval(D(expr, "Th"  )))
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    
    ######
    if(model==7)
    {
      meta <- function(kelvin,time,p,To, Ha)  ## cambian sus parametros
      {
        expr <- expression(1/(1+exp(-(log(time*(((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))))*slope))))  ## ingreso la funcion dentro de: ..time*(...))*slope..
        eval(expr)
      }
      jmeta<- function(kelvin,time,p,To, Ha)  ## cambian sus parametros
      {
        expr <- expression(1/(1+exp(-(log(time*(((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))))*slope))))  ## ingreso la funcion dentro de: ..time*(...))*slope..
        c(eval(D(expr, "p"  )),  eval(D(expr, "To" )),eval(D(expr, "Ha"  )))
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    if(model==8){
      meta <- function(kelvin,time,p,To,Ha,Tl,Hl)  ## cambian sus parametros
      {
        expr <- expression(1/(1+exp(-(log(time*(((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))))))*slope))))  ## ingreso la funcion dentro de: ..time*(...))*slope..
        eval(expr)
      }
      jmeta<- function(kelvin,time,p,To, Ha,Tl,Hl)  ## cambian sus parametros
      {
        expr <- expression(1/(1+exp(-(log(time*(((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))))))*slope))))  ## ingreso la funcion dentro de: ..time*(...))*slope..
        c(eval(D(expr, "p"  )),  eval(D(expr, "To" )),eval(D(expr, "Ha"  )),eval(D(expr, "Tl"  )),eval(D(expr, "Hl"  )))  ## segun los parametros
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    if(model==9){
      meta <- function(kelvin,time,p,To,Ha,Hh,Th)  ## cambian sus parametros
      {
        expr <- expression(1/(1+exp(-(log(time*(((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hh/1.987) * ((1/Th) - (1/kelvin))))))*slope))))  ## ingreso la funcion dentro de: ..time*(...))*slope..
        eval(expr)
      }
      jmeta<- function(kelvin,time,p,To,Ha,Hh,Th)  ## cambian sus parametros
      {
        expr <- expression(1/(1+exp(-(log(time*(((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hh/1.987) * ((1/Th) - (1/kelvin))))))*slope))))  ## ingreso la funcion dentro de: ..time*(...))*slope..
        c(eval(D(expr, "p"  )),  eval(D(expr, "To" )),eval(D(expr, "Ha"  )),eval(D(expr, "Hh"  )),eval(D(expr, "Th"  )))  ## segun los parametros
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    
    if(model==10){ #nuevooooooooooooooooooooooo
      
      meta <- function(kelvin,time,p,Ha,Hl,Tl)  ## cambian sus parametros
      {
        expr <- expression(1/(1+exp(-(log(time*((p * (kelvin/(298.16)) * exp((Ha/1.987) * ((1/298.16) - (1/kelvin))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((1/1.987) * ((1/1) - (1/kelvin))))))*slope))))  ## ingreso la funcion dentro de: ..time*(...))*slope..
        eval(expr)
      }
      jmeta<- function(kelvin,time,p,Ha,Hl,Tl)  ## cambian sus parametros
      {
        expr <- expression(1/(1+exp(-(log(time*((p * (kelvin/(298.16)) * exp((Ha/1.987) * ((1/298.16) - (1/kelvin))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((1/1.987) * ((1/1) - (1/kelvin))))))*slope))))  ## ingreso la funcion dentro de: ..time*(...))*slope..
        c(eval(D(expr, "p"  )),eval(D(expr, "Ha"  )),eval(D(expr, "Hl"  )),eval(D(expr, "Tl"  )))  ## segun los parametros
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    if(model==11){
      meta <- function(celsius,time,Tmin,b)
      {
        expr <- expression(1/(1+exp(-(log(time*b*(celsius-Tmin))*slope))))
        eval(expr)
      }
      jmeta<- function(celsius,time,Tmin,b)
      {
        expr <- expression(1/(1+exp(-(log(time*b*(celsius-Tmin))*slope))))
        c(eval(D(expr, "Tmin" )), eval(D(expr, "b" )))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,fcall = meta, jcall = jmeta,celsius = celsius,time=time,des1 = des1)
    }
    
    if(model==12)
    {
      meta <- function(celsius,time,b1,b2,b3,b4,b5)
      {
        expr <- expression(1/(1+exp(-(log(time*(b1*10^(((((celsius-b3)/(b3-b2)-(1/(1+0.28*b4+0.72*log(1+b4))))+exp(b4*((celsius-b3)/(b3-b2)-(1/(1+0.28*b4+0.72*log(1+b4))))))/((1+b4)/(1+1.5*b4+0.39*b4^2)))^2)*(1-b5+b5*((((celsius-b3)/(b3-b2)- (1/(1+0.28*b4+0.72*log(1+b4))))+exp(b4*((celsius-b3)/(b3-b2)- (1/(1+0.28*b4+0.72*log(1+b4))))))/((1+b4)/(1+1.5*b4+0.39*b4^2)))^2)))*slope))))
        eval(expr)
      }
      jmeta<- function(celsius,time,b1,b2,b3,b4,b5)
      {
        expr <- expression(1/(1+exp(-(log(time*(b1*10^(((((celsius-b3)/(b3-b2)-(1/(1+0.28*b4+0.72*log(1+b4))))+exp(b4*((celsius-b3)/(b3-b2)-(1/(1+0.28*b4+0.72*log(1+b4))))))/((1+b4)/(1+1.5*b4+0.39*b4^2)))^2)*(1-b5+b5*((((celsius-b3)/(b3-b2)- (1/(1+0.28*b4+0.72*log(1+b4))))+exp(b4*((celsius-b3)/(b3-b2)- (1/(1+0.28*b4+0.72*log(1+b4))))))/((1+b4)/(1+1.5*b4+0.39*b4^2)))^2)))*slope))))
        c(eval(D(expr, "b1")), eval(D(expr, "b2")), eval(D(expr, "b3"  )), eval(D(expr, "b4" )),
          eval(D(expr, "b5" )))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    if(model==13)
    {
      meta <- function(celsius,time,Y,Tmax, p, v)
      {
        expr <- expression(1/(1+exp(-(log(time*(Y*(exp(p*celsius)-exp(p*Tmax-(Tmax-celsius)/v))))*slope))))
        eval(expr)
      }
      jmeta<- function(celsius,time,Y,Tmax, p, v)
      {
        expr <- expression(1/(1+exp(-(log(time*(Y*(exp(p*celsius)-exp(p*Tmax-(Tmax-celsius)/v))))*slope))))
        c(eval(D(expr, "Y")),eval(D(expr, "Tmax")), eval(D(expr, "p")), eval(D(expr, "v")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    if(model==14)
    {
      meta <- function(celsius,time,alfa,k,Tmax, p, v)
      {
        expr <- expression(1/(1+exp(-(log(time*(alfa*((1/(1+k*exp(-p*celsius)))-exp(-(Tmax-celsius)/v))))*slope))))
        eval(expr)
      }
      jmeta<- function(celsius,time,alfa,k,Tmax, p, v)
      {
        expr <- expression(1/(1+exp(-(log(time*(alfa*((1/(1+k*exp(-p*celsius)))-exp(-(Tmax-celsius)/v))))*slope))))
        c(eval(D(expr, "alfa")),eval(D(expr, "k")),eval(D(expr, "Tmax")), eval(D(expr, "p")), eval(D(expr, "v")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    
    if(model==15)
    {
      meta <- function(celsius,time,aa,To,Tmax)
      {
        expr <- expression(1/(1+exp(-(log(time*(aa*celsius*(celsius-To)*(Tmax-celsius)^0.5))*slope))))
        eval(expr)
      }
      jmeta<- function(celsius,time,aa,To,Tmax)
      {
        expr <- expression(1/(1+exp(-(log(time*(aa*celsius*(celsius-To)*(Tmax-celsius)^0.5))*slope))))
        c(eval(D(expr, "aa")), eval(D(expr, "To")), eval(D(expr, "Tmax")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    if(model==16)
    {
      meta <- function(celsius,time,aa,To,Tmax,d)
      {
        expr <- expression(1/(1+exp(-(log(time*(aa*celsius*(celsius-To)*(Tmax-celsius)^d))*slope)))) #REVISAR!
        eval(expr)
      }
      jmeta<- function(celsius,time,aa,To,Tmax,d)
      {
        expr <- expression(1/(1+exp(-(log(time*(aa*celsius*(celsius-To)*(Tmax-celsius)^d))*slope)))) #REVISAR!
        c(eval(D(expr, "d")),eval(D(expr, "To")),eval(D(expr, "aa")),eval(D(expr, "Tmax")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    
    if(model==17)
    {
      meta <- function(celsius,time,Rmax,Topc,k1,k2)
      {
        expr <- expression(1/(1+exp(-(log(time*(Rmax*(1+exp(k1+k2*(Topc))))/(1+exp(k1+k2*celsius)))*slope))))
        eval(expr)
      }
      jmeta <- function(celsius,time,Rmax,Topc,k1,k2)
      {
        expr <- expression(1/(1+exp(-(log(time*(Rmax*(1+exp(k1+k2*(Topc))))/(1+exp(k1+k2*celsius)))*slope))))
        c(eval(D(expr, "Rmax")), eval(D(expr, "k1")), eval(D(expr, "k2")), eval(D(expr, "Topc")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    
    
    if(model==18)
    {
      meta <- function(celsius,time,d,Y,Tmax, v)
      {
        expr <- expression(1/(1+exp(-(log(time*(Y*((celsius)^2/((celsius)^2+d^2)-exp(-(Tmax-(celsius))/v))))*slope))))
        eval(expr)
      }
      jmeta<- function(celsius,time,To,d,Y,Tmax, v)
      {
        expr <- expression(1/(1+exp(-(log(time*(Y*((celsius)^2/((celsius)^2+d^2)-exp(-(Tmax-(celsius))/v))))*slope))))
        c(eval(D(expr, "d")),eval(D(expr, "Y")),eval(D(expr, "Tmax")),eval(D(expr, "v")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    
    if(model==19){
      meta <- function(celsius,time,Tl, p, dt,L)
      {
        expr <- expression(1/(1+exp(-(log(time*(exp(p*celsius)-exp(-(p*Tl-(celsius-Tl))/dt)+L))*slope))))
        eval(expr)
      }
      jmeta <- function(celsius,time,Tl, p, dt,L)
      {
        expr <- expression(1/(1+exp(-(log(time*(exp(p*celsius)-exp(-(p*Tl-(celsius-Tl))/dt)+L))*slope))))
        c(eval(D(expr, "Tl")),eval(D(expr, "L")), eval(D(expr, "p")), eval(D(expr, "dt")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time),est1)))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    
  }
  if(test=="probit"){
    if(model==1){
      meta <- function(kelvin,time,Ha, Hl,Tl,Hh,Th)
      {
        expr <- expression(pnorm(log(time*((coefi[1]+coefi[2]*((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))) * (kelvin/(((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th))))) * exp((Ha/1.987) * ((1/((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))) - (1/kelvin))))/
                                       (1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((Hh/1.987) * ((1/Th) - (1/kelvin)))))*slope))
        eval(expr)
      }
      jmeta<- function(kelvin,time, Ha, Hl,Tl,Hh,Th)
      {
        expr <- expression(pnorm(log(time*((coefi[1]+coefi[2]*((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))) * (kelvin/(((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th))))) * exp((Ha/1.987) * ((1/((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))) - (1/kelvin))))/
                                       (1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((Hh/1.987) * ((1/Th) - (1/kelvin)))))*slope))
        c( eval(D(expr, "Ha"  )),  eval(D(expr, "Hl" )),
           eval(D(expr, "Tl" )), eval(D(expr, "Hh" )), eval(D(expr, "Th" )))
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    if(model==2){
      meta <- function(kelvin,time,To, Ha, Hl,Tl,Hh,Th)
      {
        expr <- expression(pnorm(log(time*((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/
                                       (1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((Hh/1.987) * ((1/Th) - (1/kelvin)))))*slope))
        eval(expr)
      }
      jmeta<- function(kelvin,time,To, Ha, Hl,Tl,Hh,Th)
      {
        expr <- expression(pnorm(log(time*((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/
                                       (1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((Hh/1.987) * ((1/Th) - (1/kelvin)))))*slope))
        c(eval(D(expr, "To")),  eval(D(expr, "Ha"  )),  eval(D(expr, "Hl" )),
          eval(D(expr, "Tl" )), eval(D(expr, "Hh" )), eval(D(expr, "Th" )))
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    if(model==3){
      meta <- function(kelvin,time,p,To, Ha, Hl,Tl,Hh,Th)
      {
        expr <- expression(pnorm(log(time*((p) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/
                                       (1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((Hh/1.987) * ((1/Th) - (1/kelvin)))))*slope))
        eval(expr)
      }
      jmeta<- function(kelvin,time,p,To, Ha, Hl,Tl,Hh,Th)
      {
        expr <- expression(pnorm(log(time*((p) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/
                                       (1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((Hh/1.987) * ((1/Th) - (1/kelvin)))))*slope))
        c(eval(D(expr, "p")), eval(D(expr, "To")),  eval(D(expr, "Ha"  )),  eval(D(expr, "Hl" )),
          eval(D(expr, "Tl" )), eval(D(expr, "Hh" )), eval(D(expr, "Th" )))
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    if(model==4){
      meta <- function(kelvin,time,p,To, Ha)  ## fijate los parametros
      {
        expr <- expression(pnorm(log(time*((p * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))))*slope))   ## ingreso la funcion dentro de: ..time*(...))*slope..
        eval(expr)
      }
      jmeta<- function(kelvin,time,p,To, Ha)
      {
        expr <- expression(pnorm(log(time*((p * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))))*slope))   ## ingreso la funcion dentro de: ..time*(...))*slope..
        c(eval(D(expr, "p"  )),  eval(D(expr, "To" )),eval(D(expr, "Ha"  ))) ## igual q el de logit
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    
    ######
    if(model==5){
      meta <- function(kelvin,time,p,To, Ha, Tl, Hl)
      {
        expr <- expression(pnorm(log(time*((p * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))))))*slope))
        eval(expr)
      }
      jmeta<- function(kelvin,time,p,To, Ha, Tl, Hl)
      {
        expr <- expression(pnorm(log(time*((p * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))))))*slope))
        c(eval(D(expr, "p"  )),  eval(D(expr, "To" )),eval(D(expr, "Ha"  )),eval(D(expr, "Tl"  )),eval(D(expr, "Hl"  )))
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    ######
    if(model==6){
      meta <- function(kelvin,time,p,To, Ha, Hh, Th)
      {
        expr <- expression(pnorm(log(time*((p * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hh/1.987) * ((1/Th) - (1/kelvin))))))*slope))
        eval(expr)
      }
      jmeta<- function(kelvin,time,p,To, Ha, Hh, Th)
      {
        expr <- expression(pnorm(log(time*((p * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hh/1.987) * ((1/Th) - (1/kelvin))))))*slope))
        c(eval(D(expr, "p"  )),  eval(D(expr, "To" )),eval(D(expr, "Ha"  )),eval(D(expr, "Hh"  )),eval(D(expr, "Th"  )))
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    
    
    ######
    if(model==7){
      meta <- function(kelvin,time,p,To, Ha)  ## fijate los parametros
      {
        expr <- expression(pnorm(log(time*(((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))))*slope))   ## ingreso la funcion dentro de: ..time*(...))*slope..
        eval(expr)
      }
      jmeta<- function(kelvin,time,p,To, Ha)
      {
        expr <- expression(pnorm(log(time*(((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))))*slope))   ## ingreso la funcion dentro de: ..time*(...))*slope..
        c(eval(D(expr, "p"  )),  eval(D(expr, "To" )),eval(D(expr, "Ha"  ))) ## igual q el de logit
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    if(model==8){
      meta <- function(kelvin,time,p,To, Ha,Tl,Hl)  ## fijate los parametros
      {
        expr <- expression(pnorm(log(time*(((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))))))*slope))   ## ingreso la funcion dentro de: ..time*(...))*slope..
        eval(expr)
      }
      jmeta<- function(kelvin,time,p,To, Ha,Tl,Hl)
      {
        expr <- expression(pnorm(log(time*(((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))))))*slope))   ## ingreso la funcion dentro de: ..time*(...))*slope..
        c(eval(D(expr, "p"  )),  eval(D(expr, "To" )),eval(D(expr, "Ha"  )),eval(D(expr, "Tl"  )),eval(D(expr, "Hl"  )))  ## segun los parametros
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    if(model==9){
      meta <- function(kelvin,time,p,To,Ha,Hh,Th)  ## cambian sus parametros
      {
        expr <- expression(pnorm(log(time*(((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hh/1.987) * ((1/Th) - (1/kelvin))))))*slope))   ## ingreso la funcion dentro de: ..time*(...))*slope..
        eval(expr)
      }
      jmeta <- function(kelvin,time,p,To,Ha,Hh,Th)  ## cambian sus parametros
      {
        expr <- expression(pnorm(log(time*(((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hh/1.987) * ((1/Th) - (1/kelvin))))))*slope))   ## ingreso la funcion dentro de: ..time*(...))*slope..
        c(eval(D(expr, "p"  )),  eval(D(expr, "To" )),eval(D(expr, "Ha"  )),eval(D(expr, "Hh"  )),eval(D(expr, "Th"  )))  ## segun los parametros
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    if(model==10){
      meta <- function(kelvin,time,p,Ha,Hl,Tl)  ## cambian sus parametros
      {
        expr <- expression(pnorm(log(time*((p * (kelvin/(298.16)) * exp((Ha/1.987) * ((1/298.16) - (1/kelvin))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((1/1.987) * ((1/1) - (1/kelvin))))))*slope))   ## ingreso la funcion dentro de: ..time*(...))*slope..
        eval(expr)
      }
      jmeta <- function(kelvin,time,p,Ha,Hl,Tl)  ## cambian sus parametros
      {
        expr <- expression(pnorm(log(time*((p * (kelvin/(298.16)) * exp((Ha/1.987) * ((1/298.16) - (1/kelvin))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((1/1.987) * ((1/1) - (1/kelvin))))))*slope))   ## ingreso la funcion dentro de: ..time*(...))*slope..
        c(eval(D(expr, "p"  )),eval(D(expr, "Ha"  )),eval(D(expr, "Hl"  )),eval(D(expr, "Tl"  )))  ## segun los parametros
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    if(model==11){
      meta <- function(celsius,time,Tmin,b)
      {
        expr <- expression(pnorm(log(time*b*(celsius-Tmin))*slope))
        eval(expr)
      }
      jmeta<- function(celsius,time,Tmin,b)
      {
        expr <- expression(pnorm(log(time*b*(celsius-Tmin))*slope))
        c(eval(D(expr, "Tmin" )), eval(D(expr, "b" )))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    
    if(model==12){
      meta <- function(celsius,time,b1,b2,b3,b4,b5)
      {
        expr <- expression(pnorm(log(time*(b1*10^(((((celsius-b3)/(b3-b2)-(1/(1+0.28*b4+0.72*log(1+b4))))+exp(b4*((celsius-b3)/(b3-b2)-(1/(1+0.28*b4+0.72*log(1+b4))))))/((1+b4)/(1+1.5*b4+0.39*b4^2)))^2)*(1-b5+b5*((((celsius-b3)/(b3-b2)- (1/(1+0.28*b4+0.72*log(1+b4))))+exp(b4*((celsius-b3)/(b3-b2)- (1/(1+0.28*b4+0.72*log(1+b4))))))/((1+b4)/(1+1.5*b4+0.39*b4^2)))^2)))*slope))
        eval(expr)
      }
      jmeta<- function(celsius,time,b1,b2,b3,b4,b5)
      {
        expr <- expression(pnorm(log(time*(b1*10^(((((celsius-b3)/(b3-b2)-(1/(1+0.28*b4+0.72*log(1+b4))))+exp(b4*((celsius-b3)/(b3-b2)-(1/(1+0.28*b4+0.72*log(1+b4))))))/((1+b4)/(1+1.5*b4+0.39*b4^2)))^2)*(1-b5+b5*((((celsius-b3)/(b3-b2)- (1/(1+0.28*b4+0.72*log(1+b4))))+exp(b4*((celsius-b3)/(b3-b2)- (1/(1+0.28*b4+0.72*log(1+b4))))))/((1+b4)/(1+1.5*b4+0.39*b4^2)))^2)))*slope))
        c(eval(D(expr, "b1")), eval(D(expr, "b2")), eval(D(expr, "b3"  )), eval(D(expr, "b4" )),
          eval(D(expr, "b5" )))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    
    if(model==13){
      meta <- function(celsius,time,Y,Tmax, p, v)
      {
        expr <- expression(pnorm(log(time*(Y*(exp(p*celsius)-exp(p*Tmax-(Tmax-celsius)/v))))*slope))
        eval(expr)
      }
      jmeta<- function(celsius,time,Y,Tmax, p, v)
      {
        expr <- expression(pnorm(log(time*(Y*(exp(p*celsius)-exp(p*Tmax-(Tmax-celsius)/v))))*slope))
        c(eval(D(expr, "Y")),eval(D(expr, "Tmax")), eval(D(expr, "p")), eval(D(expr, "v")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    if(model==14){
      meta <- function(celsius,time,alfa,k,Tmax, p, v)
      {
        expr <- expression(pnorm(log(time*(alfa*((1/(1+k*exp(-p*celsius)))-exp(-(Tmax-celsius)/v))))*slope))
        eval(expr)
      }
      jmeta<- function(celsius,time,alfa,k,Tmax, p, v)
      {
        expr <- expression(pnorm(log(time*(alfa*((1/(1+k*exp(-p*celsius)))-exp(-(Tmax-celsius)/v))))*slope))
        c(eval(D(expr, "alfa")),eval(D(expr, "k")),eval(D(expr, "Tmax")), eval(D(expr, "p")), eval(D(expr, "v")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    
    if(model==15){
      meta <- function(celsius,time,aa,To,Tmax)
      {
        expr <- expression(pnorm(log(time*(aa*celsius*(celsius-To)*(Tmax-celsius)^0.5))*slope))
        eval(expr)
      }
      jmeta<- function(celsius,time,aa,To,Tmax)
      {
        expr <- expression(pnorm(log(time*(aa*celsius*(celsius-To)*(Tmax-celsius)^0.5))*slope))
        c(eval(D(expr, "aa")), eval(D(expr, "To")), eval(D(expr, "Tmax")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    if(model==16){
      meta <- function(celsius,time,aa,To,Tmax,d)
      {
        expr <- expression(pnorm(log(time*(aa*celsius*(celsius-To)*(Tmax-celsius)^d))*slope)) #REVISAR!
        eval(expr)
      }
      jmeta<- function(celsius,time,aa,To,Tmax,d)
      {
        expr <- expression(pnorm(log(time*(aa*celsius*(celsius-To)*(Tmax-celsius)^d))*slope)) #REVISAR!
        c(eval(D(expr, "d")),eval(D(expr, "To")),eval(D(expr, "aa")),eval(D(expr, "Tmax")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    
    if(model==17){
      meta <- function(celsius,time,Rmax,Topc,k1,k2)
      {
        expr <- expression(pnorm(log(time*(Rmax*(1+exp(k1+k2*(Topc))))/(1+exp(k1+k2*celsius)))*slope))
        eval(expr)
      }
      jmeta <- function(celsius,time,Rmax,Topc,k1,k2)
      {
        expr <- expression(pnorm(log(time*(Rmax*(1+exp(k1+k2*(Topc))))/(1+exp(k1+k2*celsius)))*slope))
        c(eval(D(expr, "Rmax")), eval(D(expr, "k1")), eval(D(expr, "k2")), eval(D(expr, "Topc")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    
    if(model==18){
      meta <- function(celsius,time,d,Y,Tmax, v)
      {
        expr <- expression(pnorm(log(time*(Y*((celsius)^2/((celsius)^2+d^2)-exp(-(Tmax-(celsius))/v))))*slope))
        eval(expr)
      }
      jmeta<- function(celsius,time,d,Y,Tmax, v)
      {
        expr <- expression(pnorm(log(time*(Y*((celsius)^2/((celsius)^2+d^2)-exp(-(Tmax-(celsius))/v))))*slope))
        c(eval(D(expr, "d")),eval(D(expr, "Y")),eval(D(expr, "Tmax")),eval(D(expr, "v")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    
    if(model==19){
      meta <- function(celsius,time,Tl, p, dt,L)
      {
        expr <- expression(pnorm(log(time*(exp(p*celsius)-exp(-(p*Tl-(celsius-Tl))/dt)+L))*slope))
        eval(expr)
      }
      jmeta <- function(celsius,time,Tl, p, dt,L)
      {
        expr <- expression(pnorm(log(time*(exp(p*celsius)-exp(-(p*Tl-(celsius-Tl))/dt)+L))*slope))
        c(eval(D(expr, "Tl")),eval(D(expr, "L")), eval(D(expr, "p")), eval(D(expr, "dt")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    
  }
  if(test=="cloglog"){
    if(model==1){
      meta <- function(kelvin,time,Ha, Hl,Tl,Hh,Th)
      {
        expr <- expression(1-exp(-exp(-(-log(-log(0.5)))+(log(time*((coefi[1]+coefi[2]*((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))) * (kelvin/(((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th))))) * exp((Ha/1.987) * ((1/((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))) - (1/kelvin))))/
                                                                (1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((Hh/1.987) * ((1/Th) - (1/kelvin)))))*slope))))
        eval(expr)
      }
      jmeta<- function(kelvin,time, Ha, Hl,Tl,Hh,Th)
      {
        expr <- expression(1-exp(-exp(-(-log(-log(0.5)))+(log(time*((coefi[1]+coefi[2]*((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))) * (kelvin/(((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th))))) * exp((Ha/1.987) * ((1/((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))) - (1/kelvin))))/
                                                                (1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((Hh/1.987) * ((1/Th) - (1/kelvin)))))*slope))))
        c(eval(D(expr, "Ha"  )),  eval(D(expr, "Hl" )),
          eval(D(expr, "Tl" )), eval(D(expr, "Hh" )), eval(D(expr, "Th" )))
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    if(model==2){
      meta <- function(kelvin,time,To, Ha, Hl,Tl,Hh,Th)
      {
        expr <- expression(1-exp(-exp(-(-log(-log(0.5)))+(log(time*((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/
                                                                (1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((Hh/1.987) * ((1/Th) - (1/kelvin)))))*slope))))
        eval(expr)
      }
      jmeta<- function(kelvin,time,To, Ha, Hl,Tl,Hh,Th)
      {
        expr <- expression(1-exp(-exp(-(-log(-log(0.5)))+(log(time*((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/
                                                                (1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((Hh/1.987) * ((1/Th) - (1/kelvin)))))*slope))))
        c(eval(D(expr, "To")),  eval(D(expr, "Ha"  )),  eval(D(expr, "Hl" )),
          eval(D(expr, "Tl" )), eval(D(expr, "Hh" )), eval(D(expr, "Th" )))
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    if(model==3){
      meta <- function(kelvin,time,p,To, Ha, Hl,Tl,Hh,Th)
      {
        expr <- expression(1-exp(-exp(-(-log(-log(0.5)))+(log(time*((p) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/
                                                                (1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((Hh/1.987) * ((1/Th) - (1/kelvin)))))*slope))))
        eval(expr)
      }
      jmeta<- function(kelvin,time,p,To, Ha, Hl,Tl,Hh,Th)
      {
        expr <- expression(1-exp(-exp(-(-log(-log(0.5)))+(log(time*((p) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/
                                                                (1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((Hh/1.987) * ((1/Th) - (1/kelvin)))))*slope))))
        c(eval(D(expr, "p")), eval(D(expr, "To")),  eval(D(expr, "Ha"  )),  eval(D(expr, "Hl" )),
          eval(D(expr, "Tl" )), eval(D(expr, "Hh" )), eval(D(expr, "Th" )))
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    if(model==4){
      meta <- function(kelvin,time,p,To, Ha)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*((p * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))))*slope))))   ## ingreso la funcion dentro de: ..time*(...))*slope..
        eval(expr)
      }
      jmeta<- function(kelvin,time,p,To, Ha)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*((p * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))))*slope))))   ## ingreso la funcion dentro de: ..time*(...))*slope..
        c(eval(D(expr, "p"  )),  eval(D(expr, "To" )),eval(D(expr, "Ha"  ))) 
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    if(model==5){
      meta <- function(kelvin,time,p,To, Ha, Tl, Hl)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*((p * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))))))*slope))))
        eval(expr)
      }
      jmeta<- function(kelvin,time,p,To, Ha, Tl, Hl)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*((p * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))))))*slope))))
        c(eval(D(expr, "p"  )),  eval(D(expr, "To" )),eval(D(expr, "Ha"  )),eval(D(expr, "Tl"  )),eval(D(expr, "Hl"  )))
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    if(model==6){
      meta <- function(kelvin,time,p,To, Ha, Hh, Th)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*((p * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hh/1.987) * ((1/Th) - (1/kelvin))))))*slope))))
        eval(expr)
      }
      jmeta<- function(kelvin,time,p,To, Ha, Hh, Th)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*((p * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hh/1.987) * ((1/Th) - (1/kelvin))))))*slope))))
        c(eval(D(expr, "p"  )),  eval(D(expr, "To" )),eval(D(expr, "Ha"  )),eval(D(expr, "Hh"  )),eval(D(expr, "Th"  )))
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    if(model==7){
      meta <- function(kelvin,time,p,To, Ha)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))))*slope))))   ## ingreso la funcion dentro de: ..time*(...))*slope..
        eval(expr)
      }
      jmeta<- function(kelvin,time,p,To, Ha)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))))*slope))))   ## ingreso la funcion dentro de: ..time*(...))*slope..
        c(eval(D(expr, "p"  )),  eval(D(expr, "To" )),eval(D(expr, "Ha"  ))) 
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    if(model==8){
      meta <- function(kelvin,time,p,To, Ha,Tl,Hl)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))))))*slope))))   ## ingreso la funcion dentro de: ..time*(...))*slope..
        eval(expr)
      }
      jmeta<- function(kelvin,time,p,To, Ha,Tl,Hl)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))))))*slope))))   ## ingreso la funcion dentro de: ..time*(...))*slope..
        c(eval(D(expr, "p"  )),  eval(D(expr, "To" )),eval(D(expr, "Ha"  )),eval(D(expr, "Tl"  )),eval(D(expr, "Hl"  )))  ## segun los parametros
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    if(model==9){
      meta <- function(kelvin,time,p,To,Ha,Hh,Th)  ## cambian sus parametros
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hh/1.987) * ((1/Th) - (1/kelvin))))))*slope))))   ## ingreso la funcion dentro de: ..time*(...))*slope..
        eval(expr)
      }
      jmeta <- function(kelvin,time,p,To,Ha,Hh,Th)  ## cambian sus parametros
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(((coefi[1]+coefi[2]*To) * (kelvin/(To)) * exp((Ha/1.987) * ((1/To) - (1/kelvin))))/(1 + exp((Hh/1.987) * ((1/Th) - (1/kelvin))))))*slope))))   ## ingreso la funcion dentro de: ..time*(...))*slope..
        c(eval(D(expr, "p"  )),  eval(D(expr, "To" )),eval(D(expr, "Ha"  )),eval(D(expr, "Hh"  )),eval(D(expr, "Th"  )))  ## segun los parametros
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    if(model==10){
      meta <- function(kelvin,time,p,Ha,Hl,Tl) ## cambian sus parametros
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*((p * (kelvin/(298.16)) * exp((Ha/1.987) * ((1/298.16) - (1/kelvin))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((1/1.987) * ((1/1) - (1/kelvin))))))*slope))))   ## ingreso la funcion dentro de: ..time*(...))*slope..
        eval(expr)
      }
      jmeta <- function(kelvin,time,p,Ha,Hl,Tl)  ## cambian sus parametros
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*((p * (kelvin/(298.16)) * exp((Ha/1.987) * ((1/298.16) - (1/kelvin))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/kelvin))) + exp((1/1.987) * ((1/1) - (1/kelvin))))))*slope))))   ## ingreso la funcion dentro de: ..time*(...))*slope..
        c(eval(D(expr, "p"  )),eval(D(expr, "Ha"  )),eval(D(expr, "Hl"  )),eval(D(expr, "Tl"  )))  ## segun los parametros
      }
      fcn     <- function(est1, kelvin,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(kelvin = kelvin,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          kelvin = kelvin,time=time,des1 = des1)
    }
    
    if(model==11){
      meta <- function(celsius,time,Tmin,b)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*b*(celsius-Tmin))*slope))))
        eval(expr)
      }
      jmeta<- function(celsius,time,Tmin,b)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*b*(celsius-Tmin))*slope))))
        c(eval(D(expr, "Tmin" )), eval(D(expr, "b" )))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    if(model==12){
      meta <- function(celsius,time,b1,b2,b3,b4,b5)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(b1*10^(((((celsius-b3)/(b3-b2)-(1/(1+0.28*b4+0.72*log(1+b4))))+exp(b4*((celsius-b3)/(b3-b2)-(1/(1+0.28*b4+0.72*log(1+b4))))))/((1+b4)/(1+1.5*b4+0.39*b4^2)))^2)*(1-b5+b5*((((celsius-b3)/(b3-b2)- (1/(1+0.28*b4+0.72*log(1+b4))))+exp(b4*((celsius-b3)/(b3-b2)- (1/(1+0.28*b4+0.72*log(1+b4))))))/((1+b4)/(1+1.5*b4+0.39*b4^2)))^2)))*slope))))
        eval(expr)
      }
      jmeta<- function(celsius,time,b1,b2,b3,b4,b5)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(b1*10^(((((celsius-b3)/(b3-b2)-(1/(1+0.28*b4+0.72*log(1+b4))))+exp(b4*((celsius-b3)/(b3-b2)-(1/(1+0.28*b4+0.72*log(1+b4))))))/((1+b4)/(1+1.5*b4+0.39*b4^2)))^2)*(1-b5+b5*((((celsius-b3)/(b3-b2)- (1/(1+0.28*b4+0.72*log(1+b4))))+exp(b4*((celsius-b3)/(b3-b2)- (1/(1+0.28*b4+0.72*log(1+b4))))))/((1+b4)/(1+1.5*b4+0.39*b4^2)))^2)))*slope))))
        c(eval(D(expr, "b1")), eval(D(expr, "b2")), eval(D(expr, "b3"  )), eval(D(expr, "b4" )),
          eval(D(expr, "b5" )))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    if(model==13){
      meta <- function(celsius,time,Y,Tmax, p, v)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(Y*(exp(p*celsius)-exp(p*Tmax-(Tmax-celsius)/v))))*slope))))
        eval(expr)
      }
      jmeta<- function(celsius,time,Y,Tmax, p, v)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(Y*(exp(p*celsius)-exp(p*Tmax-(Tmax-celsius)/v))))*slope))))
        c(eval(D(expr, "Y")),eval(D(expr, "Tmax")), eval(D(expr, "p")), eval(D(expr, "v")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    if(model==14){
      meta <- function(celsius,time,alfa,k,Tmax, p, v)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(alfa*((1/(1+k*exp(-p*celsius)))-exp(-(Tmax-celsius)/v))))*slope))))
        eval(expr)
      }
      jmeta<- function(celsius,time,alfa,k,Tmax, p, v)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(alfa*((1/(1+k*exp(-p*celsius)))-exp(-(Tmax-celsius)/v))))*slope))))
        c(eval(D(expr, "alfa")),eval(D(expr, "k")),eval(D(expr, "Tmax")), eval(D(expr, "p")), eval(D(expr, "v")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    
    if(model==15){
      meta <- function(celsius,time,aa,To,Tmax)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(aa*celsius*(celsius-To)*(Tmax-celsius)^0.5))*slope))))
        eval(expr)
      }
      jmeta<- function(celsius,time,aa,To,Tmax)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(aa*celsius*(celsius-To)*(Tmax-celsius)^0.5))*slope))))
        c(eval(D(expr, "aa")), eval(D(expr, "To")), eval(D(expr, "Tmax")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    if(model==16){
      meta <- function(celsius,time,aa,To,Tmax,d)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(aa*celsius*(celsius-To)*(Tmax-celsius)^d))*slope)))) #REVISAR!
        eval(expr)
      }
      jmeta<- function(celsius,time,aa,To,Tmax,d)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(aa*celsius*(celsius-To)*(Tmax-celsius)^d))*slope)))) #REVISAR!
        c(eval(D(expr, "d")),eval(D(expr, "To")),eval(D(expr, "aa")),eval(D(expr, "Tmax")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    
    if(model==17){
      meta <- function(celsius,time,Rmax,Topc,k1,k2)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(Rmax*(1+exp(k1+k2*(Topc))))/(1+exp(k1+k2*celsius)))*slope))))
        eval(expr)
      }
      jmeta <- function(celsius,time,Rmax,Topc,k1,k2)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(Rmax*(1+exp(k1+k2*(Topc))))/(1+exp(k1+k2*celsius)))*slope))))
        c(eval(D(expr, "Rmax")), eval(D(expr, "k1")), eval(D(expr, "k2")), eval(D(expr, "Topc")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    
    if(model==18)
    {
      meta <- function(celsius,time,d,Y,Tmax, v)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(Y*((celsius)^2/((celsius)^2+d^2)-exp(-(Tmax-(celsius))/v))))*slope))))
        eval(expr)
      }
      jmeta<- function(celsius,time,To,d,Y,Tmax, v)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(Y*((celsius)^2/((celsius)^2+d^2)-exp(-(Tmax-(celsius))/v))))*slope))))
        c(eval(D(expr, "d")),eval(D(expr, "Y")),eval(D(expr, "Tmax")),eval(D(expr, "v")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    
    if(model==19){
      
      meta <- function(celsius,time,Tl, p, dt,L)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(exp(p*celsius)-exp(-(p*Tl-(celsius-Tl))/dt)+L))*slope))))
        eval(expr)
      }
      jmeta <- function(celsius,time,Tl, p, dt,L)
      {
        expr <- expression(1-exp(-exp(-log(-log(0.5)))+((log(time*(exp(p*celsius)-exp(-(p*Tl-(celsius-Tl))/dt)+L))*slope))))
        c(eval(D(expr, "Tl")),eval(D(expr, "L")), eval(D(expr, "p")), eval(D(expr, "dt")))
      }
      fcn     <- function(est1, celsius,time, des1, fcall, jcall)
        (des1 - do.call("fcall", c(list(celsius = celsius,time=time), as.list(est1))))
      metamodel <- nls.lm(par = est1, fn = fcn,
                          fcall = meta, jcall = jmeta,
                          celsius = celsius,time=time,des1 = des1)
    }
    
  }
  stderro<-diag(ginv(metamodel$hessian))
  estimate<-as.data.frame(metamodel$par)
  slope<-slope
  estimados<-as.data.frame(metamodel$par)
  meta<-meta
  ind <- as.list(estimate)
  for (i in names(ind))
  {
    temp <- ind[[i]]
    storage.mode(temp) <- "double"
    assign(i, temp)
  }
  
  if(model==1) salida<-list(estimados=estimate,slope=slope,estmshape=estimados,meta=meta,p=(coefi[1]+coefi[2]*(Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th))),To=(Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)),stderro=stderro) 
  if(model==2 || model==7 || model==8 || model==9) salida<-list(estimados=estimate,slope=slope,estmshape=estimados,meta=meta,p=(coefi[1]+coefi[2]*To),To=To,stderro=stderro) 
  if(model==3 || model==4 || model==5 || model==6 || model==10 || model==11 || model==12 || model==13 || model==14 || model==15 || model==16 || model==17 || model==18 || model==19)
  {salida<-list(estimados=estimate,slope=slope,estmshape=estimados,meta=meta,stderro=stderro)}
  
  return(salida)
}
