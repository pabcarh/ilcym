#' ovi_func
#'
#' no details
#' 
#' @param modelm integer value, non linear model (position)
#' @param ovifemal data frame, temperatures and oviposition
#' @param alg a character, algorithm
#' @param inim numeric vector of the initial values for the modeling
#' @param pesos numeric vector
#' @param weights numeric vector
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
ovi_func <- function (modelm, ovifemal, alg, inim,pesos,weights)
{
  x<-ovifemal[,2]
  y<-ovifemal[,3]
  if(alg=="Newton")
  {
    if(modelm==1)
    {
      form<-y~(1-exp(-(a*x+b*x^2+c*x^3)))
      frm<-"O(E) = 1-exp(-aE-bE^2-cE^3))"
      if(pesos=="WLS") modelo<-nls(form,inim,data=ovifemal,weights=weights)
      if(pesos=="LS") modelo<-nls(form,inim,data=ovifemal)
      coefi<-as.data.frame(t(coef(modelo)))
      return(list(estimados=coefi,f=form,frm=frm))
    }
    if(modelm==2)
    {
      form<-y~pgamma(x,a,b)
      frm<-"O(E) = pgamma(E,a,b)"
      if(pesos=="WLS")  modelo<-nls(form,inim,data=ovifemal,weights=weights)
      if(pesos=="LS") modelo<-nls(form,inim,data=ovifemal)
      coefi<-as.data.frame(t(coef(modelo)))
      return(list(estimados=coefi,f=form,frm=frm))
    }
    if(modelm==3)
    {
      form<-y~1/(1+exp(a+b*x))
      frm<- "O(E) = 1/(1+exp(a+bE))"
      if(pesos=="WLS")modelo<-nls(form,inim,data=ovifemal,weights=weights)
      if(pesos=="LS")modelo<-nls(form,inim,data=ovifemal)
      coefi<-as.data.frame(t(coef(modelo)))
      return(list(estimados=coefi,frm=frm,f=form))
    }
    if(modelm==4) #para datos alrederor de cero, como heydi
    {
      form<-y~1-exp(-a*x^b)
      frm= "O(E) = 1-exp(-aE^b)"
      if(pesos=="WLS")  modelo<-nls(form,inim,data=ovifemal,weights=weights)
      if(pesos=="LS")  modelo<-nls(form,inim,data=ovifemal)
      coefi<-as.data.frame(t(coef(modelo)))
      return(list(estimados=coefi,frm=frm,f=form))
    }
    if(modelm==5) #para datos alrederor de cero, como heydi
    {
      form<-y~1-exp(-((x-a)/n)^b)
      frm= "O(E) = 1-exp(-((x-a)/n)^b)"
      if(pesos=="WLS")  modelo<-nls(form,inim,data=ovifemal,weights=weights)
      if(pesos=="LS")  modelo<-nls(form,inim,data=ovifemal)
      coefi<-as.data.frame(t(coef(modelo)))
      return(list(estimados=coefi,frm=frm,f=form))
    }
    if(modelm==6)
    {
      form<-y~pweibull(x,a,b)
      frm<-"O(E) = pweibull(E,a,b)"
      if(pesos=="WLS")  modelo<-nls(form,inim,data=ovifemal,weights=weights)
      if(pesos=="LS") modelo<-nls(form,inim,data=ovifemal)
      coefi<-as.data.frame(t(coef(modelo)))
      return(list(estimados=coefi,f=form,frm=frm))
    }
    
  }
  
  if(alg=="Marquardtr")
  {
    if(modelm==1)
    {
      x<-ovifemal[,2]
      y<-ovifemal[,3]
      ini<-as.list(inim)
      ind <- as.list(ini)
      for (i in names(ind)) {
        temp <- ini[[i]]
        storage.mode(temp) <- "double"
        assign(i, temp)
      }
      f <- function(x,a,b,c)
      {
        expr <- expression(1 - exp(-(a * x + b * x^2 + c * x^3)))
        eval(expr)
      }
      j <- function(x,a,b,c)
      {
        expr <- expression(1 - exp(-(a * x + b * x^2 + c * x^3)))
        c(eval(D(expr, "a")),eval(D(expr, "b")),eval(D(expr, "c")))
      }
      fcn     <- function(ini, x, y, fcall, jcall)
        (y - do.call("fcall", c(list(x = x), as.list(ini))))
      out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
      estimate<-as.data.frame(out$par)
      frm<-"O(E) = 1 - exp(-(a.x + b.x^2 + c.x^3)))"
    }
    
    if(modelm==2)
    {
      x<-ovifemal[,2]
      y<-ovifemal[,3]
      ini<-as.list(inim)
      ind <- as.list(ini)
      for (i in names(ind)) {
        temp <- ini[[i]]
        storage.mode(temp) <- "double"
        assign(i, temp)
      }
      f <- function(x,a,b)
      {
        expr <- expression(pgamma(x,a,b))
        eval(expr)
      }
      j <- function(x,a,b)
      {
        expr <- expression(pgamma(x,a,b))
        c(eval(D(expr, "a")),eval(D(expr, "b")))
      }
      fcn     <- function(ini, x, y, fcall, jcall)
        (y - do.call("fcall", c(list(x = x), as.list(ini))))
      out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
      estimate<-as.data.frame(out$par)
      frm<-"O(E) = pgamma(x,a,b)"
    }
    
    if(modelm==3)
    {
      x<-ovifemal[,2]
      y<-ovifemal[,3]
      ini<-as.list(inim)
      ind <- as.list(ini)
      for (i in names(ind)) {
        temp <- ini[[i]]
        storage.mode(temp) <- "double"
        assign(i, temp)
      }
      f <- function(x,a,b)
      {
        expr <- expression(1/(1+exp(a+b*x)))
        eval(expr)
      }
      j <- function(x,a,b)
      {
        expr <- expression(1/(1+exp(a+b*x)))
        c(eval(D(expr, "a")),eval(D(expr, "b")))
      }
      fcn     <- function(ini, x, y, fcall, jcall)
        (y - do.call("fcall", c(list(x = x), as.list(ini))))
      out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
      estimate<-as.data.frame(out$par)
      frm<-"O(E) = 1/(1+exp(a+b.x))"
    }
    
    
    if(modelm==4)
    {
      x<-ovifemal[,2]
      y<-ovifemal[,3]
      ini<-as.list(inim)
      ind <- as.list(ini)
      for (i in names(ind)) {
        temp <- ini[[i]]
        storage.mode(temp) <- "double"
        assign(i, temp)
      }
      f <- function(x,a,b)
      {
        expr <- expression(1-exp(-a*x^b))
        eval(expr)
      }
      j <- function(x,a,b)
      {
        expr <- expression(1-exp(-a*x^b))
        c(eval(D(expr, "a")),eval(D(expr, "b")))
      }
      fcn     <- function(ini, x, y, fcall, jcall)
        (y - do.call("fcall", c(list(x = x), as.list(ini))))
      out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
      estimate<-as.data.frame(out$par)
      frm<-"O(E) = 1-exp(-a.x^b)"
    }
    
    
    if(modelm==5)
    {
      x<-ovifemal[,2]
      y<-ovifemal[,3]
      ini<-as.list(inim)
      ind <- as.list(ini)
      for (i in names(ind)) {
        temp <- ini[[i]]
        storage.mode(temp) <- "double"
        assign(i, temp)
      }
      f <- function(x,a,b,n)
      {
        expr <- expression(1-exp(-((x-a)/n)^b))
        eval(expr)
      }
      j <- function(x,a,b,n)
      {
        expr <- expression(1-exp(-((x-a)/n)^b))
        c(eval(D(expr, "a")),eval(D(expr, "b")),eval(D(expr, "n")))
      }
      fcn     <- function(ini, x, y, fcall, jcall)
        (y - do.call("fcall", c(list(x = x), as.list(ini))))
      out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
      estimate<-as.data.frame(out$par)
      frm<-"O(E) = 1-exp(-((x-a)/n)^b)"
    }
    
    if(modelm==6)
    {
      x<-ovifemal[,2]
      y<-ovifemal[,3]
      ini<-as.list(inim)
      ind <- as.list(ini)
      for (i in names(ind)) {
        temp <- ini[[i]]
        storage.mode(temp) <- "double"
        assign(i, temp)
      }
      f <- function(x,a,b)
      {
        expr <- expression(pweibull(x,a,b))
        eval(expr)
      }
      j <- function(x,a,b)
      {
        expr <- expression(pweibull(x,a,b))
        c(eval(D(expr, "a")),eval(D(expr, "b")))
      }
      fcn     <- function(ini, x, y, fcall, jcall)
        (y - do.call("fcall", c(list(x = x), as.list(ini))))
      out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
      estimate<-as.data.frame(out$par)
      frm<-"O(E) = pweibull(x,a,b)"
    }
    
    stderro<-diag(ginv(out$hessian))
    return(list(estimados=estimate,f=f,frm=frm,stdovi=stderro))
  }
}
