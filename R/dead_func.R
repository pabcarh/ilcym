#' dead_func
#'
#' no details
#' 
#' @param proc a character, type of variable (mortal, taza)
#' @param modelm integer value, non linear model (position)
#' @param datm data frame, temperatures and mortality
#' @param alg a character, algorithm
#' @param inim numeric vector, initial values to model
#' @param pesos numeric vector
#' @param weights numeric vector
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
dead_func <- function (proc,modelm, datm, alg, inim, pesos,weights){
  x<-datm[,1]
  y<-datm[,2]
  if(alg=="Newton"){
    if(modelm==1)
    {
      form<-y~a*x^2+b*x+c
      if(proc=="mortal") frm<- "m(T) = aT?+bT+c" else frm<- "f(T) = aT?+bT+c"
    }
    if(modelm==2)
    {
      form<-y~a*sqrt(x)+b*x+c
      if(proc=="mortal") frm<- "m(T) = a?sqrt(T)+bT+c" else frm<- "f(T) = a?sqrt(T)+bT+c"
    }
    if(modelm==3)
    {
      form<-y~a*(1/sqrt(x))+b*x+c
      if(proc=="mortal") frm<- "m(T) = a/sqrt(T)+bT+c" else frm<- "f(T) = a/sqrt(T)+bT+c"
    }
    if(modelm==4)
    {
      form<-y~a*(1/x^2)+b*x+c
      if(proc=="mortal") frm<- "m(T) = a/T?+bT+c" else frm<- "f(T) = a/T?+bT+c"
    }
    if(modelm==5)
    {
      form<-y~a*(1/x)+b*x+c
      if(proc=="mortal") frm<- "m(T) = a/T+bT+c" else frm<- "f(T) = a/T+bT+c"
    }
    if(modelm==6)
    {
      form<-y~b*x+a*log(x)+c
      if(proc=="mortal") frm<- "m(T) = a?ln(T)+bT+c" else frm<- "f(T) = a?ln(T)+bT+c"
    }
    if(pesos=="WLS") modelo<-nls(form,inim,data=datm,weights=weights)
    if(pesos=="LS") modelo<-nls(form,inim,data=datm)
    coefi<-coef(modelo)
    return(list(estimados=coefi,f=form,modelo=modelo,ecua=frm))
  }
  if(alg=="Marquardtr")
  {
    ini<-as.list(inim)
    for (i in names(ini)) {
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    if(modelm==9)
    {
      f <- function(x,a,b,y0,x0)
      {
        expr <- expression(y0 + a * exp(-0.5 * ((x-x0)/b)^2))
        eval(expr)
      }
      j <- function(x,a,b,y0,x0)
      {
        expr <- expression(y0 + a * exp(-0.5 * ((x-x0)/b)^2))
        c(eval(D(expr, "a")),eval(D(expr, "b")), eval(D(expr, "x0")), eval(D(expr, "y0")))
      }
      g<-y~y0 + a * exp(-0.5 * ((x-x0)/b)^2)
      if(proc=="mortal") frm<- "m(T) = y0+a?exp(-0.5((x-x0)/b)?)" else frm<- "f(T) = y0+a?exp(-0.5((x-x0)/b)?)"
    }
    if(modelm==7)
    {
      f <- function(x,a,b,c,d)
      {
        expr <- expression(1/(1+a*exp(-b*((x-c)/d)^2)))
        eval(expr)
      }
      j <- function(x,a,b,c,d)
      {
        expr <- expression(1/(1+a*exp(-b*((x-c)/d)^2)))
        c(eval(D(expr, "a")),eval(D(expr, "b")), eval(D(expr, "c")), eval(D(expr, "d")))
      }
      g<-y~1/(1+a*exp(-b*((x-c)/d)^2))
      if(proc=="mortal") frm<- "m(T) = 1/(1+a?exp(-b((x-c)/d)?))" else frm<- "f(T) = 1/(1+a?exp(-b((x-c)/d)?))"
    }
    if(modelm==8)
    {
      f <- function(x,a,b,c,xo)
      {
        expr <- expression((a*exp(-b*((x-xo)/c)^2)))
        eval(expr)
      }
      j <- function(x,a,b,c,xo)
      {
        expr <- expression((a*exp(-b*((x-xo)/c)^2)))
        c(eval(D(expr, "a")),eval(D(expr, "b")), eval(D(expr, "c")), eval(D(expr, "xo")))
      }
      g<-y~a*exp(-b*((x-xo)/c)^2)
      if(proc=="mortal") frm<- "m(T) = a?exp(-b((x-xo)/c)?)" else frm<- "f(T) = a?exp(-b((x-xo)/c)?)"
    }
    if(modelm==10){
      f <- function(x,y0,a,b,x0)
      {
        expr <- expression(y0 + a * exp(-0.5 * (log(abs(x/x0))/b)^2))
        eval(expr)
      }
      j <- function(x,y0,a,b,x0)
      {
        expr <- expression(y0 + a * exp(-0.5 * (log(abs(x/x0))/b)^2))
        c(eval(D(expr, "y0")),eval(D(expr, "a")), eval(D(expr, "b")), eval(D(expr, "x0")))
      }
      g<-y~y0 + a * exp(-0.5 * (log(abs(x/x0))/b)^2)
      if(proc=="mortal") frm<- "m(T) = y0+a?exp(-0.5(log(abs(x/x0))/b)?)" else frm<- "f(T) = y0+a?exp(-0.5(log(abs(x/x0))/b)?)"
    }
    
    if(modelm==11){
      f <- function(x,b1,b2,b3,d)
      {
        expr <- expression(b1+b2*x+b3*x^d)
        eval(expr)
      }
      j <- function(x,b1,b2,b3,d)
      {
        expr <- expression(b1+b2*x+b3*x^d)
        c(eval(D(expr, "b1")),eval(D(expr, "b2")), eval(D(expr, "b3")), eval(D(expr, "d")))
      }
      g<-y ~ b1+b2*x+b3*x^d
      if(proc=="mortal") frm<- "m(T) = b1+b2.x+b3.x^d " else frm<- "f(T) = b1+b2.x+b3.x^d"
    }
    
    
    if(modelm==12){
      f <- function(x,b1,b2,b3)
      {
        expr <- expression(exp(b1+b2*x+b3*x^2))
        eval(expr)
      }
      j <- function(x,b1,b2,b3)
      {
        expr <- expression(exp(b1+b2*x+b3*x^2))
        c(eval(D(expr, "b1")),eval(D(expr, "b2")), eval(D(expr, "b3")))
      }
      g<-y ~ exp(b1+b2*x+b3*x^2)
      if(proc=="mortal") frm<- "m(T) = exp(b1+b2.x+b3.x?) " else frm<- "f(T) = exp(b1+b2.x+b3.x?)"
    }
    
    if(modelm==13){
      f <- function(x,b1,b2,b3,b4,b5)
      {
        expr <- expression(1-(b4/(1+b5*exp(b1+b2*x+b3*x^2))))
        eval(expr)
      }
      j <- function(x,b1,b2,b3,b4,b5)
      {
        expr <- expression(1-(b4/(1+b5*exp(b1+b2*x+b3*x^2))))
        c(eval(D(expr, "b1")),eval(D(expr, "b2")), eval(D(expr, "b3")), eval(D(expr, "b4")), eval(D(expr, "b5")))
      }
      g<-y ~ 1-(b4/(1+b5*exp(b1+b2*x+b3*x^2)))
      if(proc=="mortal") frm<- "m(T) = 1-(b4/(1+b5.exp(b1+b2.x+b3.x?)))" else frm<- "f(T) = 1-(b4/(1+b5.exp(b1+b2.x+b3.x?)))"
    }
    
    if(modelm==14){
      f <- function(x,b1,b2,b3)
      {
        expr <- expression(exp(b1+b2*x+b3*sqrt(x)))
        eval(expr)
      }
      j <- function(x,b1,b2,b3)
      {
        expr <- expression(exp(b1+b2*x+b3*sqrt(x)))
        c(eval(D(expr, "b1")),eval(D(expr, "b2")), eval(D(expr, "b3")))
      }
      g<-y ~ exp(b1+b2*x+b3*sqrt(x))
      if(proc=="mortal") frm<- "m(T) = exp(b1+b2.x+b3.(x)^0.5)" else frm<- "f(T) = exp(b1+b2.x+b3.(x)^0.5)"
    }
    
    
    
    if(modelm==15){
      f <- function(x,b1,b2,b3,b4,b5)
      {
        expr <- expression(1-(b4/(1+b5*exp(b1+b2*x+b3*sqrt(x)))))
        eval(expr)
      }
      j <- function(x,b1,b2,b3,b4,b5)
      {
        expr <- expression(1-(b4/(1+b5*exp(b1+b2*x+b3*sqrt(x)))))
        c(eval(D(expr, "b1")),eval(D(expr, "b2")), eval(D(expr, "b3")), eval(D(expr, "b4")), eval(D(expr, "b5")))
      }
      g<-y ~ exp(1-(b4/(1+b5*exp(b1+b2*x+b3*sqrt(x)))))
      if(proc=="mortal") frm<- "m(T) = 1-(b4/(1+b5.exp(b1+b2.x+b3.(x)^0.5)))" else frm<- "f(T) = 1-(b4/(1+b5.exp(b1+b2.x+b3.(x)^0.5)))"
    }
    
    
    if(modelm==16){
      f <- function(x,b1,b2,b3)
      {
        expr <- expression(exp(b1+b2*x+b3*(1/sqrt(x))))
        eval(expr)
      }
      j <- function(x,b1,b2,b3)
      {
        expr <- expression(exp(b1+b2*x+b3*(1/sqrt(x))))
        c(eval(D(expr, "b1")),eval(D(expr, "b2")), eval(D(expr, "b3")))
      }
      g<-y ~ exp(b1+b2*x+b3*(1/sqrt(x)))
      if(proc=="mortal") frm<- "m(T) = exp(b1+b2.x+b3.(1/(x)^0.5))" else frm<- "f(T) = exp(b1+b2.x+b3.(1/(x)^0.5))"
    }
    
    
    
    if(modelm==17){
      f <- function(x,b1,b2,b3,b4,b5)
      {
        expr <- expression(1-(b4/(1+b5*exp(b1+b2*x+b3*(1/sqrt(x))))))
        eval(expr)
      }
      j <- function(x,b1,b2,b3,b4,b5)
      {
        expr <- expression(1-(b4/(1+b5*exp(b1+b2*x+b3*(1/sqrt(x))))))
        c(eval(D(expr, "b1")),eval(D(expr, "b2")), eval(D(expr, "b3")), eval(D(expr, "b4")), eval(D(expr, "b5")))
      }
      g<-y ~ 1-(b4/(1+b5*exp(b1+b2*x+b3*(1/sqrt(x)))))
      if(proc=="mortal") frm<- "m(T) = 1-(b4/(1+b5.exp(b1+b2.x+b3.(1/(x)^0.5))))" else frm<- "f(T) = 1-(b4/(1+b5.exp(b1+b2.x+b3.(1/(x)^0.5))))"
    }
    
    if(modelm==18){
      f <- function(x,b1,b2,b3)
      {
        expr <- expression(exp(b1+b2*x+b3*(1/x)))
        eval(expr)
      }
      j <- function(x,b1,b2,b3)
      {
        expr <- expression(exp(b1+b2*x+b3*(1/x)))
        c(eval(D(expr, "b1")),eval(D(expr, "b2")), eval(D(expr, "b3")))
      }
      g<-y ~ exp(b1+b2*x+b3*(1/x))
      if(proc=="mortal") frm<- "m(T) = exp(b1+b2.x+b3.(1/x))" else frm<- "f(T) = exp(b1+b2.x+b3.(1/x))"
    }
    
    
    
    if(modelm==19){
      f <- function(x,b1,b2,b3,b4,b5)
      {
        expr <- expression(1-(b4/(1+b5*exp(b1+b2*x+b3*(1/x)))))
        eval(expr)
      }
      j <- function(x,b1,b2,b3,b4,b5)
      {
        expr <- expression(1-(b4/(1+b5*exp(b1+b2*x+b3*(1/x)))))
        c(eval(D(expr, "b1")),eval(D(expr, "b2")), eval(D(expr, "b3")), eval(D(expr, "b4")), eval(D(expr, "b5")))
      }
      g<-y ~ 1-(b4/(1+b5*exp(b1+b2*x+b3*(1/x))))
      if(proc=="mortal") frm<- "m(T) = 1-(b4/(1+b5.exp(b1+b2.x+b3.(1/x))))" else frm<- "f(T) = 1-(b4/(1+b5.exp(b1+b2.x+b3.(1/x))))"
    }
    
    
    if(modelm==20){
      f <- function(x,b1,b2,b3,d)
      {
        expr <- expression(exp(b1+b2*x+b3*x^d))
        eval(expr)
      }
      j <- function(x,b1,b2,b3,d)
      {
        expr <- expression(exp(b1+b2*x+b3*x^d))
        c(eval(D(expr, "b1")),eval(D(expr, "b2")), eval(D(expr, "b3")),eval(D(expr, "d")))
      }
      g<-y ~ exp(b1+b2*x+b3*x^d)
      if(proc=="mortal") frm<- "m(T) = exp(b1+b2.x+b3.x^d)" else frm<- "f(T) = exp(b1+b2.x+b3.x^d)"
    }
    
    if(modelm==21){
      f <- function(x,b1,b2,b3,b4,b5,d)
      {
        expr <- expression(1-(b4/(1+b5*exp(b1+b2*x+b3*x^d))))
        eval(expr)
      }
      j <- function(x,b1,b2,b3,b4,b5,d)
      {
        expr <- expression(1-(b4/(1+b5*exp(b1+b2*x+b3*x^d))))
        c(eval(D(expr, "b1")),eval(D(expr, "b2")), eval(D(expr, "b3")), eval(D(expr, "b4")), eval(D(expr, "b5")),eval(D(expr, "d")))
      }
      g<-y ~ 1-(b4/(1+b5*exp(b1+b2*x+b3*x^d)))
      if(proc=="mortal") frm<- "m(T) = 1-(b4/(1+b5.exp(b1+b2.x+b3.x^d)))" else frm<- "f(T) = 1-(b4/(1+b5.exp(b1+b2.x+b3.x^d)))"
    }
    
    
    if(modelm==22){
      f <- function(x,b1,b2,b3)
      {
        expr <- expression(exp(b1+b2*x+b3*log(x)))
        eval(expr)
      }
      j <- function(x,b1,b2,b3)
      {
        expr <- expression(exp(b1+b2*x+b3*log(x)))
        c(eval(D(expr, "b1")),eval(D(expr, "b2")), eval(D(expr, "b3")))
      }
      g<-y ~ exp(b1+b2*x+b3*log(x))
      if(proc=="mortal") frm<- "m(T) = exp(b1+b2.x+b3.ln(x))" else frm<- "f(T) = exp(b1+b2.x+b3.ln(x))"
    }
    
    
    
    if(modelm==23){
      f <- function(x,b1,b2,b3,b4,b5)
      {
        expr <- expression(1-(b4/(1+b5*exp(b1+b2*x+b3*log(x)))))
        eval(expr)
      }
      j <- function(x,b1,b2,b3,b4,b5)
      {
        expr <- expression(1-(b4/(1+b5*exp(b1+b2*x+b3*log(x)))))
        c(eval(D(expr, "b1")),eval(D(expr, "b2")), eval(D(expr, "b3")), eval(D(expr, "b4")), eval(D(expr, "b5")))
      }
      g<-y ~ 1-(b4/(1+b5*exp(b1+b2*x+b3*log(x))))
      if(proc=="mortal") frm<- "m(T) = 1-(b4/(1+b5.exp(b1+b2.x+b3.ln(x))))" else frm<- "f(T) = 1-(b4/(1+b5.exp(b1+b2.x+b3.ln(x))))"
    }
    
    
    
    if(modelm==24){
      f <- function(x,rm,Topt,Troh)
      {
        expr <- expression(1-rm*exp((-0.5)*(-(x-Topt)/Troh)^2))
        eval(expr)
      }
      j <- function(x,rm,Topt,Troh)
      {
        expr <- expression(1-rm*exp((-0.5)*(-(x-Topt)/Troh)^2))
        c(eval(D(expr, "rm")),eval(D(expr, "Topt")), eval(D(expr, "Troh")))
      }
      g<-y ~ 1-rm*exp((-0.5)*(-(x-Topt)/Troh)^2)
      if(proc=="mortal") frm<- "m(T) = 1-rm.exp((-0.5).(-(x-Topt)/Troh)^2)" else frm<- "f(T) = 1-rm.exp((-0.5).(-(x-Topt)/Troh)^2)"
    }
    
    
    if(modelm==25){
      f <- function(x,rm,Topt,Troh)
      {
        expr <- expression(1-rm*exp((-0.5)*(-(log(x)-log(Topt))/Troh)^2))
        eval(expr)
      }
      j <- function(x,rm,Topt,Troh)
      {
        expr <- expression(1-rm*exp((-0.5)*(-(log(x)-log(Topt))/Troh)^2))
        c(eval(D(expr, "rm")),eval(D(expr, "Topt")), eval(D(expr, "Troh")))
      }
      g<-y ~ 1-rm*exp((-0.5)*(-(log(x)-log(Topt))/Troh)^2)
      if(proc=="mortal") frm<- "m(T) = 1-rm.exp((-0.5).(-(log(x)-log(Topt))/Troh)^2)" else frm<- "f(T) = 1-rm.exp((-0.5).(-(log(x)-log(Topt))/Troh)^2)"
    }
    
    
    if(modelm==26){
      f <- function(x,Topt,B,H)
      {
        expr <- expression(1 - 1/(exp((1+exp(-(x-Topt)/B))*(1+exp(-(Topt-x)/B))*H)))
        eval(expr)
      }
      j <- function(x,Topt,B,H)
      {
        expr <- expression(1 - 1/(exp((1+exp(-(x-Topt)/B))*(1+exp(-(Topt-x)/B))*H)))
        c(eval(D(expr, "Topt")),eval(D(expr, "B")), eval(D(expr, "H")))
      }
      g<-y ~ 1 - 1/(exp((1+exp(-(x-Topt)/B))*(1+exp(-(Topt-x)/B))*H))
      if(proc=="mortal") frm<- "m(T) = 1 - 1/(exp((1+exp(-(x-Topt)/B)).(1+exp(-(Topt-x)/B)).H))" else frm<- "f(T) = 1 - 1/(exp((1+exp(-(x-Topt)/B)).(1+exp((Topt-x)/B)).H))"
    }
    
    
    if(modelm==27){
      f <- function(x,Tl,Th,B,H)
      {
        expr <- expression(1 - 1/(exp((1+exp(-(x-Tl)/B))*(1+exp(-(Th-x)/B))*H)))
        eval(expr)
      }
      j <- function(x,Tl,Th,B,H)
      {
        expr <- expression(1 - 1/(exp((1+exp(-(x-Tl)/B))*(1+exp(-(Th-x)/B))*H)))
        c(eval(D(expr, "Tl")),eval(D(expr, "Th")), eval(D(expr, "B")), eval(D(expr, "H")))
      }
      g<-y ~ 1 - 1/(exp((1+exp(-(x-Tl)/B))*(1+exp(-(Th-x)/B))*H))
      if(proc=="mortal") frm<- "m(T) = 1 - 1/(exp((1+exp(-(x-Tl)/B)).(1+exp(-(Th-x)/B)).H))" else frm<- "f(T) = 1 - 1/(exp((1+exp(-(x-Tl)/B)).(1+exp(-(Th-x)/B)).H))"
    }
    
    
    if(modelm==28){
      f <- function(x,Topt,Bl,Bh,H)
      {
        expr <- expression(1 - 1/(exp((1+exp(-(x-Topt)/Bl))*(1+exp(-(Topt-x)/Bh))*H)))
        eval(expr)
      }
      j <- function(x,Topt,Bl,Bh,H)
      {
        expr <- expression(1 - 1/(exp((1+exp(-(x-Topt)/Bl))*(1+exp(-(Topt-x)/Bh))*H)))
        c(eval(D(expr, "Topt")),eval(D(expr, "Bl")),eval(D(expr, "Bh")), eval(D(expr, "H")))
      }
      g<-y ~ 1 - 1/(exp((1+exp(-(x-Topt)/Bl))*(1+exp(-(Topt-x)/Bh))*H))
      if(proc=="mortal") frm<- "m(T) = 1 - 1/(exp((1+exp(-(x-Topt)/Bl)).(1+exp(-(Topt-x)/Bh)).H))" else frm<- "f(T) = 1 - 1/(exp((1+exp(-(x-Topt)/Bl)).(1+exp(-(Topt-x)/Bh)).H))"
    }
    
    if(modelm==29){
      f <- function(x,Tl,Th,Bl,Bh,H)
      {
        expr <- expression(1 - 1/(exp((1+exp(-(x-Tl)/Bl))*(1+exp(-(Th-x)/Bh))*H)))
        eval(expr)
      }
      j <- function(x,Tl,Th,Bl,Bh,H)
      {
        expr <- expression(1 - 1/(exp((1+exp(-(x-Tl)/Bl))*(1+exp(-(Th-x)/Bh))*H)))
        c(eval(D(expr, "Tl")),eval(D(expr, "Th")),eval(D(expr, "Bl")),eval(D(expr, "Bh")), eval(D(expr, "H")))
      }
      g<-y ~ 1 - 1/(exp((1+exp(-(x-Tl)/Bl))*(1+exp(-(Th-x)/Bh))*H))
      if(proc=="mortal") frm<- "m(T) = 1 - 1/(exp((1+exp(-(x-Tl)/Bl)).(1+exp(-(Th-x)/Bh)).H))" else frm<- "f(T) = 1 - 1/(exp((1+exp(-(x-Tl)/Bl)).(1+exp(-(Th-x)/Bh)).H))"
    }
    
    
    
    if(modelm==30){
      f <- function(x,Topt,B,H)
      {
        expr <- expression(1 - H/(exp(1+exp(-(x-Topt)/B))*(1+exp(-(Topt-x)/B))))
        eval(expr)
      }
      j <- function(x,Topt,B,H)
      {
        expr <- expression(1 - H/(exp(1+exp(-(x-Topt)/B))*(1+exp(-(Topt-x)/B))))
        c(eval(D(expr, "Topt")),eval(D(expr, "B")), eval(D(expr, "H")))
      }
      g<-y ~ 1 - H/(exp(1+exp(-(x-Topt)/B))*(1+exp(-(Topt-x)/B)))
      if(proc=="mortal") frm<- "m(T) = 1 - H/(exp(1+exp(-(x-Topt)/B)).(1+exp(-(Topt-x)/B)))" else frm<- "f(T) = 1 - H/(exp(1+exp(-(x-Topt)/B)).(1+exp(-(Topt-x)/B)))"
    }
    
    
    if(modelm==31){
      f <- function(x,Tl,Th,B,H)
      {
        expr <- expression(1 - H/(exp(1+exp(-(x-Tl)/B))*(1+exp(-(Th-x)/B))))
        eval(expr)
      }
      j <- function(x,Tl,Th,B,H)
      {
        expr <- expression(1 - H/(exp(1+exp(-(x-Tl)/B))*(1+exp(-(Th-x)/B))))
        c(eval(D(expr, "Tl")),eval(D(expr, "Th")), eval(D(expr, "B")), eval(D(expr, "H")))
      }
      g<-y ~ 1 - H/(exp(1+exp(-(x-Tl)/B))*(1+exp(-(Th-x)/B)))
      if(proc=="mortal") frm<- "m(T) = 1 - H/(exp(1+exp(-(x-Tl)/B))*(1+exp(-(Th-x)/B)))" else frm<- "f(T) = 1 - H/(exp(1+exp(-(x-Tl)/B))*(1+exp(-(Th-x)/B)))"
    }
    
    if(modelm==32){
      f <- function(x,Topt,Bl,Bh,H)
      {
        expr <- expression(1 - H/(exp(1+exp(-(x-Topt)/Bl))*(1+exp(-(Topt-x)/Bh))))
        eval(expr)
      }
      j <- function(x,Topt,Bl,Bh,H)
      {
        expr <- expression(1 - H/(exp(1+exp(-(x-Topt)/Bl))*(1+exp(-(Topt-x)/Bh))))
        c(eval(D(expr, "Topt")),eval(D(expr, "Bl")),eval(D(expr, "Bh")), eval(D(expr, "H")))
      }
      g<-y ~ 1 - H/(exp(1+exp(-(x-Topt)/Bl))*(1+exp(-(Topt-x)/Bh)))
      if(proc=="mortal") frm<- "m(T) = 1 - H/(exp(1+exp(-(x-Topt)/Bl)).(1+exp(-(Topt-x)/Bh)))" else frm<- "f(T) = 1 - H/(exp(1+exp(-(x-Topt)/Bl)).(1+exp(-(Topt-x)/Bh))"
    }
    
    
    if(modelm==33){
      f <- function(x,Tl,Th,Bl,Bh,H)
      {
        expr <- expression(1 - H/(exp(1+exp(-(x-Tl)/Bl))*(1+exp(-(Th-x)/Bh))))
        eval(expr)
      }
      j <- function(x,Tl,Th,Bl,Bh,H)
      {
        expr <- expression(1 - H/(exp(1+exp(-(x-Tl)/Bl))*(1+exp(-(Th-x)/Bh))))
        c(eval(D(expr, "Tl")),eval(D(expr, "Th")),eval(D(expr, "Bl")),eval(D(expr, "Bh")), eval(D(expr, "H")))
      }
      g<-y ~ 1 - H/(exp(1+exp(-(x-Tl)/Bl))*(1+exp(-(Th-x)/Bh)))
      if(proc=="mortal") frm<- "m(T) = 1 - H/(exp(1+exp(-(x-Tl)/Bl)).(1+exp(-(Th-x)/Bh)))" else frm<- "f(T) = 1 - H/(exp(1+exp(-(x-Tl)/Bl)).(1+exp(-(Th-x)/Bh)))"
    }
    
    
    
    if(modelm==34){
      f <- function(x,Hl,Tl,Hh,Th,Bm)
      {
        expr <- expression(Bm + ((1-1/(1+exp((Hl/1.987)*((1/Tl)-(1/(x+273.15))))+exp((Hh/1.987)*((1/Th)-(1/(x+273.15))))))*(1-Bm))) #DM: se agrego 2 parentesis "Bm + (...)*"
        eval(expr)
      }
      j <- function(x,Hl,Tl,Hh,Th,Bm)
      {
        expr <- expression((Bm + (1-1/(1+exp((Hl/1.987)*((1/Tl)-(1/(x+273.15))))+exp((Hh/1.987)*((1/Th)-(1/(x+273.15))))))*(1-Bm))) #DM: se agrego 2 parentesis "Bm + (...)*"
        c(eval(D(expr, "Hl")),eval(D(expr, "Tl")),eval(D(expr, "Hh")),eval(D(expr, "Th")), eval(D(expr, "Bm")))
      }
      g<-y ~ Bm + ((1-1/(1+exp((Hl/1.987)*((1/Tl)-(1/(x+273.15))))+exp((Hh/1.987)*((1/Th)-(1/(x+273.15))))))*(1-Bm)) #DM: se agrego 2 parentesis "Bm + (...)*"
      if(proc=="mortal") frm<- "m(T) = Bm + ((1-1/(1+exp((Hl/1.987).((1/Tl)-(1/x)))+exp((Hh/1.987)*((1/Th)-(1/x))))).(1-Bm))" else frm<- "f(T) = Bm + ((1-1/(1+exp((Hl/1.987).((1/Tl)-(1/x)))+exp((Hh/1.987)*((1/Th)-(1/x))))).(1-Bm))" #DM: se agrego 2 parentesis "Bm + (...)*"
    }
    
    
    
    if(modelm==35){
      f <- function(x,a1,b1,a2,b2)
      {
        expr <- expression((1 - exp(-(exp(a1+b1*x)))) + (1 - exp(-(exp(a2+b2*x)))))
        eval(expr)
      }
      j <- function(x,a,b)
      {
        expr <- expression((1 - exp(-(exp(a1+b1*x)))) + (1 - exp(-(exp(a2+b2*x)))))
        c(eval(D(expr, "a1")),eval(D(expr, "b1")),eval(D(expr, "a2")),eval(D(expr, "b2")))
      }
      g<-y ~ (1 - exp(-(exp(a1+b1*x)))) + (1 - exp(-(exp(a2+b2*x))))
      if(proc=="mortal") frm<- "m(T) = (1 - exp(-(exp(a1+b1.T)))) + (1 - exp(-(exp(a2+b2.T))))" else frm<- "f(T) = (1 - exp(-(exp(a1+b1.T)))) + (1 - exp(-(exp(a2+b2.T))))"
    }
    
    
    if(modelm==36){
      f <- function(x,w)
      {
        expr <- expression((w-x)^(-1))
        eval(expr)
      }
      j <- function(x,w)
      {
        expr <- expression((w-x)^(-1))
        c(eval(D(expr, "w")))
      }
      g<-y ~ (w-x)^(-1)
      if(proc=="mortal") frm<- "m(T) = (w-x)^(-1)" else frm<- "f(T) = (w-x)^(-1)"
    }
    
    
    if(modelm==37){
      f <- function(x,a1,b1,a2,b2)
      {
        expr <- expression(a1*exp(b1*x) + a2*exp(b2*x))
        eval(expr)
      }
      j <- function(x,a1,b1,a2,b2)
      {
        expr <- expression(a1*exp(b1*x) + a2*exp(b2*x))
        c(eval(D(expr, "a1")),eval(D(expr, "b1")),eval(D(expr, "a2")),eval(D(expr, "b2")))
      }
      g<-y ~ a1*exp(b1*x) + a2*exp(b2*x)
      if(proc=="mortal") frm<- "m(T) = a1.exp(b1.x) + a2.exp(b2.x)" else frm<- "f(T) = a1.exp(b1.x) + a2.exp(b2.x)"
    }
    
    if(modelm==38){
      f <- function(x,a1,b1,a2,b2,c1)
      {
        expr <- expression(a1*exp(b1*x) + a2*exp(b2*x) + c1)
        eval(expr)
      }
      j <- function(x,a1,b1,a2,b2,c1)
      {
        expr <- expression(a1*exp(b1*x) + a2*exp(b2*x) + c1)
        c(eval(D(expr, "a1")),eval(D(expr, "b1")),eval(D(expr, "a2")),eval(D(expr, "b2")),eval(D(expr, "c1")))
      }
      g<-y ~ a1*exp(b1*x)+c1
      if(proc=="mortal") frm<- "m(T) = a1.exp(b1.x) + a2.exp(b2.x) + c1" else frm<- "f(T) = a1.exp(b1.x) + a2.exp(b2.x) + c1"
    }
    
    
    if(modelm==39){
      f <- function(x,a,b,nn)
      {
        expr <- expression(a*(abs(x-b))^nn)
        eval(expr)
      }
      j <- function(x,a,b,nn)
      {
        expr <- expression(a*(abs(x-b))^nn)
        c(eval(D(expr, "a")),eval(D(expr, "b")),eval(D(expr, "nn")))
      }
      g<-y ~ a*(abs(x-b))^nn
      if(proc=="mortal") frm<- "m(T) = a.(abs(x-b))^nn" else frm<- "f(T) = a.(abs(x-b))^nn"
    }
    
    
    
    if(modelm==40){
      f <- function(x,a,To,Tl,d)
      {
        expr <- expression(a*x*(x-To)*(Tl-x)^(1/d))
        eval(expr)
      }
      j <- function(x,a,To,Tl,d) #DM: Se corrigi?, antes tenia los parametros del modelo 39
      {
        expr <- expression(a*x*(x-To)*(Tl-x)^(1/d))
        c(eval(D(expr, "a")),eval(D(expr, "To")),eval(D(expr, "Tl")),eval(D(expr, "d")))
      }
      g<-y ~ a*x*(x-To)*(Tl-x)^(1/d)
      if(proc=="mortal") frm<- "m(T) = a.x.(x-To).(Tl-x)^(1/d)" else frm<- "f(T) = a.x.(x-To).(Tl-x)^(1/d)"
    }
    
    
    if(modelm==41){
      f <- function(x,a,To,Tl,d)
      {
        expr <- expression(exp(a*x*(x-To)*(Tl-x)^(1/d)))
        eval(expr)
      }
      j <- function(x,a,To,Tl,d) #DM: Se corrigi?, antes tenia los parametros del modelo 39
      {
        expr <- expression(exp(a*x*(x-To)*(Tl-x)^(1/d)))
        c(eval(D(expr, "a")),eval(D(expr, "To")),eval(D(expr, "Tl")),eval(D(expr, "d")))
      }
      g<-y ~ exp(a*x*(x-To)*(Tl-x)^(1/d))
      if(proc=="mortal") frm<- "m(T) = exp(a.x.(x-To).(Tl-x)^(1/d))" else frm<- "f(T) = exp(a.x.(x-To).(Tl-x)^(1/d))"
    }
    
    
    if(modelm==42){
      f <- function(x,a,Tmax, Tmin,n,m)
      {
        expr <- expression(a*((x-Tmin)^n)*(Tmax-x)^m)
        eval(expr)
      }
      j <- function(x,a,Tmax, Tmin,n,m)
      {
        expr <- expression(a*((x-Tmin)^n)*(Tmax-x)^m)
        c(eval(D(expr, "a")),eval(D(expr, "Tmax")),eval(D(expr, "Tmin")),eval(D(expr, "n")),eval(D(expr, "m")))
      }
      g<-y ~ a*((x-Tmin)^n)*(Tmax-x)^m
      if(proc=="mortal") frm<- "m(T) = a.((x-Tmin)^n).(Tmax-x)^m" else frm<- "f(T) = a.((x-Tmin)^n).(Tmax-x)^m"
    }
    
    if(modelm==43){
      f <- function(x,Dmin,k,Tp,lamb)
      {
        expr <- expression(1/((Dmin/2) * (exp(k*(x-Tp)) + exp(-(x-Tp)*lamb))))
        eval(expr)
      }
      j <- function(x,Dmin,k,Tp,lamb)
      {
        expr <- expression(1/((Dmin/2) * (exp(k*(x-Tp)) + exp(-(x-Tp)*lamb))))
        c(eval(D(expr, "Dmin")),eval(D(expr, "k")),eval(D(expr, "Tp")),eval(D(expr, "lamb")))
      }
      g<-y ~ 1/((Dmin/2) * (exp(k*(x-Tp)) + exp(-(x-Tp)*lamb)))
      if(proc=="mortal") frm<- "m(T) = 1/((Dmin/2) . (exp(k.(x-Tp)) + exp(-(x-Tp).lamb)))" else frm<- "f(T) = 1/((Dmin/2) . (exp(k.(x-Tp)) + exp(-(x-Tp).lamb)))"
    }
    
    
    if(modelm==44){
      f <- function(x,a,Tl, Th,B)
      {
        expr <- expression(a*(1-exp(-(x-Tl)/B))*(1-exp(-(Th-x)/B)))
        eval(expr)
      }
      j <- function(x,a,Tl, Th,B)
      {
        expr <- expression(a*(1-exp(-(x-Tl)/B))*(1-exp(-(Th-x)/B)))
        c(eval(D(expr, "a")),eval(D(expr, "Tl")),eval(D(expr, "Th")),eval(D(expr, "B")))
      }
      g<-y ~ a*(1-exp(-(x-Tl)/B))*(1-exp(-(Th-x)/B))
      if(proc=="mortal") frm<- "m(T) = a.(1-exp(-(x-Tl)/B)).(1-exp(-(Th-x)/B))" else frm<- "f(T) = a.(1-exp(-(x-Tl)/B)).(1-exp(-(Th-x)/B))"
    }
    
    if(modelm==45){
      f <- function(x,a,Tl, Th,Bl,Bh)
      {
        expr <- expression(exp(a*(1-exp(-(x-Tl)/Bl))*(1-exp(-(Th-x)/Bh))))
        eval(expr)
      }
      j <- function(x,a,Tl, Th,Bl,Bh)
      {
        expr <- expression(exp(a*(1-exp(-(x-Tl)/Bl))*(1-exp(-(Th-x)/Bh))))
        c(eval(D(expr, "a")),eval(D(expr, "Tl")),eval(D(expr, "Th")),eval(D(expr, "Bl")),eval(D(expr, "Bh")))
      }
      g<-y ~ exp(a*(1-exp(-(x-Tl)/Bl))*(1-exp(-(Th-x)/Bh)))
      if(proc=="mortal") frm<- "m(T) = exp(a.(1-exp(-(x-Tl)/Bl)).(1-exp(-(Th-x)/Bh)))" else frm<- "f(T) = exp(a.(1-exp(-(x-Tl)/Bl)).(1-exp(-(Th-x)/Bh)))"
    }
    
    
    fcn     <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), ini)))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    estimate<-as.data.frame(out$par)
    stdmor<-diag(ginv(out$hessian))
    return(list(estimados=estimate,f=g,modelo=frm,stdmort=stdmor,ecua=frm))
  }
}