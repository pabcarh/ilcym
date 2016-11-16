#' shape
#'
#' no details
#'
#' @param model a character, dicotomic model ('logit','probit', 'cloglog')
#' @param datashap data frame, estimated rates
#' @param datao data frame
#' @param ini numeric vector of the initial values for the modeling
#' @param coefi numeric vector, linear regression parameters
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
shape<-function(model,datashap,datao,ini,coefi){
  
  ####  Shape&DeMichele sin p y To
  if(model==1) {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]+273.15
      y<-datashap[,2]
    }else{
      x<-datao[,1]+273.15
      y<-datao[,2]
    }
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,Ha, Hl,Tl,Hh,Th){
      expr <- expression(((coefi[1]+coefi[2]*((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))) * (x/(((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th))))) * exp((Ha/1.987) * ((1/((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))) - (1/x))))/
                           (1 + exp((Hl/1.987) * ((1/Tl) - (1/x))) + exp((Hh/1.987) * ((1/Th) - (1/x)))))
      eval(expr)
    }
    j <- function(x, Ha, Hl,Tl,Hh,Th){
      expr <- expression(((coefi[1]+coefi[2]*((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th))))* (x/(((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th))))) * exp((Ha/1.987) * ((1/((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))) - (1/x))))/
                           (1 + exp((Hl/1.987) * ((1/Tl) - (1/x))) + exp((Hh/1.987) * ((1/Th) - (1/x)))))
      c(eval(D(expr, "Ha"  )), eval(D(expr, "Hl" )),
        eval(D(expr, "Tl" )), eval(D(expr, "Hh" )), eval(D(expr, "Th" )))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,p=(coefi[1]+coefi[2]*((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))),stderro=stderro))
  }
  
  ####  Shape&DeMichele con To
  if(model==2) {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]+273.15
      y<-datashap[,2]
    }else{
      x<-datao[,1]+273.15
      y<-datao[,2]
    }
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,To, Ha, Hl,Tl,Hh,Th){
      expr <- expression(((coefi[1]+coefi[2]*To) * (x/(To)) * exp((Ha/1.987) * ((1/To) - (1/x))))/
                           (1 + exp((Hl/1.987) * ((1/Tl) - (1/x))) + exp((Hh/1.987) * ((1/Th) - (1/x)))))
      eval(expr)
    }
    j <- function(x, To, Ha, Hl,Tl,Hh,Th){
      expr <- expression(((coefi[1]+coefi[2]*To)* (x/(To)) * exp((Ha/1.987) * ((1/To) - (1/x))))/
                           (1 + exp((Hl/1.987) * ((1/Tl) - (1/x))) + exp((Hh/1.987) * ((1/Th) - (1/x)))))
      c(eval(D(expr, "To")), eval(D(expr, "Ha"  )), eval(D(expr, "Hl" )),
        eval(D(expr, "Tl" )), eval(D(expr, "Hh" )), eval(D(expr, "Th" )))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,p=(coefi[1]+coefi[2]*as.numeric(estimate[1])),stderro=stderro))
  }
  
  ####  Shape&DeMichele con To y p
  if(model==3) {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]+273.15
      y<-datashap[,2]
    }else{
      x<-datao[,1]+273.15
      y<-datao[,2]
    }
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,p,To, Ha, Hl,Tl,Hh,Th){
      expr <- expression((p * (x/(To)) * exp((Ha/1.987) * ((1/To) - (1/x))))/
                           (1 + exp((Hl/1.987) * ((1/Tl) - (1/x))) + exp((Hh/1.987) * ((1/Th) - (1/x)))))
      eval(expr)
    }
    j <- function(x, p,To, Ha, Hl,Tl,Hh,Th){
      expr <- expression((p* (x/(To)) * exp((Ha/1.987) * ((1/To) - (1/x))))/
                           (1 + exp((Hl/1.987) * ((1/Tl) - (1/x))) + exp((Hh/1.987) * ((1/Th) - (1/x)))))
      c(eval(D(expr, "p")),eval(D(expr, "To")), eval(D(expr, "Ha"  )), eval(D(expr, "Hl" )),
        eval(D(expr, "Tl" )), eval(D(expr, "Hh" )), eval(D(expr, "Th" )))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,p=p,stderro=stderro))
  }
  
  ####  Shape&DeMichele solo con To, p y Ha
  if(model==4) {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]+273.15
      y<-datashap[,2]
    }else{
      x<-datao[,1]+273.15
      y<-datao[,2]
    }
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,p,To, Ha){
      expr <- expression((p * (x/(To)) * exp((Ha/1.987) * ((1/To) - (1/x)))))
      eval(expr)
    }
    j <- function(x, p,To, Ha){
      expr <- expression((p* (x/(To)) * exp((Ha/1.987) * ((1/To) - (1/x)))))
      c(eval(D(expr, "p")),eval(D(expr, "To")), eval(D(expr, "Ha"  )))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,p=p,stderro=stderro))
  }
  
  ####  Shape&DeMichele con To, p, Ha, Hl, Tl
  if(model==5) {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]+273.15
      y<-datashap[,2]
    }else{
      x<-datao[,1]+273.15
      y<-datao[,2]
    }
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,p,To, Ha, Hl,Tl){
      expr <- expression((p * (x/(To)) * exp((Ha/1.987) * ((1/To) - (1/x))))/
                           (1 + exp((Hl/1.987) * ((1/Tl) - (1/x)))))
      eval(expr)
    }
    j <- function(x, p,To, Ha, Hl,Tl){
      expr <- expression((p* (x/(To)) * exp((Ha/1.987) * ((1/To) - (1/x))))/
                           (1 + exp((Hl/1.987) * ((1/Tl) - (1/x)))))
      c(eval(D(expr, "p")),eval(D(expr, "To")), eval(D(expr, "Ha"  )), eval(D(expr, "Hl" )),
        eval(D(expr, "Tl" )))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,p=p,stderro=stderro))
  }
  
  ####  Shape&DeMichele con To, p, Ha, Hh, Th
  if(model==6) {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]+273.15
      y<-datashap[,2]
    }else{
      x<-datao[,1]+273.15
      y<-datao[,2]
    }
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,p,To, Ha,Hh,Th){
      expr <- expression((p * (x/(To)) * exp((Ha/1.987) * ((1/To) - (1/x))))/
                           (1 + exp((Hh/1.987) * ((1/Th) - (1/x)))))
      eval(expr)
    }
    j <- function(x, p,To, Ha, Hh,Th){
      expr <- expression((p* (x/(To)) * exp((Ha/1.987) * ((1/To) - (1/x))))/
                           (1 + exp((Hh/1.987) * ((1/Th) - (1/x)))))
      c(eval(D(expr, "p")),eval(D(expr, "To")), eval(D(expr, "Ha"  )),
        eval(D(expr, "Hh" )), eval(D(expr, "Th" )))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,p=p,stderro=stderro))
  }
  
  ####  Shape&DeMichele con To, Ha
  
  
  if(model==7) {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]+273.15
      y<-datashap[,2]
    }else{
      x<-datao[,1]+273.15
      y<-datao[,2]
    }
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,To, Ha){
      expr <- expression(((coefi[1]+coefi[2]*To) * (x/(To)) * exp((Ha/1.987) * ((1/To) - (1/x)))))
      eval(expr)
    }
    j <- function(x, To, Ha){
      expr <- expression(((coefi[1]+coefi[2]*To)* (x/(To)) * exp((Ha/1.987) * ((1/To) - (1/x)))))
      c(eval(D(expr, "To")), eval(D(expr, "Ha"  )))
    }
    
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,p=(coefi[1]+coefi[2]*To),stderro=stderro)) ## este valor de To no es el estimado por el modelo sino del valor inicial (por eso el valor de p es el mismo que da la funcion "prueba")
  }
  
  ####  Shape&DeMichele con To, Ha, Hl, Tl
  if(model==8) {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]+273.15
      y<-datashap[,2]
    }else{
      x<-datao[,1]+273.15
      y<-datao[,2]
    }
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    
    f <- function(x,To, Ha, Hl,Tl){
      expr <- expression(((coefi[1]+coefi[2]*To) * (x/(To)) * exp((Ha/1.987) * ((1/To) - (1/x))))/
                           (1 + exp((Hl/1.987) * ((1/Tl) - (1/x)))))
      eval(expr)
    }
    j <- function(x,To, Ha, Hl,Tl){
      expr <- expression(((coefi[1]+coefi[2]*To) * (x/(To)) * exp((Ha/1.987) * ((1/To) - (1/x))))/
                           (1 + exp((Hl/1.987) * ((1/Tl) - (1/x)))))
      c(eval(D(expr, "To")), eval(D(expr, "Ha"  )), eval(D(expr, "Hl" )),
        eval(D(expr, "Tl" )))
    }
    
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,p=(coefi[1]+coefi[2]*To),stderro=stderro))
  }
  
  ####  Shape&DeMichele con To, Ha, Hh, Th
  if(model==9) {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]+273.15
      y<-datashap[,2]
    }else{
      x<-datao[,1]+273.15
      y<-datao[,2]
    }
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    
    f <- function(x,To, Ha,Hh,Th){
      expr <- expression(((coefi[1]+coefi[2]*To) * (x/(To)) * exp((Ha/1.987) * ((1/To) - (1/x))))/
                           (1 + exp((Hh/1.987) * ((1/Th) - (1/x)))))
      eval(expr)
    }
    j <- function(x,To, Ha, Hh,Th){
      expr <- expression(((coefi[1]+coefi[2]*To) * (x/(To)) * exp((Ha/1.987) * ((1/To) - (1/x))))/
                           (1 + exp((Hh/1.987) * ((1/Th) - (1/x)))))
      c(eval(D(expr, "To")), eval(D(expr, "Ha"  )),
        eval(D(expr, "Hh" )), eval(D(expr, "Th" )))
    }
    
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,p=(coefi[1]+coefi[2]*To),stderro=stderro))
  }
  
  if(model==10){ ####  Shape&DeMichele con To, Hh y Th, constantes
    if(length(datashap[,1])>(length(ini)+1))
    {
      x<-datashap[,1]+273.15
      y<-datashap[,2]
    }else
    {
      x<-datao[,1]+273.15
      y<-datao[,2]
    }
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind))
    {
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,p,Ha, Hl,Tl)
    {
      expr <- expression((p * (x/(298.16)) * exp((Ha/1.987) * ((1/298.16) - (1/x))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/x))) + exp((1/1.987) * ((1/1) - (1/x)))))
      eval(expr)
    }
    j <- function(x,p,Ha, Hl,Tl)
    {
      expr <- expression((p * (x/(298.16)) * exp((Ha/1.987) * ((1/298.16) - (1/x))))/
                           (1 + exp((Hl/1.987) * ((1/Tl) - (1/x))) + exp((1/1.987) * ((1/1) - (1/x)))))
      c(eval(D(expr, "p")), eval(D(expr, "Ha"  )), eval(D(expr, "Hl" )),
        eval(D(expr, "Tl" )))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,p=p,stderro=stderro))
  }
  
  
  ################## nuevos modelos de sharpe de michelle
  
  
  if(model==11) {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]+273.15
      y<-datashap[,2]
    }else{
      x<-datao[,1]+273.15
      y<-datao[,2]
    }
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,p,Ha, Hl,Tl,Hh,Th){
      
      expr <- expression((p * (x/298.16) * exp((Ha/1.987) * ((1/298.16) - (1/x))))/
                           (1 + exp((Hl/1.987) * ((1/Tl) - (1/x))) + exp((Hh/1.987) * ((1/Th) - (1/x)))))
      eval(expr)
    }
    j <- function(x,p,Ha, Hl,Tl,Hh,Th){
      expr <- expression((p * (x/298.16) * exp((Ha/1.987) * ((1/298.16) - (1/x))))/
                           (1 + exp((Hl/1.987) * ((1/Tl) - (1/x))) + exp((Hh/1.987) * ((1/Th) - (1/x)))))
      c(eval(D(expr, "p")),eval(D(expr, "Ha"  )), eval(D(expr, "Hl" )),
        eval(D(expr, "Tl" )), eval(D(expr, "Hh" )), eval(D(expr, "Th" )))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,p=p,stderro=stderro))
  }
  
  
  if(model==12) {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]+273.15
      y<-datashap[,2]
    }else{
      x<-datao[,1]+273.15
      y<-datao[,2]
    }
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x, Ha, Hl,Tl,Hh,Th){
      expr <- expression(((coefi[1]+coefi[2]*298.16) * (x/298.16) * exp((Ha/1.987) * ((1/298.16) - (1/x))))/
                           (1 + exp((Hl/1.987) * ((1/Tl) - (1/x))) + exp((Hh/1.987) * ((1/Th) - (1/x)))))
      eval(expr)
    }
    j <- function(x, Ha, Hl,Tl,Hh,Th){
      expr <- expression(((coefi[1]+coefi[2]*298.16) * (x/298.16) * exp((Ha/1.987) * ((1/298.16) - (1/x))))/
                           (1 + exp((Hl/1.987) * ((1/Tl) - (1/x))) + exp((Hh/1.987) * ((1/Th) - (1/x)))))
      c(eval(D(expr, "Ha"  )), eval(D(expr, "Hl" )),
        eval(D(expr, "Tl" )), eval(D(expr, "Hh" )), eval(D(expr, "Th" )))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,p=(coefi[1]+coefi[2]*298.16),stderro=stderro))
  }
  
  
  
  if(model==13) {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]+273.15
      y<-datashap[,2]
    }else{
      x<-datao[,1]+273.15
      y<-datao[,2]
    }
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,p,Ha, Hl,Tl){
      
      expr <- expression((p * (x/298.16) * exp((Ha/1.987) * ((1/298.16) - (1/x))))/
                           (1 + exp((Hl/1.987) * ((1/Tl) - (1/x)))))
      eval(expr)
    }
    j <- function(x,p,Ha, Hl,Tl){
      expr <- expression((p * (x/298.16) * exp((Ha/1.987) * ((1/298.16) - (1/x))))/
                           (1 + exp((Hl/1.987) * ((1/Tl) - (1/x)))))
      c(eval(D(expr, "p")),eval(D(expr, "Ha"  )), eval(D(expr, "Hl" )),
        eval(D(expr, "Tl" )))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,p=p,stderro=stderro))
  }
  
  if(model==14) {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]+273.15
      y<-datashap[,2]
    }else{
      x<-datao[,1]+273.15
      y<-datao[,2]
    }
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,p,Ha, Hh,Th){
      
      expr <- expression((p * (x/298.16) * exp((Ha/1.987) * ((1/298.16) - (1/x))))/
                           (1 + exp((Hh/1.987) * ((1/Th) - (1/x)))))
      eval(expr)
    }
    j <- function(x,p,Ha, Hh,Th){
      expr <- expression((p * (x/298.16) * exp((Ha/1.987) * ((1/298.16) - (1/x))))/
                           (1 + exp((Hh/1.987) * ((1/Th) - (1/x)))))
      c(eval(D(expr, "p")),eval(D(expr, "Ha"  )), eval(D(expr, "Hh" )), eval(D(expr, "Th" )))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))                                    
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,p=p,stderro=stderro))
  }
  
  
  ################## fin nuevos modelos de sharpe de michelle
  
  
  
  
  
  ####  deva 1
  if(model==15)   {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,Tmin,b){
      expr <- expression(b*(x-Tmin))
      eval(expr)
    }
    j <- function(x,Tmin,b){
      expr <- expression(b*(x-Tmin))
      c(eval(D(expr, "Tmin" )), eval(D(expr, "b" )))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  ####  deva no lineal o Higgis
  if(model==16){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,b1,b2,b3,b4,b5){
      expr <- expression(b1*10^(-((((x-b3)/(b3-b2)-(1/(1+0.28*b4+0.72*log(1+b4))))+exp(b4*((x-b3)/(b3-b2)-(1/(1+0.28*b4+0.72*log(1+b4))))))/((1+b4)/(1+1.5*b4+0.39*b4^2)))^2)*(1-b5+b5*((((x-b3)/(b3-b2)- (1/(1+0.28*b4+0.72*log(1+b4))))+exp(b4*((x-b3)/(b3-b2)- (1/(1+0.28*b4+0.72*log(1+b4))))))/((1+b4)/(1+1.5*b4+0.39*b4^2)))^2)) #DM: Se agreg? un menos "-v^2"
      eval(expr)
    }
    j <- function(x,b1,b2,b3,b4,b5){
      expr <- expression(b1*10^(-((((x-b3)/(b3-b2)-(1/(1+0.28*b4+0.72*log(1+b4))))+exp(b4*((x-b3)/(b3-b2)-(1/(1+0.28*b4+0.72*log(1+b4))))))/((1+b4)/(1+1.5*b4+0.39*b4^2)))^2)*(1-b5+b5*((((x-b3)/(b3-b2)- (1/(1+0.28*b4+0.72*log(1+b4))))+exp(b4*((x-b3)/(b3-b2)- (1/(1+0.28*b4+0.72*log(1+b4))))))/((1+b4)/(1+1.5*b4+0.39*b4^2)))^2)) #DM: Se agreg? un menos "-v^2"
      c(eval(D(expr, "b1")), eval(D(expr, "b2")), eval(D(expr, "b3"  )), eval(D(expr, "b4" )),
        eval(D(expr, "b5" )))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  ####  Logan 1
  if(model==17) {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,Y,Tmax, p, v){
      expr <- expression(Y*(exp(p*x)-exp(p*Tmax-(Tmax-x)/v)))
      eval(expr)
    }
    j <- function(x,Y,Tmax, p, v){
      expr <- expression(Y*(exp(p*x)-exp(p*Tmax-(Tmax-x)/v)))
      c(eval(D(expr, "Y")),eval(D(expr, "Tmax")), eval(D(expr, "p")), eval(D(expr, "v")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  ###  Logan 2
  if(model==18) {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,alfa,k,Tmax, p, v){
      expr <- expression(alfa*((1/(1+k*exp(-p*x)))-exp(-(Tmax-x)/v)))
      eval(expr)
    }
    j <- function(x,alfa,k,Tmax, p, v){
      expr <- expression(2*((1/(1+3*exp(-p*x)))-exp(-(Tmax-x)/v)))
      c(eval(D(expr, "alfa")),eval(D(expr, "k")),eval(D(expr, "Tmax")), eval(D(expr, "p")), eval(D(expr, "v")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  ######  Briere 1
  if(model==19) {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,aa,To,Tmax){
      expr <- expression(aa*x*(x-To)*(Tmax-x)^0.5)
      eval(expr)
    }
    j <- function(x,aa,To,Tmax){
      expr <- expression(aa*x*(x-To)*(Tmax-x)^0.5)
      c(eval(D(expr, "aa")), eval(D(expr, "To")), eval(D(expr, "Tmax")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  ##  Briere 2
  if(model==20){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,aa,To,Tmax,d){
      expr <- expression(aa*x*(x-To)*(Tmax-x)^d)
      eval(expr)
    }
    j <- function(x,aa,To,Tmax,d){
      expr <- expression(aa*x*(x-To)*(Tmax-x)^d)
      c(eval(D(expr, "d")),eval(D(expr, "To")),eval(D(expr, "aa")),eval(D(expr, "Tmax")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  ####  Stinner
  if(model==21)    {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]
      datt<-datashap}else{
        x<-datao[,1]
        y<-datao[,2]
        datt<-datao}
    maximos<-subset(datt,datt[,2]==max(datt[,2]))
    xp<-c(1)
    for(i in 1:length(datt[,1])) ifelse(x[i]<=maximos[,1],xp[i]<-x[i],xp[i]<-2*(maximos[,1])-x[i])
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(xp,Rmax,Topc,k1,k2){
      expr <- expression((Rmax*(1+exp(k1+k2*(Topc))))/(1+exp(k1+k2*xp)))
      eval(expr)
    }
    j <- function(xp,Rmax,Topc,k1,k2){
      expr <- expression((Rmax*(1+exp(k1+k2*(Topc))))/(1+exp(k1+k2*xp)))
      c(eval(D(expr, "Rmax")), eval(D(expr, "k1")), eval(D(expr, "k2")), eval(D(expr, "Topc")))
    }
    fcn <- function(ini, xp, y, fcall, jcall)
      (y - do.call("fcall", c(list(xp = xp), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,xp = xp, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #######  Hilber y logan    ---  Logan typo III
  if(model==22){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,d,Y,Tmax, v){
      expr <- expression(Y*((x)^2/((x)^2+d^2)-exp(-(Tmax-(x))/v)))
      eval(expr)
    }
    j <- function(x,d,Y,Tmax, v){
      expr <- expression(Y*((x)^2/((x)^2+d^2)-exp(-(Tmax-(x))/v)))
      c(eval(D(expr, "d")),eval(D(expr, "Y")),eval(D(expr, "Tmax")),eval(D(expr, "v")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  ####  Latin 2
  if(model==23) {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,Tl, p, dt,L){
      expr <- expression(exp(p*x)-exp(-(p*Tl-(x-Tl))/dt)+L)
      eval(expr)
    }
    j <- function(x,Tl, p, dt,L){
      expr <- expression(exp(p*x)-exp(-(p*Tl-(x-Tl))/dt)+L)
      c(eval(D(expr, "Tl")),eval(D(expr, "L")), eval(D(expr, "p")), eval(D(expr, "dt")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  ##  Lineal
  if(model==24)  {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,Inter,Slop){
      expr <- expression(Inter + Slop*x)
      eval(expr)
    }
    j <- function(x,Inter,Slop){
      expr <- expression(Inter + Slop*x)
      c(eval(D(expr, "Inter")),eval(D(expr, "Slop")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  ##  exponencial simple
  if(model==25)  {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,b1,b2){
      expr <- expression(b1*exp(b2*x))
      eval(expr)
    }
    j <- function(x,b1,b2){
      expr <- expression(b1*exp(b2*x))
      c(eval(D(expr, "b1")),eval(D(expr, "b2")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  ##  Tb Model (Logan)
  if(model==26)  {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,sy,b,Tb,DTb){
      expr <- expression(sy*exp(b*(x-Tb)-exp(b*(x-Tb)/DTb)))
      eval(expr)
    }
    j <- function(x,sy,b,Tb,DTb){
      expr <- expression(sy*exp(b*(x-Tb)-exp(b*(x-Tb)/DTb)))
      c(eval(D(expr, "sy")),eval(D(expr, "b")),eval(D(expr, "Tb")),eval(D(expr, "DTb")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  ##  Exponential Model (Logan)
  if(model==27)  {
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,sy,b,Tb){
      expr <- expression(sy*exp(b*(x-Tb)))
      eval(expr)
    }
    j <- function(x,sy,b,Tb){
      expr <- expression(sy*exp(b*(x-Tb)))
      c(eval(D(expr, "sy")),eval(D(expr, "b")),eval(D(expr, "Tb")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  ##  Exponential Tb (Logan)
  if(model==28){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,b,Tmin){
      expr <- expression(exp(b*(x-Tmin))-1)
      eval(expr)
    }
    j <- function(x,b,Tmin){
      expr <- expression(exp(b*(x-Tmin))-1)
      c(eval(D(expr, "b")),eval(D(expr, "Tmin")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  ##  Square root model of Ratkowsky
  if(model==29){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,b,Tb){
      expr <- expression(b*(x-Tb)^2)
      eval(expr)
    }
    j <- function(x,b,Tb){
      expr <- expression(b*(x-Tb)^2)
      c(eval(D(expr, "b")),eval(D(expr, "Tb")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #Davidson
  if(model==30){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,k,a,b){
      expr <- expression(k/(1+exp(a-b*x))) #DM: Se cambio "+b*x" por "-b*x"
      eval(expr)
    }
    j <- function(x,k,a,b){
      expr <- expression(k/(1+exp(a-b*x))) #DM: Se cambio "+b*x" por "-b*x"
      c(eval(D(expr, "k")),eval(D(expr, "a")),eval(D(expr, "b")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  
  #Pradham - 1
  if(model==31){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,R,Tm,To){
      expr <- expression(R*exp((-1/2)*((x-Tm)/To))^2) #DM: Se agrego "^2"
      eval(expr)
    }
    j <- function(x,R,Tm,To){
      expr <- expression(R*exp((-1/2)*((x-Tm)/To))^2) #DM: Se agrego "^2"
      c(eval(D(expr, "R")),eval(D(expr, "Tm")),eval(D(expr, "To")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #Angilletta Jr.
  if(model==32){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,a,b,c,d){
      expr <- expression(a*exp((-1/2)*(abs(x-b)/c)^d)) #DM: Se cambio "abs((x-b)/c)" por "abs(x-b)/c"
      eval(expr)
    }
    j <- function(x,a,b,c,d){
      expr <- expression(a*exp((-1/2)*(abs(x-b)/c)^d)) #DM: Se cambio "abs((x-b)/c)" por "abs(x-b)/c"
      c(eval(D(expr, "a")),eval(D(expr, "b")),eval(D(expr, "c")),eval(D(expr, "d")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #Stinner 2
  if(model==33){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,Rmax,k1,k2,Topc){
      expr <- expression(Rmax*(exp(k1+k2*Topc))/(1+exp(k1+k2*x)))
      eval(expr)
    }
    j <- function(x,Rmax,k1,k2,Topc){
      expr <- expression(Rmax*(exp(k1+k2*Topc))/(1+exp(k1+k2*x)))
      c(eval(D(expr, "Rmax")),eval(D(expr, "k1")),eval(D(expr, "k2")),eval(D(expr, "Topc")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #Hilbert
  if(model==34){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,Tb,Tmax,d,Y,v){
      expr <- expression(Y*((x-Tb)^2/((x-Tb)+d^2)-exp(-(Tmax-(x-Tb))/v))) #DM: Se quito el cuadrado en el denominador
      eval(expr)
    }
    j <- function(x,Tb,Tmax,d,Y,v){
      expr <- expression(Y*((x-Tb)^2/((x-Tb)+d^2)-exp(-(Tmax-(x-Tb))/v))) #DM: Se quito el cuadrado en el denominador
      c(eval(D(expr, "Tb")),eval(D(expr, "Tmax")),eval(D(expr, "d")),eval(D(expr, "Y")),eval(D(expr, "v")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  #Lactin 2
  if(model==35){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,Tl, p, dt){
      expr <- expression(exp(p*x)-exp(p*Tl-(Tl-x)/dt)) #DM: Se omitio un menos "-(p" por "p"  y se invirtio "(x-Tl))/dt" por "(Tl-x)/dt"
      eval(expr)
    }
    j <- function(x,Tl, p, dt){
      expr <- expression(exp(p*x)-exp(p*Tl-(Tl-x)/dt)) #DM: Se omitio un menos "-(p" por "p"  y se invirtio "(x-Tl))/dt" por "(Tl-x)/dt"
      c(eval(D(expr, "Tl")),eval(D(expr, "p")),eval(D(expr, "dt")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #Anlytis-1
  if(model==36){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,P,Tmax, Tmin,n,m){
      expr <- expression(P*(((x-Tmin)/(Tmax-Tmin))^n)*(1-(x-Tmin)/(Tmax-Tmin))^m)
      eval(expr)
    }
    j <- function(x,P,Tmax, Tmin,n,m){
      expr <- expression(P*(((x-Tmin)/(Tmax-Tmin))^n)*(1-(x-Tmin)/(Tmax-Tmin))^m)
      c(eval(D(expr, "P")),eval(D(expr, "Tmax")),eval(D(expr, "Tmin")),eval(D(expr, "n")),eval(D(expr, "m")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #Anlytis-2
  if(model==37){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,P,Tmax, Tmin,n,m){
      expr <- expression((P*(((x-Tmin)/(Tmax-Tmin))^n)*(1-(x-Tmin)/(Tmax-Tmin)))^m) #Se cambio "(P*" por "((P*"
      eval(expr)
    }
    j <- function(x,P,Tmax, Tmin,n,m){
      expr <- expression((P*(((x-Tmin)/(Tmax-Tmin))^n)*(1-(x-Tmin)/(Tmax-Tmin)))^m) #Se cambio "(P*" por "((P*"
      c(eval(D(expr, "P")),eval(D(expr, "Tmax")),eval(D(expr, "Tmin")),eval(D(expr, "n")),eval(D(expr, "m")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  
  #Anlytis-3
  if(model==38){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,a,Tmax, Tmin,n,m){
      expr <- expression(a*((x-Tmin)^n)*(Tmax-x)^m)
      eval(expr)
    }
    j <- function(x,a,Tmax, Tmin,n,m){
      expr <- expression(a*((x-Tmin)^n)*(Tmax-x)^m)
      c(eval(D(expr, "a")),eval(D(expr, "Tmax")),eval(D(expr, "Tmin")),eval(D(expr, "n")),eval(D(expr, "m")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  
  #Allahyari
  if(model==39){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,P,Tmax, Tmin,n,m){
      expr <- expression(P*(((x-Tmin)/(Tmax-Tmin))^n)*(1-((x-Tmin)/(Tmax-Tmin))^m)) #CAMBIO!
      eval(expr)
    }
    j <- function(x,P,Tmax, Tmin,n,m){
      expr <- expression(P*(((x-Tmin)/(Tmax-Tmin))^n)*(1-((x-Tmin)/(Tmax-Tmin))^m)) #CAMBIO!
      c(eval(D(expr, "P")),eval(D(expr, "Tmax")),eval(D(expr, "Tmin")),eval(D(expr, "n")),eval(D(expr, "m")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #Briere 3
  if(model==40){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,aa,To, Tmax){
      expr <- expression(aa*(x-To)*(Tmax-x)^0.5)
      eval(expr)
    }
    j <- function(x,aa,To, Tmax){
      expr <- expression(aa*(x-To)*(Tmax-x)^0.5)
      c(eval(D(expr, "aa")),eval(D(expr, "To")),eval(D(expr, "Tmax")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #Briere 4
  if(model==41){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,aa,To, Tmax,n){
      expr <- expression(aa*(x-To)*(Tmax-x)^(1/n))
      eval(expr)
    }
    j <- function(x,aa,To, Tmax,n){
      expr <- expression(aa*(x-To)*(Tmax-x)^(1/n))
      c(eval(D(expr, "aa")),eval(D(expr, "To")),eval(D(expr, "Tmax")),eval(D(expr, "n")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #Kontodimas-1
  if(model==42){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,aa,Tmin,Tmax){
      expr <- expression(aa*((x-Tmin)^2)*(Tmax-x))
      eval(expr)
    }
    j <- function(x,aa,Tmin,Tmax){
      expr <- expression(aa*((x-Tmin)^2)*(Tmax-x))
      c(eval(D(expr, "aa")),eval(D(expr, "Tmin")),eval(D(expr, "Tmax")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  
  #Kontodimas-2
  if(model==43){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,Dmin,Topt,K,lmda){
      expr <- expression(2/(Dmin*(exp(K*(x-Topt)) + exp((-lmda)*(x-Topt)))))
      eval(expr)
    }
    j <- function(x,Dmin,Topt,K,lmda){
      expr <- expression(2/(Dmin*(exp(K*(x-Topt)) + exp((-lmda)*(x-Topt)))))
      c(eval(D(expr, "Dmin")),eval(D(expr, "Topt")),eval(D(expr, "K")),eval(D(expr, "lmda")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  
  #Kontodimas-3
  if(model==44){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,a1,b1,c1,d1,f1,g1){
      expr <- expression(x*exp(a1-b1/x)/(1 + exp(c1-d1/x) + exp(f1-g1/x)))
      eval(expr)
    }
    j <- function(x,a1,b1,c1,d1,f1,g1){
      expr <- expression(x*exp(a1-b1/x)/(1 + exp(c1-d1/x) + exp(f1-g1/x)))
      c(eval(D(expr, "a1")),eval(D(expr, "b1")),eval(D(expr, "c1")),eval(D(expr, "d1")),eval(D(expr, "f1")),eval(D(expr, "g1")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  
  #Ratkowsky 2
  if(model==45){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,aa,Tmin,Tmax,b){
      expr<-expression((aa*(x-Tmin)*(1-exp((b*(Tmax-x)))))^2) #Agregar 1-exp
      eval(expr)
    }
    j <- function(x,aa,Tmin,Tmax,b){
      expr<-expression((aa*(x-Tmin)*(1-exp((b*(Tmax-x)))))^2) #Agregar 1-exp
      c(eval(D(expr, "aa")),eval(D(expr, "Tmin")),eval(D(expr, "Tmax")),eval(D(expr, "b")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  
  #Janish-1
  if(model==46){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,Dmin,Topt,K){
      expr <- expression(2/(Dmin*(exp(K*(x-Topt)) + exp((-K)*(x-Topt)))))
      eval(expr)
    }
    j <- function(x,Dmin,Topt,K){
      expr <- expression(2/(Dmin*(exp(K*(x-Topt)) + exp((-K)*(x-Topt)))))
      c(eval(D(expr, "Dmin")),eval(D(expr, "Topt")),eval(D(expr, "K")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #Janish-2
  if(model==47){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,c,a,b,Tm){
      expr <- expression(2*c/(a^(x-Tm) + b^(Tm-x)))
      eval(expr)
    }
    j <- function(x,c,a,b,Tm){
      expr <- expression(2*c/(a^(x-Tm) + b^(Tm-x)))
      c(eval(D(expr, "c")),eval(D(expr, "a")),eval(D(expr, "b")),eval(D(expr, "Tm")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #Tanigoshi
  if(model==48){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,a0,a1,a2,a3){
      expr <- expression(a0+a1*x+a2*x^2+a3*x^3)
      eval(expr)
    }
    j <- function(x,a0,a1,a2,a3){
      expr <- expression(a0+a1*x+a2*x^2+a3*x^3)
      c(eval(D(expr, "a0")),eval(D(expr, "a1")),eval(D(expr, "a2")),eval(D(expr, "a3")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #Wang-Lan-Ding
  if(model==49){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,k,a,b,c,Tmin,Tmax,r){
      expr <- expression(k*(1-exp((-a)*(x-Tmin)))*(1-exp(b*(x-Tmax)))/(1+exp((-r)*(x-c))))
      eval(expr)
    }
    j <- function(x,k,a,b,c,Tmin,Tmax,r){
      expr <- expression(k*(1-exp((-a)*(x-Tmin)))*(1-exp(b*(x-Tmax)))/(1+exp((-r)*(x-c))))
      c(eval(D(expr, "k")),eval(D(expr, "a")),eval(D(expr, "b")),eval(D(expr, "c")),eval(D(expr, "Tmin")),eval(D(expr, "Tmax")),eval(D(expr, "r")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  
  ## modelos adaptados para senescencia
  
  #Stinner-3
  if(model==50){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,c1,k1,k2){
      expr <- expression(c1/(1+exp(k1+k2*x)))
      eval(expr)
    }
    j <- function(x,c1,k1,k2){
      expr <- expression(c1/(1+exp(k1+k2*x)))
      c(eval(D(expr, "c1")),eval(D(expr, "k1")),eval(D(expr, "k2")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  
  #Stinner-4
  if(model==51){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,c1,c2,k1,k2,To){
      expr <- expression(c1/(1+exp(k1+k2*x)) + c2/(1+exp(k1+k2*(2*To-x))))
      eval(expr)
    }
    j <- function(x,c1,c2,k1,k2,To){
      expr <- expression(c1/(1+exp(k1+k2*x)) + c2/(1+exp(k1+k2*(2*To-x))))
      c(eval(D(expr, "x")),eval(D(expr, "c1")),eval(D(expr, "c2")),eval(D(expr, "k1")),eval(D(expr, "k2")),eval(D(expr, "To")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #Logan-3
  if(model==52){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,sy,b,Tmin,Tmax,DTb){
      expr <- expression(sy*exp(b*(x-Tmin)-exp(b*Tmax - (Tmax-(x-Tmin))/DTb)))
      eval(expr)
    }
    j <- function(x,sy,b,Tmin,Tmax,DTb){
      expr <- expression(sy*exp(b*(x-Tmin)-exp(b*Tmax - (Tmax-(x-Tmin))/DTb)))
      c(eval(D(expr, "sy")),eval(D(expr, "b")),eval(D(expr, "Tmin")),eval(D(expr, "Tmax")),eval(D(expr, "DTb")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #Logan-4
  if(model==53){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,alph,k,b,Tmin,Tmax,Dt){
      expr <- expression(alph*(1/(1+k*exp(-b*(x-Tmin))) - exp(-(Tmax-(x-Tmin))/Dt)))
      eval(expr)
    }
    j <- function(x,alph,k,b,Tmin,Tmax,Dt){
      expr <- expression(alph*(1/(1+k*exp(-b*(x-Tmin))) - exp(-(Tmax-(x-Tmin))/Dt)))
      c(eval(D(expr, "alph")),eval(D(expr, "k")),eval(D(expr, "b")),eval(D(expr, "Tmin")),eval(D(expr, "Tmax")),eval(D(expr, "Dt")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #Logan-5
  if(model==54){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,alph,k,b,Tmax,Dt){
      expr <- expression(alph*(1/(1+k*exp(-b*x)) - exp(-(Tmax-x)/Dt)))
      eval(expr)
    }
    j <- function(x,alph,k,b,Tmax,Dt){
      expr <- expression(alph*(1/(1+k*exp(-b*x)) - exp(-(Tmax-x)/Dt)))
      c(eval(D(expr, "alph")),eval(D(expr, "k")),eval(D(expr, "b")),eval(D(expr, "Tmax")),eval(D(expr, "Dt")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #Hilber y logan 2
  if(model==55){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,trid,D,Tmax,Dt){
      expr <- expression(trid*( (x^2)/(x^2+D)  - exp(-(Tmax-x)/Dt))) #Cambio de "x/Dt" por "x)/Dt"
      eval(expr)
    }
    j <- function(x,trid,D,Tmax,Dt){
      expr <- expression(trid*( (x^2)/(x^2+D)  - exp(-(Tmax-x)/Dt))) #Cambio de "x/Dt" por "x)/Dt"
      c(eval(D(expr, "trid")),eval(D(expr, "Tmax")),eval(D(expr, "Dt")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #Hilber y logan 3
  if(model==56){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,trid,Tmax,Tmin,D,Dt,Smin){
      expr <- expression(trid*(((x-Tmin)^2)/((x-Tmin)^2 + D) - exp(-(Tmax-(x-Tmin))/Dt)) + Smin)
      eval(expr)
    }
    j <- function(x,trid,Tmax,Tmin,D,Dt,Smin){
      expr <- expression(trid*(((x-Tmin)^2)/((x-Tmin)^2 + D) - exp(-(Tmax-(x-Tmin))/Dt)) + Smin)
      c(eval(D(expr, "trid")),eval(D(expr, "Tmax")),eval(D(expr, "Tmin")),eval(D(expr, "D")),eval(D(expr, "Dt")),eval(D(expr, "Smin")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #Taylor
  if(model==57){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,rm,Topt,Troh){ #DM: se quito Smin
      expr <- expression(rm*exp(-(0.5)*(-(x-Topt)/Troh)^2))
      eval(expr)
    }
    j <- function(x,rm,Topt,Troh){ #DM: se quito Smin
      expr <- expression(rm*exp(-(0.5)*(-(x-Topt)/Troh)^2))
      c(eval(D(expr, "rm")),eval(D(expr, "Topt")),eval(D(expr, "Troh"))) #DM: se quito Smin
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #Lactin 3
  if(model==58){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,Tl, p, dt,lamb){
      expr <- expression(exp(p*x)-exp(p*Tl-(Tl-x)/dt) + lamb) #DM: Se cambi? "(p*Tl-(Tl-x))/dt)" por "p*Tl-(Tl-x)/dt)"
      eval(expr)
    }
    j <- function(x,Tl, p, dt,lamb){
      expr <- expression(exp(p*x)-exp(p*Tl-(Tl-x)/dt) + lamb) #DM: Se cambi? "(p*Tl-(Tl-x))/dt)" por "p*Tl-(Tl-x)/dt)"
      c(eval(D(expr, "Tl")),eval(D(expr, "p")),eval(D(expr, "Dt")),eval(D(expr, "lamb")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
  
  #Sigmoid or logistic
  if(model==59){
    if(length(datashap[,1])>(length(ini)+1)){
      x<-datashap[,1]
      y<-datashap[,2]}else{
        x<-datao[,1]
        y<-datao[,2]}
    ini<-as.list(ini)
    ind <- as.list(ini)
    for (i in names(ind)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    f <- function(x,c1,a,b){
      expr <- expression(c1/(1+exp(a+b*x)))
      eval(expr)
    }
    j <- function(x,c1,a,b){
      expr <- expression(c1/(1+exp(a+b*x)))
      c(eval(D(expr, "c1")),eval(D(expr, "a")),eval(D(expr, "b")))
    }
    fcn <- function(ini, x, y, fcall, jcall)
      (y - do.call("fcall", c(list(x = x), as.list(ini))))
    out <- nls.lm(par = ini, fn = fcn,fcall = f, jcall = j,x = x, y = y)
    stderro<-diag(ginv(out$hessian))
    estimate<-as.data.frame(out$par)
    return(list(estimados=estimate,f=f,stderro=stderro))
  }
}
