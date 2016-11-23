#' prueba
#'
#' no details
#'
#' @param model integer value, non linear model (position)
#' @param datashap data frame, estimated rates
#' @param datao data frame
#' @param ini numeric vector of the initial values for the modeling
#' @param corrx numeric vector of two values, the axis range for X
#' @param corry numeric vector of two values, the axis range for y
#' @param punt integer value, extent of the plot
#' @param labx a character
#' @param laby a character
#' @param titulo a character
#' @param grises logical; default is FALSE
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
prueba<-function(model,datashap,datao,ini,corrx,corry,punt,labx,laby,titulo, grises=FALSE){
  if(grises==TRUE){ccol=c("gray10","gray20","gray30")}else{ccol=c(4,1,2)}
  
  if(model==1 || model==2  || model==7 || model==8 || model==9 || model==12){
    datlinea<-datashap
    #plot(datlinea[,1],datlinea[,2],frame=F,pch=19,col=ccol[1],cex=1.3,xlab=labx,ylab=laby,xlim=corrx,ylim=corry,axes=F,xaxt = "n",main=titulo)   ## 1
    #axis(1, xaxp=c(corrx,5))
    #axis(2,las=2)
    #arrows(datlinea[,1],(datlinea[,2]-(datlinea[,2]-datlinea[,3])), datlinea[,1],(datlinea[,2]+(datlinea[,4]-datlinea[,2])), length=0.1,angle=90, code=3,col=ccol[1]) ## 1
    punt<-punt
    #points(datlinea[,1][punt],datlinea[,2][punt],pch=19,col=ccol[2],cex=1.3) ## 2
    #arrows(datlinea[,1][punt],(datlinea[,2][punt]-(datlinea[,2][punt]-datlinea[,3][punt])), datlinea[,1][punt],(datlinea[,2][punt]+(datlinea[,4][punt]-datlinea[,2][punt])), length=0.1,angle=90, code=3,col=ccol[2])  ## 2
    
    if(length(datashap[,1])>=(length(ini)+1)){
      x<-datashap[,1]+273.15
      y<-datashap[,2]
      dataa<-data.frame(x = datlinea[punt,1]+273.15, y =datlinea[punt,2])
    }else{
      x<-datao[,1]+273.15
      y<-datao[,2]
      dataa<-data.frame(x = datlinea[punt,1]+273.15, y = c((datlinea[punt,2]),(datlinea[punt,4]),(datlinea[punt,3])))
    }
    coefi<-as.numeric(stats::coef(stats::lm(dataa[,2]~dataa[,1])))
    xlinea<-seq(0,max(dataa[,1],na.rm=TRUE),length=1000)
    ylinea<-coefi[1]+coefi[2]*(xlinea+273.15)
    #lines(xlinea,ylinea,col=ccol[2],lty=2,lwd=1)  ## 2
    if(model==1){
      ini<-as.list(ini)
      for (i in names(ini)){
        temp <- ini[[i]]
        storage.mode(temp) <- "double"
        assign(i, temp)
      }
      f <- function(x,Ha, Hl,Tl,Hh,Th){
        expr <- expression(((coefi[1]+coefi[2]*((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))) * (x/(((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th))))) * exp((Ha/1.987) * ((1/((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))) - (1/x))))/
                             (1 + exp((Hl/1.987) * ((1/Tl) - (1/x))) + exp((Hh/1.987) * ((1/Th) - (1/x)))))
        eval(expr)
      }
      p<-coefi[1]+coefi[2]*((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))
      lineashap<-seq(0,50,length=1000)
      ylineash<-f((lineashap+273.15),Ha, Hl,Tl,Hh,Th)
      salidas<-list(coefi=coefi,ini=as.data.frame(ini),p=p)
    }
    if(model==2){
      ini<-as.list(ini)
      for (i in names(ini)){
        temp <- ini[[i]]
        storage.mode(temp) <- "double"
        assign(i, temp)
      }
      f <- function(x,To,Ha, Hl,Tl,Hh,Th){
        expr <- expression(((coefi[1]+coefi[2]*To) * (x/(To)) * exp((Ha/1.987) * ((1/(To)) - (1/x))))/
                             (1 + exp((Hl/1.987) * ((1/Tl) - (1/x))) + exp((Hh/1.987) * ((1/Th) - (1/x)))))
        eval(expr)
      }
      p<-coefi[1]+coefi[2]*To
      lineashap<-seq(0,50,length=1000)
      ylineash<-f((lineashap+273.15),To,Ha, Hl,Tl,Hh,Th)
      salidas<-list(coefi=coefi,ini=as.data.frame(ini),p=p)
    }
    if(model==7)   {
      ini<-as.list(ini)
      for (i in names(ini))
      {
        temp <- ini[[i]]
        storage.mode(temp) <- "double"
        assign(i, temp)
      }
      
      
      f <- function(x,To,Ha)
      {
        expr <- expression(((coefi[1]+coefi[2]*To) * (x/(To)) * exp((Ha/1.987) * ((1/To) - (1/x)))))
        eval(expr)
      }
      
      p<-coefi[1]+coefi[2]*To   ### aqui crea un objeto para p el cual no esta incluido en la ecuacion final
      lineashap<-seq(0,50,length=1000)
      ylineash<-f((lineashap+273.15),To,Ha)
      salidas<-list(coefi=coefi,ini=as.data.frame(ini),p=p)
    }
    
    if(model==8)   {
      ini<-as.list(ini)
      for (i in names(ini))
      {
        temp <- ini[[i]]
        storage.mode(temp) <- "double"
        assign(i, temp)
      }
      
      f <- function(x,To,Ha, Hl,Tl)
      {
        expr <- expression(((coefi[1]+coefi[2]*To) * (x/(To)) * exp((Ha/1.987) * ((1/To) - (1/x))))/
                             (1 + exp((Hl/1.987) * ((1/Tl) - (1/x)))))
        eval(expr)
      }
      
      p<-coefi[1]+coefi[2]*To   ### aqui crea un objeto para p el cual no esta incluido en la ecuacion final
      lineashap<-seq(0,50,length=1000)
      ylineash<-f((lineashap+273.15),To,Ha, Hl,Tl)
      salidas<-list(coefi=coefi,ini=as.data.frame(ini),p=p)
    }
    
    if(model==9)   {
      ini<-as.list(ini)
      for (i in names(ini))
      {
        temp <- ini[[i]]
        storage.mode(temp) <- "double"
        assign(i, temp)
      }
      
      f <- function(x,To,Ha,Hh,Th)
      {
        expr <- expression(((coefi[1]+coefi[2]*To) * (x/(To)) * exp((Ha/1.987) * ((1/To) - (1/x))))/
                             (1 + exp((Hh/1.987) * ((1/Th) - (1/x)))))
        eval(expr)
      }
      
      p<-coefi[1]+coefi[2]*To   ### aqui crea un objeto para p el cual no esta incluido en la ecuacion final
      lineashap<-seq(0,50,length=1000)
      ylineash<-f((lineashap+273.15),To,Ha,Hh,Th)
      salidas<-list(coefi=coefi,ini=as.data.frame(ini),p=p)
    }
    
    ################## nuevos modelos de sharpe de michelle
    
    
    if(model==12){
      ini<-as.list(ini)
      for (i in names(ini)){
        temp <- ini[[i]]
        storage.mode(temp) <- "double"
        assign(i, temp)
      }
      f <- function(x,Ha,Hl,Tl,Hh,Th){
        expr <- expression(((coefi[1]+coefi[2]*298.16) * (x/298.16) * exp((Ha/1.987) * ((1/298.16) - (1/x))))/
                             (1 + exp((Hl/1.987) * ((1/Tl) - (1/x))) + exp((Hh/1.987) * ((1/Th) - (1/x)))))
        eval(expr)
      }
      p<-coefi[1]+coefi[2]*298.16
      lineashap<-seq(0,50,length=1000)
      ylineash<-f((lineashap+273.15),Ha, Hl,Tl,Hh,Th)
      salidas<-list(coefi=coefi,ini=as.data.frame(ini),p=p)
    }
    
    
    
    
    
    ################## nuevos modelos de sharpe de michelle
    
    
    #lines(lineashap,ylineash,col=ccol[3],lwd=2)
    return(salidas)
  }else{
    datlinea<-datashap
    #plot(datlinea[,1],datlinea[,2],frame=F,pch=19,col=ccol[1],cex=1.3,xlim=corrx,ylim=corry,axes=F,xaxt = "n",xlab=as.character(labx),ylab=laby) ## 1
    #axis(1, xaxp=c(corrx,5), main=titulo)
    #axis(2,las=2)
    #arrows(datlinea[,1],(datlinea[,2]-(datlinea[,2]-datlinea[,3])), datlinea[,1],(datlinea[,2]+(datlinea[,4]-datlinea[,2])), length=0.1,angle=90, code=3,col=ccol[1]) ## 1
    lineashap<-seq(0,50,length=1000)
    ini<-as.list(ini)
    for (i in names(ini)){
      temp <- ini[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    
    if(model==3)   {
      lineashap<-lineashap+273.15
      expre<-expression((p * (lineashap/(To)) * exp((Ha/1.987) * ((1/To) - (1/lineashap))))/
                          (1 + exp((Hl/1.987) * ((1/Tl) - (1/lineashap))) + exp((Hh/1.987) * ((1/Th) - (1/lineashap)))))
      ylineash<-eval(expre)
      coefi=NULL
    }
    
    if(model==4)   {
      lineashap<-lineashap+273.15
      expre<-expression((p * (lineashap/(To)) * exp((Ha/1.987) * ((1/To) - (1/lineashap)))))
      ylineash<-eval(expre)
      coefi=NULL
    }
    
    if(model==5){
      lineashap<-lineashap+273.15
      expre<-expression((p * (lineashap/(To)) * exp((Ha/1.987) * ((1/To) - (1/lineashap))))/
                          (1 + exp((Hl/1.987) * ((1/Tl) - (1/lineashap)))))
      ylineash<-eval(expre)
      coefi=NULL
    }
    
    if(model==6)   {
      lineashap<-lineashap+273.15
      expre<-expression((p * (lineashap/(To)) * exp((Ha/1.987) * ((1/To) - (1/lineashap))))/
                          (1 + exp((Hh/1.987) * ((1/Th) - (1/lineashap)))))
      ylineash<-eval(expre)
      coefi=NULL
    }
    
    if(model==10)   {
      lineashap<-lineashap+273.15
      expre<-expression((p * (lineashap/(298.16)) * exp((Ha/1.987) * ((1/298.16) - (1/lineashap))))/
                          (1 + exp((Hl/1.987) * ((1/Tl) - (1/lineashap))) + exp((1/1.987) * ((1/1) - (1/lineashap)))))
      ylineash<-eval(expre)
      coefi=NULL
    }
    
    ################## nuevos modelos de sharpe de michelle
    
    if(model==11)   {
      lineashap<-lineashap+273.15
      expre<-expression((p * (lineashap/298.16) * exp((Ha/1.987) * ((1/298.16) - (1/lineashap))))/
                          (1 + exp((Hl/1.987) * ((1/Tl) - (1/lineashap))) + exp((Hh/1.987) * ((1/Th) - (1/lineashap)))))
      ylineash<-eval(expre)
      coefi=NULL
    }
    
    if(model==13)   {
      lineashap<-lineashap+273.15
      expre<-expression((p * (lineashap/298.16) * exp((Ha/1.987) * ((1/298.16) - (1/lineashap))))/
                          (1 + exp((Hl/1.987) * ((1/Tl) - (1/lineashap)))))
      ylineash<-eval(expre)
      coefi=NULL
    }
    
    
    if(model==14)   {
      lineashap<-lineashap+273.15
      expre<-expression((p * (lineashap/298.16) * exp((Ha/1.987) * ((1/298.16) - (1/lineashap))))/
                          (1 + exp((Hh/1.987) * ((1/Th) - (1/lineashap)))))
      ylineash<-eval(expre)
      coefi=NULL
    }
    
    
    ################## fin nuevos modelos de sharpe de michelle
    
    
    if(model==15){
      f <- function(x,Tmin,b){
        expr <- expression(b*(x-Tmin))
        eval(expr)
      }
      ylineash<-f(lineashap,Tmin,b)
      coefi=NULL
    }
    
    if(model==16){
      f <- function(x,b1,b2,b3,b4,b5){
        expr <- expression(b1*10^(-((((x-b3)/(b3-b2)-(1/(1+0.28*b4+0.72*log(1+b4))))+exp(b4*((x-b3)/(b3-b2)-(1/(1+0.28*b4+0.72*log(1+b4))))))/((1+b4)/(1+1.5*b4+0.39*b4^2)))^2)*(1-b5+b5*((((x-b3)/(b3-b2)- (1/(1+0.28*b4+0.72*log(1+b4))))+exp(b4*((x-b3)/(b3-b2)- (1/(1+0.28*b4+0.72*log(1+b4))))))/((1+b4)/(1+1.5*b4+0.39*b4^2)))^2)) #DM: Se agreg? un menos "-v^2"
        eval(expr)
      }
      ylineash<-f(lineashap,b1,b2,b3,b4,b5)
      coefi=NULL
    }
    
    if(model==17){
      f <- function(x,Y,Tmax, p, v){
        expr <- expression(Y*(exp(p*x)-exp(p*Tmax-(Tmax-x)/v)))
        eval(expr)
      }
      ylineash<-f(lineashap,Y,Tmax, p, v)
      coefi=NULL
    }
    
    if(model==18){
      f <- function(x,alfa,k,Tmax, p, v){
        expr <- expression(alfa*((1/(1+k*exp(-p*x)))-exp(-(Tmax-x)/v)))
        eval(expr)
      }
      ylineash<-f(lineashap,alfa,k,Tmax, p, v)
      coefi=NULL
    }
    
    if(model==19){
      f <- function(x,aa,To,Tmax){
        expr <- expression(aa*x*(x-To)*(Tmax-x)^0.5)
        eval(expr)
      }
      ylineash<-f(lineashap,aa,To,Tmax)
      coefi=NULL
    }
    if(model==20){
      f <- function(x,aa,To,Tmax,d){
        expr <- expression(aa*x*(x-To)*(Tmax-x)^d)
        eval(expr)
      }
      ylineash<-f(lineashap,aa,To,Tmax,d)
      coefi=NULL
    }
    
    if(model==21){
      f <- function(xp,Rmax,Topc,k1,k2){
        expr <- expression((Rmax*(1+exp(k1+k2*(Topc))))/(1+exp(k1+k2*xp)))
        eval(expr)
      }
      ylineash<-f(lineashap,Rmax,Topc,k1,k2)
      coefi=NULL
    }
    
    if(model==22){
      f <- function(x,d,Y,Tmax, v){
        expr <- expression(Y*((x)^2/((x)^2+d^2)-exp(-(Tmax-(x))/v)))
        eval(expr)
      }
      
      ylineash<-f(lineashap,d,Y,Tmax, v)
      coefi=NULL
    }
    
    if(model==23){
      f <- function(x, Tl, p, dt, L){
        expr <- expression(exp(p*x)-exp(p*Tl-(Tl-x)/dt)+L)
        eval(expr)
      }
      ylineash<-f(lineashap, Tl, p, dt, L)
      coefi=NULL
    }
    
    if(model==24){
      f <- function(x,Inter,Slop){
        expr <- expression(Inter + Slop*x)
        eval(expr)
      }
      ylineash<-f(lineashap,Inter,Slop)
      coefi=NULL
    }
    
    
    if(model==25){
      f <- function(x,b1,b2){
        expr <- expression(b1*exp(b2*x))
        eval(expr)
      }
      ylineash<-f(lineashap,b1,b2)
      coefi=NULL
    }
    
    if(model==26){
      f <- function(x,sy,b,Tb,DTb){
        expr <- expression(sy*exp(b*(x-Tb)-exp(b*(x-Tb)/DTb)))
        eval(expr)
      }
      ylineash<-f(lineashap,sy,b,Tb,DTb)
      coefi=NULL
    }
    
    if(model==27){
      f <- function(x,sy,b,Tb){
        expr <- expression(sy*exp(b*(x-Tb)))
        eval(expr)
      }
      ylineash<-f(lineashap,sy,b,Tb)
      coefi=NULL
    }
    
    if(model==28){
      f <- function(x,b,Tmin){
        expr <- expression(exp(b*(x-Tmin))-1)
        eval(expr)
      }
      ylineash<-f(lineashap,b,Tmin)
      coefi=NULL
    }
    
    if(model==29){
      f <- function(x,b,Tb){
        expr <- expression(b*(x-Tb)^2)
        eval(expr)
      }
      ylineash<-f(lineashap,b,Tb)
      coefi=NULL
    }
    
    if(model==30){
      f <- function(x,k,a,b){
        expr <- expression(k/(1+exp(a-b*x))) #DM: Se cambio "+b*x" por "-b*x"
        eval(expr)
      }
      ylineash<-f(lineashap,k,a,b)
      coefi=NULL
    }
    
    
    if(model==31){
      f <- function(x,R,Tm,To){
        expr <- expression(R*exp((-1/2)*((x-Tm)/To))^2) #DM: Se agrego "^2"
        eval(expr)
      }
      ylineash<-f(lineashap,R,Tm,To)
      coefi=NULL
    }
    
    if(model==32){
      f <- function(x,a,b,c,d){
        expr <- expression(a*exp((-1/2)*(abs(x-b)/c)^d)) #DM: Se cambio "abs((x-b)/c)" por "abs(x-b)/c"
        eval(expr)
      }
      ylineash<-f(lineashap,a,b,c,d)
      coefi=NULL
    }
    
    if(model==33){
      f <- function(x,Rmax,k1,k2,Topc){
        expr <- expression(Rmax*(exp(k1+k2*Topc))/(1+exp(k1+k2*x)))
        eval(expr)
      }
      ylineash<-f(lineashap,Rmax,k1,k2,Topc)
      coefi=NULL
    }
    
    
    if(model==34){
      f <- function(x,Tb,Tmax,d,Y,v){
        expr <- expression(Y*((x-Tb)^2/((x-Tb)+d^2)-exp(-(Tmax-(x-Tb))/v))) #DM: Se elevo al cuadrado en el numerador y se quito en el denominador
        eval(expr)
      }
      ylineash<-f(lineashap,Tb,Tmax,d,Y,v)
      coefi=NULL
    }
    
    
    
    if(model==35){
      f <- function(x,Tl, p, dt){
        expr <- expression(exp(p*x)-exp(p*Tl-(Tl-x)/dt)) #DM: Se omitio un menos "-(p" por "p"  y se invirtio "(x-Tl))/dt" por "(Tl-x)/dt"
        eval(expr)
      }
      ylineash<-f(lineashap,Tl, p, dt)
      coefi=NULL
    }
    
    
    if(model==36){
      f <- function(x,P,Tmax, Tmin,n,m){
        expr <- expression(P*(((x-Tmin)/(Tmax-Tmin))^n)*(1-(x-Tmin)/(Tmax-Tmin))^m)
        eval(expr)
      }
      ylineash<-f(lineashap,P,Tmax, Tmin,n,m)
      coefi=NULL
    }
    
    
    
    if(model==37){
      f <- function(x,P,Tmax, Tmin,n,m){
        expr <- expression((P*(((x-Tmin)/(Tmax-Tmin))^n)*(1-(x-Tmin)/(Tmax-Tmin)))^m) #Se cambio "(P*" por "((P*"
        eval(expr)
      }
      ylineash<-f(lineashap,P,Tmax, Tmin,n,m)
      coefi=NULL
    }
    
    
    if(model==38){
      f <- function(x,a,Tmax, Tmin,n,m){
        expr <- expression(a*((x-Tmin)^n)*(Tmax-x)^m)
        eval(expr)
      }
      ylineash<-f(lineashap,a,Tmax, Tmin,n,m)
      coefi=NULL
    }
    
    
    
    if(model==39){
      f <- function(x,P,Tmax, Tmin,n,m){
        expr <- expression(P*(((x-Tmin)/(Tmax-Tmin))^n)*(1-((x-Tmin)/(Tmax-Tmin))^m)) #CAMBIO!
        eval(expr)
      }
      ylineash<-f(lineashap,P,Tmax, Tmin,n,m)
      coefi=NULL
    }
    
    
    
    if(model==40){
      f <- function(x,aa,To, Tmax){
        expr <- expression(aa*(x-To)*(Tmax-x)^0.5)
        eval(expr)
      }
      ylineash<-f(lineashap,aa,To, Tmax)
      coefi=NULL
    }
    
    
    if(model==41){
      f <- function(x,aa,To, Tmax,n){
        expr <- expression(aa*(x-To)*(Tmax-x)^(1/n))
        eval(expr)
      }
      ylineash<-f(lineashap,aa,To, Tmax,n)
      coefi=NULL
    }
    
    
    
    if(model==42){
      f <- function(x,aa,Tmin,Tmax){
        expr <- expression(aa*((x-Tmin)^2)*(Tmax-x))
        eval(expr)
      }
      ylineash<-f(lineashap,aa,Tmin,Tmax)
      coefi=NULL
    }
    
    
    
    if(model==43){
      f <- function(x,Dmin,Topt,K,lmda){
        expr <- expression(2/(Dmin*(exp(K*(x-Topt)) + exp((-lmda)*(x-Topt)))))
        eval(expr)
      }
      ylineash<-f(lineashap,Dmin,Topt,K,lmda)
      coefi=NULL
    }
    
    
    if(model==44){
      f <- function(x,a1,b1,c1,d1,f1,g1){
        expr <- expression(x*exp(a1-b1/x)/(1 + exp(c1-d1/x) + exp(f1-g1/x)))
        eval(expr)
      }
      ylineash<-f(lineashap,a1,b1,c1,d1,f1,g1)
      coefi=NULL
    }
    
    
    
    if(model==45){
      f <- function(x,aa,Tmin,Tmax,b){
        expr<-expression((aa*(x-Tmin)*(1-exp((b*(Tmax-x)))))^2) #Agregar 1-exp
        eval(expr)
      }
      ylineash<-f(lineashap,aa,Tmin,Tmax,b)
      coefi=NULL
    }
    
    
    
    if(model==46){
      f <- function(x,Dmin,Topt,K){
        expr <- expression(2/(Dmin*(exp(K*(x-Topt)) + exp((-K)*(x-Topt)))))
        eval(expr)
      }
      ylineash<-f(lineashap,Dmin,Topt,K)
      coefi=NULL
    }
    
    
    
    if(model==47){
      f <- function(x,c,a,b,Tm){
        expr <- expression(2*c/(a^(x-Tm) + b^(Tm-x)))
        eval(expr)
      }
      ylineash<-f(lineashap,c,a,b,Tm)
      coefi=NULL
    }
    
    
    if(model==48){
      f <- function(x,a0,a1,a2,a3){
        expr <- expression(a0+a1*x+a2*x^2+a3*x^3)
        eval(expr)
      }
      ylineash<-f(lineashap,a0,a1,a2,a3)
      coefi=NULL
    }
    
    
    if(model==49){
      f <- function(x,k,a,b,c,Tmin,Tmax,r){
        expr <- expression(k*(1-exp((-a)*(x-Tmin)))*(1-exp(b*(x-Tmax)))/(1+exp((-r)*(x-c))))
        eval(expr)
      }
      ylineash<-f(lineashap,k,a,b,c,Tmin,Tmax,r)
      coefi=NULL
    }
    
    ## funciones adaptadas a senescencia
    
    if(model==50){
      f <- function(x,c1,k1,k2){
        expr <- expression(c1/(1+exp(k1+k2*x)))
        eval(expr)
      }
      ylineash<-f(lineashap,c1,k1,k2)
      coefi=NULL
    }
    
    
    if(model==51){
      f <- function(x,c1,c2,k1,k2,To){
        expr <- expression(c1/(1+exp(k1+k2*x)) + c2/(1+exp(k1+k2*(2*To-x))))
        eval(expr)
      }
      ylineash<-f(lineashap,c1,c2,k1,k2,To)
      coefi=NULL
    }
    
    
    
    if(model==52){
      f <- function(x,sy,b,Tmin,Tmax,Dtb){
        expr <- expression(sy*exp(b*(x-Tmin)-exp(b*Tmax - (Tmax-(x-Tmin))/DTb)))
        eval(expr)
      }
      ylineash<-f(lineashap,sy,b,Tmin,Tmax,Dtb)
      coefi=NULL
    }
    
    if(model==53){
      f <- function(x,alph,k,b,Tmin,Tmax,Dt){
        expr <- expression(alph*(1/(1+k*exp(-b*(x-Tmin))) - exp(-(Tmax-(x-Tmin))/Dt)))
        eval(expr)
      }
      ylineash<-f(lineashap,alph,k,b,Tmin,Tmax,Dt)
      coefi=NULL
    }
    
    if(model==54){
      f <- function(x,alph,k,b,Tmax,Dt){
        expr <- expression(alph*(1/(1+k*exp(-b*x)) - exp(-(Tmax-x)/Dt)))
        eval(expr)
      }
      ylineash<-f(lineashap,alph,k,b,Tmax,Dt)
      coefi=NULL
    }
    
    if(model==55){
      f <- function(x,trid,Tmax,Dt){
        expr <- expression(trid*( (x^2)/(x^2+D)  - exp(-(Tmax-x)/Dt))) #Cambio de "x/Dt" por "x)/Dt"
        eval(expr)
      }
      ylineash<-f(lineashap,trid,Tmax,Dt)
      coefi=NULL
    }
   
    if(model==57){
      f <- function(x,rm,Topt,Troh){ #DM: se quito Smin
        expr <- expression(rm*exp(-(0.5)*(-(x-Topt)/Troh)^2))
        eval(expr)
      }
      ylineash<-f(lineashap,rm,Topt,Troh) #DM: se quito Smin
      coefi=NULL
    }
    
    if(model==58){
      f <- function(x,Tl, p, dt,lamb){
        expr <- expression(exp(p*x)-exp(p*Tl-(Tl-x)/dt) + lamb) #DM: Se cambi? "(p*Tl-(Tl-x))/dt)" por "p*Tl-(Tl-x)/dt)"  
        eval(expr)
      }
      ylineash<-f(lineashap,Tl, p, dt,lamb)
      coefi=NULL
    }
    
    if(model==59){
      f <- function(x,c1,a,b){
        expr <- expression(c1/(1+exp(a+b*x)))
        eval(expr)
      }
      ylineash<-f(lineashap,c1,a,b)
      coefi=NULL
    }
    
    
    #if(model==3 || model==4 || model==5 || model==6 || model==10 || model==11 || model==13 || model==14) {lines(lineashap-273.15,ylineash,col=ccol[3],lwd=2)}else{lines(lineashap,ylineash,col=ccol[3],lwd=2)} ##  3
    salidas<-list(ini=as.data.frame(ini),coefi=coefi)                                                  
    return(salidas)
    
  }
}
