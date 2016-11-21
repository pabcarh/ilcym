#' grafmort
#'
#' no details
#' 
#' @param proc a character, type of variable ("mortal", "taza")
#' @param modelm a character, non linear model
#' @param estimor numeric vector, estimated coefficients of the nonlinear model
#' @param g formula, equation
#' @param datm data frame, temperatures and mortality
#' @param corrx numeric vector of two values, the axis range for X
#' @param corry numeric vector of two values, the axis range for Y
#' @param mini double, initial value of line range
#' @param maxi double, final value of line range
#' @param limit a character, plot the confidence interval ("si", "no")
#' @param tam integer value, extent of the plot
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
grafmort <- function(proc, modelm, estimor, g, datm, corrx=NULL, corry=NULL, mini, maxi, limit,tam,labx, laby,titulo,grises=FALSE){
  if(grises==TRUE){ccol=c("gray20","gray30")}else{ccol=c(4,2)}
  data<-datm
  x<-data[,1]
  estshap<-as.list(estimor)
  finalshap<-estshap
  for (i in names(finalshap)) {
    temp <- finalshap[[i]]
    storage.mode(temp) <- "double"
    assign(i, temp)
  }
  xl<-seq(mini,maxi,length=2000)
  if(modelm==1 ||modelm==2 || modelm==3 || modelm==4 || modelm==5 || modelm==6)
  {
    if(modelm==1)  fom<-y~a*x^2+b*x+c
    if(modelm==2)  fom<-y~a*sqrt(x)+b*x+c
    if(modelm==3)  fom<-y~a*(1/sqrt(x))+b*x+c
    if(modelm==4)  fom<-y~a*(1/x^2)+b*x+c
    if(modelm==5)  fom<-y~a*(1/x)+b*x+c
    if(modelm==6)  fom<-y~b*x+a*log(x)+c
    forex <- fom[[length(fom)]]
    ylexpx<-as.expression(forex)
    ylx<-eval(ylexpx)
    x<-xl
    ylexp<-as.expression(forex)
    yl<-eval(ylexp)
  }
  
  if(modelm==9) {
    exprex<-expression(y0 + a * exp(-0.5 * ((x-x0)/b)^2))
    ylx<-eval(exprex)
    expre<-expression(y0 + a * exp(-0.5 * ((xl-x0)/b)^2))
    yl<-eval(expre)
  }
  if(modelm==7) {
    exprex<-expression(1/(1+a*exp(-b*((x-c)/d)^2)))
    ylx<-eval(exprex)
    expre<-expression(1/(1+a*exp(-b*((xl-c)/d)^2)))
    yl<-eval(expre)
  }
  if(modelm==8) {
    exprex<-expression(a*exp(-b*((x-xo)/c)^2))
    ylx<-eval(exprex)
    expre<-expression(a*exp(-b*((xl-xo)/c)^2))
    yl<-eval(expre)
  }
  if(modelm==10) {
    exprex<-expression(y0 + a * exp(-0.5 * (log(abs(x/x0))/b)^2))
    ylx<-eval(exprex)
    expre<-expression(y0 + a * exp(-0.5 * (log(abs(xl/x0))/b)^2))
    yl<-eval(expre)
  }
  
  if(modelm==11) {
    exprex<-expression(b1+b2*x+b3*x^d)
    ylx<-eval(exprex)
    expre<-expression(b1+b2*xl+b3*xl^d)
    yl<-eval(expre)
  }
  
  
  if(modelm==12) {
    exprex<-expression(exp(b1+b2*x+b3*x^2))
    ylx<-eval(exprex)
    expre<-expression(exp(b1+b2*xl+b3*xl^2))
    yl<-eval(expre)
  }
  
  if(modelm==13) {
    exprex<-expression(1-(b4/(1+b5*exp(b1+b2*x+b3*x^2))))
    ylx<-eval(exprex)
    expre<-expression(1-(b4/(1+b5*exp(b1+b2*xl+b3*xl^2))))
    yl<-eval(expre)
  }
  
  if(modelm==14) {
    exprex<-expression(exp(b1+b2*x+b3*sqrt(x)))
    ylx<-eval(exprex)
    expre<-expression(exp(b1+b2*xl+b3*sqrt(xl)))
    yl<-eval(expre)
  }
  
  if(modelm==15) {
    exprex<-expression(1-(b4/(1+b5*exp(b1+b2*x+b3*sqrt(x)))))
    ylx<-eval(exprex)
    expre<-expression(1-(b4/(1+b5*exp(b1+b2*xl+b3*sqrt(xl)))))
    yl<-eval(expre)
  }
  
  if(modelm==16) {
    exprex<-expression(exp(b1+b2*x+b3*(1/sqrt(x))))
    ylx<-eval(exprex)
    expre<-expression(exp(b1+b2*xl+b3*(1/sqrt(xl))))
    yl<-eval(expre)
  }
  
  if(modelm==17) {
    exprex<-expression(1-(b4/(1+b5*exp(b1+b2*x+b3*(1/sqrt(x))))))
    ylx<-eval(exprex)
    expre<-expression(1-(b4/(1+b5*exp(b1+b2*xl+b3*(1/sqrt(xl))))))
    yl<-eval(expre)
  }
  
  if(modelm==18) {
    exprex<-expression(exp(b1+b2*x+b3*(1/x)))
    ylx<-eval(exprex)
    expre<-expression(exp(b1+b2*xl+b3*(1/xl)))
    yl<-eval(expre)
  }
  
  if(modelm==19) {
    exprex<-expression(1-(b4/(1+b5*exp(b1+b2*x+b3*(1/x)))))
    ylx<-eval(exprex)
    expre<-expression(1-(b4/(1+b5*exp(b1+b2*xl+b3*(1/xl)))))
    yl<-eval(expre)
  }
  
  if(modelm==20) {
    exprex<-expression(exp(b1+b2*x+b3*x^d))
    ylx<-eval(exprex)
    expre<-expression(exp(b1+b2*xl+b3*xl^d))
    yl<-eval(expre)
  }
  
  if(modelm==21) {
    exprex<-expression(1-(b4/(1+b5*exp(b1+b2*x+b3*x^d))))
    ylx<-eval(exprex)
    expre<-expression(1-(b4/(1+b5*exp(b1+b2*xl+b3*xl^d))))
    yl<-eval(expre)
  }
  
  if(modelm==22) {
    exprex<-expression(exp(b1+b2*x+b3*log(x)))
    ylx<-eval(exprex)
    expre<-expression(exp(b1+b2*xl+b3*log(xl)))
    yl<-eval(expre)
  }
  
  if(modelm==23) {
    exprex<-expression(1-(b4/(1+b5*exp(b1+b2*x+b3*log(x)))))
    ylx<-eval(exprex)
    expre<-expression(1-(b4/(1+b5*exp(b1+b2*xl+b3*log(xl)))))
    yl<-eval(expre)
  }
  
  if(modelm==24) {
    exprex<-expression(1-rm*exp((-0.5)*(-(x-Topt)/Troh)^2))
    ylx<-eval(exprex)
    expre<-expression(1-rm*exp((-0.5)*(-(xl-Topt)/Troh)^2))
    yl<-eval(expre)
  }
  
  if(modelm==25) {
    exprex<-expression(1-rm*exp((-0.5)*(-(log(x)-log(Topt))/Troh)^2))
    ylx<-eval(exprex)
    expre<-expression(1-rm*exp((-0.5)*(-(log(xl)-log(Topt))/Troh)^2))
    yl<-eval(expre)
  }
  
  if(modelm==26) {
    exprex<-expression(1 - 1/(exp((1+exp(-(x-Topt)/B))*(1+exp(-(Topt-x)/B))*H)))
    ylx<-eval(exprex)
    expre<-expression(1 - 1/(exp((1+exp(-(xl-Topt)/B))*(1+exp(-(Topt-xl)/B))*H)))
    yl<-eval(expre)
  }
  
  if(modelm==27) {
    exprex<-expression(1 - 1/(exp((1+exp(-(x-Tl)/B))*(1+exp(-(Th-x)/B))*H)))
    ylx<-eval(exprex)
    expre<-expression(1 - 1/(exp((1+exp(-(xl-Tl)/B))*(1+exp(-(Th-xl)/B))*H)))
    yl<-eval(expre)
  }
  
  if(modelm==28) {
    exprex<-expression(1 - 1/(exp((1+exp(-(x-Topt)/Bl))*(1+exp(-(Topt-x)/Bh))*H)))
    ylx<-eval(exprex)
    expre<-expression(1 - 1/(exp((1+exp(-(xl-Topt)/Bl))*(1+exp(-(Topt-xl)/Bh))*H)))
    yl<-eval(expre)
  }
  
  if(modelm==29) {
    exprex<-expression(1 - 1/(exp((1+exp(-(x-Tl)/Bl))*(1+exp(-(Th-x)/Bh))*H)))
    ylx<-eval(exprex)
    expre<-expression(1 - 1/(exp((1+exp(-(xl-Tl)/Bl))*(1+exp(-(Th-xl)/Bh))*H)))
    yl<-eval(expre)
  }
  
  if(modelm==30) {
    exprex<-expression(1 - H/(exp(1+exp(-(x-Topt)/B))*(1+exp(-(Topt-x)/B))))
    ylx<-eval(exprex)
    expre<-expression(1 - H/(exp(1+exp(-(xl-Topt)/B))*(1+exp(-(Topt-xl)/B))))
    yl<-eval(expre)
  }
  
  if(modelm==31) {
    exprex<-expression(1 - H/(exp(1+exp(-(x-Tl)/B))*(1+exp(-(Th-x)/B))))
    ylx<-eval(exprex)
    expre<-expression(1 - H/(exp(1+exp(-(xl-Tl)/B))*(1+exp(-(Th-xl)/B))))
    yl<-eval(expre)
  }
  
  if(modelm==32) {
    exprex<-expression(1 - H/(exp(1+exp(-(x-Topt)/Bl))*(1+exp(-(Topt-x)/Bh))))
    ylx<-eval(exprex)
    expre<-expression(1 - H/(exp(1+exp(-(xl-Topt)/Bl))*(1+exp(-(Topt-xl)/Bh))))
    yl<-eval(expre)
  }
  
  if(modelm==33) {
    exprex<-expression(1 - H/(exp(1+exp(-(x-Tl)/Bl))*(1+exp(-(Th-x)/Bh))))
    ylx<-eval(exprex)
    expre<-expression(1 - H/(exp(1+exp(-(xl-Tl)/Bl))*(1+exp(-(Th-xl)/Bh))))
    yl<-eval(expre)
  }
  
  if(modelm==34) {
    exprex<-expression(Bm + ((1-1/(1+exp((Hl/1.987)*((1/Tl)-(1/(x+273.15))))+exp((Hh/1.987)*((1/Th)-(1/(x+273.15))))))*(1-Bm))) #DM: se agrego 2 parentesis "Bm + (...)*"
    ylx<-eval(exprex)
    expre<-expression(Bm + ((1-1/(1+exp((Hl/1.987)*((1/Tl)-(1/(xl+273.15))))+exp((Hh/1.987)*((1/Th)-(1/(xl+273.15))))))*(1-Bm))) #DM: se agrego 2 parentesis "Bm + (...)*"
    yl<-eval(expre)
  }
  
  if(modelm==35) {
    exprex<-expression((1 - exp(-(exp(a1+b1*x)))) + (1 - exp(-(exp(a2+b2*x)))))
    ylx<-eval(exprex)
    expre<-expression((1 - exp(-(exp(a1+b1*xl)))) + (1 - exp(-(exp(a2+b2*xl)))))
    yl<-eval(expre)
  }
  
  
  if(modelm==36) {
    exprex<-expression((w-x)^(-1))
    ylx<-eval(exprex)
    expre<-expression((w-xl)^(-1))
    yl<-eval(expre)
  }
  
  if(modelm==37) {
    exprex<-expression(a1*exp(b1*x) + a2*exp(b2*x))
    ylx<-eval(exprex)
    expre<-expression(a1*exp(b1*xl) + a2*exp(b2*xl))
    yl<-eval(expre)
  }
  
  if(modelm==38) {
    exprex<-expression(a1*exp(b1*x) + a2*exp(b2*x) + c1)
    ylx<-eval(exprex)
    expre<-expression(a1*exp(b1*xl) + a2*exp(b2*xl) + c1)
    yl<-eval(expre)
  }
  
  if(modelm==39) {
    exprex<-expression(a*(abs(x-b))^nn)
    ylx<-eval(exprex)
    expre<-expression(a*(abs(xl-b))^nn)
    yl<-eval(expre)
  }
  
  if(modelm==40) {
    exprex<-expression(a*x*(x-To)*(Tl-x)^(1/d))
    ylx<-eval(exprex)
    expre<-expression(a*xl*(xl-To)*(Tl-xl)^(1/d))
    yl<-eval(expre)
  }
  
  if(modelm==41) {
    exprex<-expression(exp(a*x*(x-To)*(Tl-x)^(1/d)))
    ylx<-eval(exprex)
    expre<-expression(exp(a*xl*(xl-To)*(Tl-xl)^(1/d)))
    yl<-eval(expre)
  }
  
  if(modelm==42) {
    exprex<-expression(a*((x-Tmin)^n)*(Tmax-x)^m)
    ylx<-eval(exprex)
    expre<-expression(a*((xl-Tmin)^n)*(Tmax-xl)^m)
    yl<-eval(expre)
  }
  
  if(modelm==43) {
    exprex<-expression(1/((Dmin/2) * (exp(k*(x-Tp)) + exp(-(x-Tp)*lamb))))
    ylx<-eval(exprex)
    expre<-expression(1/((Dmin/2) * (exp(k*(xl-Tp)) + exp(-(xl-Tp)*lamb))))
    yl<-eval(expre)
  }
  
  if(modelm==44) {
    exprex<-expression(a*(1-exp(-(x-Tl)/B))*(1-exp(-(Th-x)/B)))
    ylx<-eval(exprex)
    expre<-expression(a*(1-exp(-(xl-Tl)/B))*(1-exp(-(Th-xl)/B)))
    yl<-eval(expre)
  }
  
  if(modelm==45) {
    exprex<-expression(exp(a*(1-exp(-(x-Tl)/Bl))*(1-exp(-(Th-x)/Bh))))
    ylx<-eval(exprex)
    expre<-expression(exp(a*(1-exp(-(xl-Tl)/Bl))*(1-exp(-(Th-xl)/Bh))))
    yl<-eval(expre)
  }
  
  
  
  df<-length(data[,1])-length(estshap)
  sdli<-sqrt(sum((data[,2]-ylx)^2)/df)
  linf<-yl+sdli*qt(0.025,df)
  lsup<-yl-sdli*qt(0.025,df)
  if(limit=="yes"){
    if(proc=="mortal"){
      par(cex=tam)
      if(is.null(corrx)){corrx=c(min(data[,1],0),max(data[,1])*1.2);corrx2=seq(0,10*round(corrx[2]/10,1),5)}else{corrx2=seq(corrx[1],corrx[2],5)} ## cambio
      if(is.null(corry)){corry=c(0,100*(max(data[,2])*1.01));corry2=seq(0,round(max(corry)),10)}else{corry2=seq(corry[1],corry[2],10)}              ## cambio
      
      plot(data[,1],data[,2]*100,frame=F,col=ccol[1],xlim=corrx,ylim=corry,xlab=labx,ylab=laby,pch=19,axes=F,xaxt = "n",main=titulo) ## 4
      #axis(1, xaxp=c(corrx,5))
      #axis(2,las=2)
      axis(1, corrx2)  ## cambio
      axis(2, corry2,las=2) ## cambio
      lines(xl, yl*100, lwd=2,col=ccol[2]) ## 2
      lines(xl, linf*100, lwd=1,col=ccol[1],lty=2) ## 4
      lines(xl, lsup*100, lwd=1,col=ccol[1],lty=2) ## 4
      #arrows(data[,1],(data[,2]-(data[,2]-data[,3])), data[,1],(data[,2]+(data[,4]-data[,2])), length=0.1,angle=90, code=3,col=ccol[1]) ## 4
      valx=xl[100*yl<=100];valx1=min(valx,na.rm=TRUE);valx2=max(valx,na.rm=TRUE) # Se agrego na.rm=TRUE
      return(list(valxs=c(valx1,valx2)))
    }
    if(proc=="taza"){
      par(cex=tam)
      if(is.null(corrx)){corrx=c(min(data[,1],0),max(data[,1])*1.2);corrx2=seq(0,10*round(corrx[2]/10,1),5)}else{corrx2=seq(corrx[1],corrx[2],5)} ## cambio
      if(is.null(corry)){corry=c(0,max(data[,2])*1.01);corry2=seq(0,1,0.1)}else{corry2=seq(corry[1],corry[2],0.1)}              ## cambio
      
      plot(data[,1],data[,2],frame=F,col=ccol[1],xlim=corrx,ylim=corry,xlab=labx,ylab=laby,pch=19,axes=F,xaxt = "n",main=titulo) ## 4
      #axis(1, xaxp=c(corrx,5))
      #axis(2,las=2)
      axis(1, corrx2)  ## cambio
      axis(2, corry2,las=2) ## cambio
      lines(xl, yl, lwd=2,col=ccol[2]) ## 2
      lines(xl, linf, lwd=1,col=ccol[1],lty=2) ## 4
      lines(xl, lsup, lwd=1,col=ccol[1],lty=2) ## 4
      #arrows(data[,1],(data[,2]-(data[,2]-data[,3])), data[,1],(data[,2]+(data[,4]-data[,2])), length=0.1,angle=90, code=3,col=4)
    }
  }
  if(limit=="no"){
    if(proc=="mortal"){
      par(cex=tam)
      if(is.null(corrx)){corrx=c(min(data[,1],0),max(data[,1])*1.2);corrx2=seq(0,10*round(corrx[2]/10,1),5)}else{corrx2=seq(corrx[1],corrx[2],5)} ## cambio
      if(is.null(corry)){corry=c(0,100*(max(data[,2])*1.1));corry2=seq(0,round(max(corry)),10)}else{corry2=seq(corry[1],corry[2],10)}              ## cambio
      plot(data[,1],data[,2]*100,frame=F,col=ccol[1],xlim=corrx,ylim=corry,xlab=labx,ylab=laby,pch=19,axes=F,xaxt = "n",main=titulo) ## 4
      #axis(1, xaxp=c(corrx,5))
      #axis(2,las=2)
      axis(1, corrx2)  ## cambio
      axis(2, corry2,las=2) ## cambio
      lines(xl, yl*100, lwd=2,col=ccol[2]) ## 2
      valx=xl[100*yl<=100];valx1=min(valx,na.rm=TRUE);valx2=max(valx,na.rm=TRUE) # Se agrego na.rm=TRUE
      return(list(valxs=c(valx1,valx2)))
    }
    if(proc=="taza"){
      par(cex=tam)
      if(is.null(corrx)){corrx=c(min(data[,1],0),max(data[,1])*1.2);corrx2=seq(0,10*round(corrx[2]/10,1),5)}else{corrx2=seq(corrx[1],corrx[2],5)} ## cambio
      if(is.null(corry)){corry=c(0,max(data[,2])*1.1);corry2=seq(0,1,0.1)}else{corry2=seq(corry[1],corry[2],0.1)}              ## cambio
      plot(data[,1],data[,2],frame=F,col=ccol[1],xlim=corrx,ylim=corry,xlab=labx,ylab=laby,pch=19,axes=F,xaxt = "n",main=titulo) ## 4
      #axis(1, xaxp=c(corrx,5))
      #axis(2,las=2)
      axis(1, corrx2)  ## cambio
      axis(2, corry2,las=2) ## cambio
      lines(xl, yl, lwd=2,col=ccol[2]) ## 2
    }
  }
}