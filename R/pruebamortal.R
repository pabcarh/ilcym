#' pruebamortal
#'
#' no details
#' 
#' @param proc a character, type of variable (mortal, taza)
#' @param modelm integer value, non linear model (position)
#' @param datm data frame, temperatures and mortality
#' @param inim numeric vector, initial values to model
#' @param corrx numeric vector of two values, the axis range for X
#' @param corry numeric vector of two values, the axis range for Y
#' @param mini double, initial value of line range
#' @param maxi double, final value of line range
#' @param labx a character
#' @param laby a character
#' @param titulo a character

#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
pruebamortal<-function(proc,modelm,datm,inim,corrx,corry,mini,maxi,labx, laby,titulo)
{
  if(modelm==1)  form<-expression(a*xlinea^2+b*xlinea+c)
  if(modelm==2)  form<-expression(a*sqrt(xlinea)+b*xlinea+c)
  if(modelm==3)  form<-expression(a*(1/sqrt(xlinea))+b*xlinea+c)
  if(modelm==4)  form<-expression(a*(1/xlinea^2)+b*xlinea+c)
  if(modelm==5)  form<-expression(a*(1/xlinea)+b*xlinea+c)
  if(modelm==6)  form<-expression(b*xlinea+a*log(xlinea)+c)
  if(modelm==7)  form<-expression(1/(1+a*exp(-b*((xlinea-c)/d)^2))) # log normal, 4 parameters
  if(modelm==8) form<-expression((a*exp(-b*((xlinea-xo)/c)^2)))
  if(modelm==9)  form<-expression(y0 + a * exp(-0.5 * ((xlinea-x0)/b)^2))
  if(modelm==10) form<-expression(y0 + a * exp(-0.5 * (log(abs(xlinea/x0))/b)^2))
  
  # nuevos modelos
  if(modelm==11) form<-expression(b1+b2*xlinea+b3*xlinea^d)
  if(modelm==12) form<-expression(exp(b1+b2*xlinea+b3*xlinea^2))
  if(modelm==13) form<-expression(1-(b4/(1+b5*exp(b1+b2*xlinea+b3*xlinea^2))))
  if(modelm==14) form<-expression(exp(b1+b2*xlinea+b3*sqrt(xlinea)))
  if(modelm==15) form<-expression(1-(b4/(1+b5*exp(b1+b2*xlinea+b3*sqrt(xlinea)))))
  if(modelm==16) form<-expression(exp(b1+b2*xlinea+b3*(1/sqrt(xlinea))))
  if(modelm==17) form<-expression(1-(b4/(1+b5*exp(b1+b2*xlinea+b3*(1/sqrt(xlinea))))))
  if(modelm==18) form<-expression(exp(b1+b2*xlinea+b3*(1/xlinea)))
  if(modelm==19) form<-expression(1-(b4/(1+b5*exp(b1+b2*xlinea+b3*(1/xlinea)))))
  if(modelm==20) form<-expression(exp(b1+b2*xlinea+b3*xlinea^d))
  if(modelm==21) form<-expression(1-(b4/(1+b5*exp(b1+b2*xlinea+b3*xlinea^d))))
  if(modelm==22) form<-expression(exp(b1+b2*xlinea+b3*log(xlinea)))
  if(modelm==23) form<-expression(1-(b4/(1+b5*exp(b1+b2*xlinea+b3*log(xlinea)))))
  if(modelm==24) form<-expression(1-rm*exp((-0.5)*(-(xlinea-Topt)/Troh)^2))
  if(modelm==25) form<-expression(1-rm*exp((-0.5)*(-(log(xlinea)-log(Topt))/Troh)^2))
  if(modelm==26) form<-expression(1 - 1/(exp((1+exp(-(xlinea-Topt)/B))*(1+exp(-(Topt-xlinea)/B))*H)))
  if(modelm==27) form<-expression(1 - 1/(exp((1+exp(-(xlinea-Tl)/B))*(1+exp(-(Th-xlinea)/B))*H)))
  if(modelm==28) form<-expression(1 - 1/(exp((1+exp(-(xlinea-Topt)/Bl))*(1+exp(-(Topt-xlinea)/Bh))*H)))
  if(modelm==29) form<-expression(1 - 1/(exp((1+exp(-(xlinea-Tl)/Bl))*(1+exp(-(Th-xlinea)/Bh))*H)))
  if(modelm==30) form<-expression(1 - H/(exp(1+exp(-(xlinea-Topt)/B))*(1+exp(-(Topt-xlinea)/B))))
  if(modelm==31) form<-expression(1 - H/(exp(1+exp(-(xlinea-Tl)/B))*(1+exp(-(Th-xlinea)/B))))
  if(modelm==32) form<-expression(1 - H/(exp(1+exp(-(xlinea-Topt)/Bl))*(1+exp(-(Topt-xlinea)/Bh))))
  if(modelm==33) form<-expression(1 - H/(exp(1+exp(-(xlinea-Tl)/Bl))*(1+exp(-(Th-xlinea)/Bh))))
  if(modelm==34) form<-expression(Bm + ((1-1/(1+exp((Hl/1.987)*((1/Tl)-(1/(xlinea+273.15))))+exp((Hh/1.987)*((1/Th)-(1/(xlinea+273.15))))))*(1-Bm))) #DM: se agrego 2 parentesis "Bm + (...)*"
  if(modelm==35) form<-expression((1 - exp(-(exp(a1+b1*xlinea)))) + (1 - exp(-(exp(a2+b2*xlinea)))))
  if(modelm==36) form<-expression((w-xlinea)^(-1))
  if(modelm==37) form<-expression(a1*exp(b1*xlinea) + a2*exp(b2*xlinea))
  if(modelm==38) form<-expression(a1*exp(b1*xlinea) + a2*exp(b2*xlinea)+c1)
  if(modelm==39) form<-expression(a*(abs(xlinea-b))^nn)
  if(modelm==40) form<-expression(a*xlinea*(xlinea-To)*(Tl-xlinea)^(1/d))
  if(modelm==41) form<-expression(exp(a*xlinea*(xlinea-To)*(Tl-xlinea)^(1/d)))
  if(modelm==42) form<-expression(a*((xlinea-Tmin)^n)*(Tmax-xlinea)^m)
  if(modelm==43) form<-expression(1/((Dmin/2) * (exp(k*(xlinea-Tp)) + exp(-(xlinea-Tp)*lamb))))
  if(modelm==44) form<-expression(a*(1-exp(-(xlinea-Tl)/B))*(1-exp(-(Th-xlinea)/B)))
  if(modelm==45) form<-expression(exp(a*(1-exp(-(xlinea-Tl)/Bl))*(1-exp(-(Th-xlinea)/Bh))))
  
  # fin de nuevos modelos
  
  
  if(proc=="mortal") graphics::plot(datm[,1],datm[,2]*100,frame=F,pch=19,col=4,cex=1.3,xlim=corrx,ylim=corry,axes=F,xaxt = "n",xlab=labx,ylab=laby,main=titulo)
  if(proc=="taza") graphics::plot(datm[,1],datm[,2],frame=F,pch=19,col=4,cex=1.3,xlim=corrx,ylim=corry,axes=F,xaxt = "n",xlab=labx,ylab=laby,main=titulo)
  graphics::axis(1, xaxp=c(corrx,5))
  graphics::axis(2,las=2)
  ini<-as.list(inim)
  for (i in names(ini))
  {
    temp <- ini[[i]]
    storage.mode(temp) <- "double"
    assign(i, temp)
  }
  xlinea<-seq(mini,maxi,length=100)
  ylinea<-eval(form)
  if(proc=="mortal") graphics::lines(xlinea,ylinea*100,lwd=2,col=2)
  if(proc=="taza") graphics::lines(xlinea,ylinea,lwd=2,col=2)
  salidas<-list(ini=as.data.frame(ini))
  return(salidas)
}
