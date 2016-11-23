#' pruebaovi
#'
#' no details
#'
#' @param modelm integer value, non linear model (position)
#' @param ovifemal data frame, temperatures and oviposition
#' @param inim numeric vector of the initial values for the modeling
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
pruebaovi<-function(modelm, ovifemal,inim,corrx,corry,mini,maxi,labx,laby,titulo)
{
  x<-ovifemal[,2]
  y<-ovifemal[,3]
  if(modelm==1)  form<-expression((1-exp(-(a*xlinea+b*xlinea^2+c*xlinea^3))))
  if(modelm==2)  form<-expression(pgamma(xlinea,a,b))
  if(modelm==3)  form<-expression(1/(1+exp(a+b*xlinea)))
  if(modelm==4)  form<-expression(1-exp(-a*xlinea^b))
  if(modelm==5)  form<-expression(1-exp(-((xlinea-a)/n)^b))
  if(modelm==6)  form<-expression(pweibull(xlinea,a,b))
  
  plot(ovifemal[,2],ovifemal[,3]*100,frame=F,pch=19,col=4,cex=1.3,xlim=corrx,ylim=corry,axes=F,xaxt = "n",xlab=labx,ylab=laby,main=titulo)
  axis(1, xaxp=c(corrx,5))
  axis(2,las=2)
  ini<-as.list(inim)
  for (i in names(ini))
  {
    temp <- ini[[i]]
    storage.mode(temp) <- "double"
    assign(i, temp)
  }
  xlinea<-seq(mini,maxi,length=1000)
  ylinea<-eval(form)
  lines(xlinea,ylinea*100,lwd=2,col=2)
  salidas<-list(ini=as.data.frame(ini))
  return(salidas)
}
