#' grafovi
#'
#' no details
#' 
#' @param modelm integer value, non linear model (position)
#' @param estimor numeric vector, estimated coefficients of the nonlinear model
#' @param g formula, equation
#' @param ovifemal data frame, temperatures and oviposition
#' @param corrx numeric vector of two values, the axis range for X
#' @param mini double, initial value of line range
#' @param maxi double, final value of line range
#' @param ldx double, x position of the legend
#' @param ldy double, y position of the legend
#' @param sdli numeric vector, standar deviation estimated
#' @param qto double, quantile for a normal distribution (P=0.95, etc)
#' @param limit numeric vector, standar deviation estimated
#' @param tam integer value, extent of the plot
#' @param labx a character
#' @param laby a character
#' @param titulo a character
#' @param grises logical; default is FALSE
#' @param xi double, estimated median from female senescence
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
grafovi <- function(modelm, estimor, g, ovifemal, corrx=NULL,  mini, maxi, ldx,ldy,sdli,qto,limit,tam,labx=NULL, laby=NULL, titulo=NULL,grises=FALSE,xi){
  if(grises==TRUE)
  {colo = c("gray5","gray10","gray15","gray20","gray25","gray30","gray35","gray40","gray45","gray50","gray55","gray60","gray65")
  ccol=c("gray20","gray30")
  }else{colo <- 15:28;ccol=c(4,1)}
  x<-ovifemal[,2]
  y<-ovifemal[,3]
  estshap<-as.list(estimor)
  finalshap<-estshap
  for (i in names(finalshap)) {
    temp <- finalshap[[i]]
    storage.mode(temp) <- "double"
    assign(i, temp)
  }
  xl<-seq(mini,maxi,length=1000)
  if(modelm==1){  f<-function(x){1-exp(-(a*x+b*x^2+c*x^3))}}
  if(modelm==2){  f<-function(x){pgamma(x,a,b)}}
  if(modelm==3){  f<-function(x){1/(1+exp(a+b*x))}}
  if(modelm==4){  f<-function(x){1-exp(-a*x^b)}}
  if(modelm==5){  f<-function(x){1-exp(-((x-a)/n)^b)}}
  if(modelm==6){  f<-function(x){pweibull(x,a,b)}}
  yl<-f(xl)
  
  nag <- length(levels(factor(ovifemal[, 1])))
  not <- levels(factor(ovifemal[, 1]))
  agente <- subset(ovifemal, ovifemal[, 1] == not[1])
  par(cex=tam)
  
  if(is.null(corrx)){corrx=c(0,max(ovifemal[,2])+max(ovifemal[,2])*0.1)} ## nuevo
  corrx2=seq(corrx[1],corrx[2]+0.3,0.3) ## cambio
  plot(agente[, 2], agente[, 3]*100,frame=F,ylab=laby,
       xlab=labx,pch=19,cex=1.3,col=ccol[2],xlim=corrx,ylim=c(0,100),axes=F,xaxt = "n") ## 1
  axis(1, corrx2)  ## cambio
  #axis(1, xaxp=c(corrx,5))
  axis(2,las=2)
  lines(c(-1,xi),c(50,50),lwd=1,lty=2)
  lines(c(xi,xi),c(-10,50),lwd=1,lty=2)
  
  linf<-yl+sdli*qto
  lsup<-yl-sdli*qto
  if(limit=="yes"){
    lines(xl, linf*100, lwd=1,col=ccol[1],lty=2)  ## 4
    lines(xl, lsup*100, lwd=1,col=ccol[1],lty=2)
  } ## 4
  pcho=c(18,15,17,0,1,2,5,6,8)
  for(i in 1:nag){
    agente <- subset(ovifemal, ovifemal[, 1] == not[i])
    points(agente[, 2], agente[, 3]*100,pch=pcho[i],cex=1.3,col=colo[i]) ## colo
  }
  legend(ldx,ldy,cex=1.2,not,pch=pcho[1:nag],col=colo,lty = 3) ## colo
  lines(xl, yl*100, lwd=2,col=ccol[2]) ## 1
}
