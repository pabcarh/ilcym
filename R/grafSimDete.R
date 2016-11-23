#' grafSimDete
#'
#' no details
#' 
#' @param matrizOut matrix, life parameters estimated by temperature
#' @param estadios character vector, all states
#' @param labx a character
#' @param laby a character
#' @param titulo a character
#' @param legx double, x position of the legend
#' @param legy double, y position of the legend
#' @param corrx1 double, minimun axis valuefor X
#' @param corrx2 double, maximun axis valuefor X
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
grafSimDete<-function(matrizOut, estadios, labx,laby,titulo,legx,legy,corrx1,corrx2){
  numEstadios = length(estadios)
  names1 = colnames(matrizOut)
  names1=names1[1:(length(names1)-2)]
  matrizOut[matrizOut<1]=1
  dat = matrizOut
  cap=rep(10,15)
  tamY = cap^(0:14)
  tamsY = seq(0,30,length.out=15)
  
  x1 = 1:nrow(dat)
  #x1 = corrx
  z1 = log(dat[,numEstadios-1])
  mod1 = lm(z1~x1)
  
  plot(1:nrow(dat),log(dat[,1]),type="n", xlab=labx, ylab=laby,axes=FALSE, xlim = c(corrx1,corrx2))
  #plot(corrx,log(dat[,1]),type="n", xlab=labx, ylab=laby,axes=FALSE)
  grid(NA,15,lwd=1)
  title(titulo)
  axis(1,seq(corrx1,corrx2,round(corrx2/8)))
  axis(2,tamsY,tamY, las=2, cex.axis=0.6)
  abline(mod1, lwd=2)
  
  for(i in 1:(numEstadios-1)){
    points(1:nrow(dat),log(dat[,i]),type="l",col=(i+1))
    #points(corrx,log(dat[,i]),type="l",col=(i+1))
  }
  
  cols = c(2:(i+1),1)
  legend(legx,legy,c(names1,paste("Expon.(",names1[length(names1)],")")),col =cols,lty = 1)
  
}
