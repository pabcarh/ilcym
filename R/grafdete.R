#' grafdete
#'
#' no details
#' 
#' @param simudeter matrix, life parameters estimated by temperature
#' @param dia integer value, day evaluated
#' @param corrx numeric vector of two values, the axis range for X
#' @param corry numeric vector of two values, the axis range for y#' @param riter list, objects about phenology
#' @param tam integer value, extent of the plot
#' @param lgx double, x position of the legend
#' @param lgy double, y position of the legend
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
grafdete<-function(simudeter,dia,corrx,corry,tam,lgx,lgy)
{
  tablesimu<-matrix(0,nrow(simudeter),ncol(simudeter))
  for(i in 1:nrow(simudeter)) for(j in 1:ncol(simudeter)) tablesimu[i,j]<-ifelse(simudeter[i,j]>1,simudeter[i,j],1)
  par(cex=tam)
  plot(1:dia,log(tablesimu[,1]),type="l",col=1,frame=F,ylim=corry,xlim=corrx)
  for(i in 2:(length(estadios)-1)) lines(1:dia,log(tablesimu[,i]),col=i)
  legend(lgx,lgy,names(simudeter)[-c(length(estadios),(length(estadios)+1))],col =1:i,lty = 1)
  return(list(table=tablesimu))
}