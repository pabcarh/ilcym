#' PlotDevath
#'
#' no details
#'
#' @param matri data frame, counts of insect by days and others variables
#' @param parametros numeric vector, slopes of modeling outputs
#' @param test a character, model selected
#' @param corrx numeric vector of two values, the axis range for X
#' @param lgx double, x position of the legend
#' @param lgy double, y position of the legend
#' @param tam integer value, extent of the plot
#' @param labx a character
#' @param laby a character
#' @param titulo a character
#' @param grises logical; default is FALSE
#' @param intervalo double, scale in days
#' @param exis.val numeric vector, if there are live individuals in a temperature
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
PlotDevath <- function(matri, parametros,test,corrx,lgx,lgy,tam,labx = NULL, laby = NULL,titulo=NULL,grises=FALSE,intervalo=intervalo,exis.val)
{
  
  slope <- parametros[-1:-(length(parametros) - 1)]
  intercepto <- parametros[-(length(parametros))]
  nag <- length(levels(factor(matri[, 1])))
  not <- levels(factor(matri[, 1]))
  p50 <- li <- ls <- rep(0, nag)
  #for (i in 1:nag) agente <- subset(matri, matri[, 1] == not[i])
  agente <- subset(matri, matri[, 1] == not[nag])
  
  graphics::par(cex=tam)
  graphics::plot(agente[, 3], agente[, 5]*100/agente[,6], axes=F,xaxt = "n",xlim =corrx , ylim = c(0, 100),
                 frame = F, ylab = laby, xlab = labx,main=titulo)
  #plot(agente[, 3], agente[, 5]*100/agente[,6], axes=F,xaxt = "n",xlim =corrx , ylim = c(0, 100),
  #     frame = F, ylab = laby, xlab = labx,main=titulo, cex.lab=0.6)
  
  #axis(1, xaxp=c(corrx,5))
  graphics::axis(1, seq(corrx[1],corrx[2],0.5))
  graphics::axis(2,seq(0,100,20),paste(seq(0,100,20),"%",sep=""),las=2)
  #pcho <- 15:28
  pcho <- c(18,15,17,0,1,2,5,6,8) #nuevo cambio de puntos puedo aumentar mas (comprobar)
  if(grises==TRUE)
  {co <- c("gray5","gray10","gray15","gray20","gray25","gray30","gray35","gray40","gray45","gray50","gray55","gray60","gray65","gray70")
  }else{co <- 1:15}
  for (i in 1:nag) { ## bucle por temperatura
    agente <- subset(matri, matri[, 1] == not[i])
    
    if (test == "probit") {
      psp <- seq(0.01, 0.99, 0.01)
      pr <- stats::qnorm(psp)
      p <- stats::pnorm(pr)
      w <- (pr - intercepto[i])/slope
      graphics::lines(w, p*100, col = co[i], lwd = 1.8)
      #p50[i] <- (-intercepto[i])/slope
      if(agente[1,2]==intervalo & agente[1,5]==agente[1,4] & exis.val[i]==0){p50[i]=100}else{p50[i] <- (-intercepto[i])/slope}
    }else{
      ps <- seq(0, 10, 0.001)
      if (test == "cloglog") {
        p <- 1 - exp(-exp((intercepto[i] + slope * ps)))
        #p50[i] <- ( -(-log(-log(0.5))) - intercepto[i])/slope
        if( agente[1,2]==intervalo & agente[1,5]==agente[1,4] & exis.val[i]==0){p50[i]=100}else{p50[i] <- ( -(-log(-log(0.5))) - intercepto[i])/slope}
      }
      if (test == "logit") {
        p <- 1/(1 + exp(-(intercepto[i] + slope * ps)))
        #p50[i] <- ((log((1 - 0.5)/0.5)) - intercepto[i])/slope
        if( agente[1,2]==intervalo & agente[1,5]==agente[1,4] & exis.val[i]==0){p50[i]=100}else{p50[i] <- ((log((1 - 0.5)/0.5)) - intercepto[i])/slope}
      }
      graphics::lines(ps, p*100, col = co[i], lwd = 1.8)
    }
    graphics::points(agente[, 3], agente[, 5]*100/agente[,6], pch = pcho[i], cex = 1.5, col = co[i])
    minx <- min(matri[, 3])
    li[i] <- p50[i] - stats::qt(0.975, length(matri[, 3]) - length(intercepto) -
                                  length(slope)) * (1/slope)/sqrt(length(intercepto) + length(slope))
    ls[i] <- p50[i] + stats::qt(0.975, length(matri[, 3]) - length(intercepto) -
                                  length(slope)) * (1/slope)/sqrt(length(intercepto) + length(slope))
    
    graphics::lines(rbind(c(minx - 0.8, 50), c(p50[i], 50), c(p50[i], -100)), col = "gray50")
    graphics::lines(rbind(c(ls[i], 50), c(li[i], 50)), col = "black", lwd = 2)
    graphics::lines(rbind(c(ls[i], -0.03*100), c(li[i], -0.03*100)), col = "black", lwd = 2)
    graphics::lines(rbind(c(li[i], -0.03*100 - 0.02*100), c(li[i], -0.03*100 + 0.02*100)), col = "black", lwd = 2)
    graphics::lines(rbind(c(ls[i], -0.03*100 - 0.02*100), c(ls[i], -0.03*100 + 0.02*100)), col = "black", lwd = 2)
    graphics::lines(rbind(c(li[i], 0.5*100 - 0.02*100), c(li[i], 0.5*100 + 0.02*100)), col = "black", lwd = 2)
    graphics::lines(rbind(c(ls[i], 0.5*100 - 0.02*100), c(ls[i], 0.5*100 + 0.02*100)), col = "black", lwd = 2)
  }
  graphics::legend(lgx,lgy, not, pch = pcho, col = co, lty = 1.5)#mean(matri[,5])#max(matri[,3])
  mat<-data.frame(x = (as.numeric(not)), y = (1/exp(p50)), Lower = (1/exp(ls)), Upper = (1/exp(li)))
  datao<-data.frame(x = (as.numeric(not)), y = c((1/exp(p50)),(1/exp(ls)),(1/exp(li))))
  #if( agente[1,2]==intervalo & length(p50[p50==100])>=1 ){p50[p50==100]=-Inf}  ## el valor de cero puede afectar algunas funciones cuando se hace el ajuste
  p50[p50==100]=-Inf
  
  ############## DM: Inicio: Lo siguiente se agrego para usarlo en los graficos de Development Time y omitir las T extremas donde no se desarrolla el insecto, de modo que el grafico represente mejor la realidad
  vv=exis.val==1; mat<-mat[vv,]; #p50<-p50[vv]
  vvv=c(vv,vv,vv); datao<-datao[vvv,];
  ############## DM: Fin
  
  salidas <- list(mat=mat,shapMi=datao,median=p50)
  return(salidas)
}
