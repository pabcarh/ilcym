#' mortality
#'
#' no details
#' 
#' @param opc integer value (1: life table data, 2: Cohort data)
#' @param datos data frame
#' @param estadios character vector, all states
#' @param est a character
#' @param tp numeric vector of the constant temperature
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
mortality<-function (opc, datos, estadios, est, tp){
  if(opc == 1){
    tip<-1:length(estadios)
    for (is in 1:(length(estadios))) if (estadios[is] == est) data<-datos[[is]]
    nag <- length(levels(factor(data[, 1])))
    not <- levels(factor(data[, 1]))
    agente <- list(1)
    mortalidad <- c(1)
    for (i in 1:nag){
      agente[i] <- list(subset(data, data[, 1] == not[i]))
      if(agente[[i]][1,3] == agente[[i]][1,4]){mortalidad[i]<-1}else{
        mortalidad[i]<-round((mean(agente[[i]][,3])-sum(agente[[i]][,4]))/mean(agente[[i]][,3]),3)}
    }
  }
  if(opc == 2){
    N=sapply(datos,nrow)
    cum.mort<-list(1)
    num.est=(1:(length(estadios)))[estadios==est]
    num.temp=1:num.est
    for(k1 in num.temp)
    {
      temporal=mort.select(datos, estadios, k1)
      cum.mort[[k1]]=sapply(temporal,length)
    }
    cum.mort=do.call(rbind.data.frame, cum.mort)
    if(num.est>1)
    {
      dcum.mort=apply(cum.mort,2,descum)
      mortalidad=dcum.mort[num.est,]/(N-cum.mort[num.est-1,])
      mortalidad[is.na(mortalidad)]=1
      names(mortalidad)=tp;mortalidad=(t(mortalidad))[,1]
    }else{
      mortalidad=cum.mort[num.est,]/(N)
      mortalidad[is.na(mortalidad)]=1
      names(mortalidad)=tp;mortalidad=(t(mortalidad))[,1]
    }
  }
  data <- data.frame(x = tp, y = mortalidad)
  dat <- data.frame(T = tp, Mortality = mortalidad)
  return(list(datam = data))
}
