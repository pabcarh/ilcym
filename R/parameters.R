#' parameters
#'
#' no details
#' 
#' @param N integer value, number of simulations
#' @param estad character vector, all states
#' @param ltb matrix, life table counts
#' @param poli double
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
parameters <- function (N, estad, ltb,poli=1){
  estadios <- estad
  lifetable <- ltb
  s.x <- m.x <- l.x <- rep(0, nrow(lifetable))
  for (i in 1:nrow(lifetable)) l.x[i] <- sum(lifetable[i, 1:(length(estadios))])/N  ## porcentaje de vivos por dias
  # for (i in 1:nrow(lifetable)) if(l.x[i]==0) m.x[i] <- 0 else m.x[i] <- lifetable[,(length(estad)+1)][i]/N/l.x[i]
  for (i in 1:nrow(lifetable)) if(l.x[i]==0) m.x[i] <- 0 else m.x[i] <- poli*lifetable[,(length(estad)+1)][i]/N/l.x[i]
  ## el objeto m.x no es afectado por el efecto poliembrionico pero asumimos la mortalidad del host es la unica que afecta ala mortalidad del parasitoide
  s.x1 <- c(l.x[1], l.x[2])
  s.x2 <- rep(0, nrow(lifetable) - 2)
  for (i in 1:(nrow(lifetable) - 2)) s.x2[i] <- l.x[i + 2]/l.x[i + 1]
  s.x <- c(s.x1, s.x2) ## es el porcentaje de vivos segun la cantidad del dia anterior
  Ro <- sum(l.x * m.x)  ## este parametro depende del procentaje de vivos de los inmaduros
  Day <- 0:(nrow(lifetable)-1)
  Tg <- sum(Day * l.x * m.x)/Ro
  GRR <- sum(m.x)
  r <- log(Ro)/Tg
  euler <- sum(l.x * m.x * exp(-r * (Day + 1)))
  r1 <- r - 1e-06
  for (k in 1:100) {
    euler1 <- sum(l.x * m.x * exp(-r1 * (Day + 1)))
    r2 <- ((1 - euler) * (r1 - r)/(euler1 - euler)) + r
    r1 <- r2
  }
  rf <- r1
  Dt <- log(2)/rf
  lamda <- exp(rf)
  T <- log(Ro)/rf
  parameters <- t(data.frame(r = rf, Ro = Ro, GRR = GRR, T = T, lambda = lamda, Dt = Dt))
  colnames(parameters) <- "Parameters"
  print(parameters)
  if(poli>1)
  {
    matrizOut=data.frame(ltb[,1:(length(estadios)-2)],poli*ltb[,(ncol(ltb)-2):ncol(ltb)])
    nombres=colnames(ltb);nombres2=c(paste(nombres[1:(length(estadios)-2)],"Host",sep="_"),paste(nombres[(ncol(ltb)-2):ncol(ltb)],"Paras",sep="_"))
    colnames(matrizOut)=nombres2
  }else{matrizOut=ltb}
  return(list(ltb=matrizOut,m.x=m.x,s.x=s.x,l.x=as.matrix(l.x),parametro=parameters))
}
