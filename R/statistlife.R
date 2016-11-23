#' statistlife
#'
#' no details
#' 
#' @param data matriz, observed life table (fluctuating)
#' @param estad character vector, all states
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
statistlife<-function (data, estad){
  estadios <- estad
  esadult <- estadios[(length(estadios) - 1):length(estadios)]
  datat <- as.matrix(data)
  estadoa <- as.list(as.data.frame((datat)))
  estadob <- as.list(rep(0, ncol(data)))
  num <- matrix(as.numeric(datat), nrow(datat), ncol(datat))
  numc <- num
  for (i in 1:nrow(numc)) for (j in 1:ncol(numc)) if (is.na(numc[i,j]))
    numc[i, j] <- 0
  newegg <- matrix(apply(numc, 1, sum), nrow(numc))
  est1 <- subset(estadios, estadios != esadult[2])
  nopupa1 <- nopupa2 <- rep(0, (length(est1) - 2))
  estados3 <- huevo <- matrix(rep(0, length(nopupa1) * ncol(data)),
                              nrow = ncol(data), ncol = length(nopupa1))
  for (is in 1:(length(est1) - 1)) {
    est2 <- est1[is + 1]
    ifelse(est2 == esadult[1], caso <- "caso2", caso <- "caso1")
    if (caso == "caso1") {
      for (i in 1:ncol(data)) {
        estadoa[[i]] <- as.matrix(estadoa[[i]])
        estadob[[i]] <- length(subset(estadoa[[i]], estadoa[[i]] ==
                                        est2))
      }
      estados3[, is] <- as.numeric(as.character(as.matrix(estadob)))
      for (i in 1:length(estados3[, 1])) for (j in 1:length(nopupa1)) if (estados3[i,j] != 0)
        estados3[i, j] <- 1
      nopupa1[is] <- sum(estados3[, is])
      huevo[, is] <- Estado(t(data), estad, est1[is])
      for (i in 1:length(estados3[, 1])) for (j in 1:length(nopupa1)) ifelse(estados3[i,j] == 1, huevo[i, j] <- huevo[i, j], huevo[i,
                                                                                                                                 j] <- 0)
      for (i in 1:length(estados3[, 1])) for (j in 1:length(nopupa1)) ifelse(huevo[i,
                                                                                   j] == 0, huevo[i, j] <- NA, huevo[i, j] <- huevo[i,
                                                                                                                                    j])
      nopupa2[is] <- mean(huevo[, is], na.rm = T)
    }
  }
  if (caso == "caso2") {
    for (i in 1:ncol(data)) {
      estadoa[[i]] <- as.matrix(estadoa[[i]])
      estadob[[i]] <- length(subset(estadoa[[i]], estadoa[[i]] == esadult[2]))
    }
    estados2 <- as.numeric(as.character(as.matrix(estadob)))
    for (i in 1:length(estados2)) if (estados2[i] != 0)
      estados2[i] <- 1
    male <- sum(estados2)
    num <- matrix(as.numeric(datat), nrow(datat), ncol(datat))
    a <- as.list(as.data.frame((num)))
    b <- as.list(rep(0, ncol(data)))
    for (i in 1:ncol(data)) b[[i]] <- sum(as.matrix(table(a[[i]])))
    female1 <- as.numeric(as.character(as.matrix(b)))
    for (i in 1:length(female1)) if (female1[i] != 0)
      female1[i] <- 1
    female <- sum(female1)
    for (i in 1:ncol(data)) {
      estadoa[[i]] <- as.matrix(estadoa[[i]])
      estadob[[i]] <- length(subset(estadoa[[i]], estadoa[[i]] == (estad[is + 3])))
    }
    estados5 <- as.numeric(as.character(as.matrix(estadob)))
    for (i in 1:length(estados5)) if (estados5[i] != 0)
      estados5[i] <- 1
    dead <- sum(estados5)
    pupa1 <- male + female
    capullo <- Estado(t(data), estad, est1[is])
    capullo2 <- rep(0, length(female1))
    for (i in 1:length(female1)) ifelse(estados2[i] == 0 &
                                          female1[i] == 0, capullo2[i] <- 0, capullo2[i] <- 1)
    for (i in 1:length(female1)) ifelse(capullo2[i] == 1,
                                        capullo[i] <- capullo[i], capullo[i] <- 0)
    for (i in 1:length(female1)) ifelse(capullo[i] == 0,
                                        capullo[i] <- NA, capullo[i] <- capullo[i])
    pupa <- mean(capullo, na.rm = T)
  }
  f <- c(ncol(data), female/(female + male), male, female,
         ncol(data) - male - female, male + female + dead, sum(newegg)/female)
  observations <- data.frame(f)
  rownames(observations) <- c("Number of insect :", "Sex ratio :",
                              "Males :", "Females :", "Immature death :", "", "Eggs/Females :")
  colnames(observations) <- c("Observations")
  todo <- rep(0, length(female1))
  for (i in 1:length(female1)) todo[i] <- sum(capullo[i], huevo[i,])
  todos <- mean(todo, na.rm = T)
  time <- c(round(nopupa2, 3), round(pupa, 3), round(todos,
                                                     3))
  insect <- c(nopupa1, pupa1, "")
  d1 <- rep(0, (length(nopupa2) - 1))
  for (i in 1:(length(nopupa2) - 1)) d1[i] <- nopupa1[i + 1] *100/nopupa1[i]
  d2 <- pupa1 * 100/nopupa1[-1:-(length(nopupa2) - 1)]
  percente <- c(round((nopupa1[1]/(ncol(data))) * 100, 2),
                round(d1, 2), round(d2, 2))/100
  percent <- c(paste(percente * 100, "%"), "")
  acc <- rep(1, length(nopupa1) + 1)
  acc[1] <- percente[1]
  for (i in 2:(length(nopupa1) + 1)) acc[i] <- acc[i - 1] *percente[i]
  accume <- acc
  accum <- c(paste(round(100 * accume, 3), "%"), "")
  survival <- data.frame(cbind(time))
  rownames(survival) <- c(est1[-(length(est1))], "Total")
  colnames(survival) <- c("Time")
  inmor <- rep(0, (length(nopupa1) - 1))
  for (i in 1:(length(nopupa1) - 1)) inmor[i] <- nopupa1[i] -nopupa1[i + 1]
  insectmor <- c(ncol(data) - nopupa1[1], inmor, nopupa1[-1:-(length(nopupa2) -1)] - pupa1)
  pmor <- rep(0, (length(nopupa1) - 1))
  for (i in 1:(length(nopupa1) - 1)) pmor[i] <- (nopupa1[i] -nopupa1[i + 1]) * 100/nopupa1[i]
  ##
  contar1<-function(vec,est){n=length(vec[vec==est]);return(n)}
  contar2<-function(vec){n=length(vec[!is.na(vec)]);return(n)}
  
  myd<-function(vec,estad,est)
  {
    pos=(1:length(estad))[estad==est]
    vec=data.frame(vec,1)
    v1=rownames(vec[vec[,1]==est,])
    if(est!=estad[length(estadios) - 2]){ 
      if(length(v1) != 0){
        n=as.numeric(v1[length(v1)]) 
        if(vec[n+1,1]==estad[pos+1]){mort=0;dsr=1}else{mort=1;dsr=0}
        salida=c(Mort=mort,Desar=dsr,di=n+1)
      }else(salida=c(Mort=NA,Desar=NA,di=NA))
      return(salida)  
    }else{
      if(length(v1) != 0){
        n=as.numeric(v1[length(v1)]) 
        if(vec[n+1,1]==estad[pos+2] || is.numeric(estad[pos+1])){mort=0;dsr=1}else{mort=1;dsr=0}
        salida=c(Mort=mort,Desar=dsr,di=n+1)
      }else(salida=c(Mort=NA,Desar=NA,di=NA))
      return(salida)  
    }    
  }
  
  datos=data
  estadinm=estadios[1:(length(estadios) - 2)]
  
  mor=c(1:length(estadinm))
  
  for(i in 1:length(estadinm))
  {
    o1=apply(datos,2,myd,estadios,estadinm[i]);o1=data.frame(t(o1))
    o2=na.omit(o1)
    mort=aggregate(Mort ~ di,data=o2,sum);v1=sum(mort[,2])
    desa=aggregate(Desar ~ di,data=o2,sum);v2=sum(desa[,2])
    mor[i]=v1/sum(v1,v2)
  }
  
  ##
  permore <- c(round(((ncol(data) - nopupa1[1])/(ncol(data))) *
                       100, 2), round(pmor, 2), round((nopupa1[-1:-(length(nopupa2) -
                                                                      1)] - pupa1) * 100/nopupa1[-1:-(length(nopupa2) - 1)],2))
  permor <- paste(round(mor, 3))
  #permor <- paste(round(permore/100, 3))
  acc <- c(paste(round((1 - accume) * 100, 3), "%"))
  mortality <- data.frame(cbind(permor))
  rownames(mortality) <- est1[-(length(est1))]
  colnames(mortality) <- c("Percent")
  #print(observations)
  cat("\n")
  #print(survival)
  cat("\n")
  #print(mortality)
  salida <- list(Survival_Time = survival, Mortality = mortality)
  return(salida)
}
