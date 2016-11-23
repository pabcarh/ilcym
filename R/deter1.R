#' deter1
#'
#' no details
#' 
#' @param cuadro matrix, predictions of phenology variables by day
#' @param estadios character vector, all states
#' @param hfeno list, objects about phenology
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
deter1<-function(cuadro,estadios,hfeno){
  dia<-366 #tama?o fila
  day<-0:363
  inmad <-  estadios[-(length(estadios)-1):-(length(estadios))]
  maduros   <-  estadios[(length(estadios)-1):(length(estadios))]
  ratest<-Rate<-rest<-list(1)
  for(i in 1:length(inmad))      rest[i]<-list(cuadro[,2*i]/cuadro[,(2*i-1)])
  for(i in 1:(length(inmad)+1))  Rate[i]<-list(cuadro[,(2*i-1)])
  Rate[[i+1]]<-cuadro[,(2*i)]
  ovi<-cuadro[,(2*i)+1]
  riter<-dij<-dacc<-pij<-plij<-list(1)
  for(k in 1:length(inmad))
  {
    pd<-d<-enditer<-pl<-p<-c(0)
    i<-1
    while(i<=100)#length(day))
    {
      d[i]<-Rate[[k]][i]*day[i]
      pd[i]<-distrimodeldeim(d[i],k,hfeno)
      if(i==1)p[i]<-pd[i]
      if(i!=1)p[i]<-pd[i]-pd[i-1]
      ifelse(pd[i]==1,enditer[i]<-0,enditer[i]<-1)
      if(enditer[i]==1) i<-i+1 else break
    }
    for(h in 2:i)ifelse(pd[h]==1,pl[h]<-1,pl[h]<-p[h]/(1-pd[h-1]))
    dij[[k]]<-d
    dacc[[k]]<-pd
    pij[[k]]<-p
    plij[[k]]<-pl
    riter[[k]]<-enditer
  }
  for(kk in 1:length(maduros))
  {
    pd<-d<-enditer<-pl<-p<-c(0)
    i<-1
    while(i<=length(day))
    {
      d[i]<-Rate[[k]][i]*day[i]
      pd[i]<-distrimodeldema(d[i],kk)
      if(i==1)p[i]<-pd[i]
      if(i!=1)p[i]<-pd[i]-pd[i-1]
      ifelse(pd[i]==1,enditer[i]<-0,enditer[i]<-1)
      if(enditer[i]==1) i<-i+1 else break
    }
    for(h in 2:i)ifelse(pd[h]==1,pl[h]<-1,pl[h]<-p[h]/(1-pd[h-1]))
    dij[[kk+length(inmad)]]<-d
    dacc[[kk+length(inmad)]]<-pd
    pij[[kk+length(inmad)]]<-p
    plij[[kk+length(inmad)]]<-pl
    riter[[kk+length(inmad)]]<-enditer
  }
  return(list(dij,dacc,pij,plij,riter,rest,Rate,ovi))
}
