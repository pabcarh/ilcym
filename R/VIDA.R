#' VIDA
#'
#' no details
#' 
#' @param vidalife matriz, observed life table (fluctuating)
#' @param Tablelife matriz, counts for the simulated life table
#' @param p matrix, life parameters estimated by temperature
#' @param cuadro matrix, predictions of phenology variables by day
#' @param estadios character vector, all states
#' @param N integer value, number of insects to simulate
#' @param M integer value, number of days
#' @param hfeno list, objects about phenology
#' @param Rs numeric vector, sexual ratios by day
#' @param Rs2 numeric vector, sexual ratios by day with position
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
VIDA<-function(vidalife,Tablelife,p,cuadro,estadios,N,M,hfeno,Rs,Rs2=NULL){
  table_life<-estimado<-res<-list(1)
  NUM <- 3;ni=length(estadios)-2
  tablita<-matrix(NA,6,NUM)
  RES<-matrix(NA,ni+1,NUM)
  RESM<-matrix(NA,ni,NUM)
  
  for(k in 1:NUM){
    estimado [[k]]<-simulacion(cuadro,estadios,N,M,hfeno,Rs)$mat
    table_life[[k]] <- life.table2(estimado[[k]], estadios)$life.table
    tablita[,k]<-PARAlife(N, estadios,table_life[[k]])$parametro
    RES[,k]<-(statistlife(estimado[[k]], estadios)$Survival_Time)[,1]
    RESM[,k]<-as.numeric(as.character((statistlife0(estimado[[k]],estadios)$Mortality)[,1]))
  }
  
  tablita=t(na.omit(t(tablita)))
  MEANlife<-apply(tablita,1,median)
  SElife<-apply(tablita,1,se)
  
  MEANres<-apply(RES[1:ni,],1,mean)
  SEres<-apply(RES[1:ni,],1,se)
  
  MEANresm<-apply(RESM,1,median) ## este vector contiene el porcentaje promedio de mortalidad de cada estado inmaduro
  SEresm<-apply(RESM,1,se) ## contiene los errores estandar de los porcentajes de mortalidad de cada estado
  
  OBS1<-statistlife(vidalife, estadios)$Survival_Time
  OBS2<-statistlife(vidalife,estadios)$Mortality
  
  colnames(p)<-"Observed";t95=abs(qt(0.05/2,2))
  
  P<-c(1)
  for(s in 1:6){P[s]<-t.test(tablita[s,],mu=p[s,],alt="two.side",conf.level = 0.99)$p.value} ## p-value de la prueba de hipotesis de los parametros
  CUADRO1<-data.frame(Simulated=paste(round(MEANlife,3),"(","?",round(t95*SElife,3),")"),round(p,3),P=round(P,4))
  cat("Life-table parameters\n")
  print(CUADRO1)
  
  P<-c(1)
  for(s in 1:ni){P[s]<-t.test(RES[s,],mu=OBS1[s,],alt="two.side")$p.value}
  CUADRO2<-data.frame(Simulated=paste(round(MEANres,3),"(","?",round(t95*SEres,3),")"),Observed=OBS1[1:ni,],P=round(P,4))
  cat("\nDevelopment time (days)\n");rownames(CUADRO2)=estadios[1:ni]
  print(CUADRO2)
  
  P<-c(1)
  for(s in 1:ni){P[s]<- t.test(RESM[s,],mu=as.numeric(as.character(OBS2[s,])),alt="two.side")$p.value}
  CUADRO3<-data.frame(Simulated=paste(round(MEANresm,3),"(","?",round(t95*SEresm,3),")"),Observed=OBS2,P=round(P,4))
  cat("\nMortality (%)\n")
  print(CUADRO3)
}
