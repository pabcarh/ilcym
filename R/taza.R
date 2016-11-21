#' taza
#'
#' no details
#' 
#' @param opc integer value (1: life table data, 2: Cohort data)
#' @param datos data frame
#' @param tp numeric vector of the constant temperature
#' @param colum integer value
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
taza<-function(opc,datos,tp,colum){## Todo esto funciona si la muerte esta codificada como "NA"
  if(opc==1){## hace que se creen sub grupos por cada temperatura
    data2<-datos
    not <- levels(factor(data2[,1]))
    agente <- list(1)
    for (i in 1:length(not)) agente[i] <- list(subset(data2, data2[, 1] == not[i]))
    datos<-agente
    tp<-as.numeric(not)
  }
  ovi<-list(1)
  if(opc==1){colum=2;for(i in 1:length(datos)) ovi[i]<-list(t(datos[[i]][,colum:ncol(datos[[i]])]))}else{
    for(i in 1:length(datos)) if(ncol(datos[[i]])>colum){ovi[i]<-list(t(datos[[i]][,(colum+1):ncol(datos[[i]])]))}else{ovi[i]<-list(t(matrix(NA,nrow(datos[[i]]),2)))}}
  
  sftab<-list(1)
  for(k in 1:length(datos)) sftab[k]<-list(apply(ovi[[k]],2,summ)) ## me deberia dar el numero total de huevos por cada hembra por temperatura
  novi<-sdovi<-meaovi<-c(1)
  todos<-list(1)
  
  for(i in 1:length(tp)){
    novi[i]<-length(na.omit(sftab[[i]]))
    if(novi[i]==0){todos[[i]]<-NA}else{todos[[i]]<-na.omit(sftab[[i]])}
    meaovi[i]<-mean(na.omit(sftab[[i]]))
    sdovi[i]<-sd(na.omit(sftab[[i]]))
  }
  meaovi[is.nan(meaovi)]=0
  
  todo<-cbind(tp[1],todos[[1]])
  for(i in 2:length(tp)) todo<-rbind(todo,cbind(tp[i],todos[[i]]))
  todo<-data.frame(todo)
  colnames(todo)<-c("x","y")
  tablaovi<-data.frame(x=tp,y=meaovi,sd=sdovi,n=novi)
  female<-data.frame(x=tp,y=meaovi)
  salida<-list(femal=tablaovi,todo=todo)
  return(salida)
}
