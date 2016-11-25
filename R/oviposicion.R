#' oviposicion
#'
#' no details
#'
#' @param opc integer value (1: life table data, 2: Cohort data)
#' @param datos data frame
#' @param colum integer value
#' @param mediana numeric vector
#' @param estadios character vector, all states
#' @param tp numeric vector of the constant temperature
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
oviposicion<-function(opc,datos,colum,mediana,estadios,tp=NULL)  ## cambio
{
  if(opc==1){## ultimo cambio
    data<-datos ## ultimo cambio
    not <- levels(factor(data[,1]))
    agente <- list(1)
    for (i in 1:length(not)) agente[i] <- list(subset(data, data[, 1] == not[i]))
    datos<-agente
    tp<-as.numeric(not)
    ovi<-sumf<-dias<-frecf<-list(1)
    for(i in 1:length(datos))
    {
      if(dim(datos[[i]])[1] == 1 & is.na(datos[[1]][1,colum]) ){ovi[[i]]=NA;sumf[[i]]=0;dias[[i]]=0;frecf[[i]]=0}else{
        ovi[i]<-list(t(datos[[i]][,colum:ncol(datos[[i]])])) ## seleccionan todas las columnas 
        sumf[i]<-list(apply(ovi[[i]],1,media))  ## promedio de huevos por cada dia
        dias[[i]]<-1:nrow(ovi[[i]])
        frecf[[i]]=cumsum(sumf[[i]]/sum(sumf[[i]],na.rm=T)) ## frecuencia relativa acumulada promedio
      }
    }
  }else
  {
    ovi<-sumf<-dias<-frecf<-list(1)
    for(i in 1:length(datos))
    {
      if(ncol(datos[[i]]) == length(estadios)){ovi[[i]]=NA;sumf[[i]]=0;dias[[i]]=0;frecf[[i]]=0}else{
        ovi[i]<-list(t(datos[[i]][,colum:ncol(datos[[i]])])) ## seleccionan todas las columnas 
        sumf[i]<-list(apply(ovi[[i]],1,media))  ## promedio de huevos por cada dia
        dias[[i]]<-1:nrow(ovi[[i]])
        frecf[[i]]=cumsum(sumf[[i]]/sum(sumf[[i]],na.rm=T)) ## frecuencia relativa acumulada promedio
      }
    }
  }
  female<-cbind(T=tp[1],x=dias[[1]]/exp(mediana[1]),y=frecf[[1]])
  for(i in 2:length(datos))female<-rbind(female,cbind(T=tp[i],x=dias[[i]]/exp(mediana[i]),y=frecf[[i]]))
  rownames(female)=1:nrow(female)
  female<-data.frame(female)
  female[is.nan(female[,2]),2]<-0
  return(list(femal=female))
}
