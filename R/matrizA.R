#' matrizA
#'
#' no details
#' 
#' @param estadios character vector, all states
#' @param hfeno list, objects about phenology
#' @param Table data frame, temperatures (tmin and tmax)
#' @param steps integer value, between 1 and 48
#' @param modelim integer value, position of the nonlinear model (development rate)
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
matrizA <- function(estadios, hfeno, Table, steps, modelim){
  inmaduros <-  estadios[-(length(estadios)-1):-(length(estadios))]
  maduros   <-  estadios[(length(estadios)-1):(length(estadios))]
  nmax=nrow(Table)
  matriz<-matrix(0,ncol=2*(length(inmaduros)*2+length(maduros)+1),nrow=nrow(Table))
  Table2=cbind(id=1:nrow(Table),Table)
  
  for(K in 1:length(inmaduros))
  {
    #  Desarrollo  # extraccion de funciones y parametros
    
    parametrosc <- hfeno$pdv_dv[[K]]
    funciont <- as.expression(hfeno$fdv_dv[[K]][[3]])
    
    #   Mortalidad # extraccion de funciones y parametros
    
    parametrosm <- hfeno$pmortal[[K]]
    funcionm <- as.expression(hfeno$mortal[[K]][[3]])
    
    RM=apply(Table2,1,RateI,Table2,K,parametrosc,parametrosm,funciont,funcionm,nmax,steps,modelim=modelim) ## procesamiento de tasa de desarrollo y mortalidad por cada temperatura
    matriz[,4*K-3]=RM[1,];matriz[,4*K-2]=RM[2,]
    matriz[,4*K-1]=RM[3,];matriz[,4*K]=RM[4,]
  }
  #  Hembras
  
  parametrosc <- hfeno$pfh_h
  for (i in names(parametrosc)){temp <- parametrosc[i];storage.mode(temp) <- "double";assign(i, temp)}
  formulac <- hfeno$fh_h
  funciont <- as.expression(formulac[[3]])
  RM=apply(Table2,1,RateI,Table2,(K+1),parametrosc,funciont=funciont,nmax=nmax,steps=steps,J=NA,modelim=modelim)
  matriz[,4*(K+1)-3]=RM[1,];matriz[,4*(K+1)-2]=RM[2,]
  
  #  Machos
  
  parametrosc <- hfeno$pfm_m;
  for (i in names(parametrosc)){temp <- parametrosc[i];storage.mode(temp) <- "double";assign(i, temp)}
  formulac <- hfeno$fm_m
  funciont <- as.expression(formulac[[3]])
  RM=apply(Table2,1,RateI,Table2,(K+2),parametrosc,funciont=funciont,nmax=nmax,steps=steps,J=NA,modelim=modelim)
  matriz[,4*(K+1)-1]=RM[1,];matriz[,4*(K+1)]=RM[2,]
  matriz[matriz>1]<-1;matriz[matriz<0]<-0
  
  #  Fecundidad
  
  parametrosc <- hfeno$ptazaeh_h
  for (i in names(parametrosc)){temp <- parametrosc[i];storage.mode(temp) <- "double";assign(i, temp)}
  formulac <- hfeno$ftazaeh_h
  funciont <- as.expression(formulac[[3]])
  RM=apply(Table2,1,RateI,Table2,(K+3),parametrosc,funciont=funciont,nmax=nmax,steps=steps,J=3,modelim=modelim)
  matriz[,4*(K+1)+1]=RM[1,];matriz[,4*(K+1)+2]=RM[2,]
  
  vectorPares = (1:(ncol(matriz)/2))*2
  matrizA = matriz[,vectorPares]
  return(matrizA)
}
