#' Ratios
#'
#' no details
#' 
#' @param Table data frame, temperatures (tmin and tmax)
#' @param modelim integer value, position of the nonlinear model (development rate for inmature states)
#' @param modelm integer value, position of the nonlinear model (development rate for adult states)
#' @param estadios character vector, all states
#' @param xi double, estimated median from female senescence
#' @param steps integer value, between 1 and 48
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
Ratios<-function(Table,modelim,modelm,estadios,xi, steps){
  inmaduros <-  estadios[-(length(estadios)-1):-(length(estadios))]
  maduros   <-  estadios[(length(estadios)-1):(length(estadios))]
  matriz <- matrix(0,ncol=(length(inmaduros)*2+length(maduros)+1),nrow=nrow(Table))
  K<-1
  while(K<=length(inmaduros)){
    #  Desarrollo  # extraccion de funciones y parametros
    parametrosc <- hfeno$pdv_dv[[K]]
    parametrosc<-as.list(parametrosc)
    for (i in names(parametrosc)){
      temp <- parametrosc[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    formulac <- hfeno$fdv_dv[[K]]
    forexc <- formulac[[length(formulac)]]
    funcionc <- as.expression(forexc)
    
    #   Mortalidad # extraccion de funciones y parametros
    parametrosm <- hfeno$pmortal[[K]]
    parametrosm<-as.list(parametrosm)
    for (i in names(parametrosm)){
      temp <- parametrosm[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    formulam <- hfeno$mortal[[K]]
    forexm <- formulam[[length(formulam)]]
    funcionm <- as.expression(forexm)
    M<-R<-M2<-R2<-rep(0,length(Table[,1]))
    
    Ratesumme<-Rateacum<-Mortalitysumme<-Mortacum<-rep(0,length(Table[,1]))
    Ratetot1<-Mortality1<-Ratetot2<-Mortality2<-0
    
    for(i in 1:length(Table[,1])){#recorre cada valor de la tabla d temperaturas
      
      M<-(Table[,1][i]+Table[,2][i])/2
      R<-(Table[,2][i]-Table[,1][i])/2
      ifelse(i!=length(Table[,1]),M2<-(Table[,2][i]+Table[,1][i+1])/2,M2<-(Table[,2][i]+Table[,1][1])/2)
      ifelse(i!=length(Table[,1]),R2<-(Table[,2][i]-Table[,1][i+1])/2,R2<-(Table[,2][i]-Table[,1][1])/2)
      
      Ratetot1<-Mortality1<-Ratetot2<-Mortality2<-0
      #M <- round(M,1)
      #R <- round(R,1)
      #M2 <- round(M2,1)
      #R2 <- round(R2,1)
      
      
      Ratetot1<-Mortality1<-Ratetot2<-Mortality2<-0
      for(j in 0:(steps-1)){#por defecto 48 pasos
        T<-pi/steps*(j+.5)
        ## eval1 ##
        DEGC<-R*cos(T)+M
        if(modelim[K]==1 || modelim[K]==2 || modelim[K]==3 ||
           modelim[K]==4 || modelim[K]==5 || modelim[K]==6 ||
           modelim[K]==7 || modelim[K]==8 || modelim[K]==9 ||
           modelim[K]==10){
          DEGK<-DEGC+273.15
          x<-DEGK
        }else{ x<-DEGC }
        
        Rate<-eval(funcionc)
        Ratetot1<-Ratetot1+Rate
        x<-DEGC
        Mort<-eval(funcionm)
        if(Mort > 1){
          Mort = 1
        }else{
          Mort = Mort
        }
        Mortality1 <- Mortality1+ Mort
        
        ## eval2 ##
        T<-pi/steps*(j+.5)#borrar
        DEGC<-R2*cos(T)+M2
        if(modelim[K]==1 || modelim[K]==2 || modelim[K]==3 ||
           modelim[K]==4 || modelim[K]==5 || modelim[K]==6 ||
           modelim[K]==7  || modelim[K]==8 || modelim[K]==9 ||
           modelim[K]==10){
          DEGK<-DEGC+273.15
          x<-DEGK
        }else{x<-DEGC}
        
        Rate<-eval(funcionc)
        Ratetot2<-Ratetot2+Rate
        x<-DEGC
        Mort<-eval(funcionm)
        if(Mort > 1){
          Mort = 1
        }else{
          Mort = Mort
        }
        Mortality2 <- Mortality2+ Mort
        Mortality2
        
      }
      
      Ratetot1 <- Ratetot1/steps
      Ratetot2 <- Ratetot2/steps
      Ratetot<-(Ratetot1+Ratetot2)/2
      matriz[i,2*K-1] <- Ratetot
      
      Mortality1 <- Mortality1/steps
      Mortality2 <- Mortality2/steps
      Mortality<- (Mortality1+Mortality2)/2
      matriz[i,2*K] <- Mortality
      
    }
    K<-K+1
  }
  for(J in 1:(length(maduros)+1)){
    if(J==1)parametrosc <- hfeno$pfh_h
    if(J==2)parametrosc <- hfeno$pfm_m
    if(J==3)parametrosc <- hfeno$ptazaeh_h
    parametrosc<-as.list(parametrosc)
    for (i in names(parametrosc)){
      temp <- parametrosc[[i]]
      storage.mode(temp) <- "double"
      assign(i, temp)
    }
    if(J==1) formulac <- hfeno$fh_h
    if(J==2) formulac <- hfeno$fm_m
    if(J==3) formulac <- hfeno$ftazaeh_h
    forexc <- formulac[[length(formulac)]]
    funcionc <- as.expression(forexc)
    M<-R<-M2<-R2<-rep(0,length(Table[,1]))
    
    Ratesumme<-Rateacum<-rep(0,length(Table[,1]))
    ratetot1<-ratetot2<-ratetot<-0
    fectot1<-fectot2<-fectot<-fectotal<-0
    
    for(i in 1:length(Table[,1])){
      Ratetot1<-Ratetot2<-Ratetot<-0
      if(J==3){
        fectot1<-fectot2<-fectot<-0
      }
      M<-(Table[,1][i]+Table[,2][i])/2
      R<-(Table[,2][i]-Table[,1][i])/2
      ifelse(i!=length(Table[,1]),M2<-(Table[,2][i]+Table[,1][i+1])/2,M2<-(Table[,2][i]+Table[,1][1])/2)
      ifelse(i!=length(Table[,1]),R2<-(Table[,2][i]-Table[,1][i+1])/2,R2<-(Table[,2][i]-Table[,1][1])/2)
      
      for(j in 0:(steps-1)){
        T<-pi/steps*(j+.5)
        #### eval1 ###
        DEGC<-R*cos(T)+M
        x<-DEGC
        if(J==1 || J==2){
          if(modelm[J]==1 || modelm[J]==2 || modelm[J]==3 ||
             modelm[J]==4 || modelm[J]==5 || modelm[J]==6 ||
             modelm[J]==7 || modelm[J]==8 || modelm[J]==9 ||
             modelm[J]==10){
            DEGK<-DEGC+273.15
            x<-DEGK
          }else{x<-DEGC}
        }
        
        Rate<-eval(funcionc)
        Ratetot1<-Ratetot1+Rate
        if(J==3){
          x<-DEGC
          fec1<-eval(funcionc)
          if(fec1>0){ fec1=fec1}else{ fec1=0}
          fectot1<-fectot1+fec1
        }
        
        #### eval2 ###
        DEGC<-R2*cos(T)+M2
        if(J==1 || J==2){
          if(modelm[J]==1 || modelm[J]==2 || modelm[J]==3 ||
             modelm[J]==4 || modelm[J]==5 || modelm[J]==6 ||
             modelm[J]==7 || modelm[J]==8 || modelm[J]==9 ||
             modelm[J]==10){
            DEGK<-DEGC+273.15
            x<-DEGK
          }else{x<-DEGC}
        }
        #if(J==3) x<-DEGC
        Rate<-eval(funcionc)
        Ratetot2<-Ratetot2+Rate
        if(J==3){
          x<-DEGC
          fec2<-eval(funcionc)
          if(fec2>0){
            fec2=fec2
          }else{
            fec2=0
          }
          fectot2<-fectot2+fec2
        }
      }
      Ratetot1 <- Ratetot1/steps
      Ratetot2 <- Ratetot2/steps
      Ratetot<-(Ratetot1+Ratetot2)/2
      if(J==1){
        matriz[i, length(inmaduros)*2+1] <- Ratetot
      }
      if(J==2){
        matriz[i, length(inmaduros)*2+2] <- Ratetot
      }
      if(J==3){
        fectot1 = fectot1 / steps
        fectot2 = fectot2 / steps
        fectot = (fectot1 + fectot2) / 2
        if(fectot < 0){
          fectot = 0
        }
        matriz[i, length(inmaduros)*2+3] <- fectot
      }
    }
  }
  
  T<-1/matriz[,1]+1/matriz[,3]+1/matriz[,5]+xi/matriz[,7]
  Ro<-matriz[,9]*(1-matriz[,2])*(1-matriz[,4])*(1-matriz[,6])/2
  rm<-log(Ro)/T
  J<-exp(rm)
  Dt<-log(2)/rm
  E<-(1-matriz[,2])^matriz[,1]
  L<-(1-matriz[,4])^matriz[,3]
  P<-(1-matriz[,6])^matriz[,5]
  ERI<-(1-ifelse(length(E[E==0])==0,0,length(E[E==0])/12))*
    (1-ifelse(length(L[L==0])==0,0,length(L[L==0])/12))*
    (1-ifelse(length(P[P==0])==0,0,length(P[P==0])/12))
  
  return(list(MATRIZ=data.frame(matriz),ERI=ERI,J=J))
}
