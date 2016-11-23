#' simultemp
#'
#' no details
#' 
#' @param N integer value, number of simulations
#' @param sexratio double, between 0 to 1, sex ratio is defined in the ilcym project
#' @param isFixed logical; default is FALSE
#' @param temp numeric vector, temperatures
#' @param reps vector of integer values, repetitions for each constant temperature
#' @param xi double, estimated median from female senescence
#' @param steps integer value, between 1 and 48
#' @param poli double
#' @param params list
#' @param .. Phenology objetc or arrays
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
simultemp<-function(N,sexratio,isFixed,temp,reps,xi,steps, poli, params, ..){
  n1=length(temp)
  hfeno=params$hfeno
  estadios<-params$estadios
  modelim<-params$modelim
  modelm<-params$modelm  
  modelim=c(modelim,modelm)
  pars=data.frame(Temper=NA,r=NA,Ro=NA,GRR=NA,T=NA,lambda=NA,Dt=NA)
  
  for(i in 1:n1)
  {
    Table=data.frame(V1=rep(temp[i],365),V2=rep(temp[i],365))
    
    if(isFixed){
      Rs=rep(sexratio,365)
    }else{
      paramR=params$paramR
      if(is.na(paramR)){
        sexratio <- params$ff
        Rs=rep(sexratio,365)
      }else{
        Table2=cbind(id=1:nrow(Table),Table)
        #paramR=params$paramR
        ff=params$ff
        Rs=apply(Table2,1,RateR,Table2,nmax=nrow(Table),steps=steps,ff,paramR)
      }
    }
    
    
    for(j in 1:reps[i])
    {
      #cuadro<-Ratios(Table,modelim,modelme,estadios, xi,steps)$MATRIZ
      
      cuadro <- matrizA(estadios, hfeno, Table, steps, modelim=modelim)
      #insectos<-simulacion(cuadro,estadios,N,M,hfeno,sexratio)
      
      
      ageclases <- AgeClases(cuadro, estadios, hfeno)
      estadiosAges <- ageclases$estadiosAges
      oviFreq <- ageclases$oviFreq
      ageClassMatriz <- ageclases$ageClassMatriz
      Day=1
      severalYear=FALSE
      Steps=364
      simu2<-simulacionUnaGeneracion(Day,estadiosAges, oviFreq, ageClassMatriz, hfeno, cuadro, estadios,N,sexratio,Rs,Steps)
      matrizOut<-simu2$matrizOut
      
      pravida<-parameters(N, estadios, matrizOut,poli )
      parametros<-pravida$parametro
      pars=rbind(pars,c(Temper=temp[i],t(parametros)))
    }
  }
  pars=pars[-1,]
  return(pars)
}
