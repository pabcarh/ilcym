#' simulacionUnaGeneracion
#'
#' no details
#' 
#' @param Day integer value, day evaluated
#' @param estadiosAges list
#' @param oviFreq matrix, predictions of oviposition
#' @param ageClassMatriz matrix
#' @param hfeno list, objects about phenology
#' @param cuadro matrix, predictions of phenology variables by day
#' @param estadios character vector, all states
#' @param numIni integer value
#' @param sexratio double, between 0 to 1, sex ratio is defined in the ilcym project
#' @param Rs numeric vector, sexual ratios by day
#' @param Steps integer value, between 1 and 48
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
simulacionUnaGeneracion <- function(Day, estadiosAges, oviFreq, ageClassMatriz, hfeno, cuadro, estadios,numIni,sexratio,Rs,Steps){
  Step=0
  
  matrizA = cuadro
  
  sizeMatriz = nrow(matrizA)
  
  inmad <-  estadios[-(length(estadios)-1):-(length(estadios))]
  
  numInmaduros = length(inmad)
  
  maduros   <-  estadios[(length(estadios)-1):(length(estadios))]
  
  numMaduros = length(maduros)
  
  NEWIND = 0
  
  dias<-254#porq?
  
  
  
  #matrizOut = matrix(NA, nrow=(nrow(matrizA)-1), (numInmaduros+numMaduros+1))
  
  matrizOut = matrix(NA, nrow=nrow(matrizA), (numInmaduros + numMaduros + 1))
  
  stagesN = c(rep(0,dias))
  
  stagesN1 = c(numIni, rep(0,(dias-1)))
  
  femalesN = c(rep(NA,dias))
  
  eggsN = stagesN1
  
  stagesDev = c(rep(0,dias))
  
  
  
  ## vectorN <- list(1)
  
  ## vectorNini <- list(1)
  
  
  
  #if(!severalYear){
  
  vectorN <- list(1)
  
  vectorNini <- list(1)
  
  for(k in 1:(numInmaduros + numMaduros)){
    
    v1 = c(100,rep(0,(ageClassMatriz[1,k]-1)))
    
    v2 = c(0,rep(0,(ageClassMatriz[1,k]-1)))
    
    v3 = c(0,rep(0,(ageClassMatriz[1,k]-1)))#
    
    
    
    if(k == 1){
      
      vectorN[[k]] = list(v1,v2,v3)#
      
      vectorNini[[k]] = v1
      
    }else{
      
      vectorN[[k]] = list(v2,v2,v3)#
      
      vectorNini[[k]] = v2
      
    }
    
  }
  
  
  
  matrizOut[1,] = c(numIni,rep(0,(numInmaduros + numMaduros)))# numero inicial de individuos
  
  #}
  
  ## else{
  
  ##     vectorNini <- list(1)
  
  ##     for(k in 1:(numInmaduros + numMaduros)){
  
  ##         vectorNini[[k]] <- vectorN[[k]][[1]]
  
  ##         matrizOut[1,k] = round(numInd[[k]])
  
  ##     }
  
  ## }
  
  
  
  #while(Day < (sizeMatriz-1)){
  
  while(Step < 363){
    
    for(k in 1:(numInmaduros + numMaduros)){
      
      #Ageclass = ageClassMatriz[Day,k]
      
      Ageclass = ageClassMatriz[Day+1,k]
      
      
      
      if(k <= numInmaduros){
        
        stagesN = c(rep(0,Ageclass))
        
      }
      
      
      
      if(k == (numInmaduros+1)){
        
        femalesN = c(rep(0,Ageclass))
        
      }
      
      
      
      matrizAges = estadiosAges[[k]]
      
      NEWIND = 0
      
      
      
      E <- c(rep(0, Ageclass))
      
      EL <- c(rep(0, Ageclass))
      
      Eage <- c(rep(0, Ageclass))
      
      pacc <- c(rep(0, Ageclass))
      
      p1 <- c(rep(0, Ageclass))
      
      p2 <- c(rep(0, Ageclass))
      
      
      
      if(k <= numInmaduros){
        
        s1 <- c(rep(NA, Ageclass))
        
      }
      
      
      
      if(k == numInmaduros+1){
        
        vectorN[[k]][[2]] = vectorN[[k]][[3]]#
        
      }
      
      if(k == numInmaduros+2){
        
        vectorN[[k]][[2]] = vectorN[[k]][[3]]#
        
      }
      
      
      
      for(Ii in 1:(Ageclass-1)){
        
        Eage[Ii] = matrizAges[Day+1, Ii+1]
        
        
        
        if(k <= numInmaduros){
          
          survivalk = (1 - matrizA[Day,k*2]) ^ matrizA[Day,k*2-1]
          
        }
        
        
        
        if(k <= numInmaduros){
          
          if(Eage[Ii] > 1){
            
            s1[Ii] = 1
            
          }else{
            
            s1[Ii] = survivalk
            
          }
          
        }
        
        
        
        if(k <= numInmaduros){
          
          pacc[Ii] = distrimodeldeim(Eage[Ii], k, hfeno)
          
        }else{
          
          pacc[Ii] = distrimodeldema(Eage[Ii], k - numInmaduros, hfeno)
          
        }
        
        
        
        if(Ii == 1){
          
          p1[Ii] = pacc[Ii] - 0
          
        }else{
          
          p1[Ii] = pacc[Ii] - pacc[Ii - 1]
          
        }
        
        
        
        if(p1[Ii] == 0){
          
          p2[Ii] = 0
          
        }else{
          
          if(Ii == 1){
            
            p2[Ii] = p1[Ii] / (1 - 0)
            
          }else{
            
            p2[Ii] = p1[Ii] / (1 - pacc[Ii - 1])
            
          }
          
        }
        
        
        
        if(Day==1){
          
          E[Ii] = vectorNini[[k]][Ii]
          
        }else{
          
          E[Ii] = vectorN[[k]][[1]][Ii]#modifique el 1 por el 2
          
        }
        
        
        
        
        
        if(k <= numInmaduros){
          
          EL[Ii] = E[Ii] * p2[Ii] * s1[Ii]
          
        }else{
          
          EL[Ii] = E[Ii] * p2[Ii]
          
        }
        
        
        
        if(k <= numInmaduros){
          
          E1 = E[Ii] * s1[Ii] - EL[Ii]
          
        }else{
          
          E1 = E[Ii] - EL[Ii]
          
        }
        
        
        
        stagesN[Ii+1] = E1
        
        
        
        if(k <= numInmaduros){
          
          NEWIND = NEWIND + EL[Ii]
          
        }
        
        
        
      } # fin for clases
      
      
      
      vectorN[[k]][[2]] = c(vectorN[[k]][[2]][1], stagesN[2:length(stagesN)])
      
      vectorN[[k]][[1]] = vectorN[[k]][[2]]
      
      
      
      if(k <= (numInmaduros-1)){
        
        vectorN[[k]][[3]][1]=NEWIND
        
        vectorN[[k+1]][[2]] = vectorN[[k]][[3]]
        
      }
      
      
      
      if(k == numInmaduros){
        
        vectorN[[k+1]][[3]][1] = NEWIND * Rs[Day]
        
        vectorN[[k+2]][[3]][1] = NEWIND * (1-Rs[Day])
        
        
        
        vectorN[[k+1]][[2]] = vectorN[[k]][[3]]#
        
      }
      
      
      
      if(k == (numInmaduros+1)){
        
        femalesN = E
        
      }
      
      
      
      matrizOut[Day+1,k] = round(sum(vectorN[[k]][[2]]))
      
      #matrizOut[1:3,]
      
      
      
    }#fin for estadios
    
    
    
    
    
    #Oviposition
    
    Ageclass = ageClassMatriz[Day,numInmaduros+1]
    
    matrizAges = estadiosAges[[numInmaduros+1]]
    
    Ovi = c(rep(0,Ageclass))
    
    
    
    E <- c(rep(0, Ageclass))
    
    EL <- c(rep(0, Ageclass))
    
    Eage <- c(rep(0, Ageclass))
    
    pacc <- c(rep(0, Ageclass))
    
    p1 <- c(rep(0, Ageclass))
    
    p2 <- c(rep(0, Ageclass))
    
    
    
    for(Ii in 1:(Ageclass-2)){
      
      E[Ii] = vectorN[[k-1]][[2]][Ii]
      
      Eage[Ii] = oviFreq[Day, Ii]
      
      E1 = matrizA[Day,k*2-1]
      
      E1 = E1 * Eage[Ii] * E[Ii]
      
      Ovi[Ii+1] = E1
      
      NEWIND = NEWIND + E1
      
    }
    
    
    
    Ovi = na.omit(Ovi)
    
    matrizOut[Day+1,(numInmaduros + numMaduros +1)] = round(sum(Ovi),2)
    
    
    
    #vectorN[[1]][[2]][1] = NEWIND
    
    matrizOut[Day+1,1] = round(sum(vectorN[[1]][[2]]))
    
    vectorN[[1]][[1]] = vectorN[[1]][[2]]
    
    
    
    #                      if(round(matrizOut[Day,(numInmaduros + numMaduros)]) > 0 && round(matrizOut[Day+1,(numInmaduros + numMaduros)]) == 0){
    
    #                                  break;
    
    #                      }
    
    if(matrizOut[Day,1]==0 && matrizOut[Day,2]==0 && matrizOut[Day,3]==0 && matrizOut[Day,4]==0 && matrizOut[Day,5]==0){
      
      #matrizOut[Day:nrow(matrizOut),]="Dead"
      break;
      
    }
    
    
    
    NEWIND=0
    
    Day=Day+1
    
    Step = Step + 1
    
  }#fin bucle steps
  
  
  
  
  
  namesMatriz = c(estadios, "New Egg");
  
  colnames(matrizOut) = namesMatriz
  
  matrizOutDead = matrizOut
  matrizOutDead[is.na(matrizOutDead)]="Dead"
  
  matrizOut[is.na(matrizOut)]=0
  
  return(list(matrizOut=matrizOut, vectorN=vectorN, matrizOutDead=matrizOutDead))
}
