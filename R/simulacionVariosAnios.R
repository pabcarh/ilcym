#' simulacionVariosAnios
#'
#' no details
#' 
#' @param matrizOut matrix, life parameters estimated by temperature
#' @param vectorN list
#' @param estadios character vector, all states
#' @param estadiosAges list
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
simulacionVariosAnios <- function(matrizOut, vectorN, estadios, estadiosAges){
  
  inmad <-  estadios[-(length(estadios)-1):-(length(estadios))]
  numInmaduros = length(inmad)
  maduros   <-  estadios[(length(estadios)-1):(length(estadios))]
  numMaduros = length(maduros)
  sizeMatriz = nrow(matrizOut)
  numTotal <- sum(matrizOut[sizeMatriz,1:(numInmaduros+numMaduros)])
  numIni <- 0
  numInd <- list(1)
  
  for(k in 1:(numInmaduros+numMaduros)){
    if(k <= numInmaduros){
      slope= 2
    }else{
      slope= 1.8
    }
    
    age=0
    class=1
    matrizAges = estadiosAges[[k]]
    while(age < slope){
      age =  matrizAges[sizeMatriz,class]
      Number = vectorN[[k]][[2]][class] / numTotal * 100 
      estadiosAges[[k]][1,class] = age
      vectorN[[k]][[1]][class] = Number
      class = class+1
    }
    
    numInd[[k]] <- sum(vectorN[[k]][[2]] / numTotal * 100)
  }
  return(list(vectorN=vectorN,numInd=numInd))	
}
