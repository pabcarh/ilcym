#' AgeClases
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
AgeClases<- function(cuadro, estadios, hfeno){
  matrizA=cuadro
  sizeMatrizA = nrow(matrizA)
  inmad <-  estadios[-(length(estadios)-1):-(length(estadios))]
  numInmaduros = length(inmad)
  maduros   <-  estadios[(length(estadios)-1):(length(estadios))]
  numMaduros = length(maduros)
  estadiosAges<-list(1)
  ageClassMatriz = matrix(NA, nrow=sizeMatrizA, ncol=(numInmaduros+numMaduros))
  
  Day = 1
  Ageclass = 1
  age = 0
  
  #primera linea
  for(k in 1:(numInmaduros + numMaduros)){
    Ageclass = 1
    age = 0
    Ageclass0 = c(rep(NA, 253))# porq 253??
    matrizAges = matrix(NA, nrow=sizeMatrizA, ncol=(length(Ageclass0)+1))
    #matrizAges = matrix(NA, nrow=366, ncol=(length(Ageclass0)+1))
    matrizAges[Day,1] = 0
    
    ## if(k==1){
    ##     Ageclass = Ageclass + 1
    ## }
    
    Ageclass = Ageclass + 1
    #for(Ii in 2:length(Ageclass0)){
    for(Ii in 1:length(Ageclass0)){
      if(k == (numInmaduros + numMaduros)){
        Ageclass0[Ii] = matrizA[(sizeMatrizA+1)-Ii,k*2-2]
      }else{
        Ageclass0[Ii] = matrizA[(sizeMatrizA+1)-Ii,k*2-1]
      }
      
      age = age + Ageclass0[Ii]
      matrizAges[Day,Ii+1] = age
    }
    if(k <= numInmaduros){
      ageClassMatriz[Day,k] = length(matrizAges[1,][matrizAges[1,]<2])+1
    }else{
      ageClassMatriz[Day,k] = length(matrizAges[1,][matrizAges[1,]<1.8])+1
    }
    estadiosAges[[k]] = matrizAges
  }
  
  # a partir de la 2da linea
  for(k in 1:(numInmaduros + numMaduros)){
    Day=2
    Ageclass = 1
    age = 0
    matrizAges = estadiosAges[[k]]
    while(Day <= sizeMatrizA){#porq 367 en excel?
      matrizAges[Day,1] = 0
      if(k == (numInmaduros + numMaduros)){
        Rate = matrizA[Day,k*2-2]
      }else{
        Rate = matrizA[Day,k*2-1]
      }
      
      if(k <= numInmaduros){
        while(age < 2){ # cambiar a valor del slop
          ageacc = matrizAges[Day-1, Ageclass]
          age = ageacc + Rate
          matrizAges[Day, Ageclass+1] = age
          Ageclass = Ageclass + 1
          
          if(Ageclass == 254){ ## por que 254?????
            age = 3
          }
        }
        ageClassMatriz[Day,k] = Ageclass
      }else{
        while(age < 1.8){ # cambiar a valor del slop
          ageacc = matrizAges[Day-1, Ageclass]
          age = ageacc + Rate
          matrizAges[Day, Ageclass+1] = age
          Ageclass = Ageclass + 1
          
          
          if(Ageclass == 254){     ## antes : if(Ageclass == 255){
            age = 3
          }
        }
        ageClassMatriz[Day,k] = Ageclass
      }
      age = 0
      Day = Day + 1
      Ageclass = 1
    }
    estadiosAges[[k]] = matrizAges
  }
  
  ## Ovifreq ##
  matrizAges = estadiosAges[[numInmaduros+1]]
  oviFreq = matrix(NA, nrow=sizeMatrizA, ncol=length(Ageclass0))
  Day = 1
  Ageclass = 1
  age1 = 0
  age2 = 0
  
  parametrosc <- hfeno$povih_h
  parametrosc<-as.list(parametrosc)
  for (i in names(parametrosc)){
    temp <- parametrosc[[i]]
    storage.mode(temp) <- "double"
    assign(i, temp)
  }
  
  formulac <- hfeno$fovih_h
  forexc <- formulac[[length(formulac)]]
  funcionc <- as.expression(forexc)
  
  while(Day < sizeMatrizA){#porq 367 en excel?
    matrizAges[Day,1] = 0
    age1 = matrizAges[Day, Ageclass]
    x = age1
    OviFrec1 = eval(funcionc)
    
    while(age2 < 1.8 & Ageclass<254){
      age2 = matrizAges[Day, Ageclass]
      x= age2
      OviFrec2 = eval(funcionc)
      OviFrec = OviFrec2 - OviFrec1
      OviFrec1 = OviFrec2
      if(Ageclass == 1){
        OviFrec = 0
      }
      oviFreq[Day,Ageclass] = OviFrec
      Ageclass = Ageclass + 1
    }
    age2 = 0
    Day = Day + 1
    Ageclass = 1
  }
  
  
  return(list(estadiosAges=estadiosAges, oviFreq=oviFreq, ageClassMatriz=ageClassMatriz))
}
