#' simulacion
#'
#' no details
#' 
#' @param cuadro matrix, predictions of phenology variables by day
#' @param estadios character vector, all states
#' @param n integer value, number of insects to simulate
#' @param m integer value, number of days
#' @param hfeno list, objects about phenology
#' @param Rs numeric vector, sexual ratios by day
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
simulacion<-function(cuadro,estadios,n,m,hfeno,Rs){
  mat<-matrix("dead",ncol=n,nrow=(m))
  for(z in 1:n){
    k <- 1
    Insect_age <- 0
    survival_p <- 1
    dead <- 0
    Day <- 1
    slope<- 0
    djx <- 0
    mjx <- 0
    kj <- NA
    rand_d <- round(runif(1),3) #' random number between 0 and 1 with 3 digits
    rand_m <- round(runif(1),3)
    rand_ovi <- round(runif(1),3)
    rand_sex = 0
    
    if(rand_ovi == 0){
      rand_ovi <- 0.0001
    }
    
    sizeInmaduros = length(estadios[-(length(estadios)-1):-(length(estadios))])
    numFemale = sizeInmaduros + 1
    numMale = sizeInmaduros + 2
    
    #loop for immature life stages
    while(dead != 1 || k <= sizeInmaduros){
      if(k > sizeInmaduros){
        break
      }
      if(dead ==1){
        break
      }
      
      #read the parameter for variation and for each specific stage
      slope = hfeno$slope_dv[[k]]
      
      #read value from Matrix "A"
      djx <- cuadro[Day,k*2-1]
      mjx <- cuadro[Day,k*2]
      
      #calculate physiological age
      if(Insect_age == 0){
        Insect_age = Insect_age + (djx / 2)
      }else{
        Insect_age <- Insect_age + djx
      }
      
      dev_p <- distrimodel(Insect_age,k,hfeno)
      
      ################################################
      # Se evalua si el insecto muere en cierto estado
      ################################################
      
      #evaluate survival
      survival_p = survival_p * ((1 - mjx) ^ djx)  ## esta probabilidad es muy alta
      survival_rand = (1 - rand_m)
      
      if(survival_p < survival_rand || Day==365){
        dead = 1
        mat[Day,z] <- "dead"
        survival_p = 1
        k = 1
      }
      
      ##################################################
      # Se evalua si el insecto pasa al siguiente estado
      ##################################################
      
      #evaluate development to next stage
      if(dev_p > rand_d){   ### si dev_p es demasiado peque?o nunca va pasar al siguiente estado fenol?gico
        #read djx from Matrix "A" of the next stage
        djx <- cuadro[Day,(k+1)*2-1]
        Insect_age = djx / 2
        
        rand_d <- round(runif(1),3) #' random number for the next stage
        rand_m <- round(runif(1),3)
        
        survival_p = 1
        k = k + 1
      }
      
      #inmature dev stages
      if(k <= sizeInmaduros){
        kj = estadios[k]
      }
      
      if(k > sizeInmaduros){
        Insect_age = 0
      }
      if(dead == 1){
        kj = "dead"
      }
      mat[Day,z] <- kj
      Day = Day + 1
      dead
      
    }
    # end loop for inmatures
    
    #print(paste(Day,Rs[Day],sep="_")) ## Taza calculada en base ala temperatura que tiene seg?n el dia cuando el insecto termina de ser inmaduro
    if(k == (sizeInmaduros+1) && Insect_age == 0 && Day!=nrow(cuadro) + 1){
      rand_sex <- round(runif(1),3)
      
      if(rand_sex < Rs[Day]){ ## Aqu? se aplica
        k = numFemale
      }else{
        k = numMale
      }
    }
    
    #loop for adult males
    while(dead != 1 && k != numFemale && Day!=nrow(cuadro) + 1){
      if(dead ==1){
        break
      }
      
      #read the parameter for variation and for each specific stage
      slope = hfeno$slope_snm
      
      #read value from Matrix "A"
      djx <- cuadro[Day,k*2-2]
      
      #calculate physiological age
      if(Insect_age == 0){
        Insect_age = Insect_age + (djx / 2)
      }else{
        Insect_age <- Insect_age + djx
      }
      
      dev_p <- distriModelAdults(Insect_age,k, sizeInmaduros,hfeno)
      #dev_p = log(Insect_age)
      #dev_p <- distrimodel(dev_p,k)
      
      #evaluate survival of males
      if(dev_p > rand_d){
        dead <- 1
        Insect_age = 0
        
        rand_d <- round(runif(1),3) #' random number for the next stage
        rand_m <- round(runif(1),3)
        
        survival_p = 1
      }
      
      kj = estadios[k]
      
      if(dead == 1){
        kj = "dead"
      }
      
      mat[Day,z] <- kj
      Day = Day + 1
    }# end loop for males
    #print(Day)
    #loop for adult females
    while(dead != 1 && k != numMale && Day!= nrow(cuadro) + 1){
      if(dead ==1){
        break
      }
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
      
      slope = hfeno$slope_snh
      djx <- cuadro[Day,k*2-1]
      
      #calculate physiological age
      if(Insect_age == 0){
        djx <- djx + cuadro[Day-1,k*2-1]
      }
      
      Insect_age = Insect_age + djx
      dev_p <- distriModelAdults(Insect_age,k, sizeInmaduros,hfeno)
      #dev_p = log(Insect_age)
      #dev_p <- distrimodel(dev_p,k)
      
      #evaluate survival of females
      if(dev_p > rand_d){
        dead = 1
        Insect_age = 0
        rand_d <- round(runif(1),3) #' random number for the next stage
        rand_m <- round(runif(1),3)
        
        survival_p = 1
      }
      
      #calculate oviposition
      #read parameters
      if(dead < 1){
        alfa = hfeno$povih_h[[1]]
        beta = hfeno$povih_h[[2]]
        cv = 0.3 # Averiguar acerca de este valor
        
        #read ovitot from Matrix "A"
        Ovitot = cuadro[Day-1,k*2+1]
        
        #calculate age-dependent oviposition frequency
        x  <-  Insect_age - djx
        OviFrec1  <- eval(funcionc)
        x  <-  Insect_age
        OviFrec2  <- eval(funcionc)
        OviFrec = OviFrec2 - OviFrec1
        
        #calculate oviposition
        Ovitot = qnorm(rand_ovi, Ovitot, Ovitot * cv);Ovitot[Ovitot<0]=0
        Oviposition = Ovitot * OviFrec
        Oviposition = round(Oviposition, 0)
        kj<- Oviposition
      }
      
      if(dead == 1){
        kj = "dead"
      }
      #print("error5");print(paste(Day,z))
      
      mat[Day,z] <- kj
      Day = Day + 1
    }#end loop females and oviposition
    
    Day = 1
  }# end loop individuos
  return(list(mat=data.frame(mat)))
}
