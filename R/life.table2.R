#' life.table2
#'
#' no details
#' 
#' @param data matriz, observed life table (fluctuating)
#' @param estad character vector, all states
#' @param ovi logical, default is TRUE, if there is a table of oviposition (fluctuating)
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
life.table2<-function (data, estad,ovi=TRUE){
  if(ovi)
  {
    maduros   <-  estad[(length(estad)-1):(length(estad))]
    inmaduros <-  estad[-(length(estad)-1):-(length(estad))]
    estadios <- estad
    datat <- as.matrix(data)
    num <- matrix(as.numeric(datat), nrow(datat), ncol(datat))
    numc <- num
    for (i in 1:nrow(numc)) for (j in 1:ncol(numc)) if (is.na(numc[i, j])) numc[i, j] <- 0
    newegg <- matrix(apply(numc, 1, sum), nrow(numc))
    day <- 0:(nrow(num) - 1)
    oviposition <- apply(numc, 2, sum)
    estados1 <- matrix(rep(0, length(day)), nrow = length(day),
                       ncol = length(estadios))
    for (i in 1:day) for (j in 1:length(estadios)) estados1[, j] <- Estado(data, estad, estadios[j])
    life.table <- data.frame(estados1, newegg)
    nombres <- rep(0, length(estadios))
    for (i in 1:length(estadios)) nombres[i] <- estadios[i]
    colnames(life.table) <- c(nombres, "new.Egg")
    #trues <- rep(TRUE, nrow(life.table))
    #for (i in 1:nrow(life.table)) trues[i] <- sum(life.table[i, ]) != 0
    #life.table <- life.table[trues, ]
    #print(life.table)
    return(list(life.table=life.table))
  }else
  {
    datos=data;estad=estadios
    n=ncol(datos);m=nrow(datos)
    #n1=ncol(ovipo);m1=nrow(ovipo)
    tab1 <- matrix(0,nrow = m,ncol = (length(estad)+1))
    colnames(tab1)=c(estad,"new.Egg")
    #dat=matrix(as.numeric(as.matrix(ovipo)),m1,n1)
    for(i in 1:(length(estad)+1))
    {
      #if(i==(length(estad)+1)){tab1[,i]=apply(dat,1,sum,na.rm=T);next} ## oviposicion
      if(i==(length(estad)+1)){tab1[,i]=rep(0,m);next} ## oviposicion
      tab1[,i]=apply(datos,1,contar1,estad[i])
    }
    return(list(life.table=data.frame(tab1)))
  } 
}
