#' life.tableD
#'
#' no details
#' 
#' @param cuadro matrix, predictions of phenology variables by day
#' @param estad character vector, all states
#' @param n integer value, number of insects to simulate
#' @param m integer value, number of days
#' @param hfeno list, objects about phenology
#' @param Rs numeric vector, sexual ratios by day
#' @param Amp integer value
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
life.tableD<-function(cuadro,estad,n,m,hfeno,Rs,Amp=1)
{
  tab1 <- matrix(0,nrow = m,ncol = (length(estad)-1))
  tab2 <- matrix(0,nrow = m,ncol = 3)
  
  k=1:(length(estad)-2)
  k2=length(k)+1
  
  colnames(tab1)=c(estad[k],"nMad")
  colnames(tab2)=c(estad[-k],"Oviposicion")
  
  tab1[1,1]=n
  
  Insect_age <- (cuadro[1,k*2-1] / 2)
  Insect_age_m<-rep(0,2)
  survival_p <- 1
  Rat<-c(1);Rat[1]=0;Ov=0
  
  for(Day in 2:m)
  {
    dev_p<-dev_p_m<-c(1)
    djx <- cuadro[Day-1,k*2-1]
    mjx <- cuadro[Day-1,k*2]
    Insect_age <- Insect_age + djx
    for(i in k){dev_p[i]=distrimodel(Insect_age[i],k[i])} 
    survival_p = survival_p * ((1 - mjx/Amp) ^ (djx/Amp)) 
    nviv=round(survival_p*tab1[Day-1,-ncol(tab1)])
    ndes=round(dev_p*nviv)
    ndif=nviv-ndes
    tab1[Day,]=c(ndif,0)+c(0,ndes)
    for(i in k){if(tab1[Day-1,k[i]]==0 & tab1[Day,k[i]]>0){survival_p[i] = 1;Insect_age[i] = djx[i] / 2}}
    
    djx <- c(cuadro[Day,k2*2-1],cuadro[Day,k2*2])
    if(tab1[Day-1,k2]==0 & tab1[Day,k2]>0){Insect_age_m = c(cuadro[Day,k2*2-1]/2,cuadro[Day,k2*2]/2)}else{Insect_age_m = Insect_age_m + (djx / 2)}
    Insect_age_m = Insect_age_m + (djx / 2)
    
    for(j in 1:2){dev_p_m[j]=distriModelAdults(Insect_age_m[j],k2+j-1,length(k))}

    ndes2=c(round(Rs[Day]*tab1[Day,k2]),round((1-Rs[Day])*tab1[Day,k2]))
    nvivt=c(tab2[Day-1,1],tab2[Day-1,2])
    nviv2=nvivt-round(c(dev_p_m[1]*nvivt[1],dev_p_m[2]*nvivt[2]))
    
    ntot=ndes2+nviv2
    tab2[Day,]=c(ntot,cuadro[Day,k2*2+1]*ntot[1])
    
    if(tab2[Day,1]>0){Ov=1}

    parametrosc <- hfeno$povih_h
    for (i in names(parametrosc)){temp <- parametrosc[i];storage.mode(temp) <- "double";assign(i, temp)}
    formulac <- hfeno$fovih_h
    funcionO <- as.expression(formulac[[3]])
    
    x=Insect_age_m[1];Rat[Day]=eval(funcionO)*Ov
    
  }
  Ov2=rep(1,m);Ov2[tab2[,1]==0]=0
  Rat=descum(Rat)
  lifeT=data.frame(tab1[,-k2],tab2[,1:2],Oviposicion=round(tab2[,3]*Rat*Ov2))
  ncero=apply(lifeT,1,sum);ncero=ncero==0
  lifeT=lifeT[!ncero,]
  return(list(life.table=lifeT))
  
}
