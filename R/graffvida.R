#' graffvida
#'
#' no details
#' 
#' @param Tablelife matriz, counts for the observed life table
#' @param cuadro matrix, predictions of phenology variables by day
#' @param estadios character vector, all states
#' @param N integer value, number of insects to simulate
#' @param M integer value, number of days
#' @param hfeno list, objects about phenology
#' @param labx a character
#' @param laby a character
#' @param corrx numeric vector of two values, the axis range for X
#' @param corry numeric vector of two values, the axis range for y
#' @param lgx double, x position of the legend
#' @param lgy double, y position of the legend
#' @param Rs numeric vector, sexual ratios by day
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
graffvida<-function(Tablelife,cuadro,estadios,N,M,hfeno,labx,laby,corrx,corry,lgx,lgy, Rs){
  table_life<-estimado<-list(1)
  NUM <- 8
  TABLEmean<-matrix(0, nrow=nrow(cuadro), ncol=length(estadios)+1)
  TABLEmin<-matrix(0, nrow=nrow(cuadro), ncol=length(estadios)+1)
  TABLEmax<-matrix(0, nrow=nrow(cuadro), ncol=length(estadios)+1)
  for(k in 1:NUM)
  {
    estimado [[k]]<- simulacion(cuadro,estadios,N,M,hfeno,Rs)$mat
    table_life[[k]] <- life.table2(estimado[[k]], estadios)$life.table
  }
  
  for(i in 1:ncol(table_life[[1]]))
  {
    for(j in 1:(nrow(table_life[[1]])-1))
    {
      TABLEmin[j,i] <- min(c(table_life[[1]][j,i],table_life[[2]][j,i],table_life[[3]][j,i],table_life[[4]][j,i]))
      TABLEmax[j,i] <- max(c(table_life[[1]][j,i],table_life[[2]][j,i],table_life[[3]][j,i],table_life[[4]][j,i]))
      TABLEmean[j,i]<- mean(c(table_life[[1]][j,i],table_life[[2]][j,i],table_life[[3]][j,i],table_life[[4]][j,i]))
    }
  }
  plot(1:nrow(TABLEmean),TABLEmean[,1],ylim=corry,xlim=corrx,type="l", col=1,lwd=1,frame=F,ylab=laby,xlab=labx)
  
  for(i in 1:length(estadios))
  {
    temp1=TABLEmean[,i];temp1[temp1<0.5]=NA;lines(1:nrow(TABLEmean),temp1,lwd=1,col=i)
    temp2=TABLEmin[,i];temp2[temp2<0.5]=NA;lines(1:nrow(TABLEmin),temp2,lwd=1,col=i,lty=4)
    temp3=TABLEmax[,i];temp3[temp3<0.5]=NA;lines(1:nrow(TABLEmax),temp3,lwd=1,col=i,lty=4)
  }
  newTablelife <- matrix(NA,nrow(Tablelife), ncol(Tablelife))
  matrizPosic <- matrix(NA,nrow(Tablelife), length(estadios))
  posic1 <- 1
  for(j in 1:length(estadios))
  {
    posic1 <- 1
    temp1 = Tablelife[1,j]
    for(t in 2:length(Tablelife[,j])){
      if(Tablelife[(t-1),j] == Tablelife[t,j]){
        temp1 <- c(temp1,NA)
        posic1 <- c(posic1, t)
      }else{
        temp1 <- c(temp1,Tablelife[t,j])
      }
      
    }
    posic1 <- posic1[2:length(posic1)]
    matrizPosic[,j] = c(posic1,rep(0, (nrow(Tablelife)-length(posic1))))
    temp <- temp1;
    temp[is.na(temp)] = 0
    newTablelife[,j] <- temp
    temp[temp<0.5]=NA#
    points(1:length(temp),temp,col=j, pch=19)
  }
  legend(lgx,lgy,c(estadios),col=1:length(estadios),lty=rep(1,length(estadios)))
  
  newTablelife[,length(estadios)+1] = Tablelife[,length(estadios)+1]

  Indic.sim <- function(y,yl){sqe <- sum((y-yl)^2);DE <- sqrt(sqe);return(DE)}
  itab<-matrix(NA,length(estadios),1);colnames(itab)=c("Euclidian.Dist");rownames(itab)=estadios
  cat("\n", "\n","\n","Fitting indicator for each state","\n")
  for(i in 1:length(estadios))
  {
    y=c(Tablelife[-matrizPosic[,i],i],rep(0,100))
    yl=c(TABLEmean[-matrizPosic[,i],i],rep(0,100))
    id1=(1:length(y))[y!=0];id2=(1:length(yl))[yl!=0];id=union(id1,id2)
    y=y[id]+ 1e-12;yl=yl[id]+ 1e-12
    y=c(y,rep(1e-12,4));yl=c(yl,rep(1e-12,4)) 
    itab[i,]=Indic.sim(y,yl)
  }
  print(itab)
  
}
