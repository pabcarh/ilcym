#' ovip.time
#'
#' no details
#' 
#' @param dat data frame
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
ovip.time<-function(dat)
{
  dat <- data.frame(matrix(as.numeric(as.matrix(dat)),nrow(dat),ncol(dat)))
  tmps=unique(c(dat[,1]))
  
  mat2=data.frame(Temp=NA,Dia=NA,Frec=NA,Frec.A=NA,N=NA,Median.W.Days=NA);n<-c(1)
  for(i in 1:length(tmps))
  {
    dat1=dat[dat[,1]==tmps[i],-1];tovc=apply(dat1,1,tovip);mat1=table(tovc)
    if(length(mat1) > 1){
      if(names(mat1[1])=="0"){n[i]=nrow(dat1)-mat1[1];mat1=mat1[-1]}else{n[i]=nrow(dat1)} ## no considero los ceros
      mmd=median(rep(as.numeric(names(mat1)),mat1))
      tempo=data.frame(Temp=tmps[i],Dia=names(mat1),Frec=c(mat1),Frec.A=cumsum(mat1),N=n[i],Median.W.Days=mmd)
      mat2=rbind(mat2,tempo)
    }
    
    #                dat1=dat[dat[,1]==tmps[i],-1];tovc=apply(dat1,1,tovip);mat1=table(tovc)
    #                if(names(mat1[1])=="0"){n[i]=nrow(dat1)-mat1[1];mat1=mat1[-1]}else{n[i]=nrow(dat1)} ## no considero los ceros
    #                mmd=median(rep(as.numeric(names(mat1)),mat1))
    #                tempo=data.frame(Temp=tmps[i],Dia=names(mat1),Frec=mat1,Frec.A=cumsum(mat1),N=n[i],Median.W.Days=mmd)
    #                mat2=rbind(mat2,tempo)
  }
  mat2=mat2[-1,]
  mat2[,2]=as.numeric(mat2[,2]);mat2[,1]=as.factor(mat2[,1])
  
  colnames(mat2)=c("T","Day's","Freq","Cum.Freq","Sample","Median.W.Days")
  cat("\n")
  cat("\nTABLE OF FREQUENCY\n")
  print(mat2)
  
  return(list(mat2=mat2,tmps=tmps))
}
