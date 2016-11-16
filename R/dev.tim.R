#' dev.tim
#'
#' no details
#'
#' @param matri data frame, counts of insect by days and others variables
#' @param tp numeric vector of the constant temperature
#'
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
dev.tim<-function(matri,tp)
{
  tmed<-c(1);n1<-c(1);acum1<-c(1)
  for(i in 1:length(tp))
  {
    mat1=matri[matri[,1]==tp[i],];n1[i]=nrow(mat1)
    #if(n1[i]!=1){acum1<-c(acum1,cumsum(mat1[,5]))}else{acum1<-c(acum1,mat1[1,4])}
    acum1<-c(acum1,cumsum(mat1[,5]))
    nn=rep(as.numeric(mat1[,2]),mat1[,5])
    if(length(nn) != 0){tmed[i]=stats::median(nn)}else{tmed[i]=0}
  }
  acum1=acum1[-1]
  T=rep(tp,n1)
  return(list(acum1=acum1,T=T,tmed=tmed,n1=n1))
}