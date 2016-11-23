#' matriz
#'
#' no details
#' 
#' @param data data frame
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
matriz<-function(data)
{
  mes<-c(31,28,31,30,31,30,31,31,30,31,30,31)
  mat<-cbind(round(seq(data[,1][1],data[,1][2],length=mes[1]),1),round(seq(data[,2][1],data[,2][2],length=mes[1]),1))
  for(i in 2:nrow(data))
  {
    if(i != nrow(data)) mat<-rbind(mat,cbind(round(seq(data[,1][i],data[,1][i+1],length=mes[i]),1),round(seq(data[,2][i],data[,2][i+1],length=mes[i]),1)))
    else mat<-rbind(mat,cbind(round(seq(data[,1][i],data[,1][i-11],length=mes[i]),1),round(seq(data[,2][i],data[,2][i-11],length=mes[i]),1)))
  }
  return(list(temperaturas=mat))
}