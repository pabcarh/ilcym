#' mort.select
#'
#' no details
#'
#' @param datos data frame
#' @param estadios vector of length variable
#' @param num integer value
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
mort.select<-function(datos, estadios, num)
{
  datas<- mort<-list(1)
  mortalidad<-c(1)
  tip<-1:length(estadios)
  #for (is in 1:(length(estadios))) if (estadios[is] == est) num<-tip[is]
  for (i in 1:length(datos)) {
    datas[i]<-list(datos[[i]][,num])
    for (j in 1:length(datas[[i]])) if(is.na(datas[[i]][j])) datas[[i]][j]<-"NA"
    mort[i] <- list(subset(datas[[i]], datas[[i]] == "NA"))
  }
  return(mort)
}
