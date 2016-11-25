#' addmale
#'
#' no details
#'
#' @param datapre data frame
#' @param tp numeric vector of the constant temperature
#' @param estadios character vector, all states
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
addmale<-function(datapre,tp,estadios)
{
  Male<-data.frame(cbind(tp[1],1,sum(subset(datapre[[length(estadios)-2]][,4],datapre[[length(estadios)-2]][,1]==tp[1])),sum(subset(datapre[[length(estadios)-2]][,4],datapre[[length(estadios)-2]][,1]==tp[1]))))
  for (i in 2:length(tp))
  {
    Male<-rbind(Male,data.frame(cbind(tp[i],1,sum(subset(datapre[[length(estadios)-2]][,4],datapre[[length(estadios)-2]][,1]==tp[i])),sum(subset(datapre[[length(estadios)-2]][,4],datapre[[length(estadios)-2]][,1]==tp[i])))))
  }
  return(Male)
}