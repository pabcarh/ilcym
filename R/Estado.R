#' Estado
#'
#' no details
#' 
#' @param data matriz, observed life table (fluctuating)
#' @param estad character vector, all states
#' @param est a character
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
Estado <- function (data, estad,est){
  maduros   <-  estad[(length(estad)-1):(length(estad))]
  if(est==maduros[1]){
    datat <- as.matrix(data)
    num <- matrix(as.numeric(datat), nrow(datat), ncol(datat))
    a <- as.list(as.data.frame(t(num)))
    b <- as.list(rep(0, nrow(data)))
    for (i in 1:nrow(data)) b[[i]] <- sum(as.matrix(table(a[[i]])))
    female <- as.numeric(as.character(as.matrix(b)))
    return(female)
    
  }
  if(est!=maduros[1]){
    datat <- as.matrix(data)
    estadoa <- as.list(as.data.frame(t(datat)))
    estadob <- as.list(rep(0, nrow(data)))
    for (i in 1:nrow(data)) {
      estadoa[[i]] <- as.matrix(estadoa[[i]])
      estadob[[i]] <- length(subset(estadoa[[i]], estadoa[[i]] == est))
    }
    estados <- as.numeric(as.character(as.matrix(estadob)))
    return(estados)
  }
}