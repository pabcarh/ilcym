#' conteo
#'
#' no details
#'
#' @param data data frame
#' @param est character, current state
#'
#' @return None
#'
#' @examples
#' #building
#'
#' @export
conteo <- function (data,est) 
{
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
