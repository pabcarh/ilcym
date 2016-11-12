#' temp
#'
#' A special curve using specific parameters and also its modelling
#'
#' If no estimated values of the parameters of the model must enable "modeling = TRUE", then include the initial values
#'
#' @param datos data frame.
#' @param estadios vector of characters.
#' @param est a character.
#' @return character a string; default is 'finished!'.
#' @author Pablo Carhuapoma Ramos
#' @family example
#' @example inst/examples/example_temp.R
#' @export
temp <- function(datos,estadios,est) {
  tip<-1:length(estadios)
  for (is in 1:(length(estadios))) if (estadios[is] == est) dat<-datos[[is]]
  not <- as.numeric(levels(factor(dat[, 1])))
  return(list(temperature=not))
}
