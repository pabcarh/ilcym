#' verif.val
#'
#' development time sub function
#'
#' no details
#'
#' @param datos data frame.
#' @param cate a character.
#' 
#' @return None.
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#' 
#' @export
verif.val<-function(datos,cate)
{
  verif=rep(FALSE,length(datos))
  for(j in 1:length(datos))
  {
    verif[j]=findval(datos[[j]],cate)$exist
  }
  return(verif)
}

