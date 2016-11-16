#' verif.val.num
#'
#' no details
#'
#' @param datos data frame
#' @param cate character
#' @param estadios vector of length variable
#'
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
verif.val.num<-function(datos,cate,estadios)
{
  if(cate==estadios[length(estadios)-1])
  {
    verif=rep(FALSE,length(datos))
    for(j in 1:length(datos))
    {
      verif[j]=findval.num(datos[[j]])$exist
      #if(exist==TRUE){break} ## corregir
    }
    return(verif)
  }
}