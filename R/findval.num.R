#' findval.num
#'
#' no details
#'
#' @param dat data frame
#'
#' @return None
#'
#' @examples
#' #building
#'
#' @export
findval.num<-function(dat)
{n1=ncol(dat)
for(i in 1:n1)
{
  p1=unique(dat[,i])
  p1=as.character(p1)
  p1=as.numeric(p1)
  exist=FALSE
  for(j in 1:length(p1))
  {
    pp=p1[j]+1
    if(!is.na(pp)){exist=TRUE;break}
  }
  if(exist==TRUE){break}
}
return(list(exist=exist))
}