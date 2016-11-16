#' temp
#'
#' development time sub function
#'
#' no details
#'
#' @param dat data frame.
#' @param cate vector of characters.
#' 
#' @return None.
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#' 
#' @export
findval<-function(dat,cate)
{n1=ncol(dat)
result=rep(0,n1)
for(i in 1:n1)
{
  p1=unique(dat[,i])==cate
  if(length(p1[p1==TRUE])>=1){result[i]=1}else{result[i]=0}
}
ntot=sum(result)
if(ntot>=1){exist=TRUE}else{exist=FALSE}
return(list(ntot=ntot,exist=exist))
}
