#' PARAlife
#'
#' no details
#' 
#' @param N integer value, number of simulations
#' @param estad character vector, all states
#' @param ltb matrix, life table counts
#' @param Rs2 numeric vector, sexual ratios by day with position
#' @param ovipo matrix,table of oviposition (fluctuating)
#' @param nfem2 double, number of females obtained from the external table (fluctuating)
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
PARAlife<-function (N, estad, ltb,Rs2=NULL,ovipo=NULL,nfem2=NULL){
  if(is.null(ovipo))
  {
    estadios <- estad
    lifetable <- ltb
    s.x <- m.x <- l.x <- rep(0, nrow(lifetable))
    for (i in 1:nrow(lifetable)) l.x[i] <- sum(lifetable[i, 1:(length(estadios))])/N ## es un porcentaje de los conteos
    for (i in 1:nrow(lifetable)) if(l.x[i]==0) m.x[i] <- 0 else m.x[i] <- lifetable[,(length(estad)+1)][i]/N/l.x[i] ## este objeto d
    s.x1 <- c(l.x[1], l.x[2])
    s.x2 <- rep(0, nrow(lifetable) - 2)
    for (i in 1:(nrow(lifetable) - 2)) s.x2[i] <- l.x[i + 2]/l.x[i + 1]
    s.x <- c(s.x1, s.x2)
    Ro <- sum(l.x * m.x)
    Day <- 0:(nrow(lifetable)-1)
    Tg <- sum(Day * l.x * m.x)/Ro
    GRR <- sum(m.x)
    r <- log(Ro)/Tg
    euler <- sum(l.x * m.x * exp(-r * (Day + 1)))
    r1 <- r - 1e-06
    for (k in 1:100) {
      euler1 <- sum(l.x * m.x * exp(-r1 * (Day + 1)))
      r2 <- ((1 - euler) * (r1 - r)/(euler1 - euler)) + r
      r1 <- r2
    }
    rf <- r1
    Dt <- log(2)/rf
    lamda <- exp(rf)
    T <- log(Ro)/rf
    parameters <- t(data.frame(r = rf, Ro = Ro, GRR = GRR, T = T, lambda = lamda, Dt = Dt))
    colnames(parameters) <- "Parameters"
    #print(parameters)
    return(list(parametro=parameters))
  }else
  {
    n1=length(estadios)
    n2=min((1:nrow(ltb))[ltb[,n1-1]!=0]) ## AQUI PONGO -1 (vector de hembras) CONSIDERANDO QUE HAY 2 ADULTOS, ESTO PUEDE CAMBIAR PARA LOS NUEVOS DATOS
    W1=n2+nrow(ovipo)-1
    x=n2:W1+1.5
    n1=ncol(ovipo);m1=nrow(ovipo);c1=ncol(ovipo)
    dat=matrix(as.numeric(as.matrix(ovipo)),m1,n1)
    neggs = apply(dat,1,sum,na.rm=T) ## num of eggs
    neggs2 = round(ltb[1,1]*c1/nfem2)
    me_f = apply(dat,1,mean,na.rm=T) ## mean eggs on female
    sr= Rs2 ## sexual ratio -- constante por el momento
    mef_f = sr*me_f
    ies = rep(nfem2/ltb[1,1],length(x)) ##inmature estages survival
    nhem<-function(vec){n3=length(vec[!is.na(vec)]);return(n3)}
    nfem = apply(dat,1,nhem) ## num of females by day
    asurv = nfem/nfem[1]  ## adult survival
    lx= ies*asurv
    mx1 = mef_f*asurv ;mx2 = neggs*sr/lx/neggs2
    lx.mx1 = lx*mx1 ;lx.mx2=lx*mx2
    x.lx.mx1 = x*lx.mx1 ;x.lx.mx2 = x*lx.mx2
    S.mx1=sum(mx1)   ;S.mx2=sum(mx2)   ;S.lx.mx1=sum(lx.mx1)   ;S.lx.mx2=sum(lx.mx2)   ;S.x.lx.mx1=sum(x.lx.mx1)   ;S.x.lx.mx2= sum(x.lx.mx2)
    T01=S.x.lx.mx1/S.lx.mx1   ;T02=S.x.lx.mx2/S.lx.mx2
    
    r1=(log(S.lx.mx1))/T01   ;r2=(log(S.lx.mx2))/T02
    euler = (exp((-1)*r1*x))*lx.mx1  ## se tiene que jugar con el "r1" para que la suma de este vector sea 1
    Seuler = sum(euler)
    #r1 puede ser modificable hasta que Seuler sea 1
    T1=(log(S.lx.mx1))/r1   ;T2=(log(S.lx.mx2))/r2
    Dt1=log(2)/r1   ;Dt2=log(2)/r2
    lambda1=exp(r1)   ;lambda2=exp(r2)
    Ro1=S.lx.mx1   ;Ro2=S.lx.mx2
    GRR1=S.mx1   ;GRR2=S.mx2
    #data.frame(neggs,me_f,sr,mef_f,ies,nfem,asurv,lx,mx1,mx2,lx.mx1,lx.mx2,x.lx.mx1,x.lx.mx2)
    parameters=data.frame(c(r=r1,Ro=Ro1,GRR=GRR1,T=T1,lambda=lambda1,Dt=Dt1));colnames(parameters)="Parameters"
    #parameters=data.frame(c(r=r2,Ro=Ro2,GRR=GRR2,T=T2,lambda=lambda2,Dt=Dt2));colnames(parameters)="Parameters" ## este es el otro calculo segun la plantilla de excel
    return(list(parametro=parameters))
  }
}
