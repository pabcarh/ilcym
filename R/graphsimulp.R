#' graphsimulp
#'
#' no details
#' 
#' @param pars data frame, temperatures and their estimations of life parameters
#' @param modls vector of integer values, simple non linear model (positions)
#' @param direc a character
#' @param ecua logical; default is FALSE
#' @param nomec logical; default is FALSE
#' @param ejex double, the axis limit for X
#' @param ejey double, the axis limit for Y
#' @param tecu integer value, size of the equation
#' @param grises logical; default is FALSE
#' @param tit a character
#' @param ax integer value, size of axis font for X
#' @param ay integer value, size of axis font for Y
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
graphsimulp<-function(pars,modls=c(1,1,1,1,1,1),direc,ecua,nomec,ejex,ejey,tecu,grises,tit,ax,ay){
  
  if(grises==TRUE){ccol=c("gray20","gray50")}else{ccol=c("royalblue","red")}        

  etip=c("Rm","Ro","GRR","GL","Lambda","Dt")
  
  nfile=etip[!is.na(modls)]
  
  NMS=paste(nfile,collapse="_")
  
  pathIm = paste(direc,NMS,".jpg",sep="")
  
  mx<-mmx<-max(pars[,1],na.rm = TRUE);mnx<-min(pars[,1],na.rm = TRUE)
  if(mx<=40){mx=40}else{mx=mx+10}
  N=length(modls[!is.na(modls)]) 
  if(N==1){par(mfrow=c(1,1));  if(is.na(ejex)){ejex=-1.73};  if(is.na(ejey)){ejey=-1.61}  }
  if(N==2){par(mfrow=c(2,1));  if(is.na(ejex)){ejex=-0.69};  if(is.na(ejey)){ejey=-1.60}  }
  if(N==3 || N==4){par(mfrow=c(2,2));  if(is.na(ejex)){ejex=-0.86};  if(is.na(ejey)){ejey=-0.84}  }
  if(N==5 || N==6){par(mfrow=c(3,2));  if(is.na(ejex)){ejex=-0.66};  if(is.na(ejey)){ejey=-1.25}  }
  vars=c(2,3,4,5,6,7)
  nomv=c(Temper="Temperature",r = expression(Intrinsic~rate~group("(",r[m],")")),Ro = "Net reproduction rate (Ro)",GRR = "Gross reproduction rate (GRR)",T = "Generation length in days (GL)",lambda = expression(Finite~rate~of~increase~group("(",symbol("l"),")") ),Dt = "Doubling time (Dt)") ## nombre de los parametros
  modelss=c("cubic","quadratic","logarithmic","exponential")  ## tipo de modelos propuestos

  modsel<-matrix(NA,6,4)
  colnames(modsel)<-c("R2","R2_Adj","AIC","Deviance")
  sal2<-et<-list(1);ubix=20

  for(k in vars) ## gr?fico por cada par?metro
  { 
    if(ecua==TRUE){if(is.na(modls[k-1])){next}} ## este codigo hace que solo ejecute los parametros que tienen un modelo escogido
    my=max(pars[,k],na.rm = TRUE);mny=min(pars[,k],na.rm = TRUE)
    if(mny<(my-mny) & my>0){mny=0}
    plot(pars[,1],pars[,k],xlim=c(0,mx),ylim=c(mny-0.05*(abs(mny)),my+0.13*(abs(my))),xlab="Temperature ?C",ylab=nomv[k],pch=19,col=ccol[1],axes=FALSE)  ## 1
    axis(1,at=seq(0,mx,by=10),line=ejex,cex.axis=ax) ## falta incrementar los tama?os de los subtitulos a 12
    axis(2,las=2,line=ejey,cex.axis=ay)
    x=pars[,1]
    y=pars[,k]
    if(ecua==TRUE)
    {
      if(modls[k-1]==1){f1=as.formula(y ~ b0 + b1*x + b2*x^2 + b3*x^3);inis=list(b0=1,b1=1,b2=2,b3=2);ff <- function(x){b0 + b1*x + b2*x^2 + b3*x^3}}
      if(modls[k-1]==2){f1=as.formula(y ~ b0 + b1*x + b2*x^2);inis=list(b0=1,b1=1,b2=2);ff <- function(x){b0 + b1*x + b2*x^2}}
      if(modls[k-1]==3){f1=as.formula(y ~ b0 + b1*x + b2*log(x));inis0=coef(lm(y~x+I(log(x))));names(inis0)=NULL;inis=list(b0=inis0[1],b1=inis0[2],b2=inis0[3]);ff <- function(x){b0 + b1*x + b2*log(x)}}
      if(modls[k-1]==4){f1=as.formula(y ~ b0 + b1*x + b2*exp(x));inis0=coef(lm(y~x+I(exp(x))));names(inis0)=NULL;inis=list(b0=inis0[1],b1=inis0[2],b2=inis0[3]);ff <- function(x){b0 + b1*x + b2*exp(x)}}

      if(tit==TRUE){title(paste("Using the",modelss[modls[k-1]],"model"))}
      
      out <- nls(f1, start = inis,trace = TRUE)
      yl=fitted(out)
      sqe=sum(residuals(out)^2)
      sal <- coef(out);sal2[[k-1]]=sal

      r<-1-sqe/sum((y-mean(y))^2) # este es el R^2
      r_ajus<- 1 - ((length(x) - 1) / (length(x) - length(sal))) * (1-r)  # este es el R^2 adjus
      AC<-AIC(out) # este es el AIC 
      Dev=deviance(out) # este es el indicador  Deviance

      modsel0<-c(R2=round(r,3),R2_Adj=round(r_ajus,3),AIC=round(AC,3),Deviance=round(Dev,3))
      modsel[k-1,]<-modsel0

      if((k-1)==1)
      {
        if(modls[k-1]==1){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(bolditalic(r[m]) == list(sal1) + list(sal2)*T + list(sal3)*T^2 + list(sal4)*T^3,list(sal1=sal[1],sal2=sal[2],sal3=sal[3],sal4=sal[4])),cex=tecu))}
        if(modls[k-1]==2){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(bolditalic(r[m]) == list(sal1) + list(sal2)*T + list(sal3)*T^2,list(sal1=sal[1],sal2=sal[2],sal3=sal[3])),cex=tecu))}
        if(modls[k-1]==3){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(bolditalic(r[m]) == list(sal1) + list(sal2)*T + list(sal3)*log(T),list(sal1=sal[1],sal2=sal[2],sal3=sal[3])),cex=tecu))}
        if(modls[k-1]==4){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(bolditalic(r[m]) == list(sal1) + list(sal2)*T + list(sal3)*e^(T),list(sal1=sal[1],sal2=sal[2],sal3=sal[3])),cex=tecu))}
      }

      if((k-1)==2 || (k-1)==3 || (k-1)==4 || (k-1)==6)
      {
        if(modls[k-1]==1){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(bolditalic(list(etip1)) == list(sal1) + list(sal2)*T + list(sal3)*T^2 + list(sal4)*T^3,list(etip1=etip[k-1],sal1=sal[1],sal2=sal[2],sal3=sal[3],sal4=sal[4])),cex=tecu))}
        if(modls[k-1]==2){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(bolditalic(list(etip1)) == list(sal1) + list(sal2)*T + list(sal3)*T^2,list(etip1=etip[k-1],sal1=sal[1],sal2=sal[2],sal3=sal[3])),cex=tecu))}
        if(modls[k-1]==3){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(bolditalic(list(etip1)) == list(sal1) + list(sal2)*T + list(sal3)*log(T),list(etip1=etip[k-1],sal1=sal[1],sal2=sal[2],sal3=sal[3])),cex=tecu))}
        if(modls[k-1]==4){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(bolditalic(list(etip1)) == list(sal1) + list(sal2)*T + list(sal3)*e^(T),list(etip1=etip[k-1],sal1=sal[1],sal2=sal[2],sal3=sal[3])),cex=tecu))}
      } 
      if((k-1)==5)
      {
        if(modls[k-1]==1){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(symbol(list(etip1)) == list(sal1) + list(sal2)*T + list(sal3)*T^2 + list(sal4)*T^3,list(etip1=etip[k-1],sal1=sal[1],sal2=sal[2],sal3=sal[3],sal4=sal[4])),cex=tecu))}
        if(modls[k-1]==2){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(symbol(list(etip1)) == list(sal1) + list(sal2)*T + list(sal3)*T^2,list(etip1=etip[k-1],sal1=sal[1],sal2=sal[2],sal3=sal[3])),cex=tecu))}
        if(modls[k-1]==3){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(symbol(list(etip1)) == list(sal1) + list(sal2)*T + list(sal3)*log(T),list(etip1=etip[k-1],sal1=sal[1],sal2=sal[2],sal3=sal[3])),cex=tecu))}
        if(modls[k-1]==4){et[[k-1]]=expression(text(ubix, my+0.09*(my),substitute(symbol(list(etip1)) == list(sal1) + list(sal2)*T + list(sal3)*e^(T),list(etip1=etip[k-1],sal1=sal[1],sal2=sal[2],sal3=sal[3])),cex=tecu))}
      }
      if(nomec==TRUE){eval(et[[k-1]])}
      for (i in names(sal))
      { 
        temp <- sal[i] 
        storage.mode(temp) <- "double"
        assign(i, temp)
      }
      curve(ff,add=TRUE,from=mnx,to=mmx,col=ccol[2],lwd=2)  ## 2
    }
  }

  if(ecua==TRUE) 
  {
    plot(1:12,1:12,type="n",axes=FALSE,xlab="",ylab="")
    grid(0,6,lty=1)
    ubix=8.5;ubix2=8.5
    lines(c(4.5,4.5),c(0,13),lty=1,lwd=2,col="gray90",type = "l")
    box(lty = "solid", col = "gray30")

    if(!is.na(modls[1])){
      sal=sal2[[1]]
      my=(12/1.09);k=2;if(modls[1]==1){ubix=8.5};if(modls[1]==2){ubix=7.6};if(modls[1]==3 || modls[1]==4){ubix=6.7}
      text(1.8,12,nomv[2], col="gray40") ## para rm
      eval(et[[1]])
      text(ubix2-3, my*(1.09)-0.94,  substitute(bolditalic(R^{2})==list(r2),list(r2=modsel[1,1])))
      text(ubix2-1.2, my*(1.09)-0.97,  substitute(bolditalic(R^{2}~Adj)==list(r2a),list(r2a=modsel[1,2])))
      text(ubix2+0.6, my*(1.09)-1,  substitute(bolditalic(AIC)==list(aic),list(aic=modsel[1,3])))
      text(ubix2+2.6, my*(1.09)-1,  substitute(bolditalic(Deviance)==list(dev),list(dev=modsel[1,4])))
    }
    if(!is.na(modls[2])){
      sal=sal2[[2]]
      my=(12/1.09)-1.8;k=3;if(modls[2]==1){ubix=8.5};if(modls[2]==2){ubix=7.6};if(modls[2]==3 || modls[2]==4){ubix=6.7}
      text(2.2,10,nomv[3], col="gray40") ## para Ro
      eval(et[[2]])
      text(ubix2-3, my*(1.09)-0.94,  substitute(bolditalic(R^{2})==list(r2),list(r2=modsel[2,1])))
      text(ubix2-1.2, my*(1.09)-0.97,  substitute(bolditalic(R^{2}~Adj)==list(r2a),list(r2a=modsel[2,2])))
      text(ubix2+0.6, my*(1.09)-1,  substitute(bolditalic(AIC)==list(aic),list(aic=modsel[2,3])))
      text(ubix2+2.6, my*(1.09)-1,  substitute(bolditalic(Deviance)==list(dev),list(dev=modsel[2,4])))
    }
    if(!is.na(modls[3])){
      sal=sal2[[3]]
      my=(12/1.09)-3.6;k=4;if(modls[3]==1){ubix=8.5};if(modls[3]==2){ubix=7.6};if(modls[3]==3 || modls[3]==4){ubix=6.7}
      text(2.5,8,nomv[4], col="gray40")
      eval(et[[3]])
      text(ubix2-3, my*(1.09)-0.94,  substitute(bolditalic(R^{2})==list(r2),list(r2=modsel[3,1])))
      text(ubix2-1.2, my*(1.09)-0.97,  substitute(bolditalic(R^{2}~Adj)==list(r2a),list(r2a=modsel[3,2])))
      text(ubix2+0.6, my*(1.09)-1,  substitute(bolditalic(AIC)==list(aic),list(aic=modsel[3,3])))
      text(ubix2+2.6, my*(1.09)-1,  substitute(bolditalic(Deviance)==list(dev),list(dev=modsel[3,4])))
    }
    if(!is.na(modls[4])){
      sal=sal2[[4]]
      my=(12/1.09)-5.4;k=5;if(modls[4]==1){ubix=8.5};if(modls[4]==2){ubix=7.6};if(modls[4]==3 || modls[4]==4){ubix=6.7}
      text(2.5,6,nomv[5], col="gray40")
      eval(et[[4]])
      text(ubix2-3, my*(1.09)-0.94,  substitute(bolditalic(R^{2})==list(r2),list(r2=modsel[4,1])))
      text(ubix2-1.2, my*(1.09)-0.97,  substitute(bolditalic(R^{2}~Adj)==list(r2a),list(r2a=modsel[4,2])))
      text(ubix2+0.6, my*(1.09)-1,  substitute(bolditalic(AIC)==list(aic),list(aic=modsel[4,3])))
      text(ubix2+2.6, my*(1.09)-1,  substitute(bolditalic(Deviance)==list(dev),list(dev=modsel[4,4])))
    }
    if(!is.na(modls[5])){
      sal=sal2[[5]]
      my=(12/1.09)-7.2;k=6;if(modls[5]==1){ubix=8.5};if(modls[5]==2){ubix=7.6};if(modls[5]==3 || modls[5]==4){ubix=6.7}
      text(2.4,4,nomv[6], col="gray40") ## para lambda
      eval(et[[5]])
      text(ubix2-3, my*(1.09)-0.94,  substitute(bolditalic(R^{2})==list(r2),list(r2=modsel[5,1])))
      text(ubix2-1.2, my*(1.09)-0.97,  substitute(bolditalic(R^{2}~Adj)==list(r2a),list(r2a=modsel[5,2])))
      text(ubix2+0.6, my*(1.09)-1,  substitute(bolditalic(AIC)==list(aic),list(aic=modsel[5,3])))
      text(ubix2+2.6, my*(1.09)-1,  substitute(bolditalic(Deviance)==list(dev),list(dev=modsel[5,4])))
    }
    if(!is.na(modls[6])){
      sal=sal2[[6]]
      my=(12/1.09)-9;k=7;if(modls[6]==1){ubix=8.5};if(modls[6]==2){ubix=7.6};if(modls[6]==3 || modls[6]==4){ubix=6.7}
      text(2,2,nomv[7], col="gray40") ## para Dt
      eval(et[[6]])
      text(ubix2-3, my*(1.09)-0.94,  substitute(bolditalic(R^{2})==list(r2),list(r2=modsel[6,1])))
      text(ubix2-1.2, my*(1.09)-0.97,  substitute(bolditalic(R^{2}~Adj)==list(r2a),list(r2a=modsel[6,2])))
      text(ubix2+0.6, my*(1.09)-1,  substitute(bolditalic(AIC)==list(aic),list(aic=modsel[6,3])))
      text(ubix2+2.6, my*(1.09)-1,  substitute(bolditalic(Deviance)==list(dev),list(dev=modsel[6,4])))
    }
  }
  return(list(coefs=sal2,modsel=modsel,et=et, pathIm=pathIm))
}
