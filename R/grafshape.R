#' grafshape
#'
#' no details
#' 
#' @param model integer value, non linear model (position)
#' @param estshap list of some numeric vectors, the parameters values estimated
#' @param datashap data frame, estimated rates
#' @param qtt double, quantile of a normal distribution (p=0.95)
#' @param sdli double, standar deviation estimated
#' @param corrx numeric vector of two values, the axis range for X
#' @param corry numeric vector of two values, the axis range for y
#' @param mini double, minimun value of the temperature
#' @param maxi double, maximun value of the temperature
#' @param coefi numeric vector, linear regression parameters
#' @param limit logical; default is FALSE, when It considers the confidence lines to curve
#' @param tam integer value, extent of the plot
#' @param labx a character
#' @param laby a character
#' @param titulo a character
#' @param grises logical; default is FALSE
#' @param scaleY double, interval for the Y axis range
#' @param scaleX double, interval for the X axis range
#' @param est a character
#' @param estadios character vector, all states
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
grafshape<-function(model,estshap,datashap,qtt,sdli,corrx=NULL,corry=NULL,mini,maxi,coefi,limit,tam,labx=NULL,laby=NULL, titulo=NULL,grises=FALSE, scaleY, scaleX, est, estadios){
  if(grises==TRUE){ccol=c("gray20","gray30")}else{ccol=c(4,2)}
  
  data<-datashap
  
  estshap<-as.list(estshap)
  
  finalshap<-estshap
  
  for (i in names(finalshap)) {
    
    temp <- finalshap[[i]]
    
    storage.mode(temp) <- "double"
    
    assign(i, temp)
    
  }
  
  
  
  xl<-seq(mini,maxi,length=200)
  
  if(model==1)   {
    
    xl<-xl+273.15
    
    expre<-expression(((coefi[1]+coefi[2]*((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))) * (xl/(((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th))))) * exp((Ha/1.987) * ((1/((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))) - (1/xl))))/
                        
                        (1 + exp((Hl/1.987) * ((1/Tl) - (1/xl))) + exp((Hh/1.987) * ((1/Th) - (1/xl)))))
    
    yl<-eval(expre)
    
  }
  
  if(model==2)   {
    
    xl<-xl+273.15
    
    expre<-expression(((coefi[1]+coefi[2]*To) * (xl/(To)) * exp((Ha/1.987) * ((1/To) - (1/xl))))/
                        
                        (1 + exp((Hl/1.987) * ((1/Tl) - (1/xl))) + exp((Hh/1.987) * ((1/Th) - (1/xl)))))
    
    yl<-eval(expre)
    
  }
  
  if(model==3)   {
    
    xl<-xl+273.15
    
    expre<-expression((p * (xl/(To)) * exp((Ha/1.987) * ((1/To) - (1/xl))))/
                        
                        (1 + exp((Hl/1.987) * ((1/Tl) - (1/xl))) + exp((Hh/1.987) * ((1/Th) - (1/xl)))))
    
    yl<-eval(expre)
    
  }
  
  
  
  if(model==4){
    
    xl<-xl+273.15
    
    expre<-expression((p * (xl/(To)) * exp((Ha/1.987) * ((1/To) - (1/xl)))))
    
    yl<-eval(expre)
    
  }
  
  
  
  ######
  
  
  
  if(model==5){
    
    xl<-xl+273.15
    
    expre<-expression((p * (xl/(To)) * exp((Ha/1.987) * ((1/To) - (1/xl))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/xl)))))
    
    yl<-eval(expre)
    
  }
  
  
  
  if(model==6){
    
    xl<-xl+273.15
    
    expre<-expression((p * (xl/(To)) * exp((Ha/1.987) * ((1/To) - (1/xl))))/(1 + exp((Hh/1.987) * ((1/Th) - (1/xl)))))
    
    yl<-eval(expre)
    
  }
  
  ###### solo se aumenta la funcion considerando la transformacion a grados kelvin
  
  
  
  if(model==7){
    
    xl<-xl+273.15
    
    expre<-expression(((coefi[1]+coefi[2]*To) * (xl/(To)) * exp((Ha/1.987) * ((1/To) - (1/xl)))))
    
    yl<-eval(expre)
    
  }
  
  
  
  if(model==8){
    
    xl<-xl+273.15
    
    expre<-expression(((coefi[1]+coefi[2]*To) * (xl/(To)) * exp((Ha/1.987) * ((1/To) - (1/xl))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/xl)))))
    
    yl<-eval(expre)
    
  }
  
  
  
  if(model==9){
    
    xl<-xl+273.15
    
    expre<-expression(((coefi[1]+coefi[2]*To) * (xl/(To)) * exp((Ha/1.987) * ((1/To) - (1/xl))))/(1 + exp((Hh/1.987) * ((1/Th) - (1/xl)))))
    
    yl<-eval(expre)
    
  }
  
  
  
  if(model==10){
    
    xl<-xl+273.15
    
    expre<-expression((p * (xl/(298.16)) * exp((Ha/1.987) * ((1/298.16) - (1/xl))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/xl))) + exp((1/1.987) * ((1/1) - (1/xl)))))
    
    yl<-eval(expre)
    
  }
  
  
  
  ## nuevos modelos sharpe
  
  if(model==11){
    
    xl<-xl+273.15
    
    expre<-expression((p * (xl/298.16) * exp((Ha/1.987) * ((1/298.16) - (1/xl))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/xl))) + exp((Hh/1.987) * ((1/Th) - (1/xl)))))
    
    yl<-eval(expre)
    
  }
  
  
  if(model==12){
    
    xl<-xl+273.15
    
    expre<-expression(((coefi[1]+coefi[2]*298.16) * (xl/298.16) * exp((Ha/1.987) * ((1/298.16) - (1/xl))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/xl))) + exp((Hh/1.987) * ((1/Th) - (1/xl)))))
    
    yl<-eval(expre)
    
  }
  
  if(model==13){
    
    xl<-xl+273.15
    
    expre<-expression((p * (xl/298.16) * exp((Ha/1.987) * ((1/298.16) - (1/xl))))/(1 + exp((Hl/1.987) * ((1/Tl) - (1/xl)))))
    
    yl<-eval(expre)
    
  }
  
  if(model==14){
    
    xl<-xl+273.15
    
    expre<-expression((p * (xl/298.16) * exp((Ha/1.987) * ((1/298.16) - (1/xl))))/(1 + exp((Hh/1.987) * ((1/Th) - (1/xl)))))
    
    yl<-eval(expre)
    
  }
  
  
  
  
  ## fin de nuevos modelos sharpe
  
  
  if(model==15){
    
    expre<-expression(b*(xl-Tmin))
    
    yl<-eval(expre)
    
  }
  
  if(model==16)   {
    
    expre<-expression(b1*10^(-((((xl-b3)/(b3-b2)-(1/(1+0.28*b4+0.72*log(1+b4))))+exp(b4*((xl-b3)/(b3-b2)-(1/(1+0.28*b4+0.72*log(1+b4))))))/((1+b4)/(1+1.5*b4+0.39*b4^2)))^2)*(1-b5+b5*((((xl-b3)/(b3-b2)- (1/(1+0.28*b4+0.72*log(1+b4))))+exp(b4*((xl-b3)/(b3-b2)- (1/(1+0.28*b4+0.72*log(1+b4))))))/((1+b4)/(1+1.5*b4+0.39*b4^2)))^2)) #DM: Se agreg? un menos "-v^2"
    
    yl<-eval(expre)
    
  }
  
  
  
  if(model==17) {
    
    expre<-expression(Y*(exp(p*xl)-exp(p*Tmax-(Tmax-xl)/v)))
    
    yl<-eval(expre)
    
  }
  
  if(model==18){
    
    expre<-expression(alfa*((1/(1+k*exp(-p*xl)))-exp(-(Tmax-xl)/v)))
    
    yl<-eval(expre)
    
  }
  
  
  
  
  
  if(model==19){
    
    expre<-expression(aa*xl*(xl-To)*(Tmax-xl)^0.5)
    
    yl<-eval(expre)
    
  }
  
  if(model==20){
    
    expre<-expression(aa*xl*(xl-To)*(Tmax-xl)^d)
    
    yl<-eval(expre)
    
  }
  
  
  
  if(model==21) {
    
    expre<-expression((Rmax*(1+exp(k1+k2*(Topc))))/(1+exp(k1+k2*xl)))
    
    yl<-eval(expre)
    
  }
  
  
  
  if(model==22){
    
    expre<-expression(Y*(xl^2/(xl^2+d^2)-exp(-(Tmax-xl)/v)))
    
    yl<-eval(expre)
    
  }
  
  
  
  if(model==23) {
    
    expre<-expression(exp(p*xl)-exp(-(p*Tl-(xl-Tl))/dt)+L)
    
    yl<-eval(expre)
    
  }
  
  
  
  if(model==24) {
    
    expre<-expression(Inter + Slop*xl)
    
    yl<-eval(expre)
    
  }
  
  
  
  if(model==25) {
    
    expre<-expression(b1*exp(b2*xl))
    
    yl<-eval(expre)
    
  }
  
  
  
  if(model==26) {
    
    expre<-expression(sy*exp(b*(xl-Tb)-exp(b*(xl-Tb)/DTb)))
    
    yl<-eval(expre)
    
  }
  
  
  
  if(model==27) {
    
    expre<-expression(sy*exp(b*(xl-Tb)))
    
    yl<-eval(expre)
    
  }
  
  
  
  if(model==28) {
    
    expre<-expression(exp(b*(xl-Tmin))-1)
    
    yl<-eval(expre)
    
  }
  
  
  
  if(model==29) {
    
    expre<-expression(b*(xl-Tb)^2)
    
    yl<-eval(expre)
    
  }
  
  
  
  if(model==30) {
    
    expre<-expression(k/(1+exp(a-b*xl))) #DM: Se cambio "+b*xl" por "-b*xl"
    
    yl<-eval(expre)
    
  }
  
  if(model==31) {
    
    expre<-expression(R*exp((-1/2)*((xl-Tm)/To))^2) #DM: Se agrego "^2"
    
    yl<-eval(expre)
    
  }
  
  if(model==32) {
    
    expre<-expression(a*exp((-1/2)*(abs(xl-b)/c)^d)) #DM: Se cambio "abs((x-b)/c)" por "abs(x-b)/c"
    
    yl<-eval(expre)
    
  }
  
  if(model==33) {
    
    expre<-expression(Rmax*(exp(k1+k2*Topc))/(1+exp(k1+k2*xl)))
    
    yl<-eval(expre)
    
  }
  
  if(model==34) {
    
    expre<-expression(Y*((xl-Tb)^2/((xl-Tb)+d^2)-exp(-(Tmax-(xl-Tb))/v))) #DM: Se quito al cuadrado en el denominador
    
    yl<-eval(expre)
    
  }
  
  if(model==35) {
    
    expre<-expression(exp(p*xl)-exp(p*Tl-(Tl-xl)/dt)) #DM: Se omitio un menos "-(p" por "p"  y se invirtio "(x-Tl))/dt" por "(Tl-x)/dt"
    
    yl<-eval(expre)
    
  }
  
  if(model==36) {
    
    expre<-expression(P*(((xl-Tmin)/(Tmax-Tmin))^n)*(1-(xl-Tmin)/(Tmax-Tmin))^m)
    
    yl<-eval(expre)
    
  }
  
  if(model==37) {
    
    expre<-expression((P*(((xl-Tmin)/(Tmax-Tmin))^n)*(1-(xl-Tmin)/(Tmax-Tmin)))^m) #Se cambio "(P*" por "((P*"
    
    yl<-eval(expre)
    
  }
  
  if(model==38) {
    
    expre<-expression(a*((xl-Tmin)^n)*(Tmax-xl)^m)
    
    yl<-eval(expre)
    
  }
  
  if(model==39) {
    
    expre<-expression(P*(((xl-Tmin)/(Tmax-Tmin))^n)*(1-((xl-Tmin)/(Tmax-Tmin))^m)) #CAMBIO!
    
    yl<-eval(expre)
    
  }
  
  if(model==40) {
    
    expre<-expression(aa*(xl-To)*(Tmax-xl)^0.5)
    
    yl<-eval(expre)
    
  }
  
  if(model==41) {
    
    expre<-expression(aa*(xl-To)*(Tmax-xl)^(1/n))
    
    yl<-eval(expre)
    
  }
  
  if(model==42) {
    
    expre<-expression(aa*((xl-Tmin)^2)*(Tmax-xl))
    
    yl<-eval(expre)
    
  }
  
  if(model==43) {
    
    expre<-expression(2/(Dmin*(exp(K*(xl-Topt)) + exp((-lmda)*(xl-Topt)))))
    
    yl<-eval(expre)
    
  }
  
  if(model==44) {
    
    expre<-expression(xl*exp(a1-b1/xl)/(1 + exp(c1-d1/xl) + exp(f1-g1/xl)))
    
    yl<-eval(expre)
    
  }
  
  if(model==45) {
    
    expre<-expression((aa*(xl-Tmin)*(1-exp((b*(Tmax-xl)))))^2) #Agregar 1-exp
    
    yl<-eval(expre)
    
  }
  
  if(model==46) {
    
    expre<-expression(2/(Dmin*(exp(K*(xl-Topt)) + exp((-K)*(xl-Topt)))))
    
    yl<-eval(expre)
    
  }
  
  if(model==47) {
    
    expre<-expression(2*c/(a^(xl-Tm) + b^(Tm-xl)))
    
    yl<-eval(expre)
    
  }
  
  if(model==48) {
    
    expre<-expression(a0+a1*xl+a2*xl^2+a3*xl^3)
    
    yl<-eval(expre)
    
  }
  
  if(model==49) {
    
    expre<-expression(k*(1-exp((-a)*(xl-Tmin)))*(1-exp(b*(xl-Tmax)))/(1+exp((-r)*(xl-c))))
    
    yl<-eval(expre)
    
  }
  
  
  ## funciones adaptadas a senescencia
  
  if(model==50) {
    
    expre<-expression(c1/(1+exp(k1+k2*xl)))
    
    yl<-eval(expre)
    
  }
  
  if(model==51) {
    
    expre<-expression(c1/(1+exp(k1+k2*xl)) + c2/(1+exp(k1+k2*(2*To-xl))))
    
    yl<-eval(expre)
    
  }
  
  if(model==52) {
    
    expre<-expression(sy*exp(b*(xl-Tmin)-exp(b*Tmax - (Tmax-(xl-Tmin))/DTb)))
    
    yl<-eval(expre)
    
  }
  
  if(model==53) {
    
    expre<-expression(alph*(1/(1+k*exp(-b*(xl-Tmin))) - exp(-(Tmax-(xl-Tmin))/Dt)))
    
    yl<-eval(expre)
    
  }
  
  if(model==54) {
    
    expre<-expression(alph*(1/(1+k*exp(-b*xl)) - exp(-(Tmax-xl)/Dt)))
    
    yl<-eval(expre)
    
  }
  
  if(model==55) {
    
    expre<-expression(trid*( (xl^2)/(xl^2+D)  - exp(-(Tmax-xl)/Dt))) #Cambio de "x/Dt" por "x)/Dt"
    
    yl<-eval(expre)
    
  }
  
  if(model==56) {
    
    expre<-expression(trid*(((xl-Tmin)^2)/((xl-Tmin)^2 + D) - exp(-(Tmax-(xl-Tmin))/Dt)) + Smin)
    
    yl<-eval(expre)
    
  }
  
  if(model==57) {
    
    expre<-expression(rm*exp(-(0.5)*(-(xl-Topt)/Troh)^2))
    
    yl<-eval(expre)
    
  }
  
  if(model==58) {
    
    expre<-expression(exp(p*xl)-exp(p*Tl-(Tl-x)/dt) + lamb) #DM: Se cambi? "(p*Tl-(Tl-x))/dt)" por "p*Tl-(Tl-x)/dt)"  
    
    yl<-eval(expre)
    
  }
  
  if(model==59) {
    
    expre<-expression(c1/(1+exp(a+b*xl)))
    
    yl<-eval(expre)
    
  }
  
  
  
  par(cex=tam)
  
  if(is.null(corrx)){corrx=c(min(data[,1],0),max(data[,1])+max(data[,1])*0.2);corrx2=seq(0,10*round(corrx[2]/10,1),5)}else{corrx2=seq(corrx[1],corrx[2],scaleX)} ## cambio
  
  if(is.null(corry)){corry=c(0,max(data[,2],yl)+max(data[,2])*0.2);corry2=seq(0,round(max(corry),1),0.1)}else{corry2=seq(corry[1],corry[2],scaleY)}              ## cambio
  
  plot(data[,1],data[,2],frame=F,col=ccol[1],pch=19,xlim=corrx,ylim=corry,xlab=labx,ylab=laby,axes=F,xaxt = "n", main=titulo) ## 4
  n2= length(data[,1])
  #posic <- (1:7)[max(dataa[,2])==dataa[,2]][1] # excluye el ultimo valor para la regresion
  #dataaa <- dataa[1:posic,] # excluye el ultimo valor para la regresion
  
  posic=(1:n2)[max(data[,2])==data[,2]][1]
  modr= lm(data[1:posic,2]~data[1:posic,1])
  if((1:length(estadios))[estadios==est]<length(estadios)-1){
    abline(modr,col=ccol[2],lty=2,lwd=2)
  }
  #axis(1, xaxp=c(corrx,5))
  
  #axis(2,las=2)
  
  axis(1, corrx2)  ## cambio
  
  axis(2, corry2,las=2) ## cambio
  
  
  
  linf<-yl+sdli*qtt
  
  lsup<-yl-sdli*qtt
  
  if(limit=="yes"){
    
    if(model==3 || model==4 || model==5 || model==6 || model==7 || model==8 || model==9 || model==10 || model==11 || model==12 || model==13 || model==14){  ############ solo aumente el tipo de modelo
      
      lines(xl-273.15, yl, lwd=2,col=ccol[2]) ## 2
      
      lines(xl-273.15, linf, lwd=1,col=ccol[1],lty=2)
      
      lines(xl-273.15, lsup, lwd=1,col=ccol[1],lty=2)
      
    }else{
      
      lines(xl, yl, lwd=2,col=ccol[2]) ## 2
      
      lines(xl, linf, lwd=1,col=ccol[1],lty=2)
      
      lines(xl, lsup, lwd=1,col=ccol[1],lty=2)
      
    }}
  
  if(limit=="no") if(model==3 || model==4 || model==5  || model==6 || model==10 || model==11 || model==12 || model==13 || model==14) lines(xl-273.15, yl, lwd=2,col=ccol[2]) else lines(xl, yl, lwd=2,col=ccol[2]) ## 2 ..... solo aumente el tipo de modelo
  
  arrows(data[,1],(data[,2]-(data[,2]-data[,3])), data[,1],(data[,2]+(data[,4]-data[,2])), length=0.1,angle=90, code=3,col=ccol[1]) ## 4
  
  if(limit=="yes"){
    
    if(model==1 || model==2) {
      
      lines(xl-273.15, linf, lwd=1,col=ccol[1],lty=2) ## 4
      
      lines(xl-273.15, lsup, lwd=1,col=ccol[1],lty=2) ## 4
      
      yli<-coefi[1]+coefi[2]*(xl)
      
      lines(xl-273.15, yl, lwd=2,col=ccol[2]) ## 2
      
      #lines(xl-273.15, yli, col=ccol[2],lty=2,lwd=2) ## 2
      
    }}
  
  if(limit=="no"){
    
    if(model==1 || model==2 || model==7 || model==8 || model==9 || model==12 ) {
      
      yli<-coefi[1]+coefi[2]*(xl)
      
      lines(xl-273.15, yl, lwd=2,col=ccol[2]) ## 2
      
      #lines(xl-273.15, yli, col=ccol[2],lty=2,lwd=2) ## 2
      
    }}
  
  if(model==1) {
    
    points(((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th)))-273.15,coefi[1]+coefi[2]*((Hl-Hh)/(1.987*log(-Hl/Hh)+(Hl/Tl)-(Hh/Th))),pch=22,cex=2)
    
  }
  
  if(model==2 || model==7 || model==8 || model==9 || model==12) {
    if(model==12){To=298.16}
    points(To-273.15,coefi[1]+coefi[2]*To,pch=22,cex=2)
    
  }
  
}
