#' coef_mort
#'
#' no details
#' 
#' @param proc a character, type of variable (mortal, taza)
#' @param modelm a character, non linear model
#' @param estimor numeric vector, estimated coefficients of the nonlinear model
#' @param stdmortg numeric vector, standar deviation estimated
#' @param modelo a character, non linear model
#' @param gg a character, equation
#' @param datm data frame, temperatures and mortality
#' @param alg a character, algorithm
#' @param pesos numeric vector
#' @param weights numeric vector
#' @param g formula, equation
#' @param modelim character vector, all stages of life

#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
coef_mort<-function(proc,modelm,estimor,stdmortg,modelo,gg,datm,alg,pesos,weights,g,modelim)
{
  x<-datm[,1]
  y<-datm[,2]
  ind <- as.list(estimor)
  for (i in names(ind))
  {
    temp <- estimor[[i]]
    storage.mode(temp) <- "double"
    assign(i, temp)
  }
  if(alg=="Newton")
  {
    if(modelm==1 ||modelm==2 || modelm==3 || modelm==4 || modelm==5 || modelm==6)
    {
      yl<-predict(modelo)
      art=0.00000000000000000000000000001
      sqe <- sum(residuals(modelo)^2)+art
      qq<-summary(modelo)
      tabla<-qq$"parameters"
      stderror<-tabla[,2]  #Agrege DM
      aic<-AIC(modelo) #Quite ,k=length(estimor)  DM
      if(pesos=="WLS") r<-1-sum(residuals(modelo)^2)/sum(weights*(y-mean(y))^2)
      if(pesos=="LS")  r<-1-sum(residuals(modelo)^2)/sum((y-mean(y))^2)
      r_ajus<- 1 - ((length(x) - 1) / (length(x) - length(estimor))) * (1-r)
      if(pesos=="WLS") MSC<-log(sum(weights*(yl-mean(yl))^2)/sum(weights*(y-mean(y))^2))-2*length(estimor)/length(x)
      if(pesos=="LS") MSC<-log(sum((yl-mean(yl))^2)/sum((y-mean(y))^2))-2*length(estimor)/length(x)
      #var <- (sum((y-yl)^2)/(length(x) - length(estimor)))
      anva.1 <- c(length(coef(modelo))-1,length(x)-length(coef(modelo)),length(x)-1)
      anva.2 <- c((sum(y^2)-length(y)*mean(y)^2)-sqe,sqe,sum(y^2)-length(y)*mean(y)^2)
      anva.3 <- c(anva.2[1]/anva.1[1],anva.2[2]/anva.1[2],NA)
      anva.4 <- c(anva.3[1]/anva.3[2],NA,NA)
      anva.5 <- c(1-pf(anva.4[1],anva.1[1],anva.1[2]),NA,NA)
      anva   <- cbind(anva.1,round(anva.2,4),round(anva.3,4),round(anva.4,4),round(anva.5,4))
      rownames(anva) <- c("Model","Error","Total")
      colnames(anva) <- c("DF","SS","MS","Fc","p")
      if(proc=="mortal"){
        cat("MORTALITY FOR TEMPERATURE\n")
        dat <- data.frame(T = datm[,1], Mortality = datm[,2])}
      if(proc=="taza"){
        cat("FECUNDYTY FOR TEMPERATURE\n")
        dat <- data.frame(T = datm[,1], Fecundity = datm[,2])}
      print(dat)
      cat("\nNONLINEAR REGRESSION MODEL\n")
      cat("\nMethod:", alg)
      cat("\nFormula: ", g,"\n")
      cat("\nParameters\n")
      print(round(qq$"parameters",5))
      cat("\nAnalysis of Variance\n")
      print(anva,na.print = "")
      #n<-length(x); k<-length(estimor);  if(n/k<40) {aic<-aic+2*k*(k+1)/(n-k-1); nombre="AICc"} else nombre="AIC"  #Nuevo DM
      FRAMESEL<-data.frame(R2=round(r,4),R2_Adj=round(r_ajus,4),SSR=round(sqe,4),AIC=round(aic,4),MSC=round(MSC,4))
      #colnames(FRAMESEL)[4]<-nombre # Nuevo DM
      rownames(FRAMESEL)<-c("")
      cat("\nSelection criteria")
      print(t(FRAMESEL))
      para<-t(matrix(paste(round(estimor,5),"(","?",round(as.numeric(tabla[,2]),5),")")))
      colnames(para)<-names(estimor)
      param<-data.frame(g,para,d.f.=paste(anva.1[1],",",anva.1[2]),F=round(anva.4[1],3),P.value=round(anva.5[1],3),FRAMESEL)
    }
  }
  if(alg=="Marquardtr")
  {
    if(modelm==9)   {
      expre<-expression(y0 + a * exp(-0.5 * ((x-x0)/b)^2))
      yl<-eval(expre)
      gg<-y~y0 + a * exp(-0.5 * ((x-x0)/b)^2)
    }
    if(modelm==7){
      expre<-expression(1/(1+a*exp(-b*((x-c)/d)^2)))
      yl<-eval(expre)
      gg<-y~1/(1+a*exp(-b*((x-c)/d)^2))
    }
    if(modelm==8) {
      expre<-expression(a*exp(-b*((x-xo)/c)^2))
      yl<-eval(expre)
      gg<-y~a*exp(-b*((x-xo)/c)^2)
    }
    if(modelm==10) {
      expre<-expression(y0 + a * exp(-0.5 * (log(abs(x/x0))/b)^2))
      yl<-eval(expre)
      gg<-y~y0 + a * exp(-0.5 * (log(abs(x/x0))/b)^2)
    }
    
    if(modelm==11) {
      expre<-expression(b1+b2*x+b3*x^d)
      yl<-eval(expre)
      gg<-y~b1+b2*x+b3*x^d
    }
    
    if(modelm==12) {
      expre<-expression(exp(b1+b2*x+b3*x^2))
      yl<-eval(expre)
      gg<-y~exp(b1+b2*x+b3*x^2)
    }
    
    if(modelm==13) {
      expre<-expression(1-(b4/(1+b5*exp(b1+b2*x+b3*x^2))))
      yl<-eval(expre)
      gg<-y~1-(b4/(1+b5*exp(b1+b2*x+b3*x^2)))
    }
    
    
    if(modelm==14) {
      expre<-expression(exp(b1+b2*x+b3*sqrt(x)))
      yl<-eval(expre)
      gg<-y~exp(b1+b2*x+b3*sqrt(x))
    }
    
    
    if(modelm==15) {
      expre<-expression(1-(b4/(1+b5*exp(b1+b2*x+b3*sqrt(x)))))
      yl<-eval(expre)
      gg<-y~1-(b4/(1+b5*exp(b1+b2*x+b3*sqrt(x))))
    }
    
    if(modelm==16) {
      expre<-expression(exp(b1+b2*x+b3*(1/sqrt(x))))
      yl<-eval(expre)
      gg<-y~exp(b1+b2*x+b3*(1/sqrt(x)))
    }
    
    if(modelm==17) {
      expre<-expression(1-(b4/(1+b5*exp(b1+b2*x+b3*(1/sqrt(x))))))
      yl<-eval(expre)
      gg<-y~1-(b4/(1+b5*exp(b1+b2*x+b3*(1/sqrt(x)))))
    }
    
    if(modelm==18) {
      expre<-expression(exp(b1+b2*x+b3*(1/x)))
      yl<-eval(expre)
      gg<-y~exp(b1+b2*x+b3*(1/x))
    }
    
    if(modelm==19) {
      expre<-expression(1-(b4/(1+b5*exp(b1+b2*x+b3*(1/x)))))
      yl<-eval(expre)
      gg<-y~1-(b4/(1+b5*exp(b1+b2*x+b3*(1/x))))
    }
    
    if(modelm==20) {
      expre<-expression(exp(b1+b2*x+b3*x^d))
      yl<-eval(expre)
      gg<-y~exp(b1+b2*x+b3*x^d)
    }
    
    if(modelm==21) {
      expre<-expression(1-(b4/(1+b5*exp(b1+b2*x+b3*x^d))))
      yl<-eval(expre)
      gg<-y~1-(b4/(1+b5*exp(b1+b2*x+b3*x^d)))
    }
    
    if(modelm==22) {
      expre<-expression(exp(b1+b2*x+b3*log(x)))
      yl<-eval(expre)
      gg<-y~exp(b1+b2*x+b3*log(x))
    }
    
    if(modelm==23) {
      expre<-expression(1-(b4/(1+b5*exp(b1+b2*x+b3*log(x)))))
      yl<-eval(expre)
      gg<-y~1-(b4/(1+b5*exp(b1+b2*x+b3*log(x))))
    }
    
    if(modelm==24) {
      expre<-expression(1-rm*exp((-0.5)*(-(x-Topt)/Troh)^2))
      yl<-eval(expre)
      gg<-y~1-rm*exp((-0.5)*(-(x-Topt)/Troh)^2)
    }
    
    if(modelm==25) {
      expre<-expression(1-rm*exp((-0.5)*(-(log(x)-log(Topt))/Troh)^2))
      yl<-eval(expre)
      gg<-y~1-rm*exp((-0.5)*(-(log(x)-log(Topt))/Troh)^2)
    }
    
    if(modelm==26) {
      expre<-expression(1 - 1/(exp((1+exp(-(x-Topt)/B))*(1+exp(-(Topt-x)/B))*H)))
      yl<-eval(expre)
      gg<-y~1 - 1/(exp((1+exp(-(x-Topt)/B))*(1+exp(-(Topt-x)/B))*H))
    }
    
    if(modelm==27) {
      expre<-expression(1 - 1/(exp((1+exp(-(x-Tl)/B))*(1+exp(-(Th-x)/B))*H)))
      yl<-eval(expre)
      gg<-y~1 - 1/(exp((1+exp(-(x-Tl)/B))*(1+exp(-(Th-x)/B))*H))
    }
    
    if(modelm==28) {
      expre<-expression(1 - 1/(exp((1+exp(-(x-Topt)/Bl))*(1+exp(-(Topt-x)/Bh))*H)))
      yl<-eval(expre)
      gg<-y~1 - 1/(exp((1+exp(-(x-Topt)/Bl))*(1+exp(-(Topt-x)/Bh))*H))
    }
    
    if(modelm==29) {
      expre<-expression(1 - 1/(exp((1+exp(-(x-Tl)/Bl))*(1+exp(-(Th-x)/Bh))*H)))
      yl<-eval(expre)
      gg<-y~1 - 1/(exp((1+exp(-(x-Tl)/Bl))*(1+exp(-(Th-x)/Bh))*H))
    }
    
    
    if(modelm==30) {
      expre<-expression(1 - H/(exp(1+exp(-(x-Topt)/B))*(1+exp(-(Topt-x)/B))))
      yl<-eval(expre)
      gg<-y~1 - H/(exp(1+exp(-(x-Topt)/B))*(1+exp(-(Topt-x)/B)))
    }
    
    if(modelm==31) {
      expre<-expression(1 - H/(exp(1+exp(-(x-Tl)/B))*(1+exp(-(Th-x)/B))))
      yl<-eval(expre)
      gg<-y~1 - H/(exp(1+exp(-(x-Tl)/B))*(1+exp(-(Th-x)/B)))
    }
    
    if(modelm==32) {
      expre<-expression(1 - H/(exp(1+exp(-(x-Topt)/Bl))*(1+exp(-(Topt-x)/Bh))))
      yl<-eval(expre)
      gg<-y~1 - H/(exp(1+exp(-(x-Topt)/Bl))*(1+exp(-(Topt-x)/Bh)))
    }
    
    if(modelm==33) {
      expre<-expression(1 - H/(exp(1+exp(-(x-Tl)/Bl))*(1+exp(-(Th-x)/Bh))))
      yl<-eval(expre)
      gg<-y~1 - H/(exp(1+exp(-(x-Tl)/Bl))*(1+exp(-(Th-x)/Bh)))
    }
    
    if(modelm==34) {
      expre<-expression(Bm + ((1-1/(1+exp((Hl/1.987)*((1/Tl)-(1/(x+273.15))))+exp((Hh/1.987)*((1/Th)-(1/(x+273.15))))))*(1-Bm))) #DM: se agrego 2 parentesis "Bm + (...)*"
      yl<-eval(expre)
      gg<-y~Bm + ((1-1/(1+exp((Hl/1.987)*((1/Tl)-(1/(x+273.15))))+exp((Hh/1.987)*((1/Th)-(1/(x+273.15))))))*(1-Bm)) #DM: se agrego 2 parentesis "Bm + (...)*"
    }
    
    if(modelm==35) {
      expre<-expression((1 - exp(-(exp(a1+b1*x)))) + (1 - exp(-(exp(a2+b2*x)))))
      yl<-eval(expre)
      gg<-y~(1 - exp(-(exp(a1+b1*x)))) + (1 - exp(-(exp(a2+b2*x))))
    }
    
    if(modelm==36) {
      expre<-expression((w-x)^(-1))
      yl<-eval(expre)
      gg<-y ~ (w-x)^(-1)
    }
    
    if(modelm==37) {
      expre<-expression(a1*exp(b1*x) + a2*exp(b2*x))
      yl<-eval(expre)
      gg<-y ~ a1*exp(b1*x) + a2*exp(b2*x)
    }
    
    if(modelm==38) {
      expre<-expression(a1*exp(b1*x) + a2*exp(b2*x) + c1)
      yl<-eval(expre)
      gg<-y ~ a1*exp(b1*x) + a2*exp(b2*x) + c1
    }
    
    
    if(modelm==39) {
      expre<-expression(a*(abs(x-b))^nn)
      yl<-eval(expre)
      gg<-y ~ a*(abs(x-b))^nn
    }
    
    
    if(modelm==40) {
      expre<-expression(a*x*(x-To)*(Tl-x)^(1/d))
      yl<-eval(expre)
      gg<-y ~ a*x*(x-To)*(Tl-x)^(1/d)
    }
    
    
    if(modelm==41) {
      expre<-expression(exp(a*x*(x-To)*(Tl-x)^(1/d)))
      yl<-eval(expre)
      gg<-y ~ exp(a*x*(x-To)*(Tl-x)^(1/d))
    }
    
    if(modelm==42) {
      expre<-expression(a*((x-Tmin)^n)*(Tmax-x)^m)
      yl<-eval(expre)
      gg<-y ~ a*((x-Tmin)^n)*(Tmax-x)^m
    }
    
    if(modelm==43) {
      expre<-expression(1/((Dmin/2) * (exp(k*(x-Tp)) + exp(-(x-Tp)*lamb))))
      yl<-eval(expre)
      gg<-y ~ 1/((Dmin/2) * (exp(k*(x-Tp)) + exp(-(x-Tp)*lamb)))
    }
    
    if(modelm==44) {
      expre<-expression(a*(1-exp(-(x-Tl)/B))*(1-exp(-(Th-x)/B)))
      yl<-eval(expre)
      gg<-y ~ a*(1-exp(-(x-Tl)/B))*(1-exp(-(Th-x)/B))
    }
    
    if(modelm==45) {
      expre<-expression(exp(a*(1-exp(-(x-Tl)/Bl))*(1-exp(-(Th-x)/Bh))))
      yl<-eval(expre)
      gg<-y ~ exp(a*(1-exp(-(x-Tl)/Bl))*(1-exp(-(Th-x)/Bh)))
    }
    
    sqe <- sum((y-yl)^2)
    if(length(x) == length(estimor)) stop("number of parameters = number of temperatures") #DM se agrego!
    var <- (sum((y-yl)^2)/(length(x) - length(estimor)))
    stderror<-sqrt(var*stdmortg)
    #if(length(x) == length(estimor)) stderror<-rep(NA,length(estimor)) #JC
    names(stderror)<-names(estimor) # Agregue DM
    tvalues<-estimor/stderror
    tvalta<-qt(0.025,length(x) - length(estimor))
    pvalues<-1-pt(as.numeric(tvalues),length(x) - length(estimor))
    r<-1-sqe/sum((y-mean(y))^2)
    r_ajus<- 1 - ((length(x) - 1) / (length(x) - length(estimor))) * (1-r)
    #AC<-length(x)*log(sqe/length(x))+2*length(estimor) #Agregue /length(x) DM
    AC<-length(x)*(log(2*pi*sqe/length(x))+1)+2*(length(estimor)+1) #Cambio segun Nonlinear Regression with R, pag 105. DM
    MSC<-log(sum((yl-mean(yl))^2)/sum((y-mean(y))^2))-2*length(estimor)/length(x)
    anva.1 <- c(length(estimor)-1,length(x)-length(estimor),length(x)-1)
    anva.2 <- c((sum(y^2)-length(y)*mean(y)^2)-sqe,sqe,sum(y^2)-length(y)*mean(y)^2)
    anva.3 <- c(anva.2[1]/anva.1[1],anva.2[2]/anva.1[2],NA)
    anva.4 <- c(anva.3[1]/anva.3[2],NA,NA)
    anva.5 <- c(1-pf(anva.4[1],anva.1[1],anva.1[2]),NA,NA)
    anva   <- cbind(anva.1,round(anva.2,4),round(anva.3,4),round(anva.4,4),round(anva.5,4))
    rownames(anva) <- c("Model","Error","Total")
    colnames(anva) <- c("DF","SS","MS","Fc","p")
    if(proc=="mortal"){
      cat("MORTALITY FOR TEMPERATURE\n")
      dat <- data.frame(T = datm[,1], Mortality = datm[,2])}
    if(proc=="taza"){
      cat("FECUNDYTY FOR TEMPERATURE\n")
      dat <- data.frame(T = datm[,1], Fecundity = datm[,2])}
    print(dat)
    cat("\nNONLINEAR REGRESSION MODEL\n")
    cat("\nMethod:", alg)
    cat("\nFormula: ", g,"\n")
    cat("\nParameters\n")
    estshap2<-data.frame(matrix(round(estimor,4)),round(stderror,5),formatC(round(as.numeric(tvalues),5)),round(pvalues,5))
    colnames(estshap2)<-c("Estimate","Std.Error","t value","Pr(>|t|)")
    rownames(estshap2)<-c(colnames(estimor))
    print(estshap2)
    cat("\nAnalysis of Variance\n")
    print(anva,na.print = "")
    #n<-length(x); k<-length(estimor); if(n/k<40) {AC<-AC+2*k*(k+1)/(n-k-1); nombre="AICc"} else nombre="AIC" #Nuevo DM
    FRAMESEL<-data.frame(R2=round(r,4),R2_Adj=round(r_ajus,4),SSR=round(sqe,4),AIC=round(AC,4),MSC=round(MSC,4))
    #colnames(FRAMESEL)[4]<-nombre #Nuevo DM
    rownames(FRAMESEL)<-c("")
    cat("\nSelection criteria")
    print(t(FRAMESEL))
    para<-t(matrix(paste(round(estimor,5),"(","?",round(stderror,5),")")))
    colnames(para)<-names(estimor)
    param<-data.frame(modelim,para,d.f.=paste(anva.1[1],",",anva.1[2]),F=round(anva.4[1],3),P.value=round(anva.5[1],3),FRAMESEL)
    
  }
  return(list(parmor=param,frames=FRAMESEL, Std.Error=stderror))      # Agregue "Std.Error=stderror" DM
}