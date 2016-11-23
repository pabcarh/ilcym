#' coef_ovi
#'
#' no details
#' 
#' @param modelm integer value, non linear model (position)
#' @param estimor numeric vector, estimated coefficients of the nonlinear model
#' @param stdovi numeric vector, standar deviation estimated
#' @param g formula, equation
#' @param frm formula, non linear equation
#' @param ovifemal data frame, temperatures and oviposition
#' @param alg a character, algorithm
#' @param pesos numeric vector
#' @param weights numeric vector
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
coef_ovi<-function(modelm,estimor,stdovi,g,frm,ovifemal,alg,pesos,weights)
{
  x<-ovifemal[,2]
  y<-ovifemal[,3]
  ind <- as.list(estimor)
  for (i in names(ind))
  {
    temp <- estimor[[i]]
    storage.mode(temp) <- "double"
    assign(i, temp)
  }
  if(alg=="Newton")
  {
    if(pesos=="WLS") modelo<-nls(g,estimor,data=ovifemal,weights=weights)
    if(pesos=="LS") modelo<-nls(g,estimor,data=ovifemal)
    art=0.00000000000000000000000000001
    sqe<-sum((residuals(modelo))^2)+art
    qq<-summary(modelo)
    tabla<-qq$"parameters"
    stderror<-tabla[,2]  #Agrege DM
    aic<-AIC(modelo) #Quite ,k=length(estimor)  DM
    if(pesos=="WLS") r<-1-sum(residuals(modelo)^2)/sum(weights*(y-mean(y))^2)
    if(pesos=="LS") r<-1-sum(residuals(modelo)^2)/sum((y-mean(y))^2)
    r_ajus<- 1 - ((length(x) - 1) / (length(x) - length(estimor))) * (1-r)
    yl<-predict(modelo)
    if(pesos=="WLS") MSC<-log(sum(weights*(yl-mean(yl))^2)/sum(weights*(y-mean(y))^2))-2*length(estimor)/length(x)
    if(pesos=="LS") MSC<-log(sum((yl-mean(yl))^2)/sum((y-mean(y))^2))-2*length(estimor)/length(x)
    anva.1 <- c(length(coef(modelo))-1,length(x)-length(coef(modelo)),length(x)-1)
    anva.2 <- c((sum(y^2)-length(y)*mean(y)^2)-sqe,sqe,sum(y^2)-length(y)*mean(y)^2)
    anva.3 <- c(anva.2[1]/anva.1[1],anva.2[2]/anva.1[2],NA)
    anva.4 <- c(anva.3[1]/anva.3[2],NA,NA)
    anva.5 <- c(1-pf(anva.4[1],anva.1[1],anva.1[2]),NA,NA)
    anva   <- cbind(anva.1,round(anva.2,4),round(anva.3,4),round(anva.4,4),round(anva.5,4))
    rownames(anva) <- c("Model","Error","Total")
    colnames(anva) <- c("DF","SS","MS","Fc","p")
    cat("\nNONLINEAR REGRESSION MODEL\n")
    cat("\nMethod:", alg)
    cat("\nFormula: ", frm,"\n")
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
    param<-data.frame(frm,para,d.f.=paste(anva.1[1],",",anva.1[2]),F=round(anva.4[1],3),P.value=round(anva.5[1],3),FRAMESEL)
    
  }
  if(alg=="Marquardtr")
  {
    
    if(modelm==1)   {
      yl<-g(x,a,b,c)
      g<-y~ 1 - exp(-(a * x + b * x^2 + c * x^3))
    }
    if(modelm==2)   {
      yl<-g(x,a,b)
      g<-y~ pgamma(x,a,b)
    }
    if(modelm==3)   {
      yl<-g(x,a,b)
      g<-y~ 1/(1+exp(a+b*x))
    }
    if(modelm==4)   {
      yl<-g(x,a,b)
      g<-y~ 1-exp(-a*x^b)
    }
    if(modelm==5)   {
      yl<-g(x,a,b,n)
      g<-y~1-exp(-((x-a)/n)^b)
    }
    if(modelm==6)   {
      yl<-g(x,a,b)
      g<-y~ pweibull(x,a,b)
    }
    
    sqe <- sum((y-yl)^2)
    var <- (sum((y-yl)^2)/(length(x) - length(estimor)))
    stderror<-sqrt(var*stdovi)
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
    cat("\nNONLINEAR REGRESSION MODEL\n")
    cat("\nMethod:", alg)
    cat("\nFormula: ", frm,"\n")
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
    param<-data.frame(frm,para,d.f.=paste(anva.1[1],",",anva.1[2]),F=round(anva.4[1],3),P.value=round(anva.5[1],3),FRAMESEL)
    
  }
  df<-length(x)-length(ind)
  sdli<-sqrt(sqe/df)
  qto<-qt(0.025,df)
  salidas<-list(sdli=sdli,qto=qto,parovi=param,gg=g,frames=FRAMESEL,yl=yl, Std.Error=stderror)      # Agregue "Std.Error=stderror" DM)
  return(salidas)
}
