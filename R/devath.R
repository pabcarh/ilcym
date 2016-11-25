#' devath
#'
#' Function to obtain a table of frequencies, statistics of the data, and model parameters of presence of change of stage with respect to the temperature and the logarithm of the days, without considering the intercept, and including the adjustment indicators and the intervals confidence.
#'
#' @param opc integer value (1: life table data, 2: Cohort data)
#' @param datos data frame
#' @param estadios vector of length variable
#' @param est a character, name of the state evaluated
#' @param tp numeric vector of the constant temperature
#' @param intervalo double, scale in days
#' @param modelo a character, dicotomic model ('logit','probit', 'cloglog')
#' @param poli double, multiplier value of fecundity (only for parasitoids)
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @family example
#' @example inst/examples/example_devath.R
#' @export
devath <- function(opc, datos, estadios, est, tp,intervalo,modelo,poli=1)
{
  if (opc == 1) {
    tip<-1:length(estadios);exis.val<-c(1)
    for (is in 1:(length(estadios))){  if (estadios[is] == est){ data<-datos[[is]]; data[,4]=data[,4]*poli ; data[,3]=data[,3]*poli }}
    nag <- length(levels(factor(data[, 1])))
    not <- levels(factor(data[, 1]))
    tran <- puntos <- poracum <- por <- tol <- agente <- sampli<-list(1)
    for (i in 1:nag) {
      agente[i] <- list(subset(data, data[, 1] == not[i]))
      tol[i] <- list(agente[[i]][, 4])
      por[i] <- list(rep(0, nrow(agente[[i]])))
      for (j in 1:length(tol[[i]])) por[[i]][j] <- (tol[[i]][j])
      poracum[i] <- list(rep(0, length(tol[[i]])))
      for (j in 1:length(tol[[i]])) {
        if (j == 1) poracum[[i]][j] <- por[[i]][j]
        if (j != 1) poracum[[i]][j] <- por[[i]][j] + poracum[[i]][j - 1]
      }
      tran[i] <- list(round(poracum[[i]], 4))
      sampli[[i]]<-tran[[i]][length(tran[[i]])]
      if( nrow(agente[[i]])==1 & agente[[i]][1,3]==agente[[i]][1,4]){exis.val[i]=0}else{exis.val[i]=1}
    }
    cap=rep(TRUE,length(tp))
    matriz <- cbind(agente[[1]][, 1], agente[[1]][, 2],
                    log(agente[[1]][, 2]), agente[[1]][, 3],tran[[1]],sampli[[1]])
    for (i in 2:nag) matriz <- rbind(matriz, cbind(agente[[i]][, 1],
                                                   agente[[i]][, 2], log(agente[[i]][, 2]),
                                                   agente[[i]][, 3],tran[[i]],sampli[[i]]))
    matri <- as.data.frame(matriz)
    medobs<-sd<-se<-c(0)
    rept<-list(1)
    for (i in 1:nag)
    {
      rept[[i]]<-rep(agente[[i]][,2],agente[[i]][,4])
      medobs[i]<-stats::median(rept[[i]])
      sd[i]<- sd(rept[[i]])
      se[i]<-sd[i]/sqrt(length(rept[[i]]))
    }
    matriztotal<-as.matrix(cbind(agente[[1]][,1:2],round(log(agente[[1]][,2]),3),agente[[1]][,3:4],sampli[[1]],agente[[1]][,3]-sampli[[1]],round(agente[[1]][,4]/sampli[[1]],3),round(tran[[1]]*100/sampli[[1]],3)))
    for (i in 2:nag) matriztotal<-rbind(matriztotal,as.matrix(cbind(agente[[i]][,1:2],round(log(agente[[i]][,2]),3),agente[[i]][,3:4],sampli[[i]],agente[[i]][,3]-sampli[[i]],round(agente[[i]][,4]/sampli[[i]],3),round(tran[[i]]*100/sampli[[i]],3))))
    ntabla<-unique(matriztotal[,7])
  }
  if(opc==2){
    datas<-dia<-cont<-puntos<-poracum<-por<-fre<-tran <- n <- temp <- list(1)
    muestra<-posi<-exis.val<-c(1);posi[1]=NA
    tip<-1:length(estadios)
    num=tip[estadios==est]
    for (i in 1:length(datos)) {
      
      datas[i]<-list(datos[[i]][,num])
      vect1=datas[[i]]; if(length(vect1)==length(vect1[is.na(vect1)])){exis.val[i]<-FALSE}else{exis.val[i]<-TRUE}
      fre[i] <- list(seq(0,max(c(datas[[i]],0),na.rm=T)+intervalo,intervalo))
      
      for (j in 1:length(fre[[i]])) fre[[i]][j] <- length(subset(datas[[i]],datas[[i]] == intervalo*(j-1)))
      por[i] <- list(seq(0,max(c(datas[[i]],0),na.rm=T)+intervalo,intervalo))
      
      sfre=sum(fre[[i]]);if(sfre!=0){por[[i]] <- fre[[i]]/sfre}else{por[[i]] <- fre[[i]]}
      
      poracum[[i]] <- cumsum(por[[i]][-length(por[[i]])])
      
      puntos[i] <- list(as.character(poracum[[i]]))
      tran[i] <- list(round(as.numeric(puntos[[i]]), 2))

      cont[i] <- list(seq(0,max(c(datas[[i]],0),na.rm=T)+intervalo,intervalo))
      
      
      for (j in 1:length(poracum[[i]])) {
        cont[[i]][j] <- intervalo*(j-1)
        dia[[i]] <- subset(cont[[i]], cont[[i]] != 0)
      }
      n[i] <- list(rep(length(datas[[i]]), length(tran[[i]])))
      temp[i] <- list(rep(tp[i], length(tran[[i]])))
    }
    cap=rep(TRUE,length(tp))  
    matriz <- cbind(NA, NA, NA, NA, NA,NA)  
    for (i in 1:length(datos))
    {
      if(length(dia[[i]][is.na(dia[[i]])])==length(dia[[i]]))
      {cap[i]=FALSE
      next
      }
      matriz <- rbind(matriz, cbind(temp[[i]],dia[[i]], log(dia[[i]]), n[[i]], round(c(tran[[i]][-1],1.0)*n[[i]]),n[[i]]))
    }
    matriz1 <- as.data.frame(matriz[-1,])
    if(dim(matriz1)[1]==0){stop("Hay valores vacios en los conjuntos de datos de temperaturas")}
    
    nag <- length(levels(factor(matriz1[, 1])))
    not <- levels(factor(matriz1[, 1]))
    matriz2 <- agente <- list(1)
    for (i in 1:nag) {
      if(sum(fre[[i]])!=0)
      {
        ind3=1:length(fre[[i]][-1]);ind3=ind3[fre[[i]][-1]!=0]
        agente[i] <- list(subset(matriz1, matriz1[, 1] == not[i]))
        
        matriz2[i] <- list(agente[[i]][ind3, ])
        
        matriz2[[i]]<-rbind(agente[[i]][min(ind3)-1,],matriz2[[i]])  
      }else{ind3=1; agente[i] <- list(subset(matriz1, matriz1[, 1] == not[i])); matriz2[[i]]<-agente[[i]]}
    }
    matri <- matriz2[[1]]
    for (i in 2:nag) matri <- rbind(matri, matriz2[[i]])
    medobs<-sd<-se<-c(0)
    rept<-estdis<-list(1)
    for (i in 1:nag){
      estdis[[i]]<-subset(matri, matri[, 1] == not[i])
      for(j in nrow(estdis[[i]]):1){ if(j==1){if(estdis[[i]][,2][j]!=1){estdis[[i]][,6][j]<-0}else{estdis[[i]][,6][j]<-estdis[[i]][,5][j]} }else{ estdis[[i]][,6][j]<-estdis[[i]][,5][j]-estdis[[i]][,5][j-1]}}
      rept[[i]]<-rep(estdis[[i]][,2],estdis[[i]][,6])
      if(length(rept[[i]])==0){medobs[i]<-0;sd[i]=NA;se[i]=NA}else{
        medobs[i]<-stats::median(rept[[i]]);sd[i]<- sd(rept[[i]])
        se[i]<-sd[i]/sqrt(length(rept[[i]]))}
    }
    estdism<-estdis[[1]]
    for(i in 2:nag) estdism<-rbind(estdism,cbind(estdis[[i]]))
    matriztotal<-data.frame(estdism[,1],estdism[,2],estdism[,3],estdism[,4],estdism[,6],estdism[,6]/estdism[,4],estdism[,5])
  }
  Temperature <- as.factor(matri[,1])
  muestra<-matri[,6]
  Slope <- matri[,3]
  Develomep <- matri[,5]
  aicp <- suppressWarnings(warning(stats::AIC(modelop<-stats::glm(cbind(Develomep,muestra-Develomep) ~ Temperature - 1 + Slope, family = stats::binomial('probit')))))
  aicl <- suppressWarnings(warning(stats::AIC(modelol<-stats::glm(cbind(Develomep,muestra-Develomep) ~ Temperature - 1 + Slope, family = stats::binomial('logit')))))
  aicc <- suppressWarnings(warning(stats::AIC(modeloc<-stats::glm(cbind(Develomep,muestra-Develomep) ~ Temperature - 1 + Slope, family = stats::binomial('cloglog')))))
  rp<-1-sum(((matri[,5]/matri[,6])-round(stats::fitted.values(modelop),4))^2)/sum(((matri[,5]/matri[,6])-mean(matri[,5]/matri[,6]))^2)
  r_ajusp<- 1 - ((length(matri[,2]) - 1) / (length(matri[,2]) - length(tp)+1)) * (1-rp)
  rl<-1-sum(((matri[,5]/matri[,6])-round(stats::fitted.values(modelol),4))^2)/sum(((matri[,5]/matri[,6])-mean(matri[,5]/matri[,6]))^2)
  r_ajusl<- 1 - ((length(matri[,2]) - 1) / (length(matri[,2]) - length(tp)+1)) * (1-rl)
  rc<-1-sum(((matri[,5]/matri[,6])-round(stats::fitted.values(modeloc),4))^2)/sum(((matri[,5]/matri[,6])-mean(matri[,5]/matri[,6]))^2)
  r_ajusc<- 1 - ((length(matri[,2]) - 1) / (length(matri[,2]) - length(tp)+1)) * (1-rc)
  mscp<-log(sum((round(stats::fitted.values(modelop),4)-mean(round(stats::fitted.values(modelop),4)))^2)/sum(((matri[,5]/matri[,6])-mean((matri[,5]/matri[,6])))^2))-2*(length(tp)+1)/length(matri[,2])
  mscl<-log(sum((round(stats::fitted.values(modelol),4)-mean(round(stats::fitted.values(modelol),4)))^2)/sum(((matri[,5]/matri[,6])-mean((matri[,5]/matri[,6])))^2))-2*(length(tp)+1)/length(matri[,2])
  mscc<-log(sum((round(stats::fitted.values(modeloc),4)-mean(round(stats::fitted.values(modeloc),4)))^2)/sum(((matri[,5]/matri[,6])-mean((matri[,5]/matri[,6])))^2))-2*(length(tp)+1)/length(matri[,2])
  model<-c("probit","logit","cloglog")
  aic<-as.numeric(c(aicp,aicl,aicc))
  maic<-min(aic)
  mod=modelo
  aicv=aic[model==modelo]
  if (mod == "probit") {
    familia <- stats::family(stats::glm(cbind(Develomep,muestra-Develomep) ~ Temperature - 1 + Slope, family = stats::binomial('probit')))
    coeficientes <- stats::coef(summary(stats::glm(cbind(Develomep,muestra-Develomep) ~ Temperature - 1 + Slope, family = stats::binomial('probit'))))
    parametros <- stats::coef(stats::glm(cbind(Develomep,muestra-Develomep) ~ Temperature - 1 + Slope, family = stats::binomial('probit')))
    modelo<-stats::glm(cbind(Develomep,muestra-Develomep) ~ Temperature - 1 + Slope, family = stats::binomial('probit'))
    matri[,7]<-100*round(stats::fitted.values(modelop),4)
    medobsp<-sdp<-sep<-c(0)
    reptp<-estdisp<-list(1)
    for (i in 1:nag){
      estdisp[[i]]<-subset(matri, matri[, 1] == not[i])
      for(j in nrow(estdisp[[i]]):1) if(j==1)estdisp[[i]][,5][j]<-0 else estdisp[[i]][,5][j]<-estdisp[[i]][,7][j]-estdisp[[i]][,7][j-1]
      reptp[[i]]<-rep(estdisp[[i]][,2],estdisp[[i]][,5])
      
      if(length(reptp[[i]])==0){medobs[i]<-0;sdp[i]=NA;sep[i]=NA}else{
        medobsp[i]<-stats::median(reptp[[i]])
        sdp[i]<- sd(reptp[[i]])
        sep[i]<-sdp[i]/sqrt(length(reptp[[i]]))}
    }
    
  }
  if (mod=="logit") {
    familia <- stats::family(stats::glm(cbind(Develomep,muestra-Develomep) ~ Temperature - 1 + Slope, family = stats::binomial('logit')))
    coeficientes <- stats::coef(summary(stats::glm(cbind(Develomep,muestra-Develomep) ~ Temperature - 1 + Slope, family = stats::binomial('logit'))))
    parametros <- stats::coef(stats::glm(cbind(Develomep,muestra-Develomep) ~ Temperature - 1 + Slope, family = stats::binomial('logit')))
    modelo<-stats::glm(cbind(Develomep,muestra-Develomep) ~ Temperature -1 + Slope, family = stats::binomial('probit'))
    matri[,7]<-100*round(stats::fitted.values(modelol),4)
    medobsp<-sdp<-sep<-c(0)
    reptp<-estdisp<-list(1)
    for (i in 1:nag){
      estdisp[[i]]<-subset(matri, matri[, 1] == not[i])
      for(j in nrow(estdisp[[i]]):1) if(j==1)estdisp[[i]][,5][j]<-0 else estdisp[[i]][,5][j]<-estdisp[[i]][,7][j]-estdisp[[i]][,7][j-1]
      reptp[[i]]<-rep(estdisp[[i]][,2],estdisp[[i]][,5])
      
      if(length(reptp[[i]])==0){medobs[i]<-0;sdp[i]=NA;sep[i]=NA}else{
        medobsp[i]<-stats::median(reptp[[i]])
        sdp[i]<- sd(reptp[[i]])
        sep[i]<-sdp[i]/sqrt(length(reptp[[i]]))}
    }
    
  }
  if (mod=="cloglog") {
    familia <- stats::family(stats::glm(cbind(Develomep,muestra-Develomep) ~ Temperature - 1 + Slope, family = stats::binomial('cloglog')))
    coeficientes <- stats::coef(summary(stats::glm(cbind(Develomep,muestra-Develomep) ~ Temperature - 1 + Slope, family = stats::binomial('cloglog'))))  ###### el cambio es aqui
    parametros <- stats::coef(stats::glm(cbind(Develomep,muestra-Develomep) ~ Temperature - 1 + Slope, family = stats::binomial('cloglog')))
    modelo<-stats::glm(cbind(Develomep,muestra-Develomep) ~ Temperature - 1 + Slope, family = stats::binomial('probit'))
    matri[,7]<-100*round(stats::fitted.values(modeloc),4)
    medobsp<-sdp<-sep<-c(0)
    reptp<-estdisp<-list(1)
    for (i in 1:nag){
      estdisp[[i]]<-subset(matri, matri[, 1] == not[i])
      for(j in nrow(estdisp[[i]]):1) if(j==1)estdisp[[i]][,5][j]<-0 else estdisp[[i]][,5][j]<-estdisp[[i]][,7][j]-estdisp[[i]][,7][j-1]
      reptp[[i]]<-rep(estdisp[[i]][,2],estdisp[[i]][,5])
      
      if(length(reptp[[i]])==0){medobs[i]<-0;sdp[i]=NA;sep[i]=NA}else{
        medobsp[i]<-stats::median(reptp[[i]])
        sdp[i]<- sd(reptp[[i]])
        sep[i]<-sdp[i]/sqrt(length(reptp[[i]]))}
    }
  }
  parametros <- as.vector(parametros)
  slope <- parametros[nag + 1]
  intercepto <- parametros[1:nag]
  if(mod=="probit") cua1<-as.data.frame(t(as.matrix(c("probit",paste(round(slope,3),"(","?",round(as.numeric(coeficientes[,2][nag+1]),3),")"),round(as.numeric(coeficientes[,3][nag+1]),3),round(as.numeric(coeficientes[,4][nag+1]),3),
                                                      round(stats::deviance(modelop),3),round(aic[1],3),round(mscp,3),round(rp,3),round(r_ajusp,3)))))
  if(mod=="logit") cua1<-as.data.frame(t(as.matrix(c("logit",paste(round(slope,3),"(","?",round(as.numeric(coeficientes[,2][nag+1]),3),")"),round(as.numeric(coeficientes[,3][nag+1]),3),round(as.numeric(coeficientes[,4][nag+1]),3),
                                                     round(stats::deviance(modelol),3),round(aic[2],3),round(mscl,3),round(rl,3),round(r_ajusl,3)))))
  if(mod=="cloglog") cua1<-as.data.frame(t(as.matrix(c("cloglog",paste(round(slope,3),"(","?",round(as.numeric(coeficientes[,2][nag+1]),3),")"),round(as.numeric(coeficientes[,3][nag+1]),3),round(as.numeric(coeficientes[,4][nag+1]),3),
                                                       round(stats::deviance(modeloc),3),round(aic[3],3),round(mscc,3),round(rc,3),round(r_ajusc,3)))))
  if(mod=="logit" || mod=="probit") p50 <- (-intercepto)/slope else p50 <- (-(-log(-log(0.5)))-intercepto)/slope   ##  ERROR
  cat("\nESTIMATION OF PARAMETERS")
  print(familia)
  pval1=rep(0.0001,nrow(coeficientes))
  coeficientes[coeficientes[,4]<0.0001,4]=pval1[coeficientes[,4]<0.0001]
  coeficientes=cbind(round(coeficientes[,1:3],3),"Pr(>|z|)"=coeficientes[,4])
  print(coeficientes)
  
  crite<-data.frame(Deviance=c(stats::deviance(modelop),stats::deviance(modelol),stats::deviance(modeloc)), AIC=round(aic,3),MSC=round(c(mscp,mscl,mscc),3),R_Squared=round(c(rp,rl,rc),3),Adj_R_squared=round(c(r_ajusp,r_ajusl,r_ajusc),3))
  rownames(crite)<-model
  cat("\n")
  cat("\nSELECTION CRITERIA\n")
  print(crite)
  cat("\n")
  li<-ls<-c(1)
  for(i in 1:length(tp[cap])){
    li[i] <- p50[i] - stats::qt(0.975, length(matri[, 3]) - length(intercepto) -
                           1) * (1/slope)/sqrt(length(intercepto) + 1)
    ls[i] <- p50[i] + stats::qt(0.975, length(matri[, 3]) - length(intercepto) -
                           1) * (1/slope)/sqrt(length(intercepto) + 1)}
  tabmedi1<-data.frame(Temperature=tp[cap],Log_median=round(p50,3),Log_lower=round(li,3),Log_upper=round(ls,3),Days=round(exp(p50),3),Lower=round(exp(li),3),Upper=round(exp(ls),3),SD=round(sdp,3),SE=round(sep,3))
  tabmedi2<-data.frame(Temperature=tp[cap],Days_obs=medobs,SD_obs=round(sd,3),SE_obs=round(se,3))
  if(opc==2)colnames(matriztotal)<-c("T","Day's","Ln(Day's)","Sample","Dev_Sen","fr_Dev","Fr_Dev (%)")
  if(opc==1)colnames(matriztotal)<-c("T","Day's","Ln(Day's)","Sample","Dev_Sen","Sample_Dev","Sample_mort","fr_Dev","Fr_Dev (%)")
  cat("\nESTIMATED\n")
  print(tabmedi1)
  cat("\nOBSERVED\n")
  print(tabmedi2)
  cat("\n")
  cat("\nTABLE OF FREQUENCY\n")
  print(matriztotal)
  salidas <- list(Slope = slope, Temperaturas = as.numeric(not),
                  Interceptos = intercepto, matriz = matri, parametros = parametros,model=mod,tabest=tabmedi1,tabobs=tabmedi2,mattotal=matriztotal,cua=cua1,exis.val=exis.val, Std.Error=coeficientes[,2])
  return(salidas)
}
