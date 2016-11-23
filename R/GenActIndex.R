#' GenActIndex
#'
#' no details
#' 
#' @param posic integer value, position of the cell
#' @param TempsAll list, objetcs about geographic details of the temperatures
#' @param coords matrix, coordinates of the extent (all FLT files have the same extent)
#' @param x1 numeric vector, range of the longitude position
#' @param y1 numeric vector, range of the latitude position
#' @param modelim integer value, position of the nonlinear model (development rate for inmature states)
#' @param modelm integer value, position of the nonlinear model (development rate for adult states)
#' @param estadios character vector, all states
#' @param xi double, estimated median from female senescence
#' @param steps integer value, between 1 and 48
#' @param filtro numeric vector, range of the limit temperature
#' @param hfeno list, objects about phenology
#' 
#' @return None
#' @author Pablo Carhuapoma Ramos
#' @examples
#' #building
#'
#' @export
GenActIndex<-function(posic,TempsAll=TempsAll,coords=coords,x1=x1,y1=y1,modelim=modelim,modelm=modelm,estadios=estadios,xi=xi,steps=steps,filtro=NULL,hfeno)
{
  d1<-1:length(x1);d2<-1:length(y1);filtroin=TRUE
  plon=d1[x1==coords[posic,1]]
  plat=d2[y1==coords[posic,2]]
  
  Table=data.frame(mini=c(TempsAll$zTmin1[plon,plat],TempsAll$zTmin2[plon,plat],TempsAll$zTmin3[plon,plat],TempsAll$zTmin4[plon,plat],TempsAll$zTmin5[plon,plat],TempsAll$zTmin6[plon,plat],TempsAll$zTmin7[plon,plat],TempsAll$zTmin8[plon,plat],TempsAll$zTmin9[plon,plat],TempsAll$zTmin10[plon,plat],TempsAll$zTmin11[plon,plat],TempsAll$zTmin12[plon,plat]),
                   maxi=c(TempsAll$zTmax1[plon,plat],TempsAll$zTmax2[plon,plat],TempsAll$zTmax3[plon,plat],TempsAll$zTmax4[plon,plat],TempsAll$zTmax5[plon,plat],TempsAll$zTmax6[plon,plat],TempsAll$zTmax7[plon,plat],TempsAll$zTmax8[plon,plat],TempsAll$zTmax9[plon,plat],TempsAll$zTmax10[plon,plat],TempsAll$zTmax11[plon,plat],TempsAll$zTmax12[plon,plat]))
  
  Table=Table/10;Table2=cbind(id=1:nrow(Table),Table) ## filtrar en esta temperatura
  
  
  if(length(Table[is.na(Table)])==0)
  {
    if(!is.null(filtro)){tmm=apply(Table,2,mean,na.rm=TRUE);if(tmm[1] > filtro[1] && tmm[2] < filtro[2]){filtroin=TRUE}else{filtroin=FALSE}}
    if(filtroin)
    {
      inmaduros <-  estadios[-(length(estadios)-1):-(length(estadios))]
      maduros   <-  estadios[(length(estadios)-1):(length(estadios))]
      nmax=nrow(Table)
      matriz<-matrix(0,ncol=2*(length(inmaduros)*2+length(maduros)+1),nrow=nrow(Table))
      
      for(K in 1:length(inmaduros))
      {
        #  Desarrollo  # extraccion de funciones y parametros
        
        parametrosc <- hfeno$pdv_dv[[K]]
        funciont <- as.expression(hfeno$fdv_dv[[K]][[3]])
        
        #   Mortalidad # extraccion de funciones y parametros
        
        parametrosm <- hfeno$pmortal[[K]]
        funcionm <- as.expression(hfeno$mortal[[K]][[3]])
        
        RM=apply(Table2,1,RateI,Table2,K,parametrosc,parametrosm,funciont,funcionm,nmax,steps,modelim=modelim) ## procesamiento de tasa de desarrollo y mortalidad por cada temperatura
        matriz[,4*K-3]=RM[1,];matriz[,4*K-2]=RM[2,]
        matriz[,4*K-1]=RM[3,];matriz[,4*K]=RM[4,]
      }
      #  Hembras
      
      parametrosc <- hfeno$pfh_h
      for (i in names(parametrosc)){temp <- parametrosc[i];storage.mode(temp) <- "double";assign(i, temp)}
      formulac <- hfeno$fh_h
      funciont <- as.expression(formulac[[3]])
      RM=apply(Table2,1,RateI,Table2,(K+1),parametrosc,funciont=funciont,nmax=nmax,steps=steps,J=NA,modelim=modelim)
      matriz[,4*(K+1)-3]=RM[1,];matriz[,4*(K+1)-2]=RM[2,]
      
      #  Machos
      
      parametrosc <- hfeno$pfm_m;
      for (i in names(parametrosc)){temp <- parametrosc[i];storage.mode(temp) <- "double";assign(i, temp)}
      formulac <- hfeno$fm_m
      funciont <- as.expression(formulac[[3]])
      RM=apply(Table2,1,RateI,Table2,(K+2),parametrosc,funciont=funciont,nmax=nmax,steps=steps,J=NA,modelim=modelim)
      matriz[,4*(K+1)-1]=RM[1,];matriz[,4*(K+1)]=RM[2,]
      matriz[matriz>1]<-1;matriz[matriz<0]<-0
      
      #  Fecundidad
      
      parametrosc <- hfeno$ptazaeh_h
      for (i in names(parametrosc)){temp <- parametrosc[i];storage.mode(temp) <- "double";assign(i, temp)}
      formulac <- hfeno$ftazaeh_h
      funciont <- as.expression(formulac[[3]])
      RM=apply(Table2,1,RateI,Table2,(K+3),parametrosc,funciont=funciont,nmax=nmax,steps=steps,J=3,modelim=modelim)
      matriz[,4*(K+1)+1]=RM[1,];matriz[,4*(K+1)+2]=RM[2,]
      
      
      TF=data.frame(Ind=1:nmax)
      for(i in 1:K){TF1=data.frame(DaiSurv1=(1-matriz[,4*i-1])^matriz[,4*i-3],DaiSurv2=(1-matriz[,4*i])^matriz[,4*i-2]);TF=cbind(TF,TF1)}
      for(i in 1:K){Gl0=data.frame(Gl01=1/matriz[,4*i-3],Gl02=1/matriz[,4*i-2]);TF=cbind(TF,Gl0)}
      TF$Gl11=(1/matriz[,4*(K)+1])/(1/xi);TF$Gl22=(1/matriz[,4*(K)+2])/(1/xi);TF=TF[,-1]  ## aqui faltaba esta ultima linea pues no consideraba la tasa de desarrollo de la hembra en el ultimo mes
      temp1=TF[,c(2*K+seq(1,2*K,2),ncol(TF)-1)];temp2=TF[,c(2*K+seq(2,2*K,2),ncol(TF))]
      TF$Gl1=apply(temp1,1,sum);TF$Gl2=apply(temp2,1,sum)
      temp1<-temp2<-1;for(j in 1:K){temp1=temp1*(1-matriz[,4*j-1]);temp2=temp2*(1-matriz[,4*j])} ## aqui usaba como contador ala "i" en ves de la "j"
      TF$Is1=temp1;TF$Is2=temp2
      TF$Ro1=(TF$Is1*matriz[,4*K+5])/2;TF$Ro1[TF$Ro1<0]=0; TF$Ro2=(TF$Is2*matriz[,4*K+6])/2;TF$Ro2[TF$Ro2<0]=0
      TF$rm1=log(TF$Ro1)/TF$Gl1; TF$rm2=log(TF$Ro2)/TF$Gl2
      TF$Fr1=exp(TF$rm1); TF$Fr2=exp(TF$rm2)
      TF$Dt1=log(2)/TF$rm1; TF$Dt2=log(2)/TF$rm2
      #Ro=mean(c(TF[,"Ro1"],TF[,"Ro2"]))
      
      meses=c(31,28,31,30,31,30,31,31,30,31,30,31);r1=rep(1,nrow(TF))
      GI=sum((meses-1)/TF[,"Gl1"] + 1/TF[,"Gl2"])
      #temp.fr1=TF[,"Fr1"];ind1=temp.fr1>1;temp.fr1=temp.fr1[ind1];temp.meses=meses[ind1]  ## hay un objeto q no se utiliza
      #temp.fr2=TF[,"Fr2"];ind2=temp.fr2>1;temp.fr2=temp.fr2[ind2]
      #if(length(temp.fr2)==0){temp.fr2=1}
      #if(length(temp.fr1)==0){AI=NA}else{AI=log(prod(c(temp.fr1^(temp.meses-1),temp.fr2)),10)}
      temp.fr1=TF[,"Fr1"];ind1=temp.fr1<1;temp.fr1[ind1]=1
      temp.fr2=TF[,"Fr2"];ind2=temp.fr2<1;temp.fr2[ind2]=1
      AI=log(prod(c(temp.fr1^(meses-1),temp.fr2)),10) ## en este caso hubo el problema que no consideraba a los valores del ultimo mes si el resto valores eran menores a 1
      
      #ERI=1;for(s in 1:K){n0=(sum(meses[TF[,2*s-1]==0]-1) + sum(r1[TF[,2*s]==0]))/365;ERI=ERI*(1-n0)}
      #if(is.null(DL)){DL=0} ## La diferencia limite segun su origen siempre ha de ser positivo
      #ERI=1;for(s in 1:K){n0=(sum(meses[TF[,2*s-1]<=DL]-1) + sum(r1[TF[,2*s]<=DL]))/365;ERI=ERI*(1-n0)}
      
      #       ERI1=c(apply(TF[,seq(1,2*K,2)]*TF[,"Ro1"],1,mean,na.rm = TRUE), apply(TF[,seq(2,2*K,2)]*TF[,"Ro2"],1,mean,na.rm = TRUE))
      #       ERI=mean(ERI1)/max(ERI1);ERI[is.na(ERI)]=0;ERI[ERI>1]=1
      
      
      ERI=1;n0=(sum(meses[TF[,"Ro1"] <= 1]-1) + sum(r1[TF[,"Ro1"] <= 1]))/365;ERI=ERI*(1-n0)
      
      
      #indices=c(AI=AI,GI=GI,ERI=ERI,Ro=Ro)
      indices=c(AI=AI,GI=GI,ERI=ERI)
      return(list(indices=indices))
    }else{indices=c(AI=0,GI=0,ERI=0);return(list(indices=indices))}
  }else{
    
    #indices=c(AI=NA,GI=NA,ERI=NA,Ro=NA)
    indices=c(AI=NA,GI=NA,ERI=NA)
    return(list(indices=indices))
  }
  
}
