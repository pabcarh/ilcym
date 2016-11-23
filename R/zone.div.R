#' zone.div
#'
#' no details
#' 
#' @param dir1 character, path of the minimun temperature in FLT format
#' @param dir2 character, path of the maximun temperature in FLT format
#' @param ilon numeric vector, range of the longitude
#' @param ilat numeric vector, range of the latitude
#' @param R double, number of regions to divide
#' @param dir.out character, final path of the folder for the indexes (ASCII files of ERI,GI,AI)
#' @param name.out character, name of the ASCII files
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
zone.div<-function(dir1,dir2,ilon,ilat,R,dir.out,name.out=name.out,modelim=modelim,modelm=modelm,estadios=estadios,xi=xi,steps=steps,filtro=NULL,hfeno){
  
  lats=p.area(ilat,R)
  
  #############################################
  # Corriendo por cada area y almacenar su data
  #############################################
  
  for(j in R:1){
    
    TempsAll=GenMatriz(dir1,dir2,ilon,lats[j,])
    
    coords=TempsAll$coords;nfil=nrow(coords)
    
    x1=TempsAll$x1
    
    y1=TempsAll$y1
    
    
    RFinal=matrix(NA,nfil,3)
    
    ################################
    
    # Corriendo por punto del Area j
    
    ################################
    
    #system.time(
    
    for(i in 1:nfil){
      
      RFinal[i,]=GenActIndex(i,TempsAll=TempsAll,coords=coords,x1=x1,y1=y1,modelim=modelim,modelm=modelm,estadios=estadios,xi=xi,steps=steps,filtro,hfeno=hfeno)$indices
    }
    
    #)
    
    rm(TempsAll);rm(x1);rm(y1)
    
    Inds=data.frame(coords[1:nfil,],RFinal)
    
    rm(RFinal)
    
    ##################################################################
    
    # Generando los archivos con encabezado: Lon - Lat - GI - AI - ERI
    
    ##################################################################
    
    write.table(Inds,paste(dir.out,"file",j,".txt",sep=""),row.names = F)
    
    rm(Inds)
  }
  
  ###########################################
  
  # Corriendo por cada area y uniendo su data
  
  ###########################################
  
  Inds=read.table(paste(dir.out,"file",R,".txt",sep=""),header=T)
  
  if(R!=1)
  {
    for(j in (R-1):1){
      
      TempInds=read.table(paste(dir.out,"file",j,".txt",sep=""),header=T)
      Inds=rbind(Inds,TempInds)
    }
  }
  
  #mins<-apply(Inds[,-c(1,2)], 2,min, na.rm=T)
  #maxs<-apply(Inds[,-c(1,2)], 2,max, na.rm=T)
  
  rm(TempInds)
  gridded(Inds) = ~x+y ## Creando el objeto Grid
  
  writeAsciiGrid(Inds["X1"], na.value = -9999,paste(dir.out,"AI.asc",sep=""))
  
  writeAsciiGrid(Inds["X2"], na.value = -9999,paste(dir.out,"GI.asc",sep=""))
  
  writeAsciiGrid(Inds["X3"], na.value = -9999,paste(dir.out,"ERI.asc",sep=""))
  
  #writeAsciiGrid(Inds["X4"], na.value = -9999,paste(dir.out,"Ro.asc",sep=""))
  
  ### creacion d FLT(binario) ###
  #writeGDAL(Inds["X3"], paste(dir.out,name.out,"_ERI.flt",sep=""), drivername = "EHdr")
  #writeGDAL(Inds["X1"], paste(dir.out,name.out,"_AI.flt",sep=""), drivername = "EHdr")
  #writeGDAL(Inds["X2"], paste(dir.out,name.out,"_GI.flt",sep=""), drivername = "EHdr")
  
  #return(list(mins=mins,maxs=maxs))
}
