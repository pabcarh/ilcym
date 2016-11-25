datos<-list(1)
estadios1<-list.files(system.file("extdata",package = "ilcym"))[-5]
# order them acording sequence of life cycle
estadios1<-estadios1[c(1,3,5,2,4)]
# separate only the state names
estadios<-gsub(".tsv", "",estadios1)

# upload the cohort data about PTM (potato pest)
for(i in 1:length(estadios1)) # CAMBIO
{
  datos[[i]]<-read.table(system.file("extdata", estadios1[i], package = "ilcym"), header = FALSE)
}

# evaluate an state (egg for example)
est<-estadios[1]

estadiosINMD<-estadios[1:(length(estadios)-2)]
estadiosADLT<-estadios[(length(estadios)-1):length(estadios)]

intervalo<-1
temperatura<-temp(datos,estadios,est)
tp<-temperatura$temperature
opc=1 #cohort data
modelo<-"logit"

distri<-devath(opc,datos,estadios,est,tp,intervalo=intervalo,modelo)
matri<-distri$matriz
slope <- distri$Slope; names(slope)<-"slope"
parametros<-distri$parametros
test<-distri$model
exis.val <- distri$exis.val
corrx<-c(0,5)

plotF<-PlotDevath (matri,parametros,test,corrx,lgx=0.1,lgy=90,tam=0.8,
                  labx="Development time (Ln-days)",laby="Accumulated frequency (%)",
                  titulo=modelo,grises=FALSE,intervalo=intervalo,exis.val=exis.val)
