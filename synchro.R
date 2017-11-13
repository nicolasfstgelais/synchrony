library(lubridate)
library(plyr)
library(foreach)
library(parallel)
library(doSNOW)
library(geosphere)
library(synchrony)
library(cdlTools)
library(ggplot2)
library(fiftystater)

# Calculate the number of cores 
no_cores <- detectCores() - 2

# Initiate cluster
cl <- makeCluster(no_cores)
registerDoSNOW(cl)

#data=dataYM
#chem="TN"
#site="MonitoringLocationIdentifier"
#lim=20
#state=NA
#state=35

synchro_calc<-function(data,chem,site="MonitoringLocationIdentifier",state=NA,lim)
{
  cl <- makeCluster(no_cores)
  registerDoSNOW(cl)
  
  data=data[!is.na(data[,chem]),]
  if(!is.na(state)){data=data[data$StateCode==state,]}
  sites=unique(data[,site])
  x <-
  foreach(i=1:length(sites), .combine='rbind') %:%
  foreach(j= (i+1):length(sites), .combine='rbind') %dopar% {
    R1=data[data[,site]%in%sites[i],]
    R2=data[data[,site]%in%sites[j],]
    int=intersect(R1$ym,R2$ym)
    if(length(int)>lim){
      cor=synchrony::kendall.w (cbind(R1[R1$ym%in%int,chem],R2[R2$ym%in%int,chem]), nrands = 0 ,quiet=F)
      #cor=cor.test(R1[R1$ym%in%int,chem],R2[R2$ym%in%int,chem],method ="spear")
      dist=geosphere::distm (as.numeric(R1[1,c("LongitudeMeasure","LatitudeMeasure")]), as.numeric(R2[1,c("LongitudeMeasure","LatitudeMeasure")]), fun = geosphere::distHaversine)
      #corTemp[c,"s1"]=sites[i]
      #corTemp[c,"s2"]=sites[j]
      #corTemp[c,"pval"]=cor$p.value
      #corTemp[c,"coeff"]=cor$estimate
      #data.frame(s1=sites[i],s2=sites[j],pval=cor$p.value,coeff=cor$estimate,dist=dist)
      data.frame(s1=sites[i],s2=sites[j],pval=cor$pval,coeff=cor$w.corrected,dist=dist)
      
      #print(c)
      #c=c+1
    }
  }
  return(x)
  }


#write.csv(x,"data/corTN.csv")


#general results
corTN=read.csv("data/corTN.csv")
corTP=read.csv("data/corTP.csv")

corTP$dist=corTP$dist/1000
corTN$dist=corTN$dist/1000


plot(log(corTP$coeff[corTP$pval<1]+1)~(corTP$dist[corTP$pval<1]),xlim=c(1,10000000))
plot((corTN$coeff[corTN$pval<1])~log(corTN$dist[corTN$pval<1]),col="red")
points((corTP$coeff[corTP$pval<1])~log(corTP$dist[corTP$pval<1]))

plot((corTN$coeff[corTN$pval<1])~log(corTN$dist[corTN$pval<1]),col="red",type="n")
lines(lowess((corTN$coeff[corTN$coeff>0])~log(corTN$dist[corTN$coeff>0]), f = 0.7), col = "red",lwd=4)
lines(lowess((corTP$coeff[corTP$coeff>0])~log(corTP$dist[corTP$coeff>0]), f = 0.7), col = "chartreuse3",lwd=4)

lines(lowess((corTN$coeff[corTN$coeff<0])~log(corTN$dist[corTN$coeff<0]), f = 0.7), col = "red",lwd=4)
lines(lowess((corTP$coeff[corTP$coeff<0])~log(corTP$dist[corTP$coeff<0]), f = 0.7), col = "chartreuse3",lwd=4)



plot((corTN$coeff[corTN$pval<0.05])~corTN$dist[corTN$pval<0.05],col="red",log="x",xlim=c(1,10000000))
points((corTP$coeff[corTP$pval<0.05])~(corTP$dist[corTP$pval<0.05]),log="x")


plot(density(corTP$dist),col="chartreuse3",lwd=2)
lines(density(corTN$dist),col="red",lwd=2)


dataYM=read.csv("data/dataYM.csv")
stations=read.csv("data/station.csv")


dataYM=merge(dataYM,stations[,c("MonitoringLocationIdentifier","LatitudeMeasure","LongitudeMeasure","StateCode")],by ="MonitoringLocationIdentifier",all.x = T)

dataState=ddply(data, .(StateCode),summarize,lat=mean((LatitudeMeasure),na.rm=TRUE),long=mean((LongitudeMeasure),na.rm=TRUE))
rownames(dataState)=dataState[,1];dataState=dataState[,-1]

dataYM=dataYM[rowSums(!is.na(dataYM[,c("TN","TP")])) == 2,]

# by state
cl <- makeCluster(no_cores)
registerDoSNOW(cl)

corTNbyState=list()
corTPbyState=list()

for(i in as.character(unique(dataYM$StateCode))){

  sTN=synchro_calc(dataYM,chem="TN",state=i,lim = 20)
  sTP=synchro_calc(dataYM,chem="TP",state=i,lim = 20)
  
  if(is.null(sTN)){next}
  if(nrow(sTN)==1){next}
  
  corTNbyState[[i]]=sTN
  corTPbyState[[i]]=sTP
  

  #ppi=300
  #jpeg(paste0("figures/",i,".jpg"),width=7*ppi, height=6*ppi,bg="transparent",res=ppi)
  #plot( corTPbyState[[i]]$coeff~ log(corTPbyState[[i]]$dist),col="chartreuse3",pch=16)
  #points( corTNbyState[[i]]$coeff~ log(corTNbyState[[i]]$dist),col="red",pch=16)
  #lines(lowess(corTNbyState[[i]]$coeff~log(corTNbyState[[i]]$dist+1), f = 0.7), col = "red",lwd=4)
  #lines(lowess(corTPbyState[[i]]$coeff~log(corTPbyState[[i]]$dist+1), f = 0.7), col = "chartreuse3",lwd=4)
  #dev.off()
  print(i)
}
stopCluster(cl)

corByState=matrix(NA,nrow(dataState),4,dimnames=list(as.character(rownames(dataState)),c("TN","TP","lat","long")))
corByState=as.data.frame(corByState)


for(i in (names(corTNbyState)))
{
  corByState[i,"TN"]=mean(corTNbyState[[i]]$coeff)
  corByState[i,"TP"]=mean(corTPbyState[[i]]$coeff)
  corByState[i,c("lat","long")]=as.numeric(dataState[i,])
  }
corByState$stateName=tolower(apply(as.matrix(rownames(corByState)),1,fips,"Name"))


corByState=corByState[-which(is.na(corByState$TN)),]
# TN map


#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('blue','red'))

#This adds a column of color values
# based on the y values
col <- rbPal(10)[as.numeric(cut(corByState[,"TN"],breaks = 10))]
col <- rbPal(10)[as.numeric(cut(corByState[,"TP"],breaks = 10))]


plot(corByState$long,corByState$lat,col=col,pch=16,cex=1.5)
#text(corByState$long,corByState$lat,labels = round(corByState$TN,digits = 2))

plot(corByState$TP,corByState$TN)



#data("fifty_states") # this line is optional due to lazy data loading


# map_id creates the aesthetic mapping to the state name column in your data
p <- ggplot(corByState, aes(map_id = stateName)) + 
  # map points to the fifty_states shape data
  geom_map(aes(fill = TP), map = fifty_states) + 
  expand_limits(x = fifty_states$long, y = fifty_states$lat) +
  coord_map() +
  scale_x_continuous(breaks = NULL) + 
  scale_y_continuous(breaks = NULL) +
  labs(x = "", y = "") +
  theme(legend.position = "bottom", 
        panel.background = element_blank())

p
# add border boxes to AK/HI
p + fifty_states_inset_boxes() 


## by year

tempo=cut(dataYM$ym,breaks=c(-Inf,1970,1975,1980,1985,1990,1995,2000,2005,2010,Inf))

corTNbyPeriod=list()
corTPbyPeriod=list()

cl <- makeCluster(no_cores)
registerDoSNOW(cl)


for(i in sort(unique(tempo))){
  print(i)
  
  sTN=synchro_calc(dataYM[tempo%in%i,],chem="TN",lim = 50)
  sTP=synchro_calc(dataYM[tempo%in%i,],chem="TP",lim = 50)
  
  if(is.null(sTN)){next}
  if(nrow(sTN)==1){next}
  
  corTNbyPeriod[[i]]=sTN
  corTPbyPeriod[[i]]=sTP
  
  
  ppi=300
  jpeg(paste0("figures/",i,".jpg"),width=7*ppi, height=6*ppi,bg="transparent",res=ppi)
  plot( corTPbyPeriod[[i]]$coeff~ log(corTPbyPeriod[[i]]$dist+1),col="chartreuse3",pch=16)
  points( corTNbyPeriod[[i]]$coeff~ log(corTNbyPeriod[[i]]$dist+1),col="red",pch=16)
  lines(lowess(corTNbyPeriod[[i]]$coeff~log(corTNbyPeriod[[i]]$dist+1), f = 0.7), col = "red",lwd=4)
  lines(lowess(corTPbyPeriod[[i]]$coeff~log(corTPbyPeriod[[i]]$dist+1), f = 0.7), col = "chartreuse3",lwd=4)
  dev.off()
  print(i)
}

corByPeriod=matrix(NA,length(unique(tempo)),2,dimnames=list(unique(tempo),c("TN","TP")))

i="(2005,2010]"

c=1
col=heat.colors(length(unique(tempo)))
for(i in sort(unique(tempo),decreasing = T))
{
  corByPeriod[i,"TN"]=mean(corTNbyPeriod[[i]]$coeff)
  corByPeriod[i,"TP"]=mean(corTPbyPeriod[[i]]$coeff)
  if(is.na(corByPeriod[i,"TN"]))next
  if(c==1)plot(density(corTPbyPeriod[[i]]$coeff),col=col[c],lwd=2,ylim=c(0,5))
    if(c!=1)lines(density(corTPbyPeriod[[i]]$coeff),col=col[c],lwd=2)
  c=c+1
}

legend("topleft", legend=sort(unique(tempo),decreasing = T), col=col,pch=16)



plot(corByPeriod[,"TP"])


for(i in 1:length(corTNbyState)){
  print(median(corTNbyState[[i]]$coeff))
  print(median(corTPbyState[[i]]$coeff))
}
lapply(corTNbyState, f)

corTPbyState[[15]]


site=29

te[[site]]$coeff)
median(corTNbyState[[site]]$coeff)


lapply(corTNbyState[[]],median)
corTN=
corTP=synchro_calc(dataYM,chem="TP",state=42,lim = 10)

a=corTNbyState[[51]]


median(corTN$coeff)
median(corTP$coeff)

plot((corTN$coeff[corTN$pval<1])~log(corTN$dist[corTN$pval<1]),col="red",pch=16)
points((corTP$coeff[corTP$pval<1])~log(corTP$dist[corTP$pval<1]),col="chartreuse3",pch=16)

lines(lowess((corTN$coeff[])~log(corTN$dist[]+1), f = 0.7), col = "red",lwd=4)
lines(lowess((corTP$coeff[])~log(corTP$dist[]+1), f = 0.7), col = "chartreuse3",lwd=4)





#select rows for which we have both TP and TN



  p60=dataYM$ym<1970
  p70=dataYM$ym<1980&dataYM$ym>=1970
  p80=dataYM$ym<1990&dataYM$ym>=1980
  p90=dataYM$ym<2000&dataYM$ym>=1990
  p20=dataYM$ym>2000

  
  corTP=synchro_calc(dataYM,"TP","MonitoringLocationIdentifier")
  
  corTN=synchro_calc(dataYM,"TN","MonitoringLocationIdentifier")
  
  corTN60=synchro_calc(dataYM[p60,],"TN","MonitoringLocationIdentifier",20)
  corTN70=synchro_calc(dataYM[p70,],"TN","MonitoringLocationIdentifier",20)
  corTN80=synchro_calc(dataYM[p80,],"TN","MonitoringLocationIdentifier",20)
  corTN90=synchro_calc(dataYM[p90,],"TN","MonitoringLocationIdentifier",20)
  corTN20=synchro_calc(dataYM[p20,],"TN","MonitoringLocationIdentifier",20)
  
  plot( corTN80$coeff[corTN80$pval<0.05],pch=16)
  points(corTN90$coeff[corTN90$pval<0.05],col="red",pch=16)
  points(corTN20$coeff[corTN20$pval<0.05],col="green",pch=16)
  
  
  
  corTP60=synchro_calc(dataYM[p60,],"TP","MonitoringLocationIdentifier",20)
  corTP70=synchro_calc(dataYM[p70,],"TP","MonitoringLocationIdentifier",20)
  corTP80=synchro_calc(dataYM[p80,],"TP","MonitoringLocationIdentifier",20)
  corTP90=synchro_calc(dataYM[p90,],"TP","MonitoringLocationIdentifier",20)
  corTP20=synchro_calc(dataYM[p20,],"TP","MonitoringLocationIdentifier",20)
  
  
    
    synchro80=synchro_calc(dataYM[p80,c("BQMA","ym","TN")])
  synchro90=synchro_calc(dataYM[p90,c("BQMA","ym","TN")])
  synchro95=synchro_calc(dataYM[p95,c("BQMA","ym","TN")])
  synchro20=synchro_calc(dataYM[p20,c("BQMA","ym","TN")])
  synchro85=synchro_calc(dataYM[p85,c("BQMA","ym","TN")])
  
  synchro80$period=1980
  synchro85$period=1985
  synchro90$period=1990
  synchro95$period=1995
  synchro20$period=2000
  



synchro=synchro_calc(dataYM[,c("MonitoringLocationIdentifier","ym","TN")])
