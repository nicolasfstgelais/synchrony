library(lubridate)
library(plyr)

data=dataYM
chem="TP"
site="MonitoringLocationIdentifier"

corTP=synchro_calc(dataYM,"TP","MonitoringLocationIdentifier")
corTN=synchro_calc(dataYM,"TN","MonitoringLocationIdentifier")


synchro_calc<-function(data,chem,site,lim)
{
  data=data[!is.na(data[,chem]),]
  sites=unique(data[,site])
  corTemp=data.frame(s1=NA,s2=NA,coeff=NA,pval=NA)
  c=1
  i=1
  j=2
  for(i in 1:length(sites))
  {
    for(j in (i+1):length(sites)){
      R1=data[data[,site]%in%sites[i],]
      R2=data[data[,site]%in%sites[j],]
      int=intersect(R1$ym,R2$ym)
      if(length(int)<lim){next}
      cor=cor.test(R1[R1$ym%in%int,chem],R2[R2$ym%in%int,chem],method ="pear")
      corTemp[c,"s1"]=sites[i]
      corTemp[c,"s2"]=sites[j]
      corTemp[c,"pval"]=cor$p.value
      corTemp[c,"coeff"]=cor$estimate
      print(c)
      c=c+1
    }
  }
  return(corTemp)
}
  
  cTP=corTemp
  
  plot(cTP$coeff[cTP$pval<0.05],pch=16)
  points(cTN$coeff[cTN$pval<0.05],col="red",pch=16)
  
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
  

data=read.csv("data/resultWide.csv")

data$year=year(data$ActivityStartDate)
data$month=month(data$ActivityStartDate)
data$ym=data$year+(data$month/12)

data$Phosphorus=LtoN(data$Phosphorus)
data$Nitrogen=LtoN(data$Nitrogen)

dataYM=ddply(data, .(MonitoringLocationIdentifier,ym,year),summarize,TN=mean((Nitrogen),na.rm=TRUE),TP=mean((Phosphorus),na.rm=TRUE))

dataYM[is.nan(dataYM$TN),"TN"]=NA
dataYM[is.nan(dataYM$TP),"TP"]=NA


write.csv(dataYM,"dataYM.csv")

synchro=synchro_calc(dataYM[,c("MonitoringLocationIdentifier","ym","TN")])
