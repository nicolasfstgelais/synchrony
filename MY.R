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
