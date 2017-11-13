library(tidyr)
library(dplyr)

# long to wide

LtoN <- function(x) as.numeric(as.character(x))


data=read.csv("data/result.csv")
head(data)
colnames(data)
key="CharacteristicName"
values="ResultMeasureValue"
headers=c("ActivityStartDate","MonitoringLocationIdentifier")


data.sel=data[,c(headers,values,key)]
groupVar=list(data.sel$ActivityStartDate,data.sel$MonitoringLocationIdentifier,data.sel$CharacteristicName)

data.sel=aggregate(data.sel, by=groupVar,FUN=mean) 

# if more than one measure per site/day -> take the first one, would be better to g by depth
data.sel.dist=distinct(data.sel, ActivityStartDate,MonitoringLocationIdentifier,CharacteristicName,.keep_all = T)



dataWide<- spread(data.sel.dist,CharacteristicName, ResultMeasureValue)

dataWide$Phosphorus=LtoN(dataWide$Phosphorus)
dataWide$Nitrogen=LtoN(dataWide$Nitrogen)

write.csv(dataWide,"data/resultWide.csv")
