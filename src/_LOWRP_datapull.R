## Libraries
#devtools::install_github("SwampThingPaul/AnalystHelper")
library(AnalystHelper);
library(openxlsx)
library(plyr)
library(reshape2)
library(dssrip)
library(zoo)
library(classInt)
#
library(magrittr)
library(flextable)
library(ggplot2)


#Paths
wd="C:/Julian_LaCie/_GitHub/LOWRP_ModelEval"

paths=paste0(wd,c("/Plots/","/Export/","/Data/","/GIS","/src/","/_documents/"))
# Folder.Maker(paths);#One and done. Creates folders in working directory.
plot.path=paths[1]
export.path=paths[2]
data.path=paths[3]


LOWRP.dat=data.frame()

## LOK Data
vars=c("STAGE","LOK_NET_INFLOWS","LAKESTO")

dss_out=opendss(paste0(data.path,"RSMBN/ASR/RSMBN_output.dss"))  
for(i in 1:length(vars)){
  paths=paste0("/RSMBN/LOK/",vars[i],"/01JAN1965 - 01JAN2005/1DAY/SIMULATED/")  
  tmp=data.frame(getFullTSC(dss_out,paths))
  tmp$Date=date.fun(date.fun(rownames(tmp))-lubridate::ddays(1))
  rownames(tmp)<-NULL
  tmp=tmp[,c(2,1)]
  colnames(tmp)=c("Date","Data.Value")
  tmp$parameter=vars[i]
  tmp$SITE="LOK"
  tmp$Alt="ASR"
  LOWRP.dat=rbind(tmp,LOWRP.dat)
  print(i)
}

RSM.sites=c("S77","S77BF","S78","S79","BASIN2ECALOOS","BASIN2WCALOOS",
            "ECALOOS2BASIN","WCALOOS2BASIN","S80","S80_MAXFC","BASIN2C44",
            "C44TOBASIN","S308","S308BF","S308_QFC",
            "TMC2EST","S48","S49","NSF2EST",
            "S80_QPFCSOURCE_LAKE","S79_QPFCSOURCE_LAKE")

dss_out=opendss(paste0(data.path,"RSMBN/ASR/RSMBN_output.dss"))  
for(i in 1:length(RSM.sites)){
  paths=paste0("/RSMBN/",RSM.sites[i],"/FLOW/01JAN1965 - 01JAN2016/1DAY/SIMULATED/")  
  tmp=data.frame(getFullTSC(dss_out,paths))
  tmp$Date=date.fun(date.fun(rownames(tmp))-lubridate::ddays(1))
  rownames(tmp)<-NULL
  tmp=tmp[,c(2,1)]
  colnames(tmp)=c("Date","Data.Value")
  tmp$parameter="FLOW"
  tmp$SITE=RSM.sites[i]
  tmp$Alt="ASR"
  LOWRP.dat=rbind(tmp,LOWRP.dat)
  print(i)
}

sum(!(c(RSM.sites,"LOK")%in%unique(LOWRP.dat$SITE)))
unique(LOWRP.dat$Alt)


# write.csv(LOWRP.dat,paste0(export.path,"LOWRP_AltASR_data.csv"),row.names = F)