## 
# Title:      LOWRP Model Evaluation
# Objective:  Evalute discharges to the estuaries
#             To support Draft PIR/EIS review
# Created by: Paul Julian; pjulian@sccf.org
# Created on: 02/22/2022
## 
## 

## BAD ## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

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

## Functions
consec.startend=function(var){
  runs=rle(var)
  myruns = which(runs$values == TRUE)
  runs.lengths.cumsum = cumsum(runs$lengths)
  ends = runs.lengths.cumsum[myruns]
  newindex = ifelse(myruns>1, myruns-1, 0)
  starts = runs.lengths.cumsum[newindex] + 1
  if (0 %in% newindex) starts = c(1,starts)
  rslt=list(starts=starts,ends=ends)
  return(rslt)
}

# -------------------------------------------------------------------------
list.files(paste0(data.path,"RSMBN"))
alts=c("ECB","FWO","ALT1BWR","ASR")
n.alts=length(alts)

cols.alts=c(grey.colors(2),wesanderson::wes_palette("Zissou1",2,"continuous"))

# Discharge ---------------------------------------------------------------
RSM.sites=c("S79","S80","S80_QPFCSOURCE_LAKE","S79_QPFCSOURCE_LAKE",
            "S77","S308","S77_QFC","S308_QFC","S308BF","S77BF",
            "TMC2EST","S48","S49","NSF2EST","S2","S3","S4BP",
            "S351_QFC","S351_FC_SHIFT2_ENVTARG","S354_QFC","S354_FC_SHIFT2_ENVTARG",
            "S77_QFC","C10A_QFC")
# S271 = C10A
q.dat=data.frame()
for(j in 1:n.alts){
  dss_out=opendss(paste0(data.path,"RSMBN/",alts[j],"/RSMBN_output.dss"))  
  
  for(i in 1:length(RSM.sites)){
    paths=paste0("/RSMBN/",RSM.sites[i],"/FLOW/01JAN1965 - 01JAN2005/1DAY/SIMULATED/")  
    tmp=data.frame(getFullTSC(dss_out,paths))
    tmp$Date=date.fun(date.fun(rownames(tmp))-lubridate::ddays(1))
    rownames(tmp)<-NULL
    tmp$SITE=RSM.sites[i]
    tmp$Alt=alts[j]
    q.dat=rbind(tmp,q.dat)
    print(i)
  }
}

for(j in 2:n.alts){
  dss_out=opendss(paste0(data.path,"RSMBN/",alts[j],"/RSMBN_output.dss"))  
  
    paths=paste0("/RSMBN/STA2TMC/FLOW/01JAN1965 - 01JAN2005/1DAY/SIMULATED/")  
    tmp=data.frame(getFullTSC(dss_out,paths))
    tmp$Date=date.fun(date.fun(rownames(tmp))-lubridate::ddays(1))
    rownames(tmp)<-NULL
    tmp$SITE="STA2TMC"
    tmp$Alt=alts[j]
    q.dat=rbind(tmp,q.dat)
    print(j)
  }
unique(q.dat$Alt)
unique(q.dat$SITE)


q.dat=q.dat[order(q.dat$Alt,q.dat$SITE,q.dat$Date),]
q.dat$CY=as.numeric(format(q.dat$Date,"%Y"))
q.dat$month=as.numeric(format(q.dat$Date,"%m"))

q.dat.xtab=reshape2::dcast(q.dat,Alt+Date+CY+month~SITE,value.var="FLOW",function(x)mean(x,na.rm=T))
head(q.dat.xtab)
q.dat.xtab$S79.14d=with(q.dat.xtab,ave(S79,Alt,FUN=function(x) c(rep(NA,13),rollapply(x,width=14,FUN=function(x)mean(x,na.rm=T)))))
q.dat.xtab$S80.14d=with(q.dat.xtab,ave(S80,Alt,FUN=function(x) c(rep(NA,13),rollapply(x,width=14,FUN=function(x)mean(x,na.rm=T)))))
q.dat.xtab$S79.30d=with(q.dat.xtab,ave(S79,Alt,FUN=function(x) c(rep(NA,29),rollapply(x,width=30,FUN=function(x)mean(x,na.rm=T)))))

## SLE
sle.gw=read.table(paste0(data.path,"gw_flws.txt"),header=F,col.names = c("Date","sle_gw"))
sle.gw$Date=with(sle.gw,date.fun(Date,form="%m/%d/%Y"))
q.dat.xtab=merge(q.dat.xtab,sle.gw,"Date",all.x=T)
q.dat.xtab=q.dat.xtab[order(q.dat.xtab$Alt,q.dat.xtab$Date),]

# q.dat.xtab$SLE.S80trib=rowSums(q.dat.xtab[,c("S80","TMC2EST","S48","S49","NSF2EST","sle_gw")],na.rm=T)

q.dat.CY=ddply(q.dat.xtab,c("CY","Alt"),summarise,CRE=sum(S79),SLE=sum(SLE.S80trib))#SLE.S80trib.cfs=sum(SLE.S80trib,na.rm=T))
q.dat.CY.mean=ddply(q.dat.CY,'Alt',summarise,mean.CRE=mean(CRE),mean.SLE=mean(SLE))

q.dat.xtab.melt=reshape2::melt(q.dat.xtab,id.vars=c("Date","Alt","CY","month"))
q.dat.CY=ddply(subset(q.dat.xtab.melt,variable%in%c("S79_QPFCSOURCE_LAKE","SLE.S80trib")),
                      c("CY","Alt","variable"),summarise,TFlow=sum(cfs.to.acftd(value)))#SLE.S80trib.cfs=sum(SLE.S80trib,na.rm=T))

q.dat.CY.mean=reshape2::dcast(q.dat.CY,variable~Alt,value.var = "TFlow",mean)
q.dat.CY.mean=q.dat.CY.mean[c(2,1),c("variable",alts)]

apply(q.dat.CY.mean[,2:5],2,sum)

# q.dat.CY=ddply(q.dat,c("CY","Alt",'SITE'),summarise,TFlow.cfs=sum(FLOW,na.rm=T))
# reshape2::dcast(q.dat.CY,Alt~SITE,value.var = "TFlow.cfs",sum,na.rm=T)

# tmp=subset(q.dat.xtab,Alt=="ASR")
# vars=c("Date", "Alt","NSF2EST", "S308", "S308_QFC","S48", "S49", 
#        "S80", "S80_QPFCSOURCE_LAKE", "TMC2EST", "sle_gw")
# write.csv(tmp[,vars],paste0(export.path,"tmp_SLE_asr.csv"),row.names=F)
# 
# tmp=subset(q.dat.xtab,Alt=="FWO")
# vars=c("Date", "Alt","NSF2EST", "S308", "S308_QFC","S48", "S49",
#        "S80", "S80_QPFCSOURCE_LAKE", "TMC2EST", "sle_gw")
# write.csv(tmp[,vars],paste0(export.path,"tmp_SLE_FWOr.csv"),row.names=F)

# tmp=subset(q.dat.xtab,Alt=="ECB")
# vars=c("Date", "Alt","NSF2EST", "S308", "S308_QFC","S48", "S49",
#        "S80", "S80_QPFCSOURCE_LAKE", "TMC2EST", "sle_gw")
# write.csv(tmp[,vars],paste0(export.path,"tmp_SLE_ECB.csv"),row.names=F)

# flood control -----------------------------------------------------------
RSM.sites.region=data.frame(SITE=c("S351_QFC","S351_FC_SHIFT2_ENVTARG",
                                   "S354_QFC","S354_FC_SHIFT2_ENVTARG",
                                   "S77_QFC","S308_QFC","C10A_QFC"),
                            Region=c(rep("WCAs",4),"Cal",'StL',"LWLagoon"))

regq.dat=merge(subset(q.dat,SITE%in%RSM.sites.region$SITE),RSM.sites.region,"SITE")

regq.dat.CY=ddply(regq.dat,c("CY","Alt",'Region'),summarise,TFlow.kAcft=sum(FLOW,na.rm=T))
reshape2::dcast(regq.dat.CY,Alt~Region,value.var = "TFlow.kAcft",mean)


regq.dat.CY=ddply(regq.dat,c("CY","Alt",'Region'),summarise,TFlow.kAcft=sum(cfs.to.acftd(FLOW),na.rm=T)/1000)
regq.dat.CY.mean=reshape2::dcast(regq.dat.CY,Alt~Region,value.var = "TFlow.kAcft",mean)
regq.dat.CY.mean=regq.dat.CY.mean[match(alts,regq.dat.CY.mean$Alt),]
regq.dat.CY.mean=regq.dat.CY.mean[,c("Alt","WCAs","Cal","StL","LWLagoon")]
regq.dat.CY.mean=regq.dat.CY.mean[match(alts,regq.dat.CY.mean$Alt),]

tmp=regq.dat.CY.mean[,c("WCAs","Cal","StL","LWLagoon")]
rownames(tmp)<-alts
cols2=rev(wesanderson::wes_palette("Zissou1",4,"continuous"))

ylim.val=c(0,1200);by.y=400;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
# png(filename=paste0(plot.path,"AvgFloodControl.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,2,0.25,1),oma=c(1,2,0.5,0.25),lwd=0.5);
layout(matrix(c(1:2),1,2,byrow=T),heights=c(1,0.4))

x=barplot(t(tmp),beside=F,col=NA,border=NA,ylim=ylim.val,axes=F,ann=F,names.arg = rep(NA,n.alts))
abline(h=ymaj,lty=1,col=adjustcolor("grey",0.5))
x=barplot(t(tmp),beside=F,col=cols2,ylim=ylim.val,axes=F,ann=F,names.arg = rep(NA,n.alts),add=T)
with(regq.dat.CY.mean,text(x,WCAs/2,round(WCAs,0),cex=0.75,col="white"))
with(regq.dat.CY.mean,text(x,WCAs+(((Cal+WCAs)-WCAs)/2),round(regq.dat.CY.mean$Cal,0),cex=0.75))
with(regq.dat.CY.mean,text(x,(WCAs+Cal)+(((Cal+WCAs+StL)-(Cal+WCAs))/2),round(regq.dat.CY.mean$StL,0),cex=0.75))
with(regq.dat.CY.mean,text(x,Cal+WCAs+StL+LWLagoon,round(regq.dat.CY.mean$LWLagoon,0),pos=3,offset=0.1,cex=0.75))
axis_fun(2,ymaj,ymin,ymaj)
axis_fun(1,x,x,alts,line=-0.5,las=1,cex=0.8);box(lwd=1)
mtext(side=2,line=3,"Discharge Volume (x1000 Ac-Ft Y\u207B\u00B9)")
mtext(side=1,line=1.5,"Alternatives")

# par(family="serif",mar=c(0,2,0.25,1))
plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA,xaxs="i",yaxs="i")
legend(-0.15,0.75,legend=c("Water Conservation Areas (S351 & S354)","Caloosahatchee River (S77)","St. Lucie River (S308)","Lake Worth Lagoon (S271)"),
       pch=22,
       lty=0,lwd=0.01,
       col="black",
       pt.bg=cols2,
       pt.cex=1.25,ncol=1,cex=0.9,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0,yjust=1)
text(1,0.2,"Mean annual flood control releases\nfrom Lake Okeechobee for the\n41 year (1965 - 2005)simulation period of record.",adj=1,xpd=NA,cex=0.65,font=3)
dev.off()

# Lake Discharges ---------------------------------------------------------
q.dat.xtab$S79_GT2100=with(q.dat.xtab,ifelse(S79>=2100,1,0))
q.dat.xtab$S80_GT1400=with(q.dat.xtab,ifelse(S80>=1400,1,0))
q.dat.xtab$S79_LT2100=with(q.dat.xtab,ifelse(S79<2100,1,0))
q.dat.xtab$S80_LT1400=with(q.dat.xtab,ifelse(S80<1400,1,0))

CRE.GT2100_annual=ddply(q.dat.xtab,c("CY","Alt"),summarise,Q.Lake.GT2100=sum(cfs.to.acftd(S79_QPFCSOURCE_LAKE[S79_GT2100==1])/1000,na.rm=T))
CRE.GT2100_annual.mean=ddply(CRE.GT2100_annual,"Alt",summarise,mean.val=mean(Q.Lake.GT2100))
CRE.GT2100_annual.mean=CRE.GT2100_annual.mean[match(alts,CRE.GT2100_annual.mean$Alt),]
CRE.GT2100_annual.mean=CRE.GT2100_annual.mean[match(alts,CRE.GT2100_annual.mean$Alt),]

SLE.GT1400_annual=ddply(q.dat.xtab,c("CY","Alt"),summarise,Q.Lake.GT1400=sum(cfs.to.acftd(S80_QPFCSOURCE_LAKE[S80_GT1400==1])/1000,na.rm=T))
SLE.GT1400_annual.mean=ddply(SLE.GT1400_annual,"Alt",summarise,mean.val=mean(Q.Lake.GT1400))
SLE.GT1400_annual.mean=SLE.GT1400_annual.mean[match(alts,SLE.GT1400_annual.mean$Alt),]
SLE.GT1400_annual.mean=SLE.GT1400_annual.mean[match(alts,SLE.GT1400_annual.mean$Alt),]

# png(filename=paste0(plot.path,"Lakedischarge_Annualmean.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
layout(matrix(c(1:2),1,2,byrow=T))
par(family="serif",mar=c(2,2,0.25,1),oma=c(4,3,2,0.25),lwd=0.5);

ylim.val=c(0,600);by.y=100;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=barplot(CRE.GT2100_annual.mean$mean.val,col=NA,border=NA,ylim=ylim.val,space=c(0),axes=F,ann=F)
abline(h=ymaj,lty=1,col=adjustcolor("grey",0.5))
x=barplot(CRE.GT2100_annual.mean$mean.val,
          col=adjustcolor(cols.alts,0.5),space=c(0),ylim=ylim.val,axes=F,ann=F,add=T)
axis_fun(2,ymaj,ymin,ymaj)
axis_fun(1,x,x,alts,cex=0.8,las=2)
box(lwd=1)
mtext(side=3,adj=0,"Mean Annual Discharge\n\u22652100 cfs at S-79",cex=0.75)
text(x,CRE.GT2100_annual.mean$mean.val,
     round(CRE.GT2100_annual.mean$mean.val,0),col="black",pos=1,cex=1)
# mtext(side=3,adj=1,"CY1965 - 2016",cex=0.75)
mtext(side=1,line=4,"Alternative")
mtext(side=2,line=2.5,"Lake Discharge\n(x1000 Ac-Ft Yr\u207B\u00B9)")
mtext(side=3,line=-1.25,adj=0," CRE")
# mtext(side=1,line=3,adj=1,"Flow Tag: S79_QPFCSOURCE_LAKE",cex=0.5,col=adjustcolor("black",0.5),font=3)

ylim.val=c(0,200);by.y=50;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=barplot(SLE.GT1400_annual.mean$mean.val,col=NA,border=NA,ylim=ylim.val,space=c(0),axes=F,ann=F)
abline(h=ymaj,lty=1,col=adjustcolor("grey",0.5))
x=barplot(SLE.GT1400_annual.mean$mean.val,
          col=adjustcolor(cols.alts,0.5),space=c(0),ylim=ylim.val,axes=F,ann=F,add=T)
axis_fun(2,ymaj,ymin,ymaj)
axis_fun(1,x,x,alts,cex=0.8,las=2)
box(lwd=1)
mtext(side=3,adj=0,"Mean Annual Discharge\n\u22651400 cfs at S-80",cex=0.75)
mtext(side=3,line=-1.25,adj=0," SLE")
text(x,SLE.GT1400_annual.mean$mean.val,round(SLE.GT1400_annual.mean$mean.val,0),col="black",pos=1,cex=1)
mtext(side=3,adj=1,"CY1965 - 2005",cex=0.75)
mtext(side=1,line=4,"Alternative")
dev.off()

LOK.SLE.annQ=ddply(q.dat.xtab,c("CY","Alt"),summarise,TFlow.acft=sum(cfs.to.acftd(S80_QPFCSOURCE_LAKE),na.rm=T))
LOK.SLE.annQ.sum=ddply(LOK.SLE.annQ,"Alt",summarise,meanQ=round(mean(TFlow.acft/1000,na.rm=T),0))
LOK.SLE.annQ.sum=LOK.SLE.annQ.sum[match(alts,LOK.SLE.annQ.sum$Alt),]
LOK.SLE.annQ.sum

LOK.CRE.annQ=ddply(q.dat.xtab,c("CY","Alt"),summarise,TFlow.acft=sum(cfs.to.acftd(S79_QPFCSOURCE_LAKE),na.rm=T))
LOK.CRE.annQ.sum=ddply(LOK.CRE.annQ,"Alt",summarise,meanQ=round(mean(TFlow.acft/1000,na.rm=T),0))
LOK.CRE.annQ.sum=LOK.CRE.annQ.sum[match(alts,LOK.CRE.annQ.sum$Alt),]
LOK.CRE.annQ.sum


# RECOVER Salinity Envelope -----------------------------------------------
## CRE
q.dat.xtab$S79_QPFCSOURCE_LAKE.14d=with(q.dat.xtab,ave(S79_QPFCSOURCE_LAKE,Alt,FUN=function(x) c(rep(NA,13),rollapply(x,width=14,FUN=function(x)mean(x,na.rm=T)))))
q.dat.xtab$CRE.low=with(q.dat.xtab,ifelse(S79.14d<750,1,0)) # RECOVER Low
q.dat.xtab$CRE.low1=with(q.dat.xtab,ifelse(S79.14d<457,1,0))
q.dat.xtab$CRE.low2=with(q.dat.xtab,ifelse(S79.14d>=457&S79.14d<750,1,0))
q.dat.xtab$CRE.opt=with(q.dat.xtab,ifelse(S79.14d>=750&S79.14d<2100,1,0)) # RECOVER Optimum
q.dat.xtab$CRE.high=with(q.dat.xtab,ifelse(S79.14d>=2100&S79.14d<2600,1,0)) # RECOVER Stress
q.dat.xtab$CRE.high1=with(q.dat.xtab,ifelse(S79.14d>=2600&S79.14d<4500,1,0))
q.dat.xtab$CRE.high2=with(q.dat.xtab,ifelse(S79.14d>=4500&S79.14d<6500,1,0))
q.dat.xtab$CRE.high3=with(q.dat.xtab,ifelse(S79.14d>=6500,1,0))
q.dat.xtab$CRE.dam=with(q.dat.xtab,ifelse(S79.14d>=2600,1,0)) # RECOVER Damaging

## SLE
q.dat.xtab$SLE.S80trib=rowSums(q.dat.xtab[,c("S80","TMC2EST","S48","S49","NSF2EST","sle_gw",'STA2TMC')],na.rm=T)
q.dat.xtab$SLE.S80trib.14d=with(q.dat.xtab,ave(SLE.S80trib,Alt,FUN=function(x) c(rep(NA,13),rollapply(x,width=14,FUN=function(x)mean(x,na.rm=T)))))
q.dat.xtab$S80_QPFCSOURCE_LAKE.14d=with(q.dat.xtab,ave(S80_QPFCSOURCE_LAKE,Alt,FUN=function(x) c(rep(NA,13),rollapply(x,width=14,FUN=function(x)mean(x,na.rm=T)))))
q.dat.xtab$S308_QFC.14d=with(q.dat.xtab,ave(S308_QFC,Alt,FUN=function(x) c(rep(NA,13),rollapply(x,width=14,FUN=function(x)mean(x,na.rm=T)))))

# 2007 Sal Metric
q.dat.xtab$sltrib=rowSums(q.dat.xtab[,c("TMC2EST","S48","S49","NSF2EST","sle_gw",'STA2TMC')],na.rm=T)
q.dat.xtab$slbsn=with(q.dat.xtab,S80-S308_QFC+sltrib)
q.dat.xtab$slbsn.14d=with(q.dat.xtab,ave(slbsn,Alt,FUN=function(x) c(rep(NA,13),rollapply(x,width=14,FUN=function(x)mean(x,na.rm=T)))))
q.dat.xtab$sltotal=with(q.dat.xtab,S80+sltrib)
q.dat.xtab$sltotal.14d=with(q.dat.xtab,ave(sltotal,Alt,FUN=function(x) c(rep(NA,13),rollapply(x,width=14,FUN=function(x)mean(x,na.rm=T)))))

# 
q.dat.xtab$SLE.low=with(q.dat.xtab,ifelse(SLE.S80trib.14d<150,1,0)) # RECOVER Low
q.dat.xtab$SLE.opt=with(q.dat.xtab,ifelse(SLE.S80trib.14d>=150&SLE.S80trib.14d<1400,1,0)) # RECOVER Optimum
q.dat.xtab$SLE.high=with(q.dat.xtab,ifelse(SLE.S80trib.14d>=1400&SLE.S80trib.14d<1700,1,0)) # RECOVER stress
q.dat.xtab$SLE.dam=with(q.dat.xtab,ifelse(SLE.S80trib.14d>=1700,1,0)) # RECOVER damaging
q.dat.xtab$SLE.high1=with(q.dat.xtab,ifelse(SLE.S80trib.14d>=1700&SLE.S80trib.14d<4000,1,0))
q.dat.xtab$SLE.high2=with(q.dat.xtab,ifelse(SLE.S80trib.14d>=4000,1,0))
q.dat.xtab$SLE.2000=with(q.dat.xtab,ifelse(sltotal.14d>2000,1,0))
q.dat.xtab$SLE.2000.bsn=with(q.dat.xtab,ifelse(slbsn.14d>2000,1,0))

##
q.dat.xtab$CRE.low.count=0
q.dat.xtab$CRE.low.LOK.count=0
q.dat.xtab$CRE.low.basin.count=0
q.dat.xtab$CRE.low1.count=0
q.dat.xtab$CRE.low2.count=0
q.dat.xtab$CRE.opt.count=0
q.dat.xtab$CRE.opt.LOK.count=0
q.dat.xtab$CRE.opt.basin.count=0
q.dat.xtab$CRE.high.count=0
q.dat.xtab$CRE.high_2100.count=0
q.dat.xtab$CRE.high1.count=0
q.dat.xtab$CRE.high2.count=0
q.dat.xtab$CRE.high3.count=0
q.dat.xtab$CRE.high.LOK.count=0
q.dat.xtab$CRE.high.basin.count=0
q.dat.xtab$CRE.dam.count=0
q.dat.xtab$CRE.dam.LOK.count=0
q.dat.xtab$CRE.dam.basin.count=0
q.dat.xtab$CRE.high3.LOK.count=0
q.dat.xtab$CRE.high3.basin.count=0

q.dat.xtab$SLE.low.count=0
q.dat.xtab$SLE.low.LOK.count=0
q.dat.xtab$SLE.low.basin.count=0
q.dat.xtab$SLE.opt.count=0
q.dat.xtab$SLE.opt.LOK.count=0
q.dat.xtab$SLE.opt.basin.count=0
q.dat.xtab$SLE.high.count=0

q.dat.xtab$SLE.high1.count=0
q.dat.xtab$SLE.high2.count=0
q.dat.xtab$SLE.high.LOK.count=0
q.dat.xtab$SLE.high.basin.count=0
q.dat.xtab$SLE.dam.count=0
q.dat.xtab$SLE.dam.LOK.count=0
q.dat.xtab$SLE.dam.basin.count=0
q.dat.xtab$SLE.high2.LOK.count=0
q.dat.xtab$SLE.high2.basin.count=0

###
q.dat.xtab2=data.frame()
for(j in 1:length(alts)){
  tmp=subset(q.dat.xtab,Alt==alts[j])
  for(i in 14:nrow(tmp)){
    ## CRE
    tmp$CRE.low.count[i]=with(tmp,ifelse(CRE.low[i]==1&sum(CRE.low.count[(i-13):(i-1)],na.rm=T)==0,1,0))  
    tmp$CRE.low.LOK.count[i]=with(tmp,
                                  ifelse(CRE.low.count[i]==1,
                                         ifelse((S79.14d[i]-S79_QPFCSOURCE_LAKE.14d[i])<=750,1,0),0))
    tmp$CRE.low.basin.count[i]=with(tmp,CRE.low.count[i]-CRE.low.LOK.count[i])
    
    tmp$CRE.low1.count[i]=with(tmp,ifelse(CRE.low1[i]==1&sum(CRE.low1.count[(i-13):(i-1)],na.rm=T)==0,1,0))  
    tmp$CRE.low2.count[i]=with(tmp,ifelse(CRE.low2[i]==1&sum(CRE.low2.count[(i-13):(i-1)],na.rm=T)==0,1,0))  
    tmp$CRE.opt.count[i]=with(tmp,ifelse(CRE.opt[i]==1&sum(CRE.opt.count[(i-13):(i-1)],na.rm=T)==0,1,0))
    tmp$CRE.opt.LOK.count[i]=with(tmp,
                                  ifelse(CRE.opt.count[i]==1,
                                         ifelse((S79.14d[i]-S79_QPFCSOURCE_LAKE.14d[i])<=750,1,0),0))
    tmp$CRE.opt.basin.count[i]=with(tmp,CRE.opt.count[i]-CRE.opt.LOK.count[i])
    
    tmp$CRE.high.count[i]=with(tmp,ifelse(CRE.high[i]==1&sum(CRE.high.count[(i-13):(i-1)],na.rm=T)==0,1,0))
    tmp$CRE.high1.count[i]=with(tmp,ifelse(CRE.high1[i]==1&sum(CRE.high1.count[(i-13):(i-1)],na.rm=T)==0,1,0))
    tmp$CRE.high2.count[i]=with(tmp,ifelse(CRE.high2[i]==1&sum(CRE.high2.count[(i-13):(i-1)],na.rm=T)==0,1,0))
    tmp$CRE.high3.count[i]=with(tmp,ifelse(CRE.high3[i]==1&sum(CRE.high3.count[(i-13):(i-1)],na.rm=T)==0,1,0))
    tmp$CRE.high.LOK.count[i]=with(tmp,
                                   ifelse(CRE.high.count[i]==1,
                                          ifelse((S79.14d[i]-S79_QPFCSOURCE_LAKE.14d[i])<=2100,1,0),0))
    tmp$CRE.high.basin.count[i]=with(tmp,CRE.high.count[i]-CRE.high.LOK.count[i])
    tmp$CRE.dam.count[i]=with(tmp,ifelse(CRE.dam[i]==1&sum(CRE.dam.count[(i-13):(i-1)],na.rm=T)==0,1,0))
    tmp$CRE.dam.LOK.count[i]=with(tmp,
                                  ifelse(CRE.dam.count[i]==1,
                                         ifelse((S79.14d[i]-S79_QPFCSOURCE_LAKE.14d[i])<=2600,1,0),0))
    tmp$CRE.dam.basin.count[i]=with(tmp,CRE.dam.count[i]-CRE.dam.LOK.count[i])
    
    tmp$CRE.high3.LOK.count[i]=with(tmp,
                                    ifelse(CRE.high3.count[i]==1,
                                           ifelse((S79.14d[i]-S79_QPFCSOURCE_LAKE.14d[i])<=6500,1,0),0))
    tmp$CRE.high3.basin.count[i]=with(tmp,CRE.high3.count[i]-CRE.high3.LOK.count[i])
    
    ## SLE
    tmp$SLE.low.count[i]=with(tmp,ifelse(SLE.low[i]==1&sum(SLE.low.count[(i-13):(i-1)],na.rm=T)==0,1,0)) 
    tmp$SLE.low.LOK.count[i]=with(tmp,
                                  ifelse(SLE.low.count[i]==1,
                                         ifelse((SLE.S80trib.14d[i]-S80_QPFCSOURCE_LAKE.14d[i])<=150,1,0),0))
    tmp$SLE.low.basin.count[i]=with(tmp,SLE.low.count[i]-SLE.low.LOK.count[i])
    tmp$SLE.opt.count[i]=with(tmp,ifelse(SLE.opt[i]==1&sum(SLE.opt.count[(i-13):(i-1)],na.rm=T)==0,1,0))  
    tmp$SLE.opt.LOK.count[i]=with(tmp,
                                  ifelse(SLE.opt.count[i]==1,
                                         ifelse((SLE.S80trib.14d[i]-S80_QPFCSOURCE_LAKE.14d[i])<=150,1,0),0))
    tmp$SLE.opt.basin.count[i]=with(tmp,SLE.opt.count[i]-SLE.opt.LOK.count[i])
    tmp$SLE.high.count[i]=with(tmp,ifelse(SLE.high[i]==1&sum(SLE.high.count[(i-13):(i-1)],na.rm=T)==0,1,0))
    tmp$SLE.high.LOK.count[i]=with(tmp,
                                   ifelse(SLE.high.count[i]==1,
                                          ifelse((SLE.S80trib.14d[i]-S80_QPFCSOURCE_LAKE.14d[i])<=1400,1,0),0))
    tmp$SLE.high.basin.count[i]=with(tmp,SLE.high.count[i]-SLE.high.LOK.count[i])
    tmp$SLE.dam.count[i]=with(tmp,ifelse(SLE.dam[i]==1&sum(SLE.dam.count[(i-13):(i-1)],na.rm=T)==0,1,0))  
    tmp$SLE.dam.LOK.count[i]=with(tmp,
                                  ifelse(SLE.dam.count[i]==1,
                                         ifelse((SLE.S80trib.14d[i]-S80_QPFCSOURCE_LAKE.14d[i])<=1700,1,0),0))
    tmp$SLE.dam.basin.count[i]=with(tmp,SLE.dam.count[i]-SLE.dam.LOK.count[i])
    tmp$SLE.high1.count[i]=with(tmp,ifelse(SLE.high1[i]==1&sum(SLE.high1.count[(i-13):(i-1)],na.rm=T)==0,1,0))
    tmp$SLE.high2.count[i]=with(tmp,ifelse(SLE.high2[i]==1&sum(SLE.high2.count[(i-13):(i-1)],na.rm=T)==0,1,0))
    tmp$SLE.high2.LOK.count[i]=with(tmp,
                                    ifelse(SLE.high2.count[i]==1,
                                           ifelse((SLE.S80trib.14d[i]-S80_QPFCSOURCE_LAKE.14d[i])<=4000,1,0),0))
    tmp$SLE.high2.basin.count[i]=with(tmp,SLE.high2.count[i]-SLE.high2.LOK.count[i])
    
    # tmp$SLE.2000.count[i]=with(tmp,ifelse(SLE.2000[i]==1&sum(SLE.2000.count[(i-13):(i-1)],na.rm=T)==0,1,0))
    # tmp$SLE.2000.count.bsn[i]=with(tmp,ifelse(SLE.2000.bsn[i]==1&sum(SLE.2000.count.bsn[(i-13):(i-1)],na.rm=T)==0,1,0))
    # tmp$SLE.2000.count.LOK[i]=with(tmp,SLE.2000.count[i]-SLE.2000.count.bsn[i])
    # tmp$SLE.2000.count.LOK[i]=with(tmp,
    #                                 ifelse(SLE.2000.count[i]==1,
    #                                        ifelse((SLE.S80trib.14d[i]-S308_QFC.14d[i])<=2000,1,0),0))
    # tmp$SLE.2000.count.bsn[i]=with(tmp,SLE.2000.count[i]-SLE.2000.count.LOK[i])
  }
  q.dat.xtab2=rbind(q.dat.xtab2,tmp)
  print(j)
}


# Old RECOVER Salinity envelope based on monthly Q
# SLE
q.dat.xtab$SLE.2000.count=0
q.dat.xtab$SLE.2000.count.LOK=0
q.dat.xtab$SLE.2000.count.bsn=0
q.dat.xtab3=data.frame()
for(j in 1:length(alts)){
  tmp=subset(q.dat.xtab,Alt==alts[j])
  for(i in 14:nrow(tmp)){
    tmp$SLE.2000.count[i]=with(tmp,ifelse(sum(SLE.2000[(i-13):i],na.rm=T)==14&sum(SLE.2000.count[(i-13):(i-1)],na.rm=T)==0,1,0))
    tmp$SLE.2000.count.bsn[i]=with(tmp,ifelse(sum(SLE.2000.bsn[(i-13):i],na.rm=T)==14&sum(SLE.2000.count.bsn[(i-13):(i-1)],na.rm=T)==0,1,0))
  }
  q.dat.xtab3=rbind(q.dat.xtab3,tmp)
  print(j)
}

SLE_SalEnv2000=ddply(q.dat.xtab3,"Alt",summarise,
                     SLE.2000.total=sum(SLE.2000.count,na.rm = T),
                     SLE.2000.bsn=sum(SLE.2000.count.bsn,na.rm = T))
SLE_SalEnv2000$SLE.2000.lok=with(SLE_SalEnv2000,SLE.2000.total-SLE.2000.bsn)

SLE.sal.env=ddply(q.dat.xtab,c("Alt","CY","month"),summarise,month.Q=mean(SLE.S80trib,na.rm=T))
SLE.sal.env$lowQ=with(SLE.sal.env,ifelse(month.Q<=350,1,0))
SLE.sal.env$highQ=with(SLE.sal.env,ifelse(month.Q>=2000&month.Q<3000,1,0))
SLE.sal.env$ext.highQ=with(SLE.sal.env,ifelse(month.Q>3000,1,0))
SLE.sal.env.sum=ddply(SLE.sal.env,"Alt",summarise,
                      N.low=sum(lowQ,na.rm=T),
                      N.highQ=sum(highQ,na.rm=T),
                      N.ext.highQ=sum(ext.highQ,na.rm=T))
SLE.sal.env.sum=SLE.sal.env.sum[match(alts,SLE.sal.env.sum$Alt),]
SLE.sal.env.sum

SLE.sal.env.sum2=merge(SLE.sal.env.sum,SLE_SalEnv2000[,c("Alt","SLE.2000.bsn","SLE.2000.lok")],"Alt")
vars=c("Alt", "N.low","SLE.2000.bsn", "SLE.2000.lok", "N.highQ", "N.ext.highQ")
SLE.sal.env.sum2=SLE.sal.env.sum2[match(alts,SLE.sal.env.sum2$Alt),vars]



# CRE
CRE.sal.env=ddply(q.dat.xtab,c("Alt","CY","month"),summarise,month.Q=mean(S79,na.rm=T))
CRE.sal.env$lowQ=with(CRE.sal.env,ifelse(month.Q<450,1,0))
CRE.sal.env$highQ=with(CRE.sal.env,ifelse(month.Q>2800,1,0))
CRE.sal.env$ext.highQ=with(CRE.sal.env,ifelse(month.Q>4500,1,0))

CRE.sal.env.sum=ddply(CRE.sal.env,c("Alt"),summarise,
      N.lowQ=sum(lowQ),
      N.highQ=sum(highQ),
      N.ext.highQ=sum(ext.highQ))
CRE.sal.env.sum=CRE.sal.env.sum[match(alts,CRE.sal.env.sum$Alt),]

CRE.sal.env.consec=data.frame()
for(j in 1:length(alts)){
  tmp=subset(CRE.sal.env,Alt==alts[j])
  tmp$highQ.N=0
  for(i in 2:nrow(tmp)){
    tmp$highQ.N[i]=with(tmp,ifelse(highQ[i-1]==0&highQ[i]>0,1,
                                   ifelse(highQ[i-1]>0&highQ[i]>0,1,0)))
    
  }
  CRE.highQ=consec.startend(tmp$highQ.N>0)
  tmp$sum.highQ=0
  for(i in 1:length(CRE.highQ$ends)){
    tmp[CRE.highQ$ends[i],]$sum.highQ=with(tmp[c(CRE.highQ$starts[i]:CRE.highQ$ends[i]),],sum(highQ.N,na.rm=T))
  }
  CRE.sal.env.consec=rbind(tmp,CRE.sal.env.consec)
  print(j)
}
rslt.CREHighQ=ddply(CRE.sal.env.consec,c("Alt","sum.highQ"),summarise,count.event=N.obs(sum.highQ))
## works!! consistent with old # of consecutive months plot  

# Trying to figure out consecutive 14-day period
# length(seq(date.fun("1965-01-01"),date.fun("2005-12-31"),"1 days"))
# length(seq(date.fun("1965-01-01"),date.fun("2005-12-31"),"14 days"))
# 
# length(seq(date.fun("1965-01-01"),date.fun("2016-12-31"),"1 days"))
# length(seq(date.fun("1965-01-01"),date.fun("2016-12-31"),"14 days"))
# length(seq(date.fun("1965-01-01"),date.fun("2016-12-31"),"1 months"))
# 
# test=data.frame(date=seq(date.fun("1965-01-01"),date.fun("2016-12-31"),"1 days"))
# test$p14=as.numeric(format(test$date,"%j"))%/%14L+1L
# tmp=data.frame(Date=seq(date.fun("1965-01-01"),date.fun("2016-12-31"),"1 day"))
# # tmp$D14.period=as.numeric(format(tmp$Date,'%j'))%/%15L+1L
# tmp$D14.period=(1:nrow(tmp))%/%14L+1L
# # nrow(tmp)
# length(unique(tmp$D14.period))
# 
# test=ddply(tmp,"D14.period",summarise,N.val=N.obs(D14.period))
# tail(test)

##
vars=c(paste0("CRE.",c("low.count","opt.count","high.LOK.count","high.basin.count","dam.LOK.count","dam.basin.count","high3.count")),
       paste0("SLE.",c("low.count","opt.count","high.LOK.count","high.basin.count","dam.LOK.count","dam.basin.count","high2.count")))

SalEnv_count.melt=reshape2::melt(q.dat.xtab2[,c("Alt",vars)],id.vars = "Alt")
SalEnv_count=reshape2::dcast(SalEnv_count.melt,Alt~variable,value.var = "value",sum)
SalEnv_count=SalEnv_count[match(alts,SalEnv_count$Alt),]
SalEnv_count

SalEnv_count$SLE.stress_dam.LOK=with(SalEnv_count,SLE.high.LOK.count+SLE.dam.LOK.count)
SalEnv_count$CRE.stress_dam.LOK=with(SalEnv_count,CRE.high.LOK.count+CRE.dam.LOK.count)
SalEnv_count$per_stressdam.lok=with(SalEnv_count,((CRE.stress_dam.LOK-CRE.stress_dam.LOK[2])/CRE.stress_dam.LOK[2])*100)
SalEnv_count

SalEnv_count.LOSOM=data.frame(
  Alt=c("ECB19","NA25f","PA25"),
  CRE.stress_dam.LOK=c(184+198,183+187,66+81),
  SLE.stress_dam.LOK=c(164+159,145+140,30+38),
  CRE.opt=c(469,588,765),
  SLE.opt=c(831,865,912)
)

plot(CRE.stress_dam.LOK~SLE.stress_dam.LOK,SalEnv_count,ylim=c(50,400),xlim=c(50,400))
abline(0,1)
text(CRE.stress_dam.LOK~SLE.stress_dam.LOK,SalEnv_count,Alt,pos=3,offset=0.5)
points(CRE.stress_dam.LOK~SLE.stress_dam.LOK,SalEnv_count.LOSOM,pch=22,bg="forestgreen")
text(CRE.stress_dam.LOK~SLE.stress_dam.LOK,SalEnv_count.LOSOM,Alt,pos=3,offset=0.5)

# png(filename=paste0(plot.path,"RECOVER_CRESLE_opt.png"),width=5,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,2,1,1),oma=c(2,3,0.5,1),lwd=0.5);
ylim.val=c(0,1100);by.y=250;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=ylim.val;xmaj=ymaj;xmin=ymin
plot(CRE.opt.count~SLE.opt.count,SalEnv_count,ylim=xlim.val,xlim=xlim.val,type="n",ann=F,axes=F)
abline(h=ymaj,v=xmaj,lty=2,col=adjustcolor("grey",0.5))
abline(0,1)
points(CRE.opt.count~SLE.opt.count,SalEnv_count,pch=21,bg="indianred1",cex=1.25,lwd=0.1)
text(CRE.opt.count~SLE.opt.count,subset(SalEnv_count,!Alt%in%c("ECB","ASR","ALT1BWR")),Alt,pos=4,cex=0.75,offset=0.5,col="indianred1")
text(CRE.opt.count~SLE.opt.count,subset(SalEnv_count,Alt%in%c("ECB","ASR")),Alt,pos=2,cex=0.75,offset=0.5,col="indianred1")
text(CRE.opt.count~SLE.opt.count,subset(SalEnv_count,Alt%in%c("ALT1BWR")),Alt,pos=1,cex=0.75,offset=0.5,col="indianred1")
points(CRE.opt~SLE.opt,SalEnv_count.LOSOM,pch=22,bg="forestgreen",cex=1.25,lwd=0.1)
text(CRE.opt~SLE.opt,SalEnv_count.LOSOM,Alt,pos=2,col="forestgreen",cex=0.75,offset=0.5)
axis_fun(2,ymaj,ymin,ymaj,line=-0.5)
axis_fun(1,xmaj,xmin,xmaj);box(lwd=1)
mtext(side=2,line=3,"CRE Optimum Flow Events")
mtext(side=1,line=2,"SLE Optimum Flow Events")
mtext(side=3,adj=0,"Based on 14-day event count")
legend("topleft",legend=c("LOWRP (1965 - 2005)","LOSOM (1965 - 2016)"),
       pch=c(21,22),
       lty=c(0,0),lwd=0.1,
       col=c("black"),
       pt.bg=c("indianred1","forestgreen"),
       pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

vars.CRE=paste("CRE",c("low.count","opt.count","high.LOK.count","dam.LOK.count","high3.count"),sep=".")
vars.SLE=paste("SLE",c("low.count","opt.count","high.LOK.count","dam.LOK.count","high2.count"),sep=".")
CRE.SalEnv_count3=SalEnv_count[,c("Alt",vars.CRE)]
SLE.SalEnv_count3=SalEnv_count[,c("Alt",vars.SLE)]
CRE.labs=c("Low Flow\n(<750 cfs)","Optimum\n(750 - 2100 cfs)","Stress From LOK\n(2100 - 2600 cfs)","Damaging From LOK\n(>2600 cfs)","Extreme\n(>6500 cfs)")
SLE.labs=c("Low Flow\n(<150 cfs)","Optimum\n(150 - 1400 cfs)","Stress From LOK\n(1400 - 1700 cfs)","Damaging From LOK\n(>1700 cfs)","Extreme\n(>4000 cfs)")
# png(filename=paste0(plot.path,"RECOVER_CRE_SalEnv_.png"),width=7,height=3,units="in",res=200,type="windows",bg="white")
layout(matrix(c(1:5),1,5,byrow=T))
par(family="serif",mar=c(2,2,1,1),oma=c(3,3,2,1),lwd=0.5);

ymax=c(800,500,100,120,70)
yval=ymax/2
for(i in 2:6){
  ylim.val=c(0,ymax[i-1]);ymaj=seq(ylim.val[1],ylim.val[2],yval[i-1]);ymin=seq(ylim.val[1],ylim.val[2],yval[i-1]/2)
  x=barplot(CRE.SalEnv_count3[,i],col=adjustcolor(cols.alts,0.5),ylim=ylim.val,axes=F,ann=F,space=0)
  axis_fun(2,ymaj,ymin,ymaj)
  axis_fun(1,x,x,alts,cex=0.7,las=2)
  # abline(v=c(x[1]+(x[2]-x[1])/2,x[3]+(x[4]-x[3])/2),lwd=1)
  box(lwd=1)
  mtext(side=3,adj=0,CRE.labs[i-1],cex=0.7)
  text(x,CRE.SalEnv_count3[,i],
       round(CRE.SalEnv_count3[,i],0),font=2,col="black",pos=3,cex=0.75,offset=0.25)
  # if(i==2){mtext(side=3,adj=0,line=-1.25," CRE")}
}
mtext(side=4,line=0.5,"Caloosahatchee")
mtext(side=1,line=1.75,outer=T,"Alternative")
mtext(side=2,line=0.75,outer=T,"Count of 14-Day Periods")
dev.off()

# png(filename=paste0(plot.path,"RECOVER_SLE_SalEnv_.png"),width=7,height=3,units="in",res=200,type="windows",bg="white")
layout(matrix(c(1:5),1,5,byrow=T))
par(family="serif",mar=c(2,2,1,1),oma=c(3,3,2,1),lwd=0.5);

ymax=c(30,1000,80,120,70)
yval=ymax/2
for(i in 2:6){
  ylim.val=c(0,ymax[i-1]);ymaj=seq(ylim.val[1],ylim.val[2],yval[i-1]);ymin=seq(ylim.val[1],ylim.val[2],yval[i-1]/2)
  x=barplot(SLE.SalEnv_count3[,i],col=adjustcolor(cols.alts,0.5),ylim=ylim.val,axes=F,ann=F,space=0)
  axis_fun(2,ymaj,ymin,ymaj)
  axis_fun(1,x,x,alts,cex=0.7,las=2)
  # abline(v=c(x[1]+(x[2]-x[1])/2,x[3]+(x[4]-x[3])/2),lwd=1)
  box(lwd=1)
  mtext(side=3,adj=0,SLE.labs[i-1],cex=0.7)
  text(x,SLE.SalEnv_count3[,i],
       round(SLE.SalEnv_count3[,i],0),font=2,col="black",pos=3,cex=0.75,offset=0.25)
  # if(i==2){mtext(side=3,adj=0,line=-1.25," CRE")}
}
mtext(side=4,line=0.5,"St Lucie")
mtext(side=1,line=1.75,outer=T,"Alternative")
mtext(side=2,line=0.75,outer=T,"Count of 14-Day Periods")
dev.off()

SalEnv_count.FWO=reshape2::dcast(SalEnv_count.melt,variable~Alt,value.var = "value",sum)
SalEnv_count.FWO=SalEnv_count.FWO[,c("variable",alts)]
SalEnv_count.FWO$ASR.PerFWO=with(SalEnv_count.FWO,((ASR-FWO)/FWO)*100)
SalEnv_count.FWO$ALT1BWR.PerFWO=with(SalEnv_count.FWO,((ALT1BWR-FWO)/FWO)*100)

# png(filename=paste0(plot.path,"RECOVER_CRE_SalEnv_FWO.png"),width=6.5,height=3.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,2,1,1),oma=c(2,2,1,1),lwd=0.5);
ylim.val=c(-30,20);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=barplot(subset(SalEnv_count.FWO,variable%in%vars.CRE)$ASR.PerFWO,
        ylim=ylim.val,space=0,
        axes=F,ann=F)
abline(h=0,lwd=1)
axis_fun(1,x,x,CRE.labs,line=0,cex=0.6)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
with(subset(SalEnv_count.FWO,variable%in%vars.CRE),
     text(x,ASR.PerFWO,
          round(ASR.PerFWO,1),
          pos=ifelse(ASR.PerFWO<0,1,3),offset=0.2))
mtext(side=3,adj=0,"Caloosahatchee Estuary\nAlternative: ASR")
mtext(side=1,line=2.5,outer=F,"Salinity Envelope Category")
mtext(side=2,line=2.5,outer=F,"% Difference relative to FWO")
dev.off()

# png(filename=paste0(plot.path,"RECOVER_SLE_SalEnv_FWO.png"),width=6.5,height=3.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,2,1,1),oma=c(2,2,1,1),lwd=0.5);
ylim.val=c(-25,15);by.y=5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=barplot(subset(SalEnv_count.FWO,variable%in%vars.SLE)$ASR.PerFWO,
          ylim=ylim.val,space=0,
          axes=F,ann=F)
abline(h=0,lwd=1)
axis_fun(1,x,x,SLE.labs,line=0,cex=0.6)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
with(subset(SalEnv_count.FWO,variable%in%vars.SLE),
     text(x,ASR.PerFWO,
          round(ASR.PerFWO,1),
          pos=ifelse(ASR.PerFWO<0,1,3),offset=0.2))
mtext(side=3,adj=0,"St Lucie Estuary\nAlternative: ASR")
mtext(side=1,line=2.5,outer=F,"Salinity Envelope Category")
mtext(side=2,line=2.5,outer=F,"% Difference relative to FWO")
dev.off()

CRE.sal.env.sum

vars.CRE=paste("CRE",c("low.count","opt.count","high.LOK.count","dam.LOK.count","high3.count"),sep=".")
CRE.SalEnv_count3=SalEnv_count[,c("Alt",vars.CRE)]
CRE.labs=c("Low Flow\n(<750 cfs)","Optimum\n(750 - 2100 cfs)","Stress From LOK\n(2100 - 2600 cfs)","Damaging From LOK\n(>2600 cfs)","Extreme\n(>6500 cfs)")
# png(filename=paste0(plot.path,"RECOVER_CRE_SalEnv_all.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
layout(matrix(c(1:8),2,4,byrow=T))
par(family="serif",mar=c(2,2,1,1),oma=c(3,3,2,1),lwd=0.5);

CRE.labs.old=c("Low Flow\n(<450 cfs Basin & LOK)","High Flow\n(>2800 cfs Basin & LOK)","Extreme High Flow\n(>4500 cfs Basin & LOK)")
ymax=c(150,100,100)
yval=ymax/2
for(i in 2:4){
  ylim.val=c(0,ymax[i-1]);ymaj=seq(ylim.val[1],ylim.val[2],yval[i-1]);ymin=seq(ylim.val[1],ylim.val[2],yval[i-1]/2)
  x=barplot(CRE.sal.env.sum[,i],col=adjustcolor(cols.alts,0.5),ylim=ylim.val,axes=F,ann=F,space=0)
  axis_fun(2,ymaj,ymin,ymaj)
  axis_fun(1,x,x,NA,cex=0.7,las=2)
  # abline(v=c(x[1]+(x[2]-x[1])/2,x[3]+(x[4]-x[3])/2),lwd=1)
  box(lwd=1)
  mtext(side=3,adj=0,CRE.labs.old[i-1],cex=0.5)
  text(x,CRE.sal.env.sum[,i],
       round(CRE.sal.env.sum[,i],0),font=2,col="black",pos=3,cex=0.75,offset=0.25)
  if(i==2){mtext(side=2,line=2.5,"Number of Exceedances")}
}

plot(0:1,0:1,ann=F,axes=F,type="n")
text(0,0.5,adj=0,cex=0.75,
     "2007 RECOVER Salinity Envelope\nPerformance Metric based on\nmean monthly discharge at S-79")

ymax=c(800,500,100,150)
yval=ymax/2
for(i in 2:5){
  ylim.val=c(0,ymax[i-1]);ymaj=seq(ylim.val[1],ylim.val[2],yval[i-1]);ymin=seq(ylim.val[1],ylim.val[2],yval[i-1]/2)
  x=barplot(CRE.SalEnv_count3[,i],col=adjustcolor(cols.alts,0.5),ylim=ylim.val,axes=F,ann=F,space=0)
  axis_fun(2,ymaj,ymin,ymaj)
  axis_fun(1,x,x,alts,cex=0.7,las=2)
  # abline(v=c(x[1]+(x[2]-x[1])/2,x[3]+(x[4]-x[3])/2),lwd=1)
  box(lwd=1)
  mtext(side=3,adj=0,CRE.labs[i-1],cex=0.5)
  text(x,CRE.SalEnv_count3[,i],
       round(CRE.SalEnv_count3[,i],0),font=2,col="black",pos=3,cex=0.75,offset=0.25)
  # if(i==2){mtext(side=3,adj=0,line=-1.25," CRE")}
  if(i==2){mtext(side=2,line=2.5,"Count of 14-Day Periods")}
}
mtext(side=4,line=0.5,"2020 RECOVER PM")
mtext(side=1,line=1.75,outer=T,"Alternative")
dev.off()


# png(filename=paste0(plot.path,"RECOVER_CRE_SalEnv_all_2.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
layout(matrix(c(1:8),2,4,byrow=T))
par(family="serif",mar=c(2,2,1,1),oma=c(3,3,2,1),lwd=0.5);

CRE.labs.old=c("Low Flow\n(<450 cfs Basin & LOK)","High Flow\n(>2800 cfs Basin & LOK)","Extreme High Flow\n(>4500 cfs Basin & LOK)")
ymax=c(150,100,100)
yval=ymax/2
for(i in 2:4){
  ylim.val=c(0,ymax[i-1]);ymaj=seq(ylim.val[1],ylim.val[2],yval[i-1]);ymin=seq(ylim.val[1],ylim.val[2],yval[i-1]/2)
  x=barplot(CRE.sal.env.sum[,i],col=adjustcolor(cols.alts,0.5),ylim=ylim.val,axes=F,ann=F,space=0)
  axis_fun(2,ymaj,ymin,ymaj)
  axis_fun(1,x,x,NA,cex=0.7,las=2)
  # abline(v=c(x[1]+(x[2]-x[1])/2,x[3]+(x[4]-x[3])/2),lwd=1)
  box(lwd=1)
  mtext(side=3,adj=0,CRE.labs.old[i-1],cex=0.5)
  text(x,CRE.sal.env.sum[,i],
       round(CRE.sal.env.sum[,i],0),font=2,col="black",pos=3,cex=0.75,offset=0.25)
  if(i==2){mtext(side=2,line=2.5,"Number of Exceedances")}
}
mtext(side=4,line=0.5,"2007 RECOVER PM")
plot(0:1,0:1,ann=F,axes=F,type="n")
# text(0,0.5,adj=0,cex=0.75,
#     "2007 RECOVER Salinity Envelope\nPerformance Metric based on\nmean monthly discharge at S-79")

ymax=c(800,500,100,150)
yval=ymax/2
for(i in 2:5){
  ylim.val=c(0,ymax[i-1]);ymaj=seq(ylim.val[1],ylim.val[2],yval[i-1]);ymin=seq(ylim.val[1],ylim.val[2],yval[i-1]/2)
  x=barplot(CRE.SalEnv_count3[,i],col=adjustcolor(cols.alts,0.5),ylim=ylim.val,axes=F,ann=F,space=0)
  axis_fun(2,ymaj,ymin,ymaj)
  axis_fun(1,x,x,alts,cex=0.7,las=2)
  # abline(v=c(x[1]+(x[2]-x[1])/2,x[3]+(x[4]-x[3])/2),lwd=1)
  box(lwd=1)
  mtext(side=3,adj=0,CRE.labs[i-1],cex=0.5)
  text(x,CRE.SalEnv_count3[,i],
       round(CRE.SalEnv_count3[,i],0),font=2,col="black",pos=3,cex=0.75,offset=0.25)
  # if(i==2){mtext(side=3,adj=0,line=-1.25," CRE")}
  if(i==2){mtext(side=2,line=2.5,"Count of 14-Day Periods")}
}
mtext(side=4,line=0.5,"2020 RECOVER PM")
mtext(side=1,line=1.75,outer=T,"Alternative")
dev.off()

SLE.sal.env.sum2
# png(filename=paste0(plot.path,"RECOVER_SLE_SalEnv_all.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
layout(matrix(c(1:10),2,5,byrow=T))
par(family="serif",mar=c(2,2,1,1),oma=c(3,3,2,1),lwd=0.5);

SLE.labs.old=c(
  "Low Flow\n(<350 cfs Basin & LOK)",
  "High Flow\n(>2000 cfs Basin - 14d)",
  "High Flow\n(>2000 cfs LOK - 14d)",
  "Optimum Flow\n(2000 - 3000 cfs\nBasin & LOK)",
  "Damaging Flow\n(>3000 cfs\nBasin & LOK)")
ymax=rep(120,5)
yval=ymax/2
for(i in 2:6){
  ylim.val=c(0,ymax[i-1]);ymaj=seq(ylim.val[1],ylim.val[2],yval[i-1]);ymin=seq(ylim.val[1],ylim.val[2],yval[i-1]/2)
  x=barplot(SLE.sal.env.sum2[,i],col=adjustcolor(cols.alts,0.5),ylim=ylim.val,axes=F,ann=F,space=0)
  axis_fun(2,ymaj,ymin,ymaj)
  axis_fun(1,x,x,NA,cex=0.7,las=2)
  # abline(v=c(x[1]+(x[2]-x[1])/2,x[3]+(x[4]-x[3])/2),lwd=1)
  box(lwd=1)
  mtext(side=3,adj=0,SLE.labs.old[i-1],cex=0.5)
  text(x,SLE.sal.env.sum2[,i],
       round(SLE.sal.env.sum2[,i],0),font=2,col="black",pos=3,cex=0.75,offset=0.25)
  if(i==2){mtext(side=2,line=2.5,"Number of Exceedances")}
}
mtext(side=4,line=0.5,"2007 RECOVER PM")

ymax=c(30,1000,80,120)
yval=ymax/2
for(i in 2:5){
  ylim.val=c(0,ymax[i-1]);ymaj=seq(ylim.val[1],ylim.val[2],yval[i-1]);ymin=seq(ylim.val[1],ylim.val[2],yval[i-1]/2)
  x=barplot(SLE.SalEnv_count3[,i],col=adjustcolor(cols.alts,0.5),ylim=ylim.val,axes=F,ann=F,space=0)
  axis_fun(2,ymaj,ymin,ymaj)
  axis_fun(1,x,x,alts,cex=0.7,las=2)
  # abline(v=c(x[1]+(x[2]-x[1])/2,x[3]+(x[4]-x[3])/2),lwd=1)
  box(lwd=1)
  mtext(side=3,adj=0,SLE.labs[i-1],cex=0.5)
  text(x,SLE.SalEnv_count3[,i],
       round(SLE.SalEnv_count3[,i],0),font=2,col="black",pos=3,cex=0.75,offset=0.25)
  # if(i==2){mtext(side=3,adj=0,line=-1.25," CRE")}
  if(i==2){mtext(side=2,line=2.5,outer=F,"Count of 14-Day Periods")}
}
mtext(side=4,line=0.5,"2020 RECOVER PM")
mtext(side=1,line=1.75,outer=T,"Alternative")

plot(0:1,0:1,ann=F,axes=F,type="n")
dev.off()


SLE.sal.env.sum2
# png(filename=paste0(plot.path,"RECOVER_SLE_SalEnv_all_2.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
layout(matrix(c(1:8),2,4,byrow=T))
par(family="serif",mar=c(2,2,1,1),oma=c(3,3,2,1),lwd=0.5);

SLE.labs.old=c(
  "Low Flow\n(<350 cfs Basin & LOK)",
  "High Flow\n(>2000 cfs Basin - 14d)",
  "High Flow\n(>2000 cfs LOK - 14d)",
  "Optimum Flow\n(2000 - 3000 cfs\nBasin & LOK)",
  "Damaging Flow\n(>3000 cfs\nBasin & LOK)")
ymax=rep(120,5)
yval=ymax/2
for(i in 2:4){
  ylim.val=c(0,ymax[i-1]);ymaj=seq(ylim.val[1],ylim.val[2],yval[i-1]);ymin=seq(ylim.val[1],ylim.val[2],yval[i-1]/2)
  x=barplot(SLE.sal.env.sum2[,i],col=adjustcolor(cols.alts,0.5),ylim=ylim.val,axes=F,ann=F,space=0)
  axis_fun(2,ymaj,ymin,ymaj)
  axis_fun(1,x,x,NA,cex=0.7,las=2)
  # abline(v=c(x[1]+(x[2]-x[1])/2,x[3]+(x[4]-x[3])/2),lwd=1)
  box(lwd=1)
  mtext(side=3,adj=0,SLE.labs.old[i-1],cex=0.5)
  text(x,SLE.sal.env.sum2[,i],
       round(SLE.sal.env.sum2[,i],0),font=2,col="black",pos=3,cex=0.75,offset=0.25)
  if(i==2){mtext(side=2,line=2.5,"Number of Exceedances")}
}
mtext(side=4,line=0.5,"2007 RECOVER PM")

plot(0:1,0:1,ann=F,axes=F,type="n")

ymax=c(30,1000,80,120)
yval=ymax/2
for(i in 2:5){
  ylim.val=c(0,ymax[i-1]);ymaj=seq(ylim.val[1],ylim.val[2],yval[i-1]);ymin=seq(ylim.val[1],ylim.val[2],yval[i-1]/2)
  x=barplot(SLE.SalEnv_count3[,i],col=adjustcolor(cols.alts,0.5),ylim=ylim.val,axes=F,ann=F,space=0)
  axis_fun(2,ymaj,ymin,ymaj)
  axis_fun(1,x,x,alts,cex=0.7,las=2)
  # abline(v=c(x[1]+(x[2]-x[1])/2,x[3]+(x[4]-x[3])/2),lwd=1)
  box(lwd=1)
  mtext(side=3,adj=0,SLE.labs[i-1],cex=0.5)
  text(x,SLE.SalEnv_count3[,i],
       round(SLE.SalEnv_count3[,i],0),font=2,col="black",pos=3,cex=0.75,offset=0.25)
  # if(i==2){mtext(side=3,adj=0,line=-1.25," CRE")}
  if(i==2){mtext(side=2,line=2.5,outer=F,"Count of 14-Day Periods")}
}
mtext(side=4,line=0.5,"2020 RECOVER PM")
mtext(side=1,line=1.75,outer=T,"Alternative")
dev.off()


# CRE MFL -----------------------------------------------------------------
q.dat.xtab$mfl.exceed=with(q.dat.xtab,ifelse(is.na(S79.30d)==T,0,ifelse(S79.30d<457,1,0)))
CRE.mfl.rslt=data.frame()
q.dat1.xtab.mfl=data.frame()
for(j in 1:n.alts){
  
  tmp=subset(q.dat.xtab,Alt==alts[j])
  ## Adapted from mflst_cre_v2.py
  for(i in 2:nrow(tmp)){
    if(tmp$mfl.exceed[i-1]==1&tmp$mfl.exceed[i]==0){
      tmp$exceed_end[i-1]=1 #found the last exceedance dates 
    }else{
      tmp$exceed_end[i-1]=0
    }
  }
  
  # subset(tmp,exceed_end==1)
  
  tmp$countdown=0
  tmp$exceed2=NA
  counts=0
  exc_n=0
  for(i in 30:nrow(tmp)){
    # rest counts
    if(tmp$mfl.exceed[i-1]==0&tmp$mfl.exceed[i]==1){
      counts=1
    }
    
    if(tmp$exceed_end[i]==1){
      if(tmp$countdown[i-1]<1){
        tmp$countdown[i]=365
      }else{
        tmp$countdown[i]=tmp$countdown[i-1]-1
        if(tmp$countdown[i]==0 & tmp$mfl.exceed[i]==1){
          tmp$countdown[i]=365
        }
      }
    }else{
      tmp$countdown[i]=tmp$countdown[i-1]-1
      counts=counts+1
      
      if(counts>366 & tmp$mfl.exceed[i]==1){
        tmp$countdown[i]=365
        counts=0
      }
      if(tmp$countdown[i]==0 & tmp$mfl.exceed[i]==1){
        tmp$countdown[i]=365
      }
    }
    
    #identify yearly violations
    if(tmp$countdown[i]<0){
      if(tmp$mfl.exceed[i]==1){
        if(tmp$mfl.exceed[i-1]!=1){
          tmp$exceed2[i]=1
          exc_n=exc_n+1}else{
            tmp$exceed2[i]=0
          }
      }else{tmp$exceed2[i]=0}
    }else{
      if(tmp$countdown[i]==365 & tmp$exceed_end[i]==0){
        tmp$exceed2[i]==1
        exc_n=exc_n+1
      }else{
        tmp$exceed2[i]=0
      }
    }
  }
  
  
  counts
  exc_n
  
  CRE.mfl.rslt=rbind(CRE.mfl.rslt,data.frame(Alt=alts[j],N.exceed=exc_n))
  q.dat1.xtab.mfl=rbind(q.dat1.xtab.mfl,tmp)
  print(j)
}

q.dat1.xtab.mfl$exc=with(q.dat1.xtab.mfl,ifelse(countdown<0&mfl.exceed==1,1,0))
q.dat1.xtab.mfl$plot_exc=with(q.dat1.xtab.mfl,ifelse(countdown<0&mfl.exceed==1,S79.30d,NA))
q.dat1.xtab.mfl$plot_exc365=with(q.dat1.xtab.mfl,ifelse(countdown>0&mfl.exceed==1,S79.30d,NA))

plot(exc~Date,subset(q.dat1.xtab.mfl,Alt=="FWO"))
summary(m1 <- glm(exc~Alt, family="poisson", data=q.dat1.xtab.mfl))


for(i in 1:n.alts){
# png(filename=paste0(plot.path,"CRE_MFL_Alt_",alts[i],".png"),width=6.5,height=5,units="in",res=200,type="windows",bg="white")
layout(matrix(1:2,2,1),heights=c(1,0.2))
par(family="serif",mar=c(1,3,0.75,1),oma=c(1.5,1,2.5,0.5));

ylim.val=c(0,15000);by.y=5000;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=date.fun(c("1965-01-01","2006-01-01"));xmaj=seq(xlim.val[1],xlim.val[2],"10 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")

plot(S79.30d~Date,q.dat1.xtab.mfl,type="n",yaxs="i",xlim=xlim.val,ylim=ylim.val,axes=F,ann=F)
abline(h=ymaj,v=xmaj,lty=1,col=adjustcolor("grey",0.5))
abline(h=457,col="brown",lty=2)
with(subset(q.dat1.xtab.mfl,Alt==alts[i]),lines(Date,S79.30d,col="blue",lty=1.5))
with(subset(q.dat1.xtab.mfl,Alt==alts[i]),lines(Date,plot_exc,col="orange",lty=1.5))
with(subset(q.dat1.xtab.mfl,Alt==alts[i]),lines(Date,plot_exc365,col="grey",lty=1.5))
axis_fun(1,xmaj,xmin,format(xmaj,"%Y"),line=-0.5)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=3,"Discharge (cfs)")
mtext(side=1,line=1.5,"Year")
mtext(side=3,paste0("MFL Recovery Water Body - Caloosahatchee River 30 Day Averege Flow at S79\n",
                    alts[i]," : ",subset(CRE.mfl.rslt,Alt==alts[i])$N.exceed,
                    " exceedance in 41 years of simualtion"))
plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA)
legend(0.5,-0.5,legend=c("30-day Moving Average","MFL Criteria (457 cfs)","Exceedance","Exceedance w/in 365 Days"),
       pch=NA,
       lty=c(1,2,1,1),lwd=2,
       col=c("blue","brown","orange","grey"),
       pt.bg=NA,
       pt.cex=1.5,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()
}

# png(filename=paste0(plot.path,"CRE_MFL_Alt_all.png"),width=6.5,height=7,units="in",res=200,type="windows",bg="white")
layout(matrix(1:5,5,1),heights=c(1,1,1,1,0.4))
par(family="serif",mar=c(0.75,2,0.5,1),oma=c(1.5,3,0.75,0.5));

# ylim.val=c(0,15000);by.y=5000;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
ylim.val=c(1,15000);ymaj=log.scale.fun(ylim.val,"major");ymin=log.scale.fun(ylim.val,"minor")
xlim.val=date.fun(c("1965-01-01","2006-01-01"));xmaj=seq(xlim.val[1],xlim.val[2],"10 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
for(i in 1:n.alts){
plot(S79.30d~Date,q.dat1.xtab.mfl,type="n",yaxs="i",xlim=xlim.val,ylim=ylim.val,axes=F,ann=F,log='y')
abline(h=ymaj,v=xmaj,lty=1,col=adjustcolor("grey",0.5))
abline(h=457,col="brown",lty=2)
with(subset(q.dat1.xtab.mfl,Alt==alts[i]),lines(Date,S79.30d,col="blue",lty=1.5))
with(subset(q.dat1.xtab.mfl,Alt==alts[i]),lines(Date,plot_exc,col="orange",lty=1.5))
with(subset(q.dat1.xtab.mfl,Alt==alts[i]),lines(Date,plot_exc365,col="grey",lty=1.5))
if(i==4){axis_fun(1,xmaj,xmin,format(xmaj,"%Y"),line=-0.5)}else{axis_fun(1,xmaj,xmin,NA,line=-0.5)}
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=1,adj=0,cex=0.7,line=-1.5,
      paste0("  Alternative: ",alts[i],";\n  ",subset(CRE.mfl.rslt,Alt==alts[i])$N.exceed,' MFL Exceedances'))
if(i==4){mtext(side=1,line=1.5,"Year")}
if(i==1){mtext(side=3,adj=0,"Caloosahatchee River Estuary")}
}
mtext(side=2,line=1.5,"Discharge (cfs)",outer=T)

par(mar=c(1,2,1,1))
plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA)
mtext(side=1,adj=1,"CY 1965 - 2005",cex=0.8)
legend(0.5,0,legend=c("30-day Moving Average","MFL Criteria (457 cfs)","Exceedance","Exceedance w/in 365 Days"),
       pch=NA,
       lty=c(1,2,1,1),lwd=2,
       col=c("blue","brown","orange","grey"),
       pt.bg=NA,
       pt.cex=1.5,ncol=2,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()


CRE.mfl.rslt=CRE.mfl.rslt[match(alts,CRE.mfl.rslt$Alt),]
CRE.mfl.rslt$NA25_perchange=with(CRE.mfl.rslt,round(((N.exceed-N.exceed[2])/N.exceed[2])*100,2))
CRE.mfl.rslt$zscore=with(CRE.mfl.rslt,N.exceed-mean(N.exceed)/sd(N.exceed))



# png(filename=paste0(plot.path,"CREMFL_Alts.png"),width=6.5,height=3.5,units="in",res=200,type="windows",bg="white")
ylim.val=c(0,40);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
par(family="serif",mar=c(2,2,0.5,0.25),oma=c(1,2,0.75,1),lwd=0.5);

x=barplot(CRE.mfl.rslt$N.exceed,beside=F,ylim=ylim.val,col=NA,border=NA,axes=F,ann=F,names.arg=rep(NA,length(alts)))
abline(h=ymaj,lty=1,col=adjustcolor("grey",0.5))
barplot(CRE.mfl.rslt$N.exceed,beside=F,ylim=ylim.val,col=adjustcolor(cols.alts,0.5),axes=F,ann=F,names.arg=rep(NA,length(alts)),add=T)
text(x,CRE.mfl.rslt$N.exceed,CRE.mfl.rslt$N.exceed,pos=3,offset=0.25)
text(x,CRE.mfl.rslt$N.exceed,paste0(CRE.mfl.rslt$NA25_perchange,"%"),pos=1,font=3)
axis_fun(1,x,x,alts,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2.5,"MFL Exceedances")
mtext(side=1,line=1.75,"Alternative")
mtext(side=3,adj=0,"Caloosahatchee MFL")
mtext(side=3,adj=1,"CY 1965 - 2005")
dev.off()

density.group.all=data.frame()
for(j in 1:length(alts)){
  
    tmp=subset(q.dat.xtab,Alt==alts[j])
    
    rng.val=range(log(tmp[tmp$S79!=0,]$S79),na.rm=T)
    nx=length(tmp$S79)
    w=rep(1/nx,nx)
    n.den=512#nx; #default is 512
    dens=density(log(tmp[tmp$S79!=0,]$S79),
                 # weights=w,   
                 kernel="gaussian",
                 from=rng.val[1],
                 to=rng.val[2],
                 n = n.den)
    dens.dat=data.frame(x=exp(dens$x),
                    y=dens$y,
                    scaled =  dens$y / max(dens$y, na.rm = TRUE),
                    ndensity = dens$y / max(dens$y, na.rm = TRUE),
                    count=dens$y*nx,
                    n=nx,
                    group="S79")
  
    rng.val=range(log(tmp[tmp$SLE.S80trib!=0,]$SLE.S80trib),na.rm=T)
    nx=length(tmp$SLE.S80trib)
    w=rep(1/nx,nx)
    n.den=512#nx; #default is 512
    dens=density(log(tmp[tmp$SLE.S80trib!=0,]$SLE.S80trib),
                 # weights=w,   
                 kernel="gaussian",
                 from=rng.val[1],
                 to=rng.val[2],
                 n = n.den)
    dens.dat=rbind(dens.dat,data.frame(x=exp(dens$x),
                        y=dens$y,
                        scaled =  dens$y / max(dens$y, na.rm = TRUE),
                        ndensity = dens$y / max(dens$y, na.rm = TRUE),
                        count=dens$y*nx,
                        n=nx,
                        group="SLE.S80trib")
                   )
    dens.dat$Alt=alts[j]
    density.group.all=rbind(density.group.all,dens.dat)
    print(j)
    }
plot(y~x,density.group.all,log="x")

# png(filename=paste0(plot.path,"Estuary_DaQDesnity.png"),width=6.5,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.5,0.5),oma=c(2,1,1,0.5));
layout(matrix(c(1:8),4,2,byrow=F))

ylim.val=c(0,2);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(10,20000);xmaj=log.scale.fun(xlim.val,"major");xmin=log.scale.fun(xlim.val,"minor")

for(i in 1:length(alts)){
plot(y~x,density.group.all,log="x",xlim=xlim.val,ylim=ylim.val,yaxs="i",type="n",axes=F,ann=F)
abline(h=ymaj,v=xmaj,lty=1,col=adjustcolor("grey",0.5))
with(subset(density.group.all,Alt==alts[i]&group=="S79"),
       shaded.range(x,rep(0.001,length(x)),y,bg=cols.alts[i],lty=1,col.adj = 0.5))
axis_fun(2,ymaj,ymin,format(ymaj))
if(i==1){mtext(side=3,adj=0,"Caloosahatchee River (S79)",cex=0.8)}
if(i==4){axis_fun(1,xmaj,xmin,xmaj)}else{axis_fun(1,xmaj,xmin,NA)}
box(lwd=1)
mtext(side=3,adj=0,line=-1.25,paste0("  ",alts[i]),cex=0.75)
}

xlim.val=c(100,15000);xmaj=log.scale.fun(xlim.val,"major");xmin=log.scale.fun(xlim.val,"minor")
ylim.val=c(0,1);by.y=0.25;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
for(i in 1:length(alts)){
  plot(y~x,density.group.all,log="x",xlim=xlim.val,ylim=ylim.val,yaxs="i",type="n",axes=F,ann=F)
  abline(h=ymaj,v=xmaj,lty=1,col=adjustcolor("grey",0.5))
  with(subset(density.group.all,Alt==alts[i]&group=="SLE.S80trib"),
       shaded.range(x,rep(0.001,length(x)),y,bg=cols.alts[i],lty=1,col.adj = 0.5))
  axis_fun(2,ymaj,ymin,format(ymaj))
  if(i==1){mtext(side=3,adj=0,"St Lucie & Tributaries (S80 + North Fork)",cex=0.8)}
  if(i==4){axis_fun(1,xmaj,xmin,xmaj)}else{axis_fun(1,xmaj,xmin,NA)}
  box(lwd=1)
  mtext(side=3,adj=0,line=-1.25,paste0("  ",alts[i]),cex=0.75)
}
mtext(side=1,line=1,outer=T,"Daily Discharge (ft\u00B3 sec\u207B\u00B9)")
mtext(side=2,line=-0.5,outer=T,"Density")
dev.off()
# Lake Stage --------------------------------------------------------------

lakeO.stage=data.frame()
for(i in 1:n.alts){
  dss_out=opendss(paste0(data.path,"RSMBN/",alts[i],"/RSMBN_output.dss"))  
  paths=paste0("/RSMBN/LOK/STAGE/01JAN1965 - 01JAN2005/1DAY/SIMULATED/")
  tmp=data.frame(getFullTSC(dss_out,paths))
  tmp$Date=date.fun(date.fun(rownames(tmp))-lubridate::ddays(1))
  rownames(tmp)<-NULL
  tmp$Alt=alts[i]
  lakeO.stage=rbind(tmp,lakeO.stage)
  print(i)
}
unique(lakeO.stage$Alt)
range(lakeO.stage$Date)

lakeO.stage$recess_7day=with(lakeO.stage,ave(STAGE,Alt,FUN=function(x) c(rep(NA,6),diff(x,lag=6))))

lakeO.stage$month=as.numeric(format(lakeO.stage$Date,'%m'))
lakeO.stage$day=as.numeric(format(lakeO.stage$Date,'%d'))
lakeO.stage$CY=as.numeric(format(lakeO.stage$Date,'%Y'))
lakeO.stage$DoY=as.numeric(format(lakeO.stage$Date,'%j'))
lakeO.stage$WY=WY(lakeO.stage$Date)
lakeO.stage$low.stg=with(lakeO.stage,ifelse(STAGE<11,1,0))
lakeO.stage$vlow.stg=with(lakeO.stage,ifelse(STAGE<=10,1,0))
lakeO.stage$High.stg=with(lakeO.stage,ifelse(STAGE>16,1,0))
lakeO.stage$vHigh.stg=with(lakeO.stage,ifelse(STAGE>=17,1,0))

ann.peak=ddply(lakeO.stage,c("CY","Alt"),summarise,max.stg=max(STAGE,na.rm=T))
ann.peak$GT17=with(ann.peak,ifelse(round(max.stg,1)>=16.9,1,0))

ann.peak=merge(ann.peak,data.frame(Alt=alts,plot.y=1:4),"Alt")

# png(filename=paste0(plot.path,"LOK_GT17_timeline.png"),width=7,height=2.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,1.25,0.5),oma=c(2,2,1,1),lwd=0.1);

xlim.val=c(1965,2005);by.x=10;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
ylim.val=c(0.75,4.25);by.y=1;ymaj=1:5#seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(plot.y~CY,ann.peak,type="n",xlim=xlim.val,ylim=ylim.val,ann=F,axes=F)
abline(h=ymaj,v=xmaj,lty=1,col=adjustcolor("grey",0.5))
for(i in 1:5){
  points(plot.y~CY,subset(ann.peak,GT17==1&Alt==alts[i]),pch=21,bg=adjustcolor(cols.alts[i],0.5),cex=1.5)
}
axis_fun(2,1:4,1:4,alts)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);box(lwd=1)
mtext(side=2,line=4,"Alternative")
mtext(side=1,line=1.5,"Calendar Year")
mtext(side=3,adj=0,line=1,"Lake Okeechobee")
mtext(side=3,adj=0,cex=0.75,"Annual maximum stage \u2265 17Ft NGVD29",col="grey50")
dev.off()


xlim.val=c(0,1);by.x=0.5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],0.1)
ylim.val=c(8,18);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
# png(filename=paste0(plot.path,"LOK_stage_curve.png"),width=6.5,height=3.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,1.5,1,0.5),oma=c(2,2.5,0.5,0.5),lwd=0.5);
layout(matrix(1:2,1,2,byrow=F))

ECB.sd=ecdf_fun(subset(lakeO.stage,Alt==alts[1])$STAGE)
ECB.sd$proportion=1-ECB.sd$proportion
FWO.sd=ecdf_fun(subset(lakeO.stage,Alt==alts[2])$STAGE)
FWO.sd$proportion=1-FWO.sd$proportion

plot(value~proportion,ECB.sd,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(ECB.sd,lines(proportion,value,col=adjustcolor(cols.alts[1],0.5),lwd=1))
with(FWO.sd,lines(proportion,value,col=adjustcolor(cols.alts[1],0.5),lwd=1,lty=2))
tmp=ecdf_fun(subset(lakeO.stage,Alt==alts[3])$STAGE)
tmp$proportion=1-tmp$proportion
lines(value~proportion,tmp,col=adjustcolor(cols.alts[3],0.5),lwd=2)
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.5)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
legend("topright",legend=alts[1:3],
       lty=c(1,2,1),lwd=c(1,1,2),col=adjustcolor(c(cols.alts[1],cols.alts[1],cols.alts[3]),0.5),
       ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=1,xpd=NA,xjust=0,yjust=1)
mtext(side=2,line=2,"Stage (Ft, NGVD29)",cex=1)
mtext(side=3,adj=0,"Lake Okeechobee")

plot(value~proportion,ECB.sd,ylim=ylim.val,xlim=xlim.val,ann=F,axes=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(ECB.sd,lines(proportion,value,col=adjustcolor(cols.alts[1],0.5),lwd=1))
with(FWO.sd,lines(proportion,value,col=adjustcolor(cols.alts[1],0.5),lwd=1,lty=2))
tmp=ecdf_fun(subset(lakeO.stage,Alt==alts[4])$STAGE)
tmp$proportion=1-tmp$proportion
lines(value~proportion,tmp,col=adjustcolor(cols.alts[4],0.5),lwd=2)
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.5)
axis_fun(2,ymaj,ymin,NA);box(lwd=1)
legend("topright",legend=alts[c(1,2,4)],
       lty=c(1,2,1),lwd=c(1,1,2),col=adjustcolor(c(cols.alts[1],cols.alts[1],cols.alts[4]),0.5),
       ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=1,xpd=NA,xjust=0,yjust=1)
mtext(side=1,outer=T,line=0.5,"Proportion of Time \u2265 Stage Elevation",cex=1)
mtext(side=3,adj=1,"CY1965 - 2005")
dev.off()



days.POS=ddply(lakeO.stage,"Alt",summarise,N.days=N.obs(STAGE),
               sum.mod.high=sum(STAGE>15,na.rm=T),
               sum.low=sum(low.stg),
               sum.vlow=sum(vlow.stg),
               sum.High=sum(High.stg),
               sum.vHigh=sum(vHigh.stg))
days.POS=days.POS[match(alts,days.POS$Alt),]

labs=c("% LOK Stage >17 Ft NGVD29 (Exterme High)",
       "% LOK Stage >15 Ft NGVD29 (Moderate High)",
       "% LOK Stage <11 Ft NGVD29 (Moderate Low)",
       "% LOK Stage <10 Ft NGVD29 (Extreme Low)")
# png(filename=paste0(plot.path,"LOK_stg.png"),width=4,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,1,1),oma=c(5,2,2,1),lwd=0.5);
layout(matrix(1:4,4,1,byrow=F))

ylim.val=c(0,2);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=barplot((days.POS$sum.vHigh/14975)*100,col=adjustcolor(cols.alts,0.5),
          ylim=ylim.val,space=0,axes=F,ann=F)
axis_fun(2,ymaj,ymin,ymaj)
axis_fun(1,x,x,NA,cex=0.8,las=2)
box(lwd=1)
mtext(side=3,adj=0,line=1.25,"LOK")
mtext(side=3,adj=0,labs[1],cex=0.7)
text(x,(days.POS$sum.vHigh/14975)*100,
     round((days.POS$sum.vHigh/14975)*100,1),font=2,col="black",pos=3,cex=1,offset=0.25)

ylim.val=c(0,40);by.y=20;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=barplot((days.POS$sum.mod.high/14975)*100,col=adjustcolor(cols.alts,0.5),
          ylim=ylim.val,space=0,axes=F,ann=F)
axis_fun(2,ymaj,ymin,ymaj)
axis_fun(1,x,x,NA,cex=0.8,las=2)
box(lwd=1)
mtext(side=3,adj=0,labs[2],cex=0.7)
text(x,(days.POS$sum.mod.high/14975)*100,
     round((days.POS$sum.mod.high/14975)*100,1),font=2,col="black",pos=3,cex=1,offset=0.25)

ylim.val=c(0,14);by.y=7;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=barplot((days.POS$sum.low/14975)*100,col=adjustcolor(cols.alts,0.5),
          ylim=ylim.val,space=0,axes=F,ann=F)
axis_fun(2,ymaj,ymin,ymaj)
axis_fun(1,x,x,NA,cex=0.8,las=2)
box(lwd=1)
mtext(side=3,adj=0,labs[3],cex=0.7)
text(x,(days.POS$sum.low/14975)*100,
     round((days.POS$sum.low/14975)*100,1),font=2,col="black",pos=3,cex=1,offset=0.25)

ylim.val=c(0,6);by.y=3;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=barplot((days.POS$sum.vlow/14975)*100,col=adjustcolor(cols.alts,0.5),
          ylim=ylim.val,space=0,axes=F,ann=F)
axis_fun(2,ymaj,ymin,ymaj)
axis_fun(1,x,x,alts,cex=0.8,las=2)
box(lwd=1)
mtext(side=3,adj=0,labs[4],cex=0.7)
text(x,(days.POS$sum.vlow/14975)*100,
     format(round((days.POS$sum.vlow/14975)*100,1)),font=2,col="black",pos=3,cex=1,offset=0.25)
mtext(side=2,line=0.5,"Percent of Time",outer=T)
mtext(side=1,line=4,"Alternative",outer=F)
dev.off()


# Old RECOVER Stage Envelope ----------------------------------------------
# stg.env.vals=data.frame(
#   month=c(sort(rep(1:12,2)),12),
#   day=c(rep(c(1,15),12),31),
#   Lower=c(14.50,14.50,14.23,14.00,13.75,13.50,13.23,13.00,12.47,rep(12,7),12.27,12.50,13.03,13.50,14.05,rep(14.5,4)),
#   Upper=c(16.00,16.00,15.73,15.50,15.25,15.00,14.73,14.50,13.97,13.50,13.23,13.00,13.00,13.00,13.27,13.50,14.05,14.50,15.03,15.50,15.77,rep(16,4))
# )
# 
# plot(1:25,stg.env.vals$Lower,ylim=c(11,17),type="l")
# lines(1:25,stg.env.vals$Upper)

stg.env.vals=read.csv(paste0(wd,"/_resources/LOK_2007StgEnv_Da.csv"))
stg.env.vals$Date=date.fun(stg.env.vals$Date,form="%m/%d/%Y")
stg.env.vals.da=ddply(stg.env.vals,"Jday",summarise,Env_Low=mean(Env_Low),Env_High=mean(Env_High))

plot(Env_Low~Jday,stg.env.vals.da,ylim=c(11,17),type="l")
lines(Env_High~Jday,stg.env.vals.da)

stg.env.07=lakeO.stage
stg.env.07=merge(stg.env.07,stg.env.vals.da,by.x="DoY",by.y="Jday",all.x=T)
# stg.env.07=merge(stg.env.07,stg.env.vals[,c("Date","Env_Low","Env_High")],by.x="Date",by.y="Date",all.x=T)
stg.env.07=stg.env.07[order(stg.env.07$Alt,stg.env.07$Date),]

stg.env.07$dev_lower=with(stg.env.07,pmax(0,Env_Low-STAGE))
stg.env.07$dev_upper=with(stg.env.07,pmax(0,STAGE-Env_High))
stg.env.07$env=with(stg.env.07,ifelse(dev_upper+dev_lower>0,0,1))

stg.score.07=ddply(stg.env.07,"Alt",summarise,
                   N.env=(sum(env,na.rm=T)/14975)*100,
                   N.below.env=(sum(dev_lower>0)/14975)*100,
                   N.above.env=(sum(dev_upper>0)/14975)*100,
                   sum.dev_lower.ftdays=sum(dev_lower),
                   sum.dev_upper.ftdays=sum(dev_upper),
                   stg_env=(sum(env)/14975)*100,
                   ex.low=(sum(STAGE<10.0)/14975)*100,
                   ex.high=(sum(STAGE>17.0)/14975)*100)
stg.score.07$sum.dev_lower.ftwks=stg.score.07$sum.dev_lower/7
stg.score.07$sum.dev_upper.ftwks=stg.score.07$sum.dev_upper/7


stg.score.07$std.score.exlow=with(stg.score.07,((ex.low/100)*36*52)*-0.185+100)
stg.score.07$std.score.exhigh=with(stg.score.07,((ex.high/100)*36*52)*-0.253+100)
# 41-yr Equation
stg.score.07$std.score.lower=with(stg.score.07,ifelse(sum.dev_lower.ftwks<218.7,100,sum.dev_lower.ftwks*-0.05227+111.429))
stg.score.07$std.score.higher=with(stg.score.07,sum.dev_upper.ftwks*-0.0469+100)

# 36-yr RECOVER Equation
# stg.score.07$std.score.lower.36y=with(stg.score.07,ifelse(sum.dev_lower.ftwks<192,100,sum.dev_lower.ftwks*-0.0595+111.429))
# stg.score.07$std.score.higher.36y=with(stg.score.07,sum.dev_upper.ftwks*-0.0534+100)

stg.score.07$wgh.score.exlow=with(stg.score.07,(0.25*std.score.exlow))
stg.score.07$wgh.score.exhigh=with(stg.score.07,(0.50*std.score.exhigh))
stg.score.07$wgh.score.std.lower=with(stg.score.07,(0.10*std.score.lower))
stg.score.07$wgh.score.std.higher=with(stg.score.07,(0.15*std.score.higher))
stg.score.07$wgh.scr=rowSums(stg.score.07[,c("wgh.score.exlow",
                                             "wgh.score.exhigh",
                                             "wgh.score.std.lower",
                                             "wgh.score.std.higher")])
stg.score.07=stg.score.07[match(alts,stg.score.07$Alt),]


# png(filename=paste0(plot.path,"LOK_OldScore_total_parsed.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
ylim.val=c(0,105);by.y=20;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
par(family="serif",mar=c(1,2,0.5,1),oma=c(2.5,2,0.5,0.5),lwd=0.5);
layout(matrix(1:2,1,2),widths=c(1,0.5))

cols.stack=viridis::magma(4,0.5)
x=barplot(t(stg.score.07[,c("wgh.score.exlow",
                          "wgh.score.exhigh",
                          "wgh.score.std.lower",
                          "wgh.score.std.higher")]),
        col=cols.stack,ylim=ylim.val,
        border=adjustcolor("grey30",0.5),names.arg=rep(NA,n.alts),
        space=0.05,axes=F,ann=F)
text(x,stg.score.07$wgh.scr,round(stg.score.07$wgh.scr,2),pos=3,offset=0.2)
axis_fun(1,x,x,alts,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"Standardized Score")
mtext(side=1,line=2,"Alternative")

plot(0:1,0:1,ann=F,axes=F,type="n")
mtext(side=3,line=-2,"Lake Okeechobee\nEnvelope Penalty Score\nAll Years")
legend("center",legend=rev(c("Extreme Low", "Extreme High","Below Envelope","Above Envelope")),
       pch=22,lty=0,lwd=0.1,pt.bg=rev(cols.stack),pt.cex=2,col=adjustcolor("grey30",0.5),
       ncol=1,cex=0.8,bty="n",y.intersp=1,x.intersp=1,xpd=NA,xjust=0,yjust=1,
       title.adj = 0,title=" Performance Measure")
mtext(side=1,line=-1,adj=1,"Simulation Period of\nRecord (1965-2005). \n2007 RECOVER LOK PM",cex=0.75)
dev.off()


# png(filename=paste0(plot.path,"LOK_OldScore_total.png"),width=4,height=3,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1.5,2,0.5,0.25),oma=c(2.5,2,0.5,1),lwd=0.5);

ylim.val=c(0,100);by.y=50;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)

x=barplot(stg.score.07$wgh.scr,ylim=ylim.val,col=NA,border=NA,axes=F,ann=F,names.arg=NA)
abline(h=ymaj,lty=1,col=adjustcolor("grey",0.5))
barplot(stg.score.07$wgh.scr,beside=F,ylim=ylim.val,col=adjustcolor(cols.alts,0.5),axes=F,ann=F,names.arg=NA,add=T)
with(stg.score.07,text(x,wgh.scr,round(wgh.scr,1),pos=3,offset=0.2,cex=0.8))
axis_fun(1,x,x,alts)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,adj=0,"2007 RECOVER PM (1965 - 2005)",cex=0.75)
mtext(side=1,line=2.5,"Alternative")
mtext(side=2,outer=F,line=2,"Weighted Score (%)")
dev.off()

# png(filename=paste0(plot.path,"LOK_OldScore.png"),width=4,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1.5,2,0.5,0.25),oma=c(2.5,3,0.5,1),lwd=0.5);
layout(matrix(1:4,4,1,byrow=T),widths=c(1,0.5))
ylim.val=c(0,120);by.y=50;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)

x=barplot(stg.score.07$std.score.exlow,ylim=ylim.val,col=NA,border=NA,axes=F,ann=F,names.arg=NA)
abline(h=ymaj,lty=1,col=adjustcolor("grey",0.5))
barplot(stg.score.07$std.score.exlow,beside=F,ylim=ylim.val,col="forestgreen",axes=F,ann=F,names.arg=NA,add=T)
with(stg.score.07,text(x,std.score.exlow,round(std.score.exlow,2),pos=3,offset=0.2,cex=0.8))
axis_fun(1,x,x,NA)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,adj=0,"Stage Below 10 Feet NGVD (1965 - 2005)",cex=0.75)

x=barplot(stg.score.07$std.score.exhigh,ylim=ylim.val,col=NA,border=NA,axes=F,ann=F,names.arg=NA)
abline(h=ymaj,lty=1,col=adjustcolor("grey",0.5))
barplot(stg.score.07$std.score.exhigh,beside=F,ylim=ylim.val,col="forestgreen",axes=F,ann=F,names.arg=NA,add=T)
with(stg.score.07,text(x,std.score.exhigh,format(round(std.score.exhigh,2)),pos=3,offset=0.2,cex=0.8))
axis_fun(1,x,x,NA)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,adj=0,"Stage Above 17 Feet NGVD (1965 - 2005)",cex=0.75)

x=barplot(stg.score.07$std.score.higher,ylim=ylim.val,col=NA,border=NA,axes=F,ann=F,names.arg=NA)
abline(h=ymaj,lty=1,col=adjustcolor("grey",0.5))
barplot(stg.score.07$std.score.higher,beside=F,ylim=ylim.val,col="forestgreen",axes=F,ann=F,names.arg=NA,add=T)
with(stg.score.07,text(x,std.score.higher,format(round(std.score.higher,2)),pos=3,offset=0.2,cex=0.8))
axis_fun(1,x,x,NA)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,adj=0,"Score Above Envelope - Weekly Calc (1965 - 2005)",cex=0.75)

x=barplot(stg.score.07$std.score.lower,ylim=ylim.val,col=NA,border=NA,axes=F,ann=F,names.arg=NA)
abline(h=ymaj,lty=1,col=adjustcolor("grey",0.5))
barplot(stg.score.07$std.score.lower,beside=F,ylim=ylim.val,col="forestgreen",axes=F,ann=F,names.arg=NA,add=T)
with(stg.score.07,text(x,std.score.lower,format(round(std.score.lower,2)),pos=3,offset=0.2,cex=0.8))
axis_fun(1,x,x,alts)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=3,adj=0,"Score Below Envelope - Weekly Calc (1965 - 2005)",cex=0.75)
mtext(side=1,line=2.5,"Alternative")
mtext(side=2,outer=T,line=0.5,"Standard Score (%)")
dev.off()



stg.env.07.eom=data.frame()
for(i in 1:n.alts){
  tmp=subset(stg.env.07,Alt==alts[i])
  
  tmp$eom.dev_lower=with(tmp,ifelse(as.numeric(format(date.fun(Date+lubridate::ddays(1)),"%d"))==1,dev_lower,0))
  tmp$eom.dev_upper=with(tmp,ifelse(as.numeric(format(date.fun(Date+lubridate::ddays(1)),"%d"))==1,dev_upper,0))
  
  stg.env.07.eom=rbind(tmp,stg.env.07.eom)
  print(i)
}

stg.score.07.eom=ddply(stg.env.07.eom,"Alt",summarise,
                       sum.dev_lower=sum(eom.dev_lower), #foot-months
                       sum.dev_upper=sum(eom.dev_upper)) #foot-months
stg.score.07.eom$sum.dev_lower.ftwks=with(stg.score.07.eom,sum.dev_lower/12*52)
stg.score.07.eom$sum.dev_upper.ftwks=with(stg.score.07.eom,sum.dev_upper/12*52)
stg.score.07.eom$std.score.lower=with(stg.score.07.eom,ifelse(sum.dev_lower.ftwks<218.7,100,sum.dev_lower.ftwks*-0.05227+111.429))
stg.score.07.eom$std.score.upper=with(stg.score.07.eom,sum.dev_upper.ftwks*-0.0469+100)

# NEW RECOVER Stage Envelope ----------------------------------------------
## 
AprSep=seq(date.fun("1965-04-15"),date.fun("1965-09-15"),"1 days")
MayAug=seq(date.fun("1965-05-01"),date.fun("1965-09-01"),"1 days")

stg.env=lakeO.stage
# stg.env$STAGE=round(stg.env$STAGE,2)
stg.env$CY=as.numeric(format(stg.env$Date,"%Y"))
stg.env$month_day=as.character(format(stg.env$Date,"%m_%d"))
stg.env$AugDec.stg=with(stg.env,ifelse(as.numeric(format(Date,"%m"))%in%c(8:12),STAGE,NA))
stg.env$JuneJuly.stg=with(stg.env,ifelse(as.numeric(format(Date,"%m"))%in%c(6,7),STAGE,NA))
stg.env$AprSep.stg=with(stg.env,ifelse(month_day%in%as.character(format(AprSep,"%m_%d")),STAGE,NA))
stg.env$MayAug.stg=with(stg.env,ifelse(month_day%in%as.character(format(MayAug,"%m_%d")),STAGE,NA))

stg.env.CY=ddply(stg.env,c("Alt","CY"),summarise,
                 max.stg=max(STAGE,na.rm=T),
                 AugDec.max=max(AugDec.stg,na.rm=T),
                 JuneJuly_13=sum(JuneJuly.stg>13,na.rm=T),
                 MayAug_11.5=sum(MayAug.stg<11.5,na.rm=T),
                 AprSep_12=sum(AprSep.stg<12,na.rm=T))

# 1 = Normal
# 2 = Recovery
stg.env.CY$env=1
env.rslt=data.frame()
for(i in 1:length(alts)){
  tmp=subset(stg.env.CY,Alt==alts[i])
  for(j in 2:nrow(tmp)){
    tmp$env[j]=with(tmp,
                    ifelse(max.stg[j-1]>17|JuneJuly_13[j-1]>=30,2,
                           ifelse(AugDec.max[j-1]<=16&(MayAug_11.5[j-1]>=60|AprSep_12[j-1]>=90),1,tmp$env[j-1])))
    
  }
  env.rslt=rbind(env.rslt,tmp)
}

env.rslt$env2=with(env.rslt,ifelse(env==2,0,1))
env.rslt=merge(env.rslt,data.frame(Alt=alts,PlotOffset=rev(seq(2,8,2))),"Alt")
env.rslt$env.plt=with(env.rslt,env2+PlotOffset)
env.rslt$env.f=with(env.rslt,ifelse(env==1,"normal","recovery"))
env.rslt$Alt=factor(env.rslt$Alt,levels=alts)

#env.count=reshape2::dcast(env.rslt,Alt~env.f,value.var = "env",function(x) N.obs(x))
env.count=ddply(env.rslt,c("Alt","env.f"),summarise,N.val=N.obs(env.f))
env.count$axs.val=9:2

# png(filename=paste0(plot.path,"LOK_NewEnv.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,1.75,1),oma=c(2,3,1,2));

ylim.val=c(1.5,9.5)
xlim.val=c(1965,2005);by.x=4;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
plot(env2~CY,env.rslt,type="n",axes=F,ann=F,xlim=xlim.val,xaxs="i",yaxs="i",ylim=ylim.val)
abline(v=c(xmaj,xmin),h=c(2:19),lty=c(1,3),lwd=0.5,col=adjustcolor("grey",0.5))
abline(h=c(2:19),lty=c(3),lwd=0.5,col=c("grey"))
for(i in 1:length(alts)){
  with(subset(env.rslt,Alt==alts[i]),lines(env.plt~CY,type="s",col=cols.alts[i],lwd=3))
}
axis_fun(1,xmaj,xmin,xmaj,line=-0.75,cex=0.75)
axis_fun(2,seq(2.5,8.5,2),seq(2.5,8.5,2),rev(alts))
axis_fun(4,env.count$axs.val,env.count$axs.val,env.count$N.val,cex=0.6)
abline(h=seq(3.5,9.5,2))
box(lwd=1)
mtext(side=3,adj=1,"Upper Step = Normal Envelope\nLower Step = Recovery Envelope",font=3)
mtext(side=1,line=1.25,"Calendar Year")
mtext(side=4,line=1.25,"Normal/Recovery Envelop Count")
dev.off()


library(LORECOVER)
lakeO.stage$Data.Value=lakeO.stage$STAGE

norm.lakeO.stage.scr=data.frame()
for(i in 1:n.alts){
  tmp=subset(lakeO.stage,Alt==alts[i])[,c("Date","Data.Value")]
  rslt=LORECOVER::norm_env(tmp)
  rslt$Alt=alts[i]
  norm.lakeO.stage.scr=rbind(norm.lakeO.stage.scr,rslt)
  print(i)
}

norm.lakeO.stage.scr=rename(norm.lakeO.stage.scr,c("penalty"="norm.score"))
norm.lakeO.stage.scr$WY=WY(norm.lakeO.stage.scr$Date)

rec.lakeO.stage.scr=data.frame()
for(i in 1:n.alts){
  tmp=subset(lakeO.stage,Alt==alts[i])[,c("Date","Data.Value")]
  rslt=LORECOVER::rec_env(tmp)
  rslt$Alt=alts[i]
  rec.lakeO.stage.scr=rbind(rec.lakeO.stage.scr,rslt)
  print(i)
}
rec.lakeO.stage.scr=rename(rec.lakeO.stage.scr,c("penalty"="rec.score"))
rec.lakeO.stage.scr$WY=WY(rec.lakeO.stage.scr$Date)

lakeO.stage.scr=merge(norm.lakeO.stage.scr,rec.lakeO.stage.scr[,c("Date","Alt","rec.score")],c("Date","Alt"))
lakeO.stage.scr$CY=as.numeric(format(lakeO.stage.scr$Date,"%Y"))
lakeO.stage.scr$WY=WY(lakeO.stage.scr$Date)
lakeO.stage.scr=lakeO.stage.scr[order(lakeO.stage.scr$Alt,lakeO.stage.scr$Date),]

vars=c("Alt","CY","env")
lakeO.stage.scr=merge(lakeO.stage.scr,env.rslt[,vars],c("Alt","CY"))
lakeO.stage.scr$score=with(lakeO.stage.scr,ifelse(env==1,norm.score,rec.score))
lakeO.stage.scr$Alt_CY=with(lakeO.stage.scr,paste(Alt,CY,sep="_"))
lakeO.stage.scr$cum.abs.pen=with(lakeO.stage.scr,ave(score,Alt_CY,FUN=function(x)cumsum(abs(x))))
lakeO.stage.scr$Alt=factor(lakeO.stage.scr$Alt,levels=alts)

boxplot(cum.abs.pen~Alt,lakeO.stage.scr,outline=F)
stg.scr.sum=ddply(lakeO.stage.scr,"Alt",summarise,mean.val=mean(cum.abs.pen,na.rm=T))
stg.scr.sum$FWO.perdiff=with(stg.scr.sum,((mean.val-mean.val[1])/mean.val[1])*100)

env.pen.sum=ddply(lakeO.stage.scr,"Alt",summarise,
                  N.val=N.obs(score),
                  pen_above=sum(score[score>0],na.rm=T),
                  pen_below=sum(abs(score)[score<0],na.rm=T),
                  per_below=(sum(score<0)/N.obs(score))*100,
                  per0=(sum(score==0,na.rm=T)/N.obs(score))*100,
                  per_above=(sum(score>0)/N.obs(score))*100)
env.pen.sum$total.pen=rowSums(env.pen.sum[,c("pen_above","pen_below")])
env.pen.sum=env.pen.sum[match(alts,env.pen.sum$Alt),]
env.pen.sum$FWO_PerBelow=with(env.pen.sum,(per_below-per_below[1])/per_below[1])*100
env.pen.sum$FWO_PerWith=with(env.pen.sum,(per0-per0[1])/per0[1])*100
env.pen.sum$FWO_PerAbove=with(env.pen.sum,(per_above-per_above[1])/per_above[1])*100
env.pen.sum=env.pen.sum[match(alts,env.pen.sum$Alt),]
# env.pen.sum$FWO_totalpen=with(env.pen.sum,(total.pen-total.pen[1])/total.pen[1])*100

env.pen.sum.melt=reshape2::melt(env.pen.sum[,c("Alt","pen_above","pen_below","total.pen")],id.var="Alt")
env.pen.sum.FWO=reshape2::dcast(env.pen.sum.melt,variable~Alt,value.var = "value")
env.pen.sum.FWO$ASR.PerFWO=with(env.pen.sum.FWO,((ASR-FWO)/FWO)*100)

lakeO.stage.scr$month=as.numeric(format(lakeO.stage.scr$Date,"%m"))
env.pen.sum.maysept=ddply(subset(lakeO.stage.scr,month%in%seq(5,9,1)),"Alt",summarise,
                          N.val=N.obs(score),
                          pen_above=sum(score[score>0],na.rm=T),
                          pen_below=sum(abs(score)[score<0],na.rm=T),
                          per_below=(sum(score<0)/N.obs(score))*100,
                          per0=(sum(score==0,na.rm=T)/N.obs(score))*100,
                          per_above=(sum(score>0)/N.obs(score))*100)
env.pen.sum.maysept$total.pen=rowSums(env.pen.sum.maysept[,c("pen_above","pen_below")])
env.pen.sum.maysept=env.pen.sum.maysept[match(alts,env.pen.sum.maysept$Alt),]


env.pen.sum.rec=ddply(subset(lakeO.stage.scr,env==2),"Alt",summarise,
                          N.val=N.obs(score),
                          pen_above=sum(score[score>0],na.rm=T),
                          pen_below=sum(abs(score)[score<0],na.rm=T),
                          per_below=(sum(score<0)/N.obs(score))*100,
                          per0=(sum(score==0,na.rm=T)/N.obs(score))*100,
                          per_above=(sum(score>0)/N.obs(score))*100)
env.pen.sum.rec$total.pen=rowSums(env.pen.sum.rec[,c("pen_above","pen_below")])
env.pen.sum.rec=env.pen.sum.rec[match(alts,env.pen.sum.rec$Alt),]

cols.IMC=c(rgb(238,232,170,maxColorValue = 255),rgb(143,188,143,maxColorValue = 255))
# png(filename=paste0(plot.path,"LOK_NewEnvScore_AllYrs.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
ylim.val=c(0,30000);by.y=10000;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
par(family="serif",mar=c(1,2,0.5,0.25),oma=c(2,2,0.75,1),lwd=0.5);
layout(matrix(1:2,1,2,byrow=T),widths=c(1,0.5))

x=barplot(t(env.pen.sum[,c("pen_below","pen_above")]),beside=F,ylim=ylim.val,col=NA,border=NA,axes=F,ann=F,names.arg=rep(NA,nrow(env.pen.sum)))
abline(h=ymaj,lty=1,col=adjustcolor("grey",0.5))
barplot(t(env.pen.sum[,c("pen_below","pen_above")]),beside=F,ylim=ylim.val,col=cols.IMC,axes=F,ann=F,names.arg=rep(NA,nrow(env.pen.sum)),add=T)
with(env.pen.sum,text(x,pen_below/2,round(pen_below,0)))
with(env.pen.sum,text(x,pen_below+abs((pen_below-total.pen)/2),round(pen_above,0)))
with(env.pen.sum,text(x,total.pen,round(total.pen,0),pos=3))
axis_fun(1,x,x,alts,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=3,"Stage Envelope Penalty Score")
mtext(side=1,line=1.75,"Alternative")

plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA)
mtext(side=3,line=-2,"Lake Okeechobee\nEnvelope Penalty Score\nAll Years")
legend(0.5,0.5,legend=c("Total Score","Upper Penalty\n(Above Envelope)","Lower Penalty\n(Below Envelope)"),
       pch=22,
       lty=0,lwd=0.01,
       col=c("white","black","black"),
       pt.bg=c("white",rev(cols.IMC)),
       pt.cex=2,ncol=1,cex=1,bty="n",y.intersp=1.25,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5,
       title.adj=0,title="Penalty Score")
mtext(side=1,line=-1,adj=1,"Simulation Period of Record\n(1965-2005). Includes Normal\nand Recovery envelopes.\n2020 RECOVER LOK PM",cex=0.75)
dev.off()

# png(filename=paste0(plot.path,"LOK_NewEnvScore_MaySep.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
ylim.val=c(0,20000);by.y=5000;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
par(family="serif",mar=c(1,2,0.5,0.25),oma=c(2,2,0.75,1),lwd=0.5);
layout(matrix(1:2,1,2,byrow=T),widths=c(1,0.5))

x=barplot(t(env.pen.sum.maysept[,c("pen_below","pen_above")]),beside=F,ylim=ylim.val,col=NA,border=NA,axes=F,ann=F,names.arg=rep(NA,nrow(env.pen.sum)))
abline(h=ymaj,lty=1,col=adjustcolor("grey",0.5))
barplot(t(env.pen.sum.maysept[,c("pen_below","pen_above")]),beside=F,ylim=ylim.val,col=cols.IMC,axes=F,ann=F,names.arg=rep(NA,nrow(env.pen.sum)),add=T)
with(env.pen.sum.maysept,text(x,pen_below/2,round(pen_below,0)))
with(env.pen.sum.maysept,text(x,pen_below+abs((pen_below-total.pen)/2),round(pen_above,0)))
with(env.pen.sum.maysept,text(x,total.pen,round(total.pen,0),pos=3))
axis_fun(1,x,x,alts,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=3,"Stage Envelope Penalty Score")
mtext(side=1,line=1.75,"Alternative")

plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA)
mtext(side=3,line=-3,"Lake Okeechobee\nEnvelope Penalty Score\nMay-Sep\nAll Years")
legend(0.5,0.5,legend=c("Lower Penalty\n(Below Envelope)","Upper Penalty\n(Above Envelope)"),
       pch=22,
       lty=0,lwd=0.01,
       col="black",
       pt.bg=c(cols.IMC),
       pt.cex=2,ncol=1,cex=1,bty="n",y.intersp=1.5,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5,
       title.adj=0,title="Total Score")
mtext(side=1,line=-1,adj=1,"Simulation Period of Record\n(1965-2005). Includes Normal\nand Recovery envelope.\n2021 RECOVER LOK PM",cex=0.75)
dev.off()

# png(filename=paste0(plot.path,"LOK_NewEnvScore_RecYrs.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
ylim.val=c(0,30000);by.y=10000;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
par(family="serif",mar=c(1,2,0.5,0.25),oma=c(2,2,0.75,1),lwd=0.5);
layout(matrix(1:2,1,2,byrow=T),widths=c(1,0.5))

x=barplot(t(env.pen.sum.rec[,c("pen_below","pen_above")]),beside=F,ylim=ylim.val,col=NA,border=NA,axes=F,ann=F,names.arg=rep(NA,nrow(env.pen.sum)))
abline(h=ymaj,lty=1,col=adjustcolor("grey",0.5))
barplot(t(env.pen.sum.rec[,c("pen_below","pen_above")]),beside=F,ylim=ylim.val,col=cols.IMC,axes=F,ann=F,names.arg=rep(NA,nrow(env.pen.sum)),add=T)
with(env.pen.sum.rec,text(x,pen_below/2,round(pen_below,0)))
with(env.pen.sum.rec,text(x,pen_below+abs((pen_below-total.pen)/2),round(pen_above,0)))
with(env.pen.sum.rec,text(x,total.pen,round(total.pen,0),pos=3))
axis_fun(1,x,x,alts,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=3,"Stage Envelope Penalty Score")
mtext(side=1,line=1.75,"Alternative")

plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA)
mtext(side=3,line=-2,"Lake Okeechobee\nEnvelope Penalty Score\nRecovery Years")
legend(0.5,0.5,legend=c("Total Score","Upper Penalty\n(Above Envelope)","Lower Penalty\n(Below Envelope)"),
       pch=22,
       lty=0,lwd=0.01,
       col=c("white","black","black"),
       pt.bg=c("white",rev(cols.IMC)),
       pt.cex=2,ncol=1,cex=1,bty="n",y.intersp=1.25,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5,
       title.adj=0,title="Penalty Score")
mtext(side=1,line=-1,adj=1,"Simulation Period of Record\n(1965-2005). Includes only years\nin the recovery envelopes.\n2020 RECOVER LOK PM",cex=0.75)
dev.off()




# png(filename=paste0(plot.path,"LOK_NewEnvScore.png"),width=6.5,height=5,units="in",res=200,type="windows",bg="white")
ylim.val=c(0,30000);by.y=10000;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
par(family="serif",mar=c(1.5,2,0.5,0.25),oma=c(2,2.5,0.75,1),lwd=0.5);
layout(matrix(c(1:2,3,3),2,2,byrow=F),widths=c(1,0.5))

x=barplot(t(env.pen.sum[,c("pen_below","pen_above")]),beside=F,ylim=ylim.val,col=NA,border=NA,axes=F,ann=F,names.arg=rep(NA,nrow(env.pen.sum)))
abline(h=ymaj,lty=1,col=adjustcolor("grey",0.5))
barplot(t(env.pen.sum[,c("pen_below","pen_above")]),beside=F,ylim=ylim.val,col=cols.IMC,axes=F,ann=F,names.arg=rep(NA,nrow(env.pen.sum)),add=T)
with(env.pen.sum,text(x,pen_below/2,round(pen_below,0)))
with(env.pen.sum,text(x,pen_below+abs((pen_below-total.pen)/2),round(pen_above,0)))
with(env.pen.sum,text(x,total.pen,round(total.pen,0),pos=3))
# axis_fun(1,x,x,alts,line=-0.5)
axis_fun(1,x,x,NA,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
text(x[4]+min(diff(x))*0.7,ylim.val[2]-ylim.val[2]*0.05,"A",xpd=NA,cex=1.5)
mtext(side=3,adj=0,"All Year")

x=barplot(t(env.pen.sum.rec[,c("pen_below","pen_above")]),beside=F,ylim=ylim.val,col=NA,border=NA,axes=F,ann=F,names.arg=rep(NA,nrow(env.pen.sum)))
abline(h=ymaj,lty=1,col=adjustcolor("grey",0.5))
barplot(t(env.pen.sum.rec[,c("pen_below","pen_above")]),beside=F,ylim=ylim.val,col=cols.IMC,axes=F,ann=F,names.arg=rep(NA,nrow(env.pen.sum)),add=T)
with(env.pen.sum.rec,text(x,pen_below/2,round(pen_below,0)))
with(env.pen.sum.rec,text(x,pen_below+abs((pen_below-total.pen)/2),round(pen_above,0)))
with(env.pen.sum.rec,text(x,total.pen,round(total.pen,0),pos=3))
axis_fun(1,x,x,alts,line=-0.5)
axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
text(x[4]+min(diff(x))*0.7,ylim.val[2]-ylim.val[2]*0.05,"B",xpd=NA,cex=1.5)
mtext(side=2,line=1.25,outer=T,"Stage Envelope Penalty Score")
mtext(side=1,line=1.75,"Alternative")
mtext(side=3,adj=0,"Recovery Years Only")

plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA)
legend(0.5,0.5,legend=c("Total Score","Upper Penalty\n(Above Envelope)","Lower Penalty\n(Below Envelope)"),
       pch=22,
       lty=0,lwd=0.01,
       col=c("white","black","black"),
       pt.bg=c("white",rev(cols.IMC)),
       pt.cex=2,ncol=1,cex=1,bty="n",y.intersp=1.25,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5,
       title.adj=0,title="Penalty Score")
mtext(side=1,line=-1,adj=1,"Simulation Period of Record\n(1965-2005).\n2020 RECOVER LOK PM",cex=0.75)
dev.off()
# Consecutive Events Stage ------------------------------------------------
highstg_consec=data.frame()
for(j in 1:length(alts)){
  tmp=subset(lakeO.stage,Alt==alts[j])
  tmp$stg17=0
  tmp$stg16=0
  for(i in 2:nrow(tmp)){
    # tmp$consec[i]=with(tmp,ifelse(Q_GT2100[i-1]>0&Q_GT2100[i]>0,1,0))
    tmp$stg16[i]=with(tmp,ifelse(High.stg[i-1]==0&High.stg[i]>0,1,
                                 ifelse(High.stg[i-1]>0&High.stg[i]>0,1,0)))
    tmp$stg17[i]=with(tmp,ifelse(vHigh.stg[i-1]==0&vHigh.stg[i]>0,1,
                                 ifelse(vHigh.stg[i-1]>0&vHigh.stg[i]>0,1,0)))
    
  }
  
  highstg=consec.startend(tmp$stg16>0)
  tmp$sum.stg16=0
  for(i in 1:length(highstg$ends)){
    tmp[highstg$ends[i],]$sum.stg16=with(tmp[c(highstg$starts[i]:highstg$ends[i]),],sum(stg16,na.rm=T))
  }
  vhighstg=consec.startend(tmp$stg17>0)
  tmp$sum.stg17=0
  for(i in 1:length(vhighstg$ends)){
    tmp[vhighstg$ends[i],]$sum.stg17=with(tmp[c(vhighstg$starts[i]:vhighstg$ends[i]),],sum(stg17,na.rm=T))
  }
  
  highstg_consec=rbind(tmp,highstg_consec)
  print(j)
}

bks=c(1,15,30,60,90,180,365)
rslt.stg16=ddply(highstg_consec,c("Alt","sum.stg16"),summarise,count.event=N.obs(sum.stg16))
rslt.stg16$cat=findInterval(rslt.stg16$sum.stg16,bks,left.open = FALSE,rightmost.closed = TRUE)
rslt.stg17=ddply(highstg_consec,c("Alt","sum.stg17"),summarise,count.event=N.obs(sum.stg17))
rslt.stg17$cat=findInterval(rslt.stg17$sum.stg17,bks,left.open = FALSE,rightmost.closed = TRUE)

rslt.stg16.sum=reshape2::dcast(subset(rslt.stg16,cat>0),cat~Alt,value.var="count.event",sum,na.rm=T)
rslt.stg17.sum=reshape2::dcast(subset(rslt.stg17,cat>0),cat~Alt,value.var="count.event",sum,na.rm=T)

# fill
tmp=rslt.stg17.sum[1,]
tmp[,1:5]=tmp[,1:5]=c(4,rep(0,4))
rslt.stg17.sum=rbind(rslt.stg17.sum,tmp)
tmp=rslt.stg17.sum[1,]
tmp[,1:5]=tmp[,1:5]=c(5,rep(0,4))
rslt.stg17.sum=rbind(rslt.stg17.sum,tmp)
tmp=rslt.stg17.sum[1,]
tmp[,1:5]=tmp[,1:5]=c(6,rep(0,4))
rslt.stg17.sum=rbind(rslt.stg17.sum,tmp)


rslt.stg16.sum=rslt.stg16.sum[,c("cat",alts)]
rslt.stg17.sum=rslt.stg17.sum[,c("cat",alts)]
xlabs=c("< 15", "15 - 30","30 - 60","60 - 90","90 - 180", "> 180")

# png(filename=paste0(plot.path,"LOK_highstg_events.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.5,2,0.5,1),oma=c(4,3,1,0.25),lwd=0.5);
layout(matrix(c(1:12),6,2,byrow=F))

# ylim.max=c(150,20,20,10,4)
for(i in 1:length(xlabs)){
  ylim.val=c(0,10);by.y=5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  
  tmp=t(rslt.stg16.sum[i,2:ncol(rslt.stg16.sum)])
  x=barplot(tmp,beside=T,col=adjustcolor(cols.alts,0.5),
            ylim=ylim.val,axes=F,
            names.arg = rep(NA,length(tmp)),
            space=c(0),yaxs="i",xaxs="i")
  text(x,tmp,tmp,pos=3,offset=0.25)
  axis_fun(2,ymaj,ymin,format(ymaj))
  if(i==length(xlabs)){axis_fun(1,x,x,alts,las=2,cex=0.8,line=-0.25)}else{axis_fun(1,x,x,NA)}
  mtext(side=2,line=2.5,paste0("# of Events\n",xlabs[i]," days"),cex=0.6)
  box(lwd=1)
  if(i==1){mtext(side=3,adj=0,"LOK Daily Stage > 16 Ft NGVD",cex=0.8)}
}
for(i in 1:length(xlabs)){
  ylim.val=c(0,10);by.y=5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
  
  tmp=t(rslt.stg17.sum[i,2:ncol(rslt.stg17.sum)])
  x=barplot(tmp,beside=T,col=adjustcolor(cols.alts,0.5),
            ylim=ylim.val,axes=F,
            names.arg = rep(NA,length(tmp)),
            space=c(0),yaxs="i",xaxs="i")
  text(x,tmp,tmp,pos=3,offset=0.25)
  axis_fun(2,ymaj,ymin,format(ymaj))
  if(i==length(xlabs)){axis_fun(1,x,x,alts,las=2,cex=0.8,line=-0.25)}else{axis_fun(1,x,x,NA)}
  # mtext(side=2,line=2.5,paste0("# of Events\n",xlabs[i]," days"),cex=0.75)
  box(lwd=1)
  if(i==1){mtext(side=3,adj=0,"LOK Daily Stage \u2265 17 Ft NGVD",cex=0.8)}
}
mtext(side=1,outer=T,line=3,"Alternative")

dev.off()