#########################################################
####### Figure 2 ########################################
###########################################################
###### Species richness analyses

#####################################################
##### Data import

setwd("C:/Users/torn/Downloads")
load("Renv_OTU_GF3.RData")

require(readxl)

########################
### New polished list
### A) all
dfAll<-as.data.frame(read_xlsx("relSpecies_all_polished_v3.xlsx", sheet="All but no cf v2 (richness)"))
row.names(dfAll)<-dfAll$Species
##### Remove unnecessary columns
df<-dfAll[,-c(1:4)]
#### seperate regimes withinthe same observer period
oldspecies<-df[,ENV$Year>=2009 & ENV$Year < 2013]
newspecies<-df[,ENV$Year>=2013 & ENV$Year <= 2017]
#############################
dfsub<-df[,ENV$Year>=2009 & ENV$Year <= 2017]
ENVsub<-ENV[ENV$Year>=2009 & ENV$Year <= 2017,]

## loop for richness (number of species/genera e.g.)
alpha1<-c()
for (i in 1:ncol(df)){
  alpha1[i]<-length(df[df[,i]>0,i])
}

####### Statistics
require(changepoint)
m<-cpt.meanvar(alpha1[ENV$Year>=2009 & ENV$Year <= 2017])
ENV$Date[ENV$Year>=2009 & ENV$Year <= 2017][44]
#################

##### Time series regul and deseas
require(pastecs)
Data2<-data.frame(ENV$Date, alpha1)
d<-difftime(as.Date(Data2[,1]), as.Date(Data2[1,1]), units="days")
Data2<-data.frame(as.vector(d), as.Date(Data2[,1]), as.numeric(Data2[,2])) 
str(Data2) ### Check that the values are in the right format
Data2[1,2]
tail(Data2)
regData<-regul(Data2[,1], 
               Data2[,3], 
               frequency=12, 
               units="daystoyears",
               datemin="01/11/2005", 
               dateformat ="d/m/Y", tol=7,#tol=13
               n=170, method="linear")
ts<-tseries(regData)

##subset to years where Diana counted
alpha2<-alpha1[ENV$Year >=2009 & ENV$Year <= 2017]
ENV2<-ENV[ENV$Year >=2009 & ENV$Year <= 2017,]

par(mfrow=c(1,1))
plot(ENV$Date, alpha1, type="l", col="grey", axes=F, xlab="Date", ylab="Species richness", ylim=c(0,60))
axis(1, at=as.POSIXct(paste0(2006:2020,"-01-01")), labels=c(2006:2020), las=2)
axis(2)
points(ENV$Date, alpha1, pch=16, col="grey")
points(ENV2$Date, alpha2, col="blue", type="l")
points(ENV2$Date, alpha2, col="blue", pch=16)
#points(seasdf$Date, seasdf$SST, type="l", lty=2, col="red")
segments(x0=as.POSIXct("2009-01-01"), y0=mean(alpha1[ENV$Year >=2009 & ENV$Date <= as.POSIXct("2012-12-18")]),
         x1=as.POSIXct("2013-01-01"), y1=mean(alpha1[ENV$Year >=2009 & ENV$Date <= as.POSIXct("2012-12-18")]), lwd=3, col="red")
segments(x0=as.POSIXct("2013-01-01"), y0=mean(alpha1[ENV$Date >= as.POSIXct("2012-12-18") & ENV$Year <= 2017]),
         x1=as.POSIXct("2018-01-01"), y1=mean(alpha1[ENV$Date >= as.POSIXct("2012-12-18") & ENV$Year <= 2017]), lwd=3, col="red")

#segments(x0=as.POSIXct("2006-01-1"), y0=mean(alpha1[ENV$Year< 2009]),
#         x1=as.POSIXct("2009-01-01"), y1=mean(alpha1[ENV$Year < 2009]), lwd=3, col="black", lty=2)

segments(x0=as.POSIXct("2006-01-01"), y0=mean(alpha1[ENV$Date <= as.POSIXct("2012-12-18")]),
         x1=as.POSIXct("2013-01-01"), y1=mean(alpha1[ENV$Date <= as.POSIXct("2012-12-18")]), lwd=3, col="pink")

segments(x0=as.POSIXct("2013-01-01"), y0=mean(alpha1[ENV$Date >= as.POSIXct("2012-12-18") ]),
         x1=as.POSIXct("2020-12-01"), y1=mean(alpha1[ENV$Date >= as.POSIXct("2012-12-18") ]), lwd=3, col="pink")


#abline(lm(alpha1[ENV$Year >=2009 & ENV$Date <= as.POSIXct("2013-01-28")] ~ ENV$Date[ENV$Year >=2009 & ENV$Date <= as.POSIXct("2013-01-28")]), lty=2, lwd=3, col="green")
#abline(lm(alpha1[ENV$Date >= as.POSIXct("2013-01-28") & ENV$Year <= 2017] ~ ENV$Date[ENV$Date >= as.POSIXct("2013-01-28") & ENV$Year <= 2017]), lty=2, lwd=3, col="green")

legend("topright", 
       legend=c("Counts by same observer (DiKr)","Counts by different observers", 
                "Mean of 2 regimes (same counter DiKr)", "Mean of 2 regimes (all data)"), 
       col=c("blue","grey", "red", "pink"), 
       lty=c(1,1,1,1), pch=c(16,16,NA, NA), lwd=c(1,1,2,2), box.lty = 0)
text(as.POSIXct(paste0("2011-01-01")), 12, labels = "Changepoint analysis")
text(as.POSIXct(paste0("2011-01-01")), 10, labels = "MBIC penalty= 13.7")
text(as.POSIXct(paste0("2011-01-01")), 8, labels = "Date: Dec 2012")
#title(adj=0.02, line=-1,main="A")



##########################################################
####### Figure 3 ##########################################
###########################################################
###### Atlantic type water plots
#############################################################
