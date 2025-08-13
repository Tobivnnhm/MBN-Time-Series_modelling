

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
dfAll<-dfAll[,-c(1:4)]
df<-dfAll

#### seperate regimes withinthe same observer period
oldspecies<-df[,ENV$Year>=2009 & ENV$Year < 2013]
newspecies<-df[,ENV$Year>=2013 & ENV$Year <= 2017]

##### unique species analyses
df2<-data.frame(old=rowMeans(oldspecies), new=rowMeans(newspecies))
df2<-df2[rowSums(df2) > 0,]
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
cpt.mean(alpha1[ENV$Year>=2009 & ENV$Year <= 2017])
m<-cpt.meanvar(alpha1[ENV$Year>=2009 & ENV$Year <= 2017])
str(m)
cpt.meanvar(alpha1[ENV$Year>=2009 & ENV$Year <= 2017])
cpt.var(alpha1[ENV$Year>=2009 & ENV$Year <= 2017])
ENV$Date[ENV$Year>=2009 & ENV$Year <= 2017][44]

require(diptest)
dip.test(alpha1[ENV$Year>=2009 & ENV$Year <= 2017])
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
plot(regData)
ts<-tseries(regData)
#plot(decompose(ts))
dec<-decompose(ts)
dec$seasonal

require(lubridate)
seasdf<-data.frame(Date=date_decimal(as.numeric(time(ts))),
                   SSTd=as.numeric(ts),
                   SST=as.numeric(dec$seasonal)+mean(as.numeric(ts)))


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









#######################################################
# ######### #########  
# #             #
# ##########    #
#         ##    #     
# ##########    #      
#######################################################


################################################################
################################################################
####### Statistics on richness

## loop for richness (number of species/genera e.g.)
alpha1<-c()
for (i in 1:ncol(df)){
  alpha1[i]<-length(df[df[,i]>0,i])
}

##subset to years where Diana counted
alpha2<-alpha1[ENV$Year >=2009 & ENV$Year <= 2017]
ENV2<-ENV[ENV$Year >=2009 & ENV$Year <= 2017,]

plot(ENV2$Date, alpha2, type="b")
alpha3<-data.frame(alpha=alpha1, Date=as.character(ENV$Date))


# regularize
alphats<-data.frame(alpha=alpha2)
alphats$Date<-as.Date(ENV2$Date)
alphats$d<-difftime(as.Date(alphats$Date), as.Date(alphats$Date[1]), units="days") ### convert dates to days since first sampling
## 171 months
require(pastecs)
regData<-regul(as.vector(alphats$d), alphats$alpha, frequency=12, units="daystoyears",
               datemin="14/01/2009", 
               dateformat ="d/m/Y", tol=13,#tol=13
               n=9*12, method="spline")
plot(regData)

ts0<-tseries(regData)
plot(ts0)
plot(decompose(ts0))

#### deseasonalize
require(forecast)
require(lubridate)
decalpha<-seasadj(decompose(ts0))

decalphadf<-data.frame(alpha=as.numeric(decalpha), Date=as.numeric(time(decalpha)))
alphadf<-data.frame(alpha=as.numeric(ts0), Date=as.numeric(time(ts0)))

alphadf$Date<-seq(as.Date("2009-01-15"), as.Date("2017-12-15"), by="month")

par(mfrow=c(2,1))
plot(decalpha)
plot(alpha2, col="red", lty=2, type="l")

plot(alphadf$Date, alphadf$alpha, type="b")
points(decalphadf$Date, decalphadf$alpha, type="b", col="red")
### no deseasonalization needed
#### USE: alphadf

str(alphadf)


##################################################################
### 3) Nitrate & Co
##################################################################

MBN_physchem<-read.csv("MBN_2D_data.csv")
MBN_physchem$Month<-month(MBN_physchem$Date)
MBN_physchem$Year<-year(MBN_physchem$Date)
MBN1<-MBN_physchem

Allraw<-merge(MBN1, alpha3, by="Date", all=T)

str(MBN1)
str(alpha3)

MBN_physchemreg<-read.csv("MBN_2D_data_regul.csv")
MBN_physchemreg$Month<-month(MBN_physchemreg$Date)
MBN_physchemreg$Year<-year(MBN_physchemreg$Date)
MBN2<-MBN_physchemreg
str(MBN2)
MBN2$Date<-as.Date(MBN2$Date)
Allreg<-merge(MBN2, alpha3, by="Date", all=T)
str(Allreg)
Allreg$Date<-as.Date(paste0(substr(Allreg$Date,start=1, stop=7),"-02"))

Allreg2<-merge(Allreg, MBN2, by="Date", all=T)
View(Allreg2)

####################################################
################################################
### 1) trend analyses ###############################
#### test for a monotonous trend over time

require(EnvStats)
print(kendallSeasonalTrendTest(alpha1, season=ENV$Month, year= ENV$Year))
require(Kendall)
kendallTrendTest(alpha1)
print(kendallSeasonalTrendTest(noaa$NPP, season=noaa$Month, year= noaa$Year))
print(kendallSeasonalTrendTest(MBN1$NOx_int, season=MBN1$Month, year= MBN1$Year))
print(kendallSeasonalTrendTest(MBN1$SiOH4_int, season=MBN1$Month, year= MBN1$Year))

str(Allraw)

### all years
for (i in 3:ncol(Allraw)){
  dat<-data.frame(var=Allraw[,i], Month=Allraw$Month, Year=Allraw$Year)
  dat<-na.omit(dat)
  print(colnames(Allraw)[i])
  print(kendallSeasonalTrendTest(dat$var, season=dat$Month, year= dat$Year))
}

#### before tipping point
for (i in 3:ncol(Allraw)){
  dat<-data.frame(var=Allraw[,i], Month=Allraw$Month, Year=Allraw$Year)
  dat<-na.omit(dat)
  dat<-dat[dat$Year <= 2012,]
  print(colnames(Allraw)[i])
  print(kendallSeasonalTrendTest(dat$var, season=dat$Month, year= dat$Year))
}

#### after tipping point
for (i in 3:ncol(Allraw)){
  dat<-data.frame(var=Allraw[,i], Month=Allraw$Month, Year=Allraw$Year)
  dat<-na.omit(dat)
  dat<-dat[dat$Year > 2012,]
  print(colnames(Allraw)[i])
  print(kendallSeasonalTrendTest(dat$var, season=dat$Month, year= dat$Year))
}

print(kendallSeasonalTrendTest(noaa$NPP, season=noaa$Month, year= noaa$Year))
print(kendallSeasonalTrendTest(noaa$NPP[noaa$Year<=2012], season=noaa$Month[noaa$Year<=2012], year= noaa$Year[noaa$Year<=2012]))
print(kendallSeasonalTrendTest(noaa$NPP[noaa$Year>2012], season=noaa$Month[noaa$Year>2012], year= noaa$Year[noaa$Year>2012]))

###########################################################################
#################### Tipping point/ regime shift diagnostics #############
#### is there a new stable system? #########################################

#########################################################################
##### Bimodal distribution ################################################
###https://www.nature.com/articles/s41559-020-1273-8
#### the data before and aftwer the tipping point should have diffferent unimodal distribution -> overall bimodal distribution

require(diptest)
hist(alpha1[ENV$Year>=2009 & ENV$Year <= 2017])
dip.test(alpha1[ENV$Year>=2009 & ENV$Year <= 2017])
#### bi modal distribution (--> evidence for tipping point)

################################################################
####### changepoint analyses ######################################
####https://www.nature.com/articles/s41598-021-93843-z#Sec7
### detect sudden/ stepwise changes

require(changepoint)

cpt.mean(alpha1[ENV$Year>=2009 & ENV$Year <= 2017])
cpt.meanvar(alpha1[ENV$Year>=2009 & ENV$Year <= 2017])
cpt.var(alpha1[ENV$Year>=2009 & ENV$Year <= 2017])
ENV[ENV$Year>=2009 & ENV$Year <= 2017,][46,]
ENV[79,] ### mean changepoint in march 2013
ENV[87,] ### meanvar changepoint in dec 2013
ENV[ENV$Year>=2009 & ENV$Year <= 2017,][46,]


######################################
##### structure change/ breakpoint analyses
#############https://www.nature.com/articles/s41598-021-93843-z#Sec7
#### after a tipping point correlations between parameters should be different

require(strucchange)

require(dplyr)
merge_alpha_data <- function(df) {
  # Group by Date
  grouped_df <- df %>%
    group_by(Date) %>%
    summarise(
      Temp5 = first(na.omit(Temp5)),
      Sal5 = first(na.omit(Sal5)),
      SiOH4_int = first(na.omit(SiOH4_int)),
      NOx_int = first(na.omit(NOx_int)),
      PO4_int = first(na.omit( PO4_int)),
      alpha = first(na.omit(alpha)),
      across(where(is.numeric), ~ first(na.omit(.))) # Keep other numeric columns
    ) %>%
    ungroup() # Remove grouping
  
  return(grouped_df)
}

Allreg3 <- merge_alpha_data(Allreg)



MBN3<-data.frame(Allreg3$Date, Allreg3$alpha, Allreg3$Temp5)
MBN3<-data.frame(Allreg3$Date, Allreg3$alpha, Allreg3$Sal5)
MBN3<-data.frame(Allreg3$Date, Allreg3$alpha, Allreg3$NOx_int)
MBN3<-data.frame(Allreg3$Date, Allreg3$alpha, Allreg3$SiOH4_int)

MBN3<-na.omit(MBN3)
MBN3<-MBN3[year(MBN3$Allreg3.Date) >=2009 & year(MBN3$Allreg3.Date) <=2017,]
breakpoints(MBN3[,2] ~ MBN3[,3], breaks = 1)
summary(breakpoints(MBN3[,2] ~ MBN3[,3], breaks = 1))

i=37

MBN3$Allreg3.Date[i]

par(mfrow=c(2,2))
plot(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3], pch=16,
     xlab="reg. Surface Temperatue at 0-5m [°C]", ylab="reg. species richness",
     xlim=c(min(MBN3[,3]),max( MBN3[,3])), ylim=c(0,max( MBN3[,2])))
abline(lm(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3]))
points(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3], col="red", pch=16)
abline(lm(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3]), col="red")

summary(lm(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3]))
summary(lm(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3]))

text(2,55, labels = paste0("breakpoint at ",paste0(MBN3$Allreg3.Date[i])))

title(adj=0.02, line=-1,main="a)")




plot(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3], pch=16,
     xlab="reg. Surface Salinity at 0-5m", ylab="reg. species richness",
     xlim=c(min(MBN3[,3]),max( MBN3[,3])), ylim=c(0,max( MBN3[,2])))
abline(lm(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3]))
points(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3], col="red", pch=16)
abline(lm(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3]), col="red")

summary(lm(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3]))
summary(lm(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3]))

text(27,55, labels = paste0("breakpoint at ",paste0(MBN3$Allreg3.Date[i])))

title(adj=0.02, line=-1,main="b)")


plot(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3], pch=16,
     xlab="reg. Mean NOX 0-25 m [µmol L-1]", ylab="reg. species richness",
     xlim=c(min(MBN3[,3]),max( MBN3[,3])), ylim=c(0,max( MBN3[,2])))
abline(lm(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3]))
points(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3], col="red", pch=16)
abline(lm(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3]), col="red")

summary(lm(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3]))
summary(lm(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3]))

text(3.5,55, labels = paste0("breakpoint at ",paste0(MBN3$Allreg3.Date[i])))

title(adj=0.02, line=-1,main="c)")



plot(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3], pch=16,
     xlab="Mean Silicate 0-25 m [µmol L-1]", ylab="species richness",
     xlim=c(min(MBN3[,3]),max( MBN3[,3])), ylim=c(0,max( MBN3[,2])))
abline(lm(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3]))
points(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3], col="red", pch=16)
abline(lm(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3]), col="red")

summary(lm(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3]))
summary(lm(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3]))

text(2.5, 55, labels = paste0("breakpoint at ",paste0(MBN3$Allreg3.Date[i])))

title(adj=0.02, line=-1,main="d)")






################################################################



MBN3<-data.frame(Allreg3$Date, Allreg3$NOx_int, Allreg3$Temp5)
MBN3<-na.omit(MBN3)
breakpoints(MBN3[,2] ~ MBN3[,3], breaks = 1)
summary(breakpoints(MBN3[,2] ~ MBN3[,3], breaks = 1))
i=71

par(mfrow=c(2,2))
plot(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3], pch=16,
     xlab="reg. Surface Temperatue at 0-5m [°C]", ylab="reg. Mean NOX 0-25 m [µmol L-1]",
     xlim=c(min(MBN3[,3]),max( MBN3[,3])), ylim=c(0,max( MBN3[,2])))
abline(lm(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3]))
points(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3], col="red", pch=16)
abline(lm(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3]), col="red")

summary(lm(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3]))
summary(lm(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3]))

text(2,11, labels = paste0("breakpoint at ",paste0(MBN3$Allreg3.Date[i])))

title(adj=0.02, line=-1,main="a)")


MBN3<-data.frame(Allreg3$Date, Allreg3$NOx_int, Allreg3$Sal5)
MBN3<-na.omit(MBN3)
breakpoints(MBN3[,2] ~ MBN3[,3], breaks = 1)
summary(breakpoints(MBN3[,2] ~ MBN3[,3], breaks = 1))
i=71
plot(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3], pch=16,
     xlab="reg. Surface Salinity at 0-5m [°C]", ylab="reg. Mean NOX 0-25 m [µmol L-1]",
     xlim=c(min(MBN3[,3]),max( MBN3[,3])), ylim=c(0,max( MBN3[,2])))
abline(lm(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3]))
points(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3], col="red", pch=16)
abline(lm(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3]), col="red")

summary(lm(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3]))
summary(lm(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3]))

text(26,11, labels = paste0("breakpoint at ",paste0(MBN3$Allreg3.Date[i])))

title(adj=0.02, line=-1,main="b)")


MBN3<-data.frame(Allreg3$Date, Allreg3$NOx_int, Allreg3$SiOH4_int)
MBN3<-na.omit(MBN3)
breakpoints(MBN3[,2] ~ MBN3[,3], breaks = 1)
summary(breakpoints(MBN3[,2] ~ MBN3[,3], breaks = 1))
i=68
plot(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3], pch=16,
     xlab="reg. Mean Silicate 0-25 m [µmol L-1]", ylab="reg. Mean NOX 0-25 m [µmol L-1]",
     xlim=c(min(MBN3[,3]),max( MBN3[,3])), ylim=c(0,max( MBN3[,2])))
abline(lm(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3]))
points(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3], col="red", pch=16)
abline(lm(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3]), col="red")

summary(lm(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3]))
summary(lm(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3]))

text(2,11, labels = paste0("breakpoint at ",paste0(MBN3$Allreg3.Date[i])))

title(adj=0.02, line=-1,main="c)")

MBN3<-data.frame(Allreg3$Date, Allreg3$NOx_int, Allreg3$PO4_int)
MBN3<-na.omit(MBN3)
breakpoints(MBN3[,2] ~ MBN3[,3], breaks = 1)
summary(breakpoints(MBN3[,2] ~ MBN3[,3], breaks = 1))
i=70
plot(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3], pch=16,
     xlab="reg. Mean PO4 0-25 m [µmol L-1]", ylab="reg. Mean NOX 0-25 m [µmol L-1]",
     xlim=c(min(MBN3[,3]),max( MBN3[,3])), ylim=c(0,max( MBN3[,2])))
abline(lm(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3]))
points(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3], col="red", pch=16)
abline(lm(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3]), col="red")

summary(lm(MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date <= paste0(MBN3$Allreg3.Date[i]),3]))
summary(lm(MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),2] ~ MBN3[MBN3$Allreg3.Date > paste0(MBN3$Allreg3.Date[i]),3]))

text(0.3,11, labels = paste0("breakpoint at ",paste0(MBN3$Allreg3.Date[i])))

title(adj=0.02, line=-1,main="d)")

























################################################################################################
###############################
#### indicspecies 2 diff abund taxa


########################
### New polished list
### A) all
dfAll<-as.data.frame(read_xlsx("relSpecies_all_polished_v3.xlsx", sheet="simplified for indicator"))

row.names(dfAll)<-dfAll$Species
colnames(dfAll)
### B) specie sin less than 3 samples removed
str(dfAll)
df <- dfAll[dfAll$samples >= 3, ]
##### Remove unnecessary columns
dfAll<-dfAll[,-c(1:4)]
df<-df[,-c(1:4)]

####################
#### ALL
##### Subsets of old and new regime
oldspecies<-dfAll[,ENV$Year>=2009 & ENV$Year < 2013]
rowSums(oldspecies)
newspecies<-dfAll[,ENV$Year>=2013 & ENV$Year <= 2017]
rowSums(newspecies)



##### merge them and remove species not present between 2009 and 2017, show the Mean rel abundance
df2<-data.frame(old=rowMeans(oldspecies), new=rowMeans(newspecies), both=rowMeans(dfAll[,ENV$Year>=2009 & ENV$Year < 2013]), 
                Allmean=rowMeans(dfAll), Allmax=apply(dfAll, 1, max, na.rm = TRUE))
df22<-df2
df2<-df2[rowSums(df2) > 0,]

##### merge them and remove species not present between 2009 and 2017, but now show the sample numbers
df3 <- data.frame(
  old = apply(oldspecies, 1, function(row) sum(row > 0, na.rm = TRUE)),
  new = apply(newspecies, 1, function(row) sum(row > 0, na.rm = TRUE))
)
df3<-df3[rowSums(df3) > 0,]
  df3$Both<-rowSums(df3)

##### Now remove species in less (<=) than 2 samples
df4<-df3[df3$Both >2,]
df4<-df3[rownames(df3) %in% rownames(df),]

unique_rownames_df3 <- setdiff(rownames(df3), rownames(df4))
result_table <- df3[unique_rownames_df3, , drop = FALSE]
result_table[result_table$old == 0 | result_table$new == 0,]
nrow(result_table[result_table$old == 0 | result_table$new == 0,])

# merge dataframes to show df4 with mean rel abundance AND numbe rof sample present
merged_df <- merge(df2[rownames(df4), , drop = FALSE],
                   df4[rownames(df4), , drop = FALSE],
                   by = "row.names", all = FALSE)
rownames(merged_df) <- merged_df$Row.names
merged_df$Row.names <- NULL
nrow(merged_df)

######################################################
##### Numbers for the Results of the Manuscript
########################################################
##### Total species #### For MS text
nrow(dfAll)  ### species found in total
#####################
##### Total species only in 2009-2017#### For MS text
nrow(dfAll)-nrow(df2)  ### species were removed because npot present between 2009-2017
nrow(df2)  #### remaining species
##########################
##### Species unique in old and new regime
#####################
sum(df2$new == 0, na.rm = TRUE)  ##### species only present only present before 2013
sum(df2$new == 0, na.rm = TRUE)/nrow(df2) * 100 ### in %
sum(df2$old == 0, na.rm = TRUE) #### species only present after 2013
sum(df2$old == 0, na.rm = TRUE)/nrow(df2) * 100 ### in %
###### After removing species in <=2 samples
#sum(df3$Both <= 2, na.rm = TRUE)
#sum(df3$Both <= 2, na.rm = TRUE)/nrow(df3)
1-(nrow(df4)/nrow(df3))
nrow(df4)
######################
sum(df4$new == 0, na.rm = TRUE)  ##### species only present only present before 2013
sum(df4$new == 0, na.rm = TRUE)/nrow(df4) * 100 ### in %
sum(df4$old == 0, na.rm = TRUE) #### species only present after 2013
sum(df4$old == 0, na.rm = TRUE)/nrow(df4) * 100 ### in %
#############
nrow(merged_df[merged_df$new.x == 0,])
#nrow(merged_df[merged_df$new.x == 0 & merged_df$old.x < 0.1,])
nrow(merged_df[merged_df$new.x == 0 & merged_df$Allmean < 0.1,])
merged_df[merged_df$new.x == 0 & merged_df$Allmean > 1,]


##################################
###### Now indicator taxa

library(indicspecies)

dfsub<-dfAll[,ENV$Year>=2009 & ENV$Year <= 2017]
ENVsub<-ENV[ENV$Year>=2009 & ENV$Year <= 2017,]

ind<-as.data.frame(t(dfsub))
ind$Date<-as.Date(ENVsub$Date)
ind$group<-rep(NA, nrow(ind))
for (i in 1:nrow(ind)){
  if (as.Date(ind$Date[i]) <= as.Date("2012-12-18")){
    ind$group[i] <-"old"
  } else { ind$group[i] <-"new" }
}


abund<- as.data.frame(t(dfsub))
time<-ind$group

inv = multipatt(abund, time, func = "r.g", control = how(nperm=9999))
summary(inv)

#### Data frame with indicator species
dfinv<-na.omit(inv$sign[inv$sign$p.value <= 0.05,])

#### merge with Full list (After removing doubletons)
merged_df$a<-rownames(merged_df)
dfinv$a<-rownames(dfinv)

merged_df <- merge(merged_df, dfinv, by = "a", all = TRUE)
colnames(merged_df)[1]<-"Species"
d<-merged_df
d[is.na(d)] <- 0

####################################
##### Numbers for manuscript
############################
nrow(dfinv[dfinv$s.old==0,])
nrow(dfinv[dfinv$s.new==0,])
###################
nrow(na.omit(d[d$s.old == 1 & d$new.y == 0,]))
nrow(na.omit(d[d$s.old == 1 & d$new.y > 0 ,]))
na.omit(d[d$s.old == 1 & d$new.y > 0 & d$both <0.1,])
na.omit(d[d$s.old == 1 & d$new.y > 0 & d$both >1,])
#################################
(d[d$s.old == 1 & d$new.y > 0,])  #### Indicator only
nrow((d[d$s.old == 1 & d$new.y > 0,])) 
(d[d$s.old == 1 & d$new.y == 0,])  #### Indicator and unique
nrow((d[d$s.old == 1 & d$new.y == 0,]))
(d[d$s.old == 0 & d$new.y == 0,])  #### unique only
nrow((d[d$s.old == 0 & d$new.y == 0,]))
nrow((d[d$s.old == 1 | d$new.y == 0,]) ) ### All
#################################
(d[d$s.new == 1 & d$old.y > 0,])  #### Indicator only
(d[d$s.new == 1 & d$old.y == 0,])  #### Indicator and unique
(d[d$s.new == 0 & d$old.y == 0,])  #### unique only
#### New regime
d[d$s.new == 1 & d$Allmean <0.1,]
d[d$s.new == 1 & d$Allmean >1,]
d[d$old.x == 0 & d$Allmean <0.1,]
d[d$old.x == 0 & d$Allmean >1,]

d[d$s.old == 1 &  d$Allmean <0.1,]
d[d$new.x == 0 & d$s.old == 1 & d$Allmean <0.1,]
d[d$new.x == 0 & d$s.old ==0 & d$Allmean <0.1,]
d[d$new.x > 0 & d$s.old == 1 & d$Allmean <0.1,]

d[d$s.old == 1 & d$new.x > 0 & d$Allmean <0.1,]
d[d$s.old == 1 & d$new.x > 0 & d$Allmean >1,]
#write.csv(df2, file="regime_shift_species_polished.csv")




##################################################
#### Overview table





MTable<-data.frame(Species=rownames(dfAll),
                   Samples=as.numeric(apply(dfAll, 1, function(row) sum(row > 0, na.rm = TRUE))),
                   max=as.numeric(apply(dfAll, 1, max, na.rm = TRUE)),
                   mean=as.numeric(rowMeans(dfAll)),
                   uniquePre2013=0,
                   indPre2013=0,
                   uniquePost2013=0,
                   indPost2013=0,
                   uniqueObs=0
                   )

dd<-merge(d, MTable, by="Species", all=T)
View(dd)

dd$uniquePre2013[dd$new.y == 0]<-1

dd$uniquePost2013[dd$old.y == 0]<-1
dd$uniqueObs[is.na(dd$old.x)]<-1
dd$indPre2013[dd$s.old == 1]<-1
dd$indPost2013[dd$s.new == 1]<-1

colnames(dd)

MTable<-dd[,c(1,15:22)]

View(MTable)

write.csv(MTable, file="Indicatro_Unique_Table_2025_v3.csv")



































dfspecies<-df

df<-dfspecies
#############################


dfsub<-df[,ENV$Year>=2009 & ENV$Year <= 2017]
ENVsub<-ENV[ENV$Year>=2009 & ENV$Year <= 2017,]


library(indicspecies)

ind<-as.data.frame(t(dfsub))
ind$Date<-as.Date(ENVsub$Date)
ind$group<-rep(NA, nrow(ind))
for (i in 1:nrow(ind)){
  if (as.Date(ind$Date[i]) <= as.Date("2013-04-16")){
    ind$group[i] <-"old"
  } else { ind$group[i] <-"new" }
}


abund<- as.data.frame(t(dfsub))
time<-ind$group

inv = multipatt(abund, time, func = "r.g", control = how(nperm=9999))
summary(inv)
