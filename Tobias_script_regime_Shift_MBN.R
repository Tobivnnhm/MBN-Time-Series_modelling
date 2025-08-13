
load("Renv_OTU_GF3.RData")
nrow(relSpecies_all)
nrow(relGenus_all)

setwd("C:/Users/torn/Downloads")

############ SPecies level
#### polished species
# write.csv(relSpecies_all, file="relSpecies_all.csv")
require(readxl)
df<-as.data.frame(read_xlsx("relSpecies_all_polished.xlsx", sheet=1))
row.names(df)<-df$Species
colnames(df)
df<-df[,-c(1:4)]
### remove species which occur in less than 3 samples (ca 1% rounded up)
### only include species level taxa,excluding cf
sublist<-as.data.frame(read_xlsx("relSpecies_all_polished.xlsx", sheet=2))

df  <-df[row.names(df) %in% sublist[,1],]

########################
### New polished list
### A) all
dfAll<-as.data.frame(read_xlsx("relSpecies_all_polished_v3.xlsx", sheet=2))
row.names(dfAll)<-dfAll$Species
colnames(dfAll)
### B) specie sin less than 3 samples removed
str(dfAll)
df <- dfAll[dfAll$samples >= 3, ]
##### Remove unnecessary columns
dfAll<-dfAll[,-c(1:4)]
df<-df[,-c(1:4)]

oldspecies<-df[,ENV$Year>=2009 & ENV$Year < 2013]
rowSums(oldspecies)
newspecies<-df[,ENV$Year>=2013 & ENV$Year <= 2017]
rowSums(newspecies)

df2<-data.frame(old=rowMeans(oldspecies), new=rowMeans(newspecies))
df2<-df2[rowSums(df2) > 0,]

#write.csv(df2, file="regime_shift_species_polished.csv")

dfspecies<-df

#### Genus level
############## polsiehd genera
# write.csv(relGenus_all, file="relGenus_all.csv")

df<-as.data.frame(read_xlsx("relGenus_all_polished.xlsx", sheet=1))
row.names(df)<-df$Genus
colnames(df)
df<-df[,-c(1:4)]
### remove species which occur in less than 3 samples (ca 1% rounded up)
### only include species level taxa,excluding cf
sublist<-as.data.frame(read_xlsx("relGenus_all_polished.xlsx", sheet=2))

df  <-df[row.names(df) %in% sublist[,1],]

oldspecies<-df[,ENV$Year>=2009 & ENV$Year < 2013]
rowSums(oldspecies)
newspecies<-df[,ENV$Year>=2013 & ENV$Year <= 2017]
rowSums(newspecies)

df2<-data.frame(old=rowMeans(oldspecies), new=rowMeans(newspecies))
df2<-df2[rowSums(df2) > 0,]
#pmax(df2[23,])


#write.csv(df2, file="regime_shift_Genus_polished.csv")

dfGenus<-df

################################

summary(lm(colSums(dfGenus != 0)~colSums(dfspecies != 0)))
plot(colSums(dfGenus != 0)~colSums(dfspecies != 0))

##############
### decide for genus or species for followinfg steps

df<-dfspecies
#############################


dfsub<-df[,ENV$Year>=2009 & ENV$Year <= 2017]
ENVsub<-ENV[ENV$Year>=2009 & ENV$Year <= 2017,]

## loop for richness (number of species/genera e.g.)
alpha1<-c()
for (i in 1:ncol(df)){
  alpha1[i]<-length(df[df[,i]>0,i])
}

require(changepoint)
cpt.mean(alpha1[ENV$Year>=2009 & ENV$Year <= 2017])
m<-cpt.meanvar(alpha1[ENV$Year>=2009 & ENV$Year <= 2017])
str(m)
cpt.meanvar(alpha1[ENV$Year>=2009 & ENV$Year <= 2017])
cpt.var(alpha1[ENV$Year>=2009 & ENV$Year <= 2017])
ENV$Date[ENV$Year>=2009 & ENV$Year <= 2017][46]

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
segments(x0=as.POSIXct("2009-01-01"), y0=mean(alpha1[ENV$Year >=2009 & ENV$Date <= as.POSIXct("2013-02-26")]),
         x1=as.POSIXct("2013-01-01"), y1=mean(alpha1[ENV$Year >=2009 & ENV$Date <= as.POSIXct("2013-02-26")]), lwd=3, col="red")
segments(x0=as.POSIXct("2013-01-01"), y0=mean(alpha1[ENV$Date >= as.POSIXct("2013-02-26") & ENV$Year <= 2017]),
         x1=as.POSIXct("2018-01-01"), y1=mean(alpha1[ENV$Date >= as.POSIXct("2013-02-26") & ENV$Year <= 2017]), lwd=3, col="red")

#segments(x0=as.POSIXct("2006-01-1"), y0=mean(alpha1[ENV$Year< 2009]),
#         x1=as.POSIXct("2009-01-01"), y1=mean(alpha1[ENV$Year < 2009]), lwd=3, col="black", lty=2)

segments(x0=as.POSIXct("2006-01-01"), y0=mean(alpha1[ENV$Date <= as.POSIXct("2013-02-26")]),
         x1=as.POSIXct("2013-01-01"), y1=mean(alpha1[ENV$Date <= as.POSIXct("2013-02-26")]), lwd=3, col="pink")

segments(x0=as.POSIXct("2013-01-01"), y0=mean(alpha1[ENV$Date >= as.POSIXct("2013-02-26") ]),
         x1=as.POSIXct("2020-12-01"), y1=mean(alpha1[ENV$Date >= as.POSIXct("2013-02-26") ]), lwd=3, col="pink")


#abline(lm(alpha1[ENV$Year >=2009 & ENV$Date <= as.POSIXct("2013-01-28")] ~ ENV$Date[ENV$Year >=2009 & ENV$Date <= as.POSIXct("2013-01-28")]), lty=2, lwd=3, col="green")
#abline(lm(alpha1[ENV$Date >= as.POSIXct("2013-01-28") & ENV$Year <= 2017] ~ ENV$Date[ENV$Date >= as.POSIXct("2013-01-28") & ENV$Year <= 2017]), lty=2, lwd=3, col="green")

legend("topright", 
       legend=c("Counts by same observer (DiKr)","Counts by different observers", 
                "Mean of 2 regimes (same counter DiKr)", "Mean of 2 regimes (all data)"), 
       col=c("blue","grey", "red", "pink"), 
       lty=c(1,1,1,1), pch=c(16,16,NA, NA), lwd=c(1,1,2,2), box.lty = 0)
text(as.POSIXct(paste0("2011-01-01")), 12, labels = "Changepoint analysis")
text(as.POSIXct(paste0("2011-01-01")), 10, labels = "MBIC penalty= 13.7")
text(as.POSIXct(paste0("2011-01-01")), 8, labels = "Date: 2013-02-26")
#title(adj=0.02, line=-1,main="A")


###############################
#### indicspecies 2 diff abund taxa


########################
### New polished list
### A) all
dfAll<-as.data.frame(read_xlsx("relSpecies_all_polished_v3.xlsx", sheet=3))
row.names(dfAll)<-dfAll$Species
colnames(dfAll)
### B) specie sin less than 3 samples removed
str(dfAll)
df <- dfAll[dfAll$samples >= 3, ]
##### Remove unnecessary columns
dfAll<-dfAll[,-c(1:4)]
df<-df[,-c(1:4)]

oldspecies<-df[,ENV$Year>=2009 & ENV$Year < 2013]
rowSums(oldspecies)
newspecies<-df[,ENV$Year>=2013 & ENV$Year <= 2017]
rowSums(newspecies)

df2<-data.frame(old=rowMeans(oldspecies), new=rowMeans(newspecies))
df2<-df2[rowSums(df2) > 0,]

str(df2)
rownames(df2[df2$old == 0, ])
rownames(df2[df2$new == 0, ])

#write.csv(df2, file="regime_shift_species_polished.csv")

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
