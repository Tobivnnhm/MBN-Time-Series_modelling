####################################################
### Tipping point /regime shift analyses
####################################################

##############################
###### Pipeline###############
############################
### 1) Regularize and deseasonalize data if needed (test for raw vs deseasonalized data)
### 3) Monotonous trend analyses (seasonally adjusted Mann-Kendall test)
### 4) Tipping point/ regime shift analyses (univariate)
###    a) bimodal dsitribution of data (dip.test, histogram)
###    b) changepoint analyses
###    c) Differnces before and after the tipping point (kruskall wallis)
### 5) Tipping point/ regime shift analyses (correlations, multivariate)
###    a) breakpoint analyses (Are correlations changing?)
###    b) catastrophe theory/ cusp model (is a linear change in drivers leading to a sudden change in response variable?)


####################################################
##### Data import and standardization #############
setwd("C:/Users/torn/Documents")

###################################################################
### 1) Alpha diversity (Richness)
##################################################################

load("Renv_OTU_GF3.RData")

df<-relSpecies_all
df<-relGenus_all
df<-relFamily_all
df<-relOrder_all
df<-relClass_all

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


###################################################################
### 2) Primary production (Cpmb2 model)
##################################################################
require(lubridate)
noaa<-read.csv("dataNPP_CbPM2.csv")[,1:2]
noaa$month<-month(noaa$DATE.YMD)
noaa$year<-year(noaa$DATE.YMD)
colnames(noaa)<-c("Date","NPP","Month","Year")
noaa$Date<-as.Date(noaa$Date)

str(noaa)
str(alphadf)

Allreg<- merge(noaa, alphadf, by="Date", all=T)
Allreg<-Allreg[,1:5]

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
hist(decalpha)#, breaks=c(0:40))
hist(noaa$NPP)
hist(MBN1$NOx_int, breaks=c(0:18))
hist(MBN1$SiOH4_int, breaks=c(0:7))

dip.test(alpha1[ENV$Year>=2009 & ENV$Year <= 2017])
dip.test(decalpha)
dip.test(noaa$NPP)
dip.test(MBN1$NOx_int)
dip.test(MBN1$SiOH4_int)

for (i in c(8:ncol(Allreg))){
  print(colnames(Allreg)[i])
  print(dip.test(na.omit(Allreg[,i])))
  hist(na.omit(Allreg[,i]))
}
str(Allreg)
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

cpt.mean(na.omit(decalpha))
cpt.meanvar(na.omit(decalpha))

cpt.mean(noaa$NPP)
cpt.meanvar(noaa$NPP)
noaa[129,]

cpt.mean(na.omit(MBN1$NOx_int))
cpt.meanvar(na.omit(MBN1$NOx_int))
MBN1[79,]

for (i in 2: 11){
  print(colnames(MBN2)[i])
print(cpt.meanvar(na.omit(MBN2[,i])))
}
MBN2[182,]
plot(MBN1[,2], ylim=c(0,75), type="b")
##### statistical confirmed tipping point :)

##############################################
#### Differences before and after

### alpha
alpha3<-data.frame(alpha=alpha1[ENV$Year>=2009 & ENV$Year <= 2017])
alpha3$group<-(c(rep("old",length(alpha1[ENV$Year>=2009 & ENV$Year <= 2012])), rep("new",length(alpha1[ENV$Year>2012 & ENV$Year <= 2017]))))

kruskal.test(alpha3$alpha, alpha3$group)
aggregate(alpha3, by=list(alpha3$group), FUN=mean)
aggregate(alpha3, by=list(alpha3$group), FUN=var)

#### NoAA NPP
noaa3<-data.frame(NPP=noaa$CHLA..CbPM2.NPP..mg.C.m2.day.)
noaa3$group<-c(rep("old",nrow(noaa[noaa$year <= 2012,])), rep("new",nrow(noaa[noaa$year>2012,])))
noaa3<-na.omit(noaa3)

kruskal.test(noaa3$NPP, noaa3$group)
aggregate(noaa3, by=list(noaa3$group), FUN=mean)
aggregate(noaa3, by=list(noaa3$group), FUN=var)

# NOx
Allraw3<-data.frame(NOx=Allraw$NOx_int, Year=Allraw$Year)
Allraw3$group<-c(rep("old",143), rep("new",146))
Allraw3<-na.omit(Allraw3)

kruskal.test(Allraw3$NOx, Allraw3$group)
aggregate(Allraw3, by=list(Allraw3$group), FUN=mean)
aggregate(Allraw3, by=list(Allraw3$group), FUN=var)

Allreg3<-data.frame(NOx=MBN2$NOx_int, Year=MBN2$Year)
Allreg3$group<-c(rep("old",87), rep("new",96))
Allreg3<-na.omit(Allreg3)

kruskal.test(Allreg3$NOx, Allreg3$group)
aggregate(Allreg3, by=list(Allreg3$group), FUN=mean)
aggregate(Allreg3, by=list(Allreg3$group), FUN=var)

# euphZ
Allraw3<-data.frame(NOx=Allraw$euphz, Year=Allraw$Year)
Allraw3$group<-c(rep("old",143), rep("new",146))
Allraw3<-na.omit(Allraw3)

kruskal.test(Allraw3$NOx, Allraw3$group)
aggregate(Allraw3, by=list(Allraw3$group), FUN=mean)
aggregate(Allraw3, by=list(Allraw3$group), FUN=var)


######################################
##### structure change/ breakpoint analyses
#############https://www.nature.com/articles/s41598-021-93843-z#Sec7
#### after a tipping point correlations between parameters should be different

require(strucchange)

MBN3<-data.frame(Allraw$Date, Allraw$NOx_int, Allraw$Chl_int)
MBN3<-data.frame(Allraw$Date, Allraw$alpha, Allraw$Chl_int)
MBN3<-data.frame(Allreg$Date, Allreg$NPP, Allreg$alpha)

MBN3<-data.frame(Allraw$Date, Allraw$NOx_int, Allraw$StratI)
MBN3<-data.frame(Allraw$Date, Allraw$alpha, Allraw$StratI)

MBN3<-data.frame(Allreg2$Date, Allreg2$NOx_int, Allreg2$StratI)
MBN3<-data.frame(Allreg2$Date, Allreg2$NOx_int, Allreg2$Temp5)
MBN3<-data.frame(Allreg2$Date, Allreg2$alpha, Allreg2$StratI)


MBN3<-data.frame(Allreg$Date, Allreg$NPP, Allreg$StratI)

MBN3<-data.frame(Allraw$Date, Allraw$NOx_int, Allraw$Temp5)
MBN3<-data.frame(Allraw$Date, Allraw$alpha, Allraw$Temp5)
MBN3<-data.frame(Allreg$Date, Allreg$NPP, Allreg$Temp5)

MBN3<-data.frame(Allraw$Date, Allraw$NOx_int, Allraw$SiOH4_int)
MBN3<-data.frame(Allraw$Date, Allraw$alpha, Allraw$SiOH4_int)
MBN3<-data.frame(Allreg$Date, Allreg$NPP, Allreg$Temp5)

MBN3<-na.omit(MBN3)
MBN3$Allraw.Date<-as.Date(MBN3$Allraw.Date)
MBN3$Allraw.Date<-as.Date(MBN3$Allreg2.Date)

breakpoints(MBN3[,2] ~ MBN3[,3])
str(MBN3)

MBN3[56,]
plot(MBN3[MBN3$Allraw.Date <= "2012-12-30",2] ~ MBN3[MBN3$Allraw.Date <= "2012-12-30",3], pch=16)
abline(lm(MBN3[MBN3$Allraw.Date <= "2012-12-30",2] ~ MBN3[MBN3$Allraw.Date <= "2012-12-30",3]))
points(MBN3[MBN3$Allraw.Date > "2012-12-30",2] ~ MBN3[MBN3$Allraw.Date > "2012-12-30",3], col="red", pch=16)
abline(lm(MBN3[MBN3$Allraw.Date > "2012-12-30",2] ~ MBN3[MBN3$Allraw.Date > "2012-12-30",3]), col="red")

summary(lm(MBN3[MBN3$Allraw.Date <= "2012-12-30",2] ~ MBN3[MBN3$Allraw.Date <= "2012-12-30",3]))
summary(lm(MBN3[MBN3$Allraw.Date > "2012-12-30",2] ~ MBN3[MBN3$Allraw.Date > "2012-12-30",3]))

#########################################################################
#### stochastic cusp model (SCM) to investigate the stability ????
#############https://www.nature.com/articles/s41598-021-93843-z#Sec7
###################################################################

## We validated the fitted SCM by assessing 
### (i) the significance of SSB in the linear model of zt, 
### (ii) evidence for the existence of bimodality in zt in the cusp area, 
### (iii) the percentage of observations in the cusp area, and 
### (iv) the goodness of the SCM fit using Cobb's pseudo-R2. 
### (v) Moreover, we compared the fitted SCM to alternative linear and logistic regression models often used to confront linear and continuous dynamics with the non-linear discontinuous regime shift case.

MBN3<-data.frame(Allraw$Date, Allraw$NOx_int, Allraw$Temp5, Allraw$StratI, Allraw$SiOH4_int, Allraw$PO4_int, Allraw$Chl_int)
MBN3<-data.frame(Allreg$Date, Allreg$NOx_int, Allreg$Temp5, Allreg$StratI, Allreg$SiOH4_int, Allreg$PO4_int, Allreg$Chl_int)
MBN3<-data.frame(Allraw$Date, Allraw$alpha, Allraw$Temp5, Allraw$StratI, Allraw$SiOH4_int, Allraw$PO4_int, Allraw$Chl_int, Allraw$NOx_int)
MBN3<-MBN3[Allraw$Year >=2009 & Allraw$Year <=2017,]
MBN3<-data.frame(Allreg$Date, Allreg$alpha, Allreg$Temp5, Allreg$StratI, Allreg$SiOH4_int, Allreg$PO4_int, Allreg$Chl_int, Allreg$NOx_int)
MBN3<-data.frame(Allreg$Date, Allreg$NPP, Allreg$Temp5, Allreg$StratI)

MBN3<-na.omit(MBN3)

require(cusp)
MBN4<-MBN3[,2:ncol(MBN3)]

#### y : state variable
### alpha: nominal/control variable
#### beta: bifurcation variable
####m -1: removes intercept

### NOx  model (raw)
colnames(MBN4)
fit = cusp(y~MBN4[,1], alpha~MBN4[,4]+MBN4[,5]+MBN4[,6], beta~MBN4[,2]+MBN4[,3])
fit = cusp(y~MBN4[,1], alpha~MBN4[,4]+MBN4[,6], beta~MBN4[,2])
fit = cusp(y~MBN4[,1], alpha~MBN4[,5], beta~MBN4[,3]-1)

summary(fit)
summary(fit, logist=T)
plot(fit)
cusp3d(fit)

### NOx  model (reg)
colnames(MBN4)
fit = cusp(y~MBN4[,1], alpha~MBN4[,4]+MBN4[,5]+MBN4[,6]+MBN4[,3], beta~MBN4[,2]+MBN4[,3])
fit = cusp(y~MBN4[,1], alpha~MBN4[,4]+MBN4[,5]+MBN4[,6]+MBN4[,3], beta~MBN4[,2])
fit = cusp(y~MBN4[,1], alpha~MBN4[,4]+MBN4[,6]+MBN4[,3], beta~MBN4[,2])
fit = cusp(y~MBN4[,1], alpha~MBN4[,6]+MBN4[,3]-1, beta~MBN4[,2])
fit = cusp(y~MBN4[,1], alpha~MBN4[,6]-1, beta~MBN4[,2])

#fit = cusp(y~MBN4[,1], alpha~MBN4[,4]+MBN4[,6], beta~MBN4[,2]+MBN4[,3])
#fit = cusp(y~MBN4[,1], alpha~MBN4[,6], beta~MBN4[,2]+MBN4[,3])
#fit = cusp(y~MBN4[,1], alpha~MBN4[,4]+MBN4[,6]-1, beta~MBN4[,2]+MBN4[,3])

summary(fit)
summary(fit, logist=T)
plot(fit)
cusp3d(fit)

plot(fit$y, col="black", type="l")
points(fit$y, col="black")
points(fit$fitted.values, col="grey")
points(fit$fitted.values, col="grey", type="l")



### alpha diversity model
fit = cusp(y~MBN4[,1], alpha~MBN4[,4]+MBN4[,5]+MBN4[,6]+MBN4[,7], beta~MBN4[,2]+MBN4[,3])
fit = cusp(y~MBN4[,1], alpha~MBN4[,4]+MBN4[,5]+MBN4[,7], beta~MBN4[,2]+MBN4[,3])
fit = cusp(y~MBN4[,1], alpha~MBN4[,4]+MBN4[,5]-1, beta~MBN4[,2]) #### BEST
fit = cusp(y~MBN4[,1], alpha~MBN4[,4]+MBN4[,5]-1, beta~MBN4[,2]-1) #### BEST

#fit = cusp(y~MBN4[,1], alpha~MBN4[,5]-1, beta~MBN4[,3])
colnames(MBN4)

summary(fit)
summary(fit, logist=T)
plot(fit)
cusp3d(fit)

plot(fit$y, col="black", type="l")
points(fit$y, col="black")
points(fit$fitted.values, col="grey")
points(fit$fitted.values, col="grey", type="l")

plot(fit)

pred<-data.frame(fit$linear.predictors)
pred$cusp<-pred[,1] < -0.1 & pred[,2]< 0.5 |  pred[,1] < -0.2 & pred[,2]< 1
pred$col=as.factor(pred$cusp)
levels(pred$col)<-c("black","red")

plot(fit$y, col="black", type="l")
points(fit$y, col="green",pch=16)
points(fit$fitted.values, col=as.character(pred$col), type="p", pch=16)
points(fit$fitted.values, col="grey", type="l")




### alpha diversity model (reg)
fit = cusp(y~MBN4[,1], alpha~MBN4[,4]+MBN4[,5]+MBN4[,6]+MBN4[,7], beta~MBN4[,2]+MBN4[,3])
fit = cusp(y~MBN4[,1], alpha~MBN4[,4]+MBN4[,5]-1, beta~MBN4[,3])
fit = cusp(y~MBN4[,1], alpha~MBN4[,5]-1, beta~MBN4[,3])
colnames(MBN4)

summary(fit)
summary(fit, logist=T)
plot(fit)
cusp3d(fit)




###################################################################
## Differences between the two regimes <=2012 / > 2012
## kruskall wallis and heatmap of disappearing taxa
##################################################################

mean(alpha1[ENV$Year <= 2012 & ENV$Year >= 2009])
mean(alpha1[ENV$Year > 2012 & ENV$Year <= 2017])
alpha2<-data.frame(richness=alpha1, Year=ENV$Year, cat=rep(NA, length(alpha1)))
alpha2$cat[ENV$Year <= 2012 & ENV$Year >= 2009]<-"old"
alpha2$cat[ENV$Year > 2012 & ENV$Year <= 2017]<-"new"
kruskal.test(alpha2$richness, g=alpha2$cat)
### sign lower past 2013

pa_Genus<-relSpecies_all
pa_Genus[pa_Genus > 0] <- 1

pa_Genusold<-pa_Genus[,ENV$Year <= 2012 & ENV$Year >= 2009]
pa_Genusnew<-pa_Genus[,ENV$Year > 2012 & ENV$Year <= 2017]
pa_genusoldnew<-as.matrix(data.frame(old=rowSums(pa_Genusold), new=rowSums(pa_Genusnew)))
pa_genusoldnew[pa_genusoldnew > 0] <- 1

heatmap(pa_genusoldnew, scale="none",  Colv=NA, col=c("white","black"))

pa_Genus<-pa_Genus[order(rowSums(-pa_Genus)),]
require(vegan)
heatmap(pa_Genus, Colv=NA, scale="none", col=c("white","black"), Rowv=NA)


plot(as.Date(colnames(OTU)), alpha1, type="p", ylab="number of genera")
points(as.Date(colnames(OTU)), alpha1, type="l")
abline(v=as.Date(paste0(2006:2020,"-01-01")), col="grey")
mtext(at=as.Date(paste0(2006:2019,"-06-15")), 2006:2019)

##############################################################################
### Early warning Signals/Signs (EWS)
## tipping point warning signs? EWI/S (Early warning indicators/signals)
#https://www.sciencedirect.com/science/article/pii/S0272771414003199?casa_token=4a1My2pjlywAAAAA:ZRcdSAKycQBYZXyGNu5h-JpIVRJc7FS6mP5VVBZiutvXIPB55QZZtPAWhXKSXXjcFtEkQrvi-qo
# https://www.pnas.org/content/pnas/113/50/E8089.full.pdf

#### alpha
### increased autocorrelation
autocor<-data.frame(alpha=decalpha, Date=rep(NA, length(decalpha)), autocor=rep(NA, length(decalpha)))
for(i in 0:(length(na.omit(decalpha))-10)){
  autocor[(i+1),3]<-acf(na.omit(decalpha)[(1+i):(10+i)], lag.max=1)$acf[2]
}
plot(autocor$autocor, type="l")
abline(v=c(79,85))
#### Maybe


### increasing variance
varchange<-data.frame(alpha=alpha1, Date=ENV$Date, varchange=rep(NA, length(alpha1)))
for(i in 0:(length(alpha1)-10)){
  varchange[(i+1),3]<-var(alpha1[(1+i):(10+i)])
}
plot(varchange$varchange, type="l")
abline(v=c(79,85))
#### Yes!!!


#### NPP
### increased autocorrelation
autocor<-data.frame(alpha=noaa$NPP, Date=noaa$Date, autocor=rep(NA, length(noaa$NPP)))
for(i in 0:(length(noaa$NPP)-10)){
  autocor[(i+1),3]<-acf(noaa$NPP[(1+i):(10+i)], lag.max=1)$acf[2]
}
plot(autocor$autocor, type="l")
abline(v=c(129,138))
#### NO

### increasing variance
varchange<-data.frame(alpha=noaa$NPP, Date=noaa$Date, varchange=rep(NA, length(noaa$NPP)))
for(i in 0:(length(noaa$NPP)-10)){
  varchange[(i+1),3]<-var(noaa$NPP[(1+i):(10+i)])
}
plot(varchange$varchange, type="l")
abline(v=c(129,138))
#### YES

##########################NOX
################################
cpt.mean(na.omit(MBN2$NOx_int))
cpt.meanvar(na.omit(MBN2$NOx_int))
### increased autocorrelation
autocor<-data.frame(alpha=MBN2$NOx_int, Date=MBN2$Date, autocor=rep(NA, length(MBN2$NOx_int)))
for(i in 0:(length(MBN2$NOx_int)-10)){
  autocor[(i+1),3]<-acf(MBN2$NOx_int[(1+i):(10+i)], lag.max=1)$acf[2]
}
plot(autocor$autocor, type="l")
abline(v=c(79))
#### NO

### increasing variance
varchange<-data.frame(alpha=MBN2$NOx_int, Date=MBN2$Date, varchange=rep(NA, length(MBN2$NOx_int)))
for(i in 0:(length(MBN2$NOx_int)-10)){
  varchange[(i+1),3]<-var(MBN2$NOx_int[(1+i):(10+i)])
}
plot(varchange$varchange, type="l")
abline(v=c(79))
#### YES


################# EWS using a pacvkagae

require(earlywarnings)
bdstest_ews(na.omit(MBN2$NOx_int), ARMAoptim=T, embdim=3,epsilon=0.5,boots=200,logtransform=FALSE,interpolate=FALSE)
summary(bds.test(na.omit(MBN2$NOx_int)))
ddjnonparam_ews(na.omit(noaa$NPP))
decalphadf


