###################################################################################
#### Multivariate time series analyses ###########################################
##################################################################################
#### v1.6 Tobias Vonnahme 18.03.2024 #############################################
#################################################################################
#### Input: deseasonalized & regularized covariates (except AOI/NAOI stay as raw data) 
####        Community data, main groups, regularized and deseasonalized #########
#################################################################################
### 1. Data import and merging
### 2. Trend analyses
### 3. Autocorrelation tests
### 4. MARSS model for Nutrients using phycial drivers/covariates
### 5. Use MARSS model output to fill NAs in Nutrient data
### 6. MARSS model for phytoplankton community (main groups)
####################################################################################
### required packages: readxl, Kendall, MARSS
###################################################################################

setwd("C:/Users/torn/OneDrive - Grønlands Naturinstitut/Documents") #set working directory
#load("Renv_ts_261021.RData") #load environment if available
#load("Renv_ts_280921.RData") 

########################################################
#### 1) Read in the data and merge into one data frame ###
#######################################################

### covariates (simplified data to 2D, regularized, deseasonalized)
covdat<-read.csv("MBN_2D_data_deseas.csv")
#covdat<-read.csv("MBN_2D_data_regul.csv")
### Optional datsets:
### Before desesonalisation (but regularized): "MBN_2D_data_regul.csv"
### Raw data (before regularisation): "MBN_2D_data.csv"
if(colnames(covdat)[1] == "X"){covdat<-covdat[,-1]} #remove the first column that is sometimes added via csv conversion in R

### Zooplankton data
zoodat<-read.csv("Zooplankton_deseasonalized_MARSS.csv")
if(colnames(zoodat)[1] == "X"){zoodat<-zoodat[,-1]} #remove the first column that is sometimes added via csv conversion in R

### community data (subset of main groups, regularized, deseasonalized)
comdat<-read.csv("Phytoplankton_community_groups_deseason_v2.csv")
#comdat<-read.csv("Phytoplankton_community_groups_original.csv")
### "Phytoplankton_community_groups_original.csv"
### "Phytoplankton_community_groups_deseason.csv"
### "Phytoplankton_community_groups_clr.csv"
if(colnames(comdat)[1] == "X"){comdat<-comdat[,-1]} #remove the first column that is sometimes added via csv conversion in R

### Arctic (AO) and North Atlantic (NAO) Oscillation (Raw data)
require(readxl)
AO<-as.data.frame(read_xlsx("AO_Index.xlsx"))
colnames(AO)<-c("Date","AO","NAO")

#### merge community and covariate data frames
## check dimensions
dim(comdat)
dim(zoodat)
dim(covdat)
dim(AO)
## check if dates are the same
identical(covdat$Date, comdat$Date)
comdat$Date<-covdat$Date
zoodat$Date<-covdat$Date[1:nrow(zoodat)]
AO$Date<-covdat$Date
identical(comdat$Date, as.character(as.Date(AO$Date)))
## make dates identifcal (its a regularized datasets so the difference is negligible)
require(lubridate)
str(comdat)
str(covdat)
str(zoodat)
str(AO)
### merge
Alldat<-merge(merge(covdat, comdat, by="Date", all=TRUE),AO, by="Date", all=TRUE)
Alldat<-merge(merge(merge(covdat, comdat, by="Date", all=TRUE),AO, by="Date", all=TRUE), zoodat, by="Date", all=TRUE)
View(Alldat)
### in 2011 there was no NOx measurements (replace interpolated values with NA for now)
### Later use the MARSS model to model the values in 2011
Alldat$NOx_int[substr(Alldat$Date, start=1, stop=4)=="2011"]<-NA



###################################
#### 2) Trend analyses
#### Mann Kendall test testing for monotonous non-linear trends
#####################################

require(Kendall)

p<-data.frame(variable=rep(NA, ncol(Alldat)-1), p=rep(NA, ncol(Alldat)-1), tau=rep(NA, ncol(Alldat)-1),
              padj=rep(NA, ncol(Alldat)-1), sign=rep(NA, ncol(Alldat)-1))

for (i in 2:ncol(Alldat)){
  print(colnames(Alldat[i]))
  print(MannKendall(Alldat[,i]))
  p[i,1]<-colnames(Alldat[i])
  p[i,2]<-MannKendall(Alldat[,i])$sl
  p[i,3]<-MannKendall(Alldat[,i])$tau
}

p<-p[-1,]
p<-p[-c(15,22,26),]
p$padj<-p.adjust(p$p, method = "fdr")

for (i in 1:nrow(p)){
  if (p$padj[i]<=0.001){p$sign[i] <- "***"}
  else if (p$padj[i]<=0.01){p$sign[i] <- "**"}
  else if (p$padj[i]<=0.05){p$sign[i] <- "*"}
  else if (p$p[i]<=0.05){p$sign[i] <- "(*)"}
  else {p$sign[i] <- "ns"}
}

p

#### Results (p<0.05)
# decrease: Temp5, N, Pseudo_niotzschia  Fragi, other diat
# increase: Sal150, Si, P, Thalassionema, Phaeocystis, Cilio NAO

#################################################################
### 3) Autocorrelation tests
##################################################################

# BUT looking at the plots its not a monotonous trend, but stochastic variation
#### Check for autocorrelation (lines outside the blue dotted lines show potentially significant autocorrelation at with a time leg of x)
require(forecast)
par(mfrow=c(3,4))
for (i in 2:13){
  acf(na.omit(Alldat[,i]))}
for (i in 14:25){
  acf(na.omit(Alldat[,i]))}
par(mfrow=c(1,3))
acf(na.omit(Alldat[,18]))
acf(na.omit(Alldat[,12]))
acf(na.omit(Alldat[,21]))
#### check the orders of a potential ARIMA model
for (i in 2:ncol(Alldat)){
  print(auto.arima(Alldat[,i]))}


############################################################################
#### 4 a) Multivariate state space modelling ###################################
#### MARSS model (Nutrients as variable, physical data and Chl as covariates)
############################################################################

AlldatOld<-Alldat
Alldat<-AlldatOld 

Alldat$year<-as.numeric(substr(Alldat[,1], start=1, stop=4)) #add column with the year for later subsetting
full2<-as.matrix(Alldat[,2:ncol(Alldat)]) #convert dataframe to matrix format
years<-full2[,"year"]>=2005 & full2[,"year"]<=2020 #option to subse tthe dataset to certain years
#years<-full2[,"year"]>=2005 & full2[,"year"]<=2018 #option to subse tthe dataset to certain years. here for zoolankton


colnames(full2) # check the order of the columns
dat <- t(full2[years,colnames(full2)[c(8,9,10)]]) #subset data of the nutrients
covariates <- t(full2[years,colnames(full2)[c(2,6,7,22)]]) #simplest model of covariates
covariates <- t(full2[years,colnames(full2)[c(2,5,7,21)]]) #simplest model of covariates


dim(dat) 
#dat<-dat[,2:183] ## subset observation data to allow a lagged response to covariates up to 1 month
#row.names(covariates)
#covariates<-as.matrix(t(data.frame(StratInd=covariates[1,2:183],  T150m<-covariates[2,1:182], 
#                                   Chl<-covariates[3,1:182], AO=covariates[4,1:182]))) #define a time lag in the differen covariates
#### covariates at t-x affect observations at t (lagged response)

# z-score the response variables
the.mean <- apply(dat,1,mean,na.rm=TRUE)
the.sigma <- sqrt(apply(dat,1,var,na.rm=TRUE))
dat <- (dat-the.mean)*(1/the.sigma)

## ----msscov-z-score-covar-data-----------------------------------------------------------------
the.mean <- apply(covariates,1,mean,na.rm=TRUE)
the.sigma <- sqrt(apply(covariates,1,var,na.rm=TRUE))
covariates <- (covariates-the.mean)*(1/the.sigma)

## ----msscov-plank-plot, fig=TRUE, echo=FALSE, fig.cap='(ref:msscov-plank-dat)', warning=FALSE----
LWA <- ts(cbind(t(dat), t(covariates)), start=c(2005,11), freq=12) #bind the observations and covariates
plot(LWA, main="", yax.flip=TRUE) # plot the time series

### Model testing
require(MARSS)
A <- U <- x0 <- "zero"       #### No trend/drift 
c = covariates         #### 
C <- "unequal"#"unconstrained"   #### effects of different process errors on states ,  "equal" gives lower AIC, but doesnt make sense here!
Q <- "diagonal and unequal"#"unconstrained"  ### Test also "diagonal and unequal" "diagonal and equal" "equalvarcov" "identity" "unconstrained"
Z <- "identity"     ### only 1 time series (no different states)
B <- "diagonal and equal" ### state at t-1 affects state at t equally among both states, tried "diagonal and unequal" "identity
R <- "equalvarcov" #matrix(list(0.05,0,0,0.1),2,2) ### same effect of different observation errors (same sampling + similar measurement methods) diagonal and unequal before
tinitx <- 1              ### Needed for estimating B
y <- dat 

model.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R, c=c, C=C, x0=x0, tinitx=tinitx)
kem.1<-MARSS(y, model = model.list) #AIC: 1055 
#### First model estimate

U2 <- "unconstrained"
model.list <- list(B = B, U = U2, Q = Q, Z = Z, A = A, R = R, c=c, C=C, x0=x0, tinitx=tinitx)
kem.2<-MARSS(y, model = model.list) #AIC: 1063  -> worse fit
#### -> worse fit and low values for U estimate (discard trend (U)) -> Discard

C2<-"unconstrained"
model.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R, c=c, C=C2, x0=x0, tinitx=tinitx)
kem.3<-MARSS(y, model = model.list) #AIC: 1057 -> equal fit
#### -> same fit, unconstrained has a higher DF than unequal, also unequal effects of different covariates makes more sense
### equal is another option but makes no sense and has higher AIC  -> Discard

Q2<-"unconstrained"
model.list <- list(B = B, U = U, Q = Q2, Z = Z, A = A, R = R, c=c, C=C, x0=x0, tinitx=tinitx)
kem.4<-MARSS(y, model = model.list) #AIC: 1042 -> Better
#### Lower AIC and higher values for Q, meaninful with unequeal covariates and variance 
#### process errors have different variance for different taxa (var, diagonal) and year-to year changes covary/ process errors are dependent (inherent in compositionality data!!) -> Keep

B2 <- "diagonal and unequal"
model.list <- list(B = B2, U = U, Q = Q2, Z = Z, A = A, R = R, c=c, C=C, x0=x0, tinitx=tinitx)
kem.5<-MARSS(y, model = model.list) #AIC: 1046 -> Worse
#### state (taxa proportion) at t-1 affects taxa at t equally (diagnoal and equal) among all taxa (Makes sense here) -> Discard

R2 <- "diagonal and unequal"
model.list <- list(B = B, U = U, Q = Q2, Z = Z, A = A, R = R2, c=c, C=C, x0=x0, tinitx=tinitx)
kem.6<-MARSS(y, model = model.list) #AIC: 1044 -> better but: "Warning: the  R.(NOx_int,NOx_int)  parameter value has not converged."
#### observation errors are different for the different taxa, but year-to year changes in the obs error are independent (diagonal) -> Keep

### define the best model fit
kem.t<-kem.6
kem.t
### calculate confidence intervals (if CI do not overlap with 0 thye are considered significant)
MARSSparamCIs(kem.t)

### diagnostic plots
plot(kem.t)

#######################################
######## 5) use the modelled PO4 and NOx data to fill NA gaps in dataset
##########################################

dat <- t(full2[years,colnames(full2)[c(8,9,10)]]) #subset observation data of the nutrients
dat<-dat[,2:183]
# z-score method for the response variables (we need mean and sigma for inverse z-scaling)
the.mean <- apply(dat,1,mean,na.rm=TRUE)
the.sigma <- sqrt(apply(dat,1,var,na.rm=TRUE))
#zdat <- (dat-the.mean)*(1/the.sigma) #z score function
### --> reverse: dat <- the.mean + the.sigma * zdat

#### extract the modelled NOx value for 2011
rs<-residuals(kem.t)
rs[rs$.rownames=="NOx_int" & is.na(rs$value),]
zNOx2011<-rs$.fitted[rs$.rownames=="NOx_int" & is.na(rs$value)]
NOx2011<-c(the.mean[2] + the.sigma[2] * zNOx2011) ## reverse the z-standardization
### reverse z standardisation

### extract modelled PO4 values for 2020
rs[rs$.rownames=="PO4_int" & is.na(rs$value),]
zpo42020<-rs$.fitted[rs$.rownames=="PO4_int" & is.na(rs$value)]
po42020<-c(the.mean[3] + the.sigma[3] * zpo42020)# reverse z standardisation

#### Add modelled values to Alldat data frame
Alldat$NOx_int[Alldat$year == 2011]<-NOx2011
Alldat$PO4_int[(nrow(Alldat)-length(po42020)+1):nrow(Alldat)]<-po42020

length(Alldat$PO4_int[(nrow(Alldat)-length(po42020)+1):nrow(Alldat)])
length(po42020)

#######################################################################
#### 6) MARSS model for communities
##########################################################
##########################################################

#### Original without grazers
full<-as.matrix(na.omit(Alldat[,c(2:(ncol(Alldat)-4),ncol(Alldat) )]))
#### With grazers
full<-as.matrix(na.omit(Alldat[,2:ncol(Alldat)]))

### prepare dataset with observations (y) and covariates (c) in matrix format
### Matrix conversion for community
years <- full[,"year"]>=2005 & full[,"year"]<=2020

########################################
#### Best model fit
colnames(full)
dat <- t(full[years,colnames(full)[c(17,11,20)]]) #"most dominant"
#dat <- t(full[years,colnames(full)[c(17,11, 14, 16, 18, 19, 20, 12, 13)]])  #All algae
#dat <- t(full[years,colnames(full)[c(17,11, 14, 16, 18, 19, 20)]]) # deseasonalized groups, subset observation data of the most abundant algae groups (Chaetoceros, Thalassiosira, Phaeocystis)
covariates <- t(full[years,colnames(full)[c(2,8,9,10, 5, 21)]]) # no grazer
covariates <- t(full[years,colnames(full)[c(2,8,9,10, 5, 21, 23:25)]]) # GRazer


dim(dat)
dat<-dat[,6:ncol(dat)] ## no Grazer ## subset observation data to allow a lagged response to covariates up to 6 months 183
dat<-dat[,6:ncol(dat)] ## Grazer ## subset observation data to allow a lagged response to covariates up to 6 months 183
#covariates<-as.matrix(t(data.frame(StratInd=covariates[1,6:169], Si=covariates[2,4:167], NOx=covariates[3,4:167], PO4=covariates[4,4:167], 
#                                   S150m=covariates[5,1:164], AO=covariates[6,6:169]))) # 
#covariates<-as.matrix(t(data.frame(StratInd=covariates[1,6:169], Si=covariates[2,4:167], NOx=covariates[3,4:167], PO4=covariates[4,4:167], 
#                                   S150m=covariates[5,1:164], AO=covariates[6,6:169],
#                                   SiP=covariates[2,4:167]/covariates[4,4:167], SiN=covariates[2,4:167]/covariates[3,4:167] ,PN=covariates[4,4:167]/covariates[3,4:167]))) # include ratios
covariates<-as.matrix(t(data.frame(sal5=covariates[3,6:169], Si=covariates[2,4:167], NOx=covariates[3,4:167], PO4=covariates[4,4:167], 
                                   S150m=covariates[5,1:164], AO=covariates[6,6:169]))) # 
covariates<-as.matrix(t(data.frame(sal5=covariates[3,6:152], Si=covariates[2,4:150], NOx=covariates[3,4:150], PO4=covariates[4,4:150], 
                                   S150m=covariates[5,1:147], AO=covariates[6,6:152],
                                   Cir=covariates[7,2:148], Biv=covariates[8,2:148], Zoo=covariates[9,2:148]))) #

# z-score the response variables
the.mean <- apply(dat,1,mean,na.rm=TRUE)
the.sigma <- sqrt(apply(dat,1,var,na.rm=TRUE))
dat <- (dat-the.mean)*(1/the.sigma)
## ----msscov-z-score-covar-data-----------------------------------------------------------------
the.mean <- apply(covariates,1,mean,na.rm=TRUE)
the.sigma <- sqrt(apply(covariates,1,var,na.rm=TRUE))
covariates <- (covariates-the.mean)*(1/the.sigma)
## ----msscov-plank-plot, fig=TRUE, echo=FALSE, fig.cap='(ref:msscov-plank-dat)', warning=FALSE----

dim(dat)
dim(covariates)

LWA <- ts(cbind(t(dat), t(covariates)), start=c(2005,11), freq=12) #bind the observations and covariates
plot(LWA, main="", yax.flip=TRUE) # plot the time series

### Model testing
require(MARSS)
A <- U <- x0 <- "zero"       #### No trend/drift ### demeaned data (z scaled)
c = covariates         #### 
C <- "unequal" #"unconstrained"   #### effects of different process errors on states ,  "equal" gives lower AIC, but doesnt make sense here!
Q <- "diagonal and unequal"#"unconstrained"  ### Test also "diagonal and unequal" "diagonal and equal" "equalvarcov" "identity" "unconstrained"
Z <- "identity"     ### only 1 time series (no different states)
B <- "diagonal and equal" ### state at t-1 affects state at t equally among both states, tried "diagonal and unequal" "identity
R <- "equalvarcov" #matrix(list(0.05,0,0,0.1),2,2) ### same effect of different observation errors (same sampling + similar measurement methods) diagonal and unequal before
tinitx <- 1              ### Needed for estimating B
y <- dat 

####################################

######################

model.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R, c=c, C=C, x0=x0, tinitx=tinitx)
kem.1<-MARSS(y, model = model.list) #
kem.1$AIC #AIC: 1219  #1119.6
#### First model estimate

U2 <- "unconstrained" ### allow drift
model.list <- list(B = B, U = U2, Q = Q, Z = Z, A = A, R = R, c=c, C=C, x0=x0, tinitx=tinitx)
kem.2<-MARSS(y, model = model.list) #
kem.2$AIC #AIC: 1225 #1125.5  -> worse fit
#### -> worse fit and low values for U estimate (discard trend (U)) Also meaningless for z standardized data (somewhere said)
###  -> Discard

C2<-"unconstrained"
model.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R, c=c, C=C2, x0=x0, tinitx=tinitx)
kem.3<-MARSS(y, model = model.list) #
kem.3$AIC  #AIC: 3215 # 1119.6-> equal fit
#### -> same fit, unconstrained has a higher DF than unequal, also unequal effects of different covariates makes more sense
### equal is another option but makes no sense and has higher AIC  -> Discard

Q2<-"unconstrained"
model.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R, c=c, C=C, x0=x0, tinitx=tinitx)
kem.4<-MARSS(y, model = model.list) #
kem.4$AIC #AIC: 3117 #1113.1 -> Better
#### Lower AIC and higher values for Q, meaninful with unequeal covariates and variance 
#### process errors have different variance for different taxa (var, diagonal) and year-to year changes covary/ process errors are dependent (inherent in compositionality data!!) -> Keep

B2 <- "unconstrained" ## allow interspecies interactions  -> rejected: no significant interactions
B2 <- "diagonal and unequal" ## only intraspecies interactions
model.list <- list(B = B2, U = U, Q = Q2, Z = Z, A = A, R = R, c=c, C=C, x0=x0, tinitx=tinitx)
kem.5<-MARSS(y, model = model.list) #
kem.5$AIC  #AIC: 3098 # 1106.4-> Better
#### state (taxa proportion) at t-1 affects taxa at t equally (diagnoal and equal) among all taxa (Makes sense here) -> Discard
#### state (taxa proportion) at t-1 affects taxa at t differently (diagnoal and unequal)  -> Keep

R2 <- "diagonal and unequal"
model.list <- list(B = B2, U = U, Q = Q2, Z = Z, A = A, R = R2, c=c, C=C, x0=x0, tinitx=tinitx)
kem.6<-MARSS(y, model = model.list) #
kem.6$AIC  #AIC: 2974  # 1112.0 -> better but: "Warning: the  R.(NOx_int,NOx_int)  parameter value has not converged."
#### observation errors are different for the different taxa, but year-to year changes in the obs error are independent (diagonal) -> Keep
### worse when Grazers are included -> Discard

### define the best model fit
#kem.6 without grazer
#kem.5 with grazer
kem.t<-kem.5
#kem.t
### calculate confidence intervals (if CI do not overlap with 0 thye are considered significant)
MARSSparamCIs(kem.t)

### diagnostic plots
str(kem.t)
plot(kem.t)

kem.t$coef[abs(kem.t$coef)>=0.09]

kem.tend<-kem.t
par(mfrow=c(2,2))
rs<-residuals(kem.tend)
plot(rs$value[rs$.rownames=="Chaetoceros"], type="p", pch=15, col="grey", main="Chaeotoceros")
points(rs$value[rs$.rownames=="Chaetoceros"], type="l", col="grey", lwd=0.1)
points(rs$.fitted[rs$.rownames=="Chaetoceros"], type="l", col="red", lty=2, lwd=2)

plot(rs$value[rs$.rownames=="Phaeocystis"], type="p", pch=15, col="grey", main="Phaeocystis")
points(rs$value[rs$.rownames=="Phaeocystis"], type="l", col="grey", lwd=0.1)
points(rs$.fitted[rs$.rownames=="Phaeocystis"], type="l", col="red", lty=2, lwd=2)

plot(rs$value[rs$.rownames=="Thalassiosira"], type="p", pch=15, col="grey", main="Thalassiosira")
points(rs$value[rs$.rownames=="Thalassiosira"], type="l", col="grey", lwd=0.1)
points(rs$.fitted[rs$.rownames=="Thalassiosira"], type="l", col="red", lty=2, lwd=2)

### Model diagnostics
data.frame(Model=c("observation error", "process error with RW", "process error around mean"),
           AICc=round(c(kem.1$AICc,
                        kem.2$AICc,
                        kem.3$AICc),1))

require(forecast)
fr <- forecast(kem.4, newdata=list(y=NULL, c=kem.4$call$model$c, d=NULL), h=100)
plot(fr)

str(kem.3)

i <- 3; #process error around mean
j <- 1; #First state
par(mfrow=c(2,1))
modn <- paste("kem",i,sep=".")
resid.j <- MARSSresiduals(get(modn), type="tt1")$model.residuals[j,]
plot.ts(resid.j, 
        ylab="Residual")
abline(h=0, lty="dashed")
acf(resid.j, na.action=na.pass)


kem.3
summary(kem.3)
str(kem.3)
plot(kem.3)
str(kem.3)


kem.3$marss$data
kem.3$model$data
kem.3$states
residuals(kem.3)



pr<-predict(kem.6)
plot(pr$pred$y[1:148], type="p", pch=15, col="grey")
points(pr$pred$y[1:148], type="l", col="grey", lwd=0.1)
points(pr$pred$estimate[1:148], type="l", col="red", lty=2, lwd=2)

fitted(kem.3)
logLik(kem.2)
AIC(kem.1)
coef(kem.3)

tsSmooth(kem.3)

pr<-predict(kem.6)
plot(pr$pred$y[184:366], type="p", pch=15, col="grey")
points(pr$pred$y[184:366], type="l", col="grey", lwd=0.1)
points(pr$pred$estimate[184:366], type="l", col="red", lty=2, lwd=2)


