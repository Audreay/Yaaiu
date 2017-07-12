rm(list=objects())
library(dplyr)
library(grpreg)
library(grpSLOPE)
library(RColorBrewer)
library(ranger)
library(stackeR)
library(gamwave)
library(wavethresh)
library(ncvreg)
library(yarrr)

##############import mean forecasts from the mean model
data.mean.forecast<-readRDS("/Users/yannig/Documents/GAMwavelet/Hybrid_model_update2/Results/data.mean.forecast.RDS")

##############import individual load data
MtestHJF<-readRDS("/Users/yannig/Documents/GAMwavelet/Data/irish_data_corrected_italia.RDS")

##############import tarif options
load("/Users/yannig/Documents/Open_data/Irish/Data_CER_clean//tarif.RData")
########correction of false "b" in tarif2
tarif2$ResidentialTariffallocation[which(tarif2$ResidentialTariffallocation=="b")]<-"B"
tarif2$ResidentialTariffallocation<-factor(tarif2$ResidentialTariffallocation)
##############import survey data
survey<-read.table("/Users/yannig/Documents/Open_data/Irish/Data_CER_clean/residential_survey_data_summary.txt",sep='\t',header=T)

#######building of residential database consumption
resI<-tarif2$ID[which((tarif2$Code==1))]   ##extract residential load data
sel<-names(MtestHJF)%in%paste("I",resI,sep="")
residential<-MtestHJF[,sel]


########join tarif2, survey and the mean consumption of each customer
m<-apply(residential,2,mean)
Mat<-data.frame(resI,m)
names(Mat)<-c("ID","mean")

data_anova<-left_join(Mat, survey, by='ID')
data_anova<-left_join(data_anova, tarif2, by='ID')

Mat$ID<-paste0("I", Mat$ID)

#######correction of the data set
data_anova<-data_anova[-which(data_anova$SOCIALCLASS==0),]
data_anova$SOCIALCLASS<-droplevels(data_anova$SOCIALCLASS)
data_anova$ResidentialTariffallocation<-droplevels(data_anova$ResidentialTariffallocation)
data_anova$ResidentialStimulusallocation<-droplevels(data_anova$ResidentialStimulusallocation)
data_anova$BUILT.YEAR2<-as.numeric(data_anova$BUILT.YEAR)
data_anova$HOME.APPLIANCE..White.goods.2<-as.numeric(data_anova$HOME.APPLIANCE..White.goods.)
sel0<-which(data_anova$BUILT.YEAR2<1000)
data_anova$BUILT.YEAR2[sel0]<-mean(data_anova$BUILT.YEAR2[-sel0])
data_anova$ID<-paste0("I", data_anova$ID)
data_anova<-na.omit(data_anova)

####################################################################
#################subsampling
####################################################################
####we build the data fra me of 50 aggregated customers selected at random
set.seed(seed=110)
s<-sample(data_anova$ID, 50)
m<-residential[, which(colnames(residential)%in%s)]


mean.forecast<-data.mean.forecast[data.mean.forecast$ID%in%s,]$mean.forecast.grpreg
m.center<-m-matrix(rep(mean.forecast, each=nrow(m)), byrow=F, ncol=ncol(m))
#m.center<-m

##control if the data are demeaned
mean(colMeans(m))
mean(colMeans(m.center))


data_agg<-data.frame(MtestHJF$Temp, MtestHJF$DayType, MtestHJF$Instant, MtestHJF$Date)
names(data_agg)<-c("Temp","DayType","Instant", "Date")
data_agg$Load.center<-rowMeans(m.center)
data_agg$Load<-rowMeans(m)

######intercept
Cte<-mean(mean.forecast)
#Cte<-0

####we work on the first 1000 temporal observation 
N0<-1000
N1<-1000
data_agg0<-data_agg[1:N0,]
data_agg1<-data_agg[(N0+1):(N0+N1),]

################################################################################
################gamwave addtive model with only irregular composant (wavelets)
################################################################################
####create the design matrix: temperature and effects of the instants/daytype, for estimation (0) and prediction (1)
M_day0<-createIndicators(data_agg0, "DayType", prefix = "DT")
M_day1<-createIndicators(data_agg1, "DayType", prefix = "DT")

names(M_day0)

A0<-data_agg0$Temp
A1<-data_agg1$Temp
for(i in c(1:7))
{
  A0<-cbind(A0, M_day0[,6+i]*data_agg0$Instant)
  A1<-cbind(A1, M_day1[,6+i]*data_agg1$Instant)
}


####estimation of the model using gamwave
dim(A0)
plot(data_agg0$Load, type='l')
gmw<-GAMwv(y=data_agg0$Load.center, A0, penalty="SCAD", Hybrid=F, conv.thresh=0.05, max.iter = 10, trace=T, nfolds = 10)

plot(gmw$relative.distance) ###plot of the relative distance between the estimated effect of 2 succesive iterations

######TO DO: improve the code for the design of wevelet basis generation, e.g. when the covariate is sparse
######due to the indicators variable, maybe there is a gain also in ncvreg, e.g. not use CV at each iteration?
######clearly the CV computation is the most time consuming operation, but changeing nfolds seems to have an impact on the relative distance
######between each iteration
######use the formula structure
gmw$fitted.center<-gmw$fitted
gmw$fitted<-gmw$fitted+Cte

gmw$forecast.center<-predict(gmw, newdata=A1)$forecast
gmw$forecast<-gmw$forecast.center+Cte

######TO DO: improve the predic function such that it computes the A matrix automatically from the orginal data base

###############################################################################################
###################GAM with regular components (splines)
###############################################################################################
gsp<-gam(Load.center~ s(Instant, bs="cr", k=10, by=DayType)+s(Temp, k=10, bs='cr'), data=data_agg0, method="REML")

gsp$fitted.center<-gsp$fitted
gsp$fitted<-gsp$fitted.center+Cte

gsp$forecast.center<-predict(gsp, newdata=data_agg1)
gsp$forecast<-gsp$forecast.center+Cte


###############################################################################################
###################hybrid GAM
###############################################################################################
#the irregular part
gam_backfit_wavelet<-GAMwv(y=data_agg0$Load.center, A0[,-1], penalty="SCAD", Hybrid=F, conv.thresh=0.05, max.iter = 10, trace=T, nfolds = 10)
yres<-data_agg0$Load.center-gam_backfit_wavelet$fitted
#regular part
data_agg0$yres<-yres
gam_backfit_spline<-bam(yres ~ s(Temp, bs='cr'), data=data_agg0, method="fREML")

hybrid<-list()
hybrid$gsp<-gam_backfit_spline
hybrid$gmw<-gam_backfit_wavelet

hybrid$fitted.center<-hybrid$gsp$fitted+hybrid$gmw$fitted
hybrid$fitted<-hybrid$fitted.center+Cte

hybrid$gsp$forecast.center<-predict(hybrid$gmw, newdata=A1[,-1])$forecast
hybrid$gsp$forecast<-hybrid$gsp$forecast.center+Cte

hybrid$gmw$forecast.center<-predict(hybrid$gsp, newdata=data_agg1)
hybrid$gmw$forecast<-hybrid$gmw$forecast.center+Cte

hybrid$forecast.center<-hybrid$gsp$forecast.center+hybrid$gmw$forecast.center
hybrid$forecast<-hybrid$forecast.center+Cte



###########################################################################
####################plots
###########################################################################

##adjusment
col<-piratepal("basel", length.out = 4)

a<-1
b<-7*48
plot(data_agg0$Load[a:b], type='l')
lines(gmw$fitted[a:b], col=col[1])
lines(gsp$fitted[a:b], col=col[2])
lines(hybrid$fitted[a:b], col=col[3])
legend("top", legend=c("gmw","gsp","hybrid"), col=col, bty='n', lty=1, ncol=3)

##forecasts
a<-1
b<-7*48
plot(data_agg1$Load[a:b], type='l')
lines(gmw$forecast[a:b], col=col[1])
lines(gsp$forecast[a:b], col=col[2])
lines(hybrid$forecast[a:b], col=col[3])
legend("top", legend=c("gmw","gsp", "hybrid"), col=col, bty='n', lty=1, ncol=3)

#####fitted adjusment measures
mean((data_agg0$Load-gmw$fitted)^2)
mean((data_agg0$Load-gsp$fitted)^2)
mean((data_agg0$Load-hybrid$fitted)^2)

#####forecast adjusment measures
mean((data_agg1$Load-gmw$forecast)^2)
mean((data_agg1$Load-gsp$forecast)^2)
mean((data_agg1$Load-hybrid$forecast)^2)


#####export the results in a .RDS file
results<-list()
results$fittedMSE<-c(mean((data_agg0$Load-gmw$fitted)^2), mean((data_agg0$Load-gsp$fitted)^2),mean((data_agg0$Load-hybrid$fitted)^2))
names(results$fittedMSE)<-c("gmw", "gsp", "hybrid")
results$forecastMSE<-c(mean((data_agg1$Load-gmw$forecast)^2), mean((data_agg1$Load-gsp$forecast)^2), mean((data_agg1$Load-hybrid$forecast)^2))
names(results$forecastMSE)<-c("gmw", "gsp", "hybrid")
results$gsp<-gsp
results$gmw<-gmw
results$hybrid<-hybrid
results$data0<-data_agg0
results$data1<-data_agg1

saveRDS(results, file="/Users/yannig/Documents/GAMwavelet/Results/hybrid2.RDS")








##############################################################################
#############introducing other covariate and selection (not finished)
##############################################################################
###smooth temperature
expSmooth=function(x,alpha)
{
  xsmooth=x
  for(i in c(2:length(x)))
  {
    xsmooth[i]<-(1-alpha)*xsmooth[i-1]+alpha*x[i]
  }
  return(xsmooth)
}

data_agg$Temp_s05<-expSmooth(data_agg$Temp, 0.5)
data_agg$Temp_s09<-expSmooth(data_agg$Temp, 0.9)

###min/max daily temperature
names(data_agg)
head(data_agg$Date)
DateDay<-format(data_agg$Date, "%Y/%m/%d")
data_agg$DateDay<-DateDay
TempMin<-tapply(data_agg$Temp, INDEX=DateDay, FUN=min)
TempMax<-tapply(data_agg$Temp, INDEX=DateDay, FUN=max)

DataTemp<-data.frame(TempMin, TempMax)
DataTemp$DateDay<-names(TempMax)

data_agg<-left_join(data_agg, DataTemp, by="DateDay")

##################################################################################################################
#################################variable selection for hybrid models: not finished!!!!!!
##################################################################################################################
N<-1000
data_agg_deb<-data_agg[1:N,]

###1st step: gam model
nknots<-10
fitGAM<-gam(Load~ s(Instant, bs="cr", k=nknots, by=DayType)+s(Temp, bs='cr', k=nknots)+s(Temp_s05, bs='cr', k=nknots)+
              s(Temp_s09, bs='cr', k=nknots)+s(TempMin, k=nknots, bs='cr')+s(TempMax, k=nknots, bs='cr'), data=data_agg_deb)
summary(fitGAM)
fitGAM.terms=predict(fitGAM, type="lpmatrix")
dim(fitGAM.terms)

d<-(ncol(fitGAM.terms)-1)/(nknots-1)

grp<-NULL
for(i in c(1:d))
{
  grp<-c(grp, rep(i, nknots-1))
}
  

mod<- grpreg(X=fitGAM.terms[,-1], y=data_agg_deb$Load, group=grp, penalty = "grLasso")
cvmod <- cv.grpreg(X=fitGAM.terms[,-1], y=data_agg_deb$Load, group = grp, penalty = "grLasso")
coefgrLasso <- coef(cvmod) ## Beta at minimum CVE
coefgrLasso[!coefgrLasso==0 & names(coefgrLasso) != "(Intercept)"]

plot(mod, label=T)














