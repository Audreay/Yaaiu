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

##############import individual load data
#load("/Users/yannig/Documents/Open_data/Irish/Data_CER_clean/MtestHJF.RData")
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
####we build 1 data frame for modelling the mean part of the process, not depend on time
#set.seed(seed=110)
#s<-sample(data_anova$ID, 50) ##subsampling 50 customers
s<-data_anova$ID  ##all the population

m<-residential[, which(colnames(residential)%in%s)]
ID<-rep(colnames(m),each=nrow(m))
length(ID)
Temp<-rep(MtestHJF$Temp,length(s))
DayType<-rep(MtestHJF$DayType,length(s))
Instant<-rep(MtestHJF$Instant,length(s))

m_deply<-as.vector(as.matrix(m))
Mat<-data.frame(ID, m_deply, Temp, DayType, Instant)
dim(Mat)
names(Mat)<-c("ID","Load","Temp","DayType","Instant")
sel<-data_anova$ID%in%Mat$ID
data_sub_linear<-data_anova[sel,]
dim(data_sub_linear)
summary(data_anova)
dim(data_sub_linear)
dim(data_anova)

#saveRDS(data_sub_linear, "/Users/yannig/Documents/GAMwavelet/Hybrid_model_update2/Data/data_mean.RDS")

####################################################################
#################models of the individual mean (independant of time) on subsamples
####################################################################
intercept<-FALSE
####design matrix building
if(intercept==FALSE)
{
  mod<-lm(mean~SOCIALCLASS+OWNERSHIP+BUILT.YEAR2+HEAT.HOME+HEAT.WATER+WINDOWS.doubleglazed+HOME.APPLIANCE..White.goods.
          +ResidentialStimulusallocation-1, data=data_sub_linear)
  explX <- model.matrix(mod)
  summary(mod)
}
if(intercept==TRUE)
{
  mod<-lm(mean~SOCIALCLASS+OWNERSHIP+BUILT.YEAR2+HEAT.HOME+HEAT.WATER+WINDOWS.doubleglazed+HOME.APPLIANCE..White.goods.
          +ResidentialStimulusallocation, data=data_sub_linear)
  explX <- model.matrix(mod)[,-1]
}

########building groups for the group lasso selection, some modalities of factor could have 0 observations due to subsampling
########in that case they are deleted

covnames<-labels(terms(mod))
grp<-NULL
for(i in c(1:length(covnames)))
{
  pos<-grep(covnames[i], colnames(explX))
  grp<-c(grp, rep(i, length(pos)))
}

y <-data_sub_linear$mean

length(y)
dim(explX)

mod<- grpreg(X=explX, y=y, group=grp, penalty = "grLasso")
cvmod <- cv.grpreg(X=explX, y=y, group = grp, penalty = "grLasso")
coefgrLasso <- coef(cvmod) ## Beta at minimum CVE

coefgrLasso[!coefgrLasso==0 & names(coefgrLasso) != "(Intercept)"]
plot(mod, label=T)


########to do: adding groupSLOPE, testing the model on the overall population

######fitted mmodel on the mean
mean.forecast<-predict(cvmod, X=explX, lambda=cvmod$lambda.min)

plot(y, mean.forecast, pch=20)
lines(y, y, col='red')

######save of the grpreg model
saveRDS(cvmod, file= "/Users/yannig/Documents/GAMwavelet/Hybrid_model_update2/Results/mean_grpreg.RDS")

######save of the mean correction
data.mean.forecast<-data.frame(data_sub_linear$ID, data_sub_linear$mean, mean.forecast)
names(data.mean.forecast)<-c("ID", "mean","mean.forecast.grpreg")
head(data.mean.forecast)

saveRDS(data.mean.forecast, file= "/Users/yannig/Documents/GAMwavelet/Hybrid_model_update2/Results/data.mean.forecast.RDS")

