#####################################################################################################
#####     AQUATIC WARBLER TREND MONITORING FROM REPEATED TRANSECT SURVEYS           #################
##### 	DATA PREPARATION FOR MULTI-YEAR ANALYSIS IN JAGS  				#################
#####################################################################################################
## based on Kery et al 2005
## modified by steffen.oppel@rspb.org.uk on 24 Aug 2013
## modified by steffen.oppel@vogelwarte.ch in January 2026


rm(list=ls())
YEAR<-2025					### LAST YEAR WITH DATA FOR MONITORING
nyears<-length(c(2011:YEAR))		### number of years for analysis
nsites<-50					### number of transects




library(tidyverse)
library(data.table)
library(jagsUI)
library(readxl)
library(IPMbook)
library(reshape2)
filter<-dplyr::filter
select<-dplyr::select




############################################################################
#
# 1. LOAD DATA-------------
# 
##############################################################################
try(setwd("C:/Users/jba/OneDrive - Vogelwarte/Projects/Aquatic Warbler Steffen/Github repo/AQWA"),silent=T)
try(setwd("C:/STEFFEN/OneDrive - Vogelwarte/ExternalCollaborations/AQWA"),silent=T)
try(setwd("C:/Users/sop/OneDrive - Vogelwarte/ExternalCollaborations/AQWA"),silent=T)

AW<-fread("data/AQWA_counts_POL.csv") %>%
  mutate(Date=dmy(Date))
head(AW)
years<-unique(AW$Year)
table(AW$No_singing_males)

hab<-fread("data/AQWA_habitat_POL.csv") %>%
  mutate(Date=dmy(Date)) %>%
  dplyr::filter(!is.na(Transect))
head(hab)
dim(hab)
dim(AW)


### CONVERT DATES INTO A CONTINUOUS NUMBER STARTING ON the earliest day ever surveyed
AW$yday<-yday(AW$Date)
min(AW$yday)
AW$yday<-AW$yday-(min(AW$yday)-1)



########## SPLIT SURVEYS DATA FRAME INTO SITE AND SURVEY SPECIFIC INFORMATION ###########
# sites<-aggregate(tree_height~Transect_ID+Block+Year+Dead_vegetation+Water+Veg_height, data=surveys, FUN='unique')
# surveys<-aggregate(Date~Transect_ID+Count+Year+MeanTemperatureC+MeanWindSpeedKmh+Precipitationmm+Observer, data=surveys, FUN='unique')
# head(surveys)
# dim(surveys)
# head(sites)
# dim(sites)

sites<-hab %>%
  dplyr::filter(!is.na(Transect)) %>%
  dplyr::filter(Transect!="") %>%
  group_by(Year, Transect) %>%
  summarise(veg_height=mean(Av_veget_height),tree=mean(Trees_shrubs),tree_height=mean(Tree_shrub_av_height), reed=mean(Reed), water=mean(Water_level), litter=mean(Litter)) %>%
  ungroup()
  #mutate(site=paste(Transect,"_",Part_of_transect))
dim(sites)


habsurveys<-hab %>%
  group_by(Year, Transect) %>%
  summarise(day=first(Date)) %>%
  ungroup()
dim(habsurveys)

surveys<-AW %>%
  group_by(Year,Transect_no, Control) %>%
  summarise(N=sum(No_singing_males),day=first(yday))
dim(surveys)


### TRYING TO UNDERSTAND HOW OFTEN EACH TRANSECT WAS COUNTED

surveys %>%
  group_by(Year,Transect_no) %>%
  summarise(counts=max(Control)) %>%
  pivot_wider(names_from='Year', values_from='counts')


nsites<-length(unique(sites$Transect))
nyears<-length(unique(surveys$Year))


# ### ADD THE NUMBER OF AQWA TO EACH SURVEY, including zero counts (which are not in the imported query)
# 
# surveys$AQWA <- 0  										 # creates a blank data frame that has zero for every count
# surveys$Transect<-paste(surveys$Transect_ID,surveys$Count,surveys$year,sep = "_")		# creates a matching expression for each transect and count
# birds$Transect<-paste(birds$Transect_ID,birds$Count, birds$Year ,sep = "_")			# creates a matching expression for each transect and count
# surveys$AQWA[match(birds$Transect,surveys$Transect)]<-birds$SumOfnumber_males		# adds the data to the right transect and count based on the matching expression for each transect and count
# surveys$AQWA[is.na(surveys$AQWA)]<-0								# all NA in the AQWA field are because of zero counts
# head(surveys)
# 
# 
# #### ADD MISSING SURVEYS ###########
# 
# 
# if (dim(surveys)[1]==nyears*3*nsites){    					### this must be same dimensions as 50*3*2=300 
# surveys<-surveys[order(surveys$Transect_ID,surveys$Count, decreasing=F),]
# head(surveys)



######################################################################################################
# 
# 2. CREATE ARRAY WITH SURVEY DATA TO SAVE FOR JAGS -------------------------
# 
#######################################################################################################

# number of transects changes over time
# number of controls varies by year and transect
# needs robust loop to generate a full array


# If duplicates per cell **do not** exist:
AQWA.y <- acast(surveys, Transect_no ~ Control ~ Year, value.var = "N")
obs.cov.day <- acast(surveys, Transect_no ~ Control ~ Year, value.var = "day")
obs.cov.eff <- acast(surveys, Transect_no ~ Control ~ Year, value.var = "N")
obs.cov.eff <- ifelse(is.na(obs.cov.eff),0,1)
             




# ######################################################################################################
# # 
# # 3. MANIPULATE SITE AND OBSERVATION COVARIATE DATA ---------------------------
# # 
# #######################################################################################################
# 
# head(surveys)
# surveys$Transect<-paste(surveys$Transect_ID,surveys$year ,sep = "_")		# create new unique variable for year and transect
# surveys<-surveys[order(surveys$Transect,surveys$Count,decreasing=F),]
# surveys_red<-surveys[,c(10,2,3:8)]								# retained only columns with relevant obsCovariates
# head(surveys_red)
# 
# 
# ###### TAKE MEAN VALUES OF HABITAT DATA ACROSS BLOCKS ####################################
# head(sites)
# str(sites)
# surv<-aggregate(tree_height~Transect_ID+Year, data=sites, FUN='mean')
# surv$Transect<-paste(surv$Transect_ID,surv$Year ,sep = "_")
# sites$Transect<-paste(sites$Transect_ID,sites$Year ,sep = "_")
# 
# for (i in surv$Transect){
# surv$veg1[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Veg_height==1))[1]/5
# surv$veg2[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Veg_height==2))[1]/5
# surv$veg3[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Veg_height==3))[1]/5
# surv$veg4[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Veg_height==4))[1]/5
# 
# surv$wat1[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Water==1))[1]/5
# surv$wat2[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Water==2))[1]/5
# surv$wat3[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Water==3))[1]/5
# surv$wat4[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Water==4))[1]/5
# 
# surv$lit1[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Dead_vegetation==1))[1]/5
# surv$lit2[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Dead_vegetation==2))[1]/5
# surv$lit3[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Dead_vegetation==3))[1]/5
# surv$lit4[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Dead_vegetation==4))[1]/5
# }
# head(surv)
# dim(surv)
# 
# 
# 
# 
# ######################################################################################################
# # 
# # 3. CREATE SITE LEVEL COVARIATES FOR ABUNDANCE AND STANDARDIZE THEM -----------------
# # 
# #######################################################################################################
# 
# 
# ## SORT THE TABLES SO THEY ALL HAVE THE SAME ORDER
# 
# siteCov<-sites[order(sites$Transect, decreasing=F),] 
# head(siteCov)
# 
# ### STANDARDISE SITE COVARIATES FOR USE IN WINBUGS
# 
# meantree<-mean(siteCov$tree_height, na.rm = TRUE)
# sdtree<-sd(siteCov$tree_height, na.rm = TRUE)
# siteCov$tree_height<-(siteCov$tree_height-meantree)/sdtree
# 
# 
# #### CAST THE SITE COV DATA FRAME INTO MATRIX WITH 1 COLUMN PER YEAR
# siteCovList<-array(NA, dim=c(nsites,nyears, 13))
# 
# for (col in 3:15){
# b<-siteCov[,c(1,2,col)]
# dis<-cast(b, Transect_ID~Year)							# pivot table to create data frame with one line per transect per year, and each column reflecting the observations per distance band on each count survey
# dis<-dis[order(dis$Transect_ID,decreasing=F),]
# # fix(dis)
# # head(dis)
# siteCovList[,,(col-2)]<-as.matrix(dis[,2:(nyears+1)], dimnames=NULL)
# }
# #treeheight<-siteCovList[,,1]
# veg1<-siteCovList[,,2]
# veg2<-siteCovList[,,3]
# veg3<-siteCovList[,,4]
# veg4<-siteCovList[,,5]
# 
# wat1<-siteCovList[,,6]
# wat2<-siteCovList[,,7]
# wat3<-siteCovList[,,8]
# wat4<-siteCovList[,,9]
# 
# lit1<-siteCovList[,,10]
# lit2<-siteCovList[,,11]
# lit3<-siteCovList[,,12]
# lit4<-siteCovList[,,13]
# 
# 
# ######################################################################################################
# # 
# # 4. CREATE OBSERVATION LEVEL COVARIATE ARRAYS FOR DETECTION AND STANDARDIZE THEM ------------------
# # 
# #######################################################################################################
# 
# ## SORT THE TABLES SO THEY ALL HAVE THE SAME ORDER
# 
# obsCov<-surveys_red[order(surveys_red$Transect, decreasing=F),] 
# head(obsCov)
# 
# 
# ### STANDARDIZE COVARIATES FOR WINBUGS
# 
# meant<-mean(obsCov$MeanTemperatureC, na.rm = TRUE)
# sdt<-sd(obsCov$MeanTemperatureC, na.rm = TRUE)
# obsCov$temp<-(obsCov$MeanTemperatureC-meant)/sdt
# 
# meant<-mean(obsCov$ MeanWindSpeedKmh, na.rm = TRUE)
# sdt<-sd(obsCov$MeanWindSpeedKmh, na.rm = TRUE)
# obsCov$wind<-(obsCov$ MeanWindSpeedKmh-meant)/sdt
# 
# meant<-mean(obsCov$Precipitationmm, na.rm = TRUE)
# sdt<-sd(obsCov$Precipitationmm, na.rm = TRUE)
# obsCov$rain<-(obsCov$Precipitationmm-meant)/sdt
# 
# meant<-mean(obsCov$Day, na.rm = TRUE)
# sdt<-sd(obsCov$Day, na.rm = TRUE)
# obsCov$Day<-(obsCov$Day-meant)/sdt
# 
# 
# ### create array for each covariate
# 
# rain<-array(NA, dim=c(nsites,3,nyears))
# wind<-array(NA, dim=c(nsites,3,nyears))
# temp<-array(NA, dim=c(nsites,3,nyears))
# day<-array(NA, dim=c(nsites,3,nyears))
# 
# 
# ### fill in array for each covariate
# for (y in 2011:YEAR){
# obsC<-subset(obsCov, year==y)
# #obsC<-merge(siteCov,obsC, by="Transect",all.x=T)		## abandoned. ensured NA for counts that did not take place and ensured proper size of array
# 
# y<-match(y,c(2011:YEAR))						## translates the year (2011, 2012, etc.) into consecutive number (1,2,...) for array dimensions
# x<-cast(obsC, Transect ~ Count, value='temp')
# x2<-as.matrix(x[,2:4])
# temp[,,y]<-x2
# 
# x<-cast(obsC, Transect ~ Count, value='rain')
# rain[,,y]<-as.matrix(x[,2:4])
# 
# x<-cast(obsC, Transect ~ Count, value='wind')
# wind[,,y]<-as.matrix(x[,2:4])
# 
# x<-cast(obsC, Transect ~ Count, value='Day')
# day[,,y]<-as.matrix(x[,2:4])
# }
# 
# 
# 
# ######################################################################################################
# # 
# # 5. REPLACE ALL NA IN COVARIATES otherwise "undefined node" error ------------------
# # 
# #######################################################################################################
# 
# 
# veg1[is.na(veg1)]<-0
# veg4[is.na(veg4)]<-0
# veg2[is.na(veg2)]<-0
# veg3[is.na(veg3)]<-0
# veg23<-veg2+veg3
# 
# lit1[is.na(lit1)]<-0
# lit4[is.na(lit4)]<-0
# lit2[is.na(lit2)]<-0
# lit3[is.na(lit3)]<-0
# 
# wat1[is.na(wat1)]<-0
# wat4[is.na(wat4)]<-0
# wat2[is.na(wat2)]<-0
# wat3[is.na(wat3)]<-0
# 
# for (d in 1:nyears){							### replace missing dates with mean for each survey round in each year
# day[is.na(day[,1,d]),1,d]<-mean(day[,1,d], na.rm=T)
# day[is.na(day[,2,d]),2,d]<-mean(day[,2,d], na.rm=T)
# day[is.na(day[,3,d]),3,d]<-mean(day[,3,d], na.rm=T)
# temp[is.na(temp[,1,d]),1,d]<-mean(temp[,1,d], na.rm=T)
# temp[is.na(temp[,2,d]),2,d]<-mean(temp[,2,d], na.rm=T)
# temp[is.na(temp[,3,d]),3,d]<-mean(temp[,3,d], na.rm=T)
# wind[is.na(wind[,1,d]),1,d]<-mean(wind[,1,d], na.rm=T)
# wind[is.na(wind[,2,d]),2,d]<-mean(wind[,2,d], na.rm=T)
# wind[is.na(wind[,3,d]),3,d]<-mean(wind[,3,d], na.rm=T)
# rain[is.na(rain[,1,d]),1,d]<-mean(rain[,1,d], na.rm=T)
# rain[is.na(rain[,2,d]),2,d]<-mean(rain[,2,d], na.rm=T)
# rain[is.na(rain[,3,d]),3,d]<-mean(rain[,3,d], na.rm=T)
# }
# 
# 
# ######################################################################################################
# # 
# # 6. CREATE INPUT DATA FOR JAGS ------------------------
# # 
# #######################################################################################################
# 
# # check data dimensions
# dim(AQWA.y)
# dim(veg2)
# dim(lit2)
# dim(veg4)
# dim(lit4)
# dim(wat4)
# 
# 
# dim(day)
# dim(wind)
# dim(rain)
# dim(temp)
# 


#### SAVE WORKSPACE

try(setwd("C:/Users/jba/OneDrive - Vogelwarte/Projects/Aquatic Warbler Steffen/Github repo/AQWA"),silent=T)
try(setwd("C:/STEFFEN/OneDrive - Vogelwarte/ExternalCollaborations/AQWA"),silent=T)
try(setwd("C:/Users/sop/OneDrive - Vogelwarte/ExternalCollaborations/AQWA"),silent=T)
rm(AQWA_count,b,birds,bugs.dir,col,dis,i,meant,meantree,missing,obsC,obsCov,sdt,sdtree, siteCov, siteCovList,sites,startdate,surv,x,x2,y,YEAR)
ls()
save.image("data\\AQWA_trend_input.RData")

bugs.data<-list(M = AQWA.y, nsite=nsites, nrep=3, primocc=seq(1:nyears), nyear=nyears,veg2=veg2,veg3=veg3,veg4=veg4, wat2=wat2,wat3=wat3,wat4=wat4)
rm(AQWA.y,d,day,lit1,lit2,lit3,lit4,nsites,nyears,rain,temp,veg1,veg2,veg23,veg3,veg4,wat1,wat2,wat3,wat4,wind)
save.image("C:\\STEFFEN\\MANUSCRIPTS\\Submitted\\AQWA_census\\PeerageofScience\\AquaticWarbler_monitoring_data.RData")




######################################################################################################
# 
# 7. SPECIFY TREND MODEL IN JAGS ------------------------
# 
#######################################################################################################


model{
  
  # Priors
  loglam~dunif(-5,5)##mean abundance prior
  trend~dunif(-10,10)##    trend prior
  
  for(i in 1:nsite){
    lam.site[i]~dnorm(0,tau.site)T(-20, 20)## needs to be T(-20, 20) in JAGS
  }
  tau.site<-1/(sigma.site*sigma.site)
  sigma.site~dunif(0,10)
  for(year in 1:nyear){
    p0[year]~dunif(0,1)
    logitp0[year]<-log(p0[year]/(1-p0[year]))
  }
  tau.lp<-1/(sigma.p*sigma.p)
  sigma.p~dunif(0,10)
  
  
  # State and observation models
  for(year in 1:nyear){
    for(i in 1:nsite){
      log(lambda[i,year])<- loglam + trend*primocc[year]+lam.site[i]## 'primocc' is the primary occasion - here: years
      #log(lambda[i,year])<- loglam + lam.site[i]
      N[i,year]~dpois(lambda[i,year])
      
      for(t in 1:nrep){
        M[i,t,year]~dbin(p[i,t,year],N[i,year])
        p[i,t,year] <- exp(lp[i,t,year])/(1+exp(lp[i,t,year]))
        lp[i,t,year] ~ dnorm(mu.lp[i,t,year], tau.lp)T(-20, 20) # H.c.
        mu.lp[i,t,year]<-logitp0[year]
      }
    }
  }
  
  # Computation of fit statistic (Bayesian p-value)
  # Fit statistic for observed data
  # Also, generate replicate data and compute fit stats for them
  for(year in 1:nyear){
    totalN[year]<-sum(N[,year])
    
    for(i in 1:nsite){
      for(t in 1:nrep){
        
        # Actual data
        eval[i,t,year] <-N[i,year]*p[i,t,year] # Expected value
        sd.resi[i,t,year]<-sqrt(eval[i,t,year]*(1-p[i,t,year])) +0.5
        E[i,t,year]<-(M[i,t,year]-eval[i,t,year])/ sd.resi[i,t,year]
        E2[i,t,year] <- pow(E[i,t,year],2)
        
        # Replicate data sets
        M.new[i,t,year]~dbin(p[i,t,year],N[i,year])
        E.new[i,t,year]<-(M.new[i,t,year]-eval[i,t,year])/sd.resi[i,t,year]
        E2.new[i,t,year] <- pow(E.new[i,t,year], 2)
      }
    }
  }
  fit <- sum(E2[,,])# Sum up squared residuals for actual data set
  fit.new <- sum(E2.new[,,]) # Sum up for replicate data sets
}

