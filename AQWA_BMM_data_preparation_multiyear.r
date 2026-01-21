#####################################################################################################
#####     AQUATIC WARBLER TREND MONITORING FROM REPEATED TRANSECT SURVEYS           #################
##### 	DATA PREPARATION FOR MULTI-YEAR ANALYSIS IN JAGS  				#################
#####################################################################################################
## based on Kery et al 2005
## modified by steffen.oppel@rspb.org.uk on 24 Aug 2013
## modified by steffen.oppel@vogelwarte.ch in January 2026

YEAR<-2025					### LAST YEAR WITH DATA FOR MONITORING
nyears<-length(c(2011:YEAR))		### number of years for analysis
nsites<-50					### number of transects




rm(list=ls())
library(tidyverse)
library(data.table)
library(jagsUI)
library(readxl)
library(IPMbook)
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


hab<-fread("data/AQWA_habitat_POL.csv") %>%
  mutate(Date=dmy(Date)) %>%
  dplyr::filter(!is.na(Transect))
head(hab)





########## SPLIT SURVEYS DATA FRAME INTO SITE AND SURVEY SPECIFIC INFORMATION ###########
sites<-aggregate(tree_height~Transect_ID+Block+Year+Dead_vegetation+Water+Veg_height, data=surveys, FUN='unique')
surveys<-aggregate(Date~Transect_ID+Count+Year+MeanTemperatureC+MeanWindSpeedKmh+Precipitationmm+Observer, data=surveys, FUN='unique')
head(surveys)
dim(surveys)
head(sites)
dim(sites)


### CONVERT DATES INTO A CONTINUOUS NUMBER STARTING ON 10 MAY EACH YEAR

surveys$Date<-as.Date(surveys$Date, format="%d/%m/%Y")                 #formats the column into a Date object
surveys$year<-as.numeric(format(surveys$Date, format="%Y"))            # extracts the year from the Date column and creates a new year column
surveys$Year<-NULL							            # removes the column 3 (Year) so that year is the last column
for (y in 2011:2013){							# loops over years specified to calculate the season day in each year
startdate<-as.Date(sprintf("%d-05-10 00:00",y))             # set a startdate, here 10 May (must be smaller than first survey day)
surveys$Day[surveys$year==y]<-as.numeric(surveys$Date[surveys$year==y]-startdate)   #calculate the difference between the start date and the survey date
}
surveys$Date<-NULL
head(surveys)
head(birds)
#fix(surveys)
#fix(birds)
birds$Date<-NULL



### ADD THE NUMBER OF AQWA TO EACH SURVEY, including zero counts (which are not in the imported query)

surveys$AQWA <- 0  										 # creates a blank data frame that has zero for every count
surveys$Transect<-paste(surveys$Transect_ID,surveys$Count,surveys$year,sep = "_")		# creates a matching expression for each transect and count
birds$Transect<-paste(birds$Transect_ID,birds$Count, birds$Year ,sep = "_")			# creates a matching expression for each transect and count
surveys$AQWA[match(birds$Transect,surveys$Transect)]<-birds$SumOfnumber_males		# adds the data to the right transect and count based on the matching expression for each transect and count
surveys$AQWA[is.na(surveys$AQWA)]<-0								# all NA in the AQWA field are because of zero counts
head(surveys)


#### ADD MISSING SURVEYS ###########


if (dim(surveys)[1]==nyears*3*nsites){    					### this must be same dimensions as 50*3*2=300 
surveys<-surveys[order(surveys$Transect_ID,surveys$Count, decreasing=F),]
head(surveys)



######################################################################################################
# 
# 2. CREATE ARRAY WITH SURVEY DATA TO SAVE FOR JAGS -------------------------
# 
#######################################################################################################

### create array to be filled with data
AQWA.y<-array(NA, dim=c(nsites,3,nyears))

#### CAST THE MOLTEN DATA FRAME INTO MATRIX WITH 1 COLUMN PER COUNT and fill in array
for (y in 1:nyears){
b<-subset(surveys, year==c(2011:YEAR)[y])
b$Transect<-paste(b$Transect_ID,b$year,sep = "_")		# creates a matching expression for each transect and count
dis<-cast(b, Transect~Count, value="AQWA")							# pivot table to create data frame with one line per transect per year, and each column reflecting the observations per distance band on each count survey
dis<-dis[order(dis$Transect,decreasing=F),]
# fix(dis)
# head(dis)
AQWA.y[,,y]<-as.matrix(dis[,2:4], dimnames=NULL)
}

}else{print("something is wrong with the survey data - dimensions not correct")} ## close the IF loop that only runs when number of surveys is correct



######################################################################################################
# 
# 3. MANIPULATE SITE AND OBSERVATION COVARIATE DATA ---------------------------
# 
#######################################################################################################

head(surveys)
surveys$Transect<-paste(surveys$Transect_ID,surveys$year ,sep = "_")		# create new unique variable for year and transect
surveys<-surveys[order(surveys$Transect,surveys$Count,decreasing=F),]
surveys_red<-surveys[,c(10,2,3:8)]								# retained only columns with relevant obsCovariates
head(surveys_red)


###### TAKE MEAN VALUES OF HABITAT DATA ACROSS BLOCKS ####################################
head(sites)
str(sites)
surv<-aggregate(tree_height~Transect_ID+Year, data=sites, FUN='mean')
surv$Transect<-paste(surv$Transect_ID,surv$Year ,sep = "_")
sites$Transect<-paste(sites$Transect_ID,sites$Year ,sep = "_")

for (i in surv$Transect){
surv$veg1[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Veg_height==1))[1]/5
surv$veg2[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Veg_height==2))[1]/5
surv$veg3[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Veg_height==3))[1]/5
surv$veg4[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Veg_height==4))[1]/5

surv$wat1[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Water==1))[1]/5
surv$wat2[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Water==2))[1]/5
surv$wat3[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Water==3))[1]/5
surv$wat4[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Water==4))[1]/5

surv$lit1[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Dead_vegetation==1))[1]/5
surv$lit2[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Dead_vegetation==2))[1]/5
surv$lit3[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Dead_vegetation==3))[1]/5
surv$lit4[surv$Transect==i]<-dim(subset(sites[sites$Transect==i,], Dead_vegetation==4))[1]/5
}
head(surv)
dim(surv)


if ('B35' %in% surv$Transect_ID[surv$Year==2012]){
missing<-(surv[surv$Transect_ID=='B35' & surv$Year==2012,])
missing[,3:16]<-NA
missing$Year<-2011		### add year when the transect was skipped
surv<-rbind(surv,missing)}
surv$Transect<-paste(surv$Transect_ID,surv$Year,sep = "_")		# creates a matching expression for each transect and count

if ('B21' %in% surv$Transect_ID[surv$Year==2012]){
missing<-(surv[surv$Transect_ID=='B21' & surv$Year==2012,])
missing[,3:16]<-NA
missing$Year<-2013		### add year when the transect was skipped
surv<-rbind(surv,missing)}
surv$Transect<-paste(surv$Transect_ID,surv$Year,sep = "_")}		# creates a matching expression for each transect and count


sites<-surv[order(surv$Transect_ID, surv$Year,decreasing=F),]
for(i in 1:2){
sites[,i]<-as.factor(sites[,i])
}
sites$Transect<-NULL
dim(sites)      # this should be 150 transects
sites$yearnum<-ifelse(sites$Year=="2011",2011,ifelse(sites$Year=="2012",2012,2013)) 
head(sites)
#fix(sites)



######################################################################################################
# 
# 3. CREATE SITE LEVEL COVARIATES FOR ABUNDANCE AND STANDARDIZE THEM -----------------
# 
#######################################################################################################


## SORT THE TABLES SO THEY ALL HAVE THE SAME ORDER

siteCov<-sites[order(sites$Transect, decreasing=F),] 
head(siteCov)

### STANDARDISE SITE COVARIATES FOR USE IN WINBUGS

meantree<-mean(siteCov$tree_height, na.rm = TRUE)
sdtree<-sd(siteCov$tree_height, na.rm = TRUE)
siteCov$tree_height<-(siteCov$tree_height-meantree)/sdtree


#### CAST THE SITE COV DATA FRAME INTO MATRIX WITH 1 COLUMN PER YEAR
siteCovList<-array(NA, dim=c(nsites,nyears, 13))

for (col in 3:15){
b<-siteCov[,c(1,2,col)]
dis<-cast(b, Transect_ID~Year)							# pivot table to create data frame with one line per transect per year, and each column reflecting the observations per distance band on each count survey
dis<-dis[order(dis$Transect_ID,decreasing=F),]
# fix(dis)
# head(dis)
siteCovList[,,(col-2)]<-as.matrix(dis[,2:(nyears+1)], dimnames=NULL)
}
#treeheight<-siteCovList[,,1]
veg1<-siteCovList[,,2]
veg2<-siteCovList[,,3]
veg3<-siteCovList[,,4]
veg4<-siteCovList[,,5]

wat1<-siteCovList[,,6]
wat2<-siteCovList[,,7]
wat3<-siteCovList[,,8]
wat4<-siteCovList[,,9]

lit1<-siteCovList[,,10]
lit2<-siteCovList[,,11]
lit3<-siteCovList[,,12]
lit4<-siteCovList[,,13]


######################################################################################################
# 
# 4. CREATE OBSERVATION LEVEL COVARIATE ARRAYS FOR DETECTION AND STANDARDIZE THEM ------------------
# 
#######################################################################################################

## SORT THE TABLES SO THEY ALL HAVE THE SAME ORDER

obsCov<-surveys_red[order(surveys_red$Transect, decreasing=F),] 
head(obsCov)


### STANDARDIZE COVARIATES FOR WINBUGS

meant<-mean(obsCov$MeanTemperatureC, na.rm = TRUE)
sdt<-sd(obsCov$MeanTemperatureC, na.rm = TRUE)
obsCov$temp<-(obsCov$MeanTemperatureC-meant)/sdt

meant<-mean(obsCov$ MeanWindSpeedKmh, na.rm = TRUE)
sdt<-sd(obsCov$MeanWindSpeedKmh, na.rm = TRUE)
obsCov$wind<-(obsCov$ MeanWindSpeedKmh-meant)/sdt

meant<-mean(obsCov$Precipitationmm, na.rm = TRUE)
sdt<-sd(obsCov$Precipitationmm, na.rm = TRUE)
obsCov$rain<-(obsCov$Precipitationmm-meant)/sdt

meant<-mean(obsCov$Day, na.rm = TRUE)
sdt<-sd(obsCov$Day, na.rm = TRUE)
obsCov$Day<-(obsCov$Day-meant)/sdt


### create array for each covariate

rain<-array(NA, dim=c(nsites,3,nyears))
wind<-array(NA, dim=c(nsites,3,nyears))
temp<-array(NA, dim=c(nsites,3,nyears))
day<-array(NA, dim=c(nsites,3,nyears))


### fill in array for each covariate
for (y in 2011:YEAR){
obsC<-subset(obsCov, year==y)
#obsC<-merge(siteCov,obsC, by="Transect",all.x=T)		## abandoned. ensured NA for counts that did not take place and ensured proper size of array

y<-match(y,c(2011:YEAR))						## translates the year (2011, 2012, etc.) into consecutive number (1,2,...) for array dimensions
x<-cast(obsC, Transect ~ Count, value='temp')
x2<-as.matrix(x[,2:4])
temp[,,y]<-x2

x<-cast(obsC, Transect ~ Count, value='rain')
rain[,,y]<-as.matrix(x[,2:4])

x<-cast(obsC, Transect ~ Count, value='wind')
wind[,,y]<-as.matrix(x[,2:4])

x<-cast(obsC, Transect ~ Count, value='Day')
day[,,y]<-as.matrix(x[,2:4])
}



######################################################################################################
# 
# 5. REPLACE ALL NA IN COVARIATES otherwise "undefined node" error ------------------
# 
#######################################################################################################


veg1[is.na(veg1)]<-0
veg4[is.na(veg4)]<-0
veg2[is.na(veg2)]<-0
veg3[is.na(veg3)]<-0
veg23<-veg2+veg3

lit1[is.na(lit1)]<-0
lit4[is.na(lit4)]<-0
lit2[is.na(lit2)]<-0
lit3[is.na(lit3)]<-0

wat1[is.na(wat1)]<-0
wat4[is.na(wat4)]<-0
wat2[is.na(wat2)]<-0
wat3[is.na(wat3)]<-0

for (d in 1:nyears){							### replace missing dates with mean for each survey round in each year
day[is.na(day[,1,d]),1,d]<-mean(day[,1,d], na.rm=T)
day[is.na(day[,2,d]),2,d]<-mean(day[,2,d], na.rm=T)
day[is.na(day[,3,d]),3,d]<-mean(day[,3,d], na.rm=T)
temp[is.na(temp[,1,d]),1,d]<-mean(temp[,1,d], na.rm=T)
temp[is.na(temp[,2,d]),2,d]<-mean(temp[,2,d], na.rm=T)
temp[is.na(temp[,3,d]),3,d]<-mean(temp[,3,d], na.rm=T)
wind[is.na(wind[,1,d]),1,d]<-mean(wind[,1,d], na.rm=T)
wind[is.na(wind[,2,d]),2,d]<-mean(wind[,2,d], na.rm=T)
wind[is.na(wind[,3,d]),3,d]<-mean(wind[,3,d], na.rm=T)
rain[is.na(rain[,1,d]),1,d]<-mean(rain[,1,d], na.rm=T)
rain[is.na(rain[,2,d]),2,d]<-mean(rain[,2,d], na.rm=T)
rain[is.na(rain[,3,d]),3,d]<-mean(rain[,3,d], na.rm=T)
}


######################################################################################################
# 
# 6. CREATE INPUT DATA FOR JAGS ------------------------
# 
#######################################################################################################

# check data dimensions
dim(AQWA.y)
dim(veg2)
dim(lit2)
dim(veg4)
dim(lit4)
dim(wat4)


dim(day)
dim(wind)
dim(rain)
dim(temp)



#### SAVE WORKSPACE
rm(AQWA_count,b,birds,bugs.dir,col,dis,i,meant,meantree,missing,obsC,obsCov,sdt,sdtree, siteCov, siteCovList,sites,startdate,surv,surveys,surveys_red,x,x2,y,YEAR)
ls()
save.image("C:\\STEFFEN\\RSPB\\AquaticWarbler\\Analysis\\Trend_monitoring\\AQWA_BMM_input.RData")
save.image("A:\\RSPB\\AquaticWarbler\\Analysis\\Trend_monitoring\\AQWA_BMM_input.RData")

bugs.data<-list(M = AQWA.y, nsite=nsites, nrep=3, primocc=seq(1:nyears), nyear=nyears,veg2=veg2,veg3=veg3,veg4=veg4, wat2=wat2,wat3=wat3,wat4=wat4)
rm(AQWA.y,d,day,lit1,lit2,lit3,lit4,nsites,nyears,rain,temp,veg1,veg2,veg23,veg3,veg4,wat1,wat2,wat3,wat4,wind)
save.image("C:\\STEFFEN\\MANUSCRIPTS\\Submitted\\AQWA_census\\PeerageofScience\\AquaticWarbler_monitoring_data.RData")


