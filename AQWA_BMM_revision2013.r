#######################################################################################################################################
#################################### AQUATIC WARBLER MONITORING DATA ANALYSIS 2012 ####################################################
#######################################################################################################################################
# written by steffen.oppel@rspb.org.uk
# included 2012 data on 30 Oct 2012
# TRANSECT B35 was only surveyed in 2012!
# changed veg_height to proportions on 16 Nov 2012
# included additional models to overcome transect heterogeneity
# 17 May 2012: updated script to remove veg1 (redundant in combination with veg2+veg3+veg4)!
# 7 August 2013: reverted to BMM (rather than gdistsamp) due to reviewer criticism
# 8 October 2013 - removed Observer from all models as confounded with abundance
# 13 December 2013 - included data from 2013, CHANGED B44 Count 1 in 2013 from 13 to 3 birds.



############### LOAD THE REQUIRED PACKAGES #########################################################

library(RODBC)
library(unmarked)
library(Hmisc)
library(reshape)
library(slam)

############### SET THE WORKING DIRECTORY #########################################################

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\AquaticWarbler\\Raw_Data")
setwd("C:\\STEFFEN\\RSPB\\AquaticWarbler\\Raw_Data")


#####################################################################################################################################################
#############    LOAD RAW DATA AND MANIPULATE FOR UNMARKED FORMAT      ##############################################################################
#####################################################################################################################################################

AQWA_count <- odbcConnectAccess2007 ('AQWA_transect_monitoring.accdb')
surveys <- sqlQuery(AQWA_count, "SELECT * FROM block_surveys")
birds <- sqlQuery(AQWA_count, "SELECT * FROM transect_totals")
odbcClose(AQWA_count)
head(surveys)

setwd("S:\\ConSci\\DptShare\\SteffenOppel\\RSPB\\AquaticWarbler\\Analysis\\Trend_monitoring")
setwd("C:\\STEFFEN\\RSPB\\AquaticWarbler\\Analysis\\Trend_monitoring")

########## SPLIT SURVEYS DATA FRAME INTO SITE AND SURVEY SPECIFIC INFORMATION ###########
sites<-aggregate(tree_height~Transect_ID+Block+Year+Dead_vegetation+Water+Veg_height, data=surveys, FUN='unique')
surveys<-aggregate(Date~Transect_ID+Count+Year+MeanTemperatureC+MeanWindSpeedKmh+Precipitationmm+Observer, data=surveys, FUN='unique')
head(surveys)
dim(surveys)
head(sites)
dim(sites)


########## FIGURE OUT DISTRIBUTION OF VEGETATION HEIGHT ACROSS BIEBRZA ##################
#### changed on 6 Nov 2012 - the previous version $density somehow underestimated the lowest height category
#### changed on 12 Nov 2012 - Lawki Marsh is only 4690 ha
#### changed back on 6 Dec because not logical to use only part of Biebrza

BIEBRZA<-8280
OTOP<-4690

#VEG_HEIGHT_AREA2012<-((hist(sites$Veg_height[sites$Year==2012], plot=F, breaks=c(0,1.5,2.5,3.5,4.5))$counts)/sum(hist(sites$Veg_height[sites$Year==2012], plot=F, breaks=c(0,1.5,2.5,3.5,4.5))$counts))*OTOP
#VEG_HEIGHT_AREA2011<-((hist(sites$Veg_height[sites$Year==2011], plot=F, breaks=c(0,1.5,2.5,3.5,4.5))$counts)/sum(hist(sites$Veg_height[sites$Year==2011], plot=F, breaks=c(0,1.5,2.5,3.5,4.5))$counts))*OTOP

VEG_HEIGHT_prop2013<-((hist(sites$Veg_height[sites$Year==2013], plot=F, breaks=c(0,1.5,2.5,3.5,4.5))$counts)/sum(hist(sites$Veg_height[sites$Year==2013], plot=F, breaks=c(0,1.5,2.5,3.5,4.5))$counts))
VEG_HEIGHT_prop2012<-((hist(sites$Veg_height[sites$Year==2012], plot=F, breaks=c(0,1.5,2.5,3.5,4.5))$counts)/sum(hist(sites$Veg_height[sites$Year==2012], plot=F, breaks=c(0,1.5,2.5,3.5,4.5))$counts))
VEG_HEIGHT_prop2011<-((hist(sites$Veg_height[sites$Year==2011], plot=F, breaks=c(0,1.5,2.5,3.5,4.5))$counts)/sum(hist(sites$Veg_height[sites$Year==2011], plot=F, breaks=c(0,1.5,2.5,3.5,4.5))$counts))




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
### this is critical to have correct dimensions of the input matrix for both the bird numbers and the observation covariates
if ('B27' %in% surveys$Transect_ID){
missing<-(surveys[surveys$Transect_ID=='B27',])[1,]
missing$Count<-3
#missing[,4:10]<-NA
surveys<-rbind(surveys,missing)}
dim(surveys)   ### this must be same dimensions as 49*3*2=294+3 (B35 only surveyed in 2012) = 297
surveys<-surveys[order(surveys$Transect_ID,surveys$Count, decreasing=F),]
head(surveys)


#####################################################################################################################################################
#############    CREATE UNMARKED DATA FRAME       ###################################################################################################
#####################################################################################################################################################


#### CAST THE MOLTEN DATA FRAME INTO MATRIX WITH 1 COLUMN PER COUNT AND DISTANCE CATEGORY (9 COLUMNS)
surveys$Transect<-paste(surveys$Transect_ID,surveys$year,sep = "_")		# creates a matching expression for each transect and count
dis<-cast(surveys, Transect~Count, value="AQWA")							# pivot table to create data frame with one line per transect per year, and each column reflecting the observations per distance band on each count survey
# fix(dis)
head(dis)
AQWA.d<-as.matrix(dis[,2:4], dimnames=NULL)							      # convert to matrix which is required for 'unmarked'
head(AQWA.d)
dim(AQWA.d)													# must be 148 - 49 transects in 2011, and 50 in 2012, 49 in 2013


### OBSERVATION COVARIATE DATA FRAME for unmarked Frame class

head(surveys)
surveys$Transect<-paste(surveys$Transect_ID,surveys$year ,sep = "_")		# create new unique variable for year and transect
surveys<-surveys[order(surveys$Transect,surveys$Count,decreasing=F),]
surveys_red<-surveys[,c(3:8,10)]								# retained only columns with relevant obsCovariates
head(surveys_red)
dim(surveys_red)   ### this must be same dimensions as 3*dim(AQWA.d)



### SITE COVARIATE DATA FRAME for unmarked Frame class
### needs one value per transect for vegetation etc. variables that are measured at block level
### changed below on 31 Oct 2012 to use metric covariates for water, litter, and veg height
### changed on 17 Nov to use prop of each veg height and water




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

sites<-surv[order(surv$Transect_ID, surv$Year,decreasing=F),]
for(i in 1:2){
sites[,i]<-as.factor(sites[,i])
}
sites$Transect<-NULL
dim(sites)      # this should be 61 transects
sites$yearnum<-ifelse(sites$Year=="2011",2011,2012) 
head(sites)
#fix(sites)




##### COMBINE RESPONSE AND OBSERVATION COVARIATES TO UNMARKED FRAME AND STANDARDIZE NUMERIC COVARIATES #######

AQWAdis<-unmarkedFramePCount(y=AQWA.d, siteCovs=sites, obsCovs=surveys_red)
siteCovs(AQWAdis)[,c(3,16)] <- scale(siteCovs(AQWAdis)[,c(3,16)])
obsCovs(AQWAdis)[,c(1:3,6)] <- scale(obsCovs(AQWAdis)[,c(1:3,6)])
summary(AQWAdis)


##################################################################################################
######### CHECK FOR OVERDISPERSION ###############################################################
##################################################################################################

surveys$veg1<-rep(sites$veg1, each=3)
surveys$veg2<-rep(sites$veg2, each=3)
surveys$veg3<-rep(sites$veg3, each=3)
surveys$veg4<-rep(sites$veg4, each=3)

test<-glm(AQWA~Transect_ID+year+veg1+veg2+veg3+veg4+Observer+MeanTemperatureC, data=surveys, family=quasipoisson)
summary(test)

### DISPERSION PARAMETER = 0.9938887 (0.9619426 when all Biebrza is considered)



###################################################################################################
######## ANALYSIS OF DATA #########################################################################
###################################################################################################

### first formula is for detection, second formula is for abundance  ####

#####  BIOLOGICALLY PLAUSIBLE MODELS EVALUATED  #####

veg_temp_obs <- pcount(~ Year + MeanTemperatureC ~ 1+veg2+veg3+veg4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
veg_rain_obs <- pcount(~ Year + Precipitationmm ~ 1+veg2+veg3+veg4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
veg_day_obs <- pcount(~ Year + Day ~ 1+veg2+veg3+veg4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
veg_wind_obs <- pcount(~ Year + MeanWindSpeedKmh ~ 1+veg2+veg3+veg4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
veg_observer <- pcount(~ Year + 1 ~ 1+veg2+veg3+veg4+ yearnum, data= AQWAdis,  mixture = "P", K=70)

litter_temp_obs <- pcount(~ Year + MeanTemperatureC ~ 1+lit2+lit3+lit4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
litter_rain_obs <- pcount(~ Year + Precipitationmm ~ 1+lit2+lit3+lit4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
litter_day_obs <- pcount(~ Year + Day ~ 1+lit2+lit3+lit4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
litter_wind_obs <- pcount(~ Year + MeanWindSpeedKmh ~ 1+lit2+lit3+lit4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
litter_observer <- pcount(~ Year + 1 ~ 1+lit2+lit3+lit4+ yearnum, data= AQWAdis,  mixture = "P", K=70)

water_temp_obs <- pcount(~ Year + MeanTemperatureC ~ 1+ wat2+wat3+wat4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
water_rain_obs <- pcount(~ Year + Precipitationmm ~ 1+ wat2+wat3+wat4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
water_day_obs <- pcount(~ Year + Day ~ 1+ wat2+wat3+wat4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
water_wind_obs <- pcount(~ Year + MeanWindSpeedKmh ~ 1+ wat2+wat3+wat4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
water_observer <- pcount(~ Year + 1 ~ 1+ wat2+wat3+wat4+ yearnum, data= AQWAdis,  mixture = "P", K=70)


#####  SAME SUITE OF MODELS ACCOUNTING FOR VEGETATION HEIGHT DETECTION DIFFERENCES  #####

veg_temp_VegYear <- pcount(~ Year + MeanTemperatureC+veg2+veg3+veg4 ~ 1+veg2+veg3+veg4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
veg_rain_VegYear <- pcount(~ Year + Precipitationmm+veg2+veg3+veg4 ~ 1+veg2+veg3+veg4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
veg_day_VegYear <- pcount(~ Year + Day+veg2+veg3+veg4 ~ 1+veg2+veg3+veg4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
veg_wind_VegYear <- pcount(~ Year + MeanWindSpeedKmh+veg2+veg3+veg4 ~ 1+veg2+veg3+veg4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
veg_VegYear <- pcount(~ Year + 1+veg2+veg3+veg4 ~ 1+veg2+veg3+veg4+ yearnum, data= AQWAdis,  mixture = "P", K=70)

litter_temp_VegYear <- pcount(~ Year + MeanTemperatureC+veg2+veg3+veg4 ~ 1+lit2+lit3+lit4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
litter_rain_VegYear <- pcount(~ Year + Precipitationmm+veg2+veg3+veg4 ~ 1+lit2+lit3+lit4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
litter_day_VegYear <- pcount(~ Year + Day+veg2+veg3+veg4 ~ 1+lit2+lit3+lit4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
litter_wind_VegYear <- pcount(~ Year + MeanWindSpeedKmh+veg2+veg3+veg4 ~ 1+lit2+lit3+lit4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
litter_VegYear <- pcount(~ Year + 1+veg2+veg3+veg4 ~ 1+lit2+lit3+lit4+ yearnum, data= AQWAdis,  mixture = "P", K=70)

water_temp_VegYear <- pcount(~ Year + MeanTemperatureC+veg2+veg3+veg4 ~ 1+ wat2+wat3+wat4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
water_rain_VegYear <- pcount(~ Year + Precipitationmm+veg2+veg3+veg4 ~ 1+ wat2+wat3+wat4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
water_day_VegYear <- pcount(~ Year + Day+veg2+veg3+veg4 ~ 1+ wat2+wat3+wat4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
water_wind_VegYear <- pcount(~ Year + MeanWindSpeedKmh+veg2+veg3+veg4 ~ 1+ wat2+wat3+wat4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
water_VegYear <- pcount(~ Year + 1+veg2+veg3+veg4 ~ 1+ wat2+wat3+wat4+ yearnum, data= AQWAdis,  mixture = "P", K=70)

#####  SAME SUITE OF MODELS WITHOUT ANNUAL DETECTABILITY DIFFERENCES  #####

veg_temp_obsVeg <- pcount(~ MeanTemperatureC+veg2+veg3+veg4 ~ 1+veg2+veg3+veg4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
veg_rain_obsVeg <- pcount(~ Precipitationmm+veg2+veg3+veg4 ~ 1+veg2+veg3+veg4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
veg_day_obsVeg <- pcount(~ Day+veg2+veg3+veg4 ~ 1+veg2+veg3+veg4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
veg_wind_obsVeg <- pcount(~ MeanWindSpeedKmh+veg2+veg3+veg4 ~ 1+veg2+veg3+veg4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
veg_observerVeg <- pcount(~ 1+veg2+veg3+veg4 ~ 1+veg2+veg3+veg4+ yearnum, data= AQWAdis,  mixture = "P", K=70)

litter_temp_obsVeg <- pcount(~ MeanTemperatureC+veg2+veg3+veg4 ~ 1+lit2+lit3+lit4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
litter_rain_obsVeg <- pcount(~ Precipitationmm+veg2+veg3+veg4 ~ 1+lit2+lit3+lit4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
litter_day_obsVeg <- pcount(~ Day+veg2+veg3+veg4 ~ 1+lit2+lit3+lit4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
litter_wind_obsVeg <- pcount(~ MeanWindSpeedKmh+veg2+veg3+veg4 ~ 1+lit2+lit3+lit4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
litter_observerVeg <- pcount(~ 1+veg2+veg3+veg4 ~ 1+lit2+lit3+lit4+ yearnum, data= AQWAdis,  mixture = "P", K=70)

water_temp_obsVeg <- pcount(~ MeanTemperatureC+veg2+veg3+veg4 ~ 1+ wat2+wat3+wat4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
water_rain_obsVeg <- pcount(~ Precipitationmm+veg2+veg3+veg4 ~ 1+ wat2+wat3+wat4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
water_day_obsVeg <- pcount(~ Day+veg2+veg3+veg4 ~ 1+ wat2+wat3+wat4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
water_wind_obsVeg <- pcount(~ MeanWindSpeedKmh+veg2+veg3+veg4 ~ 1+ wat2+wat3+wat4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
water_observerVeg <- pcount(~ 1+veg2+veg3+veg4 ~ 1+ wat2+wat3+wat4+ yearnum, data= AQWAdis,  mixture = "P", K=70)

null <- pcount(~ 1 ~ 1+ yearnum, data= AQWAdis,  mixture = "P", K=70)


fl <- fitList(null=null, water_observer=water_observer,water_wind_obs=water_wind_obs,water_day_obs=water_day_obs,water_rain_obs=water_rain_obs,water_temp_obs=water_temp_obs,litter_observer=litter_observer,litter_wind_obs=litter_wind_obs,litter_day_obs=litter_day_obs, veg_temp_obs=veg_temp_obs, veg_rain_obs=veg_rain_obs, veg_day_obs=veg_day_obs, veg_wind_obs=veg_wind_obs, veg_observer=veg_observer, litter_temp_obs=litter_temp_obs, litter_rain_obs=litter_rain_obs,
water_observerVeg=water_observerVeg,water_wind_obsVeg=water_wind_obsVeg,water_day_obsVeg=water_day_obsVeg,water_rain_obsVeg=water_rain_obsVeg,water_temp_obsVeg=water_temp_obsVeg,litter_observerVeg=litter_observerVeg,litter_wind_obsVeg=litter_wind_obsVeg,litter_day_obsVeg=litter_day_obsVeg, veg_temp_obsVeg=veg_temp_obsVeg, veg_rain_obsVeg=veg_rain_obsVeg, veg_day_obsVeg=veg_day_obsVeg, veg_wind_obsVeg=veg_wind_obsVeg, veg_observerVeg=veg_observerVeg, litter_temp_obsVeg=litter_temp_obsVeg, litter_rain_obsVeg=litter_rain_obsVeg,
water_VegYear=water_VegYear,water_wind_VegYear=water_wind_VegYear,water_day_VegYear=water_day_VegYear,water_rain_VegYear=water_rain_VegYear,water_temp_VegYear=water_temp_VegYear, veg_temp_VegYear=veg_temp_VegYear, veg_rain_VegYear=veg_rain_VegYear, veg_day_VegYear=veg_day_VegYear, veg_wind_VegYear=veg_wind_VegYear, veg_VegYear=veg_VegYear,litter_VegYear=litter_VegYear,litter_wind_VegYear=litter_wind_VegYear,litter_day_VegYear=litter_day_VegYear, litter_temp_VegYear=litter_temp_VegYear, litter_rain_VegYear=litter_rain_VegYear)

ms <- modSel(fl, nullmod="null")
ms
sink("AQWA_BMM_AIC_table2013_no_obs.txt")
ms
sink()
 


###################################################################################################
######## SUMMARY OF TOP MODEL #########################################################################
###################################################################################################

### run model again without intercept so extraction of real parameter estimates is easier
top <- pcount(~ -1+Year + Precipitationmm+veg2+veg3+veg4 ~ 1+ wat2+wat3+wat4+ yearnum, data= AQWAdis,  mixture = "P", K=70)
out<-summary(top)



setwd("S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\AQWA_census\\BCI")
setwd("C:\\STEFFEN\\MANUSCRIPTS\\in_press\\AQWA_census")


out_tab<-rbind(out$state, out$det)
write.table(out_tab, "clipboard", sep="\t")


###################################################################################################
######## GOF TEST OF TOP MODEL #########################################################################
###################################################################################################

# Function adopted from Sillet et al. 2012
freeTuke <- function(fm) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    sum((sqrt(observed) - sqrt(expected))^2)
    }

# should set nsim=100 or more, but nsim=200 takes a few hours to run
pb <- parboot(top, freeTuke, nsim=20, report=1)





########################################################
#### DENSITY ESTIMATES AND PROJECTED POPULATON SIZE ####
########################################################

### after removing 'Year' as predictor in model:
new<-data.frame(veg1=rep(0,4), veg2=rep(1,4), veg3=rep(0,4), veg4=rep(0,4), yearnum=0, Precipitationmm=0,wat1=c(1,0,0,0), wat2=c(0,1,0,0), wat3=c(0,0,1,0), wat4=c(0,0,0,1))

RESULT<-predict(top,type='state',newdata=new,appendData=T)
names(RESULT)[1:4]<-c('density','se','lcl','ucl')
#RESULT$plot<-c(1,2,3,4,1.2,2.2,3.2,4.2)
RESULT$plot<-c(1,2,3,4)

RESULT


########################################################################################
### PLOTTING A GRAPH showing mean abundance across water heights   ##################
########################################################################################
pdf("Fig2.pdf", width=9,height=7)
jpeg("Fig2.jpg", width=800,height=600, units = "px", quality = 100)
tiff(filename = "Fig2.tiff",width=1600,height=1200, units = "px", compression = "none", pointsize = 28)
par(mar=c(5,6,0,0))
errbar(RESULT$plot,RESULT$density,RESULT$lcl, RESULT$ucl, type='p', cex=1.6, xlab="Water level", ylab="Aquatic Warbler abundance (males/transect)", axes=F,xlim=c(0.5,4.5), ylim=c(0,35), frame=F, cex.lab=1.6, mgp=c(3.5,1,0))
axis(1, at=c(0.5,1,2,3,4,4.5), labels=c("","dry","wet ground","< 15 cm","> 15 cm",""), cex.lab=1.6, cex.axis=1.6)
axis(2, at=seq(0,35,5), labels=T, cex.lab=1.6, cex.axis=1.6, las=1)
dev.off()

### FOR MULTIPLE YEARS:
pdf("AQWA_abundance_water_heightLawki.pdf", width=8,height=6)
par(mar=c(5,4.5,0,0))
errbar(RESULT$plot[RESULT$Year==2011],RESULT$density[RESULT$Year==2011],RESULT$lcl[RESULT$Year==2011], RESULT$ucl[RESULT$Year==2011], type='p', cex=1.6, xlab="Water depth", ylab="", axes=F,xlim=c(0.5,4.5), ylim=c(0,30), frame=F, cex.lab=1.6, mgp=c(3.5,1,0))
par(new=T)
errbar(RESULT$plot[RESULT$Year==2012],RESULT$density[RESULT$Year==2012],RESULT$lcl[RESULT$Year==2012], RESULT$ucl[RESULT$Year==2012], type='p', pch=1, cex=1.6, xlab="", ylab="Aquatic Warbler abundance (males/transect)", axes=F,xlim=c(0.5,4.5), ylim=c(0,30), frame=F, cex.lab=1.6, mgp=c(4,1,0))
axis(1, at=c(0.5,1.1,2.1,3.1,4.1,4.5), labels=c("","dry","wet ground","< 15 cm","> 15 cm",""), cex.lab=1.6, cex.axis=1.6)
axis(2, at=seq(0,30,5), labels=T, cex.lab=1.6, cex.axis=1.6, las=1)
legend(0.5,30, pch=c(16,1), cex=1.6, legend=c("2011","2012"), bty='n')
dev.off()



########################################################################################
### PLOTTING A GRAPH showing mean detection probability across vegetation heights   #####
########################################################################################
### after removing 'Year' as predictor in model:
pdf("Fig3.pdf", width=8,height=6)
jpeg("Fig3.jpg", width=800,height=600, units = "px", quality = 100)
tiff(filename = "Fig3.tiff",width=1600,height=1200, units = "px", compression = "none", pointsize = 28)

par(mar=c(5,6,0,0))

new<-data.frame(wat1=rep(0,12), wat2=rep(1,12), wat3=rep(0,12), wat4=rep(0,12), Precipitationmm=1.25,veg1=rep(c(1,0,0,0),3), veg2=rep(c(0,1,0,0),3), veg3=rep(c(0,0,1,0),3), veg4=rep(c(0,0,0,1),3), Year=rep(as.factor(c('2011','2012','2013')), each=4))
RESULT<-predict(top,type='det',newdata=new,appendData=T)
names(RESULT)[1:4]<-c('density','se','lcl','ucl')
RESULT$plot<-c(1,2,3,4,1.1,2.1,3.1,4.1,1.2,2.2,3.2,4.2)

#pdf("AQWA_detection_probability_VegHeight.pdf", width=8,height=6)
#par(mar=c(5,4.5,0,0))
errbar(RESULT$plot,RESULT$density,RESULT$lcl, RESULT$ucl, type='p', pch=rep(c(5,16,4),each=4), cex=1.6, xlab="Vegetation height", ylab="Detection probability", axes=F,xlim=c(0.5,4.5), ylim=c(0,1), frame=F, cex.lab=1.6, mgp=c(3.5,1,0))
axis(1, at=c(0.5,1.1,2.1,3.1,4.1,4.5), labels=c("","< 40 cm","41-80 cm","81-120 cm","> 120 cm",""), cex.lab=1.6, cex.axis=1.6)
axis(2, at=seq(0,1,0.2), labels=T, cex.lab=1.6, cex.axis=1.6, las=1, mgp=c(5.5,0.7,0))
#text(0.5,1,"A", cex= 2.3)
legend(3.7,1,pch=c(5,16,4), legend=c('2011','2012','2013'), bty='n',cex=1.4)
dev.off()






#####################################################################################################################################
### DISCARDED MODELS INCLUDING OBSERVER AND YEAR EFFECT    ######################
#####################################################################################################################################


#####  BIOLOGICALLY PLAUSIBLE MODELS EVALUATED  #####
temp_obs <- pcount(~ MeanTemperatureC+Observer ~ Year, data= AQWAdis,  mixture = "P", K=70)
rain_obs <- pcount(~ Precipitationmm+Observer ~ Year, data= AQWAdis,  mixture = "P", K=70)
day_obs <- pcount(~ Day+Observer ~ Year, data= AQWAdis,  mixture = "P", K=70)
wind_obs <- pcount(~ MeanWindSpeedKmh+Observer ~ Year, data= AQWAdis,  mixture = "P", K=70)
observer <- pcount(~Observer ~ Year, data= AQWAdis,  mixture = "P", K=70)
null <- pcount(~1 ~ Year, data= AQWAdis,  mixture = "P", K=70)

veg_temp_obs <- pcount(~ MeanTemperatureC+Observer ~ Year+veg2+veg3+veg4, data= AQWAdis,  mixture = "P", K=70)
veg_rain_obs <- pcount(~ Precipitationmm+Observer ~ Year+veg2+veg3+veg4, data= AQWAdis,  mixture = "P", K=70)
veg_day_obs <- pcount(~ Day+Observer ~ Year+veg2+veg3+veg4, data= AQWAdis,  mixture = "P", K=70)
veg_wind_obs <- pcount(~ MeanWindSpeedKmh+Observer ~ Year+veg2+veg3+veg4, data= AQWAdis,  mixture = "P", K=70)
veg_observer <- pcount(~ 1+Observer ~ Year+veg2+veg3+veg4, data= AQWAdis,  mixture = "P", K=70)

litter_temp_obs <- pcount(~ MeanTemperatureC+Observer ~ Year+lit2+lit3+lit4, data= AQWAdis,  mixture = "P", K=70)
litter_rain_obs <- pcount(~ Precipitationmm+Observer ~ Year+lit2+lit3+lit4, data= AQWAdis,  mixture = "P", K=70)
litter_day_obs <- pcount(~ Day+Observer ~ Year+lit2+lit3+lit4, data= AQWAdis,  mixture = "P", K=70)
litter_wind_obs <- pcount(~ MeanWindSpeedKmh+Observer ~ Year+lit2+lit3+lit4, data= AQWAdis,  mixture = "P", K=70)
litter_observer <- pcount(~ 1+Observer ~ Year+lit2+lit3+lit4, data= AQWAdis,  mixture = "P", K=70)

water_temp_obs <- pcount(~ MeanTemperatureC+Observer ~ Year+ wat2+wat3+wat4, data= AQWAdis,  mixture = "P", K=70)
water_rain_obs <- pcount(~ Precipitationmm+Observer ~ Year+ wat2+wat3+wat4, data= AQWAdis,  mixture = "P", K=70)
water_day_obs <- pcount(~ Day+Observer ~ Year+ wat2+wat3+wat4, data= AQWAdis,  mixture = "P", K=70)
water_wind_obs <- pcount(~ MeanWindSpeedKmh+Observer ~ Year+ wat2+wat3+wat4, data= AQWAdis,  mixture = "P", K=70)
water_observer <- pcount(~ 1+Observer ~ Year+ wat2+wat3+wat4, data= AQWAdis,  mixture = "P", K=70)


#####  SAME SUITE OF MODELS ACCOUNTING FOR VEGETATION HEIGHT DETECTION DIFFERENCES  #####
temp_obsVeg <- pcount(~ MeanTemperatureC+Observer+veg2+veg3+veg4 ~ Year, data= AQWAdis,  mixture = "P", K=70)
rain_obsVeg <- pcount(~ Precipitationmm+Observer+veg2+veg3+veg4 ~ Year, data= AQWAdis,  mixture = "P", K=70)
day_obsVeg <- pcount(~ Day+Observer+veg2+veg3+veg4 ~ Year, data= AQWAdis,  mixture = "P", K=70)
wind_obsVeg <- pcount(~ MeanWindSpeedKmh+Observer+veg2+veg3+veg4 ~ Year, data= AQWAdis,  mixture = "P", K=70)
observerVeg <- pcount(~ 1+Observer+veg2+veg3+veg4 ~ Year, data= AQWAdis,  mixture = "P", K=70)

veg_temp_obsVeg <- pcount(~ MeanTemperatureC+Observer+veg2+veg3+veg4 ~ Year+veg2+veg3+veg4, data= AQWAdis,  mixture = "P", K=70)
veg_rain_obsVeg <- pcount(~ Precipitationmm+Observer+veg2+veg3+veg4 ~ Year+veg2+veg3+veg4, data= AQWAdis,  mixture = "P", K=70)
veg_day_obsVeg <- pcount(~ Day+Observer+veg2+veg3+veg4 ~ Year+veg2+veg3+veg4, data= AQWAdis,  mixture = "P", K=70)
veg_wind_obsVeg <- pcount(~ MeanWindSpeedKmh+Observer+veg2+veg3+veg4 ~ Year+veg2+veg3+veg4, data= AQWAdis,  mixture = "P", K=70)
veg_observerVeg <- pcount(~ 1+Observer+veg2+veg3+veg4 ~ Year+veg2+veg3+veg4, data= AQWAdis,  mixture = "P", K=70)

litter_temp_obsVeg <- pcount(~ MeanTemperatureC+Observer+veg2+veg3+veg4 ~ Year+lit2+lit3+lit4, data= AQWAdis,  mixture = "P", K=70)
litter_rain_obsVeg <- pcount(~ Precipitationmm+Observer+veg2+veg3+veg4 ~ Year+lit2+lit3+lit4, data= AQWAdis,  mixture = "P", K=70)
litter_day_obsVeg <- pcount(~ Day+Observer+veg2+veg3+veg4 ~ Year+lit2+lit3+lit4, data= AQWAdis,  mixture = "P", K=70)
litter_wind_obsVeg <- pcount(~ MeanWindSpeedKmh+Observer+veg2+veg3+veg4 ~ Year+lit2+lit3+lit4, data= AQWAdis,  mixture = "P", K=70)
litter_observerVeg <- pcount(~ 1+Observer+veg2+veg3+veg4 ~ Year+lit2+lit3+lit4, data= AQWAdis,  mixture = "P", K=70)

water_temp_obsVeg <- pcount(~ MeanTemperatureC+Observer+veg2+veg3+veg4 ~ Year+ wat2+wat3+wat4, data= AQWAdis,  mixture = "P", K=70)
water_rain_obsVeg <- pcount(~ Precipitationmm+Observer+veg2+veg3+veg4 ~ Year+ wat2+wat3+wat4, data= AQWAdis,  mixture = "P", K=70)
water_day_obsVeg <- pcount(~ Day+Observer+veg2+veg3+veg4 ~ Year+ wat2+wat3+wat4, data= AQWAdis,  mixture = "P", K=70)
water_wind_obsVeg <- pcount(~ MeanWindSpeedKmh+Observer+veg2+veg3+veg4 ~ Year+ wat2+wat3+wat4, data= AQWAdis,  mixture = "P", K=70)
water_observerVeg <- pcount(~ 1+Observer+veg2+veg3+veg4 ~ Year+ wat2+wat3+wat4, data= AQWAdis,  mixture = "P", K=70)



########################################################################################
### PLOTTING A GRAPH showing mean density across vegetation heights   ##################
########################################################################################


#pdf("AQWA_density_veg_heightLawki.pdf", width=8,height=6)
par(mar=c(5,6,0,0))
errbar(RESULT$plot[RESULT$Year==2011],RESULT$density[RESULT$Year==2011],RESULT$lcl[RESULT$Year==2011], RESULT$ucl[RESULT$Year==2011], type='p', cex=1.6, xlab="Vegetation height", ylab="", axes=F,xlim=c(0.5,4.5), ylim=c(0,1), frame=F, cex.lab=1.6, mgp=c(3.5,1,0))
par(new=T)
errbar(RESULT$plot[RESULT$Year==2012],RESULT$density[RESULT$Year==2012],RESULT$lcl[RESULT$Year==2012], RESULT$ucl[RESULT$Year==2012], type='p', pch=1, cex=1.6, xlab="", ylab="Aquatic Warbler density (males/ha)", axes=F,xlim=c(0.5,4.5), ylim=c(0,1), frame=F, cex.lab=1.6, mgp=c(4,1,0))
axis(1, at=c(0.5,1.1,2.1,3.1,4.1,4.5), labels=c("","< 40 cm","40-80 cm","80-120 cm","> 120 cm",""), cex.lab=1.6, cex.axis=1.6)
axis(2, at=seq(0,1,0.2), labels=T, cex.lab=1.6, cex.axis=1.6, las=1)
legend(0.5,1.0, pch=c(16,1), cex=1.6, legend=c("2011","2012"), bty='n')
#dev.off()




########################################################################################
### PLOTTING A GRAPH showing mean detection probability over time   #####
########################################################################################
#new<-data.frame(wat2=1, wat3=0, wat4=0, Day=seq(9,53,1), veg2=1, veg3=0, veg4=0)
new<-data.frame(wat2=1, wat3=0, wat4=0, Precipitationmm=seq(0,99.7,1/3), veg2=1, veg3=0, veg4=0, Year=as.factor(c('2011','2012','2013')))

RESULT<-predict(top,type='det',newdata=new,appendData=T)
names(RESULT)[1:4]<-c('density','se','lcl','ucl')
plot(RESULT$density[RESULT$Year=='2013']~RESULT$Precipitationmm[RESULT$Year=='2013'], type='l', frame=F, xlim=c(0,100),ylim=c(0,1),xlab="", ylab="Detection probability", cex.lab=1.6, cex.axis=1.6, axes=F)
par (new=T)
plot(RESULT$lcl[RESULT$Year=='2013']~RESULT$Precipitationmm[RESULT$Year=='2013'], type='l', lty=2, frame=F, xlim=c(0,100),ylim=c(0,1),xlab="Daily precipitation (mm)", ylab="", cex.lab=1.6, cex.axis=1.6, axes=F)
par (new=T)
plot(RESULT$ucl[RESULT$Year=='2013']~RESULT$Precipitationmm[RESULT$Year=='2013'], type='l', lty=2, frame=F, xlim=c(0,100),ylim=c(0,1),xlab="", ylab="", cex.lab=1.6, cex.axis=1.6, axes=F)



rx <- range(surveys$Day)
startdate<-as.Date("2013-05-10") 
range<-c(startdate+rx)
daterange<-seq(range[1], range[2],by='5 days')

axis(2, at=seq(0,1,0.2), labels=T, cex.lab=1.6, cex.axis=1.6, las=1, mgp=c(5.5,0.7,0))
axis(1, at=seq(0,100,10), labels=T, cex.lab=1.6, cex.axis=1.6)		#format(daterange, "%d %b")
text(0,1,"B", cex= 2.3)

dev.off()








