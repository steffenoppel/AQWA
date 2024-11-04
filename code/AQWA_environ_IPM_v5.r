#########################################################################################################
#
# INTEGRATED POPULATION MODEL FOR AQUATIC WARBLER IN HUNGARY
#
#########################################################################################################
# written by steffen.oppel@rspb.org.uk in September 2016
# based on Kery and Schaub 2012
# triggered by request from Jochen Bellebaum and Zsolt Vegvari
# only count data available from Hungary
# demographic parameters cannot be informed by data (only priors)
# ultimate goal is to figure out whether winter (NAO, drought in Africa...) or summer (flooding, fire) explain extinction
# updated 2 Oct 2016 by Jochen Bellebaum
# updated 18 Oct 2016 by Steffen Oppel - introduced covariates and added graphical output
# updated 19 Oct 2016 by Steffen Oppel - introduced immigration and summer rain, and incorporated low count for July counts
# updated 19 Oct 2016 by Steffen Oppel - changed 'flood' from water level to binary variable (>130 cm water level)
# updated 20 Oct 2016 by Steffen Oppel - changed 'flood' to >140 cm water level, reduce prior for phij, and made immigration more flexible
# updated 20 Oct 2016 by Steffen Oppel - tried various approaches to constrain survival and fec - can all be achieved with immigration, but phij always > phia

# NO FURTHER IMPROVEMENT POSSIBLE UNLESS fec AND imm CAN BE DESCRIBED BY SOME ENVIRONMENTAL VARIABLES

# updated 10 Feb 2017 by Steffen Oppel - re-specified juvenile survival as proportion of adult survival to avoid juv surv being higher than ad surv
# updated 17 Feb vy Steffen ppel after discussion with Jochen Bellebaum - removed all co-variates to estimate correlation afterwards

library(jagsUI)
library(readxl)


############################################################################
#
# LOAD DATA
# 
##############################################################################
setwd("A:\\RSPB\\AquaticWarbler\\Analysis")
setwd("C:\\STEFFEN\\RSPB\\AquaticWarbler\\Analysis")


AW<-read_excel("IPM_inputdata.xlsx", sheet = "Tabelle1")
head(AW)
years<-AW$Year


# Population counts (from years 1976 to 2010)
y <- round(AW$AW.geomean,0)		# ENTER DATA OR READ IN AS VECTOR

# NAO data - Winter.NAO
NAO <- AW$Winter.NAO		# ENTER DATA OR READ IN AS VECTOR

# Precipitation data SAHEL - Winter.precip.anomaly
winrain <- AW$Winter.precip.anomaly		# ENTER DATA OR READ IN AS VECTOR

# Precipitation data in Hungary
rain <- AW$HU.precip.May		

# WATER LEVEL HUNGARY - HU.water.cm
flood <- AW$HU.water.cm

# Inundation data HUNGARY - IND.inundation
inund <- AW$IND.inundation

# Fire data HUNGARY - HU.burnt.area
fire <- AW$HU.burnt.area

# AREA AVAILABLE IN HUNGARY - HU.area
area <- log(AW$HU.area)				## use area on log scale for poisson draw of immigrants

# COUNT IN JUNE OR JULY
count <- ifelse(AW$AW.count=="June",1,2)


# STANDARDISE data(AVOID LARGE NUMERALS)

flood<-ifelse(flood>140,1,2)
flood[is.na(flood)]<-2			### set missing values to normal non-flood year

meanrain<-mean(rain, na.rm = TRUE)
sdrain<-sd(rain, na.rm = TRUE)
rain<-(rain-meanrain)/sdrain
rain[is.na(rain)]<-0			### set missing values to mean

meanarea<-mean(area, na.rm = TRUE)
sdarea<-sd(area, na.rm = TRUE)
area<-(area-meanarea)/sdarea
area[is.na(area)]<-0			### set missing values to mean


#jags.data <- list(nyears = length(years), y = y, NAO=NAO,flood=flood, fire=fire, inund=inund, winrain=winrain, rain=rain, count=count, area=area)
jags.data <- list(nyears = length(years), y = y, flood=flood, count=count, area=area)




##############################################################################
#
# IPM WITH FIXED ENVIRONMENTAL EFFECTS ON SURVIVAL AND FECUNDITY
# 
##############################################################################

sink("AQWA.IPM.null.jags")
cat("
model {


#----------------------------------------
# 1. Define the priors for the parameters
#----------------------------------------

# Initial population sizes
N1[1] ~ dnorm(20, 0.0001)T(0,)           # 1-year old individuals
NadSurv[1] ~ dnorm(60, 0.0001)T(0,)      # Adults >= 2 years
#imm[1]<-0					# set immigration to 0 in first year

# DEMOGRAPHIC PARAMETERS INFORMED BY DATA FROM OTHER STUDIES
mphij ~ dunif(0.3,0.5)       		# mean juvenile survival OFFSET - this is the proportion of adult survival that specifies how much lower juvenile survival is compared to adult survival
#mphij<- mphij.prop*mphia      		# mean juvenile survival as specified by multiplying offset by adult survival
mphia ~ dunif(0.45,0.55)		# mean adult survival; Dyrcz: apparent survival males 0.67 (95% CI: 0.537-0.785); females 0.42 (0.290-0.563). 
mfec ~ dunif(0,5)			# mean fecundity per breeding attempt, derived from Kubacka et al. 2013 by multiplying brood size and nest survival (3.8–4.1)*(0.36-0.87)

# CHANGE THE SCALE OF DEMOGRAPHIC PARAMETERS TO FACILITATE INCORPORATION OF COVARIATES
l.mphij<-log(mphij/(1-mphij))	# juvenile survival probability on logit scale;
l.mphia<-log(mphia/(1-mphia))			# adult survival probability on logit scale

# PRIORS FOR ENVIRONMENTAL EFFECTS
#beta.f.rain ~ dnorm(0, 0.0001)T(-10,10)		# prior for flood effect on fecundity
#beta.NAO ~ dnorm(0, 0.0001)T(-10,10)		# prior for NAO effect on survival
#beta.inund ~ dnorm(0, 0.0001)T(-10,10)		# prior for inundation effect on survival
#beta.winrain ~ dnorm(0, 0.0001)T(-10,10)		# prior for winter rain effect on survival
#beta.fire ~ dnorm(0, 0.0001)T(-10,10)		# prior for fire effect on fecundity
#beta.imm ~ dnorm(0, 0.0001)T(-10,10)		# prior for immigration that is governed by available area

# PRIORS FOR FLOOD EFFECT ON NUMBER OF BROODS
nbroods[1] ~ dunif(0,1)			# mean number of broods raised by the population in FLOOD YEARS
nbroods[2] ~ dunif(1.2,2)		# mean number of broods raised by the population in NORMAL YEARS

# PRIORS FOR COUNT EFFECT ON NUMBER OF OBSERVED BIRDS
count.p[2] ~ dunif(0,1)		# mean number of counted animals is lower when count is in July
count.p[1] <-1				# mean number of counted animals is normal

# RANDOM ANNUAL EFFECT ON PRODUCTIVITY
tau.fec<-1/pow(sigma.fec,2)
sigma.fec~dunif(0,5)

# RANDOM ANNUAL EFFECT ON SURVIVAL
tau.surv<-1/pow(sigma.surv,2)
sigma.surv~dunif(0,2)


#--------------------------------------------------
# 2. Relate parameters to environmental variables
#--------------------------------------------------

for (t in 1:(nyears-1)){
   eps.surv[t] ~ dnorm(0, tau.surv)												# random variation around fecundity
   logit(phij[t]) <- l.mphij+ eps.surv[t]			#beta.winrain*winrain[t]  # Juv. apparent survival ALSO TRY IND INUNDATION AND W.rain  + beta.NAO*NAO[t]  + beta.inund*inund[t]  
   logit(phia[t]) <- l.mphia+ eps.surv[t]			#+ beta.winrain*winrain[t]  # Adult apparent survival ALSO TRY IND INUNDATION AND W.rain  + beta.NAO*NAO[t]  + beta.inund*inund[t]  
   eps[t] ~ dnorm(0, tau.fec)												# random variation around fecundity
   fpre[t] <- mfec  + eps[t]     				# Productivity per breeding attempt beta.fire*fire[t] + beta.f.rain*rain[t] + 
   f[t] <- min(fpre[t],4.2)												# cap fledglings per brood to 4.2
   nbrood[t] <- nbroods[flood[t]]  											# FIXED TO 1 brood when flood and 2 broods when no flood
}


#-----------------------
# 3. Derived parameters
#-----------------------

# Population growth rate
for (t in 1:(nyears-1)){
   lambda[t] <- Ntot[t+1] / Ntot[t]
   }


#--------------------------------------------
# 4. Likelihood for population population count data
#--------------------------------------------
   # 4.1 System process
   for (t in 2:nyears){
      	chicks[t-1] <- 0.5 * f[t-1] * nbrood[t-1] * Ntot[t-1]		## nbroods is avg number of broods per year, f is fecundity per clutch, not per year
      	chicksrd[t-1] <- round(chicks[t-1])
		N1[t] ~ dbin(phij[t-1],chicksrd[t-1]) 
		NadSurv[t] ~ dbin(phia[t-1],round(Ntot[t-1]))
   		#log(imm.offset[t]) <- (beta.imm*area[t])*(flood[t]-1)			# IMMIGRATION SCALES TO AVAILABLE AREA and is only available in non-flood years
		#imm[t]~dpois(imm.offset[t]) 								# number of immigrants, will equal 0 in flood years
		#imm.offset[t]~dunif(0,200) 								# annually varying number of immigrant
		#imm[t]<-imm.offset[t]*(beta.imm*area[t])*(flood[t]-1)				# number of immigrants scaled to available area, will equal 0 in flood years
	}

   # 4.2 Observation process
   for (t in 1:nyears){
	capacity[t] ~ dnorm(700,0.01)     
	Ntot[t] <- min(max(1,(NadSurv[t] + N1[t])), capacity[t])			# +imm[t]
	Ntotobs[t] <- Ntot[t]*count.p[count[t]]								## includes factor <1 for counts in July		
	Ntotrd[t] <- round(Ntotobs[t])			
      y[t] ~ dpois(Ntotrd[t])
      }
}
",fill = TRUE)
sink()



# Initial values
inits <- function(){list(
nbroods= c(runif(1, 0,1),runif(1, 1.2,2)),
mphij = runif(1, 0.3,0.5),
mphia = runif(1, 0.45,0.55),
mfec = runif(1, 0,5))}

#mphia = runif(1, 0.29,0.785),
#mfec = runif(1, 1.1,3.7))}

# Parameters monitored
#parameters <- c("Ntot","phij","phia","f","nbrood","beta.f.rain","beta.fire","beta.winrain","lambda","count.p","mphij.prop")
parameters <- c("Ntot","phij","phia","f","nbrood","lambda","count.p","mphij.prop")


# MCMC settings
ni <- 150000
nt <- 2
nb <- 50000
nc <- 4

 
# Call JAGS from R
#ipm.model <- jags(jags.data, inits, parameters, "A:\\RSPB\\AquaticWarbler\\Analysis\\AQWA.IPM.imm5.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.cores=nc, parallel=T)
ipm.model <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\AquaticWarbler\\Analysis\\AQWA.IPM.null.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.cores=nc, parallel=T)




############################################################################
#
# WRITE ALL OUTPUT INTO A TEXT FILE
# 
##############################################################################
print(ipm.model, dig=3)
#write.table(ipm.model$summary, "AQWA_HU_Output_table_v5.csv", sep=",")




############################################################################
#
# MAKE A GRAPH OF THE POPULATION TRAJECTORY AND KEY DEMOGRAPHIC PARAMETERS
# 
##############################################################################
nyears = length(years)


#pdf("AQWA_Hungary_model_output_null.pdf", width=13, height=10)

par(mfrow=c(2,2), mar=c(4,5,0,0),oma=c(0,0,0,0))


### PLOT 1 - POPULATION TRAJECTORY

lower <- upper <- lowerhyp <- upperhyp <- numeric()
for (i in 1:nyears){
   lower[i] <- quantile(ipm.model $sims.list$Ntot[,i], 0.025)
   upper[i] <- quantile(ipm.model $sims.list$Ntot[,i], 0.975)
   }
plot(ipm.model $mean$Ntot, type = "b", ylim = c(0, 750), ylab = "Population size", xlab = "", las = 1, pch = 16, frame = F, cex = 1.5, axes=F, main="", cex.lab=1.5, col='blue')

axis(1, at=c(1:(nyears)), labels=c(1977:2011), cex.axis=1.3)
axis(2, at=seq(0,750,250), labels=T, las=1, cex.axis=1.3)
legend(x = 1, y = 750, legend = c("Counts", "Estimated trajectory"), pch = c(4, 16), col = c("red", "blue"), bty = "n")

segments(1:nyears, lower, 1:nyears, upper, col="black")
points(y, pch = 4, cex = 1.2, col="red")



### PLOT 2 - ADULT SURVIVAL

lower <- upper <- lowerhyp <- upperhyp <- numeric()
for (i in 1:(nyears-1)){
   lower[i] <- quantile(ipm.model $sims.list$phia[,i], 0.025)
   upper[i] <- quantile(ipm.model $sims.list$phia[,i], 0.975)
   }
plot(ipm.model $mean$phia, type = "b", ylim = c(0, 1), ylab = "Adult survival probability", xlab = "", las = 1, pch = 16, frame = F, cex = 1.5, axes=F, main="", cex.lab=1.5, col='blue')
axis(1, at=c(1:(nyears)), labels=c(1977:2011), cex.axis=1.3)
axis(2, at=seq(0,1,0.2), labels=T, las=1, cex.axis=1.3)
segments(1:(nyears-1), lower, 1:(nyears-1), upper, col="black")



### PLOT 3 - JUVENILE SURVIVAL

lower <- upper <- lowerhyp <- upperhyp <- numeric()
for (i in 1:(nyears-1)){
   lower[i] <- quantile(ipm.model $sims.list$phij[,i], 0.025)
   upper[i] <- quantile(ipm.model $sims.list$phij[,i], 0.975)
   }
plot(ipm.model $mean$phij, type = "b", ylim = c(0, 1), ylab = "Juvenile survival probability", xlab = "", las = 1, pch = 16, frame = F, cex = 1.5, axes=F, main="", cex.lab=1.5, col='blue')
axis(1, at=c(1:(nyears)), labels=c(1977:2011), cex.axis=1.3)
axis(2, at=seq(0,1,0.2), labels=T, las=1, cex.axis=1.3)
segments(1:(nyears-1), lower, 1:(nyears-1), upper, col="black")


### PLOT 4 - FECUNDITY

lower <- upper <- lowerhyp <- upperhyp <- numeric()
for (i in 1:(nyears-1)){
   lower[i] <- quantile(ipm.model $sims.list$f[,i], 0.025)*quantile(ipm.model $sims.list$nbrood[,i], 0.025)
   upper[i] <- quantile(ipm.model $sims.list$f[,i], 0.975)*quantile(ipm.model $sims.list$nbrood[,i], 0.975)
   }
plot(ipm.model $mean$f*ipm.model $mean$nbrood, type = "b", ylim = c(0, 10), ylab = "Fledglings / female", xlab = "", las = 1, pch = 16, frame = F, cex = 1.5, axes=F, main="", cex.lab=1.5, col='blue')
axis(1, at=c(1:(nyears)), labels=c(1977:2011), cex.axis=1.3)
axis(2, at=seq(0,10,2), labels=T, las=1, cex.axis=1.3)
segments(1:(nyears-1), lower, 1:(nyears-1), upper, col="black")


dev.off()




############################################################################
#
# CORRELATE OUTPUT WITH ENVIRONMENTAL EFFECTS
# 
##############################################################################


l.fitted<-l.lower<-l.upper<-ad.fitted<-ad.lower<-ad.upper<-ju.fitted<-ju.lower<-ju.upper<-pr.fitted<-pr.lower<-pr.upper<-numeric()

for (i in 1:(nyears-1)){
l.fitted[i]<-quantile(ipm.model $sims.list$lambda[,i], 0.5)
l.lower[i]<-quantile(ipm.model $sims.list$lambda[,i], 0.025)
l.upper[i]<-quantile(ipm.model $sims.list$lambda[,i], 0.975)

ad.fitted[i]<-quantile(ipm.model $sims.list$phia[,i], 0.5)
ad.lower[i]<-quantile(ipm.model $sims.list$phia[,i], 0.025)
ad.upper[i]<-quantile(ipm.model $sims.list$phia[,i], 0.975)

ju.fitted[i]<-quantile(ipm.model $sims.list$phij[,i], 0.5)
ju.lower[i]<-quantile(ipm.model $sims.list$phij[,i], 0.025)
ju.upper[i]<-quantile(ipm.model $sims.list$phij[,i], 0.975)

pr.fitted[i]<-quantile(ipm.model $sims.list$f[,i], 0.5)*quantile(ipm.model $sims.list$nbrood[,i], 0.5)
pr.lower[i]<-quantile(ipm.model $sims.list$f[,i], 0.025)*quantile(ipm.model $sims.list$nbrood[,i], 0.025)
pr.upper[i]<-quantile(ipm.model $sims.list$f[,i], 0.975)*quantile(ipm.model $sims.list$nbrood[,i], 0.975)}


############## NAO INDEX ###########

pdf("AQWA_Hungary_NAO_correlations.pdf", width=9, height=9)

par(mfrow=c(2,2), mar=c(4,5,0,0),oma=c(0,0,0,0))

plot(l.fitted~NAO[1:34], xlim=c(-2,2), ylim=c(0,2.0), xlab="NAO Index",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(NAO[1:34],l.lower,NAO[1:34],l.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(l.fitted,NAO[1:34],alternative = c("two.sided"),method = "spearman")
text(-2,2, sprintf("r = %f, p = %g",round(test$estimate,3), round(test$p.value,3)), adj=0)

plot(ad.fitted~NAO[1:34], xlim=c(-2,2), ylim=c(0,1.0), xlab="NAO Index",ylab="Adult survival probability",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(NAO[1:34],ad.lower,NAO[1:34],ad.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(ad.fitted,NAO[1:34],alternative = c("two.sided"),method = "spearman")
text(-2,1, sprintf("r = %f, p = %g",round(test$estimate,3), round(test$p.value,3)), adj=0)

plot(ju.fitted~NAO[1:34], xlim=c(-2,2), ylim=c(0,1.0), xlab="NAO Index",ylab="Juvenile survival probability",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(NAO[1:34],ju.lower,NAO[1:34],ju.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(ju.fitted,NAO[1:34],alternative = c("two.sided"),method = "spearman")
text(-2,1, sprintf("r = %f, p = %g",round(test$estimate,3), round(test$p.value,3)), adj=0)

plot(pr.fitted~NAO[1:34], xlim=c(-2,2), ylim=c(0,10), xlab="NAO Index",ylab="N fledglings per female",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(NAO[1:34],pr.lower,NAO[1:34],pr.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(pr.fitted,NAO[1:34],alternative = c("two.sided"),method = "spearman")
text(-2,10, sprintf("r = %f, p = %g",round(test$estimate,3), round(test$p.value,3)), adj=0)

dev.off()



############## WINRAIN ###########


pdf("AQWA_Hungary_winrain_correlations.pdf", width=9, height=9)

par(mfrow=c(2,2), mar=c(4,5,0,0),oma=c(0,0,0,0))

plot(l.fitted~winrain[1:34], xlim=c(-2500,500), ylim=c(0,2.0), xlab="Winter rainfall anomaly (mm)",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(winrain[1:34],l.lower,winrain[1:34],l.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(l.fitted,winrain[1:34],alternative = c("two.sided"),method = "spearman")
text(-2500,2, sprintf("r = %f, p = %g",round(test$estimate,3), round(test$p.value,3)), adj=0)

plot(ad.fitted~winrain[1:34], xlim=c(-2500,500), ylim=c(0,1.0), xlab="Winter rainfall anomaly (mm)",ylab="Adult survival probability",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(winrain[1:34],ad.lower,winrain[1:34],ad.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(ad.fitted,winrain[1:34],alternative = c("two.sided"),method = "spearman")
text(-2500,1, sprintf("r = %f, p = %g",round(test$estimate,3), round(test$p.value,3)), adj=0)

plot(ju.fitted~winrain[1:34], xlim=c(-2500,500), ylim=c(0,1.0), xlab="Winter rainfall anomaly (mm)",ylab="Juvenile survival probability",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(winrain[1:34],ju.lower,winrain[1:34],ju.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(ju.fitted,winrain[1:34],alternative = c("two.sided"),method = "spearman")
text(-2500,1, sprintf("r = %f, p = %g",round(test$estimate,3), round(test$p.value,3)), adj=0)

plot(pr.fitted~winrain[1:34], xlim=c(-2500,500), ylim=c(0,10), xlab="Winter rainfall anomaly (mm)",ylab="N fledglings per female",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(winrain[1:34],pr.lower,winrain[1:34],pr.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(pr.fitted,winrain[1:34],alternative = c("two.sided"),method = "spearman")
text(-2500,10, sprintf("r = %f, p = %g",round(test$estimate,3), round(test$p.value,3)), adj=0)

dev.off()

############## INUND ###########



pdf("AQWA_Hungary_inund_correlations.pdf", width=9, height=9)

par(mfrow=c(2,2), mar=c(4,5,0,0),oma=c(0,0,0,0))

plot(l.fitted~inund[1:34], xlim=c(7500,18000), ylim=c(0,2.0), xlab="Inner Niger Delta inundation",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(inund[1:34],l.lower,inund[1:34],l.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(l.fitted,inund[1:34],alternative = c("two.sided"),method = "spearman")
text(7500,2, sprintf("r = %f, p = %g",round(test$estimate,3), round(test$p.value,3)), adj=0)

plot(ad.fitted~inund[1:34], xlim=c(7500,18000), ylim=c(0,1.0), xlab="Inner Niger Delta inundation",ylab="Adult survival probability",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(inund[1:34],ad.lower,inund[1:34],ad.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(ad.fitted,inund[1:34],alternative = c("two.sided"),method = "spearman")
text(7500,1, sprintf("r = %f, p = %g",round(test$estimate,3), round(test$p.value,3)), adj=0)

plot(ju.fitted~inund[1:34], xlim=c(7500,18000), ylim=c(0,1.0), xlab="Inner Niger Delta inundation",ylab="Juvenile survival probability",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(inund[1:34],ju.lower,inund[1:34],ju.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(ju.fitted,inund[1:34],alternative = c("two.sided"),method = "spearman")
text(7500,1, sprintf("r = %f, p = %g",round(test$estimate,3), round(test$p.value,3)), adj=0)

plot(pr.fitted~inund[1:34], xlim=c(7500,18000), ylim=c(0,10), xlab="Inner Niger Delta inundation",ylab="N fledglings per female",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(inund[1:34],pr.lower,inund[1:34],pr.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(pr.fitted,inund[1:34],alternative = c("two.sided"),method = "spearman")
text(7500,10, sprintf("r = %f, p = %g",round(test$estimate,3), round(test$p.value,3)), adj=0)

dev.off()



############## RAIN ###########



pdf("AQWA_Hungary_rain_correlations.pdf", width=9, height=9)

par(mfrow=c(2,2), mar=c(4,5,0,0),oma=c(0,0,0,0))

plot(l.fitted~rain[1:34], xlim=c(0,200), ylim=c(0,2.0), xlab="Summer rainfall (mm)",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(rain[1:34],l.lower,rain[1:34],l.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(l.fitted,rain[1:34],alternative = c("two.sided"),method = "spearman")
text(0,2, sprintf("r = %f, p = %g",round(test$estimate,3), round(test$p.value,3)), adj=0)

plot(ad.fitted~rain[1:34], xlim=c(0,200), ylim=c(0,1.0), xlab="Summer rainfall (mm)",ylab="Adult survival probability",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(rain[1:34],ad.lower,rain[1:34],ad.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(ad.fitted,rain[1:34],alternative = c("two.sided"),method = "spearman")
text(0,1, sprintf("r = %f, p = %g",round(test$estimate,3), round(test$p.value,3)), adj=0)

plot(ju.fitted~rain[1:34], xlim=c(0,200), ylim=c(0,1.0), xlab="Summer rainfall (mm)",ylab="Juvenile survival probability",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(rain[1:34],ju.lower,rain[1:34],ju.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(ju.fitted,rain[1:34],alternative = c("two.sided"),method = "spearman")
text(0,1, sprintf("r = %f, p = %g",round(test$estimate,3), round(test$p.value,3)), adj=0)

plot(pr.fitted~rain[1:34], xlim=c(0,200), ylim=c(0,10), xlab="Summer rainfall (mm)",ylab="N fledglings per female",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(rain[1:34],pr.lower,rain[1:34],pr.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(pr.fitted,rain[1:34],alternative = c("two.sided"),method = "spearman")
text(0,10, sprintf("r = %f, p = %g",round(test$estimate,3), round(test$p.value,3)), adj=0)

dev.off()





