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





# STANDARDISE AND data(AVOID LARGE NUMERALS)
meanNAO<-mean(NAO, na.rm = TRUE)
sdNAO<-sd(NAO, na.rm = TRUE)
NAO<-(NAO-meanNAO)/sdNAO

#meanflood<-mean(flood, na.rm = TRUE)
#sdflood<-sd(flood, na.rm = TRUE)
#flood<-(flood-meanflood)/sdflood
flood<-ifelse(flood>140,1,2)
flood[is.na(flood)]<-2			### set missing values to normal non-flood year

meanwinrain<-mean(winrain, na.rm = TRUE)
sdwinrain<-sd(winrain, na.rm = TRUE)
winrain<-(winrain-meanwinrain)/sdwinrain
winrain[is.na(winrain)]<-0			### set missing values to mean

meanrain<-mean(rain, na.rm = TRUE)
sdrain<-sd(rain, na.rm = TRUE)
rain<-(rain-meanrain)/sdrain
rain[is.na(rain)]<-0			### set missing values to mean


meaninund<-mean(inund, na.rm = TRUE)
sdinund<-sd(inund, na.rm = TRUE)
inund<-(inund-meaninund)/sdinund
inund[is.na(inund)]<-0			### set missing values to mean


meanfire<-mean(fire, na.rm = TRUE)
sdfire<-sd(fire, na.rm = TRUE)
fire<-(fire-meanfire)/sdfire
fire[is.na(fire)]<-0			### set missing values to mean

meanarea<-mean(area, na.rm = TRUE)
sdarea<-sd(area, na.rm = TRUE)
area<-(area-meanarea)/sdarea
area[is.na(area)]<-0			### set missing values to mean


jags.data <- list(nyears = length(years), y = y, NAO=NAO,flood=flood, fire=fire, inund=inund, winrain=winrain, rain=rain, count=count, area=area)
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
mphij.prop ~ dunif(0,1)       		# mean juvenile survival OFFSET - this is the proportion of adult survival that specifies how much lower juvenile survival is compared to adult survival
mphij<- mphij.prop*mphia      		# mean juvenile survival as specified by multiplying offset by adult survival
mphia ~ dunif(0.2,0.9)		# mean adult survival; Dyrcz: apparent survival males 0.67 (95% CI: 0.537-0.785); females 0.42 (0.290-0.563). 
mfec ~ dunif(0.2,4)			# mean fecundity per breeding attempt, derived from Kubacka et al. 2013 by multiplying brood size and nest survival (3.8–4.1)*(0.36-0.87)

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
sigma.fec~dunif(0,2)

# RANDOM ANNUAL EFFECT ON SURVIVAL
tau.surv<-1/pow(sigma.surv,2)
sigma.surv~dunif(0,5)


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
mphij.prop = runif(1, 0,1),
mphia = runif(1, 0.45,0.75),
mfec = runif(1, 0.2,3.0))}

#mphia = runif(1, 0.29,0.785),
#mfec = runif(1, 1.1,3.7))}

# Parameters monitored
parameters <- c("Ntot","phij","phia","f","nbrood","beta.f.rain","beta.fire","beta.winrain","lambda","count.p","mphij.prop")


# MCMC settings
ni <- 15000
nt <- 1
nb <- 10000
nc <- 4

 
# Call JAGS from R
#ipm.model <- jags(jags.data, inits, parameters, "A:\\RSPB\\AquaticWarbler\\Analysis\\AQWA.IPM.imm5.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.cores=nc, parallel=T)
ipm.model <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\AquaticWarbler\\Analysis\\AQWA.IPM.imm5.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.cores=nc, parallel=T)




############################################################################
#
# WRITE ALL OUTPUT INTO A TEXT FILE
# 
##############################################################################
print(ipm.model, dig=3)
write.table(ipm.model$summary, "AQWA_HU_Output_table_v5.csv", sep=",")




############################################################################
#
# MAKE A GRAPH OF THE POPULATION TRAJECTORY AND KEY DEMOGRAPHIC PARAMETERS
# 
##############################################################################
nyears = length(years)


pdf("AQWA_Hungary_model_output_v5.pdf", width=13, height=10)

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
legend(x = 1, y = 700, legend = c("Counts", "Estimated trajectory"), pch = c(4, 16), col = c("red", "blue"), bty = "n")

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
# MAKE A GRAPH OF THE ENVIRONMENTAL EFFECTS
# 
##############################################################################

pdf("AQWA_Hungary_effect_sizes_v5.pdf", width=9, height=7)

par(mar=c(3,5,0,1),oma=c(0,0,0,0))
plot(ipm.model $summary[172:177,1], ylim = c(-3, 3), ylab = "Mean effect size", xlab = "", las = 1, pch = 16, frame = F, cex = 1.5, axes=F, main="", cex.lab=1.5, col='blue')
segments(1:6, ipm.model $summary[172:177,3], 1:6, ipm.model $summary[172:177,7], col="black")
axis(1, at=c(1:6), labels=c("NAO","immigration","rain.fec","fire","inund","win.rain"), cex.axis=1.2)
axis(2, at=seq(-3,3,1), labels=T, las=1, cex.axis=1.3)
abline(h=0, lty=2)

dev.off()





############################################################################
#
# MAKE A GRAPH OF THE CORRELATIONS WITH POP GROWTH RATE
# 
##############################################################################




l.fitted<-l.lower<-l.upper<-ad.fitted<-ad.lower<-ad.upper<-ju.fitted<-ju.lower<-ju.upper<-pr.fitted<-pr.lower<-pr.upper<-imm.fitted<-imm.lower<-imm.upper<-numeric()

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

imm.fitted[i]<-log(quantile(ipm.model $sims.list$imm[,i], 0.5))
imm.lower[i]<-log(quantile(ipm.model $sims.list$imm[,i], 0.025)+1)
imm.upper[i]<-log(quantile(ipm.model $sims.list$imm[,i], 0.975)+1)

pr.fitted[i]<-quantile(ipm.model $sims.list$f[,i], 0.5)
pr.lower[i]<-quantile(ipm.model $sims.list$f[,i], 0.025)
pr.upper[i]<-quantile(ipm.model $sims.list$f[,i], 0.975)}



pdf("AQWA_Hungary_lambda_correlations_v5.pdf", width=9, height=9)

par(mfrow=c(2,2), mar=c(4,5,0,0),oma=c(0,0,0,0))

plot(l.fitted~ad.fitted, xlim=c(0,1), ylim=c(0,2.0), xlab="Adult survival probability",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(ad.lower,l.fitted,ad.upper,l.fitted ,col="gray", lty=1, lwd=0.5)
segments(ad.fitted,l.lower,ad.fitted,l.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(l.fitted,ad.fitted,alternative = c("two.sided"),method = "spearman")
text(0.1,2.0, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)

plot(l.fitted~ju.fitted, xlim=c(0,1), ylim=c(0,2.0), xlab="Juvenile survival probability",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(ju.lower,l.fitted,ju.upper,l.fitted ,col="gray", lty=1, lwd=0.5)
segments(ju.fitted,l.lower,ju.fitted,l.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(l.fitted,ju.fitted,alternative = c("two.sided"),method = "spearman")
text(0,2.0, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)

plot(l.fitted~pr.fitted, xlim=c(0,4), ylim=c(0,2.0), xlab="Fledglings / brood",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(pr.lower,l.fitted,pr.upper,l.fitted ,col="gray", lty=1, lwd=0.5)
segments(pr.fitted,l.lower,pr.fitted,l.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(l.fitted,pr.fitted,alternative = c("two.sided"),method = "spearman")
text(0,2.0, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)

plot(l.fitted~imm.fitted, xlim=c(0,6), ylim=c(0,2.0), xlab="log(number) of immigrants",ylab="Population growth rate",las=1, type='p', pch=16, main="",frame=FALSE, cex.axis=1.3, cex=1, cex.lab=1.3)
segments(imm.lower,l.fitted,imm.upper,l.fitted ,col="gray", lty=1, lwd=0.5)
segments(imm.fitted,l.lower,imm.fitted,l.upper ,col="gray", lty=1, lwd=0.5)
test<-cor.test(l.fitted,imm.fitted,alternative = c("two.sided"),method = "spearman")
text(0,2.0, sprintf("r = %f, p = %g",test$estimate, test$p.value), adj=0)



dev.off()


cor.test(winrain, inund)





