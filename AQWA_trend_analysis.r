#######################################################################################################################################
#################################### AQUATIC WARBLER MONITORING DATA ANALYSIS 2026 ####################################################
#######################################################################################################################################
# written by steffen.oppel@rspb.org.uk
# TRANSECT B35 was only surveyed in 2012!
# 17 May 2012: updated script to remove veg1 (redundant in combination with veg2+veg3+veg4)!
# 7 August 2013: reverted to BMM (rather than gdistsamp) due to reviewer criticism
# 25 Aug 2013: removed beta-binomial mixture models, as simple BMM provided perfect fit for data even without covariates
# 27 August: ran a varieties of site-cov models, but DIC is always lowest for the one incl veg2, veg3, and veg4
# 7 October 2013: tried to implement Adam Butler's suggestions for incorporating power analysis in BUGS code -> TRAP UNDEFINED REAL RESULT
# 9 October 2013: updated model 3 to include wat and veg as abund and det level covariates
# 12 Dec 2013: updated model to include day instead of temp after repeated model selection un unmarked 
## RESULT: DIC is LOWEST for model 3 with wat and veg as covariates (model 1: DIC = 3048.9, model 2: DIC = 2873.6, model 3: DIC = 2179.3)



############### LOAD THE REQUIRED PACKAGES #########################################################

library(jagsUI)



############### SET THE WORKING DIRECTORY #########################################################

try(setwd("C:/Users/jba/OneDrive - Vogelwarte/Projects/Aquatic Warbler Steffen/Github repo/AQWA"),silent=T)
try(setwd("C:/STEFFEN/OneDrive - Vogelwarte/ExternalCollaborations/AQWA"),silent=T)
try(setwd("C:/Users/sop/OneDrive - Vogelwarte/ExternalCollaborations/AQWA"),silent=T)


######################################################################################################
# 
# 1. READ IN THE DATA THAT HAVE BEEN FORMATTED FOR JAGS
# 
#######################################################################################################

load("data\\AQWA_trend_input.RData")

######################################################################################################
# 
# DATA SCREENING AND SIMPLE SUMMARY
# 
######################################################################################################
y<-AQWA.y
head(y)


######################################################################################################
# 
# MODEL 1: Simple Binomial-mixture model for trend estimation as described by K?ry et al. 2009 - no covariates
# 
######################################################################################################
# Specify model in BUGS language
sink("models/NmixTrend.jags")
cat("
model{

# Priors
loglam~dunif(-5,5)					##	mean abundance prior
trend~dunif(-10,10)					##    trend prior

for(i in 1:nsite){
lam.site[i]~dnorm(0,tau.site)T(-20, 20)		## needs to be T(-20, 20) in JAGS
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
log(lambda[i,year])<- loglam + trend*primocc[year]+lam.site[i]		## 'primocc' is the primary occasion - here: years
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
fit <- sum(E2[,,])		# Sum up squared residuals for actual data set
fit.new <- sum(E2.new[,,]) 	# Sum up for replicate data sets
}
",fill = TRUE)
sink()


# Bundle data
R = nrow(y)
T = ncol(y)
nyears = dim(y)[3]
jags.data <- list(M = y, nsite = R, nrep = T, nyear=nyears, primocc=seq(1:nyears))

# Initial values
Nst <- apply(y, c(1, 3), max) + 1
Nst[is.na(Nst)] <- 1
inits <- function(){list(N = Nst, loglam = runif(1, -3, 3), sigma.site = runif(1, 0, 1), p0 = runif(nyears, 0, 1))}

# Parameters monitored
parameters <- c("loglam","lam.site","totalN","N", "fit", "fit.new","trend")

# MCMC settings
ni <- 500000				## fit statistic needs >400000 to converge!!
nt <- 40
nb <- 100000
nc <- 3


# Call JAGS from R
trend.model <- jags(jags.data,
                  inits,
                  parameters,
                  "models/NmixTrend.jags",
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.cores=nc, parallel=T)

# Evaluation of fit
plot(model1$sims.list$fit, model1$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE, xlim = c(100, 300), ylim = c(100, 300))
abline(0, 1, lwd = 2, col = "black")
mean(model1$sims.list$fit.new > model1$sims.list$fit)
mean(model1$mean$fit) / mean(model1$mean$fit.new)

# Summarize posteriors
print(model1, dig = 2)



######################################################################################################
# 
# MODEL 2: Binomial-mixture model for trend estimation INCORPORATING COVARIATES FOR ABUNDANCE (but not detection)
# 
######################################################################################################


# Specify model in BUGS language
sink("BinmixSiteCovariateTrend.txt")
cat("
model{

# Priors
loglam~dunif(-5,5)					##	mean abundance prior
trend~dunif(-10,10)					##    trend prior
#beta.lit2~dunif(-5,5)
#beta.lit3~dunif(-5,5)
#beta.lit4~dunif(-5,5)
#bveg1~dunif(-5,5)
bveg2~dunif(-5,5)
bveg3~dunif(-5,5)
bveg4~dunif(-5,5)

for(i in 1:nsite){
lam.site[i]~dnorm(0,tau.site)	##I(-20, 20)		## site-specific random effect
}
tau.site<-1/(sigma.site*sigma.site)
sigma.site~dunif(0,10)

for(year in 1:nyear){
p0[year]~dunif(0,1)					## detection probability
logitp0[year]<-log(p0[year]/(1-p0[year]))
}
tau.lp<-1/(sigma.p*sigma.p)
sigma.p~dunif(0,10)


# State and observation models
for(year in 1:nyear){
for(i in 1:nsite){
log(lambda[i,year])<- loglam + trend*primocc[year]+bveg2*veg2[i,year]+ bveg3*veg3[i,year]+ bveg4*veg4[i,year]+lam.site[i]		##  'primocc' is the primary occasion - here: years
N[i,year]~dpois(lambda[i,year])

for(t in 1:nrep){
M[i,t,year]~dbin(p[i,t,year],N[i,year])
p[i,t,year] <- exp(lp[i,t,year])/(1+exp(lp[i,t,year]))
lp[i,t,year] ~ dnorm(mu.lp[i,t,year], tau.lp)			##I(-20, 20) # H.c.
mu.lp[i,t,year]<-logitp0[year]
}
}

### DERIVED PARAMETER FOR EACH YEAR ###
totalN[year]<-sum(N[,year])
det.prob[year]<-mean(p[,,year])
}

# Computation of fit statistic (Bayesian p-value)
# Fit statistic for observed data
# Also, generate replicate data and compute fit stats for them
for(year in 1:nyear){


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
fit <- sum(E2[,,])		# Sum up squared residuals for actual data set
fit.new <- sum(E2.new[,,]) 	# Sum up for replicate data sets
}
",fill = TRUE)
sink()


# Bundle data
R = nrow(y)
T = ncol(y)
nyears = dim(y)[3]
win.data <- list(M = y, nsite = R, nrep = T, nyear=nyears, primocc=seq(1:nyears), veg2=veg2, veg3=veg3, veg4=veg4)

# Initial values
Nst <- apply(y, c(1, 3), max) + 1
Nst[is.na(Nst)] <- 1
inits <- function(){list(N = Nst,
loglam = runif(1, -3, 3),
sigma.site = runif(1, 0, 1),
#beta.lit2=runif(1,-5,5),
#beta.lit3=runif(1,-5,5),
#beta.lit4=runif(1,-5,5),
bveg2=runif(1,-5,5),
bveg3=runif(1,-5,5),
bveg4=runif(1,-5,5),
#btemp=runif(1,-5,5),
p0 = runif(nyears, 0, 1))}

# Parameters monitored
params <- c("lam.site","loglam", "trend","totalN","bveg2","bveg3","bveg4","fit", "fit.new","det.prob")


# MCMC settings
ni <- 250000
nt <- 5
nb <- 50000
nc <- 3

# Call WinBUGS from R
#model2 <- bugs(win.data, inits, params, "BinmixSiteCovariateTrend.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
model2 <- jags(win.data, inits, params, "BinmixSiteCovariateTrend.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
model2<-model2$BUGSoutput

# Evaluation of fit
plot(model2$sims.list$fit, model2$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE, xlim = c(100, 300), ylim = c(100, 300))
abline(0, 1, lwd = 2, col = "black")
mean(model2$sims.list$fit.new > model2$sims.list$fit)
mean(model2$mean$fit) / mean(model2$mean$fit.new)

# Summarize posteriors
print(model2, dig = 2)



######################################################################################################
# 
# MODEL 3: Binomial-mixture model for trend estimation INCORPORATING COVARIATES FOR ABUNDANCE AND DETECTION
# 
######################################################################################################


# Specify model in BUGS language
sink("BinmixAllCovariateTrend.txt")
cat("

model{

# Priors
loglam~dunif(-5,5)##mean abundance prior
trend~dunif(-10,10)##    trend prior
beta.wat2~dunif(-5,5)
beta.wat3~dunif(-5,5)
beta.wat4~dunif(-5,5)
bveg2~dunif(-5,5)
bveg3~dunif(-5,5)
bveg4~dunif(-5,5)
brain~dunif(-5,5)

for(i in 1:nsite){
lam.site[i]~dnorm(0,tau.site)##I(-20, 20)## site-specific random effect
}
tau.site<-1/(sigma.site*sigma.site)
sigma.site~dunif(0,10)

for(year in 1:nyear){
p0[year]~dunif(0,1)## detection probability
logitp0[year]<-log(p0[year]/(1-p0[year]))
}
tau.lp<-1/(sigma.p*sigma.p)
sigma.p~dunif(0,10)


# State and observation models
for(year in 1:nyear){
for(i in 1:nsite){
log(lambda[i,year])<- loglam + trend*primocc[year]+beta.wat2*wat2[i,year]+beta.wat3*wat3[i,year]+beta.wat4*wat4[i,year]+lam.site[i]## 'primocc' is the primary occasion - here: years
N[i,year]~dpois(lambda[i,year])

for(t in 1:nrep){
M[i,t,year]~dbin(p[i,t,year],N[i,year])
p[i,t,year] <- exp(lp[i,t,year])/(1+exp(lp[i,t,year]))
lp[i,t,year] ~ dnorm(mu.lp[i,t,year], tau.lp)##I(-20, 20) # H.c.
mu.lp[i,t,year]<-logitp0[year] + brain*rain[i,t,year]+ bveg2*veg2[i,year]+ bveg3*veg3[i,year]+ bveg4*veg4[i,year]
}
}

### DERIVED PARAMETER FOR EACH YEAR ###
totalN[year]<-sum(N[,year])
}

# Computation of fit statistic (Bayesian p-value)
# Fit statistic for observed data
# Also, generate replicate data and compute fit stats for them
for(year in 1:nyear){


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
",fill = TRUE)
sink()


# Bundle data
R = nrow(y)
T = ncol(y)
nyears = dim(y)[3]
win.data <- list(M = y, nsite = R, nrep = T, nyear=nyears, primocc=seq(1:nyears), wat2=wat2, wat3=wat3, wat4=wat4, rain=rain,veg2=veg2,veg3=veg3,veg4=veg4)

# Initial values
Nst <- apply(y, c(1, 3), max) + 1
Nst[is.na(Nst)] <- 1
inits <- function(){list(N = Nst,
loglam = runif(1, -3, 3),
sigma.site = runif(1, 0, 1),
beta.wat2=runif(1,-5,5),
beta.wat3=runif(1,-5,5),
beta.wat4=runif(1,-5,5),
bveg2=runif(1,-5,5),
bveg3=runif(1,-5,5),
bveg4=runif(1,-5,5),
brain=runif(1,-5,5),
p0 = runif(nyears, 0, 1))}

# Parameters monitored
params <- c("loglam","lam.site","trend","totalN","beta.wat2","beta.wat3","beta.wat4","bveg2","bveg3","bveg4","brain","fit", "fit.new")

# MCMC settings
ni <- 750000
nt <- 5
nb <- 150000
nc <- 3

# Call WinBUGS from R (BRT 6000 seconds for 350000 iterations)
#model3 <- bugs(win.data, inits, params, "BinmixAllCovariateTrend.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
model3 <- jags(win.data, inits, params, "BinmixAllCovariateTrend.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
model3<-model3$BUGSoutput

# Evaluation of fit
plot(model3$sims.list$fit, model3$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE, xlim = c(100, 300), ylim = c(100, 300))
abline(0, 1, lwd = 2, col = "black")
mean(model3$sims.list$fit.new > model3$sims.list$fit)
mean(model3$mean$fit) / mean(model3$mean$fit.new)

# Summarize posteriors
print(model3, dig = 2)
write.table(model3$summary,"S:\\ConSci\\DptShare\\SteffenOppel\\MANUSCRIPTS\\in_prep\\AQWA_census\\BCI\\final_model_out.csv", sep=",")



######################################################################################################
# 
# MODEL 4: Binomial-mixture model for trend estimation INCORPORATING COVARIATES FOR DETECTION BUT NOT FOR ABUNDANCE
# 
######################################################################################################


# Specify model in BUGS language
sink("BinmixDetCovariateTrend.txt")
cat("

model{

# Priors
loglam~dunif(-5,5)##mean abundance prior
trend~dunif(-10,10)##    trend prior
btemp~dunif(-5,5)
bveg2~dunif(-5,5)
bveg3~dunif(-5,5)
bveg4~dunif(-5,5)

for(i in 1:nsite){
lam.site[i]~dnorm(0,tau.site)##I(-20, 20)## site-specific random effect
}
tau.site<-1/(sigma.site*sigma.site)
sigma.site~dunif(0,10)

for(year in 1:nyear){
p0[year]~dunif(0,1)## detection probability
logitp0[year]<-log(p0[year]/(1-p0[year]))
}
tau.lp<-1/(sigma.p*sigma.p)
sigma.p~dunif(0,10)


# State and observation models
for(year in 1:nyear){
for(i in 1:nsite){
log(lambda[i,year])<- loglam + trend*primocc[year]+lam.site[i]## 'primocc' is the primary occasion - here: years
N[i,year]~dpois(lambda[i,year])

for(t in 1:nrep){
M[i,t,year]~dbin(p[i,t,year],N[i,year])
p[i,t,year] <- exp(lp[i,t,year])/(1+exp(lp[i,t,year]))
lp[i,t,year] ~ dnorm(mu.lp[i,t,year], tau.lp)##I(-20, 20) # H.c.
mu.lp[i,t,year]<-logitp0[year] + btemp*temp[i,t,year]+ bveg2*veg2[i,year]+ bveg3*veg3[i,year]+ bveg4*veg4[i,year]
}
}

### DERIVED PARAMETER FOR EACH YEAR ###
totalN[year]<-sum(N[,year])
}

# Computation of fit statistic (Bayesian p-value)
# Fit statistic for observed data
# Also, generate replicate data and compute fit stats for them
for(year in 1:nyear){


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
",fill = TRUE)
sink()


# Bundle data
R = nrow(y)
T = ncol(y)
nyears = dim(y)[3]
win.data <- list(M = y, nsite = R, nrep = T, nyear=nyears, primocc=seq(1:nyears), temp=day,veg2=veg2,veg3=veg3,veg4=veg4)

# Initial values
Nst <- apply(y, c(1, 3), max) + 1
Nst[is.na(Nst)] <- 1
inits <- function(){list(N = Nst,
loglam = runif(1, -3, 3),
sigma.site = runif(1, 0, 1),
bveg2=runif(1,-5,5),
bveg3=runif(1,-5,5),
bveg4=runif(1,-5,5),
btemp=runif(1,-5,5),
p0 = runif(nyears, 0, 1))}

# Parameters monitored
params <- c("lam.site","loglam", "trend","totalN","bveg2","bveg3","bveg4","btemp","fit", "fit.new")

# MCMC settings
ni <- 250000
nt <- 5
nb <- 100000
nc <- 3

# Call WinBUGS from R (BRT 6000 seconds for 350000 iterations)
model4 <- jags(win.data, inits, params, "BinmixDetCovariateTrend.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
model4<-model4$BUGSoutput

# Evaluation of fit
plot(model4$sims.list$fit, model4$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE, xlim = c(100, 300), ylim = c(100, 300))
abline(0, 1, lwd = 2, col = "black")
mean(model4$sims.list$fit.new > model4$sims.list$fit)
mean(model4$mean$fit) / mean(model4$mean$fit.new)

# Summarize posteriors
print(model4, dig = 2)





######################################################################################################
# 
# MODEL 5: Combined State-space and Binomial-mixture model for trend estimation INCORPORATING COVARIATES FOR ABUNDANCE (but not detection)
# 
######################################################################################################


# Specify model in BUGS language
sink("Binmix_SSM.txt")
cat("
model{

# Priors
loglam~dunif(-5,5)					##	mean abundance prior
#trend~dunif(-10,10)					##    trend prior
bveg2~dunif(-5,5)
bveg3~dunif(-5,5)
bveg4~dunif(-5,5)
tau.proc<-1/(sigma.proc*sigma.proc)		## site-specific temporal process variation
sigma.proc~dunif(0,10)
mean.trend~dunif(0.1,10)					## site-specific trend prior
totalN[1]~dunif(50,5000)						## prior for first year total abundance across all transects


for(i in 1:nsite){
lam.site[i]~dnorm(0,tau.site)					## site-specific random effect
}
tau.site<-1/(sigma.site*sigma.site)
sigma.site~dunif(0,10)

for(year in 1:nyear){
p0[year]~dunif(0,1)					## detection probability
logitp0[year]<-log(p0[year]/(1-p0[year]))
}
tau.lp<-1/(sigma.p*sigma.p)
sigma.p~dunif(0,10)


# State and observation models
for(year in 1:nyear){
for(i in 1:nsite){
log(lambda[i,year])<- loglam + bveg2*veg2[i,year]+ bveg3*veg3[i,year]+ bveg4*veg4[i,year]+lam.site[i]		## trend estimation moved to SSM likelihood
N[i,year]~dpois(lambda[i,year])

for(t in 1:nrep){
M[i,t,year]~dbin(p[i,t,year],N[i,year])
p[i,t,year] <- exp(lp[i,t,year])/(1+exp(lp[i,t,year]))
lp[i,t,year] ~ dnorm(mu.lp[i,t,year], tau.lp)
mu.lp[i,t,year]<-logitp0[year]
}			## close loop over T replicate counts
}			## close loop over sites

### DERIVED PARAMETER FOR EACH YEAR ###
#totalN[year]<-sum(N[,year])
detect.prob[year]<-mean(p[,,year])
}			## close loop over years



### STATE-SPACE MODEL FOR TREND ESTIMATION ###
for(year in 1:(nyear-1)){			### start loop over every year
totalN[year+1]<- sum(N[,year])*trend[year]			###
trend[year]~dnorm(mean.trend, tau.proc)
}


# Computation of fit statistic (Bayesian p-value)
# Fit statistic for observed data
# Also, generate replicate data and compute fit stats for them
for(year in 1:nyear){


for(i in 1:nsite){
for(t in 1:nrep){

# Actual data
eval[i,t,year] <-N[i,year]*p[i,t,year] # Expected value
sd.resi[i,t,year]<-sqrt(eval[i,t,year]*(1-p[i,t,year])) +0.5
E[i,t,year]<-(M[i,t,year]-eval[i,t,year])/ sd.resi[i,t,year]
E2[i,t,year] <- pow(E[i,t,year],2)

# Replicate data sets
M.new[i,t,year]~dbin(p[i,t,year],N[i,year])			### may have to truncate to I(0,
E.new[i,t,year]<-(M.new[i,t,year]-eval[i,t,year])/sd.resi[i,t,year]
E2.new[i,t,year] <- pow(E.new[i,t,year], 2)
}
}
}
fit <- sum(E2[,,])		# Sum up squared residuals for actual data set
fit.new <- sum(E2.new[,,]) 	# Sum up for replicate data sets
}
",fill = TRUE)
sink()


# Bundle data
R = nrow(y)
T = ncol(y)
nyears = dim(y)[3]
win.data <- list(M = y, nsite = R, nrep = T, nyear=nyears, veg2=veg2,veg3=veg3,veg4=veg4)


# Initial values
Nst <- apply(y, c(1, 3), max) + 1
Nst[is.na(Nst)] <- 20	### manual fix for the transect B35 that was not surveyed in 2011
inits <- function(){list(
N = Nst,
loglam = runif(1, -3, 3),
sigma.site = runif(1, 0, 1),
bveg2=runif(1,-5,5),
bveg3=runif(1,-5,5),
bveg4=runif(1,-5,5),
p0 = runif(nyears, 0, 1))}

# Parameters monitored
params <- c("mean.trend","detect.prob","totalN","bveg2","bveg3","bveg4","fit", "fit.new","loglam","tau.proc","sigma.proc")

# MCMC settings
ni <- 35000
nt <- 5
nb <- 15000
nc <- 3

# Call WinBUGS from R
#model5 <- bugs(win.data, inits, params, "Binmix_SSM.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
model5 <- jags(win.data, inits, params, "Binmix_SSM.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
model5<-model5$BUGSoutput

# Evaluation of fit
plot(model5$sims.list$fit, model5$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE, xlim = c(100, 300), ylim = c(100, 300))
abline(0, 1, lwd = 2, col = "black")
mean(model5$sims.list$fit.new > model5$sims.list$fit)
mean(model5$mean$fit) / mean(model5$mean$fit.new)

# Summarize posteriors
print(model5, dig = 2)





######################################################################################################
# 
# EXPORT SUMMARIES AND SAVE OUTPUT
# 
######################################################################################################

sink("AQWA_BMM_trend_model_summaries.txt")
print(model1,dig = 2)
print(model2,dig = 2)
print(model3, dig = 4)
print(model4, dig = 2)
sink()

save.image("AQWA_BMM_BUGS_trend_model.RData")




######################################################################################################
# 
# COMPARE RAW COUNT AND ESTIMATED ABUNDANCE 
# 
######################################################################################################
count2011<-apply(y[,,1],2,sum,na.rm=T)
count2012<-apply(y[,,2],2,sum,na.rm=T)
count2013<-apply(y[,,3],2,sum,na.rm=T)


model3$summary[62,1]/649				## for max count in 2011
model3$summary[62,3]/649				## for max count in 2011
model3$summary[62,7]/649				## for max count in 2011

model3$summary[63,1]/608				## for max count in 2012
model3$summary[63,3]/608				## for max count in 2012
model3$summary[63,7]/608				## for max count in 2012

model3$summary[64,1]/622				## for max count in 2013
model3$summary[64,3]/622				## for max count in 2013
model3$summary[64,7]/622				## for max count in 2013



########################################################################################
### ESTIMATE RATIO OF FULL COUNT VS MODEL ESTIMATE AND OVERALL DENSITY   ############
########################################################################################
mean(c(model3$summary[62,1]/649,model3$summary[63,1]/608,model3$summary[64,1]/622)) *2594
mean(c(model3$summary[62,3]/649,model3$summary[63,3]/608,model3$summary[64,3]/622)) *2594
mean(c(model3$summary[62,7]/649,model3$summary[63,7]/608,model3$summary[64,7]/622)) *2594

(model3$summary[63,1]/608) *2594
(model3$summary[63,3]/608)*2594
(model3$summary[63,7]/608)*2594


