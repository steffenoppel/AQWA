###################################################################################################
########   AQUATIC WARBLER PVA TO EXPLORE REINTRODUCTION TO GERMANY ################
###################################################################################################
# written by Steffen Oppel on 28 October 2021 (steffen.oppel@rspb.org.uk)
# simulation taken from https://naes.unr.edu/shoemaker/teaching/NRES-470/PVA1_421.html
# requested by Susanne Arbeiter to examine needs for reinforcement

# updated on 3 Nov 2024 by adopting AQWA IPM code from 2016
# IPM written without data for population in Hungary
# adopted model structure for equally data-poor situation in Germany

# some modifications from Jaume on 6 November 2024 to:
# Change priors to uniform
# Modify some standard deviations that were generating unrealistic parameter uncertainties
# Add a random effect on productivity

# updated by Jaume on 28 November 2024 to add changes suggested in meeting:
# a past scenario where mowing prevented second broods
# future scenarios where second broods happen again (and another one without mowing just for comparison)
# a modelling of either constrained or unconstrained reproduction by habitat availability
# the release of exactly 50 individuals yearly, instead of a range between 45 and 50 

#The projection scenarios are the following:

#1) Mowing (habitat unconstrained) + 5 release years
#2) Mowing (habitat unconstrained) + 10 release years

#3) No mowing (habitat unconstrained) + 5 release years
#4) No mowing (habitat unconstrained) + 10 release years

#5) No mowing + habitat constrained 180 + 5 release years
#6) No mowing + habitat constrained 180 + 10 release years

#7) No mowing + habitat constrained 360 + 5 release years
#8) No mowing + habitat constrained 360 + 10 release years

library(popbio)
library(doParallel)
library(foreach)
library(tidyverse)
library(data.table)
library(jagsUI)
library(readxl)
filter<-dplyr::filter
select<-dplyr::select


############################################################################
#
# 1. LOAD DATA-------------
# 
##############################################################################
try(setwd("C:/Users/jba/OneDrive - Vogelwarte/Projects/Aquatic Warbler Steffen/Github repo/AQWA"),silent=T)
try(setwd("C:/STEFFEN/OneDrive - Vogelwarte/AQWA"),silent=T)
try(setwd("C:/Users/sop/OneDrive - Vogelwarte/AQWA"),silent=T)

AW<-fread("data/AQWA_counts_GER.csv") %>%
  rename(Year=V1) %>%
  filter(!is.na(Year)) %>%
  filter(Year>2002) # omit the counts from the 1990s that suggest massive decreases
head(AW)
years<-AW$Year

#Number of future scenarios modelled 
nscenarios <- 8

#Number of years for the projections into the future
nprojyears <- 20

#Number of birds released, per scenario and year
releases <- matrix(nrow = nscenarios, ncol = nprojyears, data = 50)
releases[c(1,3,5,7), 6:20] <- 0
releases[c(2,4,6,8), 11:20] <- 0

#A variable to model density-dependence in the number of breeding males that can reproduce. Remember, in the first scenarios, there is no density dependence. Then, we have some with 180 males allowed (120 ha), finally, one with 360 males allowed (240 ha). But we are not really interested in the number of males, but in the number of females. So far I will calculate the number of females inside the model, so that I can consider this extra added uncertainty around the sex ratio. 

maxrepm <- numeric(length = nscenarios)
maxrepm[1:4] <- 10000 #Dummy limit
maxrepm[5:6] <- 180
maxrepm[7:8] <- 360

# Population counts (from years 2003 to 2024)
y <- AW$Nmales		# ENTER DATA OR READ IN AS VECTOR

jags.data <- list(ncountyears = length(years-1), y = y[1:21], nscenarios = nscenarios, nprojyears = nprojyears, releases = releases, maxrepm = maxrepm)

#Caution, we're missing data from year 2024, that leads to an error in the model. Let's add it manually.
jags.data$y[22] <- 3


# Let's define new priors. Perhaps I should make them slightly wider. We have to discuss this. The important thing is, they resemble the prior knowledge.

mean(runif(1e6, 0.28, 0.56)) #Mean is 0.42
pa <- rnorm(1e6, 0.42, 0.03)
hist(pa)
quantile(pa)

mean(runif(1e6, 0.2, 0.44)) #Mean is 0.32
pj <- rnorm(1e6, 0.32, 0.025)
quantile(pj) #Nice

mean(runif(1e6, 2.4, 4)) #Mean is 3.2
f1 <- rnorm(1e6, 3.2, .16)
quantile(f1)

mean(runif(1e6, 1, 2.5)) #Mean is 1.75
f2 <- rnorm(1e6, 1.75, .15)
quantile(f2)

pm <- rnorm(1e6, 0.56, 0.01)
hist(pm)

db <- rnorm(1e6, 0.25, 0.05)
quantile(db)

##############################################################################
#
# IPM WITH INFORMATIVE PRIORS FOR SURVIVAL AND FECUNDITY and releases --------
# 
##############################################################################

sink("models/AQWA.IPM.Mowing.Scenarios.v2.jags")
cat("
model {

#----------------------------------------
# 1. Define the priors for the parameters
#----------------------------------------

# Initial population sizes

N1[1,1] ~ dnorm(60, 0.0001)T(0,)           # 1-year old individuals
NadSurv[1,1] ~ dnorm(100, 0.0001)T(0,)      # Adults >= 2 years

# DEMOGRAPHIC PARAMETERS INFORMED BY DATA FROM OTHER STUDIES
mphia ~ dnorm(0.42,1/(0.03^2))T(0,1)			### survival of adult females
mphij ~ dnorm(0.32,1/(0.025^2))T(0,1)				### survival of first year females
mfec1 ~ dnorm(3.2,1/(0.16^2))			   ### fecundity = number of fledglings raised per FIRST brood
mfec2 ~ dnorm(1.75,1/(0.15^2))			   ### fecundity = number of fledglings raised per SECOND brood, should be slightly lower than first brood
prop.males ~ dnorm(0.56, 1/(0.01^2))T(0,1)  ### proportion of population that is male and can therefore be counted in population - 56%


#--------------------------------------------------
# 2. Random variation in annual survival and productivity
#--------------------------------------------------
# CHANGE THE SCALE OF DEMOGRAPHIC PARAMETERS TO FACILITATE INCORPORATION OF COVARIATES
l.mphij<-log(mphij/(1-mphij))	    # juvenile survival probability on logit scale;
l.mphia<-log(mphia/(1-mphia))			# adult survival probability on logit scale
log.mfec1 <- log(mfec1)            #first-brood productivity on log scale
log.mfec2 <- log(mfec2)           #second-brood productivity on log scale

# RANDOM ANNUAL EFFECTS ON SURVIVAL AND PRODUCTIVITY

tau.surv<-1/pow(sigma.surv,2)
sigma.surv~dunif(0,5)

tau.prod <- pow(sigma.prod, -2)
sigma.prod ~ dunif(0,2)

for (t in 1:(ncountyears+nprojyears)){
   eps.surv[t] ~ dnorm(0, tau.surv)						# random variation around annual survival
   eps.fec[t] ~ dnorm(0, tau.prod)            # random variation around productivity
   logit(phij[t]) <- l.mphij+ eps.surv[t]			# add random effect - could add environmental predictors here (if we had any)
   logit(phia[t]) <- l.mphia+ eps.surv[t]			# add random effect - could add environmental predictors here (if we had any) 
   log(fec1[t]) <- log.mfec1 + eps.fec[t]     # add random effect
}

#Fec 2, the fecundity of second broods, is only applied to future scenarios now

for(t in 1:nprojyears){

  log(fec2[t]) <- log.mfec2 + eps.fec[t+ncountyears] #To use the same eps for both
  
}

#-----------------------
# 3. Derived parameters
#-----------------------

# Population growth rate
for (t in 1:(ncountyears-1)){
   lambda[t] <- Ntot[1,t+1] / max(1,Ntot[1,t]) # prevent invalid parent error when Ntot=0
   loglambda[t]<-log(lambda[t])## for calculating geometric mean of overall population growth rate
}
    #### OVERALL POPULATION GROWTH RATE  #########
    mean.lambda<-exp((1/(ncountyears-1))*sum(loglambda[1:(ncountyears-1)]))   # Geometric mean


#--------------------------------------------
# 4. Likelihood for population count data
#--------------------------------------------
   # 4.1 System process
   for (t in 2:ncountyears){
        
      	chicks[1,t-1] <- (Ntot[1,t-1])*(1-prop.males)* fec1[t-1] # total fecundity
      	chicksrd[1,t-1] <- round(chicks[1,t-1])
		N1[1,t] ~ dbin(phij[t-1],chicksrd[1,t-1]) 
		NadSurv[1,t] ~ dbin(phia[t-1],round((Ntot[1,t-1])))
	}

   # 4.2 Observation process
   for (t in 1:ncountyears){
   	    Ntot[1,t] <- NadSurv[1,t] + N1[1,t]								## total population		
	      Ntotobs[t] <- max(1,(Ntot[1,t]))*prop.males								## only males are counted, which is 56% of population
	      Ntotrd[t] <- round(Ntotobs[t])
        y[t] ~ dpois(Ntotrd[t])
        
   }
      
  #To make calculations easy, let's replicate the calculations for the past into the 8 scenarios for all objects that are involved as well into the future projections
  
  for(ns in 2:nscenarios){
  
    Ntot[ns,1:ncountyears] <- Ntot[1,1:ncountyears]
    N1[ns,1:ncountyears] <- N1[1,1:ncountyears]
    NadSurv[ns,1:ncountyears] <- NadSurv[1,1:ncountyears]
    chicks[ns,1:(ncountyears-1)] <- chicks[1,1:(ncountyears-1)]
    chicksrd[ns,1:(ncountyears-1)] <- chicksrd[1,1:(ncountyears-1)]

  }
      
# -------------------------------------------------        
# 5. PREDICTION INTO THE FUTURE WITH ONGOING RELEASES
# -------------------------------------------------

#Scenarios to be represented

#1) Mowing (habitat unconstrained) + 5 release years
#2) Mowing (habitat unconstrained) + 10 release years

#3) No mowing (habitat unconstrained) + 5 release years
#4) No mowing (habitat unconstrained) + 10 release years

#5) No mowing + habitat constrained 180 + 5 release years
#6) No mowing + habitat constrained 180 + 10 release years

#7) No mowing + habitat constrained 360 + 5 release years
#8) No mowing + habitat constrained 360 + 10 release years


#Hence, there should be variables storing all modifications in parameters involved in these simulations. These are:

#DB and NfemDB should have two dimensions: Nscenario, and time. DB can be manipulated a priori, but not NfemDB

#This should do it for the proportion of double broods

#Scenarios 1 and 2, no double broods
for(ns in 1:2){
  for(t in 1:nprojyears){
  
  db[ns,t] <- 0

  }
}

#The rest, double broods

for(ns in 3:nscenarios){
  for(t in 1:nprojyears){
  
  db[ns,t] ~ dnorm(0.25,1/(0.05^2))T(0,1)

  }
}

#Females allowed to breed in every scenario

maxrepf[1:4] <- maxrepm[1:4]

for(ns in 5:nscenarios){

  maxrepf[ns] <- (maxrepm[ns]*(1-prop.males))/prop.males

}


#Another variable, nintroyears, specifies the number of years for which the birds are released. It is already provided in the data.

### INCLUDE SCENARIOS FOR N CAPTIVE BIRDS AND NUMBER OF YEARS
    
    ## POPULATION PROCESS

    for(ns in 1:nscenarios){
      for (t in (ncountyears+1):(ncountyears+nprojyears)){
        
            #Double broods start in the future now, adapt the index of the parameters to it.
            # changed to proportion of actually breeding females AFTER density dependence
            # Nfemdb[ns,t-ncountyears] ~ dbin(db[ns,t-ncountyears],round((Ntot[ns,t-1])*(1-prop.males)))  ## random draw of double brooding females
            
            #Apply density-dependence
            Nfem.breed1[ns,t-ncountyears] <- min(Ntot[ns,t-1]*(1-prop.males), maxrepf[ns])
            Nfem.breed2[ns,t-ncountyears] ~ dbin(db[ns,t-ncountyears],round(Nfem.breed1[ns,t-ncountyears]))  

      	    chicks[ns,t-1] <- Nfem.breed1[ns,t-ncountyears] * fec1[t] + Nfem.breed2[ns,t-ncountyears]*fec2[t-ncountyears]			# total fecundity
      	    chicksrd[ns,t-1] <- round(chicks[ns,t-1]) + releases[ns,t-ncountyears]
		        N1[ns,t] ~ dbin(phij[t-1],chicksrd[ns,t-1]) 
		        NadSurv[ns,t] ~ dbin(phia[t-1],round((Ntot[ns,t-1])))
		   	    Ntot[ns,t] <- NadSurv[ns,t] + N1[ns,t]								## total population	
  	} #Time loop
  } #Nscenarios loop
}
",fill = TRUE)
sink()



# Initial values
inits <- function(){list(
  mphij= runif(1, 0.25,0.35),
  mphia = runif(1, 0.45,0.55),
  mfec1 = runif(1, 2.5,4.0),
  mfec2 = runif(1, 1,2.5))}

# Parameters monitored
parameters <- c("Ntot","mphij","mphia","mfec1","mfec2","mean.lambda","prop.males", "phij", "phia", "fec1", "fec2", "db", "Nfemdb")


# MCMC settings
ni <- 75000
nt <- 5
nb <- 25000
nc <- 4


# Call JAGS from R
ipm.model <- jags(jags.data,
                  inits,
                  parameters,
                  "models/AQWA.IPM.Mowing.Scenarios.v2.jags",
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.cores=nc, parallel=T)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure 1: POPULATION TRAJECTORY IN PAST AND FUTURE BY SCENARIO ----------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ntotdf <- data.frame(Ntot = numeric(length = length(ipm.model$mean$Ntot)), cip = numeric(length = length(ipm.model$mean$Ntot)), cim = numeric(length = length(ipm.model$mean$Ntot)), Scenario = c(rep("Mowing + 5y", 42), rep("Mowing + 10y", 42), rep("No mowing + 5y Habitat Unconstrained", 42), rep("No mowing + 10y Habitat Unconstrained", 42), rep("No mowing + 5y Habitat Constrained 120 ha", 42), rep("No mowing + 10y Habitat Constrained 120 ha", 42), rep("No mowing + 5y Habitat Constrained 240 ha", 42), rep("No mowing + 10y Habitat Constrained 240 ha", 42)), year = rep(1:42, nscenarios))

ntotdf$Ntot <- as.numeric(t(ipm.model$mean$Ntot))
ntotdf$cim <- as.numeric(t(ipm.model$q2.5$Ntot))
ntotdf$cip <- as.numeric(t(ipm.model$q97.5$Ntot))

ntotdf$Scenario <- factor(ntotdf$Scenario, levels = c("Mowing + 5y", "Mowing + 10y", "No mowing + 5y Habitat Unconstrained", "No mowing + 10y Habitat Unconstrained", "No mowing + 5y Habitat Constrained 120 ha", "No mowing + 10y Habitat Constrained 120 ha", "No mowing + 5y Habitat Constrained 240 ha", "No mowing + 10y Habitat Constrained 240 ha"))

ntotdf$year <- factor(ntotdf$year)

#The plot

poptrends <- ggplot(data = ntotdf, aes(x = year, y = Ntot, col = Scenario, group = Scenario)) + geom_line(size = 1.1) + geom_vline(aes(xintercept = 23), linetype = 2) + geom_ribbon(aes(ymin = cim, ymax = cip), linetype = 2, fill = NA) + facet_wrap(~Scenario, nrow = 4)  + theme_bw() + theme(legend.position = "none") + scale_x_discrete(breaks = seq(from = 0, to = 42, by = 6)) + xlab("") + ylab("")

#ggsave("output/Scenario_projections.jpg", width=181,height=141, quality=100, units="mm")


#Extinction probability (could plot this as well)

w0 <- function(x){
  
  ext <- sum(x == 0) / length(x)
  return(ext)
  
}

ny <- jags.data$ncountyears

extprob <- matrix(nrow = nscenarios, ncol = nprojyears, data = 0)

#Quick plot, has to be improved in case we want to keep it

par(mfrow = c(4,2))
for(i in 1:nscenarios){
  for(t in 1:nprojyears){
    extprob[i,t] <- w0(ipm.model$sims.list$Ntot[,i,t+jags.data$ncountyears])
  }
  plot(extprob[i,], type = "l", lwd = 1.9, ylim = c(0,1), ylab = "P(ext)", xlab = "Year", main = paste0("Scenario ", i))
}

#It looks like extinction probabilities 20 years into the future are not large. What about population decline probabilities? (That is, what is the probability that the population size at the end of the projections is lower than that at the end of the releases, in each scenario?)

pdecline <- numeric(length = nscenarios)
relscenario <- rep(c(5,10),4)

for(i in 1:nscenarios){
  
  pdecline[i] = sum(ipm.model$sims.list$Ntot[,i,jags.data$ncountyears+nprojyears] < ipm.model$sims.list$Ntot[,i,(jags.data$ncountyears+relscenario[i])]) / dim(ipm.model$sims.list$Ntot)[1]
  
}

pdecline #Quite large

#On the retrospective analysis (past trends), see how inference has modified the priors and the estimates.

par(mfrow = c(2,2))

#Juvenile survival

plot(ipm.model$mean$phij, type = "l", col = "green", ylim = c(0.18, 0.45), lwd = 2, xlab = "Year", ylab = "Estimate", main = "Juvenile survival")
lines(ipm.model$q2.5$phij, lty = "dashed", lwd = 1.5, col = "green")
lines(ipm.model$q97.5$phij, lty = "dashed", lwd = 1.5, col = "green")
abline(v = ny, lty = "dashed")

#Adult survival

plot(ipm.model$mean$phia, type = "l", col = "green", ylim = c(0.25,0.6), lwd = 2, xlab = "Year", ylab = "Estimate", main = "Adult survival")
lines(ipm.model$q2.5$phia, lty = "dashed", lwd = 1.5, col = "green")
lines(ipm.model$q97.5$phia, lty = "dashed", lwd = 1.5, col = "green")
abline(v = ny, lty = "dashed")

#First-year productivity

plot(ipm.model$mean$fec1, type = "l", col = "red", ylim = c(1.5,5.6), lwd = 2, xlab = "Year", ylab = "Estimate", main = "First-year productivity")
lines(ipm.model$q2.5$fec1, lty = "dashed", lwd = 1.5, col = "red")
lines(ipm.model$q97.5$fec1, lty = "dashed", lwd = 1.5, col = "red")
abline(v = ny, lty = "dashed")

#Second-year productivity

plot(ipm.model$mean$fec2, type = "l", col = "red", ylim = c(0.7,3.2), lwd = 2, xlab = "Year", ylab = "Estimate", main = "Second-year productivity")
lines(ipm.model$q2.5$fec2, lty = "dashed", lwd = 1.5, col = "red")
lines(ipm.model$q97.5$fec2, lty = "dashed", lwd = 1.5, col = "red")
abline(v = ny, lty = "dashed")

#See how priors compare to final estimates

par(mfrow = c(2,2))

#Adult survival

adprior <- dnorm(seq(from = 0, to = 1, by = 0.01), mean = 0.42, sd = 0.03)
plot(adprior, type = "l", xaxt = "n", ylab = "", xlab = "Estimate value", main = "Adult Survival \n Red = Prior, Blue = Estimate", ylim = c(0,15), col = "red", lwd = 2)
axis(1, at = seq(from = 0, to = 100, by = 10), labels = seq(from = 0, to = 1, by = 0.1))
ipm.model$mean$mphia
ipm.model$sd$mphia
adresult <- dnorm(seq(from = 0, to = 1, by = 0.01), mean = ipm.model$mean$mphia, sd = ipm.model$sd$mphia)
lines(adresult, type = "l", lty = "dashed", col = "blue", lwd = 2)

#Juvenile survival

juvprior <- dnorm(seq(from = 0, to = 1, by = 0.01), mean = 0.32, sd = 0.025)
plot(juvprior, type = "l", xaxt = "n", ylab = "", xlab = "Estimate value", main = "Juvenile Survival \n Red = Prior, Blue = Estimate", ylim = c(0,19), col = "red", lwd = 2)
axis(1, at = seq(from = 0, to = 100, by = 10), labels = seq(from = 0, to = 1, by = 0.1))
ipm.model$mean$mphij
ipm.model$sd$mphij
juvresult <- dnorm(seq(from = 0, to = 1, by = 0.01), mean = ipm.model$mean$mphij, sd = ipm.model$sd$mphij)
lines(juvresult, type = "l", lty = "dashed", col = "blue", lwd = 2)

#First-brood productivity

p1prior <- dnorm(seq(from = 2, to = 4, by = 0.1), mean = 3.2, sd = 0.16)
plot(p1prior, type = "l", xaxt = "n", ylab = "", xlab = "Estimate value", main = "First-brood productivity \n Red = Prior, Blue = Estimate", ylim = c(0,3), col = "red", lwd = 2)
axis(1, at = seq(from = 1, to = 21, by = 1), labels = seq(from = 2, to = 4, by = 0.1))
ipm.model$mean$mfec1
ipm.model$sd$mfec1
p1result <- dnorm(seq(from = 2, to = 4, by = 0.1), mean = ipm.model$mean$mfec1, sd = ipm.model$sd$mfec1)
lines(p1result, type = "l", lty = "dashed", col = "blue", lwd = 2)

#Second-brood productivity

p2prior <- dnorm(seq(from = 0.5, to = 3.2, by = 0.1), mean = 1.75, sd = 0.15)
plot(p2prior, type = "l", xaxt = "n", ylab = "", xlab = "Estimate value", main = "Second-brood productivity \n Red = Prior, Blue = Estimate", ylim = c(0,3), col = "red", lwd = 2)
axis(1, at = seq(from = 1, to = 28, by = 1), labels = seq(from = 0.5, to = 3.2, by = 0.1))
ipm.model$mean$mfec2
ipm.model$sd$mfec2
p2result <- dnorm(seq(from = 0.5, to = 3.2, by = 0.1), mean = ipm.model$mean$mfec2, sd = ipm.model$sd$mfec2)
lines(p2result, type = "l", lty = "dashed", col = "blue", lwd = 2)

#Compared to models where we allowed second broods, here the change in the estimates is very very mild.

############################################################################
#
# WRITE ALL OUTPUT INTO A TEXT FILE ----
# 
##############################################################################
out<-as.data.frame(ipm.model$summary)  
out$parameter<-row.names(ipm.model$summary)
names(out)[c(12,5,3,7)]<-c('parm','median','lcl','ucl')
print(ipm.model, dig=3)
write.table(out, "output/AQWA_GER_model_JAUME_nscenarios_output.csv", sep=",")




