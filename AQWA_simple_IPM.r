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


library(jagsUI)


############################################################################
#
# LOAD DATA
# 
##############################################################################

years<-seq(1976,2010,1)


# Population counts (from years 1976 to 2010)
y <- c(68,60,62,88,104,122,170,208,176,188,198,214,224,220)		# ENTER DATA OR READ IN AS VECTOR

# NAO data
NAO <- c(68,60,62,88,104,122,170,208,176,188,198,214,224,220)		# ENTER DATA OR READ IN AS VECTOR

# Precipitation data
rain <- c(68,60,62,88,104,122,170,208,176,188,198,214,224,220)		# ENTER DATA OR READ IN AS VECTOR





##############################################################################
#
# IPM WITH FIXED ENVIRONMENTAL EFFECTS ON SURVIVAL AND FECUNDITY
# 
##############################################################################

sink("AQWA.IPM.simple.jags")
cat("
model {



#----------------------------------------
# 1. Define the priors for the parameters
#----------------------------------------

# Initial population sizes
N1[1] ~ dnorm(20, 0.0001)T(0,)           # 1-year old individuals
NadSurv[1] ~ dnorm(60, 0.0001)T(0,)      # Adults >= 2 years

# DEMOGRAPHIC PARAMETERS INFORMED BY DATA FROM OTHER STUDIES
mphij ~ dunif(0.05,0.85)       			# mean juvenile survival, UPDATE WITH PROPER RANGE
mphia ~ dunif(0.10,0.85)			# mean adult survival, UPDATE WITH PROPER RANGE
mfec ~ dunif(1.1,3.7)			# mean fecundity per breeding attempt, derived from Kubacka et al. 2013 by multiplying brood size and nest survival (3.8–4.1)*(0.36-0.87)
nbroods ~ dunif(1,2)			# mean number of broods raised by the population

#
l.mphij<-log(mphij/(1-mphij))			# juvenile survival probability on logit scale
l.mphia<-log(mphia/(1-mphia))			# adult survival probability on logit scale


# PRIORS FOR ENVIRONMENTAL EFFECTS

beta.f ~ dnorm(0, 0.0001)T(-10,10)		# prior for rain effect on fecundity
beta.phi ~ dnorm(0, 0.0001)T(-10,10)		# prior for NAO effect on survival
beta.brood ~ dnorm(0, 0.0001)T(-10,10)		# prior for rain effect on number of broods




#--------------------------------------------------
# 2. Relate parameters to environmental variables
#--------------------------------------------------

for (t in 1:(nyears-1)){

   logit(phij[t]) <- l.mphij + beta.phi*NAO[t]  # Juv. apparent survival
   logit(phia[t]) <- l.mphia + beta.phi*NAO[t]  # Adult apparent survival
   f[t] <- mfec  + beta.f*rain[t]         	# Productivity per breeding attempt
   nbrood[t] <- nbroods  + beta.brood*rain[t]   # Number of broods

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
      	chicks[t] <- 0.5 * f[t-1] * nbrood[t-1] * Ntot[t-1]		## nbroods is avg number of broods per year, f is fecundity per clutch, not per year
      	N1[t] ~ dbin(phij[t-1],chicks[t-1]) 
	NadSurv[t] ~ dbin(phia[t-1],Ntotrd[t-1]) 		

   # 4.2 Observation process
   for (t in 1:nyears){
	capacity[t] ~ dnorm(600,0.01)     
	Ntot[t] <- min(max(1,(NadSurv[t] + N1[t])), capacity[t])		# enforces a carrying capacity of 600 and avoids extinction
	Ntotrd[t] <- round(Ntot[t])			
      y[t] ~ dpois(Ntot[t])
      }

}
",fill = TRUE)
sink()


# Bundle data

jags.data <- list(nyears = nyears, y = y, NAO=NAO,rain=rain)


# Initial values
inits <- function(){list(
nbroods	= runif(1, 1,2),
mphij = runif(1, 0.05,0.85),
mphia = runif(1, 0.1,0.85),
mfec = runif(1, 1.1,3.7)}


# Parameters monitored
parameters <- c("Ntot","phij","phia","f","nbrood","beta.phi","beta.f","beta.brood")


# MCMC settings
ni <- 150000
nt <- 2
nb <- 50000
nc <- 3

 
# Call JAGS from R
ipm.model <- jags(jags.data, inits, parameters, "AQWA.IPM.simple.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.cores=nc, parallel=T)
print(ipm.model, dig=3)




