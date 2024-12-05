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

# updated by Steffen on 4 Dec 2024 to include resighting data for survival
# best survival model includes age effect for survival and sex effect for detection probability

#The projection scenarios are the following:

#1) Mowing (habitat unconstrained) + 5 release years
#2) Mowing (habitat unconstrained) + 10 release years

#3) No mowing (habitat unconstrained) + 5 release years
#4) No mowing (habitat unconstrained) + 10 release years

#5) No mowing + habitat constrained 180 + 5 release years
#6) No mowing + habitat constrained 180 + 10 release years

#7) No mowing + habitat constrained 360 + 5 release years
#8) No mowing + habitat constrained 360 + 10 release years

rm(list=ls())
library(popbio)
library(doParallel)
library(foreach)
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
try(setwd("C:/STEFFEN/OneDrive - Vogelwarte/AQWA"),silent=T)
try(setwd("C:/Users/sop/OneDrive - Vogelwarte/AQWA"),silent=T)

AW<-read_excel("data/AW_raw_data.xlsx", sheet="Lit.Pop")[,1:2] %>%
  filter(year<2024)
head(AW)
years<-AW$year



### load and prepare mark-resight data
## see script AQWA_survival.r for various alternative options
AW_CH<-read_excel("data/AW_raw_data.xlsx", sheet="returns") %>%
  rename(ID=`ring ID`) %>%
  mutate(sex=if_any(everything(), ~str_detect(tolower(.), "female"))) %>%
  mutate(sex=ifelse(is.na(sex),1,2)) %>%
  mutate(`2018`=ifelse(`release year`==2018,1,0)) %>%
  mutate(`2019`=ifelse(`release year`==2019,1,`2019`)) %>%
  select(ID, sex,`2018`,`2019`,`2020`,`2021`,`2022`)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(AW_CH[,3:7], 1, get.first)

# Create age matrices X indicating age classes
age.mat<-as.matrix(AW_CH[,3:7])
age.mat[,2]<-ifelse(age.mat[,1]==1,2,1)
age.mat[,3:5]<-2


#Number of birds released, per scenario and year
releases <- AW %>%
  mutate(count=ifelse(year==2018,49,0)) %>%
  mutate(count=ifelse(year==2019,50,count))

# Population counts (from years 2003 to 2024)
y <- AW$count		# ENTER DATA OR READ IN AS VECTOR

# compile data (because survival is only fraction of time series, and we cannot estimate temporal variability, we keep survival separate)
jags.data <- list(ncountyears = length(years),
                  y = y,
                  releases = releases$count,
                  ## add survival data
                  n.marked=dim(AW_CH)[1],
                  n.markocc=dim(AW_CH[,3:7])[2],
                  age=age.mat,
                  f=f,
                  sex=as.numeric(AW_CH$sex),
                  y.mark = as.matrix(AW_CH[,3:7]))


##############################################################################
#
# IPM WITH DATA FOR SURVIVAL AND INFORMATIVE PRIORSFECUNDITY and releases --------
# 
##############################################################################

sink("models/AQWA.IPM.lithuania.jags")
cat("
model {

#----------------------------------------
# 1. Define the priors for the parameters
#----------------------------------------

# Initial population sizes

N1[1] ~ dnorm(3, 0.0001)T(0,)           # 1-year old individuals
NadSurv[1] ~ dnorm(5, 0.0001)T(0,)      # Adults >= 2 years
Ntot[1] <- NadSurv[1] + N1[1]	

# DEMOGRAPHIC PARAMETERS INFORMED BY DATA FROM OTHER STUDIES
mfec1 ~ dnorm(3.2,1/(0.16^2))			   ### fecundity = number of fledglings raised per FIRST brood
mfec2 ~ dnorm(1.75,1/(0.15^2))			   ### fecundity = number of fledglings raised per SECOND brood, should be slightly lower than first brood
prop.males ~ dnorm(0.56, 1/(0.01^2))T(0,1)  ### proportion of population that is male and can therefore be counted in population - 56%


# SURVIVAL PRIORS FOR AGE AND SEX GROUPS

mphi[2] ~ dnorm(0.42,1/(0.03^2))T(0,1)			### survival of adult birds
mphi[1] ~ dnorm(0.32,1/(0.025^2))T(0,1)		### survival of first year birds
      for (i in 1:n.marked){
        for (t in f[i]:(n.markocc-1)){
          phi[i,t] <- mphi[age[i,t]]
          p[i,t] <- mean.p[sex[i]]
        } #t
        p[i,n.markocc] <- mean.p[sex[i]]
      } #i
      
      for(sx in 1:2) {
        mean.p[sx] ~ dunif(0.05, 0.75)                  # resighting probability differs between sexes
      }


#--------------------------------------------------
# 2. Random variation in annual survival and productivity
#--------------------------------------------------
# CHANGE THE SCALE OF DEMOGRAPHIC PARAMETERS TO FACILITATE INCORPORATION OF COVARIATES
l.mphij<-log(mphi[1]/(1-mphi[1]))	    # juvenile survival probability on logit scale;
l.mphia<-log(mphi[2]/(1-mphi[2]))			# adult survival probability on logit scale
log.mfec1 <- log(mfec1)            #first-brood productivity on log scale
log.mfec2 <- log(mfec2)           #second-brood productivity on log scale

# RANDOM ANNUAL EFFECTS ON SURVIVAL AND PRODUCTIVITY

tau.surv<-1/pow(sigma.surv,2)
sigma.surv~dunif(0,5)

tau.prod <- pow(sigma.prod, -2)
sigma.prod ~ dunif(0,2)

for (t in 1:(ncountyears)){
   eps.surv[t] ~ dnorm(0, tau.surv)						# random variation around annual survival
   eps.fec[t] ~ dnorm(0, tau.prod)            # random variation around productivity
   logit(phij[t]) <- l.mphij+ eps.surv[t]			# add random effect - could add environmental predictors here (if we had any)
   logit(phia[t]) <- l.mphia+ eps.surv[t]			# add random effect - could add environmental predictors here (if we had any) 
   log(fec1[t]) <- log.mfec1 + eps.fec[t]     # add random effect
   log(fec2[t]) <- log.mfec2 + eps.fec[t]     # To use the same eps for both
   db[t] ~ dnorm(0.25,1/(0.05^2))T(0,1)       # proportion of females breeding twice
}


#-----------------------
# 3. Derived parameters
#-----------------------

# Population growth rate
for (t in 1:(ncountyears-1)){
   lambda[t] <- Ntot[t+1] / max(1,Ntot[t]) # prevent invalid parent error when Ntot=0
   loglambda[t]<-log(lambda[t])## for calculating geometric mean of overall population growth rate
}
    #### OVERALL POPULATION GROWTH RATE  #########
    mean.lambda<-exp((1/(ncountyears-1))*sum(loglambda[1:(ncountyears-1)]))   # Geometric mean


#--------------------------------------------
# 4. Likelihood for population count data
#--------------------------------------------
   # 4.1 System process
   for (t in 2:ncountyears){
        
    Nfem.breed2[t-1] ~ dbin(db[t-1],round((Ntot[t-1])*(1-prop.males)))    
    chicks[t-1] <- (Ntot[t-1])*(1-prop.males)* fec1[t-1]  + Nfem.breed2[t-1]*fec2[t-1]  # total fecundity
    chicksrd[t-1] <- round(chicks[t-1]) + releases[t]
		N1[t] ~ dbin(phij[t-1],chicksrd[t-1]) 
		NadSurv[t] ~ dbin(phia[t-1],round((Ntot[t-1])))
		Ntot[t] <- NadSurv[t] + N1[t]								## total population	
		
	}

   # 4.2 Observation process
   for (t in 1:(ncountyears-4)){    ### curtail likelihood so that last few years are only projection
	
	      Ntotobs[t] <- max(1,(Ntot[t]))*prop.males								## only males are counted, which is 56% of population
	      Ntotrd[t] <- round(Ntotobs[t])
        y[t] ~ dpois(Ntotrd[t])
        
   }
   
   # 4.3 Survival estimation
   for (i in 1:n.marked){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.markocc){
      # State process
      z[i,t] ~ dbern(z[i,t-1] * phi[i,t-1])
      # Observation process
      y.marked[i,t] ~ dbern(z[i,t] * p[i,t])
    } #t
  } #i
}
",fill = TRUE)
sink()



# Initial values
inits <- function(){list(
  z= zInit(as.matrix(AW_CH[,3:7])),
  mean.p = runif(2, 0.2,0.7),
  mphi= c(runif(1, 0.25,0.35),runif(1, 0.45,0.55)),
  mfec1 = runif(1, 2.5,4.0),
  mfec2 = runif(1, 1,2.5))}

# Parameters monitored
parameters <- c("Ntot","mphi","mfec1","mfec2","mean.lambda","prop.males", "mean.p")


# MCMC settings
ni <- 75000
nt <- 5
nb <- 25000
nc <- 4


# Call JAGS from R
ipm.model <- jags(jags.data,
                  inits,
                  parameters,
                  "models/AQWA.IPM.lithuania.jags",
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.cores=nc, parallel=T)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure 1: POPULATION TRAJECTORY IN PAST AND FUTURE BY SCENARIO ----------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AW$Ntot <- as.numeric(t(ipm.model$mean$Ntot))
AW$cim <- as.numeric(t(ipm.model$q2.5$Ntot))
AW$cip <- as.numeric(t(ipm.model$q97.5$Ntot))


#The plot

ggplot(data = AW) +
  geom_point(aes(x = year, y = count/ipm.model$mean$prop.males), col = "blue") +
  geom_line(aes(x = year, y = Ntot), col = "firebrick", linewidth=1) +
  geom_vline(aes(xintercept = 2018), linetype = 2) +
  geom_ribbon(aes(x = year, ymin = cim, ymax = cip), linetype = 2, fill = "firebrick", alpha=0.2) +
  theme_bw() + theme(legend.position = "none") +
  xlab("Year") + ylab("N AQWA males in Lithuania")

ggsave("output/Lithuania_projections.jpg", width=181,height=141, quality=100, units="mm")



############################################################################
#
# WRITE ALL OUTPUT INTO A TEXT FILE ----
# 
##############################################################################
out<-as.data.frame(ipm.model$summary)  
out$parameter<-row.names(ipm.model$summary)
names(out)[c(12,5,3,7)]<-c('parm','median','lcl','ucl')
print(ipm.model, dig=3)
write.table(out, "output/AQWA_LIT_model_output.csv", sep=",")




