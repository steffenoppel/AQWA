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

# updated by Steffen on 7 Dec 2024 to include more scenarios and updated priors for fecundity

# updated by Steffen on 4 Feb 2025 to adjust projection horizon to 25 years and releases to 5,10,15 years

# updated by Steffen on 11 Feb 2025 to remove temporal random effects and fix adult survival to informed prior (few data to estimate it)


## better priors needed for survival (from other migratory Acrocephalus warblers)
# Apparent survival for males (ϕ=0.26±0.06, CI 0.18−0.41) was higher than that for females (ϕ=0.18±0.07, CI 0.08–0.36) from https://www.tandfonline.com/doi/full/10.1080/00063657.2018.1448366#d1e445
# Annual adult survival ranged from 0.69 ± 0.05 to 0.88 ± 0.03 from https://bioone.org/journals/ardea/volume-103/issue-2/arde.v103i2.a5/Climatic-Influences-on-Survival-of-Migratory-African-Reed-Warblers-Acrocephalus/10.5253/arde.v103i2.a5.full
# survival rates were 60.3 ± 6.0 (se) for males and 54.9 ± 10.3 for females. Survival rates at Site B were 32.9 ± 16.0 for males and 52.0 ± 22.4 for females from https://www.tandfonline.com/doi/abs/10.1080/03078698.2006.9674347
# changed on 10 Feb 2025 to specify juvenile survival as proportion of adult survival

## priors from AQWA studies in Poland
# survival rates were 0.42 for females (0.290-0.563).(Dyrcz from AQWA Conservation handbook)
# survival rates were 0.36±0.04 (standard error) for M and F in SE Poland (Dyrcz from AQWA Conservation handbook)
# and the detection probability was 0.57±0.09 for males and 0.26±0.10


#The projection scenarios are the following:
# RUN EACH SCENARIO WITH AND WITHOUT IMPROVED SURVIVAL

#1) Mowing (habitat unconstrained) + 5 release years
#2) Mowing (habitat unconstrained) + 10 release years
#3) Mowing (habitat unconstrained) + 15 release years

#4) No mowing (habitat unconstrained) + 5 release years
#5) No mowing (habitat unconstrained) + 10 release years
#6) No mowing (habitat unconstrained) + 15 release years

#7) No mowing + habitat constrained 200ha + 5 release years
#8) No mowing + habitat constrained 200ha + 10 release years
#9) No mowing + habitat constrained 200ha + 15 release years

#10) No mowing + habitat constrained 400ha + 5 release years
#11) No mowing + habitat constrained 400ha + 10 release years
#12) No mowing + habitat constrained 400ha + 15 release years

#13) No mowing + habitat constrained 1200ha + 5 release years
#14) No mowing + habitat constrained 1200ha + 10 release years
#15) No mowing + habitat constrained 1200ha + 15 release years

#16) No mowing + habitat constrained 2400ha + 5 release years
#17) No mowing + habitat constrained 2400ha + 10 release years
#18) No mowing + habitat constrained 2400ha + 15 release years

#19) No mowing (habitat unconstrained) no releases but starting population at N1 (retrospective hypothetical scenario)
#20) Mowing (habitat unconstrained) no releases but starting population at N1 (retrospective hypothetical scenario with improved survival)


rm(list=ls())
library(popbio)
# library(doParallel)
# library(foreach)
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

AW<-fread("data/AQWA_counts_GER.csv") %>%
  rename(Year=V1) %>%
  filter(!is.na(Year)) %>%
  filter(Year>2002) # omit the counts from the 1990s that suggest massive decreases
head(AW)
years<-AW$Year



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

#Number of future scenarios modelled 
nscenarios <- 20

#Number of years for the projections into the future
nprojyears <- 25

#Number of birds released, per scenario and year
releases <- matrix(nrow = nscenarios, ncol = nprojyears, data = 50)
releases[c(1,4,7,10,13,16), 6:25] <- 0
releases[c(2,5,8,11,14,17), 11:25] <- 0
releases[c(3,6,9,12,15,18), 16:25] <- 0
releases[c(19,20), ] <- 0
releases <- rbind(releases,releases)  ## double to include with and without survival improvement

#A variable to model density-dependence in the number of breeding males that can reproduce.
# Remember, in the first scenarios, there is no density dependence.
# Then, we have some with 180 males allowed (120 ha), finally, one with 360 males allowed (240 ha).
# But we are not really interested in the number of males, but in the number of females.
# So far I will calculate the number of females inside the model, so that I can consider this extra added uncertainty around the sex ratio.
# changed to fixed number of females based on Susanne's email (16 Dec 2024)

# maxrepm <- numeric(length = nscenarios)
# maxrepm[1:6] <- 10000 #Dummy limit
# maxrepm[c(7:9, 13:15)] <- 180
# maxrepm[10:12] <- 360

maxrepf <- numeric(length = nscenarios)
maxrepf[c(1:6,19,20)] <- 10000 #Dummy limit
maxrepf[c(7:9)] <- 24
maxrepf[10:12] <- 48
maxrepf[13:15] <- 144
maxrepf[16:18] <- 288
maxrepf <- c(maxrepf,maxrepf)  ## double to include with and without survival improvement

## IMPROVEMENT IN SURVIVAL
impsurv <- numeric(length = nscenarios*2)
impsurv[1:nscenarios] <- 1 # no change in survival for first 19 scenarios
impsurv[(nscenarios+1):(nscenarios*2)] <- 1.05 # 5% improvement in survival for second 19 scenarios

# Population counts (from years 2003 to 2024)
y <- AW$Nmales		# ENTER DATA OR READ IN AS VECTOR

# compile data (because survival is only fraction of time series, and we cannot estimate temporal variability, we keep survival separate)
jags.data <- list(ncountyears = length(years-1),
                  y = y[1:21],
                  nscenarios = nscenarios*2,
                  nprojyears = nprojyears,
                  releases = releases,
                  maxrepf = maxrepf,
                  impsurv = impsurv,
                  ## add survival data
                  n.marked=dim(AW_CH)[1],
                  n.markocc=dim(AW_CH[,3:7])[2],
                  age=age.mat,
                  f=f,
                  sex=as.numeric(AW_CH$sex),
                  y.marked = as.matrix(AW_CH[,3:7]))

#Caution, we're missing data from year 2024, that leads to an error in the model. Let's add it manually.
jags.data$y[22] <- 3


### TESTING REASONABLE DISTRIBUTION FOR INFORMATIVE PRIORS ###

## adult female survival (for projection)
# 0.42 for females (0.290-0.563).(Dyrcz from AQWA Conservation handbook, p. 67)
# survival rates were 0.36±0.04 .(Foucher from AQWA Conservation handbook, p.68)
pa <- runif(1e6, 0.29, 0.56) #Mean is 0.42
pa <- rnorm(1e6, 0.385, 0.04)  ## standard error based on n=79 individuas that were resighted
#pa <- rbeta(1e6, 38, 65) ## attempt to use beta distribution
hist(pa)
quantile(pa)

## detection probability differs by sex
# detection probability was 0.57±0.09 for males and 0.26±0.10 for females.
pm <- rnorm(1e6, 0.57, 0.04)  ## males
hist(pm)
quantile(pm)

pf <- rnorm(1e6, 0.26, 0.04)  ## males
hist(pf)
quantile(pf)

## juvenile survival
mean(runif(1e6, 0.2, 0.44)) #Mean is 0.32
pj <- rnorm(1e6, 0.32, 0.025)
quantile(pj) #Nice

## productivity first clutch
mean(runif(1e6, 2.3, 4)) #Mean is 3.2
f1 <- rnorm(1e6, 3.15, .16)
quantile(f1)

## productivity second clutch
mean(runif(1e6, 1.8, 3.7)) #Mean is 2.75
f2 <- rnorm(1e6, 2.75, .2)
quantile(f2)

## sex ratio (proportion males)
pm <- rnorm(1e6, 0.56, 0.01)
hist(pm)

## probortion of double brooding females
db <- rnorm(1e6, 0.25, 0.07)
quantile(db)
hist(db)



hist(rbeta(1e6, 50, 5)) ## proportion of adult survival that we could use for juvenile survival



##############################################################################
#
# IPM WITH DATA FOR SURVIVAL AND INFORMATIVE PRIORSFECUNDITY and releases --------
# 
##############################################################################

sink("models/AQWA.IPM.surv.Scenarios.v10.jags")
cat("
model {

#----------------------------------------
# 1. Define the priors for the parameters
#----------------------------------------

# Initial population sizes

N1[1,1] ~ dnorm(60, 0.0001)T(0,)           # 1-year old individuals
NadSurv[1,1] ~ dnorm(100, 0.0001)T(0,)      # Adults >= 2 years

# DEMOGRAPHIC PARAMETERS INFORMED BY DATA FROM OTHER STUDIES
mfec1 ~ dnorm(3.15,1/(0.16^2))			   ### fecundity = number of fledglings raised per FIRST brood
mfec2 ~ dnorm(2.75,1/(0.2^2))			   ### fecundity = number of fledglings raised per SECOND brood, should be slightly lower than first brood
prop.males ~ dnorm(0.56, 1/(0.01^2))T(0,1)  ### proportion of population that is male and can therefore be counted in population - 56%


# SURVIVAL PRIORS FOR AGE AND SEX GROUPS

phi.ad ~ dnorm(0.385,0.04)T(0.29,)   ## adult femae survival for population projections based on info from Poland
mphi[2] ~ dbeta(1.2,1.2)			### survival of adult birds for estimation - NOT USED IN POPULATION PROCESS BECAUSE BASED ON TOO FEW DATA
mphi[1] ~	dbeta(1.2,1.2)	### survival of first year birds
      for (i in 1:n.marked){
        for (t in f[i]:(n.markocc-1)){
          phi[i,t] <- mphi[age[i,t]]
          p[i,t] <- mean.p[sex[i]]
        } #t
        p[i,n.markocc] <- mean.p[sex[i]]
      } #i
      
# PRIORS FOR RESIGHTING PROBABILITY BASED ON SE POLAND
mean.p[1] ~ dnorm(0.57, 0.04)                  # resighting probability for males
mean.p[2] ~ dnorm(0.26, 0.04)T(0,)            # resighting probability for females



#--------------------------------------------------
# 2. Random variation in annual survival and productivity
#--------------------------------------------------

for (t in 1:(ncountyears+nprojyears)){

   phij[t] <- mphi[1]			# fixed effect for every occasion - no temporal variability
   phia[t] <- mphi[2]			# fixed effect for every occasion - no temporal variability
   fec1[t] <-mfec1
}

#Fec 2, the fecundity of second broods, is only applied to future scenarios now

for(t in 1:nprojyears){
  fec2[t]<-mfec2
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
		NadSurv[1,t] ~ dbin(phi.ad,round((Ntot[1,t-1])))
	}

   # 4.2 Observation process
   for (t in 1:ncountyears){
   	    Ntot[1,t] <- NadSurv[1,t] + N1[1,t]								## total population		
	      Ntotobs[t] <- max(1,(Ntot[1,t]))*prop.males								## only males are counted, which is 56% of population
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
   
      
  #To make calculations easy, let's replicate the calculations for the past into the 8 scenarios for all objects that are involved as well into the future projections
  # ## TWO SCENARIOS SHOULD PROJECT HISTORIC TREND AND THUS NEED TO START AT HISTORIC POP SIZE
  # COULD NOT WORK OUT HOW TO FIX THIS - cannot define node twice and also cannot leave gaps in matrix
  # cumbersome long manual specification
  
  for(ns in 2:18){
  
    Ntot[ns,1:ncountyears] <- Ntot[1,1:ncountyears]
    N1[ns,1:ncountyears] <- N1[1,1:ncountyears]
    NadSurv[ns,1:ncountyears] <- NadSurv[1,1:ncountyears]
    chicks[ns,1:(ncountyears-1)] <- chicks[1,1:(ncountyears-1)]
    chicksrd[ns,1:(ncountyears-1)] <- chicksrd[1,1:(ncountyears-1)]

  }
  for(t in 2:ncountyears){
    Ntot[19,t] <- Ntot[1,1]
    N1[19,t] <- N1[1,1]
    NadSurv[19,t] <- NadSurv[1,1]
    chicks[19,(t-1)] <- chicks[1,1]
    chicksrd[19,(t-1)] <- chicksrd[1,1]
    Ntot[20,t] <- Ntot[1,1]
    N1[20,t] <- N1[1,1]
    NadSurv[20,t] <- NadSurv[1,1]
    chicks[20,(t-1)] <- chicks[1,1]
    chicksrd[20,(t-1)] <- chicksrd[1,1]
  }
  
  for(ns in 21:(nscenarios-2)){
  
    Ntot[ns,1:ncountyears] <- Ntot[1,1:ncountyears]
    N1[ns,1:ncountyears] <- N1[1,1:ncountyears]
    NadSurv[ns,1:ncountyears] <- NadSurv[1,1:ncountyears]
    chicks[ns,1:(ncountyears-1)] <- chicks[1,1:(ncountyears-1)]
    chicksrd[ns,1:(ncountyears-1)] <- chicksrd[1,1:(ncountyears-1)]

  }
  
  for(t in 2:ncountyears){
    Ntot[39,t] <- Ntot[1,1]
    N1[39,t] <- N1[1,1]
    NadSurv[39,t] <- NadSurv[1,1]
    chicks[39,(t-1)] <- chicks[1,1]
    chicksrd[39,(t-1)] <- chicksrd[1,1]
    Ntot[40,t] <- Ntot[1,1]
    N1[40,t] <- N1[1,1]
    NadSurv[40,t] <- NadSurv[1,1]
    chicks[40,(t-1)] <- chicks[1,1]
    chicksrd[40,(t-1)] <- chicksrd[1,1]
  }
      
# -------------------------------------------------        
# 5. PREDICTION INTO THE FUTURE WITH ONGOING RELEASES
# -------------------------------------------------

#Scenarios to be represented

#1) Mowing (habitat unconstrained) + 5 release years
#2) Mowing (habitat unconstrained) + 10 release years
#3) Mowing (habitat unconstrained) + 15 release years

#4) No mowing (habitat unconstrained) + 5 release years
#5) No mowing (habitat unconstrained) + 10 release years
#6) No mowing (habitat unconstrained) + 15 release years

#7) No mowing + habitat constrained 200ha + 5 release years
#8) No mowing + habitat constrained 200ha + 10 release years
#9) No mowing + habitat constrained 200ha + 15 release years

#10) No mowing + habitat constrained 400ha + 5 release years
#11) No mowing + habitat constrained 400ha + 10 release years
#12) No mowing + habitat constrained 400ha + 15 release years

#13) No mowing + habitat constrained 1200ha + 5 release years
#14) No mowing + habitat constrained 1200ha + 10 release years
#15) No mowing + habitat constrained 1200ha + 15 release years

#16) No mowing + habitat constrained 2400ha + 5 release years
#17) No mowing + habitat constrained 2400ha + 10 release years
#18) No mowing + habitat constrained 2400ha + 15 release years

#19) No mowing (habitat unconstrained) no releases but starting population at N1 (retrospective hypothetical scenario)
#20) Mowing (habitat unconstrained) no releases but starting population at N1 (retrospective hypothetical scenario)


#Hence, there should be variables storing all modifications in parameters involved in these simulations. These are:

#DB and NfemDB should have two dimensions: Nscenario, and time. DB can be manipulated a priori, but not NfemDB

#This should do it for the proportion of double broods

#Scenarios 1 - 3 and 20, no double broods
for(ns in c(1:3,20:23,40)){
  for(t in 1:nprojyears){
  
  db[ns,t] <- 0

  }
}

#The rest, double broods

for(ns in c(4:19,24:39)){
  for(t in 1:nprojyears){
  
  db[ns,t] ~ dnorm(0.25,1/(0.07^2))T(0,1)

  }
}

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
		        N1[ns,t] ~ dbin(min(1,phij[t-1]*impsurv[ns]),chicksrd[ns,t-1]) 
		        NadSurv[ns,t] ~ dbin(min(1,phi.ad*impsurv[ns]),round((Ntot[ns,t-1])))
		   	    Ntot[ns,t] <- NadSurv[ns,t] + N1[ns,t]								## total population	
  	} #Time loop
  	
  	
  	## FINAL POPULATION GROWTH RATE (steady state after reintroduction ends)
  	# Population growth rate
    for (t in (ncountyears+10):(ncountyears+nprojyears-1)){
      lambda.f[ns,t] <- Ntot[ns,t+1] / max(1,Ntot[ns,t]) # prevent invalid parent error when Ntot=0
      loglambda.f[ns,t]<-log(lambda.f[ns,t])## for calculating geometric mean of overall population growth rate
    }
    #### OVERALL POPULATION GROWTH RATE  #########
    proj.lambda[ns]<-exp((1/10)*sum(loglambda.f[ns,(ncountyears+10):(ncountyears+nprojyears-1)]))   # Geometric mean
  	
  } #Nscenarios loop
}
",fill = TRUE)
sink()



# Initial values
inits <- function(){list(
  z= zInit(as.matrix(AW_CH[,3:7])),
  prop.phij = rbeta(1,50,5),
  phi.ad = rnorm(1,0.38,0.02),
  mean.p = c(runif(1, 0.4,0.6),runif(1, 0.2,0.3)),
  mphi= c(runif(1, 0.25,0.35),runif(1, 0.45,0.55)),    ## using rbeta results in crazy output and model does not converge at all
  mfec1 = runif(1, 2.8,3.8),
  mfec2 = runif(1, 1.9,3.0))}

# Parameters monitored
parameters <- c("phi.ad", "mean.p","mphi","mfec1","mfec2","mean.lambda","prop.males","db", "Nfemdb",
                "proj.lambda","Ntot")


# MCMC settings
ni <- 75000
nt <- 5
nb <- 25000
nc <- 4


# Call JAGS from R
ipm.model <- jags(jags.data,
                  inits,
                  parameters,
                  "models/AQWA.IPM.surv.Scenarios.v10.jags",
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.cores=nc, parallel=T)


############################################################################
# CHECK CONVERGENCE AND WRITE ALL OUTPUT INTO A TEXT FILE ----
##############################################################################
out<-as.data.frame(ipm.model$summary)  
out$parameter<-row.names(ipm.model$summary)
names(out)[c(12,5,3,7)]<-c('parm','median','lcl','ucl')
print(ipm.model, dig=3)
out %>% arrange(desc(Rhat)) %>% select(parm, median, lcl, ucl, Rhat, n.eff)
out %>% filter(!startsWith(parm,"Ntot")) %>% filter(!startsWith(parm,"db")) %>%select(parm, median, lcl, ucl, Rhat, n.eff)
write.table(out, "output/AQWA_GER_model_nscenarios_output.v10.csv", sep=",")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure 1: POPULATION TRAJECTORY IN PAST AND FUTURE BY SCENARIO ----------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ntotdf <- data.frame(Ntot = numeric(length = length(ipm.model$mean$Ntot)),
                     cip = numeric(length = length(ipm.model$mean$Ntot)),
                     cim = numeric(length = length(ipm.model$mean$Ntot)),
                     Scenario = c(rep("Mowing + 5y Habitat Unconstrained", 47),
                                  rep("Mowing + 10y Habitat Unconstrained", 47),
                                  rep("Mowing + 15y Habitat Unconstrained", 47),
                                  rep("No mowing + 5y Habitat Unconstrained", 47),
                                  rep("No mowing + 10y Habitat Unconstrained", 47),
                                  rep("No mowing + 15y Habitat Unconstrained", 47),
                                  rep("No mowing + 5y Habitat Constrained 200 ha", 47),
                                  rep("No mowing + 10y Habitat Constrained 200 ha", 47),
                                  rep("No mowing + 15y Habitat Constrained 200 ha", 47),
                                  rep("No mowing + 5y Habitat Constrained 400 ha", 47),
                                  rep("No mowing + 10y Habitat Constrained 400 ha", 47),
                                  rep("No mowing + 15y Habitat Constrained 400 ha", 47),
                                  rep("No mowing + 5y Habitat Constrained 1200 ha", 47),
                                  rep("No mowing + 10y Habitat Constrained 1200 ha", 47),
                                  rep("No mowing + 15y Habitat Constrained 1200 ha", 47),
                                  rep("No mowing + 5y Habitat Constrained 2400 ha", 47),
                                  rep("No mowing + 10y Habitat Constrained 2400 ha", 47),
                                  rep("No mowing + 15y Habitat Constrained 2400 ha", 47),
                                  rep("Past trajectory without mowing", 47),
                                  rep("Past trajectory with mowing", 47)),
                     survival=rep(c("no survival improvement","5% survival improvement"), each=nscenarios*47),
                     year = rep(1:47, nscenarios*2))

ntotdf$Ntot <- as.numeric(t(ipm.model$mean$Ntot))
ntotdf$cim <- as.numeric(t(ipm.model$q2.5$Ntot))
ntotdf$cip <- as.numeric(t(ipm.model$q97.5$Ntot))

ntotdf$Scenario <- factor(ntotdf$Scenario, levels = c("Mowing + 5y Habitat Unconstrained", "Mowing + 10y Habitat Unconstrained", "Mowing + 15y Habitat Unconstrained",
                                                      "No mowing + 5y Habitat Unconstrained", "No mowing + 10y Habitat Unconstrained","No mowing + 15y Habitat Unconstrained",
                                                      "No mowing + 5y Habitat Constrained 200 ha", "No mowing + 10y Habitat Constrained 200 ha","No mowing + 15y Habitat Constrained 200 ha",
                                                      "No mowing + 5y Habitat Constrained 400 ha", "No mowing + 10y Habitat Constrained 400 ha","No mowing + 15y Habitat Constrained 400 ha",
                                                      "No mowing + 5y Habitat Constrained 1200 ha", "No mowing + 10y Habitat Constrained 1200 ha","No mowing + 15y Habitat Constrained 1200 ha",
                                                      "No mowing + 5y Habitat Constrained 2400 ha", "No mowing + 10y Habitat Constrained 2400 ha","No mowing + 15y Habitat Constrained 2400 ha",
                                                      "Past trajectory without mowing", "Past trajectory with mowing"))

#ntotdf$year <- factor(ntotdf$year)

#The plot for the past (with hypothetical changes of not mowing and not mowing and survival improvement)

poptrendspast <- 
  ntotdf %>%
  filter(Scenario %in% c("Past trajectory without mowing", "Past trajectory with mowing")) %>%
  filter(year>21) %>%
  mutate(year=year-21) %>%
  # bind_rows(ntotdf %>%
  #             filter(Scenario=="Mowing + 5y Habitat Unconstrained") %>%
  #             filter(year<22) %>%
  #             mutate(Scenario="Observed past")
  # ) %>%
  
  ggplot(aes(x = year+2002, y = Ntot, col = Scenario, group = Scenario)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = cim, ymax = cip), linetype = 2, fill = NA) +
  facet_grid(survival~Scenario, scales="free_y")  +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Year") + ylab("Aquatic Warbler population size")
poptrendspast
ggsave("output/Past_population_trajectories.jpg", width=241,height=191, quality=100, units="mm")



#The plot for the future

poptrends <- 
  ntotdf %>%
  filter(!(Scenario %in% c("Past trajectory without mowing", "Past trajectory with mowing"))) %>%

  ggplot(aes(x = year+2002, y = Ntot, col = Scenario, linetype=survival, group = survival)) +
  geom_line(linewidth = 1.1) +
  geom_vline(aes(xintercept = 2025), linetype = 2) +
  geom_ribbon(aes(ymin = cim, ymax = cip,linetype=survival), fill = NA) +
  facet_wrap(~Scenario, ncol = 3, scales="free_y")  +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Year") + ylab("Aquatic Warbler population size")
poptrends
#ggsave("output/Scenario_projections.jpg", width=241,height=191, quality=100, units="mm")


# Extinction probability (could plot this as well)
## revised to show prop samples where Ntot<3

w0 <- function(x){
  ext <- sum(x < 2) / length(x)
  return(ext)
}

ny <- jags.data$ncountyears
# 
# extprob <- matrix(nrow = nscenarios, ncol = nprojyears, data = 0)
# 
# #Quick plot, has to be improved in case we want to keep it
# 
# par(mfrow = c(4,2))
# for(i in 1:nscenarios){
#   for(t in 1:nprojyears){
#     extprob[i,t] <- w0(ipm.model$sims.list$Ntot[,i,nprojyears+jags.data$ncountyears])
#   }
#   plot(extprob[i,t], type = "l", lwd = 1.9, ylim = c(0,1), ylab = "P(ext)", xlab = "Year", main = paste0("Scenario ", i))
# }

## simpler plot

fut.ext<-as_tibble(rbind(ipm.model$samples[[1]],ipm.model$samples[[2]],ipm.model$samples[[3]],ipm.model$samples[[4]])) %>%
  dplyr::select(tidyselect::starts_with("Ntot")) %>%
  dplyr::select(tidyselect::ends_with(",47]"))

names(fut.ext)<- rep(c("Mowing + 5y Habitat Unconstrained", "Mowing + 10y Habitat Unconstrained", "Mowing + 15y Habitat Unconstrained",
                   "No mowing + 5y Habitat Unconstrained", "No mowing + 10y Habitat Unconstrained","No mowing + 15y Habitat Unconstrained",
                   "No mowing + 5y Habitat Constrained 200 ha", "No mowing + 10y Habitat Constrained 200 ha","No mowing + 15y Habitat Constrained 200 ha",
                   "No mowing + 5y Habitat Constrained 400 ha", "No mowing + 10y Habitat Constrained 400 ha","No mowing + 15y Habitat Constrained 400 ha",
                   "No mowing + 5y Habitat Constrained 1200 ha", "No mowing + 10y Habitat Constrained 1200 ha","No mowing + 15y Habitat Constrained 1200 ha",
                   "No mowing + 5y Habitat Constrained 2400 ha", "No mowing + 10y Habitat Constrained 2400 ha","No mowing + 15y Habitat Constrained 2400 ha",
                   "Past trajectory without mowing","Past trajectory with mowing"), 2)


fut.ext %>%
  pivot_longer(names(fut.ext), names_to = "Scenario", values_to = "N") %>%
  #gather(key="Scenario",value="N") %>%
  mutate(survival=rep(rep(c("no survival improvement","5% survival improvement"), each=(nscenarios)),40000)) %>%
  filter(!(str_detect(pattern="15y", string=Scenario))) %>%
  filter(!(Scenario %in% c("Past trajectory without mowing", "Past trajectory with mowing"))) %>%
  group_by(Scenario,survival) %>%
  summarize(ext.prob=w0(N)) %>%
  
  ### start the plot ###
  ggplot(aes(x = Scenario, y=ext.prob, fill = Scenario)) +                       # Draw overlaying histogram
  geom_bar(stat="identity",alpha = 0.6) +
  facet_wrap(~survival, ncol = 1)  +
  labs(x="Future scenarios", y="Probability of extinction within 20 years",
       fill="Scenario") +
  scale_fill_discrete() +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x=element_blank(), 
        axis.title=element_text(size=14),
        strip.text=element_text(size=14, color="black"),
        strip.background=element_rect(fill="white", colour="black"),
        legend.text=element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="bottom",
        legend.direction="horizontal",
        panel.grid.major = element_line(size=.1, color="grey94"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"))


ggsave("output/Extinction_probability.jpg", width=351,height=191, quality=100, units="mm")






#It looks like extinction probabilities 20 years into the future are not large. What about population decline probabilities? (That is, what is the probability that the population size at the end of the projections is lower than that at the end of the releases, in each scenario?)
# 
# pdecline <- numeric(length = nscenarios)
# relscenario <- rep(c(5,10),4)
# 
# for(i in 1:nscenarios){
#   
#   pdecline[i] = sum(ipm.model$sims.list$Ntot[,i,jags.data$ncountyears+nprojyears] < ipm.model$sims.list$Ntot[,i,(jags.data$ncountyears+relscenario[i])]) / dim(ipm.model$sims.list$Ntot)[1]
#   
# }
# 
# pdecline #Quite large

## plot histograms of future growth rates

as_tibble(rbind(ipm.model$samples[[1]],ipm.model$samples[[2]],ipm.model$samples[[3]],ipm.model$samples[[4]])) %>%
  dplyr::select(tidyselect::starts_with("proj.lambda")) %>%
  rename(nochange=`proj.lambda[1]`,
         no.mowing.unlim.hab=`proj.lambda[4]`,
         no.mowing.lim.hab.200=`proj.lambda[27]`,
         no.mowing.lim.hab.1200=`proj.lambda[33]`) %>%
  dplyr::select(-tidyselect::starts_with("proj.lambda")) %>%
  gather(key="Scenario",value="lambda") %>%
  filter(lambda>0.5) %>% ## remove bizarre growth rates caused by extinct populations?
  mutate(Scenario=if_else(Scenario=="no.mowing.lim.hab.200", "improved survival, no mowing, 200 ha habitat",
                          if_else(Scenario=="nochange", "current mowing pressure and survival",
                          if_else(Scenario=="no.mowing.unlim.hab","current survival and no mowing","improved survival, no mowing, 1200 ha habitat")))) %>%
  
  ### start the plot ###
  ggplot(aes(x = lambda, fill = Scenario)) +                       # Draw overlaying histogram
  geom_histogram(position = "identity", alpha = 0.2, bins = 80, aes(y = after_stat(density)), color="black") +
  geom_density(alpha=0.5) +
  #geom_vline(aes(xintercept = 1), colour="indianred3", size=1) +
  geom_segment(aes(x = 1, y = 0, xend = 1, yend = 10), colour="gray15", linetype = "dashed", linewidth=1)+
  
  labs(x="Future population growth rate", y="Probability density",
       fill="Scenario") +
  scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) +
  scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title = element_text(size=14),
        legend.position="inside",
        legend.position.inside=c(0.14,0.88),
        panel.grid.major = element_line(size=.1, color="grey94"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"))


ggsave("output/Future_pop_growth_rates.jpg", width=235,height=181, quality=100, units="mm")



#On the retrospective analysis (past trends), see how inference has modified the priors and the estimates.
## NOT UPDATED THIS CODE WITH SCENARIOS V4
# par(mfrow = c(2,2))
# 
# #Juvenile survival
# 
# plot(ipm.model$mean$phij, type = "l", col = "green", ylim = c(0.18, 0.45), lwd = 2, xlab = "Year", ylab = "Estimate", main = "Juvenile survival")
# lines(ipm.model$q2.5$phij, lty = "dashed", lwd = 1.5, col = "green")
# lines(ipm.model$q97.5$phij, lty = "dashed", lwd = 1.5, col = "green")
# abline(v = ny, lty = "dashed")
# 
# #Adult survival
# 
# plot(ipm.model$mean$phia, type = "l", col = "green", ylim = c(0.25,0.6), lwd = 2, xlab = "Year", ylab = "Estimate", main = "Adult survival")
# lines(ipm.model$q2.5$phia, lty = "dashed", lwd = 1.5, col = "green")
# lines(ipm.model$q97.5$phia, lty = "dashed", lwd = 1.5, col = "green")
# abline(v = ny, lty = "dashed")
# 
# #First-year productivity
# 
# plot(ipm.model$mean$fec1, type = "l", col = "red", ylim = c(1.5,5.6), lwd = 2, xlab = "Year", ylab = "Estimate", main = "First-year productivity")
# lines(ipm.model$q2.5$fec1, lty = "dashed", lwd = 1.5, col = "red")
# lines(ipm.model$q97.5$fec1, lty = "dashed", lwd = 1.5, col = "red")
# abline(v = ny, lty = "dashed")
# 
# #Second-year productivity
# 
# plot(ipm.model$mean$fec2, type = "l", col = "red", ylim = c(0.7,3.2), lwd = 2, xlab = "Year", ylab = "Estimate", main = "Second-year productivity")
# lines(ipm.model$q2.5$fec2, lty = "dashed", lwd = 1.5, col = "red")
# lines(ipm.model$q97.5$fec2, lty = "dashed", lwd = 1.5, col = "red")
# abline(v = ny, lty = "dashed")
# 
# #See how priors compare to final estimates
# 
# par(mfrow = c(2,2))
# 
# #Adult survival
# 
# adprior <- dnorm(seq(from = 0, to = 1, by = 0.01), mean = 0.42, sd = 0.03)
# plot(adprior, type = "l", xaxt = "n", ylab = "", xlab = "Estimate value", main = "Adult Survival \n Red = Prior, Blue = Estimate", ylim = c(0,15), col = "red", lwd = 2)
# axis(1, at = seq(from = 0, to = 100, by = 10), labels = seq(from = 0, to = 1, by = 0.1))
# ipm.model$mean$mphi[2]
# ipm.model$sd$mphi[2]
# adresult <- dnorm(seq(from = 0, to = 1, by = 0.01), mean = ipm.model$mean$mphi[2], sd = ipm.model$sd$mphi[2])
# lines(adresult, type = "l", lty = "dashed", col = "blue", lwd = 2)
# 
# #Juvenile survival
# 
# juvprior <- dnorm(seq(from = 0, to = 1, by = 0.01), mean = 0.32, sd = 0.025)
# plot(juvprior, type = "l", xaxt = "n", ylab = "", xlab = "Estimate value", main = "Juvenile Survival \n Red = Prior, Blue = Estimate", ylim = c(0,19), col = "red", lwd = 2)
# axis(1, at = seq(from = 0, to = 100, by = 10), labels = seq(from = 0, to = 1, by = 0.1))
# ipm.model$mean$mphi[1]
# ipm.model$sd$mphi[1]
# juvresult <- dnorm(seq(from = 0, to = 1, by = 0.01), mean = ipm.model$mean$mphi[1], sd = ipm.model$sd$mphi[1])
# lines(juvresult, type = "l", lty = "dashed", col = "blue", lwd = 2)
# 
# #First-brood productivity
# 
# p1prior <- dnorm(seq(from = 2, to = 4, by = 0.1), mean = 3.2, sd = 0.16)
# plot(p1prior, type = "l", xaxt = "n", ylab = "", xlab = "Estimate value", main = "First-brood productivity \n Red = Prior, Blue = Estimate", ylim = c(0,3), col = "red", lwd = 2)
# axis(1, at = seq(from = 1, to = 21, by = 1), labels = seq(from = 2, to = 4, by = 0.1))
# ipm.model$mean$mfec1
# ipm.model$sd$mfec1
# p1result <- dnorm(seq(from = 2, to = 4, by = 0.1), mean = ipm.model$mean$mfec1, sd = ipm.model$sd$mfec1)
# lines(p1result, type = "l", lty = "dashed", col = "blue", lwd = 2)
# 
# #Second-brood productivity
# 
# p2prior <- dnorm(seq(from = 0.5, to = 3.2, by = 0.1), mean = 1.75, sd = 0.15)
# plot(p2prior, type = "l", xaxt = "n", ylab = "", xlab = "Estimate value", main = "Second-brood productivity \n Red = Prior, Blue = Estimate", ylim = c(0,3), col = "red", lwd = 2)
# axis(1, at = seq(from = 1, to = 28, by = 1), labels = seq(from = 0.5, to = 3.2, by = 0.1))
# ipm.model$mean$mfec2
# ipm.model$sd$mfec2
# p2result <- dnorm(seq(from = 0.5, to = 3.2, by = 0.1), mean = ipm.model$mean$mfec2, sd = ipm.model$sd$mfec2)
# lines(p2result, type = "l", lty = "dashed", col = "blue", lwd = 2)
# 
# #Compared to models where we allowed second broods, here the change in the estimates is very very mild.




