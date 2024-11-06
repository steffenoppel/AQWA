###################################################################################################
########   AQUATIC WARBLER PVA TO EXPLORE REINTRODUCTION TO GERMANY ################
###################################################################################################
# written by Steffen Oppel on 28 October 2021 (steffen.oppel@rspb.org.uk)
# simulation taken from https://naes.unr.edu/shoemaker/teaching/NRES-470/PVA1_421.html
# requested by Susanne Arbeiter to examine needs for reinforcement

# updated on 3 Nov 2024 by adopting AQWA IPM code from 2016
# IPM written without data for population in Hungary
# adopted model structure for equally data-poor situation in Germany


# Grundsätzlich möchte ich mit dem paper die Fragen beantworten:
# Was ist demographisch mit der Population passiert? (wenn möglich, sonst kann ich darauf am ehesten verzichten)
# Wieviel Jahre müssen wir auswildern, bis die Population stabil ist (z.B. < 0.10 extinction probability after 50 years)?
# Wie können wir während der nächsten Jahre überprüfen, ob wir auf diesem Pfad sind?


### NEED TO DO: BUILD IN DENSITY/HABITAT DEPENDENCE
## ADD SOME ENVIRONMENTAL CONSTRAINT THAT CAN BE LIFTED IN FUTURE


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
setwd("C:/Users/jba/OneDrive - Vogelwarte/Projects/Aquatic Warbler Steffen/Github repo/AQWA")

AW<-fread("data/AQWA_counts_GER.csv") %>%
  rename(Year=V1) %>%
  filter(!is.na(Year)) %>%
  filter(Year>2002) # omit the counts from the 1990s that suggest massive decreases
head(AW)
years<-AW$Year


# Population counts (from years 2003 to 2024)
y <- AW$Nmales		# ENTER DATA OR READ IN AS VECTOR
jags.data <- list(ncountyears = length(years-1), y = y[1:21], nintroyears = 30)

#Caution, we're missing data from year 2024, that leads to an error in the model. Let's add it manually.
jags.data$y[22] <- 3


#Let's define new priors

mean(runif(1e6, 0.28, 0.56)) #Mean is 0.42
pa <- rnorm(1e6, 0.42, 0.03)
hist(pa)
quantile(pa)

#pa <- (rnorm(1e6, 0.42, 0.01))
#spa <- rnorm(1e6, 0.011, 0.01)
#pasim <- (rnorm(1e6, sample(pa, size = 1e6, replace = TRUE), sample(spa, size = 1e6, replace = TRUE))) #NAs: negative sigma (no probs, just ignore)
#hist(pasim)
#quantile(pasim, na.rm = T) #Looks more or less cool

mean(runif(1e6, 0.2, 0.44)) #Mean is 0.32
pj <- rnorm(1e6, 0.32, 0.025)
quantile(pj) #Nice


#pj <- rnorm(1e6, 0.32, 0.01)
#spj <- rnorm(1e6, 0.03, 0.02)
#pjsim <- (rnorm(1e6, sample(pj, size = 1e6, replace = TRUE), sample(spj, size = 1e6, replace = TRUE))) #NAs: negative sigma (no probs, just ignore)

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

sink("models/AQWA.IPM.Jaume.jags")
cat("
model {


#----------------------------------------
# 1. Define the priors for the parameters
#----------------------------------------

# Initial population sizes
N1[1] ~ dnorm(60, 0.0001)T(0,)           # 1-year old individuals
NadSurv[1] ~ dnorm(100, 0.0001)T(0,)      # Adults >= 2 years

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

# RANDOM ANNUAL EFFECT ON SURVIVAL
tau.surv<-1/pow(sigma.surv,2)
sigma.surv~dunif(0,5)

for (t in 1:(ncountyears+nintroyears)){
   eps.surv[t] ~ dnorm(0, tau.surv)						# random variation around annual survival
   logit(phij[t]) <- l.mphij+ eps.surv[t]			# add random effect - could add environmental predictors here (if we had any)
   logit(phia[t]) <- l.mphia+ eps.surv[t]			# add random effect - could add environmental predictors here (if we had any) 
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
        db[t-1] ~ dnorm(0.25,1/(0.05^2))T(0,1)            ### proportion of double-brooding females who would then have F1+F2 fecundity												# random variation around fecundity
        Nfemdb[t-1]  ~ dbin(db[t-1],round((Ntot[t-1])*(1-prop.males)))  ## random draw of double brooding females
      	chicks[t-1] <- (Ntot[t-1])*(1-prop.males)* mfec1 + Nfemdb[t-1]*mfec2			# total fecundity
      	chicksrd[t-1] <- round(chicks[t-1])
		N1[t] ~ dbin(phij[t-1],chicksrd[t-1]) 
		NadSurv[t] ~ dbin(phia[t-1],round((Ntot[t-1])))
	}

   # 4.2 Observation process
   for (t in 1:ncountyears){
   	    Ntot[t] <- NadSurv[t] + N1[t]								## total population		
	      Ntotobs[t] <- max(1,(Ntot[t]))*prop.males								## only males are counted, which is 56% of population
	      Ntotrd[t] <- round(Ntotobs[t])
        y[t] ~ dpois(Ntotrd[t])
        #  y[t] ~ dbin(prop.males,round(Ntot[t]))                 ## changed from a Poisson to a binomial draw where 56% of the population can be recorded as males - results in inconsistent node error
   }
      
      
# -------------------------------------------------        
# 5. PREDICTION INTO THE FUTURE WITH ONGOING RELEASES
# -------------------------------------------------
### POTENTIALLY INCLUDE SCENARIOS FOR N CAPTIVE BIRDS AND NUMBER OF YEARS

    #for (ncr in 1:scen.capt.release){
    CAPT.ADD[ncountyears] <-48   ## 48 were released in first year

      # SPECIFY IMPROVEMENT OF SURVIVAL
      #for (is in 1:scen.imp.surv){ 
    
        # improve.surv ~ dunif(1,1.10)
        # fut.survival[ncr,is] <-imp.surv[is]*mean(ann.surv.terrvis[1:nyears.terrvis])

        ## POPULATION PROCESS
        
        for (t in (ncountyears+1):(ncountyears+nintroyears)){
        
            # CAPTIVE RELEASE OF JUVENILE BIRDS

            CAPT.ADD[t] ~ dunif(44.51,50.49)   ## randomly draw a number for each year. Changed to allow equal probabilities to all numbers. Otherwise we could use a categorical, it would probably be more beautiful
            db[t-1] ~ dunif(0,0.5)            ### proportion of double-brooding females who would then have F1+F2 fecundity												# random variation around fecundity
            Nfemdb[t-1]  ~ dbin(db[t-1],round((Ntot[t-1])*0.44))  ## random draw of double brooding females
      	    chicks[t-1] <- (Ntot[t-1])*0.44* mfec1 + Nfemdb[t-1]*mfec2			# total fecundity
      	    chicksrd[t-1] <- round(chicks[t-1]) + round(CAPT.ADD[t-1])
		        N1[t] ~ dbin(phij[t-1],chicksrd[t-1]) 
		        NadSurv[t] ~ dbin(phia[t-1],round((Ntot[t-1])))
		   	    Ntot[t] <- NadSurv[t] + N1[t]								## total population	
	}


    #  } # end scenario of imp.surv

    #} # end scenario of n capt released
      
      
      
      
      
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
parameters <- c("Ntot","mphij","mphia","mfec1","mfec2","mean.lambda","prop.males", "phij", "phia")


# MCMC settings
ni <- 75000
nt <- 5
nb <- 25000
nc <- 4


# Call JAGS from R
#ipm.model <- jags(jags.data, inits, parameters, "A:\\RSPB\\AquaticWarbler\\Analysis\\AQWA.IPM.imm5.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.cores=nc, parallel=T)
ipm.model <- jags(jags.data,
                  inits,
                  parameters,
                  "models/AQWA.IPM.Jaume.jags",
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.cores=nc, parallel=T)



#First thing: does inference affect the priors?

ipm.model$mean$Ntot
ipm.model$q97.5$Ntot

#First result: confidence intervals in population projections are narrower

plot(ipm.model$mean$Ntot, type = "l", lwd = 2, col = "blue", ylim = c(0,400), ylab = "Counts", xlab = "Year", main = "Population counts")
lines(ipm.model$q2.5$Ntot, lwd = 1.5, lty = "dashed", col = "blue")
lines(ipm.model$q97.5$Ntot, lwd = 1.5, lty = "dashed", col = "blue")

#Extinction probability

w0 <- function(x){
  
  ext <- sum(x == 0) / length(x)
  return(ext)
  
}

ny <- jags.data$ncountyears

extprob <- apply(ipm.model$sims.list$Ntot, 2, w0)[(ny+1):(ny+jags.data$nintroyears)] 

#Probability of extinction is 0 so far. Here we would have to play with the number of years with releases to see how the results change as a function of it.

library(ggplot2)

#Juvenile survival

plot(ipm.model$mean$phij, type = "l", col = "green", ylim = c(0.18, 0.45), lwd = 2, xlab = "Year", ylab = "Estimate", main = "Juvenile survival")
lines(ipm.model$q2.5$phij, lty = "dashed", lwd = 1.5, col = "green")
lines(ipm.model$q97.5$phij, lty = "dashed", lwd = 1.5, col = "green")
abline(v = ny, lty = "dashed")

#Adult survival

plot(ipm.model$mean$phia, type = "l", col = "red", ylim = c(0.25,0.6), lwd = 2, xlab = "Year", ylab = "Estimate", main = "Adult survival")
lines(ipm.model$q2.5$phia, lty = "dashed", lwd = 1.5, col = "red")
lines(ipm.model$q97.5$phia, lty = "dashed", lwd = 1.5, col = "red")
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

par(mfrow = c(1,1))
#How does the male proportion looks like?

ipm.model$mean$prop.males
ipm.model$q2.5$prop.males
ipm.model$q97.5$prop.males

#Nice

############################################################################
#
# WRITE ALL OUTPUT INTO A TEXT FILE ----
# 
##############################################################################
out<-as.data.frame(ipm.model$summary)  
out$parameter<-row.names(ipm.model$summary)
names(out)[c(12,5,3,7)]<-c('parm','median','lcl','ucl')
print(ipm.model, dig=3)
write.table(out, "output/AQWA_GER_model_JAUME_output.csv", sep=",")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure 1: POPULATION TRAJECTORY IN PAST AND FUTURE ----------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## retrieve the past population estimates (2003-2023)
AQWA.pred<-out[(grep("Ntot\\[",out$parm)),c(12,5,3,7)] %>%
  mutate(Year=seq(min(AW$Year),(max(AW$Year)+jags.data$nintroyears),1)) %>%
  mutate(ucl=ifelse(ucl>500,500,ucl))


## CREATE A COLOUR PALETTE FOR THE NUMBER OF CHICKS RELEASED
colfunc <- colorRampPalette(c("cornflowerblue", "firebrick"))


ggplot()+
  geom_ribbon(data=AQWA.pred,aes(x=Year, ymin=lcl,ymax=ucl),alpha=0.2, fill="firebrick")+
  geom_line(data=AQWA.pred, aes(x=Year, y=median),linewidth=1, col="firebrick")+
  geom_point(data=AW,aes(x=Year, y=Ntot), size=2,col='darkblue')+
  
  ## format axis ticks
  scale_y_continuous(name="N Aquatic Warblers", limits=c(0,400),breaks=seq(0,400,50), labels=as.character(seq(0,400,50)))+
  scale_x_continuous(name="Year", breaks=seq(2000,2055,5), labels=as.character(seq(2000,2055,5)))+
  
  ## add vertical line for when introduction starts
  geom_vline(aes(xintercept=2023), linetype=2, col="forestgreen", linewidth=2) +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=16, color="black"),
        axis.text.x=element_text(size=12, color="black",angle=45, vjust = 1, hjust=1), 
        axis.title=element_text(size=16),
        legend.position = c(0.90, 0.82))

#ggsave("output/AQWA_projection.jpg", width=180,height=131, quality=100, units="mm")


