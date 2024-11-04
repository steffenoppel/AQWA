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
setwd("C:\\Users\\sop\\OneDrive - Vogelwarte\\AQWA")


AW<-fread("AQWA_counts_GER.csv") %>%
  rename(Year=V1) %>%
  filter(!is.na(Year)) %>%
  filter(Year>2002) # omit the counts from the 1990s that suggest massive decreases
head(AW)
years<-AW$Year


# Population counts (from years 2003 to 2024)
y <- AW$Ntot		# ENTER DATA OR READ IN AS VECTOR
jags.data <- list(ncountyears = length(years-1), y = y[1:21], nintroyears = 30)




##############################################################################
#
# IPM WITH INFORMATIVE PRIORS FOR SURVIVAL AND FECUNDITY and releases --------
# 
##############################################################################

sink("AQWA.IPM.intro.jags")
cat("
model {


#----------------------------------------
# 1. Define the priors for the parameters
#----------------------------------------

# Initial population sizes
N1[1] ~ dnorm(60, 0.0001)T(0,)           # 1-year old individuals
NadSurv[1] ~ dnorm(100, 0.0001)T(0,)      # Adults >= 2 years

# DEMOGRAPHIC PARAMETERS INFORMED BY DATA FROM OTHER STUDIES
mphia ~ dunif(0.28,0.56)				### survival of adult females
mphij ~ dunif(0.20,0.44)				### survival of first year females
mfec1 ~ dunif(2.4,4.0)			   ### fecundity = number of fledglings raised per FIRST brood
mfec2 ~ dunif(1.0,2.5)			   ### fecundity = number of fledglings raised per SECOND brood, should be slightly lower than first brood



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
   lambda[t] <- Ntot[t+1] / min(1,Ntot[t]) # prevent invalid parent error when Ntot=0
   }


#--------------------------------------------
# 4. Likelihood for population count data
#--------------------------------------------
   # 4.1 System process
   for (t in 2:ncountyears){
        db[t-1] ~ dunif(0,0.5)            ### proportion of double-brooding females who would then have F1+F2 fecundity												# random variation around fecundity
        Nfemdb[t-1]  ~ dbin(db[t-1],round((Ntot[t-1])*0.44))  ## random draw of double brooding females
      	chicks[t-1] <- (Ntot[t-1])*0.44* mfec1 + Nfemdb[t-1]*mfec2			# total fecundity
      	chicksrd[t-1] <- round(chicks[t-1])
		N1[t] ~ dbin(phij[t-1],chicksrd[t-1]) 
		NadSurv[t] ~ dbin(phia[t-1],round((Ntot[t-1])))
	}

   # 4.2 Observation process
   for (t in 1:ncountyears){
   	    Ntot[t] <- NadSurv[t] + N1[t]								## total population		
	      Ntotobs[t] <- max(1,(Ntot[t]))*0.56								## only males are counted, which is 56% of population		
	      Ntotrd[t] <- round(Ntotobs[t])			
        y[t] ~ dpois(Ntotrd[t])
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

            CAPT.ADD[t] ~ dunif(45,50)   ## randomly draw a number for each year
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
parameters <- c("Ntot","phij","phia","mfec1","mfec2","db","lambda")


# MCMC settings
ni <- 15000
nt <- 1
nb <- 10000
nc <- 4


# Call JAGS from R
#ipm.model <- jags(jags.data, inits, parameters, "A:\\RSPB\\AquaticWarbler\\Analysis\\AQWA.IPM.imm5.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.cores=nc, parallel=T)
ipm.model <- jags(jags.data,
                  inits,
                  parameters,
                  "C:\\Users\\sop\\OneDrive - Vogelwarte\\AQWA\\AQWA.IPM.intro.jags",
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.cores=nc, parallel=T)




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








##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######
#################### SIMULATION OF POPULATION TRAJECTORY ACROSS RANGE OF PARAMETERS ########################################
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######

### SPECIFY RANGE OF PARAMETERS FOR DEMOGRAPHIC MODEL ###

pop.sizes<-seq(40,50,10)			### population size in individuals - ENTER THE RANGE YOU MAY RELEASE (assumed to be 50:50 males and females)
intro_years<-c(5,10,20)       ### enter the number of years in which the above number will be introduced EVERY YEAR
Sa<-seq(0.28,0.56,0.05)				### survival of adult females
Sj<-seq(0.20,0.44,0.05)				### survival of first year females
F1<-seq(2.4,4.0,0.2)			   ### fecundity = number of fledglings raised per FIRST brood
F2<-seq(1.0,2.5,0.17)			   ### fecundity = number of fledglings raised per SECOND brood, should be slightly lower than first brood
DB<-seq(0,0.5,0.1)            ### proportion of double-brooding females who would then have F1+F2 fecundity



### SPECIFY PARAMETERS FOR POPULATION VIABILITY ANALYSIS ###

nyears <- 30                    ### number of years over which simulations are run 
nreps <- 500                    ### number of stochastic simulations
K <- seq(54,643,100)            ### Carrying capacity of GER/Poland habitat for Aquatic Warbler

### SPECIFY PARAMETERS FOR STOCHASTIC COMPONENTS OF PVA ###

Catastrophe_prob <- 0.001       ### probability of a major disease, fire or other catastrophe occurring - 0.1% chance of major catastrophe
Catastrophe_severity <- 0.9     ### magnitude of catastrophe, in terms of size of population that will survive (90%)
breedfail_prob<-0.0001          ### chance that breeding season is 50% lower due to crazy weather events
extinctcutoffs<-c(1,2,5,10)     ### below which number the population is considered functionally extinct








##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######
#################### SETTING UP THE POPULATION MATRIX ########################################
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######

## SIMPLE RUN TO TEST FOR LATER INCORPORATION IN LOOP
## 56% of adult population is male

bird.matrix<-expression(
  0,  (F1*0.44)+(F2*0.44*DB),   ## proportion of males is 0.56
  Sj, Sa)

bird.vr<-list(F1=mean(F1),F2=mean(F2),Sa=mean(Sa), Sj=mean(Sj),DB=mean(DB))
A<-matrix(sapply(bird.matrix, eval,bird.vr , NULL), nrow=sqrt(length(bird.matrix)), byrow=TRUE)
projections<-pop.projection(A,n=c(50,100),iterations=50)
projections$lambda




##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######
#################### DEFINING FUNCTIONS FOR STOCHASTIC POPULATION VIABILITY ANALYSIS ########################################
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######


####  Density-dependent Ricker model to simulate population growth

Ricker <- function(prev_abund,K){       # this is a function for computing next-year abundance -- includes env stochasticity
  prev_abund * exp(log(R_max))*(1-(prev_abund/K))
}



####  Stochastic population viability analysis function
R_max <- 1
PVAdemo <- function(nreps,nyears,Init_N,Kloop,Catastrophe_prob,Catastrophe_severity,F1loop,F2loop,Sjloop,Saloop,DBloop,breedfail_prob,intro_years){
  PopArray2 <- array(0,dim=c((nyears+1),nreps))
  
  ## start looping through replicates
  
  for(rep in 1:nreps){
    
    # set initial abundance
    PopArray2[1,rep] <- Init_N      ### initial abundance of released birds in year 1
    
    ### loop through years
    for(y in 2:(nyears+1)){
      
      
      ### CREATE LESLIE MATRIX WITH SUBSET OF VITAL RATES
      ### FECUNDITY CAN SUFFER STOCHASTIC BAD SEASON
      
      F1YEAR<-ifelse(rbinom(1,1,breedfail_prob)==1,F1loop*0.5,F1loop) 			#### random draw to halve breeding success season due to bad weather
      F2YEAR<-ifelse(rbinom(1,1,breedfail_prob)==1,F2loop*0.5,F2loop) 			#### random draw to halve breeding success season due to bad weather
      bird.vr<-list(F1=F1YEAR, F2=F1YEAR, DB=DBloop, Sa=Saloop,Sj=Sjloop)
      A<-matrix(sapply(bird.matrix, eval,bird.vr , NULL), nrow=sqrt(length(bird.matrix)), byrow=TRUE)
      pop.size<-c(0,(Init_N/2))  ## number of introduced birds is equal sex ratio, and all adults
      projections<-pop.projection(A,n=pop.size,iterations=50)
      
      
      ### STOCHASTIC RICKER POPULATION MODEL
      
      R_max <- projections$lambda       # Maximum rate of growth (max lambda)
      
      
      ### stochasticity and density dependence
      nextyear <- max(0,trunc(Ricker(prev_abund=PopArray2[y-1,rep],K=Kloop)))    ### calculates abundance based on Ricker model, rounded to integer and set to a min of 0
      
      ### random catastrophic events
      if(runif(1)<Catastrophe_prob) nextyear <- nextyear*Catastrophe_severity
      
      ### continuing releases
      if(y<intro_years) nextyear <- nextyear+Init_N
      
      ### 
      PopArray2[y,rep] <- nextyear
      
    }
  }
  
  return(PopArray2)
}



#### CALCULATE PROPORTION OF SIMULATIONS WHERE SPECIES GOES EXTINCT

Extinction_bysim <- function(simdata, threshold){
  sum(apply(Default,2,function(t)  min(t)<threshold))/ncol(simdata)				### define extinction as <threshold birds
}




##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######
#################### LOOP OVER ALL COMBINATIONS OF DEMOGRAPHY AND RUN STOCHASTIC PVA ########################################
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######


### COMPREHENSIVE TABLE OF ALL COMBINATIONS OF DEMOGRAPHIC PARAMETERS ###

simul_in<-expand.grid(pop.sizes, Sa, Sj,F1,intro_years,K,DB) %>%
  mutate(F2=F2[match(Var4,F1)]) ### vary F1 and F2 simultaneously so that F2 is always lower than F1
dim(simul_in)
names(simul_in)<-c('pop.size','Sa','Sj','F1','intro_years','K','DB','F2')
SIM_OUT<-data.frame()


### setup parallel backend to use 8 processors

cl<-makeCluster(8)
registerDoParallel(cl, cores=8)



### START PARALLEL LOOP AND RETURN LAMBDA AND PROPORTION OF SIMULATIONS WHERE SPECIES GOES EXTINCT

SIM_OUT<-foreach(s=c(1:dim(simul_in)[1]), .packages='popbio',.combine=rbind) %dopar% {
  
  Init_N <- simul_in[s,1]
  Default <- PVAdemo(nreps,nyears,Init_N,
                     Kloop=simul_in[s,6],
                     Catastrophe_prob,Catastrophe_severity,
                     F1loop=simul_in[s,4],
                     F2loop=simul_in[s,8],
                     Sjloop=simul_in[s,3],
                     Saloop=simul_in[s,2],
                     DBloop=simul_in[s,7],
                     breedfail_prob,
                     intro_years=simul_in[s,5])
  
  
  ### CALCULATING MEAN POPULATION GROWTH RATE
  out<-simul_in[s,]
  bird.vr<-list(F1=simul_in[s,4], F2=simul_in[s,8], Sa=simul_in[s,2],Sj=simul_in[s,3], DB=simul_in[s,7])
  A<-matrix(sapply(bird.matrix, eval,bird.vr , NULL), nrow=sqrt(length(bird.matrix)), byrow=TRUE)
  pop.size<-c(0,(Init_N/2))
  projections<-pop.projection(A,n=pop.size,iterations=50)
  out$lambda<-projections$lambda
  
  
  ### CALCULATING EXTINCTION PROBABILITY
  OUT<-data.frame()
  for (t in extinctcutoffs){
    out$ExtThresh<-t
    out$OUTCOME<-Extinction_bysim(Default,t)
    OUT<-rbind(OUT,out)
  }
  
  
  return(OUT)
}
stopCluster(cl)

write.table(SIM_OUT,"AQWA_PVA_simulation_output_v1.csv", sep=",", row.names=FALSE)







##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######
#################### SUMMARISE OUTPUT AND PLOT EXTINCTION PROBABILITY WITH RESPECT TO HABITAT  ##############################
##########~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~######
#SIM_OUT<-fread("AQWA_PVA_simulation_output_v1.csv")
#names(SIM_OUT)<-c('pop.size','Sa','Sj','F','intro_years','lambda','ExtThresh','OUTCOME')
max(SIM_OUT$lambda)
dim(SIM_OUT)
head(SIM_OUT)
hist(SIM_OUT$lambda)
SIM_OUT %>% filter(is.na(OUTCOME))
SIM_OUT %>% filter(OUTCOME>0)
SIM_OUT %>% filter(lambda>1)

#### PLOT SIMULATION OUTPUT ####
SIM_OUT %>% 
  #filter(pop.size==50) %>%
  group_by(pop.size,intro_years, K) %>%
  summarise(prop.ext=mean(OUTCOME), lcl=quantile(OUTCOME,0.25), ucl=quantile(OUTCOME,0.75)) %>%
  #mutate(extinct=ifelse(ExtThresh==1,"<1 bird",paste("<",ExtThresh," birds"))) %>%
  
  ggplot() +
  geom_line(aes(x=K, y=prop.ext, colour=as.factor(pop.size))) +
  facet_wrap(~intro_years) +
  geom_ribbon(aes(x=K, ymin=lcl,ymax=ucl, fill=as.factor(pop.size)), alpha=0.3) +
  xlab("Carrying capacity (extent of habitat)") +
  ylab("Prob. of Aquatic Warbler extinct within 30 years") +
  labs(colour=expression(paste("Number introduced \n per year")),
       fill=expression(paste("Number introduced \n per year")),
       title="Aquatic Warbler re-introduction years") +
  theme(panel.background=element_rect(fill="white", colour="black"),
        legend.background = element_rect(),
        legend.title = element_text(size=16),
        legend.text=element_text(size=14),
        legend.position="inside",
        legend.position.inside=c(0.9,0.87),
        axis.text=element_text(size=16, color="black"), 
        axis.title=element_text(size=18), 
        strip.text.x=element_text(size=16, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())

ggsave("AquaticWarbler_PVA_ExtProb_by_habitat.jpg", width=11, height=8)










