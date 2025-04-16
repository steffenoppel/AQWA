###################################################################################################
########   AQUATIC WARBLER PVA TO EXPLORE REINTRODUCTION TO GERMANY ################
###################################################################################################
# written by Steffen Oppel and Jaume Badia in February 2025
# simulation taken from https://naes.unr.edu/shoemaker/teaching/NRES-470/PVA1_421.html
# additional IPM requested by Susanne Arbeiter to examine effects on source population

# scenarios with no mowing and unlimited habitat only for Biebrza population
# each year 10-12 nests are removed and females can only produce one secondary clutch

rm(list=ls())
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
  select(ID, sex,`2018`,`2019`,`2020`,`2021`,`2022`,`release year`)

## check raw sum of birds that return
AW_CH %>% rowwise() %>%
  #filter(`release year`==2019) %>%
  mutate(ret=ifelse(`2018`==1,
                            max(`2019`,`2020`,`2021`,`2022`),max(`2020`,`2021`,`2022`))) %>%
  group_by(sex) %>%
  summarise(mean(ret))
  

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(AW_CH[,3:7], 1, get.first)

# Create age matrices X indicating age classes
age.mat<-as.matrix(AW_CH[,3:7])
age.mat[,2]<-ifelse(age.mat[,1]==1,2,1)
age.mat[,3:5]<-2

#Number of future scenarios modelled 
nscenarios <- 16

#Number of nests removed
n_nest_removed <- seq(0,15,1)

#Number of years for the projections into the future
nprojyears <- 25


# ## CHANGE IN SURVIVAL (HIGHER OR LOWER)
# impsurv <- numeric(length = nscenarios*3)
# impsurv[1:nscenarios] <- 1 # no change in survival for first 19 scenarios
# impsurv[(nscenarios+1):(nscenarios*2)] <- 1.05 # 5% improvement in survival for second 19 scenarios
# impsurv[(nscenarios*2+1):(nscenarios*3)] <- 0.95 # 5% reduction in survival for third 19 scenarios

# Population counts (from years 2003 to 2024)
y <- AW$Nmales+2408		# MADE UP NUMBER TO MATCH COUNT IN 2024

# compile data (because survival is only fraction of time series, and we cannot estimate temporal variability, we keep survival separate)
jags.data <- list(ncountyears = length(years-1),
                  y = y,
                  nscenarios = nscenarios,
                  nprojyears = nprojyears,
                  n_nest_removed = n_nest_removed,
                  #impsurv = impsurv,
                  ## add survival data
                  n.marked=dim(AW_CH)[1],
                  n.markocc=dim(AW_CH[,3:7])[2],
                  age=age.mat,
                  f=f,
                  sex=as.numeric(AW_CH$sex),
                  y.marked = as.matrix(AW_CH[,3:7]))

#Caution, we're missing data from year 2024, that leads to an error in the model. Let's add it manually.
jags.data$y[22] <- 2411




##############################################################################
#
# IPM WITH DATA FOR SURVIVAL AND ANNUAL REMOVAL OF NESTS --------
# 
##############################################################################

sink("models/AQWA.IPM.source.pop.jags")
cat("
model {

#----------------------------------------
# 1. Define the priors for the parameters
#----------------------------------------

# Initial population sizes

N1[1,1] ~ dnorm(600, 0.0001)T(0,)           # 1-year old individuals
NadSurv[1,1] ~ dnorm(1600, 0.0001)T(0,)      # Adults >= 2 years
Nfem.breed1[1,1] <- (Ntot[1,1]*(1-prop.males))
Nfem.breed2[1,1] ~ dbin(db[1,1],round(Nfem.breed1[1,1]))
chicks[1,1] <- Nfem.breed1[1,1] * mfec1 + Nfem.breed2[1,1]*mfec2			# total fecundity
chicksrd[1,1] <- round(chicks[1,1])



# DEMOGRAPHIC PARAMETERS INFORMED BY DATA FROM OTHER STUDIES
mfec1 ~ dnorm(3.15,1/(0.16^2))			   ### fecundity = number of fledglings raised per FIRST brood
mfec2 ~ dnorm(2.75,1/(0.2^2))			   ### fecundity = number of fledglings raised per SECOND brood, should be slightly lower than first brood
prop.males ~ dnorm(0.56, 1/(0.01^2))T(0,1)  ### proportion of population that is male and can therefore be counted in population - 56%


# SURVIVAL PRIORS FOR AGE AND SEX GROUPS

phi.ad ~ dnorm(0.385,1/(0.04^2))  ## adult femae survival for population projections based on info from Poland
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
mean.p[1] ~ dnorm(0.57, 1/(0.04^2))                  # resighting probability for males
mean.p[2] ~ dnorm(0.26, 1/(0.04^2))T(0,)            # resighting probability for females



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

# 4.1. double brooding for all years

for(ns in 1:nscenarios){
  for(t in 1:(ncountyears+nprojyears)){
  
  db[ns,t] ~ dnorm(0.25,1/(0.07^2))T(0,1)

  }
}

# 4.2 System process
for (t in 2:ncountyears){
   
    Nfem.breed1[1,t] <- (Ntot[1,t]*(1-prop.males))
    Nfem.breed2[1,t] ~ dbin(db[1,t],round(Nfem.breed1[1,t]))

    chicks[1,t] <- Nfem.breed1[1,t] * fec1[t] + Nfem.breed2[1,t]*fec2[t]			# total fecundity
    chicksrd[1,t] <- round(chicks[1,t])
        
		N1[1,t] ~ dbin(phij[t-1],chicksrd[1,t-1]) 
		NadSurv[1,t] ~ dbin(phi.ad,round((Ntot[1,t-1])))
}

# 4.3 Observation process
for (t in 1:ncountyears){
   	    Ntot[1,t] <- NadSurv[1,t] + N1[1,t]								## total population		
	      Ntotobs[t] <- max(1,(Ntot[1,t]))*prop.males								## only males are counted, which is 56% of population
	      Ntotrd[t] <- round(Ntotobs[t])
        y[t] ~ dpois(Ntotrd[t])
        
}
   

# 4.4 Survival estimation
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
   

      
# -------------------------------------------------        
# 5. PREDICTION INTO THE FUTURE WITH ONGOING CAPTURE FOR TRANSLOCATION
# -------------------------------------------------


### COPY THE PAST VALUES FOR ALL SCENARIOS

for(ns in 2:nscenarios){
  for(t in 1:ncountyears){
  
    Nfem.breed1[ns,t] <- Nfem.breed1[1,t]
    Nfem.breed2[ns,t] <- Nfem.breed2[1,t]
    chicks[ns,t] <- chicks[1,t]
    chicksrd[ns,t] <- chicksrd[1,t]
    N1[ns,t] <- N1[1,t]
		NadSurv[ns,t] <- NadSurv[1,t]
		Ntot[ns,t] <-Ntot[1,t]
  }
}





    
## POPULATION PROCESS

    for(ns in 1:nscenarios){
      for (t in (ncountyears+1):(ncountyears+nprojyears)){
        
            # Apply the removal of nests
            Nfem.breed1[ns,t] <- (Ntot[ns,t]*(1-prop.males)) - n_nest_removed[ns]    ### N breeders is reduced by number of nests that are removed
            Nfem.breed2[ns,t] ~ dbin(db[ns,t],round(Nfem.breed1[ns,t] + n_nest_removed[ns]))  ### the removed nests woud still allow the females to second brood

      	    chicks[ns,t] <- Nfem.breed1[ns,t] * fec1[t] + Nfem.breed2[ns,t]*mfec2			# total fecundity
      	    chicksrd[ns,t] <- round(chicks[ns,t])
		        N1[ns,t] ~ dbin(phij[t-1],max(chicksrd[ns,t-1],1)) 
		        NadSurv[ns,t] ~ dbin(phi.ad,round((Ntot[ns,t-1])))
		   	    Ntot[ns,t] <- NadSurv[ns,t] + N1[ns,t]								## total population	
  	} #Time loop
  	
  	
  	## FINAL POPULATION GROWTH RATE (steady state after reintroduction ends)
  	# Population growth rate
    for (t in (ncountyears):(ncountyears+nprojyears-1)){
      lambda.f[ns,t] <- Ntot[ns,t+1] / max(1,Ntot[ns,t]) # prevent invalid parent error when Ntot=0
      loglambda.f[ns,t]<-log(lambda.f[ns,t])## for calculating geometric mean of overall population growth rate
    }
    #### OVERALL POPULATION GROWTH RATE  #########
    proj.lambda[ns]<-exp((1/10)*sum(loglambda.f[ns,(ncountyears):(ncountyears+nprojyears-1)]))   # Geometric mean
  	
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
parameters <- c("proj.lambda","Ntot")


# MCMC settings
ni <- 75000
nt <- 5
nb <- 25000
nc <- 4


# Call JAGS from R
ipm.model <- jags(jags.data,
                  inits,
                  parameters,
                  "models/AQWA.IPM.source.pop.jags",
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
write.table(out, "output/AQWA_source_pop_scenarios_output.csv", sep=",")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Figure 1: POPULATION TRAJECTORY IN FUTURE BY NUMBER OF NESTS REMOVED ----------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ntotdf <- data.frame(Ntot = numeric(length = length(ipm.model$mean$Ntot)),
                     cip = numeric(length = length(ipm.model$mean$Ntot)),
                     cim = numeric(length = length(ipm.model$mean$Ntot)),
                     Scenario = rep(seq(0,15,1), each=47),
                     year = rep(1:47, nscenarios)) %>%
  mutate(year=year+2002)

ntotdf$Ntot <- as.numeric(t(ipm.model$mean$Ntot))
ntotdf$cim <- as.numeric(t(ipm.model$q2.5$Ntot))
ntotdf$cip <- as.numeric(t(ipm.model$q97.5$Ntot))



## gradient shaded ribbon is too much effort: https://stackoverflow.com/questions/4913117/gradient-in-geom-ribbon

poptrends <- 
  ntotdf %>%

  ggplot(aes(x = year, y = Ntot, col = as.factor(as.character(Scenario)))) +
  geom_line(linewidth = 1.1) +
  geom_vline(aes(xintercept = 2025), linetype = 2) +
  geom_ribbon(aes(ymin = cim, ymax = cip,fill=as.factor(as.character(Scenario))), colour=NA, alpha=0.3) +
  facet_wrap(~Scenario, ncol = 4, scales="free_y")  +
  theme_bw() +
  theme(legend.position = "none",
        strip.text=element_text(size=10)) +
  # theme(legend.position="bottom",
  #     legend.direction="horizontal") +
  xlab("Year") + ylab("Aquatic Warbler population size")
poptrends
ggsave("output/Source_pop_scenario_projections.jpg", width=241,height=241, quality=100, units="mm")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CALCULATE PROBABILITY OF DECLINE ----------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## plot histograms of future growth rates

as_tibble(rbind(ipm.model$samples[[1]],ipm.model$samples[[2]],ipm.model$samples[[3]],ipm.model$samples[[4]])) %>%
  dplyr::select(tidyselect::starts_with("proj.lambda")) %>%
  gather(key="Scenario",value="lambda") %>%
  mutate(N_nests=as.numeric(str_extract(Scenario, "(\\d)+"))-1) %>%
  mutate(N_nests=factor(N_nests, levels=as.character(seq(0,15,1)))) %>%
  

  ### start the plot ###
  ggplot(aes(x = lambda, fill = N_nests)) +                       # Draw overlaying histogram
  geom_histogram(position = "identity", alpha = 0.1, bins = 80, aes(y = after_stat(density)), color="black") +
  geom_density(alpha=0.3) +
  #geom_vline(aes(xintercept = 1), colour="indianred3", size=1) +
  geom_segment(aes(x = 1, y = 0, xend = 1, yend = 17), colour="firebrick", linetype = "dashed", linewidth=1)+
  
  labs(x="Future population growth rate", y="Probability density",
       fill="Number of nests removed annually") +
  scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) +
  scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title = element_text(size=14),
        legend.position="inside",
        legend.position.inside=c(0.18,0.72),
        panel.grid.major = element_line(size=.1, color="grey94"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"))


ggsave("output/Source_pop_future_pop_growth_rates.jpg", width=235,height=181, quality=100, units="mm")



#### TABLE OF PROBABILITY OF DECLINE

prob_decline<-as_tibble(rbind(ipm.model$samples[[1]],ipm.model$samples[[2]],ipm.model$samples[[3]],ipm.model$samples[[4]])) %>%
  dplyr::select(tidyselect::starts_with("proj.lambda")) %>%
  gather(key="Scenario",value="lambda") %>%
  mutate(N_nests=as.numeric(str_extract(Scenario, "(\\d)+"))-1) %>%
  mutate(decline=ifelse(lambda<1,1,0)) %>%
  group_by(N_nests) %>%
  summarise(prob=mean(decline))

fwrite(prob_decline,"AQWA_source_pop_prob_decline_nest_removal.csv")
