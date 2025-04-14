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

N1[1,1] ~ dnorm(60, 0.0001)T(0,)           # 1-year old individuals
NadSurv[1,1] ~ dnorm(100, 0.0001)T(0,)      # Adults >= 2 years

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
   
            Nfem.breed1[1,t] <- (Ntot[1,t-1]*(1-prop.males))
            Nfem.breed2[1,t] ~ dbin(db[1,t],round(Nfem.breed1[1,t]))

      	    chicks[1,t-1] <- Nfem.breed1[1,t] * fec1[t] + Nfem.breed2[1,t]*fec2[t]			# total fecundity
      	    chicksrd[1,t-1] <- round(chicks[1,t-1])
        
      	# chicks[1,t-1] <- (Ntot[1,t-1])*(1-prop.males)* fec1[t-1] # total fecundity
      	# chicksrd[1,t-1] <- round(chicks[1,t-1])
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
  for(t in 2:(ncountyears)){
  
    Nfem.breed1[ns,t] <- Nfem.breed1[1,t]
    Nfem.breed2[ns,t] <- Nfem.breed2[1,t]
    # chicks[ns,t] <- chicks[1,t]
    # chicksrd[ns,t] <- chicksrd[1,t]
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

      	    chicks[ns,t-1] <- Nfem.breed1[ns,t-1] * fec1[t] + Nfem.breed2[ns,t-1]*mfec2			# total fecundity
      	    chicksrd[ns,t-1] <- round(chicks[ns,t-1])
		        N1[ns,t] ~ dbin(phij[t-1],chicksrd[ns,t-1]) 
		        NadSurv[ns,t] ~ dbin(phi.ad,round((Ntot[ns,t-1])))
		   	    Ntot[ns,t] <- NadSurv[ns,t] + N1[ns,t]								## total population	
  	} #Time loop
  	
  	
  	## FINAL POPULATION GROWTH RATE (steady state after reintroduction ends)
  	# Population growth rate
    for (t in (ncountyears+2):(ncountyears+nprojyears-1)){
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
write.table(out, "output/AQWA_GER_model_nscenarios_output.v13.csv", sep=",")




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
                     survival=rep(c("no survival improvement","5% survival increase","5% survival decrease"), each=nscenarios*47),
                     year = rep(1:47, nscenarios*3))

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
ggsave("output/Past_population_trajectories.jpg", width=191,height=241, quality=100, units="mm")



#The plot for the future
## gradient shaded ribbon is too much effort: https://stackoverflow.com/questions/4913117/gradient-in-geom-ribbon

poptrends <- 
  ntotdf %>%
  filter(!(Scenario %in% c("Past trajectory without mowing", "Past trajectory with mowing"))) %>%
  mutate(survival=ifelse(survival=="no survival improvement","no survival change",survival)) %>%
  

  ggplot(aes(x = year+2002, y = Ntot, col = survival, linetype=survival, group = survival)) +
  geom_line(linewidth = 1.1) +
  geom_vline(aes(xintercept = 2025), linetype = 2) +
  #geom_ribbon(aes(ymin = cim, ymax = cip,linetype=survival), fill = NA) +
  geom_ribbon(aes(ymin = cim, ymax = cip,fill=survival), colour=NA, alpha=0.3) +
  facet_wrap(~Scenario, ncol = 3, scales="free_y")  +
  theme_bw() +
  #theme(legend.position = "none") +
  theme(legend.position="bottom",
      legend.direction="horizontal") +
  xlab("Year") + ylab("Aquatic Warbler population size")
poptrends
#ggsave("output/Scenario_projections_v13.jpg", width=241,height=241, quality=100, units="mm")



##### ALTERED THIS TO PROBABILITY OF TARGET POP SIZE OF 100
# Extinction probability (could plot this as well)
## revised to show prop samples where Ntot>100

w0 <- function(x){
  ext <- sum(x > 99) / length(x)
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
                   "Past trajectory without mowing","Past trajectory with mowing"), 3)

habitat<- rep(c(10000, 10000, 10000,
                10000, 10000,10000,
                       200, 200,200,
                       400, 400,400,
                       1200, 1200,1200,
                       2400, 2400,2400,
                       9999,9999), 3)

rel.years<- rep(c(5, 10, 15,
                       5, 10,15,
                       5, 10,15,
                       5, 10,15,
                       5, 10,15,
                       5, 10,15,
                       0,0), 3)


target.prob<- fut.ext %>%
  pivot_longer(names(fut.ext), names_to = "Scenario", values_to = "N") %>%
  #gather(key="Scenario",value="N") %>%
  mutate(habitat=rep(habitat,40000)) %>%
  mutate(releases=rep(rel.years,40000)) %>%
  mutate(survival=rep(rep(c("no survival change","5% survival increase","5% survival decrease"), each=(nscenarios)),40000)) %>%
  #filter(!(str_detect(pattern="15y", string=Scenario))) %>%
  filter((str_detect(pattern="No mowing", string=Scenario))) %>%
  filter(!(Scenario %in% c("Past trajectory without mowing", "Past trajectory with mowing"))) %>%
  group_by(Scenario,survival,habitat,releases) %>%
  summarize(ext.prob=w0(N))
  # mutate(habitat=as.factor(habitat)) %>%
  # ungroup() %>%
  #mutate(habitat=ifelse(habitat==3000,"unlimited",habitat)) %>%
  #filter(!ext.prob==0) %>%
  
  ### start the plot ###
  #ggplot(aes(x = Scenario, y=ext.prob, fill = Scenario)) +                       # Draw overlaying histogram
  #geom_bar(stat="identity",alpha = 0.6) +
  ggplot(data=target.prob,aes(x = as.factor(habitat), y=ext.prob, fill = as.factor(releases), colour = as.factor(releases))) +                       # Draw overlaying histogram
  #geom_line(linewidth=1.5) +
  geom_bar(stat="identity",position=position_dodge(width=0.3), width=0.2) +
  facet_wrap(~survival, ncol = 1)  +
  labs(x="Extent of habitat (ha)", y="Probability of target population size after 25 years",
       fill="N of release years", colour="N of release years") +
  #scale_fill_discrete() +
  #guides(fill = guide_legend(nrow = 5)) +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=12, color="black"),
        #axis.text.x=element_blank(), 
        axis.title=element_text(size=14),
        strip.text=element_text(size=14, color="black"),
        strip.background=element_rect(fill="white", colour="black"),
        legend.text=element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="inside",
        legend.position.inside=c(0.12,0.9),
        legend.background=element_blank(),
        #legend.direction="horizontal",
        panel.grid.major = element_line(linewidth=.1, color="grey94"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"))

write.table(target.prob, "output/AQWA_Target100_pop_size_probability_25y.csv", sep=",")
ggsave("output/Target_pop_size_probability.jpg", width=351,height=241, quality=100, units="mm")

################################################################################
#######EXTRACTING PROBABILITY OF DECREASE 10 YEARS AFTER END OF RELEASES########
################################################################################

#On the projections: Calculations accounting for uncertainty

samples <- ipm.model$sims.list$Ntot #40000 samples, 60 scenarios, 47 years

#End of releases (releases start at t = 27 + 1).  
endrel <- rep(c(5,10,15), 18) + 22
endrel10 <- endrel + 10

#Remove past scenarios from the samples 
samplesr <- samples[,-c(19, 20, 39, 40, 59,60),]

nsc <- 54
ndraws <- 40000

calcs <- array(data = NA, dim = c(ndraws, nsc, 2)) #Nscenarios x Nsamples x Nkey (end releases vs end projections)

for(i in 1:ndraws){
  for(j in 1:nsc){
    
    calcs[i,j,1] <- samplesr[i, j, endrel[j]]
    calcs[i,j,2] <- samplesr[i, j, endrel10[j]]
    
  }
}




subcalcs <- calcs[,,2] - calcs[,,1] #This is the absolute difference, that is, values above 0 mean the population has increased 10 years after release, values below 0 mean otherwise.

#Now, let's extract probabilities
probdecline <- numeric(length = nsc)

for(i in 1:nsc)
  probdecline[i] <- 1- (sum(subcalcs[,i] >= 0) / ndraws)


rescalcs <- data.frame(Scenario = rep(levels(ntotdf$Scenario)[1:18],3),
                       survival = rep(c("no survival change","5% survival increase","5% survival decrease"), each=18),
                       probdecline = probdecline,
                       mean = apply(subcalcs, 2, mean),
                       cimin = apply(subcalcs, 2, quantile, probs = c(0.025)),
                       cimax = apply(subcalcs, 2, quantile, probs = c(0.975)))


rescalcs %>%
  #pivot_longer(names(fut.ext), names_to = "Scenario", values_to = "N") %>%
  #gather(key="Scenario",value="N") %>%
  mutate(habitat=habitat[c(1:18,21:38,41:58)]) %>%
  mutate(releases=rel.years[c(1:18,21:38,41:58)]) %>%
  filter((str_detect(pattern="No mowing", string=Scenario))) %>%

  ### start the plot ###
  ggplot(aes(x = as.factor(habitat), y=probdecline, fill = as.factor(releases), colour = as.factor(releases))) +                       # Draw overlaying histogram
  geom_bar(stat="identity",position=position_dodge(width=0.3), width=0.2) +
  facet_wrap(~survival, ncol = 1)  +
  labs(x="Extent of habitat (ha)", y="Probability of population decline 10 years after end of reinforcement",
       fill="N of release years", colour="N of release years") +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=12, color="black"),
        axis.title=element_text(size=14),
        strip.text=element_text(size=14, color="black"),
        strip.background=element_rect(fill="white", colour="black"),
        legend.text=element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="bottom",
        legend.background=element_blank(),
        legend.direction="horizontal",
        panel.grid.major = element_line(linewidth=.1, color="grey94"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"))


ggsave("output/Future_decrease_probability.jpg", width=351,height=241, quality=100, units="mm")
write.table(rescalcs, "output/AQWA_GER_decrease_probabilities_10y_post_release.csv", sep=",")




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
         no.mowing.curr.surv=`proj.lambda[4]`,
         mowing.imp.surv=`proj.lambda[21]`,
         # no.mowing.lim.hab.200=`proj.lambda[27]`,
         # no.mowing.lim.hab.1200=`proj.lambda[33]`,
         no.mowing.imp.surv=`proj.lambda[24]`,
         no.mowing.red.surv=`proj.lambda[44]`) %>%
  dplyr::select(-tidyselect::starts_with("proj.lambda")) %>%
  gather(key="Scenario",value="lambda") %>%
  filter(lambda>0.5) %>% ## remove bizarre growth rates caused by extinct populations?
  mutate(Scenario=if_else(Scenario=="mowing.imp.surv", "improved survival, no second broods",
                          if_else(Scenario=="nochange", "current survival, no second broods",
                          if_else(Scenario=="no.mowing.red.surv", "reduced survival and second broods",
                          if_else(Scenario=="no.mowing.curr.surv","current survival and second broods","improved survival and second broods"))))) %>%
  
  ### start the plot ###
  ggplot(aes(x = lambda, fill = Scenario)) +                       # Draw overlaying histogram
  geom_histogram(position = "identity", alpha = 0.2, bins = 80, aes(y = after_stat(density)), color="black") +
  geom_density(alpha=0.5) +
  #geom_vline(aes(xintercept = 1), colour="indianred3", size=1) +
  geom_segment(aes(x = 1, y = 0, xend = 1, yend = 12), colour="firebrick", linetype = "dashed", linewidth=1)+
  
  labs(x="Future population growth rate", y="Probability density",
       fill="Scenario (with unlimited habitat)") +
  scale_fill_viridis_d(alpha=0.3,begin=0,end=0.98,direction=1) +
  scale_color_viridis_d(alpha=1,begin=0,end=0.98,direction=1) +
  
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20),
        legend.text=element_text(size=12),
        legend.title = element_text(size=14),
        legend.position="inside",
        legend.position.inside=c(0.18,0.82),
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




