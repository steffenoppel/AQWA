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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. LOAD DATA (prepared in script AQWA_BMM_data_preparation_multiyear.r) ----------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

load("data\\AQWA_trend_input.RData")
y<-AQWA.y
head(y)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2. FORMULATE TREND MODEL AND ANALYSE DATA ----------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## 2.1. Specify model in BUGS language ------
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


## 2.2. Bundle data, initial values and specify monitored parameters -----------
R = nrow(y)
T = ncol(y)
nyears = dim(y)[3]
jags.data <- list(M = y, nsite = R, nrep = T, nyear=nyears, primocc=seq(1:nyears))

# Initial values
Nst <- apply(ifelse(is.na(y),0,y), c(1, 3), max) + 1
Nst[is.na(Nst)] <- 1
inits <- function(){list(N = Nst, loglam = runif(1, -3, 3), sigma.site = runif(1, 0, 1), p0 = runif(nyears, 0, 1))}

# Parameters monitored
parameters <- c("loglam","lam.site","totalN","N", "fit", "fit.new","trend")

# MCMC settings
ni <- 50000				## fit statistic needs >400000 to converge!!
nt <- 40
nb <- 10000
nc <- 3


## 2.3. RUN MODEL IN JAGS  -----------
trend.model <- jags(jags.data,
                  inits,
                  parameters,
                  "models/NmixTrend.jags",
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.cores=nc, parallel=T)



## 2.4. SUMMARISE OUTPUT AND EVALUATE MODEL  -----------------------------------

# Evaluation of fit
plot(trend.model$sims.list$fit, trend.model$sims.list$fit.new, main = "", xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", frame.plot = FALSE)
abline(0, 1, lwd = 2, col = "black")
mean(trend.model$sims.list$fit.new > trend.model$sims.list$fit)
mean(trend.model$mean$fit) / mean(trend.model$mean$fit.new)

# Summarize posteriors
print(trend.model, dig = 2)




## 2.5. CHECK CONVERGENCE AND WRITE ALL OUTPUT INTO A TEXT FILE ----------------

out<-as.data.frame(trend.model$summary)  
out$parameter<-row.names(trend.model$summary)
names(out)[c(12,5,3,7)]<-c('parm','median','lcl','ucl')
print(trend.model, dig=3)
out %>% arrange(desc(Rhat)) %>% select(parm, median, lcl, ucl, Rhat, n.eff)
write.table(out, "output/AQWA_Biebrza_trend_model_output.csv", sep=",")




## 2.6. FIGURE OF POPULATION TRAJECTORY  ---------------------------------------

data.frame(Ntot = trend.model$mean$totalN,
                     cip = trend.model$q2.5$totalN,
                     cim = trend.model$q97.5$totalN,
                     year=seq(2011,2025,1)) %>%


  
  ggplot(aes(x = year, y = Ntot)) +
  geom_line(linewidth = 1.1, col="firebrick") +
  geom_ribbon(aes(ymin = cim, ymax = cip), linetype = 2, fill="firebrick", alpha=0.2) +
  xlab("Year") + ylab("Aquatic Warbler population size") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=14, color="black"),
        axis.title=element_text(size=16),
        panel.grid.major = element_line(linewidth=.1, color="grey94"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"))



ggsave("output/AQWA_Biebrza_trend_2011_2025.jpg", width=291,height=201, quality=100, units="mm")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. COMPARE RAW DATA WITH ESTIMATES ----------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AW %>% ##dplyr::filter(Year==2015) %>% dplyr::filter(Transect_no=="TBB74")
  group_by(Year,Transect_no, Control) %>%
  summarise(N=sum(No_singing_males),day=first(yday)) %>%
  ungroup() %>%
  group_by(Year,Transect_no) %>%
  summarise(Max=mean(N)) %>% ##print(n=80)
  ungroup() %>%
  group_by(Year) %>%
  #summarise(mean=mean(Max,na.rm=T),lcl=quantile(Max,0.025, na.rm=T),ucl=quantile(Max,0.975,na.rm=T)) %>%
  summarise(mean=mean(Max,na.rm=T),sd=sd(Max,na.rm=T)) %>% mutate(lcl=mean-0.5*sd,ucl=mean+0.5*sd, Data="Raw data (Mean N per transect)") %>%
  
  bind_rows(
    data.frame(mean = trend.model$mean$totalN,
               lcl = trend.model$q2.5$totalN,
               ucl = trend.model$q97.5$totalN,
               Year=seq(2011,2025,1),
               Data="Model estimates (Total Number)")) %>%
  
  
  ggplot(aes(x = Year, y = mean, col=Data, fill=Data)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), linetype = 2, alpha=0.2) +
  scale_x_continuous(breaks=seq(2011,2025,2), labels=seq(2011,2025,2)) +
  facet_wrap(~Data, ncol=1, scales="free_y") +
  xlab("Year") + ylab("Aquatic Warbler abundance") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=14, color="black"),
        axis.title=element_text(size=16),
        strip.text=element_text(size=14, color="black"),
        strip.background=element_rect(fill="white", colour="black"),
        legend.text=element_text(size=12),
        legend.title = element_text(size=14),
        legend.position="none",
        legend.position.inside=c(0.18,0.82),
        panel.grid.major = element_line(linewidth=.1, color="grey94"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, colour = "black"))


ggsave("output/AQWA_Biebrza_trend_2011_2025.jpg", width=201,height=291, quality=100, units="mm")



