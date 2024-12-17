###################################################################################################
########   AQUATIC WARBLER SURVIVAL ANALYSIS TO CONTRIBUTE TO PVA ################
###################################################################################################
# to be added to the IPM once we get this to run

## updated on 4 Dec 2024 to explore various model formulations

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

AW_CH<-read_excel("data/AW_raw_data.xlsx", sheet="returns") %>%
  rename(ID=`ring ID`) %>%
  mutate(sex=if_any(everything(), ~str_detect(tolower(.), "female"))) %>%
  mutate(sex=ifelse(is.na(sex),1,2)) %>%
  mutate(`2018`=ifelse(`release year`==2018,1,0)) %>%
  mutate(`2019`=ifelse(`release year`==2019,1,`2019`)) %>%
  select(ID, sex,`2018`,`2019`,`2020`,`2021`,`2022`)

head(AW_CH)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
f <- apply(AW_CH[,3:7], 1, get.first)

# Create age matrices X indicating age classes
age.mat<-as.matrix(AW_CH[,3:7])
age.mat[,2]<-ifelse(age.mat[,1]==1,2,1)
age.mat[,3:5]<-2


### prepare the data
jags.data <- list(nind=dim(AW_CH)[1],
                  n.occasions=dim(AW_CH[,3:7])[2],
                  age=age.mat,
                  f=f,
                  sex=as.numeric(AW_CH$sex),
                  y = as.matrix(AW_CH[,3:7]))




##############################################################################
#
# simple CJS model to estimate survival --------
# 
##############################################################################

### tried 3 model formulations: constant, sex.det, and age.sex, and DIC is lowest for age.sex


sink("models/AQWA.surv.age.sex.det.beta.jags")
cat("
    model {
      # Priors and constraints
      for (i in 1:nind){
        for (t in f[i]:(n.occasions-1)){
          phi[i,t] <- beta[age[i,t]]
          p[i,t] <- mean.p[sex[i]]
        } #t
        p[i,n.occasions] <- mean.p[sex[i]]
      } #i
      
      for (ag in 1:2){
        beta[ag] ~ dbeta(1.2, 1.2)              # Priors for age- and sex-specific survival
        # for(sx in 1:2) {
        #   beta[ag,sx] ~ dunif(0.2, 0.6)              # Priors for age- and sex-specific survival
        # }
      }
        for(sx in 1:2) {
          mean.p[sx] ~ dbeta(1.2, 1.2)                  # resighting probability is the same for all birds
        }
      
      
      # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
      # State process
      z[i,t] ~ dbern(z[i,t-1] * phi[i,t-1])
      # Observation process
      y[i,t] ~ dbern(z[i,t] * p[i,t])
    } #t
  } #i
    }
    ",fill = TRUE)
sink()

# Initial values

inits <- function(){list(
  z= zInit(as.matrix(AW_CH[,3:7])),
  #beta= matrix(runif(4, 0.4,0.59),ncol=2),
  beta= runif(2, 0.4,0.59),
  #mean.phi= runif(1, 0.4,0.59),
  mean.p = runif(2, 0.2,0.7))}

# Parameters monitored
parameters <- c("beta","mean.p")


# MCMC settings
ni <- 7500
nt <- 5
nb <- 2500
nc <- 4


# Call JAGS from R
surv.model <- jags(jags.data,
                   inits,
                   parameters,
                   "models/AQWA.surv.age.sex.de.beta.jags",
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.cores=nc, parallel=T)




############################################################################
#
# WRITE ALL OUTPUT INTO A TEXT FILE ----
# 
##############################################################################
out<-as.data.frame(surv.model$summary)  
out$parameter<-row.names(surv.model$summary)
names(out)[c(12,5,3,7)]<-c('parm','median','lcl','ucl')
print(surv.model, dig=3)
write.table(out, "output/AQWA_surv_model_output.csv", sep=",")


#Uniform prior
beta[1]    0.391  0.089  0.234  0.388  0.572    FALSE 1 1.001  2452
beta[2]    0.309  0.076  0.206  0.296  0.495    FALSE 1 1.000  4000
mean.p[1]  0.631  0.091  0.415  0.650  0.745    FALSE 1 1.001  2524
mean.p[2]  0.231  0.116  0.074  0.209  0.534    FALSE 1 1.000  4000
deviance  59.551 11.753 38.548 58.933 82.994    FALSE 1 1.000  4000

#Beta prior
mean     sd   2.5%    50%  97.5% overlap0 f  Rhat n.eff
beta[1]    0.362  0.111  0.195  0.343  0.638    FALSE 1 1.003  2261
beta[2]    0.283  0.090  0.129  0.275  0.480    FALSE 1 1.000  4000
mean.p[1]  0.737  0.141  0.418  0.756  0.956    FALSE 1 1.002  1067
mean.p[2]  0.271  0.147  0.081  0.241  0.662    FALSE 1 1.000  3809
deviance  50.788 15.725 24.738 49.159 86.911    FALSE 1 1.002  2003

