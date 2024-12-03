###################################################################################################
########   AQUATIC WARBLER SURVIVAL ANALYSIS TO CONTRIBUTE TO PVA ################
###################################################################################################
# to be added to the IPM once we get this to run

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

sink("models/AQWA.surv.jags")
cat("
    model {
      # Priors and constraints
      for (i in 1:nind){
        for (t in f[i]:(n.occasions-1)){
          phi[i,t] <- beta[age[i,t], sex[i]]
          #p[i,t] <- mean.p
        } #t
      } #i
      
      for (ag in 1:2){
        for(sx in 1:2) {
          beta[ag,sx] ~ dunif(0.2, 0.6)              # Priors for age- and sex-specific survival
        }
      }
      
      mean.p ~ dunif(0.05, 0.75)                  # resighting probability is the same for all birds
      
      # Likelihood 
      for (i in 1:nind){
        # Define latent state at first capture
        z[i,f[i]] <- 1
        for (t in (f[i]+1):n.occasions){
          # State process
          z[i,t] ~ dbern(mu1[i,t])
          mu1[i,t] <- phi[i,t-1] * z[i,t-1]
          # Observation process
          y[i,t] ~ dbern(mu2[i,t])
          mu2[i,t] <- mean.p * z[i,t]
        } #t
      } #i
    }
    ",fill = TRUE)
sink()

# Initial values

# Function to create a matrix of initial values for latent state z
cjs.init.z <- function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2 <- max(which(ch[i,]==1))
    ch[i,f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]){
    ch[i,1:f[i]] <- NA
  }
  return(ch)
}

inits <- function(){list(
  z= cjs.init.z(as.matrix(AW_CH[,3:7]), f),
  beta= matrix(runif(4, 0.4,0.6),ncol=2),
  mean.p = runif(1, 0.2,0.7))}

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
                  "models/AQWA.surv.jags",
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




