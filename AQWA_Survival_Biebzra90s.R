###Survival estimation using the data from Dyrcz & Zdunek, 1993, Ibis, "Breeding ecology of the Aquatic Warbler Acrocephalus paludicola on the Biebrza marshes, northeast Poland####

#Remember: all ringed birds were adults

marr <- matrix(nrow = 5, ncol = 6, data = 0)
rownames(marr) <- c("1986", "1987", "1988", "1989", "1990")
colnames(marr) <- c("1987", "1988", "1989", "1990", "1991", "Missed")

marr[1,] <- c(3,1,1,1,0,2)
marr[2,] <- c(0,1,0,0,0,15)
marr[3,] <- c(0,0,10,5,2,13)
marr[4,] <- c(0,0,0,10,6,44)
marr[5,] <- c(0,0,0,0,6,74)

#Bundle data
jags.data <- list(marr=marr, rel=rowSums(marr), nyears=ncol(marr))

cat(file = "dyrcz.mult", "
    model{
    
    #Priors and lineal models
    phi.mean ~ dbeta(1.2,1.2)
    p.mean ~ dbeta(1.2,1.2)
    
    for (t in 1:(nyears-1)){
      phi[t] <- phi.mean
      p[t] <- p.mean
  }
  # Likelihood
  # Define the multinomial likelihood
  for (t in 1:(nyears-1)){
    marr[t,1:nyears] ~ dmulti(pi[t,], rel[t])
  }
  # Define the cell probabilities of the m-array
  for (t in 1:(nyears-1)){
  # Main diagonal
  q[t] <- 1 - p[t] # Probability of non-recapture
  pi[t,t] <- phi[t] * p[t]
  # Above main diagonal
  for (j in (t+1):(nyears-1)){
    pi[t,j] <- prod(phi[t:j]) * prod(q[t:(j-1)]) * p[j]
  } #j
  # Below main diagonal
  for (j in 1:(t-1)){
    pi[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(nyears-1)){
    pi[t,nyears] <- 1-sum(pi[t,1:(nyears-1)])
  }
}
")
    
# Initial values
inits <- function(){list(phi.mean=0.4)}
# Parameters monitored
parameters <- c("phi.mean", "p.mean")
# MCMC settings
ni <- 3000; nb <- 1000; nc <- 3; nt <- 1; na <- 1000
# Call JAGS (ART <1 min), check convergence and summarize posteriors
out19 <- jags(jags.data, inits, parameters, "dyrcz.mult", n.iter=ni, n.burnin=nb, n.chains=nc,
              n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out19) # Not shown

print(out19)

#           mean    sd   2.5%    50%  97.5% overlap0 f  Rhat n.eff
#phi.mean  0.775 0.097  0.586  0.775  0.958    FALSE 1 1.000  6000
#p.mean    0.217 0.050  0.135  0.211  0.331    FALSE 1 1.000  6000
#deviance 60.980 1.978 59.057 60.393 66.327    FALSE 1 1.001  4114