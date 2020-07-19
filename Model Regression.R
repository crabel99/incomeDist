# Define functions
# AIC compensated for small sample sizes
AICc <- function(k, n, lpd) {
  # k - number of coefficients
  # n - number of data points
  # lpd - log likelihood 
  f <-  2 * k * (n + k)/(n - k - 1) - 2 * lpd
  return(f)
}

# Load Data
load("./Table 2.1.RData")

# Load Libraries
library(rstan)
library(StanHeaders)
library(BH)
library(RcppEigen)
library(Rcpp)
library(inline)
library(fitdistrplus)

rstan_options(auto_write = TRUE)

expose_stan_functions("./gibbs_model.stan")
expose_stan_functions("./maxwell_model.stan")

# Configure data for processing
years <- names(totalIncome)
years <- years[-1]

# Initialization data
inits <- list(
  list ('sigma' = 100, 'T' = 10000, 'r0' = 50000, 'alpha' = 1),
  list ('sigma' = 5000, 'T' = 20000, 'r0' = 75000, 'alpha' = 1.4),
  list ('sigma' = 10000, 'T' = 40000, 'r0' = 100000, 'alpha' = 1.6),
  list ('sigma' = 100, 'T' = 60000, 'r0' = 150000, 'alpha' = 2))
nChains <- length(inits)
nCores <- nChains

# Output prototypes
fitEval <- data.frame(aicGibbs = double(),
                      aicMaxwell = double(),
                      sseGibbs = double(),
                      sseMaxwell = double())
regGibbs <- data.frame("C" = double(),
                       "T" = double(),
                       "r0" = double(),
                       "alpha" = double())
regMaxwell <- regGibbs
pdGibbs <- data.frame("sigma.shape"=double(),
                      "sigma.rate"=double(),
                      "T.shape"=double(),
                      "T.rate"=double(),
                      "r0.shape"=double(),
                      "r0.rate"=double(),
                      "alpha.shape"=double(),
                      "alpha.rate"=double())
pdMaxwell <- pdGibbs
params <- c("sigma", "T", "r0", "alpha", "lpd") # output parameters

# Loop through each of the years
for (i in 1:length(years)) {
  # The IRS added additional income brackets in 2000
  if (as.numeric(years[i]) < 2000) {
    bins <- 18
  } else {
    bins <- 22
  }
  y_hat <- vector(mode = "numeric", length = bins)
  incomeData <- list ('N' = bins,
                      'x' = totalIncome[1:bins, 1],
                      'y' = totalIncome[1:bins, i + 1])
  y <- incomeData$y/sum(incomeData$y)
  
  # Process the data
  fitGibbs <- stan(file = './gibbs_model.stan',
                   data = incomeData,
                   chains = nChains,
                   init = inits,
                   pars = params,
                   cores = nCores)
  fitMaxwell <- stan(file = './maxwell_model.stan',
                     data = incomeData,
                     chains = nChains,
                     init = inits,
                     pars = params,
                     cores = nCores)
  
  # Condition the data for the modified Gibbs distribution
  draws <- extract(fitGibbs, pars = params)
  reg <- data.frame(row.names = years[i],
                    "T" = mean(draws$T),
                    "r0" = mean(draws$r0),
                    "alpha" = mean(draws$alpha))
  for (j in 1:bins) {
    if (j < bins) {
      y_hat[j] <- gibbs_int(c(totalIncome$Min[j],
                              totalIncome$Min[j + 1]),
                            as.numeric(reg[1,]),
                            double())
    } else {
      y_hat[j] <- gibbs_int(c(totalIncome$Min[j], Inf),
                              as.numeric(reg[1,]),
                              double())
    }
  }
  reg$C <- sum(y_hat)
  y_hat <- y_hat/reg$C
  fit <- data.frame(row.names = years[i], sseGibbs = sum((y - y_hat)^2))
  fit$aicGibbs <- AICc(3, bins, mean(draws$lpd))
  
  # Determine the posterior-prior information distribution
  prior <- list(sigma = fitdist(as.vector(draws$sigma),"gamma"),
                T = fitdist(as.vector(draws$T),"gamma", method="mme"),
                r0 = fitdist(as.vector(draws$r0),"gamma", method="mme"),
                alpha = fitdist(as.vector(draws$alpha),"gamma"))
  pd <- data.frame(row.names = years[i],
                   "sigma.shape" = prior$sigma$estimate[1],
                   "sigma.rate" = prior$sigma$estimate[2],
                   "T.shape" = prior$T$estimate[1],
                   "T.rate" = prior$T$estimate[2],
                   "r0.shape" = prior$r0$estimate[1],
                   "r0.rate" = prior$r0$estimate[2],
                   "alpha.shape" = prior$alpha$estimate[1],
                   "alpha.rate" = prior$alpha$estimate[2])
  
  # export data
  regGibbs <- rbind(regGibbs, reg)
  pdGibbs <- rbind(pdGibbs, pd)
  
  # Condition the data for the modified Maxwell distribution
  draws <- extract(fitMaxwell, pars = params)
  reg <- data.frame(row.names = years[i],
                    "T" = mean(draws$T),
                    "r0" = mean(draws$r0),
                    "alpha" = mean(draws$alpha))
  for (j in 1:bins) {
    if (j < bins) {
      y_hat[j] <- maxwell_int(c(totalIncome$Min[j],
                                totalIncome$Min[j + 1]),
                              as.numeric(reg[1,]),
                              double())
    } else {
      y_hat[j] <- maxwell_int(c(totalIncome$Min[j], Inf),
                              as.numeric(reg[1,]),
                              double())
    }
  }
  reg$C <- sum(y_hat)
  y_hat <- y_hat/reg$C
  fit$sseMaxwell <- sum((y - y_hat)^2)
  fit$aicMaxwell <- AICc(3, bins, mean(draws$lpd))

  # Determine the posterior-prior information distribution
  prior <- list(sigma = fitdist(as.vector(draws$sigma),"gamma"),
                T = fitdist(as.vector(draws$T),"gamma", method="mme"),
                r0 = fitdist(as.vector(draws$r0),"gamma", method="mme"),
                alpha = fitdist(as.vector(draws$alpha),"gamma"))
  pd <- data.frame(row.names = years[i],
                   "sigma.shape" = prior$sigma$estimate[1],
                   "sigma.rate" = prior$sigma$estimate[2],
                   "T.shape" = prior$T$estimate[1],
                   "T.rate" = prior$T$estimate[2],
                   "r0.shape" = prior$r0$estimate[1],
                   "r0.rate" = prior$r0$estimate[2],
                   "alpha.shape" = prior$alpha$estimate[1],
                   "alpha.rate" = prior$alpha$estimate[2])
  
  # export data
  regMaxwell <- rbind(regMaxwell, reg)
  pdMaxwell <- rbind(pdMaxwell, pd)
  fitEval <- rbind(fitEval, fit)
  
}

# Clean up the workspace
rm(prior, fit, y_hat, reg, draws, inits, y, incomeData, years, nChains, params)
rm(pd, i, j, bins, AICc, gibbs_dist, gibbs_int, maxwell_dist, maxwell_int, nCores)

# Uncomment to delete the rstan output
# rm(fitGibbs, fitMaxwell)

# Save the workspace
save.image("./model_output.RData")
