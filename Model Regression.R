# Define functions
# AIC compensated for small sample sizes
AICc <- function(k, n, lpd) {
  # k - number of coefficients
  # n - number of data points
  # lpd - log likelihood 
  f <-  2 * k * n / (n - k - 1) - 2 * lpd
  return(f)
}

# Load Libraries
suppressPackageStartupMessages({
  library(rstan)
  library(BH)
  library(RcppEigen)
  library(Rcpp)
  library(inline)
  library(fitdistrplus)
  library(modeest)
})

# Load Data
# Select the appropriate table either "1.1", "1.1 TI", or "2.1"
table <- "1.1 TI" 
load(paste("./Table ", table,".RData", sep = ""))

# Configure data for processing
years <- names(binReturns)
binVect <- names(binMin)

# Stan's integrate_1d doesn't like the [0,1] endpoints and chokes
if (grepl(table, "1.1") | grepl(table, "1.1 TI")) {
  for (i in 1:length(years)) {
    binReturns[[years[i]]][2] <-
      binReturns[[years[i]]][2] + binReturns[[years[i]]][1]
    binReturns[[years[i]]] <- binReturns[[years[i]]][-1]
  }
  for (i in 1:length(binVect)) {
    binMin[[binVect[i]]] <- binMin[[binVect[i]]][-2]
  }
}

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

message("Compiling the Stan code modules... ", appendLF = FALSE)
suppressMessages({
  expose_stan_functions("./gibbs_model.stan")
  expose_stan_functions("./maxwell_model.stan")
})
message("COMPLETE")

# Configure the environment
nChains <- 4
nCores <- parallel::detectCores()



# Initialization data from prior distributions
initf <- function() {
  list(
    "sigma" = rgamma(1, 2, 8e-6),
    "T" = rgamma(1, 4, 8e-5),
    "r0" = rgamma(1, 4, 4e-5),
    "alpha" = rgamma(1, 2, 1))
}

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
pdGibbs <- data.frame("sigma.shape" = double(),
                      "sigma.rate" = double(),
                      "T.shape" = double(),
                      "T.rate" = double(),
                      "r0.shape" = double(),
                      "r0.rate" = double(),
                      "alpha.shape" = double(),
                      "alpha.rate" = double())
pdMaxwell <- pdGibbs
params <- c("sigma", "T", "r0", "alpha", "lpd") # output parameters

# Loop through each of the years
for (i in 1:length(years)) {
  # The IRS data bins vary from year to year
  for (j in length(binVect):1) {
    if ( years[i] >= binVect[j]) {
      bins <- length(binReturns[[i]])
      x <- binMin[[j]]
    }
  }
  y_hat <- vector(mode = "numeric", length = bins)
  y <- binReturns[[i]]
  y <- y/sum(y)
  incomeData <- list('N' = bins,
                     'x' = x,
                     'y' = y)
  
  # print(incomeData)
  # Process the data
  message(paste("Processing tax data from:",years[i]))
  suppressMessages({
    fitMaxwell <- stan(
      file = './maxwell_model.stan',
      data = incomeData,
      chains = nChains,
      init = initf,
      pars = params,
      cores = nCores
    )
    fitGibbs <- stan(
      file = './gibbs_model.stan',
      data = incomeData,
      chains = nChains,
      init = initf,
      pars = params,
      cores = nCores
    )
  })
  
  
  # Condition the data for the modified Gibbs distribution
  draws <- extract(fitGibbs, pars = params)
  reg <- data.frame(row.names = years[i],
                    "T" = mean(draws$T),
                    "r0" = mean(draws$r0),
                    "alpha" = mean(draws$alpha))
  for (j in 1:bins) {
    if (j < bins) {
      y_hat[j] <- gibbs_int(c(x[j], x[j + 1]), as.numeric(reg[1,]), double())
    } else {
      y_hat[j] <- gibbs_int(c(x[j], Inf), as.numeric(reg[1,]), double())
    }
  }
  reg$C <- sum(y_hat)
  y_hat <- y_hat/reg$C
  fit <- data.frame(row.names = years[i], 
                    sseGibbs = sum((y - y_hat)^2))
  fit$aicGibbs <-
    AICc(3, bins, mlv(draws$lpd, method = "parzen", kernel = "gaussian"))
  
  # Determine the posterior-prior information distribution
  prior <- list(sigma = fitdist(as.vector(draws$sigma),"gamma"),
                T = fitdist(as.vector(draws$T),"gamma", method = "mme"),
                r0 = fitdist(as.vector(draws$r0),"gamma", method = "mme"),
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
      y_hat[j] <- maxwell_int(c(x[j], x[j + 1]), as.numeric(reg[1,]), double())
    } else {
      y_hat[j] <- maxwell_int(c(x[j], Inf), as.numeric(reg[1,]), double())
    }
  }
  reg$C <- sum(y_hat)
  y_hat <- y_hat/reg$C
  fit$sseMaxwell <- sum((y - y_hat)^2)
  fit$aicMaxwell <-
    AICc(3, bins, mlv(draws$lpd, method = "parzen", kernel = "gaussian"))

  # Determine the posterior-prior information distribution
  prior <- list(sigma = fitdist(as.vector(draws$sigma),"gamma"),
                T = fitdist(as.vector(draws$T),"gamma", method = "mme"),
                r0 = fitdist(as.vector(draws$r0),"gamma", method = "mme"),
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
rm(prior, fit, y_hat, reg, draws, initf, y, years, nChains, params,
   pd, i, j, bins, AICc, gibbs_dist, gibbs_int, maxwell_dist, maxwell_int, 
   nCores, binVect, x)

# Uncomment to delete the rstan output
# rm(fitGibbs, fitMaxwell)

# Save the workspace
save.image(paste("./model_output_", table, ".RData", sep = ""))
