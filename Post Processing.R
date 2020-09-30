# This script is for handling all of the post processing of the data for data
# validation. It relies on shinystan and fitdistrplus to examine the data.

# Load data
load("./model_output.RData")

# Load Libraries
suppressPackageStartupMessages({
  library(rstan)
  library(shinystan)
  library(fitdistrplus)
})

message("Compiling the Stan code modules... ", appendLF = FALSE)
suppressMessages({
  expose_stan_functions("./gibbs_model.stan")
  expose_stan_functions("./maxwell_model.stan")
})
message("COMPLETE")

# Post Processing
drawsGibbs <- extract(fitGibbs, pars = c("sigma","T","r0","alpha"))
drawsMaxwell <- extract(fitMaxwell, pars = c("sigma","T","r0","alpha"))

gibbsPrior <- list(sigma = fitdist(as.vector(drawsGibbs$sigma),"gamma"),
                    T = fitdist(as.vector(drawsGibbs$T),"gamma", method = "mme"),
                    r0 = fitdist(as.vector(drawsGibbs$r0),"gamma", method = "mme"),
                    alpha = fitdist(as.vector(drawsGibbs$alpha),"gamma"))

maxwellPrior <- list(sigma = fitdist(as.vector(drawsMaxwell$sigma),"gamma"),
                      T = fitdist(as.vector(drawsMaxwell$T),"gamma", method = "mme"),
                      r0 = fitdist(as.vector(drawsMaxwell$r0),"gamma", method = "mme"),
                      alpha = fitdist(as.vector(drawsMaxwell$alpha),"gamma"))
gibbsPD <- data.frame("sigma.shape" = gibbsPrior$sigma$estimate[1],
                       "sigma.rate" = gibbsPrior$sigma$estimate[2],
                       "T.shape" = gibbsPrior$T$estimate[1],
                       "T.rate" = gibbsPrior$T$estimate[2],
                       "r0.shape" = gibbsPrior$r0$estimate[1],
                       "r0.rate" = gibbsPrior$r0$estimate[2],
                       "alpha.shape" = gibbsPrior$alpha$estimate[1],
                       "alpha.rate" = gibbsPrior$alpha$estimate[2])
maxwellPD <- data.frame("sigma.shape" = maxwellPrior$sigma$estimate[1],
                         "sigma.rate" = maxwellPrior$sigma$estimate[2],
                         "T.shape" = maxwellPrior$T$estimate[1],
                         "T.rate" = maxwellPrior$T$estimate[2],
                         "r0.shape" = maxwellPrior$r0$estimate[1],
                         "r0.rate" = maxwellPrior$r0$estimate[2],
                         "alpha.shape" = maxwellPrior$alpha$estimate[1],
                         "alpha.rate" = maxwellPrior$alpha$estimate[2])

y <- incomeData$y/sum(incomeData$y)
sseGibbs <- sum((y - y_hatGibbs)^2)
sseMaxwell <- sum((y - y_hatMaxwell)^2)

DELTA_AIC <- AICc(3, bins, sseMaxwell/bins) - AICc(3, bins, sseGibbs/bins)

plotdist(as.vector(drawsGibbs$sigma),
         "gamma",
         para = list(
           shape = gibbsPD$sigma.shape,
           rate = gibbsPD$sigma.rate))
plotdist(as.vector(drawsGibbs$alpha),
         "gamma",
         para = list(
           shape = gibbsPD$alpha.shape,
           rate = gibbsPD$alpha.rate))
plotdist(as.vector(drawsGibbs$r0),
         "gamma",
         para = list(
           shape = gibbsPD$r0.shape,
           rate = gibbsPD$r0.rate))
plotdist(as.vector(drawsGibbs$T),
         "gamma",
         para = list(
           shape = gibbsPD$T.shape,
           rate = gibbsPD$T.rate))

plotdist(as.vector(drawsMaxwell$sigma),
         "gamma",
         para = list(
           shape = maxwellPD$sigma.shape,
           rate = maxwellPD$sigma.rate))
plotdist(as.vector(drawsMaxwell$alpha),
         "gamma",
         para = list(
           shape = maxwellPD$alpha.shape,
           rate = maxwellPD$alpha.rate))
plotdist(as.vector(drawsMaxwell$r0),
         "gamma",
         para = list(
           shape = maxwellPD$r0.shape,
           rate = maxwellPD$r0.rate))
plotdist(as.vector(drawsMaxwell$T),
         "gamma",
         para = list(
           shape = maxwellPD$T.shape,
           rate = maxwellPD$T.rate))

launch_shinystan(fitGibbs, rstudio = getOption("shinystan.rstudio"))
launch_shinystan(fitMaxwell, rstudio = getOption("shinystan.rstudio"))
