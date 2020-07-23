library(EnvStats)

years <- names(binReturns)
y_maxwell <- list()
y_gibbs <- y_maxwell
binVect <- names(binMin)

for (i in 1:length(years)) {
  
  yrName = toString(years[i])
  
  for (j in length(binVect):1){
    if ( years[i] >= binVect[j]) {
      bins <- length(binReturns[[yrName]])
      x <- binMin[[j]]
    }
  }
  y_hat <- vector(mode = "numeric", length = bins)
  for (j in 1:bins) {
    if (j < bins){
      y_hat[j] <- maxwell_int(limits = c(x[j], x[j + 1]),
                              theta = as.numeric(regMaxwell[i, 1:3]), 
                              x_r = double()) / 
        regMaxwell[i, 4] * sum(binReturns[[yrName]])
    } else {
      y_hat[j] <- maxwell_int(limits = c(x[j], Inf),
                              theta = as.numeric(regMaxwell[i, 1:3]), 
                              x_r = double()) / 
        regMaxwell[i, 4] * sum(binReturns[[yrName]])
    }
  }
  y_maxwell[[toString(years[i])]] <- y_hat
  for (j in 1:bins) {
    if (j < bins){
      y_hat[j] <- gibbs_int(limits = c(x[j], x[j + 1]),
                            theta = as.numeric(regGibbs[i, 1:3]), 
                            x_r = double()) / 
        regGibbs[i, 4] * sum(binReturns[[yrName]])
    } else {
      y_hat[j] <- gibbs_int(limits = c(x[j], Inf),
                            theta = as.numeric(regGibbs[i, 1:3]), 
                            x_r = double()) / 
        regGibbs[i, 4] * sum(binReturns[[yrName]])
    }
  }
  y_gibbs[[toString(years[i])]] <- y_hat
}
rm( i, j, y_hat, bins)

i <- 25
for (j in length(binVect):1){
  if ( years[i] >= binVect[j]) {
    x <- binMin[[j]]
  }
}
yrName = toString(years[i])
plot(x, 
     binReturns[[yrName]], 
     ylim = range(c(binReturns[[yrName]], 
                    y_maxwell[[yrName]], 
                    y_gibbs[[yrName]])),
     main = paste('US Income Distribution for', yrName),
     xlab = 'Income [$ Nominal]',
     ylab = '# of People',
     log = 'xy')
lines(x, y_maxwell[[yrName]], col = "red")
lines(x, y_gibbs[[yrName]], col = "blue") 
legend("bottomleft",
       inset = 0.02,
       c("IRS Data", "Mod Maxwell", "Mod Gibbs"),
       cex=.8,
       col=c("black", "red", "blue"),
       pch=c(1, NA, NA),
       lty=c(NA, 1, 1))

# find the cross over point noted by Banerjee 2010
# With the regression done using 
fn <- function(x) crossprod(
  dpareto(x, regMaxwell$r0[i], shape = regMaxwell$alpha[i]) - 
  dgamma(x, 3/2, scale = regMaxwell$T[i]))
# solnMaxwell <- optim(regMaxwell$T[i], fn)
solnMaxwell <- optimize(fn, lower = regMaxwell$T[i], upper = 10 * regMaxwell$T[i])
fn <- function(x) crossprod(
  dpareto(x, regGibbs$r0[i], shape = regGibbs$alpha[i]) - 
  dgamma(x, 1, scale = regGibbs$T[i]))
# solnGibbs <- optim(regGibbs$T[i], fn)
solnGibbs <- optimize(fn, lower = regGibbs$T[i], upper = 10 * regGibbs$T[i])

pop <- maxwell_int(limits = c(regMaxwell[i,2], Inf),
                   theta = as.numeric(regMaxwell[i, 1:3]), 
                   x_r = double()) / regMaxwell[i, 4]
f <- 1 - regMaxwell[i,1] * (1 - pop) * sum(binReturns[[yrName]]) /
  sum(binIncome[[yrName]])/1000

# Plot income distributions
fGibbs <- function(x) gibbs_dist(x, NaN, as.numeric(regGibbs[i, 1:3]),
                                 double(), integer()) / regGibbs[i, 4]
fMaxwell <- function(x) maxwell_dist(x, NaN, as.numeric(regMaxwell[i, 1:3]),
                                     double(), integer()) / regMaxwell[i, 4]

plot(Vectorize(fMaxwell), 1, 1e8, 
     main = paste('US Income Distribution for', yrName),
     xlab = 'Income [$ Nominal]', ylab = 'Density',
     log = 'xy', col = "red")
par(new=TRUE)
plot(Vectorize(fGibbs), 1, 1e8,
     main = paste('US Income Distribution for', yrName),
     xlab = 'Income [$ Nominal]', ylab = 'Density',
     log = 'xy', col = "blue")
legend("topright",
       inset = 0.02,
       c("Mod Maxwell", "Mod Gibbs"),
       cex=.8,
       col=c("red", "blue"),
       pch=c(NA,NA),
       lty=c(1,1))