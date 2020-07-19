years <- names(totalIncome)
years <- years[-1]
y_maxwell <- list()
y_gibbs <- y_maxwell
for (i in 1:length(years)) {
  if (as.numeric(years[i]) < 2000) {
    bins <- 18
  } else {
    bins <- 22
  }
  y_hat <- vector(mode = "numeric", length = bins)
  for (j in 1:bins) {
    if (j < bins){
      y_hat[j] <- maxwell_int(limits = c(totalIncome$Min[j],
                                         totalIncome$Min[j + 1]),
                              theta = as.numeric(regMaxwell[i, 1:3]), 
                              x_r = double()) / 
        regMaxwell[i, 4] * sum(totalIncome[1:bins, i + 1])
    } else {
      y_hat[j] <- maxwell_int(limits = c(totalIncome$Min[j], Inf),
                              theta = as.numeric(regMaxwell[i, 1:3]), 
                              x_r = double()) / 
        regMaxwell[i, 4] * sum(totalIncome[1:bins, i + 1])
    }
  }
  y_maxwell[[toString(years[i])]] <- y_hat
  for (j in 1:bins) {
    if (j < bins){
      y_hat[j] <- gibbs_int(limits = c(totalIncome$Min[j],
                                       totalIncome$Min[j + 1]),
                            theta = as.numeric(regGibbs[i, 1:3]), 
                            x_r = double()) / 
        regGibbs[i, 4] * sum(totalIncome[1:bins, i + 1])
    } else {
      y_hat[j] <- gibbs_int(limits = c(totalIncome$Min[j], Inf),
                            theta = as.numeric(regGibbs[i, 1:3]), 
                            x_r = double()) / 
        regGibbs[i, 4] * sum(totalIncome[1:bins, i + 1])
    }
  }
  y_gibbs[[toString(years[i])]] <- y_hat
}
rm( i, j, y_hat, bins)

i = 9z
yrName = toString(years[i])
plot(totalIncome$Min[1:length(y_maxwell[[yrName]])], 
     totalIncome[[yrName]][1:length(y_maxwell[[yrName]])], 
     ylim = range(c(totalIncome[[yrName]][1:length(y_maxwell[[yrName]])], 
                    y_maxwell[[yrName]], 
                    y_gibbs[[yrName]])),
     main = paste('US Income Distribution for', yrName),
     xlab = 'Income [$ Nominal]',
     ylab = '# of People',
     log = 'xy')
lines(totalIncome$Min[1:length(y_maxwell[[yrName]])], 
      y_maxwell[[yrName]], col = "red")
lines(totalIncome$Min[1:length(y_maxwell[[yrName]])], 
      y_gibbs[[yrName]], col = "blue") 
