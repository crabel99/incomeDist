# Load Libraries
suppressPackageStartupMessages({
  library(readxl)
  library(rstan)
  library(BH)
  library(RcppEigen)
  library(Rcpp)
  library(inline)
  library(fitdistrplus)
})

message("Compiling the Stan code modules... ", appendLF = FALSE)
suppressMessages({
  expose_stan_functions("./gibbs_model.stan")
  expose_stan_functions("./maxwell_model.stan")
})
message("COMPLETE")

mod_maxwell <- function(x, theta) {
  T <- theta[1]
  r0 <- theta[2]
  alpha <- theta[3]
  c <- c <- maxwell_int(c(0, Inf), theta, double())
  v <- ((1 - x/r0 + (x/r0)^2) / (1 + x/r0)^2)^(r0 / (T * 6)) *
    x / (1 + (x/r0)^3)^((2 + alpha)/3) *
    exp(r0 / (sqrt(3) * T) * atan((1 - 2 * x / r0) / sqrt(3)))
  v <- -log(v / c) * v / c
  if (x == 0 || x == Inf) v <- 0
  return(v)
}

# Load and import data
load("~/Dropbox/src/statEcon/incomeDist/model_output_1.1 TI.RData")

# This data is imported from EIA Annual Energy Review Table 1.5
# The Annual Energy Review was discontinued in 2012
# https://www.eia.gov/totalenergy/data/annual/showtext.php?t=ptb0105
# The columns
# year - The year of the report
# total - Total primary energy consumption [Billion Btu]
# percap - Per Capita energy consumption [Million Btu]
# paid - Energy Expenditures [Millions Nominal Dollars]
# numbered - status column of the data preceding (R) indicates Revised data
dataNames <- c("year","total","1","percap","2","paid","3")
annualReview <-
  read_excel("Average Wage Index.xlsx",
             sheet = "Energy",
             range = "A10:G71",
             col_names = dataNames)
annualReview$paid <- suppressWarnings(as.numeric(annualReview$paid)) * 1E6
annualReview$total <- annualReview$total * 1055.05585262 # convert Billion Btu to GJ

CPI <- read_excel("CPI-U 1970-2010.xlsx", range = "A12:D620")
GDPDEF <- read_excel("GDPDEF.xls", range = "A11:B305")

EPI <- annualReview$total/annualReview$paid
plot(
  annualReview$year[!is.na(EPI)],
  EPI[!is.na(EPI)],
  ylim = range(c(EPI[!is.na(EPI)])),
  main = 'Energy Price Index 1970 to 2010',
  xlab = 'Year',
  ylab = 'EPI [GJ/$]',
  log = 'y'
)

indicies <-
  data.frame(
    year = numeric(length(1970:2010)),
    CPI = numeric(length(1970:2010)),
    GDP = numeric(length(1970:2010)),
    EPI = numeric(length(1970:2010))
  )
indiciesNorm <- indicies
for (i in 1970:2010) {
  indicies$year[i - 1969] <- i
  indicies$CPI[i - 1969] <-
    mean(CPI$`Observation Value`[grepl(as.character(i), CPI$Year, fixed = TRUE)])
  indicies$GDP[i - 1969] <-
    mean(GDPDEF$GDPDEF[grepl(as.character(i), GDPDEF$observation_date, fixed = TRUE)])
  indicies$EPI[i - 1969] <- 1 / EPI[annualReview$year == i]
  indiciesNorm[i - 1969, 1] <- i
  indiciesNorm[i - 1969, c(2, 3, 4)] <-
    indicies[i - 1969, c(2, 3, 4)] / indicies[1, c(2, 3, 4)]
}

plot(
  indiciesNorm$year,
  indiciesNorm$CPI,
  ylim = range(c(
    indiciesNorm$CPI,
    indiciesNorm$GDP,
    indiciesNorm$EPI
  )),
  main = 'Comparison of Various Deflators 1970-2010',
  xlab = 'Year',
  ylab = 'Value Relative to 1970',
  col = 'yellowgreen',
  type = 'l'
)
lines(indiciesNorm$year, indiciesNorm$GDP, col = "red")
lines(indiciesNorm$year, indiciesNorm$EPI, col = "blue")
legend(
  "topleft",
  inset = 0.02,
  c("CPI", "GDP", "EPI"),
  cex = 0.8,
  col = c("yellowgreen", "red", "blue"),
  pch = c(NA, NA, NA),
  lty = c(1, 1, 1)
)

list.reg <- as.numeric(row.names(regGibbs)) <= max(indicies$year)
list.index <- indiciesNorm$year >= min(as.numeric(row.names(regGibbs)))

Temp <- data.frame(year = indicies$year[list.index],
                   nominal = regMaxwell$T[list.reg],
                   CPI = regMaxwell$T[list.reg]/indiciesNorm$CPI[list.index],
                   GDP = regMaxwell$T[list.reg]/indiciesNorm$GDP[list.index],
                   EPI = regMaxwell$T[list.reg]/indiciesNorm$EPI[list.index])
r0 <- data.frame(year = indicies$year[list.index],
                 nominal = regMaxwell$r0[list.reg],
                 CPI = regMaxwell$r0[list.reg]/indiciesNorm$CPI[list.index],
                 GDP = regMaxwell$r0[list.reg]/indiciesNorm$GDP[list.index],
                 EPI = regMaxwell$r0[list.reg]/indiciesNorm$EPI[list.index])
s <- r0
for (i in 1:length(Temp$year)){
  s$nominal[i] <- maxwell_ent(c(Temp$nominal[i], r0$nominal[i],regMaxwell$alpha[i]), double())
  s$CPI[i] <- maxwell_ent(c(Temp$CPI[i], r0$CPI[i],regMaxwell$alpha[i]), double())
  s$GDP[i] <- maxwell_ent(c(Temp$GDP[i], r0$GDP[i],regMaxwell$alpha[i]), double())
  s$EPI[i] <- maxwell_ent(c(Temp$EPI[i], r0$EPI[i],regMaxwell$alpha[i]), double())
}
plot(
  Temp$year,
  Temp$nominal,
  ylim = range(c(
    Temp$nominal,
    Temp$CPI * Temp$nominal[1] / Temp$CPI[1],
    Temp$GDP * Temp$nominal[1] / Temp$GDP[1],
    Temp$EPI * Temp$nominal[1] / Temp$EPI[1]
  )),
  main = 'Economic Temperature Using Different Deflators 1996-2010',
  xlab = 'Year',
  ylab = '1996 Dollars',
  type = 'l'
)
lines(Temp$year, Temp$CPI * Temp$nominal[1] / Temp$CPI[1], col = "brown")
lines(Temp$year, Temp$GDP * Temp$nominal[1] / Temp$GDP[1] , col = "red")
lines(Temp$year, Temp$EPI * Temp$nominal[1] / Temp$EPI[1], col = "blue")
legend(
  "topleft",
  inset = 0.02,
  c("Nominal", "CPI", "GDP", "EPI"),
  cex = 0.8,
  col = c("black", "green", "red", "blue"),
  pch = c(NA, NA, NA, NA),
  lty = c(1, 1, 1, 1)
)

plot(
  s$year,
  s$nominal - s$nominal[1],
  ylim = range(c(
    s$nominal - s$nominal[1],
    s$CPI - s$CPI[1],
    s$GDP - s$GDP[1],
    s$EPI - s$EPI[1]
  )),
  main = 'Specific Economic Entropy Using Different Deflators 1996-2010',
  xlab = 'Year',
  ylab = 'Specific Entropy Referenced to 1996',
  type = 'l'
)
lines(s$year, s$CPI - s$CPI[1], col = "green")
lines(s$year, s$GDP - s$GDP[1] , col = "red")
lines(s$year, s$EPI - s$EPI[1], col = "blue")
legend(
  "topleft",
  inset = 0.02,
  c("Nominal", "CPI", "GDP", "EPI"),
  cex = 0.8,
  col = c("black", "green", "red", "blue"),
  pch = c(NA, NA, NA, NA),
  lty = c(1, 1, 1, 1)
)