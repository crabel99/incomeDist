library(readxl)
binReturns <- list()
binIncome <- list()
binMin <- list()
# binMin$"2013" <- c(0, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000,
#                    45000, 50000, 55000, 60000, 75000, 100000, 200000, 500000,
#                    1000000, 1500000, 2000000, 5000000, 10000000)
# binMin$"2010" <- c(0, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000,
#                  45000, 50000, 55000, 60000, 75000, 100000, 200000, 250000,
#                  500000, 1000000, 1500000, 2000000, 5000000, 10000000)
# binMin$"2000" <- binMin$"2013"
# binMin$"1993" <- c(0, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000,
#                  45000, 50000, 55000, 60000, 75000, 100000, 200000, 500000,
#                  1000000)

binMin$"2003" <- c(0, 1, 5000, 10000, 15000, 20000, 25000, 30000, 40000,
                   50000, 75000, 100000, 200000, 500000, 1000000, 1500000,
                   2000000, 5000000, 10000000)
binMin$"2000" <- c(0, 1, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
                   10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000,
                   18000, 19000, 20000, 25000, 30000, 40000, 50000, 75000,
                   100000, 200000, 500000, 1000000, 1500000, 2000000, 5000000,
                   10000000)
binMin$"1996" <- c(0, 1, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
                   10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000,
                   18000, 19000, 20000, 25000, 30000, 40000, 50000, 75000,
                   100000, 200000, 500000, 1000000)




# pathName <- "./IRS Data/Table 2.1"
pathName <- "./IRS Data/Table 1.1"

files <- list.files(path = pathName)
files <- c(files[19:length(files)],files[1:18])
# dataNames <- c("Returns", "Income")
dataNames <- c("Returns", "delete", "Income")
years <- 2017 - length(files) + 1
years <- rep(years:2017)

for (i in 1:length(files)) {
 dirName <- path.expand(paste(pathName, files[i], sep ="/"))
#   if (years[i] > 2012) {
#     tableRange <- "B11:C32"
#     incomeMin <- binMin$"00"
#   } else if(years[i] > 2009) {
#     tableRange <- "B11:C33"
#     incomeMin <- binMin$"10"
#   } else if (years[i] > 2004) {
#     tableRange <- "B11:C32"
#     incomeMin <- binMin$"00"
#   } else if (years[i] == 2004) {
#     tableRange <- "B14:C35"
#     incomeMin <- binMin$"00"
#   } else if (years[i] > 2000) {
#     tableRange <- "B9:C30"
#     incomeMin <- binMin$"00"
#   } else if (years[i] == 2000) {
#     tableRange <- "B10:C31"
#     incomeMin <- binMin$"00"
#   } else if (years[i] == 1997) {
#     tableRange <- "B10:C27"
#     incomeMin <- binMin$"93"
#   } else {
#     tableRange <- "B9:C26"
#     incomeMin <- binMin$"93"
#   }
  if (i > 7) {
    tableRange <- "B11:D29"
    incomeMin <- binMin$"2003"
  } else if (i == 7) {
    tableRange <- "B13:D47"
    incomeMin <- binMin$"2000"
  } else if (i == 5 || i == 6) {
    tableRange <- "B14:D48"
    incomeMin <- binMin$"2000"
  } else if (i == 2) {
    tableRange <- "B12:D42"
    incomeMin <- binMin$"1996"
  } else {
    tableRange <- "B13:D43"
    incomeMin <- binMin$"1996"
  }
  import <- read_xls(dirName, range = tableRange, col_names = dataNames)
  binIncome[[i]] <- import$Income
  binReturns[[i]] <- import$Returns
}
names(binIncome) <- years
names(binReturns) <- years
rm(import, dataNames, dirName, files, i, incomeMin, pathName, tableRange, years)