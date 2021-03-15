adsbGeneralConversionScript <- "C:/Users/XXXX/Documents/ADSB/CODE/adsbRead_v100.R"

masterDataFile <- "C:/Users/XXXX/Documents/ADSB/FAA_Releasable_Aircraft_Database/MASTER.txt"
mfgCodeFile <- "C:/Users/XXXX/Documents/ADSB/FAA_Releasable_Aircraft_Database/ACFTREF.txt"
dataDir <- "D:/DATA/ADSB/HALE/TSV"

setwd("C:/Users/XXXX/Documents/ADSB/HALE/")

parkName <- "HALE"

filterByBBox <- TRUE
filterByPark <- TRUE

linesInstead <- FALSE

elevationCutoff <- 18000 * 0.3048 #meters
distanceCutoff <- 300 * 1000 #km

loggerCoords <- c(-156.0452, 20.6617)

timeNewAirTour <- 900

adsbFilePattern <- "[a-zA-Z]+-\\d{8}[a-z]?.TSV.gz"

adsbFiles <- dir(dataDir, pattern = adsbFilePattern, recursive = T, full.names=T)

source(adsbGeneralConversionScript)
