library(sp)
library(rgdal)
library(splancs)
library(sf)
library(geosphere)
library(data.table)

readADSBfile <- function(fnFile) #Read in the raw ADS-B file as a data table and convert fields
{
	data <- fread(fnFile, 
								header=T, fill=T, sep = "\t")
	
	#Keep only records where there is a TSLC (time since last comm) value
	data <- data[!is.na(data$tslc)]
	
	data$TIME <- as.numeric(data$TIME)
	
	#Keep only records that have a timestamp and gps position
	data <- data[!is.na(data$TIME),]
	data <- data[!is.na(data$lat),]
	
	if(!is.numeric(data$lat)) data$lat <- as.numeric(data$lat)
	if(!is.numeric(data$lon)) data$lon <- as.numeric(data$lon)
	if(!is.numeric(data$ver_velocity)) data$ver_velocity <- as.numeric(data$ver_velocity)
	if(!is.numeric(data$hor_velocity)) data$hor_velocity <- as.numeric(data$hor_velocity)
	
	#Convert to date object, and convert the int values to the proper doubles
	data$TIME <- as.POSIXct(data$TIME,origin="1970-01-01", tz="UTC")
	data$lat <- data$lat/1e7
	data$lon <- data$lon/1e7
	data$altitude <- data$altitude/1e3
	data$heading <- data$heading/1e2
	data$hor_velocity <- data$hor_velocity/1e2
	data$ver_velocity <- data$ver_velocity/1e2
	
	#Convert the 2-byte flag field into a list of T/F values
	tempFlags <- strtoi(substr(data$validFlags,3,4), base=16)
	tempFlags <- as.logical(intToBits(tempFlags))
	tempFlags <- matrix(tempFlags,ncol = 32, byrow=T)[,1:8]
	colnames(tempFlags) <- c("valid_LATLON","valid_ALTITUDE","valid_HEADING","valid_VELOCITY","valid_CALLSIGN","valid_IDENT","SIMULATED_REPORT","valid_VERTICAL_VELOCITY")
	
	cbind(data,tempFlags)
}

mergeDailyFiles <- function(fnFiles) #Combine morning/afternoon files, combine site A/B files
{
	#Get the site names
	sites.today <- unique(gsub("(.+)-(20\\d{6}[a-z]?).TSV.gz","\\1",basename(fnFiles)))
	nSitesToday <- length(sites.today)
	
	adsbSitesData <- vector("list", length(sites.today)) 
	
	for(s in 1:length(sites.today)) #for each site, create a list object for each partial file
	{
		sitefiles.today <- grep(sites.today[s],fnFiles, value=T)
		adsbPartialData <- vector("list", length(sitefiles.today))
		for(sd in 1:length(sitefiles.today))
		{
			adsbPartialData[[sd]] <- readADSBfile(sitefiles.today[sd])
			#If there's no callsign in the datafile, create a blank field for it
			if(sum(names(adsbPartialData[[sd]]) == "callsign") == 0) adsbPartialData[[sd]]$callsign <- NA
			#Ignore all records where tlsc is > 3... these are just repeated points, not real data
			adsbPartialData[[sd]] <- adsbPartialData[[sd]][adsbPartialData[[sd]]$tslc <= 2,]
		}
		adsbSitesData[[s]] <- rbindlist(adsbPartialData)
	}
	
	if(nSitesToday > 1)
	{
		for(i in 1:nSitesToday)
		{
			#Create a field to match records from site A to site B
			#Since the loggers have different clocks, look for records
			#Where lat/lon/velocities/ICAO are the EXACT SAME
			tempMatchID <- apply(cbind(adsbSitesData[[i]]$ICAO_address,
																 adsbSitesData[[i]]$lat*1e7,
																 adsbSitesData[[i]]$lon*1e7,
																 adsbSitesData[[i]]$hor_velocity*1e2,
																 adsbSitesData[[i]]$ver_velocity*1e2),1,paste,collapse="_")
			tempMatchID <- gsub(" ","",tempMatchID)
			adsbSitesData[[i]]$matchID <- tempMatchID
		}
		
		for(i in 2:nSitesToday)
		{
			#Find the first match of exact lat/lon/velocities, and save the difference in timestamps
			#Then subtract that time diff to sync the clocks
			match.loc.master <- which(adsbSitesData[[1]]$matchID %in% adsbSitesData[[i]]$matchID)[1]
			if(!is.na(match.loc.master))
			{
				match.loc <- which(adsbSitesData[[i]]$matchID == adsbSitesData[[1]]$matchID[match.loc.master])[1]
				siteTimeDiff <- adsbSitesData[[i]][match.loc]$TIME - adsbSitesData[[1]][match.loc.master]$TIME
				adsbSitesData[[i]]$TIME <- adsbSitesData[[i]]$TIME - siteTimeDiff
			}
		}
	}
	
	fnOutput <- rbindlist(adsbSitesData, use.names=TRUE)
	
	#If there was more than one site, there will be a matchID field. Remove dupes.
	if("matchID" %in% names(fnOutput)) fnOutput <- fnOutput[!duplicated(fnOutput$matchID),]
	
	fnOutput
}

qaqcData <- function(fnData) #Make sure we're not including any weird data
{
	#Make sure there's no LAT/LON records that are invalid
	fnData <- fnData[fnData$valid_LATLON,]
	fnData <- fnData[!fnData$SIMULATED_REPORT,]
	
	fnData <- fnData[fnData$lat <= 90,]
	fnData <- fnData[fnData$lat >= -90,]
	fnData <- fnData[fnData$lon <= 180,]
	fnData <- fnData[fnData$lon >= -180,]
	
	fnData$altitude[!fnData$valid_ALTITUDE] <- -1 #Highest possible commercial flight
	fnData$altitude[is.na(fnData$altitude)] <- -1
	
	fnData
}

filterDataSpatial <- function(fnData, fnBB, fnSHP, fnCutoff, fnFilterByPark) #Keep only records that match the spatial filter
{
	tPts <- fnData[,c("lon","lat")]
	colnames(tPts) <- c("x","y")
	tPts <- as.points(tPts)
	#Rough test to see if the route is even in the bounding box
	isInBB <- tPts[,"arx"] >= fnBB['x',"min"] & 
		tPts[,"arx"] <= fnBB['x',"max"] &
		tPts[,"ary"] >= fnBB['y',"min"] &
		tPts[,"ary"] <= fnBB['y',"max"]
	
	#If there is a single point in the BB, and it never went above the upper cutoff, continue
	if(sum(isInBB) != 0 && sum(fnData$altitude[isInBB]>fnCutoff) == 0)
	{
		#If we're not filtering by the park, but just the BB
		if(!fnFilterByPark)
		{
			fnData$inStudyArea <- isInBB
			fnData
		}else{
			isInPark <- rep(FALSE,dim(fnData)[1])
			for(p in fnSHP@polygons[[1]]@Polygons)
			{
				isInPark[inout(tPts,p@coords)] <- TRUE
			}
			if(sum(isInPark) > 0 && sum(fnData$altitude[isInPark] > fnCutoff) == 0)
			{
				fnData$inStudyArea <- isInPark
				fnData
			}
		}
	}
}

getAircraftInfo <- function(fnData, fnMasterData) #Find desired fields from the FAA DB joined on ICAO/HEX
{
	fnInfo <- fnMasterData[HEXID == fnData$ICAO_address[1],,]
	tempCS <- sort(unique(fnData$callsign))[1]
	
	if(dim(fnInfo)[1]==0){
		if(is.null(fnData$callSign[1]))
		{
			fnInfo <- data.frame(t(c(NA,NA,NA,"UNKNOWN","UNKNOWN")))
		}else	fnInfo <- data.frame(t(c(NA,NA,NA,"UNKNOWN",h$callsign[1])))
		
		names(fnInfo) <- c("NNUM","MDLCODE","NAME","TYPE","CALLSIGN")
	}else{
		if(is.na(tempCS))
		{
			fnInfo <- c(fnInfo[1,1:4],"UNKNOWN")
		}else	fnInfo <- c(fnInfo[1,1:4],tempCS)
		
		fnInfo$MDLCODE <- mfgData[CODE == fnInfo$MDLCODE,,]$MODEL
		names(fnInfo) <- c("NNUM","MDLCODE","NAME","TYPE","CALLSIGN")
		#info$MDLCODE <- mfgData$MODEL[mfgData$CODE == info$MDLCODE]
	}
	fnInfo
}

convertADSB2SHP <- function(fnRecList, fnInfoList, fnIsFiltBB, fnIsFiltPark) #Convert the R object to a GeoPKG
{
	for(l in 1:length(fnRecList))
	{
		tempRec <- fnRecList[[l]]
		tempInfo <- fnInfoList[[l]]
		
		#Build the filename of the resulting output file
		#outFileNameBase <- paste(parkName,format(tempRec$TIME[1],"y%Y"),
		#												 ifelse(fnIsFiltBB,ifelse(fnIsFiltPark,"prkBnd","prkBbox"),"all"),"MP",sep="_")
		#Build the filename of the resulting output file
		outFileNameBase <- paste(parkName,
														 ifelse(fnIsFiltBB,ifelse(fnIsFiltPark,"prkBnd","prkBbox"),"all"),"MP",sep="_")
		
		outName <- paste(outDir,outFileNameBase,".gpkg",sep="")
		
		#Create a layername from the ac type
		lyrName <- gsub("\\s(.)","\\U\\1",tempInfo$TYPE,perl=T)
		
		sqlQry <- sprintf("SELECT * from \"%s\" WHERE UniqueFlt = \"%s\"", lyrName, tempRec$UniqueFlt[1])
		
		#If the Uniqueflight is not already in the GPKG, continue
		if(!exists(outName) || dim(sf::st_read(outName,query=sqlQry))[1] == 0)
		{
			testXYZ <- as.matrix(tempRec[,c("lon", "lat", "altitude")])
			
			testXYZ <- SpatialPoints(testXYZ, proj4string = CRS("+init=epsg:4326"))
			#testXYZM <- spTransform(testXYZM, CRS("+init=epsg:5070"))
			testXYZ <- as.data.frame(testXYZ)
			testXYZ[,1] <- round(testXYZ[,1],5)
			testXYZ[,2] <- round(testXYZ[,2],5)
			testXYZ[,3] <- round(testXYZ[,3],1)
			
			#Create an XYZ multipoint geometry
			outMPrec <- st_multipoint(as.matrix(testXYZ))
			#outLINErec <- st_linestring(as.matrix(testXYZ))
			
			testTemp <- st_sf(UniqueFlt = tempRec$UniqueFlt[1],
												Date = format(tempRec$TIME[1],"%Y-%m-%d"),
												StartTime = format(head(tempRec$TIME,1),"%H:%M:%S"),
												EndTime = format(tail(tempRec$TIME,1),"%H:%M:%S"),
												TimeConverted = Sys.time(),
												MedianHVelocity_ms = format(median(tempRec$hor_velocity), digits=0),
												TailNum = tempInfo$NNUM,
												CallSign = tempInfo$CALLSIGN,
												Model = tempInfo$MDLCODE,
												Type = tempInfo$TYPE,
												Name = tempInfo$NAME,
												geom = st_sfc(outMPrec), crs = 4326)
			#crs = "+proj=longlat +datum=WGS84")
			
			sf::st_write(testTemp, outName, layer=lyrName, append=T, quiet=T)
		}
	}
}

createUniqueFlightIDs <- function(fnData, fnTimeSep) #Aircraft might make multiple flights in a day, create a field for this
{			
	#Get the times where the gps points exceed the cutoff to determine if it's a new flight
	takeoffs <- NULL
	takeoffs <- which(diff(fnData$TIME) > fnTimeSep)
	fnData$KEY <- "A"
	
	#If there's more than one flight, continue
	if(length(takeoffs) > 0)
	{		
		#if(takeoffs[1] == 1)
		#{
		#	fnData <- fnData[-1,]
		#	takeoffs <- which(diff(fnData$TIME) > fnTimeSep)
		#}
		if(length(takeoffs) > 0)
		{
			takeoffs <- c(0, takeoffs)
			#Set the flight key to A-Z
			keyMods <- LETTERS[1:length(takeoffs)]
			fnData$KEY <- keyMods[length(keyMods)]
			for(l in 1:(length(takeoffs)-1))
			{
				fnData$KEY[(takeoffs[l] + 1):(takeoffs[l+1])] <- keyMods[l]
			}
		}
	}
	
	#unique flight is ICAO+YYMMDD+KEY
	paste(fnData$ICAO_address, substr(gsub('-','',fnData$TIME),3,8),fnData$KEY,sep='')
	
}

prettyTime <- function(nSecs) #Output a human-readable timestamp
{
	sprintf("%i:%02i:%02i", floor(nSecs / 3600), floor(nSecs / 60) %% 60, floor(nSecs %% 60))
	#paste(floor(nSecs / 3600), "h ", floor(nSecs / 60) %% 60, "m ", floor(nSecs %% 60), "s ", sep="")
}

if(filterByBBox) #Load the GIS filter
{
	studyArea <- readOGR(list.files("GIS\\NPATMA",pattern="*.shp$",full.names = T))
	studyArea.bb <- as.data.frame(studyArea@bbox)
}

#############################################################################

outDir <- paste("OUTPUT",format(Sys.time(),"%Y%m%d/"),sep='_')
if(!dir.exists(outDir)) dir.create(outDir, showWarnings = TRUE)

# FAA Master Data read-in ----

if(!exists('masterData'))
{
	masterData <- fread(masterDataFile, header=T, sep=',', strip.white = T, colClasses = "character", quote="")[,c(1, 3, 7, 19, 34)]
	mfgData <- fread(mfgCodeFile, header=T, sep=',', strip.white = T, colClasses = "character", quote="")[,1:3]
	colnames(mfgData)[1] <- "CODE"
	colnames(masterData) <- c("NNUM","MDLCODE","NAME","TYPE","HEXID")
	
	masterData$TYPE[masterData$TYPE=="1"] <- "Glider" 
	masterData$TYPE[masterData$TYPE=="2"] <- "Balloon"
	masterData$TYPE[masterData$TYPE=="3"] <- "Blimp/Dirigible"
	masterData$TYPE[masterData$TYPE=="4"] <- "Fixed wing single engine"
	masterData$TYPE[masterData$TYPE=="5"] <- "Fixed wing multi engine"
	masterData$TYPE[masterData$TYPE=="6"] <- "Rotorcraft"
	masterData$TYPE[masterData$TYPE=="7"] <- "Weight-shift-control"
	masterData$TYPE[masterData$TYPE=="8"] <- "Powered Parachute"
	masterData$TYPE[masterData$TYPE=="9"] <- "Gyroplane"
	masterData$TYPE[masterData$TYPE=="H"] <- "Hybrid Lift"
	masterData$TYPE[masterData$TYPE=="O"] <- "Other"
}

#----

nPtsTotal <- 0
nFtsTotal <- 0

filesProcessed <- 0
filesProcessedTimes <- NULL

numRawFiles <- length(adsbFiles)

#Set up the progressbar
pbInfoPtrn <- paste("[%0",nchar(as.character(numRawFiles)),"i:%i]",sep="")

pbInfo <- paste("Processing file from:",
								"1970-01-01",
								sprintf(pbInfoPtrn,filesProcessed,numRawFiles)
								)

pb <- winProgressBar(title="ADS-B geoprocessing...",
										 label=pbInfo,
										 min=0, max=numRawFiles, initial=0, width=600)

#Get all the dates for the ingest
uniqDates <- sort(unique(gsub(".+(20\\d{6})[a-z]?.TSV.gz","\\1",adsbFiles)))

startTime <- proc.time()

for(d in uniqDates)
{
	runTime <- proc.time()
	
	adsbRecsAll <- NULL
	adsbInfoAll <- NULL
	
	files.today <- adsbFiles[grep(d,basename(adsbFiles))]

	adsbData <- mergeDailyFiles(files.today)
	
	adsbData <- qaqcData(adsbData)
	
	#Exclude points that are Xkm away from the logger... mainly to get rid of weird data points
	distToLogger <- distHaversine(adsbData[,c("lon","lat")],loggerCoords)
	adsbData <- adsbData[distToLogger <= distanceCutoff,]
	
	adsbData <- within(adsbData, rm(matchID, SIMULATED_REPORT,tslc))

	#Split the file into a list of tables for each unique ICAO address
	adsbDataByHex <- split(adsbData,adsbData$ICAO_address)
	
	for(h in adsbDataByHex)
	{
		if(filterByBBox)
		{
			h <- filterDataSpatial(h, studyArea.bb, studyArea, elevationCutoff, filterByPark)
		}else h$inStudyArea <- TRUE
		
		if(!is.null(h))
		{
			nPtsTotal <- nPtsTotal + dim(h)[1]
			
			info <- getAircraftInfo(h, masterData)
			
			##Calculate Takeoffs
			h$UniqueFlt <- createUniqueFlightIDs(h, timeNewAirTour)
			
		  h <- h[,c("UniqueFlt", "TIME", "lat", "lon",	"altitude", "hor_velocity","inStudyArea")]
			
			hByUniqueFlt <- split(h, h$UniqueFlt)
			
			if(filterByBBox)
			{
				for(u in hByUniqueFlt)
				{
					if(sum(u$inStudyArea) > 0)
					{
						adsbRecsAll[[u$UniqueFlt[1]]] <- u
						adsbInfoAll[[u$UniqueFlt[1]]] <- info
					}
				}
			}else
			{
				for(u in hByUniqueFlt)
				{
					adsbRecsAll[[u$UniqueFlt[1]]] <- u
					adsbInfoAll[[u$UniqueFlt[1]]] <- info
				}
			}
		}
	}
	
	nFtsTotal <- nFtsTotal + length(adsbRecsAll)
	if(length(adsbRecsAll) > 0) convertADSB2SHP(adsbRecsAll, adsbInfoAll, filterByBBox, filterByPark)

	#Used for progressbar updates
	filesProcessed <- filesProcessed + length(files.today)
	
	filesProcessedTimes <- c(filesProcessedTimes, (proc.time() - runTime)[3])
	elapsedTime <- (proc.time() - startTime)[3]
	remainingTime <- mean(filesProcessedTimes * (numRawFiles - filesProcessed))
	
	pbInfo <- paste("Processing",
									basename(files.today)[length(files.today)],
									"---",
									sprintf(pbInfoPtrn,filesProcessed,numRawFiles),
									"---",
									paste("(Elapsed: ", prettyTime(elapsedTime), ", Remaining: ", prettyTime(remainingTime), ")", sep="")
									)
	
	setWinProgressBar(pb, filesProcessed, label=pbInfo)
}
close(pb)

print(paste("There were",prettyNum(nPtsTotal,big.mark=','),"points collected!"))
print(paste("There were",prettyNum(nFtsTotal,big.mark=','),"unique flights!"))
