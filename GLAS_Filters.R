# GLAS filtering from partial coverage
#
# source("X:/Rscripts/NWT_ALS_GLAS_Gradient/GLAS_Filters.R")
###############################################################################################
rm(list = objects())

#PACKAGES
require(rgdal)

#FUNCTIONS
CopyDir <- function(Src, Dst, Overwrite = TRUE, ...){
	Dirs <- list.dirs(Src, full.names = FALSE, recursive = TRUE)[-1]
	dir.create(paste0(Dst, "/", tail(strsplit(Src, "/")[[1]], 1)))
	for(i in 1:length(Dirs)){
		dir.create(paste0(Dst, "/", tail(strsplit(Src, "/")[[1]], 1), "/", Dirs[i]))
		file.names <- dir(paste0(Src, "/", Dirs[i])) 
		sapply(file.names, function(x) { file.copy(from = paste0(Src, x), to = paste0(Dst, x), overwrite = Overwrite) })
	}	
	
}

Opth <- "X:/NWT/Phase2/ALS_GLAS_Gradient/Footprint_Data"
dir.create(paste0(Opth, "/02_Height_Anomalies/"))
dir.create(paste0(Opth, "/03_Slope/"))
dir.create(paste0(Opth, "/04_Snow/"))
dir.create(paste0(Opth, "/05_Possible_Snow/"))
dir.create(paste0(Opth, "/06_Pulse_Energy/"))

F01pth <- "X:/NWT/Phase2/ALS_GLAS_Gradient/Footprint_Data/01_Partial_Intersection"
Files <- list.files(F01pth, ".csv")

for(i in 1:length(Files)){
#Filter level 1
	L1data <- read.csv(paste0(F01pth, "/", Files[i]), header = TRUE, stringsAsFactors = FALSE)
	L1data[, "Date"] <- as.Date(L1data[, "Date"], format = "%m/%d/%Y")
	for(X in 75:78)
		L1data[, X] <- as.Date(L1data[, X], format = "%d/%m/%Y")
	dir.create(paste0(Opth, "/01_Partial_Intersection/Shps"))
	sp <- SpatialPointsDataFrame(coords = L1data[, c("easting", "northing")], data = L1data, proj4string = CRS(paste0("+proj=utm +zone=", L1data[1, "zone"], " +datum=NAD83")))
	writeOGR(sp, paste0(Opth, "/01_Partial_Intersection/Shps"), substr(Files[i], 1, nchar(Files[i]) - 4), driver = "ESRI Shapefile", overwrite_layer = TRUE)
	
			Rmvs <- NULL
	if(any(is.na(L1data[, "X025"]))){
			Rmvs <- which(is.na(L1data[, c(66:71)]))
		if(length(Rmvs) > 0)
			L1data <- L1data[-Rmvs, ]
	}
	
#Filter level 2
	L2data <- L1data
	#Saturation Index
		SatRmvs <- which(L1data[, "Satcor"] > 2)
	if(length(SatRmvs) > 0)
		L2data <- L1data[-SatRmvs, ]
	#Waveform Extent
		ExtRmvs <- which(L2data[, "P100"] > 35)
	if(length(ExtRmvs) > 0)
		L2data <- L2data[-ExtRmvs, ]
	#First Quartile Range
		QuRmvs <- which(L2data[, "P025"] >= 50)
	if(length(QuRmvs) > 0)
		L2data <- L2data[-QuRmvs, ]
	#P100 - P095
		PRmvs <- which(L2data[, "P100"] - L2data[, "P095"] > 3)
	if(length(PRmvs) > 0)
		L2data <- L2data[-PRmvs, ]

	#Trmvs <- sum(c(length(Rmvs), length(SatRmvs), length(ExtRmvs), length(QuRmvs), length(PRmvs)), na.rm = TRUE)
	Trmvs <- sum(c(length(Rmvs), length(SatRmvs), length(ExtRmvs), length(QuRmvs)), na.rm = TRUE)
	print(paste(Trmvs, "records removed..."))
	
	dir.create(paste0(Opth, "/02_Height_Anomalies/Shps"))
	sp <- SpatialPointsDataFrame(coords = L2data[, c("easting", "northing")], data = L2data, proj4string = CRS(paste0("+proj=utm +zone=", L1data[1, "zone"], " +datum=NAD83")))
	writeOGR(sp, paste0(Opth, "/02_Height_Anomalies/Shps"), substr(Files[i], 1, nchar(Files[i]) - 4), driver = "ESRI Shapefile", overwrite_layer = TRUE)
	write.csv(L2data, paste0(Opth, "/02_Height_Anomalies/", Files[i]), row.names = FALSE, quote = FALSE)
#Filter level 3
	L3data <- L2data
	#Slope > 5
		SRmvs <- which(L3data[, "AvgALSslope"] > 5)
	if(length(SRmvs) > 0)
		L3data <- L3data[-SRmvs, ]
	#tan(slope) > 0.25 * (0.96 * X095 + 0.53)
		Dfp <- rep(NA, nrow(L3data))
	for(s in 1:nrow(L3data))
		Dfp[s] <- round(mean(c(L3data[s, "Major"], L3data[s, "Minor"])), 2)
	Sbias <- Dfp * tan(L3data[, "AvgALSslope"] * (pi / 180))
	ALSsht <- 0.96 * L3data[, "X095"] + 0.53
		BRmvs <- which((Sbias / ALSsht) > 0.25)	
	if(length(BRmvs) > 0)
		L3data <- L3data[-BRmvs, ]

	Trmvs <- sum(c(length(SRmvs), length(BRmvs)), na.rm = TRUE)
	print(paste(Trmvs, "records removed..."))
	
	if(nrow(L3data) > 0){
		dir.create(paste0(Opth, "/03_Slope/Shps"))
		sp <- SpatialPointsDataFrame(coords = L3data[, c("easting", "northing")], data = L3data, proj4string = CRS(paste0("+proj=utm +zone=", L1data[1, "zone"], " +datum=NAD83")))
		writeOGR(sp, paste0(Opth, "/03_Slope/Shps"), substr(Files[i], 1, nchar(Files[i]) - 4), driver = "ESRI Shapefile", overwrite_layer = TRUE)
		write.csv(L3data, paste0(Opth, "/03_Slope/", Files[i]), row.names = FALSE, quote = FALSE)
	#Filter level 4
		L4data <- L3data
		#Certainty of snow cover
			L2A <- which(L4data[, "Lc"] == "L2A")
				L2ARmvs <- which(L4data[L2A, "Date"] > L4data[L2A, "Snw2003"])
			if(length(L2ARmvs) > 0)
				L4data <- L4data[L2A[-L2ARmvs], ]
			L3A <- which(L4data[, "Lc"] == "L3A")
				L3ARmvs <- which(L4data[L3A, "Date"] > L4data[L3A, "Snw2004"])
			if(length(L3ARmvs) > 0)
				L4data <- L4data[L3A[-L3ARmvs], ]
	
		Trmvs <- sum(c(length(L2ARmvs), length(L3ARmvs)), na.rm = TRUE)
		print(paste(Trmvs, "records removed..."))
	
		if(nrow(L4data) > 0){
			dir.create(paste0(Opth, "/04_Snow/Shps"))
			sp <- SpatialPointsDataFrame(coords = L4data[, c("easting", "northing")], data = L4data, proj4string = CRS(paste0("+proj=utm +zone=", L1data[1, "zone"], " +datum=NAD83")))
			writeOGR(sp, paste0(Opth, "/04_Snow/Shps"), substr(Files[i], 1, nchar(Files[i]) - 4), driver = "ESRI Shapefile", overwrite_layer = TRUE)
			write.csv(L4data, paste0(Opth, "/04_Snow/", Files[i]), row.names = FALSE, quote = FALSE)
		#Filter level 5
			L5data <- L4data
			#Possibility of snow cover
				L2A <- which(L5data[, "Lc"] == "L2A")
					L2ARmvs <- which(L5data[L2A, "Date"] > L5data[L2A, "Pss2003"])
				if(length(L2ARmvs) > 0)
					L5data <- L5data[L2A[-L2ARmvs], ]
				L3A <- which(L5data[, "Lc"] == "L3A")
					L3ARmvs <- which(L5data[L3A, "Date"] > L5data[L3A, "Pss2004"])
				if(length(L3ARmvs) > 0)
					L5data <- L5data[L3A[-L3ARmvs], ]
		
			Trmvs <- sum(c(length(L2ARmvs), length(L3ARmvs)), na.rm = TRUE)
			print(paste(Trmvs, "records removed..."))
	
			if(nrow(L5data) > 0){
				dir.create(paste0(Opth, "/05_Possible_Snow/Shps"))
				sp <- SpatialPointsDataFrame(coords = L5data[, c("easting", "northing")], data = L5data, proj4string = CRS(paste0("+proj=utm +zone=", L1data[1, "zone"], " +datum=NAD83")))
				writeOGR(sp, paste0(Opth, "/05_Possible_Snow/Shps"), substr(Files[i], 1, nchar(Files[i]) - 4), driver = "ESRI Shapefile", overwrite_layer = TRUE)
				write.csv(L5data, paste0(Opth, "/05_Possible_Snow/", Files[i]), row.names = FALSE, quote = FALSE)
			#Filter level 5
				L6data <- L5data
				#GLAS pulse energy
					#Tnrg <- 0.071
					Tnrg <- 0.028
						NrgRmvs <- which(L6data[, "TxNrg"] < Tnrg)
					if(length(NrgRmvs) > 0)
						L6data <- L6data[-NrgRmvs, ]

				Trmvs <- sum(c(length(NrgRmvs)), na.rm = TRUE)
				print(paste(Trmvs, "records removed..."))
		
				if(nrow(L6data) > 0){
					dir.create(paste0(Opth, "/06_Pulse_Energy/Shps"))
					sp <- SpatialPointsDataFrame(coords = L6data[, c("easting", "northing")], data = L6data, proj4string = CRS(paste0("+proj=utm +zone=", L1data[1, "zone"	], " +datum=NAD83")))
					writeOGR(sp, paste0(Opth, "/06_Pulse_Energy/Shps"), substr(Files[i], 1, nchar(Files[i]) - 4), driver = "ESRI Shapefile", overwrite_layer = TRUE)
					write.csv(L6data, paste0(Opth, "/06_Pulse_Energy/", Files[i]), row.names = FALSE, quote = FALSE)
				}
			}
		}
	}
}


	#tmp <- CopyDir("E:/NWT_Phase2/ALS_LatitudeGradient/Footprint_Data/L2A_L3A_ALS_GLAS_Intersection/ALS_GLAS_Data_Filters", "E:/NWT_Phase2/ALS_LatitudeGradient/Footprint_Data")
