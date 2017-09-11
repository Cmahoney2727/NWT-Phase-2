# Stack kNN outputs and extract validation points
#
#	source("C:/Users/craig.mahoney/Dropbox/R/NWT/StackRaster4Validation.R")
#############################################################################################
rm(list = objects())

require(rgdal)
require(raster)

#Import reference raster
Ipth <- "X:/NWT/Phase1/RF_Analysis/Outputs/Predictions/Mosaics"
	Rfiles <- list.files(Ipth, ".tif$")
	
		Rstack <- stack()
	for(i in 1:length(Rfiles))
		Rstack <- stack(Rstack, raster(paste0(Ipth, "/", Rfiles[i])))
	
Spth <- "X:/NWT/Phase1/ALS_BT/ValidationPolygons"
Sfiles <- list.files(Spth, ".shp$")
Snames <- substr(Sfiles, 1, nchar(Sfiles) - 4)

#for(i in 2:length(Snames)){
	i <- 10
	ttt <- Sys.time()
	
		system.time(N <- ogrInfo(Spth, Snames[i])$nrows)
		system.time(d <- readOGR(Spth, Snames[i]))
		#define random number of sample points
		int <- sample.int(N, round(0.01 * N))		#round(0.01 * N)
		#subset polygons by random samples
		S <- d[int, ]
		spT <- spTransform(S, projection(Rstack))
	
		print(paste(round(0.01 * N), "samples..."))
		system.time(O <- extract(Rstack, spT, fun = mean, weights = TRUE, sp = TRUE))
		names(O)[c(17:32)] <- c("L1_CC", "L1_CC_SD", "L1_SH", "L1_SH_SD", "L2_CC", "L2_CC_SD", "L2_SH", "L2_SH_SD", "L3_CC", "L3_CC_SD", "L3_SH", "L3_SH_SD", "L4_CC", "L4_CC_SD", "L4_SH", "L4_SH_SD")
	writeOGR(O, "X:/NWT/Phase1/ALS_BT/ValidationSubset_RF", Snames[i], driver = "ESRI Shapefile")
	
	print(abs(ttt - Sys.time()))
#}