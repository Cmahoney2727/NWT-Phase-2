# Append which weather station is closest to each GLAS footprint for snow data purposes
#
#	source("X:/Rscripts/Append_WeatherStn_SnowRecords_to_GLAS.R")
#############################################################################################################
rm(list = objects())

require(rgdal)
require(foreign)
require(rgeos)

#read weather station & snow threshold data
Q <- readOGR("X:/NWT/Phase2/ALS_GLAS_Gradient/Snow_Data/Shapefiles", "WeatherStnLocations", stringsAsFactors = FALSE)
FID <- as.numeric(rownames(Q@data)) - 1
Q@data <- data.frame(FID, Q@data)

#read GLAS data
P <- "X:/NWT/Phase2/ALS_GLAS_Gradient/Footprint_Data/00_Partial_Intersection_Backup"
Files <- list.files(P, ".csv")

	dir.create("X:/NWT/Phase2/ALS_GLAS_Gradient/Footprint_Data/01_Partial_Intersection/Shps")
for(i in 1:length(Files)){
	print(Files[i])

	d <- read.csv(paste0(P, "/", Files[i]), header = TRUE, stringsAsFactors = FALSE)
	Sp <- SpatialPoints(d[, c("lon", "lat")], Q@proj4string)
	
	tmp <- gDistance(Q, Sp, byid = TRUE)
	minD <- apply(tmp, 1, function(x) order(x, decreasing = FALSE)[1])
	Out <- cbind(d, Q@data[minD, 7:10])
	colnames(Out) <- c(colnames(d), colnames(Q@data)[7:10])
	
	write.csv(Out, paste0("X:/NWT/Phase2/ALS_GLAS_Gradient/Footprint_Data/01_Partial_Intersection/", Files[i]), quote = FALSE, row.names = FALSE)

	sp <- SpatialPointsDataFrame(coords = Out[, c("easting", "northing")], data = Out, proj4string = CRS(paste0("+proj=utm +zone=", Out[1, "zone"], " +datum=NAD83")))
	writeOGR(sp, paste0("X:/NWT/Phase2/ALS_GLAS_Gradient/Footprint_Data/01_Partial_Intersection/Shps"), substr(Files[i], 1, nchar(Files[i]) - 4), driver = "ESRI Shapefile", overwrite_layer = TRUE)
}

	