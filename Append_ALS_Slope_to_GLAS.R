# Extract GLAS data from elevation laz file, convert to DEM and analysis for most common slope, then append to L2A_L3A GLAS data files
#
#	source("E:/NWT_Phase2/Scripts/R/Phase2GLASData/Append_ALS_Slope_to_GLAS.R")
#######################################################################################################################################
rm(list = objects())

require(landsat)
require(rgdal)
require(maptools)

#set LAStools path
	LTpth <- "X:/2016Titan/LAStools/bin/"
	setwd(LTpth)
#set site name
	Site <- "Transects_2011"
#set LAZ/LAS file path
	Ipth <- paste0("E:/NWT_Phase2/ALS_LatitudeGradient/ALS_Data/", Site, "/Outputs/Ground")
#set output file location
	Opth <- paste0("E:/NWT_Phase2/ALS_LatitudeGradient/ALS_Data/", Site, "/Outputs/Non_normalized_height")
#set step size (in meters) for analysis 
	Step <- 1.0
#set number of cores for use during analysis (1 to 16)
	Cores <- 12
	
#####START#####
#make appropriate folders
dir.create(paste0("E:/NWT_Phase2/ALS_LatitudeGradient/ALS_Data/", Site, "/Outputs/Ground_GLAS_FP"))
dir.create(paste0("E:/NWT_Phase2/ALS_LatitudeGradient/ALS_Data/", Site, "/Outputs/DEM"))
#merge non-normalized height laz files as one and index
	comm <- paste0("lasmerge -i ", Ipth, "/*.laz -o ", Opth, "/", Site, "_Merged_Elevation.laz")
		system(comm)
#index merged laz file for easy polygon extraction
	comm <- paste0("lasindex -i ", Opth, "/*.laz")
		system(comm)
	
#####EXTRACT POLYOGN AS LAZ#####
Pdir <- paste0("E:/NWT_Phase2/ALS_LatitudeGradient/Footprint_Data/Shapefiles/", Site, "/")
Pfiles <- list.files(Pdir, ".shp$")	

for(p in 1:length(Pfiles)){
	print(paste(p, "of", length(Pfiles)))
#clip merged and indexed laz file as a function of GLAS polygons
	comm <- paste0("lasclip -i ", Opth, "/*.laz -poly ", Pdir, Pfiles[p], " -o E:/NWT_Phase2/ALS_LatitudeGradient/ALS_Data/", Site, "/Outputs/Ground_GLAS_FP/", substr(Pfiles[p], 1, nchar(Pfiles[p]) - 4), ".laz")
		system(comm)
}

#make DEM and output as ASCII#####
Gpth <- paste0("E:/NWT_Phase2/ALS_LatitudeGradient/ALS_Data/", Site, "/Outputs/Ground_GLAS_FP")
Gfiles <- list.files(Gpth , ".laz")

	Slope <- data.frame(matrix(NA, ncol = 2, nrow = length(Gfiles), dimnames = list(NULL, c("ID", "ALSslope"))), stringsAsFactors = FALSE)
for(g in 1:length(Gfiles)){
	print(paste(g, "of", length(Gfiles)))
	comm <- paste0("las2dem -i ", Gpth, "/", Gfiles[g], " -keep_class 2 -o E:/NWT_Phase2/ALS_LatitudeGradient/ALS_Data/", Site, "/Outputs/DEM/", substr(Gfiles[g], 1, nchar(Gfiles[g]) - 4), ".asc")
		system(comm)
	#import footprints one-by-one and define slope	
	d <- readAsciiGrid(paste0("E:/NWT_Phase2/ALS_LatitudeGradient/ALS_Data/", Site, "/Outputs/DEM/", substr(Gfiles[g], 1, nchar(Gfiles[g]) - 4), ".asc"))
		dem <- slopeasp(d)
		h <- hist(dem$slope[[1]], seq(0, 90, 0.5), plot = FALSE)
	Slope[g, ] <- data.frame(ID = substr(Gfiles[g], 1, nchar(Gfiles[g]) - 4), ALSslope = h$breaks[which(h$counts == max(h$counts))], stringsAsFactors = FALSE)
}

#pair with extracted GLAS data
Spth <- "E:/NWT_Phase2/ALS_LatitudeGradient/Footprint_Data/L2A_L3A_ALS_GLAS_Intersection/UTM_Data/"
Sfiles <- list.files(Spth, ".csv")
File <- grep(Site, Sfiles, value = TRUE)
I <- read.csv(paste0(Spth, File), header = TRUE, stringsAsFactors = FALSE)

if(any(is.na(I$ID == Slope[match(I$ID, Slope$ID), "ID"])) == TRUE){
	rmvs <- which(is.na(I$ID == Slope[match(I$ID, Slope$ID), "ID"]))
	I <- I[-rmvs, ]
}

if(all(I$ID == Slope[match(I$ID, Slope$ID), "ID"]) == TRUE){
	Odata <- data.frame(I[, -ncol(I)], ALSslope = Slope[match(I$ID, Slope$ID), "ALSslope"], Country = I[, ncol(I)], stringsAsFactors = FALSE)
	Odata[, "site"] <- rep(Site, nrow(Odata))
	write.table(Odata, paste0("E:/NWT_Phase2/ALS_LatitudeGradient/Footprint_Data/L2A_L3A_ALS_GLAS_Intersection/ALSslopeData/", Site, ".csv"), sep = ",", row.names = FALSE, quote = FALSE)
}else{
	stop("ID's do not match...process terminated!")
}


