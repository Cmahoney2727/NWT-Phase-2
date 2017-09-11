# Script 2. Stack rasters and convert to point shapefile
#
#	source("X:/Rscripts/RasterStack2shp.R")
###############################################################################
rm(list = objects())

require(raster)
require(rgdal)

StackRasters <- function(pth, files = NULL, ...){

	Ext <- list()
for(f in 1:length(files))
	Ext[[f]] <- extent(raster(paste0(pth, files[f])))
	
	Elist <- lapply(Ext, as.matrix)
	Emat <- matrix(unlist(Elist), ncol = length(Elist), dimnames = list(c("xmin", "ymin", "xmax", "ymax"), NULL))
	#find smallest extent size for cropping
	MinExt <- c(max(Emat[1, ]), min(Emat[3, ]), max(Emat[2, ]), min(Emat[4, ]))
	
		S <- stack()
	for(i in 1:length(files))
		S <- stack(S, crop(raster(paste0(pth, files[i])), MinExt))
		
	return(S)
}
	

Pth <- "X:/2016Titan/Mosaicked/0802a/"
Files <- list.files(Pth, ".tif$")
Files <- Files[c(3:6, 8:9, 31:33, 11:30, 10, 2, 34, 56:58, 36:55, 35, 1, 7, 59:60)]

	print(Pth)

	Rstack <- StackRasters(Pth, Files)
		
	system.time(tmp <- rasterToPoints(Rstack, spatial = TRUE, progress = "text"))
	names(tmp)[c(1:4, 57:60)] <- c("Elv_avg", "Elv_max", "Elv_min", "Elv_std", "All_2m", "Fst_2m", "Tot_All", "Tot_Fst")
	writeOGR(tmp, "X:/2016Titan/Shapefiles", "0802a_YZF_NW", driver = "ESRI Shapefile")
	
	
	
	Pth <- "X:/2016Titan/Mosaicked/0802a/"
	Files <- list.files(Pth, ".tif$")
	Files <- Files[c(3:6, 8:9, 31:33, 11:30, 10, 2, 34, 56:59, 36:55, 35, 1, 7, 60:61)]
	#Files <- Files[c(34, 36, 60, 61)]
	
		Rstack <- stack()
	for(i in 1:length(Files)){
		if(Files[i] == "RGB_MB.tif"){
			Rstack <- stack(Rstack, stack(paste0(Pth, Files[i])))
		}else{	
			Rstack <- stack(Rstack, raster(paste0(Pth, Files[i])))
		}
	}
		
	system.time(tmp <- rasterToPoints(Rstack, spatial = TRUE, progress = "text"))
	names(tmp)[c(1:4, 34:36, 60:63)] <- c("Elv_avg", "Elv_max", "Elv_min", "Elv_std", "RGB_1", "RGB_2", "RGB_3", "All_2m", "Fst_2m", "Tot_All", "Tot_Fst")
	writeOGR(tmp, "X:/2016Titan/Shapefiles", "0802a_YZF_NW_RGB2", driver = "ESRI Shapefile")
		