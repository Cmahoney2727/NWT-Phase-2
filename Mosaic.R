# R mosaic
#
#
############################################################################################################
rm(list = objects())

require(raster)
require(raster)

Ipth <- "X:/2016Titan/ALS_Deliverables/0803a_YFZ_Nelson/Outputs/CanopyOutputs"
Files <- list.files(Ipth, ".tif$")

#Elevation
	Efiles <- grep("Elevation", Files, value = TRUE)
	Uelev <- paste0("Elevation_", unique(substr(Efiles, nchar(Efiles) - 6, nchar(Efiles) - 4)))
#Count
	Ucount <- c("All2m_c00", "First2m_c00", "TotalAll_c00", "TotalFirst_c00")
#Intensity
	Ifiles <- grep("int", Files, value = TRUE)
	Uint <- paste0("int_", unique(substr(Ifiles, nchar(Ifiles) - 6, nchar(Ifiles) - 4)))
#Height
	Uheight <- unique(substr(Files, nchar(Files) - 6, nchar(Files) - 4))[-1]
	
	ttt <- Sys.time()
	Usuffix <- c(Uelev, Ucount, Uint, Uheight)
	Opth <- "X:/2016Titan/Mosaicked/0803a/"
	dir.create(Opth)
for(i in 1:length(Usuffix)){
	print(paste0("(", Usuffix[i], ")  ", i, " of ", length(Usuffix)))
	Mfiles <- unique(grep(paste0("0_", Usuffix[i]), Files, value = TRUE))
	print(length(Mfiles))
		Rlist <- list()
	for(j in 1:length(Mfiles))
		Rlist[[j]] <- raster(paste0(Ipth, "/", Mfiles[j]))
	Rlist$fun <- mean
	Rlist$file <- paste0(Opth, Usuffix[i], ".tif")
	Rlist$format <- "GTiff"
	Rlist$datatype <- "FLT4S"
	Rlist$prj <- FALSE
	Rlist$options <- c("COMPRESS=LZW")
	print("Saving...")
	system.time(Mosaic <- do.call(mosaic, Rlist))
	#writeRaster(Mosaic, file = paste0("X:/2016Titan/ALS_Deliverables/0731c_T04_T05_YFS_YZF/Outputs/Mosaics/", Usuffix[i], ".tif"), format = "GTiff", datatype = "FLT4S", options = )
}
	print(abs(ttt - Sys.time()))
	