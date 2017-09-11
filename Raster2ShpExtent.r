# Convert raster to polygon extent
#
#
########################################################################################
rm(list = objects())

require(rgeos)
require(rgdal)
require(raster)

P <- "E:/CFS/ALS_Deliverables/Extents"
Files <- list.files(P, ".tif$")

for(i in c(5, 4, 7, 2, 6, 3)){		#1:length(Files)
	ttt <- Sys.time()
	print(Files[i])

	R <- raster(paste0(P, "/", Files[i]))
	print("Converting...")
	system.time(r <- R > -Inf)
	system.time(pp <- rasterToPolygons(r, na.rm = TRUE, dissolve = TRUE))
	print("Saving...")
	writeOGR(pp, paste0(P, "/Shps"), substr(Files[i], 1, nchar(Files[i]) - 4), driver = "ESRI Shapefile")
	
	print(abs(ttt - Sys.time()))
}