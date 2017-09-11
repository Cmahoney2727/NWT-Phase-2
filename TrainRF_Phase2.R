#	PREPARE (OPTIMAL) GLAS DATA AS OBSERVATIONS FOR RANDOM FOREST AND KNN PREDICTIONS 
#
#	source("E:/NWT/ModelFramework/Scripts/R/TrainRF.R")
#######################################################################################
rm(list = objects())

#PACKAGES
require(parallel)
require(rgdal)
require(raster)
require(randomForest)

Var <- c("SHt", "CC")
	Pth <- "X:/NWT/Phase2/Imputation/"
	Gdata <- readOGR(paste0(Pth, "GLAS_Data"), "GLAS_P2_L23A_Lv4_Predictors")@data

for(vvv in 1:length(Var)){
	#Out Directory
		dir.create(paste0(Pth, "TrainedModels/", sep = ""))
		OutDir <- paste0(Pth, "TrainedModels/", Var[vvv], "/", sep = "")
		dir.create(OutDir)
	
#####Random Forest
		Predictors <- c("Aspect", "DEM", "SlopDEM", "EOSD201", "Trc2000", "Band3", "Band4", "Band5", "NDVI")
	if(Var[vvv] == "SHt"){
		Response <- "Gsht"
			set.seed(123)
			system.time(rf <- randomForest(x = Gdata[, Predictors], y = Gdata[, Response], ntree = 250, mtry = 3))
		save(rf, file = paste(OutDir, "RF_Sht_ntree250_mtry_3.out", sep = ""), compress = TRUE)
	}	
	if(Var[vvv] == "CC"){
		Response <- "Gcc"
			set.seed(123)
			system.time(rf <- randomForest(x = Gdata[, Predictors], y = Gdata[, Response], ntree = 500, mtry = 2))
		save(rf, file = paste(OutDir, "RF_CC_ntree500_mtry2.out", sep = ""), compress = TRUE)
	}		
	
}# End vvv

