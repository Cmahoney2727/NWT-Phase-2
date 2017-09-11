# Train kNN models
#
# Author: craig.mahoney
# Created: 25/04/2014
# Last Revised: 31/07/2014
#
#	source("E:/NWT/ModelFramework/Scripts/R/TrainkNN.R")
##########################################################################
rm(list = objects())

#PACKAGES
require(parallel)
require(rgdal)
require(raster)
require(kknn)

Var <- c("SHt", "CC")
	Pth <- "X:/NWT/Phase2/Imputation/"
	Gdata <- readOGR(paste0(Pth, "GLAS_Data"), "GLAS_P2_L23A_Lv4_Predictors")@data

for(vvv in 1:length(Var)){
	#Out Directory
		dir.create(paste0(Pth, "TrainedModels/", sep = ""))
		OutDir <- paste0(Pth, "TrainedModels/", Var[vvv], "/", sep = "")
		dir.create(OutDir)	
	
#####kNN
		Predictors <- c("Aspect", "DEM", "SlopDEM", "EOSD201", "Trc2000", "Band3", "Band4", "Band5", "NDVI")
	if(Var[vvv] == "SHt"){
		Response <- "Gsht"
			set.seed(123)
			system.time(TunekNN <- train.kknn(paste0(Response, " ~ ", paste(Predictors, sep = " ", collapse = "+")), data = Gdata, kmax = 50, kernel = c("epanechnikov", "triangular", "biweight", "triweight", "cos", "inv", "gaussian", "optimal")))
		save(TunekNN, file = paste(OutDir, "kNN_Sht_Triangular_k_30.out", sep = ""), compress = TRUE)
		write.csv(Gdata[, c(Response, Predictors)], paste0(OutDir, "kNN_Sht_TrainingData.csv"), quote = FALSE, row.names = FALSE)
	}	
	if(Var[vvv] == "CC"){
		Response <- "Gcc"
			set.seed(123)
			system.time(TunekNN <- train.kknn(paste0(Response, " ~ ", paste(Predictors, sep = " ", collapse = "+")), data = Gdata, kmax = 30, kernel = c("epanechnikov", "triangular", "biweight", "triweight", "cos", "inv", "gaussian", "optimal")))
		save(TunekNN, file = paste(OutDir, "kNN_CC_Gaussian_k_30.out", sep = ""), compress = TRUE)
		write.csv(Gdata[, c(Response, Predictors)], paste0(OutDir, "kNN_CC_TrainingData.csv"), quote = FALSE, row.names = FALSE)
	}		
}

