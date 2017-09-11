# Predict area wide vegetation heights
#
# Author: craig.mahoney
# Created: 25/04/2014
# Last Revised: 31/07/2014
#
#	source("X:/NWT/Phase2/Imputation/Scripts/PredictionsRF_Phase2.R")
##########################################################################
rm(list = objects())

#####PACKAGES
require(raster)
require(rgdal)
require(matlab)
require(itertools)
require(parallel)
require(foreach)
require(doSNOW)
require(randomForest)

#####FUNCTIONS
RFsdev <- function(X, Predictions){
					SDev <- rep(NA, nrow(Predictions[[X]]))
				for(i in 1:nrow(Predictions[[X]]))
					SDev[i] <- sd(Predictions[[X]][i,])
				return(SDev)	
}

#####PREAMBLE
	ttt <- Sys.time()
	
	Pth <- "X:/NWT/Phase2/Imputation/Predictor_Data/NWT_p2_Predictors.tif"
		Rstack <- stack()
	for(bands in 1:9)
		Rstack <- stack(Rstack, raster(Pth, band = bands))
		Predictors <- c("Aspect", "DEM", "SlopDEM", "EOSD201", "Trc2000", "Band3", "Band4", "Band5", "NDVI")
		names(Rstack) <- Predictors
		
		Mulitplier <- 1500
		xdim <- ceiling( (extent(Rstack)@xmax - extent(Rstack)@xmin) / (res(Rstack)[1] * Mulitplier) )# - 1
		ydim <- ceiling( (extent(Rstack)@ymax - extent(Rstack)@ymin) / (res(Rstack)[1] * Mulitplier) )# - 1
	
	Grd <- GridTopology(c(extent(Rstack)@xmin, extent(Rstack)@ymin) + (0.5 * (res(Rstack) * Mulitplier)), res(Rstack) * Mulitplier, c(xdim, ydim))
	Sg <- SpatialGrid(Grd, proj4string = CRS(projection(Rstack)))
	Sp <- as(Sg, "SpatialPolygons")
	
		image(Rstack[["DEM"]], xlim = c(extent(Sg)@xmin, extent(Sg)@xmax), ylim = c(extent(Sg)@ymin, extent(Sg)@ymax))
		plot(Sp, add = TRUE)	
	
	P <- "X:/NWT/Phase2/Imputation/Outputs/kNN/CC/SDev/"
	File <- list.files(P, ".tif")
	Pno <- sort(as.numeric(substr(File, nchar(" S_"), nchar(File) - 4)))
	

for(Att in c("SHt", "CC")){	#c("SHt", "CC")
	Opth <- paste0("X:/NWT/Phase2/Imputation/Outputs/RF/")
	dir.create(Opth)
	dir.create(paste0(Opth, "/", Att))
	dir.create(paste0(Opth, "/", Att, "/Tiles"))
	dir.create(paste0(Opth, "/", Att, "/SDev"))

		if(Att == "SHt"){
			Response <- "Gsht"
			rf <- get(load("X:/NWT/Phase2/Imputation/TrainedModels/SHt/RF_Sht_ntree250_mtry_3.out"))
		}else if(Att == "CC"){
			Response <- "Gcc"
			rf <- get(load("X:/NWT/Phase2/Imputation/TrainedModels/CC/RF_CC_ntree500_mtry2.out"))
		}
	
	for(i in Pno){	#1:length(Sp@polygons)
		P <- SpatialPolygons( list( Polygons( list( Polygon(Sp@polygons[[i]]@Polygons[[1]]@coords) ), ID = i) ), proj4string = CRS(projection(Rstack)))
		system.time(T <- crop(Rstack, P))
		system.time(Zmat <- getValues(T))
		#Zmat <- round(Zmat, 4)
		#Remove non-vegetation land cover classes
			Zmat[Zmat[, "EOSD201"] == 0, ] <- NA
			Zmat[Zmat[, "EOSD201"] == 11, ] <- NA
			Zmat[Zmat[, "EOSD201"] == 12, ] <- NA
			Zmat[Zmat[, "EOSD201"] == 20, ] <- NA
			Zmat[Zmat[, "EOSD201"] == 31, ] <- NA
			Zmat[Zmat[, "EOSD201"] == 32, ] <- NA
			Zmat[Zmat[, "EOSD201"] == 33, ] <- NA
			Zmat[Zmat[, "EOSD201"] == 40, ] <- NA
			Zmat[Zmat[, "EOSD201"] == 51, ] <- NA
			Zmat[Zmat[, "EOSD201"] == 52, ] <- NA
			Zmat[Zmat[, "EOSD201"] == 100, ] <- NA
			
				F <- function(X, Zmat){which(is.na(Zmat[, X]))}
			Rmvs <- lapply(1:9, F, Zmat)
				Rmvs <- unique(unlist(Rmvs))
		
		if(length(Rmvs) != nrow(Zmat)){
			Ncores <- 4
			cl <- makeCluster(Ncores)
			registerDoSNOW(cl)	
		system.time(predRF <- foreach(d = isplitRows(Zmat[-Rmvs, Predictors], chunks = Ncores), .combine = c, .packages = c("randomForest")) %dopar% { predict(rf, newdata = d, predict.all = TRUE) })
			print("Finished RF predictions...")
		#system.time(SD <- parLapply(cl, x = 1:Ncores, fun = kNNsdev, Predictions = predkNN, Training, Name = Att)	)
			stopCluster(cl)
		SD <- sapply(X = seq(2, Ncores * 2, 2), RFsdev, Predictions = predRF)
			SDev <- rep(NA, nrow(Zmat[-Rmvs, ]))
		for(j in 1:length(SD)){
			if(j == 1)
				i1 <- j
				i2 <- i1 + (length(SD[[j]]) - 1)
				SDev[i1:i2] <- SD[[j]]
				i1 <- i2 + 1
		}
			Variable <- rep(NA, nrow(Zmat[-Rmvs, ]))
		for(k in seq(1, length(predRF) - 1, 2)){
			if(k == 1)
				i1 <- k
				i2 <- i1 + (length(predRF[[k]]) - 1)
				Variable[i1:i2] <- round(predRF[[k]], 3)
				i1 <- i2 + 1
		}
		#Sdev	
		Smat <- rep(-1, ncol(T) * nrow(T))
			rrr <- seq(1, ncol(T) * nrow(T), 1)[-Rmvs]
		Smat[rrr] <- SDev	
		
		sp <- SpatialPointsDataFrame((expand.grid(x = seq(extent(T)@xmin, extent(T)@xmax - res(T)[1], res(T)[1]) + (0.5 * res(T)[1]), y = rev(seq(extent(T)@ymin, extent(T)@ymax - res(T)[2], res(T)[2]) + (0.5 * res(T)[1])))), data = data.frame(SDev = Smat), proj4string = CRS(projection(Rstack)))
			gridded(sp) <- TRUE
				print("Saving SDev...")
				writeGDAL(sp, fname = paste0(Opth, "/", Att, "/SDev/S_", i, ".tif"), drivername = "GTiff", setStatistics = TRUE, type = "Float32", mvFlag = -1)
		#Predictions	
		Omat <- rep(-1, ncol(T) * nrow(T))
			rrr <- seq(1, ncol(T) * nrow(T), 1)[-Rmvs]
		Omat[rrr] <- Variable
				
		sp <- SpatialPointsDataFrame((expand.grid(x = seq(extent(T)@xmin, extent(T)@xmax - res(T)[1], res(T)[1]) + (0.5 * res(T)[1]), y = rev(seq(extent(T)@ymin, extent(T)@ymax - res(T)[2], res(T)[2]) + (0.5 * res(T)[1])))), data = data.frame(kNN = Omat), proj4string = CRS(projection(Rstack)))
			gridded(sp) <- TRUE
				print("Saving predictions...")
				writeGDAL(sp, fname = paste0(Opth, "/", Att, "/Tiles/T_", i, ".tif"), drivername = "GTiff", setStatistics = TRUE, type = "Float32", mvFlag = -1)
		}
	}
	
}#End Att
