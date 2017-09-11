# Predict area wide vegetation heights
#
# Author: craig.mahoney
# Created: 25/04/2014
# Last Revised: 31/07/2014
#
#	source("X:/Rscripts/PredictionsRF_V2.R")
##########################################################################
rm(list = objects())

require(raster)
require(rgdal)
require(matlab)
require(itertools)
require(parallel)
require(foreach)
require(doSNOW)
require(randomForest)
require(matrixStats)

#####FUNCTIONS
TileStack <- function(Rstack, Multiplier = 1500, Plot = TRUE, ...){
	#determine number of tiles to break the dataset in to	
		#Multiplier <- 1500
		xdim <- ceiling( (extent(Rstack)@xmax - extent(Rstack)@xmin) / (res(Rstack)[1] * Multiplier) )# - 1
		ydim <- ceiling( (extent(Rstack)@ymax - extent(Rstack)@ymin) / (res(Rstack)[1] * Multiplier) )# - 1
	#grid to create tiles
	Grd <- GridTopology(c(extent(Rstack)@xmin, extent(Rstack)@ymin) + (0.5 * (res(Rstack) * Multiplier)), res(Rstack) * Multiplier, c(xdim, ydim))
	Sg <- SpatialGrid(Grd, proj4string = CRS(projection(Rstack)))
	Sp <- as(Sg, "SpatialPolygons")
	#visualize tiles (not necessary)
	if(Plot == TRUE){
		image(Rstack[[4]], xlim = c(extent(Sg)@xmin, extent(Sg)@xmax), ylim = c(extent(Sg)@ymin, extent(Sg)@ymax), asp = 1)
		plot(Sp, add = TRUE)	
	}
	
	return(Sp)
}	

MosaicTiles <- function(Ipth, Ofile = paste0(getwd(), "/MosaickedTiles.tif"), ...){
	Tiles <- list.files(Ipth, ".tif$")
		Rlist <- list()
	for(j in 1:length(Tiles))
		Rlist[[j]] <- raster(paste0(Ipth, "/", Tiles[j]))
	Rlist$fun <- mean
	Rlist$file <- Ofile
	Rlist$format <- "GTiff"
	Rlist$datatype <- "FLT4S"
	Rlist$prj <- FALSE
	Rlist$options <- c("COMPRESS=LZW")
	print("Saving...")
	system.time(Mosaic <- do.call(mosaic, Rlist))
	return(Mosaic)
}

RFsdev <- function(X, Predictions){
		SDev <- rep(NA, nrow(Predictions[[X]]))
	for(i in 1:nrow(Predictions[[X]]))
		SDev[i] <- sd(Predictions[[X]][i,])
	return(SDev)	
}

RFsdev2 <- function(X, Predictions){
	SDev <- apply(Predictions[[X]], 1, sd)
	return(SDev)
}

#####PREAMBLE
	ttt <- Sys.time()
	
#directory of predictor rasters
	Ipth <- "X:/NWT/Phase1/RF_Analysis/Predictor_Layers/Composite/CompositeStack.tif"
	Mpth <- "X:/NWT/Phase1/RF_Analysis/Outputs/Predictions/Mosaics/"
	dir.create(Mpth)
	#Only use if one set has been predicted
	Tpth <- paste0("X:/NWT/Phase1/RF_Analysis/Outputs/Predictions/Level_1/SHt/Tiles")
	Tfiles <- list.files(Tpth, ".tif$")
	Tno <- sort(as.numeric(substr(Tfiles, 3, nchar(Tfiles) - 4)))
	
#####START	
	#stack rasters & rename if necessary
	Filters <- list.files("X:/NWT/Phase1/RF_Analysis/GLAS_Layers/", "Level_")[-5]
	
		Rstack <- stack()
	for(bands in 1:8)
		Rstack <- stack(Rstack, raster(Ipth, band = bands))
		names(Rstack) <- c("Band3", "Band4", "Band5", "CDEM", "CMI", "CTI", "EOSD", "SMI")
	#create tiles
		Sp <- TileStack(Rstack, Multiplier = 1000, Plot = TRUE)
		
for(fff in 2:4){	#1:length(Filters)
	for(Att in c("SHt", "CC")){	#c("SHt", "CC")
		#trained RF file path including file
			RFfile <- paste0("X:/NWT/Phase1/RF_Analysis/Outputs/Trained_Models/", "Trained_", Filters[fff], "_", Att, ".out")
		#output directory where sub-folders will be created
			Opth <- paste0("X:/NWT/Phase1/RF_Analysis/Outputs/Predictions/", Filters[fff])
			dir.create(Opth)
			dir.create(paste0(Opth, "/", Att))
			dir.create(paste0(Opth, "/", Att, "/Tiles"))
			dir.create(paste0(Opth, "/", Att, "/SDev"))
	
		#import trained RF model
			system.time(RF <- get(load(RFfile)))
			
			ttt <- Sys.time()
		for(s in Tno){	#1:length(Sp@polygons) or Tno
		#define tile to work with
			P <- SpatialPolygons( list( Polygons( list( Polygon(Sp@polygons[[s]]@Polygons[[1]]@coords) ), ID = s) ), proj4string = CRS(projection(Rstack)))
		#crop & import data values for tile
			system.time(T <- crop(Rstack, P))
			system.time(Zmat <- getValues(T))
			if(any(!is.na(Zmat))){
			#round all values to X decimals
				Zmat <- round(Zmat, 4)
			#Remove non-vegetation land cover classes
				Zmat[Zmat[, "EOSD"] == 0, ] <- NA
				Zmat[Zmat[, "EOSD"] == 11, ] <- NA
				Zmat[Zmat[, "EOSD"] == 12, ] <- NA
				Zmat[Zmat[, "EOSD"] == 20, ] <- NA
				Zmat[Zmat[, "EOSD"] == 31, ] <- NA
				Zmat[Zmat[, "EOSD"] == 32, ] <- NA
				Zmat[Zmat[, "EOSD"] == 33, ] <- NA
				Zmat[Zmat[, "EOSD"] == 40, ] <- NA
				Zmat[Zmat[, "EOSD"] == 51, ] <- NA
				Zmat[Zmat[, "EOSD"] == 52, ] <- NA
				Zmat[Zmat[, "EOSD"] == 100, ] <- NA
			#remove all NA values from predictors			
					F <- function(X, Zmat){which(is.na(Zmat[, X]))}
				Rmvs <- sapply(X = 1:8, FUN = F, Zmat)
					Rmvs <- unique(unlist(Rmvs))
				if(length(Rmvs) < 1)
					Rmvs <- -seq(1, nrow(Zmat), 1)
				if(length(which(!is.na(Zmat))) > 1){
					#cluster cores in preparation for RF predictions
						Ncores <- 6
					if(class(Zmat[-Rmvs, ]) != "numeric"){
						cl <- makeCluster(Ncores)
						registerDoSNOW(cl)	
						#make Rf predictions
							system.time(predRF <- foreach(d = isplitRows(Zmat[-Rmvs, ], chunks = Ncores), .combine = c, .packages = c("randomForest")) %dopar% { predict(RF, newdata = d, predict.all = TRUE) })
								print("Finished RF predictions...")
					#stop core cluster to relinquish memory
						stopCluster(cl)
						#prepare standard deviation classes for outputs
							system.time(SD <- sapply(X = seq(2, Ncores * 2, 2), function(X, predRF){rowSds(predRF[[X]])}, predRF))	
						#gather standard deviation values
							SDev <- rep(NA, nrow(Zmat[-Rmvs, ]))
						for(i in 1:length(SD)){
							if(i == 1)
								i1 <- i
							i2 <- i1 + (length(SD[[i]]) - 1)
							SDev[i1:i2] <- SD[[i]]
							i1 <- i2 + 1
						}
						#gather predicted (response) variable values
							Variable <- rep(NA, nrow(T) * ncol(T))
						for(i in seq(1, length(predRF) - 1, 2)){
							if(i == 1)
								i1 <- i
							i2 <- i1 + (length(predRF[[i]]) - 1)
							Variable[i1:i2] <- round(predRF[[i]], 3)
							i1 <- i2 + 1
						}
					}else{
						predRF <- predict(RF, newdata = Zmat[-Rmvs, ], predict.all = TRUE)
						SDev <- sd(predRF$individual)
						Variable <- predRF$aggregate
					}
							print("Finished preparing standard deviation classes...")
					
				#place standard deviation values in to a vector for output
					Smat <- rep(-1, ncol(T) * nrow(T))
						rrr <- seq(1, ncol(T) * nrow(T), 1)[-Rmvs]
					Smat[rrr] <- SDev	
		
					sp <- SpatialPointsDataFrame((expand.grid(x = seq(extent(T)@xmin, extent(T)@xmax - res(T)[1], res(T)[1]) + (0.5 * res(T)[1]), y = rev(seq(extent(T)@ymin, extent(T)@ymax - res(T)[2], res(T)[2]) + (0.5 * res(T)[1])))), data = data.frame(SDev = Smat), proj4string = CRS(projection(Rstack)))
						gridded(sp) <- TRUE
							print("Saving SDev...")
					writeGDAL(sp, fname = paste0(Opth, "/", Att, "/SDev/S_", s, ".tif"), drivername = "GTiff", setStatistics = TRUE, type = "Float32", mvFlag = -1)
				
				#place (response) variable values in to a vector for output
					Omat <- rep(-1, ncol(T) * nrow(T))
						rrr <- seq(1, ncol(T) * nrow(T), 1)[-Rmvs]
					Omat[rrr] <- Variable
				
					sp <- SpatialPointsDataFrame((expand.grid(x = seq(extent(T)@xmin, extent(T)@xmax - res(T)[1], res(T)[1]) + (0.5 * res(T)[1]), y = rev(seq(extent(T)@ymin, extent(T)@ymax - res(T)[2], res(T)[2]) + (0.5 * res(T)[1])))), data = data.frame(kNN = Omat), proj4string = CRS(projection(Rstack)))
						gridded(sp) <- TRUE
							print("Saving predictions...")
					writeGDAL(sp, fname = paste0(Opth, "/", Att, "/Tiles/T_", s, ".tif"), drivername = "GTiff", setStatistics = TRUE, type = "Float32", mvFlag = -1)
				}
			}
		}
			print("Mosaicking...")
			M <- MosaicTiles(Ipth = paste0(Opth, "/", Att, "/Tiles"), Ofile = paste0(Mpth, "RF_Mosaic_", Filters[fff], "_", Att, "_Predictions.tif"))
			M <- MosaicTiles(Ipth = paste0(Opth, "/", Att, "/SDev"), Ofile = paste0(Mpth, "RF_Mosaic_", Filters[fff], "_", Att, "_SDev.tif"))
		print(abs(ttt - Sys.time()))
	}
}

