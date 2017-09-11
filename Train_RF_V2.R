rm(list = objects())

require(raster)
require(rgdal)
require(randomForest)

RFimplementation <- function(Rdata, N, Mtry = length(predictors)){	
	Rname <- colnames(Rdata)[1]
	Pnames <- colnames(Rdata)[-1]
	set.seed(123)
	seeds <- sample(seq(1, 1e6, 1), N)
	
		Imp <- as.data.frame(matrix(NA, nrow = length(seeds) * length(Pnames), ncol = 3, dimnames = list(NULL, c("Pname", "%IncMSE", "IncNodePurity"))))
	for(s in 1:length(seeds)){
		set.seed(seeds[s])
		rf <- randomForest(x = Rdata[, Pnames], y = Rdata[, Rname], mtry = Mtry, importance = TRUE)
		
		if(s == 1){
			i1 <- s; i2 <- i1 + length(Pnames) - 1
		}else{
			i1 <- i2 + 1; i2 <- i1 + length(Pnames) - 1
		}
		Imp[i1:i2, ] <- data.frame(rownames(importance(rf)), importance(rf), stringsAsFactors = FALSE)
	}
	list(Imp = Imp, RF = rf)
}	

	Opth <- "X:/NWT/Phase1/RF_Analysis/Outputs/Trained_Models/"
	pth <- "X:/NWT/Phase1/RF_Analysis/GLAS_Layers/"
	Filters <- list.files(pth, "Level_")[-5]
	for(fff in 1:length(Filters)){
		File <- list.files(paste0(pth, Filters[fff]), ".out")
		Idata <- get(load(paste0(pth, Filters[fff], "/", File)))
		for(Att in c("Sht", "CC")){
			if(Att == "Sht")
				Name <- "ShtGLAS"
			if(Att == "CC")
				Name <- "CclGLAS"
				
			Rdata <- data.frame(Idata[, Name], Idata[, c(55:62)])
			colnames(Rdata) <- c(Name, "EOSD", "Band3", "Band4", "Band5", "CDEM", "CMI", "CTI", "SMI")
			Rdata <- round(Rdata, 4)
			
			system.time(tmp <- RFimplementation(Rdata, N = 1, Mtry = 8))
			rf <- tmp$RF
			save(rf, file = paste0(Opth, "Trained_", Filters[fff], "_", Att, ".out"), compress = TRUE)
		}
	}
	




	#set path to directory where rasters (response variable and predictors) are located
	### Note, 'pth' requires a '/' at the end ###
	pth <- "X:/NWT/Phase1/RF_Analysis/Predictor_Layers/Composite/CompositeStack.tif"
	#call code to stack rasters
		S <- stack()
	for(bands in 1:8)
		S <- stack(S, raster(pth, band = bands))
		names(S) <- c("Band3", "Band4", "Band5", "CDEM", "CMI", "CTI", "EOSD", "SMI")
	#call code to sample a number of random points from the raster stack
	### Note this takes a long time ###
	Idata <- SampleStack(S, n = 1000)
	#call code to trim out Na values and pointless cells
	Rdata <- TrimData(Idata)
	#call code to implement RF N times and group importance data
	Imp <- RFimplementation(Rdata, N = 5, mtry = (ncol(Rdata) - 1))
	
		Unames <- unique(Imp[, 1])
		OverImp <- as.data.frame(matrix(NA, nrow = length(Unames), ncol = 3, dimnames = list(NULL, c("Pname", "%IncMSE", "IncNodePurity"))))
	for(i in 1:length(Unames)){
		iii <- which(Imp[, 1] == Unames[i])
		OverImp[i, ] <- data.frame(Unames[i], mean(Imp[iii, 2]), mean(Imp[iii, 3]), stringsAsFactors = FALSE)
	}
	
	Mse <- OverImp[order(OverImp[, 2]), c(1, 2)]
	Purity <- OverImp[order(OverImp[, 3]), c(1, 3)]
	#Mean squared error importance
	par(mfrow = c(1, 2))
	par(mar = c(3.5, 6.5, 1, 1))
	plot(Mse[, 2], seq(1, 6, 1), xlab = "", ylab = "", yaxt = "n")
	axis(2, at = seq(1, 6, 1), labels = Mse[, 1], las = 1)
	title(xlab = "Inc MSE (%)", line = 2.5)
	#Node purity importance
	par(mar = c(3.5, 6.5, 1, 1))
	plot(Purity[, 2], seq(1, 6, 1), xlab = "", ylab = "", yaxt = "n")
	axis(2, at = seq(1, 6, 1), labels = Purity[, 1], las = 1)	
	title(xlab = "Inc node purity", line = 2.5)
	