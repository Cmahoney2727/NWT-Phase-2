#
#
#
##################################################################################
rm(list = objects())

require(raster)
require(rgdal)
require(rgeos)
require(hydroGOF)

#Shapefile
P <- "E:/NWT/NWT_Phase2/GLASData/ShpFiles/Phase2/Lv4_Filtered"
Shp <- "GLAS_P2_L23A_Lv4_Predictors"

	S <- readOGR(P, Shp)
	Sb <- gBuffer(S, byid = TRUE, width = (S@data$AvgDmtr / 2), quadsegs = 10)

#Raster stack
Rpth <- "E:/NWT/NWT_Phase2/Imputation/Mosaics"
Rs <- list.files(Rpth, ".tif$")
	
	S <- stack()
for(i in 1:length(Rs))
	S <- stack(S, raster(paste0(Rpth, "/", Rs[i])))
	
	#Assessment of regional predictions by GLAS data
		system.time(O <- extract(S, Sb, fun = mean, weights = TRUE, sp = TRUE))
	tmp <- writeOGR(O, "E:/NWT/NWT_Phase2/Imputation/Assessment", "GLAS_Assessment_kNN_RF_Predictions", driver = "ESRI Shapefile")
	
#Remove NAs
	data <- O@data[!is.na(O@data$RF_CC), ]
		Rmvs <- which(data[, "Gcc"] < 5)
	if(length(Rmvs) > 0)
		data <- data[-Rmvs, ]
		
		
	#####Stand heigt
	#RF
	Reg <- lm(data[, "Gsht"] ~ data[, "RF_SHt"])
		Rsq <- summary(Reg)$adj.r.squared
		Rmse <- summary(Reg)$sigma
	#kNN
	Reg <- lm(data[, "Gsht"] ~ data[, "kNN_SHt"])
		Rsq <- summary(Reg)$adj.r.squared
		Rmse <- summary(Reg)$sigma
	#####Crown closure
	#RF
	Reg <- lm(data[, "Gcc"] ~ data[, "RF_CC"])
		Rsq <- summary(Reg)$adj.r.squared
		Rmse <- summary(Reg)$sigma
	#kNN
	Reg <- lm(data[, "Gcc"] ~ data[, "kNN_CC"])
		Rsq <- summary(Reg)$adj.r.squared
		Rmse <- summary(Reg)$sigma
		
#####Plots
	#####Stand height
	#RF
	par(mar = c(3.5, 3.5, 1, 1))
	hist(data[, "Gsht"], seq(0, 35, 2), xlab = "", ylab = "", main = "", ylim = c(0, 2500))
	hist(data[, "RF_SHt"], seq(0, 35, 2), add = TRUE, lty = 2, lwd = 2, border = "grey60", xlab = "", ylab = "", main = "")
	title(xlab = "Stand height (m)", ylab = "Frequency", line = 2.5)
	box(bty = "l")
	#kNN
	par(mar = c(3.5, 3.5, 1, 1))
	hist(data[, "Gsht"], seq(0, 35, 2), xlab = "", ylab = "", main = "", ylim = c(0, 2500))
	hist(data[, "kNN_SHt"], seq(0, 35, 2), add = TRUE, lty = 2, lwd = 2, border = "grey60", xlab = "", ylab = "", main = "")
	title(xlab = "Crown closure (%)", ylab = "Frequency", line = 2.5)
	box(bty = "l")
	#####Crown closure
	#RF
	par(mar = c(3.5, 3.5, 1, 1))
	hist(data[, "Gcc"], seq(10, 85, 3), xlab = "", ylab = "", main = "", ylim = c(0, 1200))
	hist(data[, "RF_CC"], seq(10, 85, 3), add = TRUE, lty = 2, lwd = 2, border = "grey60", xlab = "", ylab = "", main = "")
	title(xlab = "Stand height (m)", ylab = "Frequency", line = 2.5)
	box(bty = "l")
	#kNN
	par(mar = c(3.5, 3.5, 1, 1))
	hist(data[, "Gcc"], seq(10, 85, 3), xlab = "", ylab = "", main = "", ylim = c(0, 1200))
	hist(data[, "kNN_CC"], seq(10, 85, 3), add = TRUE, lty = 2, lwd = 2, border = "grey60", xlab = "", ylab = "", main = "")
	title(xlab = "Crown closure (%)", ylab = "Frequency", line = 2.5)
	box(bty = "l")
	