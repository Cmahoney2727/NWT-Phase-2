#	source("E:/NWT_Phase2/Scripts/R/Post_Filter_GLAS_ALS_Comparisons/ALS-GLAS_Comparisons.R")

P <- "X:/NWT/Phase2/ALS_GLAS_Gradient/Footprint_Data/03_Slope/"
Files <- list.files(P, ".csv")

	d <- NULL
for(i in 1:length(Files))
	d <- rbind(d, read.csv(paste0(P, Files[i]), header = TRUE, stringsAsFactors = FALSE))
	
	#Define ALS stand height
	Y <-  0.96 * d[, "X095"] + 0.53
	#Candidates to model Y
	Xv <- data.frame(Y, d[, c("P085", "P090", "P095", "P100")], stringsAsFactors = FALSE)
		Rmvs <- which(Xv[, "P100"] > 10)
	if(length(Rmvs) > 0)
		Xv <- Xv[-Rmvs, ]
	
	par(mar = c(3.5, 3.5, 1, 1), mfrow = c(2, 2))
		Stats <- Pstats <- matrix(NA, ncol = 3, nrow = ncol(Xv) - 1, dimnames = list(NULL, c("Rsq", "Rmse", "MAD")))
		Xp <- data.frame(matrix(NA, ncol = 5, nrow = nrow(Xv), dimnames = list(NULL, c("Y", "GH085", "GH090", "GH095", "GH100"))), stringsAsFactors = FALSE)
		Xp[, "Y"] <- Xv[, "Y"]
	for(i in 2:ncol(Xv)){
		tmp <- glm(paste0("Y ~ ", names(Xv)[i]), data = Xv)
		#Compute statistics
		Rsq <- round(summary(tmp)$adj.r.squared, 2)
		Rmse <- round(summary(tmp)$sigma, 2)
		Mad <- round(mean(abs(Xv[, "Y"] - Xv[, names(Xv)[i]])), 2)
		Stats[i - 1, ] <- c(Rsq, Rmse, Mad)
		#Plot data
		#plot(Xv[, c(names(Xv)[i], "Y")], xlim = c(0, max(Xv[, c(names(Xv)[i], "Y")])), ylim = c(0, max(Xv[, c(names(Xv)[i], "Y")])), xlab = "", ylab = "")
		#abline(tmp, lty = 2)
		#title(xlab = paste("GLAS", names(Xv)[i], "(m)"), ylab = "ALS stand height (m)", line = 2.5)
		
		Xp[, i] <- coef(tmp) * Xv[, i] + 0
		
		plot(Xp[, c(names(Xp)[i], "Y")], xlim = c(0, max(Xp[, c(names(Xp)[i], "Y")])), ylim = c(0, max(Xp[, c(names(Xp)[i], "Y")])), xlab = "", ylab = "")
		abline(0, 1)
		title(xlab = paste("GLAS", names(Xp)[i], "(m)"), ylab = "ALS stand height (m)", line = 2.5)
		
		Reg <- lm(paste0("Y ~ ", names(Xp)[i], " + 0"), data = Xp)
		Rsq <- round(summary(Reg)$adj.r.squared, 2)
		Rmse <- round(summary(Reg)$sigma, 2)
		Mad <- round(mean(abs(Xp[, "Y"] - Xp[, names(Xp)[i]])), 2)
		Pstats[i - 1, ] <- c(Rsq, Rmse, Mad)
	}
	
	D <- data.frame(d[, c(66:71)], d[, c("P085", "P090", "P095", "P100")], stringsAsFactors = FALSE)
	pairs(D)
		
	
	
	
		Usite <- unique(d[, "site"])
		P <- 1:length(Usite)
	par(mar = c(3.5, 3.5, 1, 1))
	plot(NA, NA, xlim = c(0, 20), ylim = c(0, 20), xlab = "", ylab = "")
	for(i in 1:length(Usite)){
		iii <- which(d[, "site"] == Usite[i])
		points(d[iii, "P090"], Y[iii], pch = P)
	}
		abline(0, 1)
		abline(lm(Y ~ d[, "P090"] + 0), lty = 2)
	title(xlab = paste0("GLAS ", "P090", " (m)"), ylab = "ALS P095 (m)", line = 2.5)
	
	legend(x = 0, y = 20, Usite, pch = 1:length(Usite))
	
	