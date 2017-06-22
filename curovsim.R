# Load the original list of esets from curatedOvarianData, e.g.
# load("061417_esets.Rda")

# Load signatures.Rda
# load(signatures.Rda)

predictions <- vector("list", 14)

for(i in 1:length(list.sig.funs)){
	predictions[[i]] <- lapply(esets, function(x){
		tryCatch(predict(list.sig.funs[[i]], newdata = x, type = "lp")@lp, 
				error = function(e){rep(NA, ncol(x))})
			})
}

tmp <- vector("list", 14)

cindex <- function(eset, pred, nitr = 100){
	eval.data <- cbind(eset$y[,1], eset$y[, 2], pred)
	##Calculate Cox coefficient of scaled risk score, and interval:
	tryCatch(Inf.Cval(mydata=eval.data, tau=1460, itr = nitr)$Dhat,
			error = function(e){NA})
}

# Rows are prediction fcns, columns are datasets
cmat <- matrix(NA, 14, 16)

for(i in 1:length(predictions)){
	for(j in 1:length(esets)){
		cmat[i,j] <- cindex(esets[[j]], predictions[[i]][[j]])
	}
}

# ignore dataset 12 due to missingness and inability to compute c-index

tesets <- esets
tesets[[12]] <- NULL

tcmat <- cmat[,-12]

tpredictions <- predictions

# Turn predictions into ranks
# Drop dataset 12 predictions from each predictor

for(i in 1:length(tpredictions)){
	tpredictions[[i]] <- lapply(tpredictions[[i]], rank, na.last = "keep")
	tpredictions[[i]][[12]] <- NULL
}

# Subset to the ten datasets that were used in Waldron et. al.

esetnames <- c("Bentnik", "Crijns", NA, "Yoshihara10", "Mok", "Konstantinopolous",
		    NA, "Bonome", NA, "Yoshihara12", NA, "Tothill", "Dressman", "TCGA_RNA", "TCGA")

# which eset corresponds to which signature
eset2sig <- list(12, 3, NA, 6, 5, 8, NA, c(1,2), NA, 11, NA, NA, NA, NA, 10, 10)

# sizes of each dataset - culled from Waldron et. al. + original articles
sigsizes <- c(185, 185, 98, 80, 53, 43, 119, 42, 35, 413, 91, 127, 285, 511)

set.seed(32804)

nruns <- 250
outlist <- vector("list", nruns)

for(i in 1:nruns){
	# Datasets to hold out
	holdout <- sample(1:length(tesets), 3) # Keep 12, hold out 3
	
	# Signatures to remove from use
	rmsigs <- na.omit(unlist(eset2sigs[holdout]))

	if(!is.null(rmsigs)){
		tctmp <- tcmat[-rmsigs, -holdout]
		n_wt <- sigsizes[-rmsigs]
	} else {
		tctmp <- tcmat
		n_wt <- sigsizes
	}
	
	means <- apply(tctmp, 1, mean, na.rm = T)	

	cs_wt <- absnorm(means)

	n_wt <- absnorm(n_wt)
	
	unwt_c <- nwt_c <- wt_c <- rep(NA, ncol(tcmat))

	for(j in holdout){
		#Get all predictions
		allpreds <- do.call(rbind, lapply(tpredictions, "[[", j))

		if(!is.null(rmsigs)){
			prdtmp <- allpreds[-rmsigs,]
		} else {
			prdtmp <- allpreds
		}

		# unweighted predictions 
		unwt_preds <- colMeans(prdtmp, na.rm = T)

		# SS weighted predictions
		nwt_preds <- apply(prdtmp, 2, function(x){sum(x * n_wt, na.rm = T)})
		
		#weighted predictions
		wt_preds <- apply(prdtmp, 2, function(x){sum(x * cs_wt, na.rm = T)})
	
		unwt_c[j] <- cindex(tesets[[j]], unwt_preds, 10)
		nwt_c[j] <- cindex(tesets[[j]], nwt_preds, 10)
		wt_c[j] <- cindex(tesets[[j]], wt_preds, 10)
	}

	tc <- tcmat
	tc[rmsigs,] <- NA # Don't include the held-outs signatures 
	tc[,-holdout] <- NA # Report performance on the holdout datasets
	outlist[[i]] <-  scale(rbind(tc, unwt_c, nwt_c, wt_c))
}

# Generate boxplot (Figure 2)

tmp <- sapply(outlist, function(x){x}, simplify = "array")
boxout <- matrix(NA, nruns*ncol(outlist[[1]]), nrow(outlist[[i]]))
for(i in 1:ncol(boxout)){
	boxout[,i] <- as.vector(tmp[i,,])
}

labels <- c(signames, "Avg.", "SS", "CS")
labels[1] <- "Bon._263"
labels[2] <- "Bon._572"
labels[8] <- "Konst.10"
par(cex.axis = 1, cex.lab = 1, font.lab = 2)
bp <- boxplot(boxout - boxout[,ncol(boxout)], xaxt = "n", names = labels, ylim = c(-4, 4.5), las = 2, axis.cex = 1.5, ylab = "Difference in Validation C-indices (normalized)")
tm <- apply(boxout - boxout[,ncol(boxout)], 2, max, na.rm = T)
ats <- apply(cbind(bp$stats[5,], tm), 1, max) + .9
ats[c(4,6,12)] <- ats[c(4,6,12)] + 0.2
axis(1, at=seq_along(labels), labels=FALSE)
text(x=seq_along(labels), y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]),
labels=labels, srt=45, adj=1, xpd=TRUE, cex = 1, font = 2)
pcts <- apply(boxout, 2, function(x){sum(boxout[,ncol(boxout)] > x, na.rm = T)/sum(!is.na(x))})
ats <- rep(4.3,17)
ats[1:17 %% 2 == 1] <- ats[1:17 %% 2 == 1] - 0.5
text(x = 1:17, y = ats, labels = round(pcts, 2), cex = 1, font = 2)
abline(h=0, lty = "dashed")