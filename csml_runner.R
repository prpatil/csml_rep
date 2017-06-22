# From curatedOvarianData vignette
source(system.file("extdata", "patientselection.config",package="curatedOvarianData"))
sapply(ls(), function(x) if(!x %in% c("remove.samples", "duplicates")) print(get(x)))
source(system.file("extdata", "createEsetList.R", package = "curatedOvarianData"))

# Save the original eset list
#save(esets, file = "061417_esets.Rda")

# Remove esets with missing gene expression data
ridx <- which(unlist(lapply(esets, function(x){sum(is.na(exprs(x)))})) > 0)
esets <- esets[-ridx]

eset_orig <- vector("list", length(esets))

# Convert eset list to set of matrices
for(i in 1:length(esets)){
	eset_orig[[i]] <- t(exprs(esets[[i]]))
}

# Work with the intersection of rows
cn <- lapply(eset_orig, colnames)
cn_int <- Reduce(intersect, cn)

edat <- vector("list", length(esets))

for(i in 1:length(edat)){
	edat[[i]] <- eset_orig[[i]][,cn_int]
}

edat_orig <- edat

# Normalize the columns
for(i in 1:length(edat_orig)){
	edat_orig[[i]] <- apply(edat_orig[[i]], 2, scale)
}

cl<-makeCluster(detectCores() - 1)
registerDoParallel(cl)

lasso_all_nonl_bad <- tree_all_nonl_bad <- nnet_all_nonl_bad <- rf_all_nonl_bad <- mom_all_nonl_bad <- pan6_all_nonl_bad <-  vector("list", 4)
perturb <- c(0.25, 1, 5, 10)
nrun <- 100

system.time(for(i in 1:length(perturb)){

	set.seed(32084, kind = "L'Ecuyer-CMRG")
	pan6_all_nonl_bad[[i]] <- foreach(idx = 1:nrun, 
      	  .combine = rbind)  %dopar% {
			#multistudysim(rffit, rfpred, 0.25, perturb[i], perturb[i], edat_orig, simtype = "nonl")
			pan6sim(0.25, perturb[i], perturb[i], edat_orig, simtype = "nonl")
      	 }
})

tmp <- do.call(rbind, lapply(pan6_all, function(x){apply(sqrt(x), 2, mean)}))
rownames(tmp) <- c(0.25, 1, 2, 10)
xtable(tmp)

stopImplicitCluster()

labels <- c("LASSO", "CART", "Neural Net", "Mas-o-Menos", "Random Forest", "L-C-N-M")
names <- c("lasso_all", "tree_all", "nnet_all", "mom_all", "rf_all", "pan6_all")
names <- paste0(names, "_bad")

par(mfrow = c(6,1))
par(cex = 0.6)
par(mar = c(0, 0, 0, 0), oma = c(6, 8, 0.5, 12))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

mins <- whichmins <- vector("list", 6)
for(i in 1:6){
	tmp <- do.call(rbind, lapply(eval(parse(text=names[i])), function(x){apply(sqrt(x), 2, median)}))
	matplot(tmp[,-c(5,6,8,9,10)]/tmp[,10], type = "l", lwd = 3.5, ylim = c(0.25, 2), yaxt = "n", xaxt = "n", xlab ="", ylab = "")
	mins[[i]] <- apply(tmp[,-c(5,6,8,9)], 1, min)
	whichmins[[i]] <- apply(tmp[,-c(5,6,8,9)], 1, which.min)
	abline(h = 1, lty = "dashed")

	if(i == 3){
		legend(4.05, 2, c("Merged","Avg", "n-Avg", "CS-Avg", "Stack OLS", "Study OLS"),col=seq_len(ncol(tmp)),cex=2,fill=seq_len(ncol(tmp)), bty = "n", xpd = NA)
	}

	axis(2, at = c(0.5,1,1.5), cex.axis = 1.5)
	
	if(i == 6){
		axis(1, at = 1:4, labels = perturb, cex.axis = 1.5)
	}

	if(i == 3){
		mtext("Average RMSE Ratio", side = 2, line = 3.7, adj = 0.79, cex= 1.5)
	}

}

mbox <- do.call(cbind, mins)
wmbox <- do.call(cbind, whichmins)

matplot(mbox, type = "l", col = "black", ylim = c(min(mbox) - 0.2,max(mbox) + 0.2), cex =2, xaxt = "n", xlab = "Coefficient Perturbation Window", ylab = "", cex.axis = 1, cex.lab = 1.5, xlim = c(0.95, 4.05))
axis(1, at = 1:4, labels = perturb, cex.axis = 1)
#axis(2, at = floor(min(mbox)):ceiling(max(mbox)), labels = F, tck = -0.015)
#mtext(floor(min(mbox)):ceiling(max(mbox)), 2, line = 0.5, at = floor(min(mbox)):ceiling(max(mbox)))
mtext("Average Validation RMSE", 2, line = 2, cex = 1.5)

for(i in 1:ncol(mbox)){
	for(j in 1:4){
		shadowtext(j, mbox[j,i], labels = LETTERS[i], col = wmbox[j,i], r=0.12, cex = 2)
	}
}


##################
# Multi-ensemble #
##################

nrun <- 100
multiens_out_tree <- matrix(NA, nrun, 4)

set.seed(32084)
for(i in 1:nrun){
	multiens_out[i,] <- multiensemble(rffit, rfpred, treefit_ensemble, treepred, 0.25, 10, 0.25, 100, edat_orig, simtype = "normal")      	 
}

xtable(t(as.matrix(apply(sqrt(multiens_out), 2, mean))))