# Set up a dataset split for one simulation iteration
#
# Input:
# edat_orig - orignal (cleaned) list of curatedOvarianData datasets
# ndat - length(edat_orig)
# nvar - number of features to reduce each dataset to
# simtype - one of "normal", "slash", "nonl" representing three scenarios
#		described in the supplement
# good - Perturbation level for "good" i.e. low-perturbation studies
# bad - Perturbation level for "bad" i.e. high-perturbation studies
# val - Perturbation level for "val" i.e. validation studies
#
# Output:
# edat - a list of esets with generated Y

init_data <- function(edat_orig, ndat, nvar, simtype, good, bad, val){

	edat <- edat_orig
	edat <- edat[sample(1:ndat)] # Randomize dataset order

	idx <- sample(1:ncol(edat[[1]]), nvar)
	for(i in 1:ndat){
		edat[[i]] <- edat[[i]][,idx]
	}

	# Generate Y

	if(simtype == "nonl"){
		ncoef <- sample(3:nvar, 1) # Use this for interaction run
	} else {
		ncoef <- sample(2:nvar, 1)
	}

	coefs <- sample(c(runif(round(ncoef/2), -5, -0.5), runif(ncoef - round(ncoef/2), 0.5, 5)))
	vars <- sample(1:ncol(edat[[1]]), ncoef)

	for(i in 1:ndat){

		curcoefs <- sapply(coefs, function(x){runif(1, x - (i<=5)*good - (i > 5 & i <= 10)*bad - (i > 10)*val , x + (i<=5)*good + (i > 5 & i <= 10)*bad + (i > 10)*val)})
			

		if(simtype == "slash"){
			y <- (edat[[i]][,vars] %*% curcoefs) + 0.05*cbind(rnorm(nrow(edat[[i]]))/runif(nrow(edat[[i]]))) # Added "slash" noise
		} else if(simtype == "nonl"){
			y <- (edat[[i]][,vars] %*% curcoefs) + 2.8*edat[[i]][,vars[1]]*edat[[i]][,vars[2]] - 1.6*edat[[i]][,vars[1]]*edat[[i]][,vars[3]] + cbind(rnorm(nrow(edat[[i]]))) # Added interaction terms
		} else {
			y <- (edat[[i]][,vars] %*% curcoefs) + cbind(rnorm(nrow(edat[[i]]))) # Added noise
		}

		edat[[i]] <- cbind(y, edat[[i]])
		edat[[i]] <- as.data.frame(edat[[i]])
		colnames(edat[[i]]) <- c("y", paste0("V", 1:nvar))
	}

	edat
}

# Run main simulation (Figure 3)
#
# Input:
# modfit - Function to fit a learner to a dataset (from multistudysim_utils.R)
# modpred - Function to make predictions from a learner (from multistudysim_utils.R)
# good - Perturbation level for "good" i.e. low-perturbation studies
# bad - Perturbation level for "bad" i.e. high-perturbation studies
# val - Perturbation level for "val" i.e. validation studies
# edat_orig - orignal (cleaned) list of curatedOvarianData datasets
# simtype - one of "normal", "slash", "nonl" representing three scenarios
#		described in the supplement
#
# Output:
# colMeans(outmat) - Average MSE across validation datasets for each weighting scheme


multistudysim <- function(modfit, modpred, good, bad, val, edat_orig, simtype = "normal"){
	
	ndat <- length(edat_orig)
	ntrain <- 10
	nvar <- 20
	modfit <- modfit
	modpred <- modpred

	edat <- init_data(edat_orig, ndat, nvar, simtype, good, bad, val)

	# Learn the merged predictor

	matstack <- do.call(rbind, edat[1:ntrain])

	mod0 <- modfit(matstack)

	mods <- allpreds <- vector("list", ntrain)

	mses <- matrix(NA, ntrain, ntrain)

	for(i in 1:ntrain){

		mods[[i]] <- modfit(edat[[i]])

		mses[i,] <- unlist(lapply(edat[1:ntrain], function(x){
					preds <- modpred(mods[[i]], newdata = x[,-1])
					mean((preds - x[,"y"])^2)}
				))

		curpreds <-  lapply(edat[1:ntrain], function(x){modpred(mods[[i]], newdata = x[,-1])})
		allpreds <- mapply(cbind, allpreds, curpreds, SIMPLIFY = FALSE)
	}

	diag(mses) <- NA

	# CS Weights
	tt <- apply(mses, 1, mean, na.rm = T)
	weights <- absnorm(sqrt(tt), max.norm = TRUE)
	
	# Sample Size Weights
	nk <- unlist(lapply(edat, nrow))
	nwts <- absnorm(nk[1:ntrain])

	# Regression: stacked (intercept and no intercept)
	predstack <- do.call(rbind, allpreds)
	coefs_stack_noint <- nnls::nnls(predstack, matstack$y)$x
	coefs_stack_int <- nnls::nnls(cbind(rep(1,nrow(predstack)),predstack), matstack$y)$x

	# Just a safeguard against full collinearity, although I think we are OK with nnls now
	coefs_stack_noint[which(is.na(coefs_stack_noint))] <- 0
	coefs_stack_int[which(is.na(coefs_stack_int))] <- 0

	coefs_stack_noint_norm <- absnorm(coefs_stack_noint)

	# Regression: study-specific (intercept and no intercept)
	coefs_ss_noint <- mapply(function(x,y){nnls::nnls(y,x[,"y"])$x}, edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
	coefs_ss_noint <- colMeans(do.call(rbind, coefs_ss_noint), na.rm = T)
	coefs_ss_int <- mapply(function(x,y){nnls::nnls(cbind(rep(1,nrow(y)),y),x[,"y"])$x}, edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
	coefs_ss_int <- colMeans(do.call(rbind, coefs_ss_int), na.rm = T)

	coefs_ss_noint[which(is.na(coefs_ss_noint))] <- 0
	coefs_ss_int[which(is.na(coefs_ss_int))] <- 0

	coefs_ss_noint_norm <- absnorm(coefs_ss_noint)

	outmat <- matrix(NA, ndat - ntrain, 10)
	colnames(outmat) <- c("Merged", "Unweighted", "Sample_Weighted", "CS_Weighted",
				    "Stack_noint", "Stack_noint_norm", "Stack_int",
				    "SS_noint", "SS_noint_norm", "SS_int")

	for(i in (ntrain + 1):ndat){
		
		# Predictions from each SSL
		allmod <- do.call(rbind,lapply(mods, modpred, newdata = edat[[i]][,-1]))
	
		# Merged predictor
		merged <- modpred(mod0, newdata = edat[[i]][,-1])			

		# Unweighted average
		unweighted <- colMeans(allmod)

		# sample size weighted
		sample_wtd <- apply(allmod, 2, function(x){sum(nwts*x)})

		# cross-study weighted
		cs_wtd <- apply(allmod, 2, function(x){sum(weights*x)})

		# regression: stacked (noint, int, each normed)
		stack_noint <- apply(allmod, 2, function(x){sum(coefs_stack_noint*x)})
		stack_noint_norm <- apply(allmod, 2, function(x){sum(coefs_stack_noint_norm*x)})
		stack_int <- apply(allmod, 2, function(x){coefs_stack_int[1] + sum(coefs_stack_int[-1]*x)})

		# regression: study_specific (noint, int, noint normed)
		ss_noint <- apply(allmod, 2, function(x){sum(coefs_ss_noint*x)})
		ss_noint_norm <- apply(allmod, 2, function(x){sum(coefs_ss_noint_norm*x)})
		ss_int <- apply(allmod, 2, function(x){coefs_ss_int[1] + sum(coefs_ss_int[-1]*x)})

		cury <- edat[[i]][,"y"]

		outmat[i - ntrain,] <- c(mean((cury - merged)^2), mean((cury - unweighted)^2), mean((cury - sample_wtd)^2), mean((cury - cs_wtd)^2),
						 mean((cury - stack_noint)^2), mean((cury - stack_noint_norm)^2), mean((cury - stack_int)^2),
							 mean((cury - ss_noint)^2), mean((cury - ss_noint_norm)^2), mean((cury - ss_int)^2))	

	}

	colMeans(outmat)
}

# Example run
# system.time(tts <- multistudysim(lassofit, lassopred, 0.25, 5, 0.25, edat_orig))

# Run L-C-N-M panel simulation
#
# Note that learners are already built in and not provided as arguments
#
# Input:
# good - Perturbation level for "good" i.e. low-perturbation studies
# bad - Perturbation level for "bad" i.e. high-perturbation studies
# val - Perturbation level for "val" i.e. validation studies
# edat_orig - orignal (cleaned) list of curatedOvarianData datasets
# simtype - one of "normal", "slash", "nonl" representing three scenarios
#		described in the supplement
#
# Output:
# colMeans(outmat) - Average MSE across validation datasets for each weighting scheme

pan6sim <- function(good, bad, val, edat_orig, simtype = "normal"){
	
	ndat <- length(edat_orig)
	ntrain <- 10
	nvar <- 20
	modfit <- list(lassofit, treefit, nnetfit, momfit)
	modpred <- list(lassopred, treepred, nnetpred, mompred)


	edat <- init_data(edat_orig, nvar, simtype, good, bad, val)
		
	matstack <- do.call(rbind, edat[1:ntrain])

	mod0 <- vector("list", length(modfit))
	mod0preds <- rep(1, nrow(matstack))

	# Fit each learner on all the data
	for(i in 1:length(modfit)){
		mod0[[i]] <- modfit[[i]](matstack)
		mod0preds <- cbind(mod0preds, modpred[[i]](mod0[[i]], newdata = matstack[,-1]))
	}
		
	# Get stacking weights for the full data predictors
	mod0weights_int <- nnls(mod0preds, matstack$y)$x		

	mods <- vector("list", ntrain)
	allpreds <- vector("list", ntrain)

	mses <- matrix(NA, ntrain*length(modfit), ntrain)

	for(i in 1:ntrain){

		mods[[i]] <- lapply(modfit, function(mod){mod(edat[[i]])})

		mses[(4*i-3):(4*i),] <- do.call(cbind, lapply(edat[1:ntrain], function(x){
						preds <- mapply(function(mod, modpred){modpred(mod, newdata = x[,-1])}, mods[[i]], modpred)
						colMeans((preds - x[,"y"])^2)}
					))

		curpreds <-  lapply(edat[1:ntrain], function(x){mapply(function(mod, modpred){modpred(mod, newdata = x[,-1])}, mods[[i]], modpred)})
		allpreds <- mapply(cbind, allpreds, curpreds, SIMPLIFY = FALSE)
	}

	#diag(mses) <- NA # mses is a non-square matrix

	# CS Weights
	tt <- apply(mses, 1, mean, na.rm = T)
	weights <- absnorm(sqrt(tt), max.norm = TRUE)
	
	# Sample Size Weights
	nk <- unlist(lapply(edat, nrow))[1:ntrain]
	nk <- rep(nk, each = length(modfit))
	nwts <- absnorm(nk)

	# Regression: stacked (intercept and no intercept)
	predstack <- do.call(rbind, allpreds)
	coefs_stack_noint <- nnls::nnls(predstack, matstack$y)$x
	coefs_stack_int <- nnls::nnls(cbind(rep(1,nrow(predstack)),predstack), matstack$y)$x

	coefs_stack_noint[which(is.na(coefs_stack_noint))] <- 0
	coefs_stack_int[which(is.na(coefs_stack_int))] <- 0

	coefs_stack_noint_norm <- absnorm(coefs_stack_noint)

	# Regression: study-specific (intercept and no intercept)
		
	coefs_ss_noint <- mapply(function(x,y){nnls::nnls(y,x[,"y"])$x}, edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
	coefs_ss_noint <- colMeans(do.call(rbind, coefs_ss_noint), na.rm = T)

	coefs_ss_int <- mapply(function(x,y){nnls::nnls(cbind(rep(1,nrow(y)),y),x[,"y"])$x}, edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
	coefs_ss_int <- colMeans(do.call(rbind, coefs_ss_int), na.rm = T)

	coefs_ss_noint[which(is.na(coefs_ss_noint))] <- 0
	coefs_ss_int[which(is.na(coefs_ss_int))] <- 0

	coefs_ss_noint_norm <- absnorm(coefs_ss_noint)

	outmat <- matrix(NA, ndat - ntrain, 10)
	colnames(outmat) <- c("Merged", "Unweighted", "Sample_Weighted", "CS_Weighted",
				    "Stack_noint", "Stack_noint_norm", "Stack_int",
				    "SS_noint", "SS_noint_norm", "SS_int")


	for(i in (ntrain + 1):ndat){
		
		# Predictions from each SSL
		allmod <- t(do.call(cbind, lapply(mods, function(x){mapply(function(mod, modpred){modpred(mod, newdata = edat[[i]][,-1])}, x, modpred)})))
	
		# Merged predictor
		merged <- t(mapply(function(mod, modpred){modpred(mod, newdata = edat[[i]][,-1])}, mod0, modpred))
		merged <- apply(merged, 2, function(x){mod0weights_int[1] + sum(mod0weights_int[-1]*x)})

		# Unweighted average
		unweighted <- colMeans(allmod)

		# sample size weighted
		sample_wtd <- apply(allmod, 2, function(x){sum(nwts*x)})

		# cross-study weighted
		cs_wtd <- apply(allmod, 2, function(x){sum(weights*x)})

		# regression: stacked (noint, int, each normed)
		stack_noint <- apply(allmod, 2, function(x){sum(coefs_stack_noint*x)})
		stack_noint_norm <- apply(allmod, 2, function(x){sum(coefs_stack_noint_norm*x)})
		stack_int <- apply(allmod, 2, function(x){coefs_stack_int[1] + sum(coefs_stack_int[-1]*x)})

		# regression: study_specific (noint, int, noint normed)
		ss_noint <- apply(allmod, 2, function(x){sum(coefs_ss_noint*x)})
		ss_noint_norm <- apply(allmod, 2, function(x){sum(coefs_ss_noint_norm*x)})
		ss_int <- apply(allmod, 2, function(x){coefs_ss_int[1] + sum(coefs_ss_int[-1]*x)})

		cury <- edat[[i]][,"y"]
		outmat[i - ntrain,] <- c(mean((cury - merged)^2), mean((cury - unweighted)^2), mean((cury - sample_wtd)^2), mean((cury - cs_wtd)^2),
						 mean((cury - stack_noint)^2), mean((cury - stack_noint_norm)^2), mean((cury - stack_int)^2),
						 mean((cury - ss_noint)^2), mean((cury - ss_noint_norm)^2), mean((cury - ss_int)^2))	

	}

	colMeans(outmat)
}

# Run "ensemble of ensembles" simulation (supplemental)
#
# Input:
# modfit - Function to fit a learner to a dataset (from multistudysim_utils.R)
# modpred - Function to make predictions from a learner (from multistudysim_utils.R)
# modfit_ens - Function to fit an ensembling learner to a dataset (from multistudysim_utils.R)
# modpred_ens - Function to make predictions from an ensembling learner (from multistudysim_utils.R)
# good - Perturbation level for "good" i.e. low-perturbation studies
# bad - Perturbation level for "bad" i.e. high-perturbation studies
# val - Perturbation level for "val" i.e. validation studies
# ntrees - Number of SSLs to fit on each dataset
# edat_orig - orignal (cleaned) list of curatedOvarianData datasets
# simtype - one of "normal", "slash", "nonl" representing three scenarios
#		described in the supplement
#
# Output:
# colMeans(outmat) - Average MSE across validation datasets for each weighting scheme

multiensemble <- function(modfit, modpred, modfit_ens, modpred_ens, good, bad, val, ntrees = 100, edat_orig, simtype = "normal"){
	
	ndat <- length(edat_orig)
	ntrain <- 10
	nvar <- 20
	modfit <- modfit
	modfit_ens <- modfit_ens
	modpred <- modpred
	modpred_ens <- modpred_ens
	ntrees <- ntrees

	edat <- init_data(edat_orig, nvar, simtype, good, bad, val)


	alltrees <- vector("list", ntrain)
	for(i in 1:length(alltrees)){
		alltrees[[i]] <- vector("list", ntrees)
	}

	matstack <- do.call(rbind, edat[1:ntrain])

	mod0 <- modfit(matstack)

	mods <- allpreds <- vector("list", ntrain)

	mses <- matrix(NA, ntrain, ntrain)

	for(i in 1:ntrain){

		mods[[i]] <- modfit(edat[[i]])

		mses[i,] <- unlist(lapply(edat[1:ntrain], function(x){
					preds <- modpred(mods[[i]], newdata = x[,-1])
					mean((preds - x[,"y"])^2)}
				))

		curpreds <-  lapply(edat[1:ntrain], function(x){modpred(mods[[i]], newdata = x[,-1])})
		allpreds <- mapply(cbind, allpreds, curpreds, SIMPLIFY = FALSE)
	}

	diag(mses) <- NA

	# CS Weights
	tt <- apply(mses, 1, mean, na.rm = T)
	weights <- absnorm(sqrt(tt), max.norm = TRUE)
	

	# Learn the ensemble - resample n with replacement, take subset of p
	for(i in 1:ntrain){
		for(j in 1:ntrees){
			n_idx <- sample(1:nrow(edat[[i]]), replace = T)
			p_idx <- sample(1:nvar, sample(1:nvar, 1))
			alltrees[[i]][[j]] <- modfit_ens(edat[[i]], n_idx, p_idx)
		}
	}

	outmat <- matrix(NA, ndat - ntrain, 4)

	for(i in (ntrain + 1):ndat){

		allmod <- do.call(rbind,lapply(mods, modpred, newdata = edat[[i]][,-1]))
		
		# Merged predictor
		merged <- modpred(mod0, newdata = edat[[i]][,-1])			

		# Unweighted average
		unweighted <- colMeans(allmod)

		# cross-study weighted
		cs_wtd <- apply(allmod, 2, function(x){sum(weights*x)})

		ens_preds <- rowMeans(do.call(cbind, lapply(alltrees, function(x){do.call(cbind, lapply(x, function(y){modpred_ens(y, newdata = edat[[i]][,-1])}))})))
		cury <- edat[[i]][,"y"]

		outmat[i - ntrain,] <- c(mean((cury - merged)^2), mean((cury - unweighted)^2), mean((cury - cs_wtd)^2), mean((cury - ens_preds)^2))

	}

	colMeans(outmat)
}

# Parallelized wrapper for multistudysim/pan6sim
#
# Input:
# modfit - Function to fit a learner to a dataset (from multistudysim_utils.R)
# modpred - Function to make predictions from a learner (from multistudysim_utils.R)
# edat_orig - orignal (cleaned) list of curatedOvarianData datasets
# nrun - Number of iteration simulations
# val - "low" for low-perturbation validation sets, "high" for high-perturbation
# simtype - one of "normal", "slash", "nonl" representing three scenarios
#		described in the supplement
# pan6 - TRUE if pan6sim should be done instead of multistudysim
#
# Output:
# outlist - Average MSE across validation datasets for each weighting scheme,
#		length = number of perturbation levels tested


runsimpar <- function(modfit, modpred, edat_orig, nrun = 100, val = "low", simtype = "normal", pan6 = FALSE){

	cl<-makeCluster(detectCores() - 1)
	registerDoParallel(cl)

	perturb <- c(0.25, 1, 5, 10)

	outlist <-  vector("list", length(perturb))

	nrun <- nrun

	system.time(for(i in 1:length(perturb)){
		if(val == "high"){
			vcur <- perturb[i]
		} else {
			vcur <- 0.25
		}

		set.seed(32084, kind = "L'Ecuyer-CMRG")
		outlist[[i]] <- foreach(idx = 1:nrun, 
      		  .combine = rbind, 
			  .export = c("absnorm", "init_data", "multistudysim", "pan6sim", 
					  "lassofit", "treefit", "nnetfit", "momfit",
					  "lassopred", "treepred", "nnetpred", "mompred"))  %dopar% {
				if(pan6){
					pan6sim(0.25, perturb[i], vcur, edat_orig, simtype = simtype)
				} else {
					multistudysim(modfit, modpred, 0.25, perturb[i], vcur, edat_orig, simtype = simtype)
				}
      	 	}
	})

	stopImplicitCluster()
	outlist
}

# Example Run
# runsimpar(momfit, mompred, edat_orig, nrun = 10, val = "low", simtype = "slash", pan6 = FALSE)

# Plot results of multistudysim + pan6sim
#
# Input:
# outlist - list of output from one set of runs of runsimpar

simplots <- function(outlist){

	par(mfrow = c(length(outlist),1))
	par(cex = 0.6)
	par(mar = c(0, 0, 0, 0), oma = c(6, 8, 0.5, 12))
	par(tcl = -0.25)
	par(mgp = c(2, 0.6, 0))

	mins <- whichmins <- vector("list", 6)
	perturb <- c(0.25, 1, 5, 10)

	for(i in 1:length(outlist)){
		tmp <- do.call(rbind, lapply(outlist[[i]], function(x){apply(sqrt(x), 2, median)}))
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
	mtext("Average Validation RMSE", 2, line = 2, cex = 1.5)

	for(i in 1:ncol(mbox)){
		for(j in 1:4){
			shadowtext(j, mbox[j,i], labels = LETTERS[i], col = wmbox[j,i], r=0.12, cex = 2)
		}
	}

}
