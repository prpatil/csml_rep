

multistudysim <- function(modfit, modpred, good, bad, val, edat_orig, simtype = "normal"){
	
	ndat <- length(edat_orig)
	ntrain <- 10
	nvar <- 20
	modfit <- modfit
	modpred <- modpred

	edat <- edat_orig
	edat <- edat[sample(1:length(edat))] # Randomize dataset order

	idx <- sample(1:ncol(edat[[1]]), nvar)
	for(i in 1:length(edat)){
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
# system.time(tts <- multistudysim(treefit, treepred, 0.25, 5, 0.25, edat_orig))

multiensemble <- function(modfit, modpred, modfit_ens, modpred_ens, good, bad, val, ntrees = 100, edat_orig, simtype = "normal"){
	
	ndat <- length(edat_orig)
	ntrain <- 10
	nvar <- 20
	modfit <- modfit
	modfit_ens <- modfit_ens
	modpred <- modpred
	modpred_ens <- modpred_ens
	ntrees <- ntrees

		edat <- edat_orig
		edat <- edat[sample(1:length(edat))] # Randomize dataset order

		idx <- sample(1:ncol(edat[[1]]), nvar)
		for(i in 1:length(edat)){
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

pan6sim <- function(good, bad, val, edat_orig, simtype = "normal"){
	
	ndat <- length(edat_orig)
	ntrain <- 10
	nvar <- 20
	modfit <- list(lassofit, treefit, nnetfit, momfit)
	modpred <- list(lassopred, treepred, nnetpred, mompred)


		edat <- edat_orig
		edat <- edat[sample(1:length(edat))] # Randomize dataset order

		idx <- sample(1:ncol(edat[[1]]), nvar)
		for(i in 1:length(edat)){
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


		# Here's where things change

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

			#allpreds <- cbind(allpreds, unlist(lapply(edat[1:ntrain], function(x){modpred(mods[[i]], newdata = x[,-1])})))
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