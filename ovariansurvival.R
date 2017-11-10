# Make sure to run curovsim.R first so that esets
# are loaded. Use of C-index as perfomance metric is possible
# but commented out here.

# Build 200-gene signatures on each dataset and evaluate on
# every other dataset

sigpred <- function(sig, dat){
	subdat <- exprs(dat)[names(sig),]
	# Using this to rescale the predictions
	#coxph(dat$y ~ as.vector(sig %*% subdat))$linear.predictors

	as.vector(sig %*% subdat)
}

cindex <- function(eset, pred, nitr = 100){
	eval.data <- cbind(eset$y[,1], eset$y[, 2], pred)
	##Calculate Cox coefficient of scaled risk score, and interval:
	tryCatch(Inf.Cval(mydata=eval.data, tau=1460, itr = nitr)$Dhat,
			error = function(e){NA})
}

# Generate predictions & c-indices

# Rows are predictors, columns are datasets
cmat <- matrix(NA, length(edat), length(edat))
coxmat <- matrix(NA, length(edat), length(edat))
preds <- vector("list", length(edat))

ngenes <- 200

for(i in 1:length(edat)){
	cur <- edat[[i]]
	coxout <- rowCoxTests(exprs(cur), cur$y)
	coxout <- coxout[order(coxout[,3]),]
	tcoefs <- coxout[1:ngenes,1]
	names(tcoefs) <- rownames(coxout)[1:ngenes]

	coxmat[i,] <- unlist(lapply(edat, function(dat){
				preds <- sigpred(tcoefs, dat)
				coxph(dat$y ~ preds)$coef
				}))

# If you want to use c-index instead	
#	cmat[i,] <- unlist(lapply(edat, function(dat){
#			preds <- sigpred(coefs, dat)
#			cindex(dat, as.vector(preds))
#		}))

	preds[[i]] <- lapply(edat, function(dat){
			sigpred(tcoefs, dat)
			#cindex(dat, as.vector(preds))
		})
}

coxmat <- do.call(rbind, lapply(preds, function(predlist){
		mapply(function(pred, dat){
			coxph(dat$y ~ pred)$coef
		}, predlist, edat)
		}))

coxmat <- exp(coxmat)

# meta-analytic model
# need survHD

coxout <- lapply(edat, function(eset){
			rowCoxTests(exprs(eset), eset$y)
			})

coefs <- do.call(cbind, lapply(coxout, "[[", 1))
ses <- do.call(cbind, lapply(coxout, "[[", 2))

coefs <- lapply(seq_len(nrow(coefs)), function(i) coefs[i,])
ses <- lapply(seq_len(nrow(ses)), function(i) ses[i,])

metamodel_fit <- function(coefs, ses, rmidx, rn){

	curcoefs <- lapply(coefs, function(x){x[-rmidx]})
	curses <- lapply(ses, function(x){x[-rmidx]})
		
	cl <- makeCluster(4)
	clusterExport(cl, c("curcoefs", "curses"), envir = environment())
	tmp <- do.call(rbind, clusterMap(cl, function(a,b){
			out <- metafor::rma(yi = a, sei = b,method = "DL")
			c(out$beta, out$pval)},
		 curcoefs, curses, .scheduling = "dynamic"))
	stopCluster(cl)

	rownames(tmp) <- rn
	colnames(tmp) <- c("Estimate", "Pval")

	tmp[order(tmp[,2]),][1:200,1]
}


na_aug <- function(vec, rmidx, len){

	tmp <- rep(NA, len)
	offset <- 0
	i <- 1

	while(i <= len){
		if(!(i %in% rmidx)){
			tmp[i] <- vec[i - offset]
		} else {
			offset <- offset + 1
		}
		i <- i + 1
	}
	
	tmp
}

ss <- unlist(lapply(edat, ncol))

wtd_mean <- function(dat_train, dat_test, preds, cmat, rmidx, ss){
	# Only want the held-out datasets (rmidx)
	curpred <- lapply(preds, "[", rmidx)
	curpred <- lapply(curpred, "[", 1:length(rmidx))
	curpred <- sapply(1:length(rmidx), function(x){do.call(rbind, lapply(curpred, "[[", x))})

	# Drop rows and columns in cmat
	curcmat <- abs(cmat)
	curcmat[rmidx,] <- rep(NA, ncol(curcmat))
	curcmat[,rmidx] <- rep(NA, nrow(curcmat))
	diag(curcmat) <- NA

	# For cross-study weighted mean
	wts <- absnorm(rowMeans(curcmat, na.rm = TRUE))
	
	# For unweighted
	unwt <- wts
	unwt[which(!is.na(unwt))] <- 1/(nrow(curcmat) - length(rmidx))

	# For ss weighted
	curss <- ss
	curss[which(is.na(unwt))] <- NA
	curss <- absnorm(curss)

	# For regression weighting

	trainpred <- preds[-rmidx]
	trainpred <- lapply(trainpred, function(x){x[-rmidx]})

	# For stacked regression
	
	stackpred <- lapply(trainpred, unlist)
	stackpred <- do.call(cbind, stackpred)
	colnames(stackpred) <- paste0("V", 1:length(dat_train))

	stacky <- Surv(do.call(rbind, lapply(dat_train, "[[", "y")))
	stackpred <- data.frame(stackpred)
	stackpred$y <- stacky
	stackcox <- coxph(y ~ ., data = stackpred)

	stack_wt <- na_aug(absnorm(stackcox$coef), rmidx, length(curss))

	# For study-specific regression

	study_cf <- matrix(NA, length(dat_train), length(dat_train))

	for(i in 1:length(dat_train)){
		curpr <- t(do.call(rbind, lapply(trainpred, "[[", i)))
		study_cf[i,] <- coxph(dat_train[[i]]$y ~ curpr)$coef
	}
	
	study_wt <- na_aug(absnorm(colMeans(abs(study_cf))), rmidx, length(curss))
	
	prfinal_unwt <- lapply(curpred, function(x){apply(x, 2, function(y){sum(y*unwt, na.rm = T)})})
	prfinal_sswt <- lapply(curpred, function(x){apply(x, 2, function(y){sum(y*curss, na.rm = T)})})
	prfinal_wt <- lapply(curpred, function(x){apply(x, 2, function(y){sum(y*wts, na.rm = T)})})
	prfinal_stack_reg <- lapply(curpred, function(x){apply(x, 2, function(y){sum(y*stack_wt, na.rm =T)})})
	prfinal_study_reg <- lapply(curpred, function(x){apply(x, 2, function(y){sum(y*study_wt, na.rm =T)})})

	cind_unwt <- cind_sswt <- cind_wt <- cind_sr <- cind_stud <- vector("numeric", length(rmidx))

	for(i in 1:length(rmidx)){
		# for c-index
		#cind_unwt[i] <- cindex(dat_test[[i]], prfinal_unwt[[i]], nitr = 10)
		#cind_sswt[i] <- cindex(dat_test[[i]], prfinal_sswt[[i]], nitr = 10)
		#cind_wt[i] <- cindex(dat_test[[i]], prfinal_wt[[i]], nitr = 10)
		#cind_sr[i] <- cindex(dat_test[[i]], prfinal_stack_reg[[i]], nitr = 10)
		#cind_stud[i] <- cindex(dat_test[[i]], prfinal_study_reg[[i]], nitr = 10)

		# for cox coef
		cind_unwt[i] <- coxph(dat_test[[i]]$y ~ prfinal_unwt[[i]])$coef
		cind_sswt[i] <- coxph(dat_test[[i]]$y ~ prfinal_sswt[[i]])$coef
		cind_wt[i] <- coxph(dat_test[[i]]$y ~ prfinal_wt[[i]])$coef
		cind_sr[i] <- coxph(dat_test[[i]]$y ~ prfinal_stack_reg[[i]])$coef
		cind_stud[i] <- coxph(dat_test[[i]]$y ~ prfinal_study_reg[[i]])$coef

	}

	cbind(cind_unwt, cind_sswt, cind_wt, cind_sr, cind_stud, cmat[15,rmidx], sum(sign(stack_wt) < 0, na.rm = T)/(length(stack_wt) - length(rmidx)), sum(sign(study_wt) < 0, na.rm = T)/(length(study_wt) - length(rmidx)))
}

cl <- makeCluster(4)

set.seed(109304)
nsim <- 250

outmat <- matrix(NA, nsim*3, 7)
merge_outvec <- vector("numeric", nsim*3)
rmidx_vec <- vector("list", nsim)
stack_sgn <- study_sgn <- vector("numeric", nsim)

for(i in 1:nsim){
	rmidx <- sample(1:length(edat), 3)
	rmidx_vec[[i]] <- rmidx
	dat_train <- edat[-rmidx]
	dat_test <- edat[rmidx]

	y_all <- Surv(do.call(rbind, lapply(dat_train, function(x){x$y})))
	merge_all <- do.call(cbind, lapply(dat_train, exprs))

	# for single case
	mergeout <- rowCoxTests(merge_all, y_all)
	mergeout <- mergeout[order(mergeout[,3]),]
	mergecoefs <- mergeout[1:200,1]
	names(mergecoefs) <- rownames(mergeout)[1:200]

	mergepreds <- lapply(dat_test, function(x){sigpred(mergecoefs, x)})

	# c-index 
	#merge_outvec[(3*(i-1) + 1):(3*(i-1) + 3)] <- mapply(cindex, dat_test, mergepreds, nitr = 10)

	# cox coef
	merge_outvec[(3*(i-1) + 1):(3*(i-1) + 3)] <- mapply(function(dat, preds){coxph(dat$y ~ preds)$coef}, dat_test, mergepreds)

	# Meta-analysis

	curcoefs <- lapply(coefs, function(x){x[-rmidx]})
	curses <- lapply(ses, function(x){x[-rmidx]})

	clusterExport(cl, c("curcoefs", "curses"))
	tmp <- do.call(rbind, clusterMap(cl, function(a,b){
			out <- metafor::rma(yi = a, sei = b,method = "DL")
			c(out$beta, out$pval)},
		 curcoefs, curses))

	rownames(tmp) <- rownames(dat_train[[1]])
	colnames(tmp) <- c("Estimate", "Pval")

	metasig <- tmp[order(tmp[,2]),][1:200,1]

	metapreds <- lapply(dat_test, sigpred, sig = metasig)

	c_meta <- vector("numeric", length(metapreds))

	for(j in 1:length(dat_test)){
		# c-index
		#c_meta[j] <- cindex(dat_test[[j]], as.vector(metapreds[[j]]), nitr = 10)
		c_meta[j] <- coxph(dat_test[[j]]$y ~ as.vector(metapreds[[j]]))$coef
	}
		
	wt_out <- wtd_mean(dat_train, dat_test, preds, coxmat, rmidx, ss)
	outmat[(3*(i-1) + 1):(3*(i-1) + 3),] <- cbind(wt_out[,1:6], c_meta)
	stack_sgn[i] <- wt_out[1,7]
	study_sgn[i] <- wt_out[1,8]	
}

outmat <- cbind(merge_outvec, outmat)
colnames(outmat) <- c("Merged", "Avg", "n-Avg", "CS-Avg", "Reg-s", "Reg-a", "TCGA", "Meta")

# Compute absolute hazard relative to OLS-a
final <- exp(abs(outmat))/exp(abs(outmat))[,6]

boxplot(final[,c(1,8,7,2,3,4,5,6)], col = c(rep("orange", 2), "lightblue", rep("purple3", 2), rep("white", 3)), ylab = "Performance Ratio",
	  cex.lab = 2, cex.axis = 2, yaxt = "n")
axis(2, at = c(0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15), cex.axis = 1.5)
abline(h = 1, lty = "dashed")
abline(h = median(final[,4]), lty = "dashed")