# From curatedOvarianData vignette
source(system.file("extdata", "patientselection.config",package="curatedOvarianData"))
sapply(ls(), function(x) if(!x %in% c("remove.samples", "duplicates")) print(get(x)))
source(system.file("extdata", "createEsetList.R", package = "curatedOvarianData"))

# Save the original eset list - reuse this list for other analyses
# e.g. save(esets, file = "061417_esets.Rda")

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


labels <- c("LASSO", "CART", "Neural Net", "Mas-o-Menos", "Random Forest", "L-C-N-M")
names <- c("lasso", "tree", "nnet", "mom", "rf", "pan6")

# Note: the modfit/modpred arguments for pan6 don't matter as they are unused

# This produces figure 3 in the main text

lasso_normal_low <- runsimpar(lassofit, lassopred, edat_orig, nrun = 100, val = "low", simtype = "normal", pan6 = FALSE)
tree_normal_low <- runsimpar(treefit, treepred, edat_orig, nrun = 100, val = "low", simtype = "normal", pan6 = FALSE)
nnet_normal_low <- runsimpar(nnetfit, nnetpred, edat_orig, nrun = 100, val = "low", simtype = "normal", pan6 = FALSE)
mom_normal_low <- runsimpar(momfit, mompred, edat_orig, nrun = 100, val = "low", simtype = "normal", pan6 = FALSE)
rf_normal_low <- runsimpar(rffit, rfpred, edat_orig, nrun = 100, val = "low", simtype = "normal", pan6 = FALSE)
pan6_normal_low <- runsimpar(momfit, mompred, edat_orig, nrun = 100, val = "low", simtype = "normal", pan6 = TRUE)

simplots(list(lasso_normal_low, tree_normal_low, nnet_normal_low, mom_normal_low, rf_normal_low, pan6_normal_low))

# These tables appear as tables 1-6 in the supplement

xtable(do.call(rbind, lapply(lasso_normal_low, function(x){apply(sqrt(x), 2, mean)})))
xtable(do.call(rbind, lapply(tree_normal_low, function(x){apply(sqrt(x), 2, mean)})))
xtable(do.call(rbind, lapply(nnet_normal_low, function(x){apply(sqrt(x), 2, mean)})))
xtable(do.call(rbind, lapply(mom_normal_low, function(x){apply(sqrt(x), 2, mean)})))
xtable(do.call(rbind, lapply(rf_normal_low, function(x){apply(sqrt(x), 2, mean)})))
xtable(do.call(rbind, lapply(pan6_normal_low, function(x){apply(sqrt(x), 2, mean)})))

# The following appear in the supplement as figures 1-5

lasso_normal_high <- runsimpar(lassofit, lassopred, edat_orig, nrun = 100, val = "high", simtype = "normal", pan6 = FALSE)
tree_normal_high <- runsimpar(treefit, treepred, edat_orig, nrun = 100, val = "high", simtype = "normal", pan6 = FALSE)
nnet_normal_high <- runsimpar(nnetfit, nnetpred, edat_orig, nrun = 100, val = "high", simtype = "normal", pan6 = FALSE)
mom_normal_high <- runsimpar(momfit, mompred, edat_orig, nrun = 100, val = "high", simtype = "normal", pan6 = FALSE)
rf_normal_high <- runsimpar(rffit, rfpred, edat_orig, nrun = 100, val = "high", simtype = "normal", pan6 = FALSE)
pan6_normal_high <- runsimpar(momfit, mompred, edat_orig, nrun = 100, val = "high", simtype = "normal", pan6 = TRUE)

simplots(list(lasso_normal_high, tree_normal_high, nnet_normal_high, mom_normal_high, rf_normal_high, pan6_normal_high))

lasso_slash_low <- runsimpar(lassofit, lassopred, edat_orig, nrun = 100, val = "low", simtype = "slash", pan6 = FALSE)
tree_slash_low <- runsimpar(treefit, treepred, edat_orig, nrun = 100, val = "low", simtype = "slash", pan6 = FALSE)
nnet_slash_low <- runsimpar(nnetfit, nnetpred, edat_orig, nrun = 100, val = "low", simtype = "slash", pan6 = FALSE)
mom_slash_low <- runsimpar(momfit, mompred, edat_orig, nrun = 100, val = "low", simtype = "slash", pan6 = FALSE)
rf_slash_low <- runsimpar(rffit, rfpred, edat_orig, nrun = 100, val = "low", simtype = "slash", pan6 = FALSE)
pan6_slash_low <- runsimpar(momfit, mompred, edat_orig, nrun = 100, val = "low", simtype = "slash", pan6 = TRUE)

simplots(list(lasso_slash_low, tree_slash_low, nnet_slash_low, mom_slash_low, rf_slash_low, pan6_slash_low))

lasso_slash_high <- runsimpar(lassofit, lassopred, edat_orig, nrun = 100, val = "high", simtype = "slash", pan6 = FALSE)
tree_slash_high <- runsimpar(treefit, treepred, edat_orig, nrun = 100, val = "high", simtype = "slash", pan6 = FALSE)
nnet_slash_high <- runsimpar(nnetfit, nnetpred, edat_orig, nrun = 100, val = "high", simtype = "slash", pan6 = FALSE)
mom_slash_high <- runsimpar(momfit, mompred, edat_orig, nrun = 100, val = "high", simtype = "slash", pan6 = FALSE)
rf_slash_high <- runsimpar(rffit, rfpred, edat_orig, nrun = 100, val = "high", simtype = "slash", pan6 = FALSE)
pan6_slash_high <- runsimpar(momfit, mompred, edat_orig, nrun = 100, val = "high", simtype = "slash", pan6 = TRUE)

simplots(list(lasso_slash_high, tree_slash_high, nnet_slash_high, mom_slash_high, rf_slash_high, pan6_slash_high))

lasso_nonl_low <- runsimpar(lassofit, lassopred, edat_orig, nrun = 100, val = "low", simtype = "nonl", pan6 = FALSE)
tree_nonl_low <- runsimpar(treefit, treepred, edat_orig, nrun = 100, val = "low", simtype = "nonl", pan6 = FALSE)
nnet_nonl_low <- runsimpar(nnetfit, nnetpred, edat_orig, nrun = 100, val = "low", simtype = "nonl", pan6 = FALSE)
mom_nonl_low <- runsimpar(momfit, mompred, edat_orig, nrun = 100, val = "low", simtype = "nonl", pan6 = FALSE)
rf_nonl_low <- runsimpar(rffit, rfpred, edat_orig, nrun = 100, val = "low", simtype = "nonl", pan6 = FALSE)
pan6_nonl_low <- runsimpar(momfit, mompred, edat_orig, nrun = 100, val = "low", simtype = "nonl", pan6 = TRUE)

simplots(list(lasso_nonl_low, tree_nonl_low, nnet_nonl_low, mom_nonl_low, rf_nonl_low, pan6_nonl_low))

lasso_nonl_high <- runsimpar(lassofit, lassopred, edat_orig, nrun = 100, val = "high", simtype = "nonl", pan6 = FALSE)
tree_nonl_high <- runsimpar(treefit, treepred, edat_orig, nrun = 100, val = "high", simtype = "nonl", pan6 = FALSE)
nnet_nonl_high <- runsimpar(nnetfit, nnetpred, edat_orig, nrun = 100, val = "high", simtype = "nonl", pan6 = FALSE)
mom_nonl_high <- runsimpar(momfit, mompred, edat_orig, nrun = 100, val = "high", simtype = "nonl", pan6 = FALSE)
rf_nonl_high <- runsimpar(rffit, rfpred, edat_orig, nrun = 100, val = "high", simtype = "nonl", pan6 = FALSE)
pan6_nonl_high <- runsimpar(momfit, mompred, edat_orig, nrun = 100, val = "high", simtype = "nonl", pan6 = TRUE)

simplots(list(lasso_nonl_high, tree_nonl_high, nnet_nonl_high, mom_nonl_high, rf_nonl_high, pan6_nonl_high))

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