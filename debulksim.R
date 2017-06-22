# Load the original list of esets from curatedOvarianData, e.g.
# load("061417_esets.Rda")

# need to remove datasets 2,3,6,7,12 from esets because debulking is completely missing

eset_db <- esets[-c(2,3,6,7,12)]
genes <- Reduce(intersect, lapply(eset_db, rownames))

# Generate debulking outcome

for(i in 1:length(eset_db)){
	cury <- 1 * as.vector(pData(eset_db[[i]])[,"debulking"]=="optimal")
	
	eset_db[[i]] <- t(exprs(eset_db[[i]])[genes,])
	eset_db[[i]] <- data.frame(cbind(cury, eset_db[[i]]))
	colnames(eset_db[[i]])[1] <- "y"
	
	if(sum(is.na(cury)) > 0){
		eset_db[[i]] <- eset_db[[i]][-which(is.na(cury)),]
	}

	eset_db[[i]]$y <- as.factor(eset_db[[i]]$y)
}

# Run debulking simulation iteration

rundebulksim <- function(modfit, modpred, eset_db){

	set.seed(32804)	
	ndat <- 11
	ntrain <- 8
	modfit <- treefit
	modpred <- treepred

		edat <- eset_db

		edat <- edat[sample(1:length(edat))] # randomize order

		matstack <- do.call(rbind, edat[1:ntrain])

		mod0 <- modfit(matstack)

		mods <- vector("list", ntrain)

		mses <- matrix(NA, ntrain, ntrain)

		for(i in 1:ntrain){

			mods[[i]] <- modfit(edat[[i]])

			mses[i,] <- unlist(lapply(edat[1:ntrain], function(x){
						preds <- modpred(mods[[i]], newdata = x[,-1])
						LogLossBinary(x[,"y"], preds)}
					))

		}

		diag(mses) <- NA

		tt <- apply(mses, 1, mean, na.rm = T)
		weights <- absnorm(tt, max.norm = TRUE)

		nk <- unlist(lapply(edat[1:ntrain], nrow))
		nwts <- absnorm(nk)

		outmat <- matrix(NA, ndat - ntrain, 4)

		for(i in (ntrain + 1):ndat){
			merged <- modpred(mod0, newdata = edat[[i]][,-1])
			allmod <- do.call(rbind,lapply(mods, modpred, newdata = edat[[i]][,-1]))
			avg <- colMeans(allmod)
			cs_wt <- apply(allmod, 2, function(x){sum(weights*x)})
			n_wt <- apply(allmod, 2, function(x){sum(nwts*x)})
			cury <- edat[[i]][,"y"]

			outmat[i - ntrain,] <- c(LogLossBinary(cury, merged), LogLossBinary(cury, avg), LogLossBinary(cury, n_wt), LogLossBinary(cury, cs_wt))	
		}
	colMeans(outmat)
}

nrun <- 200
tree_debulk <- rf_debulk <- matrix(NA, nrun, 4)

set.seed(32804)
system.time(for(i in 1:nrun){
	tree_debulk[i,] <- rundebulksim(treefit_fs, treepred_fs, eset_db)
})

set.seed(32804)
system.time(for(i in 1:nrun){
	rf_debulk[i,] <- rundebulksim(rf_fs, rf_fs, eset_db)
})

# Generate results table - table 1 in main text
db_mat <- matrix(NA, 2, 3)
rownames(db_mat) <- c("CART", "RF")
colnames(db_mat) <- c("Merged", "Avg", "n_avg")

db_mat[1,] <- c(sum(tree_debulk[,4] < tree_debulk[,1]), sum(tree_debulk[,4] < tree_debulk[,2]), sum(tree_debulk[,4] < tree_debulk[,3]))
db_mat[2,] <- c(sum(rf_debulk[,4] < rf_debulk[,1]), sum(rf_debulk[,4] < rf_debulk[,2]), sum(rf_debulk[,4] < rf_debulk[,3]))
db_mat <- (db_mat/200)*100

xtable(db_mat, digits = 1)

# Optional boxplots
#par(mfrow = c(1,2))
#db_plot <- tree_debulk
#boxplot(db_plot[,1] - db_plot[,4], db_plot[,2] - db_plot[,4], db_plot[,3] - db_plot[,4],
#	   main = "Difference in CS RMSE and * RMSE, 100 iterations", names = c("Merged", "Avg.", "SS"),
#	   ylim = c(-0.05, 1))
#abline(h = 0, lty = "dashed")

#db_plot <- rf_debulk
#boxplot(db_plot[,1] - db_plot[,4], db_plot[,2] - db_plot[,4], db_plot[,3] - db_plot[,4],
#	   main = "Difference in CS RMSE and * RMSE, 100 iterations", names = c("Merged", "Avg.", "SS"),
#	   ylim = c(-0.05, 1))
#abline(h = 0, lty = "dashed")