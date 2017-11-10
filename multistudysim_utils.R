# Load packages

library(curatedOvarianData)
library(RcppArmadillo)
library(e1071)
library(ranger)
library(glmnet)
library(rpart)
library(genefilter)
library(nnet)
library(mboost)
library(nnls)
library(foreach)
library(parallel)
library(doParallel)
library(survHD)

# These are the learner/predictor pairs, written to be supplied generically
# to the simulation functions

momfit <- function(data){
	y <- data[,1]
	d <- apply(data, 2, function(x){x/sd(x)})
	d[,1] <- y

	score <- sign(RcppArmadillo::fastLmPure(as.matrix(d[,-1]), d[,1])$coef)/sqrt(ncol(d)-1)
	lp <- d[,-1] %*% score
	beta <- RcppArmadillo::fastLmPure(lp, d[,1])$coef

	list(score, beta)
}

mompred <- function(mod, newdata){
	d <- apply(newdata, 2, function(x){x/sd(x)})
	predscore <- d %*% mod[[1]]
	mod[[2]]*as.vector(predscore)
}

treefit <- function(data, ...){
	mod <- rpart::rpart(y ~ ., data = as.data.frame(data), ...)
	mod <- rpart::prune(mod, mod$cptable[which(abs(diff(mod$cptable[,"xerror"])) < 0.01)[1], "CP"])
	mod
}

treepred <- function(mod, newdata){
	as.vector(predict(mod, newdata=newdata))
}

treefit_fs <- function(data, ...){
	tout <- colttests(as.matrix(data[,-1]), data[,1], tstatOnly = TRUE)
	data <- data[, c("y",rownames(tout)[order(abs(tout[,1]), decreasing = T)[1:20]])]
	mod <- rpart::rpart(y ~ ., data = as.data.frame(data), ...)
	mod
}

treepred_fs <- function(mod, newdata){
	as.vector(predict(mod, newdata=newdata)[,2])
}

treefit_ensemble <- function(data, n_idx, p_idx, ...){
	form <- paste0("y ~ ", paste0("V", p_idx, collapse = " + "))
	mod <- rpart::rpart(form, data = as.data.frame(data), subset = n_idx)
	mod
}

lassofit <- function(data, ...){
	mod <- glmnet::cv.glmnet(x = as.matrix(data[,-1]), y = data[,1], ...)
}

lassopred <- function(mod, newdata){
	as.vector(predict(mod, newx=as.matrix(newdata), s="lambda.1se"))
}

rffit <- function(data, ...){
	ranger::ranger(y ~ ., data = data, write.forest = TRUE)
}

rfpred <- function(mod, newdata){
	predict(mod, data = newdata)$predictions
}

rffit_fs <- function(data, ...){
	tout <- colttests(as.matrix(data[,-1]), data[,1], tstatOnly = TRUE)
	data <- data[, c("y",rownames(tout)[order(abs(tout[,1]), decreasing = T)[1:20]])]
	ranger(y ~ ., data = data, write.forest = TRUE, probability = T)
}

rfpred_fs <- function(mod, newdata){
	predict(mod, data = newdata)$predictions[,2]
}

nnetfit <- function(data, ...){
	nnet::nnet(y ~ ., data = data, size = 10, MaxNWts = 2000, linout = T, trace = F)
}

nnetpred <- function(mod, newdata){
	as.vector(predict(mod, newdata = newdata))
}

boostfit <- function(data, ...){
	mboost::mboost(y ~ ., data = data, baselearner = "bbs")
}

boostpred <- function(mod, newdata){
	as.vector(predict(mod, newdata = newdata))
}

# Normalize by subtracting off the worst performer
# Takes as input a vector of weights and a flag indicating whether
# high values are bad (max.norm = TRUE) or low values are bad (max.norm = FALSE)

absnorm <- function(vec, max.norm = FALSE){
	sgn <- sign(vec)
	vec <- abs(vec)
	if(max.norm){
		mvec <- max(vec)	
	} else {
		mvec <- min(vec)
	}

	vec <- abs(vec - mvec)
	sgn*(vec/sum(vec))
}

# Ues to shadow letters for Figure 3 plot bottom panel
# credit to Greg Snow c/o 
# https://stackoverflow.com/questions/25631216/r-is-there-any-way-to-put-border-shadow-or-buffer-around-text-labels-en-r-plot
shadowtext <- function(x, y=NULL, labels, col='white', bg='black', 
                       theta= seq(0, 2*pi, length.out=50), r=0.1, ... ) {

    xy <- xy.coords(x,y)
    xo <- r*strwidth('A')
    yo <- r*strheight('A')

    # draw background text with small shift in x and y in background colour
    for (i in theta) {
        text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
    }
    # draw actual text in exact xy position in foreground colour
    text(xy$x, xy$y, labels, col=col, ... )
}

# Function to compute log-loss c/o
# https://www.r-bloggers.com/making-sense-of-logarithmic-loss/
#
# Used for debulking analysis where the outcome is binary and
# predictions are probabilities

LogLossBinary = function(actual, predicted, eps = 1e-15) {
   if(is.factor(actual)){
	actual <- as.numeric(actual) - 1
   }
   predicted = pmin(pmax(predicted, eps), 1-eps)
   - (sum(actual * log(predicted) + (1 - actual) * log(1 - predicted))) / length(actual)
}