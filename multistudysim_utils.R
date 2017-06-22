# Load packages

library(RcppArmadillo)
library(e1071)
library(ranger)
library(glmnet)
library(rpart)
library(genefilter)
library(nnet)
library(nnls)

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
	#tout <- colttests(as.matrix(data[,-1]), data[,1], tstatOnly = TRUE)
	#data <- data[, c("y",rownames(tout)[order(abs(tout[,1]), decreasing = T)[1:20]])]
	mod <- rpart::rpart(y ~ ., data = as.data.frame(data), ...)
	mod <- rpart::prune(mod, mod$cptable[which(abs(diff(mod$cptable[,"xerror"])) < 0.01)[1], "CP"])
	mod
}

treefit_ensemble <- function(data, n_idx, p_idx, ...){
	form <- paste0("y ~ ", paste0("V", p_idx, collapse = " + "))
	mod <- rpart::rpart(form, data = as.data.frame(data), subset = n_idx)
	#mod <- rpart::prune(mod, mod$cptable[which(abs(diff(mod$cptable[,"xerror"])) < 0.01)[1], "CP"])
	mod
}

treepred <- function(mod, newdata){
	as.vector(predict(mod, newdata=newdata))
	#as.vector(predict(mod, newdata=newdata)[,2])
}

lassofit <- function(data, ...){
	mod <- glmnet::cv.glmnet(x = as.matrix(data[,-1]), y = data[,1], ...)
}

lassopred <- function(mod, newdata){
	as.vector(predict(mod, newx=as.matrix(newdata), s="lambda.1se"))
}

rffit <- function(data, ...){
	#tout <- colttests(as.matrix(data[,-1]), data[,1], tstatOnly = TRUE)
	#data <- data[, c("y",rownames(tout)[order(abs(tout[,1]), decreasing = T)[1:20]])]
	#ranger(y ~ ., data = data, write.forest = TRUE, probability = T) # added probability for tumorgrade stuff
	ranger::ranger(y ~ ., data = data, write.forest = TRUE)
}

rfpred <- function(mod, newdata){
	#predict(mod, data = newdata)$predictions[,2] # second column for probability stuff
	predict(mod, data = newdata)$predictions
}

nnetfit <- function(data, ...){
#	tout <- colttests(as.matrix(data[,-1]), data[,1], tstatOnly = TRUE)
#	data <- data[, c("y",rownames(tout)[order(abs(tout[,1]), decreasing = T)[1:50]])]
	nnet::nnet(y ~ ., data = data, size = 10, linout = T, trace = F)
}

nnetpred <- function(mod, newdata){
	as.vector(predict(mod, newdata = newdata))
}

# Function to normalize by subtracting off the worst performer
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



# credit to Greg Snow
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


