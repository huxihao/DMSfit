library(methods) ## a bug in Rscript

## Functions from the mseq package
CVTrainIndexGene <- function(data, fold, seed)
{
	uniq_cat <- unique(data$index)
	train_index <- matrix(TRUE, fold, dim(data)[1])
	
	set.seed(seed)
	cv_cat <- sample(uniq_cat)
	for (j in 1 : length(uniq_cat))
	{
		row_num <- ceiling(j / floor(length(uniq_cat) / fold))
		if (row_num > fold)
		{
			row_num <- fold
		}
		train_index[row_num, data$index == cv_cat[j]] <- FALSE
	}
	
	return(train_index)
}

getNullCount <- function(data)
{
	uniq_cat <- unique(data$index)
	null_count <- rep(0, length(data$count))
	for (k in uniq_cat)
	{
		null_count[data$index == k] <- sum(data$count[data$index == k]) / sum(data$index == k)
	}
	
	return(null_count)
}

setOriOffset <- function(data)
{
	uniq_cat <- unique(data$index)
	log_expr <- rep(0, length(data$index))
	
	for (k in uniq_cat)
	{
		log_expr[data$index == k] <- log(sum(data$count[data$index == k]) / sum(data$index == k))
	}

	return(log_expr)
}

## Functions replacing the ones in the mseq package
expData <- function(oriData, choose, llen, rlen) {
	cat("Expanding the surrounding sequences...\n")
	seq <- as.character(oriData$seq)
	run_i <- 0
	num_i <- sum(choose)
	sseq <- character(num_i * (llen + rlen))
	for (i in 1 : length(seq))
		if (choose[i]) {
			run_i <- run_i + 1
			sseq[((run_i - 1) * (llen + rlen) + 1) : (run_i * (llen + rlen))] <- seq[(i - llen) : (i + rlen - 1)]
		}
	sseq <- factor(sseq)
	cat("set of characters = ", levels(sseq), "\n")
	sseq <- matrix(sseq, ncol = llen + rlen, byrow = TRUE)
	oriName <- colnames(oriData)
	fname <- oriName[!(oriName %in% c('tag', 'seq', 'gene', 'count'))]
	data <- data.frame(oriData$count[choose], sseq, oriData[choose,][fname])
	cname <- character(llen + rlen)
	for (i in 1 : (llen + rlen)) {
		j <- i - 1 - llen
		if (j < 0) {
			cname[i] <- paste("pM", -j, sep = '')
		} else {
			cname[i] <- paste("p", j, sep = '')
		}
	}
	colnames(data) <- c("count", cname, fname)
	cat("number of genes =", length(unique(data$index)), "\n")
	cat("length of surrounding sequences =", dim(data)[2] - 2, "\n")
	cat("number of counts (positions) =", dim(data)[1], "\n")
	cat("total number of reads =", sum(data$count), "\n")
	cat("total number of features =", dim(data)[2]-2, "\n")
	return(data)
}

glmPred <- function(model.glm, newdata) {
	newdata$log_expr = 0
	newdata$log_expr0 = 0
	newdata$log_expr1 = 0
	newdata$log_expr0v = 0
	newdata$log_expr1v = 0
	newdata$log_expr0s = 0
	newdata$log_expr1s = 0
	return(predict(model.glm, newdata)-model.glm$coef[1])
}

getPredCount <- function(data, model.glm) {
	uniq_cat <- unique(data$index)
	pred_pref <- exp(glmPred(model.glm, data))
	pred_count <- rep(0, length(data$count))
	for (k in uniq_cat) {
		pred_count[data$index == k] <- sum(data$count[data$index == k]) / 
			sum(pred_pref[data$index == k]) * pred_pref[data$index == k]
	}
	return(pred_count)
}

getPredCountWeight <- function(data, model.glm, rho) {
	uniq_cat <- unique(data$index)
	pred_pref <- exp(glmPred(model.glm, data))
	pred_count <- rep(0, length(data$count))
	for (k in uniq_cat) {
		pred_count[data$index == k] <- 
		    sum(rho[data$index == k] * data$count[data$index == k]) / sum(rho[data$index == k]) *
		    sum(data$index == k) * pred_pref[data$index == k] / sum(pred_pref[data$index == k])
	}
	return(pred_count)
}

setOffsetWeight <- function(data, rho) {
	uniq_cat <- unique(data$index)
	log_expr <- rep(0, length(data$index))
	
	for (k in uniq_cat) {
		log_expr[data$index == k] <- log(
			sum(rho[data$index == k] * data$count[data$index == k]) /
			sum(rho[data$index == k]))
	}
	return(log_expr)
}

updateOffsetWeight <- function(data, log_pred_count, log_expr, rho) {
	uniq_cat <- unique(data$index)
	log_pref <- log_pred_count - log_expr
	
	for (k in uniq_cat) {
		log_expr[data$index == k] <- log(
			sum(rho[data$index == k] * data$count[data$index == k]) /
			sum(rho[data$index == k] * exp(log_pref[data$index == k])))
	}
	return(log_expr)
}

meanWeight <- function(data, rho) {
	mu <- rep(0, length(data$index))
	for(i in unique(data$index)) {
		mu[data$index == i] <- weighted.mean(data$count[data$index == i], rho[data$index == i])
	}
	return(mu)
}

getDev <- function(pred_count, real_count) {
	dev <- sum(real_count[real_count != 0] * log(real_count[real_count != 0] / pred_count[real_count != 0]))
	dev <- dev - sum(real_count) + sum(pred_count)
	dev <- 2 * dev
	return(dev)
}

getDevMix <- function(tau0, lambda0, tau1, lambda1, real_count) {
	lambda0 <- lambda0[real_count != 0]
	lambda1 <- lambda1[real_count != 0]
	x <- real_count[real_count != 0]
	a <- log(tau0) + (-lambda0 + x * log(lambda0))
	b <- log(tau1) + (-lambda1 + x * log(lambda1))
	a[is.na(a)] = 0
	b[is.na(b)] = 0
	dev <- rep(0, length(x))
	for(i in 1:length(x)) {
		if(abs(a[i] - b[i]) > 200) { ## approximate
			dev[i] <- max(a[i], b[i])
		}else{
			dev[i] <- a[i] + log(1 + exp(b[i]-a[i]))
		}
	}
	dev <- -2 * sum(dev + x - x * log(x))
	return(dev)
}

##-----------------------
##-- Poisson Linear Model
train.PoissonLinear <- function(data, thrd=1e-4, max_iter=200) {
	cat("Train: Poisson Linear\n")

	## Step 1
	log_expr <- setOriOffset(data)
	for_dev <- 0
	cat("Iter", "mu", "alpha", "dev", "\n", sep="\t")
	for (k in 1:max_iter) {
		## Step 2
		model.glm <- glm(count ~ . - index, data = data, family = poisson(link = "log"), offset = log_expr)

		## Step 3
		log_expr <- updateOffset(data, model.glm$fitted.values, log_expr)

		## Step 4
		now_dev <- getDev(getPredCount(data, model.glm), data$count)
		cat(k, round(mean(exp(log_expr)),3), round(model.glm$coef[1],3), now_dev, "\n", sep = "\t")
		if (k > 1 && abs((now_dev - for_dev) / for_dev) < thrd) { break }
		for_dev <- now_dev
	}
	null_dev <- getDev(getNullCount(data), data$count)
	cat("R-squared =", 1 - now_dev/null_dev,"\n")

	return(list(model.glm, dev= now_dev, dev0= null_dev))
}

cv.PoissonLinear <- function(data, fold=5, thrd=1e-4, max_iter=200, seed=2012) {
	cat("Cross validation: Poisson Linear\n")
	train_index <- CVTrainIndexGene(data, fold, seed)
	
	dev <- rep(0, fold)
	null_dev <- rep(0, fold)
	
	for (j in 1 : fold) {
		train_data <- data[train_index[j, ],]
		test_data <- data[!train_index[j, ],]
		cat("Fold ", j, "of", fold, nrow(train_data), "->", nrow(test_data), "\n")
		## Train
		outcomes <- train.PoissonLinear(train_data, thrd, max_iter)
		## Predict
		dev[j] <- getDev(getPredCount(test_data, outcomes[[1]]), test_data$count)
		null_dev[j] <- getDev(getNullCount(test_data), test_data$count)
	}
	
	cat("deviance =", dev, "\n")
	cat("null deviance =", null_dev, "\n")
	
	R_squared <- 1 - sum(dev) / sum(null_dev)
	cat("R_squared =", R_squared, "\n")
	
	return(R_squared)
}

##----------------------------
##-- Mixture of Poisson Models 
train.MixturePoisson <- function(data, thrd=1e-4, max_iter=200, rep_seed=10) {
	for(i in 1:rep_seed) {
		if(i == 1) {
			best <- train.MixturePoisson_seed(data, thrd=thrd, max_iter=max_iter, seed=1)
		}else{
		new <- train.MixturePoisson_seed(data, thrd=thrd, max_iter=max_iter, seed=i)
		if(new$r2 > best$r2)
			best <- new
		}
	}
	return(best)
}

train.MixturePoisson_seed <- function(data, thrd=1e-4, max_iter=200, seed=2012) {
	cat("Train: Mixture of Poisson with seed =", seed, "\n")

	## Step 1
	set.seed(seed)
	rho0 <- runif(nrow(data), 0.499, 0.501)
	rho1 <- 1 - rho0

	## Step 2
	tau0 <- mean(rho0)
	tau1 <- 1 - tau0
	
	mu0 <- meanWeight(data, rho0)
	mu1 <- meanWeight(data, rho1)

	for_dev <- 0
	cat("Iter", "tau0", "tau1", "mu0", "mu1", "dev", "\n", sep="\t")
	for (k in 1:max_iter) {
		## Step 3
		log_pred_prob0 <- dpois(data$count, mu0, log=TRUE)
		log_pred_prob1 <- dpois(data$count, mu1, log=TRUE)

		rho0 <- tau0 / (tau0  + tau1 * exp(log_pred_prob1 - log_pred_prob0))
		rho1 <- 1 - rho0

		## Step 4
		mu0 <- meanWeight(data, rho0)
		mu1 <- meanWeight(data, rho1)

		## Step 5
		swap <- mu0 > mu1 ## mu0 should less than mu1, otherwise swap
		if(sum(swap) > 0) {
		cat("Swap the labels of", length(unique(data$index[swap])), "genes\n")
		tmp <- rho0[swap]; rho0[swap] <- rho1[swap]; rho1[swap] <- tmp
		tmp <- mu0[swap]; mu0[swap] <- mu1[swap]; mu1[swap] <- tmp
		}

		## Step 6
		tau0 <- mean(rho0)
		tau1 <- 1 - tau0

		## Step 7
		now_dev <- getDevMix(tau0, mu0, tau1, mu1, data$count)
		cat(k, round(tau0,3), round(tau1,3), round(mean(mu0),3), round(mean(mu1),3), now_dev, "\n", sep = "\t")
		if (k > 3 && abs((now_dev - for_dev) / for_dev) < thrd) { break }
		for_dev <- now_dev
	}
	null_dev <- getDev(getNullCount(data), data$count)
	r_squared = 1 - now_dev/null_dev
	cat("R-squared =", r_squared,"\n")

	return(list(r2= r_squared, dev= now_dev, dev0= null_dev, prob= rho1, pred= tau0*mu0 + tau1*mu1))
}

cv.MixturePoisson <- function(data, fold=5, thrd=1e-4, max_iter=200, seed=2012) {
	cat("Cross validation using Mixture of Poisson ...\n")
	train_index <- CVTrainIndexGene(data, fold, seed)
	dev <- rep(0, fold)
	null_dev <- rep(0, fold)
	for (j in 1 : fold) {
		test_data <- data[!train_index[j, ],]
		cat("Fold ", j, "of", fold, "0 ->", nrow(test_data), "\n")
		## Predict
		outcomes <- train.MixturePoisson(test_data, thrd, max_iter)
		dev[j] <- outcomes$dev
		null_dev[j] <- outcomes$dev0
	}
	
	cat("deviance =", dev, "\n")
	cat("null deviance =", null_dev, "\n")
	
	R_squared <- 1 - sum(dev) / sum(null_dev)
	cat("R_squared =", R_squared, "\n")
	
	return(R_squared)
}

##------------------------------------
##-- Mixture of Poisson Linear Models
train.MixturePoissonLinear <- function(data, thrd=1e-4, max_iter=200, rep_seed=10) {
	for(i in 1:rep_seed) {
		if(i == 1){
			best <- train.MixturePoissonLinear_seed(data, thrd=thrd, max_iter=max_iter, seed=1)
		}else{
		new <- train.MixturePoissonLinear_seed(data, thrd=thrd, max_iter=max_iter, seed=i)
		if(new$r2 > best$r2)
			best <- new
		}
	}
	return(best)
}

train.MixturePoissonLinear_seed <- function(data, thrd=1e-4, max_iter=200, seed=2012) {
	cat("Train: Mixture of Poisson Linear with seed =", seed, "\n")

	## Step 1
	set.seed(seed)
	data <- data[sample(nrow(data)),]
	rho0 <- runif(nrow(data), 0.499, 0.501)
	rho1 <- 1 - rho0

	## Step 2
	tau0 <- mean(rho0)
	tau1 <- 1 - tau0

	log_expr0 <- setOffsetWeight(data, rho0)
	log_expr1 <- setOffsetWeight(data, rho1)

	for_dev <- 0
	cat("Iter", "tau0", "tau1", "mu0", "mu1", "alpha0", "alpha1", "dev", "\n", sep="\t")
	for (k in 1:max_iter) {
		## Step 3
		model.glm0 <- glm(count ~ . - index, data = data, family = poisson(link = "log"), offset = log_expr0, weights = rho0)
		model.glm1 <- glm(count ~ . - index, data = data, family = poisson(link = "log"), offset = log_expr1, weights = rho1)

		## Step 4
		log_pref0 <- log_expr0 + model.glm0$coef[1] + glmPred(model.glm0, data)
		log_pref1 <- log_expr1 + model.glm1$coef[1] + glmPred(model.glm1, data)

		log_pred_prob0 <- dpois(data$count, exp(log_pref0), log=TRUE)
		log_pred_prob1 <- dpois(data$count, exp(log_pref1), log=TRUE)

		rho0 <- tau0 / (tau0 + tau1 * exp(log_pred_prob1 - log_pred_prob0))
		rho1 <- 1 - rho0

		## Step 5
		log_expr0 <- updateOffsetWeight(data, log_pref0, log_expr0, rho0)
		log_expr1 <- updateOffsetWeight(data, log_pref1, log_expr1, rho1)

		## Step 6
		swap <- log_expr0 > log_expr1 ## log_expr0 should less than log_expr1, otherwise swap
		if(sum(swap) > 0) {
		cat("Swap the labels of", length(unique(data$index[swap])), "genes\n")
		tmp <- rho0[swap]; rho0[swap] <- rho1[swap]; rho1[swap] <- tmp
		tmp <- log_expr0[swap]; log_expr0[swap] <- log_expr1[swap]; log_expr1[swap] <- tmp
		}

		## Step 7
		tau0 <- mean(rho0)
		tau1 <- 1 - tau0

		## Step 8
		now_dev <- getDevMix(tau0, getPredCountWeight(data, model.glm0, rho0), tau1, getPredCountWeight(data, model.glm1, rho1), data$count)
		cat(k, round(tau0,3), round(tau1,3), round(mean(exp(log_expr0)),2), round(mean(exp(log_expr1)),2), round(model.glm0$coef[1],3), round(model.glm1$coef[1],3), now_dev, "\n", sep = "\t")
		if (k > 1 && abs((now_dev - for_dev) / for_dev) < thrd) { break }
		for_dev <- now_dev
	}

	null_dev <- getDev(getNullCount(data), data$count)
	r_squared <- 1 - now_dev/null_dev
	cat("R-squared =", r_squared,"\n")
	cat("Number of coef is", length(model.glm0$coef)-1, "\n")

	return(list(model.glm0, model.glm1, tau0, r2=r_squared, dev= now_dev, dev0= null_dev))
}

predict.MixturePoissonLinear <- function(model, data, thrd=1e-4, max_iter=200, seed=2012) {
	cat("Predict: Mixture of Poisson Linear\n")

	model.glm0 <- model[[1]]
	model.glm1 <- model[[2]]
	tau0 <- model[[3]]

	## Step 1
	set.seed(seed)
	rho0 <- runif(nrow(data), 0.499, 0.501)
	rho1 <- 1 - rho0

	## Step 2
	#tau0 <- mean(rho0)
	tau1 <- 1 - tau0

	log_expr0 <- setOffsetWeight(data, rho0)
	log_expr1 <- setOffsetWeight(data, rho1)

	for_dev <- 0
	cat("Iter", "tau0", "tau1", "mu0", "mu1", "alpha0", "alpha1", "dev", "\n", sep="\t")
	for (k in 1:max_iter) {
		## Step 3
		#model.glm0 <- glm(count ~ . - index, data = data, family = poisson(link = "log"), offset = log_expr0, weights = rho0)
		#model.glm1 <- glm(count ~ . - index, data = data, family = poisson(link = "log"), offset = log_expr1, weights = rho1)

		## Step 4
		log_pref0 <- log_expr0 + model.glm0$coef[1] + glmPred(model.glm0, data)
		log_pref1 <- log_expr1 + model.glm1$coef[1] + glmPred(model.glm1, data)

		log_pred_prob0 <- dpois(data$count, exp(log_pref0), log=TRUE)
		log_pred_prob1 <- dpois(data$count, exp(log_pref1), log=TRUE)

		rho0 <- tau0 / (tau0 + tau1 * exp(log_pred_prob1 - log_pred_prob0))
		rho1 <- 1 - rho0

		## Step 5
		log_expr0 <- updateOffsetWeight(data, log_pref0, log_expr0, rho0)
		log_expr1 <- updateOffsetWeight(data, log_pref1, log_expr1, rho1)

		## Step 6
		swap <- log_expr0 > log_expr1 ## log_expr0 should less than log_expr1, otherwise swap
		if(sum(swap) > 0) {
		cat("Swap the labels of", length(unique(data$index[swap])), "genes\n")
		tmp <- rho0[swap]; rho0[swap] <- rho1[swap]; rho1[swap] <- tmp
		tmp <- log_expr0[swap]; log_expr0[swap] <- log_expr1[swap]; log_expr1[swap] <- tmp
		}

		## Step 7
		#tau0 <- mean(rho0)
		#tau1 <- 1 - tau0

		## Step 8
		now_dev <- getDevMix(tau0, getPredCountWeight(data, model.glm0, rho0), tau1, getPredCountWeight(data, model.glm1, rho1), data$count)
		cat(k, round(tau0,3), round(tau1,3), round(mean(exp(log_expr0)),2), round(mean(exp(log_expr1)),2), round(model.glm0$coef[1],3), round(model.glm1$coef[1],3), now_dev, "\n", sep = "\t")
		if (k > 1 && abs((now_dev - for_dev) / for_dev) < thrd) { break }
		for_dev <- now_dev
	}
	null_dev <- getDev(getNullCount(data), data$count)
	cat("R-squared =", 1 - now_dev/null_dev,"\n")

	pred_count <- rho0 * getPredCountWeight(data, model.glm0, rho0) + rho1 * getPredCountWeight(data, model.glm1, rho1)
	return(list(dev= now_dev, dev0= null_dev, prob= rho1, pred= pred_count))
}

cv.MixturePoissonLinear <- function(data, fold=5, thrd=1e-4, max_iter=200, seed=2012) {
	cat("Cross validation: Mixture of Poisson Linear\n")
	train_index <- CVTrainIndexGene(data, fold, seed)
	
	dev <- rep(0, fold)
	null_dev <- rep(0, fold)
	
	for (j in 1 : fold) {
		train_data <- data[train_index[j, ],]
		test_data <- data[!train_index[j, ],]
		cat("Fold ", j, "of", fold, nrow(train_data), "->", nrow(test_data), "\n")
		## Train
		model <- train.MixturePoissonLinear(train_data, thrd, max_iter)
		## Predict
		outcomes <- predict.MixturePoissonLinear(model, test_data, thrd, max_iter)
		dev[j] <- outcomes$dev
		null_dev[j] <- outcomes$dev0
	}
	
	cat("deviance =", dev, "\n")
	cat("null deviance =", null_dev, "\n")
	
	R_squared <- 1 - sum(dev) / sum(null_dev)
	cat("R_squared =", R_squared, "\n")
	
	return(R_squared)
}

##-------------------------------------------
##-- Mixture of Poisson Linear Models Combine
train.MixturePoissonLinearCombine <- function(data, same_dir=TRUE, thrd=1e-4, max_iter=200, rep_seed=10) {
	for(i in 1:rep_seed) {
		if (i == 1)
			best <- train.MixturePoissonLinearCombine_seed(data, same_dir, thrd, max_iter, seed=1)
		new <- train.MixturePoissonLinearCombine_seed(data, same_dir, thrd, max_iter, seed=i)
		if(new$r2 > best$r2)
			best <- new
	}
	return(best)
}

train.MixturePoissonLinearCombine_seed <- function(data, same_dir=TRUE, thrd=1e-4, max_iter=200, seed=2012) {
	cat("Train: Mixture of Poisson Linear Combine with seed =", seed, "\n")

	## Step 1
	set.seed(seed)
	rho0 <- runif(nrow(data$v1), 0.499, 0.501)
	rho1 <- 1 - rho0

	## Step 2
	tau0 <- mean(rho0)
	tau1 <- 1 - tau0

	log_expr0v <- setOffsetWeight(data$v1, rho0)
	log_expr1v <- setOffsetWeight(data$v1, rho1)
	log_expr0s <- setOffsetWeight(data$s1, rho0)
	log_expr1s <- setOffsetWeight(data$s1, rho1)

	for_dev <- 0
	cat("Iter", "tau0", "tau1", "mu0v", "mu1v", "mu0s", "mu1s", "dev", "\n", sep="\t")
	for (k in 1 : max_iter) {
		## Step 3
		model.glm0v <- glm(count ~ .-index, data = data$v1, family = poisson(link = "log"), offset = log_expr0v, weights = rho0)
		model.glm1v <- glm(count ~ .-index, data = data$v1, family = poisson(link = "log"), offset = log_expr1v, weights = rho1)
		model.glm0s <- glm(count ~ .-index, data = data$s1, family = poisson(link = "log"), offset = log_expr0s, weights = rho0)
		model.glm1s <- glm(count ~ .-index, data = data$s1, family = poisson(link = "log"), offset = log_expr1s, weights = rho1)

		## Step 4
		log_pref0v <- log_expr0v + model.glm0v$coef[1] + glmPred(model.glm0v, data$v1)
		log_pref1v <- log_expr1v + model.glm1v$coef[1] + glmPred(model.glm1v, data$v1)
		log_pref0s <- log_expr0s + model.glm0s$coef[1] + glmPred(model.glm0s, data$s1)
		log_pref1s <- log_expr1s + model.glm1s$coef[1] + glmPred(model.glm1s, data$s1)

		log_pred_prob0v <- dpois(data$v1$count, exp(log_pref0v), log=TRUE)
		log_pred_prob1v <- dpois(data$v1$count, exp(log_pref1v), log=TRUE)
		log_pred_prob0s <- dpois(data$s1$count, exp(log_pref0s), log=TRUE)
		log_pred_prob1s <- dpois(data$s1$count, exp(log_pref1s), log=TRUE)

		rho0v <- 1 / (1 + (tau1/tau0) * exp(log_pred_prob1v - log_pred_prob0v))
		rho0s <- 1 / (1 + (tau1/tau0) * exp(log_pred_prob1s - log_pred_prob0s))
		rho1v <- 1 - rho0v
		rho1s <- 1 - rho0s

		## Step 5
		log_expr0v <- updateOffsetWeight(data$v1, log_pref0v, log_expr0v, rho0)
		log_expr1v <- updateOffsetWeight(data$v1, log_pref1v, log_expr1v, rho1)
		log_expr0s <- updateOffsetWeight(data$s1, log_pref0s, log_expr0s, rho0)
		log_expr1s <- updateOffsetWeight(data$s1, log_pref1s, log_expr1s, rho1)

		## Step 6
		swap <- log_expr0v > log_expr1v
		if(sum(swap) > 0) {
		cat("Swap the labels of", length(unique(data$v1$index[swap])), "genes in V1\n")
		tmp <- rho0v[swap]; rho0v[swap] <- rho1v[swap]; rho1v[swap] <- tmp
		}
		if(same_dir) {
		swap <- log_expr0s > log_expr1s
		}else{ ## opposite direction
		swap <- log_expr0s < log_expr1s
		}
		if(sum(swap) > 0) {
		cat("Swap the labels of", length(unique(data$s1$index[swap])), "genes in S1\n")
		tmp <- rho0s[swap]; rho0s[swap] <- rho1s[swap]; rho1s[swap] <- tmp
		}
		
		## Step 7
		rho0 <- (rho0v+rho0s)/2
		rho1 <- 1 - rho0

		## Step 8
		tau0 <- mean(rho0)
		tau1 <- 1 - tau0

		log_expr0v <- updateOffsetWeight(data$v1, log_pref0v, log_expr0v, rho0)
		log_expr1v <- updateOffsetWeight(data$v1, log_pref1v, log_expr1v, rho1)
		log_expr0s <- updateOffsetWeight(data$s1, log_pref0s, log_expr0s, rho0)
		log_expr1s <- updateOffsetWeight(data$s1, log_pref1s, log_expr1s, rho1)

		## Step 9
		now_dev <- 
		  getDevMix(tau0, getPredCountWeight(data$v1, model.glm0v, rho0), tau1, getPredCountWeight(data$v1, model.glm1v, rho1), data$v1$count) +
		  getDevMix(tau1, getPredCountWeight(data$s1, model.glm0s, rho0), tau0, getPredCountWeight(data$s1, model.glm1s, rho1), data$s1$count)
		cat(k, round(tau0,3), round(tau1,3), round(mean(exp(log_expr0v))), round(mean(exp(log_expr1v))), round(mean(exp(log_expr0s))), round(mean(exp(log_expr1s))), now_dev, "\n", sep = "\t")
		if (k > 3 && (for_dev - now_dev) / for_dev < thrd) { break }
		for_dev <- now_dev
	}
	null_dev <- getDev(getNullCount(data$v1), data$v1$count) + getDev(getNullCount(data$s1), data$s1$count)
	r_squared <- 1 - now_dev/null_dev
	cat("R-squared =", r_squared,"\n")
	return(list(model.glm0v, model.glm1v, model.glm0s, model.glm1s, tau0, r2= r_squared, dev= now_dev, dev0= null_dev))
}

predict.MixturePoissonLinearCombine <- function(model, data, same_dir=TRUE, thrd=1e-4, max_iter=200, seed=2012) {
	cat("Predict: Mixture of Poisson Linear Combine\n")
	model.glm0v <- model[[1]]
	model.glm1v <- model[[2]]
	model.glm0s <- model[[3]]
	model.glm1s <- model[[4]]
	tau0 <- model[[5]]

	## Step 1
	set.seed(seed)
	rho0 <- runif(nrow(data$v1), 0.499, 0.501)
	rho1 <- 1 - rho0

	## Step 2
	#tau0 <- mean(rho0)
	tau1 <- 1 - tau0

	log_expr0v <- setOffsetWeight(data$v1, rho0)
	log_expr1v <- setOffsetWeight(data$v1, rho1)
	log_expr0s <- setOffsetWeight(data$s1, rho0)
	log_expr1s <- setOffsetWeight(data$s1, rho1)

	for_dev <- 0
	cat("Iter", "tau0", "tau1", "mu0v", "mu1v", "mu0s", "mu1s", "dev", "\n", sep="\t")
	for (k in 1 : max_iter) {
		## Step 3
		#model.glm0v <- glm(count ~ .-index, data = data$v1, family = poisson(link = "log"), offset = log_expr0v, weights = rho0)
		#model.glm1v <- glm(count ~ .-index, data = data$v1, family = poisson(link = "log"), offset = log_expr1v, weights = rho1)
		#model.glm0s <- glm(count ~ .-index, data = data$s1, family = poisson(link = "log"), offset = log_expr0s, weights = rho0)
		#model.glm1s <- glm(count ~ .-index, data = data$s1, family = poisson(link = "log"), offset = log_expr1s, weights = rho1)

		## Step 4
		log_pref0v <- log_expr0v + model.glm0v$coef[1] + glmPred(model.glm0v, data$v1)
		log_pref1v <- log_expr1v + model.glm1v$coef[1] + glmPred(model.glm1v, data$v1)
		log_pref0s <- log_expr0s + model.glm0s$coef[1] + glmPred(model.glm0s, data$s1)
		log_pref1s <- log_expr1s + model.glm1s$coef[1] + glmPred(model.glm1s, data$s1)

		log_pred_prob0v <- dpois(data$v1$count, exp(log_pref0v), log=TRUE)
		log_pred_prob1v <- dpois(data$v1$count, exp(log_pref1v), log=TRUE)
		log_pred_prob0s <- dpois(data$s1$count, exp(log_pref0s), log=TRUE)
		log_pred_prob1s <- dpois(data$s1$count, exp(log_pref1s), log=TRUE)

		rho0v <- 1 / (1 + (tau1/tau0) * exp(log_pred_prob1v - log_pred_prob0v))
		rho0s <- 1 / (1 + (tau1/tau0) * exp(log_pred_prob1s - log_pred_prob0s))
		rho1v <- 1 - rho0v
		rho1s <- 1 - rho0s

		## Step 5
		log_expr0v <- updateOffsetWeight(data$v1, log_pref0v, log_expr0v, rho0)
		log_expr1v <- updateOffsetWeight(data$v1, log_pref1v, log_expr1v, rho1)
		log_expr0s <- updateOffsetWeight(data$s1, log_pref0s, log_expr0s, rho0)
		log_expr1s <- updateOffsetWeight(data$s1, log_pref1s, log_expr1s, rho1)

		## Step 6
		swap <- log_expr0v > log_expr1v
		if(sum(swap) > 0) {
		cat("Swap the labels of", length(unique(data$v1$index[swap])), "genes in V1\n")
		tmp <- rho0v[swap]; rho0v[swap] <- rho1v[swap]; rho1v[swap] <- tmp
		}
		if(same_dir) {
		swap <- log_expr0s > log_expr1s
		}else{ ## opposite direction
		swap <- log_expr0s < log_expr1s
		}
		if(sum(swap) > 0) {
		cat("Swap the labels of", length(unique(data$s1$index[swap])), "genes in S1\n")
		tmp <- rho0s[swap]; rho0s[swap] <- rho1s[swap]; rho1s[swap] <- tmp
		}
		
		## Step 7
		rho0 <- (rho0v+rho0s)/2
		rho1 <- 1 - rho0

		## Step 8
		#tau0 <- mean(rho0)
		#tau1 <- 1 - tau0

		log_expr0v <- updateOffsetWeight(data$v1, log_pref0v, log_expr0v, rho0)
		log_expr1v <- updateOffsetWeight(data$v1, log_pref1v, log_expr1v, rho1)
		log_expr0s <- updateOffsetWeight(data$s1, log_pref0s, log_expr0s, rho0)
		log_expr1s <- updateOffsetWeight(data$s1, log_pref1s, log_expr1s, rho1)

		## Step 9
		now_dev_v <- getDevMix(tau0, getPredCountWeight(data$v1, model.glm0v, rho0), tau1, getPredCountWeight(data$v1, model.glm1v, rho1), data$v1$count)
		now_dev_s <- getDevMix(tau1, getPredCountWeight(data$s1, model.glm0s, rho0), tau0, getPredCountWeight(data$s1, model.glm1s, rho1), data$s1$count)
		cat(k, round(tau0,3), round(tau1,3), round(mean(exp(log_expr0v))), round(mean(exp(log_expr1v))), round(mean(exp(log_expr0s))), round(mean(exp(log_expr1s))), now_dev_v + now_dev_s, "\n", sep = "\t")
		if (k > 3 && (for_dev - now_dev_v - now_dev_s) / for_dev < thrd) { break }
		for_dev <- now_dev_v + now_dev_s
	}
	null_dev_v <- getDev(getNullCount(data$v1), data$v1$count)
	null_dev_s <- getDev(getNullCount(data$s1), data$s1$count)
	cat("R-squared =", 1 - (now_dev_v + now_dev_s) / (null_dev_v + null_dev_s),"\n")

	predCountV <- rho0 * getPredCountWeight(data$v1, model.glm0v, rho0) + rho1 * getPredCountWeight(data$v1, model.glm1v, rho1)
	predCountS <- rho0 * getPredCountWeight(data$s1, model.glm0s, rho0) + rho1 * getPredCountWeight(data$s1, model.glm1s, rho1)
	return(list(ProbV1=rho1v, ProbS1=rho1s,
	            PredV1=predCountV, PredS1=predCountS,
				DevV1=now_dev_v, DevS1=now_dev_s,
				Dev0V1=null_dev_v, Dev0S1=null_dev_s))
}

cv.MixturePoissonLinearCombine <- function(data, fold=5, same_dir=TRUE, seed=2011, thrd=1e-4, max_iter=200) {
	cat("Cross validation: Mixture of Poisson Linear Combine\n")
	train_index <- CVTrainIndexGene(data$v1, fold, seed)
	dev_v1 <- rep(0, fold); dev_s1 <- rep(0, fold);
	null_dev_v1 <- rep(0, fold); null_dev_s1 <- rep(0, fold);
	for (j in 1 : fold) {
		train_data <- list(v1= data$v1[train_index[j, ],], s1= data$s1[train_index[j, ],])
		test_data <- list(v1= data$v1[!train_index[j, ],], s1= data$s1[!train_index[j, ],])
		cat("Fold ", j, "of", fold, nrow(train_data$v1), "->", nrow(test_data$v1), "\n")
		## Train
		model <- train.MixturePoissonLinearCombine(train_data, same_dir, thrd, max_iter)
		## Predict
		outcomes <- predict.MixturePoissonLinearCombine(model, test_data, same_dir, thrd, max_iter)
		dev_v1[j] = outcomes$DevV1
		dev_s1[j] = outcomes$DevS1
		null_dev_v1[j] = outcomes$Dev0V1
		null_dev_s1[j] = outcomes$Dev0S1
	}
	return(list(r2v = 1-sum(dev_v1)/sum(null_dev_v1),
	            r2s = 1-sum(dev_s1)/sum(null_dev_s1)))
}


##############################################################
## fit seperately
train.MixturePoissonSeperate <- function(data, thrd=1e-4, max_iter=200) {
	v1_fitted <- train.MixturePoisson(data$v1, thrd, max_iter)
	s1_fitted <- train.MixturePoisson(data$s1, thrd, max_iter)
	return(list(ProbV1=v1_fitted$prob, ProbS1=s1_fitted$prob, 
				PredV1=v1_fitted$pred, PredS1=s1_fitted$pred, 
				DevV1=v1_fitted$dev, DevS1=s1_fitted$dev, 
				Dev0V1=v1_fitted$dev0, Dev0S1=s1_fitted$dev0))
}

train.MixturePoissonLinearSeperate <- function(data, thrd=1e-4, max_iter=200) {
	v1_model <- train.MixturePoissonLinear(data$v1, thrd, max_iter)
	s1_model <- train.MixturePoissonLinear(data$s1, thrd, max_iter)
	return(list(v1_model, s1_model))
}

predict.MixturePoissonLinearSeperate <- function(model, data, thrd=1e-4, max_iter=200) {
	v1_fitted <- predict.MixturePoissonLinear(model[[1]], data$v1, thrd, max_iter)
	s1_fitted <- predict.MixturePoissonLinear(model[[2]], data$s1, thrd, max_iter)
	return(list(ProbV1=v1_fitted$prob, ProbS1=s1_fitted$prob, 
	            PredV1=v1_fitted$pred, PredS1=s1_fitted$pred, 
				DevV1=v1_fitted$dev, DevS1=s1_fitted$dev,
				Dev0V1=v1_fitted$dev0, Dev0S1=s1_fitted$dev0))
}

##############################################################
## treat V1 as the control and fit S1
train.MixturePoissonLinearStepwise <- function(data, case=1, thrd=1e-4, max_iter=200) {
	if(case == 0) { ## s1 ~ seq + v1
		v1_model <- train.MixturePoissonLinear(data$s1, thrd, max_iter)
		data_s1 <- data$s1
		data_s1$pbv <- data$v1$count
	}
	if(case == 1) { ## s1 ~ seq + log(v1)
		v1_model <- train.MixturePoissonLinear(data$v1, thrd, max_iter)
		data_s1 <- data$s1
		data_s1$pbv <- log2(data$v1$count+1)
	}
	if(case == 2) { ## s1 ~ seq + prob(v1)
		v1_model <- train.MixturePoissonLinear(data$v1, thrd, max_iter)
		v1_fitted <- predict.MixturePoissonLinear(v1_model, data$v1, thrd, max_iter)
		data_s1 <- data$s1
		data_s1$pbv <- v1_fitted$prob
	}
	if(case == 3) { ## s1 ~ seq + prob(s1)
		v1_model <- train.MixturePoissonLinear(data$s1, thrd, max_iter)
		v1_fitted <- predict.MixturePoissonLinear(v1_model, data$s1, thrd, max_iter)
		data_s1 <- data$s1
		data_s1$pbv <- v1_fitted$prob
	}
	s1_model <- train.MixturePoissonLinear(data_s1, thrd, max_iter)
	return(list(v1_model, s1_model))
}

predict.MixturePoissonLinearStepwise <- function(model, data, case=1, thrd=1e-4, max_iter=200) {
	if(case == 0) { ## s1 ~ seq + v1
		v1_fitted <- predict.MixturePoissonLinear(model[[1]], data$s1, thrd, max_iter)
		data_s1 <- data$s1
		data_s1$pbv <- data$v1$count
	}
	if(case == 1) { ## s1 ~ seq + log(v1)
		v1_fitted <- predict.MixturePoissonLinear(model[[1]], data$v1, thrd, max_iter)
		data_s1 <- data$s1
		data_s1$pbv <- log2(data$v1$count+1)
	}
	if(case == 2) { ## s1 ~ seq + prob(v1)
		v1_fitted <- predict.MixturePoissonLinear(model[[1]], data$v1, thrd, max_iter)
		data_s1 <- data$s1
		data_s1$pbv <- v1_fitted$prob
	}
	if(case == 3) { ## s1 ~ seq + prob(s1)
		v1_fitted <- predict.MixturePoissonLinear(model[[1]], data$s1, thrd, max_iter)
		data_s1 <- data$s1
		data_s1$pbv <- v1_fitted$prob
	}
	s1_fitted <- predict.MixturePoissonLinear(model[[2]], data_s1, thrd, max_iter)
	return(list(ProbV1=v1_fitted$prob, ProbS1=s1_fitted$prob, 
	            PredV1=v1_fitted$pred, PredS1=s1_fitted$pred, 
				DevV1=v1_fitted$dev, DevS1=s1_fitted$dev,
				Dev0V1=v1_fitted$dev0, Dev0S1=s1_fitted$dev0))
}

cv.MixturePoissonLinearStepwise <- function(data, case=1, fold=5, seed=2011, thrd=1e-4, max_iter=200) {
	cat("Cross validation: Mixture of Poisson Linear Stepwise\n")
	train_index <- CVTrainIndexGene(data$v1, fold, seed)
	dev_v1 <- rep(0, fold); dev_s1 <- rep(0, fold);
	null_dev_v1 <- rep(0, fold); null_dev_s1 <- rep(0, fold);
	for (j in 1 : fold) {
		train_data <- list(v1= data$v1[train_index[j, ],], s1= data$s1[train_index[j, ],])
		test_data <- list(v1= data$v1[!train_index[j, ],], s1= data$s1[!train_index[j, ],])
		cat("Fold ", j, "of", fold, nrow(train_data$v1), "->", nrow(test_data$v1), "\n")
		## Train
		model <- train.MixturePoissonLinearStepwise(train_data, case, thrd, max_iter)
		## Predict
		outcomes <- predict.MixturePoissonLinearStepwise(model, test_data, case, thrd, max_iter)
		dev_v1[j] = outcomes$DevV1
		dev_s1[j] = outcomes$DevS1
		null_dev_v1[j] = outcomes$Dev0V1
		null_dev_s1[j] = outcomes$Dev0S1
	}
	return(list(r2v = 1-sum(dev_v1)/sum(null_dev_v1),
	            r2s = 1-sum(dev_s1)/sum(null_dev_s1)))
}

## format data
expDataCombine <- function(data_table, choose, left_span, right_span) {
	data_table$count <- data_table$v1
	encoded_v <- expData(data_table, choose, left_span, right_span)
	data_table$count <- data_table$s1
	encoded_s <- expData(data_table, choose, left_span, right_span)
	return(list(v1=encoded_v, s1=encoded_s))
}

## generate a new random data set
generateData <- function(choose, index, count, tau=0.7, relation=TRUE, seed=2011) {
	newCount <- rep(0, length(count)) ## a new copy
	for (i in unique(index)) {
		geneIndex <- choose & index == i
		counts <- count[geneIndex]
		num0 <- round(sum(geneIndex) * tau)
		num1 <- sum(geneIndex) - num0
		mu0 <- mean(counts) * 0.8
		mu1 <- mean(counts) * 1.2
		set.seed(seed)
		randomIndex <- sample(c(1:sum(geneIndex)))
		if (relation) { ## positive
			mixCount <- c(rpois(num0, mu0), rpois(num1, mu1))
		}else{
			mixCount <- c(rpois(num0, mu1), rpois(num1, mu0))
		}
		newCount[geneIndex] <- mixCount[randomIndex]# + round(rnorm(length(mixCount), 10, 3))
	}
	newCount[newCount < 0] <- 0 ## trim
	return(newCount)
}

fit_pars <- function(action, combine, window, data_file) {
	cat(action, combine, window, data_file, "\n")
	data_table <- read.csv(data_file)
	log_file <- paste(data_file, ".log", sep="")
	
	## generate data set
	#data_table$v1 <- generateData(data_table$tag>=0, data_table$index, data_table$v1, 0.5, TRUE)
	#data_table$s1 <- generateData(data_table$tag>=0, data_table$index, data_table$s1, 0.5, TRUE) 

	########################################################################
	left_span <- window/2
	right_span <- window/2

	########################################################################
	## Train
	if (action == "Train")
	{
		if (combine == "MPL") {
			encoded <- expData(data_table, data_table$tag >= 0, left_span, right_span)
			model <- train.MixturePoissonLinear(encoded)
			write(paste(names(model[[1]]$coef), sep=""), file=log_file, append=FALSE)
			write(paste(model[[1]]$coef, sep=""), file=log_file, append=TRUE)
			write(paste(model[[2]]$coef, sep=""), file=log_file, append=TRUE)
		}else{
		encoded <- expDataCombine(data_table, data_table$tag >= 0, left_span, right_span)
		if (combine == "MixPoiSep")
			model <- train.MixturePoissonSeperate(encoded)
		if (combine == "MixPoiLin") 
			model <- train.MixturePoissonLinearSeperate(encoded)
		if (combine == "MPLStepwise0")
			model <- train.MixturePoissonLinearStepwise(encoded, 0)
		if (combine == "MPLStepwise1")
			model <- train.MixturePoissonLinearStepwise(encoded, 1)
		if (combine == "MPLStepwise2")
			model <- train.MixturePoissonLinearStepwise(encoded, 2)
		if (combine == "MPLStepwise3")
			model <- train.MixturePoissonLinearStepwise(encoded, 3)
		if (combine == "MPLComSam")
			model <- train.MixturePoissonLinearCombine(encoded, TRUE)
		if (combine == "MPLComOpp")
			model <- train.MixturePoissonLinearCombine(encoded, FALSE)
		}

		## save models
		save(model, file=paste(data_file, ".model", sep=""))
	}

	########################################################################
	## Predict
	if (action == "Predict")
	{
		## setting default values for blind regions
		onlyv <- FALSE
		numRow <- nrow(data_table)
		probv <- rep(0.5, numRow)
		probs <- rep(0.5, numRow)
		predv <- rep(-1, numRow)
		preds <- rep(-1, numRow)
		load(file=paste(data_file, ".model", sep=""))
		pack_ids <- round(data_table$index, -2) # package size: default -2
		r2v_dev = 0; r2v_dev0 = 0;
		r2s_dev = 0; r2s_dev0 = 0;

		for(i in unique(pack_ids)) {
			cat("predict package", i, "...\n")
			choose <- pack_ids == i & data_table$tag >= 0
			if (combine == "MPL") {
				encoded <- expData(data_table, choose, left_span, right_span)
				outcome <- predict.MixturePoissonLinear(model, encoded)
				probv[choose] <- outcome$prob
				predv[choose] <- outcome$pred
				r2v_dev <- r2v_dev + outcome$dev
				r2v_dev0 <- r2v_dev0 + outcome$dev0
				onlyv <- TRUE
			}else{
			encoded <- expDataCombine(data_table, choose, left_span, right_span)
			if (combine == "MixPoiSep")
				outcome <- train.MixturePoissonSeperate(encoded)
			if (combine == "MixPoiLin") 
				outcome <- predict.MixturePoissonLinearSeperate(model, encoded)
			if (combine == "MPLStepwise0")
				outcome <- predict.MixturePoissonLinearStepwise(model, encoded, 0)
			if (combine == "MPLStepwise1")
				outcome <- predict.MixturePoissonLinearStepwise(model, encoded, 1)
			if (combine == "MPLStepwise2")
				outcome <- predict.MixturePoissonLinearStepwise(model, encoded, 2)
			if (combine == "MPLStepwise3")
				outcome <- predict.MixturePoissonLinearStepwise(model, encoded, 3)
			if (combine == "MPLComSam")
				outcome <- predict.MixturePoissonLinearCombine(model, encoded, TRUE)
			if (combine == "MPLComOpp")
				outcome <- predict.MixturePoissonLinearCombine(model, encoded, FALSE)
			probv[choose] <- outcome$ProbV1
			probs[choose] <- outcome$ProbS1
			predv[choose] <- outcome$PredV1
			preds[choose] <- outcome$PredS1
			r2v_dev <- r2v_dev + outcome$DevV1
			r2s_dev <- r2s_dev + outcome$DevS1
			r2v_dev0 <- r2v_dev0 + outcome$Dev0V1
			r2s_dev0 <- r2s_dev0 + outcome$Dev0S1
			}
		}
		## Save R-squared for V1 and S1
		r2v = 1 - r2v_dev/r2v_dev0
		r2s = 1 - r2s_dev/r2s_dev0

		if (onlyv) {
		data_table$pred <- predv
		data_table$prob <- probv
		write.csv(data_table[,!(names(data_table) %in% c("pred"))], file=data_file, row.names=FALSE, quote=FALSE)
		write(paste(data_file, action, combine, window, r2v, r2v_dev, sep="\t"), file=log_file)
		}else{
		data_table$pdv <- predv
		data_table$pds <- preds
		data_table$pbv <- probv
		data_table$pbs <- probs
		write.csv(data_table[,!(names(data_table) %in% c("pdv","pds"))], file=data_file, row.names=FALSE, quote=FALSE)
		write(paste(data_file, action, combine, window, r2v, r2s, r2v_dev, r2s_dev, sep="\t"), file=log_file)
		}
	}

	########################################################################
	## Cross Validation
	if (action == "CV") 
	{
		if (combine == "MPL") {
			encoded <- expData(data_table, data_table$tag >= 0, left_span, right_span)
			r2 = cv.MixturePoissonLinear(encoded)
			write(paste(data_file, action, combine, window, r2, sep="\t"), file=log_file)
		}else{
		r2v = 0; r2s = 0;
		encoded <- expDataCombine(data_table, data_table$tag >= 0, left_span, right_span)
		if (combine == "MixPoiSep") {
			r2v = cv.MixturePoisson(encoded$v1)
			r2s = cv.MixturePoisson(encoded$s1)
		}
		if (combine == "MixPoiLin") {
			r2v = cv.MixturePoissonLinear(encoded$v1)
			r2s = cv.MixturePoissonLinear(encoded$s1)
		}
		if (combine == "MPLStepwise0") {
			tmp = cv.MixturePoissonLinearStepwise(encoded, 0)
			r2v = tmp$r2v
			r2s = tmp$r2s
		}
		if (combine == "MPLStepwise1") {
			tmp = cv.MixturePoissonLinearStepwise(encoded, 1)
			r2v = tmp$r2v
			r2s = tmp$r2s
		}
		if (combine == "MPLStepwise2") {
			tmp = cv.MixturePoissonLinearStepwise(encoded, 2)
			r2v = tmp$r2v
			r2s = tmp$r2s
		}
		if (combine == "MPLStepwise3") {
			tmp = cv.MixturePoissonLinearStepwise(encoded, 3)
			r2v = tmp$r2v
			r2s = tmp$r2s
		}
		if (combine == "MPLComSam") {
			tmp = cv.MixturePoissonLinearCombine(encoded, same_dir=TRUE)
			r2v = tmp$r2v
			r2s = tmp$r2s
		}
		if (combine == "MPLComOpp") {
			tmp = cv.MixturePoissonLinearCombine(encoded, same_dir=FALSE)
			r2v = tmp$r2v
			r2s = tmp$r2s
		}
		write(paste(data_file, action, combine, window, r2v, r2s, sep="\t"), file=log_file)
		}
	}
	########################################################################
	## End of all actions
}

## The Main Function
main <- function(self.name) {
	initial.options <- commandArgs(trailingOnly = FALSE)
	file.arg.name <- "--file="
	script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
	script.basename <- basename(script.name)
	if(script.basename != self.name){return(0)}

	action <- "Train"
	combine <- "MPL"
	window <- 2
	data_file <- "../work/ProbRNA.csv"
	args <- commandArgs(TRUE)
	if(length(args) >= 1) {
		action <- args[1]	
	if(length(args) >= 2) {
		combine <- args[2]
	if(length(args) >= 3) {
		window <- as.numeric(args[3])
	if(length(args) >= 4) {
		data_file <- args[4]
	} }}}
	fit_pars(action, combine, window, data_file)
}

main('fit_po1.R')
