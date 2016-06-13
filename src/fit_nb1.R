library(methods) ## a bug in Rscript
library(MASS)
source("../src/fit_po1.R")

##-----------------------
##-- Negative Binomial Model
train.PoissonLinear <- function(data, thrd=1e-5, max_iter=200) {
	cat("Train: Negative Binomial\n")

	## Step 1
	log_expr <- setOriOffset(data)
	for_dev <- 0
	cat("Iter", "mu", "alpha", "dev", "\n", sep="\t")
	for (k in 1:max_iter) {
		## Step 2
		model.glm <- glm.nb(count ~ . - index + offset(log_expr), data = data)

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

train.MixturePoissonLinear <- function(data, thrd=1e-5, max_iter=200, rep_seed=10) {
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

train.MixturePoissonLinear_seed <- function(data, thrd=1e-5, max_iter=200, seed=2012) {
	cat("Train: Mixture of Negative Binomial with seed =", seed, "\n")

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
		model.glm0 <- glm.nb(count ~ . - index + offset(log_expr0), data = data, weights = rho0)
		model.glm1 <- glm.nb(count ~ . - index + offset(log_expr1), data = data, weights = rho1)

		## Step 4
		log_pref0 <- log_expr0 + model.glm0$coef[1] + glmPred(model.glm0, data)
		log_pref1 <- log_expr1 + model.glm1$coef[1] + glmPred(model.glm1, data)

		log_pred_prob0 <- dnbinom(data$count, size=model.glm0$theta, mu=exp(log_pref0), log=TRUE)
		log_pred_prob1 <- dnbinom(data$count, size=model.glm1$theta, mu=exp(log_pref1), log=TRUE)

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

predict.MixturePoissonLinear <- function(model, data, thrd=1e-5, max_iter=200, seed=2012) {
	cat("Predict: Mixture of Negative Binomial\n")

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
		#model.glm0 <- glm.nb(count ~ . - index + offset(log_expr0), data = data, weights = rho0)
		#model.glm1 <- glm.nb(count ~ . - index + offset(log_expr1), data = data, weights = rho1)

		## Step 4
		log_pref0 <- log_expr0 + model.glm0$coef[1] + glmPred(model.glm0, data)
		log_pref1 <- log_expr1 + model.glm1$coef[1] + glmPred(model.glm1, data)

		log_pred_prob0 <- dnbinom(data$count, size=model.glm0$theta, mu=exp(log_pref0), log=TRUE)
		log_pred_prob1 <- dnbinom(data$count, size=model.glm1$theta, mu=exp(log_pref1), log=TRUE)

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

main('fit_nb1.R')
