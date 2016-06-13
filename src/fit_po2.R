library(methods) ## a bug in Rscript
library(MASS)
source("../src/fit_po1.R")

## Functions replacing the ones in the mseq package
expData <- function(oriData, choose, llen, rlen) {
	cat("Expanding the surrounding sequences...\n")
	seq <- as.character(oriData$seq)
	run_i <- 0
	num_i <- sum(choose)
	index <- numeric(num_i)
	sseq <- character(num_i * (llen + rlen))
	mcount <- numeric(num_i)
	for (i in 1 : length(seq)) {
		if (choose[i]) {
			run_i <- run_i + 1
			sseq[((run_i - 1) * (llen + rlen) + 1) : (run_i * (llen + rlen))] <- seq[(i - llen) : (i + rlen - 1)]
			idxrng <- (1:length(seq))[oriData$index==oriData$index[i]]
			mcount[run_i] <- log(1+max(oriData$count[max(min(idxrng), i-25) : min(max(idxrng), i+25-1)]))
		}
	}
	sseq <- factor(sseq)
	cat("set of characters = ", levels(sseq), "\n")
	sseq <- matrix(sseq, ncol = llen + rlen, byrow = TRUE)
	oriName <- colnames(oriData)
	fname <- oriName[!(oriName %in% c('count', 'tag', 'seq', 'gene'))]
	data <- data.frame(oriData$count[choose], sseq, oriData[choose,][fname], mcount)
	cname <- character(llen + rlen)
	for (i in 1 : (llen + rlen)) {
		j <- i - 1 - llen
		if (j < 0) {
			cname[i] <- paste("pM", -j, sep = '')
		} else {
			cname[i] <- paste("p", j, sep = '')
		}
	}
	colnames(data) <- c("count", cname, fname, "max_count")
	cat(names(data))
	cat("number of genes =", length(unique(data$index)), "\n")
	cat("length of surrounding sequences =", dim(data)[2] - 2, "\n")
	cat("number of counts (positions) =", dim(data)[1], "\n")
	cat("total number of reads =", sum(data$count), "\n")
	cat("total number of features =", dim(data)[2]-2, "\n")
	return(data)
}

main('fit_po2.R')
