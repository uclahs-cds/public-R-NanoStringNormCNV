
apply.kd.cna.thresh <- function(output, headers, kd.thresh){
	which.cna <- colnames(output)[!colnames(output) %in% headers];
	which.n   <- which(colnames(output) %in% which.cna);

	# pull below into a separate object
	# pull below into a separate object
	na.counts <- apply(
		X = output[,which.n, drop = FALSE],
		MARGIN = 2,
		FUN = function(f) { all(is.na(f)) }
		);

	# need to exclude columns with all NAs (occurs when perchip=T and there are no reference samples on that chip!)
	if (any(na.counts)) {
		all.na <- which(na.counts);
		print(paste("dropping:", all.na));
		which.n <- which.n[-all.na];
		which.cna <- which.cna[which.n];
		output.round <- output[, which.cna, drop = FALSE];
		}
	else{
		output.round <- output[, which.cna, drop = FALSE];
		}

	# loop over each sample
	for (col.ind in 1:ncol(output.round)) {
		if(2 == length(kd.thresh)){
			cna.thresh.single <- get.sample.specific.cna.thresholds(method = 4, data = output.round[ , col.ind], percent = kd.thresh[1])[1:2];
			cna.thresh.multi <- get.sample.specific.cna.thresholds(method = 4, data = output.round[ , col.ind], percent = kd.thresh[2])[1:2];
			output.round[ , col.ind] <- as.vector(call.cna.states(data.frame(log2ratio=output.round[ , col.ind]), c(cna.thresh.multi[1], cna.thresh.single, cna.thresh.multi[2]))$CN) + 2;
		}else if(4 == length(kd.thresh)){
			thresh <- vector(length = 4);
			thresh[1] <- get.sample.specific.cna.thresholds(method = 4, data = output.round[ , col.ind], percent = kd.thresh[1])[1];	# hom del
			thresh[2] <- get.sample.specific.cna.thresholds(method = 4, data = output.round[ , col.ind], percent = kd.thresh[2])[1];	# het del
			thresh[3] <- get.sample.specific.cna.thresholds(method = 4, data = output.round[ , col.ind], percent = kd.thresh[3])[2];	# het gain
			thresh[4] <- get.sample.specific.cna.thresholds(method = 4, data = output.round[ , col.ind], percent = kd.thresh[4])[2];	# hom gain
			output.round[ , col.ind] <- as.vector(call.cna.states(data.frame(log2ratio=output.round[ , col.ind]), thresh)$CN) + 2;
			}
		}

	return(output.round);
	}
