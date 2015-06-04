
apply.kd.cna.thresh <- function(tmr2ref, kd.thresh){
	# assign header names
	headers <- c('Code.Class', 'CodeClass', 'Name', 'Accession');

	# define sample columns
	which.cna <- colnames(tmr2ref)[!colnames(tmr2ref) %in% headers];
	which.n   <- which(colnames(tmr2ref) %in% which.cna);

	# pull below into a separate object
	na.counts <- apply(
		X = tmr2ref[,which.n, drop = FALSE],
		MARGIN = 2,
		FUN = function(f) { all(is.na(f)) }
		);

	# need to exclude columns with all NAs (occurs when perchip=T and there are no reference samples on that chip!)
	if (any(na.counts)) {
		all.na <- which(na.counts);
		print(paste("dropping:", all.na));
		which.n <- which.n[-all.na];
		which.cna <- which.cna[which.n];
		cna.output <- tmr2ref[, which.cna, drop = FALSE];
		}
	else{
		cna.output <- tmr2ref[, which.cna, drop = FALSE];
		}

	# Determine the thresholds based on all patients combined (shown to be more stable if only considering small subset of patients)
	if(2 == length(kd.thresh)){
		cna.thresh.single <- get.sample.specific.cna.thresholds(method = 4, data = unlist(cna.output), percent = kd.thresh[1])[1:2];
		cna.thresh.multi <- get.sample.specific.cna.thresholds(method = 4, data = unlist(cna.output), percent = kd.thresh[2])[1:2];
		# loop over each sample
		for (col.ind in 1:ncol(cna.output)) {
			cna.output[ , col.ind] <- as.vector(call.cna.states(data.frame(log2ratio=cna.output[ , col.ind]), c(cna.thresh.multi[1], cna.thresh.single, cna.thresh.multi[2]))$CN) + 2;
			}
	}else if(4 == length(kd.thresh)){
		thresh <- vector(length = 4);
		thresh[1] <- get.sample.specific.cna.thresholds(method = 4, data = unlist(cna.output), percent = kd.thresh[1])[1];	# hom del
		thresh[2] <- get.sample.specific.cna.thresholds(method = 4, data = unlist(cna.output), percent = kd.thresh[2])[1];	# het del
		thresh[3] <- get.sample.specific.cna.thresholds(method = 4, data = unlist(cna.output), percent = kd.thresh[3])[2];	# het gain
		thresh[4] <- get.sample.specific.cna.thresholds(method = 4, data = unlist(cna.output), percent = kd.thresh[4])[2];	# hom gain
		# loop over each sample
		for (col.ind in 1:ncol(cna.output)) {
			cna.output[ , col.ind] <- as.vector(BoutrosLab.utilities.copynumber::call.cna.states(data.frame(log2ratio=cna.output[ , col.ind]), unlist(thresh))$CN) + 2;
			}
		}
	

	return(cna.output);
	}
