apply.kd.cna.thresh <- function(ratio.data, kd.values) {
	# assign header names
	headers <- c('Code.Class', 'CodeClass', 'Name', 'Accession');

	# define sample columns
	which.cna <- colnames(ratio.data)[!colnames(ratio.data) %in% headers];
	which.n   <- which(colnames(ratio.data) %in% which.cna);

	# need to exclude columns with all NAs:
	# occurs when perchip = T and there are no reference samples on that chip!
	na.counts <- apply(
		X = ratio.data[,which.n, drop = FALSE],
		MARGIN = 2,
		FUN = function(f) { all(is.na(f)) }
		);

	if (any(na.counts)) {
		all.na <- which(na.counts);
		which.n <- which.n[-all.na];
		which.cna <- which.cna[which.n];
		}
	cna.output <- ratio.data[, which.cna, drop = FALSE];
		
	# determine the thresholds based on all patients combined
	# shown to be more stable if only considering small subset of patients
	if (2 == length(kd.values)) {
		cna.thresh.single <- NanoStringNormCNV::get.cna.thresholds(ratios = unlist(cna.output), percent = kd.values[1]); # het
		cna.thresh.multi  <- NanoStringNormCNV::get.cna.thresholds(ratios = unlist(cna.output), percent = kd.values[2]); # hom

		# loop over each sample to call CNAs
		for (col.ind in 1:ncol(cna.output)) {
			cna.output[, col.ind] <- NanoStringNormCNV::tumour.normal.ratio.to.cn.state(
				ratios = cna.output[, col.ind],
				thresholds = c(cna.thresh.multi[1], cna.thresh.single, cna.thresh.multi[2])
				);
			}
	} else if (4 == length(kd.values)) {
		thresh <- vector(length = 4);

		thresh[1] <- NanoStringNormCNV::get.cna.thresholds(ratios = unlist(cna.output), percent = kd.values[1])[1]; # hom del
		thresh[2] <- NanoStringNormCNV::get.cna.thresholds(ratios = unlist(cna.output), percent = kd.values[2])[1]; # het del
		thresh[3] <- NanoStringNormCNV::get.cna.thresholds(ratios = unlist(cna.output), percent = kd.values[3])[2]; # het gain
		thresh[4] <- NanoStringNormCNV::get.cna.thresholds(ratios = unlist(cna.output), percent = kd.values[4])[2]; # hom gain

		# loop over each sample to call CNAs
		for (col.ind in 1:ncol(cna.output)) {
			cna.output[, col.ind] <- NanoStringNormCNV::tumour.normal.ratio.to.cn.state(
				ratios = cna.output[, col.ind],
				thresholds = unlist(thresh)
				);
			}
		}

	return(cna.output);
	}
