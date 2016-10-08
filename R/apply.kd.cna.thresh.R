
apply.kd.cna.thresh <- function(tmr2ref, kd.thresh, sex.info, sex.probes = NULL) {
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

	# need to exclude columns with all NAs:
	# occurs when perchip = T and there are no reference samples on that chip!
	if (any(na.counts)) {
		all.na <- which(na.counts);
		print(paste("dropping:", all.na));
		which.n <- which.n[-all.na];
		which.cna <- which.cna[which.n];
		}
	cna.output <- tmr2ref[, which.cna, drop = FALSE];

####WIP
	# if male, exclude sex chrom probes from threshold calculations
	cna.output.tmp <- cna.output;
	for (i in sex.info[sex.info$sex == "M",]$SampleID) {
		if (is.null(sex.probes)) {
			stop("Must provide names of sex chromosome probes if sample is male!");
		} else {
			cna.output.tmp[rownames(cna.output.tmp) %in% sex.probes, i] <- NA;
			}
		}
####
		
	# determine the thresholds based on all patients combined
	# shown to be more stable if only considering small subset of patients
	if (2 == length(kd.thresh)) {
		cna.thresh.single <- NanoStringNormCNV::get.sample.specific.cna.thresholds(cna.data = as.vector(cna.output.tmp), percent = kd.thresh[1]);
		cna.thresh.multi  <- NanoStringNormCNV::get.sample.specific.cna.thresholds(cna.data = as.vector(cna.output.tmp), percent = kd.thresh[2]);

		# loop over each sample
		for (col.ind in 1:ncol(cna.output)) {
			cna.output[, col.ind] <- NanoStringNormCNV::tumour.normal.ratio.to.cn.state(
				ratio.data = data.frame(ratio = cna.output[ , col.ind]),
				thresholds = c(cna.thresh.multi[1], cna.thresh.single, cna.thresh.multi[2])
				);
			}
	} else if (4 == length(kd.thresh)) {
		thresh <- vector(length = 4);

		thresh[1] <- NanoStringNormCNV::get.sample.specific.cna.thresholds(cna.data = as.vector(cna.output.tmp), percent = kd.thresh[1])[1]; # hom del
		thresh[2] <- NanoStringNormCNV::get.sample.specific.cna.thresholds(cna.data = as.vector(cna.output.tmp), percent = kd.thresh[2])[1]; # het del
		thresh[3] <- NanoStringNormCNV::get.sample.specific.cna.thresholds(cna.data = as.vector(cna.output.tmp), percent = kd.thresh[3])[2]; # het gain
		thresh[4] <- NanoStringNormCNV::get.sample.specific.cna.thresholds(cna.data = as.vector(cna.output.tmp), percent = kd.thresh[4])[2]; # hom gain

		# loop over each sample
		for (col.ind in 1:ncol(cna.output)) {
			cna.output[, col.ind] <- NanoStringNormCNV::tumour.normal.ratio.to.cn.state(
				ratio.data = data.frame(ratio = cna.output[ , col.ind]),
				thresholds = unlist(thresh)
				);
			}
		}

####WIP
	# if male, boost sex chrom probes
	cna.output.tmp <- cna.output;
	for (i in sex.info[sex.info$sex == "M",]$SampleID) {
		cna.output[rownames(cna.output) %in% sex.probes, i] <- cna.output[rownames(cna.output) %in% sex.probes, i] + 1;
		}
####

	return(cna.output);
	}
