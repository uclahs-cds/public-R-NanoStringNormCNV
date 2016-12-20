process.xy.probes <- function(ns.data, sex.info) {
	# check sex info is provided for each sample
	if (!all(colnames(ns.data[, -(1:3)]) %in% sex.info$SampleID)) {
		stop("Must provide sex information (M, F, or NA) for every sample in dataset!");
		}

	# identify xy probes by pattern
	x.genes <- grep(x = tolower(ns.data$Name), pattern = 'chrx');
	y.genes <- grep(x = tolower(ns.data$Name), pattern = 'chry');

	# notify user
	if (length(x.genes) > 0 | length(y.genes) > 0) {
		sex.probes <- c(as.vector(ns.data$Name[x.genes]), as.vector(ns.data$Name[y.genes]));
		flog.info("Identified the following as sex chromosome probes:");
		cat(paste(c("\t", sex.probes, "\n"), collapse = "\n\t"));
	} else {
		sex.probes <- NULL;
		flog.warn("Identified no sex chromosome probes!");
		}

	# extract male sex probe info and remove from 'ns.data' for separate processing
	ns.data.XY <- cbind(
		ns.data[c(x.genes, y.genes), 1:3, drop = FALSE],
		ns.data[
			c(x.genes, y.genes),
			colnames(ns.data) %in% sex.info$SampleID[sex.info$Sex %in% 'M', drop = FALSE]
			]
		);

	for (i in c(x.genes, y.genes)) {
		for (j in which(colnames(ns.data) %in% sex.info$SampleID[sex.info$Sex %in% 'M'])) {
			ns.data[i, j] <- NA;
			}
		}

	# remove chrY probes from female samples
	for (i in sex.info[sex.info$Sex %in% "F",]$SampleID) {
		ns.data[y.genes, i] <- NA;
		}
	
	# remove chrX and chrY probes where sex is not provided
	if (any(is.na(sex.info$Sex)) & ! is.null(sex.probes)) {
		flog.info("Removing XY probes where sample's sex is not available:");
		for (i in sex.info[is.na(sex.info$Sex),]$SampleID) {
			cat(paste(c("\t", i, "\n")));
			ns.data[c(x.genes, y.genes), i] <- NA;
			}
		cat("\n");
		}

	# do not return an empty data-frame
	if (ncol(ns.data.XY[, -(1:3)]) == 0) {
		ns.data.XY <- NULL;
		}

	output <- list(
		ns.data.without.maleXY = ns.data,
		ns.data.maleXY.only = ns.data.XY,
		sex.probes = sex.probes
		);
	
	return(output);
	}