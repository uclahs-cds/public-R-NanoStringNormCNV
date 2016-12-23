process.xy.probes <- function(normalized.data, sex.info) {
	# check sex info is provided for each sample
	if (!all(colnames(normalized.data[, -(1:3)]) %in% sex.info$SampleID)) {
		stop("Must provide sex information (M, F, or NA) for every sample in dataset!");
		}

	# identify xy probes by pattern
	x.genes <- grep(x = tolower(normalized.data$Name), pattern = 'chrx');
	y.genes <- grep(x = tolower(normalized.data$Name), pattern = 'chry');

	# notify user
	if (length(x.genes) > 0 | length(y.genes) > 0) {
		sex.probes <- c(as.vector(normalized.data$Name[x.genes]), as.vector(normalized.data$Name[y.genes]));
		flog.info("Identified the following as sex chromosome probes:");
		cat(paste(c("\t", sex.probes, "\n"), collapse = "\n\t"));
	} else {
		sex.probes <- NULL;
		flog.warn("Identified no sex chromosome probes!");
		}

	# extract male sex probe info and remove from 'normalized.data' for separate processing
	normalized.data.XY <- cbind(
		normalized.data[c(x.genes, y.genes), 1:3, drop = FALSE],
		normalized.data[
			c(x.genes, y.genes),
			colnames(normalized.data) %in% sex.info$SampleID[sex.info$Sex %in% 'M', drop = FALSE]
			]
		);

	for (i in c(x.genes, y.genes)) {
		for (j in which(colnames(normalized.data) %in% sex.info$SampleID[sex.info$Sex %in% 'M'])) {
			normalized.data[i, j] <- NA;
			}
		}

	# remove chrY probes from female samples
	for (i in sex.info[sex.info$Sex %in% "F",]$SampleID) {
		normalized.data[y.genes, i] <- NA;
		}
	
	# remove chrX and chrY probes where sex is not provided
	if (any(is.na(sex.info$Sex)) & ! is.null(sex.probes)) {
		flog.info("Removing XY probes where sample's sex is not available:");
		for (i in sex.info[is.na(sex.info$Sex),]$SampleID) {
			cat(paste(c("\t", i, "\n")));
			normalized.data[c(x.genes, y.genes), i] <- NA;
			}
		cat("\n");
		}

	# do not return an empty data-frame
	if (ncol(normalized.data.XY[, -(1:3)]) == 0) {
		normalized.data.XY <- NULL;
		}

	output <- list(
		normalized.data.without.maleXY = normalized.data,
		normalized.data.maleXY.only = normalized.data.XY,
		sex.probes = sex.probes
		);
	
	return(output);
	}