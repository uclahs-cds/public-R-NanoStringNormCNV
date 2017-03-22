call.copy.number.values <- function(
	normalized.data,
	reference,
	per.chip = FALSE,
	chip.info = NULL,
	thresh.method = 'round',
	cna.thresh = c(0.4, 1.5, 2.5, 3.5),
	kd.values = c(0.85, 0.95),
	multi.factor = 2,
	adjust = FALSE
	) {

	# check input
	if (! thresh.method %in% (unlist(strsplit("round KD kd none","\\s")))) {
		stop("Sorry method isn't currently supported. Please try one of 'round', 'KD', or 'none'.");
		}

	if (toupper(thresh.method) == 'KD' & (2 != length(kd.values) & 4 != length(kd.values))) {
		stop(paste(
			"Please specify two or four values for KD thresholds.",
			"The first should be for heterozygous and the second for homozygous if length is 2.",
			"If length is 4, the order should be hom deletion, het deletion, het gain, hom gain."
			));
		}

	# make sure kd values make sense
	if (toupper(thresh.method) == 'KD' & 2 == length(kd.values) & kd.values[1] > kd.values[2]) {
		stop("Invalid KD thresholds -- the first should be for heterozygous and the second for homozygous.");
		}
	if (toupper(thresh.method) == 'KD' & 4 == length(kd.values) & (kd.values[1] < kd.values[2] | kd.values[3] > kd.values[4])) {
		stop("Invalid KD thresholds -- the order should be hom deletion, het deletion, het gain, hom gain.");
		}

	# get tumour/normal ratios
	out.cna <- NanoStringNormCNV::get.tumour.normal.ratio(
		normalized.data = normalized.data,
		reference = reference,
		chip.info = chip.info,
		per.chip = per.chip
		);

	# boost probes
	out.cna <- out.cna * multi.factor;

	# if specified to make the median CN = multi.factor, adjust the values
	if (adjust) {
		out.cna <- apply(out.cna, 2, function(f) { f - (median(f, na.rm = TRUE) - multi.factor) });
		out.cna[out.cna < 0] <- 0;
		}

	# round if specified (based on NS recommendataions)
	if (thresh.method == 'round') {

		# segment using set thresholds
		out.cna.final <- NanoStringNormCNV::apply.ns.cna.thresh(
			ratio.data = out.cna,
			cna.thresh = cna.thresh
			);

	} else if (thresh.method == 'KD') {

		# segment using thresholds obtained through kernel density approach
		out.cna.final <- NanoStringNormCNV::apply.kd.cna.thresh(
			ratio.data = out.cna,
			kd.values = kd.values,
			neutral.cn = multi.factor
			);

	} else {

		# else return as is
		out.cna.final <- out.cna;

		}

	# add the probe information back to output
	header.names <- c('Code.Class', 'CodeClass', 'Name', 'Accession');
	rownames(normalized.data) <- normalized.data[, colnames(normalized.data) == 'Name'];
	out.cna.final <- cbind(
		normalized.data[,colnames(normalized.data)[colnames(normalized.data) %in% header.names], drop = FALSE],
		out.cna.final
		);

	return(out.cna.final);
	}
