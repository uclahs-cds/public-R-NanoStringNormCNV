call.copy.number.state <- function (
	normalized.data,
	reference,
	per.chip = FALSE,
	chip.info = NULL,
	thresh.method = 'round',
	multi.factor = 2,
	kd.vals = c(0.85, 0.95),
	adjust = FALSE,
	cna.thresh = c(0.4, 1.5, 2.5, 3.5)
	) {

	# Check input
	if (! thresh.method %in% (unlist(strsplit("round KD kd none","\\s")))) {
		stop("Sorry method isn't currently supported. Please try one of round, KD, or none.");
		}

	if (toupper(thresh.method) == 'KD' & (2 != length(kd.vals) & 4 != length(kd.vals))) {
		stop("Please specify two or four values for KD thresholds.  The first should be for heterozygous and the second for homozygous if length 2. If length 4, the order should be hom deletion, het deletion, het gain, hom gain.");
		}

	# make sure kd values make sense
	if (toupper(thresh.method) == 'KD' & 2 == length(kd.vals) & kd.vals[1] > kd.vals[2]) {
		stop("Invalid KD thresholds -- the first should be for heterozygous and the second for homozygous.");
		}
	if (toupper(thresh.method) == 'KD' & 4 == length(kd.vals) & (kd.vals[1] < kd.vals[2] | kd.vals[3] > kd.vals[4])) {
		print(kd.vals);
		stop("Invalid KD thresholds -- the order should be hom deletion, het deletion, het gain, hom gain.");
		}

	# get tumour ratios
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

		# segment using kernel density
		out.cna.final <- NanoStringNormCNV::apply.kd.cna.thresh(
			ratio.data = out.cna,
			kd.values = kd.vals
			);

	} else {

		# else return as is
		out.cna.final <- out.cna;

		}

	# add the probe information back to out.cna.round
	header.names <- c('Code.Class', 'CodeClass', 'Name', 'Accession');
	rownames(normalized.data) <- normalized.data[, colnames(normalized.data) == 'Name'];
	out.cna.final <- cbind(
		normalized.data[,colnames(normalized.data)[colnames(normalized.data) %in% header.names], drop = FALSE],
		out.cna.final
		);

	return(out.cna.final);
	}
