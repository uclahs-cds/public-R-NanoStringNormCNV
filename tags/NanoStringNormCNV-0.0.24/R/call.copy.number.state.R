
call.copy.number.state <- function (input, reference, sex.info, per.chip = FALSE, chip.info = NULL, thresh.method = 'round', multi.factor = 2, kd.vals = c(0.85, 0.95), adjust = FALSE, cna.thresh = c(0.4, 1.5, 2.5, 3.5)) {

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

	# grep X and Y-chromosome genes
	x.genes <- grep(x = tolower(input$Name), pattern = 'chrx');
	y.genes <- grep(x = tolower(input$Name), pattern = 'chry');

	if (length(x.genes) > 0 | length(y.genes) > 0) {
		flog.info("Identified the following as sex chromosome probes:");
		cat(paste(c("\t", input$Name[x.genes], input$Name[y.genes], "\n"), collapse = "\n\t"));
		}

	# removing chrY probes from female samples
	for (i in sex.info[sex.info$sex %in% "F",]$SampleID) {
		input[y.genes, i] <- NA;
		}
	
	# removing chrX and chrY probes where sex is not provided
	for (i in sex.info[is.na(sex.info$sex),]$SampleID) {
		input[c(x.genes, y.genes), i] <- NA;	
		}

	### Analysis
	# get tumour ratios
	out.cna <- NanoStringNormCNV::get.tumour.normal.ratio(
		ns.counts = input,
		ref = reference,
		chips.info = chip.info,
		per.chip = per.chip
		);

	# boosting probes
	out.cna <- out.cna * multi.factor;

	# if specified to make the median CN = multi.factor, adjust the values
	if (adjust) {
		out.cna <- apply(out.cna, 2, function(f) f - (median(f, na.rm = TRUE) - multi.factor));
		}

	# undo boosting of chrX/Y probe values in male samples
	for (i in sex.info[sex.info$sex %in% 'M',]$SampleID) {
		out.cna[x.genes, i] <- (out.cna[x.genes, i] / multi.factor);
		out.cna[y.genes, i] <- (out.cna[y.genes, i] / multi.factor);
		}

	# round if specified (based on NS recommendataions)
	if (thresh.method == 'round') {

		# segment using NS pre-defined thresholds
		out.cna.final <- NanoStringNormCNV::apply.ns.cna.thresh(
			tmr2ref = out.cna,
			thresh = cna.thresh
			);

	} else if (thresh.method == 'KD') {

		# segment using kernel density
		out.cna.final <- NanoStringNormCNV::apply.kd.cna.thresh(
			tmr2ref = out.cna,
			kd.thresh = kd.vals,
			sex.info = sex.info,
			sex.probes = rownames(out.cna)[c(x.genes, y.genes)]
			);

	} else {

		# else return as is
		out.cna.final <- out.cna;

		}

	# add the probe information back to out.cna.round
	header.names <- c('Code.Class', 'CodeClass', 'Name', 'Accession');
	out.cna.final <- cbind(
		input[,colnames(input)[colnames(input) %in% header.names], drop = FALSE],
		out.cna.final
		);

	return(out.cna.final);
	}
