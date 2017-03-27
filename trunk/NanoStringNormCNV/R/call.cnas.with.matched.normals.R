call.cnas.with.matched.normals <- function(
	normalized.data,
	phenodata,
	per.chip = FALSE,
	call.method = 1,
	kd.values = c(0.85, 0.95),
	use.sex.info = TRUE
	) {

	# use non-control probes only
	use.codeclass <- c("Endogenous", "Housekeeping", "Invariant");

	# ensure sample order matches
	phenodata <- phenodata[match(colnames(normalized.data)[-(1:3)], phenodata$SampleID),];
	
	sex.probes <- NULL;
	if (use.sex.info) {
		if (! 'Sex' %in% colnames(phenodata)) {
			stop("Sex information must be provided in phenodata if use.sex.info = TRUE!");
			}

		# identify and process XY probes separately
		xy.processed.data <- process.xy.probes(
			normalized.data = normalized.data,
			sex.info = phenodata[, c("SampleID", "Sex")]
			);

		sex.probes <- xy.processed.data$sex.probes;
		if (! is.null(sex.probes)) {
			normalized.data    <- xy.processed.data$normalized.data.without.maleXY;
			normalized.data.XY <- xy.processed.data$normalized.data.maleXY.only;

			# if there are few chrXY probes, inform user
			if (!is.null(normalized.data.XY) && nrow(normalized.data.XY) < 40) {
				flog.warn("Low chrX/chrY probe number! Consider using pooled normals when calling CNAs in male sex chromosomes");
				}
			
			use.genes.XY <- which(normalized.data.XY$CodeClass %in% use.codeclass);
			has.ref.XY   <- phenodata$SampleID[!(phenodata$ReferenceID %in% 'missing') & phenodata$Type == 'Tumour' & phenodata$Sex %in% 'M'];
			}
		}

	use.genes <- which(normalized.data$CodeClass %in% use.codeclass);
	has.ref   <- phenodata$SampleID[!(phenodata$ReferenceID %in% 'missing') & phenodata$Type == 'Tumour'];

	# set up output variables
	cna.raw	<- as.data.frame(matrix(
		nrow = length(use.genes),
		ncol = length(has.ref),
		dimnames = list(normalized.data$Name[use.genes], has.ref)
		));
	cna.rounded <- as.data.frame(matrix(
		nrow = length(use.genes),
		ncol = length(has.ref),
		dimnames = list(normalized.data$Name[use.genes], has.ref)
		));

	# iterate through each sample here
	for (tmr in has.ref) {
		ref <- phenodata[phenodata$SampleID == tmr,]$ReferenceID;

		tmr.ind <- which(colnames(normalized.data) == tmr);
		ref.ind <- which(colnames(normalized.data) == ref);
		input.data <- normalized.data[use.genes, c(1:3, tmr.ind, ref.ind)];
		
		if (!is.null(sex.probes)) {
			if (!is.null(normalized.data.XY)) {
				tmr.ind.XY <- which(colnames(normalized.data.XY) == tmr);
				ref.ind.XY <- which(colnames(normalized.data.XY) == ref);
				input.data.XY <- normalized.data.XY[use.genes.XY, c(1:3, tmr.ind.XY, ref.ind.XY)];
				}
			}

		chip.info <- phenodata[
			phenodata$SampleID %in% c(tmr, ref),
			c("SampleID", "Cartridge")
			];

		raw.ratios <- NanoStringNormCNV::call.copy.number.values(
			normalized.data = input.data,
			reference = ref,
			thresh.method = 'none',
			chip.info = chip.info,
			multi.factor = 2
			);
		if (tmr %in% names(raw.ratios)) {
			cna.raw[, tmr] <- raw.ratios[, tmr];
			}

		# calculate tumour/normal ratios (for male sex chrom probes, if any)
		if (!is.null(sex.probes) && tmr %in% has.ref.XY && length(use.genes.XY) > 0) {
			raw.ratios.xy <- NanoStringNormCNV::call.copy.number.values(
				normalized.data = input.data.XY,
				reference = ref,
				thresh.method = 'none',
				chip.info = chip.info,
				multi.factor = 1
				);
			if (tmr %in% names(raw.ratios.xy)) {
				cna.raw[sex.probes, tmr] <- raw.ratios.xy[, tmr];
				}
			}		

		if (call.method == 1) {
			### Naive thresholds
			thresh <- c(0.4, 1.5, 2.5, 3.5);

			# call CNAs in tumours (for autosome and female sex chrom probes)
			round.ratios <- NanoStringNormCNV::call.copy.number.values(
				normalized.data = input.data,
				reference = ref,
				per.chip = per.chip,
				chip.info = chip.info,
				cna.thresh = thresh,
				multi.factor = 2
				);
			if (tmr %in% names(round.ratios)) {
				cna.rounded[, tmr] <- round.ratios[, tmr];
				}

			# call CNAs in tumours (for male sex chrom probes)
			if (!is.null(sex.probes) && tmr %in% has.ref.XY && length(use.genes.XY) > 0) {
				round.ratios.xy <- NanoStringNormCNV::call.copy.number.values(
					normalized.data = input.data.XY,
					reference = ref,
					per.chip = per.chip,
					chip.info = chip.info,
					multi.factor = 1,
					cna.thresh = thresh
					);
				if (tmr %in% names(round.ratios.xy)) {
					cna.rounded[sex.probes, tmr] <- round.ratios.xy[, tmr];
					}
				}			
		} else if (call.method == 2) {
			### Call copy number states using kernel density values
			if ((length(kd.values) != 4 & length(kd.values) != 2) | !is.numeric(kd.values)) {
				flog.warn(paste0(
					"For 'call.method' 2, user must provide 4 or 2 kernel density values!\n",
					"Switching to default KD values: 0.85, 0.95."
					));
				kd.values <- c(0.85, 0.95);
				}

			# call CNAs in tumours (for autosome and female sex chrom probes)
			round.ratios <- NanoStringNormCNV::call.copy.number.values(
				normalized.data = input.data,
				reference = ref,
				per.chip = per.chip,
				chip.info = chip.info,
				thresh.method = 'KD',
				kd.values = kd.values,
				multi.factor = 2
				);
			if (tmr %in% names(round.ratios)) {
				cna.rounded[, tmr] <- round.ratios[, tmr];
				}

			# call CNAs in tumours (for male sex chrom probes)
			if (!is.null(sex.probes) && tmr %in% has.ref.XY && length(use.genes.XY) > 0) {
				round.ratios.xy <- NanoStringNormCNV::call.copy.number.values(
					normalized.data = input.data.XY,
					reference = ref,
					per.chip = per.chip,
					chip.info = chip.info,
					multi.factor = 1,
					thresh.method = 'KD',
					kd.values = kd.values
					);
				if (tmr %in% names(round.ratios.xy)) {
					cna.rounded[sex.probes, tmr] <- round.ratios.xy[, tmr];
					}
				}
		} else {
			stop("Argument 'call.method' accepts only values 1 or 2! Please consult documentation!");
			}
		}

	# remove all NA columns (if any)
	all.na <- which(apply(cna.rounded, 2, function(x) all(is.na(x))));
	if (length(all.na) > 0) { cna.rounded <- cna.rounded[, -all.na]; }

	# combine output
	cna.all <- list(rounded = cna.rounded, raw = cna.raw);

	return(cna.all);
	
	}
