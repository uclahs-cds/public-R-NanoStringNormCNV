call.cnas.with.pooled.normals <- function(
	normalized.data,
	phenodata,
	per.chip = FALSE,
	call.method = 0,
	kd.values = NULL,
	use.sex.info = TRUE
	) {
	
	# use non-control probes only
	use.codeclass <- c("Endogenous", "Housekeeping", "Invariant");

	# ensure sample order matches
	phenodata <- phenodata[match(colnames(normalized.data)[-(1:3)], phenodata$SampleID),];

	is.tmr <- which(phenodata$Type == 'Tumour');
	is.ref <- which(phenodata$Type == 'Reference');

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
			
			is.ref.XY <- which(phenodata[phenodata$Sex %in% 'M',]$Type == 'Reference');

			use.genes.XY <- which(normalized.data.XY$CodeClass %in% use.codeclass);
			}
		}

	use.genes <- which(normalized.data$CodeClass %in% use.codeclass);

	# calculate tumour/normal ratios (for autosome and female sex chrom probes)
	cna.raw <- NanoStringNormCNV::call.copy.number.state(
		normalized.data = normalized.data[use.genes,],
		reference = phenodata$SampleID[is.ref],
		per.chip = per.chip,
		chip.info = phenodata[, c('SampleID', 'Cartridge')],
		thresh.method = 'none',
		multi.factor = 2,
		adjust = TRUE
		);

	# calculate tumour/normal ratios (for male sex chrom probes)
	if (!is.null(sex.probes) && ncol(normalized.data.XY) > 0 && length(use.genes.XY) > 0) {
		cna.raw.XY <- NanoStringNormCNV::call.copy.number.state(
			normalized.data = normalized.data.XY[use.genes.XY,],
			reference = phenodata[phenodata$Sex %in% 'M',]$SampleID[is.ref.XY],
			per.chip = per.chip,
			chip.info = phenodata[, c('SampleID', 'Cartridge')],
			thresh.method = 'none',
			multi.factor = 1,
			adjust = TRUE
			)[, -(1:3)];

		# add to main data variable
		for (i in colnames(cna.raw.XY)) {
			cna.raw[sex.probes, i] <- cna.raw.XY[, i];
			}
		}

	# call CNAs in normals (for autosome and female sex chrom probes)
	cna.normals.unadj <- NanoStringNormCNV::call.copy.number.state(
		normalized.data = normalized.data[, c(1:3, (is.ref + 3))],
		reference = phenodata$SampleID[is.ref],
		per.chip = FALSE,
		chip.info = phenodata[, c('SampleID', 'Cartridge')],
		thresh.method = 'none',
		multi.factor = 2,
		adjust = TRUE
		)[, -(1:3)];

	# call CNAs in tumours
	if (call.method <= 1) {
		if (call.method == 0) {
			### Naive thresholds
			thresh <- c(0.4, 1.5, 2.5, 3.5);
		} else {
			### Thresholds from normal sample max/min values (excluding male sex chrom probes)
			# find median of sample minima and maxima
			minimum.median <- median(apply(cna.normals.unadj, 2, min, na.rm = TRUE));
			maximum.median <- median(apply(cna.normals.unadj, 2, max, na.rm = TRUE));

			# calculate threshold values
			thresh.offset <- (maximum.median - minimum.median) * 0.15;
			thresh <- c(
				minimum.median,
				minimum.median + thresh.offset,
				maximum.median - thresh.offset,
				maximum.median
				);
			}

		# call CNAs in tumours (for autosome and female sex chrom probes)
		cna.rounded <- NanoStringNormCNV::call.copy.number.state(
			normalized.data = normalized.data[use.genes,],
			reference = phenodata$SampleID[is.ref],
			per.chip = per.chip,
			chip.info = phenodata[, c('SampleID', 'Cartridge')],
			multi.factor = 2,
			adjust = TRUE,
			cna.thresh = thresh
			);

		# call CNAs in tumours (for male sex chrom probes)
		if (!is.null(sex.probes) && ncol(normalized.data.XY) > 0 && length(use.genes.XY) > 0) {
			cna.rounded.XY <- NanoStringNormCNV::call.copy.number.state(
				normalized.data = normalized.data.XY[use.genes.XY,],
				reference = phenodata[phenodata$Sex %in% 'M',]$SampleID[is.ref.XY],
				per.chip = per.chip,
				chip.info = phenodata[, c('SampleID', 'Cartridge')],
				multi.factor = 1,
				adjust = TRUE,
				cna.thresh = thresh
				)[, -(1:3)];

			# add to main data variable
			for (i in colnames(cna.rounded.XY)) {
				cna.rounded[sex.probes, i] <- cna.rounded.XY[, i];
				}
			}
	} else {
		### Call copy number states using kernel density values
		if (call.method == 3) {
			if ((length(kd.values) != 4 & length(kd.values) != 2) | !is.numeric(kd.values)) {
				flog.warn(paste0(
					"For 'call.method' 3, user must provide 4 or 2 kernel density values!\n",
					"Switching to default KD values (setting 'call.method' to 2)."
					));
				call.method <- 2;
				}
			}
		if (call.method == 2) { kd.values <- c(0.85, 0.95); }

		# call CNAs in tumours (for autosome and female sex chrom probes)
		cna.rounded <- NanoStringNormCNV::call.copy.number.state(
			normalized.data = normalized.data[use.genes,],
			reference = phenodata$SampleID[is.ref],
			per.chip = per.chip,
			chip.info = phenodata[, c('SampleID', 'Cartridge')],
			multi.factor = 2,
			thresh.method = 'KD',
			kd.values = kd.values,
			adjust = TRUE
			);

		# call CNAs in tumours (for male sex chrom probes)
		if (!is.null(sex.probes) && ncol(normalized.data.XY) > 0 && length(use.genes.XY) > 0) {
			cna.rounded.XY <- NanoStringNormCNV::call.copy.number.state(
				normalized.data = normalized.data.XY[use.genes.XY,],
				reference = phenodata[phenodata$Sex %in% 'M',]$SampleID[is.ref.XY],
				per.chip = per.chip,
				chip.info = phenodata[, c('SampleID', 'Cartridge')],
				multi.factor = 1,
				thresh.method = 'KD',
				kd.values = kd.values,
				adjust = TRUE
				)[, -(1:3)];

			# add to main data variable
			for (i in colnames(cna.rounded.XY)) {
				cna.rounded[sex.probes, i] <- cna.rounded.XY[, i];
				}
			}
		}

	# call CNAs in normals (for male sex chrom probes)
	if (!is.null(sex.probes) && ncol(normalized.data.XY) > 0 && length(use.genes.XY) > 0) {
		cna.normals.unadj.XY <- NanoStringNormCNV::call.copy.number.state(
			normalized.data = normalized.data.XY[, c(1:3, (is.ref.XY + 3))],
			reference = phenodata[phenodata$Sex %in% 'M',]$SampleID[is.ref.XY],
			per.chip = FALSE,
			chip.info = phenodata[, c('SampleID', 'Cartridge')],
			thresh.method = 'none',
			multi.factor = 1,
			adjust = TRUE
			)[, -(1:3)];

		# add to main data variable
		for (i in colnames(cna.normals.unadj.XY)) {
			cna.normals.unadj[sex.probes, i] <- cna.normals.unadj.XY[, i];
			}
		}

	# collect normal sample CNAs
	cna.normals <- cna.rounded[, phenodata$SampleID[is.ref]];

	cna.raw 	<- as.matrix(cna.raw[, colnames(cna.raw) %in% phenodata$SampleID[is.tmr]]);
	cna.rounded <- as.matrix(cna.rounded[, colnames(cna.rounded) %in% phenodata$SampleID[is.tmr]]);
	has.ref 	<- is.tmr[colnames(normalized.data)[is.tmr + 3] %in% colnames(cna.rounded)];

	colnames(cna.rounded) <- phenodata$SampleID[has.ref];
	colnames(cna.raw)     <- phenodata$SampleID[has.ref];

	# combine output
	cna.all <- list(
		rounded = cna.rounded,
		raw = cna.raw,
		normals = cna.normals,
		normals.unadj = cna.normals.unadj
		);

	return(cna.all);

	}
