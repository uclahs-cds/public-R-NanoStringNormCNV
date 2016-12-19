call.cnas.with.matched.normals <- function(
	normalized.data,
	phenodata,
	per.chip = FALSE,
	call.method = 1,
	kd.values = NULL,
	use.sex.info = TRUE
	) {

	# use non-control probes only
	use.codeclass <- c("Endogenous", "Housekeeping", "Invariant");

	# ensure sample order matches
	phenodata <- phenodata[match(colnames(normalized.data)[-(1:3)], phenodata$SampleID),];
	
	sex.probes <- NULL;
	if (use.sex.info) {
		# identify and process XY probes separately
		xy.processed.data <- process.xy.probes(
			ns.data = normalized.data,
			sex.info = phenodata[, c("SampleID", "sex")]
			);

		sex.probes <- xy.processed.data$sex.probes;
		if (! is.null(sex.probes)) {
			normalized.data    <- xy.processed.data$ns.data.without.maleXY;
			normalized.data.XY <- xy.processed.data$ns.data.maleXY.only;

			# if there are few chrXY probes, inform user
			if (!is.null(normalized.data.XY) && nrow(normalized.data.XY) < 40) {
				flog.warn("Low chrX/chrY probe number! Consider using pooled normals when calling CNAs in male sex chromosomes");
				}
			
			use.genes.XY <- which(normalized.data.XY$CodeClass %in% use.codeclass);
			has.ref.XY   <- phenodata$SampleID[!(phenodata$ref.name %in% 'missing') & phenodata$type == 'Tumour' & phenodata$sex %in% 'M'];
			}
		}

	use.genes <- which(normalized.data$CodeClass %in% use.codeclass);
	has.ref   <- phenodata$SampleID[!(phenodata$ref.name %in% 'missing') & phenodata$type == 'Tumour'];

	# set up output variables
	cna.raw 	<- matrix(nrow = length(use.genes), ncol = length(has.ref));
	cna.rounded <- matrix(nrow = length(use.genes), ncol = length(has.ref));

	colnames(cna.rounded) <- colnames(cna.raw) <- has.ref;
	rownames(cna.rounded) <- rownames(cna.raw) <- normalized.data$Name[use.genes];

	cna.raw 	<- as.data.frame(cna.raw);
	cna.rounded <- as.data.frame(cna.rounded);

	# iterate through each sample here
	for (tmr in has.ref) {
		ref <- phenodata[phenodata$SampleID == tmr,]$ref.name;

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

		cna.raw[, tmr] <- NanoStringNormCNV::call.copy.number.state(
			input = input.data,
			reference = ref,
			thresh.method = 'none',
			multi.factor = 2
			)[, 4];

		# calculate tumour-normal ratios (for male sex chrom probes, if any)
		if (!is.null(sex.probes) && tmr %in% has.ref.XY && length(use.genes.XY) > 0) {
			cna.raw[sex.probes, tmr] <- NanoStringNormCNV::call.copy.number.state(
				input = input.data.XY,
				reference = ref,
				thresh.method = 'none',
				multi.factor = 1
				)[, 4];
			}		

		chip.info <- phenodata[c(tmr.ind, ref.ind), c("SampleID", "cartridge")];
		
		if (call.method <= 1) {
			### NanoString recommended thresholds
			thresh <- c(0.4, 1.5, 2.5, 3.5);

			# call CNAs in tumours (for autosome and female sex chrom probes)
			cna.rounded[, tmr] <- NanoStringNormCNV::call.copy.number.state(
				input = input.data,
				reference = ref,
				per.chip = per.chip,
				chip.info = chip.info,
				cna.thresh = thresh,
				multi.factor = 2
				)[, 4];

			# call CNAs in tumours (for male sex chrom probes)
			if (!is.null(sex.probes) && tmr %in% has.ref.XY && length(use.genes.XY) > 0) {
				cna.rounded[sex.probes, tmr] <- NanoStringNormCNV::call.copy.number.state(
					input = input.data.XY,
					reference = ref,
					per.chip = per.chip,
					chip.info = chip.info,
					multi.factor = 1,
					cna.thresh = thresh
					)[, 4];
				}			
		} else {
			### Call copy number states using kernel density values
			if (call.method == 3) { 
				if ((length(kd.values) != 4 & length(kd.values) != 2) | !is.numeric(kd.values)) {
					flog.warn(paste0(
						"For 'call.method' 3, user must provide 4 or 2 kernel density values!\n",
						"Switching to default values (setting 'call.method' to 2)."
						));
					call.method <- 2;
					}
				}
			if (call.method == 2) { kd.values <- c(0.85, 0.95); }

			# call CNAs in tumours (for autosome and female sex chrom probes)
			cna.rounded[,tmr] <- NanoStringNormCNV::call.copy.number.state(
				input = input.data,
				reference = ref,
				per.chip = per.chip,
				chip.info = chip.info,
				thresh.method = 'KD',
				kd.vals = kd.values,
				multi.factor = 2
				)[, 4];

			# call CNAs in tumours (for male sex chrom probes)
			if (!is.null(sex.probes) && tmr %in% has.ref.XY && length(use.genes.XY) > 0) {
				cna.rounded[sex.probes, tmr] <- NanoStringNormCNV::call.copy.number.state(
					input = input.data.XY,
					reference = ref,
					per.chip = per.chip,
					chip.info = chip.info,
					multi.factor = 1,
					thresh.method = 'KD',
					kd.vals = kd.values
					)[, 4];
				}
			}
		}

	# combine output
	cna.all <- list(rounded = cna.rounded, raw = cna.raw);

	return(cna.all);
	
	}
