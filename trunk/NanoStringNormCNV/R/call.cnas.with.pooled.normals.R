call.cnas.with.pooled.normals <- function(normalized.data, phenodata, per.chip = FALSE, kd.option = 0, kd.values = NULL) {
	flog.warn("Currently, cannot call on chromosomes X and Y!");

	# use non-control probes (from autosomes only)
	use.genes <- which(normalized.data$CodeClass %in% c("Endogenous", "Housekeeping", "Invariant"));
	use.genes <- use.genes[!(use.genes %in% grep("chr[XY]", normalized.data$Name))];

	cna.normals <- matrix(
		nrow = length(use.genes),
		ncol = length(which(phenodata$type == 'Reference'))
		);
	cna.normals.unadj <- matrix(
		nrow = length(use.genes),
		ncol = length(which(phenodata$type == 'Reference'))
		);

	is.tmr  <- which(phenodata$type == 'Tumour');
	is.ref  <- which(phenodata$type == 'Reference');

	cna.raw <- NanoStringNormCNV::call.copy.number.state(
		input = normalized.data[use.genes,],
		reference = phenodata$SampleID[is.ref],
		per.chip = per.chip,
		chip.info = phenodata,
		thresh.method = 'none',
		adjust = TRUE
		);
	
	# make an average ref sample to use to call CNAs in normals
	norm.data.normals.only <- cbind(
		normalized.data[, c(1:3, (is.ref + 3))],
		avg.ref = apply(X = normalized.data[, (is.ref + 3)], MARGIN = 1, FUN = mean)
		);

	cna.normals.unadj <- NanoStringNormCNV::call.copy.number.state(
		input = norm.data.normals.only[use.genes,],
		reference = 'avg.ref',
		per.chip = per.chip,
		chip.info = phenodata[is.ref,],
		thresh.method = 'none',
		adjust = TRUE
		);
	cna.normals.unadj <- cna.normals.unadj[, -c(1:3)];
	
	if (kd.option <= 1) {
		if (kd.option == 0) {
			# NanoString recommended thresholds
			thresh <- c(0.4, 1.5, 2.5, 3.5);
		} else {
			# thresh.offset <- diff(range(cna.normals.unadj) * 0.15);

			thresh <- c(
				min(cna.normals.unadj),
				quantile(
					x = unlist(cna.normals.unadj),
					probs = c(0.1, 0.9),
					names = FALSE
					),
				# min(cna.normals.unadj) + thresh.offset,
				# max(cna.normals.unadj) - thresh.offset,
				max(cna.normals.unadj)
				);
			}

		cna.rounded <- NanoStringNormCNV::call.copy.number.state(
			input = normalized.data[use.genes,],
			reference = phenodata$SampleID[is.ref],
			per.chip = per.chip,
			chip.info = phenodata,
			adjust = TRUE,
			cna.thresh = thresh
			);
	} else {
		# call copy number states using kernel density values
		if (kd.option == 3) { 
			if (length(kd.values) != 4 & !is.numeric(kd.values)) {
				flog.warn(paste0(
					"For 'kd.option' 3, user must provide 4 kernel density values!\n",
					"Switching to default values ('kd.option' 2)."
					));
				kd.option <- 1;
				}
			}
		if (kd.option == 2) { kd.values <- c(0.9, 0.87, 0.93, 0.96); }# put whatever ends up being the default in apply.kd.cna.thresh here!!

		cna.rounded <- NanoStringNormCNV::call.copy.number.state(
			input = normalized.data[use.genes,],
			reference = phenodata$SampleID[is.ref],
			per.chip = per.chip,
			chip.info = phenodata,
			thresh.method = 'KD',
			kd.vals = kd.values,
			adjust = TRUE
			);
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
