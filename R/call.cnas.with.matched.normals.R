call.cnas.with.matched.normals <- function(
	normalized.data,
	phenodata,
	per.chip = FALSE,
	call.method = 0,
	kd.values = NULL
	) {
	
	# use non-control probes
	use.genes <- which(normalized.data$CodeClass %in% c("Endogenous", "Housekeeping", "Invariant"));

	has.ref 	<- which(phenodata$ref.name != 'missing' & phenodata$type == 'Tumour');
	cna.raw 	<- matrix(nrow = length(use.genes), ncol = length(has.ref));
	cna.rounded <- matrix(nrow = length(use.genes), ncol = length(has.ref));

	# iterate through each sample here
	for (tmr in 1:length(has.ref)) {
		tmr.ind <- which(colnames(normalized.data) == phenodata$SampleID[has.ref[tmr]]);
		ref.ind <- which(colnames(normalized.data) == phenodata$ref.name[has.ref[tmr]]);

		sample.sex <- phenodata[
				c(has.ref[tmr], which(phenodata$SampleID == phenodata$ref.name[has.ref[tmr]])),
				c("SampleID", "sex")
				];

		cna.raw[,tmr] <- NanoStringNormCNV::call.copy.number.state(
			input = normalized.data[use.genes, c(1:3, tmr.ind, ref.ind), drop = FALSE],
			reference = phenodata$ref.name[has.ref[tmr]],
			sex.info = sample.sex,
			thresh.method = 'none'
			)[,4];

		if (call.method <= 1) {
			if (call.method == 0) {
				# NanoString recommended thresholds
				thresh <- c(0.4, 1.5, 2.5, 3.5);
			} else {
				# Thresholds from max/min values
				thresh.offset <- diff(range(normalized.data[use.genes, ref.ind], na.rm = TRUE) * 0.15);

				thresh <- c(
					min(normalized.data[use.genes, ref.ind], na.rm = TRUE),
					# # using quantiles seems more robust to outliers!
					# quantile(
					# 	x = unlist(normalized.data[use.genes, ref.ind]),
					# 	probs = c(0.1, 0.9),
					# 	names = FALSE,
					# 	na.rm = TRUE
					# 	),
					min(normalized.data[use.genes, ref.ind], na.rm = TRUE) + thresh.offset,
					max(normalized.data[use.genes, ref.ind], na.rm = TRUE) - thresh.offset,
					max(normalized.data[use.genes, ref.ind], na.rm = TRUE)
					);
				}

			cna.rounded[,tmr] <- NanoStringNormCNV::call.copy.number.state(
				input = normalized.data[use.genes, c(1:3, tmr.ind, ref.ind), drop = FALSE],
				reference = phenodata$ref.name[has.ref[tmr]],
				sex.info = sample.sex,	
				per.chip = per.chip,
				chip.info = phenodata,
				cna.thresh = thresh
				)[,4];
		} else {
			# call copy number states using kernel density values
			if (call.method == 3) { 
				if ((length(kd.values) != 4 & length(kd.values) != 2) | !is.numeric(kd.values)) {
					flog.warn(paste0(
						"For 'call.method' 3, user must provide 4 or 2 kernel density values!\n",
						"Switching to default values (setting 'call.method' to 2)."
						));
					call.method <- 2;
					}
				}
			if (call.method == 2) { kd.values <- c(0.85, 0.95); }# put whatever ends up being the default in apply.kd.cna.thresh here!!

			cna.rounded[,tmr] <- NanoStringNormCNV::call.copy.number.state(
				input = normalized.data[use.genes, c(1:3, tmr.ind, ref.ind), drop = FALSE],
				reference = phenodata$ref.name[has.ref[tmr]],
				sex.info = sample.sex,
				per.chip = per.chip,
				chip.info = phenodata,
				thresh.method = 'KD',
				kd.vals = kd.values
				)[,4];
			}
		}

	colnames(cna.rounded) <- phenodata$SampleID[has.ref];
	colnames(cna.raw)     <- phenodata$SampleID[has.ref];

	# combine output
	cna.all <- list(rounded = cna.rounded, raw = cna.raw);

	return(cna.all);
	
	}
