call.cnas.with.matched.normals <- function(normalized.data, phenodata, per.chip = FALSE, kd.option = 0) {
	flog.warn("Currently, cannot call on chromosomes X and Y!");

	# use non-control probes (from autosomes only)
	use.genes <- which(normalized.data$CodeClass %in% c("Endogenous", "Housekeeping", "Invariant"));
	use.genes <- use.genes[!(use.genes %in% grep("chr[XY]", normalized.data$Name))];

	has.ref 	<- which(phenodata$ref.name != 'missing' & phenodata$type == 'Tumour');
	cna.raw 	<- matrix(nrow = length(use.genes), ncol = length(has.ref));
	cna.rounded <- matrix(nrow = length(use.genes), ncol = length(has.ref));

	# iterate through each sample here
	for (tmr in 1:length(has.ref)) {
		tmr.ind <- which(colnames(normalized.data) == phenodata$SampleID[has.ref[tmr]]);
		ref.ind <- which(colnames(normalized.data) == phenodata$ref.name[has.ref[tmr]]);

		cna.raw[,tmr] <- call.copy.number.state(
			input = normalized.data[use.genes, c(1:3, tmr.ind, ref.ind), drop = FALSE],
			reference = phenodata$ref.name[has.ref[tmr]],
			thresh.method = 'none',
			multi.factor = 2
			)[,4];

		if (kd.option <= 1) {
			cna.rounded[,tmr] <- NanoStringNormCNV::call.copy.number.state(
				input = normalized.data[use.genes, c(1:3, tmr.ind, ref.ind), drop = FALSE],
				reference = phenodata$ref.name[has.ref[tmr]],
				per.chip = per.chip,
				chip.info = phenodata
				)[,4];
		} else {
			cna.rounded[,tmr] <- NanoStringNormCNV::call.copy.number.state(
				input = normalized.data[use.genes, c(1:3, tmr.ind, ref.ind), drop = FALSE],
				reference = phenodata$ref.name[has.ref[tmr]],
				per.chip = per.chip,
				chip.info = phenodata,
				thresh.method = 'KD'
				)[,4];
			}
		}

	colnames(cna.rounded) <- phenodata$SampleID[has.ref];
	colnames(cna.raw)     <- phenodata$SampleID[has.ref];

	cna.all <- list(rounded = cna.rounded, raw = cna.raw);

	return(cna.all);
	}
