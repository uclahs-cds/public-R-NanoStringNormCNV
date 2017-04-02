score.runs2 <- function(normalized.data, cna.rounded, phenodata, cna.normals = NULL) {
	# remove unnecessary probes
	normalized.data <- normalized.data[normalized.data$CodeClass %in% c("Endogenous", "Housekeeping", "Invariant"),];

	# log normalized counts, if possible
	if (all(unlist(normalized.data) >= 0)) {
		normalized.data <- log10(normalized.data[, -c(1:3)] + 1);
	} else {
		normalized.data <- normalized.data[, -c(1:3)];
		}

	# set up data
	phenodata.alu1 <- phenodata[phenodata$Fragmentation == 'Alu1',];
	phenodata.soni <- phenodata[phenodata$Fragmentation == 'Sonicate',];

	norm.data.alu1 <- normalized.data[, phenodata.alu1$SampleID];
	norm.data.soni <- normalized.data[, phenodata.soni$SampleID];

	if (is.null(cna.normals)) {
		cnas <- cna.rounded;
	} else {
		cnas <- cbind(cna.rounded, cna.normals);
		}

	cnas.alu1 <- cnas[, colnames(cnas) %in% phenodata.alu1$SampleID];
	cnas.soni <- cnas[, colnames(cnas) %in% phenodata.soni$SampleID];

	# order everything
	phenodata.alu1 <- phenodata.alu1[order(phenodata.alu1$SampleID),];
	norm.data.alu1 <- norm.data.alu1[, order(colnames(norm.data.alu1))];
	cnas.alu1 <- cnas.alu1[, order(colnames(cnas.alu1))];

	phenodata.soni <- phenodata.soni[order(phenodata.soni$SampleID),];
	norm.data.soni <- norm.data.soni[, order(colnames(norm.data.soni))];
	cnas.soni <- cnas.soni[, order(colnames(cnas.soni))];

	# initialize output variables
	scores <- list(
		counts.chip = NA,
		cnas.chip = NA,
		counts.type = NA,
		cnas.type = NA
		);
	scores.alu1 <- scores;
	scores.soni <- scores;

	### Check the adjusted rand index (ARI) of multiple parameters
	## using normalized counts
	# evaluate using cartridge information
	scores.alu1$counts.chip <- NanoStringNormCNV::get.ari(
		data.to.cluster = norm.data.alu1,
		feature = phenodata.alu1$Cartridge,
		is.discrete = FALSE
		);
	scores.soni$counts.chip <- NanoStringNormCNV::get.ari(
		data.to.cluster = norm.data.soni,
		feature = phenodata.soni$Cartridge,
		is.discrete = FALSE
		);

	# evaluate using tissue type information
	scores.alu1$counts.type <- NanoStringNormCNV::get.ari(
		data.to.cluster = norm.data.alu1,
		feature = phenodata.alu1$Type,
		is.discrete = FALSE
		);
	scores.soni$counts.type <- NanoStringNormCNV::get.ari(
		data.to.cluster = norm.data.soni,
		feature = phenodata.soni$Type,
		is.discrete = FALSE
		);

	## using copy number calls
	# evaluate using cartridge information
	scores.alu1$cnas.chip <- NanoStringNormCNV::get.ari(
		data.to.cluster = cnas.alu1,
		feature = phenodata.alu1$Cartridge,
		is.discrete = FALSE
		);
	scores.soni$cnas.chip <- NanoStringNormCNV::get.ari(
		data.to.cluster = cnas.soni,
		feature = phenodata.soni$Cartridge,
		is.discrete = FALSE
		);

	# evaluate using tissue type information
	if (!is.null(cna.normals)) {
		scores.alu1$cnas.type <- NanoStringNormCNV::get.ari(
			data.to.cluster = cnas.alu1,
			feature = phenodata.alu1$Type,
			is.discrete = FALSE
			);
		scores.soni$cnas.type <- NanoStringNormCNV::get.ari(
			data.to.cluster = cnas.soni,
			feature = phenodata.soni$Type,
			is.discrete = FALSE
			);
		}

	return(list(scores.alu1 = scores.alu1, scores.soni = scores.soni));
	}
