score.runs2 <- function(normalized.data, cna.rounded, phenodata, cna.normals = NULL) {
	# remove unnecessary probes
	probe.class <- c("Endogenous", "Housekeeping", "Invariant");
	normalized.data <- normalized.data[normalized.data$CodeClass %in% probe.class,];

	# log normalized counts, if possible
	if (all(unlist(normalized.data) >= 0)) {
		normalized.data <- log10(normalized.data[, -c(1:3)] + 1);
	} else {
		normalized.data <- normalized.data[, -c(1:3)];
		}

	# set up data
	if (is.null(cna.normals))  cnas <- cna.rounded;
	if (!is.null(cna.normals)) cnas <- cbind(cna.rounded, cna.normals);
	
	pheno.alui <- phenodata[phenodata$Fragmentation == 'Alui',];
	pheno.soni <- phenodata[phenodata$Fragmentation == 'Sonicate',];

	norm.alui <- normalized.data[, pheno.alui$SampleID, drop = FALSE];
	norm.soni <- normalized.data[, pheno.soni$SampleID, drop = FALSE];
	cnas.alui <- cnas[, colnames(cnas) %in% pheno.alui$SampleID, drop = FALSE];
	cnas.soni <- cnas[, colnames(cnas) %in% pheno.soni$SampleID, drop = FALSE];

	pheno.alui.norm <- pheno.alui[(pheno.alui$SampleID %in% names(norm.alui)),];
	pheno.soni.norm <- pheno.soni[(pheno.soni$SampleID %in% names(norm.soni)),];
	pheno.alui.cnas <- pheno.alui[(pheno.alui$SampleID %in% colnames(cnas.alui)),];
	pheno.soni.cnas <- pheno.soni[(pheno.soni$SampleID %in% colnames(cnas.soni)),];

	# order everything
	pheno.alui.norm <- pheno.alui.norm[match(pheno.alui.norm$SampleID, names(norm.alui)),];
	pheno.soni.norm <- pheno.soni.norm[match(pheno.soni.norm$SampleID, names(norm.soni)),];
	pheno.alui.cnas <- pheno.alui.cnas[match(pheno.alui.cnas$SampleID, colnames(cnas.alui)),];
	pheno.soni.cnas <- pheno.soni.cnas[match(pheno.soni.cnas$SampleID, colnames(cnas.soni)),];

	# initialize output variables
	scores <- list(
		counts.chip = NA,
		cnas.chip = NA,
		counts.type = NA,
		cnas.type = NA
		);
	scores.alui <- scores;
	scores.soni <- scores;
	scores.frag <- list(counts = NA, cnas = NA);

	### Check the adjusted rand index (ARI) of multiple parameters
	## using normalized counts
	# evaluate using cartridge information
	scores.alui$counts.chip <- NanoStringNormCNV::get.ari(
		data.to.cluster = norm.alui,
		feature = pheno.alui.norm$Cartridge,
		is.discrete = FALSE
		);
	scores.soni$counts.chip <- NanoStringNormCNV::get.ari(
		data.to.cluster = norm.soni,
		feature = pheno.soni.norm$Cartridge,
		is.discrete = FALSE
		);

	# evaluate using tissue type information
	scores.alui$counts.type <- NanoStringNormCNV::get.ari(
		data.to.cluster = norm.alui,
		feature = pheno.alui.norm$Type,
		is.discrete = FALSE
		);
	scores.soni$counts.type <- NanoStringNormCNV::get.ari(
		data.to.cluster = norm.soni,
		feature = pheno.soni.norm$Type,
		is.discrete = FALSE
		);

	## using copy number calls
	# evaluate using cartridge information
	scores.alui$cnas.chip <- NanoStringNormCNV::get.ari(
		data.to.cluster = cnas.alui,
		feature = pheno.alui.cnas$Cartridge,
		is.discrete = TRUE
		);
	scores.soni$cnas.chip <- NanoStringNormCNV::get.ari(
		data.to.cluster = cnas.soni,
		feature = pheno.soni.cnas$Cartridge,
		is.discrete = TRUE
		);

	# evaluate using tissue type information
	if (!is.null(cna.normals)) {
		scores.alui$cnas.type <- NanoStringNormCNV::get.ari(
			data.to.cluster = cnas.alui,
			feature = pheno.alui.cnas$Type,
			is.discrete = TRUE
			);
		scores.soni$cnas.type <- NanoStringNormCNV::get.ari(
			data.to.cluster = cnas.soni,
			feature = pheno.soni.cnas$Type,
			is.discrete = TRUE
			);
		}

	# fragmentation method!
	norm.frag <- normalized.data[, order(colnames(normalized.data)), drop = FALSE];
	pheno.frag.norm <- phenodata[phenodata$SampleID %in% colnames(norm.frag),];
	pheno.frag.norm <- pheno.frag.norm[order(pheno.frag.norm$SampleID),];
	scores.frag$counts <- NanoStringNormCNV::get.ari(
		data.to.cluster = norm.frag,
		feature = pheno.frag.norm$Fragmentation,
		is.discrete = FALSE
		);

	cnas.frag <- cnas[, order(colnames(cnas)), drop = FALSE];
	pheno.frag.cnas <- phenodata[phenodata$SampleID %in% colnames(cnas.frag),];
	pheno.frag.cnas <- pheno.frag.cnas[order(pheno.frag.cnas$SampleID),];
	scores.frag$cnas <- NanoStringNormCNV::get.ari(
		data.to.cluster = cnas.frag,
		feature = pheno.frag.cnas$Fragmentation,
		is.discrete = TRUE
		);

	return(list(scores.alui = scores.alui, scores.soni = scores.soni, scores.frag = scores.frag));
	}
