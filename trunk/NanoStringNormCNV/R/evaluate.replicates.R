evaluate.replicates <- function(normalized.data, phenodata, cnas) {
	# get samples with replicates
	pheno.reps <- phenodata[which(phenodata$has.repl == 1),];

	nano.reps <- normalized.data[, colnames(normalized.data) %in% pheno.reps$SampleID];
	cna.reps  <- normalized.data[, colnames(cnas) %in% pheno.reps$SampleID];

	pheno.reps.nano <- pheno.reps[pheno.reps$SampleID %in% colnames(nano.reps),];
	pheno.reps.cnas <- pheno.reps[pheno.reps$SampleID %in% colnames(cna.reps),];

	# checks
	if (!(ncol(nano.reps) > 0)) {
		stop("Cannot proceed with evaluation as no replicates were identified in normalized count data!");
		}
	if (!(ncol(cna.reps) > 0)) {
		stop("Cannot proceed with evaluation as no replicates were identified in CNA call data!");
		}

	# order according to samples
	pheno.reps.nano <- pheno.reps.nano[order(pheno.reps.nano$SampleID),];
	nano.reps <- nano.reps[, pheno.reps.nano$SampleID];

	pheno.reps.cnas <- pheno.reps.cnas[order(pheno.reps.cnas$SampleID),];
	cna.reps <- cna.reps[, pheno.reps.cnas$SampleID];

	# check count variance
	var.matrix <- NanoStringNormCNV::calculate.replicate.variance(
		norm.data.reps = nano.reps,
		phenodata.reps = pheno.reps.nano
		);

	# check CNA concordance
	conc.matrix <- NanoStringNormCNV::calculate.replicate.concordance(
		cnas.reps = cna.reps,
		phenodata.reps = pheno.reps.cnas
		);
	conc.summary <- apply(conc.matrix, 2, sum)/nrow(conc.matrix);

	return(
		list(
			variance = var.matrix,
			concordance = conc.matrix,
			conc.summary = conc.summary,
			norm.data = nano.reps,
			norm.pheno = pheno.reps.nano,
			cna.data = cna.reps,
			cna.pheno = pheno.reps.cnas
			)
		);
	}
