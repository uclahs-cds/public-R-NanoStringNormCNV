evaluate.replicates <- function(phenodata, normalized.data = NULL, cna.rounded = NULL) {
	# get samples with replicates
	pheno.reps <- phenodata[which(phenodata$HasReplicate == 1),];

	count.reps <- normalized.data[, colnames(normalized.data) %in% pheno.reps$SampleID];
	cna.reps   <- cna.rounded[, colnames(cna.rounded) %in% pheno.reps$SampleID];

	pheno.reps.count <- pheno.reps[pheno.reps$SampleID %in% colnames(count.reps),];
	pheno.reps.cnas  <- pheno.reps[pheno.reps$SampleID %in% colnames(cna.reps),];

	# check for replicates; if none, set output to NULL
	if (!is.null(count.reps) && !(ncol(count.reps) > 0)) {
		flog.warn("Cannot evaluate using normalized count data as no replicates were identified in dataset!");
		var.matrix <- count.reps <- pheno.reps.count <- NULL;
	} else if (is.null(count.reps)) {
		var.matrix <- count.reps <- pheno.reps.count <- NULL;
		}

	if (!is.null(cna.reps) && !(ncol(cna.reps) > 0)) {
		flog.warn("Cannot evaluate using CNA data as no replicates were identified in dataset!");
		conc.matrix <- conc.summary <- cna.reps <- pheno.reps.cnas <- NULL;
	} else if (is.null(cna.reps)) {
		conc.matrix <- conc.summary <- cna.reps <- pheno.reps.cnas <- NULL;
		}

	# order samples
	if (!is.null(count.reps)) {
		pheno.reps.count <- pheno.reps.count[order(pheno.reps.count$SampleID),];
		count.reps <- count.reps[, pheno.reps.count$SampleID];
		}

	if (!is.null(cna.reps)) {
		pheno.reps.cnas <- pheno.reps.cnas[order(pheno.reps.cnas$SampleID),];
		cna.reps <- cna.reps[, pheno.reps.cnas$SampleID];
		}

	# calculate count variance
	if (!is.null(count.reps)) {
		var.matrix <- NanoStringNormCNV::calculate.replicate.variance(
			norm.data.reps = count.reps,
			phenodata.reps = pheno.reps.count
			);
		}

	# calculate CNA concordance
	if (!is.null(cna.reps)) {
		conc.matrix <- NanoStringNormCNV::calculate.replicate.concordance(
			cnas.reps = cna.reps,
			phenodata.reps = pheno.reps.cnas
			);
		conc.summary <- apply(conc.matrix, 2, sum)/nrow(conc.matrix);
		}
		
	return(
		list(
			variance = var.matrix,
			concordance = conc.matrix,
			conc.summary = conc.summary,
			norm.counts = count.reps,
			count.pheno = pheno.reps.count,
			cna.calls = cna.reps,
			cna.pheno = pheno.reps.cnas
			)
		);
	}
