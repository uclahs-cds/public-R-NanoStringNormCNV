evaluate.replicates <- function(phenodata, normalized.data = NULL, cna.rounded = NULL) {
	# set up output variables
	var.matrix <- conc.matrix <- conc.summary <- NULL;

	# extract samples with replicates
	pheno.reps <- phenodata[which(phenodata$HasReplicate == 1),];

	count.reps <- normalized.data[, colnames(normalized.data) %in% pheno.reps$SampleID];
	cna.reps   <- cna.rounded[, colnames(cna.rounded) %in% pheno.reps$SampleID];

	pheno.reps.count <- pheno.reps[pheno.reps$SampleID %in% colnames(count.reps),];
	pheno.reps.cna   <- pheno.reps[pheno.reps$SampleID %in% colnames(cna.reps),];

	# check for replicates
	if (!(nrow(pheno.reps.count) > 0)) {
		if (!is.null(normalized.data)) {
			flog.warn("Cannot evaluate replicates using normalized count data: no replicates identified!");
			}
		count.reps <- pheno.reps.count <- NULL;
		}

	if (!(nrow(pheno.reps.cna) > 0)) {
		if (!is.null(cna.rounded)) {
			flog.warn("Cannot evaluate replicates using CNA data: no replicates identified!");
			}
		cna.reps <- pheno.reps.cna <- NULL;
		}

	# order samples
	if (!is.null(count.reps)) {
		pheno.reps.count <- pheno.reps.count[order(pheno.reps.count$SampleID),];
		count.reps <- count.reps[, pheno.reps.count$SampleID];
		}

	if (!is.null(cna.reps)) {
		pheno.reps.cna <- pheno.reps.cna[order(pheno.reps.cna$SampleID),];
		cna.reps <- cna.reps[, pheno.reps.cna$SampleID];
		}

	# calculate count variance
	if (!is.null(count.reps)) {
		var.matrix <- NanoStringNormCNV::calculate.replicate.variance(
			normalized.data.reps = count.reps,
			phenodata.reps = pheno.reps.count
			);
		}

	# calculate CNA concordance
	if (!is.null(cna.reps)) {
		conc.matrix <- NanoStringNormCNV::calculate.replicate.concordance(
			cna.rounded.reps = cna.reps,
			phenodata.reps = pheno.reps.cna
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
			cna.pheno = pheno.reps.cna
			)
		);
	}
