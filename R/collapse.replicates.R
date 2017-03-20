collapse.replicates <- function(normalized.data, phenodata) {
	normalized.data.aggr <- normalized.data[, c('CodeClass', 'Name', 'Accession')];

	# for each sample name, merge the replicate counts to get a single value per probe per sample
	for (sample.name in unique(phenodata$Name)) {

		sample.id <- phenodata$SampleID[sample.name == phenodata$Name];

		# aggregate by median value
		normalized.data.aggr[, sample.name] <- apply(
			X = normalized.data[, sample.id, drop = FALSE],
			MARGIN = 1,
			FUN = median
			);
		}

	return(normalized.data.aggr);
	}

