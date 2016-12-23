collapse.duplicate.samples <- function(normalized.data, annotation) {
	# for each sample that is duplicated, merge the counts so you get a single count
	normalized.data.aggr <- normalized.data[, c('CodeClass', 'Name', 'Accession')];

	# loop over each sample
	for (sample.name in unique(annotation$Name)) {

		# get the colnames equivalent of sample.name
		sample.id <- annotation$SampleID[sample.name == annotation$Name];

		# aggregate this by median value
		normalized.data.aggr[, sample.name] <- apply(
			X = normalized.data[, sample.id, drop = FALSE],
			MARGIN = 1,
			FUN = median
			);
		}

	return(normalized.data.aggr);
	}

