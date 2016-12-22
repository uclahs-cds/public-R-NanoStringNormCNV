calculate.replicate.concordance <- function(cna.rounded.reps, phenodata.reps) {
	# set up output variable
	out.table <- as.data.frame(
		matrix(
			nrow = nrow(cna.rounded.reps),
			ncol = length(unique(phenodata.reps$Name)),
			dimnames = list(
				row.names(cna.rounded.reps),
				unique(phenodata.reps$Name)
				)
			)
		);

	# loop over each unique sample name and create a difference matrix
	for (this.sample in unique(phenodata.reps$Name)) {
		# get the samples pertaining to the name
		this.ID <- phenodata.reps[phenodata.reps$Name == this.sample, 'SampleID'];

		# calculate the per-gene CN concordance of replicates
		per.gene.conc <- apply(
			X = cna.rounded.reps[,which(colnames(cna.rounded.reps) == this.ID), drop = FALSE],
			MARGIN = 1,
			FUN = function(f) { ifelse(length(unique(as.numeric(f))) == 1, 1, 0); }
			);

		# add results to output table
		out.table[,this.sample] <- per.gene.conc;
		}

	# return output
	return (out.table);
	}

