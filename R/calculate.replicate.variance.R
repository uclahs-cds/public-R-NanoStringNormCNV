calculate.replicate.variance <- function (normalized.data.reps, phenodata.reps, var.function = 'sd') {
	# set up output variable
	out.table <- as.data.frame(
		matrix(
			nrow = nrow(normalized.data.reps),
			ncol = length(unique(phenodata.reps$Name)),
			dimnames = list(
				row.names(normalized.data.reps),
				unique(phenodata.reps$Name)
				)
			)
		);

	# loop over each unique sample name and create a difference matrix
	for (this.sample in unique(phenodata.reps$Name)) {
		# get the samples pertaining to the name
		this.ID <- phenodata.reps[phenodata.reps$Name == this.sample, 'SampleID'];

		# calculate the per-gene variance of replicates
		per.gene.var <- apply(
			X = normalized.data.reps[,which(colnames(normalized.data.reps) %in% this.ID), drop = FALSE],
			MARGIN = 1,
			FUN = var.function
			);

		# add results to output table
		out.table[,this.sample] <- per.gene.var;
		}

	# return output
	return (out.table);
	}
