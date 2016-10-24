calculate.replicate.variance <- function (norm.data.reps, phenodata.reps, var.function = 'sd') {

	# initiate an out.table file
	out.table <- as.data.frame(
		matrix(
			nrow = nrow(norm.data.reps),
			ncol = length(unique(phenodata.reps$Name)),
			dimnames = list(
				row.names(norm.data.reps),
				unique(phenodata.reps$Name)
				)
			)
		);

	# loop over each unique sample from Name and create a difference matrix
	for (this.sample in unique(phenodata.reps$Name)) {
		# get the two samples pertaining to the name
		this.ID <- phenodata.reps[phenodata.reps$Name == this.sample, 'SampleID'];

		# calculate the per-gene variance of replicates
		per.gene.var <- apply(
			X = norm.data.reps[,which(colnames(norm.data.reps) == this.ID), drop = FALSE],
			MARGIN = 1,
			FUN = var.function
			);

		# add that to the out.table
		out.table[,this.sample] <- per.gene.var;
		}

	# return out.table
	return (out.table);
	}
