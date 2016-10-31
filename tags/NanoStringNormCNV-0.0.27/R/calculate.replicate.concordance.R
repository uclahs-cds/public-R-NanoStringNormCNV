calculate.replicate.concordance <- function (cnas.reps, phenodata.reps) {

	# initiate an out.table file
	out.table <- as.data.frame(
		matrix(
			nrow = nrow(cnas.reps),
			ncol = length(unique(phenodata.reps$Name)),
			dimnames = list(
				row.names(cnas.reps),
				unique(phenodata.reps$Name)
				)
			)
		);

	# loop over each unique sample from Name and create a difference matrix
	for (this.sample in unique(phenodata.reps$Name)) {
		# get the samples pertaining to the name
		this.ID <- phenodata.reps[phenodata.reps$Name == this.sample, 'SampleID'];

		# calculate the per-gene CN concordance of replicates
		per.gene.conc <- apply(
			X = cnas.reps[,which(colnames(cnas.reps) == this.ID), drop = FALSE],
			MARGIN = 1,
			FUN = function(f) { ifelse(length(unique(as.numeric(f))) == 1, 1, 0); }
			);

		# add that to the out.table
		out.table[,this.sample] <- per.gene.conc;
		}

	# return out.table
	return (out.table);
	}

