

calculate.replicate.variance <- function (input, pheno, gene.names, var.function = 'sd') {

	# initiate an out.table file
	out.table <- as.data.frame(
		matrix(
			nrow = nrow(input),
			ncol = length(unique(pheno$Name)),
			dimnames = list(
				gene.names,
				unique(pheno$Name)
				)
			)
		);
	print(dim(out.table));

	# loop over each unique sample from Name and create a difference matrix
	for (this.sample in unique(pheno$Name)) {
		# get the two samples pertaining to the name
		this.ID <- pheno[pheno$Name == this.sample, 'SampleID'];

		# calculate the per-gene variance of replicates
		per.gene.var <- apply(
			X = input[,which(colnames(input) == this.ID), drop = FALSE],
			MARGIN = 1,
			FUN = var.function
			);

		# add that to the out.table
		out.table[,this.sample] <- per.gene.var;
		}

	# return out.table
	return (out.table);
	}
